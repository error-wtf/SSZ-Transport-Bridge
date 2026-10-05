"""SSZ closure backend — READ-ONLY adapter to SSZ_FULL_CLOSURE.

CORE RULE: this adapter contains ZERO SSZ physics.  Every quantity is
computed by the SSZ_FULL_CLOSURE transport core (postclosure/transport.py:
k^nu nabla_nu k^mu = 0, null congruence, eikonal phase S_r, redshift);
the adapter only wraps results into the bridge contract types.

Requires the closure checkout; raises ImportError (skip-able) when it
is unavailable, so CI without the checkout skips cleanly.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    TransportArchitecture,
    TransportOperator,
)

_CANDIDATES = [
    Path("/home/error/physics/clones/SSZ_FULL_CLOSURE/src"),
    Path("/home/error/ssz-full-closure/src"),
]
_REF = Path(os.environ.get(
    "SSZ_FULL_CLOSURE_SRC",
    str(next((c for c in _CANDIDATES if (c / "ssz_p5").exists()),
             _CANDIDATES[0]))))

import sys  # noqa: E402

if str(_REF) not in sys.path:
    sys.path.insert(0, str(_REF))

try:
    from ssz_p5.postclosure.transport import (
        eikonal_phase_and_redshift,
        load_member_metric,
        null_geodesic_transport,
    )
except ModuleNotFoundError as _exc:  # pragma: no cover - env-dependent
    # ImportError so pytest.importorskip can skip cleanly on CI.
    raise ImportError(
        f"ssz_closure backend requires the SSZ_FULL_CLOSURE checkout "
        f"(looked for {_REF}); set SSZ_FULL_CLOSURE_SRC.  "
        f"Original error: {_exc}") from _exc

ARCHITECTURE = TransportArchitecture(
    name="ssz_postclosure_transport",
    direction_parameter="inward",
    operator_class="null_geodesic",
)

OBSERVABLES = ObservableSet(
    directional=(),          # static spherical geometry: no sign-swap law
    scalar=("shapiro_dt_per_E", "reduced_phase_per_E",
            "redshift_a_to_b", "max_dE", "max_dL", "max_norm_res"),
    inverse="",              # inverse problem deferred (contract 5)
)


def _metric():
    return load_member_metric(Path(_REF).parent)


def make_operator(law: SignalLaw,
                  r_a: float = 1.5, r_b: float = 1.6) -> TransportOperator:
    """Defaults inside the locked member domain (1.409-1.638 r/r_s)."""
    """law.parameters: {"inward": +1/-1 or 0, "b": impact parameter}.

    The closure's own code does all physics: null geodesic transport
    (k^nu nabla_nu k^mu = 0 conservation check), eikonal phase and
    redshift quadrature.  The adapter only wraps.
    """
    import numpy as np  # local: only for finiteness checks

    def _solve(law_: SignalLaw) -> dict[str, Any]:
        m = _metric()
        inward = bool(law_.parameters.get("inward", 1.0) >= 0)
        b = float(law_.parameters.get("b", 2.5))
        r0 = float(law_.parameters.get("r0", r_a))
        # Domain exit is a NORMAL end for an open ray; the closure code
        # raises on it.  On domain exit we mark the geodesic and keep the
        # eikonal (exact quadrature) observables — all finite.
        tc = None
        try:
            tc = null_geodesic_transport(m, b=b, r0=r0, inward=inward)
        except RuntimeError as exc:
            if "termination event" not in str(exc):
                raise
        eik = eikonal_phase_and_redshift(m, r_a, r_b)
        out = {
            "shapiro_dt_per_E": eik["shapiro_dt_per_E"],
            "reduced_phase_per_E": eik["reduced_phase_per_E"],
            "redshift_a_to_b": eik["redshift_a_to_b"],
        }
        if tc is not None:
            out["max_dE"] = tc.max_dE
            out["max_dL"] = tc.max_dL
            out["max_norm_res"] = tc.max_norm_res
        else:
            out["max_dE"] = float("nan")
            out["max_dL"] = float("nan")
            out["max_norm_res"] = float("nan")
            out["geodesic_domain_exit"] = True
        if not all(np.isfinite(float(v)) for k, v in out.items()
                   if k != "geodesic_domain_exit"):
            raise ValueError("non-finite observable from closure transport")
        return out

    return TransportOperator(ARCHITECTURE, law, _solve)


def provenance() -> dict[str, str]:
    """Bind to the closure state: git commit + checkout path."""
    root = Path(_REF).parent
    head = "UNKNOWN"
    head_file = root / ".git/HEAD"
    if head_file.exists():
        txt = head_file.read_text().strip()
        if txt.startswith("ref: "):
            ref = root / ".git" / txt[5:]
            head = ref.read_text().strip()[:12] if ref.exists() else "UNKNOWN"
        else:
            head = txt[:12]
    return {"closure_git_head": head, "closure_src": str(_REF)}
