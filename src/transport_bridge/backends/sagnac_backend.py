"""Sagnac backend — READ-ONLY adapter to the validated reference repo.

This file contains NO Sagnac formulas.  It imports the validated
reference implementation (Sagnac-Reference-Transport) and exposes it
through the bridge contracts.  The reference repo stays the single
source of truth.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    TransportArchitecture,
    TransportOperator,
)

# The reference repo location may be provided via env; default is the
# sibling checkout in the physics workspace.
_CANDIDATES = [
    Path(__file__).resolve().parents[4] / "Sagnac-Reference-Transport" / "src",
    Path("/home/error/Sagnac-Reference-Transport/src"),
    Path(__file__).resolve().parents[3] / "Sagnac-Reference-Transport" / "src",
]
_REF = Path(__import__("os").environ.get(
    "SAGNAC_REFERENCE_SRC",
    str(next((c for c in _CANDIDATES if (c / "sagnac_reference").exists()),
             _CANDIDATES[0]))))
if str(_REF) not in sys.path:
    sys.path.insert(0, str(_REF))

from sagnac_reference.analytic import (  # noqa: E402  (reference repo)
    SagnacConfig,
    delta_tau_detector,
    directed_times,
)
from sagnac_reference.inversion import velocity_from_times  # noqa: E402
from sagnac_reference.true_chain import return_time_true_chain  # noqa: E402

ARCHITECTURE = TransportArchitecture(
    name="sagnac_ring",
    direction_parameter="v",
    operator_class="wave",
)

OBSERVABLES = ObservableSet(
    # delta_tau is SIGNED (arrival-order dependent): it flips with the
    # direction swap.  |delta_tau| is the invariant magnitude.
    directional=("dt", "delta_phi", "delta_tau_signed"),
    scalar=("total", "abs_delta_tau"),
    inverse="v",
)


def _cfg(law: SignalLaw) -> SagnacConfig:
    p = law.parameters
    return SagnacConfig(L=p["L"], c=p["c"], v=p["v"])


def solve(law: SignalLaw) -> dict[str, Any]:
    """Forward solve: analytic times + independent true-chain route."""
    cfg = _cfg(law)
    tp, tm = directed_times(cfg)
    dt = tp - tm
    total = tp + tm
    # the independent route (true chain) validates the analytic one:
    tp_chain = return_time_true_chain(cfg, +1, steps=32000).return_time
    assert abs(tp_chain - tp) < 1e-10, "true chain vs analytic mismatch"
    delta_tau = delta_tau_detector(cfg)
    return {
        "t_plus": tp,
        "t_minus": tm,
        "dt": dt,
        "total": total,
        "delta_tau": delta_tau,
        "delta_tau_signed": delta_tau,
        "abs_delta_tau": abs(delta_tau),
        "delta_phi_placeholder": dt,  # phase = omega * tau; omega is readout-side
    }


def make_operator(law: SignalLaw, omega_detector: float = 2.0) -> TransportOperator:
    def _solve(law_: SignalLaw) -> dict[str, Any]:
        out = solve(law_)
        out["delta_phi"] = omega_detector * out["delta_tau"]
        return out

    return TransportOperator(ARCHITECTURE, law, _solve)


def inverse_v(forward: dict[str, float], c: float = 1.0) -> float:
    """v from (t_plus, t_minus) — the reference repo's inversion."""
    return velocity_from_times(forward["t_plus"], forward["t_minus"], c)
