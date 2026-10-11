"""Sagnac backend -- READ-ONLY adapter to the validated reference repo.

This file contains NO Sagnac formulas.  It imports the validated
reference implementation (Sagnac-Reference-Transport) and exposes it
through the bridge contracts.  The reference repo stays the single
source of truth.
"""
from __future__ import annotations

import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    TransportArchitecture,
    TransportOperator,
)
from transport_bridge.phase_contract import PhaseContractSpec, PhaseDecomposition

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

try:
    from sagnac_reference.analytic import (
        SagnacConfig,
        delta_tau_detector,
        directed_times,
    )
    from sagnac_reference.inversion import velocity_from_times
    from sagnac_reference.true_chain import return_time_true_chain
except ModuleNotFoundError as _exc:  # pragma: no cover - env-dependent
    # ImportError (not RuntimeError) so pytest.importorskip can skip the
    # backend cleanly on machines without the reference checkout (CI).
    raise ImportError(
        f"sagnac backend requires the validated reference checkout "
        f"(looked for {_REF}); set SAGNAC_REFERENCE_SRC.  Original "
        f"error: {_exc}") from _exc

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
    """v from (t_plus, t_minus) -- the reference repo's inversion."""
    return velocity_from_times(forward["t_plus"], forward["t_minus"], c)


# --- Phase & Observable Contract V1 adapter --------------------------------
# The reference repo stays untouched: directed_times / delta_tau_detector
# are IMPORTED from it; this adapter only declares how its times map into
# the bridge's phase decomposition.
#
#   t_geom          = L/c          (v-independent reference time)
#   phi_full(d,ω)   = ω · L/(c - d·v)      (full directed eikonal phase)
#   phi_extra(d,ω) := phi_full - ω·t_geom  (defined residue; gamma share of the
#                     detector time enters via delta_tau_detector elsewhere)
#
# Double counting is excluded structurally: consumers take phi_full OR
# ω·t_geom + phi_extra (identity enforced by PhaseDecomposition).

def phase_spec() -> PhaseContractSpec:
    def decompose(law: Mapping, omega: float, direction: int) -> PhaseDecomposition:
        if direction not in (+1, -1):
            raise ValueError("direction must be +1 or -1")
        params = law.parameters if isinstance(law, SignalLaw) else law
        cfg = SagnacConfig(L=float(params["L"]), c=float(params["c"]),
                           v=float(params["v"]))
        d = float(directed_times(cfg)[0 if direction > 0 else 1])
        return PhaseDecomposition(
            omega=float(omega), t_geom=cfg.T0, phi_full=float(omega) * d,
            direction=direction, time_unit="dimensionless")

    return PhaseContractSpec(
        backend_name="sagnac_ring",
        time_unit="dimensionless",
        decompose=decompose,
        required_parameters=("L", "c", "v"),
        direction_parameter="v",
    )
