"""SSZ backend — SYNTHETIC INVERSION ONLY (for now).

Contract: this backend demonstrates the bridge's synthetic forward/
inverse loop against a KNOWN frozen geometry.  It deliberately does NOT
touch SSZ_FULL_CLOSURE runtime code — the geometry proxy here is a
transparent toy that shares only the SHAPE of the later problem:

    known frozen parameters -> synthetic observables -> recover them.

Real SSZ geometry (f, h profiles from the locked member) plugs into the
same contracts once the healthy-operator decision (T) is made.
"""
from __future__ import annotations

from typing import Any

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    TransportArchitecture,
    TransportOperator,
)

ARCHITECTURE = TransportArchitecture(
    name="ssz_geometry_proxy",
    direction_parameter="frame_sign",
    operator_class="first_order",
)

OBSERVABLES = ObservableSet(
    directional=("delta_phi", "delta_t"),
    scalar=("z_mean",),
    inverse="frame_sign",
)


def solve(law: SignalLaw) -> dict[str, Any]:
    """Synthetic transport through a Schwarzschild-like proxy.

    Parameters: r_s, r_obs, frame_sign (+1/-1).  Observables:
    coordinate-time deltas for pro/retgrade circular orbits scaled by
    the frame sign, mean redshift (scalar), delta_t (directional).

    The proxy formulas are TOY geometry (linearized Lense-Thirring-like
    frame dragging) — used ONLY to exercise the bridge's inverse loop.
    """
    p = law.parameters
    r_s, r_obs = p["r_s"], p["r_obs"]
    sign = p["frame_sign"]
    c = p.get("c", 1.0)
    # toy frame-dragging-like splitting, first order in r_s/r:
    splitting = 0.05 * (r_s / r_obs)  # magnitude, direction-independent
    t_plus = 2.0 * 3.141592653589793 * r_obs / c * (1.0 + sign * splitting)
    t_minus = 2.0 * 3.141592653589793 * r_obs / c * (1.0 - sign * splitting)
    delta_t = t_plus - t_minus
    z_mean = 0.5 * r_s / r_obs  # scalar, sign-independent
    return {
        "t_plus": t_plus,
        "t_minus": t_minus,
        "delta_t": delta_t,
        "delta_phi": sign * splitting * 10.0,  # flips with frame sign
        "z_mean": z_mean,
    }


def make_operator(law: SignalLaw) -> TransportOperator:
    return TransportOperator(ARCHITECTURE, law, solve)


def inverse_frame_sign(delta_phi: float, epsilon: float = 1e-9) -> float:
    """Recover the frame_sign from the directional phase observable."""
    if abs(delta_phi) < epsilon:
        raise ValueError("delta_phi too small to determine frame sign")
    return 1.0 if delta_phi > 0 else -1.0
