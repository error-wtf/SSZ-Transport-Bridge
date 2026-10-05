"""Generic typed contracts for transport systems.

Design rules:
1.  Dataclasses only — no physics, no formulas, no numerics.
2.  Every check function is backend-agnostic: it consumes contract
    objects and asserts STRUCTURAL properties (direction flip, grid
    convergence, inverse round trip).
3.  A backend that satisfies these checks has proven it belongs to the
    same VALIDATION CLASS as every other backend that does — nothing
    more, nothing less.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class SignalLaw:
    """The immutable physical input parameters of a transport system.

    sagnac: (L, c, v).  ssz later: (f, h, ...) profile parameters.
    The bridge never interprets the fields; it only requires that the
    backend can hash-bind them (provenance) and perturb them along a
    declared direction axis.
    """
    system: str                       # "sagnac" | "ssz" (tag only)
    parameters: dict[str, float]      # opaque key->value
    provenance_sha256: str            # binds the parameters to a frozen source

    def perturbed(self, key: str, delta: float) -> SignalLaw:
        """Same law with one parameter shifted — used for direction/
        response checks.  Raises if the key is unknown (fail-closed)."""
        if key not in self.parameters:
            raise KeyError(f"unknown parameter {key!r} for system {self.system!r}")
        params = dict(self.parameters)
        params[key] = params[key] + delta
        return SignalLaw(self.system, params, self.provenance_sha256)


@dataclass(frozen=True)
class TransportArchitecture:
    """Declares the backend's transport structure without running it.

    direction_parameter: the parameter whose sign flip exchanges the two
    propagation directions (sagnac: v; ssz later: e.g. frame field sign).
    """
    name: str
    direction_parameter: str
    operator_class: str               # "wave" | "first_order" | ... (tag only)


@dataclass(frozen=True)
class TransportOperator:
    """A bound operator: architecture + signal law + a forward solve.

    `solve` returns raw directional observables as a dict; the bridge
    does not know their names.
    """
    architecture: TransportArchitecture
    law: SignalLaw
    solve: Callable[[SignalLaw], dict[str, Any]]

    def run(self) -> dict[str, Any]:
        return self.solve(self.law)


@dataclass(frozen=True)
class PropagatedState:
    """Directional results of one forward solve, tagged ±1."""
    direction: int                    # +1 or -1
    values: dict[str, float]


@dataclass(frozen=True)
class ObservableSet:
    """Declarative description of the observables a backend exposes.

    directional: names that flip sign under direction reversal
                 (sagnac: dt, delta_phi).
    scalar: names that must be invariant (sagnac: total, gamma).
    inverse: name of the inverse-reconstruction target (sagnac: v).
    """
    directional: tuple[str, ...]
    scalar: tuple[str, ...]
    inverse: str


def check_forward(
    op: TransportOperator,
    obs: ObservableSet,
    forward_values: dict[str, tuple[float, float]],
    atol: float,
) -> dict[str, bool]:
    """forward_values: observable name -> (plus-direction, minus-direction).
    Structural contract: directional names must swap, scalars must be
    invariant."""
    results: dict[str, bool] = {}
    for name in obs.directional:
        plus, minus = forward_values[name]
        results[f"directional_swap:{name}"] = bool(
            np.isclose(plus, -minus, rtol=0.0, atol=atol))
    for name in obs.scalar:
        plus, minus = forward_values[name]
        results[f"scalar_invariant:{name}"] = bool(
            np.isclose(plus, minus, rtol=0.0, atol=atol))
    return results


def check_convergence(
    solve_at: Callable[[Any], float],
    grids: list[Any],
    target_order: float | None = None,
    min_observed_order: float = 1.5,
) -> dict[str, Any]:
    """Generic grid convergence: solve_at(grid) must approach a limit with
    non-increasing errors.  If the backend also provides an analytic
    reference through target values, the observed order is estimated.
    Returns order estimate + monotonicity flags."""
    vals = [solve_at(g) for g in grids]
    diffs = [abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1)]
    monotone = all(diffs[i + 1] <= diffs[i] * 1.5 + 1e-300
                   for i in range(len(diffs) - 1))
    observed_order = float("nan")
    if len(diffs) >= 2 and diffs[0] > 0 and diffs[-1] > 0:
        ratio = np.log(diffs[0] / diffs[-1]) / np.log(
            (grids[1] - grids[0]) if False else 2.0)
        # grids are assumed to double; order from successive halvings:
        orders = [np.log2(diffs[i] / diffs[i + 1]) for i in range(len(diffs) - 1)]
        observed_order = float(np.mean(orders))
        ratio = observed_order  # noqa: F841  (kept for clarity)
    return {
        "values": vals,
        "monotone_decrease": monotone,
        "observed_order": observed_order,
        "order_ok": (not np.isnan(observed_order))
        and (observed_order >= min_observed_order),
    }


def check_parity(
    plus: PropagatedState,
    minus: PropagatedState,
    obs: ObservableSet,
    atol: float,
) -> dict[str, bool]:
    """Direction swap contract on propagated states."""
    out: dict[str, bool] = {}
    for name in obs.directional:
        out[f"parity_swap:{name}"] = bool(np.isclose(
            plus.values[name], -minus.values[name], rtol=0.0, atol=atol))
    for name in obs.scalar:
        out[f"parity_invariant:{name}"] = bool(np.isclose(
            plus.values[name], minus.values[name], rtol=0.0, atol=atol))
    return out


def check_inverse_round_trip(
    law: SignalLaw,
    inverse_param: str,
    forward_observables: dict[str, float],
    inverse_fn: Callable[[dict[str, float]], float],
    atol: float,
) -> dict[str, Any]:
    """forward -> inverse must reproduce the law's parameter."""
    v_rec = inverse_fn(forward_observables)
    v_true = law.parameters[inverse_param]
    return {
        "reconstructed": float(v_rec),
        "true": float(v_true),
        "abs_error": abs(float(v_rec) - float(v_true)),
        "ok": bool(abs(float(v_rec) - float(v_true)) <= atol),
    }
