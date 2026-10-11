"""Phase & Observable Contract V1 -- backend-agnostic phase decomposition.

SCIENTIFIC CONTEXT (Integration order, Lino 2026-10-11):
Full eikonal phase vs geometric time share vs extra phase residue are
DISTINCT quantities.  The decomposition

    phi_full(omega) = omega * t_geom + phi_extra(omega)

is an IDENTITY enforced by construction (phi_extra is DEFINED as the
residue), so double counting is structurally excluded: a consumer must
take phi_full OR omega*t_geom + phi_extra -- never both.  The identity
gate `check_phase_decomposition` re-verifies it numerically per backend.

Units convention (V1): omega in rad per unit time, t in the backend's
time unit; phi is dimensionless (radians).  Every phase-producing
backend MUST declare its time unit and expose the three quantities
plus parity data.

Readout honesty (two-channel test, see tests/test_phase_two_channel.py):
a COMMON phase factor e^{i phi(omega)} leaves any auto-PSD invariant,
|e^{i phi} A|^2 = |A|^2.  Relative phase between channels appears only
in the cross-spectrum S_ab ~ e^{i(phi_a - phi_b)} H_a H_b*.  No
astrophysical detection claims are made anywhere in this module.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .contracts import SignalLaw

# Registered time-unit tokens; fail-closed against silent unit mixing.
TIME_UNITS = ("s", "geometric_m", "dimensionless")


@dataclass(frozen=True)
class PhaseDecomposition:
    """One direction's phase decomposition at fixed omega.

    Fields:
      omega: angular frequency (rad / time_unit), finite, >= 0.
      t_geom: geometric reference time (no effect share), time_unit.
      phi_full: full eikonal phase [rad] for this direction.
      phi_extra: DEFINED residue phi_full - omega*t_geom [rad].
      time_unit: one of TIME_UNITS.
      direction: +1 / -1 direction tag.
    """
    omega: float
    t_geom: float
    phi_full: float
    direction: int
    time_unit: str = "dimensionless"
    # computed, not user-supplied:
    phi_extra: float = field(init=False)

    def __post_init__(self) -> None:
        if self.time_unit not in TIME_UNITS:
            raise ValueError(
                f"unknown time_unit {self.time_unit!r}; must be one of {TIME_UNITS}")
        vals = (self.omega, self.t_geom, self.phi_full)
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("omega, t_geom, phi_full must be finite")
        if self.omega < 0:
            raise ValueError("omega must be >= 0 (rad/time_unit)")
        if self.t_geom < 0:
            raise ValueError("t_geom must be >= 0")
        if self.direction not in (+1, -1):
            raise ValueError("direction must be +1 or -1")
        object.__setattr__(self, "phi_extra",
                           self.phi_full - self.omega * self.t_geom)

    @property
    def identity_residual(self) -> float:
        """|phi_full - (omega*t_geom + phi_extra)| -- must be ~0 (identity)."""
        return abs(self.phi_full
                   - (self.omega * self.t_geom + self.phi_extra))


@dataclass(frozen=True)
class PhaseContractSpec:
    """Backend-declared phase surface.

    decompose(law, omega, direction) -> PhaseDecomposition.
    parity_group: observable name group that must flip under direction
    reversal (checked with the bridge's own rule engine, not re-implemented).
    """
    backend_name: str
    time_unit: str
    decompose: Callable[[Mapping, float, int], PhaseDecomposition]
    required_parameters: tuple[str, ...] = ()
    parity_observables: tuple[str, ...] = ("phi_full", "phi_extra")
    direction_parameter: str | None = None

    def __post_init__(self) -> None:
        if self.time_unit not in TIME_UNITS:
            raise ValueError(f"unknown time_unit {self.time_unit!r}")
        if not callable(self.decompose):
            raise ValueError("decompose must be callable")


def check_phase_decomposition(
    spec: PhaseContractSpec,
    law: Mapping,
    omega: float,
    direction: int = +1,
    resid_tol: float = 1e-12,
) -> dict[str, Any]:
    """Gate 1 (double counting): the decomposition identity holds, and the
    direction-reversal SEMANTICS hold.

    Parity semantics (important, physics not convention): reversing the
    propagation direction at fixed law (v -> -v for Sagnac, reversed ray
    for SSZ) gives t_geom invariant, |phi_full| positive in both
    directions (a directed time never changes sign), and the DIFFERENCE
    phase phi_full(+) - phi_full(-) is ODD under the law's
    direction_parameter flip.  A signed phase that flips with the ray
    direction is exposed separately via the backend's directional
    observables (delta_tau_signed etc.), not by relabeling phi_full.
    """
    dec = spec.decompose(law, omega, direction)
    dec_m = spec.decompose(law, omega, -direction)
    if not isinstance(dec, PhaseDecomposition):
        return {"identity_ok": False, "reason": "decompose did not return PhaseDecomposition"}
    identity_ok = dec.identity_residual <= resid_tol
    t_geom_inv = math.isclose(dec.t_geom, dec_m.t_geom,
                              rel_tol=0.0, abs_tol=resid_tol)
    positive_directed = dec.phi_full >= -resid_tol and dec_m.phi_full >= -resid_tol
    # difference phase odd under direction_parameter flip: for Sagnac this
    # is checked by flipping the direction_parameter (v -> -v) and requiring
    # Phi_diff -> -Phi_diff; generic backends expose direction_parameter via
    # the law, so we test with the flipped parameter when it exists.
    diff = dec.phi_full - dec_m.phi_full
    flipped = dict(law) if not isinstance(law, SignalLaw) else None
    diff_flipped = diff
    if flipped is not None and getattr(spec, "direction_parameter", None):
        key = spec.direction_parameter
        if key in flipped:
            flipped[key] = -flipped[key]
            dec_f = spec.decompose(flipped, omega, +1)
            dec_f_m = spec.decompose(flipped, omega, -1)
            diff_flipped = dec_f.phi_full - dec_f_m.phi_full
    diff_odd = math.isclose(diff, -diff_flipped, rel_tol=1e-9, abs_tol=resid_tol)
    return {
        "identity_ok": bool(identity_ok),
        "identity_residual": dec.identity_residual,
        "t_geom_direction_invariant": bool(t_geom_inv),
        "directed_phase_nonnegative": bool(positive_directed),
        "difference_phase_odd": bool(diff_odd),
        "pass": bool(identity_ok and t_geom_inv and positive_directed and diff_odd),
        "time_unit": dec.time_unit,
        "backend": spec.backend_name,
    }


def dimension_check(
    dec: PhaseDecomposition,
    omega_unit: str = "rad/time_unit",
) -> dict[str, Any]:
    """Gate 2 (units): phi dimensionless, omega*t_geom dimensionless in
    matching units.  V1 is a consistency gate on magnitudes and finiteness
    plus the declared unit token."""
    phi = dec.omega * dec.t_geom
    checks = {
        "phi_full_finite_dimensionless": math.isfinite(dec.phi_full),
        "phi_extra_finite": math.isfinite(dec.phi_extra),
        "omega_t_geom_finite": math.isfinite(phi),
        "time_unit_registered": dec.time_unit in TIME_UNITS,
        "omega_nonneg": dec.omega >= 0.0,
    }
    return {"checks": checks, "pass": all(checks.values()),
            "omega_unit": omega_unit}


# --- Two-channel readout physics (the measurement side) -------------------

def apply_common_phase(
    amp: np.ndarray, phase_fn: Callable[[np.ndarray], np.ndarray],
) -> np.ndarray:
    """A(omega) -> e^{i phi(omega)} A(omega)."""
    return np.exp(1j * np.asarray(phase_fn(amp))) * np.asarray(amp)


def auto_psd(a: np.ndarray) -> np.ndarray:
    """|A(omega)|^2 -- invariant under common phase by construction."""
    return np.abs(np.asarray(a)) ** 2


def cross_spectrum(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """S_ab = A * conj(B) -- carries the RELATIVE phase phi_a - phi_b."""
    return np.asarray(a) * np.conj(np.asarray(b))


def check_two_channel_readout(
    amps: tuple[np.ndarray, np.ndarray],
    phases: tuple[np.ndarray, np.ndarray],
    psd_atol: float = 1e-12,
    phase_gate_rad: float = 1e-9,
) -> dict[str, Any]:
    """Gate 3 (measurement honesty):
    (a) common phase leaves BOTH auto-PSDs invariant (allclose),
    (b) the cross-spectrum S_ab = e^{i(phi_a-phi_b)} A B* carries the
        RELATIVE phase ON TOP of the data term A B*: dividing out the
        data term must recover exactly e^{i(phi_a-phi_b)}.  This is the
        calibrated-data operation; skipping it would wrongly attribute
        the random data-term phase to the transport phase.
    """
    a, b = (np.asarray(x) for x in amps)
    ph_a, ph_b = (np.asarray(x) for x in phases)
    a_p = apply_common_phase(a, lambda _: ph_a)
    b_p = apply_common_phase(b, lambda _: ph_b)
    psd_before, psd_after = auto_psd(a), auto_psd(a_p)
    psd_b_before, psd_b_after = auto_psd(b), auto_psd(b_p)
    psd_invariant = bool(np.allclose(psd_before, psd_after, atol=psd_atol)
                         and np.allclose(psd_b_before, psd_b_after, atol=psd_atol))
    sab = cross_spectrum(a_p, b_p)
    data_term = cross_spectrum(a, b)
    rel = ph_a - ph_b
    recovered = sab / (data_term + 1e-300)
    expected_phase = np.exp(1j * rel)
    err = np.abs(recovered - expected_phase)
    cross_ok = bool(np.all(err < phase_gate_rad))
    return {
        "auto_psd_invariant": psd_invariant,
        "cross_spectrum_carries_relative_phase": cross_ok,
        "max_rel_phase_err": float(np.max(err)),
        "pass": bool(psd_invariant and cross_ok),
    }
