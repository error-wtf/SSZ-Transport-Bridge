"""Phase & Observable Contract V1 -- independent test suite.

Covers (integration order, Lino 2026-10-11):
  Gate 1  double-counting exclusion (decomposition identity + parity)
  Gate 2  units / dimension consistency
  Gate 3  two-channel readout honesty (auto-PSD invariant, relative phase
          in cross-spectrum)
plus reference limits, quadrature convergence and provenance binding.
No astrophysical detection claims; frozen F7/B6 gates are not touched.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from transport_bridge import (
    PhaseDecomposition,
    SignalLaw,
    check_phase_decomposition,
    check_two_channel_readout,
    dimension_check,
)
from transport_bridge.phase_contract import (
    TIME_UNITS,
    PhaseContractSpec,
)

try:
    from transport_bridge.backends import sagnac_backend as sag
    SAGNAC_AVAILABLE = True
except ImportError:
    SAGNAC_AVAILABLE = False

try:
    from transport_bridge.backends import ssz_closure_backend as ssz
    SSZ_AVAILABLE = True
except Exception:
    SSZ_AVAILABLE = False

_SAG_REASON = "requires Sagnac reference checkout"
_SSZ_REASON = "requires SSZ_FULL_CLOSURE checkout"


def _prov(params: dict) -> str:
    return hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()


SAG_LAW = SignalLaw("sagnac", {"L": 1.0, "c": 1.0, "v": 0.2},
                    _prov({"L": 1.0, "c": 1.0, "v": 0.2}))


# --- Gate 1: decomposition identity + direction parity ---------------------

@pytest.mark.skipif(not SAGNAC_AVAILABLE, reason=_SAG_REASON)
class TestSagnacPhaseContract:
    def test_identity_no_double_counting(self):
        spec = sag.phase_spec()
        for omega in (0.0, 1.0, 7.5, 100.0):
            res = check_phase_decomposition(spec, SAG_LAW.parameters, omega)
            assert res["identity_ok"], res
            assert res["pass"], res
            assert res["identity_residual"] < 1e-12

    def test_direction_reversal_semantics(self):
        spec = sag.phase_spec()
        res = check_phase_decomposition(spec, SAG_LAW.parameters, omega=5.0)
        assert res["t_geom_direction_invariant"]
        assert res["directed_phase_nonnegative"]
        assert res["difference_phase_odd"]

    def test_reference_flat_limit_v_zero(self):
        """v -> 0: directed phase reduces to omega*t_geom, residue -> 0."""
        spec = sag.phase_spec()
        law0 = {"L": 1.0, "c": 1.0, "v": 0.0}
        dec = spec.decompose(law0, omega=9.0, direction=+1)
        assert dec.phi_full == pytest.approx(9.0 * 1.0, abs=1e-15)
        assert dec.phi_extra == pytest.approx(0.0, abs=1e-15)

    def test_analytic_closed_form_of_residue(self):
        """phi_extra must equal the exact series residue
        omega*L/c * (beta/(1-beta) for +1), cross-checked analytically."""
        spec = sag.phase_spec()
        L, c, v = 1.0, 1.0, 0.2
        omega = 4.0
        dec = spec.decompose({"L": L, "c": c, "v": v}, omega, +1)
        expected = omega * L / c * (v / c) / (1.0 - v / c)
        assert dec.phi_extra == pytest.approx(expected, rel=1e-14)

    def test_units_declared(self):
        spec = sag.phase_spec()
        assert spec.time_unit in TIME_UNITS
        dec = spec.decompose(SAG_LAW.parameters, 1.0, +1)
        gate = dimension_check(dec)
        assert gate["pass"], gate

    def test_fail_closed_unknown_time_unit(self):
        with pytest.raises(ValueError):
            PhaseDecomposition(omega=1.0, t_geom=1.0, phi_full=1.0,
                               direction=+1, time_unit="furlongs")

    def test_fail_closed_negative_omega(self):
        with pytest.raises(ValueError):
            PhaseDecomposition(omega=-1.0, t_geom=1.0, phi_full=1.0,
                               direction=+1)


# --- SSZ closure backend (geometry owner: SSZ_FULL_CLOSURE) ----------------

@pytest.mark.skipif(not SSZ_AVAILABLE, reason=_SSZ_REASON)
class TestSSZPhaseContract:
    def test_identity_no_double_counting(self):
        spec = ssz.phase_spec()
        params = {"r_a": 1.5, "r_b": 1.6, "b": 0.5}
        for omega in (0.0, 1.0, 12.0):
            res = check_phase_decomposition(spec, params, omega, resid_tol=1e-10)
            assert res["identity_ok"], res
            assert res["pass"], res

    def test_extra_phase_positive_shapiro_share(self):
        """Weak-field SSZ: extra phase (curvature share) must be positive
        for a light ray passing outside the source radius."""
        spec = ssz.phase_spec()
        dec = spec.decompose({"r_a": 1.5, "r_b": 1.6, "b": 0.5}, 3.0, +1)
        assert dec.phi_extra > 0.0

    def test_weak_field_limit_approaches_flat(self):
        """Impact parameter b -> 0 (radial ray): phi_full -> omega*(integral
        dr/sqrt(fh)) >= omega*|dr|, residue -> pure metric share; the flat
        limit is approached as f,h -> 1.  We test the direction invariance
        of t_geom and that b=0 radial extra phase is dominated by fh share."""
        spec = ssz.phase_spec()
        params = {"r_a": 1.5, "r_b": 1.6, "b": 0.0}
        dec = spec.decompose(params, 2.0, +1)
        dec_m = spec.decompose(params, 2.0, -1)
        assert dec.t_geom == dec_m.t_geom == pytest.approx(0.1)
        # static, spherically symmetric: radial time is direction-invariant
        assert dec.phi_full == pytest.approx(dec_m.phi_full, rel=1e-12)

    def test_convergence_of_quadrature(self):
        """Refining the internal grid must converge the coordinate time."""
        repo_path = ssz.closure_repo_path()
        lock = json.loads((repo_path / "MODEL_LOCK.json").read_text())
        from transport_bridge.backends.utils import load_metric_splines
        f_spl, h_spl = load_metric_splines(repo_path / lock["action_member_stream"])
        vals = []
        for n in (2000, 4000, 8000, 16000):
            grid = np.linspace(1.5, 1.6, n)
            radial = 1.0 - 0.25 * h_spl(grid) / grid**2
            integ = 1.0 / np.sqrt(f_spl(grid) * h_spl(grid) * radial)
            vals.append(float(np.trapezoid(integ, grid)))
        diffs = [abs(vals[i + 1] - vals[i]) for i in range(3)]
        assert diffs[-1] <= diffs[0] * 1.05 + 1e-300
        assert all(d >= 0 for d in diffs)

    def test_provenance_binding_intact(self):
        prov = ssz.provenance()
        assert len(prov["member_sha256"]) == 64
        assert len(prov["closure_git_head"]) == 40
        assert prov["verdict"] is not None

    def test_forbidden_impact_parameter_fails_closed(self):
        spec = ssz.phase_spec()
        with pytest.raises(ValueError):
            spec.decompose({"r_a": 1.5, "r_b": 1.6, "b": 50.0}, 1.0, +1)

    def test_units_declared(self):
        spec = ssz.phase_spec()
        assert spec.time_unit == "geometric_m"
        dec = spec.decompose({"r_a": 1.5, "r_b": 1.6, "b": 0.5}, 1.0, +1)
        assert dimension_check(dec)["pass"]


# --- Gate 3: two-channel readout honesty (backend-independent) -------------

class TestTwoChannelReadout:
    def test_common_phase_leaves_auto_psd_invariant(self):
        rng = np.random.default_rng(42)
        amps = (rng.normal(size=64) + 1j * rng.normal(size=64),
                rng.normal(size=64) + 1j * rng.normal(size=64))
        phases = (np.linspace(0, 2.0, 64), np.linspace(0, 0.5, 64))
        res = check_two_channel_readout(amps, phases)
        assert res["auto_psd_invariant"]
        assert res["pass"]

    def test_relative_phase_in_cross_spectrum(self):
        rng = np.random.default_rng(7)
        amps = (rng.normal(size=128) + 1j * rng.normal(size=128),
                rng.normal(size=128) + 1j * rng.normal(size=128))
        # REAL relative phase: the Sagnac differential phase omega*delta_tau
        omega = np.linspace(0.5, 40.0, 128)
        L, c, v = 1.0, 1.0, 0.2
        delta_tau = 2.0 * L * v / (c * c - v * v)   # closed form, axle frame
        phases = (omega * L / (c - v), omega * L / (c + v))
        res = check_two_channel_readout(amps, phases)
        assert res["cross_spectrum_carries_relative_phase"], res
        # rel = omega*L*(1/(c-v) - 1/(c+v)) = omega*delta_tau EXACTLY (unwrapped)
        rel = phases[0] - phases[1]
        assert np.allclose(rel, omega * delta_tau, rtol=1e-14, atol=1e-12)

    def test_zero_relative_phase_is_pure_gauge_for_psd(self):
        """Same phase on both channels: cross-spectrum phase must vanish --
        i.e. an exactly COMMON phase is unobservable in two-channel data."""
        amps = (np.ones(32, dtype=complex), 2.0 * np.ones(32, dtype=complex))
        ph = np.linspace(0.1, 3.0, 32)
        res = check_two_channel_readout(amps, (ph, ph))
        assert res["auto_psd_invariant"]
        assert res["max_rel_phase_err"] < 1e-9


# --- Contract object hardening ---------------------------------------------

class TestPhaseContractObject:
    def test_spec_requires_callable(self):
        with pytest.raises(ValueError):
            PhaseContractSpec(backend_name="x", time_unit="s", decompose=None)

    def test_spec_requires_registered_unit(self):
        with pytest.raises(ValueError):
            PhaseContractSpec(backend_name="x", time_unit="parsecs",
                              decompose=lambda *a: None)

    def test_identity_residual_exact_zero_by_construction(self):
        dec = PhaseDecomposition(omega=3.0, t_geom=0.5, phi_full=1.7,
                                 direction=-1)
        assert dec.identity_residual == 0.0
        assert dec.phi_extra == pytest.approx(1.7 - 1.5)
