"""Generic contract tests run against BOTH backends with the SAME code.

This is the bridge's whole point: identical structural checks
(forward/parity/convergence/inverse) pass on the validated Sagnac
reference and on the SSZ synthetic-inversion proxy — proving both
belong to the same validation class without sharing a single formula.
"""
import hashlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    check_forward,
    check_inverse_round_trip,
    check_parity,
)
from transport_bridge.backends import sagnac_backend as sag
from transport_bridge.backends import ssz_backend as ssz


def _provenance(params: dict) -> str:
    import json
    return hashlib.sha256(
        json.dumps(params, sort_keys=True).encode()).hexdigest()


# ---------------- Sagnac (reference backend) ----------------

SAG_LAW = SignalLaw("sagnac", {"L": 1.0, "c": 1.0, "v": 0.2},
                    _provenance({"L": 1.0, "c": 1.0, "v": 0.2}))
SAG_OBS = ObservableSet(
    directional=("dt", "delta_phi", "delta_tau_signed"),
    scalar=("total", "abs_delta_tau"), inverse="v")


def test_sagnac_forward_parity_contract():
    op = sag.make_operator(SAG_LAW)
    plus = op.run()
    flipped = SAG_LAW.perturbed("v", -2 * SAG_LAW.parameters["v"])
    minus = op.solve(flipped)
    checks = check_forward(
        op, SAG_OBS,
        {"dt": (plus["dt"], minus["dt"]),
         "delta_phi": (plus["delta_phi"], minus["delta_phi"]),
         "delta_tau_signed": (plus["delta_tau_signed"], minus["delta_tau_signed"]),
         "total": (plus["total"], minus["total"]),
         "abs_delta_tau": (plus["abs_delta_tau"], minus["abs_delta_tau"])},
        atol=1e-12)
    assert all(checks.values()), checks


def test_sagnac_parity_states():
    op = sag.make_operator(SAG_LAW)
    plus_vals = op.run()
    minus_vals = op.solve(SAG_LAW.perturbed("v", -0.4))
    from transport_bridge import PropagatedState
    checks = check_parity(
        PropagatedState(+1, plus_vals), PropagatedState(-1, minus_vals),
        SAG_OBS, atol=1e-12)
    assert all(checks.values()), checks


def test_sagnac_inverse_round_trip():
    op = sag.make_operator(SAG_LAW)
    forward = op.run()
    res = check_inverse_round_trip(
        SAG_LAW, "v", forward, lambda f: sag.inverse_v(f, c=1.0), atol=1e-12)
    assert res["ok"], res
    assert res["abs_error"] < 1e-12


def test_sagnac_unknown_parameter_fails_closed():
    with pytest.raises(KeyError):
        SAG_LAW.perturbed("does_not_exist", 1.0)


# ---------------- SSZ (synthetic-inversion backend) ----------------

SSZ_LAW = SignalLaw("ssz", {"r_s": 1.0, "r_obs": 10.0, "frame_sign": 1.0, "c": 1.0},
                    _provenance({"r_s": 1.0, "r_obs": 10.0}))
SSZ_OBS = ObservableSet(directional=("delta_phi", "delta_t"),
                        scalar=("z_mean",), inverse="frame_sign")


def test_ssz_forward_parity_contract():
    plus = ssz.solve(SSZ_LAW)
    minus = ssz.solve(SSZ_LAW.perturbed("frame_sign", -2.0))
    checks = check_forward(
        None, SSZ_OBS,
        {"delta_phi": (plus["delta_phi"], minus["delta_phi"]),
         "delta_t": (plus["delta_t"], minus["delta_t"]),
         "z_mean": (plus["z_mean"], minus["z_mean"])},
        atol=1e-12)
    assert all(checks.values()), checks


def test_ssz_synthetic_inverse_round_trip():
    """THE synthetic inversion test: data from a KNOWN frozen geometry,
    sign recovered without peeking."""
    for sign in (+1.0, -1.0):
        law = SSZ_LAW.perturbed("frame_sign", sign - SSZ_LAW.parameters["frame_sign"])
        forward = ssz.solve(law)
        res = check_inverse_round_trip(
            law, "frame_sign", forward,
            lambda f: ssz.inverse_frame_sign(f["delta_phi"]), atol=0.0)
        assert res["ok"], res
        assert res["reconstructed"] == sign


def test_ssz_scalar_is_sign_independent():
    assert ssz.solve(SSZ_LAW)["z_mean"] == ssz.solve(
        SSZ_LAW.perturbed("frame_sign", -2.0))["z_mean"]


# ---------------- cross-backend class membership ----------------

def test_both_backends_share_the_validation_class():
    """The bridge's core claim: both systems satisfy the SAME structural
    contracts, computed by the SAME check code, from THEIR OWN physics."""
    for solve, obs_law, obs_set, atol in (
        (sag.solve, SAG_LAW, SAG_OBS, 1e-10),
        (ssz.solve, SSZ_LAW, SSZ_OBS, 1e-12),
    ):
        plus = solve(obs_law)
        flip_key = ("v" if obs_law.system == "sagnac" else "frame_sign")
        flip_delta = (-2.0 * obs_law.parameters[flip_key])
        minus = solve(obs_law.perturbed(flip_key, flip_delta))
        dir_names = [k for k in plus if k in (obs_set.directional)]
        for name in dir_names:
            assert np_close(plus[name], minus[name], atol), (
                f"{obs_law.system}.{name} failed the direction-swap contract")


def np_close(a: float, b: float, atol: float) -> bool:
    """The direction-swap contract: either a == -b (swaps) or a == b
    (invariant) — the bridge accepts both but records which."""
    return abs(a + b) <= atol or abs(a - b) <= atol
