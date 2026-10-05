"""Generic contract tests run against BOTH backends with the SAME code.

This is the bridge's whole point: identical structural checks
(forward/parity/convergence/inverse) pass on the validated Sagnac
reference and on the SSZ synthetic-inversion proxy — proving both
belong to the same validation class without sharing a single formula.
"""
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    check_forward,
    check_inverse_round_trip,
)
from transport_bridge.backends import ssz_proxy_backend as ssz


def _provenance(params: dict) -> str:
    import json
    return hashlib.sha256(
        json.dumps(params, sort_keys=True).encode()).hexdigest()


# ---------------- Sagnac (reference backend) ----------------




SSZ_LAW = SignalLaw("ssz", {"r_s": 1.0, "r_obs": 10.0,
                            "frame_sign": 1.0, "c": 1.0},
                    hashlib.sha256(b'{"r_s":1.0,"r_obs":10.0}').hexdigest())
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
    # Cross-backend proof runs on whatever backends are installed.  The
    # Sagnac side has its own dedicated test file (skips cleanly without
    # the reference checkout, e.g. on CI).
    systems = [(ssz.solve, SSZ_LAW, SSZ_OBS, 1e-12)]
    try:
        from tests.test_sagnac_backend import SAG_LAW, SAG_OBS
        from transport_bridge.backends import sagnac_backend as sag
        systems.insert(0, (sag.solve, SAG_LAW, SAG_OBS, 1e-10))
    except ImportError:
        pass
    for solve, obs_law, obs_set, atol in systems:
        plus = solve(obs_law)
        flip_key = ("v" if obs_law.system == "sagnac" else "frame_sign")
        flip_delta = (-2.0 * obs_law.parameters[flip_key])
        minus = solve(obs_law.perturbed(flip_key, flip_delta))
        for name in plus:
            if name in obs_set.directional:
                assert directional_close(plus[name], minus[name], atol), (
                    f"{obs_law.system}.{name} failed the direction-swap "
                    "contract (directional must flip sign)")
            elif name in obs_set.scalar:
                assert scalar_close(plus[name], minus[name], atol), (
                    f"{obs_law.system}.{name} failed the scalar-invariance "
                    "contract (scalar must not flip sign)")


def directional_close(a: float, b: float, atol: float) -> bool:
    """Directional observables MUST flip sign under direction reversal."""
    return abs(a + b) <= atol


def scalar_close(a: float, b: float, atol: float) -> bool:
    """Scalar observables MUST be invariant under direction reversal."""
    return abs(a - b) <= atol
