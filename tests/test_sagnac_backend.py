"""Sagnac reference backend contracts (skip without the reference checkout)."""
import hashlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from transport_bridge import (
        ObservableSet,
        PropagatedState,
        SignalLaw,
        check_forward,
        check_inverse_round_trip,
        check_parity,
    )
    from transport_bridge.backends import sagnac_backend as sag
    SAGNAC_AVAILABLE = True
except ImportError:
    SAGNAC_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not SAGNAC_AVAILABLE,
    reason="requires the validated Sagnac reference checkout")

if SAGNAC_AVAILABLE:
    def _provenance(params: dict) -> str:
        return hashlib.sha256(
            __import__("json").dumps(params, sort_keys=True).encode()).hexdigest()

    SAG_LAW = SignalLaw("sagnac", {"L": 1.0, "c": 1.0, "v": 0.2},
                        _provenance({"L": 1.0, "c": 1.0, "v": 0.2}))
    SAG_OBS = ObservableSet(
        directional=("dt", "delta_phi", "delta_tau_signed"),
        scalar=("total", "abs_delta_tau"), inverse="v")

    def test_forward_parity_contract():
        op = sag.make_operator(SAG_LAW)
        plus = op.run()
        minus = op.solve(SAG_LAW.perturbed("v", -0.4))
        checks = check_forward(
            op, SAG_OBS,
            {"dt": (plus["dt"], minus["dt"]),
             "delta_phi": (plus["delta_phi"], minus["delta_phi"]),
             "delta_tau_signed": (plus["delta_tau_signed"], minus["delta_tau_signed"]),
             "total": (plus["total"], minus["total"]),
             "abs_delta_tau": (plus["abs_delta_tau"], minus["abs_delta_tau"])},
            atol=1e-12)
        assert all(checks.values()), checks

    def test_parity_states():
        op = sag.make_operator(SAG_LAW)
        plus_vals = op.run()
        minus_vals = op.solve(SAG_LAW.perturbed("v", -0.4))
        checks = check_parity(
            PropagatedState(+1, plus_vals), PropagatedState(-1, minus_vals),
            SAG_OBS, atol=1e-12)
        assert all(checks.values()), checks

    def test_inverse_round_trip():
        op = sag.make_operator(SAG_LAW)
        forward = op.run()
        res = check_inverse_round_trip(
            SAG_LAW, "v", forward, lambda f: sag.inverse_v(f, c=1.0), atol=1e-12)
        assert res["ok"], res
