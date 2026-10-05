"""Real SSZ closure backend contracts (skips without the closure checkout).

Proves the SAME bridge contract structure wraps the REAL SSZ_FULL_CLOSURE
transport — with zero physics in the bridge.  Per the direction-contract
rule, static spherical SSZ declares NO sign-swap involution; the checks
here assert observable validity, reproducibility and provenance binding.
"""
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from transport_bridge import SignalLaw

try:
    from transport_bridge.backends import ssz_closure_backend as sag_cb
    SSZ_CLOSURE_AVAILABLE = True
except ImportError:
    SSZ_CLOSURE_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not SSZ_CLOSURE_AVAILABLE,
    reason="requires the SSZ_FULL_CLOSURE checkout")


def _provenance(params: dict) -> str:
    return hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()


# b must sit INSIDE the photon cone of the frozen member at r0:
# b_crit(r0=1.55) = r0/sqrt(f) ~= 2.489 on ELECTRIC_PRODUCTION_MEMBER_CURRENT.
# b=2.5 is outside -> correct ValueError guard.  Tests use b=2.0.
LAW = SignalLaw("ssz_closure", {"inward": 0.0, "b": 2.0, "r0": 1.55},
                _provenance({"inward": 0.0, "b": 2.0, "r0": 1.55}))


def test_forward_run_all_observables_finite():
    op = sag_cb.make_operator(LAW)
    out = op.run()
    for name in sag_cb.OBSERVABLES.scalar:
        assert name in out, f"missing declared observable {name}"
        assert out[name] == out[name], f"{name} is NaN"
    # geodesic conservation (when the ray stayed in domain):
    if "geodesic_domain_exit" not in out:
        assert out["max_dE"] < 1e-8
        assert out["max_dL"] < 1e-8
        assert out["max_norm_res"] < 1e-8


def test_forward_is_reproducible():
    op = sag_cb.make_operator(LAW)
    a = op.run()
    b = op.run()
    assert a == b


def test_provenance_binds_to_closure():
    p = sag_cb.provenance()
    assert p["closure_git_head"] != "UNKNOWN"
    assert "SSZ_FULL_CLOSURE" in p["closure_src"]


def test_no_sign_swap_law_declared():
    """Direction-contract rule: static spherical SSZ declares NO sign-swap
    involution (unlike Sagnac's v -> -v).  The backend must NOT invent
    one."""
    assert sag_cb.OBSERVABLES.directional == (), (
        "static spherical SSZ must not declare directional sign-flip "
        "observables — the Sagnac symmetry is not universal physics")


def test_inverse_contract_deferred():
    """The inverse problem is explicitly deferred: inverse='' in the
    observable set.  Fabricating an inverse here would violate the
    contract — asserted structurally instead."""
    assert sag_cb.OBSERVABLES.inverse == ""
