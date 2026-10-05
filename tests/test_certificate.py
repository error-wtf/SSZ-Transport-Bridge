"""Certificate regressions: no vacuous success or ambiguous provenance."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_cert_tool = ROOT / "tools/build_bridge_certificate.py"
spec = importlib.util.spec_from_file_location("certificate", _cert_tool)
certificate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(certificate)


def test_git_head_is_full():
    assert len(certificate.git_head(ROOT)) == 40


def test_incomplete_contracts_cannot_pass():
    assert not certificate.contracts_pass({"forward": True})
    assert not certificate.contracts_pass({})
    flags = dict.fromkeys(certificate.CONTRACTS, True)
    assert certificate.contracts_pass(flags)
    for name in flags:
        assert not certificate.contracts_pass({**flags, name: False})


def test_nonfinite_observables_rejected():
    for value in (float("inf"), -float("inf"), float("nan")):
        assert not certificate.finite_observables({"x": value}, ("x",))
    assert not certificate.finite_observables({}, ("x",))
    assert not certificate.finite_observables({}, ())
    assert certificate.finite_observables({"x": 1.0}, ("x",))
