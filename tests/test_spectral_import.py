"""SSZ-B1/B2 tests: frozen-export provenance and independent reconstruction.

The bridge consumes FROZEN_SPECTRAL_OPERATOR_EXPORT_V1 from the closure
checkout WITHOUT importing any SSZ module.  Skips (with reason) only when
the local closure checkout or export artifact is absent (CI).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from transport_bridge import spectral_import as si

    _EXPORT_AVAILABLE = si.EXPORT_NPZ.exists() and si.EXPORT_JSON.exists()
except ImportError:  # pragma: no cover
    si = None
    _EXPORT_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not _EXPORT_AVAILABLE, reason="frozen spectral export not on this machine")


def test_b1_hash_binding():
    export = si.load_export()
    assert export["npz_sha256"] == export["prov"]["npz_sha256"]
    assert export["prov"]["export"] == "FROZEN_SPECTRAL_OPERATOR_EXPORT_V1"
    assert export["prov"]["known_limitations"]  # limitations declared, not hidden


def test_b1_member_window_declared():
    prov = si.load_export()["prov"]
    assert prov["healthy_window_u"] == [0.62, 0.70]
    assert prov["L_values"] == [6, 12, 20, 42, 110, 420, 1000]
    # every L passed the fail-closed health diagnostics at production time
    for lval, diag in prov["health_L"].items():
        assert diag["pass"], f"L={lval} health failed"
        assert diag["min_eig_K"] > 0


def test_b2_operator_reconstruction():
    export = si.load_export()
    for L in (6, 12, 20):
        mats = si.load_operator(export, L)
        K, G, S = mats["K"], mats["G"], mats["S"]
        # symmetry conventions match the golden producer (declared contract):
        assert np.max(np.abs(K - np.swapaxes(K, 1, 2))) / max(
            1.0, np.max(np.abs(K))) < 1e-7
        assert np.max(np.abs(G - np.swapaxes(G, 1, 2))) / max(
            1.0, np.max(np.abs(G))) < 1e-7
        assert np.max(np.abs(S + np.swapaxes(S, 1, 2))) / max(
            1.0, np.max(np.abs(S))) < 1e-7
        # finite, positive pivots:
        for p in ("Dh1", "DeltaV", "pivotA0"):
            assert np.all(np.isfinite(mats[p])) and np.min(np.abs(mats[p])) > 0


def test_b2_independent_local_spectra_finite():
    export = si.load_export()
    mats = si.load_operator(export, 6)
    res = si.independent_local_spectra(mats, n_check=8)
    assert res["eigen_traces"].shape[1] == 3  # 3-field coupled sector
    assert np.all(np.isfinite(res["eigen_traces"]))
    # bridge-side health agrees with the declared production-side health:
    assert float(np.min(res["eigen_traces"])) > 0
