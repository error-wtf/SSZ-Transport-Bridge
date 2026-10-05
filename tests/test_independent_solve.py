"""SSZ-B3/B4 tests: independent global box solve, convergence ladder,
negative controls.  All fail-closed.

Skips only when the frozen export (closure checkout) is absent (CI).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from transport_bridge import independent_solve as ix
    from transport_bridge import spectral_import as si

    _EXPORT_AVAILABLE = si.EXPORT_NPZ.exists() and si.EXPORT_JSON.exists()
except ImportError:  # pragma: no cover
    ix = None
    si = None
    _EXPORT_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not _EXPORT_AVAILABLE, reason="frozen spectral export not on this machine")

NATIVE_DIR = si.CLOSURE_REPO / "data/generated/spectral/native_window" if si else None


@pytest.fixture(scope="module")
def export():
    return si.load_export()


@pytest.fixture(scope="module")
def mats6(export):
    return si.load_operator(export, 6)


# ---------------------------------------------------------------- B3 solve

def test_b3_pencil_symmetric_and_fail_closed(mats6):
    a, b = ix.assemble_box_pencil(
        mats6["r"], mats6["K"], mats6["G"], mats6["S"], mats6["M"])
    assert abs(a - a.T).max() == 0.0
    assert abs(b - b.T).max() == 0.0
    # Dirichlet truncation removed exactly 2 nodes x 3 fields
    assert a.shape[0] == (len(mats6["r"]) - 2) * 3


def test_b3_symmetry_contract_violation_refused(mats6):
    k_bad = np.array(mats6["K"], copy=True)
    k_bad[10, 0, 1] += 1e5  # break symmetry
    with pytest.raises(ValueError, match="symmetry contract"):
        ix.assemble_box_pencil(mats6["r"], k_bad, mats6["G"], mats6["S"], mats6["M"])


def test_b3_mass_guard_positive(mats6):
    _a, b = ix.assemble_box_pencil(
        mats6["r"], mats6["K"], mats6["G"], mats6["S"], mats6["M"])
    assert ix.mass_min_eigenvalue(b) > 0


def test_b3_solve_finite_positive_sector(mats6):
    a, b = ix.assemble_box_pencil(
        mats6["r"], mats6["K"], mats6["G"], mats6["S"], mats6["M"])
    lam, vec = ix.solve_near_zero(a, b, k=6)
    assert np.all(np.isfinite(lam))
    # first box modes in the healthy window are positive (omega^2 > 0)
    assert np.all(lam > 0)
    # residuals at solver precision
    res = ix.residual_norms(a, b, lam, vec)
    assert np.max(res) < 1e-10


def test_b3_no_qnm_language_in_artifact():
    """The disclaimer must be explicit; no QNM claim may appear as a
    certification."""
    art = Path(__file__).resolve().parents[1] / (
        "artifacts/SSZ_B3_B4_INDEPENDENT_SOLVE_V1.json")
    if not art.exists():  # artifact built locally, not on CI
        pytest.skip("artifact not built on this machine")
    data = json.loads(art.read_text())
    assert "box-modes-not-QNMs" in data["disclaimer"]
    assert data["boundary_treatment_detail"]["name"] == "TRUNCATED_BOX"
    assert "physical BC layer is deferred" in data[
        "boundary_treatment_detail"]["why_not_qnm_bc"]


# ------------------------------------------------------------ B4 ladder

def test_b4_ladder_converges(mats6):
    lad = ix.ladder_for_L(mats6, levels=(1, 2, 4), k=6)
    # successive refinement must reduce matched shifts monotonically for
    # the low modes (1->2 coarser than 2->4)
    s12 = np.array(lad["matched_relative_shifts"]["1to2"][:4])
    s24 = np.array(lad["matched_relative_shifts"]["2to4"][:4])
    assert np.nanmax(s24) < np.nanmax(s12)
    # finest-level shifts below 1e-6 for the low modes
    assert np.nanmax(s24) < 1e-6
    # all reported modes converged at the finest level
    assert sum(lad["convergence"]["per_mode"]) == 6


def test_b4_refinement_retains_original_nodes(mats6):
    r = np.asarray(mats6["r"], float)
    rn, mm = ix.spline_refine(r, (mats6["K"],), 2)
    assert np.allclose(rn[0::2], r, atol=1e-12)
    assert len(mm[0]) == len(rn)


# -------------------------------------------------- negative controls (B10)

def test_control_corrupted_operator_fires(mats6):
    out = ix.control_corrupted_operator(mats6)
    assert out["fired"] and out["passed"]
    assert out["max_relative_shift"] > 1e-9


def test_control_sign_flipped_K_refused(mats6):
    out = ix.control_sign_flipped_K(mats6)
    assert out["passed"]
    assert out["solver_refused"], (
        "sign-flipped K must FAIL the solve (mass guard) — silent spurious "
        "eigenvalues are unacceptable")


def test_control_wrong_boundary_fires(mats6):
    out = ix.control_wrong_boundary(mats6)
    assert out["fired"] and out["passed"]
    assert out["max_relative_shift"] > 1e-2


def test_control_wrong_member_stretch_fires(mats6):
    npz_p = NATIVE_DIR / "L6_native_window_spectroscopy.npz"
    if not npz_p.exists():
        pytest.skip("producer native npz not on this machine")
    out = ix.control_wrong_member_cross_check(mats6, np.load(npz_p))
    assert out["fired"] and out["passed"]
    assert out["max_relative_diff_vs_producer"] > 1e-6


def test_control_translation_is_noop_documented(mats6):
    """A uniform radial translation is an exact no-op for the box pencil
    (A, B depend only on spacings) — this is measured physics of the
    discretization, documented so nobody mistakes it for insensitivity."""
    r_shift = np.asarray(mats6["r"], float) + 0.01
    a0, b0 = ix.assemble_box_pencil(
        mats6["r"], mats6["K"], mats6["G"], mats6["S"], mats6["M"])
    a1, b1 = ix.assemble_box_pencil(
        r_shift, mats6["K"], mats6["G"], mats6["S"], mats6["M"])
    lam0, _ = ix.solve_near_zero(a0, b0, k=4)
    lam1, _ = ix.solve_near_zero(a1, b1, k=4)
    assert np.max(np.abs(lam1 - lam0) / np.abs(lam0)) < 1e-9


def test_loader_hash_mismatch_refused(tmp_path, mats6):
    """Tampered npz (bit flipped) must be refused by the loader — the B1
    fail-closed hash binding."""
    repo = tmp_path / "repo"
    (repo / "data/generated/spectral").mkdir(parents=True)
    npz_bytes = bytearray(si.EXPORT_NPZ.read_bytes())
    npz_bytes[-64] ^= 0xFF  # corrupt one byte inside the arrays payload
    bad_npz = repo / "data/generated/spectral/FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.npz"
    bad_npz.write_bytes(bytes(npz_bytes))
    prov = json.loads(si.EXPORT_JSON.read_text())
    (repo / "data/generated/spectral/FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.json"
     ).write_text(json.dumps(prov))
    with pytest.raises(ValueError, match="hash mismatch"):
        si.load_export(repo)


# ------------------------------------------------------------ B7 cross

def test_cross_solver_vs_producer(mats6):
    npz_p = NATIVE_DIR / "L6_native_window_spectroscopy.npz"
    if not npz_p.exists():
        pytest.skip("producer native npz not on this machine")
    cc = ix.cross_check_producer(mats6, np.load(npz_p), k=8)
    assert cc["max_relative_diff"] < 1e-6
