"""SSZ-B3/B4: independent global box solve on the frozen spectral export.

Everything here is derived from the exported coefficient matrices only
(no SSZ module import).  The reduced quadratic form declared by the
producer contract is

    L = Ydot^T K Ydot - Y'^T G Y' + Y'^T S Y - Y^T M Y ,
    K = K^T, G = G^T, S^T = -S .

Eigenproblem in lambda = omega^2 form:  A psi = lambda B psi with

    B(phi, psi) = integral phi^T K psi          (mass)
    A(phi, psi) = -[ -integral phi'^T G psi'
                     + 1/2 integral (phi'^T S psi - phi^T S psi')
                     - integral phi^T M psi ]

The S bilinear form is symmetric because S is antisymmetric; the M
antisymmetric part is invisible to the quadratic form.  Assembling the
WEAK form therefore never differentiates S — differentiating the
exported S amplifies grid noise (measured: antisymmetric defect of
M + S'/2 is 16-40% of its scale on this export, see artifact), so the
strong-form route was rejected and documented as a negative result.

P1 (piecewise-linear) Galerkin element blocks on node i..i+1, h, with
midpoint Gbar (symmetrized) and Sbar (antisymmetrized):

    A_ii      += w_i sym(M)_i            B_ii += w_i sym(K)_i
    elem: A_ii += Gbar/h,  A_(i+1,i+1) += Gbar/h
          A_(i,i+1) = -Gbar/h + Sbar/2,  A_(i+1,i) = -Gbar/h - Sbar/2

Boundary treatment: TRUNCATED_BOX — Dirichlet (psi = 0) at both ends of
the declared healthy window r in [1.42861, 1.61289].  This is NOT a QNM
boundary condition; box modes are NOT quasinormal modes.  Physical BCs
are deferred until the closure-side health decision (see
QNM_SPECTROSCOPY_PATH_STATUS artifact).

This module is importable without scipy (lazy imports) so CI can collect
it; all numerics need numpy + scipy at call time.
"""
from __future__ import annotations

import itertools
from typing import Any

import numpy as np

FIELDS = 3
BOUNDARY_TREATMENT = "TRUNCATED_BOX_DIRICHLET_BOTH_WINDOW_ENDS_NOT_QNM_BC"
DISCLAIMER = (
    "box-modes-not-QNMs: all eigenvalues are TRUNCATED_BOX Dirichlet modes of "
    "the exported operator on the healthy window; no quasinormal-mode claim is "
    "made or implied"
)


def sym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a + np.swapaxes(a, 1, 2))


def antisym(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a - np.swapaxes(a, 1, 2))


def assemble_box_pencil(
    r: np.ndarray,
    K: np.ndarray,
    G: np.ndarray,
    S: np.ndarray,
    M: np.ndarray,
    boundary: str = "dirichlet",
) -> tuple[Any, Any]:
    """Symmetric P1-Galerkin pencil A psi = lambda B psi.

    boundary="dirichlet": TRUNCATED_BOX (interior nodes only).
    boundary="natural": all nodes kept (negative-control variant).
    """
    from scipy.sparse import coo_matrix

    if boundary not in ("dirichlet", "natural"):
        raise ValueError(f"unknown boundary treatment: {boundary!r}")
    r = np.asarray(r, float)
    K, G, S, M = (np.asarray(a, float) for a in (K, G, S, M))
    n = len(r)
    if any(a.shape != (n, FIELDS, FIELDS) for a in (K, G, S, M)):
        raise ValueError("operator matrices must have shape (N, 3, 3)")
    if n < 8 or np.any(np.diff(r) <= 0):
        raise ValueError("r must be strictly increasing with at least 8 nodes")
    for name, a, sgn in (("K", K, 1), ("G", G, 1), ("S", S, -1)):
        err = np.max(np.abs(a - sgn * np.swapaxes(a, 1, 2))) / max(1.0, np.max(np.abs(a)))
        if err > 1e-7:
            raise ValueError(f"symmetry contract violated: {name} (scaled err {err:.2e})")

    h = np.diff(r)
    w = np.empty(n)
    w[0] = h[0] / 2
    w[-1] = h[-1] / 2
    w[1:-1] = (h[:-1] + h[1:]) / 2

    rows: list[int] = []
    cols: list[int] = []
    av: list[float] = []
    bv: list[float] = []

    def add(i: int, j: int, blk_a: np.ndarray, blk_b: np.ndarray | None = None) -> None:
        for a in range(FIELDS):
            for b in range(FIELDS):
                rows.append(i * FIELDS + a)
                cols.append(j * FIELDS + b)
                av.append(float(blk_a[a, b]))
                bv.append(0.0 if blk_b is None else float(blk_b[a, b]))

    for i in range(n):
        add(i, i, w[i] * sym(M[i : i + 1])[0], w[i] * sym(K[i : i + 1])[0])
    for i in range(n - 1):
        gbar = sym(0.5 * (G[i] + G[i + 1])[None])[0]
        sbar = antisym(0.5 * (S[i] + S[i + 1])[None])[0]
        hh = h[i]
        add(i, i, gbar / hh)
        add(i + 1, i + 1, gbar / hh)
        add(i, i + 1, -gbar / hh + 0.5 * sbar)
        add(i + 1, i, -gbar / hh - 0.5 * sbar)

    nd = n * FIELDS
    a_mat = coo_matrix((av, (rows, cols)), shape=(nd, nd)).tocsr()
    b_mat = coo_matrix((bv, (rows, cols)), shape=(nd, nd)).tocsr()
    if boundary == "dirichlet":
        keep = np.arange(FIELDS, (n - 1) * FIELDS)
        a_mat = a_mat[keep][:, keep].tocsr()
        b_mat = b_mat[keep][:, keep].tocsr()
    return a_mat, b_mat


def spline_refine(r: np.ndarray, mats: tuple[np.ndarray, ...], factor: int
                  ) -> tuple[np.ndarray, list[np.ndarray]]:
    """Refine the grid by an integer factor; coefficients are cubic-spline
    resampled per (field, field) component.  Original nodes are retained."""
    from scipy.interpolate import CubicSpline

    if int(factor) != factor or factor < 2:
        raise ValueError("factor must be an integer >= 2")
    parts = [r[i] + (r[i + 1] - r[i]) * np.arange(factor) / factor
             for i in range(len(r) - 1)]
    rn = np.concatenate([*parts, r[-1:]])
    out = []
    for a in mats:
        an = np.empty((len(rn), *a.shape[1:]))
        for fa in range(a.shape[1]):
            for fb in range(a.shape[2]):
                an[:, fa, fb] = CubicSpline(r, a[:, fa, fb])(rn)
        out.append(an)
    return rn, out


def mass_min_eigenvalue(b_mat: Any) -> float:
    """Exact minimum eigenvalue of the block-diagonal mass matrix B
    (B only has per-node blocks w_i sym(K_i), so this is O(n) 3x3 problems)."""
    bdiag = b_mat.tocoo()
    # drop explicit zeros (mass rows assembled as 0.0 in stiffness adds)
    nz = bdiag.data != 0.0
    rows_nz, cols_nz, vals_nz = bdiag.row[nz], bdiag.col[nz], bdiag.data[nz]
    # group entries by node block
    per_node: dict[int, dict[tuple[int, int], float]] = {}
    for rr, cc, vv in zip(rows_nz, cols_nz, vals_nz, strict=False):
        if rr // FIELDS != cc // FIELDS:
            # entry crossing node blocks would mean B is NOT block-diagonal
            raise ValueError("mass matrix is not block-diagonal; guard invalid")
        a, b = rr % FIELDS, cc % FIELDS
        per_node.setdefault(rr // FIELDS, {})[(min(a, b), max(a, b))] = float(vv)
    mn = float("inf")
    for _node, entries in per_node.items():
        blk = np.zeros((FIELDS, FIELDS))
        for (a, b), v in entries.items():
            blk[a, b] = v
            blk[b, a] = v
        mn = min(mn, float(np.linalg.eigvalsh(blk)[0]))
    return mn


def solve_near_zero(a_mat: Any, b_mat: Any, k: int = 12,
                    require_mass_positive: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """The k eigenpairs closest to lambda = 0 (shift-invert), sorted by |lambda|.

    Fail-closed: the symmetric-definite shift-invert path requires a positive
    definite mass matrix; with require_mass_positive=True the exact block
    minimum eigenvalue is checked first (ARPACK on an indefinite M silently
    returns spurious values — measured: sign-flipped K produced 0.17 instead
    of the true -lambda)."""
    from scipy.sparse.linalg import eigsh

    if k < 2 or k > a_mat.shape[0] - 2:
        raise ValueError("k out of range for shift-invert")
    if require_mass_positive:
        mn = mass_min_eigenvalue(b_mat)
        if not (mn > 0):
            raise RuntimeError(
                f"mass matrix B is not positive definite (min eig {mn:.6e}); "
                "refusing shift-invert solve — operator input is corrupt or "
                "sign-flipped")
    lam, vec = eigsh(a_mat, k=k, M=b_mat, sigma=0.0, which="LM", tol=1e-12)
    order = np.argsort(np.abs(lam))
    return lam[order], vec[:, order]


def dense_full_spectrum(a_mat: Any, b_mat: Any) -> np.ndarray:
    """Full dense generalized spectrum (coarse-grid probe of the whole
    sector, including any negative lambda)."""
    from scipy.linalg import eigh

    return eigh(a_mat.toarray(), b_mat.toarray(), eigvals_only=True)


def residual_norms(a_mat: Any, b_mat: Any, lam: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Normalized generalized-eigenproblem residuals per mode:
    ||A psi - lambda B psi|| / ((||A||_F + |lambda| ||B||_F) ||psi||_2)."""
    scale_a = max(float(np.sqrt((a_mat.multiply(a_mat)).sum())), 1e-300)
    scale_b = max(float(np.sqrt((b_mat.multiply(b_mat)).sum())), 1e-300)
    out = []
    for j in range(len(lam)):
        res = a_mat @ vec[:, j] - lam[j] * (b_mat @ vec[:, j])
        denom = (scale_a + abs(lam[j]) * scale_b) * np.linalg.norm(vec[:, j])
        out.append(float(np.linalg.norm(res) / denom))
    return np.asarray(out)


def field_fractions(r: np.ndarray, K: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Diagonal-B proxy of the per-field content of each mode:
    f_a = sum_i w_i psi_{i,a}^2 K_i[a,a] / sum_b sum_i w_i psi_{i,b}^2 K_i[b,b].
    Declared proxy: B couples fields within a node, so this is a mixing
    indicator, not an observable."""
    n = len(r)
    interior = vec.reshape(n - 2, FIELDS, -1) if vec.shape[0] == (n - 2) * FIELDS \
        else vec.reshape(n, FIELDS, -1)
    ks = np.stack([sym(K[i : i + 1])[0] for i in range(len(interior))])
    h = np.diff(r)
    w = np.empty(len(interior))
    if vec.shape[0] == (n - 2) * FIELDS:
        w = (h[:-1] + h[1:]) / 2  # interior trapezoid weights
    else:
        w = np.concatenate([[h[0] / 2], (h[:-1] + h[1:]) / 2, [h[-1] / 2]])[-len(interior):]
    diag = np.stack([np.diagonal(ks[i]) for i in range(len(interior))])  # nodes x FIELDS
    p = interior**2 * w[:, None, None] * diag[:, :, None]  # nodes x FIELDS x modes
    num = p.sum(axis=0)  # FIELDS x modes
    return num / np.maximum(num.sum(axis=0, keepdims=True), 1e-300)  # FIELDS x modes


def _common_interior_index(r_fine: np.ndarray, r_coarse: np.ndarray) -> np.ndarray:
    """Positions of the coarse INTERIOR nodes inside the fine node list."""
    idx = []
    for x in r_coarse[1:-1]:
        j = int(np.argmin(np.abs(r_fine - x)))
        if abs(r_fine[j] - x) > 1e-10:
            raise ValueError("coarse interior nodes are not fine-grid nodes")
        idx.append(j)
    return np.asarray(idx)


def match_modes(
    lam_c: np.ndarray, vec_c: np.ndarray, r_c: np.ndarray, K_c: np.ndarray,
    lam_f: np.ndarray, vec_f: np.ndarray, r_f: np.ndarray, K_f: np.ndarray,
) -> tuple[dict[int, int], np.ndarray]:
    """Match coarse modes to fine modes by B-inner products evaluated on the
    common interior nodes (coarse interior nodes are fine-grid nodes for the
    integer-factor refinement family).  Returns (coarse->fine map, overlap
    matrix).  Fail-closed: raises if the grids are incompatible."""
    from scipy.optimize import linear_sum_assignment

    idx_f = _common_interior_index(r_f, r_c)
    nc = len(r_c) - 2
    nf = len(r_f) - 2
    vc = vec_c.reshape(nc, FIELDS, -1)
    vf = vec_f.reshape(nf, FIELDS, -1)[idx_f]
    h = np.diff(r_c)
    w = (h[:-1] + h[1:]) / 2
    kint = sym(K_c[1:-1])
    nc_nodes = len(idx_f)
    weight = np.zeros((nc_nodes * FIELDS, nc_nodes * FIELDS))
    for i in range(nc_nodes):
        weight[i * FIELDS:(i + 1) * FIELDS, i * FIELDS:(i + 1) * FIELDS] = w[i] * kint[i]
    overlap = vc.reshape(nc_nodes * FIELDS, -1).T @ (
        weight @ vf.reshape(nc_nodes * FIELDS, -1))
    ci, fi = linear_sum_assignment(-np.abs(overlap))
    return {int(a): int(b) for a, b in zip(ci, fi, strict=False)}, overlap


def ladder_for_L(
    mats: dict[str, np.ndarray],
    levels: tuple[int, ...] = (1, 2, 4, 8),
    k: int = 12,
) -> dict[str, Any]:
    """Resolution ladder on one L sector.  Level 1 = native export grid;
    higher levels = integer-factor cubic-spline refinement of the exported
    coefficients.  Returns per-level spectra, matched relative shifts
    between successive levels, convergence flags (criterion: relative
    lambda shift between the two finest levels < 1e-6 with matched overlap
    > 0.5), and eigenproblem residuals."""
    r = np.asarray(mats["r"], float)
    K, G, S, M = (np.asarray(mats[x], float) for x in "KGSM")
    per_level: dict[str, dict[str, Any]] = {}
    grid_cache: dict[int, tuple[np.ndarray, tuple[np.ndarray, ...]]] = {}
    for fac in levels:
        if fac == 1:
            rr, mm = r, (K, G, S, M)
        else:
            rr, mm = spline_refine(r, (K, G, S, M), fac)
        grid_cache[fac] = (rr, mm)
        a_mat, b_mat = assemble_box_pencil(rr, *mm)
        lam, vec = solve_near_zero(a_mat, b_mat, k=k)
        res = residual_norms(a_mat, b_mat, lam, vec)
        per_level[str(fac)] = {
            "n_nodes": len(rr),
            "dof": int(a_mat.shape[0]),
            "lambda": lam.tolist(),
            "max_residual": float(np.max(res)),
            "min_lambda": float(np.min(lam)),
        }
        if fac == levels[-1]:
            per_level[str(fac)]["field_fractions"] = field_fractions(
                rr, mm[0], vec).tolist()
    shifts: dict[str, list[float]] = {}
    overlaps: dict[str, list[float]] = {}
    flags: dict[str, bool] = {}
    for lo, hi in itertools.pairwise(levels):
        r_c, m_c = grid_cache[lo]
        r_f, m_f = grid_cache[hi]
        a_c, b_c = assemble_box_pencil(r_c, *m_c)
        lam_c, vec_c = solve_near_zero(a_c, b_c, k=k)
        a_f, b_f = assemble_box_pencil(r_f, *m_f)
        lam_f, vec_f = solve_near_zero(a_f, b_f, k=k)
        pairs, overlap = match_modes(
            lam_c, vec_c, r_c, m_c[0], lam_f, vec_f, r_f, m_f[0])
        sh = []
        ov = []
        for a in range(k):
            b = pairs.get(a)
            if b is None:
                sh.append(float("nan"))
                ov.append(0.0)
                continue
            sh.append(float(abs(lam_f[b] - lam_c[a]) / max(abs(lam_f[b]), 1e-300)))
            ov.append(float(abs(overlap[a, b])))
        key = f"{lo}to{hi}"
        shifts[key] = sh
        overlaps[key] = ov
    final_key = f"{levels[-2]}to{levels[-1]}"
    flags = {
        "criterion": "relative lambda shift between two finest levels < 1e-6",
        "per_mode": [
            bool(shifts[final_key][i] < 1e-6 and overlaps[final_key][i] > 0.5)
            for i in range(k)
        ],
        "min_match_overlap_finest": float(np.min(overlaps[final_key])),
    }
    return {
        "boundary_treatment": BOUNDARY_TREATMENT,
        "levels": [int(x) for x in levels],
        "per_level": per_level,
        "matched_relative_shifts": shifts,
        "match_overlaps": overlaps,
        "convergence": flags,
        "modes_reported": k,
        "disclaimer": DISCLAIMER,
    }


def cross_check_producer(
    mats: dict[str, np.ndarray],
    native_npz: Any,
    k: int = 8,
) -> dict[str, Any]:
    """Independent re-assembly vs the producer's frozen native-window npz
    (a separate implementation in the closure repo).  Returns relative
    differences of the first k near-zero eigenvalues.  Fail-closed if the
    producer artifact is missing."""
    r = np.asarray(mats["r"], float)
    K, G, S, M = (np.asarray(mats[x], float) for x in "KGSM")
    a_mat, b_mat = assemble_box_pencil(r, K, G, S, M)
    lam, _ = solve_near_zero(a_mat, b_mat, k=k)
    prod = np.asarray(native_npz["omega2"], float)[:k]
    if len(prod) < k:
        raise ValueError("producer npz has fewer modes than requested")
    if not np.all(np.isfinite(prod)):
        raise ValueError("producer npz contains nonfinite eigenvalues")
    rel = np.abs(lam - prod) / np.abs(prod)
    return {
        "comparison": "bridge P1 assembly vs producer native_window npz (independent"
                      " implementation, same discretization family)",
        "lambda_bridge": lam.tolist(),
        "lambda_producer": prod.tolist(),
        "relative_diff": rel.tolist(),
        "max_relative_diff": float(np.max(rel)),
    }


# ---------------------------------------------------------------------------
# Negative controls (B10).  Each returns a dict with passed/fired semantics.
# ---------------------------------------------------------------------------

def control_corrupted_operator(mats: dict[str, np.ndarray], rel_pert: float = 1e-3,
                               k: int = 6) -> dict[str, Any]:
    """Perturb ONE entry of K by rel_pert (relative).  The eigenvalues must
    move measurably (fail-closed threshold 1e-9 relative)."""
    r = np.asarray(mats["r"], float)
    K = np.array(mats["K"], float, copy=True)
    i = len(r) // 2
    entry_old = float(K[i, 0, 1])
    K[i, 0, 1] = entry_old * (1.0 + rel_pert)
    K[i, 1, 0] = K[i, 0, 1]  # keep the symmetry contract
    a0, b0 = assemble_box_pencil(r, mats["K"], mats["G"], mats["S"], mats["M"])
    lam0, _ = solve_near_zero(a0, b0, k=k)
    a1, b1 = assemble_box_pencil(r, K, mats["G"], mats["S"], mats["M"])
    lam1, _ = solve_near_zero(a1, b1, k=k)
    shifts = np.abs(lam1 - lam0) / np.abs(lam0)
    fired = float(np.max(shifts)) > 1e-9
    return {
        "control": "corrupted_operator_single_entry_1e-3",
        "perturbed_entry": f"K[{i},0,1] *= {1.0 + rel_pert}",
        "max_relative_shift": float(np.max(shifts)),
        "per_mode_shifts": shifts.tolist(),
        "fired": bool(fired),
        "passed": bool(fired),
    }


def control_sign_flipped_K(mats: dict[str, np.ndarray], k: int = 4) -> dict[str, Any]:
    """K -> -K makes the mass matrix B indefinite; the symmetric solver must
    REFUSE (fail-closed) or the spectrum must shift massively."""
    r = np.asarray(mats["r"], float)
    a_mat, b_mat = assemble_box_pencil(
        r, -np.asarray(mats["K"], float), mats["G"], mats["S"], mats["M"])
    refused = False
    max_shift = float("nan")
    try:
        lam1, _ = solve_near_zero(a_mat, b_mat, k=k)
        a0, b0 = assemble_box_pencil(r, mats["K"], mats["G"], mats["S"], mats["M"])
        lam0, _ = solve_near_zero(a0, b0, k=k)
        max_shift = float(np.max(np.abs(lam1[:k] + lam0[:k]) / np.abs(lam0[:k])))
    except RuntimeError as exc:
        refused = True
        max_shift = float("inf")
        refusal = str(exc)
    else:
        refusal = ""
    passed = refused or (np.isfinite(max_shift) and max_shift > 0.5)
    out = {
        "control": "sign_flipped_K",
        "solver_refused": bool(refused),
        "refusal_message": refusal if refused else "",
        "max_relative_shift_vs_original": max_shift if not refused else None,
        "passed": bool(passed),
    }
    return out


def control_wrong_boundary(mats: dict[str, np.ndarray], k: int = 6) -> dict[str, Any]:
    """Natural (no-Dirichlet) variant must differ measurably from the
    TRUNCATED_BOX spectrum — proves the boundary treatment matters and the
    declared treatment is not vacuous."""
    r = np.asarray(mats["r"], float)
    a0, b0 = assemble_box_pencil(r, mats["K"], mats["G"], mats["S"], mats["M"])
    lam0, _ = solve_near_zero(a0, b0, k=k)
    a1, b1 = assemble_box_pencil(
        r, mats["K"], mats["G"], mats["S"], mats["M"], boundary="natural")
    lam1, _ = solve_near_zero(a1, b1, k=k)
    shift = float(np.max(np.abs(lam1[:k] - lam0[:k]) / np.abs(lam0[:k])))
    return {
        "control": "wrong_boundary_natural_instead_of_truncated_box",
        "max_relative_shift": shift,
        "fired": bool(shift > 1e-2),
        "passed": bool(shift > 1e-2),
    }


def control_wrong_member_cross_check(
    mats: dict[str, np.ndarray], native_npz: Any, r_scale: float = 1.01, k: int = 6
) -> dict[str, Any]:
    """A self-consistently-hashed export whose radial grid is STRETCHED (wrong
    member geometry) passes the loader (the hash binds file integrity, not
    physics) but MUST fail the cross-solver comparison against the frozen
    producer artifact.  Note: a uniform radial TRANSLATION is an exact no-op
    for the box pencil (A, B depend only on node spacings) — measured
    7.1e-10 — and is therefore NOT a valid discriminator; the stretch is."""
    r = np.asarray(mats["r"], float) * r_scale
    a_mat, b_mat = assemble_box_pencil(r, mats["K"], mats["G"], mats["S"], mats["M"])
    lam, _ = solve_near_zero(a_mat, b_mat, k=k)
    prod = np.asarray(native_npz["omega2"], float)[:k]
    rel = np.abs(lam - prod) / np.abs(prod)
    fired = float(np.max(rel)) > 1e-6
    return {
        "control": f"wrong_member_radial_stretch_x{r_scale}_detected_by_cross_solver",
        "max_relative_diff_vs_producer": float(np.max(rel)),
        "per_mode_diffs": rel.tolist(),
        "fired": bool(fired),
        "passed": bool(fired),
    }
