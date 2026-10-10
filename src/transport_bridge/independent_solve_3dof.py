"""Independent FULL 3-DOF solve on the research-branch production operator.

Implements the SUPERSEDED_V4_2CHANNEL_VERDICTS_V1 directive: the
independent validator consumes the FULL 3-field operator
Psi = (psi, dphi, V)^T from closure spectroscopy-real-data-20261004
(V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF), NOT the superseded 2-channel
reduction behind FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.

Physics (declared by the producer contract, re-derived here):
    L = Ydot^T K Ydot - Y'^T G Y' + Y'^T S Y - Y^T M Y  (+ R Y, R = 0)
    K = K^T, G = G^T, S^T = -S, R = -P11/2 proven identically zero.

Time-harmonic ansatz gives the symmetric generalized eigenproblem
A y = lambda B y with lambda = omega_t^2:

    B(phi, y) = integral phi^T K y                      (kinetic/mass)
    A(phi, y) = -[ -integral phi'^T G y'
                   + 1/2 integral (phi'^T S y - phi^T S' ... )
                   - integral phi^T M y ]

Assembled in the WEAK form (S never differentiated; the S bilinear form
1/2 integral (phi'^T S y - phi^T S y') is symmetric because S^T = -S).

Boundary treatment: TRUNCATED_BOX Dirichlet at both window ends
(r = 0.05 and r = 1.0, the declared integration window).  This is NOT a
QNM boundary condition: the physical outgoing-wave layer (B6) remains an
open gate.  Box modes constrain the operator; they are NOT QNM claims.

Independence: no SSZ closure import; only the frozen npz + its hash.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

RESEARCH_REPO = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE-research")
EXPORT_REL = Path("data/generated/spectral/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz")
EXPECTED_SHA256 = "ad2bee7d3994097479933c1eea787f3aa381648d098ffd36aab32a4c436db2ef"

BOUNDARY_TREATMENT = "TRUNCATED_BOX_DIRICHLET_WINDOW_ENDS_NOT_QNM_BC"
DISCLAIMER = (
    "box-modes-not-QNMs on the FULL 3-DOF production operator; no "
    "quasinormal-mode claim is made or implied (B6 outgoing-wave layer OPEN)"
)


def _sha256_file(p: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_export_3dof(repo: Path = RESEARCH_REPO) -> dict[str, Any]:
    """Hash-bound load of the V3 3-DOF production export."""
    npz_p = repo / EXPORT_REL
    if not npz_p.exists():
        raise FileNotFoundError(f"export missing: {npz_p}")
    sha = _sha256_file(npz_p)
    if sha != EXPECTED_SHA256:
        raise ValueError(
            f"export hash mismatch: got {sha}, expected {EXPECTED_SHA256}"
        )
    z = np.load(npz_p)
    r = z["r"]
    L = int(z["L"][0])
    win_mask = r <= 1.0
    return {
        "repo": str(repo),
        "sha256": sha,
        "r": r,
        "L": L,
        "K": z["K_phys"],
        "G": z["G_phys"],
        "S": z["S_phys"],
        "M": z["M_phys"],
        "R": z["R_phys"],
        "win_mask": win_mask,
        "r_win": r[win_mask],
    }


def structural_gates(ex: dict[str, Any]) -> dict[str, Any]:
    """Elementwise structural checks on the FULL 3-channel operators."""
    K, G, S, M, R = ex["K"], ex["G"], ex["S"], ex["M"], ex["R"]
    Ksym = float(np.max(np.abs(K - np.swapaxes(K, 1, 2))))
    Gsym = float(np.max(np.abs(G - np.swapaxes(G, 1, 2))))
    Sdef = float(np.max(np.abs(S + np.swapaxes(S, 1, 2))))
    Sscale = float(np.max(np.abs(S)))
    Rmax = float(np.max(np.abs(R)))
    # mass guard: min eigenvalue of K on the window (fail-closed)
    Kw = K[ex["win_mask"]]
    eigK = np.linalg.eigvalsh(0.5 * (Kw + np.swapaxes(Kw, 1, 2)))
    min_eigK = float(eigK.min())
    scaleK = float(np.max(np.abs(Kw)))
    gates = {
        "K_symmetric_max_defect": Ksym,
        "G_symmetric_max_defect": Gsym,
        "S_antisym_max_defect": Sdef,
        "S_scale": Sscale,
        "R_identically_zero": Rmax == 0.0,
        "R_max_abs": Rmax,
        "min_eig_K_window": min_eigK,
        "mass_guard_positive": min_eigK > 0.0,
        "K_scale": scaleK,
    }
    gates["all_pass"] = bool(
        Ksym < 1e-6 * scaleK
        and Gsym < 1e-6 * scaleK
        and Sdef < 1e-6 * Sscale
        and Rmax == 0.0
        and min_eigK > 0.0
    )
    return gates


def _assemble_window(r: np.ndarray, K: np.ndarray, G: np.ndarray,
                     S: np.ndarray, M: np.ndarray) -> tuple:
    """P1 Galerkin assembly on the given window mesh.

    Dirichlet at both window ends: node 0 and node n-1 are eliminated
    (all three channels each).
    """
    n = r.size
    h = np.diff(r)
    A = np.zeros((3 * n, 3 * n))
    B = np.zeros((3 * n, 3 * n))
    for i in range(n - 1):
        he = h[i]
        Gbar = 0.5 * (G[i] + G[i + 1]); Gbar = 0.5 * (Gbar + Gbar.T)
        Sbar = 0.5 * (S[i] + S[i + 1]); Sbar = 0.5 * (Sbar - Sbar.T)
        Mbar = 0.5 * (M[i] + M[i + 1]); Mbar = 0.5 * (Mbar + Mbar.T)
        Kb = 0.5 * (K[i] + K[i + 1]);   Kb = 0.5 * (Kb + Kb.T)
        w1 = he / 2.0
        for a in (i, i + 1):
            sl = slice(3 * a, 3 * a + 3)
            B[sl, sl] += w1 * Kb
            A[sl, sl] += w1 * Mbar   # +int phi M y enters A with plus sign
        gs = Gbar / he
        s2 = Sbar / 2.0
        s11 = slice(3 * i, 3 * i + 3)
        s22 = slice(3 * (i + 1), 3 * (i + 1) + 3)
        A[s11, s11] += gs
        A[s22, s22] += gs
        A[s11, s22] += -gs + s2
        A[s22, s11] += -gs - s2
    drop = np.concatenate([np.arange(3), np.arange(3 * (n - 1), 3 * n)])
    keep = np.setdiff1d(np.arange(3 * n), drop)
    return A[np.ix_(keep, keep)], B[np.ix_(keep, keep)]


def _assemble_window_sparse(r: np.ndarray, K: np.ndarray, G: np.ndarray,
                            S: np.ndarray, M: np.ndarray):
    """Sparse P1 Galerkin assembly (COO triplets), otherwise identical to
    _assemble_window.  Returns (A_csc, B_csc, keep) with keep = kept dofs."""
    from scipy.sparse import coo_matrix
    n = r.size
    h = np.diff(r)
    rows, cols, va, vb = [], [], [], []
    for i in range(n - 1):
        he = h[i]
        Gbar = 0.5 * (G[i] + G[i + 1]); Gbar = 0.5 * (Gbar + Gbar.T)
        Sbar = 0.5 * (S[i] + S[i + 1]); Sbar = 0.5 * (Sbar - Sbar.T)
        Mbar = 0.5 * (M[i] + M[i + 1]); Mbar = 0.5 * (Mbar + Mbar.T)
        Kb = 0.5 * (K[i] + K[i + 1]);   Kb = 0.5 * (Kb + Kb.T)
        w1 = he / 2.0
        for a in (i, i + 1):
            for p in range(3):
                for q in range(3):
                    rows.append(3 * a + p); cols.append(3 * a + q)
                    va.append(w1 * Mbar[p, q]); vb.append(w1 * Kb[p, q])
        gs = Gbar / he
        s2 = Sbar / 2.0
        for p in range(3):
            for q in range(3):
                for (a, b, v) in ((i, i, gs[p, q]),
                                  (i + 1, i + 1, gs[p, q]),
                                  (i, i + 1, -gs[p, q] + s2[p, q]),
                                  (i + 1, i, -gs[p, q] - s2[p, q])):
                    rows.append(3 * a + p); cols.append(3 * b + q)
                    va.append(v); vb.append(0.0)  # B: no stencil contribution
    N = 3 * n
    A = coo_matrix((va, (rows, cols)), shape=(N, N)).tocsc()
    B = coo_matrix((vb, (rows, cols)), shape=(N, N)).tocsc()
    drop = np.concatenate([np.arange(3), np.arange(3 * (n - 1), 3 * n)])
    keep = np.setdiff1d(np.arange(N), drop)
    return A[np.ix_(keep, keep)].tocsc(), B[np.ix_(keep, keep)].tocsc(), keep


def box_solve_ladder(ex: dict[str, Any], ladder=(400, 800, 1600, 3200),
                     k: int = 12) -> dict[str, Any]:
    """B4-style resolution ladder with nested-mesh mode overlaps.

    Overlap: fine-grid eigenvector restricted to nodes shared with the
    coarse mesh, dotted against every coarse eigenvector (meshes are
    nested by construction of the strided subsampling).  Converged iff
    max rel |lambda shift| between the two finest levels < 1e-6 AND
    best-match overlap > 0.5 for every mode.
    """
    from scipy.sparse.linalg import eigsh
    r_full = ex["r_win"]
    m = r_full.size
    levels = {}
    for n in ladder:
        stride = max(1, (m - 1) // (n - 1))
        idx = np.arange(0, m, stride)
        if idx[-1] != m - 1:
            idx = np.append(idx, m - 1)
        r = r_full[idx]
        A, B, keep = _assemble_window_sparse(
            r, ex["K"][ex["win_mask"]][idx], ex["G"][ex["win_mask"]][idx],
            ex["S"][ex["win_mask"]][idx], ex["M"][ex["win_mask"]][idx])
        vals, vecs = eigsh(A, k=k, M=B, sigma=0.0, which="LM")
        order = np.argsort(vals)
        levels[n] = {
            "idx": idx, "keep": keep,
            "lambdas": vals[order], "vecs": vecs[:, order],
        }
    n_fine, n_coarse = ladder[-1], ladder[-2]
    lf, lc = levels[n_fine], levels[n_coarse]
    shifts = np.abs(lf["lambdas"] - lc["lambdas"]) / np.abs(lc["lambdas"])
    fine_global = lf["idx"][lf["keep"] // 3]   # dof -> node index
    coarse_global = lc["idx"][lc["keep"] // 3]
    _, f_pos, c_pos = np.intersect1d(fine_global, coarse_global,
                                     return_indices=True)
    overlaps = []
    for j in range(k):
        v_f = lf["vecs"][f_pos, j]
        best = 0.0
        for i in range(k):
            v_c = lc["vecs"][c_pos, i]
            denom = np.linalg.norm(v_f) * np.linalg.norm(v_c)
            if denom > 0:
                best = max(best, abs(float(v_f @ v_c)) / denom)
        overlaps.append(best)
    return {
        "levels": {str(n): levels[n]["lambdas"].tolist() for n in ladder},
        "finest_pair": [n_coarse, n_fine],
        "max_rel_shift_finest": float(shifts.max()),
        "min_overlap_finest": float(min(overlaps)),
        "overlaps": overlaps,
        "shifts": shifts.tolist(),
        "converged": bool(shifts.max() < 1e-6 and min(overlaps) > 0.5),
    }


def solve_box_modes(ex: dict[str, Any], n_nodes: int, k: int = 12,
                    sigma: float | None = None):
    """Dirichlet box solve on a coarsened window mesh (dense reference path).

    n_nodes: number of window nodes to use (strided subsample of the
    declared 3500-node window).  Returns (lambdas, max_residual).
    """
    from scipy.sparse import csc_matrix
    from scipy.sparse.linalg import eigsh

    r_full = ex["r_win"]
    m = r_full.size
    stride = max(1, (m - 1) // (n_nodes - 1))
    idx = np.arange(0, m, stride)
    if idx[-1] != m - 1:
        idx = np.append(idx, m - 1)
    r = r_full[idx]
    K = ex["K"][ex["win_mask"]][idx]
    G = ex["G"][ex["win_mask"]][idx]
    S = ex["S"][ex["win_mask"]][idx]
    M = ex["M"][ex["win_mask"]][idx]
    A, B = _assemble_window(r, K, G, S, M)
    A_s = csc_matrix(A)
    B_s = csc_matrix(B)
    # shift-invert at sigma=0 (spectrum anchored near zero; plain 'SA'
    # does not converge on this pencil — measured, see B3 artifact note)
    vals, vecs = eigsh(A_s, k=k, M=B_s, sigma=0.0, which="LM")
    # normalized residual ||A v - lambda B v|| / (||A v|| + |lambda| ||B v||)
    Av = A @ vecs
    Bv = B @ vecs
    res = np.linalg.norm(Av - vals[np.newaxis, :] * Bv, axis=0) / (
        np.linalg.norm(Av, axis=0)
        + np.abs(vals) * np.linalg.norm(Bv, axis=0)
    )
    return np.sort(vals), float(res.max())
