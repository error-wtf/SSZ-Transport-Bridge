#!/usr/bin/env python3
"""High-precision Evans control solver (B6 Step 1) — mpmath rewrite.

Purpose (Carmen/Lino, 12.10.26): classify the five Re(w)~0.05 scan
"GATED_POLE_CANDIDATE"s.  Double precision saturates at |D|~1e-16 with
machine resolution ~3e-9 (proved at the certified pilot).  This solver
repeats the SAME physics in 30-80 significant digits:

  * ODE  Y'' = -G^-1 (A Y' + (B - w^2 K) Y),  A = dG/dr - S,  B = M + S'/2
    (identical to FrozenOperator.coeffs / StableEvansV3._rhs).
  * Coefficients: frozen npz tables (hash-verified), same not-a-knot
    cubic spline as the double solver, reconstructed EXACTLY per cell as
    cubic Hermite from nodal values + nodal spline derivatives, then
    evaluated in mpmath.  Spline error O(h^4) is identical in both
    solvers and precision-independent (recorded caveat).
  * Interior: Frobenius seed (same REGULAR_EXPONENTS), r_in -> r_match.
  * Exterior: incoming selection = 3 smallest-Re(mu) pencil branches at
    r_seed (same rule as LocalPencil.incoming_selection), with mu and
    eigenvectors refined by shifted inverse iteration + Rayleigh
    quotient in mpmath; alpha from the same balance formula
    (left eigenvector of L = mu^2 G + B - w^2 K for its smallest
    eigenvalue, obtained by shift-0 inverse iteration on L^H):
        alpha = -u.(R_j a0) / (mu * u.(G a0)),
        R_j = mu G + mu^2 G' + A mu + B - w^2 K.
    G' at r_seed is the EXACT spline derivative (the double solver used
    a central FD with h=1e-4 r; recorded difference, O(h^2)=1e-8 rel).
  * Propagation: fixed-step RK4 with periodic column renormalization;
    per-segment log scale factors tracked analytically.
  * Evans object: raw M(w) = [stack(Qi) | stack(Qe)] (6x6).  Diagnostics
      - det M (mp.det)  -> the HOLOMORPHIC object for a later
        argument-principle rebuild,
      - smallest eigenvalue of M^H M -> sigma_min^2 (mp.eigsy),
      - sigma_min of the GRAM-SCHMIDT-orthonormalized M for direct
        comparison with the double solver's D (QR-equivalent),
      - full singular value spectrum.
    No QR sign fix, no phase resets anywhere: M(w) is analytic in w.

Regions: pilot (certified zero, must stay a zero and now resolve BELOW
the double floor 3e-9), candidate (Re w ~ 0.05, zero or plateau?),
negative (control, must stay clean).

Usage:
  b6_mpmath_evans_v2.py pilot|candidate|negative --dps 30
        [--refine 1] [--center-only]
Output: artifacts/SSZ_B6_STEP1_MPMATH_EVANS_<REGION>_DPS<dps>.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

from transport_bridge.independent_qnm_3dof_v1 import (  # noqa: E402
    REGULAR_EXPONENTS, FrozenOperator)

NPZ = BRIDGE / "data" / "frozen_specs" / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"

PROBES = {
    "pilot":     {"center": (0.22642588643752676, 0.08455094304111954),
                  "deltas": [3e-2, 1e-2, 3e-3, 1e-3]},
    "candidate": {"center": (0.051875, 0.0285),
                  "deltas": [3e-2, 1e-2, 3e-3, 1e-3]},
    "negative":  {"center": (0.06, 0.045),
                  "deltas": [3e-2, 1e-2, 3e-3, 1e-3]},
}

STEPS_INT, STEPS_EXT = 720, 560      # RK4 steps (x --refine)
RENORM_EVERY = 12                    # RK4 steps between column renorms


# ==========================================================================
# coefficients: exact mpmath re-evaluation of the double solver's splines
# ==========================================================================
class MPCoeffs:
    """Per-cell cubic Hermite reconstruction of the FrozenOperator
    not-a-knot cubic splines: value+derivative at both nodes determine
    the unique cell polynomial, so this IS the double solver's
    interpolant, evaluated in arbitrary precision."""

    def __init__(self, op: FrozenOperator):
        from scipy.interpolate import CubicSpline
        assert np.abs(np.imag(op.G)).max() == 0, "operator must be real"
        self.r = np.asarray(op.r, float)
        self.val = {}
        self.der = {}
        for name, arr in (("G", op.G), ("S", op.S), ("M", op.M), ("K", op.K)):
            real = np.real(arr)
            self.val[name] = real
            self.der[name] = CubicSpline(self.r, real, axis=0).derivative()(self.r)
        # mpmath-cached node values
        self.mp_val = {}
        self.mp_der = {}
        for name in self.val:
            v = self.val[name]
            d = self.der[name]
            self.mp_val[name] = [[[mp.mpf(float(v[j, a, b])) for b in range(3)]
                                  for a in range(3)] for j in range(len(self.r))]
            self.mp_der[name] = [[[mp.mpf(float(d[j, a, b])) for b in range(3)]
                                  for a in range(3)] for j in range(len(self.r))]

    def _entry(self, x: mp.mpf, name: str, want_deriv: bool) -> mp.matrix:
        r = self.r
        j = int(np.searchsorted(r, float(x))) - 1
        j = max(0, min(j, len(r) - 2))
        x0 = mp.mpf(float(r[j]))
        h = mp.mpf(float(r[j + 1])) - x0
        t = (x - x0) / h
        t2 = t * t
        vj = self.mp_val[name]
        dj = self.mp_der[name]
        out = mp.matrix(3, 3)
        if want_deriv:
            c0 = (6 * t2 - 6 * t) / h
            c1 = 3 * t2 - 4 * t + 1
            c2 = (-6 * t2 + 6 * t) / h
            c3 = 3 * t2 - 2 * t
        else:
            t3 = t2 * t
            c0 = 2 * t3 - 3 * t2 + 1
            c1 = (t3 - 2 * t2 + t) * h
            c2 = -2 * t3 + 3 * t2
            c3 = (t3 - t2) * h
        for a in range(3):
            for b in range(3):
                out[a, b] = (c0 * vj[j][a][b] + c1 * dj[j][a][b]
                             + c2 * vj[j + 1][a][b] + c3 * dj[j + 1][a][b])
        return out

    def GABK(self, x: mp.mpf):
        G = self._entry(x, "G", False)
        S = self._entry(x, "S", False)
        A = self._entry(x, "G", True) - S
        B = self._entry(x, "M", False) + 0.5 * self._entry(x, "S", True)
        K = self._entry(x, "K", False)
        return G, A, B, K


# ==========================================================================
# propagation
# ==========================================================================
def lu_solve_mat(G, R):
    """Solve G X = R column-wise (mp.lu_solve is single-RHS only)."""
    n, k = R.rows, R.cols
    X = mp.matrix(n, k)
    for c in range(k):
        col = mp.matrix(n, 1)
        for i in range(n):
            col[i, 0] = R[i, c]
        sol = mp.lu_solve(G, col)
        for i in range(n):
            X[i, c] = sol[i, 0]
    return X


def renorm(Y, dY):
    """Unit-normalize each stacked column; returns log10 of scale products."""
    logs = mp.mpf(0)
    for c in range(3):
        n2 = mp.mpf(0)
        for i in range(3):
            n2 += Y[i, c] ** 2 + dY[i, c] ** 2
        nrm = mp.sqrt(n2)
        if nrm == 0:
            continue
        logs += mp.log(nrm)
        for i in range(3):
            Y[i, c] /= nrm
            dY[i, c] /= nrm
    return logs


def propagate(coeffs: MPCoeffs, Y, dY, r0: float, r1: float,
              n_steps: int, w2):
    """Fixed-step RK4 r0 -> r1; periodic renormalization, scales tracked."""
    r = mp.mpf(repr(r0))
    h = (mp.mpf(repr(r1)) - r) / n_steps
    log_scale = mp.mpf(0)

    def rhs(Yv, dYv, G, A, B, K):
        return dYv, -lu_solve_mat(G, A @ dYv + (B - w2 * K) @ Yv)

    for s in range(n_steps):
        G, A, B, K = coeffs.GABK(r)
        k1y, k1d = rhs(Y, dY, G, A, B, K)
        k2y, k2d = rhs(Y + h / 2 * k1y, dY + h / 2 * k1d, G, A, B, K)
        k3y, k3d = rhs(Y + h / 2 * k2y, dY + h / 2 * k2d, G, A, B, K)
        k4y, k4d = rhs(Y + h * k3y, dY + h * k3d, G, A, B, K)
        Y = Y + h / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
        dY = dY + h / 6 * (k1d + 2 * k2d + 2 * k3d + k4d)
        r += h
        if (s + 1) % RENORM_EVERY == 0 and s + 1 < n_steps:
            log_scale += renorm(Y, dY)
    return Y, dY, log_scale


# ==========================================================================
# exterior incoming seed at r_seed
# ==========================================================================
def refine_eigenpair(C: mp.matrix, mu0: mp.mpc, iters: int = 3):
    """Shifted inverse iteration + Rayleigh quotient on companion C."""
    n = C.rows
    v = mp.matrix(n, 1)
    for i in range(n):
        v[i, 0] = mp.mpf(1) / mp.mpf(i + 1)
    v /= mp.sqrt(mp.norm(v))
    mu = mu0
    for _ in range(iters):
        A = mp.matrix(C)
        for i in range(n):
            A[i, i] -= mu
        try:
            with mp.extraprec(mp.mp.dps * 2):
                w = mp.lu_solve(A, v)
        except (ZeroDivisionError, mp.libmp.NoConvergence):
            break
        v = w / mp.sqrt(mp.norm(w))
        Mv = C @ v
        with mp.extraprec(mp.mp.dps * 2):
            mu = mp.fdot(v, Mv)
    return mu, v


def smallest_left_eigvec(L: mp.matrix, iters: int = 8):
    """Eigenvector of L for the smallest-|eigenvalue| mode via shift-0
    inverse iteration on L^H (same object as numpy eig of L.conj().T)."""
    n = L.rows
    u = mp.matrix(n, 1)
    for i in range(n):
        u[i, 0] = mp.mpf(1) / mp.mpf(i + 2)
    u /= mp.sqrt(mp.norm(u))
    LH = L.H
    for _ in range(iters):
        with mp.extraprec(mp.mp.dps * 2):
            u = mp.lu_solve(LH, u)
        u /= mp.sqrt(mp.norm(u))
    return u


def exterior_seed(coeffs: MPCoeffs, op: FrozenOperator, w: mp.mpc,
                  r_seed: float):
    """3 incoming branches (smallest Re mu) at r_seed, refined in mpmath.
    Mirrors LocalPencil.branches + incoming_selection."""
    R = mp.mpf(repr(r_seed))
    G, A, B, K = coeffs.GABK(R)
    Gp = coeffs._entry(R, "G", True)
    w2 = w * w
    # companion: [[0, I], [G^-1 (w2 K - B), 0]]
    Low = lu_solve_mat(G, w2 * K - B)
    C = mp.matrix(6, 6)
    for i in range(3):
        for j in range(3):
            C[i, j + 3] = 1
            C[i + 3, j] = Low[i, j]
    # double-precision seeds (only starting guesses) from the SAME tables
    Gn, An, Bn, Kn = op.coeffs(r_seed)
    w2n = complex(w.real if hasattr(w, 'real') else float(mp.re(w)),
                  float(mp.im(w)))
    w2n = w2n * w2n
    C0n = w2n * Kn - Bn
    M6 = np.block([[np.zeros((3, 3), dtype=complex), np.eye(3, dtype=complex)],
                   [np.linalg.solve(np.asarray(Gn, dtype=complex),
                                    np.asarray(C0n, dtype=complex)),
                    np.zeros((3, 3), dtype=complex)]])
    ev = np.linalg.eig(M6)[0]
    order = np.argsort(ev.real)
    out = []
    for j in order[:3]:
        mu, vec6 = refine_eigenpair(C, mp.mpc(float(ev[j].real),
                                              float(ev[j].imag)))
        a0 = mp.matrix(3, 1)
        for i in range(3):
            a0[i, 0] = vec6[i, 0]
        a0 /= mp.sqrt(mp.norm(a0))
        # alpha balance formula
        L = mu * mu * G + B - w2 * K
        u = smallest_left_eigvec(L)
        Rj = mu * G + mu * mu * Gp + A * mu + B - w2 * K
        num = mp.fdot(u, Rj @ a0)
        den = mu * mp.fdot(u, G @ a0)
        alpha = -num / den if abs(den) > mp.mpf(10) ** (-mp.mp.dps + 5) else mp.mpf(0)
        out.append((mu, alpha, a0))
    return out


# ==========================================================================
# matching + diagnostics
# ==========================================================================
def gram_schmidt_cols(M: mp.matrix):
    """Orthonormalize columns (modified GS); returns log det |R| factor."""
    n, k = M.rows, M.cols
    Q = mp.matrix(n, k)
    logdet = mp.mpf(0)
    for j in range(k):
        v = mp.matrix(n, 1)
        for i in range(n):
            v[i, 0] = M[i, j]
        for i in range(j):
            r = mp.fdot(Q_T_col(Q, i), v)
            for t in range(n):
                v[t, 0] -= r * Q[t, i]
        rn = mp.sqrt(mp.norm(v))
        logdet += mp.log(rn) if rn != 0 else mp.mpf(0)
        for t in range(n):
            Q[t, j] = v[t, 0] / rn
    return Q, logdet


def Q_T_col(Q, i):
    n = Q.rows
    c = mp.matrix(n, 1)
    for t in range(n):
        c[t, 0] = mp.conj(Q[t, i])
    return c


def svals_of(M: mp.matrix):
    """Full singular value spectrum via mp.eig of G = M^H M.  Verified:
    product(svals) == |det M| to 15+ digits at the pilot (12.10.)."""
    H = M.H @ M
    vals, _ = mp.eig(H)
    vals = sorted(vals, key=lambda z: mp.re(z))
    return [max(mp.re(v), mp.mpf(0)) for v in vals]


def svals_of_sqrt(M_eigs):
    return [mp.sqrt(v) for v in M_eigs]


def evaluate_point(coeffs, op, w: mp.mpc, r_match=30.0, r_seed=55.0,
                   refine=1) -> dict:
    w2 = w * w
    # ---- interior: Frobenius seed -------------------------------------
    r_in = float(op.r_in)
    Y = mp.matrix(3, 3)
    dY = mp.matrix(3, 3)
    for c, s in enumerate(REGULAR_EXPONENTS):
        rr = mp.power(mp.mpf(repr(r_in)), mp.mpf(repr(float(s))))
        Y[c, c] = rr
        dY[c, c] = mp.mpf(repr(float(s))) * rr / mp.mpf(repr(r_in))
    Qi_Y, Qi_dY, log_si = propagate(coeffs, Y, dY, r_in, r_match,
                                    STEPS_INT * refine, w2)
    # ---- exterior: incoming from r_seed --------------------------------
    branches = exterior_seed(coeffs, op, w, r_seed)
    Y0 = mp.matrix(3, 3)
    dY0 = mp.matrix(3, 3)
    for c, (mu, alpha, a0) in enumerate(branches):
        for i in range(3):
            Y0[i, c] = a0[i, 0]
            dY0[i, c] = (mu + alpha / mp.mpf(repr(r_seed))) * a0[i, 0]
    Qe_Y, Qe_dY, log_se = propagate(coeffs, Y0, dY0, r_seed, r_match,
                                    STEPS_EXT * refine, w2)
    # ---- matching matrix (raw) -----------------------------------------
    M = mp.matrix(6, 6)
    for c in range(3):
        for i in range(3):
            M[i, c] = Qi_Y[i, c]
            M[i + 3, c] = Qi_dY[i, c]
            M[i, c + 3] = Qe_Y[i, c]
            M[i + 3, c + 3] = Qe_dY[i, c]
    with mp.extraprec(mp.mp.dps * 2):
        det_raw = mp.det(M)
        Qo, logdQ = gram_schmidt_cols(M)
    sv = svals_of(M)          # eigenvalues of M^H M, ascending
    sv_o = svals_of(Qo)       # cross-check path (orthonormalized M)
    svals = [mp.sqrt(v) for v in reversed(sv)]    # report descending
    svals_o = [mp.sqrt(v) for v in reversed(sv_o)]
    return {
        "w": [mp.nstr(mp.re(w), 18), mp.nstr(mp.im(w), 18)],
        "det_abs": mp.nstr(abs(det_raw), 12),
        "log_abs_det": mp.nstr(mp.log(abs(det_raw)), 12),
        "sigma_min": mp.nstr(svals[-1], 12),
        "sigma_min_crosscheck": mp.nstr(svals_o[-1], 12),
        "svals": [mp.nstr(s, 8) for s in svals],
        "svals_orthonormalized": [mp.nstr(s, 8) for s in svals_o],
        "log_scale_int": mp.nstr(log_si, 10),
        "log_scale_ext": mp.nstr(log_se, 10),
        "logdet_orthonormalization": mp.nstr(logdQ, 10),
        "incoming_mus": [mp.nstr(mu, 12) for mu, _, _ in branches],
    }


# ==========================================================================
# main
# ==========================================================================
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("region", choices=list(PROBES))
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--refine", type=int, default=1)
    ap.add_argument("--center-only", action="store_true")
    args = ap.parse_args()
    mp.mp.dps = args.dps

    op = FrozenOperator(str(NPZ))
    t0 = time.time()
    coeffs = MPCoeffs(op)
    print(f"coeffs prepared in {time.time()-t0:.1f}s", flush=True)

    reg = PROBES[args.region]
    cr, ci = reg["center"]
    w0 = mp.mpc(mp.mpf(repr(cr)), mp.mpf(repr(ci)))
    deltas = [0.0] if args.center_only else [0.0] + reg["deltas"]

    out = {"region": args.region, "dps": args.dps, "refine": args.refine,
           "steps": [STEPS_INT * args.refine, STEPS_EXT * args.refine],
           "renorm_every": RENORM_EVERY,
           "center": [cr, ci], "r_match": 30.0, "r_seed": 55.0,
           "notes": [
               "same frozen npz (hash-verified) + same spline as double solver",
               "G'(r_seed) = exact spline derivative (double solver: central FD h=1e-4 r)",
               "det M is the holomorphic Evans object; svals via eigsy(M^H M)",
           ],
           "probes": []}
    for delta in deltas:
        w = w0 if delta == 0 else w0 + mp.mpc(repr(delta), 0)
        t1 = time.time()
        res = evaluate_point(coeffs, op, w, refine=args.refine)
        res["delta"] = delta
        res["wall_s"] = round(time.time() - t1, 1)
        out["probes"].append(res)
        print(f"delta={delta:<8g} log|det E|={res['log_abs_det']:>12}  "
              f"sigma_min={res['sigma_min']:>12}  ({res['wall_s']}s)",
              flush=True)

    art = (BRIDGE / "artifacts"
           / f"SSZ_B6_STEP1_MPMATH_EVANS_{args.region.upper()}_"
             f"DPS{args.dps}_R{args.refine}.json")
    art.write_text(json.dumps(out, indent=2) + "\n")
    print("artifact:", art, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
