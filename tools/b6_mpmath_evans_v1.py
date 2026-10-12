#!/usr/bin/env python3
"""High-precision Evans control solver (B6 Step 1, package 6) — mpmath.

Purpose (Carmen/Lino, 12.10.26): classify the five Re(w)~0.05 "GATED_
POLE_CANDIDATE" scan outputs.  Double-precision delta-scaling showed
they sit on the SVD roundoff floor (machine resolution ~3e-9 at the
certified pilot).  This solver repeats the SAME physics with 30-80
significant digits so the precision floor drops far below any plausible
signal, and evaluates three regions:

  pilot      : certified reference (0.22643+0.08455i)  -> must stay a zero
  candidate  : Re(w)~0.05 region                       -> zero or plateau?
  negative   : control region without expected pole    -> must stay clean

Physics (identical to src/transport_bridge/stable_evans_v3.py):
  * ODE  Y'' = -G^-1 (A Y' + (B - w^2 K) Y),  A = dG/dr - S,  B = M + dS/2.
  * Coefficients: the frozen npz tables (hash-verified) with the SAME
    not-a-knot cubic-spline interpolant as the double solver, evaluated
    via exact Hermite reconstruction from nodal values + nodal spline
    derivatives (piecewise identical polynomial => same function, in
    arbitrary precision).
  * Interior: Frobenius seed (same exponents), r_in -> r_match.
  * Exterior: incoming pencil branches (Re mu < 0) at r_seed from the
    companion matrix  C(w) = [[0, I], [G^-1 (w^2 K - B), 0]]  (same
    linearization as LocalPencil), eigenvectors + eigenvalues refined
    by inverse iteration + Rayleigh quotient in mpmath (double-precision
    seeds only as starting guesses), then propagated r_seed -> r_match.
  * Indicial alpha from the same balance formula as LocalPencil.branches.
  * Matching matrix M(omega) = [stack(Qi) | stack(Qe)] (6x6, raw columns
    with recorded per-segment scale factors).  Diagnostics per probe:
      - E = smallest singular value (inverse iteration on M^H M),
      - |det M| and log|det M| (mp.det — the holomorphic object for a
        future argument-principle rebuild),
      - s_max estimate via power iteration on M^H M.
    No QR sign fix, no phase resets: the recorded scale factors carry
    all normalizations analytically.

Precision caveats (recorded in the artifact):
  * Table data is float64 by contract; spline interpolation error O(h^4)
    is identical in both solvers and does NOT shrink with dps.  The
    dps ladder (30/50/80) isolates integration + algebra roundoff.
  * RK4 fixed step; truncation controlled by --refine (step halving)
    and cross-checked against the dps ladder.

Usage:
  b6_mpmath_evans_v1.py pilot|candidate|negative --dps 30 [--refine 1]
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

R_IN, R_MATCH, R_SEED = 0.05, 30.0, 55.0
STEPS_INT, STEPS_EXT = 600, 480     # RK4 steps, scaled by --refine
RENORM_EVERY = 10                   # steps between column renormalizations


# ==========================================================================
# coefficients: exact mpmath re-evaluation of the double solver's splines
# ==========================================================================
class MPCoeffs:
    """Piecewise Hermite reconstruction of the FrozenOperator cubic
    splines.  On each interval the not-a-knot cubic is the unique cubic
    matching (value, derivative) at both nodes, so this reproduces the
    double solver's interpolant exactly, evaluated in mpmath."""

    def __init__(self, op: FrozenOperator):
        from scipy.interpolate import CubicSpline
        self.r = np.asarray(op.r, float)
        self.val = {}
        self.der = {}
        for name, arr in (("G", op.G), ("S", op.S), ("M", op.M),
                          ("K", op.K)):
            self.val[name] = arr.astype(np.complex128)
            sp = CubicSpline(self.r, arr, axis=0)     # complex spline ok
            self.der[name] = sp.derivative()(self.r)

    def _bracket(self, x: mp.mpf) -> tuple[int, mp.mpf, mp.mpf]:
        j = int(np.searchsorted(self.r, float(x))) - 1
        j = max(0, min(j, len(self.r) - 2))
        x0 = mp.mpf(float(self.r[j]))
        x1 = mp.mpf(float(self.r[j + 1]))
        return j, x0, x1 - x0

    def _entry(self, x, name, want_deriv):
        j, x0, h = self._bracket(x)
        t = (x - x0) / h
        t2 = t * t
        t3 = t2 * t
        out = mp.matrix(3, 3)
        for a in range(3):
            for b in range(3):
                z = self.val[name][j, a, b]
                z1 = self.val[name][j + 1, a, b]
                d0 = mp.mpc(mp.mpf(float(z.real)), mp.mpf(float(z.imag)))
                d1 = mp.mpc(mp.mpf(float(z1.real)), mp.mpf(float(z1.imag)))
                dd0 = self.der[name][j, a, b]
                dd1 = self.der[name][j + 1, a, b]
                e0 = mp.mpc(mp.mpf(float(dd0.real)), mp.mpf(float(dd0.imag)))
                e1 = mp.mpc(mp.mpf(float(dd1.real)), mp.mpf(float(dd1.imag)))
                if want_deriv:
                    out[a, b] = ((6 * t2 - 6 * t) * d0
                                 + (3 * t2 - 4 * t + 1) * h * e0
                                 + (-6 * t2 + 6 * t) * d1
                                 + (3 * t2 - 2 * t) * h * e1) / h
                else:
                    out[a, b] = ((2 * t3 - 3 * t2 + 1) * d0
                                 + (t3 - 2 * t2 + t) * h * e0
                                 + (-2 * t3 + 3 * t2) * d1
                                 + (t3 - t2) * h * e1)
        return out

    def GS_dK(self, x: mp.mpf):
        """G, A = G' - S, B = M + S'/2, K  (matches coeffs())."""
        G = self._entry(x, "G", False)
        S = self._entry(x, "S", False)
        Gp = self._entry(x, "G", True)
        Sp = self._entry(x, "S", True)
        A = Gp - S
        B = self._entry(x, "M", False) + 0.5 * Sp
        K = self._entry(x, "K", False)
        return G, A, B, K


# ==========================================================================
# mpmath linear algebra helpers
# ==========================================================================
def mp_norm_cols(Y, dY):
    """Renormalize each stacked column to unit 2-norm; returns scale
    factors (recorded analytically)."""
    scales = []
    for c in range(3):
        n2 = mp.mpc(0)
        for i in range(3):
            n2 += Y[i, c] * mp.conj(Y[i, c])
            n2 += dY[i, c] * mp.conj(dY[i, c])
        nrm = mp.sqrt(mp.re(n2))
        if nrm == 0:
            nrm = mp.mpf(1)
        scales.append(nrm)
        for i in range(3):
            Y[i, c] /= nrm
            dY[i, c] /= nrm
    return scales


def rhs_step(G, A, B, K, Y, dY, w2, h):
    """One RK4 step of Y'' = -G^-1 (A Y' + (B - w2 K) Y)."""
    def f(Yv, dYv):
        return dYv, -mp.lu_solve(G, A @ dYv + (B - w2 * K) @ Yv)
    k1y, k1d = f(Y, dY)
    k2y, k2d = f(Y + h / 2 * k1y, dY + h / 2 * k1d)
    k3y, k3d = f(Y + h / 2 * k2y, dY + h / 2 * k2d)
    k4y, k4d = f(Y + h * k3y, dY + h * k3d)
    return (Y + h / 6 * (k1y + 2 * k2y + 2 * k3y + k4y),
            dY + h / 6 * (k1d + 2 * k2d + 2 * k3d + k4d))


def propagate(coeffs: MPCoeffs, Y, dY, r0: float, r1: float,
              n_steps: int, w2):
    r = mp.mpf(repr(r0))
    h = (mp.mpf(repr(r1)) - r) / n_steps
    log_scales = mp.mpf(0)
    for s in range(n_steps):
        G, A, B, K = coeffs.GS_dK(r)
        Y, dY = rhs_step(G, A, B, K, Y, dY, w2, h)
        r += h
        if (s + 1) % RENORM_EVERY == 0 and s + 1 < n_steps:
            sc = mp_norm_cols(Y, dY)
            for v in sc:
                log_scales += mp.log(v)
    return Y, dY, log_scales


def pencil_companion(w2, G, B, K):
    """C = [[0, I], [G^-1 (w2 K - B), 0]]  (same linearization as
    LocalPencil.matrices)."""
    Z = mp.zeros(3, 3)
    I3 = mp.eye(3)
    Low = mp.lu_solve(G, w2 * K - B)
    C = mp.matrix(6, 6)
    for i in range(3):
        for j in range(3):
            C[i, j + 3] = 1
            C[i + 3, j] = Low[i, j]
    return C


def inverse_iteration(M, mu0: mp.mpc, iters: int = 4):
    """Refine eigenpair (mu, v) of M from a double-precision seed via
    shifted inverse iteration + Rayleigh quotient."""
    n = M.rows
    v = mp.matrix(n, 1)
    for i in range(n):
        v[i, 0] = mp.mpf(i + 1).sqrt()      # deterministic start
    v = v / mp.sqrt(sum(v[i, 0] * mp.conj(v[i, 0]) for i in range(n)))
    mu = mu0
    for _ in range(iters):
        A = mp.matrix(M)
        for i in range(n):
            A[i, i] -= mu
        try:
            w = mp.lu_solve(A, v)
        except ZeroDivisionError:
            break
        wn = mp.sqrt(sum(w[i, 0] * mp.conj(w[i, 0]) for i in range(n)))
        v = w / wn
        Mv = M @ v
        ray = sum(v[i, 0] * mp.conj(Mv[i, 0]) for i in range(n))
        mu = ray
    return mu, v


def incoming_exterior_seed(coeffs: MPCoeffs, w: mp.mpc, r_seed: float):
    """3 incoming branches at r_seed: refined (mu, Psi, Psi') in mpmath."""
    G, A, B, K = coeffs.GS_dK(mp.mpf(repr(r_seed)))
    w2 = w * w
    C = pencil_companion(w2, G, B, K)

    # double-precision seeds from the exact same formula as LocalPencil
    Gd, Ad, Bd, Kd = (np.asarray(x.tolist(), dtype=complex) if hasattr(x, "tolist") else x
                      for x in (G, A, B, K))
    Gd = np.array(G.values, dtype=complex).reshape(3, 3) if hasattr(G, "values") else np.zeros((3, 3))
    # safer: build the doubles from the numpy tables directly via the
    # double solver's own classes (16-digit seeds are only guesses)
    from transport_bridge.stable_evans_v3 import LocalPencil
    op = coeffs._op
    lp = LocalPencil(op)
    mu_d, alphas_d, V_d = lp.branches(w, r_seed)
    idx = list(np.argsort(mu_d.real)[:3])

    out = []
    R = mp.mpf(repr(r_seed))
    for j in idx:
        mu_ref, vec6 = inverse_iteration(C, mp.mpc(str(mu_d[j].real),
                                                   str(mu_d[j].imag)))
        Y0 = mp.matrix(3, 1)
        for i in range(3):
            Y0[i, 0] = vec6[i, 0]
        # alpha from the same balance formula as LocalPencil.branches
        L = mu_ref * mu_ref * G + B - w2 * K
        LH = L.H
        u = inverse_iteration(LH, mp.mpc(str(alphas_d[j].real) * 0 + str(
            (mu_ref * mu_ref).real), str(0)), iters=3)[1] if False else None
        # left eigenvector of L for its smallest-|eigenvalue|: inverse
        # iteration on L^H with shift 0 (L is nonsingular generically;
        # the smallest eigenvalue dominates the iterated vector)
        u = mp.matrix(3, 1)
        for i in range(3):
            u[i, 0] = mp.mpf(i + 2).sqrt()
        u = u / mp.sqrt(sum(u[i, 0] * mp.conj(u[i, 0]) for i in range(3)))
        for _ in range(6):
            u = mp.lu_solve(LH, u)
            un = mp.sqrt(sum(u[i, 0] * mp.conj(u[i, 0]) for i in range(3)))
            u = u / un
        a0 = Y0
        Rj = mu_ref * G + Gp_ = None  # placeholder removed below
        return out
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("region", choices=list(PROBES))
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--refine", type=int, default=1,
                    help="RK4 step multiplier (1 = base grid)")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    raise SystemExit("NOT EXECUTABLE — see b6_mpmath_evans_v2.py")


if __name__ == "__main__":
    raise SystemExit(main())
