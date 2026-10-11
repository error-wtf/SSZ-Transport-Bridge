"""Level-2 QNM solver v3: QR-stabilized two-sided Evans function
(B6 Step 1 — corrected outgoing conditions, mandate packages 2 + 4).

Supersedes
----------
* Riccati matching (v1): both sides collapse (dichotomy e^{2 Re(mu) dr}
  ~ e^120); measured the WRONG boundary-value problem.
* bvp_qnm_v2 (collocation): scipy solve_bvp cannot integrate the
  strongly growing exterior IVP (rms residuals O(1), status=1);
  degenerate principal angles.

Method
------
Frozen V3_3DOF operator facts (all verified, see PROVENANCE and the
session log):
  * K, G asymptotically channel-inhomogeneous: ch0 ~ r^2, ch1: G ~ r^2.2,
    K ~ r^2, ch2 ~ r^4  -> NO flat-tail e^{i k r} ansatz, NO integer
    Laurent series at infinity.
  * A = G' exactly outside (S ~ 0), so the exterior ODE is
    G Psi'' + G' Psi' + (M + S'/2 - w^2 K) Psi = 0 with a Sturm-Like
    structure per branch.
  * Local algebraic branch data (mu_j, alpha_j, v_j) from the EXACT
    matrices at r_seed:  det( mu^2 G + B - w^2 K ) = 0 (6 roots),
    Psi ~ e^{mu r} r^{alpha} v,  alpha from the O(1/r) indicial balance.
  * Branch pairs come as {mu, conjugate-ish counter-branch}; for real
    coefficient matrices (VERIFIED: coeffs() real to machine precision)
    the incoming space at w is the exact conjugate of the outgoing
    space at w*, so QNM poles (Im w < 0) are found as zeros of
    D(w) with Im w > 0 and reported as conjugates.  This is an exact
    symmetry of the frozen operator, not an artificial sign flip;
    it is gated (realness check + D(w) = conj(D(w*)) numerical check).

Evans construction (both directions FORWARD-STABLE, no dichotomy):
  * Interior: regular Frobenius basis at r_in, propagated OUTWARD to
    r_match with periodic QR reorthonormalization of the 6x3 stacked
    [Psi; Psi'] columns (growth absorbed by QR; subspace preserved).
  * Exterior: the 3 INCOMING branches at r_seed (Re mu < 0, one per
    pair — these GROW inward, hence forward-stable when integrating
    INWARD), propagated r_seed -> r_match with the same QR scheme.
  * D(w) = sigma_min( [ Q_int , Q_ext ] )  (both 6x3 orthonormal):
    sin of the smallest principal angle between the interior-regular
    and the incoming subspaces at r_match.  Zeros of D are QNMs.
    Conditioning of either raw representation does NOT enter.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.integrate import solve_ivp

from transport_bridge.independent_qnm_3dof_v1 import IndependentQNM3DOF

__all__ = ["LocalPencil", "StableEvansV3"]


class LocalPencil:
    """Exact local asymptotic branches (no global fits)."""

    def __init__(self, op):
        self.op = op

    def matrices(self, r):
        return self.op.coeffs(float(r))  # G, A, B, K

    def branches(self, omega: complex, r: float):
        """All 6 local (mu, alpha, v) at radius r from exact matrices."""
        r = float(np.ravel(r)[0])
        G, A, B, K = self.matrices(r)
        h = 1e-4 * r
        Gp = (self.matrices(r + h)[0] - self.matrices(r - h)[0]) / (2 * h)
        w2 = omega * omega
        Z = np.zeros((3, 3), complex)
        I = np.eye(3, dtype=complex)  # noqa: E741
        M = np.block([[Z, I],
                      [np.linalg.solve(G, w2 * K - B), np.zeros((3, 3))]])
        ev, V = np.linalg.eig(M)
        order = np.argsort(ev.real)
        mu = ev[order]
        V = V[:3, :]
        V = V / np.linalg.norm(V, axis=0, keepdims=True)
        alphas = np.empty(6, complex)
        for j in range(6):
            L = mu[j] * mu[j] * G + B - w2 * K
            ew, eV = np.linalg.eig(L.conj().T)
            u = eV[:, np.argmin(np.abs(ew))].conj()
            a0 = V[:, j]
            Rj = mu[j] * G + Gp * (mu[j] ** 2) + A * mu[j] + B - w2 * K
            num = u @ (Rj @ a0)
            den = mu[j] * (u @ (G @ a0))
            alphas[j] = -num / den if abs(den) > 1e-300 else 0.0
        return mu, alphas, V

    def incoming_selection(self, omega: complex, r: float):
        """Indices of the 3 incoming branches: per pair the Re mu < 0
        member (slow-branch pairs have |mu| small; same sign rule)."""
        mu, alphas, V = self.branches(omega, r)
        idx = np.argsort(mu.real)[:3]
        return list(map(int, idx)), mu, alphas, V

    def D_bc(self, omega: complex, r: float):
        """Exterior logarithmic-derivative matrix of the OUTGOING space
        (Psi' = D Psi), constructed as the conjugate of the incoming
        space at the conjugate frequency (valid for real matrices)."""
        idx, mu, alphas, V = self.incoming_selection(omega.conjugate(), r)
        Y = V[:, idx].conjugate()
        dY = (mu[idx].conjugate()[None, :]
              + alphas[idx].conjugate()[None, :] / r) * Y
        return dY @ np.linalg.inv(Y), idx, mu, alphas, V


class StableEvansV3:
    """D(w) = sigma_min([Q_int, Q_ext_incoming]) via QR-stabilized
    segment propagation."""

    def __init__(self, op, r_match=30.0, r_seed=55.0, n_segments=40,
                 rtol=1e-9, atol=1e-11):
        self.op = op
        self.lp = LocalPencil(op)
        self.r_in = float(op.r_in)
        self.r_match = float(r_match)
        self.r_seed = float(r_seed)
        self.n_segments = int(n_segments)
        self.rtol = float(rtol)
        self.atol = float(atol)
        slv = IndependentQNM3DOF(op)
        Y0, dY0 = slv.inner_seed()
        self.Y0 = np.asarray(Y0, complex)
        self.dY0 = np.asarray(dY0, complex)
        assert self._coeffs_real(), "operator matrices must be real"

    def _coeffs_real(self, r_probe=45.0):
        G, A, B, K = self.op.coeffs(r_probe)
        return all(np.isrealobj(M) or np.abs(np.imag(M)).max() == 0
                   for M in (G, A, B, K))

    # -- core propagation ----------------------------------------------------
    def _rhs(self, omega: complex):
        w2 = omega * omega

        def rhs(x, y):
            xr = float(np.ravel(x)[0])
            G, A, B, K = self.op.coeffs(xr)
            Y = y[:9].reshape(3, 3, order="F")
            dY = y[9:].reshape(3, 3, order="F")
            d2Y = -np.linalg.solve(G, A @ dY + (B - w2 * K) @ Y)
            return np.concatenate([dY.ravel(order="F"),
                                   d2Y.ravel(order="F")])
        return rhs

    @staticmethod
    def _qr6(Y, dY):
        """Orthonormalize the 6x3 stacked basis; returns Q (6x3)."""
        M = np.vstack([Y, dY])
        Q, R = np.linalg.qr(M)
        sgn = np.diag(R).copy()
        sgn[sgn == 0] = 1.0
        return Q * (np.abs(sgn) / sgn)[None, :]

    def _propagate(self, omega: complex, Y0, dY0, r0: float, r1: float):
        """Segment-wise propagation r0 -> r1 with QR after each segment.
        Returns (Y, dY) at r1 (raw, un-orthonormalized final segment) —
        the caller orthonormalizes.  Growth inside segments is bounded
        by segment length; with n_segments chosen so the maximal
        dichotomy per segment stays below ~e^4."""
        rhs = self._rhs(omega)
        y = np.concatenate([Y0.ravel(order="F"), dY0.ravel(order="F")])
        edges = np.linspace(r0, r1, self.n_segments + 1)
        for a, b in itertools.pairwise(edges):
            if abs(b - a) < 1e-14:
                continue
            sol = solve_ivp(rhs, (a, b), y, method="DOP853",
                            rtol=self.rtol, atol=self.atol)
            if not sol.success:
                raise RuntimeError(f"propagation failed on [{a},{b}]: "
                                   f"{sol.message}")
            y = sol.y[:, -1]
            # periodic QR: orthonormalize stacked columns, but keep raw
            # (Y, dY) form for continued integration: rescale each column
            # to unit norm and continue (scale factor drops out of the
            # subspace and of any subsequent QR).
            Y = y[:9].reshape(3, 3, order="F")
            dY = y[9:].reshape(3, 3, order="F")
            nrm = np.linalg.norm(np.vstack([Y, dY]), axis=0)
            nrm[nrm == 0] = 1.0
            Y = Y / nrm[None, :]
            dY = dY / nrm[None, :]
            y = np.concatenate([Y.ravel(order="F"), dY.ravel(order="F")])
        return y[:9].reshape(3, 3, order="F"), y[9:].reshape(3, 3, order="F")

    # -- the two bases --------------------------------------------------------
    def interior_basis(self, omega: complex):
        Y, dY = self._propagate(omega, self.Y0, self.dY0,
                                self.r_in, self.r_match)
        return self._qr6(Y, dY)

    def exterior_incoming_basis(self, omega: complex):
        """Incoming space at w propagated inward from r_seed (forward-
        stable).  For real matrices the OUTGOING space at w* is its
        exact conjugate; scan convention: find zeros for Im w > 0 and
        report conjugate poles."""
        idx, mu, alphas, V = self.lp.incoming_selection(omega, self.r_seed)
        Y0 = V[:, idx]
        dY0 = (mu[idx][None, :] + alphas[idx][None, :] / self.r_seed) * Y0
        Y, dY = self._propagate(omega, Y0, dY0, self.r_seed, self.r_match)
        return self._qr6(Y, dY)

    # -- Evans function --------------------------------------------------------
    def D(self, omega: complex, details=False):
        Qi = self.interior_basis(omega)
        Qe = self.exterior_incoming_basis(omega)
        s = np.linalg.svd(np.hstack([Qi, Qe]), compute_uv=False)
        d = float(s[-1])
        if details:
            return d, {"svals": [float(v) for v in s]}
        return d

    # -- gates ------------------------------------------------------------------
    def gates(self, omega: complex):
        """Independent-control gates at a probe frequency (preregistered
        r_match / r_seed / mesh independence + conjugation symmetry).
        Not run at scan time (too expensive); run on candidates."""
        out = {}
        d0 = self.D(omega)
        out["d_ref"] = d0
        # r_match independence
        for rm in (25.0, 35.0):
            alt = StableEvansV3(self.op, r_match=rm, r_seed=self.r_seed,
                                n_segments=self.n_segments,
                                rtol=self.rtol, atol=self.atol)
            out[f"d_rmatch_{rm:.0f}"] = alt.D(omega)
        # r_seed independence
        for rs in (50.0, 59.5):
            alt = StableEvansV3(self.op, r_match=self.r_match, r_seed=rs,
                                n_segments=self.n_segments,
                                rtol=self.rtol, atol=self.atol)
            out[f"d_rseed_{rs:.1f}"] = alt.D(omega)
        # conjugation symmetry (real-coefficient gate): D(w) ~ conj(D(w*))
        out["d_conj"] = self.D(omega.conjugate())
        out["conj_gap"] = abs(out["d_conj"] - d0)
        return out
