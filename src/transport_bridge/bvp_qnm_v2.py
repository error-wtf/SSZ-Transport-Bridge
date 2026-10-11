"""Level-2 QNM solver v2: collocation BVP + principal-angle Evans
function (corrected outgoing conditions, B6 Step 1).

Method
------
The frozen V3_3DOF operator has channel-inhomogeneous asymptotics
(G,K ~ r^2 for ch0; G ~ r^2.2 / K ~ r^2 for ch1; G,K ~ r^4 for ch2;
A = G' exactly, S ~ 0 outside).  Non-integer powers rule out an
integer-Laurent Frobenius series at infinity, and both Riccati and
direct basis integration collapse at the matching radius
(e^{2 Re(mu) (R - r_m)} ~ e^120 dichotomy).

v2 avoids integration pathologies entirely:

  * INTERIOR basis: 3 solutions from the regular Frobenius center
    seed, propagated as collocation solves (solve_bvp with boundary
    data only at r_in — global collocation has no growth-direction
    instability), orthonormalized at r_match:  Q_int (6x3).
  * EXTERIOR basis: 3 solutions seeded at r_seed with the LOCAL
    algebraic outgoing pencil data at r_seed —
        Psi(r_seed) = v_j,  Psi'(r_seed) = (mu_j + alpha_j/r_seed) v_j
    with (mu_j, alpha_j, v_j) from the EXACT coeff matrices at
    r_seed (no series fits) — collocated outward to r_match,
    orthonormalized:  Q_ext (6x3).
  * Evans function:  D(omega) = sigma_min([Q_int, Q_ext])
    = sin of the smallest principal angle between the two 3-dim
    subspaces in the 6-dim (Psi, Psi') space at r_match.
    Coordinate conditioning of either basis does NOT enter.
    D -> 0 at QNM frequencies.

Local algebraic pencil at radius r (exact matrices):
  L0(mu) = mu^2 G(r) + B(r) - w^2 K(r);  roots of det L0 = 0 via the
  6x6 companion;  alpha_j from the O(1/r) indicial balance with the
  exact G'(r) (A = G' outside):
    u_j^* [ mu G + G' mu^2 + A mu + B - w^2 K ] a0
      + mu alpha u_j^* G a0 = 0.
Outgoing selection: oscillatory branches (|Re mu| < 3 |Im mu|):
  Re mu > 0;  evanescent branches: Re mu < 0.  Sign-structured, no
  magnitude ranking, no artificial sign flips.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp

__all__ = ["LocalPencil", "PrincipalAngleEvans"]


class LocalPencil:
    """Exact local asymptotic branches of the frozen operator."""

    def __init__(self, op):
        self.op = op

    def matrices(self, r):
        return self.op.coeffs(float(r))  # G, A, B, K

    def branches(self, omega: complex, r: float, dlog=False):
        """All 6 local (mu, alpha, v) at radius r from exact matrices."""
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

    def outgoing(self, omega: complex, r: float):
        mu, alphas, V = self.branches(omega, r)
        osc = np.abs(mu.real) < 3.0 * np.abs(mu.imag)
        mask = np.where(osc, mu.real > 0, mu.real < 0)
        return mask, mu, alphas, V

    def D_bc(self, omega: complex, r: float):
        """Exterior logarithmic-derivative matrix D_r: Psi'(r) = D Psi(r)
        for the 3-dim outgoing subspace."""
        mask, mu, alphas, V = self.outgoing(omega, r)
        assert int(mask.sum()) == 3
        idx = np.where(mask)[0]
        Y = V[:, idx]
        dY = (mu[idx][None, :] + alphas[idx][None, :] / r) * Y
        return dY @ np.linalg.inv(Y), mask, mu, alphas, V


class PrincipalAngleEvans:
    """D(omega) = sigma_min([Q_int, Q_ext]) via global collocation."""

    def __init__(self, op, r_in=None, r_match=30.0, r_seed=None,
                 n_mesh=260, rtol=1e-8, atol=1e-10):
        self.op = op
        self.lp = LocalPencil(op)
        self.r_in = float(op.r_in if r_in is None else r_in)
        self.r_match = float(r_match)
        self.r_seed = float(op.r[-1] if r_seed is None else r_seed)
        self.n_mesh = int(n_mesh)
        self.rtol = float(rtol)
        self.atol = float(atol)
        # frobenius center seed (regular 3-dim basis at r_in)
        from transport_bridge.independent_qnm_3dof_v1 import IndependentQNM3DOF
        slv = IndependentQNM3DOF(op)
        Y0, dY0 = slv.inner_seed()
        self.Y0 = np.asarray(Y0, complex)
        self.dY0 = np.asarray(dY0, complex)
        self._int_sol_cache = None

    # -- first order system --------------------------------------------------
    def _sys(self, omega: complex):
        w2 = omega * omega

        def rhs(x, y):
            xr = x.item() if isinstance(x, np.ndarray) and x.ndim == 0 else (
                float(np.ravel(x)[0]))
            G, A, B, K = self.op.coeffs(float(xr))
            Y = y[:3]
            dY = y[3:]
            return np.concatenate([dY, -np.linalg.solve(G, A @ dY
                                   + (B - w2 * K) @ Y)])
        return rhs

    def _collocate(self, omega: complex, seedY, seeddY, r_bc, r_other):
        """Collocate the 3 solution columns on [min,max] (increasing x);
        the seed data constrains the end at r_bc ('left' or 'right')."""
        n = 3
        at_left = r_bc < r_other
        x = (np.linspace(r_bc, r_other, self.n_mesh) if at_left
             else np.linspace(r_other, r_bc, self.n_mesh))
        sols = []
        f = self._sys(omega)
        for c in range(n):
            y0 = np.zeros((6, x.size), complex)
            y0[:3, :] = seedY[:, c][:, None]
            y0[3:, :] = seeddY[:, c][:, None]
            target = y0[:, 0]
            if at_left:
                def bcc(ya, yb, t=target):
                    return ya - t
            else:
                def bcc(ya, yb, t=target):
                    return yb - t
            sol = solve_bvp(f, bcc, x, y0, tol=self.rtol,
                            max_nodes=self.n_mesh * 6)
            sols.append(sol)
        return sols

    def interior_basis(self, omega: complex):
        """(Y, dY) of the 3-dim regular interior basis at r_match."""
        sols = self._collocate(omega, self.Y0, self.dY0,
                               self.r_in, self.r_match)
        Y = np.column_stack([s.y[:3, -1] for s in sols])
        dY = np.column_stack([s.y[3:, -1] for s in sols])
        return Y, dY

    def exterior_basis(self, omega: complex):
        """(Y, dY) of the 3-dim outgoing exterior basis at r_match."""
        mask, mu, alphas, V = self.lp.outgoing(omega, self.r_seed)
        idx = np.where(mask)[0]
        Y0e = V[:, idx]
        dY0e = (mu[idx][None, :] + alphas[idx][None, :] / self.r_seed) * Y0e
        sols = self._collocate(omega, Y0e, dY0e,
                               self.r_seed, self.r_match)
        Y = np.column_stack([s.y[:3, -1] for s in sols])
        dY = np.column_stack([s.y[3:, -1] for s in sols])
        return Y, dY

    @staticmethod
    def _orth(self, Y, dY):
        """Orthonormalize the 3-dim subspace in C^6."""
        M = np.vstack([Y, dY])
        Q, R = np.linalg.qr(M)
        sgn = np.diag(R).copy()
        sgn[sgn == 0] = 1.0
        return Q * (np.abs(sgn) / sgn)[None, :]

    def D(self, omega: complex, details=False):
        Qi = self._orth(self, *self.interior_basis(omega))
        Qe = self._orth(self, *self.exterior_basis(omega))
        s = np.linalg.svd(np.hstack([Qi, Qe]), compute_uv=False)
        d = float(s[-1])
        if details:
            return d, {"svals": s.tolist()}
        return d

    def solver_health(self, omega: complex):
        """Collocation residual/status per solve (gate diagnostics)."""
        out = {}
        sols_int = self._collocate(omega, self.Y0, self.dY0,
                                   self.r_in, self.r_match)
        out["int_status"] = [s.status for s in sols_int]
        out["int_rms_res"] = [float(s.rms_residuals.max()) for s in sols_int]
        mask, mu, al, V = self.lp.outgoing(omega, self.r_seed)
        idx = np.where(mask)[0]
        Y0e = V[:, idx]
        dY0e = (mu[idx][None, :] + al[idx][None, :] / self.r_seed) * Y0e
        sols_ext = self._collocate(omega, Y0e, dY0e,
                                   self.r_seed, self.r_match)
        out["ext_status"] = [s.status for s in sols_ext]
        out["ext_rms_res"] = [float(s.rms_residuals.max()) for s in sols_ext]
        out["out_mask"] = mask.astype(int).tolist()
        return out
