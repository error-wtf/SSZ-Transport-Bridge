"""Asymptotic tail basis for the coupled 3-DOF SSZ resonance problem
(B6 Step 1, Level 2 — corrected outgoing conditions, mandate packages
2 and 4).

Why this module exists
----------------------
The frozen V3_3DOF operator has K, G ~ r^2 * (constant matrices) as
r -> infinity (verified: 7-term Laurent fits of coeffs() over r in
[40, 60] with relative residuals ~ 1e-5..1e-6).  Consequently:

  * a flat-tail ansatz Psi_j ~ e^{i k r} with CONSTANT k (the old
    seeded outgoing condition in independent_qnm_3dof_v1) is NOT an
    asymptotic solution of the actual ODE;
  * inward Riccati integration of the outgoing side collapses onto the
    incoming conjugate branch (all outgoing acoustic branches GROW
    outward for Im w < 0; growth contrast e^{2 Re(mu) dr} ~ e^120
    across r_out - r_match = 30);
  * the dphi channel carries a LONG-RANGE 1/r piece (B_1[1,1] ~ 27),
    so the tail series MUST keep the full Laurent tower, not a
    constant-coefficient truncation.

Correct construction — Frobenius at infinity:
  Psi_j = e^{mu_j r} r^{alpha_j} sum_{n>=0} a_n^(j) r^{-n}.

Leading balances (with G = sum_m G_m r^{2-m}, A likewise with A_0 = 0
since A ~ r, B likewise with B_0 = B_1 = 0 since B ~ r^0, K likewise):

  O(r^2):  L0(mu) a0 = 0,  L0 = mu^2 G0 + B0 - w^2 K0   (6 roots mu_j)
  O(r^1):  solvability with left null vector u of L0(mu_j):
           u^* [ mu(alpha+1) G0 + G1 mu^2 + A1 mu + B1 - w^2 K1 ] a0 = 0
           -> exact indicial equation for alpha_j
  O(r^{2-n}), n >= 2:  L0 a_n = -rhs(a_{<n}), solved with the
           orthogonality gauge a_n perp a0 (removes the 1-dim null
           space of L0, rank 2).

Truncation error of the series at n_max: relative ~ (r_ref/r)^{n_max}
at radius r_ref <= r_fit range; validated numerically in selftest()
against direct ODE residuals.

Conventions (frozen V3_3DOF, PROVENANCE_V2_B6_STEP1.json):
  ODE:  G Psi'' + A Psi' + (B - w^2 K) Psi = 0,
  A = G' - S, B = M + S'/2 (as exported by Closure).
"""

from __future__ import annotations

import numpy as np

__all__ = ["LaurentTail"]


class LaurentTail:
    """Laurent-tail model of the frozen operator + outgoing basis."""

    def __init__(self, op, r_fit=(40.0, 60.0), n_terms=10, n_samples=31):
        self.op = op
        self.r1, self.r2 = float(r_fit[0]), float(r_fit[1])
        self.n_terms = int(n_terms)
        rs = np.linspace(self.r1, self.r2, int(n_samples))
        self.rs = rs
        # uniform exponent grid r^{2-m}; vanishing leading coefficients
        # (A_0 = 0, B_0 = B_1 = 0) absorb the term-order offsets.
        P = np.vstack([rs ** (2.0 - m) for m in range(self.n_terms)]).T
        self.C = {}
        self.resid = {}
        for name, idx in (("G", 0), ("A", 1), ("B", 2), ("K", 3)):
            Y = np.array([op.coeffs(float(r))[idx] for r in rs]).reshape(len(rs), 9)
            X, *_ = np.linalg.lstsq(P, Y, rcond=None)
            self.C[name] = X.reshape(self.n_terms, 3, 3)
            pred = (P @ X).reshape(len(rs), 3, 3)
            scale = np.abs(Y.reshape(len(rs), 3, 3)).max()
            self.resid[name] = float(np.abs(pred - Y.reshape(len(rs), 3, 3)).max() / scale)

    # -- leading-order data -------------------------------------------------
    @property
    def G0(self):
        return self.C["G"][0]

    def L0(self, mu: complex, omega: complex):
        """Leading O(r^2) pencil (rank 2 on each branch)."""
        return (mu * mu * self.G0 + self.C["B"][0]
                - omega * omega * self.C["K"][0])

    def mu_branches(self, omega: complex):
        """All 6 asymptotic modal exponents: det(mu^2 G0 + B0 - w^2 K0) = 0.

        Returns (mu sorted by Re, V with unit branch vectors a0 = v_j).
        """
        B0, K0 = self.C["B"][0], self.C["K"][0]
        Z = np.zeros((3, 3), complex)
        I = np.eye(3, dtype=complex)
        M = np.block([[Z, I],
                      [np.linalg.solve(self.G0, omega * omega * K0 - B0),
                       np.zeros((3, 3), complex)]])
        ev, V = np.linalg.eig(M)
        order = np.argsort(ev.real)
        mu = ev[order]
        V = V[:, order]
        V = V[:3, :]  # companion vectors are [Psi; mu Psi]: keep Psi block
        V = V / np.linalg.norm(V, axis=0, keepdims=True)
        return mu, V

    def alpha_branches(self, omega: complex):
        """Exact O(r^1) indicial equation per branch.

        Solvability of  L0(mu_j) a1 = -[mu(alpha+1) G0 + G1 mu^2 + A1 mu
        + B1 - w^2 K1] a0  with u_j the left null vector of L0(mu_j):
          mu (u^* G0 a0) alpha = -u^* [mu G0 + G1 mu^2 + A1 mu + B1
                                       - w^2 K1] a0.
        Returns (mu, alpha, V).
        """
        mu, V = self.mu_branches(omega)
        G0, G1 = self.G0, self.C["G"][1]
        A1 = self.C["A"][1]
        B1, K1 = self.C["B"][1], self.C["K"][1]
        alphas = np.empty(6, complex)
        for j in range(6):
            L = self.L0(mu[j], omega)
            # left null vector: eigenvector of L^* with smallest |eig|
            ew, eV = np.linalg.eig(L.conj().T)
            u = eV[:, np.argmin(np.abs(ew))].conj()
            a0 = V[:, j]
            Rj = (mu[j] * G0 + G1 * (mu[j] ** 2) + A1 * mu[j]
                  + B1 - (omega * omega) * K1)
            num = float(np.real_if_close(0)) if False else (u @ (Rj @ a0))
            den = mu[j] * (u @ (G0 @ a0))
            alphas[j] = -num / den if abs(den) > 1e-300 else 0.0
        return mu, alphas, V

    # -- outgoing selection (radiation condition) ---------------------------
    def outgoing_mask(self, omega: complex):
        """Physical radiation condition per branch.

        Oscillatory acoustic branches (|Re mu| < 3|Im mu|): OUTGOING
        iff Re(mu) > 0 (outward phase velocity + outward growth for
        Im w < 0 — standard QNM outgoing).  Evanescent/overdamped
        branches: physical tail is the DECAYING one, Re(mu) < 0.
        The 6 branches come in +/- pairs, so the selection is sign-
        structured, not magnitude-ranked (no path memory, no artificial
        sign flips of k without an eigenpair — mandate item).
        """
        mu, alphas, V = self.alpha_branches(omega)
        osc = np.abs(mu.real) < 3.0 * np.abs(mu.imag)
        mask = np.where(osc, mu.real > 0, mu.real < 0)
        return mask, mu, alphas, V

    # -- Frobenius series at infinity ---------------------------------------
    def series(self, omega: complex, j: int, n_max: int = 24):
        """Downward recursion for a_n of branch j (see module docstring).

        For n >= 2:
          L0 a_n = -rhs_n,  rhs_n = sum over m>=1 and shifts of
          a_{n'}, a_{n'-1}, a_{n'-2} (explicit below), solved with the
          gauge <a_n, a0> = 0 via the bordered system
          [ L0  c ; a0^* 0 ] [a_n; g] = [-rhs; 0],  c = mu G0 a0...
        Actually the null direction of L0 is spanned by a0 itself
        (Frobenius structure), so the bordered system below is exact.

        Coefficient bookkeeping (coefficient of r^{2-n}):
          d_n = mu^2 a_n + mu(al-n+1) a_{n-1} + (al-n+2)(al-n+1) a_{n-2}
          c_n = mu a_n + (al-n+1) a_{n-1}
          0 = sum_m [ G_m d_{n-m} + A_m c_{n-m} + E_m a_{n-m} ],
              E = B - w^2 K   (same exponent grid).
        a_n unknown appears with m = 0 only (d_n, c_n have no r^{-m});
        everything else is rhs from n' < n.
        """
        Gc, Ac, Ec = self.C["G"], self.C["A"], self.C["B"] - omega**2 * self.C["K"]
        mu, alphas, V = self.alpha_branches(omega)
        mu_j, al, a0 = mu[j], alphas[j], V[:, j]
        L0 = self.L0(mu_j, omega)
        a = np.zeros((n_max, 3), complex)
        a[0] = a0
        # bordered solve: null(L0) = span(a0) (rank-2 pencil per branch)
        for n in range(1, n_max):
            rhs = np.zeros(3, complex)
            for m in range(self.n_terms):
                # a_{n-m} parts (d, c at level n-m), only for m >= 1
                # (m = 0 terms hold the a_n unknown and a_{n-1} in d_n):
                pass
            # explicit: gather for each p < n the coefficient of a_p in
            # the r^{2-n} equation, then invert analytically:
            # a_p appears via:
            #   m = n - p   : G_m d_p,  A_m c_p,  E_m a_p   (m >= 1)
            #   m = n - p+1 : G_m * [mu(al-p+1) a_p]  in d_{p-1}? no...
            # Cleaner: d and c contain a_p only in d_p and c_p (their
            # own level).  G_m d_p lands at r^{2-p-m}; want 2-n =>
            # m = n-p.  Same for A_m c_p, E_m a_p.  Plus the m=0
            # higher-derivative shift terms:
            #   d_n contains mu(al-n+1) a_{n-1} and
            #   (al-n+2)(al-n+1) a_{n-2};  c_n contains (al-n+1) a_{n-1}.
            #   These land at r^{2-n} with m = 0:  G_0 * [...], A_0 * [...]
            #   (A_0 = 0), E_0 * none.
            rhs = np.zeros(3, complex)
            # shift terms (m = 0):
            s1 = mu_j * (al - n + 1)
            s2 = (al - n + 2) * (al - n + 1)
            if n >= 2:
                rhs = rhs + Gc[0] @ (s1 * a[n - 1] + s2 * a[n - 2])
            else:
                rhs = rhs + Gc[0] @ (s1 * a[n - 1])
            # Laurent terms m >= 1:
            for p in range(n):
                m = n - p
                if m >= self.n_terms:
                    continue
                dp = mu_j * mu_j * a[p]
                if p >= 1:
                    al_p1 = al - p + 1
                    dp = dp + (mu_j * al_p1 * a[p - 1]
                               + (al_p1 + 1) * al_p1 * a[p - 2] if p >= 2
                               else mu_j * al_p1 * a[p - 1])
                cp = mu_j * a[p]
                if p >= 1:
                    cp = cp + (al - p + 1) * a[p - 1]
                rhs = rhs + Gc[m] @ dp + Ac[m] @ cp + Ec[m] @ a[p]
            # solve [L0, a0*muG0a0-col; a0^*, 0] with gauge <a_n, a0> = 0:
            # L0 a_n = -rhs, null(L0) = span(a0) (assumed; verify below)
            g = Gc[0] @ a0
            M = np.zeros((4, 4), complex)
            M[:3, :3] = L0
            M[:3, 3] = g
            M[3, :3] = a0.conj()
            b = np.zeros(4, complex)
            b[:3] = -rhs
            sol = np.linalg.solve(M, b)
            a[n] = sol[:3]
        return a

    def tail_basis(self, omega: complex, r_eval: float, n_max: int = 24):
        """3 outgoing solutions evaluated at r_eval: returns (Y, dY).

        Psi_j(r) = e^{mu_j r} r^{al_j} sum_n a_n^(j) r^{-n};
        columns equilibrated (scale-invariant for the 6x6 matching
        determinant up to a global constant).
        """
        mask, mu, alphas, _ = self.outgoing_mask(omega)
        assert int(mask.sum()) == 3, f"outgoing mask picks {int(mask.sum())} of 6"
        idx = np.where(mask)[0]
        Y = np.empty((3, 3), complex)
        dY = np.empty((3, 3), complex)
        for c, j in enumerate(idx):
            a = self.series(omega, int(j), n_max=n_max)
            # Horner in 1/r
            s = a[n_max - 1].copy()
            for n in range(n_max - 2, -1, -1):
                s = s / r_eval + a[n]
            sds = a[n_max - 1] * (-(n_max - 1))
            for n in range(n_max - 2, -1, -1):
                sds = sds / r_eval + a[n] * (-n)
            psi = np.exp(mu[j] * r_eval) * r_eval ** alphas[j] * s
            dpsi = (mu[j] + alphas[j] / r_eval) * psi + (
                np.exp(mu[j] * r_eval) * r_eval ** alphas[j] * sds / r_eval
            )
            Y[:, c] = psi
            dY[:, c] = dpsi
        nrm = np.linalg.norm(Y, axis=0)
        Y = Y / nrm[None, :]
        dY = dY / nrm[None, :]
        return Y, dY

    # -- validation ----------------------------------------------------------
    def selftest(self, omega_test=1.0 - 0.5j, n_max_list=(8, 16, 24)):
        """ODE-residual convergence of the series at r = r2 (fit edge)
        and r = 30 (matching radius); must fall like truncation order.
        Returns dict with per-(n_max, r) relative residuals and the
        Laurent fit residuals."""
        from scipy.integrate import quad  # noqa: F401  (placeholder import)

        out = {"laurent_resid": dict(self.resid), "series": {}}
        w = complex(omega_test)
        mask, mu, alphas, _ = self.outgoing_mask(w)
        op = self.op
        for n_max in n_max_list:
            for r_eval in (self.r2, 30.0):
                worst = 0.0
                for j in np.where(mask)[0]:
                    a = self.series(w, int(j), n_max=n_max)
                    # build Psi, Psi', Psi'' via term-wise differentiation
                    def fact(nn, der):
                        if der == 0:
                            return 1.0 + 0j
                        if der == 1:
                            return -float(nn)
                        return float(nn * (nn + 1))  # second 1/r deriv
                    # coefficient c_n of r^{-n} with derivative order:
                    # Psi = e^{mu r} r^al sum a_n r^{-n}
                    # term_n = a_n e^{mu r} r^{al-n}
                    # term_n' = (mu + (al-n)/r) term_n
                    # term_n'' = (mu + (al-n)/r)^2 term_n
                    #            - (al-n)/r^2 term_n
                    mu_j, al = mu[j], alphas[j]
                    S0 = S1 = S2 = 0j
                    for n in range(n_max):
                        t = (mu_j + (al - n) / r_eval)
                        w_n = (mu_j - (al - n) / r_eval)  # unused
                        base = a[n] * r_eval ** (-n)
                        S0 = S0 + base
                        S1 = S1 + t * base
                        S2 = S2 + (t * t - (al - n) / r_eval**2) * base
                    e = np.exp(mu_j * r_eval) * r_eval ** al
                    psi, d1, d2 = e * S0, e * S1, e * S2
                    G, A, B, K = op.coeffs(float(r_eval))
                    res = G @ d2 + A @ d1 + (B - w * w * K) @ psi
                    worst = max(worst, float(np.abs(res).max()
                                             / max(np.abs(psi).max(), 1e-300)))
                out["series"][f"n{n_max}_r{r_eval:.0f}"] = worst
        return out
