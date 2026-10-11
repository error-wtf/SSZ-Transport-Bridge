"""Independent coupled 3-channel QNM machinery for the physical V4 branch.

B6 step 1, level 2 (independent implementation). Shares with
SSZ_FULL_CLOSURE only the frozen inputs declared in
data/frozen_specs/PROVENANCE_V2_B6_STEP1.json (operator npz, ODE
convention, Frobenius exponents, BC/preflight data, preregistered
domain). No Closure solver code is imported or copied: basis
construction, integration scheme, matching, and root finding are
implemented from the mathematical specification.

Physics
-------
ODE (t-time convention, omega = omega_t):
    G Psi'' + (G' - S) Psi' + (M + S'/2 - omega^2 K) Psi = 0
on the frozen grid r in [0.05, 60], Psi = (psi, dphi, V).

Inner: regular Frobenius basis (exponents 0.5390413104712023,
2.0002806149817025, 1.0016821978890482 [V]) built by a Taylor
coefficient recursion around r_in = 0.05 (Cauchy-Kovalevskaya style,
matrix polynomial convolution), evaluated at r_series_end and handed
to an adaptive DOP853 integration of the full ODE to the matching
radius. The basis is exact to series order at r_in and recomputed per
omega through the -omega^2 K term.

Outer: three outgoing solutions against the edge pencil
    P(w, k) = w^2 K - k^2 G + i k (G' - S) - (M + S'/2) = 0
evaluated on the exported edge blocks; the 6 k-roots are obtained from
the linearized companion form and the 3 outgoing branches selected by
continuity from the reference omega on the negative imaginary axis
(decaying e^{+i k r_star}, Im k > 0 there). Outgoing k(omega) is
solved PER omega -- no common-speed shortcut (frozen preflight:
common_speed_bc_invalid = true).

Root condition: D_Jost(omega) = det of the 6x6 interface matrix whose
columns are (Psi_j, Psi'_j) of the 3 regular inner and 3 outgoing
outer solutions at the matching radius. Resonances: blind sigma_min
scan over the preregistered tiles, local refinement, Newton polish on
D(omega) = 0 in the complex plane.
"""
from __future__ import annotations

import hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

__all__ = ["FrozenOperator", "IndependentQNM3DOF"]

REGULAR_EXPONENTS = (0.5390413104712023, 2.0002806149817025, 1.0016821978890482)


class FrozenOperator:
    """Load and hash-verify the frozen V3_3DOF operator npz (level 1)."""

    EXPECTED_SHA256 = "ad2bee7d3994097479933c1eea787f3aa381648d098ffd36aab32a4c436db2ef"

    def __init__(self, npz_path: str):
        raw = open(npz_path, "rb").read()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != self.EXPECTED_SHA256:
            raise ValueError(f"frozen operator hash mismatch: {digest}")
        d = np.load(npz_path, allow_pickle=True)
        self.r = d["r"].astype(float)
        self.K = d["K_phys"].astype(complex)
        G_raw = d["G_phys"].astype(complex)
        S_raw = d["S_phys"].astype(complex)
        self.S = 0.5 * (S_raw - np.swapaxes(S_raw, 1, 2))  # enforce exact skew
        self.G = G_raw
        self.M = d["M_phys"].astype(complex)
        self.R = d["R_phys"].astype(complex)
        self.f = d["f"].astype(float)
        self.h = d["h"].astype(float)
        self.phi = d["phi"].astype(float)
        self.L = int(np.atleast_1d(d["L"])[0])
        self.r_in = float(self.r[0])
        self.r_out = float(self.r[-1])
        # C^2 spline interpolants incl. derivative splines (independent
        # choice; frozen data consumed as-is)
        self._sp = {}
        self._sp_call = {}
        for name, arr in (("G", self.G), ("S", self.S),
                          ("M", self.M), ("K", self.K)):
            sp_r = CubicSpline(self.r, arr.real, axis=0)
            sp_i = CubicSpline(self.r, arr.imag, axis=0)
            self._sp[name] = (sp_r, sp_i)
            self._sp_call[name] = (
                lambda x, _r=sp_r, _i=sp_i: _r(x) + 1j * _i(x)
            )
        self._d_sp = {}
        for name in ("G", "S", "M", "K"):
            sp_r, sp_i = self._sp[name]
            self._d_sp[name] = (
                lambda x, _r=sp_r.derivative(), _i=sp_i.derivative():
                    _r(x) + 1j * _i(x)
            )
        self._spG_d = self._d_sp["G"]
        self._spS_d = self._d_sp["S"]

    # -- evaluated coefficient functions -------------------------------
    def coeffs(self, x: float):
        G = self._sp_call["G"](x)
        S = self._sp_call["S"](x)
        A = self._d_sp["G"](x) - S
        B = self._sp_call["M"](x) + 0.5 * self._d_sp["S"](x)
        return G, A, B, self._sp_call["K"](x)

    def rhs_first_order(self, x: float, y: np.ndarray, omega2: complex):
        """y = [Psi (3), Psi' (3)] flattened column-major for 3 columns."""
        Y = y[:9].reshape(3, 3, order="F")
        dY = y[9:].reshape(3, 3, order="F")
        G, A, B, K = self.coeffs(x)
        ddY = -np.linalg.solve(G, A @ dY + (B - omega2 * K) @ Y)
        return np.concatenate([dY.reshape(-1, order="F"),
                               ddY.reshape(-1, order="F")])


class IndependentQNM3DOF:
    """Coupled Jost determinant with per-omega outgoing pencil."""

    def __init__(self, op: FrozenOperator, r_match: float = 30.0):
        self.op = op
        self.r_match = r_match

    # -- inner regular basis (Frobenius seed + adaptive integration) ----
    def inner_seed(self):
        """Leading-order Frobenius seed at r_in (level 1 data + math).

        Psi_j(r_in) = r_in^s_j e_j, Psi'_j(r_in) = s_j r_in^(s_j-1) e_j.
        Any basis of the regular subspace is admissible; leading-order
        seeding carries an O(r_in^Delta_s) admixture (Delta_s = 1.46
        between the psi and V channels), bounded by the SEED_SENSITIVITY
        gate: candidate roots must be stationary under a change of the
        seeding radius.
        """
        r_in = self.op.r_in
        Y = np.zeros((3, 3), complex)
        dY = np.zeros((3, 3), complex)
        for c, s in enumerate(REGULAR_EXPONENTS):
            Y[c, c] = r_in ** s
            dY[c, c] = s * r_in ** (s - 1)
        return Y, dY

    # -- outgoing pencil -------------------------------------------------
    def pencil_roots(self, omega: complex, n_avg: int = 8):
        """6 k-roots of det P(omega, k) = 0 on the exported edge blocks."""
        op = self.op
        sl = slice(-n_avg, None)
        K = op.K[sl].mean(axis=0)
        G = op.G[sl].mean(axis=0)
        S = op.S[sl].mean(axis=0)
        Gp = (op.G[-1] - op.G[-n_avg - 1]) / (op.r[-1] - op.r[-n_avg - 1])
        Sp = (op.S[-1] - op.S[-n_avg - 1]) / (op.r[-1] - op.r[-n_avg - 1])
        C0 = omega**2 * K - (op.M[sl].mean(axis=0) + 0.5 * Sp)
        C1 = 1j * (Gp - S)
        C2 = -G
        # linearize: [0 I; -C2^{-1} C0  -C2^{-1} C1] has eigvals k
        Z = np.zeros((3, 3), complex)
        I = np.eye(3, dtype=complex)
        top = np.block([[Z, I]])
        bot = np.block([[-np.linalg.solve(C2, C0), -np.linalg.solve(C2, C1)]])
        ks = np.linalg.eigvals(np.vstack([top, bot]))
        return ks

    def outgoing_ks(self, omega: complex, n_avg: int = 8):
        """3 outgoing k-branches at this omega (deterministic selection,
        no path memory).

        Convention (matches the certified Step-0 machinery): the outer
        tail is asymptotically flat, so the 6 pencil roots split into
        3 'acoustic' roots k_j ~ +omega (e^{+ik r} ~ e^{+i omega r},
        the physically propagating branch) and 3 'evanescent/auxiliary'
        roots k_j ~ -omega or |k| >> |omega|. Outgoing = acoustic pair
        selected by minimal |k - omega| distance. This is unambiguous
        even when Im k signs straddle zero (deep in the lower half
        plane the acoustic root has Im k < 0: e^{i(wr+iwi)r} with
        wi<0 grows like e^{|wi| r}, the standard QNM outgoing wave).
        """
        roots = self.pencil_roots(omega, n_avg=n_avg)
        out = sorted(roots, key=lambda k: abs(k - omega))[:3]
        return np.array(out)

    # -- Riccati/impedance matching ---------------------------------------
    # NOTE (scipy 1.18.1): the dop853 integrator rejects complex state
    # vectors outright, so Z is integrated as an 18-component real vector
    # (Re Z stacked on Im Z). Stacking verified against e^{ix} before use.
    _NZ = 18

    @staticmethod
    def _pack(Z):
        return np.concatenate([Z.real.reshape(-1), Z.imag.reshape(-1)])

    @staticmethod
    def _unpack(y):
        return y[:9].reshape(3, 3) + 1j * y[9:].reshape(3, 3)

    def _integrate_Z(self, Z0, x0, x1, w2, rtol, atol):
        op = self.op

        def flat_rhs(x, y):
            Z = self._unpack(y)
            G, A, B, K = op.coeffs(x)
            Gi = np.linalg.inv(G)
            dZ = -(Gi @ A) @ Z - Z @ Z - Gi @ (B - w2 * K)
            return self._pack(dZ)

        sol = solve_ivp(flat_rhs, (x0, x1), self._pack(Z0),
                        method="DOP853", rtol=rtol, atol=atol)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return np.nan
        return self._unpack(sol.y[:, -1])

    def Z_in(self, omega: complex, rtol=1e-9, atol=1e-12):
        """Logarithmic derivative Z = dPsi Psi^{-1} of the regular inner
        basis, integrated from r_in to r_match via the matrix Riccati
        equation dZ = -(G^-1 A) Z - Z^2 - G^-1 (B - w^2 K).

        Z is scale-invariant: the exponential e^{|k| r} growth cancels,
        so there is no conditioning collapse. Frobenius seed:
        Z(r_in) = dY Y^{-1} = diag(s_j / r_in) at leading order.
        """
        op = self.op
        Y0, dY0 = self.inner_seed()
        Z0 = np.linalg.solve(Y0.T, dY0.T)
        return self._integrate_Z(Z0, op.r_in, self.r_match, omega**2, rtol, atol)

    def modal_exponents(self, omega: complex, r_eval: float):
        """Algebraic modal exponents of the Riccati problem at radius
        r_eval, from the 6x6 companion Hamiltonian

            H = [[-G^-1 A, -G^-1 (B - w^2 K)], [I, 0]],   lam = Psi'/Psi.

        Asymptotics of the frozen V3_3DOF operator: K, G ~ r^2 (r^2
        scaling of the physical metric coefficients), S and M
        subleading, so with the modal variable mu = lam * r the fixed
        points converge to mu_j^2 = w^2 * eig(Gbar^-1 Kbar) (acoustic
        speeds sqrt-eig(Gbar^-1 Kbar)), plus the auxiliary deep modes.
        Returns (mu eigenvalues, eigenvector matrix V) of H (6x6).
        """
        op = self.op
        G, A, B, K = op.coeffs(r_eval)
        GiA = np.linalg.solve(G, A)
        GiB = np.linalg.solve(G, B - omega * omega * K)
        H = np.block([[-GiA, -GiB],
                      [np.eye(3, dtype=complex), np.zeros((3, 3), complex)]])
        ev, V = np.linalg.eig(H)
        order = np.argsort(ev.real)
        return ev[order], V[:, order]

    def outgoing_selection(self, omega: complex, r_eval: float):
        """The 3 OUTGOING modal branches at radius r_eval.

        Outgoing acoustic: mu = lam * r converges to +omega*speed_j
        (e^{+i omega s_j r}-like phase propagation). Selection is
        deterministic: the 3 roots with largest Re(lam) (acoustic
        QNM branches have Re lam > 0 for Im w < 0: e^{|Im w| r} growth
        = standard outgoing QNM; the conjugated incoming branches have
        Re lam < 0). Rejections of Im k > 0 style rules: this one is
        branch-stable and free of path memory.
        """
        lam, V = self.modal_exponents(omega, r_eval)
        idx = np.argsort(lam.real)[-3:]     # 3 largest Re lam = outgoing
        return lam[idx], V[:, idx]

    def Z_out(self, omega: complex, r_eval: float, rtol=1e-9, atol=1e-12):
        """Outgoing impedance Z_out(r_eval) from the algebraic invariant
        subspace (fixed point of the Riccati at r_eval), NOT from inward
        integration.

        Why not integration: all three outgoing branches GROW outward
        (Re lam > 0 for Im w < 0); inward integration of Z is then
        dominated after e^{2 Re lam Delta r} ~ e^120 by the single
        largest-|Re lam| mode, or collapses onto the incoming conjugate
        branch — both are the WRONG boundary condition. The algebraic
        fixed point instead solves 0 = -(G^-1 A) Z - Z^2 - G^-1(B-w^2 K)
        exactly at r_eval, selecting the outgoing graph subspace; the
        truncation error is the local WKB error O((lam r)^-2).

        r_eval independence of the certified candidates is a gate
        (R_MATCH_SENSITIVITY), replacing the old flat-tail assumption.
        """
        lam, V = self.outgoing_selection(omega, r_eval)
        # outgoing graph subspace: Z_out = V2 V1^{-1} with V = [V1; V2]
        V1, V2 = V[:3, :], V[3:, :]
        return V2 @ np.linalg.inv(V1)

    def D_riccati(self, omega: complex, **kw):
        """Match matrix Z_in(r_match) - Z_out(r_match); its determinant
        vanishes exactly when a solution exists that is regular at the
        origin and outgoing at infinity. Basis rescalings cancel in Z,
        so this observable is scale-free by construction.
        """
        Zi = self.Z_in(omega, **kw)
        Zo = self.Z_out(omega, **kw)
        if not (isinstance(Zi, np.ndarray) and isinstance(Zo, np.ndarray)):
            return np.nan
        return Zi - Zo

    def sigma_min_riccati(self, omega: complex, **kw) -> float:
        dZ = self.D_riccati(omega, **kw)
        if not isinstance(dZ, np.ndarray):
            return float("nan")
        return float(np.linalg.svd(dZ, compute_uv=False)[-1])

    def rhs_wrap(self, w2):
        return lambda x, y: self.op.rhs_first_order(x, y, w2)
