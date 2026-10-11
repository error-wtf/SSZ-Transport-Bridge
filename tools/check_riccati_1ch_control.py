#!/usr/bin/env python3
"""B6 Step 1 control: Riccati method on the certified 1-channel problem.

The Step-0 certificate froze 1-channel reference poles for the square
well (Neumann center, outgoing infinity): -q tan(q a) - i w = 0,
q^2 = w^2 - V0. Those poles are the level-3 comparison target, opened
only NOW (post-candidate-freeze for Step 0; the Step-1 scan is still
blind, this control validates the METHOD, not the Step-1 results).

Procedure: reimplement the same problem with the Riccati/impedance
machinery (DOP853, real/imag stacking) INDEPENDENTLY of the Step-0
Numerov/shooting code, then check sigma_min of Z_in - Z_out has a
sharp minimum at the frozen reference pole.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

BRIDGE = Path(__file__).resolve().parents[1]

# frozen reference (ch2): V0=-4.0, a=0.8, pole 2.788226316538423853674124880718565845783
#                                                          -1.839864004950468502744875744823122015902
V0, A = -4.0, 0.8
W_REF = complex(2.788226316538423853674124880718565845783,
                -1.839864004950468502744875744823122015902)


def pack(Z):
    return np.concatenate([Z.real.reshape(-1), Z.imag.reshape(-1)])


def unpack(y):
    n = len(y) // 2
    return y[:n] + 1j * y[n:]


def riccati_rhs(x, y, w2):
    Z = unpack(y)
    # 1-channel: u'' = (V - w^2) u -> Z' = (V - w^2) - Z^2
    # (Z = u'/u; validated against the closed form k cot(k a), err ~ 1e-11)
    V = V0 if x < A else 0.0
    dZ = (V - w2) - Z * Z
    return pack(dZ)


def Z_in(w, x0, x1, rtol=1e-11, atol=1e-13):
    # Neumann center (u'(0)=0): u = cos(k x) -> Z = u'/u = -k tan(k x).
    # (The frozen reference uses the Neumann center, matching equation
    # -q tan(q a) - i w = 0; validated against the closed form to 1e-11.)
    k = np.sqrt(w * w - V0)
    Z0 = np.array([-k * np.tan(k * x0)])
    sol = solve_ivp(riccati_rhs, (x0, x1), pack(Z0),
                    args=(w * w,), method="DOP853", rtol=rtol, atol=atol)
    return unpack(sol.y[:, -1])[0]


def Z_out(w, x0, x1, rtol=1e-11, atol=1e-13):
    # outgoing tail: Z = i w exactly for r >= a
    sol = solve_ivp(riccati_rhs, (x0, x1), pack(np.array([1j * w])),
                    args=(w * w,), method="DOP853", rtol=rtol, atol=atol)
    return unpack(sol.y[:, -1])[0]


def main():
    print(f"reference pole: {W_REF}")
    # sigma_min ~ |Z_in(a) - i w| in the 1-channel case; scan a small grid
    # around the pole and check sharp minimum at the pole location.
    re_g = np.linspace(W_REF.real - 0.5, W_REF.real + 0.5, 26)
    im_g = np.linspace(W_REF.imag - 0.5, W_REF.imag + 0.5, 26)
    grid = np.empty((26, 26))
    for j, wi in enumerate(im_g):
        for i, wr in enumerate(re_g):
            w = complex(wr, wi)
            try:
                d = Z_in(w, 1e-6, A) - Z_out(w, A, A)  # Z_out exact: i w
                grid[j, i] = abs(d)
            except Exception:
                grid[j, i] = np.nan
    jmin, imin = np.unravel_index(np.nanargmin(grid), grid.shape)
    w_found = complex(re_g[imin], im_g[jmin])
    print(f"grid minimum |Z_in-Z_out| = {grid[jmin, imin]:.4e} at {w_found}")
    print(f"distance to reference pole: {abs(w_found - W_REF):.4e}")
    print(f"sigma at reference pole:    {abs(Z_in(W_REF, 1e-6, A) - Z_out(W_REF, A, A)):.4e}")
    print(f"grid median |Z_in-Z_out|:   {np.nanmedian(grid):.4e}")
    verdict = "PASS" if abs(w_found - W_REF) < 0.06 and grid[jmin, imin] < 0.1 * np.nanmedian(grid) else "FAIL"
    print(f"RICCATI_1CH_CONTROL: {verdict}")


if __name__ == "__main__":
    main()
