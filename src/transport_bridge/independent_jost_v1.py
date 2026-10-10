"""B6 step 0 — INDEPENDENT analytic-resonance machinery (Bridge implementation).

Design (independent of any Closure solver code, which has NOT been read):
  * Square-well known-answer controls on the half-line, Neumann center,
    outgoing infinity:  u'' + (omega^2 - V(r)) u = 0,
    V = V0 for r < a, 0 beyond.
  * Shooting from r = 0 outward with u(0) = 1, u'(0) = 0; the exact exterior
    outgoing condition is imposed as a Robin condition AT r = a:
        F(w) = u'(a) - i w u(a) = 0.
    No truncated-far-field, no contour deformation for these controls.
  * BLIND search: dense |F| grid scan (real/imaginary plane), cluster the
    minima, polish each with mpmath Newton on F.  The frozen reference
    poles from Closure are NEVER used to seed or steer the search; they are
    read only AFTER the candidate list is frozen to disk (enforced by the
    runner ordering + artifact hash).

Tier labeling (three-tier agreement):
  * Level 1 (shared, copied): potential parameters and the frozen reference
    files (hash-bound via data/frozen_specs/PROVENANCE_V1.json).
  * Level 2 (independent): everything numerical in this module.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

BRIDGE = Path(__file__).resolve().parents[2]
SPEC_DIR = BRIDGE / "data" / "frozen_specs"
ART_DIR = BRIDGE / "artifacts"


# ---------------------------------------------------------------- potential
def square_well(v0: float, a: float, n: int = 4001):
    """Grid + piecewise-constant potential for one channel."""
    r = np.linspace(0.0, a * 2.0, n)
    v = np.where(r < a, v0, 0.0)
    return r, v


# ----------------------------------------------------------------- Numerov
def shoot_f(w: complex, v0: float, a: float, n: int = 4001) -> complex:
    """F(w) = u'(a) - i w u(a) with u(0)=1, u'(0)=0, Numerov outward.
    Fast (few ms) but only ~1e-4 accurate near poles: used for the BLIND
    SCAN only, never for final polish."""
    r, v = square_well(v0, a, n)
    h = r[1] - r[0]
    q = v - w * w  # u'' = q u
    ia = int(np.argmin(np.abs(r - a)))
    # force evaluation exactly at the step
    if abs(r[ia] - a) > 1e-12 and ia > 0 and (r[ia] > a):
        ia -= 1
    u = np.zeros(ia + 1, dtype=complex)
    u[0], u[1] = 1.0, 1.0 + h * h * q[1] / 2.0  # Taylor start (u'(0)=0)
    c = h * h / 12.0
    for i in range(1, ia):
        u[i + 1] = (2.0 * u[i] * (1.0 + 5.0 * c * q[i])
                    - u[i - 1] * (1.0 - c * q[i - 1])) / (1.0 - c * q[i + 1])
    up = (3.0 * u[ia] - 4.0 * u[ia - 1] + u[ia - 2]) / (2.0 * h)  # 2nd order BD
    return up - 1j * w * u[ia]


def shoot_f_ivp(w: complex, v0: float, a: float,
                rtol: float = 1e-13, atol: float = 1e-14) -> complex:
    """F(w) via adaptive DOP853 integration on [0, a] (same Robin condition).
    Accurate near poles (|F| scales ~linearly with rtol); used for polish."""
    from scipy.integrate import solve_ivp

    def rhs(r, y):
        return (y[1], (v0 - w * w) * y[0])

    sol = solve_ivp(rhs, (0.0, a), (1.0 + 0j, 0.0 + 0j), method="DOP853",
                    rtol=rtol, atol=atol)
    ua, upa = sol.y[0, -1], sol.y[1, -1]
    return upa - 1j * w * ua


def blind_scan(v0: float, a: float, wr_max: float = 32.0,
               wi_min: float = -7.0, n_wr: int = 160, n_wi: int = 140,
               top: int = 24):
    """Grid scan for |F| minima; returns clustered seed guesses."""
    wrs = np.linspace(0.25, wr_max, n_wr)
    wis = np.linspace(wi_min, -0.05, n_wi)
    grid = np.empty((n_wi, n_wr))
    for j, wi in enumerate(wis):
        for i, wr in enumerate(wrs):
            grid[j, i] = abs(shoot_f(complex(wr, wi), v0, a))
    idx = np.dstack(np.unravel_index(np.argsort(grid, axis=None)[:top * 3],
                                     grid.shape))[0]
    seeds, taken = [], []
    for j, i in idx:
        w = complex(wrs[i], wis[j])
        if any(abs(w - t) < 0.6 for t in taken):
            continue
        seeds.append(w)
        taken.append(w)
        if len(seeds) >= top:
            break
    return seeds


def polish(w0: complex, v0: float, a: float, iters: int = 30,
           tol: float = 1e-9) -> tuple[complex, float]:
    """Newton on F(w) = shoot_f_ivp (DOP853, rtol=1e-13) from a scan seed.
    Converges to |F(w)| ~ 1e-11..1e-12, i.e. ~8-10 digits in w."""
    w = complex(w0)
    for _ in range(iters):
        fw = shoot_f_ivp(w, v0, a)
        if abs(fw) < tol * 1e-3:
            break
        dw = 1e-7 * max(1.0, abs(w))
        df = (shoot_f_ivp(w + dw, v0, a) - fw) / dw
        if df == 0:
            break
        w = w - fw / df
    res = abs(shoot_f_ivp(w, v0, a))
    return w, res


# ------------------------------------------------------------------ gates
def spec_provenance() -> dict:
    prov = json.loads((SPEC_DIR / "PROVENANCE_V1.json").read_text())
    for name, meta in prov["specs"].items():
        h = hashlib.sha256((SPEC_DIR / name).read_bytes()).hexdigest()
        assert h == meta["sha256"], f"spec hash drift: {name}"
    return prov


def search_channel(v0: float, a: float, wi_min: float = -7.0):
    seeds = blind_scan(v0, a, wi_min=wi_min)
    polished = []
    for s in seeds:
        try:
            w, res = polish(s, v0, a)
        except Exception:
            continue
        if res < 1e-8 and w.real > 0 and -8.0 < w.imag < 0:
            if not any(abs(w - p[0]) < 1e-6 for p in polished):
                polished.append((w, res))
    polished.sort(key=lambda t: t[0].real)
    return polished
