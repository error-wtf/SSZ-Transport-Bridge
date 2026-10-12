#!/usr/bin/env python3
"""B6 Step 1, package 4: local delta-scaling + singular-value control (V1).

Reviewer request (Carmen, 12.10.26): classify the three refined Evans
candidates INDEPENDENT of |D| magnitude and of the argument principle.

For each center w0 (3 candidates + certified pilot pole as positive
control) measure

    |D(w0 + delta)|   vs   |delta|

on a log-spaced delta ladder along two independent directions, and fit
the local power law  |D| ~ C |delta|^m:

    m ~ 1  -> simple zero (or m ~ 2: double zero),
    m ~ 0 with constant floor -> numerical plateau (no resolution).

Additionally, at w0 and at one mid-ladder delta, repeat D with tightened
tolerance (rtol 1e-10 / 1e-11, same segment structure) and record the
FULL singular-value spectrum of the 6x6 matching matrix.  If sigma_min
changes systematically under rtol, the floor is discretization-driven.

Scan convention: candidates are Im w > 0 (scan form); physics unchanged
(D is evaluated in the same convention as the scan).

Strict-JSON artifact; every D value recorded with its parameters.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

from transport_bridge.independent_qnm_3dof_v1 import FrozenOperator  # noqa: E402
from transport_bridge.stable_evans_v3 import StableEvansV3  # noqa: E402

NPZ = BRIDGE / "data" / "frozen_specs" / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
ARTIFACT = BRIDGE / "artifacts" / "SSZ_B6_STEP1_DELTA_SCALING_V1.json"

CANDIDATES = {
    "cand1": complex(0.051875, 0.0285),
    "cand2": complex(0.0499768066, 0.2101666260),
    "cand3": complex(0.0521972656, 0.1479492188),
    "pilot_pos_control": complex(0.22642588643752676, 0.08455094304111954),
}

# delta ladder: decades from 1e-2 down to 3e-8 (10 points), two directions
DELTAS = np.logspace(-2, np.log10(3e-8), 10)
DIRECTIONS = [1.0 + 0.0j, 0.0 + 1.0j]   # radial + vertical, |dir| = 1
RTOL_REF = 1e-9
RTOL_FINE = (1e-10, 1e-11)
NSEG = 40


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fit_slope(logx: np.ndarray, logy: np.ndarray):
    """Least-squares slope on the FINITE, positive subset; returns
    (slope, intercept, n_used)."""
    mask = np.isfinite(logx) & np.isfinite(logy) & (logy > -np.inf)
    if mask.sum() < 3:
        return None, None, int(mask.sum())
    x, y = logx[mask], logy[mask]
    A = np.vstack([x, np.ones_like(x)]).T
    (m, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(m), float(b), int(mask.sum())


def main() -> int:
    t0 = time.time()
    op = FrozenOperator(str(NPZ))
    ev = StableEvansV3(op, r_match=30.0, r_seed=55.0,
                       n_segments=NSEG, rtol=RTOL_REF)

    results = {}
    for name, w0 in CANDIDATES.items():
        t_c = time.time()
        entry = {"w0": [w0.real, w0.imag]}

        # ---- delta ladder, two directions -----------------------------------
        ladders = {}
        for k, direction in enumerate(DIRECTIONS):
            pts = []
            for delta in DELTAS:
                w = w0 + direction * delta
                t_e = time.time()
                d, det = ev.D(w, details=True)
                pts.append({"delta": float(delta),
                            "D": float(d),
                            "svals": det["svals"],
                            "eval_s": round(time.time() - t_e, 2)})
            ladders[f"dir{k}"] = {
                "direction_re_im": [direction.real, direction.imag],
                "points": pts,
                "slope_loglog": None,   # filled below
            }
        # fits: use only points above a plausible noise knee is CHEATING;
        # report the full-ladder fit AND the top-half-delta fit separately
        for k in ladders:
            pts = ladders[k]["points"]
            logd = np.log10([p["delta"] for p in pts])
            logD = np.log10([p["D"] for p in pts])
            m_full, b_full, n_full = fit_slope(logd, logD)
            m_top, b_top, n_top = fit_slope(logd[:5], logD[:5])   # larger deltas
            ladders[k]["slope_loglog"] = {
                "full": {"m": m_full, "b": b_full, "n": n_full},
                "large_delta_half": {"m": m_top, "b": b_top, "n": n_top},
            }
        entry["ladders"] = ladders

        # ---- floor reference: D at the center itself -------------------------
        d_center, det_center = ev.D(w0, details=True)
        entry["D_center"] = {"D": float(d_center), "svals": det_center["svals"]}

        # ---- tolerance sweep at center + one mid delta -----------------------
        tol_sweep = {}
        mid_delta = float(DELTAS[4])
        probes = [("center", w0),
                  ("mid_delta", w0 + DIRECTIONS[0] * mid_delta)]
        for label, w_probe in probes:
            for rtol in RTOL_FINE:
                ev_f = StableEvansV3(op, r_match=30.0, r_seed=55.0,
                                     n_segments=NSEG, rtol=rtol)
                t_e = time.time()
                d, det = ev_f.D(w_probe, details=True)
                tol_sweep[f"{label}_rtol{rtol:.0e}"] = {
                    "w": [w_probe.real, w_probe.imag],
                    "rtol": rtol,
                    "D": float(d),
                    "svals": det["svals"],
                    "eval_s": round(time.time() - t_e, 2)}
        entry["tolerance_sweep"] = tol_sweep
        entry["wall_s"] = round(time.time() - t_c, 1)
        results[name] = entry
        print(f"[{name}] done in {entry['wall_s']}s "
              f"(D_center={d_center:.3e}, "
              f"slope_full={ladders['dir0']['slope_loglog']['full']['m']})",
              flush=True)

    artifact = {
        "artifact": "SSZ_B6_STEP1_DELTA_SCALING_V1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "npz_sha256": sha256(NPZ),
        "method": ("local power-law fit |D(w0+delta)| vs |delta|, "
                   "2 directions, plus rtol sweep with full sval spectrum"),
        "parameters": {"r_match": 30.0, "r_seed": 55.0,
                       "n_segments": NSEG, "rtol_ref": RTOL_REF,
                       "deltas": [float(d) for d in DELTAS]},
        "results": results,
        "wall_total_s": round(time.time() - t0, 1),
    }
    ARTIFACT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"artifact written: {ARTIFACT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
