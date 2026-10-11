#!/usr/bin/env python3
"""B6 Step 1: blind sigma_min scan over the preregistered tiles D1/D2/D3.

Replicates the preregistered no-resonance scan of the Closure with the
INDEPENDENT Bridge machinery (Riccati/impedance matching, DOP853,
deterministic outgoing selection). Blind: no Closure verdicts, poles,
or candidate lists are read anywhere in this pipeline.

Tiles come from data/frozen_specs/SSZ_DOMAIN_EXPANSION_V1_PREREGISTRATION.json.
Scan resolution: 61 x 45 per tile (sigma_min is smooth; clusters matter,
not single points). rtol=1e-7 (sigma_min stable to 3 significant digits,
verified); polish at rtol=1e-11 happens later, only on candidates.
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

NWR, NWI = 61, 45
RTOL_SCAN = 1e-7
N_WORKERS = 14

# (re0, re1, im0, im1, name)
TILE_BOUNDS = {
    "D1_HIGH_RE": (1.5, 3.0, -0.8, -0.05),
    "D2_DEEP_IM": (0.05, 1.5, -2.0, -0.8),
    "D3_HIGH_RE_DEEP_IM": (1.5, 3.0, -2.0, -0.8),
}


def scan_row(args):
    import multiprocessing as _mp
    tile, j, re_grid, wi = args
    from transport_bridge.independent_qnm_3dof_v1 import (
        FrozenOperator, IndependentQNM3DOF)

    op = FrozenOperator(str(BRIDGE / "data/frozen_specs/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"))
    slv = IndependentQNM3DOF(op)
    row = np.empty(len(re_grid))
    for i, wr in enumerate(re_grid):
        try:
            row[i] = slv.sigma_min_riccati(complex(wr, wi), rtol=RTOL_SCAN, atol=1e-10)
        except Exception:
            row[i] = np.nan
    return tile, j, row


def main():
    t_start = time.time()
    re_grids = {t: np.linspace(b[0], b[1], NWR) for t, b in TILE_BOUNDS.items()}
    im_grids = {t: np.linspace(b[2], b[3], NWI) for t, b in TILE_BOUNDS.items()}

    jobs = []
    for tile, im_g in im_grids.items():
        for j, wi in enumerate(im_g):
            jobs.append((tile, j, re_grids[tile], float(wi)))
    print(f"total rows: {len(jobs)} over {len(TILE_BOUNDS)} tiles "
          f"({NWR}x{NWI} each); {N_WORKERS} workers", flush=True)

    grids = {t: np.full((NWI, NWR), np.nan) for t in TILE_BOUNDS}
    done = 0
    ctx = __import__("multiprocessing").get_context("spawn")
    with ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=ctx) as ex:
        for tile, j, row in ex.map(scan_row, jobs):
            grids[tile][j] = row
            done += 1
            if done % 100 == 0:
                el = time.time() - t_start
                eta = el / done * (len(jobs) - done)
                print(f"  {done}/{len(jobs)} rows  [{el:.0f}s elapsed, "
                      f"{eta:.0f}s ETA]", flush=True)

    out = {
        "artifact": "SSZ_B6_STEP1_BLIND_SCAN_V1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "frozen_spec": "data/frozen_specs/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz",
        "npz_sha256": "ad2bee7d3994097479933c1eea787f3aa381648d098ffd36aab32a4c436db2ef",
        "method": "independent Riccati/impedance matching (DOP853, rtol=1e-7 scan),"
                  " deterministic outgoing selection by |k - omega|",
        "preregistration": "data/frozen_specs/SSZ_DOMAIN_EXPANSION_V1_PREREGISTRATION.json",
        "grid": {"nwr": NWR, "nwi": NWI, "rtol": RTOL_SCAN,
                 "atol": 1e-10, "workers": N_WORKERS},
        "tiles": {},
    }
    art_dir = BRIDGE / "artifacts"
    art_dir.mkdir(exist_ok=True)
    summary = {}
    for tile, g in grids.items():
        b = TILE_BOUNDS[tile]
        valid = g[np.isfinite(g)]
        jmin, imin = np.unravel_index(np.nanargmin(g), g.shape)
        out["tiles"][tile] = {
            "bounds": {"Re": [b[0], b[1]], "Im": [b[2], b[3]]},
            "grid_sigma_min": g.tolist(),
            "sigma_min_global": float(valid.min()),
            "sigma_min_median": float(np.median(valid)),
            "argmin_omega": [float(re_grids[tile][imin]), float(im_grids[tile][jmin])],
            "n_nan": int(g.size - valid.size),
        }
        summary[tile] = (float(valid.min()), float(np.median(valid)),
                         [float(re_grids[tile][imin]), float(im_grids[tile][jmin])])
    path = art_dir / "SSZ_B6_STEP1_BLIND_SCAN_V1.json"
    path.write_text(json.dumps(out, indent=1))
    print("\n=== per-tile summary ===")
    for t, (mn, md, am) in summary.items():
        print(f"{t}: sigma_min={mn:.4e}  median={md:.4e}  argmin={am}  "
              f"ratio min/median={mn/md:.3f}")
    print(f"\nwritten: {path}  [{time.time()-t_start:.0f}s total]")


if __name__ == "__main__":
    main()
