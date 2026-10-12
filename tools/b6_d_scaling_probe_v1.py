#!/usr/bin/env python3
"""B6 Step 1, package 4: |D| scaling probe on rings (V1).

Decides zero-vs-plateau WITHOUT relying on phase (double precision
cannot carry the argument principle at |W| ~ 1e-18).  For a true zero
w0 of order m:  |D(w0 + delta)| ~ C * |delta|^m.
For a numerical plateau, |D| is flat/noisy at the floor level.

Rings around two centers:
  * candidate w = 0.051875 + 0.0285j  (scan convention)
  * certified pilot pole w = 0.22642588643752676 + 0.08455094304111954j
    (scan form, known true zero -> positive control)
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
CENTERS = {
    "candidate": complex(0.051875, 0.0285),
    "pilot_anchor": complex(0.22642588643752676, 0.08455094304111954),
}
RADII = [1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2]
DIRS = [0.0, np.pi / 2, np.pi, 3 * np.pi / 2]


def main() -> int:
    t0 = time.time()
    op = FrozenOperator(str(NPZ))
    ev = StableEvansV3(op, r_match=30.0, r_seed=55.0,
                       n_segments=40, rtol=1e-9)

    results = {}
    for label, w0 in CENTERS.items():
        rows = []
        d0 = float(ev.D(w0))
        rows.append({"delta": 0.0, "absD": d0})
        print(f"{label}: |D(w0)| = {d0:.3e}", flush=True)
        for r in RADII:
            vals = []
            for th in DIRS:
                vals.append(float(abs(ev.D(w0 + r * np.exp(1j * th)))))
            m = float(np.median(vals))
            rows.append({"delta": r, "absD_median": m,
                         "absD_min": float(min(vals)),
                         "absD_max": float(max(vals))})
            print(f"  r={r:<8.0e} |D|med={m:.3e} "
                  f"[{min(vals):.2e}, {max(vals):.2e}]", flush=True)
        # power-law fit on the growing part (skip r=0 and any floor pts)
        xs = np.log([row["delta"] for row in rows[1:]])
        ys = np.log([max(row.get("absD_median", row["absD"]), 1e-300)
                     for row in rows[1:]])
        slope, intercept = np.polyfit(xs, ys, 1)
        results[label] = {"w0": [w0.real, w0.imag], "d0": d0,
                          "rows": rows,
                          "loglog_slope": round(float(slope), 3),
                          "fit_intercept": round(float(np.exp(intercept)), 6)}
        print(f"  -> loglog slope = {slope:.2f} "
              f"(true zero: ~1 simple / ~2 double; plateau: ~0)",
              flush=True)

    verdict = {
        "pilot_slope": results["pilot_anchor"]["loglog_slope"],
        "candidate_slope": results["candidate"]["loglog_slope"],
        "note": ("slope ~m at a true order-m zero; ~0 = plateau/artifact"),
    }
    artifact = {
        "artifact": "SSZ_B6_STEP1_D_SCALING_PROBE_V1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "method": "|D| on rings r=1e-6..1e-2, 4 directions, median",
        "npz_sha256": hashlib.sha256(NPZ.read_bytes()).hexdigest(),
        "results": results,
        "verdict": verdict,
    }
    out = BRIDGE / "artifacts" / "SSZ_B6_STEP1_D_SCALING_PROBE_V1.json"
    out.write_text(json.dumps(artifact, indent=1) + "\n")
    print(f"done -> {out}  (wall {time.time() - t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
