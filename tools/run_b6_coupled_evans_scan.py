#!/usr/bin/env python3
"""B6 Step 1, package 2: COUPLED 3-channel QR-Evans-Wronskian pole hunt.

Scans D(w) = sigma_min([Q_int, Q_ext_incoming]) of the FULL 3-channel
frozen operator (StableEvansV3) over preregistered tiles in the upper
half w-plane (scan convention: poles reported as conjugates, exact
real-coefficient symmetry, gated).

Validation anchor (preregistered BEFORE the scan): the certified ch2
pilot pole w_p = 0.22642588643752676 - 0.08455094304111954j (QNM sign)
must appear as a zero of the full D in its conjugated scan form
w = 0.22642588643752676 + 0.08455094304111954j, because ch2 is (G,K)-
decoupled from (ch0,ch1).  If the scan + refinement reproduces it, the
coupled machinery is validated against the 1x1 pilot; additional refined
zeros are NEW coupled candidates (ch0-ch1 block couples via B01).
  SUPERSEDED 12.10.26 (SSZ_B6_PILOT_SUPERSESSION_V1): the pilot "pole"
  was disproved by the mpmath control solver (no zero, no delta scaling,
  no dip). It is NO LONGER a valid validation anchor or reference; the
  coupled-scan reproduction of it is a machinery cross-check only, not
  evidence of a physical resonance.

Acceptance gates per candidate (mirroring the pilot's certification):
  RETIRED (SSZ_B6_GATE_RECALIBRATION_V1): the old absolute < 1e-12 floor
  and relative negative control passed pure noise vacuously -- all five
  2026-10-12 GATED_POLE_CANDIDATEs failed mpmath confirmation. Gating is
  now relative (measured local floor, dip >= 2 OD, visible |E|~|delta|
  scaling) + mandatory mpmath control-solver confirmation; see the
  recalibrated gate() below.
  * r_match independence (25/30/35): |D| stays at floor level,
  * r_seed independence (50/55/59.5): same,
  * conjugation symmetry: |D(w*) - D(w)| small relative to floor,
  * negative control: |D| at w + 0.01 (fixed offset) well above floor.
Classification: GATED_POLE_CANDIDATE / REJECTED_<reason>.

Strict-JSON artifact; preregistration written to the artifact BEFORE
results are appended (two-phase file write, phase marker kept).
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

from transport_bridge.independent_qnm_3dof_v1 import FrozenOperator  # noqa: E402
from transport_bridge.stable_evans_v3 import StableEvansV3  # noqa: E402

NPZ = BRIDGE / "data" / "frozen_specs" / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
PILOT_W = complex(0.22642588643752676, 0.08455094304111954)  # scan (Im>0) form

# ---- preregistered tiles (Re, Im bounds + step) -----------------------------
TILES = {
    "T0_PILOT_ANCHOR": {"re": (0.15, 0.31, 0.02), "im": (0.04, 0.13, 0.015)},
    "T1_LOWRE":        {"re": (0.05, 1.20, 0.08), "im": (0.03, 0.60, 0.06)},
    "T2_HIGHRE":       {"re": (1.20, 3.00, 0.15), "im": (0.05, 0.80, 0.10)},
}

FLOOR_TARGET = 1e-12     # RETIRED (SSZ_B6_GATE_RECALIBRATION_V1): kept only as a
                         # coarse pre-filter; certification now uses the relative
                         # gates G1-G7 below, never absolute |D| thresholds.
REFINE_MAX_D = 1e-6      # beyond this a "minimum" is not pursued
EVAN_SEGS = 40
RTOL = 1e-9


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    t0 = time.time()
    op = FrozenOperator(str(NPZ))
    ev = StableEvansV3(op, r_match=30.0, r_seed=55.0,
                       n_segments=EVAN_SEGS, rtol=RTOL)

    # ---- timing probe (sizes nothing else) ----------------------------------
    tp = time.time()
    d_probe = ev.D(PILOT_W)
    t_eval = time.time() - tp
    n_grid = sum(
        len(np.arange(*tile["re"][:2], tile["re"][2]) + 0j)
        * len(np.arange(*tile["im"][:2], tile["im"][2]))
        for tile in TILES.values()
    )
    est_h = n_grid * t_eval / 3600.0
    print(f"D(w) eval: {t_eval:.2f}s | grid points: {n_grid} "
          f"| est scan {est_h:.2f}h", flush=True)
    print(f"D(pilot_anchor) = {d_probe:.3e}  (validation anchor, "
          f"must be far below generic O(0.1-1) level)", flush=True)

    prereg = {
        "artifact": "SSZ_B6_STEP1_COUPLED_EVANS_SCAN_V1",
        "preregistered_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "npz_sha256": sha256(NPZ),
        "method": "full 3-channel StableEvansV3 D(w)=sigma_min([Q_int,Q_ext])",
        "tiles_preregistered": TILES,
        "validation_anchor": {"w_scan_form": [PILOT_W.real, PILOT_W.imag],
                              "origin": "SSZ_B6_STEP1_CH2_PILOT_POLE_V1"},
        "acceptance": {
            "refined_floor_target": FLOOR_TARGET,
            "gates": ["r_match 25/30/35", "r_seed 50/55/59.5",
                      "conjugation symmetry", "negative control +0.01 offset"],
        },
        "d_eval_seconds": t_eval,
        "d_at_anchor": d_probe,
    }
    out_path = BRIDGE / "artifacts" / "SSZ_B6_STEP1_COUPLED_EVANS_SCAN_V1.json"
    out_path.write_text(json.dumps({"preregistration": prereg,
                                    "phase": "SCAN_RUNNING"},
                                   indent=1) + "\n")

    # ---- scan ---------------------------------------------------------------
    rows = []
    minima = []  # (d, w, tile); sort/min always with key=lambda x: x[0]
    n_done = 0
    for tname, tile in TILES.items():
        re_a, re_b, re_s = tile["re"]
        im_a, im_b, im_s = tile["im"]
        for re_w in np.arange(re_a, re_b + 0.5 * re_s, re_s):
            for im_w in np.arange(im_a, im_b + 0.5 * im_s, im_s):
                w = complex(float(re_w), float(im_w))
                try:
                    d = ev.D(w)
                except RuntimeError as exc:
                    rows.append({"re": w.real, "im": w.imag, "D": None,
                                 "tile": tname, "error": str(exc)[:120]})
                    continue
                rows.append({"re": w.real, "im": w.imag, "D": d,
                             "tile": tname})
                minima.append((d, w, tname))
                n_done += 1
                if n_done % 25 == 0:
                    best = min(minima, key=lambda x: x[0])[0]
                    print(f"  [{n_done}/{n_grid}] last tile={tname} "
                          f"D={d:.3e} best_so_far={best:.3e} "
                          f"({time.time()-t0:.0f}s)", flush=True)

    grid_ok = [r for r in rows if r["D"] is not None]
    finite = [(r["D"], complex(r["re"], r["im"]), r["tile"])
              for r in grid_ok]
    finite.sort(key=lambda x: x[0])

    # ---- Lino's rule (12.10.26): persist the FULL grid BEFORE refinement ----
    # A crash in the optimizer must not cost the 294 expensive D(w) evals.
    grid_npz = BRIDGE / "artifacts" / "SSZ_B6_STEP1_COUPLED_EVANS_GRID_V1.npz"
    np.savez_compressed(
        grid_npz,
        re=np.array([r["re"] for r in grid_ok], dtype=np.float64),
        im=np.array([r["im"] for r in grid_ok], dtype=np.float64),
        D=np.array([r["D"] for r in grid_ok], dtype=np.float64),
        tile=np.array([r["tile"] for r in grid_ok]),
        grid_preregistered_count=n_grid,
        n_errors=len(rows) - len(grid_ok),
        d_eval_seconds=t_eval,
        d_at_anchor=d_probe,
        npz_sha256=sha256(NPZ),
    )
    print(f"  grid persisted: {grid_npz.name} "
          f"({len(grid_ok)} pts, sha256={sha256(grid_npz)[:16]}...)",
          flush=True)

    # ---- refinement: top-k distinct minima ----------------------------------
    def logD(x):
        w = complex(float(x[0]), float(x[1]))  # NM passes a 2-vector
        d = ev.D(w)
        return float(np.log(max(d, 1e-300)))

    candidates = []
    for d0, w0, tile in finite[:12]:
        if any(abs(w0 - c["w_complex"]) < 0.03 for c in candidates):
            continue
        res = minimize(logD, np.array([w0.real, w0.imag]),
                       method="Nelder-Mead",
                       options={"xatol": 1e-10, "fatol": 1e-12,
                                "maxiter": 300})
        w_ref = complex(*res.x)
        d_ref = float(np.exp(min(res.fun, 0.0)))
        candidates.append({"w_complex": w_ref, "D_refined": d_ref,
                           "D_grid": d0, "tile": tile,
                           "nm_success": bool(res.success)})
        print(f"  refined: w={w_ref:.10f} D={d_ref:.3e} "
              f"(grid {d0:.3e})", flush=True)
        if len(candidates) >= 5 or d_ref > REFINE_MAX_D:
            break

    # ---- gating per candidate ----------------------------------------------
    # Recalibrated 12.10.26 per SSZ_B6_GATE_RECALIBRATION_V1: the old absolute
    # FLOOR_TARGET + relative negative control passed pure noise vacuously
    # (all five 2026-10-12 candidates failed mpmath confirmation). Now:
    # relative dip vs the MEASURED floor + |E|~|delta| scaling are mandatory.
    def gate(w: complex) -> dict:
        g = {}
        g["d_ref"] = ev.D(w)
        for rm in (25.0, 35.0):
            alt = StableEvansV3(op, r_match=rm, r_seed=55.0,
                                n_segments=EVAN_SEGS, rtol=RTOL)
            g[f"d_rmatch_{rm:.0f}"] = alt.D(w)
        for rs in (50.0, 59.5):
            alt = StableEvansV3(op, r_match=30.0, r_seed=rs,
                                n_segments=EVAN_SEGS, rtol=RTOL)
            g[f"d_rseed_{rs:.1f}"] = alt.D(w)
        g["d_conj"] = ev.D(w.conjugate())
        g["conj_gap"] = abs(g["d_conj"] - g["d_ref"])
        g["d_negctrl_offset0.01"] = ev.D(w + 0.01)
        vals = [g[k] for k in g if k.startswith(("d_rmatch", "d_rseed"))]
        # --- recalibrated gates (SSZ_B6_GATE_RECALIBRATION_V1) ---
        # measured local floor: rosette |D| at 4 offsets around w
        _ros = [abs(ev.D(w + complex(dx, dy)))
                for dx, dy in ((1e-3, 0.0), (-1e-3, 0.0),
                               (0.0, 1e-3), (0.0, -1e-3))]
        _floor = sorted(_ros)[1]          # second-smallest: robust median-ish
        _dip = _floor / max(abs(g["d_ref"]), 1e-300)
        _scal = [abs(ev.D(w + complex(10.0 ** e, 0.0)))
                 for e in (-3.0, -4.0, -5.0)]
        _ok_lin = all(s2 > 50.0 * s1 or s1 == 0.0
                      for s1, s2 in zip(_scal, _scal[1:]))
        g["measured_floor"] = _floor
        g["dip_orders"] = math.log10(_dip) if _dip > 0 else float("inf")
        g["scaling_monotone_50x_per_decade"] = bool(_ok_lin)
        # CLASSIFICATION: 'GATED' now means dip >= 2 OD AND scaling visible.
        # Final confirmation ALWAYS requires the mpmath control solver
        # (tools/b6_mpmath_evans_v2.py) rosette + dps ladder per G5.
        g["gate_pass"] = bool(
            g["dip_orders"] >= 2.0
            and g["scaling_monotone_50x_per_decade"]
            and max(vals) < 1e-9
            and g["conj_gap"] < 0.5 * g["d_ref"] + 1e-15
            and g["d_negctrl_offset0.01"] > 100 * _floor)
        g["gate_version"] = "RECALIBRATED_V1"
        return g

    final = []
    for c in candidates:
        if c["D_refined"] >= REFINE_MAX_D:
            c["classification"] = "REJECTED_NOT_A_ZERO"
            final.append(c)
            continue
        g = gate(c["w_complex"])
        c["gates"] = g
        c["classification"] = ("GATED_POLE_CANDIDATE" if g["gate_pass"]
                               else "REJECTED_GATE_FAILURE")
        # report in QNM sign (Im<0)
        c["pole_qnm_sign"] = {"re": c["w_complex"].real,
                              "im": -c["w_complex"].imag}
        final.append(c)
        print(f"  gated: w={c['w_complex']:.10f} -> {c['classification']}",
              flush=True)

    # anchor check: was the pilot pole among the refined candidates?
    anchor_hit = min((abs(c["w_complex"] - PILOT_W) for c in candidates),
                     default=None)

    artifact = {
        "preregistration": prereg,
        "phase": "COMPLETE",
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "wall_seconds": round(time.time() - t0, 1),
        "n_grid_points": n_grid,
        "n_evaluated": len(grid_ok),
        "n_errors": len(rows) - len(grid_ok),
        "validation_anchor_result": {
            "d_at_anchor_before_scan": d_probe,
            "nearest_refined_candidate_distance": anchor_hit,
            "anchor_reproduced": bool(anchor_hit is not None
                                      and anchor_hit < 1e-6),
        },
        "grid_top20": [{k: r[k] for k in ("re", "im", "D", "tile")}
                       for r in grid_ok[:20]],
        "refined_candidates": [
            {k: (v if k != "w_complex" else [v.real, v.imag])
             for k, v in c.items()} for c in final],
        "provenance": {
            "operator": NPZ.name,
            "npz_sha256": prereg["npz_sha256"],
            "evans_segments": EVAN_SEGS,
            "rtol": RTOL,
            "supersedes": "SSZ_B6_STEP1_BLIND_SCAN_V1 (wrong BVP, Riccati)",
        },
    }
    out_path.write_text(json.dumps(artifact, indent=1, allow_nan=False) + "\n")
    print("ARTIFACT:", out_path)
    print("anchor_reproduced:",
          artifact["validation_anchor_result"]["anchor_reproduced"])
    for c in final:
        print(f"  {c['classification']}: "
              f"pole_qnm={c['pole_qnm_sign']['re']:.10f}"
              f"{c['pole_qnm_sign']['im']:+.10f}j  D={c['D_refined']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
