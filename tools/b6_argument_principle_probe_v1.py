#!/usr/bin/env python3
"""B6 Step 1, package 3: Argument-principle pole counting (V1).

Counts zeros of the COMPLEX Evans-Wronskian

    W(w) = det([Q_int(w) | Q_ext_incoming(w)])      (6x6, complex)

inside a small square contour around a candidate frequency, via the
argument principle  N = (1/2pi i) oint dW/W = total unwrapped phase
windings of W along the contour.  W and D = sigma_min(...) share zeros
because the stacked matrix is singular exactly when det vanishes.

Key requirement (reviewer, 12.10.26): classify the candidate
w = 0.051875 + 0.0285j (scan convention, Im w > 0) by an
algorithm-independent test -- not by the magnitude of D alone.

Controls run in the same artifact:
  * N at the certified pilot pole (w_p = 0.22642588643752676
    + 0.08455094304111954j, scan form)  -> must be >= 1,
  * N at an offset contour (candidate + offset with no expected zero)
    -> must be 0,
  * N around the candidate for several contour half-widths h
    (robustness against a zero sitting just outside).
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
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
CANDIDATE = complex(0.051875, 0.0285)          # scan convention (Im > 0)
PILOT_W = complex(0.22642588643752676, 0.08455094304111954)  # scan form
OFFSET_CTRL = complex(0.06, 0.045)             # no expected zero nearby


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_W(ev: StableEvansV3):
    def W(omega: complex) -> complex:
        qi = ev.interior_basis(omega)
        qe = ev.exterior_incoming_basis(omega)
        return complex(np.linalg.det(np.hstack([qi, qe])))

    return W


def winding(W, center: complex, half: float, per_side: int = 8) -> dict:
    """Unwrapped-phase winding number of W on a square contour."""
    corners = [center + complex(-half, -half), center + complex(half, -half),
               center + complex(half, half), center + complex(-half, half),
               center + complex(-half, -half)]
    ts = []
    for a, b in itertools.pairwise(corners):
        s = np.linspace(0.0, 1.0, per_side, endpoint=False)
        for t in s:
            ts.append(a + (b - a) * t)
    vals = []
    for w in ts:
        vals.append(W(w))
    vals = np.asarray(vals, complex)
    if not np.all(np.isfinite(vals)):
        return {"N": None, "n_finite": int(np.isfinite(vals).sum()),
                "n_evals": len(vals)}
    dphi = np.unwrap(np.angle(vals))
    # closure term: phase from last point back to first (shortest arc).
    # Sign convention: total phase = (last - first) + angle(v0 / v_last).
    # (v0*conj(v_last) would halve positive windings and can flip the
    # sign for N near -1; verified against W(w)=w-w0 which must give +1.)
    total = float(dphi[-1] - dphi[0] + np.angle(vals[0] * np.conj(vals[-1])))
    N = total / (2.0 * np.pi)
    return {"N": round(N, 4), "n_evals": len(vals),
            "W_abs_min": float(np.abs(vals).min()),
            "W_abs_max": float(np.abs(vals).max())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--half", type=float, nargs="+",
                    default=[0.004, 0.008, 0.015],
                    help="contour half-widths around each center")
    ap.add_argument("--per-side", type=int, default=8)
    args = ap.parse_args()

    t0 = time.time()
    op = FrozenOperator(str(NPZ))
    ev = StableEvansV3(op, r_match=30.0, r_seed=55.0,
                       n_segments=40, rtol=1e-9)
    W = make_W(ev)
    t1 = time.time()
    _ = W(PILOT_W)  # warm-up + timing probe
    t_eval = time.time() - t1

    results = {"centers": {}}
    for center, label in [(CANDIDATE, "candidate"),
                          (PILOT_W, "pilot_anchor"),
                          (OFFSET_CTRL, "offset_control")]:
        entry = {"w": [center.real, center.imag], "contours": {}}
        halves = args.half if label == "candidate" else [0.008]
        for h in halves:
            r = winding(W, center, h, args.per_side)
            entry["contours"][f"h={h}"] = r
            print(f"{label:15s} h={h:<6} N={r['N']} "
                  f"(|W|min={r.get('W_abs_min', float('nan')):.3e})",
                  flush=True)
        results["centers"][label] = entry

    n_cand = [v["N"] for v in
              results["centers"]["candidate"]["contours"].values()]
    n_pilot = results["centers"]["pilot_anchor"]["contours"]["h=0.008"]["N"]
    n_ctrl = results["centers"]["offset_control"]["contours"]["h=0.008"]["N"]
    verdict = {
        "candidate_N_values": n_cand,
        "pilot_N": n_pilot,
        "control_N": n_ctrl,
        "candidate_is_zero": bool(all(n is not None and abs(n) >= 0.5
                                      for n in n_cand)),
        "pilot_reproduced": bool(n_pilot is not None and abs(n_pilot) >= 0.5),
        "control_clean": bool(n_ctrl is not None and abs(n_ctrl) < 0.5),
    }
    verdict["overall"] = ("POLE_CONFIRMED" if (verdict["candidate_is_zero"]
                           and verdict["pilot_reproduced"]
                           and verdict["control_clean"])
                          else "NOT_CONFIRMED")

    artifact = {
        "artifact": "SSZ_B6_STEP1_ARGUMENT_PRINCIPLE_V1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "method": "argument principle on W(w)=det([Q_int|Q_ext]), 6x6",
        "npz_sha256": sha256(NPZ),
        "d_eval_seconds": round(t_eval, 2),
        "contour_per_side": args.per_side,
        "results": results,
        "verdict": verdict,
    }
    out = BRIDGE / "artifacts" / "SSZ_B6_STEP1_ARGUMENT_PRINCIPLE_V1.json"
    out.write_text(json.dumps(artifact, indent=1) + "\n")
    print(f"VERDICT: {verdict['overall']}  -> {out}", flush=True)
    print(f"wall {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
