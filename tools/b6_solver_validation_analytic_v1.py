#!/usr/bin/env python3
"""B6 Solver validation against ANALYTICALLY KNOWN resonances (Carmen, 12.10.).

Before the SSZ pole search resumes (SSZ_B6_GATE_RECALIBRATION_V1, step 3),
the recalibrated gate machinery must be proven on cases where the answer is
known in closed form:

  * simple zero   : half-line square well (V0=-2.5, a=1), Neumann family,
                    matching equation  F(w) = -q tan(q a) - i w = 0,
                    q = sqrt(w^2 - V0)  (roots: 100-digit reference file).
  * double zero   : (w0, a) solved from F = 0 AND dF/dw = 0 simultaneously
                    (exceptional point); E ~ C (w - w0)^2, slope m = 2.
  * negative ctrl : off-root point -> gates MUST reject.
  * retro-test    : the supersedet SSZ pilot (SSZ_B6_PILOT_SUPERSESSION_V1)
                    fed through the SAME new gate logic -> must be REJECTED.

E(w) is built with the SAME construction chain as the certified control
solver tools/b6_mpmath_evans_v2.py: interior basis RK4-propagated (real
ODE), matched against the analytic exterior branch at r = a, raw det
(holomorphic by construction), singular values from eig(M^H M).

Artifact: artifacts/SSZ_B6_STEP1_SOLVER_VALIDATION_V1.json
Verdict   VALIDATION_PASSED only if: simple root passes G1-G5 (slope ~ 1),
          double root passes with declared order 2 (slope ~ 2),
          negative control fails, retro-test rejects the old pilot.
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

V0 = -2.5
ART = BRIDGE / "artifacts"
REF = BRIDGE / "data" / "frozen_specs" / "HALF_LINE_SQUARE_WELL_REFERENCE_ROOTS_V2.json"

mp.mp.dps = 50


# --------------------------------------------------------------------------
# Evans construction (same chain as b6_mpmath_evans_v2)
# --------------------------------------------------------------------------
def _evans(w: mp.mpc, a: mp.mpf, n_steps: int = 2000):
    """E(w) = det [[u(a), u'(a)], [e^{iwa}, i w e^{iwa}]]

    u solves u'' + (w^2 - V0) u = 0 on [0, a] with Neumann data
    u(0)=1, u'(0)=0  (center r=0, half line).  Zero of E  <=>  F(w)=0
    (whenever cos(q a) != 0; the determinant is holomorphic everywhere).
    """
    w2 = w * w
    k2 = w2 - V0
    # RK4 propagation of [u, u']
    h = a / n_steps
    u, du = mp.mpf(1), mp.mpf(0)
    r = mp.mpf(0)

    def f(_r, y):
        return mp.matrix([y[1], -k2 * y[0]])

    y = mp.matrix([u, du])
    for _ in range(n_steps):
        k1 = f(r, y)
        k2_ = f(r + h / 2, y + h / 2 * k1)
        k3 = f(r + h / 2, y + h / 2 * k2_)
        k4 = f(r + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2_ + 2 * k3 + k4)
        r += h
    ua, dua = y[0], y[1]
    ea = mp.e**(1j * w * a)
    M = mp.matrix([[ua, dua], [ea, 1j * w * ea]])
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    # singular values via eig(M^H M) (same verified recipe as the control solver)
    H = M.H @ M
    vals, _ = mp.eig(H)
    svals = sorted((mp.sqrt(max(mp.re(v), 0)) for v in vals),
                   key=lambda z: mp.re(z), reverse=True)
    return det, svals


def evaluate(w_c: complex, a: float, n_steps: int = 2000) -> dict:
    wm, am = mp.mpc(w_c), mp.mpf(repr(a))
    det, svals = _evans(wm, am, n_steps)
    return {"abs_E": abs(det), "sigma_min": svals[-1],
            "svals": [float(mp.nstr(s, 12)) for s in svals]}


# --------------------------------------------------------------------------
# Recalibrated gate logic (SSZ_B6_GATE_RECALIBRATION_V1, G1-G5 portion;
# G6/G7 are construction-time invariants, verified once below)
# --------------------------------------------------------------------------
def rosette(w_c: complex, a: float):
    offs = [(1e-3, 0.0), (-1e-3, 0.0), (0.0, 1e-3), (0.0, -1e-3)]
    return [evaluate(w_c + complex(dx, dy), a)["abs_E"] for dx, dy in offs]


def gates(w_c: complex, a: float, declared_order: int = 1) -> dict:
    g = {}
    center = evaluate(w_c, a)
    g["abs_E_center"] = mp.nstr(center["abs_E"], 8)
    g["sigma_min_center"] = mp.nstr(center["sigma_min"], 8)

    # G1 relative dip vs measured floor
    ros = rosette(w_c, a)
    floor = sorted(ros)[1]
    g["measured_floor"] = mp.nstr(floor, 8)
    g["dip_orders"] = round(math.log10(floor / float(center["abs_E"])), 3) \
        if center["abs_E"] > 0 else float("inf")
    g["G1_relative_dip_ge2OD"] = g["dip_orders"] >= 2.0

    # G2 scaling |E(w0+delta)| ~ C |delta|^m over the delta decades
    # (center excluded from the fit: log10(0) undefined; the center only
    # feeds G1/G3 via the dip)
    deltas = [1e-3, 1e-4, 1e-5, 1e-6]
    xs = [math.log10(d) for d in deltas]
    ys = [float(mp.log10(evaluate(w_c + complex(d, 0.0), a)["abs_E"]))
          for d in deltas]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    icpt = my - slope * mx
    ss_res = sum((y - (slope * x + icpt)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    g["scaling_slope"] = round(slope, 4)
    g["scaling_R2"] = round(r2, 6)
    exp_lo, exp_hi = declared_order - 0.25, declared_order + 0.25
    g["G2_scaling_slope_matches_order"] = exp_lo <= slope <= exp_hi and r2 >= 0.99

    # G3 absolute floor: center far below measured floor
    g["G3_below_measured_floor"] = float(center["abs_E"]) < 1e-2 * floor

    # G4 negative control (fixed offset) above measured floor
    nc = evaluate(w_c + complex(0.3, 0.3), a)["abs_E"]
    g["negctrl_abs_E"] = mp.nstr(nc, 8)
    g["G4_negctrl_above_floor"] = nc > 100 * floor

    # G5 dps stability: recompute center at different dps
    old = mp.mp.dps
    mp.mp.dps = old + 20
    alt = evaluate(w_c, a)["abs_E"]
    mp.mp.dps = old
    rel = abs(alt - center["abs_E"]) / center["abs_E"]
    g["dps_stability_rel_change"] = mp.nstr(rel, 4)
    g["G5_dps_stable"] = rel < mp.mpf("0.5")

    g["gate_pass"] = all(g[k] for k in (
        "G1_relative_dip_ge2OD", "G2_scaling_slope_matches_order",
        "G3_below_measured_floor", "G4_negctrl_above_floor", "G5_dps_stable"))
    return g


def main() -> int:
    now = datetime.now(timezone.utc).isoformat()
    ref = json.loads(REF.read_text())
    a = float(ref["params"]["a"])

    # ---- construction invariants (G6/G7 analogues) ------------------------
    w_n1 = complex(*ref["N_roots"][0][:2])
    # product of svals == |det| (verified identity of the control solver)
    wm = mp.mpc(w_n1)
    det, svals = _evans(wm, mp.mpf(repr(a)))
    prod = mp.fprod(svals)
    g7_err = abs(prod - abs(det)) / abs(det)

    # residual of the reference root against the closed form (independent)
    q = mp.sqrt(wm * wm - V0)
    res_ref = abs(-q * mp.tan(q * a) - 1j * wm)

    results = {"audit": "SSZ_B6_STEP1_SOLVER_VALIDATION_V1", "date": now,
               "construction": "RK4-propagated interior + analytic exterior, raw holomorphic det (b6_mpmath_evans_v2 chain)",
               "invariants": {
                   "product_svals_over_absdet_err": mp.nstr(g7_err, 4),
                   "reference_root_closed_form_residual": mp.nstr(res_ref, 4)},
               "cases": {}}

    ok = g7_err < mp.mpf("1e-20") and res_ref < mp.mpf("1e-12")

    # ---- case 1: simple zero ----------------------------------------------
    g1 = gates(w_n1, a, declared_order=1)
    results["cases"]["simple_zero_N1"] = g1
    ok = ok and g1["gate_pass"]

    # ---- case 2: double zero (F = F_w = 0 solved for (w, a)) ---------------
    F = lambda w, av: -mp.sqrt(w * w - V0) * mp.tan(mp.sqrt(w * w - V0) * av) - 1j * w
    Fw = lambda w, av: mp.diff(lambda x: F(x, av), w)
    # ANALYTIC reduction (verified): F=0 with F'=0  =>  tan(qa) = q/iw  and
    # w = -i/a  (eliminate the tangent: differentiate F=0 w.r.t. w and use
    # tan(qa) = -q*... ) => substitute: q^2 = w^2 - V0 = 2.5 - 1/a^2, and the
    # remaining REAL condition is  q*tan(q*a) = -1/a, i.e.
    #   h(a) = sqrt(2.5 - 1/a^2) * tan(a*sqrt(2.5 - 1/a^2)) + 1/a = 0
    q_of = lambda av: mp.sqrt(mp.mpf(2.5) - 1 / (av * av))
    h = lambda av: q_of(av) * mp.tan(av * q_of(av)) + 1 / av
    a0 = mp.findroot(h, (mp.mpf("1.7"), mp.mpf("2.0")))    # secant, real
    w0 = -1j / a0                                          # exact: F=F'=0
    res_d = abs(F(w0, a0))
    res_dw = abs(Fw(w0, a0))
    g2 = gates(complex(w0), float(a0), declared_order=2)
    g2["double_point_w"] = mp.nstr(w0, 20)
    g2["double_point_a"] = mp.nstr(a0, 20)
    g2["residual_F"] = mp.nstr(res_d, 4)
    g2["residual_dF_dw"] = mp.nstr(res_dw, 4)
    results["cases"]["double_zero_exceptional"] = g2
    ok = ok and res_d < mp.mpf("1e-30") and res_dw < mp.mpf("1e-18") \
        and g2["gate_pass"]

    # ---- case 3: negative control point ------------------------------------
    off = complex(w_n1.real + 0.3, w_n1.imag + 0.3)
    g3 = gates(off, a, declared_order=1)
    results["cases"]["negative_control_must_fail"] = {
        k: g3[k] for k in ("dip_orders", "scaling_slope", "gate_pass")}
    ok = ok and not g3["gate_pass"]

    # ---- case 4: retro-test — supersedet SSZ pilot through new gates -------
    # rosette sigma_min (column-normalized, dps 50) from CONTROLSOLVER_VERDICT:
    # pilot center 7.2e-9 with neighbours of the same order (NO dip), from the
    # 5x5 region grid.  Fed through the SAME G1/G2 logic -> must be REJECTED.
    pilot_center = 7.2e-9
    pilot_rosette = [5.601e-9, 3.931e-8, 2.459e-7, 2.252e-11]  # grid neighbours
    pilot_floor = sorted(pilot_rosette)[1]
    pilot_dip = round(math.log10(pilot_floor / pilot_center), 3)
    results["cases"]["retro_old_pilot"] = {
        "dip_orders": pilot_dip,
        "G1_relative_dip_ge2OD": pilot_dip >= 2.0,
        "verdict": "REJECTED (correct: pilot is not a zero)" if pilot_dip < 2.0
                   else "PASSED (INCORRECT - would re-certify the artifact!)"}
    ok = ok and pilot_dip < 2.0

    results["verdict"] = "VALIDATION_PASSED" if ok else "VALIDATION_FAILED"
    results["rule"] = ("Recalibrated gates (G1-G5) PROVEN on analytic known "
                       "zeros incl. a double zero; negative control rejected; "
                       "supersedet pilot rejected.  Solver search may resume "
                       "per SSZ_B6_GATE_RECALIBRATION_V1 step 4.")

    out = ART / "SSZ_B6_STEP1_SOLVER_VALIDATION_V1.json"
    out.write_text(json.dumps(results, indent=2))
    print("verdict:", results["verdict"], "->", out.name)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
