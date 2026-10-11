#!/usr/bin/env python3
"""demo_b3b4_ladder.py — Kernergebnis: independent 3-DOF box solve (B3/B4).

Zeigt live aus artifacts/SSZ_B3_B4_FULL_3DOF_INDEPENDENT_SOLVE_V1.json:
  - Richardson-Ladder 438→3500 Knoten (12 Box-Moden, 6 doppelt)
  - max. Verschiebung zwischen den zwei feinsten Leveln
  - mediane beobachtete Konvergenzordnung p
  - Richardson-extrapolierte Eigenwerte
  - Gates S1..S5 und ehrliches Verdict (1e-6-Kriterium NICHT erfüllt)

Laufzeit < 5 s (nur JSON-Lesung). Aufruf vom Repo-Root:
    PYTHONPATH=src python demos/demo_b3b4_ladder.py
"""
from __future__ import annotations

import json
from pathlib import Path

BRIDGE = Path(__file__).resolve().parents[1]
ART = BRIDGE / "artifacts"


def main() -> int:
    art = json.loads((ART / "SSZ_B3_B4_FULL_3DOF_INDEPENDENT_SOLVE_V1.json")
                     .read_text(encoding="utf-8"))

    print("=" * 72)
    print("B3/B4 — INDEPENDENT 3-DOF BOX SOLVE (Certification Summary)")
    print("=" * 72)

    lad = art["s4_ladder"]
    print(f"\nRichardson-Ladder (Knotenzahlen): {lad['ladder_node_counts']}")
    print(f"Feinstes Paar: {lad['finest_pair']}")
    print(f"  max_rel_shift_finest      = {lad['max_rel_shift_finest']:.3e}"
          f"   (Schwelle 1e-6: NICHT erreicht)")
    print(f"  min_overlap_finest        = {lad['min_overlap_finest']:.6f}")
    print(f"  median_observed_order p   = {lad['median_observed_order']:.4f}")
    print(f"  strict_1e-6 criterion met = {lad['strict_1e-6_criterion_met']}")

    print("\nRichardson-extrapolierte Eigenwerte (erste 6):")
    lam = [x for x in lad["richardson_lambda_inf"] if x is not None]
    for i, v in enumerate(lam[:6]):
        print(f"  lambda_inf[{i:2d}] = {v:.9f}")

    print("\nGates:")
    for g, ok in art["gates"].items():
        print(f"  {g}: {'PASS' if ok else 'FAIL'}")

    print(f"\nWall time: {art['wall_seconds']:.1f} s")
    print("\nVerdict (wörtlich aus dem Artefakt):")
    print(f"  {art['verdict']}")
    print("\nWichtig: 'Box modes only' — das sind KEINE QNM (B6 offen).")
    print("Quelle: artifacts/SSZ_B3_B4_FULL_3DOF_INDEPENDENT_SOLVE_V1.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
