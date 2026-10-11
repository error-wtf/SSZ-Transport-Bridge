#!/usr/bin/env python3
"""demo_blind_scan.py — Kernergebnis: B6 Step-1 Blind Scan (ch2-Pilot).

Zeigt live aus artifacts/SSZ_B6_STEP1_BLIND_SCAN_V1.json und
artifacts/SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json:
  - alle 3 Tiles des deklarierten Scans mit ihren argmin-Punkten
  - sigma_min_global je Tile
  - den zertifizierten ch2-Pilot-Pol (Näherung) aus dem Pilot-Artefakt

Laufzeit < 5 s. Aufruf vom Repo-Root:
    PYTHONPATH=src python demos/demo_blind_scan.py
"""
from __future__ import annotations

import json
from pathlib import Path

BRIDGE = Path(__file__).resolve().parents[1]
ART = BRIDGE / "artifacts"


def main() -> int:
    scan = json.loads((ART / "SSZ_B6_STEP1_BLIND_SCAN_V1.json")
                      .read_text(encoding="utf-8"))
    pilot = json.loads((ART / "SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json")
                       .read_text(encoding="utf-8"))

    print("=" * 72)
    print("B6 STEP-1 — BLIND SCAN tiles + ch2 PILOT POLE")
    print("=" * 72)

    print(f"\nScan-Grid: {scan['grid']}")

    print(f"\n{'Tile':<22}{'argmin (Re, Im)':<28}{'sigma_min_global'}")
    for name, t in scan["tiles"].items():
        am = t["argmin_omega"]
        print(f"{name:<22}({am[0]:>6.2f}, {am[1]:>7.2f})   "
              f"{t['sigma_min_global']:.6f}")

    pole = pilot["pole_qnm_sign"]
    scan_sign = pilot["pole_scan_sign"]
    print("\nch2 Pilot-Pol (aus Pilot-Artefakt, supersedes Blind-Scan-Niveau):")
    print(f"  status          : {pilot['status']}")
    print(f"  pole_qnm_sign   : re={pole['re']:.9f}  im={pole['im']:.9f}")
    print(f"  pole_scan_sign  : re={scan_sign['re']:.9f}  im={scan_sign['im']:.9f}")
    print(f"  certification_gates:")
    for k, v in pilot["certification_gates"].items():
        print(f"    {k}: {v}")

    print("\nHonest caveats (wörtlich):")
    for c in pilot.get("honest_caveats", []):
        print(f"  - {c}")

    print("\nQuellen:")
    print("  artifacts/SSZ_B6_STEP1_BLIND_SCAN_V1.json")
    print("  artifacts/SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
