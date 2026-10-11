#!/usr/bin/env python3
"""demo_backend_audit.py — Kernergebnis: Anti-Zirkularitäts-Audit (ADR-001).

Führt die Deklarationsprüfung live aus (dieselbe Logik wie
tools/audit_backend_independence.py, als Bibliotheksaufruf) und zeigt
zusätzlich die Verdict-Kette aus REAL_SSZ_BRIDGE_VERDICT.json
(B9_INDEPENDENCE etc.).

Laufzeit < 10 s. Aufruf vom Repo-Root:
    PYTHONPATH=src python demos/demo_backend_audit.py
"""
from __future__ import annotations

import json
from pathlib import Path

BRIDGE = Path(__file__).resolve().parents[1]
ART = BRIDGE / "artifacts"

# dieselbe Registry wie tools/audit_backend_independence.py (nach Patch)
DECLARED_PHYSICS_FILES = [
    "src/transport_bridge/independent_jost_v1.py",
    "src/transport_bridge/bvp_qnm_v2.py",
    "src/transport_bridge/stable_evans_v3.py",
]


def main() -> int:
    print("=" * 72)
    print("ANTI-CIRCULARITY / BACKEND-INDEPENDENCE (ADR-001)")
    print("=" * 72)

    import subprocess, sys
    r = subprocess.run(
        [sys.executable, str(BRIDGE / "tools" / "audit_backend_independence.py")],
        capture_output=True, text=True, cwd=BRIDGE,
    )
    print("\n[LIVE-RUN tools/audit_backend_independence.py]")
    print(r.stdout.strip() or r.stderr.strip())
    print(f"exit code: {r.returncode}")

    print("\n[Verdict-Kette aus REAL_SSZ_BRIDGE_VERDICT.json]")
    v = json.loads((ART / "REAL_SSZ_BRIDGE_VERDICT.json").read_text())
    print(f"  VERDICT              : {v['VERDICT']}")
    print(f"  BRIDGE_CLASS_MEMBERSHIP: {v['BRIDGE_CLASS_MEMBERSHIP']}")
    b9 = v["B_RESULTS"]["B9_INDEPENDENCE"]
    print(f"  B_RESULTS.B9_INDEPENDENCE : {b9}")
    print(f"  B_RESULTS.B8_INVERSE     : "
          f"{v['B_RESULTS'].get('B8_INVERSE', 'n/a')}")
    for k in ("G119", "G121", "G122"):
        if k in v.get("B_RESULTS", {}):
            print(f"  {k}: {v['B_RESULTS'][k]}")

    print("\nDeklarierte unabhängige Physics-Backends (Registry):")
    for f in DECLARED_PHYSICS_FILES:
        p = BRIDGE / f
        print(f"  {'OK ' if p.exists() else 'MISSING'} {f}")
    return 0 if r.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
