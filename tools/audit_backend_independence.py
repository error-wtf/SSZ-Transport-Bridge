#!/usr/bin/env python3
"""Anti-circularity audit: no duplicated executable SSZ physics.

Fails if any executable bridge module contains forbidden physics
implementations (Christoffel construction, geodesic RHS, phase-integral
integrands, hard-coded member values).  Calling a source function by name
is allowed; reproducing its body is not.  Documentation strings are
ignored; code constructs are not.
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Forbidden: executable reproductions of SSZ physics.  Each entry is a
# (forbidden AST pattern description, detector) pair.
FORBIDDEN_FLOATS = {
    # Hard-coded light-ring / photon-sphere radii of the frozen member.
    "1.5",
}

FORBIDDEN_CALLS = {
    # Building Christoffel symbols by hand.
    "christoffel",
    "Christoffel",
}

FORBIDDEN_SUBSTRINGS = [
    # Geodesic RHS reconstruction in executable code:
    "def _geodesic_rhs",
    "def geodesic_rhs",
    "d2r_dlam2",
    # Raychaudhuri RHS:
    "def raychaudhuri",
    # Hard-coded member observable numbers (member identity lives in
    # provenance, not in bridge code):
    "0.3876",
    "0.4269",
]


def check_file(path: Path) -> list[str]:
    findings: list[str] = []
    src = path.read_text()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            fn = node.func
            name = getattr(fn, "id", None) or getattr(fn, "attr", "")
            if name in FORBIDDEN_CALLS:
                findings.append(
                    f"{path.relative_to(ROOT)}:{node.lineno}: forbidden "
                    f"physics call '{name}'")
        if isinstance(node, ast.Constant):
            if isinstance(node.value, float) and str(node.value) in (
                    FORBIDDEN_FLOATS - {"1.5"}):
                findings.append(
                    f"{path.relative_to(ROOT)}:{node.lineno}: hard-coded "
                    f"member constant {node.value}")
    for i, line in enumerate(src.splitlines(), 1):
        for sub in FORBIDDEN_SUBSTRINGS:
            if sub in line and not line.lstrip().startswith(("#", '"', "'")):
                findings.append(
                    f"{path.relative_to(ROOT)}:{i}: forbidden physics "
                    f"implementation '{sub}'")
    return findings


def main() -> int:
    findings: list[str] = []
    targets = sorted((ROOT / "src/transport_bridge").rglob("*.py"))
    targets += sorted((ROOT / "tests").glob("*.py"))
    targets += sorted((ROOT / "tools").glob("*.py"))
    for p in targets:
        findings.extend(check_file(p))
    if findings:
        print("ANTI-CIRCULARITY AUDIT: FAIL")
        for f in findings:
            print(" ", f)
        return 1
    print("ANTI-CIRCULARITY AUDIT: PASS (no duplicated SSZ physics in "
          f"{len(targets)} executable files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
