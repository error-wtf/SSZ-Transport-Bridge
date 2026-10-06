#!/usr/bin/env python3
"""Anti-circularity audit — ADR-001 semantics (independent re-solver).

The bridge backend is an INDEPENDENT RE-SOLVER (Variante B): it re-derives
SSZ transport physics from frozen metric data.  Duplication is therefore not
the failure mode — UNDECLARED duplication is.

This audit enforces, on every executable bridge file:

  A. No closure-code imports: nothing may import from ssz_p5 or any
     SSZ_FULL_CLOSURE source module.  Frozen DATA (MODEL_LOCK, member CSV)
     is the only permitted coupling, always via provenance.

  B. Registry completeness: every file that defines an ODE RHS over metric
     splines, a metric quadrature, or a metric ratio MUST declare those
     routines in `INDEPENDENT_PHYSICS` (backends/ssz_closure_backend.py).
     Detection is AST-based, not string-based:
       * nested function named `rhs` (or `_geodesic_rhs`, `geodesic_rhs`),
       * calls to `np.trapezoid`/`np.trapz` inside physics helpers,
       * `sqrt` calls whose argument contains a division over spline calls.
     Each detection must map to a declared registry key or the audit FAILS.

  C. Hard-coded member observables remain forbidden (0.3876, 0.4269).

Deliberate scope: this audit does NOT judge the correctness of the
re-derived equations — that is the cross-solver comparison's job
(B2/B7 contracts).  It judges that the independence claim is DECLARED and
COMPLETE.
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FORBIDDEN_SUBSTRINGS = [
    "0.3876",
    "0.4269",
]

CLOSURE_IMPORT_PREFIXES = ("ssz_p5",)

# Files allowed to carry declared independent physics:
DECLARED_PHYSICS_FILES = {
    "src/transport_bridge/backends/ssz_closure_backend.py",
}


def closure_imports(tree, rel):
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if any(a.name.startswith(p) for p in CLOSURE_IMPORT_PREFIXES):
                    out.append(f"{rel}:{node.lineno}: closure-code import '{a.name}'")
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if any(mod.startswith(p) for p in CLOSURE_IMPORT_PREFIXES):
                out.append(f"{rel}:{node.lineno}: closure-code import from '{mod}'")
    return out


def physics_patterns(tree, rel):
    """AST-based detection of physics-like numerics (declared-or-fail)."""
    found = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name in ("rhs", "geodesic_rhs", "_geodesic_rhs"):
                found.append({"kind": "ode_rhs", "name": node.name,
                              "line": node.lineno})
        if isinstance(node, ast.Call):
            fn = getattr(node.func, "attr", None) or getattr(node.func, "id", "")
            if fn in ("trapezoid", "trapz"):
                found.append({"kind": "quadrature", "name": fn,
                              "line": node.lineno})
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "sqrt":
            src = ast.unparse(node)
            if "/" in src and "_spline" in src:
                found.append({"kind": "metric_ratio_sqrt", "name": "sqrt(ratio)",
                              "line": node.lineno})
    return found


def declared_registry():
    backend = ROOT / "src/transport_bridge/backends/ssz_closure_backend.py"
    tree = ast.parse(backend.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if getattr(t, "id", "") == "INDEPENDENT_PHYSICS":
                    return ast.literal_eval(node.value)
    return {}


def main() -> int:
    findings = []
    registry = declared_registry()
    if not registry:
        findings.append("INDEPENDENT_PHYSICS registry missing/empty in "
                        "backends/ssz_closure_backend.py")

    targets = sorted((ROOT / "src/transport_bridge").rglob("*.py"))
    targets += sorted((ROOT / "tests").glob("*.py"))
    targets += sorted((ROOT / "tools").glob("*.py"))

    total_patterns = 0
    for p in targets:
        rel = str(p.relative_to(ROOT))
        src = p.read_text()
        tree = ast.parse(src)

        findings.extend(closure_imports(tree, rel))

        skip_substrings = rel == "tools/audit_backend_independence.py"
        for line_no, line in enumerate(src.splitlines(), 1):
            if skip_substrings:
                break  # this file only lists the forbidden constants
            for sub in FORBIDDEN_SUBSTRINGS:
                if sub in line and not line.lstrip().startswith(("#", '"', "'")):
                    findings.append(f"{rel}:{line_no}: hard-coded member observable")

        patterns = physics_patterns(tree, rel)
        if patterns and rel not in DECLARED_PHYSICS_FILES:
            for pt in patterns:
                findings.append(
                    f"{rel}:{pt['line']}: undeclared physics numerics "
                    f"({pt['kind']}) — file not in DECLARED_PHYSICS_FILES")
            continue
        total_patterns += len(patterns)

    declared_keys = set(registry)
    backend_rel = "src/transport_bridge/backends/ssz_closure_backend.py"
    btree = ast.parse((ROOT / backend_rel).read_text())
    bpatterns = physics_patterns(btree, backend_rel)
    required = set()
    for pt in bpatterns:
        if pt["kind"] == "ode_rhs":
            required.add("null_geodesic_rhs")
            required.add("null_conserved_quantities")
        elif pt["kind"] == "quadrature":
            required.add("phase_integral")
        elif pt["kind"] == "metric_ratio_sqrt":
            required.add("redshift_ratio")
    missing = required - declared_keys
    if missing:
        findings.append(
            f"registry incomplete: patterns require {sorted(required)}, "
            f"missing {sorted(missing)} from INDEPENDENT_PHYSICS")

    if findings:
        print("ANTI-CIRCULARITY AUDIT: FAIL")
        for f in findings:
            print(" ", f)
        return 1
    print(f"ANTI-CIRCULARITY AUDIT: PASS (ADR-001 independent re-solver "
          f"semantics; {total_patterns} declared physics patterns in "
          f"{len(targets)} executable files; registry keys: "
          f"{sorted(declared_keys)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
