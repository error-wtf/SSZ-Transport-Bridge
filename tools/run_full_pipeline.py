"""Full pipeline: tests, gates, audits, demo — one command, one report.

Runs, in order:
  1. pytest suite (counts PASS/FAIL per module)
  2. tools/audit_backend_independence.py
  3. tools/build_bridge_certificate.py
  4. tools/audit_outgoing_pencil.py
  5. tools/check_riccati_1ch_control.py
  6. tools/run_b6_step0_analytic_controls.py (8/8 gates)
  7. demos/animate_core_results.py (regenerates the three GIFs)

Then writes FULL_PIPELINE_REPORT.md at the repo root with the date,
HEAD hash, test/gate/audit results, artifact numbers, and an honest
OPEN/BLOCKED list.  Exit 0 only if every stage is green.

Usage:
    python tools/run_full_pipeline.py [--skip-demo]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
import time
from pathlib import Path

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

REPORT = BRIDGE / "FULL_PIPELINE_REPORT.md"

TOOL_ORDER = [
    ("audit_backend_independence", "tools/audit_backend_independence.py"),
    ("build_bridge_certificate", "tools/build_bridge_certificate.py"),
    ("audit_outgoing_pencil", "tools/audit_outgoing_pencil.py"),
    ("check_riccati_1ch_control", "tools/check_riccati_1ch_control.py"),
    ("run_b6_step0_analytic_controls", "tools/run_b6_step0_analytic_controls.py"),
]

PASS_RE = re.compile(r"(\d+) passed")
ARTIFACT_NUMBERS = {
    "SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json": [
        ("status", None, ""),
        ("g2_analytic_residual_max", "channels.ch1", " [ch1]"),
        ("g2_analytic_residual_max", "channels.ch2", " [ch2]"),
        ("g2_analytic_residual_max", "channels.ch3", " [ch3]"),
    ],
}


def _run(cmd: list[str], cwd: Path = BRIDGE) -> tuple[int, str, float]:
    t0 = time.monotonic()
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return proc.returncode, (proc.stdout + proc.stderr).strip(), time.monotonic() - t0


def _dotted(obj, path: str):
    cur = obj
    for part in path.split("."):
        cur = cur[part]
    return cur


def collect_artifact_numbers() -> list[tuple[str, str, str]]:
    rows = []
    for fname, keys in ARTIFACT_NUMBERS.items():
        path = BRIDGE / "artifacts" / fname
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            rows.append((fname, "<missing/unreadable>", ""))
            continue
        for key, pathstr, suffix in keys:
            try:
                val = data[key] if pathstr is None else _dotted(data, pathstr)[key]
            except (KeyError, TypeError):
                continue
            rows.append((fname, key, f"{val}{suffix}"))
    pilot = BRIDGE / "artifacts" / "SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json"
    if pilot.exists():
        data = json.loads(pilot.read_text(encoding="utf-8"))
        qnm = data["pole_qnm_sign"]
        rows.append((pilot.name, "pole_qnm_sign",
                     f"re={qnm['re']:.5f}, im={qnm['im']:.5f}"))
        rows.append((pilot.name, "status", data["status"]))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-demo", action="store_true",
                        help="skip GIF regeneration (report notes it)")
    args = parser.parse_args()

    py = sys.executable
    stage_results: list[tuple[str, bool, str, float]] = []

    # ---- stage 1: pytest ------------------------------------------------
    print("[1/7] pytest suite", flush=True)
    code, out, secs = _run([py, "-m", "pytest", "tests/", "-q"])
    m = PASS_RE.search(out)
    n_pass = int(m.group(1)) if m else 0
    fail_m = re.search(r"(\d+) failed", out)
    n_fail = int(fail_m.group(1)) if fail_m else 0
    tests_green = code == 0 and n_fail == 0
    tests_line = f"{n_pass} passed, {n_fail} failed ({secs:.1f}s)"
    print(f"      {tests_line}", flush=True)
    stage_results.append(("pytest", tests_green, tests_line, secs))

    # ---- stages 2..6: tools ---------------------------------------------
    for i, (name, rel) in enumerate(TOOL_ORDER, start=2):
        print(f"[{i}/7] {name}", flush=True)
        code, out, secs = _run([py, str(BRIDGE / rel)])
        green = code == 0
        last = out.splitlines()[-1] if out else "(no output)"
        print(f"      exit={code}  {last}", flush=True)
        stage_results.append((name, green, f"exit {code}: {last}", secs))

    # ---- stage 7: demo ----------------------------------------------------
    if args.skip_demo:
        print("[7/7] demo skipped (--skip-demo)", flush=True)
        stage_results.append(("demo", True, "skipped via --skip-demo", 0.0))
        demo_files: list[Path] = sorted((BRIDGE / "demos" / "output").glob("*.gif")) \
            if (BRIDGE / "demos" / "output").exists() else []
        demo_note = "skipped (--skip-demo)"
    else:
        print("[7/7] demos/animate_core_results.py", flush=True)
        code, out, secs = _run([py, str(BRIDGE / "demos" / "animate_core_results.py")])
        green = code == 0
        print(f"      exit={code}", flush=True)
        stage_results.append(("demo", green, f"exit {code}", secs))
        demo_files = sorted((BRIDGE / "demos" / "output").glob("*.gif"))
        demo_note = "regenerated"

    # ---- HEAD + artifact numbers -----------------------------------------
    _, head, _ = _run(["git", "rev-parse", "--short", "HEAD"])
    head = head or "(git unavailable)"
    art_rows = collect_artifact_numbers()

    all_green = all(g for _, g, _, _ in stage_results)

    # ---- report ------------------------------------------------------------
    now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# Full pipeline report",
        "",
        f"* generated: **{now}**",
        f"* HEAD: `{head}`",
        f"* interpreter: `{py}`",
        f"* overall: **{'ALL GREEN' if all_green else 'NOT GREEN — see stages'}**",
        "",
        "## Stages",
        "",
        "| # | stage | green | result | time |",
        "|---|-------|-------|--------|------|",
    ]
    for i, (name, green, res, secs) in enumerate(stage_results, start=1):
        lines.append(f"| {i} | {name} | {'yes' if green else 'NO'} | {res} | {secs:.1f}s |")
    lines += [
        "",
        "## Artifact numbers (loaded live)",
        "",
        "| artifact | key | value |",
        "|----------|-----|-------|",
    ]
    for fname, key, val in art_rows:
        lines.append(f"| `{fname}` | `{key}` | {val} |")
    lines += [
        "",
        "## Demo GIFs (demos/output/)",
        "",
    ]
    if demo_files:
        for f in demo_files:
            lines.append(f"* `{f.name}` — {f.stat().st_size} bytes ({demo_note})")
    else:
        lines.append("* (none — demo was skipped or failed)")
    lines += [
        "",
        "## OPEN / BLOCKED",
        "",
        "### OPEN (honest list, next work)",
        "",
        "* **B6 Step 1 — physical 3-DOF QNM problem is NOT certified.**",
        "  Step 0 certifies the independent machinery (square-well",
        "  known-answer controls) only. The ch2 pilot pole",
        "  (`SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json`) is a single-channel",
        "  pilot, not the 3-DOF certification; ch1/ch3 and the coupled",
        "  problem remain open.",
        "* **Inverse contract for SSZ remains deferred** (`inverse=\"\"`):",
        "  no SSZ inverse round-trip claim is made anywhere in this repo.",
        "* **SSZ synthetic backend is a development control**, not a",
        "  physical claim; it is intentionally kept as negative control.",
        "",
        "### BLOCKED",
        "",
        "* Nothing in this pipeline is blocked; all runnable stages ran.",
        "",
        "### Not claimed",
        "",
        "* No statement about real astrophysical spectroscopy is made or",
        "  implied. The bridge proves validation-class membership between",
        "  backends under the five structural contracts, nothing more.",
        "",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nreport -> {REPORT.name}  (overall: "
          f"{'ALL GREEN' if all_green else 'NOT GREEN'})")
    return 0 if all_green else 1


if __name__ == "__main__":
    raise SystemExit(main())
