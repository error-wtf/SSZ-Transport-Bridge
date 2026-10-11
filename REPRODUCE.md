# Reproducing this repo's results

Step-by-step from a clean checkout. All commands were executed on
2026-10-11 on the repository at HEAD `517ef2a` (working tree includes
the 2026-10-11 deliverables: demos, pipeline, figures, README, this
file, `PEER_REVIEW.md`, `docs/FULL_PIPELINE_REPORT_2026-10-11.md`)
with the interpreter and versions below — every timing and file size
quoted here was measured, not estimated.

## 1. Environment

```
python                    3.14.7
numpy                     2.5.3
scipy                     1.18.1
matplotlib                3.11.2
pillow                    12.3.0
mpmath                    1.3.0
pytest                    9.1.1   (venv, dev extra)
ruff                      0.17.0   (lint-clean claim)
```

## 2. Install

```
python -m venv <venv>
<venv>/bin/pip install -e ".[dev]"
<venv>/bin/pip install matplotlib pillow mpmath   # demo + B6 tools
```

The core package needs only `numpy>=1.24` (`pyproject.toml`);
`matplotlib`/`pillow` for the demos/figures, `mpmath` for the B6
Step-0 40-digit analytic controls. `requirements-ssz-backend.txt`
additionally lists `pandas`, `sympy` for the closure backend import
path.

## 3. Tests

```
cd SSZ-Transport-Bridge
<venv>/bin/python -m pytest tests/ -q
```

Expected tail of output:

```
54 passed in ~10 s
```

(measured `54 passed in 8.4s` standalone on 2026-10-11; ~10 s inside
the pipeline where pytest imports are cold).

## 4. Lint

```
<venv>/bin/python -m ruff check .
```

Expected: `All checks passed!` (ruff 0.17.0, 2026-10-11).

## 5. Demos

```
<venv>/bin/python demos/demo_b3b4_ladder.py        # < 5 s
<venv>/bin/python demos/demo_blind_scan.py         # < 5 s
<venv>/bin/python demos/demo_backend_audit.py      # < 15 s (runs audit live)
<venv>/bin/python demos/animate_core_results.py    # ~19 s, 3 GIFs
```

Expected outputs of the animation demo:

```
demos/output/b6_step0_square_well_pole_reveal.gif   292123 bytes
demos/output/bridge_contract_chain.gif              463238 bytes
demos/output/b6_step1_ch2_pilot_pole.gif            125016 bytes
```

(≈ 0.9 MB total; sizes measured 2026-10-11, Pillow GIF encoding —
byte sizes may vary by a few percent across pillow versions.)

## 6. Figures

```
<venv>/bin/python tools/make_figures_v1.py            # ~4 min (GIF encoding)
<venv>/bin/python tools/make_figures_v1.py --no-anim  # PNGs only, ~15 s
```

Writes `figures/fig1_blind_scan_tiles.png`, `fig2_b3b4_ladder_convergence.png`,
`fig3_step0_pole_reproduction.png`, `fig4_metric_coefficients.png` and
`figures/anim1_blind_scan_tiles.gif` (6 frames), `figures/anim2_b3b4_ladder.gif` (4 frames).
Every displayed number is loaded live from `artifacts/*.json` and
`data/frozen_specs/*.npz` — no handwritten physics numbers in the
plotting code.

## 7. Full pipeline

```
<venv>/bin/python tools/run_full_pipeline.py
```

Expected: exit code 0 after ~2.5–3 min with all 7 stages green and
`FULL_PIPELINE_REPORT.md` written to the repo root. `--skip-demo`
omits the GIF stage. Stage order: pytest → anti-circularity audit →
bridge certificate → outgoing-pencil audit → 1-channel Riccati
control → B6 Step-0 gate tool → demos. (Reference run 2026-10-11:
see `docs/FULL_PIPELINE_REPORT_2026-10-11.md`.)

## 8. Optional: individual gates

Each stage's tool can be run standalone; examples:

```
PYTHONPATH=src <venv>/bin/python tools/evaluate_bridge.py                 # REAL_SSZ_ADAPTER_PASS
<venv>/bin/python tools/audit_backend_independence.py       # ANTI-CIRCULARITY AUDIT: PASS
<venv>/bin/python tools/run_b6_step0_analytic_controls.py   # 8/8 gates PASS
<venv>/bin/python tools/run_b3b4_full_3dof_solve.py         # S1–S5 True (box modes only)
```

## 9. What will NOT reproduce

- B6 Step 1 poles: no certified pole exists yet (status OPEN in
  `artifacts/QNM_SPECTROSCOPY_PATH_STATUS_V1.json`); the blind-scan
  and pilot-pole artifacts are records of work in progress, not
  reproducible pass/fail gates.
- Wall-clock times vary by machine; quoted numbers are from the
  2026-10-11 reference run.
- GIF byte sizes depend on the pillow build.
