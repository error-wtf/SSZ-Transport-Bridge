# Peer review guide

What an independent reviewer should check, where the load-bearing
claims live, and what this repo does and does NOT claim. Written
2026-10-11; all file references are against the working tree at HEAD
`517ef2a` plus the 2026-10-11 deliverables (demos/, figures/,
tools/make_figures_v1.py, tools/run_full_pipeline.py, README.md,
REPRODUCE.md, docs/FULL_PIPELINE_REPORT_2026-10-11.md, this file).

## Entry points (read in this order)

1. `README.md` — the validation-class idea, the status banner, and
   embedded figures with captions citing their artifact sources.
2. `docs/REAL_SSZ_BRIDGE_SPEC.md` — ADR-001 and the closure-backend
   contract (independent re-solver semantics).
3. `src/transport_bridge/contracts.py` — the five typed contracts;
   everything else is checked against this vocabulary.
4. `tools/run_full_pipeline.py` + `docs/FULL_PIPELINE_REPORT_2026-10-11.md`
   — one-shot reproduction of every green claim (tests, audit, gates).
5. `REPRODUCE.md` — clean-checkout reproduction steps.

## Where the core claims live

| Claim | Primary evidence | Cross-check |
|---|---|---|
| Two backends pass the SAME structural contracts | `tests/test_bridge_contracts.py` | `artifacts/BRIDGE_CERTIFICATE_V1.json` |
| SSZ physics re-derived independently (no circularity) | `tools/audit_backend_independence.py` (AST-based) | `src/transport_bridge/backends/ssz_closure_backend.py` `INDEPENDENT_PHYSICS` registry |
| B6 Step 0 machinery certified (8/8 gates) | `tools/run_b6_step0_analytic_controls.py` | `artifacts/SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json` (residuals ~1e-13) |
| Blind reproduction of square-well reference poles | same artifact | `data/frozen_specs/HALF_LINE_SQUARE_WELL_REFERENCE_ROOTS_V2.json` |
| 54/54 tests green | `python -m pytest tests/ -q` | `docs/FULL_PIPELINE_REPORT_2026-10-11.md` §1 |
| Lint clean | `python -m ruff check .` | ruff 0.17.0, 2026-10-11 |
| Figures show real artifact data only | `tools/make_figures_v1.py` loads every number from `artifacts/*.json` / `data/frozen_specs/*.npz` | figure captions in README |

## Negative controls to inspect (the repo fails closed — verify that)

- `tools/audit_backend_independence.py` hard-codes two forbidden member
  observables (`0.3876`, `0.4269`) and fails if any executable file
  contains them.
- Unknown-parameter perturbation must fail the inverse contract
  (`tests/test_bridge_contracts.py`).
- B6 Step-0 negative controls separate the Neumann/Dirichlet root
  families — see the `negative_controls` block in
  `SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json`.
- The B6 audit was run against the TRUE state on 2026-10-11: the B6/B2
  solver files were added (`independent_jost_v1.py` in commit `8aa4cff`,
  `bvp_qnm_v2.py`/`stable_evans_v3.py` in commit `46253bc`) without
  being registered in `DECLARED_PHYSICS_FILES`, so the anti-circularity
  audit exited 1 at HEAD `517ef2a`.  The fix applied 2026-10-11 was
  ADR-001-conform: declare the three files in the registry (they import
  no closure code — check A passes; they exist precisely as independent
  re-solvers for cross-validation).
  A reviewer should re-run the audit and confirm no closure imports.

## Known limitations / blockers

- **B6 Step 1 is OPEN.**  The blind scan
  (`SSZ_B6_STEP1_BLIND_SCAN_V1.json`) found no certified 3-DOF QNM
  pole.  The CH2 pilot pole
  (`SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json`) carries status
  `CANDIDATE_POLE_CERTIFIED_NUMERICAL`, but only as a property of the
  **decoupled ch2 subsystem** of the frozen 3x3 operator — read its
  `honest_caveats` block (channel labeling not yet cross-checked,
  full 3-channel coupled certification outstanding).  No pole of the
  physical 3-DOF problem is claimed anywhere in this repo.
- The anti-circularity audit is AST-heuristic (nested `rhs` functions,
  `np.trapezoid` calls, `sqrt`-of-spline-division patterns).  It
  enforces *declared-ness*, not *correctness* of the re-derived
  equations — that is the B2/B7 cross-solver comparison's job.
- The synthetic SSZ backend is a development control (toy
  frame-dragging proxy), not physics.
- Real-astronomy observables are out of scope: the SSZ inverse problem
  is deferred (`inverse=""`).
- GIF/PNG "readability" was verified programmatically (PIL content-band
  statistics, `identify` frame counts), not by human eye.
- B3/B4 box modes: the strict 1e-6 ladder criterion is NOT met
  (`max_rel_shift_finest = 1.14e-5`, `SSZ_B3_B4_FULL_3DOF_INDEPENDENT_SOLVE_V1.json`);
  the artifact's own verdict states box-modes-only scope.  These are
  NOT quasi-normal modes.

## What this repo does NOT claim

- No claim of any SSZ/QNM resonance frequency of a physical object.
- No claim that the closure backend's numerics are correct — only that
  they are independent (declared, audited) and cross-checked by
  structural contracts.
- No claim that the validation-class membership extends beyond the
  three backends listed.
- No performance claims beyond the measured pipeline runtime.
