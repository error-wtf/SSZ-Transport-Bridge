# Full pipeline report

* generated: **2026-10-11 15:27 UTC**
* HEAD: `517ef2a`
* interpreter: `/home/error/physics/venv/bin/python`
* overall: **ALL GREEN**

## Stages

| # | stage | green | result | time |
|---|-------|-------|--------|------|
| 1 | pytest | yes | 54 passed, 0 failed (6.9s) | 6.9s |
| 2 | audit_backend_independence | yes | exit 0: ANTI-CIRCULARITY AUDIT: PASS (ADR-001 independent re-solver semantics; 5 declared physics patterns in 38 executable files; registry keys: ['null_conserved_quantities', 'null_geodesic_rhs', 'phase_integral', 'redshift_ratio']) | 0.1s |
| 3 | build_bridge_certificate | yes | exit 0: verdict: BRIDGE_CERTIFICATE_PASS | 0.5s |
| 4 | audit_outgoing_pencil | yes | exit 0: written: /home/error/physics/SSZ-Transport-Bridge/artifacts/SSZ_B6_STEP1_PENCIL_AUDIT_V1.json  (5 frequency samples audited) | 0.3s |
| 5 | check_riccati_1ch_control | yes | exit 0: RICCATI_1CH_CONTROL: PASS | 7.4s |
| 6 | run_b6_step0_analytic_controls | yes | exit 0: verdict: B6_STEP0_CERTIFIED | 122.8s |
| 7 | demo | yes | exit 0 | 18.7s |

## Artifact numbers (loaded live)

| artifact | key | value |
|----------|-----|-------|
| `SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json` | `status` | B6_STEP0_CERTIFIED |
| `SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json` | `g2_analytic_residual_max` | 4.0983210015919556e-14 [ch1] |
| `SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json` | `g2_analytic_residual_max` | 4.511063047890465e-14 [ch2] |
| `SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json` | `g2_analytic_residual_max` | 8.648883562053639e-14 [ch3] |
| `SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json` | `pole_qnm_sign` | re=0.22643, im=-0.08455 |
| `SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json` | `status` | CANDIDATE_POLE_CERTIFIED_NUMERICAL |

## Demo GIFs (demos/output/)

* `b6_step0_square_well_pole_reveal.gif` — 292123 bytes (regenerated)
* `b6_step1_ch2_pilot_pole.gif` — 125016 bytes (regenerated)
* `bridge_contract_chain.gif` — 463238 bytes (regenerated)

## OPEN / BLOCKED

### OPEN (honest list, next work)

* **B6 Step 1 — physical 3-DOF QNM problem is NOT certified.**
  Step 0 certifies the independent machinery (square-well
  known-answer controls) only. The ch2 pilot pole
  (`SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json`) is a single-channel
  pilot, not the 3-DOF certification; ch1/ch3 and the coupled
  problem remain open.
* **Inverse contract for SSZ remains deferred** (`inverse=""`):
  no SSZ inverse round-trip claim is made anywhere in this repo.
* **SSZ synthetic backend is a development control**, not a
  physical claim; it is intentionally kept as negative control.

### BLOCKED

* Nothing in this pipeline is blocked; all runnable stages ran.

### Not claimed

* No statement about real astrophysical spectroscopy is made or
  implied. The bridge proves validation-class membership between
  backends under the five structural contracts, nothing more.

