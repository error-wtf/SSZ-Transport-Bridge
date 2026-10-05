# REAL SSZ BRIDGE SPEC

## Purpose

Extend `SSZ-Transport-Bridge` from its current synthetic SSZ proxy to a **read-only adapter of the real `SSZ_FULL_CLOSURE` transport layer**, while preserving the bridge's core law:

> Shared contracts and types; zero shared physics formulas.

This bridge must not rederive or duplicate SSZ equations. The source of truth for SSZ physics remains `SSZ_FULL_CLOSURE`.

## Current three-layer architecture

1. `Sagnac-Reference-Transport`
   - independently validated analytic/numeric reference transport
   - `SAGNAC_REFERENCE_CLOSURE_PASS`
2. `SSZ-Transport-Bridge`
   - generic validation vocabulary and structural contracts
   - current SSZ backend is synthetic only
3. `SSZ_FULL_CLOSURE`
   - real frozen-member geometry and source-free transport
   - `TRUE_FULL_CLOSURE_PASS`

Target:

```text
Sagnac reference
    -> generic validation contracts
    -> read-only real SSZ adapter
    -> bridge-level structural verdict
```

No step may reinterpret Sagnac as SSZ physics or use Sagnac formulas inside the SSZ adapter.

---

## Critical contract correction before real SSZ hookup

The current bridge contract is too Sagnac-shaped in one place:

```python
ObservableSet.directional
```

is interpreted by `check_forward` / `check_parity` as

```text
plus == -minus
```

This is valid for antisymmetric observables such as a signed Sagnac delay, but it is **not a universal direction law**.

For real SSZ, inward/outward reversal can produce any of:

- sign reversal,
- exchange of two directional channels,
- invariance,
- reciprocal transformation,
- no symmetry contract at all for that observable.

Therefore replace the binary `directional/scalar` assumption with an explicit per-observable transformation policy.

Recommended type:

```python
@dataclass(frozen=True)
class ObservableRule:
    name: str
    transform: str   # "odd" | "even" | "swap" | "reciprocal" | "custom" | "none"
    atol: float
    rtol: float = 0.0
```

or an equivalent enum/callable model.

The bridge must never infer a parity law from the mere existence of a `direction` tag.

---

## Backends

Keep all three:

```text
src/transport_bridge/backends/
    sagnac_backend.py
    ssz_synthetic_backend.py
    ssz_closure_backend.py   # NEW
```

Do not delete the synthetic backend. It remains a development/negative-control layer.

### `ssz_closure_backend.py` law

The new backend is an adapter only.

It may:

- import documented public functions from a checked-out `SSZ_FULL_CLOSURE`, or
- consume frozen machine-readable artifacts produced by that repository.

It must not:

- copy the geodesic equations,
- copy Christoffels,
- copy the phase integral,
- reconstruct `f,h` independently,
- fit new parameters,
- change the frozen member,
- silently fall back to a toy formula.

Fail closed if the real closure checkout/artifacts are unavailable or provenance does not match.

---

## Provenance binding

The adapter must bind every run to:

- `SSZ_FULL_CLOSURE` commit SHA,
- frozen member hash,
- `MODEL_LOCK.json` identity,
- `TRUE_FULL_CLOSURE_VERDICT.json` or equivalent verdict artifact,
- relevant dependency-graph node IDs.

Recommended adapter provenance object:

```python
@dataclass(frozen=True)
class BackendProvenance:
    repository: str
    commit_sha: str
    member_sha256: str
    verdict: str
    gate_ids: tuple[str, ...]
```

If hashes disagree, bridge status must be FAIL/CANNOT_EVALUATE, never PASS.

---

# Real SSZ contract set

The real adapter should initially expose a deliberately small subset of `G110-G122` rather than trying to wrap the whole repository at once.

## RB1 — Null forward transport

Source: `G111`, `null_geodesic_transport(...)`

Physics remains in `SSZ_FULL_CLOSURE`:

\[
k^\nu \nabla_\nu k^\mu = 0.
\]

Bridge consumes only returned diagnostics, e.g.

```text
max_dE
max_dL
max_norm_res
trajectory metadata
```

Bridge contract:

- results finite,
- declared conserved quantities stay within source-owned tolerances,
- adapter provenance matches frozen member.

The bridge does **not** recompute the geodesic residual.

## RB2 — Matter forward transport

Source: `G110`.

Same rule: consume source-owned diagnostics only.

This provides a second real transport family and prevents the bridge from being optical-only.

## RB3 — Phase/readout chain

Source: `G115`, `eikonal_phase_and_redshift(...)`.

Expose separately:

```text
reduced_phase_per_E
shapiro_dt_per_E
redshift_a_to_b
```

Do not collapse them into a single generic scalar.

Structural bridge statement:

```text
transport state -> declared phase/timing/redshift observables
```

No Sagnac formula is applicable here.

## RB4 — Direction reversal / path reversal

For the current static spherical SSZ slice, do **not** invent a `frame_sign` analogue.

Use physically meaningful path reversal only where supported.

Candidate first contract:

Radial phase integral between two radii has magnitude symmetry under endpoint reversal:

\[
S_r(a\to b) = -S_r(b\to a)
\]

if the adapter/source exposes an oriented integral, while redshift obeys reciprocal behavior:

\[
z_{a\to b}\, z_{b\to a} = 1
\]

for the source's static-observer frequency ratio convention.

Important: if the current `eikonal_phase_and_redshift` implementation canonicalizes endpoints using `min/max`, then the adapter must **not pretend** it has orientation information. Either:

1. add an orientation-preserving public result in `SSZ_FULL_CLOSURE`, or
2. mark oriented phase parity `OPEN/NOT_EXPOSED` in the bridge.

Do not reconstruct the sign in the bridge.

## RB5 — Convergence / robustness

Source: `G122` and source-owned robustness artifacts.

The bridge should consume convergence evidence rather than rerun a copied SSZ solver.

Allowed checks:

- verify monotonic/refinement metadata declared by the source,
- verify the artifact belongs to the same frozen member and commit,
- run the source repository's public command if available and parse its result.

Do not impose Sagnac's observed-order threshold on SSZ unless the SSZ source itself declares that numerical method/order.

## RB6 — Negative-control discrimination

Source: `G119`.

Expose at least:

- artificial four-force detection,
- wrong phase integrand detection,
- wrong redshift law detection.

Bridge-level statement:

> the real backend carries source-owned falsifiers and the adapter preserves their PASS/FIRE status without reinterpretation.

## RB7 — Known-limit anchors

Source: `G121`.

Use source-owned Schwarzschild / Newtonian / gamma-limit evidence as a structural known-limit contract.

Again: no equations duplicated in the bridge.

## RB8 — Inverse contract

Do **not** force an inverse parameter just because Sagnac has one.

For the first real SSZ adapter, set:

```text
inverse_status = NOT_YET_DEFINED
```

unless there is a genuinely source-defined inverse problem with a frozen target parameter.

The synthetic backend may continue testing generic inversion mechanics.

Real spectroscopy inversion belongs to a later layer and must not be faked with `frame_sign`.

---

# Proposed bridge verdicts

Do not create `G150` yet.

Use a separate bridge namespace:

```text
BR-R1 provenance binding
BR-R2 real null forward transport
BR-R3 real matter forward transport
BR-R4 real phase/readout exposure
BR-R5 supported reversal symmetry
BR-R6 source-owned convergence evidence
BR-R7 source-owned falsifier preservation
BR-R8 known-limit anchors
BR-R9 no duplicated SSZ physics
BR-R10 inverse explicitly OPEN or real
```

Suggested statuses:

```text
REAL_SSZ_ADAPTER_PASS
BRIDGE_CLASS_MEMBERSHIP_PASS
REAL_SSZ_INVERSE_OPEN
```

A final bridge summary may be:

```text
REAL_SSZ_TRANSPORT_BRIDGE_PASS_WITH_INVERSE_OPEN
```

only if all mandatory forward/provenance/falsifier gates pass.

This verdict must not alter `TRUE_FULL_CLOSURE_PASS` and must not imply empirical validation of SSZ.

---

# File layout

```text
src/transport_bridge/
    contracts.py
    provenance.py
    rules.py
    backends/
        sagnac_backend.py
        ssz_synthetic_backend.py
        ssz_closure_backend.py

tests/
    test_contract_rules.py
    test_sagnac_backend.py
    test_ssz_synthetic_backend.py
    test_ssz_closure_backend.py
    test_real_ssz_provenance.py
    test_real_ssz_no_formula_duplication.py
    test_cross_backend_validation_class.py

tools/
    evaluate_bridge.py
    audit_backend_independence.py

artifacts/
    REAL_SSZ_BRIDGE_VERDICT.json
    backend_provenance.json
    bridge_dependency_graph.json

docs/
    REAL_SSZ_BRIDGE_SPEC.md
    VALIDATION_CLASS_SEMANTICS.md
```

---

# Anti-circularity / independence audit

Add a dedicated audit that fails if `SSZ-Transport-Bridge` contains forbidden SSZ physics strings/formulas beyond documentation or symbol names.

Examples to flag in executable bridge code:

```text
Christoffel construction
Gamma_tensor implementation
W(u)=u^2 f(u) implementation
1/sqrt(f*h) numerical integration
geodesic RHS equations
Raychaudhuri RHS
hard-coded light-ring values
hard-coded member observables
```

The bridge may call a source function by name; it may not reproduce its body.

Also audit that:

- Sagnac backend imports only Sagnac reference package/artifacts,
- SSZ backend imports only SSZ closure package/artifacts,
- neither backend imports the other's formulas,
- shared `contracts.py` contains no backend-specific formulas.

---

# Important semantic correction: validation class

Passing the same structural contract means only:

> the systems can be evaluated by the same abstract validation vocabulary.

It does **not** mean:

- same equations,
- same causal mechanism,
- same symmetry group,
- same observable set,
- same physical interpretation,
- Sagnac validates SSZ physics.

Document this prominently.

---

# First implementation milestone

Implement only:

1. contract-rule refactor (odd/even/swap/reciprocal/none),
2. provenance loader for `SSZ_FULL_CLOSURE`,
3. real read-only `G111` null transport adapter,
4. real read-only `G115` phase/redshift adapter,
5. source-owned `G119/G121/G122` evidence ingestion,
6. independence audit,
7. machine-readable bridge verdict.

Leave:

- real inverse reconstruction,
- spectroscopy,
- QNM mode fitting,
- astronomy data,
- JIF inversion

explicitly OPEN.

---

# Acceptance criteria for milestone 1

A milestone PASS requires all of the following:

- `Sagnac-Reference-Transport` remains untouched/read-only.
- `SSZ_FULL_CLOSURE` remains the sole owner of SSZ equations.
- real SSZ adapter binds to the actual frozen member and repo commit.
- at least G111 and G115 are consumed from the real closure path.
- source-owned G119 negative controls remain visible and correctly classified.
- source-owned G121/G122 evidence is provenance-bound.
- no synthetic `frame_sign` result is reported as a real SSZ result.
- unsupported symmetry/inverse contracts report OPEN, not fabricated PASS.
- shared contract tests do not require Sagnac-specific odd parity for all backends.
- independence audit finds zero duplicated executable SSZ physics.

Only then may the bridge claim:

```text
REAL_SSZ_ADAPTER_PASS
```

and

```text
BRIDGE_CLASS_MEMBERSHIP_PASS
```

with inverse/spectroscopy still OPEN.
