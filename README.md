# SSZ-Transport-Bridge

Generic typed-contract bridge between transport systems.  Two read-only
backends, one small shared vocabulary — **zero shared formulas**.

## The core rule

The bridge shares TYPES and CHECK STRUCTURE between backends, never
physics.  Each backend owns its physics completely; the bridge owns
only the vocabulary in which five structural contracts are stated:

```
forward  →  convergence  →  direction/parity  →  observable  →  inverse
```

A backend that survives these checks belongs to the same **validation
class** as every other backend that does — nothing more, nothing less.

## Backends

| Backend | Status | Source of truth |
|---|---|---|
| `sagnac_backend` | validated (reference layer 0) | [Sagnac-Reference-Transport](https://github.com/error-wtf/Sagnac-Reference-Transport) — read-only adapter, no formulas here |
| `ssz_synthetic_backend` | development control (toy frame-dragging proxy) | intentionally kept as negative/development control |
| `ssz_closure_backend` | **REAL SSZ transport** — read-only adapter | [SSZ_FULL_CLOSURE](https://github.com/error-wtf/SSZ_FULL_CLOSURE) `postclosure/transport.py`: k^nu nabla_nu k^mu = 0, null congruence, eikonal phase S_r, redshift — zero physics in the bridge |

Direction-contract rule: static spherical SSZ declares **no** sign-swap
involution (the Sagnac `v -> -v` symmetry is NOT universal physics).
The inverse problem for SSZ is explicitly deferred (`inverse=""`).

## Shared vocabulary (contracts.py)

```
SignalLaw               immutable parameters + provenance hash
TransportArchitecture   declares structure + direction parameter
TransportOperator       bound architecture + law + forward solve
PropagatedState         directional results, tagged ±1
ObservableSet           directional / scalar / inverse names
```

Check functions: `check_forward`, `check_parity`, `check_convergence`,
`check_inverse_round_trip`.

## What the tests prove

`tests/test_bridge_contracts.py` runs the SAME structural checks against
BOTH backends with the SAME code:

- direction swap flips directional observables, preserves scalars
- inverse round trip reproduces the frozen parameter
- unknown-parameter perturbation fails closed

Sagnac recovers `v` from `(t_+, t_-)`.  SSZ recovers the frame sign
from synthetic observables generated from a **known frozen geometry** —
the synthetic-inversion test that must pass before any real astronomy
is touched.

## Status

Both backends green (8/8 contract tests).  The reference repo is
validated to `SAGNAC_REFERENCE_CLOSURE_PASS` (SAG-S1..S10, 10/10
negative controls); this bridge inherits that validation through the
read-only adapter and adds the class-membership proof.

## Position in the overall chain

```
Sagnac analytic truth
  ↓
reference transport validated          (this chain's layer 0)
  ↓
generic bridge extracted               ← THIS REPO
  ↓
SSZ synthetic forward/inverse test     (ssz_backend, done first)
  ↓
healthy SSZ spectral operator
  ↓
certified modes → real spectroscopy
```
