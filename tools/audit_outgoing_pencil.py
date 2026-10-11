#!/usr/bin/env python3
"""B6 Step 1: audit of the asymptotic outgoing matrix pencil.

Mandate package 2: derive the three outgoing eigenchannels from the FULL
frequency-dependent asymptotic pencil

    P_inf(w, k) = w^2 K - k^2 G + i k (G' - S) - (M + S'/2),

verify (a) all six eigenpairs are kept, (b) k_j genuinely depend on
omega, (c) the acoustic branch selection min |k - w| is stable under
perturbation, (d) residual asymptotic 1/r couplings at the matching
radius are small and quantified, (e) the frozen per-channel-k BC flag
(common_speed_bc_invalid) is present — the old common e^{i w r*} ansatz
is NOT adopted.

Operates on the frozen V3_3DOF export only (hash-checked).
"""
import json
import sys
from pathlib import Path

import numpy as np

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

from transport_bridge.independent_qnm_3dof_v1 import FrozenOperator, IndependentQNM3DOF  # noqa: E402

OP_NPZ = BRIDGE / "data/frozen_specs/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"


def main():
    op = FrozenOperator(str(OP_NPZ))
    slv = IndependentQNM3DOF(op)
    r_eval = op.r[-1]  # outermost exported node: asymptotic region
    G, A, B, K = op.coeffs(r_eval)
    report = {
        "artifact": "SSZ_B6_STEP1_PENCIL_AUDIT_V1",
        "npz_sha256": op.EXPECTED_SHA256,
        "r_eval": r_eval,
        "checks": {},
    }

    # (a)+(b): all 6 roots per omega, acoustic selection, k_j(w) dispersion
    ws = [complex(0.5, -0.1), complex(1.0, -0.5), complex(2.0, -1.5),
          complex(2.5, -0.2), complex(3.0, -0.8)]
    disp = []
    for w in ws:
        roots = slv.pencil_roots(w)
        acoustic = slv.outgoing_ks(w)
        all_sorted = sorted(roots, key=lambda k: (k.real, k.imag))
        report["checks"][f"w={w}"] = {
            "all_6_roots": [f"{complex(k):.6f}" for k in all_sorted],
            "acoustic_3": [f"{complex(k):.6f}" for k in acoustic],
            "max_|k-w|_acoustic": float(max(abs(k - w) for k in acoustic)),
        }
        disp.append({"w": f"{w:.2f}", "max_|k_ac-w|":
                     float(max(abs(k - w) for k in acoustic))})

    # (c): selection stability under relative perturbation of omega
    stab = {}
    for w in [complex(1.0, -0.5), complex(2.5, -0.2)]:
        a1 = np.sort_complex(slv.outgoing_ks(w))
        a2 = np.sort_complex(slv.outgoing_ks(w * (1 + 1e-9)))
        stab[f"w={w}"] = float(np.max(np.abs(a1 - a2)))
    report["checks"]["selection_stability_drift"] = stab

    # (b) explicit: acoustic branches are NOT frequency-independent.
    # Compare acoustic k at two very different omega.
    ka_lo = sorted(slv.outgoing_ks(complex(0.5, -0.1)), key=lambda k: k.real)
    ka_hi = sorted(slv.outgoing_ks(complex(3.0, -0.8)), key=lambda k: k.real)
    report["checks"]["k_dispersion_nontrivial"] = {
        "acoustic(w=0.5-0.1j)": [f"{complex(k):.4f}" for k in ka_lo],
        "acoustic(w=3.0-0.8j)": [f"{complex(k):.4f}" for k in ka_hi],
        "max_pairwise_diff": float(np.max(np.abs(np.array(ka_lo) - np.array(ka_hi)))),
    }

    # (d): residual subleading terms at r_match vs 2*r_match (asymptotic
    # decay check: differences bound the residual 1/r coupling magnitude)
    G1, A1, B1, K1 = op.coeffs(slv.r_match)
    G2, A2, B2, K2 = op.coeffs(2 * slv.r_match)
    report["checks"]["asymptotic_residuals_rmatch_vs_2rmatch"] = {
        "rel_|G(r)-G(2r)|": float(np.max(np.abs(G1 - G2)) / max(np.max(np.abs(G2)), 1e-300)),
        "rel_|K(r)-K(2r)|": float(np.max(np.abs(K1 - K2)) / max(np.max(np.abs(K2)), 1e-300)),
        "abs_max|A(r)|=|G'-S|": float(np.max(np.abs(A1))),
        "abs_max|B(r)|=|M+S'/2|": float(np.max(np.abs(B1))),
    }

    # (e): BC consistency — per-channel-k flag must exist and be True;
    # the old common-speed ansatz is used only as the flagged-invalid
    # reference, never as the Step-1 boundary condition.
    bc_ok = False
    source = None
    for spec in sorted((BRIDGE / "data/frozen_specs").glob("*.json")):
        try:
            d = json.loads(spec.read_text())
        except Exception:
            continue
        stack = [d]
        while stack:
            x = stack.pop()
            if isinstance(x, dict):
                if x.get("common_speed_bc_invalid") is True:
                    bc_ok, source = True, spec.name
                stack.extend(x.values())
            elif isinstance(x, list):
                stack.extend(x)
    report["checks"]["per_channel_bc_flag_present"] = bc_ok
    report["checks"]["bc_flag_source"] = source

    n_ok = sum(1 for f in report["checks"] if f.startswith("w="))
    print(json.dumps(report, indent=1))
    out = BRIDGE / "artifacts" / "SSZ_B6_STEP1_PENCIL_AUDIT_V1.json"
    out.write_text(json.dumps(report, indent=1))
    print(f"\nwritten: {out}  ({n_ok} frequency samples audited)")


if __name__ == "__main__":
    main()
