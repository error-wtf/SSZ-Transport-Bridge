#!/usr/bin/env python3
"""Write SUPERSEDED_V4_2CHANNEL_VERDICTS_V1.json into SSZ-Transport-Bridge:
marks the old 2-channel constraint-artifact verdicts as superseded by the
Closure DOF re-audit (3 physical DOF, dphi dynamical), without deleting
the old evidence.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

BRIDGE = Path("/home/error/SSZ-Transport-Bridge")
CLOSURE = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main() -> int:
    sup = {
        "audit": "SUPERSEDED_V4_2CHANNEL_VERDICTS_V1",
        "date": "2026-10-07",
        "supersedes": [
            {"artifact": "B8_SECTOR_TRACKING_V4_V1.json",
             "old_claim": ("negative full-3x3 kinetic direction = 100% chi "
                            "(dphi) constraint artifact; 2-channel (psi,V) "
                            "reduction valid"),
             "why_superseded": ("Closure DOF re-audit (V4_DOF_ARCHITECTURE_"
                                 "REAUDIT_V2 + V4_UNREDUCED_CONSTRAINT_RANK_V1) "
                                 "proved from the unreduced 8x8 Euler "
                                 "descriptor that dphi carries kinetic terms "
                                 "(e1*(dphi_dot)^2, max coeff 3.8e+06) and "
                                 "NO primary constraint removes it. Sec-20.4 "
                                 "physical basis is y=(psi,dphi,V) — the "
                                 "auxiliary eliminations are h1/H2/H1/dA0/"
                                 "dA1. The 2x2 (psi,V) strike was a "
                                 "projection WITHOUT constraint justification. "
                                 "The negative dphi direction is a PHYSICAL "
                                 "kinetic pathology of the +a1 branches "
                                 "(H3-ghost)."),
             "status": "SUPERSEDED — retained as evidence"},
        ],
        "current_truth": {
            "physical_basis": ["psi", "dphi", "V"],
            "n_dof": 3,
            "evidence": [
                "V4_UNREDUCED_CONSTRAINT_RANK_V1.json",
                "V4_FULL_3DOF_BRANCH_HEALTH_V1.json",
                "V4_FULL_3DOF_OPERATOR_V2.npz/.json",
                "V4_FULL_3DOF_KINETIC_HEALTH_V1.json",
            ],
            "production_branch": {"a1": -0.5, "eps": -0.3},
            "healthy_branches": [["-0.5", "+0.3"], ["-0.5", "-0.3"]],
            "ghost_branches": [["+0.5", "+0.3"], ["+0.5", "-0.3"]],
            "closure_head": "ce1c1569b5c8ec3f28f0a69834de297e27a516f8",
        },
        "bridge_implication": (
            "The independent validator (this repository) must consume the "
            "FULL 3-field operator (psi,dphi,V) from closure commit "
            "ce1c156+, NOT the superseded 2-channel reduction. The old "
            "B8 tracking numbers remain diagnostic material only."),
        "notation": ("'chi' in older artifacts names the three-component "
                      "physical vector (psi,dphi,V); it is NOT dphi alone."),
        "closure_reaudit_sha256": sha(
            CLOSURE / "data/generated/spectral/V4_DOF_ARCHITECTURE_REAUDIT_V2.json"),
        "unreduced_rank_sha256": sha(
            CLOSURE / "data/generated/spectral/V4_UNREDUCED_CONSTRAINT_RANK_V1.json"),
    }
    out = BRIDGE / "artifacts/SUPERSEDED_V4_2CHANNEL_VERDICTS_V1.json"
    out.write_text(json.dumps(sup, indent=1) + "\n")
    print("written:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
