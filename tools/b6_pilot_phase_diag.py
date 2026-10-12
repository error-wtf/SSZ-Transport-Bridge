#!/usr/bin/env python3
"""Pilot-pole phase-continuity diagnostic for the argument-principle probe.

Question: the certified pilot pole sits at |W| ~ 1e-10..1e-7 (5+ orders
above the det roundoff floor ~1.3e-17), yet v2 counted N=0.0 on the
h=0.008 contour.  If the true zero of the underlying analytic W is
inside, the phase must wind once.  N=0 means the phase of the NUMERICAL
W (det of non-holomorphically normalized bases) does NOT wind -- the
column normalization / QR sign-fix introduces path-dependent phase
jumps that destroy the analytic winding.

Output: per-point arg(W) and |W| along the contour + total winding,
written to stdout (log).
"""
import itertools
import sys
from pathlib import Path

import numpy as np

BRIDGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRIDGE / "src"))

from transport_bridge.independent_qnm_3dof_v1 import FrozenOperator  # noqa: E402
from transport_bridge.stable_evans_v3 import StableEvansV3  # noqa: E402

op = FrozenOperator(str(BRIDGE / "data" / "frozen_specs"
                        / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"))
ev = StableEvansV3(op, r_match=30.0, r_seed=55.0, n_segments=40, rtol=1e-9)


def W(om):
    qi = ev.interior_basis(om)
    qe = ev.exterior_incoming_basis(om)
    return complex(np.linalg.det(np.hstack([qi, qe])))


c, h = 0.22642588643752676 + 0.08455094304111954j, 0.008
corners = [c + complex(-h, -h), c + complex(h, -h), c + complex(h, h),
           c + complex(-h, h), c + complex(-h, -h)]
vals = []
for a, b in itertools.pairwise(corners):
    for t in np.linspace(0, 1, 8, endpoint=False):
        vals.append(W(a + (b - a) * t))
vals = np.array(vals)

print("Pilot contour 8/side: |W| in [%.2e, %.2e]"
      % (np.abs(vals).min(), np.abs(vals).max()), flush=True)
ph = np.unwrap(np.angle(vals))
total = ph[-1] - ph[0] + np.angle(vals[0] * np.conj(vals[-1]))
print("total phase: %.3f rad => N=%.3f" % (total, total / (2 * np.pi)),
      flush=True)
for i, v in enumerate(vals):
    print("  %2d: arg=%+.4f  |W|=%.3e" % (i, np.angle(v), abs(v)), flush=True)
