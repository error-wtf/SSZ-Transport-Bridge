"""SSZ closure backend — INDEPENDENT RE-SOLVER (ADR-001, Variante B).

DECISION (2026-10-06, ADR-001): this backend is an independent JUDGE that
RE-DERIVES the SSZ transport equations from frozen metric data and solves
them numerically with its own integrators.  It does NOT import physics
functions from SSZ_FULL_CLOSURE — only frozen member data (MODEL_LOCK ->
member CSV) via hash-bound provenance.

Semantics: `independent reimplementation from frozen metric data` — NOT
`zero shared physics formulas`.  The earlier adapter wording was wrong for
this backend and has been retired (docs/REAL_SSZ_BRIDGE_SPEC.md ADR-001).

Every independently re-derived physics routine MUST be declared in
INDEPENDENT_PHYSICS below.  The anti-circularity audit
(tools/audit_backend_independence.py) enforces: (a) no closure-code
imports, (b) no undeclared physics-like numerics, (c) registry entries
present for every detected pattern.  This keeps the independence claim
auditable instead of aspirational.

Symmetry: static spherical geometry.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from transport_bridge import (
    ObservableSet,
    SignalLaw,
    TransportArchitecture,
    TransportOperator,
)
from transport_bridge.provenance import load_closure_provenance

# --- Declared independent physics (audited registry) ---------------------
# Every re-derived physics routine below MUST have an entry here.  The
# anti-circularity audit cross-checks this registry against the AST.
INDEPENDENT_PHYSICS = {
    "null_geodesic_rhs": (
        "Null geodesic RHS for static spherical metrics, independently "
        "re-derived from ds^2 = -f dt^2 + dr^2/h + r^2 dphi^2: "
        "dk^t/dlam = -(f'/f) k^t k^r; dk^r/dlam = -(h f'/2)(k^t)^2 "
        "+ (h'/(2h))(k^r)^2 + h r (k^phi)^2; dk^phi/dlam = -(2/r) k^r k^phi. "
        "Implemented in _independent_null_transport.rhs."
    ),
    "null_conserved_quantities": (
        "E = f k^t, L = r^2 k^phi, norm = -f (k^t)^2 + (k^r)^2/h "
        "+ r^2 (k^phi)^2; checked as conservation residuals along the "
        "independent solution."
    ),
    "phase_integral": (
        "Eikonal phase / Shapiro delay quadrature "
        "integral dr/sqrt(f h) (trapezoid on the spline, 20000 points)."
    ),
    "redshift_ratio": (
        "Static redshift sqrt(f_a / f_b) between endpoints."
    ),
}

# --- Independent Numerical Solvers (The Judge's Tools) ---

def _independent_null_transport(f_spline, h_spline, b, r0, inward=True, lam_span=(0.0, 30.0)):
    """Independent check of null geodesic conservation."""
    f0 = float(f_spline(r0))
    h0 = float(h_spline(r0))
    
    val = h0 * (1.0 / f0 - b**2 / r0**2)
    if val < 0:
        raise ValueError("Forbidden b at r0")
    
    kr0 = -np.sqrt(val) if inward else np.sqrt(val)
    kt0 = 1.0 / f0
    kph0 = b / r0**2

    def rhs(lam, y):
        _, r, _, kt, kr, kph = y
        f = float(f_spline(r))
        fp = float(f_spline.derivative()(r))
        h = float(h_spline(r))
        hp = float(h_spline.derivative()(r))
        return [
            kt, kr, kph,
            -(fp / f) * kt * kr,
            -(0.5 * h * fp) * kt**2 + (hp / (2.0 * h)) * kr**2 + h * r * kph**2,
            -(2.0 / r) * kr * kph,
        ]

    sol = solve_ivp(rhs, lam_span, [0.0, r0, 0.0, kt0, kr0, kph0], 
                    method="DOP853", rtol=1e-12, atol=1e-12)
    
    if not sol.success:
        raise RuntimeError("Judge solve failed")

    _t, r, _phi, kt, kr, kph = sol.y
    E_arr = f_spline(r) * kt
    L_arr = r**2 * kph
    nrm = -f_spline(r) * kt**2 + kr**2 / h_spline(r) + r**2 * kph**2
    
    return {
        "max_dE": float(np.max(np.abs(E_arr - 1.0))),
        "max_dL": float(np.max(np.abs(L_arr - b))),
        "max_norm_res": float(np.max(np.abs(nrm))),
    }

def _independent_phase_redshift(f_spline, h_spline, r_a, r_b):
    """Independent quadrature for eikonal phase and redshift."""
    r_grid = np.linspace(min(r_a, r_b), max(r_a, r_b), 20000)
    integ = float(np.trapezoid(1.0 / np.sqrt(f_spline(r_grid) * h_spline(r_grid)), r_grid))
    z = float(np.sqrt(f_spline(r_a) / f_spline(r_b)))
    
    return {
        "shapiro_dt_per_E": integ,
        "reduced_phase_per_E": integ,
        "redshift_a_to_b": z,
    }

# --- Backend Adapter ---

ARCHITECTURE = TransportArchitecture(
    name="ssz_independent_judge",
    direction_parameter="inward",
    operator_class="null_geodesic",
)

OBSERVABLES = ObservableSet(
    directional=(), 
    scalar=("shapiro_dt_per_E", "reduced_phase_per_E", 
            "redshift_a_to_b", "max_dE", "max_dL", "max_norm_res"),
    inverse="", 
)

def closure_repo_path() -> Path:
    """Fail-closed resolution of the closure checkout.

    Priority: SSZ_FULL_CLOSURE_PATH env var -> default local layout.  Raises
    with an actionable message when the checkout is unavailable — the judge
    never guesses.
    """
    import os
    p = os.environ.get("SSZ_FULL_CLOSURE_PATH",
                       "/home/error/physics/clones/SSZ_FULL_CLOSURE")
    path = Path(p)
    if not (path / "MODEL_LOCK.json").exists():
        raise RuntimeError(
            f"SSZ_FULL_CLOSURE checkout not usable at {path} — set "
            "SSZ_FULL_CLOSURE_PATH to a checkout containing MODEL_LOCK.json")
    return path


def make_operator(law: SignalLaw, r_a: float = 1.5, r_b: float = 1.6) -> TransportOperator:
    repo_path = closure_repo_path()
    
    lock = json.loads((repo_path / "MODEL_LOCK.json").read_text())
    member_p = repo_path / lock["action_member_stream"]
    
    from .utils import load_metric_splines
    f_spl, h_spl = load_metric_splines(member_p)

    def _solve(law_: SignalLaw) -> dict[str, Any]:
        inward = bool(law_.parameters.get("inward", 1.0) >= 0)
        b = float(law_.parameters.get("b", 2.5))
        r0 = float(law_.parameters.get("r0", r_a))
        
        transport = _independent_null_transport(f_spl, h_spl, b, r0, inward)
        phase = _independent_phase_redshift(f_spl, h_spl, r_a, r_b)
        
        return {**transport, **phase}

    return TransportOperator(ARCHITECTURE, law, _solve)

def provenance() -> dict[str, Any]:
    repo_path = closure_repo_path()
    p = load_closure_provenance(repo_path)
    # Expose the full provenance binding, with the historical key names the
    # contract tests assert on:
    return {
        **p,
        "closure_git_head": p["commit_sha"],
        "closure_member_sha256": p["member_sha256"],
        "closure_src": str(repo_path),
    }
