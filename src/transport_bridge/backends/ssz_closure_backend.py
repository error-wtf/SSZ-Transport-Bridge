"""SSZ closure backend — INDEPENDENT VALIDATOR.

CORE RULE: This backend is a JUDGE. It does NOT import physics functions from
SSZ_FULL_CLOSURE. It consumes frozen metric data and performs independent
numerical solves to certify the transport results.

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

def make_operator(law: SignalLaw, r_a: float = 1.5, r_b: float = 1.6) -> TransportOperator:
    repo_path = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
    
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
    repo_path = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
    p = load_closure_provenance(repo_path)
    # Expose the full provenance binding, with the historical key names the
    # contract tests assert on:
    return {
        **p,
        "closure_git_head": p["commit_sha"],
        "closure_member_sha256": p["member_sha256"],
        "closure_src": str(repo_path),
    }
