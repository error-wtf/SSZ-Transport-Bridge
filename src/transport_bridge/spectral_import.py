"""Frozen spectral export consumer — the INDEPENDENT spectral judge.

Consumes FROZEN_SPECTRAL_OPERATOR_EXPORT_V1 (npz + provenance JSON) from the
SSZ_FULL_CLOSURE checkout.  No SSZ python module is imported here: the only
inputs are the hash-verified coefficient matrices and this file's own generic
numerics (numpy/scipy).  Producer != Validator.

Implements the data side of SSZ-B1 (provenance) and SSZ-B2 (operator
reconstruction).  SSZ-B3..B10 (independent solve/certification) build on
this loader.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

CLOSURE_REPO = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
EXPORT_NPZ = CLOSURE_REPO / (
    "data/generated/spectral/FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.npz")
EXPORT_JSON = CLOSURE_REPO / (
    "data/generated/spectral/FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.json")


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_export(repo: Path = CLOSURE_REPO) -> dict[str, Any]:
    """Fail-closed load: hash binding + health window declared."""
    npz_p = repo / EXPORT_NPZ.relative_to(CLOSURE_REPO)
    json_p = repo / EXPORT_JSON.relative_to(CLOSURE_REPO)
    if not npz_p.exists():
        raise FileNotFoundError(
            f"frozen spectral export npz missing: {npz_p}")
    if not json_p.exists():
        raise FileNotFoundError(
            f"frozen spectral export provenance missing: {json_p}")
    prov = json.loads(json_p.read_text())
    digest = sha256_file(npz_p)
    if digest != prov.get("npz_sha256"):
        raise ValueError(
            f"export hash mismatch: {digest} != {prov.get('npz_sha256')}")
    if prov.get("export") != "FROZEN_SPECTRAL_OPERATOR_EXPORT_V1":
        raise ValueError("unknown export type")
    return {
        "prov": prov,
        "arrays": dict(np.load(npz_p, allow_pickle=False)),
        "npz_sha256": digest,
    }


def load_operator(export: dict[str, Any], L: int) -> dict[str, np.ndarray]:
    """Return the reduced operator matrices for multipole L."""
    key = f"{int(L)}_"
    arrays = export["arrays"]
    if key + "K" not in arrays:
        raise KeyError(f"L={L} not in frozen export")
    out = {}
    for name in ("r", "K", "R", "G", "S", "M", "Dh1", "DeltaV", "pivotA0"):
        out[name] = arrays[key + name]
    return out


def independent_local_spectra(
    mats: dict[str, np.ndarray], n_check: int = 24
) -> dict[str, Any]:
    """SSZ-B2 first numerical step, fully independent: symmetrize K and G,
    solve the local generalized eigenproblem K psi = lambda G psi at a fixed
    subset of radial nodes, and report the eigenvalue traces.

    This is NOT a QNM claim — no boundary conditions, no complex poles.  It
    is the operator-reconstruction cross-check the bridge needs before any
    spectral certification layer is attempted.
    """
    from scipy.linalg import eigh

    r = mats["r"]
    n = len(r)
    idx = np.unique(np.linspace(0, n - 1, min(n_check, n)).round().astype(int))
    traces: list[np.ndarray] = []
    for i in idx:
        K = mats["K"][i]
        G = mats["G"][i]
        Ks = (K + K.T) / 2.0
        Gs = (G + G.T) / 2.0
        # symmetric-definite generalized problem K psi = lambda G psi —
        # the pinv route is numerically wrong for indefinite G:
        w = eigh(Ks, Gs, eigvals_only=True)
        traces.append(np.sort(w))
    return {
        "r": r[idx],
        "eigen_traces": np.asarray(traces),
    }
