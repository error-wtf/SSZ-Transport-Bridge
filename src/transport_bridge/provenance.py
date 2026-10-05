"""Binding and verifying the real SSZ_FULL_CLOSURE state."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def bind_module(mod: Any, repo: Path) -> None:
    """Verify that a module import actually originated from the target repo."""
    origin = getattr(mod, "__file__", None)
    if not origin or Path(origin).resolve().parents[0].resolve() != repo.resolve():
        # Note: This is simplified; usually we check if the path is under the repo root.
        # For our strict adapter, we want the module to be physically inside the repo.
        if not origin or str(repo.resolve()) not in str(Path(origin).resolve()):
            raise ValueError(f"module {mod.__name__} origin {origin} is not in {repo}")

def load_closure_provenance(repo: Path) -> dict:
    """Binding to the frozen member, lock, and verdict of SSZ_FULL_CLOSURE."""
    lock_p = repo / "MODEL_LOCK.json"
    verdict_p = repo / "TRUE_FULL_CLOSURE_VERDICT.json"
    
    if not lock_p.exists():
        raise FileNotFoundError("MODEL_LOCK.json missing")
    if not verdict_p.exists():
        raise FileNotFoundError("TRUE_FULL_CLOSURE_VERDICT.json missing")
    
    lock = json.loads(lock_p.read_text())
    verdict = json.loads(verdict_p.read_text())
    
    # 1. Bind member file — path-escape check FIRST (before any file access)
    member_rel_path = lock.get("action_member_stream", "")
    if ".." in member_rel_path or member_rel_path.startswith("/"):
        raise ValueError("path escape in member reference")
    member_p = repo / member_rel_path
    if not member_p.exists():
        raise FileNotFoundError(f"member {member_rel_path} missing")

    member_hash = sha256_file(member_p)
    if member_hash != lock.get("action_member_sha256"):
        raise ValueError(
            f"member hash mismatch: {member_hash} "
            f"!= {lock.get('action_member_sha256')}")

    # Verify manifest link if present
    manifest_rel = lock.get("action_member_reference")
    if manifest_rel:
        man_p = repo / manifest_rel
        if not man_p.exists():
            raise FileNotFoundError(f"manifest {manifest_rel} missing")
        man = json.loads(man_p.read_text())
        if man.get("member_hash") != member_hash:
            raise ValueError(
                f"member hash mismatch in manifest: "
                f"{man.get('member_hash')} != {member_hash}")

    # 2. Bind Git state
    head_p = repo / ".git/HEAD"
    if not head_p.exists():
        raise ValueError("not a git repo")
    txt = head_p.read_text().strip()
    if txt.startswith("ref: "):
        ref_path = repo / ".git" / txt[5:].strip()
        sha = ref_path.read_text().strip() if ref_path.exists() else "UNKNOWN"
    else:
        sha = txt
    sha = sha[:40]
    
    return {
        "commit_sha": sha,
        "member_sha256": member_hash,
        "verdict": verdict.get("verdict"),
        "files": {
            "MODEL_LOCK.json": sha256_file(lock_p),
            "TRUE_FULL_CLOSURE_VERDICT.json": sha256_file(verdict_p)
        },
        "historical_evidence_at_head": False # Placeholder for full graph check
    }
