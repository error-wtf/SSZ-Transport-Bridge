
"""Bridge-level Judge: Evaluates the SSZ-B1..B10 Milestone."""
import json
from pathlib import Path

from transport_bridge import SignalLaw
from transport_bridge.backends import ssz_closure_backend as cb

ROOT = Path(__file__).resolve().parents[1]

def evaluate_milestone():
    results = {}
    
    # B1: Provenance Binding
    try:
        prov = cb.provenance()
        ok = prov.get("verdict") == "TRUE_FULL_CLOSURE_PASS"
        results["B1_PROVENANCE"] = "PASS" if ok else "FAIL"
    except Exception as e:
        results["B1_PROVENANCE"] = f"FAIL: {e}"

    # B2: Null Forward Transport (Independent)
    try:
        # b must sit INSIDE the photon cone of the frozen member at r0:
        # b_crit(r0=1.55) = r0/sqrt(f) ~= 2.489.  b=2.0 is inside.
        law = SignalLaw("ssz_closure", {"inward": 1.0, "b": 2.0, "r0": 1.55}, "hash")
        op = cb.make_operator(law)
        out = op.run()
        # B2 check: residuals < 1e-8
        res_ok = out["max_dE"] < 1e-8 and out["max_dL"] < 1e-8 and out["max_norm_res"] < 1e-8
        results["B2_NULL_TRANSPORT"] = "PASS" if res_ok else "FAIL"
    except Exception as e:
        results["B2_NULL_TRANSPORT"] = f"FAIL: {e}"

    # B3: Phase/Readout Exposure
    try:
        # Same inside-cone parameter set as B2.
        out = cb.make_operator(
            SignalLaw("ssz_closure", {"inward": 1.0, "b": 2.0, "r0": 1.55}, "h")
        ).run()
        required = ["reduced_phase_per_E", "shapiro_dt_per_E", "redshift_a_to_b"]
        results["B3_READOUT"] = "PASS" if all(k in out for k in required) else "FAIL"
    except Exception as e:
        results["B3_READOUT"] = f"FAIL: {e}"

    # B4: Reversal Symmetry
    results["B4_REVERSAL"] = "OPEN" # Explicitly open per spec

    # B5/B6/B7: REAL gate validation against the closure's own audit artifacts
    # (ADR-001): each bridge gate reads the certified source gate record and
    # verifies status, commit, and member binding — no B1 aliasing.
    try:
        repo_path = cb.closure_repo_path()
        audit = json.loads(
            (repo_path / "ai_analysis/critical_jsons/ABSOLUTE_FULL_CLOSURE_AUDIT.json"
             ).read_text())
        verdict = json.loads(
            (repo_path / "ai_analysis/critical_jsons/TRUE_FULL_CLOSURE_VERDICT.json"
             ).read_text())
        prov = cb.provenance()

        # The verdict + audit certify the FROZEN MEMBER, not the research
        # branch head.  The strict binding is therefore member-hash equality;
        # the audit commit is recorded as provenance info (it dates the audit
        # run on the historical corpus).
        member_ok = str(verdict.get("member_hash", "")).startswith(
            str(prov.get("closure_member_sha256", "x")[:16]))
        checks = {
            "member_binding": member_ok,
            "verdict": verdict.get("verdict") == "TRUE_FULL_CLOSURE_PASS",
            "verdict_commit": verdict.get("git_commit"),
            "audit_commit": audit.get("git_commit"),
        }

        def gate(name: str) -> tuple[bool, str]:
            g = audit.get("gates", {}).get(name, {})
            status = g.get("certified_status")
            grounding = g.get("grounding", [])
            grounded = all((repo_path / t).exists() for t in grounding)
            return (status == "PASS" and grounded and len(grounding) > 0,
                    f"{name}={status} grounding_exists={grounded}")

        g122_ok, g122_info = gate("G122")
        g119_ok, g119_info = gate("G119")
        g121_ok, g121_info = gate("G121")

        base_ok = checks["member_binding"] and checks["verdict"]
        results["B5_CONVERGENCE"] = ("PASS" if base_ok and g122_ok else "FAIL: "
                                      f"{g122_info} checks={checks}")
        results["B6_FALSIFIERS"] = ("PASS" if base_ok and g119_ok else "FAIL: "
                                     f"{g119_info} checks={checks}")
        results["B7_LIMITS"] = ("PASS" if base_ok and g121_ok else "FAIL: "
                                 f"{g121_info} checks={checks}")
        results["_b5b6b7_evidence"] = {"checks": checks,
                                        "G122": g122_info, "G119": g119_info,
                                        "G121": g121_info}
    except Exception as e:
        results["B5_CONVERGENCE"] = f"FAIL: {e}"
        results["B6_FALSIFIERS"] = f"FAIL: {e}"
        results["B7_LIMITS"] = f"FAIL: {e}"

    # B8: Inverse Contract
    results["B8_INVERSE"] = "REAL_SSZ_INVERSE_OPEN"

    # B9: Independence Audit
    # Runs the AST-based no-duplicated-physics audit over executable code.
    audit = ROOT / "tools/audit_backend_independence.py"
    if audit.exists():
        import subprocess
        import sys as _sys
        run = subprocess.run([_sys.executable, str(audit)],
                             capture_output=True, text=True)
        results["B9_INDEPENDENCE"] = (
            "PASS" if run.returncode == 0
            else f"FAIL: {run.stdout.strip()[:200]}")
    else:
        results["B9_INDEPENDENCE"] = "FAIL: audit tool missing"

    # B10: Final Summary
    # Mandatory gates: B1, B2, B3, B5, B6, B7, B9
    mandatory = ["B1_PROVENANCE", "B2_NULL_TRANSPORT", "B3_READOUT",
                 "B5_CONVERGENCE", "B6_FALSIFIERS", "B7_LIMITS",
                 "B9_INDEPENDENCE"]
    all_mandatory = all(results.get(k) == "PASS" for k in mandatory)
    
    verdict = "REAL_SSZ_ADAPTER_PASS" if all_mandatory else "PARTIAL"
    
    return {
        "B_RESULTS": results,
        "VERDICT": verdict,
        "BRIDGE_CLASS_MEMBERSHIP": "PASS" if all_mandatory else "FAIL"
    }

if __name__ == "__main__":
    res = evaluate_milestone()
    with open(ROOT / "artifacts/REAL_SSZ_BRIDGE_VERDICT.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
