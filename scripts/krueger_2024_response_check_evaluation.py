#!/usr/bin/env python3
"""Evaluate the one R1.9 Krueger response check against its pre-run manifest.

This scorer does not fit parameters, open held-out observations, or select another candidate.  It
binds the actual audit to the committed target, two response predictions, model-error envelope,
and trust-ratio rule.  R1.9 closes the 10 nm candidate sequence for every returned classification.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


M0 = 2.616705
M_PRED = 0.447528


def _sha(path: Path | str) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _canonical_sha(payload: dict, field: str) -> str:
    content = dict(payload)
    content.pop(field, None)
    return sha256(json.dumps(
        content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _maximum(history: list[dict], key: str) -> float:
    values = [float(item[key]) for item in history if key in item]
    return max(values, default=0.0)


def evaluate(launch_path: Path | str, audit_path: Path | str) -> dict:
    launch_path = Path(launch_path)
    audit_path = Path(audit_path)
    launch = json.loads(launch_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))

    if (launch.get("schema") != "petch.krueger-2024.r19-response-check-launch.v1"
            or launch.get("protocol_id") != "K24-PETCH-R1.9"
            or launch.get("authority") is not False
            or launch.get("held_out_profile_data_read") is not False
            or launch.get("ten_nm_sequence_closes_after_this_evaluation") is not True):
        raise ValueError("launch manifest is not the sealed R1.9 response check")

    config = audit.get("configuration", {})
    candidate = launch["candidate"]
    required = {
        "boundary_case": "base",
        "duration_s": 60.0,
        "dx_um": 0.01,
        "ion_azimuthal_closure": "axisymmetric_uniform",
        "ion_azimuthal_order": 16,
        "topology_change_policy": "continue_gas_cavity",
        "effective_mask_crosslinked_growth_fraction": candidate[
            "effective_mask_crosslinked_growth_fraction"],
        "oxide_etch_yield_scale": candidate["oxide_etch_yield_scale"],
    }
    for key, expected in required.items():
        actual = config.get(key)
        if isinstance(expected, float):
            if not np.isclose(float(actual), expected, rtol=0.0, atol=2e-14):
                raise ValueError(f"audit configuration changed sealed field {key}")
        elif actual != expected:
            raise ValueError(f"audit configuration changed sealed field {key}")

    history = list(audit.get("history", ()))
    if not history:
        raise ValueError("response audit has no trajectory history")
    final = audit.get("final_metrics", history[-1].get("metrics", {}))
    actual = np.asarray([
        float(final["mask_opening_nm"]),
        float(final["etch_depth_nm"]),
    ])
    target = np.asarray([
        float(launch["target_nm"]["mask_opening"]),
        float(launch["target_nm"]["etch_depth"]),
    ])
    local_prediction = np.asarray([
        float(launch["committed_predictions_nm"]["same_operator_local_response"][
            "mask_opening"]),
        float(launch["committed_predictions_nm"]["same_operator_local_response"][
            "etch_depth"]),
    ])
    model_gate = np.asarray([
        float(launch["committed_model_error_gate_nm"]["mask_opening"]),
        float(launch["committed_model_error_gate_nm"]["etch_depth"]),
    ])
    scaled_merit = float(np.linalg.norm((actual - target) / 5.0))
    rho = float((M0 - scaled_merit) / (M0 - M_PRED))

    reached_time = float(history[-1].get("physical_time_s", np.nan))
    status_complete = (
        audit.get("status") == "complete"
        and np.isclose(reached_time, 60.0, rtol=0.0, atol=1e-12)
        and audit.get("terminal_event") is None)
    ledger_maximum = _maximum(
        history, "maximum_material_ledger_residual_units_m2")
    radiosity_maximum = _maximum(
        history, "maximum_neutral_radiosity_relative_balance_error")
    validity_pass = all(
        item.get("validity", {}).get("within_declared_scope", True)
        for item in history)
    numerical_pass = bool(status_complete and ledger_maximum == 0.0 and validity_pass)
    target_cell_pass = bool(np.all(np.abs(actual - target) <= 5.0))
    model_envelope_pass = bool(np.all(np.abs(actual - local_prediction) <= model_gate))
    strong_pass = bool(numerical_pass and target_cell_pass and model_envelope_pass)

    if strong_pass:
        decision = "strong_response_model_pass"
    elif not numerical_pass:
        decision = "reject_numerical_contract"
    elif rho >= 0.25 and scaled_merit < M0:
        decision = "accept_hold_radius"
    elif rho >= 0.10 and scaled_merit < M0:
        decision = "marginal_shrink_and_stop"
    else:
        decision = "reject_response_model"

    payload = {
        "schema": "petch.krueger-2024.r19-response-check-evaluation.v1",
        "protocol_id": "K24-PETCH-R1.9",
        "scientific_status": (
            "10 nm response-model evidence only; no fine-grid authority or held-out validation"),
        "authority": False,
        "held_out_profile_data_read": False,
        "ten_nm_sequence_closed": True,
        "inputs": {
            "launch_manifest": {
                "path_name": launch_path.name,
                "sha256": _sha(launch_path),
            },
            "run_audit": {
                "path_name": audit_path.name,
                "sha256": _sha(audit_path),
                "config_hash": audit.get("config_hash"),
            },
        },
        "actual_nm": {
            "mask_opening": float(actual[0]),
            "etch_depth": float(actual[1]),
        },
        "target_error_nm": {
            "mask_opening": float(actual[0] - target[0]),
            "etch_depth": float(actual[1] - target[1]),
        },
        "local_prediction_error_nm": {
            "mask_opening": float(actual[0] - local_prediction[0]),
            "etch_depth": float(actual[1] - local_prediction[1]),
        },
        "trajectory_contract": {
            "status": audit.get("status"),
            "physical_time_s": reached_time,
            "maximum_material_ledger_residual_units_m2": ledger_maximum,
            "maximum_neutral_radiosity_relative_balance_error": radiosity_maximum,
            "within_declared_scope": validity_pass,
            "pass": numerical_pass,
        },
        "response_gates": {
            "target_within_one_fine_cell": target_cell_pass,
            "inside_committed_model_error_envelope": model_envelope_pass,
            "scaled_l2_merit": scaled_merit,
            "prior_r17_merit": M0,
            "committed_predicted_merit": M_PRED,
            "trust_ratio_rho": rho,
        },
        "decision": decision,
        "next_action": (
            "choose the predeclared uniform-5-nm versus certified-AMR authority path"
            if strong_pass else
            "do not run another 10 nm candidate; follow the R1.9 discrepancy branch"),
    }
    payload["evaluation_sha256"] = _canonical_sha(payload, "evaluation_sha256")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--launch-manifest", required=True)
    parser.add_argument("--run-audit", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = evaluate(args.launch_manifest, args.run_audit)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps({
        "decision": payload["decision"],
        "actual_nm": payload["actual_nm"],
        "trust_ratio_rho": payload["response_gates"]["trust_ratio_rho"],
        "next_action": payload["next_action"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
