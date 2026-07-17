#!/usr/bin/env python3
"""Build the bounded, development-only Krueger R1.9 calibration proposal.

This program consumes only the two base-case calibration observables.  It deliberately refuses to
open the transfer-observation table or any held-out profile.  The local response matrix comes from
the checksum-bound final 10 nm axisymmetric secant, while the residual comes from the completed
mixed-operator 5 nm development trajectory.  Consequently, its output may select one *development*
run, but it cannot freeze parameters or support a validation claim.

R1.9 authorizes exactly the precommitted bounded response-model check and no second 10 nm candidate.
The preceding R1.8/WP5 analysis called for a trust-region proposal but did not preregister a
numerical radius.
Rather than disguise a new radius as physics, the trust box is scaled by the absolute coordinate
changes of the immediately preceding calibration proposal that was actually evaluated by the 5 nm
development trajectory.  A unit radius in that scaled infinity norm means that no coordinate moves
farther than the last evaluated move.  If the full Newton/secant direction exceeds the box, one
scalar factor clips the entire direction, preserving its response-model direction.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md"
CAMPAIGN = ROOT / "VALIDATION_FIRST_SUPERSET_CAMPAIGN_2026-07-17.md"
BASE_TARGETS = ROOT / "data" / "experimental" / "krueger_2024" / "base_case_metrics.csv"
PARAMETERS = (
    "effective_mask_crosslinked_growth_fraction",
    "oxide_etch_yield_scale",
)
RESPONSES = ("mask_opening_nm", "etch_depth_nm")
FINE_CELL_TOLERANCE_NM = np.asarray([5.0, 5.0], dtype=float)


def _sha(path: Path | str) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _canonical_sha(payload: dict, field: str) -> str:
    canonical = dict(payload)
    canonical.pop(field, None)
    return sha256(json.dumps(
        canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _verify_embedded_sha(payload: dict, field: str) -> None:
    if payload.get(field) != _canonical_sha(payload, field):
        raise ValueError(f"{field} does not match its artifact content")


def _read_json(path: Path | str) -> tuple[Path, dict]:
    path = Path(path)
    return path, json.loads(path.read_text(encoding="utf-8"))


def _targets() -> tuple[np.ndarray, dict]:
    with BASE_TARGETS.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {
        row["metric"]: float(row["value"])
        for row in rows
        if row["split"] == "calibration"
        and row["evidence_type"] == "experiment"
        and row["metric"] in {"mask_opening", "etch_depth"}
    }
    if set(selected) != {"mask_opening", "etch_depth"}:
        raise ValueError("the two declared base calibration targets are not uniquely available")
    target = np.asarray([selected["mask_opening"], selected["etch_depth"]], dtype=float)
    return target, {
        "path_name": BASE_TARGETS.name,
        "sha256": _sha(BASE_TARGETS),
        "selection": "experiment rows with split=calibration; held-out table not opened",
        "values": dict(zip(RESPONSES, target.tolist())),
    }


def _response_model(path: Path | str) -> tuple[dict, dict]:
    path, payload = _read_json(path)
    if (payload.get("schema") != "petch.krueger-2024.base-axisymmetric-secant.v1"
            or payload.get("protocol_id") != "K24-PETCH-R1.6"
            or payload.get("held_out_profile_data_read") is not False):
        raise ValueError("response model is not the sealed base-only R1.6 axisymmetric secant")
    _verify_embedded_sha(payload, "proposal_sha256")
    if payload.get("base_target_table_sha256") != _sha(BASE_TARGETS):
        raise ValueError("response model is not bound to the current base calibration targets")
    derivation = payload.get("derivation", {})
    if (tuple(derivation.get("parameter_order", ())) != PARAMETERS
            or tuple(derivation.get("response_order", ())) != RESPONSES):
        raise ValueError("response model order is incompatible with this campaign")
    jacobian = np.asarray(derivation.get("jacobian"), dtype=float)
    if jacobian.shape != (2, 2) or not np.all(np.isfinite(jacobian)):
        raise ValueError("response model Jacobian is invalid")
    condition = float(np.linalg.cond(jacobian))
    if not np.isfinite(condition) or condition > 1e6:
        raise ValueError("response model Jacobian is singular or ill-conditioned")

    current_10 = payload.get("current_axisymmetric_endpoint", {})
    previous_parameters = np.asarray([
        float(current_10["fraction"]),
        float(current_10["oxide_etch_yield_scale"]),
    ])
    proposed_parameters = np.asarray([
        float(payload["proposed_configuration"][name]) for name in PARAMETERS], dtype=float)
    prior_update = np.asarray(derivation.get("parameter_update"), dtype=float)
    if (prior_update.shape != (2,)
            or not np.allclose(previous_parameters + prior_update, proposed_parameters,
                               rtol=0.0, atol=2e-14)):
        raise ValueError("response model's preceding proposal step is internally inconsistent")
    coordinate_scale = np.abs(prior_update)
    if np.any(~np.isfinite(coordinate_scale)) or np.any(coordinate_scale <= 0.0):
        raise ValueError("preceding evaluated step cannot define both trust coordinates")
    return payload, {
        "path_name": path.name,
        "sha256": _sha(path),
        "proposal_sha256": payload["proposal_sha256"],
        "source_protocol_id": payload["protocol_id"],
        "source_protocol_sha256": payload["protocol_sha256"],
        "jacobian": jacobian,
        "condition_number": condition,
        "preceding_proposal_update": prior_update,
        "coordinate_scale": coordinate_scale,
        "proposed_parameters": proposed_parameters,
    }


def _development_endpoint(
        endpoint_path: Path | str, receipt_path: Path | str,
        expected_parameters: np.ndarray) -> tuple[dict, dict]:
    endpoint_path, endpoint = _read_json(endpoint_path)
    receipt_path, receipt = _read_json(receipt_path)
    copied = receipt.get("copied_file_sha256", {})
    endpoint_sha = _sha(endpoint_path)
    if copied.get(endpoint_path.name) != endpoint_sha:
        raise ValueError("continuation receipt does not bind the development endpoint")
    if (receipt.get("schema")
            != "petch.krueger-2024.mixed-operator-topology-continuation.v1"
            or receipt.get("held_out_profile_data_read") is not False
            or "development evidence only" not in receipt.get("scientific_status", "")):
        raise ValueError("endpoint receipt does not classify this trajectory as development-only")

    config = endpoint.get("configuration", {})
    if (endpoint.get("status") != "complete"
            or config.get("boundary_case") != "base"
            or not np.isclose(float(config.get("duration_s", np.nan)), 60.0)
            or not np.isclose(float(config.get("dx_um", np.nan)), 0.005)
            or config.get("ion_azimuthal_closure") != "axisymmetric_uniform"
            or int(config.get("ion_azimuthal_order", -1)) != 16
            or config.get("topology_change_policy") != "continue_gas_cavity"):
        raise ValueError("not the complete R1.8 5 nm base development endpoint")
    parameters = np.asarray([float(config[name]) for name in PARAMETERS], dtype=float)
    if not np.allclose(parameters, expected_parameters, rtol=0.0, atol=2e-14):
        raise ValueError("development endpoint did not evaluate the sealed preceding proposal")
    response = np.asarray([
        float(endpoint["final_metrics"][name]) for name in RESPONSES], dtype=float)
    receipt_response = np.asarray([
        float(receipt["final_metrics"][name]) for name in RESPONSES], dtype=float)
    if not np.allclose(response, receipt_response, rtol=0.0, atol=1e-12):
        raise ValueError("receipt and endpoint responses disagree")
    return endpoint, {
        "path_name": endpoint_path.name,
        "sha256": endpoint_sha,
        "config_hash": endpoint.get("config_hash"),
        "receipt_path_name": receipt_path.name,
        "receipt_sha256": _sha(receipt_path),
        "parameters": parameters,
        "response": response,
        "operator_class": "mixed-operator 5 nm development trajectory",
    }


def _physical_direction_scale(current: np.ndarray, step: np.ndarray) -> float:
    """Largest [0, 1] scale retaining f in [0, 1] and positive oxide yield."""
    limits = [1.0]
    if step[0] > 0.0:
        limits.append((1.0 - current[0]) / step[0])
    elif step[0] < 0.0:
        limits.append(current[0] / -step[0])
    if step[1] < 0.0:
        limits.append(float(np.nextafter(current[1] / -step[1], 0.0)))
    scale = float(min(limits))
    return max(0.0, scale)


def derive(
        response_model_path: Path | str, development_endpoint_path: Path | str,
        continuation_receipt_path: Path | str) -> dict:
    _model, model_info = _response_model(response_model_path)
    _endpoint, endpoint_info = _development_endpoint(
        development_endpoint_path, continuation_receipt_path,
        model_info["proposed_parameters"])
    target, target_info = _targets()

    jacobian = model_info["jacobian"]
    current_parameters = endpoint_info["parameters"]
    current_response = endpoint_info["response"]
    current_error = current_response - target
    if np.all(np.abs(current_error) <= FINE_CELL_TOLERANCE_NM):
        raise ValueError("development endpoint already lies within one fine cell; no proposal earned")
    full_step = np.linalg.solve(jacobian, -current_error)
    full_candidate = current_parameters + full_step

    coordinate_scale = model_info["coordinate_scale"]
    normalized_full_step = np.abs(full_step) / coordinate_scale
    full_scaled_linf = float(np.max(normalized_full_step))
    if not np.isfinite(full_scaled_linf) or full_scaled_linf <= 0.0:
        raise ValueError("response model did not produce a finite nonzero proposal step")
    trust_direction_scale = min(1.0, 1.0 / full_scaled_linf)
    physical_direction_scale = _physical_direction_scale(current_parameters, full_step)
    direction_scale = float(min(trust_direction_scale, physical_direction_scale))
    safeguarded_step = direction_scale * full_step
    safeguarded_candidate = current_parameters + safeguarded_step
    if (not np.all(np.isfinite(safeguarded_candidate))
            or not 0.0 <= safeguarded_candidate[0] <= 1.0
            or safeguarded_candidate[1] <= 0.0):
        raise ValueError("safeguarded proposal lies outside physical parameter bounds")

    full_prediction = current_response + jacobian @ full_step
    safeguarded_prediction = current_response + jacobian @ safeguarded_step
    current_scaled_error = current_error / FINE_CELL_TOLERANCE_NM
    safeguarded_scaled_error = (safeguarded_prediction - target) / FINE_CELL_TOLERANCE_NM
    current_merit = float(np.linalg.norm(current_scaled_error))
    safeguarded_merit = float(np.linalg.norm(safeguarded_scaled_error))

    payload = {
        "schema": "petch.krueger-2024.development-trust-region-proposal.v1",
        "protocol_id": "K24-PETCH-R1.9",
        "protocol_sha256": _sha(PROTOCOL),
        "campaign_sha256": _sha(CAMPAIGN),
        "scientific_status": (
            "development-only safeguarded proposal from a mixed-operator 5 nm endpoint; "
            "not a parameter freeze, authoritative calibration, or validation result"),
        "authority": False,
        "held_out_profile_data_read": False,
        "generator": {
            "path_name": Path(__file__).name,
            "sha256": _sha(Path(__file__)),
        },
        "calibration_targets": target_info,
        "inputs": {
            "response_model": {
                key: value for key, value in model_info.items()
                if key not in {"jacobian", "preceding_proposal_update", "coordinate_scale",
                               "proposed_parameters"}
            },
            "development_endpoint": {
                key: value for key, value in endpoint_info.items()
                if key not in {"parameters", "response"}
            },
        },
        "orders": {
            "parameters": list(PARAMETERS),
            "responses": list(RESPONSES),
        },
        "parameter_bounds": {
            PARAMETERS[0]: {"lower": 0.0, "upper": 1.0, "source": "R1.9 protocol"},
            PARAMETERS[1]: {"lower_exclusive": 0.0, "upper": None,
                            "source": "R1.9 protocol"},
        },
        "trust_region": {
            "shape": "previous-evaluated-step-scaled-linf",
            "normalized_radius": 1.0,
            "coordinate_scale": dict(zip(PARAMETERS, coordinate_scale.tolist())),
            "scale_source": (
                "absolute parameter update of the sealed R1.6 proposal, verified to be the pair "
                "evaluated by the completed 5 nm development endpoint"),
            "policy_basis": (
                "R1.9/WP5 require a bounded local step and one proposal at a time but declare no "
                "numeric radius; using the last evaluated move is an explicit numerical safeguard, "
                "not a new physical closure"),
            "full_step_normalized_coordinates": dict(zip(
                PARAMETERS, normalized_full_step.tolist())),
            "full_step_scaled_linf_norm": full_scaled_linf,
            "trust_direction_scale": float(trust_direction_scale),
            "physical_bound_direction_scale": float(physical_direction_scale),
            "applied_direction_scale": direction_scale,
            "clipped": bool(direction_scale < 1.0),
            "direction_preserved": True,
        },
        "derivation": {
            "method": "frozen-10nm-response/full-step-plus-scaled-linf-safeguard",
            "frozen_jacobian": jacobian.tolist(),
            "condition_number": model_info["condition_number"],
            "current_parameters": dict(zip(PARAMETERS, current_parameters.tolist())),
            "current_response": dict(zip(RESPONSES, current_response.tolist())),
            "target_error_response_minus_target": dict(zip(RESPONSES, current_error.tolist())),
            "preceding_evaluated_parameter_update": dict(zip(
                PARAMETERS, model_info["preceding_proposal_update"].tolist())),
            "full_newton_secant_step": dict(zip(PARAMETERS, full_step.tolist())),
            "full_unclipped_candidate": dict(zip(PARAMETERS, full_candidate.tolist())),
            "full_predicted_response": dict(zip(RESPONSES, full_prediction.tolist())),
            "safeguarded_step": dict(zip(PARAMETERS, safeguarded_step.tolist())),
            "safeguarded_candidate": dict(zip(PARAMETERS, safeguarded_candidate.tolist())),
            "safeguarded_predicted_response": dict(zip(
                RESPONSES, safeguarded_prediction.tolist())),
            "fine_cell_tolerance_nm": dict(zip(RESPONSES, FINE_CELL_TOLERANCE_NM.tolist())),
            "current_scaled_l2_merit": current_merit,
            "safeguarded_predicted_scaled_l2_merit": safeguarded_merit,
            "predicted_scaled_l2_reduction": current_merit - safeguarded_merit,
        },
        "proposed_configuration": {
            PARAMETERS[0]: float(safeguarded_candidate[0]),
            PARAMETERS[1]: float(safeguarded_candidate[1]),
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "topology_change_policy": "continue_gas_cavity",
        },
        "next_evaluation_contract": {
            "count": 1,
            "ten_nm_sequence_closes_after_this_evaluation": True,
            "operator_status": "development proposal only until WP3/WP4 numerical gates pass",
            "acceptance_state": "pending actual response",
            "required_comparison": (
                "record actual versus predicted base-only merit reduction; accept, shrink, or grow "
                "only under an explicitly recorded policy; do not inspect held-out outcomes"),
            "authority_path": (
                "any eventual final pair requires a clean t=0 authoritative 5 nm or certified-AMR "
                "base confirmation before the sealed held-out reveal"),
        },
    }
    # Ensure every numpy value was converted before creating the embedded content checksum.
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["proposal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--response-model", required=True)
    parser.add_argument("--development-endpoint", required=True)
    parser.add_argument("--continuation-receipt", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = derive(
        args.response_model, args.development_endpoint, args.continuation_receipt)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps({
        "scientific_status": payload["scientific_status"],
        "full_unclipped_candidate": payload["derivation"]["full_unclipped_candidate"],
        "proposed_configuration": payload["proposed_configuration"],
        "applied_direction_scale": payload["trust_region"]["applied_direction_scale"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
