#!/usr/bin/env python3
"""Summarize the frozen Jeong 2023 calibration and untouched transfer predictions.

The summary keeps the calibration point out of held-out error metrics.  The paper does not report
statistical measurement uncertainty for Figure 7, so the digitization interval is treated as a
necessary image-reading check, never as the complete experimental uncertainty.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from petch.experimental_data import load_jeong_2023_etch_depths


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "jeong_2023"
DEFAULT_RESULTS = ROOT / "results" / "jeong_2023_predictive_validation"
AUTHORITATIVE_IMPLEMENTATION_FILES = (
    "scripts/jeong_2023_transfer.py",
    "src/petch/experimental_boundary.py",
    "src/petch/surface_kinetics.py",
    "src/petch/feature_step_3d.py",
)


def _atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _operator_replay_status(payload):
    stored = payload.get("provenance", {}).get(
        "implementation_checksums_sha256", {})
    current = {
        name: _sha256(ROOT / name) for name in AUTHORITATIVE_IMPLEMENTATION_FILES}
    mismatch = {
        name: {
            "stored_sha256": stored.get(name),
            "current_sha256": current[name],
        }
        for name in AUTHORITATIVE_IMPLEMENTATION_FILES
        if stored.get(name) != current[name]
    }
    return {
        "current_operator": not mismatch,
        "mismatches": mismatch,
    }


def _key(control_mode, control_value, width_nm):
    return control_mode, float(control_value), float(width_nm)


def _target_key(target):
    value = (
        target.self_bias_magnitude_v
        if target.control_mode == "ion_energy" else target.electron_density_m3)
    return _key(target.control_mode, value, target.trench_width_nm)


def _run_key(run):
    target = run["target"]
    value = (
        target["self_bias_magnitude_v"]
        if target["control_mode"] == "ion_energy" else target["electron_density_m3"])
    return _key(target["control_mode"], value, target["trench_width_nm"])


def _trend_summary(records, control_mode):
    selected = sorted(
        (item for item in records if item["control_mode"] == control_mode),
        key=lambda item: item["control_value"])
    if len(selected) < 2:
        return {
            "available": False,
            "reason": "fewer_than_two_completed_conditions",
        }
    experiment = np.asarray([item["experimental_depth_nm"] for item in selected])
    prediction = np.asarray([item["predicted_depth_nm"] for item in selected])
    experiment_delta = np.diff(experiment)
    prediction_delta = np.diff(prediction)
    return {
        "available": True,
        "control_values": [item["control_value"] for item in selected],
        "experiment_monotone_increasing": bool(np.all(experiment_delta > 0.0)),
        "prediction_monotone_increasing": bool(np.all(prediction_delta > 0.0)),
        "all_adjacent_trend_signs_match": bool(np.all(
            np.sign(experiment_delta) == np.sign(prediction_delta))),
        "experimental_endpoint_gain_nm": float(experiment[-1] - experiment[0]),
        "predicted_endpoint_gain_nm": float(prediction[-1] - prediction[0]),
    }


def build_summary(results_dir):
    results_dir = Path(results_dir)
    manifest = json.loads((results_dir / "calibration_manifest.json").read_text())
    targets = load_jeong_2023_etch_depths(DATA / "digitized_figure7_depths.csv")
    width_targets = [item for item in targets if np.isclose(item.trench_width_nm, 200.0)]

    completed = {}
    artifacts = {}
    replay_status = {}
    for path in sorted(results_dir.glob("*_width200_medium.json")):
        payload = json.loads(path.read_text())
        if payload.get("status") != "complete" or len(payload.get("runs", ())) != 1:
            continue
        run = payload["runs"][0]
        key = _run_key(run)
        if key in completed:
            raise ValueError(f"duplicate Jeong prediction for {key}")
        completed[key] = run
        replay = _operator_replay_status(payload)
        replay_status[key] = replay
        artifacts[key] = {
            "file": path.name,
            "config_hash_sha256": payload["provenance"]["config_hash_sha256"],
            "git_revision": payload["provenance"]["git_revision"],
            "operator_current": replay["current_operator"],
            "implementation_mismatches": replay["mismatches"],
        }

    anchor = manifest["calibration_anchor"]
    anchor_key = _key("ion_energy", anchor["self_bias_magnitude_v"], 200.0)
    records = []
    missing = []
    for target in sorted(
            width_targets,
            key=lambda item: (item.control_mode, _target_key(item)[1])):
        key = _target_key(target)
        control_value = key[1]
        if key == anchor_key:
            prediction_nm = manifest["refinement_evidence"]["predicted_depth_nm"]
            artifact = {
                "file": "calibration_manifest.json",
                "artifact_sha256": manifest["refinement_evidence"]["artifact_sha256"],
            }
        elif key in completed:
            prediction_nm = completed[key]["prediction"]["etch_depth_nm"]
            artifact = artifacts[key]
        else:
            missing.append({
                "control_mode": target.control_mode,
                "control_value": control_value,
                "trench_width_nm": target.trench_width_nm,
            })
            continue
        residual_nm = float(prediction_nm - target.etch_depth_nm)
        records.append({
            "control_mode": target.control_mode,
            "control_value": control_value,
            "trench_width_nm": target.trench_width_nm,
            "split": target.split,
            "role": target.role,
            "experimental_depth_nm": target.etch_depth_nm,
            "digitization_uncertainty_nm": target.digitization_uncertainty_nm,
            "predicted_depth_nm": prediction_nm,
            "residual_nm": residual_nm,
            "absolute_error_nm": abs(residual_nm),
            "within_digitization_interval": bool(
                abs(residual_nm) <= target.digitization_uncertainty_nm),
            "artifact": artifact,
        })

    held_out = [item for item in records if item["split"] != "calibration"]
    if held_out:
        residual = np.asarray([item["residual_nm"] for item in held_out])
        metrics = {
            "completed_held_out_points": len(held_out),
            "required_held_out_points_at_200nm": 5,
            "mean_absolute_error_nm": float(np.mean(np.abs(residual))),
            "root_mean_square_error_nm": float(np.sqrt(np.mean(residual ** 2))),
            "maximum_absolute_error_nm": float(np.max(np.abs(residual))),
            "within_digitization_interval_count": int(sum(
                item["within_digitization_interval"] for item in held_out)),
        }
    else:
        metrics = {
            "completed_held_out_points": 0,
            "required_held_out_points_at_200nm": 5,
        }

    evidence_complete = len(held_out) == 5 and not missing
    campaign_operator_current = bool(replay_status) and all(
        item["current_operator"] for item in replay_status.values())
    trends = {
        mode: _trend_summary(records, mode) for mode in ("ion_energy", "ion_flux")}
    if not evidence_complete:
        verdict = "incomplete"
    elif not campaign_operator_current:
        verdict = "historical_results_not_current_validation"
    elif all(item["within_digitization_interval"] for item in held_out):
        verdict = "all_held_out_points_inside_digitization_intervals"
    elif all(item["all_adjacent_trend_signs_match"] for item in trends.values()):
        verdict = "trend_transfer_passes_but_absolute_depth_misses_remain"
    else:
        verdict = "held_out_transfer_fails"

    return {
        "campaign": "jeong_2023_predictive_validation",
        "scope": "200_nm_charge_off_no_reflection_transfer_stage",
        "status": (
            "complete" if evidence_complete and campaign_operator_current
            else "historical_stale_operator" if evidence_complete
            else "partial"),
        "verdict": verdict,
        "operator_replay": {
            "authoritative_files": list(AUTHORITATIVE_IMPLEMENTATION_FILES),
            "current_operator": campaign_operator_current,
            "per_artifact": {
                "/".join(map(str, key)): value
                for key, value in replay_status.items()},
            "claim_rule": (
                "Archived scores remain historical evidence when any authoritative "
                "implementation checksum differs; they cannot certify the current operator."),
        },
        "calibration_parameter_frozen": manifest["frozen_parameter"],
        "measurement_caveat": (
            "Figure 7 statistical experimental uncertainty is not reported; the 35 nm interval "
            "quantifies digitization only and is not a total validation uncertainty."),
        "records": records,
        "missing_conditions": missing,
        "held_out_metrics": metrics,
        "trend_checks": trends,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.results_dir / "summary.json"
    payload = build_summary(args.results_dir)
    _atomic_json(output, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
