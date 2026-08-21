#!/usr/bin/env python3
"""Score the already sealed Bosch v8 prediction on revealed heldout wafers."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.bosch_process_data import load_bosch_wafer_measurements_89pt  # noqa: E402
from scripts.audit_bosch_dynamic_wall_calibration import _Prediction  # noqa: E402
from scripts.audit_bosch_wall_conditioning_calibration import _metrics  # noqa: E402
from scripts.extract_bosch_heldout_measurements_after_seal import (  # noqa: E402
    MANIFEST as REVEAL_MANIFEST,
    OUTPUT as REVEALED_OUTCOMES,
)
from scripts.seal_bosch_wafer_boundary_map_heldout_prediction import (  # noqa: E402
    OUTPUT as PREDICTION,
    SEAL as PREDICTION_SEAL,
    V1_PREREGISTRATION,
)


DIRECTORY = PREDICTION.parent
OUTPUT = DIRECTORY / "heldout_score.json"
BOOTSTRAP_SEED = 20260821
BOOTSTRAP_REPLICATES = 20000


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _prediction_objects(rows):
    return tuple(
        _Prediction(
            silicon_depth_m=np.asarray(
                row["predicted_silicon_depth_um"], dtype=float) * 1.0e-6,
            oxide_loss_m=np.asarray(
                row["predicted_oxide_loss_um"], dtype=float) * 1.0e-6,
        )
        for row in rows
    )


def _silicon_baseline_metrics(measurements, prediction, baselines):
    observed = np.stack([item.silicon_depth_um for item in measurements])
    physics = np.stack([
        item.silicon_depth_m for item in prediction]) * 1.0e6
    observed_mean = np.mean(observed, axis=1)
    physics_mean = np.mean(physics, axis=1)
    global_mean = float(baselines["global_mean_silicon_depth_um"])
    mean_map = np.asarray(baselines["mean_silicon_depth_map_um"], dtype=float)
    observed_shape = observed / observed_mean[:, None]
    physics_shape = physics / physics_mean[:, None]
    mean_map_shape = mean_map / np.mean(mean_map)
    global_mean_mae = float(np.mean(np.abs(global_mean - observed_mean)))
    global_mean_mape = float(100.0 * np.mean(np.abs(
        global_mean / observed_mean - 1.0)))
    mean_map_point = float(np.sqrt(np.mean((mean_map - observed) ** 2)))
    mean_map_shape = float(100.0 * np.sqrt(np.mean(
        (mean_map_shape - observed_shape) ** 2)))
    physics_mean_mae = float(np.mean(np.abs(physics_mean - observed_mean)))
    physics_point = float(np.sqrt(np.mean((physics - observed) ** 2)))
    physics_shape_error = float(100.0 * np.sqrt(np.mean(
        (physics_shape - observed_shape) ** 2)))
    return {
        "global_mean_depth": {
            "silicon_mean_mae_um": global_mean_mae,
            "silicon_mean_mape_percent": global_mean_mape,
        },
        "calibration_mean_map": {
            "silicon_point_rmse_um": mean_map_point,
            "normalized_shape_rmse_percent": mean_map_shape,
        },
        "physics_improvement": {
            "mean_mae_um": global_mean_mae - physics_mean_mae,
            "point_rmse_um": mean_map_point - physics_point,
            "shape_rmse_percentage_points": (
                mean_map_shape - physics_shape_error),
        },
        "gates": {
            "physics_mean_beats_global_mean": (
                physics_mean_mae < global_mean_mae),
            "physics_point_beats_calibration_mean_map": (
                physics_point < mean_map_point),
            "physics_shape_beats_calibration_mean_map": (
                physics_shape_error < mean_map_shape),
        },
    }


def _drift(measurements, predictions):
    rows = []
    for process_date in sorted({item.experiment_key[:10] for item in measurements}):
        selected = [
            index for index, item in enumerate(measurements)
            if item.experiment_key.startswith(process_date)
        ]
        wafer = np.asarray([
            measurements[index].wafer_number for index in selected], dtype=float)
        observed = np.asarray([
            measurements[index].wafer_mean_silicon_depth_um
            for index in selected])
        predicted = np.asarray([
            np.mean(predictions[index].silicon_depth_m) * 1.0e6
            for index in selected])
        design = np.column_stack((wafer, np.ones_like(wafer)))
        observed_slope = float(np.linalg.lstsq(
            design, observed, rcond=None)[0][0])
        predicted_slope = float(np.linalg.lstsq(
            design, predicted, rcond=None)[0][0])
        rows.append({
            "process_date": process_date,
            "measured_wafer_count": len(selected),
            "observed_um_per_wafer": observed_slope,
            "predicted_um_per_wafer": predicted_slope,
            "residual_um_per_wafer": predicted_slope - observed_slope,
        })
    physics = float(np.mean([
        abs(row["residual_um_per_wafer"]) for row in rows]))
    zero = float(np.mean([
        abs(row["observed_um_per_wafer"]) for row in rows]))
    return {
        "lots": rows,
        "physics_slope_mae_um_per_wafer": physics,
        "zero_slope_baseline_mae_um_per_wafer": zero,
        "physics_beats_zero_slope_baseline": physics < zero,
    }


def _bootstrap(measurements, predictions, baselines):
    random = np.random.default_rng(BOOTSTRAP_SEED)
    count = len(measurements)
    names = (
        "silicon_mean_mae_um",
        "silicon_mean_mape_percent",
        "silicon_point_rmse_um",
        "normalized_shape_rmse_percent",
        "oxide_mean_mae_um",
        "selectivity_mape_percent",
        "mean_baseline_improvement_um",
        "point_baseline_improvement_um",
        "shape_baseline_improvement_percentage_points",
    )
    samples = {name: np.empty(BOOTSTRAP_REPLICATES) for name in names}
    for replicate in range(BOOTSTRAP_REPLICATES):
        selected = random.integers(0, count, size=count)
        selected_measurements = tuple(measurements[index] for index in selected)
        selected_predictions = tuple(predictions[index] for index in selected)
        metrics = _metrics(selected_measurements, selected_predictions)
        baseline = _silicon_baseline_metrics(
            selected_measurements, selected_predictions, baselines)
        for name in names[:6]:
            samples[name][replicate] = metrics[name]
        improvement = baseline["physics_improvement"]
        samples["mean_baseline_improvement_um"][replicate] = improvement[
            "mean_mae_um"]
        samples["point_baseline_improvement_um"][replicate] = improvement[
            "point_rmse_um"]
        samples[
            "shape_baseline_improvement_percentage_points"][replicate] = (
                improvement["shape_rmse_percentage_points"])
    return {
        "unit": "wafer",
        "seed": BOOTSTRAP_SEED,
        "replicates": BOOTSTRAP_REPLICATES,
        "percentile_interval": [2.5, 97.5],
        "intervals": {
            name: {
                "p2p5": float(np.percentile(values, 2.5)),
                "median": float(np.percentile(values, 50.0)),
                "p97p5": float(np.percentile(values, 97.5)),
            }
            for name, values in samples.items()
        },
    }


def build_score():
    preregistration = _load(V1_PREREGISTRATION)
    prediction = _load(PREDICTION)
    seal = _load(PREDICTION_SEAL)
    reveal = _load(REVEAL_MANIFEST)
    prediction_hash = _hash(PREDICTION)
    if (
        seal["prediction_sha256"] != prediction_hash
        or seal["eligible_for_separate_outcome_score_after_commit"] is not True
        or prediction["target_firewall"]["heldout_outcomes_read"] is not False
        or reveal["prediction_sha256"] != prediction_hash
        or reveal["output_sha256"] != _hash(REVEALED_OUTCOMES)
        or reveal["numeric_outcome_fields_parsed_by_broker"] is not False
    ):
        raise RuntimeError("Bosch heldout score lacks a valid prior seal")
    revealed_keys = tuple(reveal["revealed_experiment_keys"])
    measurements = load_bosch_wafer_measurements_89pt(
        REVEALED_OUTCOMES, allowed_experiment_keys=revealed_keys)
    row_by_key = {
        row["experiment_key"]: row for row in prediction["predictions"]
    }
    if not set(revealed_keys).issubset(row_by_key):
        raise RuntimeError("sealed prediction is missing a measured heldout wafer")
    rows = tuple(row_by_key[item.experiment_key] for item in measurements)
    official = prediction["official_map"]
    if any(
        not np.array_equal(item.x_um, np.asarray(official["x_um"]))
        or not np.array_equal(item.y_um, np.asarray(official["y_um"]))
        for item in measurements
    ):
        raise RuntimeError("heldout measurement map differs from sealed coordinates")
    predictions = _prediction_objects(rows)
    metrics = _metrics(measurements, predictions)
    gates = preregistration["heldout_score"]["absolute_acceptance"]
    metric_to_gate = {
        "silicon_mean_mae_um": "wafer_mean_si_depth_mae_um_max",
        "silicon_mean_mape_percent": "wafer_mean_si_depth_mape_percent_max",
        "silicon_point_rmse_um": "pointwise_si_depth_rmse_um_max",
        "normalized_shape_rmse_percent": (
            "normalized_radial_shape_rmse_percent_max"),
        "oxide_mean_mae_um": "wafer_mean_oxide_loss_mae_um_max",
        "selectivity_mape_percent": "selectivity_mape_percent_max",
    }
    absolute_gates = {
        name: metrics[name] <= gates[gate]
        for name, gate in metric_to_gate.items()
    }
    baselines = _silicon_baseline_metrics(
        measurements, predictions, prediction["calibration_only_baselines"])
    drift = _drift(measurements, predictions)
    bootstrap = _bootstrap(
        measurements, predictions, prediction["calibration_only_baselines"])
    passes = all(absolute_gates.values()) and all(baselines["gates"].values())
    return {
        "schema": "petch-zenodo-bosch-v8-heldout-score-v1",
        "status": (
            "PASS: sealed v8 physics beats every absolute and empirical heldout gate"
            if passes else
            "FAIL: sealed v8 physics misses at least one frozen heldout gate"
        ),
        "measured_heldout_wafer_count": len(measurements),
        "unmeasured_heldout_process_record_count": (
            prediction["heldout_process_record_count"] - len(measurements)),
        "revealed_experiment_keys": list(revealed_keys),
        "missing_heldout_process_keys": reveal[
            "missing_heldout_process_keys"],
        "physics_metrics": metrics,
        "frozen_absolute_thresholds": gates,
        "absolute_gates": absolute_gates,
        "all_absolute_gates_pass": all(absolute_gates.values()),
        "empirical_baselines": baselines,
        "all_empirical_baseline_gates_pass": all(
            baselines["gates"].values()),
        "within_lot_depth_drift": drift,
        "wafer_bootstrap": bootstrap,
        "all_frozen_heldout_gates_pass": passes,
        "target_firewall": {
            "prediction_committed_and_pushed_before_numeric_reveal": True,
            "prediction_sha256": prediction_hash,
            "heldout_prediction_changed_after_reveal": False,
            "heldout_outcomes_read": True,
            "execution_heldout_not_preexposure_blind": True,
        },
        "claim_boundary": preregistration["claim_boundary"],
        "input_hashes": {
            "v1_preregistration": _hash(V1_PREREGISTRATION),
            "heldout_prediction": prediction_hash,
            "heldout_prediction_seal": _hash(PREDICTION_SEAL),
            "revealed_heldout_outcomes": _hash(REVEALED_OUTCOMES),
            "reveal_manifest": _hash(REVEAL_MANIFEST),
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build_score())
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("Bosch v8 heldout score is stale")
        print("Bosch v8 heldout score is current")
        return 0
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(OUTPUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
