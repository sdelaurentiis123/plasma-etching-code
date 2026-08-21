#!/usr/bin/env python3
"""Explore measured machine features associated with Bosch calibration residuals.

Only brokered calibration outcomes are opened.  Every candidate feature comes
from the independently extracted, outcome-free process summary.  The audit is
for model-form discovery: feature ranking over the same calibration set incurs
selection bias and cannot itself authorize a heldout reveal.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_bosch_wall_conditioning_calibration import (  # noqa: E402
    OUTPUT as CONDITIONING_FIT,
    _arrays,
    _inputs,
    _predict,
)


SUMMARY = (
    ROOT / "data" / "experimental" / "zenodo_17122442"
    / "process_wafer_summary.csv"
)
OUTPUT = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_wall_conditioning_depth_extension_v5"
    / "calibration_residual_feature_audit.json"
)
IDENTITY_COLUMNS = {"experiment_key", "source_group", "process_date"}


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _pearson(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if np.std(left) <= 0.0 or np.std(right) <= 0.0:
        return 0.0
    return float(np.corrcoef(left, right)[0, 1])


def _center_within_lot(values, lots):
    values = np.asarray(values, dtype=float)
    output = np.empty_like(values)
    for lot in sorted(set(lots)):
        selected = lots == lot
        output[selected] = values[selected] - np.mean(values[selected])
    return output


def _detrend_within_lot(values, wafer_numbers, lots):
    """Remove an independent intercept and linear sequence slope in every lot."""
    values = np.asarray(values, dtype=float)
    wafer_numbers = np.asarray(wafer_numbers, dtype=float)
    output = np.empty_like(values)
    for lot in sorted(set(lots)):
        selected = lots == lot
        design = np.column_stack((
            np.ones(np.count_nonzero(selected)),
            wafer_numbers[selected],
        ))
        coefficient = np.linalg.lstsq(
            design, values[selected], rcond=None)[0]
        output[selected] = values[selected] - design @ coefficient
    return output


def _per_lot_sequence_slopes(values, wafer_numbers, lots):
    values = np.asarray(values, dtype=float)
    wafer_numbers = np.asarray(wafer_numbers, dtype=float)
    output = {}
    for lot in sorted(set(lots)):
        selected = lots == lot
        design = np.column_stack((
            np.ones(np.count_nonzero(selected)),
            wafer_numbers[selected],
        ))
        coefficient = np.linalg.lstsq(
            design, values[selected], rcond=None)[0]
        output[str(int(lot))] = float(coefficient[1])
    return output


def _load_features(keys):
    with SUMMARY.open("r", newline="", encoding="utf-8") as stream:
        rows = {row["experiment_key"]: row for row in csv.DictReader(stream)}
    if not set(keys).issubset(rows):
        raise RuntimeError("process summary is missing a calibration key")
    first = rows[keys[0]]
    names = tuple(name for name in first if name not in IDENTITY_COLUMNS)
    features = {
        name: np.asarray([float(rows[key][name]) for key in keys], dtype=float)
        for name in names
    }
    return features


def _lolo_residual_correction(feature, residual, prediction, observation, lots):
    corrected = np.zeros_like(prediction)
    for lot in sorted(set(lots)):
        train = lots != lot
        test = lots == lot
        x_train = feature[train]
        scale = float(np.std(x_train))
        if scale <= 0.0:
            fitted = np.full(np.count_nonzero(test), np.mean(residual[train]))
        else:
            center = float(np.mean(x_train))
            design = np.column_stack((
                np.ones(np.count_nonzero(train)),
                (x_train - center) / scale,
            ))
            coefficient = np.linalg.lstsq(design, residual[train], rcond=None)[0]
            fitted = (
                coefficient[0]
                + coefficient[1] * (feature[test] - center) / scale)
        corrected[test] = prediction[test] - fitted
    return {
        "corrected_mae": float(np.mean(np.abs(corrected - observation))),
        "corrected_mape_percent": float(100.0 * np.mean(np.abs(
            corrected / observation - 1.0))),
    }


def build():
    fit = json.loads(CONDITIONING_FIT.read_text(encoding="utf-8"))
    coefficient = fit["conditioning_coefficients"]
    coefficients = np.asarray([
        coefficient["log_carbon_cycle"],
        coefficient["silicon_precondition"],
        coefficient["silicon_oxide_precondition"],
    ])
    measurements, traces, lot_type_by_key = _inputs()
    predictions, _law = _predict(
        coefficients, measurements, traces, lot_type_by_key)
    observed_si, observed_oxide, predicted_si, predicted_oxide = _arrays(
        measurements, predictions)
    observed_si_mean = np.mean(observed_si, axis=1)
    predicted_si_mean = np.mean(predicted_si, axis=1)
    observed_oxide_mean = np.mean(observed_oxide, axis=1)
    predicted_oxide_mean = np.mean(predicted_oxide, axis=1)
    observed_shape = observed_si / observed_si_mean[:, None]
    predicted_shape = predicted_si / predicted_si_mean[:, None]
    targets = {
        "silicon_mean_residual_um": predicted_si_mean - observed_si_mean,
        "oxide_mean_residual_um": predicted_oxide_mean - observed_oxide_mean,
        "normalized_shape_rmse_percent": 100.0 * np.sqrt(np.mean(
            (predicted_shape - observed_shape) ** 2, axis=1)),
    }
    keys = tuple(item.experiment_key for item in measurements)
    lots = np.asarray([item.lot_number for item in measurements])
    features = _load_features(keys)
    wafer_numbers = features["wafer_number"]

    rankings = {}
    sequence_detrended_rankings = {}
    for target_name, target in targets.items():
        centered_target = _center_within_lot(target, lots)
        detrended_target = _detrend_within_lot(
            target, wafer_numbers, lots)
        rows = []
        detrended_rows = []
        for feature_name, feature in features.items():
            if np.std(feature) <= 0.0:
                continue
            centered_feature = _center_within_lot(feature, lots)
            rows.append({
                "feature": feature_name,
                "raw_pearson": _pearson(feature, target),
                "raw_spearman": _pearson(rankdata(feature), rankdata(target)),
                "within_lot_pearson": _pearson(
                    centered_feature, centered_target),
                "feature_minimum": float(np.min(feature)),
                "feature_maximum": float(np.max(feature)),
            })
            detrended_feature = _detrend_within_lot(
                feature, wafer_numbers, lots)
            if np.std(detrended_feature) > 0.0:
                detrended_rows.append({
                    "feature": feature_name,
                    "sequence_partial_pearson": _pearson(
                        detrended_feature, detrended_target),
                    "feature_minimum": float(np.min(feature)),
                    "feature_maximum": float(np.max(feature)),
                })
        rows.sort(key=lambda row: abs(row["within_lot_pearson"]), reverse=True)
        detrended_rows.sort(
            key=lambda row: abs(row["sequence_partial_pearson"]),
            reverse=True)
        rankings[target_name] = rows[:20]
        sequence_detrended_rankings[target_name] = detrended_rows[:20]

    correction = {}
    for outcome_name, residual, prediction, observation in (
        ("silicon_mean", targets["silicon_mean_residual_um"],
         predicted_si_mean, observed_si_mean),
        ("oxide_mean", targets["oxide_mean_residual_um"],
         predicted_oxide_mean, observed_oxide_mean),
    ):
        rows = []
        for feature_name, feature in features.items():
            if np.std(feature) <= 0.0:
                continue
            rows.append({
                "feature": feature_name,
                **_lolo_residual_correction(
                    feature, residual, prediction, observation, lots),
            })
        rows.sort(key=lambda row: row["corrected_mae"])
        correction[outcome_name] = rows[:20]

    return {
        "schema": "petch-zenodo-bosch-calibration-residual-feature-audit-v2",
        "status": "calibration-only exploratory model-form discovery; not a prediction seal",
        "calibration_wafer_count": len(measurements),
        "candidate_feature_count": len(features),
        "feature_source_contains_outcomes": False,
        "conditioning_fit_used": str(CONDITIONING_FIT.relative_to(ROOT)),
        "residual_feature_rankings": rankings,
        "sequence_detrended_residual_feature_rankings": (
            sequence_detrended_rankings),
        "per_lot_sequence_slopes_um_per_wafer": {
            "observed_silicon_mean": _per_lot_sequence_slopes(
                observed_si_mean, wafer_numbers, lots),
            "predicted_silicon_mean": _per_lot_sequence_slopes(
                predicted_si_mean, wafer_numbers, lots),
            "silicon_mean_residual": _per_lot_sequence_slopes(
                targets["silicon_mean_residual_um"], wafer_numbers, lots),
            "observed_oxide_mean": _per_lot_sequence_slopes(
                observed_oxide_mean, wafer_numbers, lots),
            "predicted_oxide_mean": _per_lot_sequence_slopes(
                predicted_oxide_mean, wafer_numbers, lots),
        },
        "exploratory_univariate_leave_one_lot_out_corrections": correction,
        "selection_bias_warning": (
            "the best feature was selected after comparing all extracted machine "
            "summaries on calibration outcomes; a physical feature family must be "
            "frozen before any formal refit or heldout reveal"),
        "heldout_outcomes_read": False,
        "eligible_for_prediction_seal": False,
        "input_hashes": {
            "process_summary": _hash(SUMMARY),
            "conditioning_fit": _hash(CONDITIONING_FIT),
        },
    }


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    rendered = _render(build())
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
