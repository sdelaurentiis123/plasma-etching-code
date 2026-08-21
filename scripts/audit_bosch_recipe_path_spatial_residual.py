#!/usr/bin/env python3
"""Calibration-only spatial residual discovery after Bosch recipe memory v7.

This audit asks whether the remaining 89-point map error is primarily a fixed
tool fingerprint or wafer-varying process response.  It uses outcome-space
linear proxies only to select the next physical boundary model.  None of the
proxies is a reactor correction, prediction seal, or heldout forecast.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_bosch_dynamic_wall_calibration import (  # noqa: E402
    BoschExactWallResponseTable,
)
from scripts.audit_bosch_recipe_path_memory_calibration import (  # noqa: E402
    OUTPUT as V7_FIT,
    RESPONSE_TABLE as V7_RESPONSE_TABLE,
    V7_PREREGISTRATION,
    _inputs,
    _predict,
)


SUMMARY = (
    ROOT / "data" / "experimental" / "zenodo_17122442"
    / "process_wafer_summary.csv"
)
OUTPUT = V7_FIT.parent / "spatial_residual_discovery.json"
IDENTITY_COLUMNS = {"experiment_key", "source_group", "process_date"}
PHYSICAL_FEATURE_FAMILY = (
    "c4f8_platen_peak_to_peak_rms",
    "c4f8_platen_peak_to_peak_q50",
    "c4f8_platen_peak_to_peak_q90",
    "c4f8_platen_peak_to_peak_mean",
    "sf6_platen_peak_to_peak_q50",
    "sf6_platen_reflected_power_mean",
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load_features(keys):
    with SUMMARY.open("r", newline="", encoding="utf-8") as stream:
        rows = {row["experiment_key"]: row for row in csv.DictReader(stream)}
    if not set(keys).issubset(rows):
        raise RuntimeError("process summary is missing a calibration key")
    names = tuple(
        name for name in rows[keys[0]] if name not in IDENTITY_COLUMNS)
    return {
        name: np.asarray([float(rows[key][name]) for key in keys])
        for name in names
    }


def _metrics(corrected_shape, observed_shape, predicted_mean_um,
             observed_depth_um):
    corrected_depth = corrected_shape * predicted_mean_um[:, None]
    return {
        "normalized_shape_rmse_percent": float(
            100.0 * np.sqrt(np.mean((corrected_shape - observed_shape) ** 2))),
        "silicon_point_rmse_um": float(np.sqrt(np.mean(
            (corrected_depth - observed_depth_um) ** 2))),
    }


def _whole_lot_spatial_proxy(residual_shape, predicted_shape, observed_shape,
                             predicted_mean_um, observed_depth_um, lots,
                             feature=None):
    corrected = np.zeros_like(predicted_shape)
    coefficient_maps = []
    for lot in sorted(set(lots)):
        train = lots != lot
        test = lots == lot
        if feature is None:
            coefficient = np.mean(residual_shape[train], axis=0)[None, :]
            predicted_residual = np.repeat(
                coefficient, np.count_nonzero(test), axis=0)
        else:
            center = float(np.mean(feature[train]))
            scale = float(np.std(feature[train]))
            if scale <= 0.0:
                raise RuntimeError("spatial proxy feature has zero fold variance")
            standardized = (feature - center) / scale
            design = np.column_stack((
                np.ones(np.count_nonzero(train)), standardized[train]))
            coefficient = np.linalg.lstsq(
                design, residual_shape[train], rcond=None)[0]
            predicted_residual = np.column_stack((
                np.ones(np.count_nonzero(test)), standardized[test])) @ coefficient
        coefficient_maps.append(coefficient)
        corrected[test] = predicted_shape[test] - predicted_residual
    return (
        _metrics(corrected, observed_shape, predicted_mean_um, observed_depth_um),
        np.stack(coefficient_maps),
    )


def _zernike_radial(radial_order, azimuthal_order, radius):
    output = np.zeros_like(radius)
    for index in range((radial_order - azimuthal_order) // 2 + 1):
        output += (
            (-1) ** index * math.factorial(radial_order - index)
            / (
                math.factorial(index)
                * math.factorial(
                    (radial_order + azimuthal_order) // 2 - index)
                * math.factorial(
                    (radial_order - azimuthal_order) // 2 - index)
            )
            * radius ** (radial_order - 2 * index)
        )
    return output


def _zernike_design(maximum_order, radius, phi):
    columns = []
    names = []
    for radial_order in range(1, int(maximum_order) + 1):
        for azimuthal_order in range(
                radial_order % 2, radial_order + 1, 2):
            radial = _zernike_radial(
                radial_order, azimuthal_order, radius)
            if azimuthal_order == 0:
                values = (("cos", radial),)
            else:
                values = (
                    ("cos", radial * np.cos(azimuthal_order * phi)),
                    ("sin", radial * np.sin(azimuthal_order * phi)),
                )
            for phase, value in values:
                columns.append(value - np.mean(value))
                names.append({
                    "radial_order": radial_order,
                    "azimuthal_order": azimuthal_order,
                    "phase": phase,
                })
    return np.column_stack(columns), names


def _whole_lot_zernike_proxy(maximum_order, residual_shape, predicted_shape,
                             observed_shape, predicted_mean_um,
                             observed_depth_um, lots, radius, phi):
    design, names = _zernike_design(maximum_order, radius, phi)
    corrected = np.zeros_like(predicted_shape)
    coefficients = []
    condition_numbers = []
    for lot in sorted(set(lots)):
        train = lots != lot
        test = lots == lot
        repeated = np.tile(design, (np.count_nonzero(train), 1))
        coefficient = np.linalg.lstsq(
            repeated, residual_shape[train].reshape(-1), rcond=None)[0]
        corrected[test] = predicted_shape[test] - design @ coefficient
        coefficients.append(coefficient)
        condition_numbers.append(float(np.linalg.cond(repeated)))
    coefficient_array = np.stack(coefficients)
    return {
        "maximum_order": int(maximum_order),
        "coefficient_count": len(names),
        "metrics": _metrics(
            corrected, observed_shape, predicted_mean_um, observed_depth_um),
        "maximum_design_condition_number": max(condition_numbers),
        "maximum_absolute_coefficient": float(np.max(np.abs(coefficient_array))),
        "median_fold_coefficient_standard_deviation": float(np.median(
            np.std(coefficient_array, axis=0))),
    }


def _pearson(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if np.std(left) == 0.0 or np.std(right) == 0.0:
        return 0.0
    return float(np.corrcoef(left, right)[0, 1])


def build():
    fit = json.loads(V7_FIT.read_text(encoding="utf-8"))
    if fit["heldout_outcomes_read"] is not False:
        raise RuntimeError("Bosch v7 heldout firewall is not closed")
    parameters = fit["parameters"]
    coefficients = np.asarray([
        parameters["conditioning_repeat_log_coefficient"],
        parameters["silicon_precondition_coefficient"],
        parameters["silicon_oxide_precondition_coefficient"],
        parameters["log_wall_loss_per_reference_wafer"],
    ])
    measurements, process_traces, lot_type_by_date = _inputs()
    table = BoschExactWallResponseTable.load(V7_RESPONSE_TABLE)
    predictions, *_ = _predict(
        coefficients, measurements, process_traces, lot_type_by_date, table)
    observed_depth_um = np.stack([
        item.silicon_depth_um for item in measurements])
    predicted_depth_um = np.stack([
        item.silicon_depth_m for item in predictions]) * 1.0e6
    observed_mean_um = np.mean(observed_depth_um, axis=1)
    predicted_mean_um = np.mean(predicted_depth_um, axis=1)
    observed_shape = observed_depth_um / observed_mean_um[:, None]
    predicted_shape = predicted_depth_um / predicted_mean_um[:, None]
    residual_shape = predicted_shape - observed_shape
    mean_residual = np.mean(residual_shape, axis=0)
    centered_residual = residual_shape - mean_residual
    singular_values = np.linalg.svd(centered_residual, compute_uv=False)
    singular_fraction = singular_values ** 2 / np.sum(singular_values ** 2)
    lots = np.asarray([item.lot_number for item in measurements])
    keys = tuple(item.experiment_key for item in measurements)
    features = _load_features(keys)

    static_metrics, static_coefficients = _whole_lot_spatial_proxy(
        residual_shape, predicted_shape, observed_shape, predicted_mean_um,
        observed_depth_um, lots)
    feature_rows = []
    for name, feature in features.items():
        if np.std(feature) <= 0.0:
            continue
        metrics, maps = _whole_lot_spatial_proxy(
            residual_shape, predicted_shape, observed_shape,
            predicted_mean_um, observed_depth_um, lots, feature=feature)
        feature_rows.append({
            "feature": name,
            "metrics": metrics,
            "feature_minimum": float(np.min(feature)),
            "feature_maximum": float(np.max(feature)),
            "maximum_absolute_standardized_slope": float(np.max(np.abs(
                maps[:, 1]))),
        })
    feature_rows.sort(
        key=lambda row: row["metrics"]["normalized_shape_rmse_percent"])

    selected_feature = features["c4f8_platen_peak_to_peak_rms"]
    standardized = (
        selected_feature - np.mean(selected_feature)
    ) / np.std(selected_feature)
    feature_design = np.column_stack((
        np.ones(len(selected_feature)), standardized))
    feature_maps = np.linalg.lstsq(
        feature_design, residual_shape, rcond=None)[0]
    slope_map = feature_maps[1]
    radius = np.hypot(measurements[0].x_um, measurements[0].y_um) / 1.0e5
    phi = np.arctan2(measurements[0].y_um, measurements[0].x_um)
    edge_basis = 1.0 / (1.0 + np.exp(-(radius - 0.90) / 0.02))
    edge_basis -= np.mean(edge_basis)
    inner = radius <= 0.75
    outer = radius >= 0.95

    zernike_rows = [
        _whole_lot_zernike_proxy(
            order, residual_shape, predicted_shape, observed_shape,
            predicted_mean_um, observed_depth_um, lots, radius, phi)
        for order in range(1, 11)
    ]
    fixed_fraction = float(
        len(measurements) * np.sum(mean_residual ** 2)
        / np.sum(residual_shape ** 2)
    )
    return {
        "schema": "petch-zenodo-bosch-v7-spatial-residual-discovery-v1",
        "status": "calibration-only model-form discovery; no physical boundary refit or prediction seal",
        "calibration_wafer_count": len(measurements),
        "points_per_wafer": residual_shape.shape[1],
        "decomposition": {
            "v7_raw_normalized_shape_rmse_percent": float(
                100.0 * np.sqrt(np.mean(residual_shape ** 2))),
            "fixed_mean_residual_rmse_percent": float(
                100.0 * np.sqrt(np.mean(mean_residual ** 2))),
            "within_wafer_set_residual_rmse_after_all_data_fixed_map_percent": (
                float(100.0 * np.sqrt(np.mean(centered_residual ** 2)))),
            "fraction_of_squared_residual_in_fixed_mean_map": fixed_fraction,
            "centered_residual_singular_variance_fractions_first_five": (
                singular_fraction[:5].tolist()),
            "interpretation": (
                "the dominant missing term is one shared spatial equipment "
                "boundary map, not another scalar depth normalization"),
        },
        "whole_lot_output_space_proxies": {
            "shared_89_point_intercept": {
                "independent_parameter_count": 88,
                "metrics": static_metrics,
                "maximum_fold_map_difference": float(np.max(
                    np.ptp(static_coefficients[:, 0, :], axis=0))),
            },
            "intercept_plus_one_standardized_feature_slope_rankings": (
                feature_rows[:25]),
            "selection_bias_warning": (
                "131 outcome-free machine features were ranked after examining "
                "calibration residuals; the physical feature family and operator "
                "must be frozen before formal refitting"),
        },
        "selected_physical_family_diagnostic": {
            "features": list(PHYSICAL_FEATURE_FAMILY),
            "representative_feature": "c4f8_platen_peak_to_peak_rms",
            "representative_feature_mean_V": float(np.mean(selected_feature)),
            "representative_feature_standard_deviation_V": float(
                np.std(selected_feature)),
            "standardized_slope_map_rmse_percent": float(
                100.0 * np.sqrt(np.mean(slope_map ** 2))),
            "standardized_slope_map_edge_basis_pearson": _pearson(
                slope_map, edge_basis),
            "standardized_slope_outer_mean_percent": float(
                100.0 * np.mean(slope_map[outer])),
            "standardized_slope_inner_mean_percent": float(
                100.0 * np.mean(slope_map[inner])),
            "interpretation": (
                "higher C4F8-phase platen Vpp produces an edge-positive, "
                "center-negative residual mode, supporting a voltage-dependent "
                "current-conserving sheath/focus-ring transmission closure"),
        },
        "complete_zernike_output_space_capacity": zernike_rows,
        "zernike_warning": (
            "these are outcome-space capacity probes, not ion-flux fits; high "
            "orders approaching the 89-point mean map are not automatically "
            "acceptable physics"),
        "heldout_outcomes_read": False,
        "heldout_prediction_written": False,
        "eligible_for_prediction_seal": False,
        "surface_laws_changed": False,
        "positive_ion_boundary_changed": False,
        "input_hashes": {
            "v7_preregistration": _hash(V7_PREREGISTRATION),
            "v7_fit": _hash(V7_FIT),
            "v7_response_table": _hash(V7_RESPONSE_TABLE),
            "process_summary": _hash(SUMMARY),
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
