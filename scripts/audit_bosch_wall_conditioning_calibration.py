#!/usr/bin/env python3
"""Calibration-only audit of the frozen Bosch wall-conditioning closure.

This driver can read only the brokered calibration asset.  It evaluates one
shared, physical three-coefficient wall-state law through the measured-waveform
reactor, cylindrical wafer transfer, and unchanged fused Belen/La Magna
surface recurrence.  Its output is an exploratory calibration receipt, never a
heldout seal: true leave-one-lot-out refits and certification-grid refinement
remain separate required gates.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.bosch_process_data import (  # noqa: E402
    load_bosch_process_traces,
    load_bosch_wafer_measurements_89pt,
)
from petch.bosch_wafer_depth import (  # noqa: E402
    build_bosch_reference_surface_mechanisms,
)
from petch.bosch_wafer_depth_fast import (  # noqa: E402
    predict_bosch_wafer_point_depth_batch_fast,
)
from petch.reactor_global.bosch_spts_cylindrical import (  # noqa: E402
    BoschSPTSCylindricalParameters,
    DeterministicBoschSPTSCylindricalReactorToWafer,
)
from petch.reactor_global.bosch_spts_reduced import (  # noqa: E402
    BoschSPTSReducedParameters,
    BoschSPTSWallConditioningLaw,
    conditioned_bosch_spts_parameters,
)


DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
PROCESS = DATA / "Process_data.nc"
DICTIONARY = DATA / "Dictionary_process.nc"
CALIBRATION = DATA / "calibration_Si_Oxide_etch_89_points.csv"
V1_PREREGISTRATION = (
    ROOT / "results" / "curated" / "zenodo_bosch_reactor_depth_holdout_v1"
    / "preregistration.json"
)
V2_PREREGISTRATION = (
    ROOT / "results" / "curated" / "zenodo_bosch_cylindrical_depth_extension_v2"
    / "preregistration.json"
)
V5_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_wall_conditioning_depth_extension_v5"
    / "preregistration.json"
)
OUTPUT = V5_PREREGISTRATION.parent / "calibration_fit.json"
REDUCED_IMPLEMENTATION = (
    ROOT / "src" / "petch" / "reactor_global" / "bosch_spts_reduced.py")
CYLINDRICAL_IMPLEMENTATION = (
    ROOT / "src" / "petch" / "reactor_global" / "bosch_spts_cylindrical.py")
SURFACE_IMPLEMENTATION = ROOT / "src" / "petch" / "bosch_wafer_depth_fast.py"


BASE_REDUCED_PARAMETERS = BoschSPTSReducedParameters(
    f_reference_lifetime_s=0.412333985,
    film_precursor_reference_lifetime_s=2.99382399e-5,
    positive_ion_reference_lifetime_s=1.53294139e-4,
    radial_cell_count=24,
    axial_cell_count=24,
    source_ring_radius_m=0.12,
    source_radial_width_m=0.02,
    source_central_fraction=(0.0, 0.0, 0.07954333),
    diffusion_coefficient_m2_s=(2.0, 1.0, 10.0),
)
BASE_CYLINDRICAL_KEYWORDS = {
    "azimuthal_cell_count": 16,
    "ion_edge_focus_amplitude": 0.5,
    "ion_edge_focus_onset_radius_m": 0.0907157698,
    "ion_edge_focus_width_m": 0.00705256433,
    "source_cosine_coefficients": (
        (-0.0356084449,), (-0.273223495,), (0.000805789873,)),
    "source_sine_coefficients": (
        (-0.405495811,), (-1.49999895,), (-0.0191321177,)),
}
FROZEN_GATES = {
    "silicon_mean_mae_um": 1.0,
    "silicon_mean_mape_percent": 3.0,
    "silicon_point_rmse_um": 1.5,
    "normalized_shape_rmse_percent": 2.0,
    "oxide_mean_mae_um": 0.08,
    "selectivity_mape_percent": 12.0,
}


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _inputs():
    v1 = _load_json(V1_PREREGISTRATION)
    v2 = _load_json(V2_PREREGISTRATION)
    v5 = _load_json(V5_PREREGISTRATION)
    allowed = tuple(v1["split_rule"]["calibration_experiment_keys"])
    measurements = load_bosch_wafer_measurements_89pt(
        CALIBRATION, allowed_experiment_keys=allowed)
    measurement_by_key = {
        measurement.experiment_key: measurement for measurement in measurements}
    traces = {
        trace.experiment_key: trace
        for trace in load_bosch_process_traces(PROCESS, DICTIONARY)
        if trace.experiment_key in measurement_by_key
    }
    if set(traces) != set(measurement_by_key):
        raise RuntimeError("calibration process/measurement key mismatch")
    lot_type_by_date = v2["conditioning_state"]["declared_lot_types"]
    lot_type_by_key = {
        key: lot_type_by_date[key.rsplit("_", 1)[0]]
        for key in measurement_by_key
    }
    if v5["target_firewall"]["heldout_outcomes_read_at_freeze"] is not False:
        raise RuntimeError("wall-conditioning target firewall is not closed")
    ordered = tuple(measurement_by_key[key] for key in sorted(measurement_by_key))
    return ordered, traces, lot_type_by_key


def _predict(coefficients, measurements, traces, lot_type_by_key):
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (3,) or np.any(~np.isfinite(coefficients)):
        raise ValueError("conditioning coefficients must be one finite triplet")
    law = BoschSPTSWallConditioningLaw(
        log_carbon_cycle_coefficient=float(coefficients[0]),
        silicon_precondition_coefficient=float(coefficients[1]),
        silicon_oxide_precondition_coefficient=float(coefficients[2]),
    )
    first = measurements[0]
    x_m = first.x_um * 1.0e-6
    y_m = first.y_um * 1.0e-6
    if any(
        not np.array_equal(measurement.x_um, first.x_um)
        or not np.array_equal(measurement.y_um, first.y_um)
        for measurement in measurements
    ):
        raise RuntimeError("calibration wafers do not share the official 89-point map")

    lot_types = tuple(sorted(set(lot_type_by_key.values())))
    indexed_measurements_by_type = {
        lot_type: tuple(
            (index, measurement)
            for index, measurement in enumerate(measurements)
            if lot_type_by_key[measurement.experiment_key] == lot_type
        )
        for lot_type in lot_types
    }

    def solve_lot_type(lot_type):
        reduced = conditioned_bosch_spts_parameters(
            BASE_REDUCED_PARAMETERS, law, lot_type)
        model = DeterministicBoschSPTSCylindricalReactorToWafer(
            BoschSPTSCylindricalParameters(
                reduced=reduced, **BASE_CYLINDRICAL_KEYWORDS))
        response = model.source_response(x_m=x_m, y_m=y_m)
        return tuple(
            (index, model.solve(
                traces[measurement.experiment_key], x_m=x_m, y_m=y_m,
                source_response=response))
            for index, measurement in indexed_measurements_by_type[lot_type]
        )

    # Each declared conditioning type owns an independent sparse transport
    # factorization.  Run those exact operators concurrently; no physics or
    # numerical tolerance changes, and results are restored to key order.
    with ThreadPoolExecutor(max_workers=len(lot_types)) as pool:
        grouped_boundaries = tuple(pool.map(solve_lot_type, lot_types))
    boundaries = [None] * len(measurements)
    for group in grouped_boundaries:
        for index, boundary in group:
            boundaries[index] = boundary
    if any(boundary is None for boundary in boundaries):
        raise RuntimeError("conditioning evaluation missed a calibration wafer")
    silicon, oxide = build_bosch_reference_surface_mechanisms()
    predictions = predict_bosch_wafer_point_depth_batch_fast(
        tuple(boundaries), silicon, oxide)
    return predictions, law


def _arrays(measurements, predictions):
    observed_si = np.stack([item.silicon_depth_um for item in measurements])
    observed_oxide = np.stack([item.oxide_loss_um for item in measurements])
    predicted_si = np.stack([
        item.silicon_depth_m * 1.0e6 for item in predictions])
    predicted_oxide = np.stack([
        item.oxide_loss_m * 1.0e6 for item in predictions])
    return observed_si, observed_oxide, predicted_si, predicted_oxide


def _metrics(measurements, predictions):
    observed_si, observed_oxide, predicted_si, predicted_oxide = _arrays(
        measurements, predictions)
    observed_si_mean = np.mean(observed_si, axis=1)
    predicted_si_mean = np.mean(predicted_si, axis=1)
    observed_oxide_mean = np.mean(observed_oxide, axis=1)
    predicted_oxide_mean = np.mean(predicted_oxide, axis=1)
    observed_selectivity = observed_si_mean / observed_oxide_mean
    predicted_selectivity = predicted_si_mean / predicted_oxide_mean
    observed_shape = observed_si / observed_si_mean[:, None]
    predicted_shape = predicted_si / predicted_si_mean[:, None]
    return {
        "silicon_mean_mae_um": float(np.mean(np.abs(
            predicted_si_mean - observed_si_mean))),
        "silicon_mean_mape_percent": float(100.0 * np.mean(np.abs(
            predicted_si_mean / observed_si_mean - 1.0))),
        "silicon_point_rmse_um": float(np.sqrt(np.mean(
            (predicted_si - observed_si) ** 2))),
        "normalized_shape_rmse_percent": float(100.0 * np.sqrt(np.mean(
            (predicted_shape - observed_shape) ** 2))),
        "oxide_mean_mae_um": float(np.mean(np.abs(
            predicted_oxide_mean - observed_oxide_mean))),
        "oxide_mean_mape_percent": float(100.0 * np.mean(np.abs(
            predicted_oxide_mean / observed_oxide_mean - 1.0))),
        "selectivity_mape_percent": float(100.0 * np.mean(np.abs(
            predicted_selectivity / observed_selectivity - 1.0))),
        "silicon_mean_correlation": float(np.corrcoef(
            predicted_si_mean, observed_si_mean)[0, 1]),
        "oxide_mean_correlation": float(np.corrcoef(
            predicted_oxide_mean, observed_oxide_mean)[0, 1]),
        "predicted_silicon_mean_range_um": [
            float(np.min(predicted_si_mean)), float(np.max(predicted_si_mean))],
        "predicted_oxide_mean_range_um": [
            float(np.min(predicted_oxide_mean)), float(np.max(predicted_oxide_mean))],
    }


def _objective(coefficients, measurements, traces, lot_type_by_key):
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
    observed_selectivity = observed_si_mean / observed_oxide_mean
    predicted_selectivity = predicted_si_mean / predicted_oxide_mean
    count = float(len(measurements))
    return np.concatenate((
        (predicted_si_mean - observed_si_mean) / np.sqrt(count),
        ((predicted_oxide_mean - observed_oxide_mean) / 0.08) / np.sqrt(count),
        ((predicted_selectivity / observed_selectivity - 1.0) / 0.12)
        / np.sqrt(count),
        ((predicted_shape - observed_shape) / 0.02).ravel()
        / np.sqrt(count * observed_shape.shape[1]),
    ))


def _leave_one_lot_out_baselines(measurements):
    observed_si = np.stack([item.silicon_depth_um for item in measurements])
    lots = np.asarray([item.lot_number for item in measurements])
    predicted_global = np.zeros_like(observed_si)
    predicted_map = np.zeros_like(observed_si)
    for lot in sorted(set(lots)):
        train = lots != lot
        test = lots == lot
        predicted_global[test] = np.mean(observed_si[train])
        predicted_map[test] = np.mean(observed_si[train], axis=0)
    observed_mean = np.mean(observed_si, axis=1)
    predicted_global_mean = np.mean(predicted_global, axis=1)
    observed_shape = observed_si / observed_mean[:, None]
    predicted_map_shape = predicted_map / np.mean(predicted_map, axis=1)[:, None]
    return {
        "global_mean": {
            "silicon_mean_mae_um": float(np.mean(np.abs(
                predicted_global_mean - observed_mean))),
            "silicon_mean_mape_percent": float(100.0 * np.mean(np.abs(
                predicted_global_mean / observed_mean - 1.0))),
        },
        "mean_map": {
            "silicon_point_rmse_um": float(np.sqrt(np.mean(
                (predicted_map - observed_si) ** 2))),
            "normalized_shape_rmse_percent": float(100.0 * np.sqrt(np.mean(
                (predicted_map_shape - observed_shape) ** 2))),
        },
    }


def build(max_nfev=40):
    measurements, traces, lot_type_by_key = _inputs()
    fit = least_squares(
        _objective, np.zeros(3), bounds=(-1.5, 1.5),
        args=(measurements, traces, lot_type_by_key),
        max_nfev=int(max_nfev), xtol=2.0e-5, ftol=2.0e-5, gtol=2.0e-5,
        verbose=1,
    )
    predictions, law = _predict(
        fit.x, measurements, traces, lot_type_by_key)
    metrics = _metrics(measurements, predictions)
    absolute_pass = {
        name: metrics[name] <= threshold
        for name, threshold in FROZEN_GATES.items()
    }
    baselines = _leave_one_lot_out_baselines(measurements)
    beats_baselines = {
        "global_mean_depth_mae": (
            metrics["silicon_mean_mae_um"]
            < baselines["global_mean"]["silicon_mean_mae_um"]),
        "mean_map_point_rmse": (
            metrics["silicon_point_rmse_um"]
            < baselines["mean_map"]["silicon_point_rmse_um"]),
        "mean_map_normalized_shape_rmse": (
            metrics["normalized_shape_rmse_percent"]
            < baselines["mean_map"]["normalized_shape_rmse_percent"]),
    }
    return {
        "schema": "petch-zenodo-bosch-wall-conditioning-calibration-fit-v1",
        "status": "calibration-only exploratory fit; not a heldout prediction seal",
        "calibration_wafer_count": len(measurements),
        "points_per_wafer": int(measurements[0].x_um.size),
        "conditioning_coefficients": {
            "log_carbon_cycle": float(fit.x[0]),
            "silicon_precondition": float(fit.x[1]),
            "silicon_oxide_precondition": float(fit.x[2]),
        },
        "conditioning_multipliers_by_declared_type": {
            lot_type: law.multiplier(lot_type)
            for lot_type in sorted(set(lot_type_by_key.values()))
        },
        "optimizer": {
            "success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "nfev": int(fit.nfev),
            "cost": float(fit.cost),
            "coefficient_bounds": [-1.5, 1.5],
        },
        "calibration_metrics": metrics,
        "frozen_absolute_gates": FROZEN_GATES,
        "absolute_gate_pass": absolute_pass,
        "all_absolute_gates_pass": all(absolute_pass.values()),
        "leave_one_lot_out_empirical_baselines": baselines,
        "in_sample_physics_beats_empirical_baseline": beats_baselines,
        "eligible_for_prediction_seal": (
            all(absolute_pass.values()) and all(beats_baselines.values())),
        "leave_one_lot_out_physics_refits_completed": False,
        "certification_grid_refinement_completed": False,
        "heldout_prediction_written": False,
        "heldout_outcomes_read": False,
        "surface_laws_changed": False,
        "per_lot_or_per_wafer_depth_offsets": False,
        "input_hashes": {
            "process_data": _hash(PROCESS),
            "process_dictionary": _hash(DICTIONARY),
            "calibration_measurements": _hash(CALIBRATION),
            "v1_preregistration": _hash(V1_PREREGISTRATION),
            "v2_preregistration": _hash(V2_PREREGISTRATION),
            "v5_preregistration": _hash(V5_PREREGISTRATION),
        },
        "implementation_hashes": {
            "reduced_reactor": _hash(REDUCED_IMPLEMENTATION),
            "cylindrical_transfer": _hash(CYLINDRICAL_IMPLEMENTATION),
            "fused_surface_recurrence": _hash(SURFACE_IMPLEMENTATION),
            "calibration_driver": _hash(Path(__file__)),
        },
        "base_reduced_parameters": BASE_REDUCED_PARAMETERS.manifest(),
        "base_cylindrical_parameters": {
            key: value for key, value in BASE_CYLINDRICAL_KEYWORDS.items()
        },
    }


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--max-nfev", type=int, default=40)
    args = parser.parse_args(argv)
    payload = build(max_nfev=args.max_nfev)
    rendered = _render(payload)
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
