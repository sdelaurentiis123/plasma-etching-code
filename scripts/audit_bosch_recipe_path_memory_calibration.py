#!/usr/bin/env python3
"""Calibration-only audit of the preregistered Bosch recipe-path memory.

The expensive deterministic reactor/cylindrical/surface response is tabulated
at thirteen frozen log wall-loss nodes.  A PCHIP interpolant accelerates the
four-parameter search and whole-lot refits.  The chronological heldout outcome
asset is never loaded by this module.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.reactor_global.bosch_spts_reduced import (  # noqa: E402
    BoschSPTSRecipePathWallLaw,
    BoschSPTSRecipePathWallState,
    BoschSPTSWallConditioningLaw,
    advance_bosch_spts_recipe_path_wall,
)
from scripts.audit_bosch_dynamic_wall_calibration import (  # noqa: E402
    BoschExactWallResponseTable,
    FROZEN_GATES,
    V5_RESIDUAL,
    _build_response_table,
    _inputs as _v6_inputs,
    _interpolation_validation,
    _load_json,
    _lot_slopes,
    _objective_from_arrays,
    _wall_nodes,
)
from scripts.audit_bosch_wall_conditioning_calibration import (  # noqa: E402
    _leave_one_lot_out_baselines,
    _metrics,
)


V7_DIRECTORY = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_recipe_path_memory_depth_extension_v7"
)
V7_PREREGISTRATION = V7_DIRECTORY / "preregistration.json"
RESPONSE_TABLE = V7_DIRECTORY / "exact_response_table.npz"
OUTPUT = V7_DIRECTORY / "calibration_fit.json"
V6_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_dynamic_wall_depth_extension_v6"
    / "preregistration.json"
)
V6_FIT = V6_PREREGISTRATION.parent / "calibration_fit.json"
V6_RESPONSE_TABLE = V6_PREREGISTRATION.parent / "exact_response_table.npz"
REDUCED_IMPLEMENTATION = (
    ROOT / "src" / "petch" / "reactor_global" / "bosch_spts_reduced.py")
CYLINDRICAL_IMPLEMENTATION = (
    ROOT / "src" / "petch" / "reactor_global"
    / "bosch_spts_cylindrical.py")
SURFACE_IMPLEMENTATION = ROOT / "src" / "petch" / "bosch_wafer_depth_fast.py"

_LOG_FOUR = math.log(4.0)
_DYNAMIC_BOUND = _LOG_FOUR / 10.0
_PARAMETER_LOWER = np.array([-1.5, -1.5, -1.5, -_DYNAMIC_BOUND])
_PARAMETER_UPPER = np.array([1.5, 1.5, 1.5, _DYNAMIC_BOUND])
_FIXED_STARTS = (
    np.array([0.0, 0.0, 0.0, 0.0]),
    np.array([0.0, 0.0, 0.0, -0.05]),
    np.array([0.0, 0.0, 0.0, 0.08]),
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _inputs():
    preregistration = _load_json(V7_PREREGISTRATION)
    exposure = preregistration["calibration_exposure"]
    firewall = preregistration["target_firewall"]
    if (
        exposure["heldout_outcomes_examined"] is not False
        or firewall["heldout_outcomes_read_at_freeze"] is not False
        or exposure["v7_operator_implemented_before_this_freeze"] is not False
        or exposure["v7_fit_started_before_this_freeze"] is not False
    ):
        raise RuntimeError("Bosch recipe-path target firewall is not closed")
    return _v6_inputs()


def _load_or_build_response_table(measurements, process_traces, *, workers=8,
                                  node_count=13, write=False):
    expected_keys = tuple(item.experiment_key for item in measurements)
    expected_nodes = _wall_nodes(node_count)
    if RESPONSE_TABLE.exists():
        table = BoschExactWallResponseTable.load(RESPONSE_TABLE)
        if (
            table.experiment_keys != expected_keys
            or not np.array_equal(table.log_wall_multiplier_nodes, expected_nodes)
        ):
            raise RuntimeError("cached Bosch v7 wall-response table is stale")
        return table
    table = _build_response_table(
        measurements, process_traces, workers=workers, node_count=node_count)
    if write:
        table.save(RESPONSE_TABLE)
    return table


def _laws(coefficients):
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (4,) or np.any(~np.isfinite(coefficients)):
        raise ValueError("Bosch recipe-path coefficients must be one finite quartet")
    static = BoschSPTSWallConditioningLaw(
        log_carbon_cycle_coefficient=float(coefficients[0]),
        silicon_precondition_coefficient=float(coefficients[1]),
        silicon_oxide_precondition_coefficient=float(coefficients[2]),
    )
    path = BoschSPTSRecipePathWallLaw(
        log_wall_loss_per_reference_wafer=float(coefficients[3]))
    return static, path


def _recipe_path_steps(coefficients, process_traces, lot_type_by_date):
    static, path = _laws(coefficients)
    state_by_date = {}
    steps = {}
    for trace in process_traces:
        state = state_by_date.setdefault(
            trace.process_date, BoschSPTSRecipePathWallState())
        static_multiplier = static.multiplier(lot_type_by_date[trace.process_date])
        step = advance_bosch_spts_recipe_path_wall(
            trace, path, state,
            static_wall_loss_multiplier=static_multiplier,
        )
        steps[trace.experiment_key] = step
        state_by_date[trace.process_date] = step.end_state
    if len(steps) != len(process_traces):
        raise RuntimeError("recipe-path Bosch wall history lost a process trace")
    return steps, static, path


def _predict(coefficients, measurements, process_traces, lot_type_by_date, table):
    steps, static, path = _recipe_path_steps(
        coefficients, process_traces, lot_type_by_date)
    multipliers = np.asarray([
        steps[item.experiment_key].combined_wall_loss_multiplier
        for item in measurements
    ])
    return table.interpolate(multipliers), steps, multipliers, static, path


def _objective(coefficients, selected_indices, measurements, process_traces,
               lot_type_by_date, table):
    predictions, *_ = _predict(
        coefficients, measurements, process_traces, lot_type_by_date, table)
    return _objective_from_arrays(
        tuple(measurements[index] for index in selected_indices),
        tuple(predictions[index] for index in selected_indices),
    )


def _fit(selected_indices, measurements, process_traces, lot_type_by_date,
         table, *, max_nfev=120):
    attempts = []
    for start_index, start in enumerate(_FIXED_STARTS):
        result = least_squares(
            _objective,
            start,
            bounds=(_PARAMETER_LOWER, _PARAMETER_UPPER),
            args=(selected_indices, measurements, process_traces,
                  lot_type_by_date, table),
            max_nfev=int(max_nfev),
            xtol=1.0e-8,
            ftol=1.0e-8,
            gtol=1.0e-8,
            verbose=0,
        )
        attempts.append((float(result.cost), start_index, result))
    _cost, start_index, result = min(
        attempts, key=lambda item: (item[0], item[1]))
    return result, start_index, tuple(
        {
            "start_index": index,
            "cost": float(item.cost),
            "success": bool(item.success),
            "nfev": int(item.nfev),
        }
        for _cost, index, item in attempts
    )


def _parameter_manifest(coefficients):
    static, path = _laws(coefficients)
    return {
        "conditioning_repeat_log_coefficient": float(coefficients[0]),
        "silicon_precondition_coefficient": float(coefficients[1]),
        "silicon_oxide_precondition_coefficient": float(coefficients[2]),
        "log_wall_loss_per_reference_wafer": float(coefficients[3]),
        "static_law_manifest": static.manifest(),
        "recipe_path_law_manifest": path.manifest(),
    }


def _identifiability(result):
    jacobian = np.asarray(result.jac, dtype=float)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    tolerance = np.finfo(float).eps * max(jacobian.shape) * singular[0]
    rank = int(np.count_nonzero(singular > tolerance))
    condition = (
        float(singular[0] / singular[-1])
        if singular[-1] > 0.0 else math.inf
    )
    covariance_shape = np.linalg.pinv(jacobian.T @ jacobian, rcond=1.0e-12)
    standard = np.sqrt(np.maximum(np.diag(covariance_shape), 0.0))
    denominator = np.outer(standard, standard)
    correlation = np.divide(
        covariance_shape,
        denominator,
        out=np.zeros_like(covariance_shape),
        where=denominator > 0.0,
    )
    np.fill_diagonal(correlation, 0.0)
    maximum_correlation = float(np.max(np.abs(correlation)))
    span = _PARAMETER_UPPER - _PARAMETER_LOWER
    normalized_bound_distance = np.minimum(
        (result.x - _PARAMETER_LOWER) / span,
        (_PARAMETER_UPPER - result.x) / span,
    )
    dynamic_bound_contact = bool(normalized_bound_distance[3] <= 0.001)
    gates = {
        "full_column_rank": rank == 4,
        "condition_number_at_most_1e6": condition <= 1.0e6,
        "maximum_parameter_correlation_below_0p995": (
            maximum_correlation < 0.995),
        "no_unresolved_dynamic_bound_contact": not dynamic_bound_contact,
    }
    return {
        "jacobian_shape": list(jacobian.shape),
        "rank": rank,
        "singular_values": singular.tolist(),
        "condition_number": condition,
        "maximum_pairwise_parameter_correlation_magnitude": maximum_correlation,
        "normalized_distance_to_nearest_bound": normalized_bound_distance.tolist(),
        "dynamic_bound_contact": dynamic_bound_contact,
        "gates": gates,
        "passes": all(gates.values()),
    }


def _whole_lot_leave_one_out(measurements, process_traces, lot_type_by_date,
                             table, *, max_nfev):
    lots = sorted({item.lot_number for item in measurements})
    predictions = [None] * len(measurements)
    folds = []
    for lot in lots:
        train = tuple(
            index for index, item in enumerate(measurements)
            if item.lot_number != lot)
        test = tuple(
            index for index, item in enumerate(measurements)
            if item.lot_number == lot)
        fit, start_index, attempts = _fit(
            train, measurements, process_traces, lot_type_by_date, table,
            max_nfev=max_nfev)
        all_predictions, *_ = _predict(
            fit.x, measurements, process_traces, lot_type_by_date, table)
        for index in test:
            predictions[index] = all_predictions[index]
        folds.append({
            "heldout_lot": lot,
            "training_wafer_count": len(train),
            "test_wafer_count": len(test),
            "selected_start_index": start_index,
            "optimizer_success": bool(fit.success),
            "optimizer_nfev": int(fit.nfev),
            "optimizer_cost": float(fit.cost),
            "parameters": _parameter_manifest(fit.x),
            "start_attempts": attempts,
            "test_metrics": _metrics(
                tuple(measurements[index] for index in test),
                tuple(all_predictions[index] for index in test)),
        })
    if any(item is None for item in predictions):
        raise RuntimeError("whole-lot validation missed a calibration wafer")
    return tuple(predictions), folds


def build(*, workers=8, max_nfev=120, node_count=13,
          validate_interpolation=True, write_table=False):
    measurements, process_traces, lot_type_by_date = _inputs()
    if int(node_count) != 13:
        raise ValueError("v7 preregistration freezes exactly thirteen nodes")
    table = _load_or_build_response_table(
        measurements, process_traces, workers=workers,
        node_count=node_count, write=write_table)
    fit, start_index, attempts = _fit(
        tuple(range(len(measurements))), measurements, process_traces,
        lot_type_by_date, table, max_nfev=max_nfev)
    predictions, steps, multipliers, static, path = _predict(
        fit.x, measurements, process_traces, lot_type_by_date, table)
    calibration_metrics = _metrics(measurements, predictions)
    absolute_pass = {
        name: calibration_metrics[name] <= threshold
        for name, threshold in FROZEN_GATES.items()
    }
    loo_predictions, folds = _whole_lot_leave_one_out(
        measurements, process_traces, lot_type_by_date, table,
        max_nfev=max_nfev)
    loo_metrics = _metrics(measurements, loo_predictions)
    baselines = _leave_one_lot_out_baselines(measurements)
    beats_baselines = {
        "global_mean_depth_mae": (
            loo_metrics["silicon_mean_mae_um"]
            < baselines["global_mean"]["silicon_mean_mae_um"]),
        "mean_map_point_rmse": (
            loo_metrics["silicon_point_rmse_um"]
            < baselines["mean_map"]["silicon_point_rmse_um"]),
        "mean_map_normalized_shape_rmse": (
            loo_metrics["normalized_shape_rmse_percent"]
            < baselines["mean_map"]["normalized_shape_rmse_percent"]),
    }
    interpolation = (
        _interpolation_validation(
            table, measurements, process_traces, workers=workers)
        if validate_interpolation
        else {"passes": False, "status": "not run", "frozen_maximum": 0.05}
    )
    lot_slopes = _lot_slopes(measurements, predictions)
    v5_slopes = _load_json(V5_RESIDUAL)[
        "per_lot_sequence_slopes_um_per_wafer"]
    v5_slope_mae = float(np.mean([
        abs(float(v5_slopes["predicted_silicon_mean"][lot])
            - float(v5_slopes["observed_silicon_mean"][lot]))
        for lot in sorted(v5_slopes["observed_silicon_mean"])
    ]))
    path_slope_mae = float(np.mean([
        abs(row["residual_um_per_wafer"]) for row in lot_slopes.values()
    ]))
    identification = _identifiability(fit)
    response_hash = _hash(RESPONSE_TABLE) if RESPONSE_TABLE.exists() else None
    seal_prerequisites = {
        "all_absolute_calibration_gates_pass": all(absolute_pass.values()),
        "whole_lot_physics_beats_all_empirical_baselines": all(
            beats_baselines.values()),
        "within_lot_slope_error_improves_over_v5": (
            path_slope_mae < v5_slope_mae),
        "identifiability_gate_passes": identification["passes"],
        "interpolation_validation_passes": interpolation["passes"],
        "exact_selected_parameter_replay_completed": False,
        "certification_grid_refinement_completed": False,
        "heldout_prediction_hash_sealed": False,
    }
    measured_steps = [steps[item.experiment_key] for item in measurements]
    return {
        "schema": "petch-zenodo-bosch-recipe-path-memory-calibration-fit-v1",
        "status": "calibration-only surrogate search; heldout firewall remains closed",
        "calibration_wafer_count": len(measurements),
        "calibration_process_trace_count": len(process_traces),
        "unmeasured_process_traces_carried_through_wall_history": sorted(
            set(trace.experiment_key for trace in process_traces)
            - set(item.experiment_key for item in measurements)),
        "points_per_wafer": 89,
        "parameters": _parameter_manifest(fit.x),
        "optimizer": {
            "success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "nfev": int(fit.nfev),
            "cost": float(fit.cost),
            "selected_fixed_start_index": start_index,
            "fixed_start_attempts": attempts,
            "parameter_lower_bounds": _PARAMETER_LOWER.tolist(),
            "parameter_upper_bounds": _PARAMETER_UPPER.tolist(),
        },
        "identifiability": identification,
        "calibration_metrics": calibration_metrics,
        "frozen_absolute_gates": FROZEN_GATES,
        "absolute_gate_pass": absolute_pass,
        "all_absolute_gates_pass": all(absolute_pass.values()),
        "whole_lot_leave_one_out": {
            "completed": True,
            "folds": folds,
            "aggregate_metrics": loo_metrics,
        },
        "leave_one_lot_out_empirical_baselines": baselines,
        "whole_lot_physics_beats_empirical_baseline": beats_baselines,
        "within_lot_sequence": {
            "recipe_path_lot_slopes": lot_slopes,
            "v5_static_slope_mae_um_per_wafer": v5_slope_mae,
            "recipe_path_slope_mae_um_per_wafer": path_slope_mae,
            "improves_over_v5": path_slope_mae < v5_slope_mae,
        },
        "wall_history": {
            "minimum_applied_multiplier": float(np.min(multipliers)),
            "maximum_applied_multiplier": float(np.max(multipliers)),
            "minimum_mean_reference_wafer_exposure": float(min(
                step.mean_reference_wafer_exposure for step in measured_steps)),
            "maximum_mean_reference_wafer_exposure": float(max(
                step.mean_reference_wafer_exposure for step in measured_steps)),
            "minimum_sf6_to_c4f8_dose_ratio": float(min(
                step.sf6_to_c4f8_dose_ratio for step in measured_steps)),
            "maximum_sf6_to_c4f8_dose_ratio": float(max(
                step.sf6_to_c4f8_dose_ratio for step in measured_steps)),
            "static_law": static.manifest(),
            "recipe_path_law": path.manifest(),
        },
        "deterministic_acceleration": {
            "method": "PCHIP in log wall-loss multiplier over thirteen exact reactor-to-surface nodes",
            "node_count": int(table.log_wall_multiplier_nodes.size),
            "log_wall_multiplier_nodes": table.log_wall_multiplier_nodes.tolist(),
            "response_table": str(RESPONSE_TABLE.relative_to(ROOT)),
            "response_table_sha256": response_hash,
            "validation": interpolation,
            "used_only_for_calibration_search": True,
            "exact_selected_parameter_replay_completed": False,
        },
        "seal_prerequisites": seal_prerequisites,
        "eligible_for_prediction_seal": all(seal_prerequisites.values()),
        "heldout_prediction_written": False,
        "heldout_outcomes_read": False,
        "surface_laws_changed": False,
        "positive_ion_path_changed": False,
        "per_lot_or_per_wafer_depth_offsets": False,
        "wafer_number_used_as_model_input": False,
        "recipe_path_extrapolation_used": False,
        "input_hashes": {
            "v6_preregistration": _hash(V6_PREREGISTRATION),
            "v6_calibration_fit": _hash(V6_FIT),
            "v6_response_table": _hash(V6_RESPONSE_TABLE),
            "v7_preregistration": _hash(V7_PREREGISTRATION),
        },
        "implementation_hashes": {
            "reduced_reactor": _hash(REDUCED_IMPLEMENTATION),
            "cylindrical_transfer": _hash(CYLINDRICAL_IMPLEMENTATION),
            "surface_recurrence": _hash(SURFACE_IMPLEMENTATION),
            "recipe_path_calibration_driver": _hash(Path(__file__)),
        },
    }


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--max-nfev", type=int, default=120)
    parser.add_argument("--node-count", type=int, default=13)
    parser.add_argument("--skip-interpolation-validation", action="store_true")
    args = parser.parse_args(argv)
    payload = build(
        workers=args.workers,
        max_nfev=args.max_nfev,
        node_count=args.node_count,
        validate_interpolation=not args.skip_interpolation_validation,
        write_table=args.write,
    )
    rendered = _render(payload)
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
