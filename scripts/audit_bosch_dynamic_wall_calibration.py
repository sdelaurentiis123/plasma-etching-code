#!/usr/bin/env python3
"""Calibration-only audit of the preregistered Bosch dynamic wall closure.

The expensive deterministic reactor/cylindrical/surface response is tabulated
at fixed log wall-loss nodes.  A shape-preserving interpolant accelerates the
six-parameter search and whole-lot refits.  Interpolation is never a prediction
seal: independent midpoint checks and an exact selected-parameter replay remain
mandatory before heldout outcomes can be opened.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
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
    BoschSPTSDynamicWallLaw,
    BoschSPTSDynamicWallState,
    BoschSPTSWallConditioningLaw,
    advance_bosch_spts_dynamic_wall,
)
from scripts.audit_bosch_wall_conditioning_calibration import (  # noqa: E402
    BASE_CYLINDRICAL_KEYWORDS,
    BASE_REDUCED_PARAMETERS,
    FROZEN_GATES,
    _arrays,
    _leave_one_lot_out_baselines,
    _metrics,
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
V5_FIT = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_wall_conditioning_depth_extension_v5"
    / "calibration_fit.json"
)
V5_RESIDUAL = V5_FIT.parent / "calibration_residual_feature_audit.json"
V6_DIRECTORY = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_dynamic_wall_depth_extension_v6"
)
V6_PREREGISTRATION = V6_DIRECTORY / "preregistration.json"
RESPONSE_TABLE = V6_DIRECTORY / "exact_response_table.npz"
OUTPUT = V6_DIRECTORY / "calibration_fit.json"
REDUCED_IMPLEMENTATION = (
    ROOT / "src" / "petch" / "reactor_global" / "bosch_spts_reduced.py")
CYLINDRICAL_IMPLEMENTATION = (
    ROOT / "src" / "petch" / "reactor_global" / "bosch_spts_cylindrical.py")
SURFACE_IMPLEMENTATION = ROOT / "src" / "petch" / "bosch_wafer_depth_fast.py"

_LOG_FOUR = math.log(4.0)
_PARAMETER_LOWER = np.array([
    -1.5, -1.5, -1.5, math.log(0.001), math.log(0.001), -_LOG_FOUR,
])
_PARAMETER_UPPER = np.array([
    1.5, 1.5, 1.5, math.log(3.0), math.log(3.0), _LOG_FOUR,
])
_FIXED_STARTS = (
    np.array([0.0, 0.0, 0.0, math.log(0.20), math.log(0.05), 0.20]),
    np.array([0.0, 0.0, 0.0, math.log(0.05), math.log(0.20), -0.20]),
    np.array([0.0, 0.0, 0.0, math.log(0.50), math.log(0.50), 0.70]),
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


@dataclass(frozen=True)
class _Prediction:
    silicon_depth_m: np.ndarray
    oxide_loss_m: np.ndarray


@dataclass(frozen=True)
class BoschExactWallResponseTable:
    """Exact node solutions with per-wafer shape-preserving interpolation."""

    experiment_keys: tuple[str, ...]
    log_wall_multiplier_nodes: np.ndarray
    silicon_depth_m: np.ndarray
    oxide_loss_m: np.ndarray

    def __post_init__(self):
        keys = tuple(str(key) for key in self.experiment_keys)
        nodes = np.asarray(self.log_wall_multiplier_nodes, dtype=float).copy()
        silicon = np.asarray(self.silicon_depth_m, dtype=float).copy()
        oxide = np.asarray(self.oxide_loss_m, dtype=float).copy()
        expected = (nodes.size, len(keys), 89)
        if (
            len(keys) == 0
            or len(keys) != len(set(keys))
            or nodes.ndim != 1
            or nodes.size < 9
            or np.any(~np.isfinite(nodes))
            or np.any(np.diff(nodes) <= 0.0)
            or nodes[0] != -_LOG_FOUR
            or nodes[-1] != _LOG_FOUR
            or silicon.shape != expected
            or oxide.shape != expected
            or np.any(~np.isfinite(silicon))
            or np.any(~np.isfinite(oxide))
            or np.any(silicon < 0.0)
            or np.any(oxide < 0.0)
        ):
            raise ValueError("invalid Bosch exact wall-response table")
        for value in (nodes, silicon, oxide):
            value.setflags(write=False)
        object.__setattr__(self, "experiment_keys", keys)
        object.__setattr__(self, "log_wall_multiplier_nodes", nodes)
        object.__setattr__(self, "silicon_depth_m", silicon)
        object.__setattr__(self, "oxide_loss_m", oxide)

    def interpolate(self, wall_multipliers):
        multipliers = np.asarray(wall_multipliers, dtype=float)
        if (
            multipliers.shape != (len(self.experiment_keys),)
            or np.any(~np.isfinite(multipliers))
            or np.any(multipliers < 0.25)
            or np.any(multipliers > 4.0)
        ):
            raise ValueError("invalid Bosch wall multipliers for interpolation")
        query = np.log(multipliers)
        output = []
        for wafer, value in enumerate(query):
            silicon = PchipInterpolator(
                self.log_wall_multiplier_nodes,
                self.silicon_depth_m[:, wafer],
                axis=0,
                extrapolate=False,
            )(value)
            oxide = PchipInterpolator(
                self.log_wall_multiplier_nodes,
                self.oxide_loss_m[:, wafer],
                axis=0,
                extrapolate=False,
            )(value)
            output.append(_Prediction(
                silicon_depth_m=np.asarray(silicon, dtype=float),
                oxide_loss_m=np.asarray(oxide, dtype=float),
            ))
        return tuple(output)

    def save(self, path):
        np.savez_compressed(
            path,
            experiment_keys=np.asarray(self.experiment_keys, dtype="U32"),
            log_wall_multiplier_nodes=self.log_wall_multiplier_nodes,
            silicon_depth_m=self.silicon_depth_m,
            oxide_loss_m=self.oxide_loss_m,
        )

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as payload:
            return cls(
                experiment_keys=tuple(payload["experiment_keys"].tolist()),
                log_wall_multiplier_nodes=payload["log_wall_multiplier_nodes"],
                silicon_depth_m=payload["silicon_depth_m"],
                oxide_loss_m=payload["oxide_loss_m"],
            )


def _inputs():
    v1 = _load_json(V1_PREREGISTRATION)
    v2 = _load_json(V2_PREREGISTRATION)
    v6 = _load_json(V6_PREREGISTRATION)
    if (
        v6["calibration_exposure"]["heldout_outcomes_examined"] is not False
        or v6["target_firewall"]["heldout_outcomes_read_at_freeze"] is not False
    ):
        raise RuntimeError("Bosch dynamic-wall target firewall is not closed")
    allowed = tuple(v1["split_rule"]["calibration_experiment_keys"])
    measurements = load_bosch_wafer_measurements_89pt(
        CALIBRATION, allowed_experiment_keys=allowed)
    measurement_keys = {item.experiment_key for item in measurements}
    process_traces = tuple(
        trace for trace in load_bosch_process_traces(PROCESS, DICTIONARY)
        if trace.experiment_key in allowed
    )
    trace_keys = {trace.experiment_key for trace in process_traces}
    if len(process_traces) != len(allowed) or trace_keys != set(allowed):
        raise RuntimeError("calibration process allowlist is incomplete")
    if not measurement_keys.issubset(trace_keys):
        raise RuntimeError("a calibration measurement lacks a process trace")
    ordered_measurements = tuple(sorted(
        measurements, key=lambda item: item.experiment_key))
    ordered_traces = tuple(sorted(
        process_traces,
        key=lambda item: (item.process_date, item.wafer_number),
    ))
    lot_type_by_date = v2["conditioning_state"]["declared_lot_types"]
    if set(trace.process_date for trace in ordered_traces) - set(lot_type_by_date):
        raise RuntimeError("a calibration lot lacks declared conditioning metadata")
    return ordered_measurements, ordered_traces, lot_type_by_date


def _wall_nodes(count=9):
    count = int(count)
    if count < 9 or count % 2 == 0:
        raise ValueError("Bosch wall response requires an odd node count >= 9")
    unit = -np.cos(np.linspace(0.0, math.pi, count))
    nodes = _LOG_FOUR * unit
    nodes[0] = -_LOG_FOUR
    nodes[-1] = _LOG_FOUR
    return nodes


def _exact_predictions_at_multiplier(
        multiplier, measurements, trace_by_key, *, workers=8, batch_size=15):
    multiplier = float(multiplier)
    reduced = replace(
        BASE_REDUCED_PARAMETERS,
        neutral_wall_loss_multiplier=multiplier,
    )
    model = DeterministicBoschSPTSCylindricalReactorToWafer(
        BoschSPTSCylindricalParameters(
            reduced=reduced, **BASE_CYLINDRICAL_KEYWORDS))
    first = measurements[0]
    x_m = first.x_um * 1.0e-6
    y_m = first.y_um * 1.0e-6
    if any(
        not np.array_equal(item.x_um, first.x_um)
        or not np.array_equal(item.y_um, first.y_um)
        for item in measurements
    ):
        raise RuntimeError("calibration wafers do not share the official 89-point map")
    response = model.source_response(x_m=x_m, y_m=y_m)
    silicon, oxide = build_bosch_reference_surface_mechanisms()
    outputs = []
    with ThreadPoolExecutor(max_workers=int(workers)) as pool:
        for start in range(0, len(measurements), int(batch_size)):
            group = measurements[start:start + int(batch_size)]

            def solve(item):
                return model.solve(
                    trace_by_key[item.experiment_key],
                    x_m=x_m,
                    y_m=y_m,
                    source_response=response,
                )

            boundaries = tuple(pool.map(solve, group))
            outputs.extend(predict_bosch_wafer_point_depth_batch_fast(
                boundaries, silicon, oxide))
    return tuple(outputs)


def _build_response_table(measurements, process_traces, *, workers=8,
                          node_count=9):
    trace_by_key = {trace.experiment_key: trace for trace in process_traces}
    nodes = _wall_nodes(node_count)
    silicon = []
    oxide = []
    for index, log_multiplier in enumerate(nodes, start=1):
        multiplier = math.exp(float(log_multiplier))
        print(
            f"exact wall-response node {index}/{len(nodes)} "
            f"multiplier={multiplier:.9g}",
            flush=True,
        )
        predictions = _exact_predictions_at_multiplier(
            multiplier, measurements, trace_by_key, workers=workers)
        silicon.append(np.stack([item.silicon_depth_m for item in predictions]))
        oxide.append(np.stack([item.oxide_loss_m for item in predictions]))
    return BoschExactWallResponseTable(
        experiment_keys=tuple(item.experiment_key for item in measurements),
        log_wall_multiplier_nodes=nodes,
        silicon_depth_m=np.stack(silicon),
        oxide_loss_m=np.stack(oxide),
    )


def _load_or_build_response_table(measurements, process_traces, *, workers=8,
                                  node_count=9, write=False):
    expected_keys = tuple(item.experiment_key for item in measurements)
    expected_nodes = _wall_nodes(node_count)
    if RESPONSE_TABLE.exists():
        table = BoschExactWallResponseTable.load(RESPONSE_TABLE)
        if (
            table.experiment_keys != expected_keys
            or not np.array_equal(table.log_wall_multiplier_nodes, expected_nodes)
        ):
            raise RuntimeError("cached Bosch wall-response table is stale")
        return table
    table = _build_response_table(
        measurements, process_traces, workers=workers, node_count=node_count)
    if write:
        table.save(RESPONSE_TABLE)
    return table


def _laws(coefficients):
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (6,) or np.any(~np.isfinite(coefficients)):
        raise ValueError("Bosch dynamic wall coefficients must be one finite sextet")
    static = BoschSPTSWallConditioningLaw(
        log_carbon_cycle_coefficient=float(coefficients[0]),
        silicon_precondition_coefficient=float(coefficients[1]),
        silicon_oxide_precondition_coefficient=float(coefficients[2]),
    )
    dynamic = BoschSPTSDynamicWallLaw(
        deposition_rate_per_reference_wafer=math.exp(float(coefficients[3])),
        cleaning_rate_per_reference_wafer=math.exp(float(coefficients[4])),
        log_wall_loss_response=float(coefficients[5]),
    )
    return static, dynamic


def _dynamic_steps(coefficients, process_traces, lot_type_by_date):
    static, dynamic = _laws(coefficients)
    state_by_date = {}
    steps = {}
    for trace in process_traces:
        state = state_by_date.setdefault(
            trace.process_date, BoschSPTSDynamicWallState())
        static_multiplier = static.multiplier(lot_type_by_date[trace.process_date])
        step = advance_bosch_spts_dynamic_wall(
            trace, dynamic, state,
            static_wall_loss_multiplier=static_multiplier,
        )
        steps[trace.experiment_key] = step
        state_by_date[trace.process_date] = step.end_state
    if len(steps) != len(process_traces):
        raise RuntimeError("dynamic Bosch wall history lost a process trace")
    return steps, static, dynamic


def _predict(coefficients, measurements, process_traces, lot_type_by_date, table):
    steps, static, dynamic = _dynamic_steps(
        coefficients, process_traces, lot_type_by_date)
    multipliers = np.asarray([
        steps[item.experiment_key].combined_wall_loss_multiplier
        for item in measurements
    ])
    predictions = table.interpolate(multipliers)
    return predictions, steps, multipliers, static, dynamic


def _objective_from_arrays(measurements, predictions):
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


def _objective(coefficients, selected_indices, measurements, process_traces,
               lot_type_by_date, table):
    predictions, *_ = _predict(
        coefficients, measurements, process_traces, lot_type_by_date, table)
    subset_measurements = tuple(measurements[index] for index in selected_indices)
    subset_predictions = tuple(predictions[index] for index in selected_indices)
    return _objective_from_arrays(subset_measurements, subset_predictions)


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
    _cost, start_index, result = min(attempts, key=lambda item: (item[0], item[1]))
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
    static, dynamic = _laws(coefficients)
    return {
        "conditioning_repeat_log_coefficient": float(coefficients[0]),
        "silicon_precondition_coefficient": float(coefficients[1]),
        "silicon_oxide_precondition_coefficient": float(coefficients[2]),
        "deposition_rate_per_reference_wafer": (
            dynamic.deposition_rate_per_reference_wafer),
        "cleaning_rate_per_reference_wafer": (
            dynamic.cleaning_rate_per_reference_wafer),
        "log_wall_loss_response": dynamic.log_wall_loss_response,
        "static_law_manifest": static.manifest(),
        "dynamic_law_manifest": dynamic.manifest(),
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
    dynamic_bound_contact = bool(np.any(normalized_bound_distance[3:] <= 0.001))
    gates = {
        "full_column_rank": rank == 6,
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


def _lot_slopes(measurements, predictions):
    rows = {}
    for lot in sorted({item.lot_number for item in measurements}):
        indices = [
            index for index, item in enumerate(measurements)
            if item.lot_number == lot
        ]
        wafer = np.asarray([measurements[index].wafer_number for index in indices])
        observed = np.asarray([
            np.mean(measurements[index].silicon_depth_um) for index in indices])
        predicted = np.asarray([
            np.mean(predictions[index].silicon_depth_m) * 1.0e6
            for index in indices
        ])
        rows[str(lot)] = {
            "observed_um_per_wafer": float(np.polyfit(wafer, observed, 1)[0]),
            "predicted_um_per_wafer": float(np.polyfit(wafer, predicted, 1)[0]),
        }
        rows[str(lot)]["residual_um_per_wafer"] = (
            rows[str(lot)]["predicted_um_per_wafer"]
            - rows[str(lot)]["observed_um_per_wafer"])
    return rows


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


def _interpolation_validation(table, measurements, process_traces, *, workers):
    trace_by_key = {trace.experiment_key: trace for trace in process_traces}
    midpoint_nodes = 0.5 * (
        table.log_wall_multiplier_nodes[:-1]
        + table.log_wall_multiplier_nodes[1:]
    )
    node_rows = []
    maximum_fraction = 0.0
    for index, log_multiplier in enumerate(midpoint_nodes, start=1):
        multiplier = math.exp(float(log_multiplier))
        print(
            f"interpolation validation {index}/{len(midpoint_nodes)} "
            f"multiplier={multiplier:.9g}",
            flush=True,
        )
        exact = _exact_predictions_at_multiplier(
            multiplier, measurements, trace_by_key, workers=workers)
        interpolated = table.interpolate(np.full(len(measurements), multiplier))
        exact_si = np.stack([item.silicon_depth_m for item in exact]) * 1.0e6
        interp_si = np.stack([item.silicon_depth_m for item in interpolated]) * 1.0e6
        exact_oxide = np.stack([item.oxide_loss_m for item in exact]) * 1.0e6
        interp_oxide = np.stack([item.oxide_loss_m for item in interpolated]) * 1.0e6
        exact_si_mean = np.mean(exact_si, axis=1)
        interp_si_mean = np.mean(interp_si, axis=1)
        exact_oxide_mean = np.mean(exact_oxide, axis=1)
        interp_oxide_mean = np.mean(interp_oxide, axis=1)
        exact_shape = exact_si / exact_si_mean[:, None]
        interp_shape = interp_si / interp_si_mean[:, None]
        exact_selectivity = exact_si_mean / exact_oxide_mean
        interp_selectivity = interp_si_mean / interp_oxide_mean
        errors = {
            "silicon_mean_mae_um": float(np.mean(np.abs(
                interp_si_mean - exact_si_mean))),
            "silicon_point_rmse_um": float(np.sqrt(np.mean(
                (interp_si - exact_si) ** 2))),
            "normalized_shape_rmse_percent": float(100.0 * np.sqrt(np.mean(
                (interp_shape - exact_shape) ** 2))),
            "oxide_mean_mae_um": float(np.mean(np.abs(
                interp_oxide_mean - exact_oxide_mean))),
            "selectivity_mape_percent": float(100.0 * np.mean(np.abs(
                interp_selectivity / exact_selectivity - 1.0))),
        }
        fractions = {
            name: value / FROZEN_GATES[name]
            for name, value in errors.items()
        }
        maximum_fraction = max(maximum_fraction, max(fractions.values()))
        node_rows.append({
            "wall_loss_multiplier": multiplier,
            "errors": errors,
            "fractions_of_frozen_gates": fractions,
        })
    return {
        "independent_midpoint_node_count": len(midpoint_nodes),
        "nodes": node_rows,
        "maximum_error_fraction_of_any_frozen_gate": maximum_fraction,
        "frozen_maximum": 0.05,
        "passes": maximum_fraction <= 0.05,
    }


def build(*, workers=8, max_nfev=120, node_count=9,
          validate_interpolation=True, write_table=False):
    measurements, process_traces, lot_type_by_date = _inputs()
    table = _load_or_build_response_table(
        measurements, process_traces, workers=workers,
        node_count=node_count, write=write_table)
    all_indices = tuple(range(len(measurements)))
    fit, start_index, attempts = _fit(
        all_indices, measurements, process_traces, lot_type_by_date, table,
        max_nfev=max_nfev)
    predictions, steps, multipliers, static, dynamic = _predict(
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
        else {
            "passes": False,
            "status": "not run",
            "frozen_maximum": 0.05,
        }
    )
    lot_slopes = _lot_slopes(measurements, predictions)
    v5_slopes = _load_json(V5_RESIDUAL)[
        "per_lot_sequence_slopes_um_per_wafer"]
    v5_slope_mae = float(np.mean([
        abs(float(v5_slopes["predicted_silicon_mean"][lot])
            - float(v5_slopes["observed_silicon_mean"][lot]))
        for lot in sorted(v5_slopes["observed_silicon_mean"])
    ]))
    dynamic_slope_mae = float(np.mean([
        abs(row["residual_um_per_wafer"]) for row in lot_slopes.values()
    ]))
    identification = _identifiability(fit)
    response_hash = _hash(RESPONSE_TABLE) if RESPONSE_TABLE.exists() else None
    seal_prerequisites = {
        "all_absolute_calibration_gates_pass": all(absolute_pass.values()),
        "whole_lot_physics_beats_all_empirical_baselines": all(
            beats_baselines.values()),
        "within_lot_slope_error_improves_over_v5": (
            dynamic_slope_mae < v5_slope_mae),
        "identifiability_gate_passes": identification["passes"],
        "interpolation_validation_passes": interpolation["passes"],
        "exact_selected_parameter_replay_completed": False,
        "certification_grid_refinement_completed": False,
        "heldout_prediction_hash_sealed": False,
    }
    return {
        "schema": "petch-zenodo-bosch-dynamic-wall-calibration-fit-v1",
        "status": "calibration-only surrogate search; exact replay and heldout seal remain closed",
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
            "dynamic_lot_slopes": lot_slopes,
            "v5_static_slope_mae_um_per_wafer": v5_slope_mae,
            "dynamic_slope_mae_um_per_wafer": dynamic_slope_mae,
            "improves_over_v5": dynamic_slope_mae < v5_slope_mae,
        },
        "wall_history": {
            "minimum_applied_multiplier": float(np.min(multipliers)),
            "maximum_applied_multiplier": float(np.max(multipliers)),
            "minimum_mean_occupancy": float(min(
                steps[item.experiment_key].mean_occupancy
                for item in measurements)),
            "maximum_mean_occupancy": float(max(
                steps[item.experiment_key].mean_occupancy
                for item in measurements)),
            "static_law": static.manifest(),
            "dynamic_law": dynamic.manifest(),
        },
        "deterministic_acceleration": {
            "method": "PCHIP in log wall-loss multiplier over exact reactor-to-surface nodes",
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
        "input_hashes": {
            "process_data": _hash(PROCESS),
            "process_dictionary": _hash(DICTIONARY),
            "calibration_measurements": _hash(CALIBRATION),
            "v1_preregistration": _hash(V1_PREREGISTRATION),
            "v2_preregistration": _hash(V2_PREREGISTRATION),
            "v5_fit": _hash(V5_FIT),
            "v5_residual_audit": _hash(V5_RESIDUAL),
            "v6_preregistration": _hash(V6_PREREGISTRATION),
        },
        "implementation_hashes": {
            "reduced_reactor": _hash(REDUCED_IMPLEMENTATION),
            "cylindrical_transfer": _hash(CYLINDRICAL_IMPLEMENTATION),
            "surface_recurrence": _hash(SURFACE_IMPLEMENTATION),
            "dynamic_calibration_driver": _hash(Path(__file__)),
        },
    }


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--max-nfev", type=int, default=120)
    parser.add_argument("--node-count", type=int, default=9)
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
