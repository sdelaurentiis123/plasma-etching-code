"""Preregistered, intended-use validation gates for product-facing demonstrations.

This module scores predictions; it does not calibrate chemistry.  Keeping those operations separate
prevents held-out observations from silently becoming model inputs and keeps numerical, parameter-
evidence, experimental-digitization, and model-form failures distinguishable.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log, sqrt

from .experimental_data import Jeon2022DimensionlessTarget


JEON_2022_DEMO_VERSION = "jeon_2022_depth_transfer_v1"


@dataclass(frozen=True)
class Jeon2022Prediction:
    observable: str
    source_figure: str
    condition_family: str
    c4f8_fraction: float
    pulse_off_ms: float
    trench_width_nm: float
    value: float
    numerical_log_standard_error: float
    within_declared_scope: bool
    parameter_evidence_supports_prediction: bool


@dataclass(frozen=True)
class Jeon2022DemoThresholds:
    """Thresholds frozen before the first unified-engine Jeon prediction run."""

    calibration_log_rmse: float = 0.10
    held_out_width_log_rmse: float = 0.15
    held_out_pulse_log_rmse: float = 0.15
    held_out_digitization_interval_coverage: float = 0.50
    numerical_log_standard_error: float = 0.03
    maximum_calibrated_parameters: int = 3
    maximum_wall_time_s: float = 7200.0
    maximum_accelerated_compute_cost_usd: float = 2.0


@dataclass(frozen=True)
class Jeon2022DemoScore:
    version: str
    calibration_log_rmse: float
    held_out_width_log_rmse: float
    held_out_pulse_log_rmse: float
    held_out_digitization_interval_coverage: float
    pulse_reversal_correct: int
    pulse_reversal_total: int
    maximum_numerical_log_standard_error: float
    calibrated_parameter_count: int
    wall_time_s: float
    accelerated_compute_cost_usd: float
    passed: bool
    failures: tuple[str, ...]


def _target_key(item):
    return (
        item.observable,
        item.source_figure,
        item.condition_family,
        float(item.c4f8_fraction),
        float(item.pulse_off_ms),
        float(item.trench_width_nm),
    )


def _log_rmse(pairs):
    errors = [(log(prediction.value) - log(target.value)) ** 2
              for target, prediction in pairs]
    return sqrt(sum(errors) / len(errors))


def score_jeon_2022_demo(
        targets, predictions, *, calibrated_parameter_count, wall_time_s,
        accelerated_compute_cost_usd=0.0, thresholds=Jeon2022DemoThresholds()):
    """Score one complete prediction set without changing model or experimental state.

    Width ratios whose numerator is the 200 nm denominator are identities and are required for
    coverage but excluded from RMSE.  Published error-bar semantics are unknown, so interval
    coverage uses only the independently recorded digitization bounds and is named accordingly.
    """
    targets = tuple(targets)
    predictions = tuple(predictions)
    if not targets or any(not isinstance(item, Jeon2022DimensionlessTarget) for item in targets):
        raise TypeError("Jeon score requires Jeon2022DimensionlessTarget observations")
    if any(not isinstance(item, Jeon2022Prediction) for item in predictions):
        raise TypeError("Jeon score requires Jeon2022Prediction values")
    if calibrated_parameter_count < 0 or wall_time_s < 0.0 or accelerated_compute_cost_usd < 0.0:
        raise ValueError("campaign counts, time, and cost must be nonnegative")

    target_by_key = {_target_key(item): item for item in targets}
    prediction_by_key = {_target_key(item): item for item in predictions}
    if len(target_by_key) != len(targets) or len(prediction_by_key) != len(predictions):
        raise ValueError("duplicate Jeon target or prediction key")
    missing = sorted(set(target_by_key) - set(prediction_by_key))
    extra = sorted(set(prediction_by_key) - set(target_by_key))
    if missing or extra:
        raise ValueError(f"prediction coverage mismatch: missing={len(missing)}, extra={len(extra)}")
    if any(item.value <= 0.0 or item.numerical_log_standard_error < 0.0
           for item in predictions):
        raise ValueError("predictions must be positive with nonnegative numerical uncertainty")

    paired = [(target, prediction_by_key[key]) for key, target in target_by_key.items()]
    nontrivial_width = [pair for pair in paired
                        if pair[0].observable == "width_shape_depth_over_200nm"
                        and pair[0].trench_width_nm != 200.0]
    calibration = [pair for pair in nontrivial_width if pair[0].split == "calibration"]
    held_out_width = [pair for pair in nontrivial_width
                      if pair[0].split == "held_out_transfer"]
    held_out_pulse = [pair for pair in paired if pair[0].observable == "pulse_depth_over_cw"]
    if (len(calibration), len(held_out_width), len(held_out_pulse)) != (5, 40, 24):
        raise ValueError("Jeon v1 split contract changed; define a new scorecard version")

    held_out = held_out_width + held_out_pulse
    coverage = sum(
        target.digitization_lower <= prediction.value <= target.digitization_upper
        for target, prediction in held_out) / len(held_out)
    reversal = [(target, prediction) for target, prediction in held_out_pulse
                if target.pulse_off_ms == 1.0]
    reversal_correct = sum(
        (prediction.value > 1.0) if target.c4f8_fraction == 0.2
        else (prediction.value < 1.0)
        for target, prediction in reversal)
    maximum_numerical_error = max(
        prediction.numerical_log_standard_error for prediction in predictions)

    calibration_rmse = _log_rmse(calibration)
    held_out_width_rmse = _log_rmse(held_out_width)
    held_out_pulse_rmse = _log_rmse(held_out_pulse)
    failures = []
    checks = (
        (calibration_rmse <= thresholds.calibration_log_rmse, "calibration_log_rmse"),
        (held_out_width_rmse <= thresholds.held_out_width_log_rmse,
         "held_out_width_log_rmse"),
        (held_out_pulse_rmse <= thresholds.held_out_pulse_log_rmse,
         "held_out_pulse_log_rmse"),
        (coverage >= thresholds.held_out_digitization_interval_coverage,
         "held_out_digitization_interval_coverage"),
        (reversal_correct == len(reversal) == 12, "one_ms_pulse_reversal"),
        (maximum_numerical_error <= thresholds.numerical_log_standard_error,
         "numerical_log_standard_error"),
        (calibrated_parameter_count <= thresholds.maximum_calibrated_parameters,
         "calibrated_parameter_count"),
        (wall_time_s <= thresholds.maximum_wall_time_s, "wall_time_s"),
        (accelerated_compute_cost_usd <= thresholds.maximum_accelerated_compute_cost_usd,
         "accelerated_compute_cost_usd"),
        (all(prediction.within_declared_scope for prediction in predictions),
         "declared_validity_scope"),
        (all(prediction.parameter_evidence_supports_prediction for prediction in predictions),
         "parameter_evidence"),
    )
    failures.extend(name for passed, name in checks if not passed)
    return Jeon2022DemoScore(
        version=JEON_2022_DEMO_VERSION,
        calibration_log_rmse=calibration_rmse,
        held_out_width_log_rmse=held_out_width_rmse,
        held_out_pulse_log_rmse=held_out_pulse_rmse,
        held_out_digitization_interval_coverage=coverage,
        pulse_reversal_correct=reversal_correct,
        pulse_reversal_total=len(reversal),
        maximum_numerical_log_standard_error=maximum_numerical_error,
        calibrated_parameter_count=int(calibrated_parameter_count),
        wall_time_s=float(wall_time_s),
        accelerated_compute_cost_usd=float(accelerated_compute_cost_usd),
        passed=not failures,
        failures=tuple(failures),
    )
