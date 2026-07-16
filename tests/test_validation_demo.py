from pathlib import Path

import pytest

from petch.experimental_data import (
    build_jeon_2022_dimensionless_targets, load_jeon_2022_trench_depths,
)
from petch.validation_demo import Jeon2022Prediction, score_jeon_2022_demo


DATA = (Path(__file__).parents[1] / "data" / "experimental" / "jeon_2022"
        / "digitized_trench_depths.csv")


def _targets():
    # The v1 scorecard includes cross-exposure pulse/CW ratios.  Their exposure basis is not
    # reported by Jeon, so production target construction withholds them by default.  Exercise the
    # legacy v1 scorecard only under its explicit common-wall-time development hypothesis.
    return build_jeon_2022_dimensionless_targets(
        load_jeon_2022_trench_depths(DATA), pulse_exposure_basis="common_wall_time")


def _predictions(targets, *, values=None, numerical_error=0.0, evidence=True, scope=True):
    values = values or {}
    return tuple(Jeon2022Prediction(
        observable=item.observable,
        source_figure=item.source_figure,
        condition_family=item.condition_family,
        c4f8_fraction=item.c4f8_fraction,
        pulse_off_ms=item.pulse_off_ms,
        trench_width_nm=item.trench_width_nm,
        value=values.get((item.observable, item.source_figure, item.pulse_off_ms,
                          item.trench_width_nm), item.value),
        numerical_log_standard_error=numerical_error,
        within_declared_scope=scope,
        parameter_evidence_supports_prediction=evidence,
    ) for item in targets)


def test_jeon_v1_scorecard_accepts_complete_exact_predictive_transfer():
    targets = _targets()
    result = score_jeon_2022_demo(
        targets, _predictions(targets), calibrated_parameter_count=3, wall_time_s=100.0)

    assert result.passed
    assert result.failures == ()
    assert result.pulse_reversal_correct == result.pulse_reversal_total == 12
    assert result.held_out_digitization_interval_coverage == 1.0


def test_jeon_v1_scorecard_refuses_missing_prediction_coverage():
    targets = _targets()
    with pytest.raises(ValueError, match="prediction coverage mismatch"):
        score_jeon_2022_demo(
            targets, _predictions(targets)[:-1], calibrated_parameter_count=0, wall_time_s=1.0)


def test_jeon_v1_scorecard_separates_physics_evidence_numerics_and_campaign_failures():
    targets = _targets()
    values = {
        (item.observable, item.source_figure, item.pulse_off_ms, item.trench_width_nm):
            (0.95 if item.observable == "pulse_depth_over_cw" else item.value)
        for item in targets if item.observable == "pulse_depth_over_cw"
    }
    result = score_jeon_2022_demo(
        targets,
        _predictions(targets, values=values, numerical_error=0.04,
                     evidence=False, scope=False),
        calibrated_parameter_count=4,
        wall_time_s=7201.0,
        accelerated_compute_cost_usd=2.01,
    )

    assert not result.passed
    assert set(result.failures) >= {
        "held_out_pulse_log_rmse",
        "one_ms_pulse_reversal",
        "numerical_log_standard_error",
        "calibrated_parameter_count",
        "wall_time_s",
        "accelerated_compute_cost_usd",
        "declared_validity_scope",
        "parameter_evidence",
    }
