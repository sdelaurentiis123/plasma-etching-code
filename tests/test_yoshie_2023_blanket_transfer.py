from pathlib import Path

import pytest

from scripts.audit_yoshie_2023_blanket_transfer import (
    RESULT,
    build_report,
    canonical_payload,
)


ROOT = Path(__file__).parents[1]


def _condition(report, cycle, timing):
    return next(
        item for item in report["conditions"]
        if item["cycle_duration_s"] == cycle and item["timing"] == timing
    )


def test_audit_freezes_the_super_blanket_condition_without_fitting():
    report = build_report()
    timing_i = _condition(report, 8, "I")

    assert timing_i["feature_to_blanket_ratio_min"] == pytest.approx(
        2.615489)
    assert timing_i[
        "ratio_lower_bound_with_digitization_allowance"
    ] > 2.5
    assert timing_i["feature_points_above_blanket"] == 7
    assert not report["identifiability_verdict"][
        "direct_scale_transfer_certified"]
    assert report["operation"].startswith("read-only comparison")


def test_eight_second_timing_order_does_not_transfer_from_blanket():
    report = build_report()
    ranks = report["rank_order_by_cycle"]["8"]

    assert ranks["blanket_descending"] == ["II", "III", "IV", "I"]
    assert ranks["feature_median_descending"] == ["II", "I", "III", "IV"]


def test_scale_only_null_is_reported_as_a_failure_not_a_depth_prediction():
    report = build_report()
    null = report["scale_only_null"]

    assert null["feature_points"] == 49
    assert null["mean_absolute_percentage_error"] == pytest.approx(
        0.3786430484)
    assert null["rmse_nm_per_bias_min"] == pytest.approx(141.9743329)


def test_committed_audit_is_exactly_replayable():
    assert RESULT.is_file()
    assert RESULT.read_text() == canonical_payload(build_report())
    assert RESULT.is_relative_to(ROOT)
