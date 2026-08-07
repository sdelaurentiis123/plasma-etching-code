import pytest

from scripts.audit_guo_krueger_matched_time_gate import (
    _history_through_gate,
    interpolate_metrics,
)


def test_matched_time_interpolation_is_exact_at_nodes_and_linear_between():
    history = [
        {
            "physical_time_s": 0.5,
            "metrics": {"etch_depth_nm": 5.0, "mask_opening_nm": 80.0},
        },
        {
            "physical_time_s": 0.8,
            "metrics": {"etch_depth_nm": 8.0, "mask_opening_nm": 74.0},
        },
    ]
    assert interpolate_metrics(history, 0.5) == {
        "etch_depth_nm": 5.0,
        "mask_opening_nm": 80.0,
    }
    assert interpolate_metrics(history, 0.65) == pytest.approx(
        {"etch_depth_nm": 6.5, "mask_opening_nm": 77.0}
    )


def test_matched_time_interpolation_refuses_extrapolation():
    history = [
        {
            "physical_time_s": 0.5,
            "metrics": {"etch_depth_nm": 5.0, "mask_opening_nm": 80.0},
        },
        {
            "physical_time_s": 0.8,
            "metrics": {"etch_depth_nm": 8.0, "mask_opening_nm": 74.0},
        },
    ]
    with pytest.raises(ValueError, match="leaves"):
        interpolate_metrics(history, 0.4)


def test_gate_diagnostics_ignore_later_endpoint_history():
    audit = {
        "history": [
            {"physical_time_s": 0.0},
            {"physical_time_s": 4.0},
            {"physical_time_s": 4.1},
            {"physical_time_s": 60.0},
        ]
    }
    assert [row["physical_time_s"] for row in _history_through_gate(audit)] == [
        0.0,
        4.0,
    ]
