from pathlib import Path
import sys

import pytest


SCRIPTS = Path(__file__).parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import viennals_topology_comparator as comparator


def test_parse_viennals_receipt_separates_independent_reverse_branch(tmp_path):
    receipt = tmp_path / "probe.csv"
    receipt.write_text(
        "phase,phase_time_s,physical_time_s,components,void_points,active_points\n"
        "initial,0,0,2,0,10\n"
        "coat,1.75,1.75,3,20,11\n"
        "sealed_reference,0,0,3,30,12\n"
        "etch,0.5,0.5,3,30,13\n"
        "etch,1.0,1.0,2,0,14\n")

    parsed = comparator.parse_viennals_csv(receipt, 0.05)

    assert parsed["closure_time_s"] == 1.75
    assert parsed["reopening_time_s"] == 1.0
    assert parsed["sealed_reference_void_points"] == 30


def test_paired_gate_uses_dimensional_localization_budget():
    petch = [
        {"dx_um": 0.025, "dt_s": 0.125, "closure_time_s": 2.875,
         "reopening_time_s": 1.0},
        {"dx_um": 0.0125, "dt_s": 0.0625, "closure_time_s": 2.5625,
         "reopening_time_s": 1.0},
    ]
    vienna = [
        {"dx_um": 0.025, "dt_s": 0.125, "closure_time_s": 1.875,
         "reopening_time_s": 1.0},
        {"dx_um": 0.0125, "dt_s": 0.0625, "closure_time_s": 1.9375,
         "reopening_time_s": 1.0},
    ]

    result = comparator.compare_levels(petch, vienna)

    assert result["passed"] is True
    assert result["authoritative_level"]["petch_closure_error_s"] == pytest.approx(0.5625)
    assert result["authoritative_level"]["single_engine_closure_bound_s"] == pytest.approx(0.5625)
    assert result["authoritative_level"]["paired_closure_difference_s"] == pytest.approx(0.625)
    assert result["petch_refinement"]["passed"] is True
    assert result["viennals_refinement"]["passed"] is True
