import csv
import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DIGITIZED = (
    ROOT / "research_sources" / "digitized"
    / "mahoney_1994_table1_argon_100W.csv")
MANIFEST = DIGITIZED.with_name(
    "mahoney_1994_table1_argon_100W_manifest.md")
GRADER = ROOT / "scripts" / "grade_reactor_global_argon_mahoney_1994.py"
POWER_DIAGNOSTIC = (
    ROOT / "scripts" / "diagnose_reactor_global_argon_mahoney_power.py")


def _load_grader():
    specification = importlib.util.spec_from_file_location(
        "grade_reactor_global_argon_mahoney_1994", GRADER)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _load_power_diagnostic():
    specification = importlib.util.spec_from_file_location(
        "diagnose_reactor_global_argon_mahoney_power", POWER_DIAGNOSTIC)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_mahoney_table_is_visually_audited_and_retains_duplicate_pressure():
    manifest = MANIFEST.read_text(encoding="utf-8")
    with DIGITIZED.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert "native-resolution visual inspection" in manifest
    assert len(rows) == 5
    pressures = np.asarray([float(row["pressure_mTorr"]) for row in rows])
    assert np.count_nonzero(pressures == 20.0) == 2
    assert {
        row["pump_mode"] for row in rows if float(row["pressure_mTorr"]) == 20
    } == {"throttled_cryogenic", "throttled_mechanical"}


def test_mahoney_independent_board_preserves_frozen_failure():
    _, grade = _load_grader().run_board()
    assert grade["coefficient_selection_target"] is None
    assert not grade["passed"]
    assert len(grade["members"]) == 4
    assert not any(member["passed"] for member in grade["members"])
    assert all(
        member["maximum_normalized_residual"] <= 1.0e-8
        for member in grade["members"]
    )
    assert all(
        member["normalized_density_shape_log_rmse"] <= np.log(2.0)
        for member in grade["members"]
    )
    assert any(
        member["maximum_model_to_measured_density_ratio"] > 5.0
        for member in grade["members"]
    )
    assert "upper-bound" in grade["power_boundary"]


def test_committed_mahoney_grade_matches_frozen_failure():
    grade_path = (
        ROOT / "results" / "curated" / "reactor_global_argon"
        / "mahoney_1994_upper_bound_grade.json")
    grade = json.loads(grade_path.read_text(encoding="utf-8"))
    assert not grade["passed"]
    assert grade["reference_rows"] == 5


def test_target_inverted_power_is_quarantined_and_misses_hopwood_context():
    _, summary = _load_power_diagnostic().run_diagnostic()
    assert summary["claim_class"].startswith("target-informed diagnostic")
    assert summary["coefficient_selection_target"] is not None
    assert len(summary["members"]) == 4
    assert all(
        member["constant_transfer_fraction_exists"]
        for member in summary["members"]
    )
    assert not any(
        member["overlaps_hopwood_70_90_percent_context"]
        for member in summary["members"]
    )
    intervals = [
        member["constant_transfer_fraction_intersection"]
        for member in summary["members"]
    ]
    assert all(0.0 < interval[0] <= interval[1] <= 1.0 for interval in intervals)
    # This target-informed diagnostic spans materially different gas/wall
    # assumptions, but every intersection remains below the independent
    # Hopwood 70--90% context.  Avoid pinning a rounding-adjacent 25% value:
    # the Phelps equal-mass energy-frame correction legitimately moved the
    # lowest endpoint from 21.9% to 25.1% without changing that conclusion.
    assert min(interval[0] for interval in intervals) < 0.30
    assert max(interval[1] for interval in intervals) > 0.65
