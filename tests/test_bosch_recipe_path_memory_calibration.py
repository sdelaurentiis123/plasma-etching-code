from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from scripts.audit_bosch_recipe_path_memory_calibration import (
    _inputs,
    _recipe_path_steps,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_recipe_path_memory_depth_extension_v7"
    / "calibration_fit.json"
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def test_recipe_path_history_carries_unmeasured_trace_and_resets_only_by_lot():
    measurements, process_traces, lot_type_by_date = _inputs()
    steps, _static, path = _recipe_path_steps(
        np.array([0.0, 0.0, 0.0, 0.02]),
        process_traces,
        lot_type_by_date,
    )
    measurement_keys = {item.experiment_key for item in measurements}

    assert len(process_traces) == 76
    assert len(measurements) == 75
    assert set(steps) - measurement_keys == {"2024-07-02_07"}
    assert steps["2024-07-02_08"].start_state == (
        steps["2024-07-02_07"].end_state)
    assert (
        steps["2024-07-05_01"].start_state.cumulative_reference_wafer_exposure
        == 0.0
    )
    assert (
        steps["2024-07-02_10"].end_state.cumulative_reference_wafer_exposure
        > steps["2024-07-02_01"].end_state.cumulative_reference_wafer_exposure
    )
    assert path.manifest()["out_of_domain_extrapolation_allowed"] is False


def test_recipe_path_calibration_receipt_preserves_firewall_and_hashes():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert payload["heldout_outcomes_read"] is False
    assert payload["heldout_prediction_written"] is False
    assert payload["eligible_for_prediction_seal"] is False
    assert payload["all_absolute_gates_pass"] is True
    assert payload["identifiability"]["passes"] is True
    assert payload["deterministic_acceleration"]["validation"]["passes"] is True
    assert payload["whole_lot_physics_beats_empirical_baseline"] == {
        "global_mean_depth_mae": True,
        "mean_map_normalized_shape_rmse": False,
        "mean_map_point_rmse": False,
    }
    assert payload["recipe_path_extrapolation_used"] is False
    assert payload["unmeasured_process_traces_carried_through_wall_history"] == [
        "2024-07-02_07"
    ]
    for name, relative in {
        "v6_preregistration": (
            "results/curated/zenodo_bosch_dynamic_wall_depth_extension_v6/"
            "preregistration.json"),
        "v6_calibration_fit": (
            "results/curated/zenodo_bosch_dynamic_wall_depth_extension_v6/"
            "calibration_fit.json"),
        "v6_response_table": (
            "results/curated/zenodo_bosch_dynamic_wall_depth_extension_v6/"
            "exact_response_table.npz"),
        "v7_preregistration": (
            "results/curated/zenodo_bosch_recipe_path_memory_depth_extension_v7/"
            "preregistration.json"),
    }.items():
        assert payload["input_hashes"][name] == _hash(ROOT / relative)


def test_recipe_path_exact_table_is_thirteen_node_and_content_addressed():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    acceleration = payload["deterministic_acceleration"]
    table = ROOT / acceleration["response_table"]

    assert acceleration["node_count"] == 13
    assert acceleration["validation"]["independent_midpoint_node_count"] == 12
    assert acceleration["validation"][
        "maximum_error_fraction_of_any_frozen_gate"] <= 0.05
    assert acceleration["response_table_sha256"] == _hash(table)
