from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.audit_bosch_wafer_boundary_map_calibration import (
    BoschExactWallIonResponseTable,
    _law_from_coefficients,
    _load_preregistration,
    _nodes,
    _replay_gate_metrics,
)
from scripts import seal_bosch_wafer_boundary_map_heldout_prediction as seal_v8


ROOT = Path(__file__).resolve().parents[1]
V8 = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_wafer_boundary_map_depth_extension_v8"
)
REVEAL = ROOT / "data" / "experimental" / "zenodo_17122442"


def _json(name):
    return json.loads((V8 / name).read_text(encoding="utf-8"))


def _hash(name):
    return sha256((V8 / name).read_bytes()).hexdigest()


def _synthetic_table():
    wall = _nodes(13)
    ion = _nodes(13)
    wafer = np.arange(3, dtype=float)[None, None, :, None]
    point = np.arange(89, dtype=float)[None, None, None, :]
    wall_value = wall[:, None, None, None]
    ion_value = ion[None, :, None, None]
    silicon = 1.0e-6 * (
        40.0 + 2.0 * wall_value + 3.0 * ion_value
        + 0.1 * wafer + 0.001 * point)
    oxide = 1.0e-6 * (
        0.6 + 0.02 * wall_value + 0.03 * ion_value
        + 0.001 * wafer + 0.00001 * point)
    return BoschExactWallIonResponseTable(
        experiment_keys=("a", "b", "c"),
        log_wall_multiplier_nodes=wall,
        log_ion_factor_nodes=ion,
        silicon_depth_m=silicon,
        oxide_loss_m=oxide,
    )


def test_v8_preregistration_hash_firewall_is_closed():
    payload = _load_preregistration()

    assert payload["target_firewall"]["heldout_outcomes_read_at_freeze"] is False
    assert payload["calibration_exposure"]["v8_fit_started_before_this_freeze"] is False


def test_wall_ion_response_tensor_interpolates_independent_local_queries():
    table = _synthetic_table()
    wall_query = np.exp(np.array([-0.7, 0.0, 0.8]))
    conditioned = table.condition_on_wall(wall_query)
    log_ion_query = np.linspace(-1.0, 1.0, 3 * 89).reshape(3, 89)
    predictions, silicon_derivative, oxide_derivative = conditioned.evaluate(
        np.exp(log_ion_query))
    silicon = np.stack([item.silicon_depth_m for item in predictions])
    oxide = np.stack([item.oxide_loss_m for item in predictions])
    wafer = np.arange(3, dtype=float)[:, None]
    point = np.arange(89, dtype=float)[None]
    expected_silicon = 1.0e-6 * (
        40.0 + 2.0 * np.log(wall_query)[:, None]
        + 3.0 * log_ion_query + 0.1 * wafer + 0.001 * point)
    expected_oxide = 1.0e-6 * (
        0.6 + 0.02 * np.log(wall_query)[:, None]
        + 0.03 * log_ion_query + 0.001 * wafer + 0.00001 * point)

    assert silicon == pytest.approx(expected_silicon, rel=2.0e-13, abs=1e-18)
    assert oxide == pytest.approx(expected_oxide, rel=2.0e-13, abs=1e-18)
    assert silicon_derivative == pytest.approx(
        np.full((3, 89), 3.0e-6), rel=2.0e-13, abs=1e-18)
    assert oxide_derivative == pytest.approx(
        np.full((3, 89), 0.03e-6), rel=2.0e-13, abs=1e-18)


def test_wall_ion_response_rejects_extrapolation():
    table = _synthetic_table()
    with pytest.raises(ValueError, match="wall-multiplier"):
        table.condition_on_wall(np.array([0.2, 1.0, 1.0]))
    conditioned = table.condition_on_wall(np.ones(3))
    with pytest.raises(ValueError, match="ion-factor"):
        conditioned.evaluate(np.full((3, 89), math.exp(1.5)))


def test_exact_replay_law_requires_the_complete_frozen_coefficient_count():
    law = _law_from_coefficients(3, 2, np.zeros(14))
    assert law.static_maximum_order == 3
    assert law.dynamic_maximum_order == 2
    assert len(law.static_coefficients) == 9
    assert len(law.dynamic_coefficients) == 5
    with pytest.raises(ValueError, match="dynamic_coefficients"):
        _law_from_coefficients(3, 2, np.zeros(13))


def test_exact_replay_baseline_gates_remain_strict_and_slope_is_nonstrict():
    passing = {
        "silicon_mean_mae_um": np.nextafter(0.338486, 0.0),
        "silicon_point_rmse_um": np.nextafter(0.486585, 0.0),
        "normalized_shape_rmse_percent": np.nextafter(0.636619, 0.0),
        "silicon_mean_mape_percent": 0.0,
        "oxide_mean_mae_um": 0.0,
        "selectivity_mape_percent": 0.0,
    }
    assert all(_replay_gate_metrics(passing, 0.082903).values())
    for name, boundary in (
        ("silicon_mean_mae_um", 0.338486),
        ("silicon_point_rmse_um", 0.486585),
        ("normalized_shape_rmse_percent", 0.636619),
    ):
        failing = {**passing, name: boundary}
        assert not all(_replay_gate_metrics(failing, 0.082903).values())
    assert not all(_replay_gate_metrics(
        passing, np.nextafter(0.082903, 1.0)).values())


def test_committed_v8_selection_is_the_lowest_passing_frozen_candidate():
    capacity = _json("full_calibration_capacity.json")
    fit = _json("calibration_fit.json")
    assert capacity["heldout_outcomes_read"] is False
    assert fit["heldout_outcomes_read"] is False
    assert fit["heldout_prediction_written"] is False
    assert fit["eligible_for_prediction_seal"] is False
    assert len(capacity["candidates"]) == 20
    passing = sorted(
        (row for row in fit["candidates"]
         if row["all_pre_replay_gates_pass"]),
        key=lambda row: (
            row["coefficient_count"],
            row["static_maximum_order"],
            row["dynamic_maximum_order"],
        ),
    )
    selected = fit["selected_candidate"]
    assert passing
    assert selected["coefficient_count"] == passing[0]["coefficient_count"]
    assert selected["static_maximum_order"] == 9
    assert selected["dynamic_maximum_order"] == 2
    assert selected["coefficient_count"] == 59


def test_committed_v8_exact_replay_passes_without_opening_heldout():
    interpolation = _json("interpolation_validation.json")
    replay = _json("exact_replay.json")
    assert interpolation["passes"] is True
    assert (
        interpolation["maximum_error_fraction_of_any_frozen_gate"]
        <= interpolation["frozen_maximum_error_fraction"]
    )
    assert replay["all_replay_and_refinement_gates_pass"] is True
    assert replay["full_calibration_exact"]["passes"] is True
    assert replay["whole_lot_exact"]["passes"] is True
    assert replay["refinement"]["passes"] is True
    assert replay["heldout_outcomes_read"] is False
    assert replay["heldout_prediction_written"] is False
    assert replay["eligible_for_prediction_seal"] is False
    whole_lot = replay["whole_lot_exact"]
    assert whole_lot["metrics"]["silicon_mean_mae_um"] < 0.338486
    assert whole_lot["metrics"]["silicon_point_rmse_um"] < 0.486585
    assert (
        whole_lot["metrics"]["normalized_shape_rmse_percent"] < 0.636619
    )
    assert whole_lot["within_lot_slope_mae_um_per_wafer"] <= 0.082903
    assert (
        replay["refinement"]["maximum_observable_difference_gate_fraction"]
        <= 0.25
    )


def test_committed_v8_replay_input_hashes_are_exact():
    replay = _json("exact_replay.json")
    assert replay["input_hashes"]["preregistration"] == _hash(
        "preregistration.json")
    assert replay["input_hashes"]["response_table"] == _hash(
        "exact_wall_ion_response_table.npz")
    assert replay["input_hashes"]["calibration_fit"] == _hash(
        "calibration_fit.json")
    assert replay["input_hashes"]["interpolation_validation"] == _hash(
        "interpolation_validation.json")


def test_v8_prediction_builder_cannot_open_the_mixed_outcome_asset(
        monkeypatch):
    original_open = Path.open

    def guarded_open(path, *args, **kwargs):
        if path.name == "Si_Oxide_etch_89_points.csv":
            raise AssertionError("prediction writer opened heldout outcomes")
        return original_open(path, *args, **kwargs)

    def synthetic_exact(queries, *_args, **_kwargs):
        return tuple(
            type("Prediction", (), {
                "silicon_depth_m": np.full(89, 43.0e-6),
                "oxide_loss_m": np.full(89, 0.65e-6),
            })()
            for _query in queries
        )

    monkeypatch.setattr(Path, "open", guarded_open)
    monkeypatch.setattr(seal_v8, "_exact_map_predictions", synthetic_exact)
    payload = seal_v8.build_prediction(workers=1)
    firewall = payload["target_firewall"]
    assert payload["heldout_process_record_count"] == 20
    assert len(payload["predictions"]) == 20
    assert firewall["mixed_outcome_csv_opened"] is False
    assert firewall["heldout_numeric_outcomes_read"] is False
    assert firewall["eligible_for_separate_outcome_score_after_hash_commit"] is True


def test_committed_v8_heldout_prediction_is_hash_sealed_and_unrevealed():
    prediction = _json("heldout_prediction.json")
    seal = _json("heldout_prediction_seal.json")
    prediction_hash = _hash("heldout_prediction.json")
    firewall = prediction["target_firewall"]
    assert prediction["heldout_process_record_count"] == 20
    assert len(prediction["predictions"]) == 20
    assert len(set(prediction["heldout_experiment_keys"])) == 20
    assert firewall["mixed_outcome_csv_opened"] is False
    assert firewall["heldout_numeric_outcomes_read"] is False
    assert firewall["heldout_outcomes_read"] is False
    assert firewall["heldout_prediction_written"] is True
    assert firewall["target_depth_used"] is False
    assert firewall["target_radial_map_used"] is False
    assert seal["prediction_sha256"] == prediction_hash
    assert seal["heldout_outcomes_read"] is False
    assert seal["heldout_prediction_written"] is True
    assert seal["eligible_for_separate_outcome_score_after_commit"] is True
    assert seal["code_hashes"]["prediction_writer"] == sha256(
        seal_v8.SCRIPT.read_bytes()).hexdigest()


def test_committed_v8_postseal_reveal_and_score_preserve_the_firewall():
    score = _json("heldout_score.json")
    prediction_hash = _hash("heldout_prediction.json")
    reveal_csv = REVEAL / "revealed_heldout_Si_Oxide_etch_89_points.csv"
    reveal_manifest = REVEAL / (
        "revealed_heldout_Si_Oxide_etch_89_points_manifest.json")

    assert prediction_hash == (
        "56ed2429832fe77280762fbca86cb6ffa4de3fd9687aa84f3b5cfd4ca99a3b1a")
    assert score["input_hashes"]["heldout_prediction"] == prediction_hash
    assert score["input_hashes"]["revealed_heldout_outcomes"] == sha256(
        reveal_csv.read_bytes()).hexdigest()
    assert score["input_hashes"]["reveal_manifest"] == sha256(
        reveal_manifest.read_bytes()).hexdigest()
    firewall = score["target_firewall"]
    assert firewall["prediction_committed_and_pushed_before_numeric_reveal"] is True
    assert firewall["heldout_prediction_changed_after_reveal"] is False
    assert firewall["prediction_sha256"] == prediction_hash


def test_committed_v8_postseal_heldout_score_passes_every_frozen_gate():
    score = _json("heldout_score.json")

    assert score["measured_heldout_wafer_count"] == 13
    assert score["unmeasured_heldout_process_record_count"] == 7
    assert score["all_absolute_gates_pass"] is True
    assert score["all_empirical_baseline_gates_pass"] is True
    assert score["all_frozen_heldout_gates_pass"] is True
    assert all(score["absolute_gates"].values())
    assert all(score["empirical_baselines"]["gates"].values())
    metrics = score["physics_metrics"]
    assert metrics["silicon_mean_mae_um"] < 1.0
    assert metrics["silicon_mean_mape_percent"] < 3.0
    assert metrics["silicon_point_rmse_um"] < 1.5
    assert metrics["normalized_shape_rmse_percent"] < 2.0
    assert metrics["oxide_mean_mae_um"] < 0.08
    assert metrics["selectivity_mape_percent"] < 12.0
    assert score["within_lot_depth_drift"][
        "physics_beats_zero_slope_baseline"] is True
    assert score["wafer_bootstrap"]["replicates"] == 20000
