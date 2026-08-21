from hashlib import md5, sha256
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("h5py")

from petch.bosch_process_data import (
    PROCESS_DATA_MD5,
    PROCESS_DICTIONARY_MD5,
    WAFER_MEASUREMENT_89_POINT_MD5,
    load_bosch_wafer_measurements_89pt,
    load_bosch_process_traces,
    process_ingestion_manifest,
    summarize_bosch_process_traces,
)
from scripts.extract_bosch_calibration_measurements import (
    MANIFEST as CALIBRATION_MANIFEST,
    OUTPUT as CALIBRATION_MEASUREMENTS,
    main as extract_calibration_main,
)
from scripts.audit_bosch_axisymmetry_model_form import (
    OUTPUT as AXISYMMETRY_AUDIT,
    build_audit as build_axisymmetry_audit,
    main as audit_axisymmetry_main,
)
from scripts.extract_zenodo_bosch_process_features import main as extract_main
from scripts.preregister_zenodo_bosch_reactor_depth_holdout import (
    HELDOUT_DATES,
    build_preregistration,
    main as preregister_main,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
PROCESS = DATA / "Process_data.nc"
DICTIONARY = DATA / "Dictionary_process.nc"
SUMMARY = DATA / "process_wafer_summary.csv"
SUMMARY_MANIFEST = DATA / "process_wafer_summary_manifest.json"
PREREGISTRATION = (
    ROOT / "results" / "curated" / "zenodo_bosch_reactor_depth_holdout_v1"
    / "preregistration.json"
)
CYLINDRICAL_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_cylindrical_depth_extension_v2"
    / "preregistration.json"
)
SPECIES_RADIAL_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_species_radial_depth_extension_v3"
    / "preregistration.json"
)
EDGE_SHEATH_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_edge_sheath_depth_extension_v4"
    / "preregistration.json"
)
WALL_CONDITIONING_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_wall_conditioning_depth_extension_v5"
    / "preregistration.json"
)
WALL_CONDITIONING_CALIBRATION_FIT = (
    WALL_CONDITIONING_PREREGISTRATION.parent / "calibration_fit.json"
)


@pytest.fixture(scope="module")
def traces():
    return load_bosch_process_traces(PROCESS, DICTIONARY)


def test_official_process_files_are_checksum_pinned():
    assert md5(PROCESS.read_bytes()).hexdigest() == PROCESS_DATA_MD5
    assert md5(DICTIONARY.read_bytes()).hexdigest() == PROCESS_DICTIONARY_MD5
    manifest = process_ingestion_manifest()
    assert manifest["experimental_outcomes_read"] is False
    assert manifest["expected_wafer_records"] == 96


def test_process_decoder_recovers_every_wafer_without_target_labels(traces):
    assert len(traces) == 96
    assert len({trace.experiment_key for trace in traces}) == 96
    assert {len(trace.channels) for trace in traces} == {31, 44}
    assert all(0.19 <= np.median(np.diff(trace.elapsed_s)) <= 0.21 for trace in traces)
    forbidden = {"si_etch", "oxide_etch", "stepheight", "target_depth"}
    assert all(
        name.lower() not in forbidden
        for trace in traces for name in trace.channels)


def test_phase_summary_recovers_declared_hundred_cycle_structure(traces):
    summaries = summarize_bosch_process_traces(traces)
    assert len(summaries) == 96
    assert {summary.metrics["c4f8_episode_count"] for summary in summaries} == {100}
    assert all(
        100 <= summary.metrics["sf6_episode_count"] <= 102
        for summary in summaries)
    assert all(
        420.0 < summary.metrics["sf6_above_threshold_s"] < 450.0
        for summary in summaries)
    assert all(
        summary.metrics["sf6_source_load_power_q50"] > 2700.0
        for summary in summaries)
    assert all(
        summary.metrics["c4f8_platen_peak_to_peak_q50"] >
        summary.metrics["sf6_platen_peak_to_peak_q50"]
        for summary in summaries)


def test_committed_process_extraction_is_exactly_replayable():
    assert extract_main(["--check"]) == 0
    manifest = json.loads(SUMMARY_MANIFEST.read_text())
    assert manifest["summary_row_count"] == 96
    assert manifest["calculated_without_measurement_csv"] is True


def test_preregistration_never_opens_measurement_outcomes(monkeypatch):
    original_open = Path.open

    def guarded_open(path, *args, **kwargs):
        if path.name.startswith("Si_Oxide_etch_"):
            raise AssertionError("preregistration opened an experimental outcome")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded_open)
    payload = build_preregistration(SUMMARY)
    assert payload["target_firewall"]["outcome_files_opened_during_preregistration"] is False
    assert payload["target_firewall"]["preexposure_blind"] is False
    assert payload["split_rule"]["calibration_process_record_count"] == 76
    assert payload["split_rule"]["heldout_process_record_count"] == 20
    assert set(HELDOUT_DATES) == {"2024-08-21", "2024-08-22"}


def test_committed_preregistration_is_current():
    assert preregister_main(["--check"]) == 0
    payload = json.loads(PREREGISTRATION.read_text())
    assert payload["forbidden_shortcuts"]
    assert payload["heldout_score"]["absolute_acceptance"]


def test_calibration_measurement_broker_is_replayable_and_excludes_holdout():
    assert extract_calibration_main(["--check"]) == 0
    manifest = json.loads(CALIBRATION_MANIFEST.read_text())
    preregistration = json.loads(PREREGISTRATION.read_text())

    assert manifest["source_md5"] == WAFER_MEASUREMENT_89_POINT_MD5
    assert manifest["splitter_numeric_outcome_fields_parsed"] is False
    assert manifest["heldout_rows_copied_to_fit_asset"] is False
    assert manifest["output_experiment_key_count"] == 75
    assert manifest["output_row_count"] == 75 * 89
    assert manifest["excluded_heldout_experiment_key_count"] == 13
    allowed = preregistration["split_rule"]["calibration_experiment_keys"]
    measurements = load_bosch_wafer_measurements_89pt(
        CALIBRATION_MEASUREMENTS, allowed_experiment_keys=allowed)
    assert len(measurements) == 75
    assert {measurement.lot_number for measurement in measurements} == set(range(1, 9))
    assert all(measurement.x_um.size == 89 for measurement in measurements)
    assert 40.0 < measurements[0].wafer_mean_silicon_depth_um < 60.0


def test_measurement_loader_refuses_a_disallowed_key_before_numeric_conversion(
        tmp_path):
    mixed = tmp_path / "mixed.csv"
    mixed.write_text(
        "experiment_key,lot_number,wafer_number,X,Y,preox_thickness,"
        "postox_thickness,postox_thickness_nan,stepheight,oxide_etch,si_etch\n"
        "2024-08-21_01,not-an-int,also-bad,bad,bad,bad,bad,N/A,bad,bad,bad\n")
    with pytest.raises(ValueError, match="outside the allowed key set"):
        load_bosch_wafer_measurements_89pt(
            mixed, allowed_experiment_keys={"2024-07-02_01"})


def test_axisymmetry_model_form_failure_is_replayable_without_heldout():
    assert audit_axisymmetry_main(["--check"]) == 0
    audit = build_axisymmetry_audit()

    assert audit["heldout_outcomes_read"] is False
    assert audit["calibration_wafer_count"] == 75
    assert audit["points_per_wafer"] == 89
    assert audit["shared_axisymmetric_model_can_pass_frozen_gate"] is False
    assert (
        audit["oracle_lower_bounds"]
        ["shared_axisymmetric_normalized_rmse_percent"] > 2.0)
    assert (
        audit["oracle_lower_bounds"]
        ["shared_unconstrained_coordinate_map_normalized_rmse_percent"] < 1.0)


def test_cylindrical_v2_is_frozen_before_fit_and_keeps_v1_surface_and_gates():
    payload = json.loads(CYLINDRICAL_PREREGISTRATION.read_text())
    audit_hash = sha256(AXISYMMETRY_AUDIT.read_bytes()).hexdigest()

    assert payload["calibration_exposure"]["cylindrical_fit_started_before_this_freeze"] is False
    assert payload["calibration_exposure"]["heldout_outcomes_examined"] is False
    assert payload["frozen_from_v1"]["acceptance_gates_unchanged"] is True
    assert payload["frozen_from_v1"]["surface_laws_unchanged"] == [
        "Belen SF6 silicon removal",
        "La Magna/Garozzo C4F8 film and SiO2 removal",
    ]
    assert payload["model_form_audit"]["sha256"] == audit_hash
    assert payload["model_extension"]["harmonic_source"]["maximum_order"] == 4
    assert payload["selection_and_seal"]["sealed_prediction_required_before_heldout_reveal"] is True
    assert payload["target_firewall"]["heldout_outcomes_read_at_freeze"] is False


def test_species_radial_v3_is_narrowly_frozen_before_fit():
    payload = json.loads(SPECIES_RADIAL_PREREGISTRATION.read_text())

    assert payload["calibration_exposure"]["heldout_outcomes_examined"] is False
    assert payload["calibration_exposure"]["species_resolved_radial_fit_started_before_this_freeze"] is False
    assert payload["frozen_from_v1_and_v2"]["acceptance_gates_unchanged"] is True
    assert payload["frozen_from_v1_and_v2"]["surface_laws_unchanged"] == [
        "Belen SF6 silicon removal",
        "La Magna/Garozzo C4F8 film and SiO2 removal",
    ]
    assert payload["model_extension"]["new_free_parameter_count"] == 4
    assert payload["frozen_code_parent"]["species_resolved_radial_source_response_implemented"] is False
    assert payload["target_firewall"]["heldout_outcomes_read_at_freeze"] is False


def test_edge_sheath_v4_is_frozen_before_fit_and_conserves_ion_current():
    payload = json.loads(EDGE_SHEATH_PREREGISTRATION.read_text())

    assert payload["calibration_exposure"]["edge_sheath_fit_started_before_this_freeze"] is False
    assert payload["calibration_exposure"]["heldout_outcomes_examined"] is False
    assert payload["frozen_from_prior_versions"]["surface_laws_unchanged"]
    assert payload["model_extension"]["affected_channel"] == "positive_ion wafer flux only"
    assert payload["model_extension"]["new_free_parameter_count"] == 3
    assert payload["model_extension"]["total_wafer_ion_current_conserved"] is True
    assert payload["frozen_code_parent"]["edge_sheath_operator_implemented"] is False
    assert payload["target_firewall"]["heldout_outcomes_read_at_freeze"] is False


def test_wall_conditioning_v5_is_frozen_before_fit_without_depth_offsets():
    payload = json.loads(WALL_CONDITIONING_PREREGISTRATION.read_text())

    assert payload["calibration_exposure"][
        "conditioning_fit_started_before_this_freeze"] is False
    assert payload["calibration_exposure"]["heldout_outcomes_examined"] is False
    assert payload["frozen_from_prior_versions"]["surface_laws_unchanged"] == [
        "Belen SF6 silicon removal",
        "La Magna/Garozzo C4F8 film and SiO2 removal",
    ]
    law = payload["conditioning_law"]
    assert law["maximum_free_coefficients"] == 3
    assert law["wall_loss_multiplier_bounds"] == [0.25, 4.0]
    assert law["physical_action"]["lower_wafer_collection"] == "unchanged"
    assert law["physical_action"]["positive_ion_state_and_transfer"] == "unchanged"
    assert any(
        "predicted depth" in shortcut
        for shortcut in payload["forbidden_shortcuts"]
    )
    assert payload["frozen_code_parent"]["conditioning_operator_implemented"] is False
    assert payload["target_firewall"]["heldout_outcomes_read_at_freeze"] is False


def test_wall_conditioning_calibration_receipt_remains_unsealed_and_target_free():
    payload = json.loads(WALL_CONDITIONING_CALIBRATION_FIT.read_text())

    assert payload["heldout_outcomes_read"] is False
    assert payload["heldout_prediction_written"] is False
    assert payload["surface_laws_changed"] is False
    assert payload["per_lot_or_per_wafer_depth_offsets"] is False
    assert payload["all_absolute_gates_pass"] is True
    assert payload["leave_one_lot_out_physics_refits_completed"] is False
    assert payload["certification_grid_refinement_completed"] is False
    assert payload["eligible_for_prediction_seal"] is False
    assert not all(payload["in_sample_physics_beats_empirical_baseline"].values())
    assert payload["input_hashes"]["calibration_measurements"] == sha256(
        CALIBRATION_MEASUREMENTS.read_bytes()).hexdigest()
