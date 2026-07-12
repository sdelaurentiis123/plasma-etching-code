from pathlib import Path

import numpy as np
import pytest

from petch.experimental_data import (
    build_jeon_2022_dimensionless_targets,
    load_bosch_wafer_measurements,
    load_bosch_wafer_measurements_89pt,
    load_jeon_2022_plasma_controls,
    load_jeon_2022_trench_depths,
    load_krueger_2024_evidence,
)


DATA = (
    Path(__file__).parents[1]
    / "data"
    / "experimental"
    / "zenodo_17122442"
    / "Si_Oxide_etch_9_points.csv"
)
KRUEGER_DATA = Path(__file__).parents[1] / "data" / "experimental" / "krueger_2024"
DATA_89 = DATA.with_name("Si_Oxide_etch_89_points.csv")
JEON_DATA = (
    Path(__file__).parents[1] / "data" / "experimental" / "jeon_2022"
    / "digitized_trench_depths.csv")
JEON_CONTROLS = JEON_DATA.with_name("digitized_plasma_controls.csv")


def test_bosch_wafer_measurements_have_verified_provenance_and_units():
    rows = load_bosch_wafer_measurements(DATA)

    assert len(rows) == 684
    identified = [row for row in rows if row.experiment_key is not None]
    unidentified = [row for row in rows if row.experiment_key is None]
    assert len({(row.experiment_key, row.wafer_number) for row in identified}) == 75
    assert len(unidentified) == 9
    assert np.isclose(min(row.silicon_etch_um for row in rows), 38.2659)
    assert np.isclose(max(row.silicon_etch_um for row in rows), 42.8646)
    assert np.isclose(min(row.oxide_etch_um for row in rows), 0.5351)
    assert np.isclose(max(row.oxide_etch_um for row in rows), 0.7417)


def test_bosch_wafer_measurements_reject_unverified_content(tmp_path):
    altered = tmp_path / "measurements.csv"
    altered.write_bytes(DATA.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="checksum mismatch"):
        load_bosch_wafer_measurements(altered)


def test_bosch_89_point_measurements_preserve_source_missingness_and_distinct_identities():
    rows = load_bosch_wafer_measurements_89pt(DATA_89)

    assert len(rows) == 7832
    assert len({(row.experiment_key, row.lot_number, row.wafer_number) for row in rows}) == 88
    assert all(row.sampling_grid == "89_point" and row.location_id is None for row in rows)
    assert sum(row.post_oxide_original_um is None for row in rows) == 157
    assert np.isclose(min(row.silicon_etch_um for row in rows), 28.653)
    assert np.isclose(max(row.silicon_etch_um for row in rows), 52.7809)
    assert np.isclose(min(row.oxide_etch_um for row in rows), 0.4248)
    assert np.isclose(max(row.oxide_etch_um for row in rows), 0.838006352)


def test_bosch_89_point_measurements_reject_unverified_content(tmp_path):
    altered = tmp_path / "measurements_89.csv"
    altered.write_bytes(DATA_89.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_bosch_wafer_measurements_89pt(altered)


def test_krueger_2024_keeps_calibration_measurements_and_simulated_inputs_separate():
    evidence = load_krueger_2024_evidence(KRUEGER_DATA)

    assert {metric.symbol for metric in evidence.calibration_metrics} == {
        "wm", "wt", "wf", "hf", "hm", "ah",
    }
    assert all(metric.evidence_type == "experiment" for metric in evidence.calibration_metrics)
    assert all(flux.evidence_type == "HPEM_simulation" for flux in evidence.boundary_fluxes)
    assert np.isclose(sum(flux.value_cm2_s for flux in evidence.boundary_fluxes), 3.984e17)
    assert all(item.split != "calibration" for item in evidence.transfer_observations)


def test_krueger_2024_exposes_held_out_experimental_transfer_trends():
    evidence = load_krueger_2024_evidence(KRUEGER_DATA)
    experimental = {
        (item.family, item.control, item.observable): item.value
        for item in evidence.transfer_observations
        if item.evidence_type in {"experiment", "experiment_and_simulation"}
    }

    assert experimental[("oxygen_ratio", "0.5", "feature_clogged")] == "true"
    assert experimental[("oxygen_ratio", "1.5", "etch_depth_rank")] == "maximum"
    assert experimental[("oxygen_ratio", "1.5_to_2.5", "etch_depth_increase")] == "false"
    assert experimental[("low_frequency_power_kw", "4_to_8", "final_profile_difference")] == (
        "few_differences")


def test_krueger_2024_rejects_unverified_transcription(tmp_path):
    target = tmp_path / "krueger_2024"
    target.mkdir()
    for source in KRUEGER_DATA.glob("*.csv"):
        (target / source.name).write_bytes(source.read_bytes())
    with (target / "transfer_observations.csv").open("ab") as stream:
        stream.write(b"\n")

    with pytest.raises(ValueError, match="checksum mismatch"):
        load_krueger_2024_evidence(target)


def test_jeon_2022_preserves_digitization_error_and_preregistered_transfer_split():
    rows = load_jeon_2022_trench_depths(JEON_DATA)

    assert len(rows) == 54
    assert {item.trench_width_nm for item in rows} == {60, 80, 100, 150, 180, 200}
    calibration = [item for item in rows if item.split == "calibration"]
    assert len(calibration) == 6
    assert {(item.source_figure, item.c4f8_fraction, item.pulse_off_ms)
            for item in calibration} == {("4b", 0.2, 0.0)}
    assert all(item.digitization_uncertainty_nm == 35.0 for item in rows)
    assert all(item.published_errorbar_semantics == "not_specified" for item in rows)


def test_jeon_2022_held_out_pulse_response_reverses_between_radical_regimes():
    rows = load_jeon_2022_trench_depths(JEON_DATA)
    by_key = {
        (item.c4f8_fraction, item.pulse_off_ms, item.trench_width_nm): item.depth_nm
        for item in rows if item.condition_family.startswith("pulse_off")}
    for width in (60, 80, 100, 150, 180, 200):
        assert by_key[(0.2, 1.0, width)] > by_key[(0.2, 0.0, width)]
        assert by_key[(0.8, 1.0, width)] < by_key[(0.8, 0.0, width)]


def test_jeon_2022_independently_digitized_duplicate_controls_close_within_budget():
    rows = load_jeon_2022_trench_depths(JEON_DATA)
    by_panel = {
        (item.source_figure, item.c4f8_fraction, item.pulse_off_ms,
         item.trench_width_nm): item
        for item in rows}
    for fraction, duplicate_figure in ((0.2, "7b"), (0.8, "9b")):
        for width in (60, 80, 100, 150, 180, 200):
            reference = by_panel[("4b", fraction, 0.0, width)]
            duplicate = by_panel[(duplicate_figure, fraction, 0.0, width)]
            assert abs(reference.depth_nm - duplicate.depth_nm) <= max(
                reference.digitization_uncertainty_nm,
                duplicate.digitization_uncertainty_nm)


def test_jeon_2022_rejects_unverified_digitization(tmp_path):
    altered = tmp_path / "digitized.csv"
    altered.write_bytes(JEON_DATA.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_jeon_2022_trench_depths(altered)


def test_jeon_2022_controls_replay_and_remain_diagnostic_derived_inputs():
    controls = load_jeon_2022_plasma_controls(JEON_CONTROLS)
    assert len(controls) == 12
    assert all(item.evidence_type == "diagnostic_derived_digitized" for item in controls)
    assert all(item.role == "physical_boundary_input" for item in controls)
    by_key = {
        (item.condition_family, item.c4f8_fraction, item.pulse_off_ms): item
        for item in controls}
    for fraction, family in ((0.2, "pulse_off_20pct"), (0.8, "pulse_off_80pct")):
        assert (by_key[(family, fraction, 0.1)].neutral_to_ion_flux_ratio
                < by_key[(family, fraction, 0.0)].neutral_to_ion_flux_ratio)
        assert (by_key[(family, fraction, 1.0)].neutral_to_ion_flux_ratio
                > by_key[(family, fraction, 0.0)].neutral_to_ion_flux_ratio)


def test_every_jeon_depth_condition_has_a_published_physical_control_ratio():
    depths = load_jeon_2022_trench_depths(JEON_DATA)
    controls = load_jeon_2022_plasma_controls(JEON_CONTROLS)
    control_conditions = {
        (item.condition_family, item.c4f8_fraction, item.pulse_off_ms)
        for item in controls}
    assert all(
        (item.condition_family, item.c4f8_fraction, item.pulse_off_ms) in control_conditions
        for item in depths)


def test_jeon_2022_controls_reject_unverified_digitization(tmp_path):
    altered = tmp_path / "controls.csv"
    altered.write_bytes(JEON_CONTROLS.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_jeon_2022_plasma_controls(altered)


def test_jeon_dimensionless_targets_remove_rate_scale_without_split_leakage():
    targets = build_jeon_2022_dimensionless_targets(
        load_jeon_2022_trench_depths(JEON_DATA))
    shapes = [item for item in targets
              if item.observable == "width_shape_depth_over_200nm"]
    pulse = [item for item in targets if item.observable == "pulse_depth_over_cw"]

    assert len(shapes) == 54
    assert len(pulse) == 24
    assert sum(item.split == "calibration" for item in targets) == 6
    assert all(item.observable == "width_shape_depth_over_200nm"
               for item in targets if item.split == "calibration")
    assert all(item.value == 1.0 for item in shapes if item.trench_width_nm == 200.0)
    assert all(item.digitization_lower == item.digitization_upper == 1.0
               for item in shapes if item.trench_width_nm == 200.0)
    assert all(item.digitization_lower <= item.value <= item.digitization_upper
               for item in targets)


def test_jeon_held_out_dimensionless_pulse_gate_preserves_regime_reversal():
    targets = build_jeon_2022_dimensionless_targets(
        load_jeon_2022_trench_depths(JEON_DATA))
    pulse_1ms = [item for item in targets
                 if item.observable == "pulse_depth_over_cw" and item.pulse_off_ms == 1.0]
    low_radical = [item for item in pulse_1ms if item.c4f8_fraction == 0.2]
    high_radical = [item for item in pulse_1ms if item.c4f8_fraction == 0.8]

    assert len(low_radical) == len(high_radical) == 6
    assert all(item.value > 1.0 for item in low_radical)
    assert all(item.value < 1.0 for item in high_radical)
    assert all(item.cancellation_assumption == "common_etch_duration_within_pulse_series"
               for item in pulse_1ms)
