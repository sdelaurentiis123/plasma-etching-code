import csv
from hashlib import sha256
from pathlib import Path

import pytest

from petch.notching_validation_3d import (
    NOZAWA_1995_NOTCH_CURVES_SHA256,
    NOZAWA_1995_PDF_SHA256,
    NotchingBenchmarkProtocol3D,
    NotchingCalibrationReveal3D,
    NotchingHeldOutPrediction3D,
    load_nozawa_1995_notch_observations,
    score_notching_benchmark_3d,
)


NOZAWA_DATA = (Path(__file__).resolve().parents[1] / "data" / "experimental"
                / "nozawa_1995" / "digitized_notch_curves.csv")


FIELDS = [
    "condition_id", "target_family", "series_label", "control_label", "control_value",
    "control_unit", "control_pixel_x", "control_axis_transform",
    "control_axis_slope_per_pixel", "control_axis_intercept",
    "control_digitization_uncertainty", "notch_depth_um", "notch_pixel_y",
    "notch_axis_slope_um_per_pixel", "notch_axis_intercept_um",
    "digitization_uncertainty_um", "measurement_uncertainty_um",
    "measurement_uncertainty_semantics", "source_figure",
    "source_pdf_sha256", "source_pdf_page", "source_print_page",
    "source_image_sha256", "evidence_type", "split", "source_location",
]


def _write_evidence(path):
    rows = [
        ("width_1", "open_area_width", "1 um", "1", "um", "0.10", "90"),
        ("shared", "shared_pad_space_width", "shared pad", "2", "um", "0.20", "80"),
        ("individual", "individual_pad_space_width", "individual pad", "3", "um", "0.30", "70"),
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for index, (condition, family, label, control, unit, depth, pixel) in enumerate(rows):
            writer.writerow(dict(
                condition_id=condition, target_family=family, control_label=label,
                series_label=family, control_value=control, control_unit=unit,
                control_pixel_x=control, control_axis_transform="linear",
                control_axis_slope_per_pixel="1.0", control_axis_intercept="0.0",
                control_digitization_uncertainty="0.01", notch_depth_um=depth,
                notch_pixel_y=pixel, notch_axis_slope_um_per_pixel="-0.01",
                notch_axis_intercept_um="1.0", digitization_uncertainty_um="0.01",
                measurement_uncertainty_um="0.02",
                measurement_uncertainty_semantics="reported_bound",
                source_figure=f"Fig. {index + 1}",
                source_pdf_sha256=NOZAWA_1995_PDF_SHA256,
                source_pdf_page=index + 1, source_print_page=2107 + index,
                source_image_sha256=str(index + 1) * 64,
                evidence_type="experiment_digitized",
                split=("calibration" if index == 0 else "held_out_transfer"),
                source_location="Nozawa et al. 1995"))
    return sha256(path.read_bytes()).hexdigest()


def _protocol(tmp_path, *, approved=True):
    path = tmp_path / "nozawa.csv"
    digest = _write_evidence(path)
    rows = load_nozawa_1995_notch_observations(path, expected_sha256=digest)
    return NotchingBenchmarkProtocol3D(
        rows, digest, {"ion_reflection_probability": (0.0, 0.5)},
        stationarity_contract_revision="CCA-profile-stationary-R3",
        stationarity_contract_approved=approved)


def test_nozawa_loader_replays_pixels_and_protocol_commits_split(tmp_path):
    protocol = _protocol(tmp_path)

    assert len(protocol.calibration_observations) == 1
    assert {item.target_family for item in protocol.held_out_observations} == {
        "shared_pad_space_width", "individual_pad_space_width"}
    assert len(protocol.commit_sha256) == 64

    altered = tmp_path / "altered.csv"
    altered.write_bytes((tmp_path / "nozawa.csv").read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_nozawa_1995_notch_observations(
            altered, expected_sha256=protocol.source_csv_sha256)


def test_primary_nozawa_curves_replay_and_preserve_opposite_topology_trends():
    rows = load_nozawa_1995_notch_observations(
        NOZAWA_DATA, expected_sha256=NOZAWA_1995_NOTCH_CURVES_SHA256)
    by_family = {
        family: sorted(
            (item for item in rows if item.target_family == family),
            key=lambda item: item.control_value)
        for family in {item.target_family for item in rows}}

    assert len(rows) == 17
    assert len(by_family["open_area_width"]) == 7
    assert len(by_family["shared_pad_space_width"]) == 5
    assert len(by_family["individual_pad_space_width"]) == 5
    assert (by_family["open_area_width"][-1].notch_depth_um
            > by_family["open_area_width"][0].notch_depth_um)
    assert (by_family["shared_pad_space_width"][-1].notch_depth_um
            < by_family["shared_pad_space_width"][0].notch_depth_um)
    assert (by_family["individual_pad_space_width"][-1].notch_depth_um
            > by_family["individual_pad_space_width"][0].notch_depth_um)
    assert all(item.measurement_uncertainty_um is None for item in rows)


def test_notching_protocol_refuses_more_than_two_calibration_parameters(tmp_path):
    protocol = _protocol(tmp_path)
    with pytest.raises(ValueError, match="<=2 bounded parameters"):
        NotchingBenchmarkProtocol3D(
            protocol.observations, protocol.source_csv_sha256,
            {"a": (0, 1), "b": (0, 1), "c": (0, 1)},
            stationarity_contract_revision="R3", stationarity_contract_approved=True)


def test_notching_score_uses_only_held_out_rows_and_all_uncertainty_terms(tmp_path):
    protocol = _protocol(tmp_path, approved=True)
    reveal = NotchingCalibrationReveal3D.from_protocol(
        protocol, {"ion_reflection_probability": 0.2})
    predictions = tuple(NotchingHeldOutPrediction3D(
        condition_id=item.condition_id,
        predicted_notch_depth_um=item.notch_depth_um,
        charging_off_notch_depth_um=0.0,
        grid_uncertainty_um=0.001, timestep_uncertainty_um=0.001,
        sample_uncertainty_um=0.001, parameter_uncertainty_um=0.001,
        exact_hard_visibility=True, physical_validity_supports_prediction=True,
        run_manifest_sha256="a" * 64,
        calibration_reveal_sha256=reveal.reveal_sha256)
        for item in protocol.held_out_observations)

    score = score_notching_benchmark_3d(protocol, reveal, predictions)

    assert score.numerical_gate_passed
    assert score.validated_notch_prediction
    assert all(value == pytest.approx(0.034) for value in
               score.combined_uncertainty_bound_um.values())

    unsigned = _protocol(tmp_path, approved=False)
    unsigned_reveal = NotchingCalibrationReveal3D.from_protocol(
        unsigned, {"ion_reflection_probability": 0.2})
    unsigned_predictions = tuple(NotchingHeldOutPrediction3D(
        condition_id=item.condition_id,
        predicted_notch_depth_um=item.notch_depth_um,
        charging_off_notch_depth_um=0.0,
        grid_uncertainty_um=0.001, timestep_uncertainty_um=0.001,
        sample_uncertainty_um=0.001, parameter_uncertainty_um=0.001,
        exact_hard_visibility=True, physical_validity_supports_prediction=True,
        run_manifest_sha256="b" * 64,
        calibration_reveal_sha256=unsigned_reveal.reveal_sha256)
        for item in unsigned.held_out_observations)
    unsigned_score = score_notching_benchmark_3d(
        unsigned, unsigned_reveal, unsigned_predictions)
    assert unsigned_score.numerical_gate_passed
    assert not unsigned_score.validated_notch_prediction
    assert "not approved" in " ".join(unsigned_score.reasons)
