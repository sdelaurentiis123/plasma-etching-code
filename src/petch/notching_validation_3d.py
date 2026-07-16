"""Preregistered Nozawa/Hwang notching evidence and held-out scoring contracts.

This module deliberately contains no historical hard-coded plot reads.  Experimental rows enter only
through a checksum-verified, pixel-replayable table.  Calibration metadata are committed before a
held-out score is formed, and the score never mutates engine parameters.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np


NOTCHING_VALIDATION_PROTOCOL = "NOZAWA-C4-PRIMARY-DIGITIZATION-2026-07-15-R1"
NOZAWA_1995_PDF_SHA256 = (
    "87500f53f0286aae0597b14168b0991791db73138ee8d67a5e5e3bc56b329c67")
NOZAWA_1995_NOTCH_CURVES_SHA256 = (
    "2e472385e002aebf94f2f0ec299f877180786b6095791ec3079a80fc48a22ec2")
# These are the three quantitative curves selected before an engine score.  Figures 8 and 9
# deliberately use opposite electrical-connectivity topologies, so matching both is a stronger
# transfer test than matching one generic "space width" curve.
_TARGET_FAMILIES = frozenset({
    "open_area_width", "shared_pad_space_width", "individual_pad_space_width"})
_SPLITS = frozenset({"calibration", "held_out_transfer"})
_UNCERTAINTY_SEMANTICS = frozenset({
    "reported_bound", "reported_standard_deviation", "not_reported"})
_AXIS_TRANSFORMS = frozenset({"linear", "log10"})


def _is_sha256(value):
    return (isinstance(value, str) and len(value) == 64
            and all(character in "0123456789abcdef" for character in value))


@dataclass(frozen=True)
class Nozawa1995NotchObservation3D:
    condition_id: str
    target_family: str
    series_label: str
    control_label: str
    control_value: float
    control_unit: str
    control_pixel_x: float
    control_axis_transform: str
    control_axis_slope_per_pixel: float
    control_axis_intercept: float
    control_digitization_uncertainty: float
    notch_depth_um: float
    notch_pixel_y: float
    notch_axis_slope_um_per_pixel: float
    notch_axis_intercept_um: float
    digitization_uncertainty_um: float
    measurement_uncertainty_um: float | None
    measurement_uncertainty_semantics: str
    source_figure: str
    source_pdf_sha256: str
    source_pdf_page: int
    source_print_page: int
    source_image_sha256: str
    evidence_type: str
    split: str
    source_location: str

    def __post_init__(self):
        numeric = np.asarray([
            self.control_value, self.control_pixel_x,
            self.control_axis_slope_per_pixel, self.control_axis_intercept,
            self.control_digitization_uncertainty,
            self.notch_depth_um, self.notch_pixel_y,
            self.notch_axis_slope_um_per_pixel, self.notch_axis_intercept_um,
            self.digitization_uncertainty_um], dtype=float)
        if (not self.condition_id or self.target_family not in _TARGET_FAMILIES
                or not self.series_label or not self.control_label or not self.control_unit
                or np.any(~np.isfinite(numeric)) or self.control_value <= 0.0
                or self.control_digitization_uncertainty <= 0.0
                or self.control_axis_transform not in _AXIS_TRANSFORMS
                or self.notch_depth_um < 0.0
                or self.digitization_uncertainty_um <= 0.0
                or (self.measurement_uncertainty_um is not None
                    and (not np.isfinite(self.measurement_uncertainty_um)
                         or self.measurement_uncertainty_um <= 0.0))
                or self.measurement_uncertainty_semantics not in _UNCERTAINTY_SEMANTICS
                or ((self.measurement_uncertainty_um is None)
                    != (self.measurement_uncertainty_semantics == "not_reported"))
                or not self.source_figure
                or self.source_pdf_sha256 != NOZAWA_1995_PDF_SHA256
                or self.source_pdf_page <= 0 or self.source_print_page <= 0
                or not _is_sha256(self.source_image_sha256)
                or self.evidence_type != "experiment_digitized"
                or self.split not in _SPLITS or not self.source_location):
            raise ValueError("invalid Nozawa 1995 notch observation")

    @property
    def replayed_control_value(self):
        transformed = (self.control_axis_slope_per_pixel * self.control_pixel_x
                       + self.control_axis_intercept)
        return 10.0 ** transformed if self.control_axis_transform == "log10" else transformed


def load_nozawa_1995_notch_observations(
        path, *, expected_sha256, verify_checksum=True):
    """Load a figure transcription only when its content and pixel mapping are replayable."""
    path = Path(path)
    if not _is_sha256(expected_sha256):
        raise ValueError("a preregistered lowercase source-table SHA-256 is required")
    payload = path.read_bytes()
    if verify_checksum and sha256(payload).hexdigest() != expected_sha256:
        raise ValueError(f"checksum mismatch for Nozawa notch evidence: {path}")
    expected = [
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
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected:
            raise ValueError("unexpected Nozawa notch-evidence schema")
        raw = list(reader)
    if not raw:
        raise ValueError("Nozawa notch evidence is empty")
    observations = tuple(Nozawa1995NotchObservation3D(
        condition_id=row["condition_id"], target_family=row["target_family"],
        series_label=row["series_label"],
        control_label=row["control_label"],
        control_value=float(row["control_value"]), control_unit=row["control_unit"],
        control_pixel_x=float(row["control_pixel_x"]),
        control_axis_transform=row["control_axis_transform"],
        control_axis_slope_per_pixel=float(row["control_axis_slope_per_pixel"]),
        control_axis_intercept=float(row["control_axis_intercept"]),
        control_digitization_uncertainty=float(
            row["control_digitization_uncertainty"]),
        notch_depth_um=float(row["notch_depth_um"]),
        notch_pixel_y=float(row["notch_pixel_y"]),
        notch_axis_slope_um_per_pixel=float(row["notch_axis_slope_um_per_pixel"]),
        notch_axis_intercept_um=float(row["notch_axis_intercept_um"]),
        digitization_uncertainty_um=float(row["digitization_uncertainty_um"]),
        measurement_uncertainty_um=(
            None if row["measurement_uncertainty_um"] == ""
            else float(row["measurement_uncertainty_um"])),
        measurement_uncertainty_semantics=row["measurement_uncertainty_semantics"],
        source_figure=row["source_figure"],
        source_pdf_sha256=row["source_pdf_sha256"],
        source_pdf_page=int(row["source_pdf_page"]),
        source_print_page=int(row["source_print_page"]),
        source_image_sha256=row["source_image_sha256"],
        evidence_type=row["evidence_type"], split=row["split"],
        source_location=row["source_location"]) for row in raw)
    depth_replay = np.asarray([
        item.notch_axis_slope_um_per_pixel * item.notch_pixel_y
        + item.notch_axis_intercept_um for item in observations])
    depth_reported = np.asarray([item.notch_depth_um for item in observations])
    control_replay = np.asarray([item.replayed_control_value for item in observations])
    control_reported = np.asarray([item.control_value for item in observations])
    if (len({item.condition_id for item in observations}) != len(observations)
            or np.max(np.abs(depth_replay - depth_reported)) > 5e-6
            or np.max(np.abs(control_replay - control_reported)) > 5e-6
            or {item.target_family for item in observations} != _TARGET_FAMILIES
            or not any(item.split == "calibration" for item in observations)
            or not any(item.split == "held_out_transfer" for item in observations)):
        raise ValueError("Nozawa evidence violates pixel replay, uniqueness, or target coverage")
    return observations


@dataclass(frozen=True)
class NotchingBenchmarkProtocol3D:
    observations: tuple[Nozawa1995NotchObservation3D, ...]
    source_csv_sha256: str
    calibration_parameter_bounds: Mapping[str, tuple[float, float]]
    stationarity_contract_revision: str
    stationarity_contract_approved: bool = False
    revision: str = NOTCHING_VALIDATION_PROTOCOL

    def __post_init__(self):
        observations = tuple(self.observations)
        bounds = {
            str(name): tuple(float(value) for value in interval)
            for name, interval in dict(self.calibration_parameter_bounds).items()}
        calibration = [item for item in observations if item.split == "calibration"]
        held_out = [item for item in observations if item.split == "held_out_transfer"]
        if (not observations
                or any(not isinstance(item, Nozawa1995NotchObservation3D)
                       for item in observations)
                or not _is_sha256(self.source_csv_sha256)
                or not 1 <= len(bounds) <= 2
                or any(not name or len(interval) != 2 or not np.all(np.isfinite(interval))
                       or interval[1] <= interval[0] for name, interval in bounds.items())
                or not calibration or not held_out
                or len({item.target_family for item in calibration}) != 1
                or {item.target_family for item in held_out}
                    != (_TARGET_FAMILIES - {calibration[0].target_family})
                or not self.stationarity_contract_revision
                or not isinstance(self.stationarity_contract_approved, (bool, np.bool_))
                or self.revision != NOTCHING_VALIDATION_PROTOCOL):
            raise ValueError(
                "C4 requires <=2 bounded parameters on one curve and all other families held out")
        object.__setattr__(self, "observations", observations)
        object.__setattr__(self, "calibration_parameter_bounds", MappingProxyType(bounds))
        object.__setattr__(
            self, "stationarity_contract_approved", bool(self.stationarity_contract_approved))

    @property
    def calibration_observations(self):
        return tuple(item for item in self.observations if item.split == "calibration")

    @property
    def held_out_observations(self):
        return tuple(item for item in self.observations if item.split == "held_out_transfer")

    @property
    def commit_sha256(self):
        payload = dict(
            revision=self.revision,
            source_csv_sha256=self.source_csv_sha256,
            calibration_condition_ids=sorted(
                item.condition_id for item in self.calibration_observations),
            held_out_condition_ids=sorted(
                item.condition_id for item in self.held_out_observations),
            calibration_parameter_bounds={
                key: list(self.calibration_parameter_bounds[key])
                for key in sorted(self.calibration_parameter_bounds)},
            stationarity_contract_revision=self.stationarity_contract_revision,
            stationarity_contract_approved=self.stationarity_contract_approved)
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        return sha256(encoded).hexdigest()


@dataclass(frozen=True)
class NotchingCalibrationReveal3D:
    protocol_commit_sha256: str
    parameter_values: Mapping[str, float]
    reveal_sha256: str

    @classmethod
    def from_protocol(cls, protocol: NotchingBenchmarkProtocol3D, parameter_values):
        if not isinstance(protocol, NotchingBenchmarkProtocol3D):
            raise TypeError("calibration reveal requires a notching protocol")
        values = {str(key): float(value) for key, value in dict(parameter_values).items()}
        if set(values) != set(protocol.calibration_parameter_bounds):
            raise ValueError("calibration reveal must contain exactly the preregistered parameters")
        if any(not np.isfinite(value)
               or not protocol.calibration_parameter_bounds[name][0] <= value <= (
                   protocol.calibration_parameter_bounds[name][1])
               for name, value in values.items()):
            raise ValueError("calibration reveal lies outside a preregistered parameter bound")
        encoded = json.dumps(dict(
            protocol_commit_sha256=protocol.commit_sha256,
            parameter_values={key: values[key] for key in sorted(values)}),
            sort_keys=True, separators=(",", ":")).encode()
        return cls(protocol.commit_sha256, values, sha256(encoded).hexdigest())

    def __post_init__(self):
        values = {str(key): float(value) for key, value in dict(self.parameter_values).items()}
        if (not _is_sha256(self.protocol_commit_sha256)
                or not values or any(not np.isfinite(value) for value in values.values())
                or not _is_sha256(self.reveal_sha256)):
            raise ValueError("invalid calibration reveal")
        object.__setattr__(self, "parameter_values", MappingProxyType(values))


@dataclass(frozen=True)
class NotchingHeldOutPrediction3D:
    condition_id: str
    predicted_notch_depth_um: float
    charging_off_notch_depth_um: float
    grid_uncertainty_um: float
    timestep_uncertainty_um: float
    sample_uncertainty_um: float
    parameter_uncertainty_um: float
    exact_hard_visibility: bool
    physical_validity_supports_prediction: bool
    run_manifest_sha256: str
    calibration_reveal_sha256: str

    def __post_init__(self):
        values = np.asarray([
            self.predicted_notch_depth_um, self.charging_off_notch_depth_um,
            self.grid_uncertainty_um, self.timestep_uncertainty_um,
            self.sample_uncertainty_um, self.parameter_uncertainty_um], dtype=float)
        if (not self.condition_id or np.any(~np.isfinite(values)) or np.any(values < 0.0)
                or not isinstance(self.exact_hard_visibility, (bool, np.bool_))
                or not isinstance(self.physical_validity_supports_prediction, (bool, np.bool_))
                or not _is_sha256(self.run_manifest_sha256)
                or not _is_sha256(self.calibration_reveal_sha256)):
            raise ValueError("invalid held-out notching prediction")

    @property
    def numerical_uncertainty_um(self):
        # These are declared error bounds, not independent standard deviations.
        return (self.grid_uncertainty_um + self.timestep_uncertainty_um
                + self.sample_uncertainty_um + self.parameter_uncertainty_um)


@dataclass(frozen=True)
class NotchingBenchmarkScore3D:
    numerical_gate_passed: bool
    validated_notch_prediction: bool
    absolute_error_um: Mapping[str, float]
    combined_uncertainty_bound_um: Mapping[str, float]
    reasons: tuple[str, ...]
    protocol_commit_sha256: str
    calibration_reveal_sha256: str

    def __post_init__(self):
        error = {str(key): float(value) for key, value in dict(self.absolute_error_um).items()}
        uncertainty = {
            str(key): float(value)
            for key, value in dict(self.combined_uncertainty_bound_um).items()}
        if (set(error) != set(uncertainty) or not error
                or any(not np.isfinite(value) or value < 0.0
                       for value in (*error.values(), *uncertainty.values()))
                or not _is_sha256(self.protocol_commit_sha256)
                or not _is_sha256(self.calibration_reveal_sha256)):
            raise ValueError("invalid notching benchmark score")
        object.__setattr__(self, "numerical_gate_passed", bool(self.numerical_gate_passed))
        object.__setattr__(
            self, "validated_notch_prediction", bool(self.validated_notch_prediction))
        object.__setattr__(self, "absolute_error_um", MappingProxyType(error))
        object.__setattr__(
            self, "combined_uncertainty_bound_um", MappingProxyType(uncertainty))
        object.__setattr__(self, "reasons", tuple(self.reasons))


def score_notching_benchmark_3d(
        protocol: NotchingBenchmarkProtocol3D,
        reveal: NotchingCalibrationReveal3D, predictions):
    """Score held-out rows only; calibration rows cannot enter this function."""
    if (not isinstance(protocol, NotchingBenchmarkProtocol3D)
            or not isinstance(reveal, NotchingCalibrationReveal3D)
            or reveal.protocol_commit_sha256 != protocol.commit_sha256):
        raise ValueError("notching score requires the matching committed protocol and reveal")
    predictions = tuple(predictions)
    if any(not isinstance(item, NotchingHeldOutPrediction3D) for item in predictions):
        raise TypeError("notching predictions have the wrong type")
    by_id = {item.condition_id: item for item in predictions}
    expected = {item.condition_id for item in protocol.held_out_observations}
    if (len(by_id) != len(predictions) or set(by_id) != expected
            or any(item.calibration_reveal_sha256 != reveal.reveal_sha256
                   for item in predictions)):
        raise ValueError("predictions must cover every and only held-out condition after reveal")
    errors = {}
    bounds = {}
    reasons = []
    measurement_complete = True
    for observation in protocol.held_out_observations:
        prediction = by_id[observation.condition_id]
        measurement = observation.measurement_uncertainty_um
        if measurement is None:
            measurement_complete = False
            measurement = 0.0
        bound = (observation.digitization_uncertainty_um + measurement
                 + prediction.numerical_uncertainty_um)
        error = abs(prediction.predicted_notch_depth_um - observation.notch_depth_um)
        errors[observation.condition_id] = error
        bounds[observation.condition_id] = bound
        if error > bound:
            reasons.append(f"{observation.condition_id}: held-out error exceeds uncertainty")
        if (prediction.predicted_notch_depth_um - prediction.charging_off_notch_depth_um
                <= prediction.numerical_uncertainty_um):
            reasons.append(f"{observation.condition_id}: charging-specific notch is unresolved")
        if not prediction.exact_hard_visibility:
            reasons.append(f"{observation.condition_id}: final operator is not hard visibility")
        if not prediction.physical_validity_supports_prediction:
            reasons.append(f"{observation.condition_id}: physical parameter evidence is nonpredictive")
    numerical_pass = not reasons
    validation_reasons = list(reasons)
    if not measurement_complete:
        validation_reasons.append("published measurement uncertainty is not quantified")
    if not protocol.stationarity_contract_approved:
        validation_reasons.append("charging stationarity contract revision is not approved")
    return NotchingBenchmarkScore3D(
        numerical_pass, not validation_reasons, errors, bounds,
        tuple(validation_reasons), protocol.commit_sha256, reveal.reveal_sha256)
