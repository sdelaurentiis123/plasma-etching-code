"""Checksum-verified scoring for the Hwang--Giapis Figure 13 notch profile.

The open circles in Figure 13 are the Nozawa experimental contour.  They are
useful as a strict source-reproduction target, but the paper reports no
measurement uncertainty and this contour was inspected during implementation.
Consequently this module reports quantitative errors and digitization-bound
coverage while refusing to relabel the result as held-out experimental
validation.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np


HWANG_GIAPIS_1997_FIG13_PROFILE_SHA256 = (
    "6cec65637bfe85dd619b82c4a17a1ab82730fe23ef6d5b35b567bb252bf6f462")
HWANG_GIAPIS_1997_PDF_SHA256 = (
    "30a6871d6416f27e8dbbb45e9eabbca79cddf7632872f8ed185a9e193832f63d")
HWANG_GIAPIS_1997_FIG13_PAGE_IMAGE_SHA256 = (
    "94e76479f33c754d43332067bb332cc0f7bdb473407efd18d0b23f12347350b3")
HWANG_GIAPIS_1997_FIG13_CROP_SHA256 = (
    "abe6d851259c4b6f27bb8632ddfa86dd7127e586373462509ffa182fcaf82167")
HWANG_GIAPIS_1997_FIG13_REPLAY_TOLERANCE_UM = 5e-6


@dataclass(frozen=True)
class HwangGiapisFig13Observation2D:
    point_id: str
    pixel_x: float
    pixel_y: float
    notch_depth_um: float
    height_above_oxide_um: float
    x_axis_slope_um_per_pixel: float
    x_axis_intercept_um: float
    y_axis_slope_um_per_pixel: float
    y_axis_intercept_um: float
    digitization_uncertainty_um: float
    measurement_uncertainty_um: float | None
    measurement_uncertainty_semantics: str

    def __post_init__(self):
        values = np.asarray([
            self.pixel_x, self.pixel_y, self.notch_depth_um,
            self.height_above_oxide_um, self.x_axis_slope_um_per_pixel,
            self.x_axis_intercept_um, self.y_axis_slope_um_per_pixel,
            self.y_axis_intercept_um, self.digitization_uncertainty_um],
            dtype=float)
        if (not self.point_id or np.any(~np.isfinite(values))
                or self.notch_depth_um < 0.0
                or not 0.0 <= self.height_above_oxide_um <= 0.3
                or self.digitization_uncertainty_um <= 0.0
                or (self.measurement_uncertainty_um is not None
                    and (not np.isfinite(self.measurement_uncertainty_um)
                         or self.measurement_uncertainty_um <= 0.0))
                or self.measurement_uncertainty_semantics != "not_reported"
                or self.measurement_uncertainty_um is not None):
            raise ValueError("invalid Hwang--Giapis Figure 13 observation")

    @property
    def replayed_notch_depth_um(self):
        return (
            self.x_axis_slope_um_per_pixel * self.pixel_x
            + self.x_axis_intercept_um)

    @property
    def replayed_height_above_oxide_um(self):
        return (
            self.y_axis_slope_um_per_pixel * self.pixel_y
            + self.y_axis_intercept_um)


def load_hwang_giapis_1997_fig13_profile(
        path, *, expected_sha256=HWANG_GIAPIS_1997_FIG13_PROFILE_SHA256,
        verify_checksum=True):
    """Load the Figure 13 open-circle contour with pixel replay checks."""
    path = Path(path)
    payload = path.read_bytes()
    digest = sha256(payload).hexdigest()
    if verify_checksum and digest != expected_sha256:
        raise ValueError("Hwang--Giapis Figure 13 profile checksum mismatch")
    expected_fields = [
        "point_id", "series_label", "pixel_x", "pixel_y", "notch_depth_um",
        "height_above_oxide_um", "x_axis_slope_um_per_pixel",
        "x_axis_intercept_um", "y_axis_slope_um_per_pixel",
        "y_axis_intercept_um", "digitization_uncertainty_um",
        "measurement_uncertainty_um", "measurement_uncertainty_semantics",
        "source_pdf_sha256", "source_pdf_page", "source_print_page",
        "source_figure", "source_page_image_sha256", "source_crop_sha256",
        "evidence_type", "source_location",
    ]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_fields:
            raise ValueError("unexpected Hwang--Giapis Figure 13 profile schema")
        rows = list(reader)
    if len(rows) != 22:
        raise ValueError("Figure 13 evidence must contain the 22 open-circle centers")
    observations = tuple(HwangGiapisFig13Observation2D(
        point_id=row["point_id"], pixel_x=float(row["pixel_x"]),
        pixel_y=float(row["pixel_y"]),
        notch_depth_um=float(row["notch_depth_um"]),
        height_above_oxide_um=float(row["height_above_oxide_um"]),
        x_axis_slope_um_per_pixel=float(row["x_axis_slope_um_per_pixel"]),
        x_axis_intercept_um=float(row["x_axis_intercept_um"]),
        y_axis_slope_um_per_pixel=float(row["y_axis_slope_um_per_pixel"]),
        y_axis_intercept_um=float(row["y_axis_intercept_um"]),
        digitization_uncertainty_um=float(row["digitization_uncertainty_um"]),
        measurement_uncertainty_um=(
            None if row["measurement_uncertainty_um"] == ""
            else float(row["measurement_uncertainty_um"])),
        measurement_uncertainty_semantics=row["measurement_uncertainty_semantics"])
        for row in rows)
    replay_depth = np.asarray(
        [item.replayed_notch_depth_um for item in observations])
    replay_height = np.asarray(
        [item.replayed_height_above_oxide_um for item in observations])
    reported_depth = np.asarray([item.notch_depth_um for item in observations])
    reported_height = np.asarray(
        [item.height_above_oxide_um for item in observations])
    if (len({item.point_id for item in observations}) != len(observations)
            or np.max(np.abs(replay_depth - reported_depth))
            > HWANG_GIAPIS_1997_FIG13_REPLAY_TOLERANCE_UM
            or np.max(np.abs(replay_height - reported_height))
            > HWANG_GIAPIS_1997_FIG13_REPLAY_TOLERANCE_UM
            or {row["series_label"] for row in rows}
            != {"open circles (Nozawa experiment)"}
            or {row["source_pdf_sha256"] for row in rows}
            != {HWANG_GIAPIS_1997_PDF_SHA256}
            or {row["source_pdf_page"] for row in rows} != {"12"}
            or {row["source_print_page"] for row in rows} != {"81"}
            or {row["source_figure"] for row in rows} != {"Fig. 13"}
            or {row["source_page_image_sha256"] for row in rows}
            != {HWANG_GIAPIS_1997_FIG13_PAGE_IMAGE_SHA256}
            or {row["source_crop_sha256"] for row in rows}
            != {HWANG_GIAPIS_1997_FIG13_CROP_SHA256}
            or {row["evidence_type"] for row in rows}
            != {"experiment_digitized"}):
        raise ValueError("Figure 13 evidence violates pixel replay or source provenance")
    return observations


@dataclass(frozen=True)
class HwangGiapisFig13Score2D:
    observation_height_um: np.ndarray
    experimental_depth_um: np.ndarray
    predicted_depth_um: np.ndarray
    error_um: np.ndarray
    rmse_um: float
    mean_absolute_error_um: float
    maximum_absolute_error_um: float
    maximum_depth_error_um: float
    digitization_bound_coverage_fraction: float
    strict_validation_pass: bool
    claim_status: str

    def __post_init__(self):
        height = np.asarray(self.observation_height_um, dtype=float).copy()
        experiment = np.asarray(self.experimental_depth_um, dtype=float).copy()
        prediction = np.asarray(self.predicted_depth_um, dtype=float).copy()
        error = np.asarray(self.error_um, dtype=float).copy()
        values = np.asarray([
            self.rmse_um, self.mean_absolute_error_um,
            self.maximum_absolute_error_um, self.maximum_depth_error_um,
            self.digitization_bound_coverage_fraction], dtype=float)
        if (height.ndim != 1 or not height.size
                or experiment.shape != height.shape
                or prediction.shape != height.shape or error.shape != height.shape
                or np.any(~np.isfinite(height)) or np.any(~np.isfinite(experiment))
                or np.any(~np.isfinite(prediction)) or np.any(~np.isfinite(error))
                or np.any(~np.isfinite(values)) or np.any(values[:4] < 0.0)
                or not 0.0 <= self.digitization_bound_coverage_fraction <= 1.0
                or self.strict_validation_pass
                or not self.claim_status):
            raise ValueError("invalid Hwang--Giapis Figure 13 score")
        for item in (height, experiment, prediction, error):
            item.setflags(write=False)
        object.__setattr__(self, "observation_height_um", height)
        object.__setattr__(self, "experimental_depth_um", experiment)
        object.__setattr__(self, "predicted_depth_um", prediction)
        object.__setattr__(self, "error_um", error)


def score_hwang_giapis_1997_fig13_profile(
        observations, notch_depth_by_height_um, *, cell_size_um=0.005):
    """Interpolate an engine contour to the experimental marker heights.

    The prediction is sampled at cell centers.  No fit, translation, scaling,
    or contour registration is performed.
    """
    observations = tuple(observations)
    depth = np.asarray(notch_depth_by_height_um, dtype=float)
    dx = float(cell_size_um)
    if (not observations
            or any(not isinstance(item, HwangGiapisFig13Observation2D)
                   for item in observations)
            or depth.ndim != 1 or not depth.size
            or np.any(~np.isfinite(depth)) or np.any(depth < 0.0)
            or not np.isfinite(dx) or dx <= 0.0):
        raise ValueError("invalid Figure 13 prediction score inputs")
    height = (np.arange(depth.size) + 0.5) * dx
    observation_height = np.asarray(
        [item.height_above_oxide_um for item in observations])
    experiment = np.asarray([item.notch_depth_um for item in observations])
    uncertainty = np.asarray(
        [item.digitization_uncertainty_um for item in observations])
    prediction = np.interp(observation_height, height, depth)
    error = prediction - experiment
    return HwangGiapisFig13Score2D(
        observation_height, experiment, prediction, error,
        float(np.sqrt(np.mean(error * error))),
        float(np.mean(np.abs(error))),
        float(np.max(np.abs(error))),
        float(abs(np.max(depth) - np.max(experiment))),
        float(np.mean(np.abs(error) <= uncertainty)),
        False,
        (
            "QUANTITATIVE_DEVELOPMENT_REPLAY_ONLY: measurement uncertainty "
            "is unreported and the contour is not held out; no model fitting "
            "or contour registration is permitted"))
