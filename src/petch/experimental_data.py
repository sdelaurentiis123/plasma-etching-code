"""Load public experimental validation data with provenance and schema checks.

Experimental data are deliberately kept separate from solver parameters.  Loading a dataset does
not calibrate a model; it only produces units-explicit observations suitable for a validation or
calibration harness.
"""

from dataclasses import dataclass
from hashlib import md5, sha256
import csv
from pathlib import Path
from typing import Optional

import numpy as np


ZENODO_17122442_9PT_MD5 = "78515caf25e29e558e1859b92f8a4827"
ZENODO_17122442_89PT_MD5 = "446e75b040eea37b634eeb8f763a62fc"

KRUEGER_2024_SHA256 = {
    "base_case_metrics.csv": "5d51d124a93e1f942a9b999649b8adcf217662967d9ea2a40089f72940992351",
    "base_case_boundary_fluxes.csv": "ad50b6099a52d2c2cc00eb4eade496b9d75c41d19881c5fec9e905f9dfd3808b",
    "transfer_observations.csv": "85cef607f20ab5e56e606666aa7e0e6241abb546d0369277b21833542e04d425",
}

JEON_2022_DEPTH_SHA256 = (
    "0f737d0d44866684513e8f16fd3a4feab8618a81866ea02357a6b3f5da98310f")
JEON_2022_CONTROL_SHA256 = (
    "2c4e28cd4b3cbf34f356a5a7dd292a3a93ecdcbfbad291c44bba1b5f91c4ee8a")


@dataclass(frozen=True)
class BoschWaferMeasurement:
    experiment_key: Optional[str]
    lot_number: Optional[int]
    wafer_number: Optional[int]
    location_id: Optional[str]
    x_um: float
    y_um: float
    pre_oxide_um: float
    post_oxide_um: float
    step_height_um: float
    oxide_etch_um: float
    silicon_etch_um: float
    post_oxide_original_um: Optional[float] = None
    sampling_grid: str = "9_point"


@dataclass(frozen=True)
class ProfileTargetMetric:
    metric: str
    symbol: str
    value: float
    unit: str
    evidence_type: str
    split: str
    source_location: str


@dataclass(frozen=True)
class BoundaryFluxReference:
    species: str
    value_cm2_s: float
    evidence_type: str
    split: str
    source_location: str


@dataclass(frozen=True)
class TransferObservation:
    family: str
    control: str
    observable: str
    value: str
    unit: str
    evidence_type: str
    split: str
    source_location: str


@dataclass(frozen=True)
class Krueger2024Evidence:
    calibration_metrics: tuple[ProfileTargetMetric, ...]
    boundary_fluxes: tuple[BoundaryFluxReference, ...]
    transfer_observations: tuple[TransferObservation, ...]


@dataclass(frozen=True)
class Jeon2022TrenchDepth:
    source_figure: str
    condition_family: str
    c4f8_fraction: float
    pulse_off_ms: float
    trench_width_nm: float
    depth_nm: float
    pixel_y: float
    axis_slope_nm_per_pixel: float
    axis_intercept_nm: float
    digitization_uncertainty_nm: float
    published_errorbar_semantics: str
    evidence_type: str
    split: str
    source_location: str


@dataclass(frozen=True)
class Jeon2022PlasmaControl:
    source_figure: str
    condition_family: str
    c4f8_fraction: float
    pulse_off_ms: float
    neutral_to_ion_flux_ratio: float
    pixel_y: float
    axis_slope_ratio_per_pixel: float
    axis_intercept_ratio: float
    digitization_uncertainty_ratio: float
    published_errorbar_semantics: str
    evidence_type: str
    role: str
    source_location: str


def _verified_csv_rows(path, expected_fields, expected_sha256, verify_checksum):
    path = Path(path)
    payload = path.read_bytes()
    if verify_checksum and sha256(payload).hexdigest() != expected_sha256:
        raise ValueError(f"checksum mismatch for experimental evidence: {path}")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_fields:
            raise ValueError(f"unexpected experimental-evidence schema: {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"experimental evidence is empty: {path}")
    return rows


def load_krueger_2024_evidence(directory, *, verify_checksum=True):
    """Load the Krüger 2024 calibration/transfer facts without conflating evidence types.

    Table-I boundary fluxes are HPEM outputs rather than measurements.  Transfer observations retain
    the source's evidence label, and MCFPM-only values are explicitly reference-only.
    """
    directory = Path(directory)
    metric_rows = _verified_csv_rows(
        directory / "base_case_metrics.csv",
        ["metric", "symbol", "value", "unit", "evidence_type", "split", "source_location"],
        KRUEGER_2024_SHA256["base_case_metrics.csv"], verify_checksum)
    flux_rows = _verified_csv_rows(
        directory / "base_case_boundary_fluxes.csv",
        ["species", "value", "unit", "evidence_type", "split", "source_location"],
        KRUEGER_2024_SHA256["base_case_boundary_fluxes.csv"], verify_checksum)
    transfer_rows = _verified_csv_rows(
        directory / "transfer_observations.csv",
        ["family", "control", "observable", "value", "unit", "evidence_type", "split",
         "source_location"],
        KRUEGER_2024_SHA256["transfer_observations.csv"], verify_checksum)

    metrics = tuple(ProfileTargetMetric(
        metric=row["metric"], symbol=row["symbol"], value=float(row["value"]), unit=row["unit"],
        evidence_type=row["evidence_type"], split=row["split"],
        source_location=row["source_location"])
        for row in metric_rows)
    fluxes = tuple(BoundaryFluxReference(
        species=row["species"], value_cm2_s=float(row["value"]),
        evidence_type=row["evidence_type"], split=row["split"],
        source_location=row["source_location"])
        for row in flux_rows)
    observations = tuple(TransferObservation(
        family=row["family"], control=row["control"], observable=row["observable"],
        value=row["value"], unit=row["unit"], evidence_type=row["evidence_type"],
        split=row["split"], source_location=row["source_location"])
        for row in transfer_rows)

    if any(item.evidence_type != "experiment" or item.split != "calibration" for item in metrics):
        raise ValueError("base-case metrics must be experimental calibration evidence")
    if any(item.evidence_type != "HPEM_simulation" for item in fluxes):
        raise ValueError("base-case boundary fluxes must remain labeled as HPEM simulation outputs")
    if any(item.split == "calibration" for item in observations):
        raise ValueError("transfer observations must not leak into the calibration split")
    return Krueger2024Evidence(metrics, fluxes, observations)


def load_jeon_2022_trench_depths(path, *, verify_checksum=True):
    """Load preregistered SiO2 depth-transfer targets digitized from Jeon et al. (2022).

    The plotted error bars are retained as semantically unspecified rather than being mislabeled as
    standard deviations. Pixel coordinates and axis maps are checked against every reported depth,
    keeping transcription error separate from the experiment's unreported measurement uncertainty.
    """
    expected = [
        "source_figure", "condition_family", "c4f8_fraction", "pulse_off_ms",
        "trench_width_nm", "depth_nm", "pixel_y", "axis_slope_nm_per_pixel",
        "axis_intercept_nm", "digitization_uncertainty_nm",
        "published_errorbar_semantics", "evidence_type", "split", "source_location",
    ]
    raw = _verified_csv_rows(path, expected, JEON_2022_DEPTH_SHA256, verify_checksum)
    rows = tuple(Jeon2022TrenchDepth(
        source_figure=row["source_figure"], condition_family=row["condition_family"],
        c4f8_fraction=float(row["c4f8_fraction"]), pulse_off_ms=float(row["pulse_off_ms"]),
        trench_width_nm=float(row["trench_width_nm"]), depth_nm=float(row["depth_nm"]),
        pixel_y=float(row["pixel_y"]),
        axis_slope_nm_per_pixel=float(row["axis_slope_nm_per_pixel"]),
        axis_intercept_nm=float(row["axis_intercept_nm"]),
        digitization_uncertainty_nm=float(row["digitization_uncertainty_nm"]),
        published_errorbar_semantics=row["published_errorbar_semantics"],
        evidence_type=row["evidence_type"], split=row["split"],
        source_location=row["source_location"]) for row in raw)
    replay_error = np.asarray([
        item.depth_nm - (
            item.axis_slope_nm_per_pixel * item.pixel_y + item.axis_intercept_nm)
        for item in rows])
    keys = {
        (item.source_figure, item.c4f8_fraction, item.pulse_off_ms, item.trench_width_nm)
        for item in rows}
    calibration = [item for item in rows if item.split == "calibration"]
    if (len(keys) != len(rows) or np.max(np.abs(replay_error)) > 0.051
            or any(item.depth_nm <= 0.0 or item.trench_width_nm <= 0.0
                   or item.digitization_uncertainty_nm <= 0.0 for item in rows)
            or any(item.evidence_type != "experiment_digitized" for item in rows)
            or any(item.published_errorbar_semantics != "not_specified" for item in rows)
            or any(item.split not in {"calibration", "held_out_transfer"} for item in rows)
            or any(item.source_figure != "4b" or item.c4f8_fraction != 0.2
                   or item.pulse_off_ms != 0.0 for item in calibration)):
        raise ValueError("Jeon 2022 evidence violates its digitization or split contract")
    return rows


def load_jeon_2022_plasma_controls(path, *, verify_checksum=True):
    """Load diagnostic-derived neutral/ion ratios without promoting them to direct measurements."""
    expected = [
        "source_figure", "condition_family", "c4f8_fraction", "pulse_off_ms",
        "neutral_to_ion_flux_ratio", "pixel_y", "axis_slope_ratio_per_pixel",
        "axis_intercept_ratio", "digitization_uncertainty_ratio",
        "published_errorbar_semantics", "evidence_type", "role", "source_location",
    ]
    raw = _verified_csv_rows(path, expected, JEON_2022_CONTROL_SHA256, verify_checksum)
    rows = tuple(Jeon2022PlasmaControl(
        source_figure=row["source_figure"], condition_family=row["condition_family"],
        c4f8_fraction=float(row["c4f8_fraction"]), pulse_off_ms=float(row["pulse_off_ms"]),
        neutral_to_ion_flux_ratio=float(row["neutral_to_ion_flux_ratio"]),
        pixel_y=float(row["pixel_y"]),
        axis_slope_ratio_per_pixel=float(row["axis_slope_ratio_per_pixel"]),
        axis_intercept_ratio=float(row["axis_intercept_ratio"]),
        digitization_uncertainty_ratio=float(row["digitization_uncertainty_ratio"]),
        published_errorbar_semantics=row["published_errorbar_semantics"],
        evidence_type=row["evidence_type"], role=row["role"],
        source_location=row["source_location"]) for row in raw)
    replay_error = np.asarray([
        item.neutral_to_ion_flux_ratio - (
            item.axis_slope_ratio_per_pixel * item.pixel_y + item.axis_intercept_ratio)
        for item in rows])
    keys = {
        (item.source_figure, item.c4f8_fraction, item.pulse_off_ms) for item in rows}
    if (len(keys) != len(rows) or np.max(np.abs(replay_error)) > 0.051
            or any(item.neutral_to_ion_flux_ratio <= 0.0
                   or item.digitization_uncertainty_ratio <= 0.0 for item in rows)
            or any(item.evidence_type != "diagnostic_derived_digitized" for item in rows)
            or any(item.role != "physical_boundary_input" for item in rows)
            or any(item.published_errorbar_semantics != "not_specified" for item in rows)):
        raise ValueError("Jeon 2022 plasma controls violate their digitization/evidence contract")
    return rows


def load_bosch_wafer_measurements(path, *, verify_checksum=True):
    """Load the nine-point wafer measurements from Zenodo record 17122442.

    The source README states that all measurement columns are in micrometres.  The two derived etch
    columns are checked against their defining measurement identities so corrupted or misinterpreted
    files fail loudly.
    """
    path = Path(path)
    payload = path.read_bytes()
    if verify_checksum and md5(payload).hexdigest() != ZENODO_17122442_9PT_MD5:
        raise ValueError(f"checksum mismatch for experimental dataset: {path}")

    expected = [
        "experiment_key", "lot_number", "wafer_number", "loc_id", "X", "Y",
        "preox_thickness", "postox_thickness", "stepheight", "oxide_etch", "si_etch",
    ]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected:
            raise ValueError(f"unexpected experimental-data schema: {reader.fieldnames}")
        rows = [
            BoschWaferMeasurement(
                experiment_key=row["experiment_key"] or None,
                lot_number=int(row["lot_number"]) if row["lot_number"] else None,
                wafer_number=int(row["wafer_number"]) if row["wafer_number"] else None,
                location_id=row["loc_id"],
                x_um=float(row["X"]),
                y_um=float(row["Y"]),
                pre_oxide_um=float(row["preox_thickness"]),
                post_oxide_um=float(row["postox_thickness"]),
                step_height_um=float(row["stepheight"]),
                oxide_etch_um=float(row["oxide_etch"]),
                silicon_etch_um=float(row["si_etch"]),
                post_oxide_original_um=float(row["postox_thickness"]),
                sampling_grid="9_point",
            )
            for row in reader
        ]

    if not rows:
        raise ValueError(f"experimental dataset is empty: {path}")
    oxide_error = np.array([
        row.oxide_etch_um - (row.pre_oxide_um - row.post_oxide_um) for row in rows
    ])
    silicon_error = np.array([
        row.silicon_etch_um - (row.step_height_um - row.oxide_etch_um) for row in rows
    ])
    if np.max(np.abs(oxide_error)) > 1e-10 or np.max(np.abs(silicon_error)) > 1e-10:
        raise ValueError("experimental dataset violates its etch-depth measurement identities")
    return rows


def load_bosch_wafer_measurements_89pt(path, *, verify_checksum=True):
    """Load the high-spatial-resolution wafer table from Zenodo record 17122442.

    This source has a distinct schema and defining identities from the nine-point table. Its
    ``postox_thickness`` column is the processed thickness used for derived values, while
    ``postox_thickness_nan`` preserves ``N/A`` at 157 originally unavailable measurements. Here
    ``oxide_etch = preox - postox`` and ``si_etch = stepheight - postox``. The distinction is kept
    explicit rather than forcing both source tables through one accidental convention.
    """
    path = Path(path)
    payload = path.read_bytes()
    if verify_checksum and md5(payload).hexdigest() != ZENODO_17122442_89PT_MD5:
        raise ValueError(f"checksum mismatch for experimental dataset: {path}")
    expected = [
        "experiment_key", "lot_number", "wafer_number", "X", "Y",
        "preox_thickness", "postox_thickness", "postox_thickness_nan",
        "stepheight", "oxide_etch", "si_etch",
    ]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected:
            raise ValueError(f"unexpected experimental-data schema: {reader.fieldnames}")
        rows = []
        for row in reader:
            original = row["postox_thickness_nan"]
            rows.append(BoschWaferMeasurement(
                experiment_key=row["experiment_key"] or None,
                lot_number=int(row["lot_number"]) if row["lot_number"] else None,
                wafer_number=int(row["wafer_number"]) if row["wafer_number"] else None,
                location_id=None, x_um=float(row["X"]), y_um=float(row["Y"]),
                pre_oxide_um=float(row["preox_thickness"]),
                post_oxide_um=float(row["postox_thickness"]),
                step_height_um=float(row["stepheight"]),
                oxide_etch_um=float(row["oxide_etch"]),
                silicon_etch_um=float(row["si_etch"]),
                post_oxide_original_um=(None if original == "N/A" else float(original)),
                sampling_grid="89_point"))
    if not rows:
        raise ValueError(f"experimental dataset is empty: {path}")
    oxide_error = np.array([
        row.oxide_etch_um - (row.pre_oxide_um - row.post_oxide_um) for row in rows])
    silicon_error = np.array([
        row.silicon_etch_um - (row.step_height_um - row.post_oxide_um) for row in rows])
    original_error = np.array([
        row.post_oxide_original_um - row.post_oxide_um
        for row in rows if row.post_oxide_original_um is not None])
    if (np.max(np.abs(oxide_error)) > 1e-8 or np.max(np.abs(silicon_error)) > 1e-8
            or np.max(np.abs(original_error)) > 1e-12):
        raise ValueError("high-resolution dataset violates its measurement identities")
    return rows
