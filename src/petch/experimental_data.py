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
JEON_2022_ELECTRON_BIAS_SHA256 = (
    "f775e240914a8f874990596876e504e28fde331793453a405b45f4ac2945bae8")
JEONG_2023_DEPTH_SHA256 = (
    "27c170c6e2ccd2ef1c2c12f2fe641a8310c1757919d62d4e79475c39b271642d")
JEONG_2023_RADICAL_SHA256 = (
    "f923dc387070f72273817ef3302b585ae0d6005a6a73b5a5ae76aed654789d41")
JEONG_2023_XML_SHA256 = (
    "249045f4e77a47fd4e01fe77e7beb05b413d3468051f780600d4a8dbef86507c")
JEONG_2023_FIGURE6_SHA256 = (
    "3e4ea56418343dbc13bf3109e778a852181f2a473d169d35ad5ccfef0baf6d53")
JEONG_2023_FIGURE7_SHA256 = (
    "9c8acd0e9a7219ea5f99e097f7977fcb3ed490635fca3a0f9f69cb0a15a6508c")
DEBOER_2002_FIGURE9_SHA256 = (
    "ed0b72235887df70552356838e376540b26234a9084ce5c886fab45ed40d7b1b")
DEBOER_2002_PDF_SHA256 = (
    "45c245a9b19671f532945155dc16c3e00d35464eb8e49480a09f90a90498ff6c")
DEBOER_2002_FIGURE9_IMAGE_SHA256 = (
    "0f78ae30e5cc2e128f4fdb84217551fe350bd7696966c6ea40233f70a9a765c4")


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


@dataclass(frozen=True)
class Jeon2022ElectronBiasControl:
    source_figure: str
    condition_family: str
    c4f8_fraction: float
    pulse_off_ms: float
    electron_density_m3: float
    electron_pixel_y: float
    electron_axis_transform: str
    electron_axis_slope_per_pixel: float
    electron_axis_intercept: float
    electron_digitization_uncertainty_m3: float
    self_bias_magnitude_v: float
    self_bias_pixel_y: float
    self_bias_axis_slope_v_per_pixel: float
    self_bias_axis_intercept_v: float
    self_bias_digitization_uncertainty_v: float
    published_errorbar_semantics: str
    evidence_type: str
    role: str
    source_location: str


@dataclass(frozen=True)
class Jeon2022DimensionlessTarget:
    """Exposure-time-cancelling target with a digitization-only uncertainty interval."""

    observable: str
    source_figure: str
    condition_family: str
    c4f8_fraction: float
    pulse_off_ms: float
    trench_width_nm: float
    value: float
    digitization_lower: float
    digitization_upper: float
    split: str
    denominator: str
    cancellation_assumption: str


@dataclass(frozen=True)
class Jeong2023EtchDepth:
    """One fixed-duration experimental marker from Jeong et al. Figure 7."""

    source_figure: str
    control_mode: str
    trench_width_nm: float
    self_bias_magnitude_v: float
    electron_density_m3: float
    variation_percent: float
    etch_depth_nm: float
    etch_duration_s: float
    marker_pixel_x: float
    marker_pixel_y: float
    x_axis_slope_percent_per_pixel: float
    x_axis_intercept_percent: float
    y_axis_slope_nm_per_pixel: float
    y_axis_intercept_nm: float
    digitization_uncertainty_nm: float
    measurement_uncertainty_semantics: str
    evidence_type: str
    split: str
    role: str
    source_xml_sha256: str
    source_image_sha256: str
    source_location: str


@dataclass(frozen=True)
class Jeong2023RadicalDensity:
    """One nonexperimental radical-density bar from the source's plasma model."""

    source_figure: str
    species: str
    radical_class: str
    electron_density_m3: float
    particle_density_cm3: float
    bar_top_pixel_y: float
    axis_transform: str
    axis_slope_per_pixel: float
    axis_intercept: float
    digitization_uncertainty_log10: float
    evidence_type: str
    role: str
    source_xml_sha256: str
    source_image_sha256: str
    source_location: str


@dataclass(frozen=True)
class DeBoer2002Figure9Depth:
    """One replayable marker from the directly digitized de Boer Figure 9."""

    series_time_min: float
    mask_opening_um: float
    etch_depth_um: float
    marker_pixel_x: float
    marker_pixel_y: float
    x_axis_slope_um_per_pixel: float
    x_axis_intercept_um: float
    y_axis_slope_um_per_pixel: float
    y_axis_intercept_um: float
    digitization_uncertainty_x_um: float
    digitization_uncertainty_y_um: float
    measurement_uncertainty_um: float | None
    measurement_uncertainty_semantics: str
    evidence_type: str
    split: str
    role: str
    source_figure: str
    source_pdf_sha256: str
    source_image_sha256: str
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


def load_jeong_2023_etch_depths(path, *, verify_checksum=True):
    """Load the fixed-20-minute energy/flux transfer matrix from Jeong Figure 7.

    Only one marker is a magnitude-calibration anchor.  Energy response, ion-flux response,
    width transfer, and the 60 nm etch-stop behavior remain held out.  The publication does not
    define a statistical measurement uncertainty, so the stored 35 nm bound remains explicitly a
    digitization interval rather than a substitute error bar.
    """
    expected = [
        "source_figure", "control_mode", "trench_width_nm", "self_bias_magnitude_v",
        "electron_density_m3", "variation_percent", "etch_depth_nm", "etch_duration_s",
        "marker_pixel_x", "marker_pixel_y", "x_axis_slope_percent_per_pixel",
        "x_axis_intercept_percent", "y_axis_slope_nm_per_pixel", "y_axis_intercept_nm",
        "digitization_uncertainty_nm", "measurement_uncertainty_semantics", "evidence_type",
        "split", "role", "source_xml_sha256", "source_image_sha256", "source_location",
    ]
    raw = _verified_csv_rows(path, expected, JEONG_2023_DEPTH_SHA256, verify_checksum)
    rows = tuple(Jeong2023EtchDepth(
        source_figure=row["source_figure"], control_mode=row["control_mode"],
        trench_width_nm=float(row["trench_width_nm"]),
        self_bias_magnitude_v=float(row["self_bias_magnitude_v"]),
        electron_density_m3=float(row["electron_density_m3"]),
        variation_percent=float(row["variation_percent"]),
        etch_depth_nm=float(row["etch_depth_nm"]),
        etch_duration_s=float(row["etch_duration_s"]),
        marker_pixel_x=float(row["marker_pixel_x"]),
        marker_pixel_y=float(row["marker_pixel_y"]),
        x_axis_slope_percent_per_pixel=float(row["x_axis_slope_percent_per_pixel"]),
        x_axis_intercept_percent=float(row["x_axis_intercept_percent"]),
        y_axis_slope_nm_per_pixel=float(row["y_axis_slope_nm_per_pixel"]),
        y_axis_intercept_nm=float(row["y_axis_intercept_nm"]),
        digitization_uncertainty_nm=float(row["digitization_uncertainty_nm"]),
        measurement_uncertainty_semantics=row["measurement_uncertainty_semantics"],
        evidence_type=row["evidence_type"], split=row["split"], role=row["role"],
        source_xml_sha256=row["source_xml_sha256"],
        source_image_sha256=row["source_image_sha256"],
        source_location=row["source_location"]) for row in raw)
    x_replay = np.asarray([
        item.x_axis_slope_percent_per_pixel * item.marker_pixel_x
        + item.x_axis_intercept_percent for item in rows])
    y_replay = np.asarray([
        item.y_axis_slope_nm_per_pixel * item.marker_pixel_y
        + item.y_axis_intercept_nm for item in rows])
    calibration = [item for item in rows if item.split == "calibration"]
    keys = {(item.control_mode, item.trench_width_nm, item.self_bias_magnitude_v,
             item.electron_density_m3) for item in rows}
    energy = [item for item in rows if item.control_mode == "ion_energy"]
    flux = [item for item in rows if item.control_mode == "ion_flux"]
    if (len(rows) != 18 or len(keys) != len(rows)
            or np.max(np.abs(x_replay - np.asarray(
                [item.variation_percent for item in rows]))) > 2e-6
            or np.max(np.abs(y_replay - np.asarray(
                [item.etch_depth_nm for item in rows]))) > 1e-3
            or {item.trench_width_nm for item in rows} != {60.0, 100.0, 200.0}
            or len(energy) != 9 or len(flux) != 9 or len(calibration) != 1
            or (calibration[0].control_mode, calibration[0].trench_width_nm,
                calibration[0].self_bias_magnitude_v) != ("ion_energy", 200.0, 890.0)
            or calibration[0].role != "magnitude_calibration"
            or any(item.role != "held_out_prediction" for item in rows
                   if item.split == "held_out_transfer")
            or {item.self_bias_magnitude_v for item in energy} != {450.0, 890.0, 1270.0}
            or {item.electron_density_m3 for item in energy} != {2.0e15}
            or {item.self_bias_magnitude_v for item in flux} != {740.0}
            or {item.electron_density_m3 for item in flux} != {1.1e15, 1.9e15, 3.1e15}
            or any(item.etch_duration_s != 1200.0 or item.etch_depth_nm <= 0.0
                   or item.trench_width_nm <= 0.0 or item.digitization_uncertainty_nm != 35.0
                   or item.measurement_uncertainty_semantics != "not_reported"
                   or item.evidence_type != "experiment_digitized"
                   or item.source_xml_sha256 != JEONG_2023_XML_SHA256
                   or item.source_image_sha256 != JEONG_2023_FIGURE7_SHA256 for item in rows)):
        raise ValueError("Jeong 2023 depth evidence violates pixel replay or frozen split")
    return rows


def load_jeong_2023_radical_densities(path, *, verify_checksum=True):
    """Load Figure-6 plasma-model outputs without promoting them to measurements."""
    expected = [
        "source_figure", "species", "radical_class", "electron_density_m3",
        "particle_density_cm3", "bar_top_pixel_y", "axis_transform",
        "axis_slope_per_pixel", "axis_intercept", "digitization_uncertainty_log10",
        "evidence_type", "role", "source_xml_sha256", "source_image_sha256",
        "source_location",
    ]
    raw = _verified_csv_rows(path, expected, JEONG_2023_RADICAL_SHA256, verify_checksum)
    rows = tuple(Jeong2023RadicalDensity(
        source_figure=row["source_figure"], species=row["species"],
        radical_class=row["radical_class"],
        electron_density_m3=float(row["electron_density_m3"]),
        particle_density_cm3=float(row["particle_density_cm3"]),
        bar_top_pixel_y=float(row["bar_top_pixel_y"]),
        axis_transform=row["axis_transform"],
        axis_slope_per_pixel=float(row["axis_slope_per_pixel"]),
        axis_intercept=float(row["axis_intercept"]),
        digitization_uncertainty_log10=float(row["digitization_uncertainty_log10"]),
        evidence_type=row["evidence_type"], role=row["role"],
        source_xml_sha256=row["source_xml_sha256"],
        source_image_sha256=row["source_image_sha256"],
        source_location=row["source_location"]) for row in raw)
    replay_log10 = np.asarray([
        item.axis_slope_per_pixel * item.bar_top_pixel_y + item.axis_intercept
        for item in rows])
    expected_classes = {
        "C4F7": "heavy", "C3F6": "heavy", "C2F4": "heavy",
        "CF3": "light", "CF2": "light", "CF": "light",
    }
    keys = {(item.species, item.electron_density_m3) for item in rows}
    if (len(rows) != 18 or len(keys) != len(rows)
            or set(item.species for item in rows) != set(expected_classes)
            or {item.electron_density_m3 for item in rows} != {1.1e15, 1.9e15, 3.1e15}
            or np.max(np.abs(replay_log10 - np.log10(np.asarray(
        [item.particle_density_cm3 for item in rows])))) > 1e-12
            or any(item.radical_class != expected_classes[item.species]
                   or item.particle_density_cm3 <= 0.0
                   or item.axis_transform != "log10"
                   or item.digitization_uncertainty_log10 != 0.05
                   or item.evidence_type != "source_plasma_model_digitized"
                   or item.role != "nonexperimental_boundary_input"
                   or item.source_xml_sha256 != JEONG_2023_XML_SHA256
                   or item.source_image_sha256 != JEONG_2023_FIGURE6_SHA256 for item in rows)):
        raise ValueError("Jeong 2023 radical evidence violates replay or evidence class")
    return rows


def load_deboer_2002_figure9_depths(path, *, verify_checksum=True):
    """Load the direct de Boer Figure-9 pixels without substituting a transport-model curve."""
    expected = [
        "series_time_min", "mask_opening_um", "etch_depth_um",
        "marker_pixel_x", "marker_pixel_y", "x_axis_slope_um_per_pixel",
        "x_axis_intercept_um", "y_axis_slope_um_per_pixel", "y_axis_intercept_um",
        "digitization_uncertainty_x_um", "digitization_uncertainty_y_um",
        "measurement_uncertainty_um", "measurement_uncertainty_semantics",
        "evidence_type", "split", "role", "source_figure", "source_pdf_sha256",
        "source_image_sha256", "source_location",
    ]
    raw = _verified_csv_rows(
        path, expected, DEBOER_2002_FIGURE9_SHA256, verify_checksum)
    rows = tuple(DeBoer2002Figure9Depth(
        series_time_min=float(row["series_time_min"]),
        mask_opening_um=float(row["mask_opening_um"]),
        etch_depth_um=float(row["etch_depth_um"]),
        marker_pixel_x=float(row["marker_pixel_x"]),
        marker_pixel_y=float(row["marker_pixel_y"]),
        x_axis_slope_um_per_pixel=float(row["x_axis_slope_um_per_pixel"]),
        x_axis_intercept_um=float(row["x_axis_intercept_um"]),
        y_axis_slope_um_per_pixel=float(row["y_axis_slope_um_per_pixel"]),
        y_axis_intercept_um=float(row["y_axis_intercept_um"]),
        digitization_uncertainty_x_um=float(row["digitization_uncertainty_x_um"]),
        digitization_uncertainty_y_um=float(row["digitization_uncertainty_y_um"]),
        measurement_uncertainty_um=(
            None if row["measurement_uncertainty_um"] == ""
            else float(row["measurement_uncertainty_um"])),
        measurement_uncertainty_semantics=row["measurement_uncertainty_semantics"],
        evidence_type=row["evidence_type"], split=row["split"], role=row["role"],
        source_figure=row["source_figure"], source_pdf_sha256=row["source_pdf_sha256"],
        source_image_sha256=row["source_image_sha256"],
        source_location=row["source_location"]) for row in raw)
    x_replay = np.asarray([
        item.x_axis_slope_um_per_pixel * item.marker_pixel_x
        + item.x_axis_intercept_um for item in rows])
    y_replay = np.asarray([
        item.y_axis_slope_um_per_pixel * item.marker_pixel_y
        + item.y_axis_intercept_um for item in rows])
    keys = {(item.series_time_min, item.mask_opening_um) for item in rows}
    boundary = [item for item in rows if item.split == "boundary_input"]
    calibration = [item for item in rows if item.split == "calibration"]
    held_out = [item for item in rows if item.split == "held_out_transfer"]
    if (len(rows) != 16 or len(keys) != len(rows)
            # Reported coordinates are rounded to 1e-6 um while the linear axis
            # maps retain full precision, so replay closes to 0.002 um rather
            # than machine precision.  This remains over 100x tighter than the
            # declared digitization bounds (0.30/0.50 um).
            or np.max(np.abs(x_replay - np.asarray(
                [item.mask_opening_um for item in rows]))) > 2e-3
            or np.max(np.abs(y_replay - np.asarray(
                [item.etch_depth_um for item in rows]))) > 2e-3
            or {item.series_time_min for item in rows} != {5.5, 12.5, 25.0}
            or len(boundary) != 3 or len(calibration) != 1 or len(held_out) != 12
            or calibration[0].series_time_min != 12.5
            or any(item.role != "open_rate_anchor" for item in boundary)
            or calibration[0].role != "sticking_calibration"
            or any(item.role != "held_out_prediction" for item in held_out)
            or any(item.series_time_min <= 0.0 or item.mask_opening_um <= 0.0
                   or item.etch_depth_um <= 0.0
                   or item.digitization_uncertainty_x_um <= 0.0
                   or item.digitization_uncertainty_y_um <= 0.0 for item in rows)
            or any(item.measurement_uncertainty_um is not None
                   or item.measurement_uncertainty_semantics != "not_reported"
                   or item.evidence_type != "experiment_digitized"
                   or item.source_figure != "Fig. 9"
                   or item.source_pdf_sha256 != DEBOER_2002_PDF_SHA256
                   or item.source_image_sha256 != DEBOER_2002_FIGURE9_IMAGE_SHA256
                   for item in rows)):
        raise ValueError("de Boer Figure 9 evidence violates pixel replay or split provenance")
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


def load_jeon_2022_electron_bias_controls(path, *, verify_checksum=True):
    """Load digitized electron density and self-bias diagnostics as boundary evidence.

    The self-bias magnitude is an experimentally constrained sheath-energy scale, not an ion energy
    distribution. Electron densities can reproduce the paper's Bohm-flux estimate through
    :func:`jeon_2022_bohm_ion_flux_m2_s`; neither diagnostic supplies ion composition.
    """
    expected = [
        "source_figure", "condition_family", "c4f8_fraction", "pulse_off_ms",
        "electron_density_m3", "electron_pixel_y", "electron_axis_transform",
        "electron_axis_slope_per_pixel", "electron_axis_intercept",
        "electron_digitization_uncertainty_m3", "self_bias_magnitude_v",
        "self_bias_pixel_y", "self_bias_axis_slope_v_per_pixel",
        "self_bias_axis_intercept_v", "self_bias_digitization_uncertainty_v",
        "published_errorbar_semantics", "evidence_type", "role", "source_location",
    ]
    raw = _verified_csv_rows(
        path, expected, JEON_2022_ELECTRON_BIAS_SHA256, verify_checksum)
    rows = tuple(Jeon2022ElectronBiasControl(
        source_figure=row["source_figure"], condition_family=row["condition_family"],
        c4f8_fraction=float(row["c4f8_fraction"]), pulse_off_ms=float(row["pulse_off_ms"]),
        electron_density_m3=float(row["electron_density_m3"]),
        electron_pixel_y=float(row["electron_pixel_y"]),
        electron_axis_transform=row["electron_axis_transform"],
        electron_axis_slope_per_pixel=float(row["electron_axis_slope_per_pixel"]),
        electron_axis_intercept=float(row["electron_axis_intercept"]),
        electron_digitization_uncertainty_m3=float(
            row["electron_digitization_uncertainty_m3"]),
        self_bias_magnitude_v=float(row["self_bias_magnitude_v"]),
        self_bias_pixel_y=float(row["self_bias_pixel_y"]),
        self_bias_axis_slope_v_per_pixel=float(row["self_bias_axis_slope_v_per_pixel"]),
        self_bias_axis_intercept_v=float(row["self_bias_axis_intercept_v"]),
        self_bias_digitization_uncertainty_v=float(
            row["self_bias_digitization_uncertainty_v"]),
        published_errorbar_semantics=row["published_errorbar_semantics"],
        evidence_type=row["evidence_type"], role=row["role"],
        source_location=row["source_location"]) for row in raw)
    electron_replay = np.asarray([
        (10.0 ** (item.electron_axis_slope_per_pixel * item.electron_pixel_y
                  + item.electron_axis_intercept)
         if item.electron_axis_transform == "log10"
         else item.electron_axis_slope_per_pixel * item.electron_pixel_y
         + item.electron_axis_intercept)
        for item in rows])
    bias_replay = np.asarray([
        item.self_bias_axis_slope_v_per_pixel * item.self_bias_pixel_y
        + item.self_bias_axis_intercept_v for item in rows])
    keys = {(item.condition_family, item.c4f8_fraction, item.pulse_off_ms) for item in rows}
    if (len(keys) != len(rows)
            or np.max(np.abs(electron_replay / np.asarray(
                [item.electron_density_m3 for item in rows]) - 1.0)) > 3e-6
            or np.max(np.abs(bias_replay - np.asarray(
                [item.self_bias_magnitude_v for item in rows]))) > 0.051
            or any(item.electron_axis_transform not in {"linear", "log10"}
                   or item.electron_density_m3 <= 0.0
                   or item.electron_digitization_uncertainty_m3 <= 0.0
                   or item.self_bias_magnitude_v <= 0.0
                   or item.self_bias_digitization_uncertainty_v <= 0.0 for item in rows)
            or any(item.evidence_type != "diagnostic_digitized" for item in rows)
            or any(item.role != "physical_boundary_input" for item in rows)
            or any(item.published_errorbar_semantics != "not_specified" for item in rows)):
        raise ValueError("Jeon 2022 electron/bias controls violate their evidence contract")
    return rows


def jeon_2022_bohm_ion_flux_m2_s(
        control, *, electron_temperature_eV=3.0, ion_mass_u=39.948):
    """Replay Jeong et al.'s assumed Bohm ion flux from its measured electron density.

    This is a diagnostic-derived boundary estimate, not a direct ion-flux measurement. The defaults
    are the paper's assumed 3 eV electron temperature and Ar ion mass; ion composition was unmeasured.
    """
    if not isinstance(control, Jeon2022ElectronBiasControl):
        raise TypeError("Bohm replay requires a Jeon2022ElectronBiasControl")
    if electron_temperature_eV <= 0.0 or ion_mass_u <= 0.0:
        raise ValueError("electron temperature and ion mass must be positive")
    elementary_charge_c = 1.602176634e-19
    atomic_mass_kg = 1.66053906660e-27
    bohm_velocity_m_s = np.sqrt(
        electron_temperature_eV * elementary_charge_c / (ion_mass_u * atomic_mass_kg))
    return float(control.electron_density_m3 * bohm_velocity_m_s)


def jeon_2022_condition_wall_duration_s(
        reference_duration_s, wall_time_duty, exposure_basis):
    """Resolve wall duration without hiding Jeon's unreported pulse-exposure protocol.

    ``reference_duration_s`` is either a directly declared wall duration or a cumulative RF-on
    duration. Jeon et al. did not state which quantity was held fixed across their pulse sweep, so a
    nontrivial duty factor requires an explicit hypothesis. Continuous-wave runs are identical under
    both hypotheses.
    """
    duration = float(reference_duration_s)
    duty = float(wall_time_duty)
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("reference duration must be positive and finite")
    if not np.isfinite(duty) or not 0.0 < duty <= 1.0:
        raise ValueError("wall-time duty factor must lie in (0, 1]")
    if exposure_basis not in {"unspecified", "wall_time", "rf_on_time"}:
        raise ValueError("unknown pulse exposure basis")
    if duty < 1.0 and exposure_basis == "unspecified":
        raise ValueError(
            "pulsed Jeon runs require --pulse-exposure-basis wall_time or rf_on_time; "
            "the source does not report which exposure was held fixed")
    return duration / duty if exposure_basis == "rf_on_time" else duration


def _positive_ratio_interval(numerator, numerator_budget, denominator, denominator_budget):
    """Worst-case positive ratio interval; budgets are bounds, not standard deviations."""
    denominator_lower = denominator - denominator_budget
    if denominator_lower <= 0.0:
        raise ValueError("digitization budget does not bound a positive denominator")
    return (
        max(numerator - numerator_budget, 0.0) / (denominator + denominator_budget),
        (numerator + numerator_budget) / denominator_lower,
    )


def build_jeon_2022_dimensionless_targets(rows, *, pulse_exposure_basis=None):
    """Build ARDE-shape targets and, only by opt-in, cross-exposure pulse ratios.

    Width-shape targets divide each depth by the 200 nm trench depth under the same plasma
    condition, so their unknown exposure cancels.  Jeon et al. do not report whether pulse coupons
    shared wall-clock duration or cumulative RF-on duration.  Those choices differ by the inverse
    duty factor and therefore do *not* cancel in a pulsed/CW depth ratio.  Cross-exposure targets are
    consequently omitted unless the caller explicitly chooses ``common_wall_time`` or
    ``common_rf_on_time``; the unverified choice is retained in every returned target.

    Returned intervals propagate only the stored digitization budgets as worst-case bounds. They do
    not stand in for the publication's statistically unspecified experimental error bars.
    """
    rows = tuple(rows)
    if not rows or any(not isinstance(item, Jeon2022TrenchDepth) for item in rows):
        raise TypeError("dimensionless Jeon targets require Jeon2022TrenchDepth rows")
    allowed_exposure_basis = {None, "common_wall_time", "common_rf_on_time"}
    if pulse_exposure_basis not in allowed_exposure_basis:
        raise ValueError(
            "pulse_exposure_basis must be None, common_wall_time, or common_rf_on_time")
    by_condition = {}
    for item in rows:
        key = (item.source_figure, item.condition_family,
               item.c4f8_fraction, item.pulse_off_ms)
        by_width = by_condition.setdefault(key, {})
        if item.trench_width_nm in by_width:
            raise ValueError(f"duplicate Jeon condition/width row: {key}, {item.trench_width_nm}")
        by_width[item.trench_width_nm] = item

    targets = []
    for key, by_width in sorted(by_condition.items()):
        if 200.0 not in by_width:
            raise ValueError(f"Jeon condition lacks its 200 nm shape reference: {key}")
        reference = by_width[200.0]
        for width, item in sorted(by_width.items()):
            if item is reference:
                lower = upper = 1.0
            else:
                lower, upper = _positive_ratio_interval(
                    item.depth_nm, item.digitization_uncertainty_nm,
                    reference.depth_nm, reference.digitization_uncertainty_nm)
            targets.append(Jeon2022DimensionlessTarget(
                observable="width_shape_depth_over_200nm",
                source_figure=item.source_figure,
                condition_family=item.condition_family,
                c4f8_fraction=item.c4f8_fraction,
                pulse_off_ms=item.pulse_off_ms,
                trench_width_nm=width,
                value=item.depth_nm / reference.depth_nm,
                digitization_lower=lower,
                digitization_upper=upper,
                split=item.split,
                denominator="same_condition_200nm_trench_depth",
                cancellation_assumption="same_coupon_exposure_within_width_series"))

    if pulse_exposure_basis is None:
        return tuple(targets)

    pulse_assumption = (
        "explicit_common_wall_time_hypothesis_not_reported_by_source"
        if pulse_exposure_basis == "common_wall_time"
        else "explicit_common_rf_on_time_hypothesis_not_reported_by_source")
    for family in ("pulse_off_20pct", "pulse_off_80pct"):
        family_conditions = {
            key: by_width for key, by_width in by_condition.items() if key[1] == family}
        if not family_conditions:
            raise ValueError(f"Jeon pulse family is missing: {family}")
        for key, by_width in sorted(family_conditions.items()):
            source_figure, _, fraction, pulse_off = key
            if pulse_off == 0.0:
                continue
            cw_key = (source_figure, family, fraction, 0.0)
            if cw_key not in family_conditions:
                raise ValueError(f"Jeon pulse family lacks its continuous-wave reference: {key}")
            cw_by_width = family_conditions[cw_key]
            if set(by_width) != set(cw_by_width):
                raise ValueError(f"Jeon pulse/CW widths do not match: {key}")
            for width, item in sorted(by_width.items()):
                reference = cw_by_width[width]
                lower, upper = _positive_ratio_interval(
                    item.depth_nm, item.digitization_uncertainty_nm,
                    reference.depth_nm, reference.digitization_uncertainty_nm)
                targets.append(Jeon2022DimensionlessTarget(
                    observable="pulse_depth_over_cw",
                    source_figure=source_figure,
                    condition_family=family,
                    c4f8_fraction=fraction,
                    pulse_off_ms=pulse_off,
                    trench_width_nm=width,
                    value=item.depth_nm / reference.depth_nm,
                    digitization_lower=lower,
                    digitization_upper=upper,
                    split="held_out_transfer",
                    denominator="same_width_continuous_wave_depth",
                    cancellation_assumption=pulse_assumption))
    return tuple(targets)


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
