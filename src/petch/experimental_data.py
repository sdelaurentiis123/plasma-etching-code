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
