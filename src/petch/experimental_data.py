"""Load public experimental validation data with provenance and schema checks.

Experimental data are deliberately kept separate from solver parameters.  Loading a dataset does
not calibrate a model; it only produces units-explicit observations suitable for a validation or
calibration harness.
"""

from dataclasses import dataclass
from hashlib import md5
import csv
from pathlib import Path
from typing import Optional

import numpy as np


ZENODO_17122442_9PT_MD5 = "78515caf25e29e558e1859b92f8a4827"


@dataclass(frozen=True)
class BoschWaferMeasurement:
    experiment_key: Optional[str]
    lot_number: Optional[int]
    wafer_number: Optional[int]
    location_id: str
    x_um: float
    y_um: float
    pre_oxide_um: float
    post_oxide_um: float
    step_height_um: float
    oxide_etch_um: float
    silicon_etch_um: float


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
