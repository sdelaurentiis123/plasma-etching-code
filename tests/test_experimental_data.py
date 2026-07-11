from pathlib import Path

import numpy as np
import pytest

from petch.experimental_data import load_bosch_wafer_measurements


DATA = (
    Path(__file__).parents[1]
    / "data"
    / "experimental"
    / "zenodo_17122442"
    / "Si_Oxide_etch_9_points.csv"
)


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
