import csv
import hashlib
import io
import json
from pathlib import Path
import runpy

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "digitize_wise_1996_figure3.py"
CSV = (
    ROOT / "data" / "experimental" / "wise_1996_gec_icp"
    / "figure3_radial_measurements.csv"
)
MANIFEST = CSV.with_name("figure3_digitization_manifest.json")


def _module():
    return runpy.run_path(str(SCRIPT), run_name="wise_digitization_test")


def test_wise_figure3_committed_outputs_replay_exactly():
    module = _module()
    expected = module["csv_text"]()
    assert CSV.read_text(encoding="utf-8") == expected
    digest = hashlib.sha256(expected.encode()).hexdigest()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest == module["manifest"](digest)
    assert manifest["digitization"]["marker_count"] == 21
    assert manifest["digitization"]["axis_overlapping_marker_count"] == 3


def test_wise_figure3_physical_marker_board_is_consistent():
    rows = list(csv.DictReader(io.StringIO(_module()["csv_text"]())))
    by_observable = {
        name: [row for row in rows if row["observable"] == name]
        for name in ("electron_density", "electron_temperature", "plasma_potential")
    }
    assert all(len(values) == 7 for values in by_observable.values())
    assert [float(row["radial_distance_m"])
            for row in by_observable["electron_density"]] == pytest.approx(
        [0.0, 0.0127, 0.0254, 0.0381, 0.0508, 0.0635, 0.0762])
    electron_density = [
        float(row["value"]) for row in by_observable["electron_density"]]
    assert electron_density[0] == pytest.approx(1.4097561e17, rel=2.0e-8)
    assert electron_density[-1] == pytest.approx(2.0e16, rel=2.0e-9)
    temperature = [
        float(row["value"]) for row in by_observable["electron_temperature"]]
    assert max(temperature) == pytest.approx(3.88379205, rel=2.0e-8)
    potential = [
        float(row["value"]) for row in by_observable["plasma_potential"]]
    assert all(left > right for left, right in zip(potential, potential[1:]))
