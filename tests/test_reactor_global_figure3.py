import csv
import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DIGITIZED = (
    ROOT / "research_sources" / "digitized"
    / "lee_lieberman_1994_figure3_argon.csv")
MANIFEST = DIGITIZED.with_name(
    "lee_lieberman_1994_figure3_argon_manifest.json")
GRADER = ROOT / "scripts" / "grade_reactor_global_argon_figure3.py"


def _load_grader():
    specification = importlib.util.spec_from_file_location(
        "grade_reactor_global_argon_figure3", GRADER)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_figure3_digitization_is_visual_audited_and_monotonic():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    with DIGITIZED.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert manifest["visual_inspection"]["status"] == "passed"
    assert manifest["source"]["source_raster_committed"] is False
    assert manifest["row_count"] == len(rows) == 18
    pressures = np.asarray([float(row["pressure_mTorr"]) for row in rows])
    temperatures = np.asarray([
        float(row["electron_temperature_eV"]) for row in rows])
    assert np.all(np.diff(pressures) > 0.0)
    assert np.all(np.diff(temperatures) < 0.0)


def test_no_fit_figure3_reproduction_passes_frozen_gate():
    _, grade = _load_grader().run_reproduction()
    assert grade["coefficient_selection_target"] is None
    assert grade["passed"]
    assert len(grade["members"]) == 2
    assert all(member["passed"] for member in grade["members"])
    assert grade["worst_member_mean_absolute_percent_error"] <= 10.0
    assert grade["worst_member_maximum_absolute_percent_error"] <= 20.0
