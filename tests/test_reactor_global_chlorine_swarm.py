import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest

from petch.reactor_global.chlorine_swarm import (
    GONZALEZ_MAGANA_2018_PURE_CL2_SWARM_CSV_SHA256,
    GonzalezMaganaPureChlorineSwarmBoard,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "digitize_gonzalez_magana_2018_cl2_swarm.py"
DIGITIZED_CSV = (
    ROOT / "research_sources" / "digitized"
    / "gonzalez_magana_2018_pure_cl2_swarm.csv"
)
PACKAGE_CSV = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "gonzalez_magana_2018_pure_cl2_swarm.csv"
)
MANIFEST = (
    ROOT / "research_sources" / "digitized"
    / "gonzalez_magana_2018_pure_cl2_swarm_manifest.json"
)


def _digitizer():
    spec = importlib.util.spec_from_file_location("gonzalez_digitizer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_swarm_board_is_exact_hash_locked_script_replay():
    module = _digitizer()
    payload = module.csv_text()
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    assert DIGITIZED_CSV.read_text(encoding="utf-8") == payload
    assert PACKAGE_CSV.read_text(encoding="utf-8") == payload
    assert digest == GONZALEZ_MAGANA_2018_PURE_CL2_SWARM_CSV_SHA256
    assert json.loads(MANIFEST.read_text(encoding="utf-8")) == (
        module.manifest(digest))


def test_swarm_board_preserves_all_three_printed_pure_cl2_observables():
    board = GonzalezMaganaPureChlorineSwarmBoard.from_package_data()
    expected = {
        "electron_drift_velocity": 23,
        "effective_ionization_coefficient": 21,
        "density_normalized_longitudinal_diffusion": 8,
    }
    assert {
        observable: len(board.for_observable(observable))
        for observable in expected
    } == expected

    drift_100 = board.for_observable("electron_drift_velocity")[0]
    assert drift_100.reduced_field_Td == 100.0
    assert drift_100.value_si == pytest.approx(8.93e4)
    assert drift_100.relative_uncertainty_min == 0.02
    assert drift_100.relative_uncertainty_max == 0.02

    ionization_100 = board.for_observable(
        "effective_ionization_coefficient")[0]
    assert ionization_100.value_si == pytest.approx(-2.540e-21)
    assert ionization_100.relative_uncertainty_min == 0.05
    assert ionization_100.relative_uncertainty_max == 0.09

    diffusion_240 = board.for_observable(
        "density_normalized_longitudinal_diffusion")[0]
    assert diffusion_240.value_si == pytest.approx(1.38e24)
    assert diffusion_240.relative_uncertainty_min == 0.10
    assert diffusion_240.relative_uncertainty_max == 0.15


def test_swarm_board_is_a_collision_gate_not_a_reactor_or_depth_boundary():
    board = GonzalezMaganaPureChlorineSwarmBoard.from_package_data()
    assert board.supports_cross_section_validation
    assert not board.supports_reactor_state_prediction
    assert not board.supports_wafer_flux
    assert not board.supports_feature_depth
    assert all(item.supports_cross_section_validation for item in board.measurements)
    assert all(not item.supports_wafer_flux for item in board.measurements)
    assert all(not item.supports_feature_depth for item in board.measurements)


def test_source_range_conflict_and_visual_audit_remain_explicit():
    audit = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert audit["visual_audit"]["status"] == "passed"
    assert [page["pdf_page"] for page in audit["visual_audit"]["pages"]] == [
        6, 7, 8,
    ]
    assert audit["transcription"]["row_count"] == 52
    assert "100-420 Td" in audit[
        "source_conflicts_preserved"]["pure_cl2_range"]
    assert "through 460 Td" in audit[
        "source_conflicts_preserved"]["pure_cl2_range"]

    records = list(csv.DictReader(io.StringIO(PACKAGE_CSV.read_text())))
    drift_fields = [
        int(row["reduced_field_Td"])
        for row in records
        if row["observable"] == "electron_drift_velocity"
    ]
    assert drift_fields[-2:] == [440, 460]
