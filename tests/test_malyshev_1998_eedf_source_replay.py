import json
from pathlib import Path

from petch.reactor_global import (
    COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256,
    LEGACY_SIGLO_CL2_2013_SHA256,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"


def _load(name):
    return json.loads((OUTPUT / name).read_text(encoding="utf-8"))


def test_source_replays_are_conservative_evidence_bounded_fixed_boards():
    molecular = _load("malyshev_1998_eedf_source_replay.json")
    atomic = _load("malyshev_1998_eedf_atomic_cl_source_replay.json")
    for board in (molecular, atomic):
        assert board["raw_collision_payload_sha256"] == (
            LEGACY_SIGLO_CL2_2013_SHA256)
        assert board["raw_collision_bytes_committed"] is False
        assert board["condition"]["source_frequency_Hz"] == 13.56e6
        assert board["energy_grid"]["family"] == (
            "threshold_aligned_piecewise_linear_v1")
        assert board["energy_grid"]["nominal_cell_count"] == 400
        assert board["energy_grid"]["breakpoints_eV"][:3] == [0.0, 0.5, 5.0]
        assert len(board["rows"]) == 6
        assert max(
            row["maximum_normalized_residual"] for row in board["rows"]
        ) < 1.0e-9
        assert board["supports_reactor_state_prediction"] is False
        assert board["supports_wafer_flux"] is False
        assert board["supports_feature_depth"] is False
        assert all(
            row["supports_feature_depth"] is False for row in board["rows"])
        assert board["preregistration"]["fraction_selected_in_this_run"] is None
        assert board["preregistration"]["feature_depth_used_for_selection"] is False
        assert "time-periodic" in board["comparison_boundaries"]["rf_heating"]

    assert molecular["atomic_momentum_payload_sha256"] is None
    assert atomic["atomic_momentum_payload_sha256"] == (
        COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256)
    assert 12.967633 in atomic["energy_grid"][
        "inserted_collision_thresholds_eV"]
