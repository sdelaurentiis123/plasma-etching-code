import json
from pathlib import Path

import pytest

from petch.reactor_global import (
    COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256,
    HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256,
    LEGACY_SIGLO_CL2_2013_SHA256,
)
from scripts.audit_malyshev_1998_detachment_importance import audit


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"


def _load(name):
    return json.loads((OUTPUT / name).read_text(encoding="utf-8"))


def test_source_replays_are_conservative_evidence_bounded_fixed_boards():
    molecular = _load("malyshev_1998_eedf_source_replay.json")
    atomic = _load("malyshev_1998_eedf_atomic_cl_source_replay.json")
    hamilton = _load(
        "malyshev_1998_eedf_hamilton_atomic_cl_source_replay.json")
    for board in (molecular, atomic, hamilton):
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
    assert hamilton["atomic_momentum_payload_sha256"] == (
        COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256)
    assert hamilton["model_variant"] == (
        "legacy_siglo_hamilton_plus_comsol_nist_atomic_cl")
    assert hamilton["hamilton_state_cross_sections_sha256"] == (
        HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256)
    assert 12.967633 in atomic["energy_grid"][
        "inserted_collision_thresholds_eV"]
    for _, threshold in (
        ("a", 3.252), ("A", 4.348), ("b", 6.498), ("B", 7.537),
        ("C", 7.790), ("c", 7.257), ("D", 8.228), ("e", 9.219),
    ):
        assert threshold in hamilton["energy_grid"][
            "inserted_collision_thresholds_eV"]
    held_out = next(
        row for row in hamilton["rows"]
        if row["absorbed_fraction_sensitivity"] == 0.30
        and row["source_power_W"] == 500.0
    )
    assert abs(held_out["temperature_proxy_percent_error"]) < 2.0
    assert held_out["electron_density_percent_error"] < -10.0
    assert held_out["relative_cl2_proxy_error_percentage_point"] < -15.0


def test_detachment_timescale_receipt_is_low_leverage_and_reproducible():
    committed = _load("malyshev_1998_detachment_importance.json")
    assert committed == audit()
    assert committed["maximum_detachment_to_neutralization_loss_ratio"] < 0.04
    assert committed["coupled_response_is_formal_bound"] is False
    assert committed["feature_depth_used"] is False
    for row in committed["rows"]:
        assert row["double_detachment_rate_coefficient_m3_s"] < (
            2.0e-4 * row["single_detachment_rate_coefficient_m3_s"])
        expected = (
            row["electron_detachment_loss_frequency_s_inv"]
            / row["mutual_neutralization_loss_frequency_s_inv"]
        )
        assert row["detachment_to_neutralization_loss_ratio"] == (
            pytest.approx(expected, rel=2.0e-15))


def test_hamilton_held_out_grid_convergence_is_below_half_per_mille():
    receipt = _load("malyshev_1998_eedf_hamilton_grid_convergence.json")
    assert receipt["coarse_grid"]["actual_cell_count"] == 415
    assert receipt["fine_grid"]["actual_cell_count"] == 813
    assert receipt["collision_identity"][
        "hamilton_state_cross_sections_sha256"
    ] == HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256
    assert receipt["condition"] == {
        "absorbed_fraction_sensitivity": 0.30,
        "source_power_W": 500.0,
        "validation_role": "held_out_reactor_diagnostic",
    }
    assert receipt["maximum_change_metric"] == "electron_density_m3"
    assert receipt["maximum_absolute_relative_change"] < 5.0e-4
    assert receipt["numerical_convergence_passed"] is True
    assert receipt["feature_depth_used"] is False
    assert receipt["supports_feature_depth"] is False
