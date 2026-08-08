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
    assert molecular["schema"] == (
        "petch.malyshev_1998_eedf_source_replay.v2")
    assert molecular["rows"] == []
    assert molecular["all_requested_conditions_closed"] is False
    assert len(molecular["solver_failures"]) == 6
    assert all(
        "atomic-Cl electron collisions are mandatory"
        in failure["physical_interpretation"]
        for failure in molecular["solver_failures"]
    )
    assert molecular["atomic_momentum_payload_sha256"] is None

    for board in (atomic, hamilton):
        assert board["schema"] == (
            "petch.malyshev_1998_eedf_source_replay.v2")
        assert board["raw_collision_payload_sha256"] == (
            LEGACY_SIGLO_CL2_2013_SHA256)
        assert board["raw_collision_bytes_committed"] is False
        assert board["condition"]["source_frequency_Hz"] == 13.56e6
        assert board["energy_grid"]["family"] == (
            "threshold_aligned_piecewise_linear_v1")
        assert board["energy_grid"]["nominal_cell_count"] == 400
        assert board["energy_grid"]["breakpoints_eV"][:3] == [0.0, 0.5, 5.0]
        assert len(board["rows"]) == 6
        assert board["all_requested_conditions_closed"] is True
        assert board["solver_failures"] == []
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
    assert abs(held_out["relative_cl2_eq11_error_percentage_point"]) < 5.0
    assert held_out["within_reported_cl2_density_accuracy"] is True
    assert held_out[
        "chlorine_nuclei_equivalent_molecule_density_m3"
    ] == pytest.approx(
        held_out[
            "target_chlorine_nuclei_equivalent_molecule_density_m3"],
        rel=2.0e-10,
    )
    assert held_out["chlorine_particle_density_multiplier_vs_gauge_equivalent"] > 1.0


def test_isotropic_coulomb_board_is_bounded_negative_sensitivity():
    baseline = _load(
        "malyshev_1998_eedf_hamilton_atomic_cl_source_replay.json")
    coulomb = _load(
        "malyshev_1998_eedf_hamilton_atomic_cl_ee_source_replay.json")
    assert coulomb["raw_collision_payload_sha256"] == (
        LEGACY_SIGLO_CL2_2013_SHA256)
    assert coulomb["atomic_momentum_payload_sha256"] == (
        COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256)
    assert coulomb["hamilton_state_cross_sections_sha256"] == (
        HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256)
    assert coulomb["electron_electron_coulomb_model"] == (
        "isotropic_classical_debye")
    assert "electron-ion" in coulomb["comparison_boundaries"]["coulomb"]
    assert len(coulomb["rows"]) == 6
    assert max(
        row["maximum_normalized_residual"] for row in coulomb["rows"]
    ) < 2.0e-7
    assert all(
        12.0 < row["coulomb_logarithm"] < 14.0
        and row["coulomb_nonlinear_iterations"] >= 1
        and row["electron_growth_root_evaluations"] < 50
        and row["supports_feature_depth"] is False
        for row in coulomb["rows"]
    )
    baseline_held_out = next(
        row for row in baseline["rows"]
        if row["absorbed_fraction_sensitivity"] == 0.30
        and row["source_power_W"] == 500.0
    )
    coulomb_held_out = next(
        row for row in coulomb["rows"]
        if row["absorbed_fraction_sensitivity"] == 0.30
        and row["source_power_W"] == 500.0
    )
    assert coulomb_held_out["electron_density_percent_error"] < (
        baseline_held_out["electron_density_percent_error"])
    assert abs(coulomb_held_out["temperature_proxy_percent_error"]) > abs(
        baseline_held_out["temperature_proxy_percent_error"])
    assert coulomb_held_out["total_positive_ion_axial_flux_m2_s"] < (
        baseline_held_out["total_positive_ion_axial_flux_m2_s"])
    assert coulomb["supports_reactor_state_prediction"] is False
    assert coulomb["supports_wafer_flux"] is False
    assert coulomb["supports_feature_depth"] is False


def test_300W_conditioned_power_fraction_transfers_to_held_out_500W_density():
    receipt = _load("malyshev_1998_power_fraction_transfer.json")
    assert receipt["schema"] == (
        "petch.malyshev_1998_power_fraction_transfer.v2")
    assert receipt["collision_identity"][
        "raw_collision_payload_sha256"] == LEGACY_SIGLO_CL2_2013_SHA256
    assert receipt["collision_identity"][
        "atomic_momentum_payload_sha256"] == (
            COMSOL_64_ATOMIC_CL_MOMENTUM_SHA256)
    assert receipt["collision_identity"][
        "hamilton_state_cross_sections_sha256"] == (
            HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256)
    calibration = receipt["calibration"]
    assert calibration["root_converged"] is True
    assert calibration["fraction_search_bracket"] == [0.30, 0.50]
    assert calibration["fitted_absorbed_fraction"] == pytest.approx(
        0.36464679405196565, rel=2.0e-10)
    assert calibration["held_out_500W_used_for_selection"] is False
    assert calibration["temperature_used_for_selection"] is False
    assert calibration["dissociation_used_for_selection"] is False
    assert calibration["feature_depth_used_for_selection"] is False
    assert receipt["transfer"]["formal_pass_threshold"] is None

    training, held_out = receipt["rows"]
    assert training["validation_role"] == "calibration_training"
    assert held_out["validation_role"] == (
        "held_out_reactor_diagnostic_forecast")
    assert training["absorbed_fraction"] == held_out["absorbed_fraction"]
    assert abs(training["electron_density_percent_error"]) < 1.0e-4
    assert abs(held_out["electron_density_percent_error"]) < 2.0
    assert held_out["temperature_proxy_percent_error"] > 4.0
    assert held_out["relative_cl2_eq11_error_percent"] == pytest.approx(
        -12.282325573685593, rel=2.0e-9)
    assert held_out["within_reported_cl2_density_accuracy"] is True
    cl2_gate = receipt["independent_held_out_cl2_gate"]
    assert cl2_gate["passed_reported_accuracy"] is True
    assert cl2_gate["used_for_power_fraction_selection"] is False
    assert cl2_gate["reported_accuracy_percent"] == 25.0
    assert cl2_gate["relative_error_percent"] == pytest.approx(
        -12.282325573685593, rel=2.0e-9)
    assert "not a statistical sigma" in cl2_gate["accuracy_semantics"]
    assert held_out["total_positive_ion_axial_flux_m2_s"] == pytest.approx(
        5.77948763620014e19, rel=2.0e-9)
    assert receipt["supports_absorbed_power_measurement"] is False
    assert receipt["supports_reactor_state_prediction"] is False
    assert receipt["supports_wafer_flux"] is False
    assert receipt["supports_feature_depth"] is False


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
