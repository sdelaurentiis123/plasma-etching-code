import json

import numpy as np
import pytest

from petch.nozawa_replay_3d import (
    NOZAWA_POLY_SI_MATERIAL_ID,
    _charging_surface_mesh_fingerprint_3d,
    load_nozawa_1995_replay_condition,
    make_nozawa_1995_open_area_geometry_3d,
    make_nozawa_1995_poisson_system_3d,
    make_nozawa_1995_replay_setup,
    main as nozawa_replay_main,
    measure_nozawa_1995_edge_notch_3d,
    run_nozawa_1995_replay,
)
from petch.threed import extract_mesh_3d


def test_nozawa_open_area_condition_is_checksum_bound_and_pad_families_refuse():
    condition = load_nozawa_1995_replay_condition("fig10_l06s06_04")

    assert condition.line_width_um == 0.6
    assert condition.space_width_um == 0.6
    assert np.isclose(condition.open_area_width_um, 1.484704283)
    assert np.isclose(condition.measured_notch_depth_um, 0.188826367)
    with pytest.raises(NotImplementedError, match="physical pad collector"):
        load_nozawa_1995_replay_condition("fig8_shared_01")


def test_nozawa_full_mirror_cell_preserves_stack_periodicity_and_four_lines():
    condition = load_nozawa_1995_replay_condition("fig10_l06s06_04")
    contract = make_nozawa_1995_open_area_geometry_3d(condition, dx_um=0.1)
    geometry = contract.geometry

    assert geometry.phi.shape == (58, 3, 39)
    assert abs(contract.realized_open_area_width_um - condition.open_area_width_um) <= 0.05
    assert np.array_equal(geometry.material_id[0], geometry.material_id[-1])
    assert np.array_equal(geometry.material_id[:, 0], geometry.material_id[:, -1])
    assert set(np.unique(geometry.material_id)) == {0, 1, 2, 3}
    assert np.count_nonzero(geometry.material_id == NOZAWA_POLY_SI_MATERIAL_ID) > 0
    assert contract.conductor_component_count == 4
    line_width = condition.line_width_um
    ordinary_space = condition.space_width_um
    centers = np.asarray(contract.line_centers_um)
    assert np.allclose(
        np.diff(centers) - line_width,
        (ordinary_space, contract.realized_open_area_width_um, ordinary_space))
    periodic_gap = (
        contract.realized_cell_width_um - centers[-1] + centers[0] - line_width)
    assert np.isclose(periodic_gap, ordinary_space)


def test_nozawa_charge_checkpoint_mesh_fingerprint_is_exact_and_condition_specific():
    first = make_nozawa_1995_open_area_geometry_3d(
        load_nozawa_1995_replay_condition("fig10_l06s06_04"), dx_um=0.1).geometry
    second = make_nozawa_1995_open_area_geometry_3d(
        load_nozawa_1995_replay_condition("fig10_l06s06_05"), dx_um=0.1).geometry

    assert (
        _charging_surface_mesh_fingerprint_3d(first)
        == _charging_surface_mesh_fingerprint_3d(first))
    assert (
        _charging_surface_mesh_fingerprint_3d(first)
        != _charging_surface_mesh_fingerprint_3d(second))


def test_nozawa_q1_system_makes_each_poly_line_equipotential_without_connecting_lines():
    condition = load_nozawa_1995_replay_condition("fig10_l06s06_04")
    geometry = make_nozawa_1995_open_area_geometry_3d(
        condition, dx_um=0.1).geometry
    system = make_nozawa_1995_poisson_system_3d(geometry)
    charge = np.zeros(system.shape)
    for component, total in zip(
            system.floating_conductor_ids, (1e-18, -2e-18, 3e-18, -4e-18)):
        selected = system.floating_conductor_node_ids == component
        charge[selected] = total / np.count_nonzero(selected)
    potential, diagnostics = system.solve(charge)

    assert system.periodic_axes == (0, 1)
    assert system.floating_conductor_ids == (1, 2, 3, 4)
    assert diagnostics.maximum_floating_conductor_voltage_spread_v < 1e-12
    conductor_voltage = [
        np.mean(potential[system.floating_conductor_node_ids == component])
        for component in system.floating_conductor_ids]
    assert len({round(value, 9) for value in conductor_voltage}) > 1


def test_nozawa_user_setup_is_public_engine_and_precommits_evidence_split():
    setup = make_nozawa_1995_replay_setup(n_position=4)

    assert setup.process.engine == "feature-charging-coevolution-3d-v1"
    assert setup.process.solver_options["periodic_lateral"]
    assert setup.process.solver_options["neutral_forward_scatter"].material_id == 3
    assert setup.process.solver_options["profile_motion_enabled"] is False
    assert setup.process.charging_options["scramble_mode"] == "fresh"
    assert setup.process.charging_options["compatible_q1_charge_state"]
    assert setup.preflight_manifest["calibration_split_frozen_before_run"]
    assert setup.preflight_manifest["protocol_commit_sha256"] == setup.protocol.commit_sha256
    assert setup.preflight_manifest["claim_blockers"]
    json.dumps(dict(setup.preflight_manifest))


def test_nozawa_experiment_uses_terminal_window_and_unattended_horizon():
    setup = make_nozawa_1995_replay_setup(
        mode="experiment", n_position=4, terminal_window_s=2e-6,
        trajectory_emergency_max_steps=65536)

    assert setup.process.charging_options["terminal_window_s"] == 2e-6
    assert setup.process.solver_options["trajectory_adaptive_horizon"]
    assert setup.process.solver_options["trajectory_emergency_max_steps"] == 65536
    assert setup.process.solver_options["profile_motion_enabled"] is True


def test_nozawa_decreasing_gain_tail_is_bounded_and_cannot_self_certify():
    setup = make_nozawa_1995_replay_setup(
        mode="charge_audit", n_position=4,
        charging_timestep_s=1.25e-7, maximum_charging_steps=20,
        charging_timestep_policy="decreasing_gain",
        stochastic_gain_exponent=0.75, stochastic_gain_offset_steps=64)
    options = setup.process.charging_options

    assert options["timestep_policy"] == "decreasing_gain"
    assert options["stochastic_gain_exponent"] == 0.75
    assert options["stochastic_gain_offset_steps"] == 64
    assert "terminal_window_s" not in options
    assert setup.preflight_manifest["numerics"]["terminal_window_s"] is None


def test_nozawa_target_observable_is_the_inner_foot_of_grating_edge_line_y():
    setup = make_nozawa_1995_replay_setup(
        "fig10_l06s06_05", n_position=4)
    observable = measure_nozawa_1995_edge_notch_3d(
        setup, setup.geometry_contract.geometry)

    assert observable.maximum_left_notch_depth_m == 0.0
    assert observable.maximum_right_notch_depth_m == 0.0
    assert np.isclose(observable.reference_left_boundary_m, 0.9e-6)
    assert np.isclose(observable.reference_right_boundary_m, 1.5e-6)


def test_nozawa_operational_smoke_runs_the_public_common_engine_end_to_end():
    setup = make_nozawa_1995_replay_setup(n_position=4, seed=1701)
    result = setup.process.run()
    step = result.steps[-1]

    assert step.charging.accepted_steps == 1
    assert not step.charging.converged  # one-nanosecond smoke never claims saturation
    assert step.diagnostics["profile_motion_enabled"] is False
    assert step.diagnostics["profile_duration_s"] == 0.0
    assert np.array_equal(result.geometry.phi, setup.geometry_contract.geometry.phi)
    assert np.array_equal(
        result.geometry.material_id, setup.geometry_contract.geometry.material_id)
    assert step.charging.history[-1]["charge_conservation_relative_error"] < 5e-13
    assert step.charge_remap.relative_charge_balance_error < 5e-13
    assert step.feature.neutral_forward_scatter is not None
    assert step.diagnostics["neutral_forward_scatter_particle_balance_error"] < 5e-13
    assert step.charging.final_step.poisson_before.maximum_floating_conductor_voltage_spread_v < 1e-12


def test_nozawa_smoke_artifact_runner_records_history_and_target_measurement(tmp_path):
    setup = make_nozawa_1995_replay_setup(n_position=4, seed=1701)
    result, summary = run_nozawa_1995_replay(setup, tmp_path / "run")

    assert result is not None
    assert summary["status"] == "pass"
    assert summary["predicted_target_notch_depth_um"] == 0.0
    assert (tmp_path / "run" / "charging_history.json").is_file()
    assert (tmp_path / "run" / "heartbeat.json").is_file()


def test_nozawa_restart_refuses_a_source_preflight_from_the_old_geometry(tmp_path, capsys):
    setup = make_nozawa_1995_replay_setup(
        "fig10_l06s06_04", mode="charge_audit", n_position=4)
    face_count = len(extract_mesh_3d(
        setup.geometry_contract.geometry.phi,
        setup.geometry_contract.geometry.dx)[1])
    source = tmp_path / "old-three-line"
    source.mkdir()
    checkpoint = source / "checkpoint.npz"
    np.savez_compressed(
        checkpoint,
        sigma_c_per_m2=np.zeros(face_count),
        condition_id=np.asarray(setup.condition.condition_id))
    (source / "preflight_manifest.json").write_text(json.dumps({
        "condition": {"condition_id": setup.condition.condition_id},
        "geometry": {"electrical_topology": "three_line_periodic_cell"},
    }))

    with pytest.raises(SystemExit) as info:
        nozawa_replay_main([
            "--condition", setup.condition.condition_id,
            "--mode", "charge_audit",
            "--n-position", "4",
            "--restart-checkpoint", str(checkpoint),
            "--restart-sampling-epoch", "0",
            "--output", str(tmp_path / "must-not-run"),
        ])

    assert info.value.code == 2
    assert "restart source geometry disagrees" in capsys.readouterr().err
