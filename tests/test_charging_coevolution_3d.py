import json

import numpy as np
import pytest
import petch.charging_coevolution_3d as charging_coevolution_3d

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.charged_surface_response_3d import GrazingSpecularIonReflection3D
from petch.charging_coevolution_3d import (
    ExperimentalObservableTolerance3D,
    ResolvedBiasSegment3D,
    SurfaceChargingSaturationError,
    _ser_candidate_acceptance,
    integrate_surface_charging_to_saturation_3d,
    physical_surface_patch_groups_3d,
    solve_charging_coevolution_3d,
)
from petch.charging_poisson import EPS0
from petch.charging_poisson_3d import NodalPoissonSystem3D, lump_triangle_sheet_charge_3d
from petch.feature_step_3d import FeatureGeometry3D
from petch.physical_sputtering import PhysicalSputterMechanism, PhysicalSputterParameters
from petch.sheath import ECHARGE
from petch.surface_kinetics import EnergeticYield, ParameterEvidence


def _species(name, charge_number, flux_m2_s, energy_eV=20.0):
    return SpeciesBoundaryState(
        name, charge_number, 40.0 if charge_number > 0 else 5.4858e-4,
        flux_m2_s, [[0.0, 0.0, np.sqrt(energy_eV)]], [1.0])


def _flat_problem(species):
    cell_shape = (1, 1, 10)
    spacing_m = np.full(3, 0.1e-6)
    fixed = np.zeros((2, 2, 11), dtype=bool)
    fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(np.ones(cell_shape), spacing_m, fixed)
    vertices = np.array([
        [0.0, 0.0, 0.0], [0.1, 0.0, 0.0],
        [0.1, 0.1, 0.0], [0.0, 0.1, 0.0]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    areas = np.full(2, 0.005)
    centroids = vertices[faces].mean(axis=1)
    normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    boundary = PlasmaBoundaryState(tuple(species), reference_plane_m=1e-6)
    arguments = dict(
        poisson_system=system, initial_sigma_c_per_m2=np.zeros(2),
        boundary=boundary, verts=vertices, faces=faces, areas=areas,
        face_centroids=centroids, face_gas_normals=normals,
        face_material_id=np.ones(2, dtype=int),
        source_bounds=(0.0, 0.1, 0.0, 0.1), source_z=1.0,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=0.1,
        patch_scales_m=(0.05e-6, 0.2e-6),
        potential_rate_tolerance_v_s=1e-6,
        timestep_s=1e-6, maximum_steps=1,
        mesh_length_unit_m=1e-6, n_position=16, seed=43,
        trajectory_fixed_dt=0.0025, trajectory_max_steps=2000,
        transport_device="cpu")
    return system, arguments


def _feature_plane():
    dx = 0.25
    shape = (4, 4, 8)
    z = np.arange(shape[2]) * dx
    phi = np.broadcast_to(0.95 - z, shape).copy()
    return FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), dx, 1e-6)


def _plane_poisson_system(geometry):
    fixed = np.zeros(geometry.phi.shape, dtype=bool)
    fixed[:, :, -1] = True
    phi_center = sum(
        geometry.phi[i:i + geometry.phi.shape[0] - 1,
                     j:j + geometry.phi.shape[1] - 1,
                     k:k + geometry.phi.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0
    return NodalPoissonSystem3D(
        np.where(phi_center > 0.0, 3.9, 1.0),
        geometry.dx * geometry.mesh_length_unit_m, fixed)


def _physical_sputter_mechanism():
    names = (
        "bulk_material_unit_density_m3", "sputter_yield",
        "emitted_product_mass_amu", "emission_angular_model", "emission_energy_model")
    evidence = {
        name: ParameterEvidence(
            "manufactured C3 co-evolution gate", "analytic",
            supports_prediction_within_declared_domain=True)
        for name in names}
    return PhysicalSputterMechanism(PhysicalSputterParameters(
        material_name="SiO2", material_inventory_name="SiO2_formula_unit",
        projectile_species=("Ar+",), bulk_material_unit_density_m3=2.2e28,
        sputter_yield=EnergeticYield(0.2, 20.0, 100.0),
        emitted_product_name="sputtered_SiO2_unit", emitted_product_mass_amu=60.084,
        emitted_material_units_per_particle=1.0,
        emission_angular_model="diffuse_cosine", emission_energy_model="thompson",
        emission_energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0},
        evidence=evidence))


def _balanced_feature_boundary(flux=2.2e21):
    return PlasmaBoundaryState((
        _species("Ar+", 1, flux, 100.0),
        _species("electron", -1, flux, 100.0),
    ), reference_plane_m=1.75e-6, provenance={"gate": "balanced planar C3"})


def _driver_options(maximum_steps=0):
    return dict(
        patch_scales_m=(0.25e-6, 1.0e-6),
        potential_rate_tolerance_v_s=1e-5,
        timestep_s=1e-9, maximum_steps=maximum_steps,
        current_balance_tolerance=0.08, timestep_policy="fixed")


def test_face_charge_is_authoritative_and_matches_the_q1_nodal_update():
    flux = 2.0e15
    system, arguments = _flat_problem((_species("ion", 1, flux),))
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    expected_sigma = ECHARGE * flux * arguments["timestep_s"]
    projected = lump_triangle_sheet_charge_3d(
        system.shape, arguments["verts"], arguments["faces"],
        result.sigma_c_per_m2, grid_origin=arguments["potential_origin"],
        grid_spacing=arguments["potential_spacing"], coordinate_length_unit_m=1e-6)
    expected_voltage = expected_sigma * 1e-6 / EPS0
    assert not result.converged
    assert result.accepted_steps == 1
    assert np.allclose(result.sigma_c_per_m2, expected_sigma, rtol=2e-15)
    assert np.allclose(result.charge_node_c, projected, rtol=0.0, atol=0.0)
    assert np.isclose(result.potential_v[:, :, 0].mean(), expected_voltage, rtol=2e-12)
    assert result.history[0]["face_to_node_update_relative_error"] < 3e-16
    assert result.history[0]["charge_conservation_relative_error"] < 3e-16
    assert len(result.patch_balance) == 2
    assert result.diagnostics["retained_node_max_relative_current_imbalance"] == 1.0


def test_equal_currents_pass_b1_b2_without_an_unnecessary_update():
    flux = 3.0e15
    system, arguments = _flat_problem((
        _species("ion", 1, flux), _species("electron", -1, flux)))
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    assert result.converged
    assert result.accepted_steps == 0
    assert result.rejected_steps == 0
    assert np.array_equal(result.sigma_c_per_m2, np.zeros(2))
    assert np.array_equal(result.charge_node_c, np.zeros(system.shape))
    assert result.history[0]["absolute_incident_charge_c"] > 0.0
    assert result.history[0]["charge_conservation_relative_error"] < 3e-16
    assert result.diagnostics["final_potential_rate_max_v_s"] < 1e-6
    assert all(item.maximum_relative_imbalance < 1e-6
               for item in result.patch_balance)


def test_ser_uses_the_same_conservative_ode_and_records_pseudo_time():
    flux = 3.0e15
    _system, arguments = _flat_problem((
        _species("ion", 1, flux), _species("electron", -1, flux)))
    arguments.update(
        maximum_steps=2, timestep_policy="ser", maximum_timestep_s=4e-6,
        stop_on_saturation=False)
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    assert result.converged
    assert result.accepted_steps == 2
    assert result.physical_time_s == 0.0
    assert np.isclose(result.pseudo_time_s, 2e-6)
    assert np.array_equal(result.sigma_c_per_m2, np.zeros(2))
    assert all(item["charge_conservation_relative_error"] < 3e-16
               for item in result.history)


def test_ser_safeguard_uses_absolute_ode_residual_not_b2_denominator():
    accepted, reason = _ser_candidate_acceptance(True, 0.9, 1.0, 0.005)
    rejected, rejection_reason = _ser_candidate_acceptance(True, 1.006, 1.0, 0.005)

    assert accepted
    assert reason is None
    assert not rejected
    assert rejection_reason == "absolute_current_residual_growth"


def test_fresh_scrambles_advance_reproducible_seed_epochs_in_fixed_physical_time():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(
        maximum_steps=2, stop_on_saturation=False, scramble_mode="fresh",
        sampling_seed_stride=101)
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    assert result.accepted_steps == 2
    assert [item["sampling_epoch"] for item in result.history] == [0, 1, 2]
    assert [item["sampling_seed"] for item in result.history] == [43, 144, 245]
    assert all(item["scramble_mode"] == "fresh" for item in result.history)
    assert result.diagnostics["scramble_mode"] == "fresh"
    assert result.diagnostics["initial_sampling_epoch"] == 0
    assert result.diagnostics["resume_sampling_epoch"] == 2
    assert all(item["charge_conservation_relative_error"] < 3e-16
               for item in result.history)


def test_progress_callback_receives_each_certified_replayable_state_read_only():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    observed = []

    def record(**progress):
        assert not progress["sigma_c_per_m2"].flags.writeable
        assert not progress["charge_node_c"].flags.writeable
        assert not progress["potential_v"].flags.writeable
        observed.append(dict(
            accepted=progress["accepted_steps"],
            epoch=progress["resume_sampling_epoch"],
            time=progress["physical_time_s"],
            sigma=progress["sigma_c_per_m2"].copy(),
            history=dict(progress["history_item"])))

    arguments.update(
        maximum_steps=2, stop_on_saturation=False, scramble_mode="fresh",
        sampling_seed_stride=101, progress_callback=record)
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    assert [item["accepted"] for item in observed] == [0, 1, 2]
    assert [item["epoch"] for item in observed] == [0, 1, 2]
    assert np.allclose(
        [item["time"] for item in observed], [0.0, 1e-6, 2e-6], rtol=0.0, atol=0.0)
    assert np.array_equal(observed[-1]["sigma"], result.sigma_c_per_m2)
    assert all("saturation_gates_satisfied" in item["history"] for item in observed)


def test_progress_callback_failure_returns_the_same_replayable_state():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))

    def fail(**_progress):
        raise OSError("manufactured durable-storage failure")

    arguments.update(progress_callback=fail)
    with pytest.raises(
            SurfaceChargingSaturationError, match="progress persistence failed") as info:
        integrate_surface_charging_to_saturation_3d(**arguments)

    assert info.value.accepted_steps == 0
    assert info.value.resume_sampling_epoch == 0
    assert np.array_equal(info.value.sigma_c_per_m2, np.zeros(2))


def test_fresh_scramble_restart_matches_an_uninterrupted_seed_sequence_bitwise():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(
        maximum_steps=2, stop_on_saturation=False, scramble_mode="fresh",
        sampling_seed_stride=101)
    uninterrupted = integrate_surface_charging_to_saturation_3d(**arguments)

    first_arguments = dict(arguments, maximum_steps=1)
    first = integrate_surface_charging_to_saturation_3d(**first_arguments)
    resumed_arguments = dict(
        first_arguments,
        initial_sigma_c_per_m2=first.sigma_c_per_m2,
        initial_sampling_epoch=first.diagnostics["resume_sampling_epoch"])
    resumed = integrate_surface_charging_to_saturation_3d(**resumed_arguments)

    assert [item["sampling_seed"] for item in uninterrupted.history] == [43, 144, 245]
    assert [item["sampling_seed"] for item in first.history] == [43, 144]
    assert [item["sampling_seed"] for item in resumed.history] == [144, 245]
    assert resumed.diagnostics["resume_sampling_epoch"] == 2
    assert np.array_equal(resumed.sigma_c_per_m2, uninterrupted.sigma_c_per_m2)
    assert np.array_equal(resumed.charge_node_c, uninterrupted.charge_node_c)
    assert np.array_equal(resumed.potential_v, uninterrupted.potential_v)


def test_failed_fresh_lookahead_checkpoint_preserves_updated_state_and_next_epoch(monkeypatch):
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(
        maximum_steps=2, stop_on_saturation=False, scramble_mode="fresh",
        sampling_seed_stride=101)
    uninterrupted = integrate_surface_charging_to_saturation_3d(**arguments)
    original = charging_coevolution_3d.advance_dielectric_charging_3d
    calls = 0

    def fail_second_evaluation(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("manufactured lookahead failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        charging_coevolution_3d, "advance_dielectric_charging_3d", fail_second_evaluation)
    with pytest.raises(SurfaceChargingSaturationError, match="manufactured lookahead") as info:
        integrate_surface_charging_to_saturation_3d(**arguments)
    error = info.value

    assert error.accepted_steps == 0
    assert error.state_updates == 1
    assert error.resume_sampling_epoch == 1
    assert error.physical_time_s == arguments["timestep_s"]
    monkeypatch.setattr(charging_coevolution_3d, "advance_dielectric_charging_3d", original)
    resumed = integrate_surface_charging_to_saturation_3d(**dict(
        arguments, maximum_steps=1,
        initial_sigma_c_per_m2=error.sigma_c_per_m2,
        initial_sampling_epoch=error.resume_sampling_epoch))
    assert np.array_equal(resumed.sigma_c_per_m2, uninterrupted.sigma_c_per_m2)
    assert np.array_equal(resumed.charge_node_c, uninterrupted.charge_node_c)
    assert np.array_equal(resumed.potential_v, uninterrupted.potential_v)


def test_fresh_scrambles_refuse_ser_and_stale_adjoint_proposals():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(scramble_mode="fresh", timestep_policy="ser")
    with pytest.raises(ValueError, match="fresh-scramble"):
        integrate_surface_charging_to_saturation_3d(**arguments)

    arguments.update(timestep_policy="fixed", adjoint_proposals={"ion": object()})
    with pytest.raises(ValueError, match="fresh-scramble"):
        integrate_surface_charging_to_saturation_3d(**arguments)


def test_saturation_failure_checkpoint_carries_its_exact_clocks():
    error = SurfaceChargingSaturationError(
        "manufactured failure", np.zeros(2), (), 3, 1, 2.5e-6, 4.0e-6)

    assert error.accepted_steps == 3
    assert error.rejected_steps == 1
    assert error.physical_time_s == 2.5e-6
    assert error.pseudo_time_s == 4.0e-6
    assert error.state_updates == 3
    assert error.resume_sampling_epoch == 0


def test_physical_patch_groups_separate_wall_and_floor_at_a_shared_corner():
    centroid = np.array([
        [0.01, 0.02, 0.01], [0.01, 0.02, 0.03],
        [0.02, 0.02, 0.01], [0.04, 0.02, 0.01]])
    normal = np.array([
        [1.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
    group = physical_surface_patch_groups_3d(
        centroid, normal, np.ones(4, dtype=int), 0.1e-6,
        mesh_length_unit_m=1e-6)

    assert group[0] == group[1]
    assert group[2] == group[3]
    assert group[0] != group[2]


def test_quasi_static_driver_reuses_exact_transport_then_remaps_signed_charge():
    geometry = _feature_plane()
    boundary = _balanced_feature_boundary()
    response = GrazingSpecularIonReflection3D.literature_bounded_sensitivity(1, "Ar+")
    result = solve_charging_coevolution_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        _physical_sputter_mechanism(), charging_system_builder=_plane_poisson_system,
        etchable_material_ids=(1,), duration_s=0.1, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=_driver_options(), charged_surface_response=response,
        n_position=256, seed=47, trajectory_fixed_dt=0.005,
        trajectory_max_steps=1000, reinitialize=False, transport_device="cpu")

    step = result.steps[0]
    assert step.charging.converged
    assert step.charging.accepted_steps == 0
    assert step.feature.transport is step.charging.final_step.transport
    assert not np.array_equal(result.geometry.phi, geometry.phi)
    assert np.array_equal(result.sigma_c_per_m2, np.zeros_like(result.sigma_c_per_m2))
    assert step.charge_remap.relative_charge_balance_error == 0.0
    assert step.diagnostics["retained_positive_charge_c"] == 0.0
    assert step.diagnostics["retained_negative_charge_c"] == 0.0
    assert len(step.charging.patch_balance) == 2
    assert "retained_node_rms_relative_current_imbalance" in step.diagnostics
    assert "retained_node_max_relative_current_imbalance" in step.diagnostics
    assert result.run_manifest["mode"] == "quasi_static"
    assert result.run_manifest["convergence_contract_revision"] == "CCA-2026-07-13-R2"
    assert result.run_manifest["charged_surface_response"]["parameters"][
        "grazing_reflection_probability"] == 0.95
    json.dumps(dict(result.run_manifest))


def test_waveform_mode_advances_each_physical_segment_without_saturation_assumption():
    geometry = _feature_plane()
    boundary = _balanced_feature_boundary()
    high_ion = PlasmaBoundaryState((
        _species("Ar+", 1, 2.2e21, 100.0),
        _species("electron", -1, 1.1e21, 100.0),
    ), boundary.reference_plane_m, provenance={"bias_phase": "ion_rich"})
    high_electron = PlasmaBoundaryState((
        _species("Ar+", 1, 1.1e21, 100.0),
        _species("electron", -1, 2.2e21, 100.0),
    ), boundary.reference_plane_m, provenance={"bias_phase": "electron_rich"})
    waveform = (
        ResolvedBiasSegment3D(1e-9, high_ion),
        ResolvedBiasSegment3D(1e-9, high_electron),
    )
    result = solve_charging_coevolution_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        _physical_sputter_mechanism(), charging_system_builder=_plane_poisson_system,
        etchable_material_ids=(1,), duration_s=2e-9, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=_driver_options(), bias_mode="waveform_resolved",
        bias_waveform=waveform, n_position=128, seed=53,
        trajectory_fixed_dt=0.005, trajectory_max_steps=1000,
        reinitialize=False, transport_device="cpu")

    assert result.run_manifest["mode"] == "waveform_resolved"
    assert [item["duration_s"] for item in result.run_manifest["waveform"]] == [1e-9, 1e-9]
    assert all(step.charging.accepted_steps == 1 for step in result.steps)
    assert all(not step.diagnostics["saturation_required"] for step in result.steps)
    assert all(np.isclose(step.charging.physical_time_s, 1e-9) for step in result.steps)
    assert np.allclose(result.steps[0].diagnostics["patch_max_relative_imbalance"], 0.5)
    assert np.allclose(result.steps[1].diagnostics["patch_max_relative_imbalance"], 1.0)
    assert np.allclose(
        result.steps[0].diagnostics["patch_symmetric_max_relative_imbalance"], 1.0 / 3.0)
    assert "not experimentally validated" in " ".join(result.validity.known_limitations)


def test_quasi_static_refuses_waveform_and_b3_refuses_unanchored_claims():
    geometry = _feature_plane()
    boundary = _balanced_feature_boundary()
    common = dict(
        geometry=geometry, boundary=boundary,
        species_role={"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism=_physical_sputter_mechanism(),
        charging_system_builder=_plane_poisson_system, etchable_material_ids=(1,),
        duration_s=0.0, n_steps=1, source_bounds=(0.0, 0.75, 0.0, 0.75),
        source_z=1.75, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=geometry.dx, charging_options=_driver_options(),
        n_position=16, trajectory_fixed_dt=0.005, trajectory_max_steps=1000,
        reinitialize=False, transport_device="cpu")
    with pytest.raises(ValueError, match="refuses pulsed bias"):
        solve_charging_coevolution_3d(
            **common, bias_waveform=(ResolvedBiasSegment3D(1e-9, boundary),))
    with pytest.raises(ValueError, match="requires at least one"):
        solve_charging_coevolution_3d(**common, experimental_claim=True)
    with pytest.raises(ValueError, match="must not exceed"):
        ExperimentalObservableTolerance3D("notch_depth", 2.0, 1.0)
    too_small = ExperimentalObservableTolerance3D(
        "notch_depth", 0.8, 1.0, feature_extent_m=0.1e-6)
    with pytest.raises(ValueError, match="patch scale"):
        solve_charging_coevolution_3d(
            **common, experimental_claim=True, observable_tolerances=(too_small,))
    anchored = ExperimentalObservableTolerance3D(
        "notch_depth", 0.8, 1.0, feature_extent_m=0.3e-6)
    claimed = solve_charging_coevolution_3d(
        **common, experimental_claim=True, observable_tolerances=(anchored,))
    assert claimed.run_manifest["observable_tolerances"][0]["tolerance"] == 0.8
