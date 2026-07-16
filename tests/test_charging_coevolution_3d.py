import json
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
import petch.charging_coevolution_3d as charging_coevolution_3d

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.charged_surface_response_3d import GrazingSpecularIonReflection3D
from petch.charging_coupled_3d import DielectricChargingStep3DResult
from petch.charging_coevolution_3d import (
    ExperimentalObservableTolerance3D,
    ResolvedBiasSegment3D,
    SurfaceChargingSaturationError,
    _ser_candidate_acceptance,
    _terminal_window_current_mean,
    integrate_surface_charging_to_saturation_3d,
    physical_surface_patch_groups_3d,
    propose_compatible_q1_pseudo_time_step_3d,
    solve_charging_coevolution_3d,
)
from petch.charging_poisson import EPS0
from petch.charging_checkpoint_3d import PhysicalChargingCheckpoint3D
from petch.charging_poisson_3d import (
    CompatibleQ1SurfaceChargeProjector3D,
    NodalPoissonSystem3D,
    lump_triangle_sheet_charge_3d,
)
from petch.charging_stationarity_3d import ProfileChargingStationarityContract3D
from petch.feature_step_3d import FeatureGeometry3D
from petch.physical_sputtering import PhysicalSputterMechanism, PhysicalSputterParameters
from petch.physical_api import (
    COMMON_CHARGING_ENGINE,
    COMMON_CHARGING_ENSEMBLE_ENGINE,
    PhysicalChargingEnsembleProcess,
    PhysicalChargingProcess,
)
from petch.sheath import ECHARGE
from petch.surface_kinetics import EnergeticYield, ParameterEvidence
from petch.twist_campaign_3d import TwistEnsembleRefinementContract3D


def _species(name, charge_number, flux_m2_s, energy_eV=20.0):
    return SpeciesBoundaryState(
        name, charge_number, 40.0 if charge_number > 0 else 5.4858e-4,
        flux_m2_s, [[0.0, 0.0, np.sqrt(energy_eV)]], [1.0])


def _flat_problem(species, *, periodic=False):
    cell_shape = (1, 1, 10)
    spacing_m = np.full(3, 0.1e-6)
    fixed = np.zeros((2, 2, 11), dtype=bool)
    fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(
        np.ones(cell_shape), spacing_m, fixed,
        periodic_axes=((0, 1) if periodic else ()))
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
        periodic_lateral=periodic, transport_device="cpu")
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


def test_compatible_q1_charge_state_preserves_an_injective_physical_step():
    flux = 2.0e15
    system, arguments = _flat_problem((_species("ion", 1, flux),))
    legacy = integrate_surface_charging_to_saturation_3d(**arguments)
    compatible = integrate_surface_charging_to_saturation_3d(
        **dict(arguments, compatible_q1_charge_state=True))

    assert np.allclose(
        compatible.sigma_c_per_m2, legacy.sigma_c_per_m2,
        rtol=2e-15, atol=2e-30)
    assert np.array_equal(compatible.charge_node_c, legacy.charge_node_c)
    assert np.array_equal(compatible.potential_v, legacy.potential_v)
    assert compatible.diagnostics["compatible_q1_charge_state"]
    assert compatible.diagnostics["q1_face_coupling_rank"] == 2
    assert compatible.diagnostics["q1_face_coupling_nullity"] == 0
    assert compatible.diagnostics["initial_unresolved_face_charge_fraction"] == 0.0
    assert compatible.diagnostics["maximum_unresolved_face_current_fraction"] < 2e-15
    assert compatible.history[0]["unresolved_face_current_projection_l1_c"] < 2e-30
    assert compatible.history[0]["face_to_node_update_relative_error"] < 3e-16
    assert np.allclose(
        compatible.history[0]["patch_q1_resolved_max_ion_normalized_imbalance"],
        compatible.history[0]["patch_max_relative_imbalance"],
        rtol=2e-15, atol=2e-15)
    assert np.allclose(
        compatible.history[0]["patch_q1_unresolved_max_ion_normalized_imbalance"],
        0.0, rtol=0.0, atol=2e-15)
    projected = lump_triangle_sheet_charge_3d(
        system.shape, arguments["verts"], arguments["faces"],
        compatible.sigma_c_per_m2, grid_origin=arguments["potential_origin"],
        grid_spacing=arguments["potential_spacing"], coordinate_length_unit_m=1e-6)
    assert np.array_equal(compatible.charge_node_c, projected)


def test_compatible_charge_state_pools_floating_conductor_faces_without_field_loss():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    fixed = np.zeros((2, 2, 11), dtype=bool)
    fixed[:, :, -1] = True
    conductor = np.zeros(fixed.shape, dtype=int)
    conductor[:, :, 0] = 4
    arguments["poisson_system"] = NodalPoissonSystem3D(
        np.ones((1, 1, 10)), np.full(3, 0.1e-6), fixed,
        floating_conductor_node_ids=conductor)
    result = integrate_surface_charging_to_saturation_3d(
        **dict(arguments, compatible_q1_charge_state=True))

    assert result.diagnostics["compatible_q1_charge_state"]
    assert result.diagnostics["q1_face_coupling_nullity"] == 1
    assert np.allclose(
        result.sigma_c_per_m2, np.mean(result.sigma_c_per_m2),
        rtol=2e-14, atol=2e-30)
    assert result.history[0]["face_to_node_update_relative_error"] < 5e-13
    assert result.history[0]["unresolved_face_current_fraction"] >= 0.0
    assert (
        result.final_step.poisson_after.maximum_floating_conductor_voltage_spread_v
        < 2e-12)


def test_periodic_compatible_charge_step_has_no_field_seam():
    system, arguments = _flat_problem(
        (_species("ion", 1, 2e15),), periodic=True)
    result = integrate_surface_charging_to_saturation_3d(
        **dict(arguments, compatible_q1_charge_state=True))

    assert result.diagnostics["poisson_periodic_axes"] == (0, 1)
    assert result.diagnostics["poisson_independent_node_shape"] == (1, 1, 11)
    assert np.array_equal(result.potential_v[0, :, :], result.potential_v[-1, :, :])
    assert np.array_equal(result.potential_v[:, 0, :], result.potential_v[:, -1, :])
    assert np.array_equal(
        result.charge_node_c, system.canonicalize_charge(result.charge_node_c))
    assert result.history[0]["face_to_node_update_relative_error"] < 3e-15


def test_compatible_q1_charge_state_filters_a_manufactured_grid_invisible_current(
        monkeypatch):
    vertices = np.asarray([
        [float(i), float(j), 0.0]
        for i in range(4) for j in range(4)])

    def vertex(i, j):
        return 4 * i + j

    faces = []
    for i in range(3):
        for j in range(3):
            faces.extend((
                [vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1)],
                [vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1)]))
    faces = np.asarray(faces, dtype=int)
    areas = np.full(len(faces), 0.5)
    fixed = np.zeros((4, 4, 2), dtype=bool)
    fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(np.ones((3, 3, 1)), 1e-6, fixed)
    projector = CompatibleQ1SurfaceChargeProjector3D.from_triangles(
        system.shape, vertices, faces, coordinate_length_unit_m=1e-6)
    trial = np.random.default_rng(916).normal(size=len(faces))
    null_current = trial - projector.project_face_charge(trial)
    null_current *= 2e-15 / np.max(np.abs(null_current))
    resolved_current = projector.project_face_charge(
        np.random.default_rng(917).normal(size=len(faces)))
    resolved_current *= 2e-16 / np.max(np.abs(resolved_current))
    face_net_current = null_current + resolved_current
    positive_face_current = np.maximum(face_net_current, 0.0)
    negative_face_current = np.maximum(-face_net_current, 0.0)
    positive_node_current = projector.node_charge_from_face_charge(positive_face_current)
    negative_node_current = projector.node_charge_from_face_charge(negative_face_current)
    physical_area = areas * 1e-12

    def manufactured_step(*, charge_node_c, duration_s, **_options):
        potential, _diagnostics = system.solve(charge_node_c)
        absolute_incident = float(
            (positive_face_current.sum() + negative_face_current.sum()) * duration_s)
        return DielectricChargingStep3DResult(
            charge_node_c=charge_node_c, charge_increment_node_c=np.zeros_like(charge_node_c),
            positive_face_current_density_a_m2=positive_face_current / physical_area,
            negative_face_current_density_a_m2=negative_face_current / physical_area,
            face_current_density_a_m2=face_net_current / physical_area,
            positive_current_node_a=positive_node_current,
            negative_current_node_a=negative_node_current,
            potential_before_v=potential,
            potential_after_v=potential,
            surface_transfer=SimpleNamespace(relative_charge_balance_error=0.0),
            transport=None, poisson_before=None, poisson_after=None,
            bidirectional_method_hint={}, bidirectional_sampling_provenance={},
            diagnostics=dict(
                incident_charge_c=0.0, positive_incident_charge_c=absolute_incident / 2.0,
                negative_incident_charge_c=-absolute_incident / 2.0,
                absolute_incident_charge_c=absolute_incident, deposited_charge_c=0.0,
                charge_conservation_residual_c=0.0,
                response_initial_bounce_budget=0, response_final_bounce_budget=0,
                response_emergency_bounce_limit=0,
                response_bounce_budget_extension_count=0,
                response_derived_bounce_budget=0,
                transport_lineage_replay_count=0,
                transport_lineage_replay_eligible_count=0,
                transport_lineage_replay_fraction=0.0,
                transport_edge_launch_inset_count=0,
                transport_trajectory_horizon_extension_count=0,
                transport_trajectory_initial_max_steps=1,
                transport_trajectory_final_max_steps=1,
                transport_trajectory_emergency_max_steps=1),
            known_limitations=())

    monkeypatch.setattr(
        charging_coevolution_3d, "advance_dielectric_charging_3d", manufactured_step)
    arguments = dict(
        poisson_system=system, initial_sigma_c_per_m2=np.zeros(len(faces)),
        boundary=PlasmaBoundaryState((_species("ion", 1, 1e15),), 2e-6),
        verts=vertices, faces=faces, areas=areas,
        face_centroids=vertices[faces].mean(axis=1),
        face_gas_normals=np.tile([0.0, 0.0, 1.0], (len(faces), 1)),
        face_material_id=np.ones(len(faces), dtype=int),
        source_bounds=(0.0, 3.0, 0.0, 3.0), source_z=1.0,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=1.0,
        patch_scales_m=(1e-6, 3e-6), potential_rate_tolerance_v_s=1e-6,
        timestep_s=1e-6, maximum_steps=1, mesh_length_unit_m=1e-6,
        trajectory_max_steps=1, transport_device="cpu")
    legacy = integrate_surface_charging_to_saturation_3d(**arguments)
    compatible = integrate_surface_charging_to_saturation_3d(
        **dict(arguments, compatible_q1_charge_state=True))

    assert projector.nullity == 4
    assert np.linalg.norm(legacy.face_charge_c) > 1e-21
    assert np.linalg.norm(compatible.face_charge_c) > 1e-23
    assert np.linalg.norm(compatible.face_charge_c) < 0.2 * np.linalg.norm(
        legacy.face_charge_c)
    assert np.allclose(
        compatible.charge_node_c, legacy.charge_node_c,
        rtol=0.0, atol=2e-36)
    assert np.allclose(
        compatible.potential_v, legacy.potential_v,
        rtol=2e-13, atol=2e-23)
    assert compatible.history[0]["unresolved_face_current_fraction"] > 0.9
    assert compatible.diagnostics["q1_face_coupling_nullity"] == 4
    assert compatible.history[0]["unresolved_face_current_projection_l1_c"] > 1e-21
    assert compatible.history[0]["face_to_node_update_relative_error"] < 5e-13
    assert max(compatible.history[0][
        "patch_q1_functional_null_sensitivity_max"]) > 0.1
    assert max(compatible.history[0][
        "patch_q1_unresolved_max_ion_normalized_imbalance"]) > 0.0
    assert len(compatible.history[0][
        "gate_patch_q1_resolved_max_ion_normalized_imbalance"]) == 2


def test_compatible_q1_pseudo_time_proposal_filters_null_current_and_caps_voltage():
    vertices = np.asarray([
        [float(i), float(j), 0.0]
        for i in range(4) for j in range(4)])

    def vertex(i, j):
        return 4 * i + j

    faces = []
    for i in range(3):
        for j in range(3):
            faces.extend((
                [vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1)],
                [vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1)]))
    faces = np.asarray(faces, dtype=int)
    areas = np.full(len(faces), 0.5)
    fixed = np.zeros((4, 4, 2), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(np.ones((3, 3, 1)), 1e-6, fixed)
    projector = CompatibleQ1SurfaceChargeProjector3D.from_triangles(
        system.shape, vertices, faces, coordinate_length_unit_m=1e-6)
    trial = np.random.default_rng(918).normal(size=len(faces))
    null_current = trial - projector.project_face_charge(trial)
    null_current *= 2e-15 / np.max(np.abs(null_current))
    resolved_current = projector.project_face_charge(
        np.random.default_rng(919).normal(size=len(faces)))
    resolved_current *= 2e-16 / np.max(np.abs(resolved_current))
    physical_area = areas * 1e-12
    current_density = (null_current + resolved_current) / physical_area
    proposal = propose_compatible_q1_pseudo_time_step_3d(
        system, vertices, faces, areas, np.zeros(len(faces)), current_density, 1e-6,
        mesh_length_unit_m=1e-6, maximum_potential_jump_v=1e6)

    expected_face_charge = projector.project_face_charge(
        (null_current + resolved_current) * 1e-6)
    assert proposal.diagnostics["current_q1_invisible_fraction"] > 0.9
    assert np.allclose(
        proposal.face_charge_c, expected_face_charge,
        rtol=2e-14, atol=2e-35)
    assert proposal.diagnostics["q1_node_update_relative_l1_error"] < 5e-13
    assert abs(proposal.diagnostics["global_charge_error_c"]) < 2e-33
    assert proposal.diagnostics["physical_time_advanced_s"] == 0.0
    assert proposal.diagnostics["exact_audit_required"]
    jump = proposal.diagnostics["maximum_potential_jump_v"]
    assert jump > 0.0
    with pytest.raises(ValueError, match="potential jump"):
        propose_compatible_q1_pseudo_time_step_3d(
            system, vertices, faces, areas, np.zeros(len(faces)), current_density, 1e-6,
            mesh_length_unit_m=1e-6, maximum_potential_jump_v=0.5 * jump)


def test_c3_records_inline_trajectory_horizon_contract_without_changing_the_step():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(
        trajectory_adaptive_horizon=True,
        trajectory_emergency_max_steps=4000)
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    item = result.history[0]
    assert item["transport_trajectory_horizon_extension_count"] == 0
    assert item["transport_trajectory_initial_max_steps"] == 2000
    assert item["transport_trajectory_final_max_steps"] == 2000
    assert item["transport_trajectory_emergency_max_steps"] == 4000
    assert result.diagnostics["trajectory_adaptive_horizon"]
    assert result.diagnostics["trajectory_emergency_max_steps"] == 4000
    assert result.diagnostics[
        "maximum_transport_trajectory_horizon_extension_count"] == 0

    with pytest.raises(ValueError, match="invalid C3"):
        integrate_surface_charging_to_saturation_3d(**dict(
            arguments, trajectory_emergency_max_steps=None))


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


def test_terminal_window_requires_elapsed_physical_time_and_gates_integrated_currents():
    flux = 3.0e15
    _system, arguments = _flat_problem((
        _species("ion", 1, flux), _species("electron", -1, flux)))
    arguments.update(
        maximum_steps=2, terminal_window_s=2e-6,
        stop_on_saturation=True)
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    assert result.converged
    assert result.accepted_steps == 2
    assert [item["terminal_window_ready"] for item in result.history] == [False, False, True]
    assert [item["saturation_gates_satisfied"] for item in result.history] == [
        False, False, True]
    assert result.diagnostics["gate_evaluation_mode"] == "terminal_window"
    assert result.diagnostics["terminal_window_s"] == 2e-6
    assert result.diagnostics["terminal_window_steps"] == 2
    assert result.diagnostics["final_potential_rate_max_v_s"] == 0.0
    assert result.diagnostics["final_instantaneous_potential_rate_max_v_s"] == 0.0
    assert np.array_equal(
        result.terminal_window_positive_face_current_density_a_m2,
        result.final_step.positive_face_current_density_a_m2)
    assert np.array_equal(
        result.terminal_window_negative_face_current_density_a_m2,
        result.final_step.negative_face_current_density_a_m2)
    assert all(item.b2_maximum_ion_normalized_imbalance == 0.0
               for item in result.patch_balance)


def test_terminal_window_does_not_average_ratios_or_hide_systematic_drift():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(
        maximum_steps=2, terminal_window_s=2e-6,
        stop_on_saturation=False)
    result = integrate_surface_charging_to_saturation_3d(**arguments)

    assert not result.converged
    assert result.history[-1]["terminal_window_ready"]
    assert result.diagnostics["final_potential_rate_max_v_s"] > 0.0
    assert result.diagnostics["final_maximum_patch_relative_imbalance"] == 1.0
    assert all(item.b2_maximum_ion_normalized_imbalance == 1.0
               for item in result.patch_balance)
    assert np.all(
        result.terminal_window_negative_face_current_density_a_m2 == 0.0)


def test_terminal_window_sparse_sliding_average_never_creates_negative_current():
    # This cancellation pattern made the old add/subtract accumulator return -1 on the first
    # face after a large sparse sample expired: (2**54 + 1) rounds to 2**54, then subtracting
    # 2**54 loses the positive unit current.  Direct reduction of the live interval is exact.
    samples = [
        {"positive": np.array([2.0 ** 54, 0.0])},
        {"positive": np.array([1.0, 3.0])},
        {"positive": np.array([0.0, 0.0])},
    ]

    mean = _terminal_window_current_mean(samples[1:], "positive")

    assert np.array_equal(mean, np.array([1.0, 3.0]))


def test_terminal_window_rejects_pseudo_time_and_nonintegral_fixed_windows():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))

    with pytest.raises(ValueError, match="fresh-scramble proposal or timestep"):
        integrate_surface_charging_to_saturation_3d(**dict(
            arguments, terminal_window_s=2e-6, timestep_policy="ser"))
    with pytest.raises(ValueError, match="integer multiple"):
        integrate_surface_charging_to_saturation_3d(**dict(
            arguments, terminal_window_s=1.5e-6))


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


def test_physical_poisson_arrivals_require_fresh_unselfcertified_time_evolution():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    with pytest.raises(ValueError, match="fresh-scramble proposal or timestep"):
        integrate_surface_charging_to_saturation_3d(
            **dict(arguments, physical_arrival_statistics="poisson"))
    result = integrate_surface_charging_to_saturation_3d(**dict(
        arguments, physical_arrival_statistics="poisson", scramble_mode="fresh",
        stop_on_saturation=False, maximum_steps=1))

    assert not result.converged
    assert result.diagnostics["physical_arrival_statistics"] == "poisson"
    assert all(item["physical_arrival_statistics"] == "poisson"
               for item in result.history)
    assert all(item["charge_conservation_relative_error"] < 3e-16
               for item in result.history)


def test_decreasing_gain_tail_is_conservative_replayable_and_never_self_certifies():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    schedule = dict(
        maximum_steps=3, stop_on_saturation=False, scramble_mode="fresh",
        sampling_seed_stride=101, timestep_policy="decreasing_gain",
        terminal_window_s=None, stochastic_gain_exponent=0.75,
        stochastic_gain_offset_steps=4)
    uninterrupted = integrate_surface_charging_to_saturation_3d(
        **dict(arguments, **schedule))

    gains = np.asarray([
        arguments["timestep_s"] * (4.0 / (4.0 + age)) ** 0.75
        for age in range(3)])
    expected_sigma = ECHARGE * flux * np.sum(gains)
    assert not uninterrupted.converged
    assert uninterrupted.physical_time_s == 0.0
    assert np.isclose(uninterrupted.pseudo_time_s, np.sum(gains))
    assert np.allclose(
        uninterrupted.sigma_c_per_m2, expected_sigma, rtol=3e-15)
    assert [item["stochastic_gain_age_steps"]
            for item in uninterrupted.history] == [0, 1, 2, 3]
    assert np.allclose(
        [item["timestep_s"] for item in uninterrupted.history],
        [gains[0], gains[0], gains[1], gains[2]], rtol=2e-15)
    assert uninterrupted.diagnostics["resume_stochastic_gain_age_steps"] == 3
    assert uninterrupted.diagnostics["stochastic_gain_exponent"] == 0.75

    first = integrate_surface_charging_to_saturation_3d(
        **{**arguments, **schedule, "maximum_steps": 1})
    resumed = integrate_surface_charging_to_saturation_3d(**{
        **arguments, **schedule, "maximum_steps": 2,
        "initial_sigma_c_per_m2": first.sigma_c_per_m2,
        "initial_sampling_epoch": first.diagnostics["resume_sampling_epoch"],
        "initial_stochastic_gain_age_steps": first.diagnostics[
            "resume_stochastic_gain_age_steps"]})
    assert np.array_equal(
        resumed.sigma_c_per_m2, uninterrupted.sigma_c_per_m2)
    assert np.array_equal(resumed.charge_node_c, uninterrupted.charge_node_c)
    assert np.array_equal(resumed.potential_v, uninterrupted.potential_v)
    assert resumed.diagnostics["resume_stochastic_gain_age_steps"] == 3


def test_decreasing_gain_refuses_frozen_samples_windows_and_self_certification():
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    base = dict(
        arguments, timestep_policy="decreasing_gain", terminal_window_s=None,
        stop_on_saturation=False)
    with pytest.raises(ValueError, match="fresh-scramble"):
        integrate_surface_charging_to_saturation_3d(**base)
    with pytest.raises(ValueError, match="fresh-scramble"):
        integrate_surface_charging_to_saturation_3d(**dict(
            base, scramble_mode="fresh", terminal_window_s=1e-6))
    with pytest.raises(ValueError, match="fresh-scramble"):
        integrate_surface_charging_to_saturation_3d(**dict(
            base, scramble_mode="fresh", stop_on_saturation=True))


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
    assert error.resume_stochastic_gain_age_steps == 0


def test_decreasing_gain_failure_checkpoint_preserves_exact_schedule_age(monkeypatch):
    flux = 2.0e15
    _system, arguments = _flat_problem((_species("ion", 1, flux),))
    arguments.update(
        maximum_steps=2, stop_on_saturation=False, scramble_mode="fresh",
        sampling_seed_stride=101, timestep_policy="decreasing_gain",
        terminal_window_s=None, stochastic_gain_exponent=0.75,
        stochastic_gain_offset_steps=4, initial_stochastic_gain_age_steps=7)
    original = charging_coevolution_3d.advance_dielectric_charging_3d
    calls = 0

    def fail_second_evaluation(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("manufactured decreasing-gain lookahead failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        charging_coevolution_3d, "advance_dielectric_charging_3d", fail_second_evaluation)
    with pytest.raises(
            SurfaceChargingSaturationError,
            match="manufactured decreasing-gain lookahead failure") as info:
        integrate_surface_charging_to_saturation_3d(**arguments)

    assert info.value.state_updates == 1
    assert info.value.resume_sampling_epoch == 1
    assert info.value.resume_stochastic_gain_age_steps == 8


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
    assert result.run_manifest["initial_geometry"]["phi"]["sha256"]
    assert result.run_manifest["species_role"]["Ar+"] == "energetic_bombardment"
    assert result.run_manifest["etchable_material_ids"] == [1]
    assert result.run_manifest["convergence_contract_revision"] == "CCA-2026-07-13-R2"
    assert result.run_manifest["surface_mechanism"]["type"]
    assert "machine_readable_provenance" in result.run_manifest["surface_mechanism"]
    assert result.run_manifest["charged_surface_response"]["parameters"][
        "grazing_reflection_probability"] == 0.95
    assert result.run_manifest["schema"] == "petch-charging-run-manifest-3d-v1"
    budget = result.run_manifest["recovery_and_error_budget"]
    assert budget["conservation"]["maximum_charge_remap_relative_balance_error"] == 0.0
    assert budget["inline_recovery"]["float64_lineage_replay_count"] >= 0
    json.dumps(dict(result.run_manifest))


def test_quasi_static_failure_preserves_the_exact_fresh_sampling_restart_epoch():
    geometry = _feature_plane()
    boundary = PlasmaBoundaryState((
        _species("Ar+", 1, 2.2e21, 100.0),
        _species("electron", -1, 1.1e21, 100.0),
    ), reference_plane_m=1.75e-6, provenance={"gate": "unbalanced restart"})
    options = _driver_options(maximum_steps=2)
    options.update(scramble_mode="fresh")

    with pytest.raises(
            SurfaceChargingSaturationError,
            match="failed signed B1/B2 saturation gates") as info:
        solve_charging_coevolution_3d(
            geometry, boundary,
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _physical_sputter_mechanism(), charging_system_builder=_plane_poisson_system,
            etchable_material_ids=(1,), duration_s=0.0, n_steps=1,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            charging_options=options, n_position=16, seed=47,
            trajectory_fixed_dt=0.005, trajectory_max_steps=1000,
            reinitialize=False, transport_device="cpu")

    error = info.value
    assert error.accepted_steps == 2
    assert error.state_updates == 2
    assert error.resume_sampling_epoch == 2
    assert error.history[-1]["sampling_epoch"] == 2


def test_public_charging_process_runs_the_same_unified_c3_engine(tmp_path):
    geometry = _feature_plane()
    process = PhysicalChargingProcess(
        geometry=geometry, boundary=_balanced_feature_boundary(),
        species_role={"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism=_physical_sputter_mechanism(),
        charging_system_builder=_plane_poisson_system,
        etchable_material_ids=(1,), duration_s=0.01, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=_driver_options(),
        solver_options=dict(
            n_position=32, seed=48, trajectory_fixed_dt=0.005,
            trajectory_max_steps=1000, reinitialize=False, transport_device="cpu"))

    result = process.run()

    assert process.engine == COMMON_CHARGING_ENGINE
    assert result.engine == COMMON_CHARGING_ENGINE
    assert result.solve.run_manifest["exact_operator"].startswith("hard visibility")
    assert result.surface_charge_c_per_m2.shape == result.solve.sigma_c_per_m2.shape
    assert result.steps[0].charging.converged

    continuation = process.continue_from(
        result, duration_s=0.01, n_steps=1, continuation_seed_stride=101)
    continued = continuation.run()
    assert continuation.geometry is result.geometry
    assert continuation.solver_options["seed"] == 149
    assert np.array_equal(
        continuation.solver_options["initial_sigma_c_per_m2"],
        result.surface_charge_c_per_m2)
    assert (continuation.solver_options["initial_surface_state"]
            is result.surface_state)
    assert continued.run_manifest["initial_surface_state_supplied"]
    assert continued.run_manifest["initial_surface_state_mesh_fingerprint"] == (
        result.solve.surface_state_mesh_fingerprint)

    checkpoint = PhysicalChargingCheckpoint3D.from_result(result.solve)
    path = tmp_path / "charged-step.npz"
    checkpoint.save(path)
    restored = PhysicalChargingCheckpoint3D.load(path)
    disk_continuation = process.continue_from_checkpoint(
        restored, duration_s=0.01, n_steps=1, continuation_seed_stride=101)
    assert disk_continuation.geometry.phi.flags.writeable is False
    assert disk_continuation.solver_options["restart_source_manifest_sha256"] == (
        checkpoint.source_manifest_sha256)
    assert np.array_equal(
        disk_continuation.solver_options["initial_sigma_c_per_m2"],
        result.surface_charge_c_per_m2)


def test_public_charging_process_refuses_required_field_override():
    geometry = _feature_plane()
    with pytest.raises(ValueError, match="cannot override required charged case fields"):
        PhysicalChargingProcess(
            geometry=geometry, boundary=_balanced_feature_boundary(),
            species_role={"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            mechanism=_physical_sputter_mechanism(),
            charging_system_builder=_plane_poisson_system,
            etchable_material_ids=(1,), duration_s=0.01, n_steps=1,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            charging_options=_driver_options(), solver_options={"charging_options": {}})


def test_profile_stationary_driver_uses_disjoint_ensembles_and_keeps_r2_diagnostics():
    geometry = _feature_plane()
    boundary = _balanced_feature_boundary()
    contract = ProfileChargingStationarityContract3D(
        potential_drift_tolerance_v=1e6,
        current_relative_l1_tolerance=10.0,
        transported_flux_relative_l1_tolerance=10.0,
        profile_velocity_relative_l1_tolerance=10.0,
        profile_increment_tolerance_m=1.0,
        minimum_independent_replicates=2,
        confidence_multiplier=0.0)
    result = solve_charging_coevolution_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        _physical_sputter_mechanism(), charging_system_builder=_plane_poisson_system,
        etchable_material_ids=(1,), duration_s=0.01, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=_driver_options(maximum_steps=1),
        charging_acceptance="profile_stationary", stationarity_contract=contract,
        stationarity_block_steps=1, stationarity_scoring_replicates=2,
        n_position=16, seed=49, trajectory_fixed_dt=0.005,
        trajectory_max_steps=1000, reinitialize=False, transport_device="cpu")

    step = result.steps[0]
    assert step.profile_stationarity.passed
    assert step.diagnostics["profile_stationarity_satisfied"]
    assert "retained_node_rms_relative_current_imbalance" in step.diagnostics
    assert "retained_node_max_relative_current_imbalance" in step.diagnostics
    assert len(set(
        step.profile_stationarity.diagnostics["first_scoring_sampling_epochs"]
        + step.profile_stationarity.diagnostics["second_scoring_sampling_epochs"])) == 4
    assert result.run_manifest["charging_acceptance"] == "profile_stationary"
    assert result.run_manifest["convergence_contract_revision"] == contract.revision
    assert not result.run_manifest["stationarity_contract"][
        "experimental_claim_authorized"]


def test_profile_stationary_draft_refuses_experimental_claims():
    geometry = _feature_plane()
    boundary = _balanced_feature_boundary()
    contract = ProfileChargingStationarityContract3D(
        1.0, 0.1, 0.1, 0.1, 1e-9)
    with pytest.raises(ValueError, match="cannot authorize experimental claims"):
        solve_charging_coevolution_3d(
            geometry, boundary,
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _physical_sputter_mechanism(), charging_system_builder=_plane_poisson_system,
            etchable_material_ids=(1,), duration_s=0.0, n_steps=1,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            charging_options=_driver_options(), charging_acceptance="profile_stationary",
            stationarity_contract=contract, stationarity_block_steps=1,
            experimental_claim=True, observable_tolerances=(
                ExperimentalObservableTolerance3D("notch", 0.1, 0.1),))


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


def test_physical_time_resolved_mode_coevolves_finite_arrivals_without_quasistatic_gate():
    geometry = _feature_plane()
    result = solve_charging_coevolution_3d(
        geometry, _balanced_feature_boundary(),
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        _physical_sputter_mechanism(), charging_system_builder=_plane_poisson_system,
        etchable_material_ids=(1,), duration_s=2e-9, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=dict(
            _driver_options(), physical_arrival_statistics="poisson"),
        bias_mode="physical_time_resolved", n_position=32, seed=59,
        trajectory_fixed_dt=0.005, trajectory_max_steps=1000,
        reinitialize=False, transport_device="cpu")

    assert result.run_manifest["mode"] == "physical_time_resolved"
    assert result.run_manifest["waveform"] is None
    assert all(step.charging.accepted_steps == 1 for step in result.steps)
    assert all(not step.diagnostics["saturation_required"] for step in result.steps)
    assert all(step.charging.diagnostics["physical_arrival_statistics"] == "poisson"
               for step in result.steps)
    assert "one realization" in " ".join(result.validity.known_limitations)


def test_public_finite_arrival_ensemble_runs_distinct_reproducible_seeds():
    geometry = _feature_plane()
    process = PhysicalChargingProcess(
        geometry=geometry, boundary=_balanced_feature_boundary(),
        species_role={"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism=_physical_sputter_mechanism(),
        charging_system_builder=_plane_poisson_system,
        etchable_material_ids=(1,), duration_s=1e-9, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=dict(
            _driver_options(), physical_arrival_statistics="poisson"),
        solver_options=dict(
            bias_mode="physical_time_resolved", n_position=16, seed=61,
            trajectory_fixed_dt=0.005, trajectory_max_steps=1000,
            reinitialize=False, transport_device="cpu"))
    ensemble = PhysicalChargingEnsembleProcess(
        process, realization_count=2, seed_stride=101).run()

    assert ensemble.engine == COMMON_CHARGING_ENSEMBLE_ENGINE
    assert ensemble.seeds == (61, 162)
    assert ensemble.realization_count == 2
    assert ensemble.mean_levelset.shape == geometry.phi.shape
    assert np.all(ensemble.standard_deviation_levelset >= 0.0)
    assert not ensemble.statistical_claim_ready

    refinement_contract = TwistEnsembleRefinementContract3D(
        minimum_realizations=4, mean_displacement_tolerance_m=1e-9,
        standard_deviation_tolerance_m=1e-9,
        onset_probability_tolerance=0.1)
    base_campaign = PhysicalChargingEnsembleProcess(
        process, realization_count=8, seed_stride=101)
    invalid_refined = PhysicalChargingEnsembleProcess(
        replace(process, solver_options=dict(process.solver_options, n_position=24)),
        realization_count=8, seed_stride=101)
    with pytest.raises(ValueError, match="only double n_position"):
        base_campaign.run_twist_refinement(
            invalid_refined, aspect_ratio=4.0, measurement_contract={},
            refinement_contract=refinement_contract)


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
