import numpy as np
import pytest

from petch.boundary_state import (
    IonEnergyTransverseMaxwellianDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    maxwellian_electron_boundary_state,
)
from petch.boundary_transport_3d import BidirectionalSamplingProvenance3D
from petch.charging_coupled_3d import (
    CurrentBalanceMetrics3D,
    DielectricChargingConvergenceError,
    _freeze_certified_bidirectional_options,
    advance_dielectric_charging_3d,
    current_balance_metrics_3d,
    integrate_dielectric_charging_transient_3d,
    solve_dielectric_charging_steady_3d,
)
from petch.charging_poisson import EPS0
from petch.charging_poisson_3d import NodalPoissonSystem3D
from petch.charged_surface_cascade_3d import ChargedSurfaceCascade3DResult
from petch.charged_surface_response_3d import PerfectAbsorberChargedSurfaceResponse3D
from petch.sheath import ECHARGE


def _flat_dielectric_problem(species):
    cell_shape = (1, 1, 10); spacing_m = np.full(3, 0.1e-6)
    fixed = np.zeros((2, 2, 11), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(np.ones(cell_shape), spacing_m, fixed)
    vertices = np.array([
        [0.0, 0.0, 0.0], [0.1, 0.0, 0.0],
        [0.1, 0.1, 0.0], [0.0, 0.1, 0.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    areas = np.full(2, 0.005)
    boundary = PlasmaBoundaryState(tuple(species), reference_plane_m=1e-6)
    arguments = dict(
        poisson_system=system, boundary=boundary, verts=vertices, faces=faces, areas=areas,
        source_bounds=(0.0, 0.1, 0.0, 0.1), source_z=1.0,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=0.1,
        mesh_length_unit_m=1e-6, n_position=16, seed=43,
        trajectory_fixed_dt=0.0025, trajectory_max_steps=2000,
        transport_device="cpu")
    return system, arguments


def _species(name, charge_number, flux_m2_s, energy_eV=20.0):
    return SpeciesBoundaryState(
        name, charge_number, 40.0 if charge_number > 0 else 5.4858e-4,
        flux_m2_s, [[0.0, 0.0, np.sqrt(energy_eV)]], [1.0])


def test_current_balance_metrics_separate_local_low_throughput_and_patch_balance():
    positive = np.array([1.0, 9.0, 4.0, 4.0])
    negative = np.array([0.0, 9.0, 5.0, 3.0])
    raw = current_balance_metrics_3d(positive, negative)
    patches = current_balance_metrics_3d(
        positive, negative, group=np.array([0, 0, 1, 1]))

    assert isinstance(raw, CurrentBalanceMetrics3D)
    assert raw.maximum_relative_imbalance == 1.0
    assert raw.throughput_weighted_rms_relative_imbalance < raw.rms_relative_imbalance
    assert np.array_equal(patches.group, [0, 1])
    assert np.allclose(patches.positive_current_a, [10.0, 8.0])
    assert np.allclose(patches.negative_current_a, [9.0, 8.0])
    assert np.isclose(patches.maximum_relative_imbalance, 1.0 / 19.0)
    assert np.isclose(patches.global_relative_imbalance, 1.0 / 35.0)
    assert patches.positive_current_a.flags.writeable is False


def test_current_balance_metrics_validate_integrated_current_patch_inputs():
    with pytest.raises(ValueError, match="matching nonempty"):
        current_balance_metrics_3d([1.0], [-1.0])
    with pytest.raises(ValueError, match="integer array"):
        current_balance_metrics_3d([1.0], [1.0], group=[0.0])
    with pytest.raises(ValueError, match="retain"):
        current_balance_metrics_3d([1.0], [1.0], group=np.array([-1]))


def test_freezing_discovered_bidirectional_map_retains_certification_ceilings():
    method = {"ion": np.array(["forward", "adjoint"])}
    frozen = _freeze_certified_bidirectional_options(dict(
        forward_log2_samples=8, max_forward_log2_samples=12,
        adjoint_log2_samples=7, max_adjoint_log2_samples=11,
        face_quadrature_points=3, max_face_quadrature_points=15), method)

    assert frozen["forward_log2_samples"] == 12
    assert frozen["adjoint_log2_samples"] == 11
    assert frozen["face_quadrature_points"] == 15
    assert not frozen["require_certification"]
    assert np.array_equal(frozen["method_hint"]["ion"], method["ion"])


def test_externally_frozen_bidirectional_map_keeps_explicit_scoring_levels():
    method = {"ion": np.array(["forward", "adjoint"])}
    frozen = _freeze_certified_bidirectional_options(dict(
        method_hint=method, forward_log2_samples=9, max_forward_log2_samples=13,
        adjoint_log2_samples=9, max_adjoint_log2_samples=13), method)

    assert frozen["forward_log2_samples"] == 9
    assert frozen["adjoint_log2_samples"] == 9


def test_freezing_certified_map_uses_measured_sampling_levels_before_ceilings():
    method = {"ion": np.array(["forward", "adjoint", "adjoint"])}
    sampling = {"ion": BidirectionalSamplingProvenance3D(
        10, np.array([7, 8, 9]), np.array([3, 6, 12]), np.arange(4))}
    frozen = _freeze_certified_bidirectional_options(dict(
        forward_log2_samples=8, max_forward_log2_samples=14,
        adjoint_log2_samples=7, max_adjoint_log2_samples=13,
        face_quadrature_points=3, max_face_quadrature_points=24), method, sampling)

    assert frozen["forward_log2_samples"] == 10
    assert frozen["adjoint_log2_samples"] == 9
    assert frozen["face_quadrature_points"] == 12
    assert not frozen["require_certification"]


def test_physical_3d_charging_step_conserves_incident_charge_and_capacitance():
    flux = 2.0e15; duration = 1.0e-3
    system, arguments = _flat_dielectric_problem((_species("ion", 1, flux),))
    result = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration, **arguments)

    area_m2 = 0.01 * 1e-12
    expected_charge = ECHARGE * flux * area_m2 * duration
    expected_sigma = ECHARGE * flux * duration
    expected_voltage = expected_sigma * 1e-6 / EPS0
    assert np.isclose(result.charge_increment_node_c.sum(), expected_charge, rtol=1e-14)
    assert abs(result.diagnostics["charge_conservation_residual_c"]) < 1e-30
    assert result.surface_transfer.outgoing == ()
    assert result.diagnostics["surface_transfer_charge_balance_residual_c"] == 0.0
    assert result.diagnostics["surface_transfer_relative_charge_balance_error"] == 0.0
    assert np.allclose(result.potential_before_v, 0.0, atol=1e-14)
    assert np.isclose(result.potential_after_v[:, :, 0].mean(), expected_voltage, rtol=2e-12)

    refined_arguments = dict(arguments); refined_arguments["n_position"] = 256
    refined = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration, **refined_arguments)
    assert (np.std(refined.potential_after_v[:, :, 0])
            < np.std(result.potential_after_v[:, :, 0]))


def test_explicit_perfect_absorber_uses_unified_response_path_without_changing_engine_result():
    flux = 2.0e15; duration = 1.0e-3
    system, arguments = _flat_dielectric_problem((_species("ion", 1, flux),))
    baseline = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration, **arguments)
    explicit = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration,
        charged_surface_response=PerfectAbsorberChargedSurfaceResponse3D(),
        face_material_id=np.array(["SiO2", "SiO2"]),
        face_gas_normals=np.tile([0.0, 0.0, 1.0], (2, 1)), **arguments)

    assert np.array_equal(explicit.face_current_density_a_m2,
                          baseline.face_current_density_a_m2)
    assert np.array_equal(explicit.charge_increment_node_c,
                          baseline.charge_increment_node_c)
    assert np.array_equal(explicit.potential_after_v, baseline.potential_after_v)
    assert explicit.surface_transfer.completed
    assert explicit.surface_transfer.relative_charge_balance_error < 3e-16


def test_physical_time_final_audit_uses_the_same_charged_surface_response_operator():
    flux = 2.0e15
    system, arguments = _flat_dielectric_problem((_species("ion", 1, flux),))
    result = integrate_dielectric_charging_transient_3d(
        initial_charge_node_c=np.zeros(system.shape), timestep_s=1e-6, n_steps=0,
        charged_surface_response=PerfectAbsorberChargedSurfaceResponse3D(),
        face_material_id=np.array(["SiO2", "SiO2"]),
        face_gas_normals=np.tile([0.0, 0.0, 1.0], (2, 1)), **arguments)

    assert isinstance(result.surface_transfer, ChargedSurfaceCascade3DResult)
    assert result.surface_transfer.completed
    assert result.surface_transfer.relative_charge_balance_error < 3e-16


def test_equal_positive_and_negative_incident_currents_leave_dielectric_uncharged():
    flux = 3.0e15
    species = (_species("ion", 1, flux), _species("electron", -1, flux))
    system, arguments = _flat_dielectric_problem(species)
    result = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=2e-3, **arguments)

    assert np.allclose(result.face_current_density_a_m2, 0.0, atol=1e-20)
    assert np.allclose(result.charge_increment_node_c, 0.0, atol=1e-32)
    assert np.allclose(result.potential_after_v, 0.0, atol=1e-14)


def test_charging_step_routes_directional_ions_forward_and_maxwellian_electrons_adjoint():
    flux = 3.0e15
    ion = _species("ion", 1, flux)
    electron = maxwellian_electron_boundary_state(
        4.0, flux, n_transverse=3, n_normal=4,
        reference_plane_m=1e-6).species[0]
    system, arguments = _flat_dielectric_problem((ion, electron))
    faces = arguments["faces"]
    centroids = arguments["verts"][faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)

    result = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=2e-3,
        transport_estimator={"ion": "forward", "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        periodic_lateral=True, **arguments)

    assert set(result.transport.hit_probability) == {"ion", "electron"}
    assert "field_adjoint_gather_3d" in result.transport.transport_model
    assert "fixed_step_nodal_field_3d" in result.transport.transport_model
    one_species_charge = ECHARGE * flux * 0.01e-12 * 2e-3
    assert abs(result.charge_increment_node_c.sum()) < 1e-6 * one_species_charge


def test_charging_step_consumes_a_certified_bidirectional_face_event_measure():
    flux = 3.0e15
    ion = SpeciesBoundaryState(
        "ion", 1, 40.0, flux, [[0.0, 0.0, np.sqrt(20.0)]], [1.0],
        density_model=IonEnergyTransverseMaxwellianDensity(
            np.array([19.0, 21.0]), np.array([1.0]), 0.1))
    electron = maxwellian_electron_boundary_state(
        4.0, flux, n_transverse=3, n_normal=4,
        reference_plane_m=1e-6).species[0]
    system, arguments = _flat_dielectric_problem((ion, electron))
    arguments = dict(arguments); arguments["trajectory_max_steps"] = 100000
    faces = arguments["faces"]
    centroids = arguments["verts"][faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)

    result = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=1e-6,
        transport_estimator="bidirectional", face_centroids=centroids,
        face_gas_normals=normals, periodic_lateral=True,
        bidirectional_options=dict(
            forward_log2_samples=8, adjoint_log2_samples=6, n_replicates=4,
            element_absolute_tolerance=0.08, element_relative_tolerance=0.1,
            face_quadrature_points=4), **arguments)

    assert result.transport.transport_model.endswith("bidirectional_3d_periodic_cell")
    assert set(result.transport.hit_probability) == {"ion", "electron"}
    assert set(result.bidirectional_sampling_provenance) == {"ion", "electron"}
    assert all(population.event_energy_eV.size > 0
               for population in result.transport.surface_fluxes.energetic_fluxes)


def test_second_charging_step_uses_first_steps_self_consistent_field():
    flux = 1.0e15
    system, arguments = _flat_dielectric_problem((_species("ion", 1, flux),))
    # A uniform positive sheet reaches +10 V after this physical interval.
    duration = 10.0 * EPS0 / (ECHARGE * flux * 1e-6)
    first = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration, **arguments)
    second = advance_dielectric_charging_3d(
        charge_node_c=first.charge_node_c, duration_s=duration, **arguments)

    impact_energy = second.transport.surface_fluxes.energetic_fluxes[0].event_energy_eV
    assert np.isclose(first.potential_after_v[:, :, 0].mean(), 10.0, rtol=2e-12)
    assert np.allclose(second.potential_before_v, first.potential_after_v, rtol=1e-13)
    assert np.isclose(impact_energy.mean(), 10.0, atol=3e-4)


def test_physical_time_driver_matches_manual_steps_and_audits_final_current_state():
    flux = 1.0e15
    system, arguments = _flat_dielectric_problem((_species("ion", 1, flux),))
    duration = 2.0 * EPS0 / (ECHARGE * flux * 1e-6)
    first = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration, **arguments)
    second = advance_dielectric_charging_3d(
        charge_node_c=first.charge_node_c, duration_s=duration, **arguments)

    transient = integrate_dielectric_charging_transient_3d(
        initial_charge_node_c=np.zeros(system.shape), timestep_s=duration, n_steps=2,
        **arguments)

    assert np.allclose(transient.charge_node_c, second.charge_node_c, rtol=1e-13)
    assert np.allclose(transient.potential_v, second.potential_after_v, rtol=1e-13)
    assert len(transient.history) == 3
    assert transient.charge_history_node_c.shape == (3, *system.shape)
    assert np.allclose(transient.charge_history_node_c[-1], transient.charge_node_c)
    assert transient.diagnostics["updates_completed"] == 2
    assert transient.history[-1]["physical_time_s"] == 2.0 * duration
    assert transient.history[-1]["max_relative_current_imbalance_node"] == 1.0
    assert abs(transient.diagnostics["cumulative_charge_conservation_residual_c"]) < 1e-30
    assert transient.surface_transfer.outgoing == ()
    assert transient.diagnostics["final_surface_transfer_charge_balance_residual_c_s"] == 0.0
    assert transient.positive_current_node_a.flags.writeable is False


def _manufactured_floating_boundary():
    ion = _species("ion", 1, 1e15, energy_eV=100.0)
    # Ten times the electron flux. Below -1 V the 90% one-eV population reflects, leaving the
    # 10% twenty-eV tail: its landing current then exactly equals the ion current.
    electron = SpeciesBoundaryState(
        "electron", -1, 5.4858e-4, 1e16,
        [[0.0, 0.0, 1.0], [0.0, 0.0, np.sqrt(20.0)]], [0.9, 0.1])
    return ion, electron


def test_steady_3d_solver_converges_the_physical_local_current_equation():
    system, arguments = _flat_dielectric_problem(_manufactured_floating_boundary())
    result = solve_dielectric_charging_steady_3d(
        initial_charge_node_c=np.zeros(system.shape), max_iter=10, min_iter=2,
        current_balance_tol=1e-12, beta=0.5, response_energy_eV=4.0,
        **arguments)

    support = (result.positive_current_node_a + result.negative_current_node_a) > 0.0
    surface_voltage = result.potential_v[:, :, 0]
    assert result.converged
    assert result.history[0]["max_relative_current_imbalance"] > 0.8
    assert result.history[-1]["max_relative_current_imbalance"] == 0.0
    assert np.allclose(
        result.positive_current_node_a[support], result.negative_current_node_a[support],
        rtol=1e-14)
    assert np.all((-20.0 < surface_voltage) & (surface_voltage < -1.0))


def test_steady_3d_solver_refuses_to_label_an_unevaluated_proposal_converged():
    system, arguments = _flat_dielectric_problem(_manufactured_floating_boundary())
    with pytest.raises(DielectricChargingConvergenceError) as caught:
        solve_dielectric_charging_steady_3d(
            initial_charge_node_c=np.zeros(system.shape), max_iter=1, min_iter=1,
            current_balance_tol=1e-12, **arguments)

    result = caught.value.result
    assert not result.converged
    assert len(result.history) == 1
    assert np.allclose(result.charge_node_c, 0.0)
    assert np.allclose(result.potential_v, 0.0)


def test_steady_3d_solver_rejects_a_current_balance_worsening_trial():
    ion = _species("ion", 1, 1e15, energy_eV=100.0)
    electron = SpeciesBoundaryState(
        "electron", -1, 5.4858e-4, 1e16,
        [[0.0, 0.0, 1.0], [0.0, 0.0, np.sqrt(5.0)],
         [0.0, 0.0, np.sqrt(20.0)]],
        [0.8, 0.1, 0.1])
    system, arguments = _flat_dielectric_problem((ion, electron))
    result = solve_dielectric_charging_steady_3d(
        initial_charge_node_c=np.zeros(system.shape), max_iter=10, min_iter=2,
        current_balance_tol=1e-12, beta=4.0, response_energy_eV=4.0,
        maximum_voltage_step=30.0, **arguments)

    assert result.converged
    assert result.rejected_steps > 0
    assert result.history[-1]["rms_relative_current_imbalance"] == 0.0


def test_anderson_update_converges_same_physical_floating_root():
    system, arguments = _flat_dielectric_problem(_manufactured_floating_boundary())
    result = solve_dielectric_charging_steady_3d(
        initial_charge_node_c=np.zeros(system.shape), max_iter=8, min_iter=2,
        current_balance_tol=1e-12, beta=1.0, nonlinear_update="anderson",
        anderson_depth=3, **arguments)

    assert result.converged
    assert result.history[-1]["max_relative_current_imbalance"] == 0.0


def _continuous_maxwellian_floating_problem():
    temperature = 4.0; ion_flux = 1e15
    ion = SpeciesBoundaryState(
        "ion", 1, 40.0, ion_flux, [[0.0, 0.0, 10.0]], [1.0],
        density_model=IonEnergyTransverseMaxwellianDensity(
            np.array([99.0, 101.0]), np.array([1.0]), 0.01))
    electron = maxwellian_electron_boundary_state(
        temperature, 10.0 * ion_flux, n_transverse=3, n_normal=4,
        reference_plane_m=1e-6).species[0]
    boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=1e-6)
    spacing_m = np.array([200e-6, 200e-6, 0.1e-6])
    fixed = np.zeros((2, 2, 11), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(np.ones((1, 1, 10)), spacing_m, fixed)
    vertices = np.array([
        [-100.0, -100.0, 0.0], [100.0, -100.0, 0.0],
        [100.0, 100.0, 0.0], [-100.0, 100.0, 0.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3]]); areas = np.full(2, 20000.0)
    barrier = temperature * np.log(10.0)
    sigma = -EPS0 * barrier / 1e-6
    initial_charge = np.zeros(system.shape)
    initial_charge[:, :, 0] = sigma * spacing_m[0] * spacing_m[1] / 4.0
    arguments = dict(
        poisson_system=system, initial_charge_node_c=initial_charge, boundary=boundary,
        verts=vertices, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        potential_origin=(-100.0, -100.0, 0.0), potential_spacing=(200.0, 200.0, 0.1),
        mesh_length_unit_m=1e-6, n_position=16, seed=59,
        trajectory_fixed_dt=0.01, trajectory_max_steps=20000,
        transport_device="cpu", max_iter=1, min_iter=1,
        current_balance_tol=1e-8, phase_space_replicates=8,
        current_confidence_sigma=2.0, require_converged=False)
    return arguments


def test_replicated_joint_phase_space_narrows_current_confidence_envelope():
    arguments = _continuous_maxwellian_floating_problem()
    coarse = solve_dielectric_charging_steady_3d(
        **arguments, phase_space_log2_samples=6)
    fine = solve_dielectric_charging_steady_3d(
        **arguments, phase_space_log2_samples=10)

    coarse_envelope = coarse.history[0]["confidence_envelope_max_relative_current_imbalance"]
    fine_envelope = fine.history[0]["confidence_envelope_max_relative_current_imbalance"]
    assert not coarse.converged and not fine.converged
    assert np.max(fine.net_current_stderr_node_a) > 0.0
    assert fine_envelope < coarse_envelope


def test_current_estimator_raises_nested_sobol_level_until_uncertainty_is_resolved():
    arguments = _continuous_maxwellian_floating_problem()
    result = solve_dielectric_charging_steady_3d(
        **arguments, phase_space_log2_samples=6,
        phase_space_max_log2_samples=10, current_estimator_relative_tol=0.03)

    state = result.history[0]
    assert state["current_estimator_converged"]
    assert 6 < state["phase_space_log2_samples"] <= 10
    assert state["current_estimator_max_relative_uncertainty"] <= 0.03


def test_current_replicates_require_full_continuous_phase_space_sampling():
    system, arguments = _flat_dielectric_problem(_manufactured_floating_boundary())
    with pytest.raises(ValueError, match="joint continuous-density"):
        solve_dielectric_charging_steady_3d(
            initial_charge_node_c=np.zeros(system.shape), max_iter=2, min_iter=1,
            phase_space_replicates=2, **arguments)
