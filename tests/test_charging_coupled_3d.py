import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.charging_coupled_3d import (
    DielectricChargingConvergenceError,
    advance_dielectric_charging_3d,
    solve_dielectric_charging_steady_3d,
)
from petch.charging_poisson import EPS0
from petch.charging_poisson_3d import NodalPoissonSystem3D
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
    assert np.allclose(result.potential_before_v, 0.0, atol=1e-14)
    assert np.isclose(result.potential_after_v[:, :, 0].mean(), expected_voltage, rtol=2e-12)

    refined_arguments = dict(arguments); refined_arguments["n_position"] = 256
    refined = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=duration, **refined_arguments)
    assert (np.std(refined.potential_after_v[:, :, 0])
            < np.std(result.potential_after_v[:, :, 0]))


def test_equal_positive_and_negative_incident_currents_leave_dielectric_uncharged():
    flux = 3.0e15
    species = (_species("ion", 1, flux), _species("electron", -1, flux))
    system, arguments = _flat_dielectric_problem(species)
    result = advance_dielectric_charging_3d(
        charge_node_c=np.zeros(system.shape), duration_s=2e-3, **arguments)

    assert np.allclose(result.face_current_density_a_m2, 0.0, atol=1e-20)
    assert np.allclose(result.charge_increment_node_c, 0.0, atol=1e-32)
    assert np.allclose(result.potential_after_v, 0.0, atol=1e-14)


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
    system, arguments = _flat_dielectric_problem(_manufactured_floating_boundary())
    result = solve_dielectric_charging_steady_3d(
        initial_charge_node_c=np.zeros(system.shape), max_iter=10, min_iter=2,
        current_balance_tol=1e-12, beta=4.0, response_energy_eV=4.0,
        maximum_voltage_step=8.0, **arguments)

    assert result.converged
    assert result.rejected_steps > 0
    assert result.history[-1]["rms_relative_current_imbalance"] == 0.0
