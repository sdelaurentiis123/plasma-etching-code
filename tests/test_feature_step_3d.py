import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.charging_poisson_3d import NodalPoissonSystem3D
from petch.feature_step_3d import (
    FeatureGeometry3D,
    advance_feature_step_3d,
    conservative_remap_surface_state,
    solve_feature_3d,
)
from petch.surface_kinetics import (
    EnergeticYield,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
)
from petch.threed import extract_mesh_3d


def _evidence():
    names = {
        "site_density_m2", "bulk_formula_density_m3", "polymer_monolayer_density_m2",
        "complex_formation_probability", "polymer_deposition_probability_on_substrate",
        "polymer_deposition_probability_on_polymer", "oxygen_polymer_etch_probability",
        "bare_sio2_yield", "complex_sio2_yield", "polymer_sputter_yield",
    }
    return {name: ParameterEvidence("manufactured moving-plane gate", "analytic") for name in names}


def _mechanism():
    yield_law = EnergeticYield(0.2, 20.0, 100.0)
    return ReducedSiO2FluorocarbonMechanism(ReducedSiO2FluorocarbonParameters(
        site_density_m2=5e18, bulk_formula_density_m3=2.2e28,
        polymer_monolayer_density_m2=4e18,
        complex_formation_probability={"CF2": 0.0},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={}, oxygen_species="O",
        oxygen_polymer_etch_probability=0.0,
        bare_sio2_yield=yield_law, complex_sio2_yield=yield_law,
        polymer_sputter_yield=yield_law, evidence=_evidence()))


def _plane_geometry():
    dx = 0.25; shape = (4, 4, 8); top = 0.95
    z = np.arange(shape[2]) * dx
    phi = np.broadcast_to(top - z, shape).copy()
    material = np.where(phi > 0.0, 1, 0)
    return FeatureGeometry3D(phi, material, dx, 1e-6), top


def _boundary():
    # Y=0.2 at 100 eV and Gamma=2.2e21 m^-2 s^-1 gives V=2e-8 m/s = 0.02 um/s.
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 0.0, [[0.0, 0.0, 1.0]], [1.0])
    return PlasmaBoundaryState((ion, neutral), reference_plane_m=1.75e-6)


def _charging_boundary():
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    electron = SpeciesBoundaryState(
        "electron", -1, 5.4858e-4, 2.2e22,
        [[0.0, 0.0, 1.0], [0.0, 0.0, np.sqrt(20.0)]], [0.9, 0.1])
    return PlasmaBoundaryState((ion, electron), reference_plane_m=1.75e-6)


def _plane_poisson_system(geometry):
    fixed = np.zeros(geometry.phi.shape, dtype=bool); fixed[:, :, -1] = True
    return NodalPoissonSystem3D(
        np.ones(tuple(np.asarray(geometry.phi.shape) - 1)),
        geometry.dx * geometry.mesh_length_unit_m, fixed)


def _area_weighted_height(phi, dx):
    _, _, centroids, areas = extract_mesh_3d(phi, dx)
    return float(np.dot(centroids[:, 2], areas) / areas.sum())


def test_one_physical_3d_step_moves_a_uniform_sio2_plane_by_flux_yield_over_density():
    geometry, initial_height = _plane_geometry()
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16384, seed=3, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.02, atol=0.002)
    assert np.isclose(result.diagnostics["max_velocity_m_s"], 2e-8, rtol=0.08)
    assert result.diagnostics["cfl_substeps"] == 1
    assert result.validity.within_declared_scope
    assert "conservative surface-state remap" in " ".join(result.validity.known_limitations)
    assert result.state_remap_diagnostics["old_topology"] == (1, 1)
    assert result.state_remap_diagnostics["new_topology"] == (1, 1)


def test_feature_step_refuses_surface_history_without_matching_mesh_fingerprint():
    geometry, _ = _plane_geometry()
    from petch.surface_kinetics import SiO2SurfaceState
    with pytest.raises(ValueError, match="fingerprint mismatch"):
        advance_feature_step_3d(
            geometry, _boundary(),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            surface_state=SiO2SurfaceState.bare((18,)), n_position=8,
            surface_state_mesh_fingerprint="not-the-current-mesh",
            reinitialize=False, transport_device="cpu")


def test_feature_step_accepts_surface_history_with_exact_current_mesh_fingerprint():
    geometry, _ = _plane_geometry()
    common = dict(
        geometry=geometry, boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism=_mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=8, reinitialize=False, transport_device="cpu")
    first = advance_feature_step_3d(**common)
    replay = advance_feature_step_3d(
        **common, surface_state=first.surface.state,
        surface_state_mesh_fingerprint=first.surface_state_mesh_fingerprint)
    assert replay.surface_state_mesh_fingerprint == first.surface_state_mesh_fingerprint
    assert np.array_equal(replay.surface.state.complex_fraction, first.surface.state.complex_fraction)


def test_surface_remap_preserves_material_integrals_and_coverage_bounds():
    from petch.surface_kinetics import SiO2SurfaceState
    old_centroid = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
    new_centroid = old_centroid + [0.0, 0.0, -0.1]
    old_area = np.array([1.0, 2.0, 1.5, 0.5])
    new_area = np.array([0.8, 2.2, 1.4, 0.6])
    material = np.array([1, 1, 2, 2])
    state = SiO2SurfaceState(
        [0.1, 0.8, 0.3, 0.6], [1e18, 2e18, 3e18, 4e18], [2e17, 4e17, 6e17, 8e17])
    remapped, diagnostics = conservative_remap_surface_state(
        state, old_centroid, old_area, material, new_centroid, new_area, material,
        dx=1.0, mesh_length_unit_m=1e-6)

    for material_id in (1, 2):
        old = material == material_id; new = material == material_id
        for before, after in (
                (state.complex_fraction, remapped.complex_fraction),
                (state.polymer_units_m2, remapped.polymer_units_m2),
                (state.removed_formula_units_m2, remapped.removed_formula_units_m2)):
            assert np.isclose(np.dot(before[old], old_area[old]),
                              np.dot(after[new], new_area[new]), rtol=2e-13)
    assert np.all((remapped.complex_fraction >= 0.0) & (remapped.complex_fraction <= 1.0))
    assert diagnostics["maximum_nearest_distance"] <= 0.1 + 1e-12


def test_multistep_solver_carries_remapped_state_and_matches_planar_total_motion():
    geometry, initial_height = _plane_geometry()
    result = solve_feature_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=2.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=16384, seed=13, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.04, atol=0.004)
    assert len(result.steps) == 2
    assert result.steps[0].next_surface_state_mesh_fingerprint == (
        result.steps[1].surface_state_mesh_fingerprint)
    assert result.surface_state_mesh_fingerprint == (
        result.steps[-1].next_surface_state_mesh_fingerprint)
    assert result.validity.within_declared_scope


def test_feature_step_uses_supplied_3d_potential_for_ion_energy_and_surface_velocity():
    geometry, initial_height = _plane_geometry()
    z = np.arange(geometry.phi.shape[2]) * geometry.dx
    potential = np.broadcast_to(10.0 * z / 1.75, geometry.phi.shape).copy()
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=geometry.dx, trajectory_fixed_dt=0.005,
        trajectory_max_steps=1000, n_position=16384, seed=37,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    expected_energy = 100.0 + 10.0 * (1.0 - initial_height / 1.75)
    expected_yield = 0.2 * (expected_energy - 20.0) / (100.0 - 20.0)
    expected_velocity = 2.2e21 * expected_yield / 2.2e28
    average_velocity = np.dot(
        result.surface.etch_velocity_m_s, result.active_face_area) / result.active_face_area.sum()
    assert result.transport.transport_model == "collisionless_fixed_step_nodal_field_3d"
    assert np.isclose(average_velocity, expected_velocity, rtol=3e-4)
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.02


def test_feature_step_solves_charge_reuses_ion_events_and_excludes_electron_from_chemistry():
    geometry, initial_height = _plane_geometry()
    mechanism = _mechanism()
    result = advance_feature_step_3d(
        geometry, _charging_boundary(),
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism, etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        trajectory_fixed_dt=0.005, trajectory_max_steps=2000,
        charging_poisson_system=_plane_poisson_system(geometry),
        charging_options=dict(
            max_iter=30, min_iter=2, current_balance_tol=1e-12,
            beta=0.5, response_energy_eV=4.0),
        n_position=64, seed=61, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    assert result.charging is not None and result.charging.converged
    assert result.transport is result.charging.transport
    support = (result.charging.positive_current_node_a
               + result.charging.negative_current_node_a) > 0.0
    assert np.allclose(
        result.charging.positive_current_node_a[support],
        result.charging.negative_current_node_a[support], rtol=1e-14)
    surface_voltage = result.charging.potential_v[:, :, 4]
    assert np.all((-20.0 < surface_voltage) & (surface_voltage < -1.0))

    populations = {item.name: item for item in result.transport.surface_fluxes.energetic_fluxes}
    assert set(populations) == {"Ar+", "electron"}
    ion_only_velocity = (
        populations["Ar+"].yield_rate_m2_s(mechanism.parameters.bare_sio2_yield)
        [result.active_face_index] / mechanism.parameters.bulk_formula_density_m3)
    assert np.allclose(result.surface.etch_velocity_m_s, ion_only_velocity, rtol=2e-13)
    assert result.diagnostics["self_consistent_charging"]
    assert result.diagnostics["charging_converged"]
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.02


def test_multistep_charged_profile_refuses_missing_charge_and_permittivity_remap():
    geometry, _ = _plane_geometry()
    with pytest.raises(ValueError, match="permittivity rebuild"):
        solve_feature_3d(
            geometry, _charging_boundary(),
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _mechanism(), etchable_material_ids=(1,), duration_s=2.0, n_steps=2,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            trajectory_fixed_dt=0.005,
            charging_poisson_system=_plane_poisson_system(geometry))
