from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from petch.boundary_state import (
    PlasmaBoundaryState, SpeciesBoundaryState, maxwellian_electron_boundary_state,
)
from petch.charging_poisson_3d import NodalPoissonSystem3D
from petch.feature_step_3d import (
    FeatureGeometry3D,
    _physical_volume_topology_signature,
    _remove_unresolved_subcell_solid_components,
    advance_feature_step_3d,
    conservative_remap_surface_state,
    make_rectangular_trench_geometry_3d,
    solve_feature_3d,
)
from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_kinetics import (
    EnergeticYield,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
)
from petch.tabulated_chemistry import TabulatedSiClArMechanism, TabulatedSiSurfaceState
from petch.threed import extract_mesh_3d


INTERACTION_DATA = (
    Path(__file__).parents[1] / "data" / "surface_interactions" / "kounis_melas_2024")
SI_ATOM_DENSITY_M3 = 8.0 / (5.43e-10) ** 3


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


def _si_cl_ar_boundary():
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2e21, [[0.0, 0.0, 10.0]], [1.0])
    chlorine = SpeciesBoundaryState(
        "Cl2", 0, 70.906, 2e22, [[0.0, 0.0, 1.0]], [1.0])
    return PlasmaBoundaryState((ion, chlorine), reference_plane_m=1.75e-6)


def _si_cl_ar_mechanism():
    table = load_kounis_melas_2024_tables(INTERACTION_DATA).reactive_ion_etch
    return TabulatedSiClArMechanism(
        table, SI_ATOM_DENSITY_M3,
        ParameterEvidence(
            "Kounis-Melas OSTI 2589032 RIE in.lammps: diamond-Si lattice a=5.43 angstrom",
            "source_derived", supports_prediction_within_declared_domain=True))


def _plane_poisson_system(geometry):
    fixed = np.zeros(geometry.phi.shape, dtype=bool); fixed[:, :, -1] = True
    # Q1 value at the cell centre; this manufactured planar gate has one dielectric material.
    phi_center = sum(
        geometry.phi[i:i + geometry.phi.shape[0] - 1,
                     j:j + geometry.phi.shape[1] - 1,
                     k:k + geometry.phi.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0
    epsilon_r = np.where(phi_center > 0.0, 3.9, 1.0)
    return NodalPoissonSystem3D(
        epsilon_r,
        geometry.dx * geometry.mesh_length_unit_m, fixed)


def _area_weighted_height(phi, dx):
    _, _, centroids, areas = extract_mesh_3d(phi, dx)
    return float(np.dot(centroids[:, 2], areas) / areas.sum())


@pytest.mark.parametrize("dx", (0.02, 0.01, 0.005))
@pytest.mark.parametrize("opening_width", (0.06, 0.08, 0.10, 0.15, 0.18, 0.20))
def test_rectangular_trench_is_one_connected_substrate_at_jeon_widths(dx, opening_width):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=3.0 * dx, domain_height=2.35, dx=dx,
        opening_width=opening_width, mask_thickness=0.7,
        substrate_top=1.4, etched_depth=3.0 * dx)

    assert _physical_volume_topology_signature(geometry, (1,)) == (1, 1)


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
    assert not result.validity.parameter_evidence_supports_prediction
    assert "bare_sio2_yield" in result.validity.nonpredictive_parameters
    assert "conservative surface-state remap" in " ".join(result.validity.known_limitations)
    assert result.state_remap_diagnostics["old_topology"] == (1, 1)
    assert result.state_remap_diagnostics["new_topology"] == (1, 1)


def test_deterministic_face_gather_feature_step_is_independent_of_forward_particle_budget():
    geometry, _ = _plane_geometry()
    common = dict(
        geometry=geometry, boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism=_mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    one = advance_feature_step_3d(**common, n_position=1, seed=3)
    many = advance_feature_step_3d(**common, n_position=1024, seed=99)

    assert one.transport.transport_model == "collisionless_deterministic_face_gather_3d"
    assert np.array_equal(one.face_velocity_mesh_units_s, many.face_velocity_mesh_units_s)
    assert np.array_equal(one.geometry.phi, many.geometry.phi)


def test_deterministic_face_gather_moves_plane_by_same_dimensional_flux_law():
    geometry, initial_height = _plane_geometry()
    result = advance_feature_step_3d(
        geometry, _boundary(),
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
        seed=3, cfl_number=0.3, reinitialize=False, transport_device="cpu")

    final_height = _area_weighted_height(result.geometry.phi, geometry.dx)
    assert np.isclose(initial_height - final_height, 0.02, atol=0.002)
    assert np.isclose(result.diagnostics["max_velocity_m_s"], 2e-8, rtol=1e-12)
    assert result.transport.transport_model == "collisionless_deterministic_face_gather_3d"


def test_feature_step_diffusely_reemits_unreacted_neutrals_with_global_balance():
    geometry, _ = _plane_geometry()
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0])
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 3e20, [[0.0, 0.0, 1.0]], [1.0])
    boundary = PlasmaBoundaryState((ion, neutral), reference_plane_m=1.75e-6)
    result = advance_feature_step_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=64, seed=3, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={"rays_per_face": 16, "seed": 5})

    audit = result.diagnostics["neutral_radiosity"]["CF2"]
    assert audit["source_rate_s"] > 0.0
    assert audit["reacted_rate_s"] == 0.0
    assert np.isclose(audit["source_rate_s"], audit["escaped_rate_s"], rtol=1e-12)
    assert audit["relative_balance_error"] < 1e-12
    assert "flux_conservative_diffuse_radiosity" in result.transport.transport_model


def test_feature_step_radiosity_requires_explicit_probability_for_every_pinned_material():
    geometry, _ = _plane_geometry()
    material = np.array(geometry.material_id, copy=True)
    material[2:, :, :] = np.where(material[2:, :, :] > 0, 2, 0)
    mixed = FeatureGeometry3D(
        geometry.phi, material, geometry.dx, geometry.mesh_length_unit_m)
    with pytest.raises(ValueError, match="missing neutral reaction probability for material 2"):
        advance_feature_step_3d(
            mixed, _boundary(),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            n_position=16, seed=3, reinitialize=False, transport_device="cpu",
            neutral_radiosity_options={"rays_per_face": 8, "seed": 5})


def _static_trench_floor_neutral_flux(opening_width, *, rays_per_face=32):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.8, cell_length=0.15, domain_height=1.4, dx=0.05,
        opening_width=opening_width, mask_thickness=0.2,
        substrate_top=0.9, etched_depth=0.5)
    source_z = geometry.phi.shape[2] * geometry.dx - geometry.dx
    quadrature = maxwellian_electron_boundary_state(
        0.026, 3e20, n_transverse=3, n_normal=4,
        reference_plane_m=source_z * geometry.mesh_length_unit_m).species[0]
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, quadrature.flux_m2_s,
        quadrature.velocity_sqrt_eV, quadrature.weight)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[0.0, 0.0, 30.0]], [1.0])
    boundary = PlasmaBoundaryState(
        (ion, neutral), reference_plane_m=source_z * geometry.mesh_length_unit_m)
    parameters = replace(
        _mechanism().parameters, complex_formation_probability={"CF2": 0.2})
    mechanism = ReducedSiO2FluorocarbonMechanism(parameters)
    result = advance_feature_step_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.8, 0.0, 0.15), source_z=source_z,
        n_position=32, seed=7, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": rays_per_face, "seed": 11, "periodic_lateral": True,
            "domain_size": np.asarray(geometry.phi.shape) * geometry.dx,
            "nonetchable_reaction_probability_by_material": {2: {"CF2": 0.2}},
        })
    floor = result.active_face_centroid[:, 2] < 0.45
    flux = result.transport.surface_fluxes.neutral_flux_m2_s["CF2"][
        result.active_face_index[floor]]
    return float(np.average(flux, weights=result.active_face_area[floor]))


def test_diffuse_neutral_transport_widens_from_local_plane_to_trench_width_ring():
    narrow = _static_trench_floor_neutral_flux(0.2)
    wide = _static_trench_floor_neutral_flux(0.4)

    assert 0.0 < narrow < wide < 3e20


def test_unetched_rectangular_trench_mesh_drops_only_zero_measure_csg_faces():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=0.06, domain_height=2.35, dx=0.02,
        opening_width=0.2, mask_thickness=0.7, substrate_top=1.4, etched_depth=0.0)
    _, _, _, areas = extract_mesh_3d(geometry.phi, geometry.dx)

    assert areas.size > 0
    assert np.all(areas > 0.0)


def test_subcell_solid_component_is_removed_but_one_cell_support_is_preserved():
    phi = -np.ones((5, 5, 5))
    material = np.ones_like(phi, dtype=int)
    phi[0:2, 0:2, 0:2] = 1.0
    phi[4, 4, 4] = 0.1

    cleaned, removed = _remove_unresolved_subcell_solid_components(phi, material, (1,), 1.0)

    assert removed == 1
    assert np.all(cleaned[0:2, 0:2, 0:2] > 0.0)
    assert cleaned[4, 4, 4] < 0.0


def test_diffuse_neutral_trench_floor_flux_converges_with_form_factor_rule():
    coarse = _static_trench_floor_neutral_flux(0.2, rays_per_face=16)
    fine = _static_trench_floor_neutral_flux(0.2, rays_per_face=32)

    assert abs(coarse - fine) / fine < 0.15


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


def test_multistep_charged_profile_refuses_a_fixed_geometry_poisson_operator():
    geometry, _ = _plane_geometry()
    with pytest.raises(ValueError, match="geometry-dependent Poisson builder"):
        solve_feature_3d(
            geometry, _charging_boundary(),
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _mechanism(), etchable_material_ids=(1,), duration_s=2.0, n_steps=2,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            trajectory_fixed_dt=0.005,
            charging_poisson_system=_plane_poisson_system(geometry))


def test_multistep_quasistatic_charging_rebuilds_material_operator_and_reconverges():
    geometry, initial_height = _plane_geometry(); systems = []

    def build(current_geometry):
        system = _plane_poisson_system(current_geometry); systems.append(system)
        return system

    result = solve_feature_3d(
        geometry, _charging_boundary(),
        {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        _mechanism(), etchable_material_ids=(1,), duration_s=10.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        trajectory_fixed_dt=0.005, trajectory_max_steps=2000,
        charging_system_builder=build,
        charging_options=dict(
            max_iter=30, min_iter=2, current_balance_tol=1e-12,
            beta=0.5, response_energy_eV=4.0),
        n_position=64, seed=67, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    assert len(systems) == 2 and len(result.steps) == 2
    assert all(step.charging is not None and step.charging.converged for step in result.steps)
    assert np.count_nonzero(systems[1].epsilon_r == 3.9) < np.count_nonzero(
        systems[0].epsilon_r == 3.9)
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.15
    assert "quasi-static charging" in " ".join(result.validity.known_limitations)


def test_feature_step_refuses_nonconverged_charging_for_profile_motion():
    geometry, _ = _plane_geometry()
    with pytest.raises(ValueError, match="requires a converged"):
        advance_feature_step_3d(
            geometry, _charging_boundary(),
            {"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
            _mechanism(), etchable_material_ids=(1,), duration_s=1.0,
            source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            trajectory_fixed_dt=0.005,
            charging_poisson_system=_plane_poisson_system(geometry),
            charging_options={"require_converged": False})


def test_second_chemistry_runs_through_unchanged_transport_remap_and_interface_engine():
    geometry, initial_height = _plane_geometry(); mechanism = _si_cl_ar_mechanism()
    result = solve_feature_3d(
        geometry, _si_cl_ar_boundary(),
        {"Ar+": "energetic_bombardment", "Cl2": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=4.0, n_steps=2,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=4096, seed=73, cfl_number=0.3, reinitialize=False,
        transport_device="cpu")

    expected_velocity_m_s = 2e21 * 0.24182079610957588 / SI_ATOM_DENSITY_M3
    mean_velocity = np.mean(result.steps[0].surface.etch_velocity_m_s)
    assert isinstance(result.surface_state, TabulatedSiSurfaceState)
    assert np.isclose(mean_velocity, expected_velocity_m_s, rtol=0.01)
    assert all(step.surface.table_fingerprint == mechanism.table.fingerprint
               for step in result.steps)
    assert result.validity.parameter_evidence_supports_prediction
    assert result.validity.nonpredictive_parameters == ()
    assert result.surface_state.removed_atoms_m2.size == result.steps[-1].next_active_face_area.size
    assert _area_weighted_height(result.geometry.phi, geometry.dx) < initial_height - 0.03
