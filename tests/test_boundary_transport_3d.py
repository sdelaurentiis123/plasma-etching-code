import numpy as np
import pytest
import warp as wp

from petch.boundary_state import (
    PlasmaBoundaryState, SpeciesBoundaryState, folded_normal_tangential_proposal,
    maxwellian_electron_boundary_state, qmc_boundary_proposal,
)
from petch.boundary_transport_3d import (
    estimate_diffuse_form_factors_3d,
    gather_boundary_state_ballistic_3d,
    gather_boundary_state_field_adjoint_3d,
    merge_boundary_transport_results_3d,
    trace_boundary_state_bidirectional_field_3d,
    trace_boundary_state_field_3d,
    trace_boundary_state_first_hit_3d,
)
from petch.surface_kinetics import (
    EnergeticYield,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
)


DEVICES = ["cpu"] + (["cuda:0"] if wp.is_cuda_available() else [])


def _flat_unit_plane():
    verts = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    areas = np.array([0.5, 0.5])
    return verts, faces, areas


def _boundary(position_m=None):
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2e19,
        velocity_sqrt_eV=[[0.0, 0.0, 10.0], [0.0, 0.0, np.sqrt(20.0)]],
        weight=[0.25, 0.75], position_m=position_m,
        provenance={"source": "manufactured"})
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 3e19, velocity_sqrt_eV=[[0.0, 0.0, 1.0]], weight=[1.0],
        provenance={"source": "manufactured"})
    return PlasmaBoundaryState(
        (ion, neutral), reference_plane_m=1e-6,
        provenance={"source": "manufactured flat-plane gate"})


@pytest.mark.parametrize("device", DEVICES)
def test_first_hit_3d_preserves_dimensional_species_flux_and_exact_energy_angle_events(device):
    verts, faces, areas = _flat_unit_plane()
    result = trace_boundary_state_first_hit_3d(
        _boundary(), {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts, faces, areas, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=256, seed=7, device=device)

    assert result.transport_model == "collisionless_absorbing_first_hit_3d"
    assert result.hit_probability == {"Ar+": 1.0, "CF2": 1.0}
    assert result.escape_probability == {"Ar+": 0.0, "CF2": 0.0}
    assert result.truncation_probability == {"Ar+": 0.0, "CF2": 0.0}
    neutral = result.surface_fluxes.neutral_flux_m2_s["CF2"]
    assert np.isclose(np.dot(neutral, areas), 3e19, rtol=1e-12)

    energetic = result.surface_fluxes.energetic_fluxes[0]
    assert energetic.name == "Ar+"
    assert np.isclose(np.dot(energetic.flux_m2_s, areas), 2e19, rtol=1e-12)
    assert set(np.round(energetic.event_energy_eV, 12)) == {20.0, 100.0}
    assert np.allclose(energetic.event_cosine_incidence, 1.0)
    law = EnergeticYield(0.2, 20.0, 100.0, energy_exponent=2.0)
    integrated_yield_rate = np.dot(energetic.yield_rate_m2_s(law), areas)
    assert np.isclose(integrated_yield_rate, 2e19 * 0.25 * 0.2, rtol=1e-12)


def test_deterministic_face_gather_reproduces_open_plane_flux_without_particle_tallies():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    result = gather_boundary_state_ballistic_3d(
        _boundary(), {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts, faces, areas, centroids, np.tile([0.0, 0.0, 1.0], (2, 1)),
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3, device="cpu")

    assert result.transport_model == "collisionless_deterministic_face_gather_3d"
    assert np.isclose(result.hit_probability["Ar+"], 1.0, rtol=1e-15)
    assert np.isclose(result.hit_probability["CF2"], 1.0, rtol=1e-15)
    assert np.allclose(result.surface_fluxes.neutral_flux_m2_s["CF2"], 3e19)
    energetic = result.surface_fluxes.energetic_fluxes[0]
    assert np.allclose(energetic.flux_m2_s, 2e19)
    assert set(np.round(energetic.event_energy_eV, 12)) == {20.0, 100.0}
    assert np.allclose(energetic.event_cosine_incidence, 1.0)


def test_deterministic_face_gather_uses_first_visible_surface_only():
    bottom, plane_faces, _ = _flat_unit_plane()
    top = bottom + [0.0, 0.0, 0.5]
    verts = np.vstack((bottom, top))
    faces = np.vstack((plane_faces, plane_faces + 4))
    areas = np.full(4, 0.5)
    centroids = verts[faces].mean(axis=1)
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 3e19, velocity_sqrt_eV=[[0.0, 0.0, 1.0]], weight=[1.0])
    boundary = PlasmaBoundaryState((neutral,), reference_plane_m=1e-6)
    result = gather_boundary_state_ballistic_3d(
        boundary, {"CF2": "neutral_reactant"}, verts, faces, areas, centroids,
        np.tile([0.0, 0.0, 1.0], (4, 1)),
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=1, device="cpu")

    flux = result.surface_fluxes.neutral_flux_m2_s["CF2"]
    assert np.array_equal(flux[:2], np.zeros(2))
    assert np.allclose(flux[2:], 3e19)
    assert np.isclose(result.hit_probability["CF2"], 1.0, rtol=1e-15)


def test_periodic_face_gather_projects_triangle_quadrature_to_exact_boundary_flux():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 3e19, velocity_sqrt_eV=[[0.2, 0.0, 1.0]], weight=[1.0])
    boundary = PlasmaBoundaryState((neutral,), reference_plane_m=1e-6)
    result = gather_boundary_state_ballistic_3d(
        boundary, {"CF2": "neutral_reactant"}, verts, faces, areas, centroids,
        np.tile([0.0, 0.0, 1.0], (2, 1)),
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3,
        periodic_lateral=True, domain_size=(1.0, 1.0, 1.0), device="cpu")

    flux = result.surface_fluxes.neutral_flux_m2_s["CF2"]
    assert result.hit_probability["CF2"] == 1.0
    assert np.isclose(np.dot(flux, areas), 3e19, rtol=2e-15)


def test_periodic_face_gather_conserves_each_energy_angle_atom_not_only_total_flux():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 2e19,
        velocity_sqrt_eV=[
            [0.1, 0.0, 10.0], [-0.2, 0.0, np.sqrt(20.0)], [0.0, 0.0, 2.0]],
        weight=[0.25, 0.75, 0.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1e-6)
    result = gather_boundary_state_ballistic_3d(
        boundary, {"Ar+": "energetic_bombardment"}, verts, faces, areas, centroids,
        np.tile([0.0, 0.0, 1.0], (2, 1)),
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3,
        periodic_lateral=True, domain_size=(1.0, 1.0, 1.0), device="cpu")

    energetic = result.surface_fluxes.energetic_fluxes[0]
    energies = np.unique(energetic.event_energy_eV)
    integrated = [
        np.dot(
            energetic.event_flux_m2_s[energetic.event_energy_eV == energy],
            areas[energetic.event_face[energetic.event_energy_eV == energy]])
        for energy in energies]
    assert np.allclose(energies, [20.04, 100.01])
    assert np.isclose(integrated[0], 0.75 * 2e19, rtol=2e-6)
    assert np.isclose(integrated[1], 0.25 * 2e19, rtol=2e-6)


def test_boundary_to_surface_chain_conserves_dimensional_formula_unit_removal():
    verts, faces, areas = _flat_unit_plane()
    transport = trace_boundary_state_first_hit_3d(
        _boundary(), {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts, faces, areas, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=256, seed=11, device="cpu")
    law = EnergeticYield(0.2, 20.0, 100.0, energy_exponent=2.0)
    evidence_names = {
        "site_density_m2", "bulk_formula_density_m3", "polymer_monolayer_density_m2",
        "complex_formation_probability", "polymer_deposition_probability_on_substrate",
        "polymer_deposition_probability_on_polymer", "oxygen_polymer_etch_probability",
        "bare_sio2_yield", "complex_sio2_yield", "polymer_sputter_yield",
    }
    evidence = {
        name: ParameterEvidence("manufactured flat-plane chain", "analytic")
        for name in evidence_names}
    mechanism = ReducedSiO2FluorocarbonMechanism(ReducedSiO2FluorocarbonParameters(
        site_density_m2=5e18, bulk_formula_density_m3=2.2e28,
        polymer_monolayer_density_m2=4e18,
        complex_formation_probability={"CF2": 0.0},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={}, oxygen_species="O",
        oxygen_polymer_etch_probability=0.0, bare_sio2_yield=law,
        complex_sio2_yield=law, polymer_sputter_yield=law, evidence=evidence))
    duration_s = 2.0
    surface = mechanism.advance(
        SiO2SurfaceState.bare((2,)), transport.surface_fluxes, duration_s)

    removed_per_source_area = np.dot(surface.state.removed_formula_units_m2, areas)
    expected = 2e19 * (0.25 * 0.2 + 0.75 * 0.0) * duration_s
    assert np.isclose(removed_per_source_area, expected, rtol=1e-12)


def test_disjoint_transport_merge_preserves_exact_event_objects_and_probabilities():
    verts, faces, areas = _flat_unit_plane(); boundary = _boundary()
    ion_boundary = PlasmaBoundaryState((boundary.get("Ar+"),), boundary.reference_plane_m)
    neutral_boundary = PlasmaBoundaryState((boundary.get("CF2"),), boundary.reference_plane_m)
    common = dict(
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=16, seed=19, device="cpu")
    ion = trace_boundary_state_first_hit_3d(
        ion_boundary, {"Ar+": "energetic_bombardment"}, **common)
    neutral = trace_boundary_state_first_hit_3d(
        neutral_boundary, {"CF2": "neutral_reactant"}, **common)
    merged = merge_boundary_transport_results_3d(ion, neutral)

    assert merged.surface_fluxes.energetic_fluxes[0] is ion.surface_fluxes.energetic_fluxes[0]
    assert np.array_equal(
        merged.surface_fluxes.neutral_flux_m2_s["CF2"],
        neutral.surface_fluxes.neutral_flux_m2_s["CF2"])
    assert merged.hit_probability == {"Ar+": 1.0, "CF2": 1.0}
    with pytest.raises(ValueError, match="disjoint"):
        merge_boundary_transport_results_3d(ion, ion)


def test_zero_nodal_field_reproduces_ballistic_3d_event_measure():
    verts, faces, areas = _flat_unit_plane()
    common = dict(
        boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=64, seed=29, device="cpu")
    ballistic = trace_boundary_state_first_hit_3d(**common)
    field = trace_boundary_state_field_3d(
        **common, nodal_potential_v=np.zeros((2, 2, 2)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, fixed_dt=0.01, max_steps=200)

    ballistic_events = ballistic.surface_fluxes.energetic_fluxes[0]
    field_events = field.surface_fluxes.energetic_fluxes[0]
    assert np.array_equal(field_events.event_face, ballistic_events.event_face)
    assert np.allclose(field_events.event_energy_eV, ballistic_events.event_energy_eV, atol=2e-5)
    assert np.allclose(field_events.event_cosine_incidence,
                       ballistic_events.event_cosine_incidence, atol=2e-7)
    assert np.allclose(field.surface_fluxes.neutral_flux_m2_s["CF2"],
                       ballistic.surface_fluxes.neutral_flux_m2_s["CF2"])


def test_periodic_field_transport_wraps_lateral_crossings_without_energy_change():
    verts, faces, areas = _flat_unit_plane()
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[2.0, 0.0, 1.0]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1e-6)
    arguments = dict(
        boundary=boundary, species_role={"Ar+": "energetic_bombardment"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=np.zeros((2, 2, 2)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        n_position=64, seed=31, fixed_dt=0.01, max_steps=300, device="cpu")

    open_cell = trace_boundary_state_field_3d(**arguments)
    periodic = trace_boundary_state_field_3d(**arguments, periodic_lateral=True)

    events = periodic.surface_fluxes.energetic_fluxes[0]
    assert open_cell.escape_probability["Ar+"] > 0.0
    assert periodic.hit_probability["Ar+"] == 1.0
    assert periodic.escape_probability["Ar+"] == 0.0
    assert np.allclose(events.event_energy_eV, 5.0, atol=2e-5)


def test_periodic_adjoint_field_gather_reproduces_flat_maxwellian_flux_and_energy():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)
    boundary = maxwellian_electron_boundary_state(
        4.0, 2e19, n_transverse=3, n_normal=4, reference_plane_m=1e-6)

    result = gather_boundary_state_field_adjoint_3d(
        boundary, {"electron": "charge_carrier"}, verts, faces, areas, centroids, normals,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=np.zeros((2, 2, 2)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        face_quadrature_points=3, fixed_dt=0.005, max_steps=1000,
        periodic_lateral=True, device="cpu")

    population = result.surface_fluxes.energetic_fluxes[0]
    integrated_flux = np.dot(population.flux_m2_s, areas)
    integrated_energy_rate = np.dot(
        np.bincount(population.event_face,
                    weights=population.event_flux_m2_s * population.event_energy_eV,
                    minlength=len(faces)), areas)
    assert np.isclose(result.hit_probability["electron"], 1.0, rtol=2e-6)
    assert np.isclose(integrated_flux, 2e19, rtol=2e-6)
    assert np.isclose(integrated_energy_rate / integrated_flux, 8.0, rtol=2e-6)


def test_adjoint_and_forward_field_transport_reproduce_maxwellian_barrier_tail():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)
    boundary = maxwellian_electron_boundary_state(
        4.0, 2e19, n_transverse=5, n_normal=8, reference_plane_m=1e-6)
    proposal = qmc_boundary_proposal(boundary.species[0], 12, seed=37)
    potential = np.zeros((2, 2, 2)); potential[:, :, 0] = -4.0
    common = dict(
        boundary=boundary, species_role={"electron": "charge_carrier"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        fixed_dt=0.0025, max_steps=4000, periodic_lateral=True, device="cpu")

    forward = trace_boundary_state_field_3d(
        **common, phase_space_log2_samples=14)
    adjoint = gather_boundary_state_field_adjoint_3d(
        **common, centroids=centroids, gas_normals=normals,
        face_quadrature_points=3, proposal_by_species={"electron": proposal})
    expected = np.exp(-1.0)

    assert np.isclose(forward.hit_probability["electron"], expected, atol=2.0 / 2 ** 14)
    assert np.isclose(adjoint.hit_probability["electron"], expected, rtol=0.015)
    assert np.isclose(
        adjoint.hit_probability["electron"], forward.hit_probability["electron"], rtol=0.015)


def test_bidirectional_field_transport_certifies_and_preserves_barrier_event_measure():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)
    boundary = maxwellian_electron_boundary_state(
        4.0, 2e19, n_transverse=3, n_normal=4, reference_plane_m=1e-6)
    potential = np.zeros((2, 2, 2)); potential[:, :, 0] = -4.0

    result = trace_boundary_state_bidirectional_field_3d(
        boundary, {"electron": "charge_carrier"}, verts, faces, areas, centroids, normals,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        forward_log2_samples=10, adjoint_log2_samples=8, n_replicates=4,
        element_absolute_tolerance=0.035, element_relative_tolerance=0.05,
        face_quadrature_points=3, fixed_dt=0.0025, max_steps=4000,
        periodic_lateral=True, seed=97, device="cpu")

    selection = result.selection_by_species["electron"]
    expected = np.exp(-1.0)
    assert selection.converged
    assert np.all(selection.estimator_consistent)
    assert np.isclose(result.transport.hit_probability["electron"], expected, atol=0.035)
    assert result.transport.surface_fluxes.energetic_fluxes[0].event_energy_eV.size > 0

    frozen = trace_boundary_state_bidirectional_field_3d(
        boundary, {"electron": "charge_carrier"}, verts, faces, areas, centroids, normals,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        forward_log2_samples=10, adjoint_log2_samples=8, n_replicates=4,
        element_absolute_tolerance=0.035, element_relative_tolerance=0.05,
        face_quadrature_points=3, fixed_dt=0.0025, max_steps=4000,
        periodic_lateral=True, seed=97, device="cpu",
        method_hint={"electron": selection.method}, require_certification=False)
    assert np.array_equal(frozen.selection_by_species["electron"].method, selection.method)
    assert np.allclose(
        frozen.transport.surface_fluxes.energetic_fluxes[0].flux_m2_s,
        result.transport.surface_fluxes.energetic_fluxes[0].flux_m2_s, rtol=0.0, atol=0.0)


def test_source_aligned_adjoint_rejects_surface_local_folded_proposal():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)
    boundary = maxwellian_electron_boundary_state(
        4.0, 2e19, n_transverse=3, n_normal=4, reference_plane_m=1e-6)
    folded = folded_normal_tangential_proposal(
        boundary.get("electron"), +1, name="electron")
    with pytest.raises(ValueError, match="surface-local coordinates"):
        gather_boundary_state_field_adjoint_3d(
            boundary, {"electron": "charge_carrier"}, verts, faces, areas,
            centroids, normals, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
            nodal_potential_v=np.zeros((2, 2, 2)), potential_origin=(0.0, 0.0, 0.0),
            potential_spacing=1.0, mesh_length_unit_m=1e-6,
            proposal_by_species={"electron": folded},
            proposal_frame_by_species={"electron": "source_aligned"},
            face_quadrature_points=1, fixed_dt=0.0025, max_steps=4000,
            periodic_lateral=True, device="cpu")


def test_subset_adjoint_gather_matches_full_collision_mesh_measure_on_selected_face():
    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.broadcast_to([0.0, 0.0, 1.0], centroids.shape)
    boundary = maxwellian_electron_boundary_state(
        4.0, 2e19, n_transverse=3, n_normal=4, reference_plane_m=1e-6)
    common = dict(
        boundary=boundary, species_role={"electron": "charge_carrier"},
        verts=verts, faces=faces, areas=areas, centroids=centroids, gas_normals=normals,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=np.zeros((2, 2, 2)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        face_quadrature_points=4, face_position_seed=113,
        fixed_dt=0.005, max_steps=1000, periodic_lateral=True, device="cpu")

    full = gather_boundary_state_field_adjoint_3d(**common)
    subset = gather_boundary_state_field_adjoint_3d(**common, gather_face_indices=[0])
    full_flux = full.surface_fluxes.energetic_fluxes[0].flux_m2_s
    subset_flux = subset.surface_fluxes.energetic_fluxes[0].flux_m2_s

    assert subset_flux[0] == full_flux[0]
    assert subset_flux[1] == 0.0


def test_linear_nodal_potential_gives_exact_electrostatic_energy_gain_under_refinement():
    verts, faces, areas = _flat_unit_plane()
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[0.0, 0.0, np.sqrt(20.0)]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1e-6)
    potential = np.zeros((2, 2, 2)); potential[:, :, 1] = 10.0
    arguments = dict(
        boundary=boundary, species_role={"Ar+": "energetic_bombardment"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        n_position=16, seed=31, max_steps=1000, device="cpu")
    coarse = trace_boundary_state_field_3d(**arguments, fixed_dt=0.02)
    fine = trace_boundary_state_field_3d(**arguments, fixed_dt=0.005)
    coarse_energy = coarse.surface_fluxes.energetic_fluxes[0].event_energy_eV
    fine_energy = fine.surface_fluxes.energetic_fluxes[0].event_energy_eV

    assert np.allclose(fine_energy, 30.0, atol=2e-3)
    assert np.max(np.abs(fine_energy - 30.0)) <= np.max(np.abs(coarse_energy - 30.0)) + 1e-6


def test_joint_phase_space_qmc_resolves_analytic_maxwellian_barrier_tail():
    # A wide collector removes lateral finite-target escape from this one-dimensional analytic gate.
    verts = np.array([
        [-100.0, -100.0, 0.0], [100.0, -100.0, 0.0],
        [100.0, 100.0, 0.0], [-100.0, 100.0, 0.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    areas = np.full(2, 20000.0)
    temperature = 4.0; barrier = temperature * np.log(10.0)
    boundary = maxwellian_electron_boundary_state(
        temperature, 1e19, n_transverse=3, n_normal=4,
        reference_plane_m=1e-6)
    potential = np.zeros((2, 2, 2)); potential[:, :, 0] = -barrier
    arguments = dict(
        boundary=boundary, species_role={"electron": "energetic_bombardment"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(-100.0, -100.0, 0.0),
        potential_spacing=(200.0, 200.0, 1.0), mesh_length_unit_m=1e-6,
        seed=47, fixed_dt=0.0025, max_steps=2000, device="cpu")
    coarse = trace_boundary_state_field_3d(
        **arguments, phase_space_log2_samples=8)
    fine = trace_boundary_state_field_3d(
        **arguments, phase_space_log2_samples=12)

    expected = np.exp(-barrier / temperature)
    assert coarse.transport_model.endswith("joint_qmc_3d")
    assert abs(fine.hit_probability["electron"] - expected) <= 1.0 / 2 ** 12
    assert (abs(fine.hit_probability["electron"] - expected)
            <= abs(coarse.hit_probability["electron"] - expected))
    landed_flux = np.dot(
        fine.surface_fluxes.energetic_fluxes[0].flux_m2_s, areas)
    assert np.isclose(
        landed_flux, 1e19 * fine.hit_probability["electron"], rtol=1e-14)


def test_field_3d_uses_each_physical_grid_spacing_in_electric_gradient():
    verts, faces, areas = _flat_unit_plane()
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[0.0, 0.0, np.sqrt(20.0)]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1e-6)
    # Ten volts along z over one mesh unit. Deliberately choose unequal x/y spacing so
    # treating the grid as isotropic would produce the wrong z extent or gradient.
    potential = np.broadcast_to(
        np.linspace(10.0, 0.0, 11), (3, 4, 11)).copy()
    result = trace_boundary_state_field_3d(
        boundary, {"Ar+": "energetic_bombardment"}, verts, faces, areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=(0.5, 1.0 / 3.0, 0.1), mesh_length_unit_m=1e-6,
        n_position=16, seed=41, fixed_dt=0.0025, max_steps=2000, device="cpu")

    impact_energy = result.surface_fluxes.energetic_fluxes[0].event_energy_eV
    assert np.allclose(impact_energy, 10.0, atol=3e-3)


def test_field_3d_separates_time_horizon_truncation_from_physical_escape():
    verts, faces, areas = _flat_unit_plane()
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[0.0, 0.0, np.sqrt(20.0)]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1e-6)
    arguments = dict(
        boundary=boundary, species_role={"Ar+": "energetic_bombardment"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        nodal_potential_v=np.zeros((2, 2, 2)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=1.0, mesh_length_unit_m=1e-6,
        n_position=16, seed=41, fixed_dt=0.005, max_steps=1, device="cpu")
    with pytest.raises(RuntimeError, match="exhausted max_steps"):
        trace_boundary_state_field_3d(**arguments)

    diagnostic = trace_boundary_state_field_3d(**arguments, allow_truncation=True)
    assert diagnostic.hit_probability["Ar+"] == 0.0
    assert diagnostic.escape_probability["Ar+"] == 0.0
    assert diagnostic.truncation_probability["Ar+"] == 1.0


@pytest.mark.parametrize("device", DEVICES)
def test_first_hit_3d_reports_geometric_oblique_incidence_without_angle_fit(device):
    verts = np.array([
        [-2.0, -2.0, 0.0], [2.0, -2.0, 0.0],
        [2.0, 2.0, 0.0], [-2.0, 2.0, 0.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    areas = np.array([8.0, 8.0])
    species = SpeciesBoundaryState(
        "ion", 1, 40.0, 1e19,
        velocity_sqrt_eV=[[0.5, 0.0, 1.0]], weight=[1.0])
    boundary = PlasmaBoundaryState((species,), reference_plane_m=1e-6)
    result = trace_boundary_state_first_hit_3d(
        boundary, {"ion": "energetic_bombardment"}, verts, faces, areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=16, seed=5, device=device)

    events = result.surface_fluxes.energetic_fluxes[0]
    assert result.hit_probability["ion"] == 1.0
    assert np.allclose(events.event_energy_eV, 1.25)
    assert np.allclose(events.event_cosine_incidence, 1.0 / np.sqrt(1.25), atol=2e-7)
    assert np.isclose(np.dot(events.flux_m2_s, areas), 1e19, rtol=1e-12)


def test_first_hit_3d_requires_complete_role_and_reference_plane_contracts():
    verts, faces, areas = _flat_unit_plane()
    common = dict(
        boundary=_boundary(), verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=8, device="cpu")
    with pytest.raises(ValueError, match="classify every"):
        trace_boundary_state_first_hit_3d(species_role={"Ar+": "energetic_bombardment"}, **common)
    with pytest.raises(ValueError, match="reference_plane"):
        trace_boundary_state_first_hit_3d(
            species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            mesh_origin_m=(0.0, 0.0, 1e-6), **common)
    with pytest.raises(ValueError, match="power of two"):
        trace_boundary_state_first_hit_3d(
            species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            **{**common, "n_position": 7})
    with pytest.raises(ValueError, match="areas must match"):
        trace_boundary_state_first_hit_3d(
            species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            areas=np.array([0.5, 0.4]),
            **{name: value for name, value in common.items() if name != "areas"})


def test_first_hit_3d_default_ray_distance_includes_source_plane_offset():
    verts, faces, areas = _flat_unit_plane()
    base = _boundary()
    far_boundary = PlasmaBoundaryState(
        base.species, reference_plane_m=100e-6, provenance=base.provenance)
    result = trace_boundary_state_first_hit_3d(
        far_boundary, {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts, faces, areas, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=100.0,
        mesh_length_unit_m=1e-6, n_position=8, seed=3, device="cpu")

    assert result.hit_probability == {"Ar+": 1.0, "CF2": 1.0}


def test_first_hit_3d_refuses_unimplemented_spatial_boundary_density():
    verts, faces, areas = _flat_unit_plane()
    positions = np.array([[0.25e-6, 0.25e-6], [0.75e-6, 0.75e-6]])
    with pytest.raises(ValueError, match="spatially uniform"):
        trace_boundary_state_first_hit_3d(
            _boundary(position_m=positions),
            {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
            verts, faces, areas, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
            mesh_length_unit_m=1e-6, n_position=8, device="cpu")


def test_diffuse_form_factor_estimator_classifies_open_plane_as_escape():
    verts, faces, _ = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    factors = estimate_diffuse_form_factors_3d(
        verts, faces, centroids, np.tile([0.0, 0.0, 1.0], (2, 1)),
        rays_per_face=16, seed=4, domain_size=(1.0, 1.0, 1.0), device="cpu")

    assert factors.source_face.size == factors.target_face.size == 0
    assert np.array_equal(factors.escape_fraction, [1.0, 1.0])


def test_diffuse_form_factor_estimator_replays_cavity_exchange_deterministically():
    bottom, faces, _ = _flat_unit_plane()
    top = bottom + [0.0, 0.0, 1.0]
    verts = np.vstack((bottom, top))
    faces = np.vstack((faces, faces + 4))
    centroids = verts[faces].mean(axis=1)
    normals = np.vstack((
        np.tile([0.0, 0.0, 1.0], (2, 1)),
        np.tile([0.0, 0.0, -1.0], (2, 1))))
    arguments = dict(
        verts=verts, faces=faces, centroids=centroids, gas_normals=normals,
        rays_per_face=64, seed=9, domain_size=(1.0, 1.0, 2.0), device="cpu")
    first = estimate_diffuse_form_factors_3d(**arguments)
    replay = estimate_diffuse_form_factors_3d(**arguments)

    assert first.transfer_fraction.size > 0
    assert np.array_equal(first.source_face, replay.source_face)
    assert np.array_equal(first.target_face, replay.target_face)
    assert np.array_equal(first.transfer_fraction, replay.transfer_fraction)
    assert np.array_equal(first.escape_fraction, replay.escape_fraction)
    outgoing = first.escape_fraction + np.bincount(
        first.source_face, weights=first.transfer_fraction, minlength=4)
    assert np.array_equal(outgoing, np.ones(4))


@pytest.mark.skipif(not wp.is_cuda_available(), reason="CUDA device unavailable")
def test_first_hit_3d_cpu_cuda_event_measure_parity():
    verts, faces, areas = _flat_unit_plane()
    arguments = dict(
        boundary=_boundary(),
        species_role={"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=256, seed=19)
    cpu = trace_boundary_state_first_hit_3d(**arguments, device="cpu")
    cuda = trace_boundary_state_first_hit_3d(**arguments, device="cuda:0")

    cpu_events = cpu.surface_fluxes.energetic_fluxes[0]
    cuda_events = cuda.surface_fluxes.energetic_fluxes[0]
    assert np.array_equal(cpu_events.event_face, cuda_events.event_face)
    assert np.array_equal(cpu_events.event_energy_eV, cuda_events.event_energy_eV)
    assert np.allclose(
        cpu_events.event_cosine_incidence, cuda_events.event_cosine_incidence,
        rtol=0.0, atol=2e-7)
    assert np.allclose(cpu_events.event_flux_m2_s, cuda_events.event_flux_m2_s)
    assert np.allclose(
        cpu.surface_fluxes.neutral_flux_m2_s["CF2"],
        cuda.surface_fluxes.neutral_flux_m2_s["CF2"])
