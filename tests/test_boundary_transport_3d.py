import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.boundary_transport_3d import trace_boundary_state_first_hit_3d
from petch.surface_kinetics import (
    EnergeticYield,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
)


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


def test_first_hit_3d_preserves_dimensional_species_flux_and_exact_energy_angle_events():
    verts, faces, areas = _flat_unit_plane()
    result = trace_boundary_state_first_hit_3d(
        _boundary(), {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts, faces, areas, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, n_position=256, seed=7, device="cpu")

    assert result.transport_model == "collisionless_absorbing_first_hit_3d"
    assert result.hit_probability == {"Ar+": 1.0, "CF2": 1.0}
    assert result.escape_probability == {"Ar+": 0.0, "CF2": 0.0}
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


def test_first_hit_3d_reports_geometric_oblique_incidence_without_angle_fit():
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
        mesh_length_unit_m=1e-6, n_position=16, seed=5, device="cpu")

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
