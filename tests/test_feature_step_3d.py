import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import FeatureGeometry3D, advance_feature_step_3d
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
    assert "no conservative remap yet" in " ".join(result.validity.known_limitations)


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
