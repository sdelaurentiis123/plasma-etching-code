import numpy as np
import pytest
import petch.feature_step_3d as feature_step_module

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import (
    SurfaceTopologyChangeError, advance_feature_step_3d,
    _retire_strictly_extinguished_material_levelsets,
    make_rectangular_trench_geometry_3d,
)
from petch.material_mechanism_3d import MaterialMechanismRouter3D, MaterialSurfaceState3D
from petch.physical_sputtering import PhysicalSputterMechanism, PhysicalSputterParameters
from petch.surface_kinetics import EnergeticYield, FaceResolvedEnergeticFlux, ParameterEvidence, SurfaceFluxes
from petch.surface_product_redeposition_3d import (
    SurfaceProductRedepositionContract3D, SurfaceProductRedepositionLaw3D,
)


def _sputter(material, inventory, product, reference_yield):
    evidence = {
        name: ParameterEvidence(
            "manufactured material-router gate", "analytic",
            supports_prediction_within_declared_domain=True)
        for name in (
            "bulk_material_unit_density_m3", "sputter_yield",
            "emitted_product_mass_amu", "emission_angular_model",
            "emission_energy_model")}
    return PhysicalSputterMechanism(PhysicalSputterParameters(
        material_name=material, material_inventory_name=inventory,
        projectile_species=("Ar+",), bulk_material_unit_density_m3=1e28,
        sputter_yield=EnergeticYield(reference_yield, 20.0, 100.0),
        emitted_product_name=product, emitted_product_mass_amu=28.0,
        emitted_material_units_per_particle=1.0,
        emission_angular_model="diffuse_cosine", emission_energy_model="thompson",
        emission_energy_parameters={
            "surface_binding_energy_eV": 4.0, "maximum_energy_eV": 100.0},
        evidence=evidence))


def _router():
    return MaterialMechanismRouter3D(
        {1: _sputter("substrate", "substrate_units", "substrate_product", 0.2),
         2: _sputter("mask", "mask_units", "mask_product", 0.05)},
        provenance={1: "manufactured substrate law", 2: "manufactured mask law"})


def test_material_router_keeps_state_velocity_and_product_ledgers_separate():
    material = np.array([1, 1, 2, 2])
    fluxes = SurfaceFluxes({}, (FaceResolvedEnergeticFlux(
        "Ar+", 4, np.arange(4), np.full(4, 1e20),
        np.full(4, 100.0), np.ones(4)),))
    router = _router()
    state = router.initial_state_by_material(material)

    result = router.advance_by_material(state, fluxes, 1.0, material)

    assert isinstance(result.state, MaterialSurfaceState3D)
    assert set(result.material_results) == {1, 2}
    assert np.allclose(result.etch_velocity_m_s[:2], 2e-9)
    assert np.allclose(result.etch_velocity_m_s[2:], 0.5e-9)
    assert result.material_exchange.product_routing_complete
    assert len(result.product_populations) == 2
    assert set(result.material_exchange.removed_units_m2) == {
        "substrate_units", "mask_units"}
    assert result.validity.parameter_evidence_supports_prediction
    assert router.provenance["materials"]["1"]["evidence"] == (
        "manufactured substrate law")

    with pytest.raises(ValueError, match="machine-readable"):
        MaterialMechanismRouter3D(
            {1: _sputter("substrate", "units", "product", 0.2)},
            provenance={1: object()})


def test_common_feature_engine_moves_mask_and_substrate_with_their_own_laws():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.2, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.3,
        substrate_top=1.0, etched_depth=0.2)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 39.948, 1e21, [[0.0, 0.0, 10.0]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1.8e-6)
    before = {
        material_id: np.asarray(field).copy()
        for material_id, field in geometry.material_levelsets.items()}

    result = advance_feature_step_3d(
        geometry, boundary, {"Ar+": "energetic_bombardment"}, _router(),
        etchable_material_ids=(1, 2), duration_s=1.0,
        source_bounds=(-0.1, 1.1, -0.1, 0.3), source_z=1.8,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    assert isinstance(result.next_surface_state, MaterialSurfaceState3D)
    assert set(np.unique(result.face_material_id[result.active_face_index])) == {1, 2}
    assert result.diagnostics["product_population_count"] == 2
    assert result.diagnostics["product_routing_complete"] is True
    assert not np.array_equal(result.geometry.material_levelsets[1], before[1])
    assert not np.array_equal(result.geometry.material_levelsets[2], before[2])

    with pytest.raises(ValueError, match="material-resolved mechanism router"):
        advance_feature_step_3d(
            geometry, boundary, {"Ar+": "energetic_bombardment"},
            _sputter("wrong-for-mask", "units", "product", 0.2),
            etchable_material_ids=(1, 2), duration_s=0.0,
            source_bounds=(-0.1, 1.1, -0.1, 0.3), source_z=1.8,
            ballistic_transport="face_gather", ballistic_face_quadrature_points=1,
            reinitialize=False, transport_device="cpu")


def test_material_extinction_requires_explicit_policy_and_retires_state(monkeypatch):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.2, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.1,
        substrate_top=1.0, etched_depth=0.2)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 39.948, 1e21, [[0.0, 0.0, 10.0]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1.8e-6)
    extinct = {
        1: np.asarray(geometry.material_levelsets[1]).copy(),
        2: np.full(geometry.phi.shape, -0.1),
    }
    extinct[2][1, 1, 1] = 2.0e-14

    monkeypatch.setattr(
        feature_step_module, "_advect_exposed_material_levelsets",
        lambda *args, **kwargs: {
            material_id: value.copy() for material_id, value in extinct.items()})
    monkeypatch.setattr(
        feature_step_module, "_suppress_interior_gas_nucleation",
        lambda previous_phi, previous_layers, phi, layers, **kwargs: (
            phi, layers, 0))

    arguments = dict(
        geometry=geometry, boundary=boundary,
        species_role={"Ar+": "energetic_bombardment"}, mechanism=_router(),
        etchable_material_ids=(1, 2), duration_s=0.1,
        source_bounds=(0.0, 1.0, 0.0, 0.2), source_z=1.8,
        ballistic_transport="face_gather", ballistic_periodic_lateral=True,
        ballistic_face_quadrature_points=1, profile_periodic_lateral=True,
        surface_state_remap_backend="indexed_knn", reinitialize=False,
        transport_device="cpu",
    )
    with pytest.raises(SurfaceTopologyChangeError) as info:
        advance_feature_step_3d(**arguments)
    assert info.value.event_kind == "material_component_change"

    result = advance_feature_step_3d(
        **arguments,
        topology_change_policy="continue_gas_cavity_and_material_extinction")

    assert set(result.geometry.material_levelsets) == {1}
    assert set(np.unique(result.geometry.material_id)) == {0, 1}
    assert all(name.startswith("m1__") for name in result.next_surface_state.fields)
    event = result.diagnostics["topology_event"]
    assert event["kind"] == "material_extinction"
    assert event["retired_material_ids"] == (2,)
    assert event["material_extinction_geometry"][2][
        "classification"] == "roundoff_projected_material_extinction"
    assert event["material_extinction_geometry"][2][
        "resolved_nonnegative_cell_count"] == 0
    retirement = result.state_remap_diagnostics["retired_materials"][2]
    assert retirement["old_face_count"] > 0
    assert retirement["retired_area_m2"] > 0.0
    assert retirement["lifecycle"] == "strict_levelset_extinction"
    assert result.state_remap_diagnostics[
        "requested_surface_state_remap_backend"] == "indexed_knn"
    assert result.state_remap_diagnostics[
        "surface_state_remap_backend"] == "common_refinement"

    continued = advance_feature_step_3d(
        result.geometry, boundary, {"Ar+": "energetic_bombardment"}, _router(),
        etchable_material_ids=(1, 2), duration_s=0.0,
        source_bounds=(0.0, 1.0, 0.0, 0.2), source_z=1.8,
        surface_state=result.next_surface_state,
        surface_state_mesh_fingerprint=result.next_surface_state_mesh_fingerprint,
        ballistic_transport="face_gather", ballistic_periodic_lateral=True,
        ballistic_face_quadrature_points=1, profile_periodic_lateral=True,
        surface_state_remap_backend="indexed_knn", reinitialize=False,
        transport_device="cpu",
        topology_change_policy="continue_gas_cavity_and_material_extinction")
    assert set(continued.next_surface_state.fields) == set(
        result.next_surface_state.fields)


def test_material_extinction_refuses_above_roundoff_hidden_support():
    layers = {1: np.ones((3, 3, 3)), 2: -np.ones((3, 3, 3))}
    layers[2][1, 1, 1] = 1.0e-8
    owner = np.ones((3, 3, 3), dtype=int)

    with pytest.raises(ValueError, match="above-roundoff nonnegative support"):
        _retire_strictly_extinguished_material_levelsets(layers, owner)


def test_common_engine_routes_material_specific_products_back_to_same_material():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.2, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.3,
        substrate_top=1.0, etched_depth=0.2)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 39.948, 1e21, [[0.0, 0.0, 10.0]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1.8e-6)

    def law(name, material_id):
        return SurfaceProductRedepositionLaw3D(
            name, material_id, {1: float(material_id == 1), 2: float(material_id == 2)},
            1e28,
            parameter_sources={
                "sticking_probability_by_material": "manufactured same-material gate",
                "bulk_material_unit_density_m3": "manufactured same-material gate"},
            parameter_bounds={
                "sticking_probability_by_material": (0.0, 1.0),
                "bulk_material_unit_density_m3": (0.9e28, 1.1e28)})

    result = advance_feature_step_3d(
        geometry, boundary, {"Ar+": "energetic_bombardment"}, _router(),
        etchable_material_ids=(1, 2), duration_s=1.0,
        source_bounds=(-0.1, 1.1, -0.1, 0.3), source_z=1.8,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
        surface_product_redeposition_options={
            "contract": SurfaceProductRedepositionContract3D((
                law("material_1:substrate_product", 1),
                law("material_2:mask_product", 2))),
            "rays_per_face": 8, "seed": 11},
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    assert result.surface_product_redeposition is not None
    assert result.diagnostics["product_redeposition_enabled"] is True
    assert result.diagnostics["product_redeposition_relative_balance_error"] < 1e-10
    assert result.diagnostics["max_growth_velocity_m_s"] >= 0.0
