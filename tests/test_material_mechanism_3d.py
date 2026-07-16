import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import advance_feature_step_3d, make_rectangular_trench_geometry_3d
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
