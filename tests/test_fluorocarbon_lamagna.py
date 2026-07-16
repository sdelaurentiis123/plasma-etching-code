import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import (
    advance_feature_step_3d, conservative_remap_surface_state,
    make_rectangular_trench_geometry_3d,
)
from petch.fluorocarbon_lamagna import (
    LaMagnaFluorocarbonParameters, LaMagnaFluorocarbonState,
    LaMagnaGarozzoFluorocarbonMechanism,
)
from petch.material_mechanism_3d import MaterialMechanismRouter3D
from petch.surface_kinetics import (
    EnergeticFlux, ParameterEvidence, SteinbruchelYield, SurfaceFluxes,
)


_INPUTS = {
    "bulk_formula_density_m3", "polymer_unit_density_m3",
    "substrate_etchant_sticking_probability",
    "substrate_polymer_sticking_probability",
    "polymer_etchant_sticking_probability",
    "polymer_polymer_sticking_probability", "reference_etchant_flux_m2_s",
    "polymer_loss_rate_m2_s", "temperature_K", "chemical_rate_coefficient",
    "chemical_activation_energy_eV", "ion_enhanced_coverage_loss_factor",
    "chemical_coverage_loss_factor", "physical_sputter_yield",
    "ion_enhanced_yield", "polymer_removal_yield",
}


def _mechanism(*, transport_mode="species_specific", predictive=True):
    evidence = {
        name: ParameterEvidence(
            "manufactured La Magna gate", "analytic",
            supports_prediction_within_declared_domain=predictive)
        for name in _INPUTS}
    return LaMagnaGarozzoFluorocarbonMechanism(LaMagnaFluorocarbonParameters(
        material_name="SiO2", material_inventory_name="SiO2_formula_unit",
        etchant_species=("etchant",), polymer_species=("polymer",),
        projectile_species=("ion",), bulk_formula_density_m3=2.0e28,
        polymer_unit_density_m3=1.0e28,
        substrate_etchant_sticking_probability=0.5,
        substrate_polymer_sticking_probability=0.25,
        polymer_etchant_sticking_probability=0.4,
        polymer_polymer_sticking_probability=0.2,
        reference_etchant_flux_m2_s=1.0e20, polymer_loss_rate_m2_s=1.0e18,
        temperature_K=300.0, chemical_rate_coefficient=1.0,
        chemical_activation_energy_eV=0.0,
        ion_enhanced_coverage_loss_factor=2.0,
        chemical_coverage_loss_factor=3.0,
        physical_sputter_yield=SteinbruchelYield(
            0.1, 0.0, angular_model="kress_1999", angular_parameter=0.0),
        ion_enhanced_yield=SteinbruchelYield(
            0.2, 0.0, angular_model="kress_1999", angular_parameter=0.0),
        polymer_removal_yield=SteinbruchelYield(
            0.05, 0.0, angular_model="kress_1999", angular_parameter=0.0),
        evidence=evidence, neutral_transport_mode=transport_mode,
    ))


def _ions(flux=1.0e18):
    return EnergeticFlux("ion", flux, [100.0], [1.0], [1.0])


def test_three_coverages_and_rates_match_the_published_algebraic_equations():
    mechanism = _mechanism()
    fluxes = SurfaceFluxes(
        {"etchant": 2.0e20, "polymer": 1.0e18}, (_ions(),))

    result = mechanism.advance(LaMagnaFluorocarbonState.bare(), fluxes, 2.0)

    ion_sp = 1.0e18 * 1.0
    ion_ie = 1.0e18 * 2.0
    ion_polymer = 1.0e18 * 0.5
    expected_pe = 2.0e20 * 0.4 / (2.0e20 * 0.4 + ion_polymer)
    expected_p = 1.0e18 * 0.25 / (ion_polymer * expected_pe + 1.0e18)
    expected_e = (
        2.0e20 * 0.5 * (1.0 - expected_p)
        / (2.0 * ion_ie + 3.0 * 1.0e20 + 2.0e20 * 0.5))
    expected_chemical = 1.0e20 * expected_e
    expected_enhanced = ion_ie * expected_e
    expected_sputter = ion_sp * (1.0 - expected_e)
    expected_rate = expected_chemical + expected_enhanced + expected_sputter

    assert np.isclose(result.etchant_on_polymer_coverage, expected_pe)
    assert np.isclose(result.polymer_coverage, expected_p)
    assert np.isclose(result.etchant_coverage, expected_e)
    assert np.isclose(result.chemical_removal_rate_m2_s, expected_chemical)
    assert np.isclose(result.ion_enhanced_removal_rate_m2_s, expected_enhanced)
    assert np.isclose(result.physical_sputter_rate_m2_s, expected_sputter)
    assert np.isclose(result.removed_formula_units_m2, 2.0 * expected_rate)
    assert np.isclose(result.etch_velocity_m_s, expected_rate / 2.0e28)
    assert abs(result.etchant_site_balance_residual_m2_s) / 2.0e20 < 2e-16
    assert abs(result.polymer_site_balance_residual_m2_s) / 1.0e18 < 2e-16
    assert abs(result.etchant_on_polymer_balance_residual_m2_s) / 1.0e20 < 2e-16
    assert result.material_exchange.residual_units_m2("SiO2_formula_unit") == 0.0


def test_ion_only_limit_retains_physical_sputtering_without_reactive_channels():
    result = _mechanism().advance(
        LaMagnaFluorocarbonState.bare(), SurfaceFluxes({}, (_ions(),)), 1.0)

    assert result.etchant_coverage == 0.0
    assert result.polymer_coverage == 0.0
    assert result.chemical_removal_rate_m2_s == 0.0
    assert result.ion_enhanced_removal_rate_m2_s == 0.0
    assert result.physical_sputter_rate_m2_s == 1.0e18
    assert result.etch_velocity_m_s > 0.0


def test_polymer_saturation_grows_a_finite_film_and_closes_the_deposition_ledger():
    mechanism = _mechanism()
    result = mechanism.advance(
        LaMagnaFluorocarbonState.bare(),
        SurfaceFluxes({"etchant": 0.0, "polymer": 5.0e18}), 4.0)

    expected_rate = 5.0e18 * 0.2
    assert result.polymer_coverage == 1.0
    assert result.etch_velocity_m_s == 0.0
    assert np.isclose(result.polymer_deposition_rate_m2_s, expected_rate)
    assert np.isclose(result.deposited_polymer_units_m2, 4.0 * expected_rate)
    assert np.isclose(result.state.polymer_film_units_m2, 4.0 * expected_rate)
    assert np.isclose(result.normal_growth_velocity_m_s, expected_rate / 1.0e28)
    assert result.material_exchange.residual_units_m2("fluorocarbon_film_unit") == 0.0


def test_existing_polymer_is_removed_before_substrate_recession():
    mechanism = _mechanism()
    initial = LaMagnaFluorocarbonState(
        0.0, 0.0, 0.0, polymer_film_units_m2=1.0e17)
    result = mechanism.advance(
        initial, SurfaceFluxes({"etchant": 2.0e20, "polymer": 0.0}, (_ions(),)), 1.0)

    assert result.state.polymer_film_units_m2 == 0.0
    assert result.removed_polymer_units_m2 == initial.polymer_film_units_m2
    assert result.removed_formula_units_m2 > 0.0
    assert result.etch_velocity_m_s > (
        result.removed_polymer_units_m2 / mechanism.parameters.polymer_unit_density_m3)
    assert result.material_exchange.residual_units_m2("fluorocarbon_film_unit") == 0.0


def test_neutral_transport_can_reproduce_viennaps_or_use_species_specific_sticking():
    state = LaMagnaFluorocarbonState(0.1, 0.2, 0.3, 0.0, 0.0)
    species_specific = _mechanism(transport_mode="species_specific")
    source_compatible = _mechanism(transport_mode="viennaps_4_6_1")

    specific = species_specific.neutral_reaction_probability(state)
    compatible = source_compatible.neutral_reaction_probability(state)

    assert np.isclose(specific["etchant"], 0.5 * 0.7)
    assert np.isclose(specific["polymer"], 0.25 * 0.7)
    assert np.isclose(compatible["etchant"], 0.5 * 0.7)
    assert np.isclose(compatible["polymer"], 0.5 * 0.7)


def test_reference_factory_records_exact_source_commit_and_does_not_claim_prediction():
    mechanism = LaMagnaGarozzoFluorocarbonMechanism(
        LaMagnaFluorocarbonParameters.viennaps_4_6_1_reference(
            reference_etchant_flux_m2_s=5.0e21))

    assert mechanism.provenance["viennaps_source_commit"] == (
        "2956ed587984c6dc38be24c6e2390e10c9b2f0a7")
    assert mechanism.provenance["neutral_transport_mode"] == "viennaps_4_6_1"
    parameters = mechanism.parameters
    energy = 100.0; cosine = 0.5
    expected_sputter = (
        0.0139 * (np.sqrt(energy) - np.sqrt(18.0))
        * (1.0 + 9.3 * (1.0 - cosine ** 2)) * cosine)
    expected_enhanced = 0.0361 * (np.sqrt(energy) - np.sqrt(4.0)) * cosine
    expected_polymer = 0.0722 * (np.sqrt(energy) - np.sqrt(4.0)) * cosine
    assert np.isclose(parameters.physical_sputter_yield.evaluate(energy, cosine), expected_sputter)
    assert np.isclose(parameters.ion_enhanced_yield.evaluate(energy, cosine), expected_enhanced)
    assert np.isclose(parameters.polymer_removal_yield.evaluate(energy, cosine), expected_polymer)
    result = mechanism.advance(
        mechanism.initial_state(), SurfaceFluxes({"FC_etchant": 1.0e20}), 0.0)
    assert not result.validity.parameter_evidence_supports_prediction
    assert set(result.validity.nonpredictive_parameters) == _INPUTS


def test_mechanism_refuses_positive_undeclared_species():
    with pytest.raises(ValueError, match="no declared fluorocarbon reaction channel"):
        _mechanism().advance(
            LaMagnaFluorocarbonState.bare(), SurfaceFluxes({"photon": 1.0e18}), 1.0)


def test_material_router_preserves_mechanism_owned_polymer_growth_velocity():
    router = MaterialMechanismRouter3D(
        {1: _mechanism(), 2: _mechanism()},
        provenance={1: "manufactured oxide", 2: "manufactured oxide copy"})
    material = np.array([1, 2])
    state = router.initial_state_by_material(material)
    assert state.surface_field_remap_modes()["m1__polymer_coverage"] == "intensive"

    result = router.advance_by_material(
        state, SurfaceFluxes({"polymer": np.full(2, 5.0e18)}), 1.0, material)

    assert np.all(result.etch_velocity_m_s == 0.0)
    assert np.allclose(result.normal_growth_velocity_m_s, 1.0e-10)
    assert result.state.surface_field_remap_modes()[
        "m2__polymer_film_units_m2"] == "conservative"
    assert np.all(result.material_exchange.deposited_units_m2[
        "fluorocarbon_film_unit"] > 0.0)


def test_common_feature_engine_advects_mechanism_growth_on_the_shared_level_set():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.2, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.3,
        substrate_top=1.0, etched_depth=0.2)
    polymer = SpeciesBoundaryState(
        "polymer", 0, 100.0, 1.0e22, [[0.0, 0.0, 1.0]], [1.0])
    boundary = PlasmaBoundaryState((polymer,), reference_plane_m=1.8e-6)

    result = advance_feature_step_3d(
        geometry, boundary, {"polymer": "neutral_reactant"}, _mechanism(),
        etchable_material_ids=(1,), duration_s=0.1,
        source_bounds=(-0.1, 1.1, -0.1, 0.3), source_z=1.8,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=1,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    active_velocity = result.face_velocity_mesh_units_s[result.active_face_index]
    assert np.any(active_velocity < 0.0)
    assert result.diagnostics["max_surface_mechanism_growth_velocity_m_s"] > 0.0
    assert result.diagnostics["max_growth_velocity_m_s"] > 0.0
    assert np.any(result.surface.state.polymer_film_units_m2 > 0.0)
    assert not np.array_equal(result.geometry.phi, geometry.phi)
    assert np.all(result.surface.material_exchange.residual_units_m2(
        "fluorocarbon_film_unit") == 0.0)


def test_moving_surface_interpolates_coverages_but_conserves_film_inventory():
    state = LaMagnaFluorocarbonState(
        [1.0, 1.0], [1.0, 1.0], [0.2, 0.4], [2.0, 4.0], [5.0, 7.0])
    old_centroid = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    new_centroid = np.array([[0.1, 0.0, 0.0], [0.9, 0.0, 0.0]])
    old_area = np.array([2.0, 2.0])
    new_area = np.array([1.0, 1.0])

    remapped, diagnostics = conservative_remap_surface_state(
        state, old_centroid, old_area, np.ones(2, dtype=int),
        new_centroid, new_area, np.ones(2, dtype=int),
        dx=1.0, mesh_length_unit_m=1.0)

    assert np.all(remapped.etchant_coverage == 1.0)
    assert np.all(remapped.polymer_coverage == 1.0)
    assert np.isclose(
        np.dot(remapped.polymer_film_units_m2, new_area),
        np.dot(state.polymer_film_units_m2, old_area))
    material = diagnostics["materials"][1]
    assert material["field_remap_modes"]["polymer_coverage"] == "intensive"
    assert material["field_remap_modes"]["polymer_film_units_m2"] == "conservative"
