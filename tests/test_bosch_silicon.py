import numpy as np
import pytest

from petch.bosch_silicon import (
    BoschSiliconFluorocarbonMechanism, BoschSiliconFluorocarbonState,
)
from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import (
    advance_feature_step_3d, make_rectangular_trench_geometry_3d,
)
from petch.fluorocarbon_lamagna import (
    LaMagnaFluorocarbonParameters, LaMagnaGarozzoFluorocarbonMechanism,
)
from petch.silicon_sf6o2 import (
    BelenSiliconParameters, BelenSiliconSF6O2Mechanism,
)
from petch.surface_kinetics import (
    EnergeticFlux, ParameterEvidence, SteinbruchelYield, SurfaceFluxes,
)


_BELEN_INPUTS = {
    "site_density_m2", "bulk_si_atom_density_m3",
    "fluorine_sticking_probability", "oxygen_sticking_probability",
    "spontaneous_fluorine_removal_rate_m2_s",
    "oxygen_desorption_rate_m2_s", "physical_sputter_yield",
    "ion_enhanced_yield", "oxygen_sputter_yield",
    "fluorine_atoms_per_removed_si",
    "ion_enhanced_fluorine_release_per_si",
}
_FILM_INPUTS = {
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


def _silicon(*, predictive=True):
    evidence = {
        name: ParameterEvidence(
            "manufactured Bosch/Belen gate", "analytic",
            supports_prediction_within_declared_domain=predictive)
        for name in _BELEN_INPUTS}
    bounds = {name: (0.0, 1.0e30) for name in _BELEN_INPUTS}
    bounds.update({
        "physical_sputter_yield": {
            "prefactor_per_sqrt_eV": (0.0, 1.0),
            "threshold_energy_eV": (0.0, 100.0)},
        "ion_enhanced_yield": {
            "prefactor_per_sqrt_eV": (0.0, 10.0),
            "threshold_energy_eV": (0.0, 100.0)},
        "oxygen_sputter_yield": {
            "prefactor_per_sqrt_eV": (0.0, 10.0),
            "threshold_energy_eV": (0.0, 100.0)},
    })
    return BelenSiliconSF6O2Mechanism(BelenSiliconParameters(
        material_name="Si", material_inventory_name="Si_atom",
        fluorine_species="F", oxygen_species="O", projectile_species=("ion",),
        site_density_m2=5.0e18, bulk_si_atom_density_m3=5.0e28,
        fluorine_sticking_probability=0.5, oxygen_sticking_probability=0.25,
        spontaneous_fluorine_removal_rate_m2_s=4.0e19,
        oxygen_desorption_rate_m2_s=2.0e19,
        physical_sputter_yield=SteinbruchelYield(0.1, 25.0),
        ion_enhanced_yield=SteinbruchelYield(0.2, 25.0),
        oxygen_sputter_yield=SteinbruchelYield(0.3, 25.0),
        fluorine_atoms_per_removed_si=4.0,
        ion_enhanced_fluorine_release_per_si=2.0,
        evidence=evidence, parameter_bounds=bounds))


def _film(*, predictive=True):
    evidence = {
        name: ParameterEvidence(
            "manufactured Bosch/film gate", "analytic",
            supports_prediction_within_declared_domain=predictive)
        for name in _FILM_INPUTS}
    return LaMagnaGarozzoFluorocarbonMechanism(
        LaMagnaFluorocarbonParameters(
            material_name="Si film-clock dummy",
            material_inventory_name="discarded_dummy_substrate_unit",
            etchant_species=("F",), polymer_species=("C4F8",),
            projectile_species=("ion",), bulk_formula_density_m3=5.0e28,
            polymer_unit_density_m3=1.0e28,
            substrate_etchant_sticking_probability=0.5,
            substrate_polymer_sticking_probability=0.25,
            polymer_etchant_sticking_probability=0.4,
            polymer_polymer_sticking_probability=0.2,
            reference_etchant_flux_m2_s=1.0e20,
            polymer_loss_rate_m2_s=1.0e18, temperature_K=300.0,
            chemical_rate_coefficient=1.0,
            chemical_activation_energy_eV=0.0,
            ion_enhanced_coverage_loss_factor=2.0,
            chemical_coverage_loss_factor=3.0,
            physical_sputter_yield=SteinbruchelYield(0.1, 0.0),
            ion_enhanced_yield=SteinbruchelYield(0.2, 0.0),
            polymer_removal_yield=SteinbruchelYield(0.05, 0.0),
            evidence=evidence, neutral_transport_mode="species_specific"))


def _mechanism(*, predictive=True):
    return BoschSiliconFluorocarbonMechanism(
        _silicon(predictive=predictive), _film(predictive=predictive))


def _ions(flux=1.0e18):
    return EnergeticFlux("ion", flux, [100.0], [1.0], [1.0])


def test_film_free_no_precursor_limit_is_exactly_the_unchanged_belen_law():
    silicon = _silicon()
    mechanism = BoschSiliconFluorocarbonMechanism(silicon, _film())
    fluxes = SurfaceFluxes({"F": 2.0e20, "O": 1.0e19}, (_ions(),))

    expected = silicon.advance(silicon.initial_state(), fluxes, 2.0)
    actual = mechanism.advance(mechanism.initial_state(), fluxes, 2.0)

    assert actual.substrate_exposure_fraction == 1.0
    assert actual.etch_velocity_m_s == expected.etch_velocity_m_s
    assert actual.removed_si_atoms_m2 == expected.removed_si_atoms_m2
    assert actual.fluorine_coverage == expected.fluorine_coverage
    assert actual.oxygen_coverage == expected.oxygen_coverage
    assert actual.state.removed_si_atoms_m2 == expected.state.removed_si_atoms_m2
    assert actual.state.polymer_film_units_m2 == 0.0


def test_c4f8_only_deposits_finite_film_without_removing_silicon():
    mechanism = _mechanism()
    result = mechanism.advance(
        mechanism.initial_state(), SurfaceFluxes({"C4F8": 5.0e18}), 4.0)

    assert result.substrate_exposure_fraction == 0.0
    assert result.removed_si_atoms_m2 == 0.0
    assert result.etch_velocity_m_s == 0.0
    assert result.deposited_polymer_units_m2 > 0.0
    assert result.normal_growth_velocity_m_s > 0.0
    assert result.state.polymer_film_units_m2 == result.deposited_polymer_units_m2
    assert result.material_exchange.residual_units_m2("Si_atom") == 0.0
    assert result.material_exchange.residual_units_m2(
        "fluorocarbon_film_unit") == 0.0


def test_existing_film_is_removed_before_belen_silicon_recession():
    mechanism = _mechanism()
    initial = BoschSiliconFluorocarbonState(
        polymer_film_units_m2=1.0e17)
    fluxes = SurfaceFluxes({"F": 2.0e20}, (_ions(),))

    result = mechanism.advance(initial, fluxes, 1.0)
    bare = _silicon().advance(_silicon().initial_state(), fluxes, 1.0)

    assert result.state.polymer_film_units_m2 == 0.0
    assert result.removed_polymer_units_m2 == initial.polymer_film_units_m2
    assert 0.0 < result.substrate_exposure_fraction < 1.0
    assert result.removed_si_atoms_m2 == pytest.approx(
        bare.removed_si_atoms_m2 * result.substrate_exposure_fraction)
    assert result.removed_si_atoms_m2 > 0.0


def test_thick_film_blocks_silicon_until_a_later_step():
    mechanism = _mechanism()
    initial = BoschSiliconFluorocarbonState(
        polymer_film_units_m2=1.0e20)
    result = mechanism.advance(
        initial, SurfaceFluxes({"F": 2.0e20}, (_ions(),)), 0.01)

    assert result.substrate_exposure_fraction == 0.0
    assert result.removed_si_atoms_m2 == 0.0
    assert result.state.polymer_film_units_m2 > 0.0


def test_transport_sink_switches_from_belen_silicon_to_film_surface():
    mechanism = _mechanism()
    bare = mechanism.neutral_reaction_probability(mechanism.initial_state())
    coated = mechanism.neutral_reaction_probability(
        BoschSiliconFluorocarbonState(polymer_film_units_m2=1.0e18))

    assert bare["F"] == 0.5
    assert coated["F"] == 0.4
    assert bare["C4F8"] == 0.25
    assert coated["C4F8"] == 0.2
    assert coated["O"] == 0.0


def test_state_declares_intensive_coverages_and_conserved_inventories():
    state = BoschSiliconFluorocarbonState.bare((2,))
    modes = state.surface_field_remap_modes()

    assert modes["available_site_fraction"] == "intensive"
    assert modes["polymer_coverage"] == "intensive"
    assert modes["removed_si_atoms_m2"] == "conservative"
    assert modes["polymer_film_units_m2"] == "conservative"
    restored = state.with_conservative_surface_fields(
        state.conservative_surface_fields())
    for name, expected in state.conservative_surface_fields().items():
        assert np.array_equal(restored.conservative_surface_fields()[name], expected)


def test_composite_validity_namespaces_parameter_evidence_and_refuses_unknown_flux():
    mechanism = _mechanism(predictive=False)
    validity = mechanism.advance(
        mechanism.initial_state(), SurfaceFluxes({"F": 1.0e18}), 0.0
    ).validity

    assert not validity.parameter_evidence_supports_prediction
    assert any(name.startswith("silicon.") for name in validity.nonpredictive_parameters)
    assert any(name.startswith("film.") for name in validity.nonpredictive_parameters)
    with pytest.raises(ValueError, match="no declared Bosch reaction channel"):
        mechanism.advance(
            mechanism.initial_state(), SurfaceFluxes({"photon": 1.0e18}), 1.0)


def test_common_feature_engine_advects_bosch_film_growth_without_a_second_solver():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.2, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.3,
        substrate_top=1.0, etched_depth=0.2)
    precursor = SpeciesBoundaryState(
        "C4F8", 0, 200.0, 1.0e22, [[0.0, 0.0, 1.0]], [1.0])
    boundary = PlasmaBoundaryState((precursor,), reference_plane_m=1.8e-6)

    result = advance_feature_step_3d(
        geometry, boundary, {"C4F8": "neutral_reactant"}, _mechanism(),
        etchable_material_ids=(1,), duration_s=0.1,
        source_bounds=(-0.1, 1.1, -0.1, 0.3), source_z=1.8,
        ballistic_transport="face_gather", ballistic_face_quadrature_points=1,
        cfl_number=0.3, reinitialize=False, transport_device="cpu")

    assert result.diagnostics["max_surface_mechanism_growth_velocity_m_s"] > 0.0
    assert np.any(result.surface.state.polymer_film_units_m2 > 0.0)
    assert not np.array_equal(result.geometry.phi, geometry.phi)
    assert np.all(result.surface.material_exchange.residual_units_m2(
        "fluorocarbon_film_unit") == 0.0)
