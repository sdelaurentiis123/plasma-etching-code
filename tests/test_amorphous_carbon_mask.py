import numpy as np
import pytest

from petch.amorphous_carbon_mask import (
    AmorphousCarbonMaskMechanism,
    AmorphousCarbonMaskParameters,
    AmorphousCarbonMaskState,
    build_krueger_2024_material_router_3d,
)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


def _ions(flux=1.2e20, energy=500.0):
    return EnergeticFlux("ions", flux, [energy], [1.0], [1.0])


def _mechanism():
    return AmorphousCarbonMaskMechanism(
        AmorphousCarbonMaskParameters.krueger_2024_reduced_projection())


def test_krueger_mask_projection_preserves_calibrated_and_appendix_b_values():
    parameters = AmorphousCarbonMaskParameters.krueger_2024_reduced_projection()
    mechanism = AmorphousCarbonMaskMechanism(parameters)

    assert parameters.polymer_deposition_probability_on_carbon == {
        "CF": 0.0842, "CF2": 0.0842, "CF3": 0.0842, "C2F3": 0.0842}
    assert parameters.polymer_deposition_probability_on_polymer["CF"] == 0.1
    assert parameters.polymer_deposition_probability_on_polymer["C2F3"] == 0.03
    assert parameters.polymer_deposition_probability_on_crosslinked_polymer == {
        "CF": 0.02, "CF2": 0.02, "CF3": 0.02, "C2F3": 0.02}
    assert parameters.effective_crosslinked_growth_fraction == 0.0
    assert parameters.oxygen_polymer_etch_probability == 0.0628
    assert parameters.oxygen_carbon_etch_probability == 1.0e-5
    assert parameters.polymer_sputter_yield.reference_yield == 0.9
    assert parameters.polymer_sputter_yield.threshold_energy_eV == 20.0
    assert parameters.polymer_sputter_yield.energy_exponent == 0.5
    assert parameters.carbon_sputter_yield.reference_yield == 0.001
    assert parameters.carbon_sputter_yield.threshold_energy_eV == 200.0
    assert parameters.carbon_sputter_yield.energy_exponent == 0.4
    assert parameters.declared_inert_neutral_species == ("C3F4",)
    assert mechanism.provenance["parameters"][
        "polymer_deposition_probability_on_carbon"]["CF2"] == 0.0842
    assert mechanism.neutral_reaction_probability(
        AmorphousCarbonMaskState.bare())["C3F4"] == 0.0
    assert not mechanism.validity(SurfaceFluxes({})).parameter_evidence_supports_prediction


def test_effective_crosslinked_growth_blends_published_collision_laws_only():
    fresh = AmorphousCarbonMaskMechanism(
        AmorphousCarbonMaskParameters.krueger_2024_reduced_projection(
            effective_crosslinked_growth_fraction=0.0))
    crosslinked = AmorphousCarbonMaskMechanism(
        AmorphousCarbonMaskParameters.krueger_2024_reduced_projection(
            effective_crosslinked_growth_fraction=1.0))
    monolayer = fresh.parameters.polymer_monolayer_density_m2
    state = AmorphousCarbonMaskState(5.0 * monolayer, 0.0)
    fluxes = SurfaceFluxes({"CF": 1.0e20})

    fresh_step = fresh.advance(state, fluxes, 0.0)
    crosslinked_step = crosslinked.advance(state, fluxes, 0.0)

    assert crosslinked_step.normal_growth_velocity_m_s < (
        fresh_step.normal_growth_velocity_m_s)
    assert np.isclose(
        crosslinked._effective_polymer_growth_probability["CF"], 0.02)
    assert crosslinked.provenance["parameters"][
        "effective_crosslinked_growth_fraction"] == 1.0
    with pytest.raises(ValueError, match="invalid amorphous-carbon mask parameters"):
        AmorphousCarbonMaskParameters.krueger_2024_reduced_projection(
            effective_crosslinked_growth_fraction=1.01)


def test_c3f4_is_explicitly_nonreactive_without_being_dropped_or_refused():
    mechanism = _mechanism()
    initial = AmorphousCarbonMaskState.bare()
    fluxes = SurfaceFluxes({"C3F4": 9.5e20})

    probability = mechanism.neutral_reaction_probability(initial)
    result = mechanism.advance(initial, fluxes, 60.0)

    assert probability["C3F4"] == 0.0
    assert result.validity.within_declared_scope
    assert result.validity.unsupported_neutral_species == ()
    assert result.state.polymer_units_m2 == 0.0
    assert result.state.removed_carbon_atoms_m2 == 0.0


def test_mask_film_inventory_and_material_ledgers_close_exactly():
    mechanism = _mechanism()
    monolayer = mechanism.parameters.polymer_monolayer_density_m2
    initial = AmorphousCarbonMaskState(0.7 * monolayer, 3.0e16)
    fluxes = SurfaceFluxes(
        {"CF2": 9.4e20, "O": 7.7e20}, (_ions(1.2e20, 500.0),))

    result = mechanism.advance(initial, fluxes, 0.2)

    film_change = result.state.polymer_units_m2 - initial.polymer_units_m2
    assert np.isclose(
        film_change,
        result.deposited_polymer_units_m2 - result.removed_polymer_units_m2,
        rtol=5e-13, atol=128.0)
    assert np.isclose(
        result.state.removed_carbon_atoms_m2 - initial.removed_carbon_atoms_m2,
        result.removed_carbon_atoms_m2,
        rtol=2e-15, atol=16.0)
    assert np.all(result.material_exchange.residual_units_m2(
        "amorphous_carbon_atom") == 0.0)
    assert np.all(result.material_exchange.residual_units_m2(
        "fluorocarbon_film_unit") == 0.0)
    assert np.array_equal(
        result.material_exchange.deposited_units_m2["fluorocarbon_film_unit"],
        result.deposited_polymer_units_m2)


def test_mask_causality_oxygen_cleans_film_and_energy_exposes_carbon_erosion():
    mechanism = _mechanism()
    monolayer = mechanism.parameters.polymer_monolayer_density_m2
    film = AmorphousCarbonMaskState(2.0 * monolayer, 0.0)
    no_oxygen = mechanism.advance(
        film, SurfaceFluxes({}, (_ions(2.0e19, 100.0),)), 0.05)
    oxygen = mechanism.advance(
        film, SurfaceFluxes({"O": 2.0e21}, (_ions(2.0e19, 100.0),)), 0.05)
    low_energy = mechanism.advance(
        AmorphousCarbonMaskState.bare(),
        SurfaceFluxes({}, (_ions(2.0e20, 100.0),)), 0.05)
    high_energy = mechanism.advance(
        AmorphousCarbonMaskState.bare(),
        SurfaceFluxes({}, (_ions(2.0e20, 500.0),)), 0.05)

    assert oxygen.state.polymer_units_m2 < no_oxygen.state.polymer_units_m2
    assert oxygen.removed_polymer_units_m2 > no_oxygen.removed_polymer_units_m2
    assert low_energy.removed_carbon_atoms_m2 == 0.0
    assert high_energy.removed_carbon_atoms_m2 > 0.0
    assert high_energy.etch_velocity_m_s > low_energy.etch_velocity_m_s


def test_mask_refuses_undeclared_positive_incident_species():
    mechanism = _mechanism()
    invalid = SurfaceFluxes(
        {"mystery": 1e18},
        (EnergeticFlux("mystery+", 1e18, [500.0], [1.0], [1.0]),))

    validity = mechanism.validity(invalid)

    assert not validity.within_declared_scope
    assert validity.unsupported_neutral_species == ("mystery",)
    assert len(validity.reasons) == 2
    with pytest.raises(ValueError, match="outside declared scope"):
        mechanism.advance(AmorphousCarbonMaskState.bare(), invalid, 1.0)


def test_krueger_router_dispatches_oxide_and_mask_without_forking_the_engine():
    router = build_krueger_2024_material_router_3d()
    material = np.array([1, 1, 2, 2])
    state = router.initial_state_by_material(material)
    fluxes = SurfaceFluxes(
        {
            "C3F4": np.full(4, 9.5e20),
            "C2F3": np.full(4, 6.8e20),
            "CF": np.full(4, 4.4e20),
            "CF2": np.full(4, 9.4e20),
            "CF3": np.full(4, 8.4e19),
            "O": np.full(4, 7.7e20),
        },
        (EnergeticFlux(
            "ions", np.full(4, 1.2e20), [500.0], [1.0], [1.0]),))

    result = router.advance_by_material(state, fluxes, 0.01, material)

    assert set(result.material_results) == {1, 2}
    assert result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert np.all(result.etch_velocity_m_s[:2] > 0.0)
    assert np.all(result.etch_velocity_m_s[2:] > 0.0)
    assert np.all(result.normal_growth_velocity_m_s[2:] > 0.0)
    assert set(result.material_exchange.removed_units_m2) == {
        "SiO2_formula_unit", "fluorocarbon_film_unit", "amorphous_carbon_atom"}
    assert router.provenance["materials"]["1"]["evidence"]["role"] == "SiO2 substrate"
    assert router.provenance["materials"]["2"]["evidence"]["role"] == (
        "amorphous-carbon mask")

    scaled = build_krueger_2024_material_router_3d(
        oxide_etch_yield_scale=0.5)
    assert np.isclose(
        scaled.mechanisms[1].parameters.bare_sio2_yield.reference_yield,
        0.5 * router.mechanisms[1].parameters.bare_sio2_yield.reference_yield)
    assert scaled.mechanisms[2].parameters.carbon_sputter_yield == (
        router.mechanisms[2].parameters.carbon_sputter_yield)


def test_guo_tml_router_reuses_transport_and_mask_but_carries_transfer_limits():
    router = build_krueger_2024_material_router_3d(
        surface_model="guo_tml")
    material = np.array([1, 2])
    state = router.initial_state_by_material(material)
    fluxes = SurfaceFluxes(
        {
            "C3F4": np.full(2, 9.5e20),
            "C2F3": np.full(2, 6.8e20),
            "CF": np.full(2, 4.4e20),
            "CF2": np.full(2, 9.4e20),
            "CF3": np.full(2, 8.4e19),
            "O": np.full(2, 7.7e20),
        },
        (EnergeticFlux(
            "ions", np.full(2, 1.2e20),
            [3500.0], [1.0], [1.0]),),
    )

    result = router.advance_by_material(state, fluxes, 0.01, material)

    assert result.etch_velocity_m_s[0] > 0.0
    assert result.etch_velocity_m_s[1] > 0.0
    assert not result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert router.provenance["materials"]["1"]["evidence"][
        "claim_status"] == "out_of_board_transfer_audit_not_prediction"
    with pytest.raises(ValueError, match="no oxide yield scale"):
        build_krueger_2024_material_router_3d(
            surface_model="guo_tml", oxide_etch_yield_scale=1.01)

    cf3_endpoint = build_krueger_2024_material_router_3d(
        surface_model="guo_tml",
        guo_aggregate_ion_formula="CF3",
    )
    assert cf3_endpoint.mechanisms[1].ion_species_mapping == {"ions": "CF3"}
    assert cf3_endpoint.provenance["materials"]["1"]["evidence"]["model"][
        "ion_species_mapping"] == {"ions": "CF3"}

    finite_layer = build_krueger_2024_material_router_3d(
        surface_model="guo_tml",
        guo_translating_layer_thickness_nm=1.2,
    )
    assert finite_layer.mechanisms[
        1].translating_layer_thickness_nm == pytest.approx(1.2)
    with pytest.raises(ValueError, match="applies only"):
        build_krueger_2024_material_router_3d(
            guo_translating_layer_thickness_nm=1.2)
