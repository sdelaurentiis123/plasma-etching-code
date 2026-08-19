import pytest

from petch import (
    EnergeticFlux,
    EnergeticYield,
    ParameterEvidence,
    SurfaceFluxes,
    TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE,
    Tio2ReducedSurfaceDeck,
)
from petch.tio2_ion_dose import tio2_formula_unit_density_m3


def _evidence(*, supports=False):
    return {
        name: ParameterEvidence(
            "manufactured TiO2 deck plumbing test",
            "manufactured",
            supports_prediction_within_declared_domain=supports,
        )
        for name in TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE
    }


def _deck(*, evidence=None):
    return Tio2ReducedSurfaceDeck(
        mass_density_kg_m3=3600.0,
        site_density_m2=1.0e19,
        passivation_monolayer_density_m2=1.2e19,
        passivation_bulk_unit_density_m3=7.0e28,
        fluorination_probability={"F": 0.2, "CF2": 0.1},
        passivation_deposition_probability_on_tio2={"CF2": 0.01},
        passivation_deposition_probability_on_passivation={"CF2": 0.1},
        oxygen_species="O",
        oxygen_passivation_removal_probability=0.05,
        oxygen_blocking_probability=0.04,
        oxygen_blocker_ion_removal_yield=EnergeticYield(0.1, 10.0, 150.0),
        bare_tio2_yield=EnergeticYield(0.2, 30.0, 150.0),
        fluorinated_tio2_yield=EnergeticYield(0.8, 15.0, 150.0),
        passivation_sputter_yield=EnergeticYield(0.4, 20.0, 150.0),
        evidence=_evidence() if evidence is None else evidence,
    )


def test_tio2_deck_refuses_missing_evidence_and_unresolved_model_form():
    missing = _evidence()
    del missing["complex_sio2_yield"]
    deck = _deck(evidence=missing)
    with pytest.raises(ValueError, match="complex_sio2_yield"):
        deck.build_parameters(allow_reduced_sensitivity=True)

    deck = _deck()
    assert deck.readiness().supports_reduced_sensitivity is True
    assert deck.readiness().supports_absolute_target_prediction is False
    with pytest.raises(ValueError, match="roughness"):
        deck.build_parameters()


def test_tio2_reduced_sensitivity_deck_is_material_specific_and_nonpredictive():
    deck = _deck()
    mechanism = deck.build_mechanism(allow_reduced_sensitivity=True)
    parameters = mechanism.parameters

    assert parameters.material_name == "ALD TiO2"
    assert parameters.material_inventory_name == "TiO2_formula_unit"
    assert parameters.bulk_formula_density_m3 == pytest.approx(
        tio2_formula_unit_density_m3(3600.0)
    )
    assert "chemistry_dependent_roughness_evolution" in parameters.known_omissions

    result = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes(
            {"F": 1.0e19, "CF2": 1.0e19, "O": 1.0e19},
            (EnergeticFlux("Ar+", 1.0e18, [150.0], [1.0], [1.0]),),
        ),
        0.01,
    )
    assert "TiO2_formula_unit" in result.material_exchange.removed_units_m2
    assert result.validity.parameter_evidence_supports_prediction is False


def test_even_predictive_parameter_evidence_cannot_override_model_form_gap():
    deck = _deck(evidence=_evidence(supports=True))

    assert deck.readiness().nonpredictive_parameter_evidence == ()
    assert deck.readiness().unresolved_model_form
    assert deck.readiness().supports_absolute_target_prediction is False


def test_oxygen_competition_is_bounded_conservative_and_suppresses_fluorination():
    mechanism = _deck().build_mechanism(allow_reduced_sensitivity=True)
    ion = EnergeticFlux("Ar+", 2.0e17, [150.0], [1.0], [1.0])
    without_oxygen = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes({"F": 2.0e19, "CF2": 1.0e19}, (ion,)),
        1.0,
    )
    with_oxygen = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes(
            {"F": 2.0e19, "CF2": 1.0e19, "O": 5.0e20}, (ion,)
        ),
        1.0,
        max_step_s=0.02,
    )

    assert with_oxygen.state.oxygen_blocked_fraction > 0.0
    assert with_oxygen.state.complex_fraction < without_oxygen.state.complex_fraction
    assert (
        with_oxygen.state.complex_fraction
        + with_oxygen.state.oxygen_blocked_fraction
        <= 1.0
    )
    blocker_change = (
        with_oxygen.state.oxygen_blocked_fraction
        * mechanism.parameters.site_density_m2
    )
    assert blocker_change == pytest.approx(
        with_oxygen.formed_oxygen_blocked_sites_m2
        - with_oxygen.removed_oxygen_blocked_sites_m2,
        rel=2.0e-11,
    )
