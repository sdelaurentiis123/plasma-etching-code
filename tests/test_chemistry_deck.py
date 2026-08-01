"""Gates for the chemistry-deck framework (Deck 1: krueger_2024).

The deck must be a FAITHFUL extraction: mechanisms built from the deck have to
reproduce the legacy hardcoded-constant construction bit for bit, and the
validator has to reject incomplete, unprovenanced, or out-of-bounds decks.
"""

import copy

import numpy as np
import pytest

from petch.chemistry_deck import (
    KRUEGER_2024_DECK,
    DeckValidationError,
    build_mixed_layer_mechanisms_from_deck,
    engine_resident_drift,
    validate_deck,
)
from petch.mixed_layer import MixedLayerParams
from petch.mixed_layer_mechanism import (
    KRUEGER_2024_CHEMISORPTION_ACTIVATED,
    KRUEGER_2024_CHEMISORPTION_PROBABILITY,
    KRUEGER_2024_DEPOSITION_ON_CROSSLINKED,
    KRUEGER_2024_DEPOSITION_ON_MASK,
    KRUEGER_2024_DEPOSITION_ON_POLYMER,
    KRUEGER_2024_DEPOSITION_ON_SUBSTRATE,
    KRUEGER_2024_PRECURSOR_STOICHIOMETRY,
    MixedLayerMechanism,
)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


def _legacy_mechanisms(volatilization_yield=1.0):
    """The pre-deck construction, verbatim from the module-level maps."""
    oxide_parameters = MixedLayerParams(
        substrate="sio2", volatilization_yield=float(volatilization_yield))
    mask_parameters = MixedLayerParams(
        substrate="carbon", volatilization_yield=float(volatilization_yield))
    common = dict(
        precursor_species=dict(KRUEGER_2024_PRECURSOR_STOICHIOMETRY),
        fluorine_species=(),
        oxygen_species=("O",),
        inert_species=("C3F4",),
        chemisorption_probability=dict(KRUEGER_2024_CHEMISORPTION_PROBABILITY),
        chemisorption_activated_probability=dict(
            KRUEGER_2024_CHEMISORPTION_ACTIVATED),
        deposition_probability_on_film=dict(KRUEGER_2024_DEPOSITION_ON_POLYMER),
        deposition_probability_on_substrate=dict(
            KRUEGER_2024_DEPOSITION_ON_SUBSTRATE),
        deposition_probability_on_crosslinked=dict(
            KRUEGER_2024_DEPOSITION_ON_CROSSLINKED))
    mask_kwargs = dict(common)
    mask_kwargs["deposition_probability_on_substrate"] = dict(
        KRUEGER_2024_DEPOSITION_ON_MASK)
    return (MixedLayerMechanism(oxide_parameters, **common),
            MixedLayerMechanism(mask_parameters, **mask_kwargs))


def _mixed_condition_fluxes():
    """Three faces with genuinely different chemistry, ions as a per-event
    spectrum (so the atom-resolved path runs, not the scalar shortcut)."""
    ion = EnergeticFlux(
        name="Ar+",
        flux_m2_s=np.array([9.6e19, 4.1e19, 1.2e19]),
        energy_eV=np.array([1500.0, 600.0, 90.0]),
        cosine_incidence=np.array([1.0, 0.55, 0.08]),
        weight=np.array([0.5, 0.35, 0.15]))
    return SurfaceFluxes(
        neutral_flux_m2_s={
            "CF": np.array([4.4e20, 1.1e20, 0.0]),
            "CF2": np.array([9.4e20, 3.0e20, 5.0e19]),
            "CF3": np.array([8.4e19, 8.4e19, 8.4e19]),
            "C2F3": np.array([6.8e20, 0.0, 2.2e20]),
            "O": np.array([7.7e20, 2.0e20, 9.0e20]),
            "C3F4": np.array([9.5e20, 9.5e20, 9.5e20]),
        },
        energetic_fluxes=(ion,))


def _advance_signature(mechanism, duration_s=0.05):
    """Everything advance() returns that a caller can observe, as flat floats."""
    state = mechanism.initial_state((3,))
    fluxes = _mixed_condition_fluxes()
    result = mechanism.advance(state, fluxes, duration_s)
    signature = {}
    for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f",
                 "n_xl_film", "n_act", "removed_formula_units_m2"):
        signature[name] = np.asarray(getattr(result.state, name), dtype=float)
    for name in ("etch_velocity_m_s", "normal_growth_velocity_m_s",
                 "removed_bare_formula_units_m2", "deposited_polymer_units_m2",
                 "removed_polymer_units_m2"):
        signature[name] = np.asarray(getattr(result, name), dtype=float)
    return signature


def test_krueger_deck_validates():
    validate_deck(KRUEGER_2024_DECK)


def test_deck_rejects_missing_section():
    deck = copy.deepcopy(KRUEGER_2024_DECK)
    del deck["chemisorption"]
    with pytest.raises(DeckValidationError, match="missing sections"):
        validate_deck(deck)


def test_deck_rejects_empty_provenance():
    deck = copy.deepcopy(KRUEGER_2024_DECK)
    deck["deposition"]["provenance"] = "   "
    with pytest.raises(DeckValidationError, match="provenance"):
        validate_deck(deck)


def test_deck_rejects_out_of_bounds_probability():
    deck = copy.deepcopy(KRUEGER_2024_DECK)
    deck["chemisorption"]["bare"]["CF"] = 1.4
    with pytest.raises(DeckValidationError, match=r"out of \[0, 1\]"):
        validate_deck(deck)

    deck = copy.deepcopy(KRUEGER_2024_DECK)
    deck["sputter_laws"]["film_fresh"] = [0.9, 600.0, 500.0, 0.5]  # Eth > E0
    with pytest.raises(DeckValidationError, match="out of bounds"):
        validate_deck(deck)


def test_deck_rejects_undeclared_species():
    deck = copy.deepcopy(KRUEGER_2024_DECK)
    deck["deposition"]["on_polymer"]["C4F8"] = 0.1
    with pytest.raises(DeckValidationError, match="undeclared species"):
        validate_deck(deck)


def test_deck_build_matches_legacy_constants_bitwise():
    """The extraction is faithful: deck-built and legacy-built mechanisms
    produce identical advance() results to the last bit, on a 3-face vector
    of mixed conditions with a multi-row energetic spectrum."""
    deck_oxide, deck_mask = build_mixed_layer_mechanisms_from_deck()
    legacy_oxide, legacy_mask = _legacy_mechanisms()

    for deck_mechanism, legacy_mechanism in ((deck_oxide, legacy_oxide),
                                             (deck_mask, legacy_mask)):
        assert deck_mechanism.parameters == legacy_mechanism.parameters
        deck_signature = _advance_signature(deck_mechanism)
        legacy_signature = _advance_signature(legacy_mechanism)
        assert set(deck_signature) == set(legacy_signature)
        for name, expected in legacy_signature.items():
            actual = deck_signature[name]
            assert actual.shape == expected.shape, name
            assert np.array_equal(actual, expected), (
                f"{name}: {actual!r} != {expected!r}")


def test_deck_build_matches_legacy_probabilities():
    """Species maps land on the mechanism objects unchanged."""
    deck_oxide, deck_mask = build_mixed_layer_mechanisms_from_deck()
    legacy_oxide, legacy_mask = _legacy_mechanisms()
    for deck_mechanism, legacy_mechanism in ((deck_oxide, legacy_oxide),
                                             (deck_mask, legacy_mask)):
        for attribute in ("precursor_stoichiometry", "chemisorption_probability",
                          "chemisorption_activated_probability",
                          "deposition_probability_on_film",
                          "deposition_probability_on_substrate",
                          "deposition_probability_on_crosslinked"):
            assert (getattr(deck_mechanism, attribute)
                    == getattr(legacy_mechanism, attribute)), attribute
        for attribute in ("fluorine_species", "oxygen_species", "inert_species"):
            assert (getattr(deck_mechanism, attribute)
                    == getattr(legacy_mechanism, attribute)), attribute


def test_volatilization_override_flows_through():
    oxide, mask = build_mixed_layer_mechanisms_from_deck(
        volatilization_yield=0.5)
    assert oxide.parameters.volatilization_yield == 0.5
    assert mask.parameters.volatilization_yield == 0.5
    # And the deck itself is not mutated by the override.
    assert KRUEGER_2024_DECK["layer"]["volatilization_yield"] == 1.0


def test_explicit_parameters_override_deck_layer():
    custom = MixedLayerParams(substrate="sio2", oxidation_probability=0.01)
    oxide, _ = build_mixed_layer_mechanisms_from_deck(oxide_parameters=custom)
    assert oxide.parameters is custom


def test_public_wrapper_delegates_to_deck():
    from petch.mixed_layer_mechanism import (
        build_krueger_2024_mixed_layer_mechanisms)

    wrapper_oxide, wrapper_mask = build_krueger_2024_mixed_layer_mechanisms()
    legacy_oxide, legacy_mask = _legacy_mechanisms()
    for wrapper, legacy in ((wrapper_oxide, legacy_oxide),
                            (wrapper_mask, legacy_mask)):
        for name, expected in _advance_signature(legacy).items():
            assert np.array_equal(_advance_signature(wrapper)[name], expected), name


def test_engine_resident_constants_match_the_kernels():
    """Engine-resident deck entries are a declaration, not a fiction: they
    must agree with what the engine actually compiles in."""
    assert engine_resident_drift(KRUEGER_2024_DECK) == []
