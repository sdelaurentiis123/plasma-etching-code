import numpy as np

from petch.reactor_global import (
    ARGON_4S_METASTABLE_ENERGY_EV,
    ARGON_IONIZATION_ENERGY_EV,
    ARGON_METASTABLE_IONIZATION_ENERGY_EV,
    RateContext,
    build_lee_lieberman_argon_volume_network,
)


def test_argon_deck_imports_table3_rate_coefficients_in_si():
    network = build_lee_lieberman_argon_volume_network()
    context = RateContext(electron_temperature_eV=3.5)
    coefficients = np.array([
        reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
    ])
    expected = np.array([
        1.23e-7 * np.exp(-18.68 / 3.5),
        3.71e-8 * np.exp(-15.06 / 3.5),
        2.05e-7 * np.exp(-4.95 / 3.5),
        2.0e-7,
        6.2e-10,
    ]) * 1.0e-6
    np.testing.assert_allclose(coefficients, expected, rtol=2e-16, atol=0.0)


def test_argon_deck_separates_physical_energies_from_rate_fit_exponents():
    assert ARGON_IONIZATION_ENERGY_EV == 15.7596119
    assert ARGON_4S_METASTABLE_ENERGY_EV == 11.54835442
    assert ARGON_METASTABLE_IONIZATION_ENERGY_EV == (
        15.7596119 - 11.54835442)
    network = build_lee_lieberman_argon_volume_network()
    losses = {
        reaction.name: reaction.electron_energy_loss_eV
        for reaction in network.reactions
    }
    assert losses["e_Ar_ground_ionization"] != 18.68
    assert losses["e_Ar_metastable_excitation"] != 15.06
    assert losses["e_Ar_metastable_step_ionization"] != 4.95
    assert losses["e_Ar_metastable_superelastic_quench"] < 0.0


def test_argon_volume_deck_is_atom_and_charge_closed():
    network = build_lee_lieberman_argon_volume_network()
    for residual in network.reaction_conservation_residuals().values():
        assert residual["elements"] == {"Ar": 0.0}
        assert residual["charge_number"] == 0.0
    report = network.source_conservation_report(
        {
            "e": 2.0e16,
            "Ar": 1.0e20,
            "Ar*": 3.0e15,
            "Ar+": 2.0e16,
        },
        RateContext(electron_temperature_eV=3.0),
    )
    assert report["normalized_maximum_residual"] <= 1.0e-16


def test_pooling_event_emits_one_electron():
    network = build_lee_lieberman_argon_volume_network()
    pooling = next(
        reaction for reaction in network.reactions
        if reaction.name == "Ar_metastable_pooling_associative_ionization")
    assert pooling.products == {"Ar": 1.0, "Ar+": 1.0, "e": 1.0}
