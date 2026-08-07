import numpy as np
import pytest

from petch.reactor_global import (
    CM3_TO_M3,
    RateContext,
    build_lee_lieberman_chlorine_particle_network,
)


def test_chlorine_particle_deck_replays_table2_rates_in_si():
    network = build_lee_lieberman_chlorine_particle_network()
    temperature = 3.5
    context = RateContext(electron_temperature_eV=temperature)
    coefficients = {
        reaction.name: reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
    }
    assert coefficients["e_Cl2_nondissociative_ionization"] == pytest.approx(
        9.21e-8 * np.exp(-12.9 / temperature) * CM3_TO_M3,
        rel=2.0e-16,
    )
    assert coefficients["e_Cl2_dissociative_ionization"] == pytest.approx(
        3.88e-9 * np.exp(-15.5 / temperature) * CM3_TO_M3,
        rel=2.0e-16,
    )
    assert coefficients["e_Cl2_ion_pair_production"] == pytest.approx(
        8.55e-10 * np.exp(-12.65 / temperature) * CM3_TO_M3,
        rel=2.0e-16,
    )
    assert coefficients["e_Cl2_dissociation"] == pytest.approx(
        3.80e-8 * np.exp(-3.824 / temperature) * CM3_TO_M3,
        rel=2.0e-16,
    )
    ratio = temperature / 12.96
    atomic_ionization_coefficients = (
        1.419e-7,
        -1.864e-8,
        -5.439e-8,
        3.306e-8,
        -3.54e-9,
        -2.915e-8,
    )
    assert coefficients["e_Cl_ionization"] == pytest.approx(
        ratio ** 0.5
        * np.exp(-12.96 / temperature)
        * sum(
            coefficient * np.log10(ratio) ** order
            for order, coefficient in enumerate(
                atomic_ionization_coefficients)
        )
        * CM3_TO_M3,
        rel=2.0e-16,
    )
    assert coefficients["Clminus_Cl2plus_neutralization"] == pytest.approx(
        5.0e-14,
        rel=2.0e-16,
    )
    assert coefficients["Clminus_Clplus_neutralization"] == pytest.approx(
        5.0e-14,
        rel=2.0e-16,
    )
    assert coefficients["e_Clminus_detachment"] == pytest.approx(
        2.63e-8 * np.exp(-5.37 / temperature) * CM3_TO_M3,
        rel=2.0e-16,
    )


def test_chlorine_particle_deck_is_atom_and_charge_closed():
    network = build_lee_lieberman_chlorine_particle_network()
    assert len(network.reactions) == 9
    for residual in network.reaction_conservation_residuals().values():
        assert residual["elements"] == {"Cl": 0.0}
        assert residual["charge_number"] == 0.0
    densities = {
        "e": 2.0e16,
        "Cl2": 8.0e19,
        "Cl": 1.0e20,
        "Cl2+": 4.0e15,
        "Cl+": 1.8e16,
        "Cl-": 3.0e15,
    }
    report = network.source_conservation_report(
        densities,
        RateContext(electron_temperature_eV=3.0),
    )
    assert report["normalized_maximum_residual"] <= 2.0e-16


def test_chlorine_network_is_particle_ready_but_power_fail_closed():
    network = build_lee_lieberman_chlorine_particle_network()
    assert not network.has_complete_electron_energy_ledger
    with pytest.raises(ValueError, match="electron-energy ledger is incomplete"):
        network.electron_power_loss_density_W_m3(
            {
                "e": 2.0e16,
                "Cl2": 8.0e19,
                "Cl": 1.0e20,
                "Cl2+": 4.0e15,
                "Cl+": 1.8e16,
                "Cl-": 3.0e15,
            },
            RateContext(electron_temperature_eV=3.0),
        )


def test_attachment_and_ion_pair_rows_do_not_hide_charge_creation():
    network = build_lee_lieberman_chlorine_particle_network()
    attachment = next(
        reaction for reaction in network.reactions
        if reaction.name == "e_Cl2_dissociative_attachment")
    ion_pair = next(
        reaction for reaction in network.reactions
        if reaction.name == "e_Cl2_ion_pair_production")
    assert attachment.reactants == {"e": 1.0, "Cl2": 1.0}
    assert attachment.products == {"Cl": 1.0, "Cl-": 1.0}
    assert ion_pair.products == {
        "Cl+": 1.0,
        "Cl-": 1.0,
        "e": 1.0,
    }
