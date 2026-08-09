import math

import pytest

from petch.reactor_global.chlorine_etch_products import (
    build_lee_graves_lieberman_etch_product_network,
)
from petch.reactor_global.network import RateContext


def _densities(network):
    return {name: 1.0e16 for name in network.species_names}


def test_table4_network_is_atom_and_charge_closed():
    network = build_lee_graves_lieberman_etch_product_network()

    assert len(network.species) == 13
    assert len(network.reactions) == 19
    assert network.elements == ("Cl", "Si")
    report = network.source_conservation_report(
        _densities(network),
        RateContext(electron_temperature_eV=3.0),
    )
    assert report["normalized_maximum_residual"] < 1.0e-14


def test_table4_printed_rates_replay_at_three_eV():
    network = build_lee_graves_lieberman_etch_product_network()
    context = RateContext(electron_temperature_eV=3.0)
    coefficients = {
        reaction.name: reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
    }

    assert coefficients["k1_sicl4_ionization"] == pytest.approx(
        7.03e-8 * 1.0e-6 * math.exp(-12.44 / 3.0)
    )
    assert coefficients[
        "k6_sicl3_dissociative_ionization_sicl2"
    ] == pytest.approx(
        4.90e-8 * 1.0e-6
        * math.exp(-13.9 / 3.0 + 6.89 / 9.0 - 1.45 / 27.0)
    )
    assert coefficients[
        "k7_sicl3_dissociative_ionization_sicl"
    ] == pytest.approx(
        2.41e-8 * 1.0e-6
        * math.exp(-14.75 / 3.0 + 2.504 / 9.0)
    )
    assert coefficients["k15_si+_mutual_neutralization"] == pytest.approx(
        5.0e-8 * 1.0e-6
    )


def test_table4_particle_deck_fails_closed_for_electron_power():
    network = build_lee_graves_lieberman_etch_product_network()

    assert not network.has_complete_electron_energy_ledger
    with pytest.raises(ValueError, match="electron-energy ledger is incomplete"):
        network.electron_power_loss_density_W_m3(
            _densities(network),
            RateContext(electron_temperature_eV=3.0),
        )
