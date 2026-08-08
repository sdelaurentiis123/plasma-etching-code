import numpy as np
import pytest

from petch.reactor_global import (
    RateContext,
    build_kemaneci_2014_comsol_nonelastic_chlorine_network,
    build_kemaneci_2014_forward_chlorine_network,
    kemaneci_2014_chlorine_species,
)


def _coefficient(network, label, context):
    reaction = next(
        item for item in network.reactions
        if item.name.startswith(f"k{label:02d}_"))
    return reaction.rate_coefficient.coefficient_si(context)


def test_kemaneci_forward_species_and_reactions_are_closed_and_explicit():
    species = kemaneci_2014_chlorine_species()
    assert len(species) == 11
    assert {item.name for item in species} == {
        "e", "Cl2", "Cl2(v=1)", "Cl2(v=2)", "Cl2(v=3)",
        "Cl", "Cl(2P1/2)", "Cl(1P5/2)", "Cl2+", "Cl+", "Cl-",
    }
    network = build_kemaneci_2014_forward_chlorine_network()
    assert len(network.reactions) == 36
    labels = {int(item.name[1:3]) for item in network.reactions}
    assert labels == set(range(1, 38)) - {32}
    network.assert_closed_conservation()
    assert not network.has_complete_electron_energy_ledger


def test_kemaneci_table4_electron_fits_replay_every_forward_row():
    network = build_kemaneci_2014_forward_chlorine_network()
    temperature = 2.3
    context = RateContext(temperature, gas_temperature_K=600.0)
    t = temperature

    def attachment(a, b):
        return (
            a * t ** -1.18 * np.exp(-3.98 / t)
            + b * t ** -1.33 * np.exp(-0.11 / (t + 0.014))
        )

    def log_term(a, shift, width):
        return a * np.exp(
            -(np.log(t) + shift) ** 2 / (2.0 * width ** 2))
    expected = {
        1: 1.04e-13 * t ** -0.29 * np.exp(-8.84 / t),
        2: 5.12e-14 * t ** 0.48 * np.exp(-12.34 / t),
        3: 2.14e-13 * t ** -0.07 * np.exp(-25.26 / t),
        4: 2.27e-16 * t ** 1.92 * np.exp(-21.26 / t),
        5: attachment(3.43e-15, 3.05e-16),
        6: attachment(14.06e-15, 12.51e-16),
        7: attachment(30.18e-15, 26.84e-16),
        8: attachment(46.31e-15, 41.18e-16),
        9: 2.94e-16 * t ** 0.19 * np.exp(-18.79 / t),
        10: 3.99e-12 * t ** -1.5 * np.exp(-7.51 / t - 0.0001 / t ** 2),
        11: (3.28e-17 * t ** -1.12 * np.exp(-0.37 / t)
             + log_term(2.86e-17, 0.99, 1.06)),
        12: (1.30e-17 * t ** -1.24 * np.exp(-0.41 / t)
             + log_term(6.08e-18, 0.94, 1.02)),
        13: (3.00e-16 * t ** -1.00 * np.exp(-0.37 / t)
             + log_term(4.61e-16, 1.04, 1.10)),
        14: (3.00e-16 * t ** -1.00 * np.exp(-0.37 / t)
             + log_term(4.61e-16, 1.04, 1.10)),
        15: (1.25e-16 * t ** -1.13 * np.exp(-0.36 / t)
             + log_term(1.06e-16, 1.01, 1.06)),
        16: 9.00e-14 * t ** -0.50,
        # Independent 50-digit evaluation of the visually audited Table-4
        # row.  Keeping a literal here prevents the source exponent sign from
        # being copied from the implementation into its own regression gate.
        17: 1.2941925109736606e-14,
        18: 7.03e-17 * t ** 0.55 * np.exp(
            -2.15 / t - 1.5 / t ** 2 - 2.05 / t ** 3),
        19: 3.17e-14 * t ** 0.53 * np.exp(-13.29 / t),
        20: 3.17e-14 * t ** 0.53 * np.exp(-13.19 / t),
        21: 4.33e-14 * t ** 0.55 * np.exp(-0.15 / t - 0.85 / t ** 2),
        23: 9.02e-15 * t ** 0.92 * np.exp(-4.88 / t),
        24: 3.62e-15 * t ** 0.72 * np.exp(-25.38 / t),
    }
    for label, value in expected.items():
        assert _coefficient(network, label, context) == pytest.approx(
            value, rel=5.0e-15, abs=0.0)
    assert _coefficient(network, 22, context) == 1.0e5


def test_kemaneci_table4_heavy_rates_replay_at_declared_gas_temperature():
    network = build_kemaneci_2014_forward_chlorine_network()
    temperature = 600.0
    context = RateContext(2.0, gas_temperature_K=temperature)
    assert _coefficient(network, 25, context) == pytest.approx(
        5.0e-14 * (300.0 / temperature) ** 0.5,
        rel=3.0e-15, abs=0.0)
    assert _coefficient(network, 26, context) == 5.0e-14
    assert _coefficient(network, 27, context) == pytest.approx(
        5.0e-14 * (300.0 / temperature) ** 0.5,
        rel=3.0e-15, abs=0.0)
    for label in range(28, 32):
        assert _coefficient(network, label, context) == 5.40e-16
    assert _coefficient(network, 33, context) == pytest.approx(
        3.50e-45 * np.exp(810.0 / temperature),
        rel=3.0e-15, abs=0.0)
    assert _coefficient(network, 34, context) == pytest.approx(
        8.75e-46 * np.exp(810.0 / temperature),
        rel=3.0e-15, abs=0.0)
    for label in range(35, 38):
        assert _coefficient(network, label, context) == pytest.approx(
            1.30e-17 * (temperature / 300.0) ** 0.5,
            rel=3.0e-15, abs=0.0)


@pytest.mark.parametrize("temperature", [0.499, 10.001])
def test_kemaneci_electron_fits_refuse_temperature_extrapolation(temperature):
    network = build_kemaneci_2014_forward_chlorine_network()
    context = RateContext(temperature, gas_temperature_K=600.0)
    with pytest.raises(ValueError, match="published fit domain"):
        _coefficient(network, 1, context)


def test_kemaneci_forward_replay_cannot_masquerade_as_power_model():
    network = build_kemaneci_2014_forward_chlorine_network()
    densities = {name: 1.0e15 for name in network.species_names}
    with pytest.raises(ValueError, match="electron-energy ledger is incomplete"):
        network.electron_power_loss_density_W_m3(
            densities, RateContext(2.0, gas_temperature_K=600.0))


def test_kemaneci_comsol_nonelastic_replay_is_explicitly_44_rows():
    network = build_kemaneci_2014_comsol_nonelastic_chlorine_network()
    assert len(network.reactions) == 44
    assert sum("_comsol" in item.name for item in network.reactions) == 8
    network.assert_closed_conservation()
    assert not network.has_complete_electron_energy_ledger

    context = RateContext(2.3, gas_temperature_K=600.0)
    # COMSOL uses 13.29 in both ground- and fine-state atomic ionization fits;
    # the primary paper prints 13.19 for row 20.
    assert _coefficient(network, 20, context) == _coefficient(
        network, 19, context)


def test_kemaneci_comsol_reverse_rows_replay_raw_model_expressions():
    network = build_kemaneci_2014_comsol_nonelastic_chlorine_network()
    context = RateContext(2.3, gas_temperature_K=600.0)
    reverse = {
        int(item.name[1:3]): item for item in network.reactions
        if "_comsol" in item.name
    }
    assert set(reverse) == {10, 11, 12, 13, 14, 15, 17, 18}

    # Independent high-precision evaluations of the raw COMSOL kf*exp(dE/Te)
    # expressions.  Literals prevent a copied sign or gap from self-passing.
    expected = {
        10: 4.5032358682758177e-14,
        11: 1.8604247950516934e-17,
        12: 5.713331032574398e-18,
        13: 2.260188387782102e-16,
        14: 2.260188387782102e-16,
        15: 6.917410928737043e-17,
        17: 2.3276134651034147e-14,
        18: 2.311809726948781e-15,
    }
    for label, value in expected.items():
        assert reverse[label].rate_coefficient.coefficient_si(
            context) == pytest.approx(value, rel=5.0e-15, abs=0.0)

    with pytest.raises(ValueError, match="published fit domain"):
        reverse[10].rate_coefficient.coefficient_si(
            RateContext(0.499, gas_temperature_K=600.0))

    assert reverse[17].electron_energy_loss_eV == -1.35
    assert reverse[18].electron_energy_loss_eV == -10.17
