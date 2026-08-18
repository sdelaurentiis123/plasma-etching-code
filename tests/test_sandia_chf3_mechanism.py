import numpy as np

from petch.reactor_global.network import RateContext
from petch.reactor_global.sandia_chf3_mechanism import (
    KELVIN_PER_EV,
    SANDIA_TABLE9_ROWS,
    build_sandia_2001_chf3_table9_network,
    sandia_table9_evidence_counts,
)


def test_table9_is_complete_conserved_and_keeps_evidence_classes_visible():
    network = build_sandia_2001_chf3_table9_network()
    assert len(SANDIA_TABLE9_ROWS) == len(network.reactions) == 38
    assert network.reaction_conservation_residuals()
    assert sandia_table9_evidence_counts() == {
        "estimated": 14,
        "regressed": 24,
    }
    assert network.has_complete_electron_energy_ledger


def test_kelvin_to_ev_rate_conversion_reproduces_source_expression():
    network = build_sandia_2001_chf3_table9_network()
    temperature_eV = 2.2
    context = RateContext(
        electron_temperature_eV=temperature_eV,
        gas_temperature_K=350.0,
    )
    for row, reaction in zip(SANDIA_TABLE9_ROWS[:33], network.reactions[:33]):
        if "e" not in row.reactants:
            continue
        source_temperature_K = temperature_eV * KELVIN_PER_EV
        expected_cm3_s = (
            row.A_cm3_s_K_minus_B
            * source_temperature_K ** row.B
            * np.exp(-row.C_K / source_temperature_K)
        )
        assert np.isclose(
            reaction.rate_coefficient.coefficient_si(context),
            expected_cm3_s * 1.0e-6,
            rtol=2.0e-14,
            atol=0.0,
        )


def test_estimated_ion_ion_temperature_law_reproduces_table9():
    network = build_sandia_2001_chf3_table9_network()
    context = RateContext(
        electron_temperature_eV=2.2,
        gas_temperature_K=350.0,
    )
    row = SANDIA_TABLE9_ROWS[32]
    reaction = network.reactions[32]
    expected = row.A_cm3_s_K_minus_B * 350.0 ** row.B * 1.0e-6
    assert np.isclose(
        reaction.rate_coefficient.coefficient_si(context), expected,
        rtol=2.0e-14, atol=0.0,
    )
