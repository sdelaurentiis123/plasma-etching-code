from pathlib import Path

import numpy as np
import pytest

from petch.experimental_data import (
    KARAHASHI_2007_FIGURE4_SHA256,
    KARAHASHI_2007_FIGURE10_SHA256,
    load_karahashi_2007_cf3_product_fractions,
    load_karahashi_2007_reactive_ion_yields,
)
from petch.ion_energy_deposition import FLUOROCARBON_FILM
from petch.reactive_ion_beam import Karahashi2007ReactiveIonYieldTable
from petch.reactive_ion_event import (
    Karahashi2007CF3ProductBranchTable,
    Karahashi2007ReactiveIonEventKernel,
    mass_partitioned_projectile_fragments,
    transmit_fragmented_projectile_through_layer,
)


DATA = (
    Path(__file__).parents[1] / "data" / "experimental" / "karahashi_2007")


@pytest.fixture(scope="module")
def kernel():
    yield_table = Karahashi2007ReactiveIonYieldTable.from_observations(
        load_karahashi_2007_reactive_ion_yields(
            DATA / "figure4_reactive_ion_yields.csv"),
        source_table_sha256=KARAHASHI_2007_FIGURE4_SHA256,
    )
    product_table = Karahashi2007CF3ProductBranchTable(
        load_karahashi_2007_cf3_product_fractions(
            DATA / "figure10_cf3_product_fractions.csv"),
        source_table_sha256=KARAHASHI_2007_FIGURE10_SHA256,
    )
    return Karahashi2007ReactiveIonEventKernel(
        yield_table, product_table)


def test_cf3_fragment_energy_partitions_by_mass_and_closes_exactly():
    energy = np.array([500.0, 1000.0, 2000.0])
    fragments = mass_partitioned_projectile_fragments("CF3+", energy)
    assert [item.element for item in fragments] == ["C", "F"]
    carbon, fluorine = fragments
    assert carbon.multiplicity == 1
    assert fluorine.multiplicity == 3
    reconstructed = (
        carbon.energy_per_fragment_eV
        + 3.0 * fluorine.energy_per_fragment_eV)
    assert np.allclose(reconstructed, energy, rtol=0.0, atol=5e-13)
    assert np.all(
        fluorine.energy_per_fragment_eV
        > carbon.energy_per_fragment_eV)
    assert fluorine.energy_per_fragment_eV[1] == pytest.approx(
        1000.0 * 18.998 / (12.011 + 3.0 * 18.998))


def test_fragmented_csda_transmission_has_identity_and_finite_range():
    energy = np.array([500.0, 1000.0])
    identity = transmit_fragmented_projectile_through_layer(
        "CF3+", energy, 1.0, 0.0, FLUOROCARBON_FILM)
    assert np.allclose(identity.total_residual_energy_eV, energy)
    assert np.array_equal(
        identity.total_deposited_energy_eV, np.zeros(2))

    stopped = transmit_fragmented_projectile_through_layer(
        "CF3+", energy, 0.2, 100.0, FLUOROCARBON_FILM)
    assert np.array_equal(
        stopped.total_residual_energy_eV, np.zeros(2))
    assert np.allclose(stopped.total_deposited_energy_eV, energy)


def test_conditional_projection_uses_no_target_fit_and_stays_nonproduction(
        kernel):
    expected_yield = {500.0: 0.8638, 1000.0: 1.4703, 2000.0: 1.7549}
    for energy, expected in expected_yield.items():
        outcome = kernel.evaluate_conditional_cf3(
            energy, acknowledge_unresolved_incidence_angle=True)
        assert outcome.removed_sio2_formula_per_ion == pytest.approx(expected)
        assert sum(
            outcome.conditional_sifx_particles_per_ion.values()
        ) == pytest.approx(expected, rel=3e-16)
        assert sum(
            outcome.normalized_product_fraction.values()
        ) == pytest.approx(1.0, rel=3e-16)
        assert outcome.incident_atoms_per_ion == {"C": 1.0, "F": 3.0}
        assert not outcome.production_eligible
        assert len(outcome.unresolved_inventories) == 7


def test_event_fluorine_ledger_exposes_the_1000ev_tight_balance(kernel):
    outcome = kernel.evaluate_conditional_cf3(
        1000.0, acknowledge_unresolved_incidence_angle=True)
    assert outcome.required_f_atoms_per_ion > 3.0
    assert outcome.required_f_atoms_lower < 3.0
    assert outcome.required_f_atoms_upper > 3.0
    assert outcome.unresolved_f_balance_lower < 0.0
    assert outcome.unresolved_f_balance_upper > 0.0
    assert outcome.normalized_product_fraction["SiF2"] > (
        outcome.normalized_product_fraction["SiF"])
    assert outcome.normalized_product_fraction["SiF2"] > (
        outcome.normalized_product_fraction["SiF4"])


def test_product_kernel_refuses_interpolation_species_and_angle(kernel):
    with pytest.raises(ValueError, match="incidence angle is unreported"):
        kernel.evaluate_conditional_cf3(1000.0)
    with pytest.raises(ValueError, match="exact"):
        kernel.evaluate_conditional_cf3(
            750.0, acknowledge_unresolved_incidence_angle=True)
    with pytest.raises(ValueError, match="only for CF3"):
        kernel.product_table.evaluate_exact("CF2+", 1000.0)
    with pytest.raises(ValueError, match="resolved atomic formula"):
        mass_partitioned_projectile_fragments("ions", 1000.0)
