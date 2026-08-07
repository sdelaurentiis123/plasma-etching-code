import csv
from pathlib import Path

import numpy as np
import pytest

from petch.reactor_global import (
    ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV,
    HAMILTON_2018_CL2_DISSOCIATION_STATES,
    MOLECULAR_CHLORINE_TOTAL_IONIZATION_THRESHOLD_EV,
    RateContext,
    ReactionNetwork,
    hamilton_2018_cl2_state_dissociation_rates,
    hamilton_2018_cl2_state_dissociation_reactions,
    lee_lieberman_chlorine_species,
    nist_hayes_atomic_chlorine_ionization_rate,
    nist_molecular_chlorine_total_ionization_rate,
)

ROOT = Path(__file__).resolve().parents[1]
TABLE25 = (
    ROOT / "research_sources" / "digitized"
    / "christophorou_olthoff_1999_table25_atomic_cl_ionization.csv"
)
TABLE25_MANIFEST = TABLE25.with_name(
    "christophorou_olthoff_1999_table25_manifest.md")
TABLE12 = TABLE25.with_name(
    "christophorou_olthoff_1999_table12_cl2_total_ionization.csv")
TABLE12_MANIFEST = TABLE12.with_name(
    "christophorou_olthoff_1999_table12_manifest.md")
HAMILTON_RATES = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "hamilton_2018_cl2_state_maxwellian_rates.csv"
)


def test_nist_hayes_table25_transcription_and_evidence():
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    assert len(coefficient.electron_energy_eV) == 48
    assert coefficient.electron_energy_eV[:4] == (11.0, 12.0, 13.0, 14.0)
    assert coefficient.electron_energy_eV[-4:] == (
        170.0, 180.0, 190.0, 200.0)
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2[:4]) / 1.0e-20,
        [0.00, 0.01, 0.02, 0.24],
        rtol=0.0,
        atol=2.0e-16,
    )
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2[-4:]) / 1.0e-20,
        [2.81, 2.72, 2.68, 2.63],
        rtol=0.0,
        atol=5.0e-16,
    )
    assert coefficient.threshold_eV == (
        ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV)
    assert coefficient.relative_uncertainty == 0.14
    assert coefficient.evidence_kind == "measured"


def test_nist_hayes_executable_table_matches_pixel_audited_csv():
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    with TABLE25.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    np.testing.assert_allclose(
        coefficient.electron_energy_eV,
        [float(row["electron_energy_eV"]) for row in rows],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2) / 1.0e-20,
        [float(row["cross_section_1e_minus_20_m2"]) for row in rows],
        rtol=0.0,
        atol=5.0e-16,
    )
    manifest = TABLE25_MANIFEST.read_text(encoding="utf-8")
    assert "original-resolution visual review of all 48" in manifest
    assert (
        "6a01d03172e2d49619998e0593d14e9b547ad01803f45a6677068768cf599c25"
        in manifest
    )


@pytest.mark.parametrize("temperature", [2.0, 3.0, 5.0, 8.0, 10.0])
def test_nist_hayes_rate_is_positive_on_reactor_temperature_domain(
        temperature):
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    assert coefficient.coefficient_si(RateContext(temperature)) > 0.0
    assert coefficient.maxwellian_kernel_tail_fraction(temperature) <= 1.0e-6


def test_nist_hayes_rate_rejects_temperature_with_material_unknown_tail():
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    with pytest.raises(ValueError, match="unmeasured cross-section tail"):
        coefficient.coefficient_si(RateContext(15.0))


def test_nist_molecular_chlorine_table12_transcription_and_boundary():
    coefficient = nist_molecular_chlorine_total_ionization_rate()
    assert len(coefficient.electron_energy_eV) == 29
    assert coefficient.electron_energy_eV[:4] == (11.5, 12.0, 13.0, 14.0)
    assert coefficient.electron_energy_eV[-4:] == (
        85.0, 90.0, 95.0, 100.0)
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2[:4]) / 1.0e-20,
        [0.03, 0.11, 0.25, 0.43],
        rtol=0.0,
        atol=2.0e-16,
    )
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2[-4:]) / 1.0e-20,
        [6.28, 6.25, 6.22, 6.19],
        rtol=0.0,
        atol=2.0e-15,
    )
    assert coefficient.threshold_eV == (
        MOLECULAR_CHLORINE_TOTAL_IONIZATION_THRESHOLD_EV)
    assert coefficient.relative_uncertainty is None
    assert coefficient.evidence_kind == "published_compilation"
    assert "total" in coefficient.source


def test_nist_molecular_chlorine_executable_table_matches_pixel_audited_csv():
    coefficient = nist_molecular_chlorine_total_ionization_rate()
    with TABLE12.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    np.testing.assert_allclose(
        coefficient.electron_energy_eV,
        [float(row["electron_energy_eV"]) for row in rows],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2) / 1.0e-20,
        [float(row["cross_section_1e_minus_20_m2"]) for row in rows],
        rtol=0.0,
        atol=2.0e-15,
    )
    manifest = TABLE12_MANIFEST.read_text(encoding="utf-8")
    assert "original-resolution visual review of all 29" in manifest
    assert (
        "46a3b5a8b9aa41ec69279cbf45112161139be9ec49a7bd0e964b048a4099d5e3"
        in manifest
    )
    assert "relative production" in manifest
    assert "`Cl2+` and `Cl+` is unknown" in manifest


@pytest.mark.parametrize("temperature", [0.3, 1.0, 2.0, 3.0, 5.0])
def test_nist_molecular_chlorine_rate_covers_industrial_temperature_band(
        temperature):
    coefficient = nist_molecular_chlorine_total_ionization_rate()
    assert coefficient.coefficient_si(RateContext(temperature)) > 0.0
    assert coefficient.maxwellian_kernel_tail_fraction(temperature) <= 1.0e-6


def test_nist_molecular_chlorine_rate_rejects_unsupported_hot_eedf():
    coefficient = nist_molecular_chlorine_total_ionization_rate()
    with pytest.raises(ValueError, match="unmeasured cross-section tail"):
        coefficient.coefficient_si(RateContext(8.0))


def test_hamilton_state_rates_match_compact_table_at_every_node():
    providers = {
        state: provider
        for state, _, provider
        in hamilton_2018_cl2_state_dissociation_rates()
    }
    with HAMILTON_RATES.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 237
    assert float(rows[0]["electron_temperature_eV"]) == 0.3
    assert float(rows[-1]["electron_temperature_eV"]) == 5.0
    for row in rows:
        context = RateContext(float(row["electron_temperature_eV"]))
        actual_sum = 0.0
        for state, _ in HAMILTON_2018_CL2_DISSOCIATION_STATES:
            actual = providers[state].coefficient_si(context)
            expected = float(row[f"{state}_m3_s"])
            assert actual == pytest.approx(expected, rel=4.0e-15)
            actual_sum += actual
        assert actual_sum == pytest.approx(
            float(row["summed_state_rate_m3_s"]), rel=4.0e-15)


@pytest.mark.parametrize("temperature", [0.299, 5.001])
def test_hamilton_state_rates_fail_outside_declared_domain(temperature):
    for _, _, provider in hamilton_2018_cl2_state_dissociation_rates():
        with pytest.raises(
                ValueError, match="outside the tabulated rate domain"):
            provider.coefficient_si(RateContext(temperature))


def test_hamilton_state_reactions_close_atoms_charge_and_energy():
    reactions = hamilton_2018_cl2_state_dissociation_reactions()
    assert len(reactions) == len(HAMILTON_2018_CL2_DISSOCIATION_STATES) == 8
    network = ReactionNetwork(
        species=lee_lieberman_chlorine_species(),
        reactions=reactions,
    )
    network.assert_closed_conservation()
    assert network.has_complete_electron_energy_ledger
    assert {
        reaction.name: reaction.electron_energy_loss_eV
        for reaction in reactions
    } == {
        f"e_Cl2_dissociation_{state}": excitation_eV
        for state, excitation_eV
        in HAMILTON_2018_CL2_DISSOCIATION_STATES
    }
    densities = {
        "e": 2.0e16,
        "Cl2": 8.0e19,
        "Cl": 1.0e19,
        "Cl2+": 1.0e15,
        "Cl+": 1.0e15,
        "Cl-": 1.0e15,
    }
    assert network.electron_power_loss_density_W_m3(
        densities, RateContext(3.0)) > 0.0
