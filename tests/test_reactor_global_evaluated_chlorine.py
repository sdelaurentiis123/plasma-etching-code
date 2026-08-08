import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from petch.reactor_global import (
    ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV,
    HAMILTON_2018_CL2_DISSOCIATION_STATES,
    HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256,
    MOLECULAR_CHLORINE_TOTAL_IONIZATION_THRESHOLD_EV,
    ElectronCollisionMomentKernel,
    ElectronEnergyDistribution,
    ElectronEnergyGrid,
    RateContext,
    ReactionNetwork,
    build_hamilton_dissociation_chlorine_particle_network,
    build_lee_lieberman_chlorine_particle_network,
    hamilton_2018_cl2_state_dissociation_collision_processes,
    hamilton_2018_cl2_state_dissociation_rates,
    hamilton_2018_cl2_state_dissociation_reactions,
    lee_lieberman_chlorine_species,
    nist_cl2_dissociative_attachment_cross_section_support,
    nist_hayes_atomic_chlorine_ionization_collision_process,
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
TABLE16 = TABLE25.with_name(
    "christophorou_olthoff_1999_table16_cl2_attachment.csv")
TABLE16_MANIFEST = TABLE16.with_name(
    "christophorou_olthoff_1999_table16_manifest.md")
HAMILTON_RATES = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "hamilton_2018_cl2_state_maxwellian_rates.csv"
)
HAMILTON_CROSS_SECTIONS = (
    ROOT / "research_sources" / "digitized"
    / "hamilton_2018_cl2_state_cross_sections.csv"
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


def test_nist_hayes_nonmaxwellian_row_enforces_physical_threshold():
    process = nist_hayes_atomic_chlorine_ionization_collision_process()
    assert process.kind == "IONIZATION"
    assert process.target == "Cl"
    assert process.product == "Cl+"
    assert process.energy_loss_eV == ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV
    assert process.electron_number_change == 1
    threshold_index = process.electron_energy_eV.index(
        ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV)
    assert all(
        value == 0.0
        for value in process.cross_section_m2[:threshold_index + 1]
    )
    assert process.electron_energy_eV[-1] == 200.0
    assert process.cross_section_m2[-1] == pytest.approx(2.63e-20)


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


def test_nist_cl2_attachment_table16_transcription_and_boundary():
    support = nist_cl2_dissociative_attachment_cross_section_support()
    assert len(support.electron_energy_eV) == 42
    assert support.electron_energy_eV[:4] == (0.05, 0.10, 0.20, 0.30)
    assert support.electron_energy_eV[-4:] == (11.0, 11.2, 11.6, 11.8)
    np.testing.assert_allclose(
        np.asarray(support.cross_section_m2[:8]) / 1.0e-20,
        [1.83, 1.04, 0.32, 0.081, 0.026, 0.013, 0.0088, 0.0065],
        rtol=0.0,
        atol=5.0e-16,
    )
    np.testing.assert_allclose(
        np.asarray(support.cross_section_m2[-8:]) / 1.0e-20,
        [
            0.0051, 0.0049, 0.0048, 0.0046,
            0.0045, 0.0042, 0.0041, 0.0043,
        ],
        rtol=0.0,
        atol=5.0e-18,
    )
    assert support.relative_uncertainty is None
    assert support.evidence_kind == "published_compilation"
    assert "Table 16" in support.source


def test_nist_cl2_attachment_executable_table_matches_pixel_audited_csv():
    support = nist_cl2_dissociative_attachment_cross_section_support()
    with TABLE16.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    np.testing.assert_allclose(
        support.electron_energy_eV,
        [float(row["electron_energy_eV"]) for row in rows],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(support.cross_section_m2) / 1.0e-20,
        [float(row["cross_section_1e_minus_20_m2"]) for row in rows],
        rtol=0.0,
        atol=5.0e-16,
    )
    manifest = TABLE16_MANIFEST.read_text(encoding="utf-8")
    assert "original-resolution visual review of all 42" in manifest
    assert (
        "97f73d5fcb067bd86a1415a2ff8c4aa097b51da279b32f4a4e2d19cbb3274164"
        in manifest
    )
    assert "does not set the" in manifest
    assert "cross section to zero outside the table" in manifest


@pytest.mark.parametrize("temperature", [0.3, 1.0, 3.0, 5.0])
def test_nist_cl2_attachment_support_reports_particle_and_energy_moments(
        temperature):
    support = nist_cl2_dissociative_attachment_cross_section_support()
    context = RateContext(temperature)
    rate = support.tabulated_rate_coefficient_si(context)
    energy = support.tabulated_incident_energy_moment_eV_m3_s(context)
    assert rate > 0.0
    assert energy > 0.0
    assert 0.05 < energy / rate < 11.8


def test_nist_cl2_attachment_support_exposes_lam_temperature_tail_gap():
    support = nist_cl2_dissociative_attachment_cross_section_support()
    rate_low, rate_high = support.rate_kernel_missing_fractions(3.0)
    energy_low, energy_high = (
        support.incident_energy_kernel_missing_fractions(3.0)
    )
    assert rate_low < 0.001
    assert rate_high > 0.09
    assert energy_low < 1.0e-5
    assert energy_high > 0.24


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
            # Exact table nodes round once through exp(log(k)); the measured
            # maximum over all 1,896 state/node pairs is 6.50e-15 relative.
            assert actual == pytest.approx(
                expected, rel=8.0e-15, abs=0.0)
            actual_sum += actual
        assert actual_sum == pytest.approx(
            float(row["summed_state_rate_m3_s"]),
            rel=4.0e-15,
            abs=0.0,
        )


def test_hamilton_nonmaxwellian_rows_are_hash_gated_exact_state_arrays():
    processes = hamilton_2018_cl2_state_dissociation_collision_processes()
    assert len(processes) == len(HAMILTON_2018_CL2_DISSOCIATION_STATES) == 8
    assert HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256 == (
        "7328d289542e23f2d12b4b172a19271120a3c5b62dc0dcd22a831569365dd288"
    )
    assert hashlib.sha256(HAMILTON_CROSS_SECTIONS.read_bytes()).hexdigest() == (
        HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256)
    for process, (state, threshold) in zip(
        processes, HAMILTON_2018_CL2_DISSOCIATION_STATES
    ):
        assert process.kind == "EXCITATION"
        assert process.target == "Cl2"
        assert process.product == f"2Cl via {state}"
        assert process.energy_loss_eV == threshold
        assert len(process.electron_energy_eV) == 10_001
        assert process.electron_energy_eV[:2] == (0.0, 0.02)
        assert process.electron_energy_eV[-1] == 200.0
        assert all(
            cross_section == 0.0
            for energy, cross_section in zip(
                process.electron_energy_eV, process.cross_section_m2)
            if energy < threshold
        )


def test_hamilton_nonmaxwellian_moments_reproduce_compact_maxwellian_total():
    temperature = 3.0
    processes = hamilton_2018_cl2_state_dissociation_collision_processes()
    thresholds = tuple(process.energy_loss_eV for process in processes)
    grid = ElectronEnergyGrid.piecewise_linear(
        (0.0, 0.5, 5.0, 20.0, 80.0, 200.0),
        (200, 900, 900, 600, 400),
        inserted_boundaries_eV=thresholds,
    )
    distribution = ElectronEnergyDistribution.maxwellian(grid, temperature)
    actual = sum(
        ElectronCollisionMomentKernel.from_process(grid, process).evaluate(
            distribution,
            maximum_unresolved_population_fraction=1.0e-7,
        ).rate_coefficient_m3_s
        for process in processes
    )
    expected = sum(
        provider.coefficient_si(RateContext(temperature))
        for _, _, provider in hamilton_2018_cl2_state_dissociation_rates()
    )
    # The compact table integrates the exact same 50,000-point arrays against
    # an analytic Maxwellian kernel. The residual here is solely the declared
    # piecewise-constant EEPF discretization used by the production solver.
    assert actual == pytest.approx(expected, rel=2.0e-4)


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


def test_hamilton_particle_deck_replaces_only_legacy_neutral_dissociation():
    legacy = build_lee_lieberman_chlorine_particle_network()
    upgraded = build_hamilton_dissociation_chlorine_particle_network()
    legacy_names = {reaction.name for reaction in legacy.reactions}
    upgraded_names = {reaction.name for reaction in upgraded.reactions}
    hamilton_names = {
        f"e_Cl2_dissociation_{state}"
        for state, _ in HAMILTON_2018_CL2_DISSOCIATION_STATES
    }

    assert len(upgraded.reactions) == len(legacy.reactions) - 1 + 8
    assert "e_Cl2_dissociation" in legacy_names
    assert "e_Cl2_dissociation" not in upgraded_names
    assert upgraded_names == (
        legacy_names - {"e_Cl2_dissociation"}
    ) | hamilton_names
    upgraded.assert_closed_conservation()
    assert not upgraded.has_complete_electron_energy_ledger


@pytest.mark.parametrize("temperature", [0.3, 1.0, 3.0, 5.0])
def test_hamilton_particle_deck_sums_exact_state_resolved_rate(temperature):
    context = RateContext(temperature)
    upgraded = build_hamilton_dissociation_chlorine_particle_network()
    expected = sum(
        provider.coefficient_si(context)
        for _, _, provider in hamilton_2018_cl2_state_dissociation_rates()
    )
    actual = sum(
        reaction.rate_coefficient.coefficient_si(context)
        for reaction in upgraded.reactions
        if reaction.name.startswith("e_Cl2_dissociation_")
    )
    assert actual == pytest.approx(expected, rel=4.0e-15, abs=0.0)
