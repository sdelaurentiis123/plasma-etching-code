import numpy as np
import pytest
from pathlib import Path

from petch.boundary_state import SpeciesBoundaryState
from petch.reactor_boundary import (
    PlasmaDiagnosticState,
    ReactorSpeciesFlux,
    TabulatedReactorFluxDeck,
    build_diagnostic_virtual_sheath_boundary,
    build_tabulated_reactor_boundary,
    load_krueger_2024_reactor_flux_deck,
)
from petch.sheath import (
    CollisionlessRFSheath,
    CollisionlessWaveformSheath,
    PeriodicSheathVoltage,
)

KRUEGER_DATA = (
    Path(__file__).parents[1] / "data" / "experimental" / "krueger_2024")


def _diagnostic(**overrides):
    values = dict(
        electron_density_m3=2.0e15,
        electron_temperature_eV=3.0,
        ion_name="Ar+",
        ion_mass_amu=39.948,
        source="manufactured diagnostic gate",
        density_evidence_kind="measured",
        temperature_evidence_kind="measured",
        electropositive_bohm_flux_closure=True,
    )
    values.update(overrides)
    return PlasmaDiagnosticState(**values)


def _waveform(*, evidence_kind="assumed", dc_v=80.0, amplitude_v=20.0):
    return PeriodicSheathVoltage.sinusoidal(
        dc_v=dc_v,
        amplitude_v=amplitude_v,
        frequency_hz=4.0e5,
        source="manufactured waveform gate",
        evidence_kind=evidence_kind,
    )


def test_waveform_sheath_preserves_legacy_sinusoidal_operator():
    phase = 2.0 * np.pi * (np.arange(32) + 0.5) / 32.0
    legacy = CollisionlessRFSheath(
        V_dc=80.0, V_rf=20.0, frequency_hz=4.0e5,
        Te_eV=3.0, ion_mass_amu=39.948, thickness_m=8.0e-4)
    waveform = CollisionlessWaveformSheath(
        waveform=_waveform(), Te_eV=3.0, ion_mass_amu=39.948,
        thickness_m=8.0e-4)

    assert np.allclose(
        waveform.ion_impact_energies(phase),
        legacy.ion_impact_energies(phase),
        rtol=0.0,
        atol=1e-12,
    )


def test_multiharmonic_waveform_changes_iedf_without_changing_mean_drop():
    phase = 2.0 * np.pi * (np.arange(64) + 0.5) / 64.0
    single = CollisionlessWaveformSheath(
        waveform=_waveform(), Te_eV=3.0, ion_mass_amu=39.948,
        thickness_m=8.0e-4)
    shaped_voltage = PeriodicSheathVoltage(
        fundamental_frequency_hz=4.0e5,
        dc_v=80.0,
        harmonic_number=np.array([1, 2]),
        sine_v=np.array([20.0, 8.0]),
        cosine_v=np.array([0.0, -3.0]),
        source="manufactured dual-harmonic waveform",
        evidence_kind="assumed",
    )
    shaped = CollisionlessWaveformSheath(
        waveform=shaped_voltage, Te_eV=3.0, ion_mass_amu=39.948,
        thickness_m=8.0e-4)

    single_energy = single.ion_impact_energies(phase, steps_per_period=256)
    shaped_energy = shaped.ion_impact_energies(phase, steps_per_period=256)
    assert not np.allclose(shaped_energy, single_energy, rtol=1e-4, atol=1e-4)
    assert abs(shaped_energy.mean() - single_energy.mean()) < 2.0


def test_development_boundary_is_current_closed_and_carries_continuous_densities():
    state = build_diagnostic_virtual_sheath_boundary(
        _diagnostic(),
        _waveform(),
        reference_plane_m=2.0e-6,
        collisionless_justification="manufactured low-pressure gate",
        n_phase=32,
        normal_energy_bins=8,
        density_phase_count=512,
    )

    ion = state.get("Ar+")
    electron = state.get("electron")
    assert ion.density_model is not None
    assert electron.density_model is not None
    assert np.isclose(state.current_density_A_m2, 0.0, atol=1e-14)
    assert np.isclose(ion.flux_m2_s, electron.flux_m2_s)
    assert state.provenance["supports_prediction"] is False
    assert ion.provenance["ion_flux_closure"] == "electropositive_bohm_flux"
    assert state.provenance["volume_boltzmann_electron_term"] is False


def test_predictive_mode_refuses_assumed_waveform_and_missing_bohm_authorization():
    with pytest.raises(ValueError, match="predictive mode requires"):
        build_diagnostic_virtual_sheath_boundary(
            _diagnostic(),
            _waveform(evidence_kind="assumed"),
            reference_plane_m=2.0e-6,
            collisionless_justification="manufactured low-pressure gate",
            claim_mode="predictive",
            model_validation_reference="manufactured validation reference",
            n_phase=16,
            normal_energy_bins=8,
            density_phase_count=256,
        )

    with pytest.raises(ValueError, match="explicitly authorize"):
        build_diagnostic_virtual_sheath_boundary(
            _diagnostic(electropositive_bohm_flux_closure=False),
            _waveform(),
            reference_plane_m=2.0e-6,
            collisionless_justification="manufactured low-pressure gate",
            n_phase=16,
            normal_energy_bins=8,
            density_phase_count=256,
        )


def test_predictive_mode_accepts_evidenced_nonnegative_full_waveform():
    state = build_diagnostic_virtual_sheath_boundary(
        _diagnostic(),
        _waveform(evidence_kind="measured_sheath_voltage"),
        reference_plane_m=2.0e-6,
        collisionless_justification="mean free path exceeds the modeled sheath",
        claim_mode="predictive",
        model_validation_reference="NIST finite-transit sheath validation gate",
        n_phase=24,
        normal_energy_bins=8,
        density_phase_count=512,
    )

    assert state.provenance["supports_prediction"] is True
    assert state.get("Ar+").provenance["waveform_evidence_kind"] == (
        "measured_sheath_voltage")
    assert state.get("Ar+").mean_energy_eV > 0.0


def test_predictive_mode_refuses_a_sign_reversing_sheath_drop():
    with pytest.raises(ValueError, match="cannot reverse sign"):
        build_diagnostic_virtual_sheath_boundary(
            _diagnostic(),
            _waveform(
                evidence_kind="measured_sheath_voltage", dc_v=10.0, amplitude_v=20.0),
            reference_plane_m=2.0e-6,
            collisionless_justification="manufactured gate",
            claim_mode="predictive",
            model_validation_reference="manufactured validation reference",
            n_phase=16,
            normal_energy_bins=8,
            density_phase_count=256,
        )


def test_krueger_flux_deck_preserves_hpep_output_and_unresolved_ion_mixture():
    deck = load_krueger_2024_reactor_flux_deck(KRUEGER_DATA)

    assert deck.source_sha256 == (
        "ad50b6099a52d2c2cc00eb4eade496b9d75c41d19881c5fec9e905f9dfd3808b")
    assert deck.unresolved_species == ("ions",)
    assert not deck.supports_predictive_boundary
    assert deck.get("C3F4").flux_m2_s == 9.5e20
    assert deck.get("ions").flux_m2_s == 1.2e20
    assert deck.get("ions").charge_number is None


def test_krueger_neutral_subset_builds_without_laundering_missing_ions():
    deck = load_krueger_2024_reactor_flux_deck(KRUEGER_DATA)
    neutral_names = tuple(
        item.name for item in deck.species_fluxes if item.role == "neutral")
    boundary = build_tabulated_reactor_boundary(
        deck,
        reference_plane_m=1.0e-6,
        included_species=neutral_names,
        neutral_temperature_K=350.0,
    )

    assert {item.name for item in boundary.species} == set(neutral_names)
    assert boundary.current_density_A_m2 == 0.0
    assert boundary.provenance["complete_flux_deck_used"] is False
    assert boundary.provenance["omitted_species"] == ("ions",)
    assert boundary.provenance["unresolved_species_in_complete_deck"] == ("ions",)
    assert all(item.density_model is not None for item in boundary.species)


def test_tabulated_reactor_boundary_refuses_aggregate_ions_and_missing_kinetics():
    deck = load_krueger_2024_reactor_flux_deck(KRUEGER_DATA)
    with pytest.raises(ValueError, match="unresolved mixture"):
        build_tabulated_reactor_boundary(deck, reference_plane_m=1.0e-6)

    resolved = TabulatedReactorFluxDeck(
        species_fluxes=(
            ReactorSpeciesFlux(
                "Ar+", 1.2e20, "positive_ion", "published_distribution",
                "manufactured reactor output", charge_number=1, mass_amu=39.948),
        ),
        source="manufactured complete reactor deck",
        source_sha256="a" * 64,
    )
    with pytest.raises(ValueError, match="kinetic template"):
        build_tabulated_reactor_boundary(resolved, reference_plane_m=1.0e-6)


def test_tabulated_reactor_boundary_uses_template_shape_but_reactor_flux():
    record = ReactorSpeciesFlux(
        "Ar+", 4.0e20, "positive_ion", "published_distribution",
        "manufactured IEAD table", charge_number=1, mass_amu=39.948)
    deck = TabulatedReactorFluxDeck(
        species_fluxes=(record,),
        source="manufactured complete reactor deck",
        source_sha256="b" * 64,
    )
    template = SpeciesBoundaryState(
        name="Ar+", charge_number=1, mass_amu=39.948, flux_m2_s=1.0,
        velocity_sqrt_eV=np.array([[0.0, 0.0, 10.0], [1.0, 0.0, 9.0]]),
        weight=np.array([0.25, 0.75]),
        provenance={"distribution_source": "manufactured IEAD"},
    )
    boundary = build_tabulated_reactor_boundary(
        deck, reference_plane_m=1.0e-6, kinetic_templates={"Ar+": template})

    ion = boundary.get("Ar+")
    assert ion.flux_m2_s == 4.0e20
    np.testing.assert_array_equal(ion.velocity_sqrt_eV, template.velocity_sqrt_eV)
    np.testing.assert_array_equal(ion.weight, template.weight)
    assert ion.provenance["reactor_flux_deck_sha256"] == "b" * 64


def test_predictive_tabulated_boundary_requires_complete_predictive_evidence():
    krueger = load_krueger_2024_reactor_flux_deck(KRUEGER_DATA)
    neutral_names = tuple(
        item.name for item in krueger.species_fluxes if item.role == "neutral")
    with pytest.raises(ValueError, match="complete resolved deck"):
        build_tabulated_reactor_boundary(
            krueger,
            reference_plane_m=1.0e-6,
            included_species=neutral_names,
            claim_mode="predictive",
        )
