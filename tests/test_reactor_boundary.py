import numpy as np
import pytest
from pathlib import Path
import shutil

from petch.boundary_state import (
    DiscreteEnergyAngleDensity2D,
    DiscreteEnergyPolarAzimuthDensity3D,
    SpeciesBoundaryState,
)
from petch.amorphous_carbon_mask import AmorphousCarbonMaskParameters
from petch.reactor_boundary import (
    KRUEGER_2024_IEAD_CSV_SHA256,
    KRUEGER_2024_IEAD_METADATA_SHA256,
    KRUEGER_2024_TRANSFER_FLUX_SHA256,
    KRUEGER_2024_TRANSFER_IEAD_SHA256,
    KRUEGER_2024_TRANSFER_METADATA_SHA256,
    PlasmaDiagnosticState,
    ReactorSpeciesFlux,
    TabulatedReactorFluxDeck,
    append_global_current_balance_maxwellian_electrons,
    build_diagnostic_virtual_sheath_boundary,
    build_krueger_2024_development_boundary,
    build_krueger_2024_transfer_boundary,
    build_tabulated_reactor_boundary,
    load_krueger_2024_digitized_iead,
    load_krueger_2024_transfer_boundary_data,
    load_krueger_2024_reactor_flux_deck,
)
from petch.sheath import (
    CollisionlessRFSheath,
    CollisionlessWaveformSheath,
    PeriodicSheathVoltage,
)
from petch.surface_kinetics import (
    EnergeticFlux,
    ReducedSiO2FluorocarbonParameters,
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


def test_krueger_base_boundary_cannot_require_held_out_observations(tmp_path):
    """A base run succeeds from base boundary inputs while held-out tables are absent."""
    for name in (
            "base_case_boundary_fluxes.csv",
            "digitized_figure4_iead.csv",
            "digitized_figure4_iead_metadata.json"):
        shutil.copy2(KRUEGER_DATA / name, tmp_path / name)

    boundary = build_krueger_2024_development_boundary(
        tmp_path, reference_plane_m=2.8e-6,
        neutral_direction_polar_order=8,
        neutral_direction_azimuthal_order=16,
        ion_energy_bin_eV=250.0,
        ion_angle_bin_deg=0.25,
        ion_azimuthal_closure="axisymmetric_uniform",
        ion_azimuthal_order=16)

    assert boundary.provenance["provider"] == (
        "krueger_2024_published_flux_and_digitized_iead")
    assert not (tmp_path / "transfer_observations.csv").exists()


def test_krueger_digitized_iead_is_checksum_bound_and_preserves_joint_distribution():
    iead = load_krueger_2024_digitized_iead(KRUEGER_DATA)

    assert iead.table_sha256 == KRUEGER_2024_IEAD_CSV_SHA256
    assert iead.metadata_sha256 == KRUEGER_2024_IEAD_METADATA_SHA256
    assert np.isclose(iead.mean_energy_eV, 3465.110787140924)
    assert np.isclose(iead.probability_weight.sum(), 1.0)
    assert iead.energy_eV.min() > 400.0
    assert iead.energy_eV.max() < 4900.0
    assert np.max(np.abs(iead.signed_angle_deg)) < 4.0

    energetic = iead.energetic_flux(1.2e20)
    assert energetic.name == "ions"
    assert np.isclose(energetic.flux_m2_s, 1.2e20)
    assert np.isclose(
        np.dot(energetic.weight, energetic.energy_eV), iead.mean_energy_eV)


def test_krueger_transfer_boundary_tables_are_checksum_bound_and_complete():
    transfer = load_krueger_2024_transfer_boundary_data(KRUEGER_DATA)

    assert transfer.flux_table_sha256 == KRUEGER_2024_TRANSFER_FLUX_SHA256
    assert transfer.iead_table_sha256 == KRUEGER_2024_TRANSFER_IEAD_SHA256
    assert transfer.metadata_sha256 == KRUEGER_2024_TRANSFER_METADATA_SHA256
    assert set(transfer.flux_m2_s_by_oxygen_ratio) == {0.5, 1.0, 1.5, 2.5}
    assert set(transfer.iead_by_low_frequency_power_kw) == {0.0, 4.0, 6.0, 8.0}
    assert np.isclose(transfer.flux_m2_s_by_oxygen_ratio[0.5]["O"], 4.1e20)
    assert np.isclose(transfer.flux_m2_s_by_oxygen_ratio[2.5]["O"], 1.5e21)
    means = [
        transfer.iead_by_low_frequency_power_kw[power].mean_energy_eV
        for power in (0.0, 4.0, 6.0, 8.0)
    ]
    assert np.all(np.diff(means) > 0.0)
    assert np.allclose(means, [519.0, 2551.0, 2998.0, 3593.0], rtol=0.015)


def test_krueger_transfer_boundary_uses_only_published_inputs_for_oxygen_sweep():
    boundary = build_krueger_2024_transfer_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6,
        low_frequency_power_kw=6.0, oxygen_to_fluorocarbon_ratio=1.5,
        neutral_direction_polar_order=4, neutral_direction_azimuthal_order=8)
    transfer = load_krueger_2024_transfer_boundary_data(KRUEGER_DATA)

    assert boundary.provenance["held_out_profiles_used_to_construct_boundary"] is False
    assert boundary.provenance["power_sweep_constant_flux_closure"] is False
    assert boundary.get("O").flux_m2_s == (
        transfer.flux_m2_s_by_oxygen_ratio[1.5]["O"])
    assert boundary.get("ions").flux_m2_s == (
        transfer.flux_m2_s_by_oxygen_ratio[1.5]["ions"])
    assert all(
        item.weight.size == 32 for item in boundary.species if item.name != "ions")
    assert np.isclose(
        boundary.get("ions").mean_energy_eV,
        transfer.iead_by_low_frequency_power_kw[6.0].mean_energy_eV)


def test_krueger_power_transfer_declares_unpublished_flux_closure_and_refuses_mixing():
    boundary = build_krueger_2024_transfer_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6,
        low_frequency_power_kw=0.0,
        neutral_direction_polar_order=4, neutral_direction_azimuthal_order=8)
    base = load_krueger_2024_reactor_flux_deck(KRUEGER_DATA)

    assert boundary.provenance["power_sweep_constant_flux_closure"] is True
    assert boundary.get("ions").flux_m2_s == base.get("ions").flux_m2_s
    assert boundary.get("ions").mean_energy_eV < 600.0
    with pytest.raises(ValueError, match="belong to the 6 kW process"):
        build_krueger_2024_transfer_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            low_frequency_power_kw=4.0, oxygen_to_fluorocarbon_ratio=1.0)
    with pytest.raises(ValueError, match="must be one of"):
        build_krueger_2024_transfer_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            low_frequency_power_kw=7.0)


def test_krueger_development_boundary_uses_all_published_fluxes_without_claim_upgrade():
    boundary = build_krueger_2024_development_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6)
    ion = boundary.get("ions")

    assert {item.name for item in boundary.species} == {
        "C3F4", "C2F3", "CF", "CF2", "CF3", "O", "ions"}
    assert boundary.provenance["complete_published_flux_table_used"] is True
    assert boundary.provenance["supports_prediction"] is False
    assert boundary.provenance["aggregate_ion_mixture_unresolved"] is True
    assert ion.charge_number == 1
    assert ion.mass_amu == 39.948
    assert ion.flux_m2_s == 1.2e20
    assert np.isclose(ion.mean_energy_eV, 3465.110787140924)
    assert ion.density_model is None
    assert isinstance(ion.density_model_2d, DiscreteEnergyAngleDensity2D)
    sampled = ion.sample_flux_velocity_2d(np.array([[0.1], [0.9]]))
    assert sampled.shape == (2, 2)
    assert np.all(sampled[:, 1] > 0.0)
    assert "species-resolved" in boundary.provenance["predictive_blockers"][0]


def test_krueger_axisymmetric_lift_is_explicit_and_forward_sampleable_in_3d():
    boundary = build_krueger_2024_development_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6,
        ion_energy_bin_eV=500.0, ion_angle_bin_deg=0.5,
        ion_azimuthal_closure="axisymmetric_uniform")
    ion = boundary.get("ions")

    assert isinstance(ion.density_model, DiscreteEnergyPolarAzimuthDensity3D)
    assert boundary.provenance["ion_azimuthal_closure"] == "axisymmetric_uniform"
    assert boundary.provenance["ion_azimuthal_order"] == 16
    assert ion.weight.size == 16 * ion.density_model.probability_weight.size
    assert np.allclose(np.sum(ion.velocity_sqrt_eV ** 2, axis=1),
                       np.repeat(ion.density_model.energy_eV, 16))
    # A complete uniform azimuthal ring carries no preferred transverse direction.
    assert np.allclose(
        np.sum(ion.velocity_sqrt_eV[:, :2] * ion.weight[:, None], axis=0),
        0.0, atol=1e-13)
    sampled = ion.sample_flux_velocity(np.array([
        [0.1, 0.0], [0.1, 0.25], [0.9, 0.5], [0.9, 0.75]]))
    assert sampled.shape == (4, 3)
    assert np.all(sampled[:, 2] > 0.0)
    assert np.all(np.isfinite(sampled))
    with pytest.raises(ValueError, match="azimuthal_closure"):
        build_krueger_2024_development_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            ion_azimuthal_closure="invented")
    with pytest.raises(ValueError, match="azimuthal_order"):
        build_krueger_2024_development_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            ion_azimuthal_closure="axisymmetric_uniform",
            ion_azimuthal_order=0)


def test_krueger_compressed_quadrature_preserves_flux_moments_and_surface_yields():
    exact = build_krueger_2024_development_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6)
    compressed = build_krueger_2024_development_boundary(
        KRUEGER_DATA,
        reference_plane_m=1.0e-6,
        n_transverse_neutral=5,
        n_normal_neutral=2,
        ion_energy_bin_eV=500.0,
        ion_angle_bin_deg=0.5,
    )
    exact_ion = exact.get("ions")
    reduced_ion = compressed.get("ions")

    assert 64 <= reduced_ion.weight.size <= 120
    assert reduced_ion.weight.size < exact_ion.weight.size / 7
    assert all(
        item.weight.size == 50
        for item in compressed.species if item.charge_number == 0)
    assert compressed.provenance["total_boundary_quadrature_nodes"] < 450
    assert np.isclose(
        reduced_ion.mean_energy_eV, exact_ion.mean_energy_eV,
        rtol=2e-15, atol=0.0)
    compression = reduced_ion.provenance["numerical_quadrature"]
    assert compression["source_node_count"] == exact_ion.weight.size
    assert compression["node_count"] == reduced_ion.weight.size
    assert compression["second_energy_moment_relative_error"] < 0.002
    assert compression["mean_direction_maximum_absolute_error"] < 1e-5

    exact_flux = EnergeticFlux(
        "ions", 1.0, exact_ion.kinetic_energy_eV,
        exact_ion.velocity_sqrt_eV[:, 2]
        / np.linalg.norm(exact_ion.velocity_sqrt_eV, axis=1),
        exact_ion.weight)
    reduced_flux = EnergeticFlux(
        "ions", 1.0, reduced_ion.kinetic_energy_eV,
        reduced_ion.velocity_sqrt_eV[:, 2]
        / np.linalg.norm(reduced_ion.velocity_sqrt_eV, axis=1),
        reduced_ion.weight)
    oxide = ReducedSiO2FluorocarbonParameters.krueger_2024_reduced_projection()
    mask = AmorphousCarbonMaskParameters.krueger_2024_reduced_projection()
    laws = (
        oxide.bare_sio2_yield,
        oxide.complex_sio2_yield,
        oxide.polymer_sputter_yield,
        mask.polymer_sputter_yield,
        mask.carbon_sputter_yield,
    )
    for law in laws:
        reference = exact_flux.mean_yield(law)
        assert abs(reduced_flux.mean_yield(law) / reference - 1.0) < 6e-4


def test_krueger_direction_marginalized_neutral_rule_preserves_physical_moments():
    boundary = build_krueger_2024_development_boundary(
        KRUEGER_DATA,
        reference_plane_m=1.0e-6,
        neutral_direction_polar_order=4,
        neutral_direction_azimuthal_order=8,
    )
    neutrals = [item for item in boundary.species if item.charge_number == 0]

    assert len(neutrals) == 6
    assert all(item.weight.size == 32 for item in neutrals)
    assert boundary.provenance["neutral_quadrature"] == {
        "method": "analytic_speed_marginal_plus_angular_quadrature",
        "polar_order": 4,
        "azimuthal_order": 8,
        "nodes_per_species": 32,
        "validity_domain": (
            "field-free neutral transport with energy-independent neutral surface laws"),
    }
    for neutral in neutrals:
        direction = neutral.velocity_sqrt_eV / np.linalg.norm(
            neutral.velocity_sqrt_eV, axis=1)[:, None]
        first = np.einsum("s,sc->c", neutral.weight, direction)
        second = np.einsum("s,si,sj->ij", neutral.weight, direction, direction)
        assert np.allclose(first, [0.0, 0.0, 2.0 / 3.0], atol=2e-15, rtol=0.0)
        assert np.allclose(second, np.diag([0.25, 0.25, 0.5]), atol=2e-15, rtol=0.0)
        assert np.isclose(
            neutral.mean_energy_eV,
            2.0 * 300.0 * 8.617333262145e-5,
            atol=2e-17, rtol=0.0)
        numerical = neutral.provenance["numerical_quadrature"]
        assert numerical["analytically_marginalized_speed"] is True
        assert numerical["maximum_first_direction_moment_error"] < 2e-15
        assert numerical["maximum_second_direction_moment_error"] < 2e-15


def test_krueger_direction_marginalization_requires_both_valid_orders():
    with pytest.raises(ValueError, match="requires both"):
        build_krueger_2024_development_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            neutral_direction_polar_order=4)
    with pytest.raises(ValueError, match="polar order"):
        build_krueger_2024_development_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            neutral_direction_polar_order=1,
            neutral_direction_azimuthal_order=8)


def test_krueger_compression_refuses_implicit_or_invalid_bin_controls():
    iead = load_krueger_2024_digitized_iead(KRUEGER_DATA)
    with pytest.raises(ValueError, match="requires both"):
        iead.development_quadrature(energy_bin_eV=500.0)
    with pytest.raises(ValueError, match="positive and finite"):
        iead.development_quadrature(energy_bin_eV=0.0, angle_bin_deg=0.5)


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


def test_global_current_balance_electrons_close_krueger_boundary_without_local_law():
    heavy = build_krueger_2024_development_boundary(
        KRUEGER_DATA, reference_plane_m=2.65e-6,
        neutral_direction_polar_order=4,
        neutral_direction_azimuthal_order=8,
    )
    charged = append_global_current_balance_maxwellian_electrons(
        heavy,
        electron_temperature_eV=3.6,
        temperature_source="Krueger 2024 thesis, Fig. 6.3 and Sec. 2.2.2",
        temperature_evidence_kind="published_HPEM_output",
        n_transverse=3,
        n_normal=4,
    )

    ion = charged.get("ions")
    electron = charged.get("electron")
    assert electron.flux_m2_s == ion.flux_m2_s
    assert electron.density_model is not None
    assert electron.weight.size == 36
    assert abs(charged.current_density_A_m2) < 1e-12
    assert all(left is right for left, right in zip(charged.species[:-1], heavy.species))
    closure = charged.provenance["global_current_balance_electron_closure"]
    assert closure["local_balance_is_not_imposed"] is True
    assert closure["volume_boltzmann_electron_term"] is False
    assert electron.provenance["angular_distribution"] == (
        "Lambertian_cosine_flux_marginal")

    with pytest.raises(ValueError, match="already contains"):
        append_global_current_balance_maxwellian_electrons(
            charged,
            electron_temperature_eV=3.6,
            temperature_source="duplicate closure gate",
            temperature_evidence_kind="manufactured",
        )
