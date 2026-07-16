from pathlib import Path

import numpy as np
import pytest

from petch.experimental_boundary import (
    Jeong2023IonBoundaryClosure,
    build_jeong_2023_boundary_state,
)
from petch.experimental_data import (
    load_jeong_2023_etch_depths,
    load_jeong_2023_radical_densities,
)


DATA = Path(__file__).parents[1] / "data" / "experimental" / "jeong_2023"


def _evidence():
    depths = load_jeong_2023_etch_depths(DATA / "digitized_figure7_depths.csv")
    radicals = load_jeong_2023_radical_densities(
        DATA / "digitized_figure6_radicals.csv")
    anchor = next(item for item in depths if item.split == "calibration")
    support = tuple(item for item in radicals if item.electron_density_m3 == 1.9e15)
    return anchor, support, radicals


def test_jeong_2023_boundary_preserves_diagnostic_and_model_evidence_classes():
    anchor, support, _ = _evidence()
    boundary = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6)

    assert {item.name for item in boundary.species} == {"Ar+", "FC_total"}
    ion = boundary.get("Ar+")
    neutral = boundary.get("FC_total")
    assert np.allclose(ion.velocity_sqrt_eV[:, 2] ** 2, 890.0)
    assert ion.flux_m2_s > 0.0 and neutral.flux_m2_s > ion.flux_m2_s
    assert ion.provenance["self_bias_is_not_iedf"] is True
    assert neutral.provenance["source_evidence_type"] == "source_plasma_model_digitized"
    assert neutral.provenance["aggregate_species_collapse"] is True
    assert np.isclose(
        neutral.flux_m2_s,
        sum(neutral.provenance["species_thermal_flux_m2_s"].values()),
        rtol=1e-14)
    assert boundary.provenance["closure_supports_prediction"] is False
    assert boundary.provenance["ion_flux_mode"] == "species_density_fraction_bohm"
    assert boundary.provenance["positive_ion_species"] == ["Ar+"]


def test_jeong_2023_boundary_rejects_mixed_radical_conditions():
    anchor, support, radicals = _evidence()
    mixed = support[:-1] + (next(
        item for item in radicals
        if item.species == support[-1].species and item.electron_density_m3 == 3.1e15),)
    with pytest.raises(ValueError, match="share one source condition"):
        build_jeong_2023_boundary_state(anchor, mixed, reference_plane_m=3.0e-6)


def test_jeong_2023_boundary_partitions_etchant_and_polymer_without_losing_flux():
    anchor, support, _ = _evidence()
    aggregate = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6)
    partitioned = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6,
        radical_channel_mode="heavy_light")

    assert {item.name for item in partitioned.species} == {
        "Ar+", "FC_etchant", "FC_polymer"}
    etchant = partitioned.get("FC_etchant")
    polymer = partitioned.get("FC_polymer")
    assert etchant.provenance["grouped_species"] == ["CF", "CF2", "CF3"]
    assert polymer.provenance["grouped_species"] == ["C2F4", "C3F6", "C4F7"]
    assert np.isclose(
        etchant.flux_m2_s + polymer.flux_m2_s,
        aggregate.get("FC_total").flux_m2_s, rtol=2e-16)
    assert partitioned.provenance["radical_channel_mode"] == "heavy_light"


def test_jeong_2023_boundary_can_retain_every_plotted_radical_without_claiming_completeness():
    anchor, support, _ = _evidence()
    aggregate = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6)
    resolved = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6,
        radical_channel_mode="species_resolved")

    assert {item.name for item in resolved.species if item.charge_number == 0} == {
        "CF", "CF2", "CF3", "C2F4", "C3F6", "C4F7"}
    assert np.isclose(
        sum(item.flux_m2_s for item in resolved.species if item.charge_number == 0),
        aggregate.get("FC_total").flux_m2_s, rtol=2e-16)
    assert resolved.provenance["figure6_is_complete_reactive_inventory"] is False
    assert "atomic_F" in resolved.provenance["unreported_boundary_channels"]
    for species in resolved.species:
        if species.charge_number == 0:
            assert species.provenance["grouped_species"] == [species.name]


def test_huang_reduced_boundary_only_groups_reaction_equivalent_radicals():
    anchor, support, _ = _evidence()
    aggregate = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6)
    projected = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6,
        radical_channel_mode="huang_2019_reduced")

    assert {item.name for item in projected.species if item.charge_number == 0} == {
        "CF", "CF2", "FC_complex_02", "FC_polymer_heavy"}
    assert projected.get("FC_complex_02").provenance["grouped_species"] == [
        "CF3", "C3F6"]
    assert projected.get("FC_polymer_heavy").provenance["grouped_species"] == [
        "C2F4", "C4F7"]
    assert np.isclose(
        sum(item.flux_m2_s for item in projected.species if item.charge_number == 0),
        aggregate.get("FC_total").flux_m2_s, rtol=2e-16)


def test_jeong_2023_boundary_refuses_unknown_radical_partition():
    anchor, support, _ = _evidence()
    with pytest.raises(ValueError, match="radical_channel_mode"):
        build_jeong_2023_boundary_state(
            anchor, support, reference_plane_m=3.0e-6,
            radical_channel_mode="magic")


def test_jeong_2023_boundary_accepts_explicit_multispecies_bohm_mixture():
    anchor, support, _ = _evidence()
    closure = Jeong2023IonBoundaryClosure(
        species_mass_amu={
            "Ar+": 39.948,
            "CF+": 12.011 + 18.998403163,
            "C3F5+": 3.0 * 12.011 + 5.0 * 18.998403163,
        },
        species_density_fraction={"Ar+": 0.2, "CF+": 0.5, "C3F5+": 0.3},
        normal_energy_fraction_of_self_bias={
            "Ar+": 1.0, "CF+": 0.98, "C3F5+": 0.95},
        provenance={
            "source": "manufactured multispecies boundary gate",
            "evidence_kind": "assumed",
        },
    )
    boundary = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6, ion_closure=closure)

    ions = tuple(item for item in boundary.species if item.charge_number > 0)
    assert {item.name for item in ions} == {"Ar+", "CF+", "C3F5+"}
    electron_density = anchor.electron_density_m3
    electron_temperature_eV = 3.0
    elementary_charge_c = 1.602176634e-19
    atomic_mass_kg = 1.66053906660e-27
    for ion in ions:
        expected_flux = (
            electron_density * closure.species_density_fraction[ion.name]
            * np.sqrt(
                electron_temperature_eV * elementary_charge_c
                / (ion.mass_amu * atomic_mass_kg)))
        assert np.isclose(ion.flux_m2_s, expected_flux, rtol=2e-15)
        assert np.allclose(
            ion.velocity_sqrt_eV[:, 2] ** 2,
            anchor.self_bias_magnitude_v
            * closure.normal_energy_fraction_of_self_bias[ion.name])
        assert ion.provenance["supports_prediction"] is False
    assert np.isclose(
        boundary.provenance["positive_ion_total_flux_m2_s"],
        sum(item.flux_m2_s for item in ions), rtol=0.0, atol=0.0)


def test_jeong_2023_boundary_accepts_reactor_supplied_species_fluxes():
    anchor, support, _ = _evidence()
    closure = Jeong2023IonBoundaryClosure(
        species_mass_amu={"CF+": 31.009403163, "CF3+": 69.006209489},
        explicit_species_flux_m2_s={"CF+": 2.0e18, "CF3+": 7.5e17},
        normal_energy_fraction_of_self_bias={"CF+": 1.0, "CF3+": 1.0},
        provenance={
            "source": "manufactured validated reactor-model gate",
            "evidence_kind": "validated_reactor_model",
        },
        supports_prediction_within_declared_domain=True,
    )
    boundary = build_jeong_2023_boundary_state(
        anchor, support, reference_plane_m=3.0e-6, ion_closure=closure)

    assert boundary.get("CF+").flux_m2_s == 2.0e18
    assert boundary.get("CF3+").flux_m2_s == 7.5e17
    assert boundary.provenance["ion_flux_mode"] == "explicit_species_flux"
    assert boundary.provenance["ion_closure_supports_prediction"] is True
    # The complete boundary still cannot claim prediction because Figure 6 is a
    # selected source-model radical inventory rather than measured wall fluxes.
    assert boundary.provenance["closure_supports_prediction"] is False


@pytest.mark.parametrize(
    "kwargs",
    (
        {
            "species_mass_amu": {"Ar+": 39.948, "CF+": 31.0},
            "species_density_fraction": {"Ar+": 0.7, "CF+": 0.2},
            "normal_energy_fraction_of_self_bias": {"Ar+": 1.0, "CF+": 1.0},
        },
        {
            "species_mass_amu": {"Ar+": 39.948},
            "explicit_species_flux_m2_s": {"CF+": 1.0e18},
            "normal_energy_fraction_of_self_bias": {"Ar+": 1.0},
        },
    ),
)
def test_jeong_2023_ion_closure_refuses_incomplete_or_unnormalized_inputs(kwargs):
    with pytest.raises(ValueError, match="Jeong"):
        Jeong2023IonBoundaryClosure(**kwargs)
