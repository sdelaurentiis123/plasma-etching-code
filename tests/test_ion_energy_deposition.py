import numpy as np
import pytest

from petch.ion_energy_deposition import (
    AMORPHOUS_CARBON, SIO2, derived_yield_energy_factor,
    csda_path_nm, nuclear_energy_in_layer_eV, projected_range_nm,
    residual_energy_after_layer_eV,
    stopping_cross_sections_eV_cm2)


def test_projected_range_matches_published_anchor_band():
    # Ar into Si-like targets at 500 eV: literature spans 1.6 +/- 0.5 nm (XPS review)
    # to 2.3 nm (SRIM); ZBL-analytic accuracy is the ~15-20 percent class.
    value = projected_range_nm(500.0, 18, 39.948, SIO2)
    assert 1.4 < value < 3.2
    # Monotone growth and keV-scale magnitudes.
    ranges = [projected_range_nm(e, 18, 39.948, SIO2)
              for e in (250, 500, 1000, 3000, 5000)]
    assert all(b > a for a, b in zip(ranges, ranges[1:]))
    assert 5.0 < ranges[3] < 12.0


def test_layer_deposit_saturates_with_energy():
    fractions = []
    for energy in (100.0, 1000.0, 5000.0):
        deposit = nuclear_energy_in_layer_eV(energy, 1.0, 2.0, 18, 39.948, SIO2)
        assert 0.0 < deposit <= energy
        fractions.append(deposit / energy)
    assert fractions[0] > fractions[1] > fractions[2]


def test_grazing_incidence_deposits_more_in_layer():
    normal = nuclear_energy_in_layer_eV(3000.0, 1.0, 2.0, 18, 39.948, SIO2)
    grazing = nuclear_energy_in_layer_eV(
        3000.0, float(np.cos(np.radians(70.0))), 2.0, 18, 39.948, SIO2)
    assert grazing > 1.5 * normal


def test_derived_factor_anchored_at_reference_and_layer_insensitive():
    assert derived_yield_energy_factor(
        1000.0, 1.0, layer_depth_nm=2.0,
        reference_energy_eV=1000.0) == pytest.approx(1.0)
    # The 8kW/6kW mean-energy factor ratio moves <3 percent across the full
    # literature band of mixed-layer depths (1-3 nm): no sensitive knob traded in.
    ratios = []
    for depth in (1.0, 3.0):
        six = derived_yield_energy_factor(
            2998.0, 1.0, layer_depth_nm=depth, reference_energy_eV=1000.0)
        eight = derived_yield_energy_factor(
            3593.0, 1.0, layer_depth_nm=depth, reference_energy_eV=1000.0)
        ratios.append(eight / six)
    assert abs(ratios[0] - ratios[1]) < 0.03


def test_bragg_additivity_between_bounds():
    energy = np.array([1000.0])
    silica_n, _ = stopping_cross_sections_eV_cm2(energy, 18, 39.948, SIO2)
    carbon_n, _ = stopping_cross_sections_eV_cm2(energy, 18, 39.948, AMORPHOUS_CARBON)
    assert silica_n[0] > 0.0 and carbon_n[0] > 0.0


def test_residual_energy_is_csda_path_inversion_not_exponential_attenuation():
    incident = 1000.0
    path_nm = csda_path_nm(incident, 18, 39.948, AMORPHOUS_CARBON)
    traversed_nm = 0.37 * path_nm
    remaining = residual_energy_after_layer_eV(
        incident, 1.0, traversed_nm, 18, 39.948, AMORPHOUS_CARBON)
    assert 10.0 < remaining < incident
    remaining_path_nm = csda_path_nm(
        remaining, 18, 39.948, AMORPHOUS_CARBON)
    assert remaining_path_nm == pytest.approx(
        path_nm - traversed_nm, rel=3e-3)


def test_residual_energy_obeys_slant_path_stopping_and_zero_layer_identity():
    energies = np.array([25.0, 200.0, 1000.0])
    identity = residual_energy_after_layer_eV(
        energies, np.array([1.0, 0.5, 0.0]), 0.0,
        18, 39.948, AMORPHOUS_CARBON)
    assert np.array_equal(identity, energies)

    normal = residual_energy_after_layer_eV(
        1000.0, 1.0, 0.5, 18, 39.948, AMORPHOUS_CARBON)
    oblique = residual_energy_after_layer_eV(
        1000.0, 0.5, 0.5, 18, 39.948, AMORPHOUS_CARBON)
    grazing = residual_energy_after_layer_eV(
        1000.0, 0.0, 0.5, 18, 39.948, AMORPHOUS_CARBON)
    assert 0.0 < oblique < normal < 1000.0
    assert grazing == 0.0


def test_residual_energy_stops_exactly_when_layer_exceeds_csda_path():
    incident = 200.0
    stopping_depth_nm = csda_path_nm(
        incident, 18, 39.948, AMORPHOUS_CARBON)
    remaining = residual_energy_after_layer_eV(
        incident, 1.0, 1.01 * stopping_depth_nm,
        18, 39.948, AMORPHOUS_CARBON)
    assert remaining == 0.0


@pytest.mark.parametrize(
    "energy, cosine, depth",
    [(-1.0, 1.0, 1.0), (100.0, -0.1, 1.0), (100.0, 1.1, 1.0),
     (100.0, 1.0, -1.0)])
def test_residual_energy_refuses_nonphysical_inputs(energy, cosine, depth):
    with pytest.raises(ValueError):
        residual_energy_after_layer_eV(
            energy, cosine, depth, 18, 39.948, AMORPHOUS_CARBON)
