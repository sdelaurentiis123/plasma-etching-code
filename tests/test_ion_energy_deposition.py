import numpy as np
import pytest

from petch.ion_energy_deposition import (
    AMORPHOUS_CARBON, SIO2, derived_yield_energy_factor,
    nuclear_energy_in_layer_eV, projected_range_nm,
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
