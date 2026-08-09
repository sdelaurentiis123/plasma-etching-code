import math

import pytest

from petch.reactor_global import (
    Adf04Level,
    ChlorineVuvLine,
    CylindricalReactor,
    chlorine_ground_population_fraction,
    deterministic_uniform_chlorine_line_wafer_boundary,
)


LEVELS = (
    Adf04Level(1, "ground", 2, "P", 1.5, 0.0),
    Adf04Level(2, "ground", 2, "P", 0.5, 882.35),
    Adf04Level(12, "excited", 2, "D", 2.5, 84120.0),
)
LINE = ChlorineVuvLine(
    upper_observed_index=12,
    lower_observed_index=1,
    wavelength_nm=118.877,
    transition_probability_s_inv=5.21e8,
    upper_total_radiative_probability_s_inv=6.0e8,
    direct_excitation_rate_coefficient_cm3_s=2.0e-11,
    photon_rate_coefficient_cm3_s=1.8e-11,
)


def test_ground_fine_structure_population_is_thermal_and_normalized():
    first = chlorine_ground_population_fraction(
        LEVELS, 1, gas_temperature_K=333.0)
    second = chlorine_ground_population_fraction(
        LEVELS, 2, gas_temperature_K=333.0)
    assert first + second == pytest.approx(1.0)
    assert first > second
    with pytest.raises(ValueError, match="ground-terminating"):
        chlorine_ground_population_fraction(
            LEVELS, 12, gas_temperature_K=333.0)


def test_uniform_chlorine_line_boundary_conserves_source_to_partial_wafer():
    geometry = CylindricalReactor(radius_m=0.215, length_m=0.065)
    result = deterministic_uniform_chlorine_line_wafer_boundary(
        geometry,
        LINE,
        LEVELS,
        wafer_radius_m=0.1,
        electron_density_m3=2.0e17,
        chlorine_atom_density_m3=3.0e19,
        gas_temperature_K=333.0,
        surface_quadrature_order=6,
        direction_quadrature_order=6,
        frequency_quadrature_order=20,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    expected_emissivity = 2.0e17 * 3.0e19 * 1.8e-11 * 1.0e-6
    assert result.primary_line_emissivity_m3_s == pytest.approx(
        expected_emissivity)
    assert result.primary_line_emission_rate_s == pytest.approx(
        expected_emissivity * geometry.volume_m3)
    assert result.alternate_branch_loss_frequency_s_inv == pytest.approx(
        7.9e7)
    assert result.wafer_photon_flux_m2_s == pytest.approx(
        result.primary_line_emission_rate_s
        * result.radiation.partial_redistribution_wafer_escape_probability
        / (math.pi * 0.1 ** 2)
    )
    assert result.radiation.escape_boundary_labels[1] == (
        "lower_endcap_outside_wafer")
    assert (
        result.radiation.partial_redistribution_escape_boundary_probability.sum()
        + result.radiation.partial_redistribution_quench_probability
    ) == pytest.approx(1.0, abs=2.0e-6)
    assert result.radiation.partial_redistribution_quench_probability > 0.0
    assert result.wafer_photon_flux_m2_s > 0.0
    assert result.prediction_supported is False

    refined = deterministic_uniform_chlorine_line_wafer_boundary(
        geometry,
        LINE,
        LEVELS,
        wafer_radius_m=0.1,
        electron_density_m3=2.0e17,
        chlorine_atom_density_m3=3.0e19,
        gas_temperature_K=333.0,
        surface_quadrature_order=10,
        direction_quadrature_order=10,
        frequency_quadrature_order=20,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    assert refined.radiation.partial_redistribution_wafer_escape_probability == (
        pytest.approx(
            result.radiation.partial_redistribution_wafer_escape_probability,
            rel=3.0e-2,
        )
    )


def test_zero_coronal_source_returns_zero_wafer_flux_without_breaking_transport():
    result = deterministic_uniform_chlorine_line_wafer_boundary(
        CylindricalReactor(radius_m=0.1, length_m=0.1),
        LINE,
        LEVELS,
        wafer_radius_m=0.1,
        electron_density_m3=0.0,
        chlorine_atom_density_m3=0.0,
        gas_temperature_K=333.0,
        surface_quadrature_order=6,
        direction_quadrature_order=6,
        frequency_quadrature_order=20,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    assert result.primary_line_emission_rate_s == 0.0
    assert result.wafer_photon_flux_m2_s == 0.0
