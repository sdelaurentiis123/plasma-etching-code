import numpy as np
import pytest
from scipy.constants import atomic_mass

from petch.reactor_global import (
    AxisymmetricRadiationZoneField,
    CylindricalReactor,
    ResonanceLineData,
    deterministic_cylinder_partial_redistribution,
    deterministic_zonal_partial_redistribution,
)


LINE = ResonanceLineData(
    wavelength_nm=104.821987,
    transition_probability_s_inv=5.32e8,
    lower_statistical_weight=1.0,
    upper_statistical_weight=3.0,
    absorber_mass_kg=39.948 * atomic_mass,
    source="nist-argon-persistent-lines",
)


def test_one_zone_reduces_to_homogeneous_partial_redistribution():
    field = AxisymmetricRadiationZoneField(
        radial_edges_m=np.array([0.0, 0.1125]),
        axial_edges_m=np.array([0.0, 0.12]),
        cell_zone_index=np.array([[0]]),
        gas_temperature_K=np.array([500.0]),
        absorber_density_m3=np.array([3.0e20]),
        emitter_density_m3=np.array([1.0]),
        source="one-zone reduction test",
    )
    zonal = deterministic_zonal_partial_redistribution(
        field,
        LINE,
        surface_quadrature_order=6,
        direction_quadrature_order=6,
        frequency_quadrature_order=24,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    homogeneous = deterministic_cylinder_partial_redistribution(
        CylindricalReactor(radius_m=0.1125, length_m=0.12),
        LINE,
        absorber_density_m3=3.0e20,
        gas_temperature_K=500.0,
        velocity_changing_collision_frequency_s_inv=0.0,
        geometry_quadrature_order=6,
        frequency_quadrature_order=24,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    assert zonal.trapping_factor == pytest.approx(
        homogeneous.partial_redistribution_trapping_factor, rel=1.0e-11)
    assert zonal.complete_frequency_redistribution_trapping_factor == (
        pytest.approx(
            homogeneous.complete_redistribution_trapping_factor,
            rel=1.0e-11,
        ))
    assert zonal.transition_probability_conservation_error_maximum < 1.0e-12
    assert zonal.escape_boundary_labels == (
        "lower_endcap_wafer_plane",
        "upper_endcap",
        "cylindrical_sidewall",
    )
    assert zonal.partial_redistribution_escape_boundary_probability.sum() == (
        pytest.approx(1.0, abs=2.0e-7))
    assert zonal.complete_redistribution_escape_boundary_probability.sum() == (
        pytest.approx(1.0, abs=2.0e-7))
    assert zonal.partial_redistribution_quench_probability == pytest.approx(0.0)
    assert zonal.complete_redistribution_quench_probability == pytest.approx(0.0)
    assert (
        zonal.partial_redistribution_wafer_escape_probability
        == pytest.approx(
            zonal.partial_redistribution_escape_boundary_probability[1],
            rel=2.0e-13,
        )
    )
    assert zonal.partial_redistribution_wafer_flux_m2_s(
        total_line_emission_rate_s=2.0e12,
        wafer_area_m2=0.02,
    ) == pytest.approx(
        1.0e14 * zonal.partial_redistribution_wafer_escape_probability)


def test_cold_shell_space_frequency_solution_is_deterministic_and_conservative():
    field = AxisymmetricRadiationZoneField(
        radial_edges_m=np.array([0.0, 0.075, 0.1125]),
        axial_edges_m=np.array([0.0, 0.03, 0.09, 0.12]),
        cell_zone_index=np.array([
            [1, 0, 1],
            [1, 1, 1],
        ]),
        gas_temperature_K=np.array([800.0, 400.0]),
        absorber_density_m3=np.array([1.5e20, 3.0e20]),
        emitter_density_m3=np.array([1.0, 0.02]),
        source="synthetic hot-core/cold-shell test",
    )
    kwargs = dict(
        surface_quadrature_order=6,
        direction_quadrature_order=6,
        frequency_quadrature_order=20,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    first = deterministic_zonal_partial_redistribution(field, LINE, **kwargs)
    second = deterministic_zonal_partial_redistribution(field, LINE, **kwargs)
    assert first.trapping_factor == second.trapping_factor
    assert first.complete_frequency_redistribution_trapping_factor == (
        second.complete_frequency_redistribution_trapping_factor)
    assert np.array_equal(
        first.initial_emission_zone_probability,
        second.initial_emission_zone_probability,
    )
    assert first.trapping_factor > 1.0
    assert first.zone_source_measure_relative_volume_error_maximum < 0.06
    assert first.transition_probability_conservation_error_maximum < 1.0e-12
    assert first.terminal_probability_conservation_error_maximum < 2.0e-7
    assert first.partial_redistribution_escape_boundary_probability.sum() == (
        pytest.approx(1.0, abs=2.0e-7))
    assert first.linear_solver_relative_residual < 1.0e-7
    assert first.initial_emission_zone_probability[0] > 0.9


def test_zone_field_is_immutable_and_rejects_missing_zone():
    field = AxisymmetricRadiationZoneField(
        radial_edges_m=np.array([0.0, 0.1]),
        axial_edges_m=np.array([0.0, 0.1]),
        cell_zone_index=np.array([[0]]),
        gas_temperature_K=np.array([500.0]),
        absorber_density_m3=np.array([1.0e20]),
        emitter_density_m3=np.array([1.0]),
        source="immutability test",
    )
    with pytest.raises(ValueError):
        field.cell_zone_index[0, 0] = 1
    with pytest.raises(ValueError):
        AxisymmetricRadiationZoneField(
            radial_edges_m=np.array([0.0, 0.1]),
            axial_edges_m=np.array([0.0, 0.1]),
            cell_zone_index=np.array([[0]]),
            gas_temperature_K=np.array([500.0, 600.0]),
            absorber_density_m3=np.array([1.0e20, 1.0e20]),
            emitter_density_m3=np.array([1.0, 1.0]),
            source="missing-zone test",
        )
    with pytest.raises(ValueError):
        AxisymmetricRadiationZoneField(
            radial_edges_m=np.array([0.0, 0.1]),
            axial_edges_m=np.array([0.0, 0.1]),
            cell_zone_index=np.array([[0]]),
            gas_temperature_K=np.array([[500.0]]),
            absorber_density_m3=np.array([1.0e20]),
            emitter_density_m3=np.array([1.0]),
            source="nonscalar-zone-vector test",
        )


def test_quenching_is_a_conservative_terminal_channel():
    field = AxisymmetricRadiationZoneField(
        radial_edges_m=np.array([0.0, 0.1125]),
        axial_edges_m=np.array([0.0, 0.12]),
        cell_zone_index=np.array([[0]]),
        gas_temperature_K=np.array([500.0]),
        absorber_density_m3=np.array([3.0e20]),
        emitter_density_m3=np.array([1.0]),
        source="one-zone quenching ledger test",
    )
    result = deterministic_zonal_partial_redistribution(
        field,
        LINE,
        quenching_collision_frequency_s_inv=LINE.transition_probability_s_inv,
        surface_quadrature_order=6,
        direction_quadrature_order=6,
        frequency_quadrature_order=20,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    assert result.partial_redistribution_quench_probability > 0.0
    assert result.complete_redistribution_quench_probability > 0.0
    assert (
        result.partial_redistribution_escape_boundary_probability.sum()
        + result.partial_redistribution_quench_probability
    ) == pytest.approx(1.0, abs=2.0e-7)
    assert (
        result.complete_redistribution_escape_boundary_probability.sum()
        + result.complete_redistribution_quench_probability
    ) == pytest.approx(1.0, abs=2.0e-7)
    with pytest.raises(ValueError):
        result.partial_redistribution_wafer_flux_m2_s(
            total_line_emission_rate_s=1.0,
            wafer_area_m2=0.0,
        )


def test_partial_wafer_disk_is_separated_from_lower_endcap_annulus():
    field = AxisymmetricRadiationZoneField(
        radial_edges_m=np.array([0.0, 0.1125]),
        axial_edges_m=np.array([0.0, 0.12]),
        cell_zone_index=np.array([[0]]),
        gas_temperature_K=np.array([500.0]),
        absorber_density_m3=np.array([0.0]),
        emitter_density_m3=np.array([1.0]),
        source="partial-wafer terminal ledger test",
    )
    result = deterministic_zonal_partial_redistribution(
        field,
        LINE,
        wafer_radius_m=0.05,
        surface_quadrature_order=8,
        direction_quadrature_order=6,
        frequency_quadrature_order=20,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    assert result.escape_boundary_labels == (
        "lower_endcap_wafer_plane",
        "lower_endcap_outside_wafer",
        "upper_endcap",
        "cylindrical_sidewall",
    )
    terminal = dict(zip(
        result.escape_boundary_labels,
        result.partial_redistribution_escape_boundary_probability,
    ))
    assert terminal["lower_endcap_wafer_plane"] > 0.0
    assert terminal["lower_endcap_outside_wafer"] > 0.0
    assert terminal["lower_endcap_wafer_plane"] < (
        terminal["lower_endcap_outside_wafer"])
    assert sum(terminal.values()) == pytest.approx(1.0, abs=2.0e-7)
