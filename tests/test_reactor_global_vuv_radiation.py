import numpy as np
import pytest
from scipy.constants import atomic_mass

from petch.reactor_global import (
    CylindricalReactor,
    KemaneciCl139nmEmissionSensitivity,
    ResonanceLineData,
    deterministic_cylinder_partial_redistribution,
    deterministic_cylinder_resonance_escape,
    uniform_isotropic_cylinder_to_disk_transfer,
    uniform_isotropic_cylinder_direct_lamellar_floor_transfer,
)


AR_1048 = ResonanceLineData(
    wavelength_nm=104.821987,
    transition_probability_s_inv=5.32e8,
    lower_statistical_weight=1.0,
    upper_statistical_weight=3.0,
    absorber_mass_kg=39.948 * atomic_mass,
    source="nist-argon-persistent-lines 1048.21987 A",
)
AR_1067 = ResonanceLineData(
    wavelength_nm=106.665980,
    transition_probability_s_inv=1.32e8,
    lower_statistical_weight=1.0,
    upper_statistical_weight=3.0,
    absorber_mass_kg=39.948 * atomic_mass,
    source="nist-argon-persistent-lines 1066.65980 A",
)


def test_cylinder_disk_transfer_is_linear_and_probability_bounded():
    geometry = CylindricalReactor(radius_m=0.2, length_m=0.08)
    first = uniform_isotropic_cylinder_to_disk_transfer(
        geometry, wafer_radius_m=0.1, volume_emissivity_m3_s=2.0e20,
        quadrature_order=24,
    )
    second = uniform_isotropic_cylinder_to_disk_transfer(
        geometry, wafer_radius_m=0.1, volume_emissivity_m3_s=6.0e20,
        quadrature_order=24,
    )
    assert second.wafer_photon_flux_m2_s == pytest.approx(
        3.0 * first.wafer_photon_flux_m2_s, rel=2e-15)
    assert second.wafer_intercept_probability == pytest.approx(
        first.wafer_intercept_probability, rel=2e-15)
    assert 0.0 < first.wafer_intercept_probability < 0.5
    assert first.wafer_photon_flux_m2_s == pytest.approx(
        first.volume_emissivity_m3_s * first.geometry_flux_length_m)


def test_absorption_is_monotonic_and_quadrature_converges():
    geometry = CylindricalReactor(radius_m=0.215, length_m=0.065)
    clear20 = uniform_isotropic_cylinder_to_disk_transfer(
        geometry, wafer_radius_m=0.1, volume_emissivity_m3_s=1.0,
        quadrature_order=20,
    )
    clear28 = uniform_isotropic_cylinder_to_disk_transfer(
        geometry, wafer_radius_m=0.1, volume_emissivity_m3_s=1.0,
        quadrature_order=28,
    )
    absorbed = uniform_isotropic_cylinder_to_disk_transfer(
        geometry, wafer_radius_m=0.1, volume_emissivity_m3_s=1.0,
        extinction_coefficient_m_inv=20.0, quadrature_order=28,
    )
    assert clear28.geometry_flux_length_m == pytest.approx(
        clear20.geometry_flux_length_m, rel=2.5e-3)
    assert absorbed.wafer_photon_flux_m2_s < clear28.wafer_photon_flux_m2_s


def test_lamellar_floor_direct_receipt_is_conserved_linear_and_depth_monotone():
    geometry = CylindricalReactor(radius_m=0.215, length_m=0.065)
    opening = 310.0e-9
    open_floor = uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
        geometry,
        wafer_radius_m=0.1,
        opening_width_m=opening,
        optical_path_depth_m=0.0,
        volume_emissivity_m3_s=2.0e20,
        radial_quadrature_order=24,
        axial_quadrature_order=24,
    )
    shallow = uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
        geometry,
        wafer_radius_m=0.1,
        opening_width_m=opening,
        optical_path_depth_m=200.0e-9,
        volume_emissivity_m3_s=2.0e20,
        radial_quadrature_order=24,
        axial_quadrature_order=24,
    )
    deep = uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
        geometry,
        wafer_radius_m=0.1,
        opening_width_m=opening,
        optical_path_depth_m=660.0e-9,
        volume_emissivity_m3_s=6.0e20,
        radial_quadrature_order=24,
        axial_quadrature_order=24,
    )
    assert open_floor.direct_floor_fraction == 1.0
    assert 0.0 < deep.direct_floor_fraction < shallow.direct_floor_fraction < 1.0
    assert deep.wafer_plane_photon_flux_m2_s == pytest.approx(
        3.0 * open_floor.wafer_plane_photon_flux_m2_s, rel=2.0e-15
    )
    assert deep.trench_floor_photon_flux_m2_s == pytest.approx(
        deep.wafer_plane_photon_flux_m2_s * deep.direct_floor_fraction
    )


def test_lamellar_floor_quadrature_converges_and_zero_source_retains_geometry():
    geometry = CylindricalReactor(radius_m=0.215, length_m=0.065)
    common = dict(
        wafer_radius_m=0.1,
        opening_width_m=310.0e-9,
        optical_path_depth_m=500.0e-9,
        volume_emissivity_m3_s=1.0,
    )
    low = uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
        geometry,
        **common,
        radial_quadrature_order=24,
        axial_quadrature_order=24,
    )
    high = uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
        geometry,
        **common,
        radial_quadrature_order=40,
        axial_quadrature_order=40,
    )
    zero = uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
        geometry,
        **{**common, "volume_emissivity_m3_s": 0.0},
        radial_quadrature_order=24,
        axial_quadrature_order=24,
    )
    assert high.direct_floor_fraction == pytest.approx(
        low.direct_floor_fraction, rel=7.0e-3
    )
    assert zero.direct_floor_fraction == pytest.approx(low.direct_floor_fraction)
    assert zero.trench_floor_photon_flux_m2_s == 0.0


def test_atomic_chlorine_excitation_card_uses_printed_fit_and_refuses_extrapolation():
    model = KemaneciCl139nmEmissionSensitivity(
        radiative_survival_fraction=0.01)
    rate = model.primary_excitation_rate_m3_s(
        electron_density_m3=1.0e17,
        chlorine_atom_density_m3=2.0e19,
        electron_temperature_eV=3.0,
    )
    assert rate > 0.0
    assert model.escaping_emissivity_m3_s(
        electron_density_m3=1.0e17,
        chlorine_atom_density_m3=2.0e19,
        electron_temperature_eV=3.0,
    ) == pytest.approx(0.01 * rate)
    assert not model.supports_prediction
    with pytest.raises(ValueError, match="0.5 <= Te <= 10"):
        model.excitation_rate_coefficient_m3_s(11.0)


def test_resonance_atomic_data_reproduce_nist_argon_oscillator_strengths():
    assert AR_1048.absorption_oscillator_strength == pytest.approx(
        0.262903, rel=2.0e-6)
    assert AR_1067.absorption_oscillator_strength == pytest.approx(
        0.0675468, rel=2.0e-6)
    assert AR_1048.natural_lorentz_hwhm_hz > (
        AR_1067.natural_lorentz_hwhm_hz)


def test_deterministic_resonance_escape_has_thin_limit_and_density_monotonicity():
    geometry = CylindricalReactor(radius_m=0.1125, length_m=0.12)
    common = dict(
        gas_temperature_K=400.0,
        geometry_quadrature_order=8,
        frequency_quadrature_order=24,
    )
    vacuum = deterministic_cylinder_resonance_escape(
        geometry, AR_1067, absorber_density_m3=0.0, **common)
    dilute = deterministic_cylinder_resonance_escape(
        geometry, AR_1067, absorber_density_m3=3.0e19, **common)
    dense = deterministic_cylinder_resonance_escape(
        geometry, AR_1067, absorber_density_m3=3.0e20, **common)
    assert vacuum.escape_probability == 1.0
    assert vacuum.trapping_factor == 1.0
    assert 1.0 < dilute.trapping_factor < dense.trapping_factor
    assert dense.frequency_profile_normalization == pytest.approx(
        1.0, abs=2.0e-9)
    assert dense.line_center_mean_path_optical_depth > 100.0


def test_top_localized_resonance_source_reduces_trapping_without_randomness():
    geometry = CylindricalReactor(radius_m=0.1125, length_m=0.12)
    common = dict(
        absorber_density_m3=3.0e20,
        gas_temperature_K=400.0,
        geometry_quadrature_order=8,
        frequency_quadrature_order=24,
    )
    uniform = deterministic_cylinder_resonance_escape(
        geometry, AR_1048, **common)
    localized = deterministic_cylinder_resonance_escape(
        geometry, AR_1048, source_axial_scale_length_m=0.005, **common)
    repeated = deterministic_cylinder_resonance_escape(
        geometry, AR_1048, source_axial_scale_length_m=0.005, **common)
    assert localized.source_weighted_mean_exit_path_m < (
        uniform.source_weighted_mean_exit_path_m)
    assert localized.trapping_factor < uniform.trapping_factor
    assert repeated == localized


def test_partial_redistribution_recovers_complete_and_coherent_limits():
    geometry = CylindricalReactor(radius_m=0.1125, length_m=0.12)
    common = dict(
        absorber_density_m3=3.0e20,
        gas_temperature_K=500.0,
        geometry_quadrature_order=8,
        frequency_quadrature_order=24,
    )
    coherent = deterministic_cylinder_partial_redistribution(
        geometry,
        AR_1067,
        velocity_changing_collision_frequency_s_inv=0.0,
        **common,
    )
    redistributed = deterministic_cylinder_partial_redistribution(
        geometry,
        AR_1067,
        velocity_changing_collision_frequency_s_inv=1.0e18,
        **common,
    )
    quenched = deterministic_cylinder_partial_redistribution(
        geometry,
        AR_1067,
        velocity_changing_collision_frequency_s_inv=1.0e18,
        quenching_collision_frequency_s_inv=1.0e18,
        **common,
    )
    assert coherent.partial_redistribution_trapping_factor >= 1.0
    assert coherent.partial_redistribution_trapping_factor != pytest.approx(
        coherent.complete_redistribution_trapping_factor, rel=1.0e-3)
    assert coherent.linear_solver_relative_residual < 1.0e-7
    assert coherent.coherent_grid_points_per_lorentz_hwhm >= 12.0
    assert redistributed.partial_redistribution_trapping_factor == pytest.approx(
        redistributed.complete_redistribution_trapping_factor, rel=2.0e-9)
    assert quenched.partial_redistribution_trapping_factor == pytest.approx(
        1.0, rel=2.0e-9)


def test_partial_redistribution_thin_limit_is_one_and_is_deterministic():
    geometry = CylindricalReactor(radius_m=0.1125, length_m=0.12)
    kwargs = dict(
        absorber_density_m3=0.0,
        gas_temperature_K=500.0,
        velocity_changing_collision_frequency_s_inv=2.0e8,
        quenching_collision_frequency_s_inv=3.0e7,
        geometry_quadrature_order=8,
        frequency_quadrature_order=24,
    )
    first = deterministic_cylinder_partial_redistribution(
        geometry, AR_1048, **kwargs)
    second = deterministic_cylinder_partial_redistribution(
        geometry, AR_1048, **kwargs)
    assert first == second
    assert first.complete_redistribution_trapping_factor == pytest.approx(1.0)
    assert first.partial_redistribution_trapping_factor == pytest.approx(1.0)


def test_invalid_radiation_conditions_fail_closed():
    geometry = CylindricalReactor(radius_m=0.2, length_m=0.1)
    with pytest.raises(ValueError):
        uniform_isotropic_cylinder_to_disk_transfer(
            geometry, wafer_radius_m=0.21, volume_emissivity_m3_s=1.0)
    with pytest.raises(ValueError):
        uniform_isotropic_cylinder_to_disk_transfer(
            geometry, wafer_radius_m=0.1, volume_emissivity_m3_s=1.0,
            extinction_coefficient_m_inv=-1.0)
