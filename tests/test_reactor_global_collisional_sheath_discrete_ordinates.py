import numpy as np
import pytest

from petch.reactor_global import (
    ArgonBornMayerPhelpsCollisionModel,
    DeterministicCollisionalRFSheath,
    DeterministicDiscreteOrdinatesRFSheath,
    certify_discrete_ordinates_convergence,
)
from petch.sheath import CollisionlessWaveformSheath, PeriodicSheathVoltage


BOLTZMANN_J_K = 1.380649e-23


def _sheath(voltage_v=1000.0):
    return CollisionlessWaveformSheath(
        waveform=PeriodicSheathVoltage.sinusoidal(
            dc_v=voltage_v,
            amplitude_v=0.0,
            frequency_hz=2.0e6,
            source="manufactured discrete-ordinates sheath gate",
        ),
        Te_eV=4.0,
        ion_mass_amu=39.948,
        thickness_m=0.01,
    )


def _model(
    density,
    *,
    potential_node_count=7,
    total_energy_node_count=7,
    transverse_fraction_node_count=10,
):
    return DeterministicDiscreteOrdinatesRFSheath(
        sheath=_sheath(),
        collision_model=ArgonBornMayerPhelpsCollisionModel(),
        gas_number_density_m3=density,
        neutral_gas_temperature_K=500.0,
        source_ion_flux_m2_s=1.0e19,
        phase_count=4,
        initial_thermal_radial_order=1,
        initial_thermal_azimuth_order=2,
        potential_node_count=potential_node_count,
        total_energy_node_count=total_energy_node_count,
        transverse_fraction_node_count=transverse_fraction_node_count,
        position_quadrature_order=3,
        hazard_quadrature_order=4,
        impact_quadrature_order=2,
        collision_azimuth_order=2,
        steps_per_period=64,
        steps_per_transit=128,
    )


def test_collisionless_limit_is_exact_and_deterministic():
    first = _model(0.0).solve()
    second = _model(0.0).solve()
    assert first.ion_arrival_probability == pytest.approx(1.0, abs=2e-12)
    assert first.unresolved_probability == pytest.approx(0.0, abs=2e-12)
    assert first.expected_collision_count_lower_bound == pytest.approx(0.0)
    np.testing.assert_array_equal(
        first.distribution.velocity_sqrt_eV,
        second.distribution.velocity_sqrt_eV,
    )
    np.testing.assert_array_equal(
        first.distribution.weight, second.distribution.weight)
    assert first.provenance["ion_collision_order_closed"] is True
    assert first.provenance["linear_solve_relative_residual"] < 2.0e-12


def test_high_optical_depth_closes_infinite_ion_collision_orders():
    density = 1.33322 / (BOLTZMANN_J_K * 500.0)
    result = _model(density).solve()
    assert result.mean_total_optical_depth > 1.0
    assert result.unresolved_probability < 2.0e-11
    assert (
        result.ion_arrival_probability + result.escaped_probability
        == pytest.approx(1.0, abs=2.0e-11)
    )
    assert result.expected_collision_count_lower_bound > 1.0
    assert result.provenance["ion_backscatter_turning_resolved"] is True
    assert result.provenance["maximum_row_probability_residual"] < 2.0e-11
    assert not result.supports_fast_neutral_wafer_flux
    assert not result.supports_feature_depth


def test_implicit_density_jvp_matches_centered_finite_difference():
    density = 0.2 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    direction = 0.3 * density
    result, tangent = _model(
        density,
        potential_node_count=5,
        total_energy_node_count=5,
        transverse_fraction_node_count=7,
    ).density_jvp(direction)
    step = 1.0e-5
    plus = _model(
        density + step * direction,
        potential_node_count=5,
        total_energy_node_count=5,
        transverse_fraction_node_count=7,
    ).solve()
    minus = _model(
        density - step * direction,
        potential_node_count=5,
        total_energy_node_count=5,
        transverse_fraction_node_count=7,
    ).solve()

    def finite(name):
        return (getattr(plus, name) - getattr(minus, name)) / (2.0 * step)

    assert tangent.expected_collision_count_tangent == pytest.approx(
        finite("expected_collision_count_lower_bound"),
        rel=3.0e-5,
        abs=2.0e-8,
    )
    assert tangent.ion_arrival_probability_tangent == pytest.approx(
        finite("ion_arrival_probability"), rel=4.0e-5, abs=2.0e-8)
    finite_energy = (
        plus.distribution.mean_energy_eV
        - minus.distribution.mean_energy_eV
    ) / (2.0 * step)
    assert tangent.mean_impact_energy_tangent_eV == pytest.approx(
        finite_energy, rel=5.0e-5, abs=2.0e-5)
    assert tangent.distribution_weight_tangent.shape == result.distribution.weight.shape


def test_low_density_agrees_with_explicit_collision_order_reference():
    density = 0.03 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    implicit = _model(density).solve()
    explicit = DeterministicCollisionalRFSheath(
        sheath=_sheath(),
        collision_model=ArgonBornMayerPhelpsCollisionModel(),
        gas_number_density_m3=density,
        neutral_gas_temperature_K=500.0,
        source_ion_flux_m2_s=1.0e19,
        phase_count=4,
        initial_thermal_radial_order=1,
        initial_thermal_azimuth_order=2,
        position_quadrature_order=3,
        hazard_quadrature_order=4,
        impact_quadrature_order=2,
        collision_azimuth_order=2,
        maximum_collision_order=3,
        steps_per_period=64,
        steps_per_transit=128,
    ).solve()
    assert explicit.unresolved_probability < 1.0e-6
    assert implicit.distribution.mean_energy_eV == pytest.approx(
        explicit.distribution.mean_energy_eV, rel=2.0e-3)
    assert implicit.expected_collision_count_lower_bound == pytest.approx(
        explicit.expected_collision_count_lower_bound, rel=2.0e-2)


def test_convergence_receipt_grades_all_transport_moments():
    density = 1.33322 / (BOLTZMANN_J_K * 500.0)
    coarse = _model(
        density,
        potential_node_count=7,
        total_energy_node_count=7,
        transverse_fraction_node_count=10,
    ).solve()
    fine = _model(
        density,
        potential_node_count=9,
        total_energy_node_count=9,
        transverse_fraction_node_count=13,
    ).solve()
    receipt = certify_discrete_ordinates_convergence(
        coarse,
        fine,
        mean_energy_relative_limit=1.0e-2,
        rms_angle_relative_limit=2.0e-2,
        collision_count_relative_limit=1.0e-2,
        neutral_arrival_relative_limit=1.0e-2,
    )
    assert receipt.passed
    assert receipt.ion_transport_closed
    assert receipt.mean_energy_relative_change < 1.0e-2
    assert receipt.ion_probability_ledger_residual < 2.0e-11


def test_boundary_export_stays_explicit_about_neutral_lower_bound():
    density = 0.2 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    result = _model(density).solve()
    boundary = result.to_boundary_state(ion_name="Ar+", ion_mass_amu=39.948)
    assert boundary.get("Ar+").flux_m2_s == pytest.approx(
        result.arriving_ion_flux_m2_s)
    assert boundary.provenance["supports_feature_depth"] is False
    assert boundary.provenance["fast_neutral_boundary_is_lower_bound"] is True
