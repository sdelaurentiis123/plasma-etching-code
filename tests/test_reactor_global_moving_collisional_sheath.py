import numpy as np
import pytest

from petch.reactor_global import (
    ArgonBornMayerPhelpsCollisionModel,
    DeterministicMovingCollisionalRFSheath,
    PeriodicCurrentDensity,
    TurnerChabertCurrentDrivenSheath,
    certify_moving_sheath_convergence,
)


BOLTZMANN_J_K = 1.380649e-23


def _moving_sheath():
    return TurnerChabertCurrentDrivenSheath(
        current=PeriodicCurrentDensity(
            fundamental_frequency_hz=2.0e6,
            harmonic_number=np.array([1]),
            sine_A_m2=np.array([-4.0]),
            cosine_A_m2=np.array([0.0]),
            source="manufactured moving collisional-sheath gate",
        ),
        electron_temperature_eV=4.0,
        ion_mass_amu=39.948,
        sheath_edge_density_m3=1.0e16,
        phase_quadrature_count=1024,
    )


def _model(
    density,
    *,
    phase_node_count=4,
    position_node_count=3,
    total_energy_node_count=3,
    transverse_fraction_node_count=3,
    steps_per_period=32,
    steps_per_transit=32,
):
    return DeterministicMovingCollisionalRFSheath(
        sheath=_moving_sheath(),
        collision_model=ArgonBornMayerPhelpsCollisionModel(),
        gas_number_density_m3=density,
        neutral_gas_temperature_K=500.0,
        source_ion_flux_m2_s=1.0e19,
        phase_node_count=phase_node_count,
        position_node_count=position_node_count,
        total_energy_node_count=total_energy_node_count,
        transverse_fraction_node_count=transverse_fraction_node_count,
        initial_thermal_radial_order=1,
        output_azimuth_order=2,
        impact_quadrature_order=1,
        collision_azimuth_order=2,
        collision_event_quadrature_order=2,
        steps_per_period=steps_per_period,
        steps_per_transit=steps_per_transit,
        maximum_transit_periods=12.0,
    )


def test_collisionless_limit_recovers_moving_sheath_reference():
    result = _model(0.0).solve()
    assert result.ion_arrival_probability == pytest.approx(1.0, abs=2.0e-12)
    assert result.escaped_probability == pytest.approx(0.0, abs=2.0e-12)
    assert result.unresolved_probability == pytest.approx(0.0, abs=2.0e-12)
    assert result.expected_collision_count_lower_bound == pytest.approx(0.0)
    assert result.distribution.mean_energy_eV == pytest.approx(
        result.collisionless_reference_mean_normal_energy_eV
        + BOLTZMANN_J_K * 500.0 / 1.602176634e-19,
        rel=5.0e-2,
        abs=0.15,
    )
    assert result.provenance["RF_phase_is_kinetic_state"] is True
    assert result.provenance["moving_sheath_self_consistency_closed"] is True
    assert result.provenance["charge_exchange_birth_phase_resolved"] is True
    assert result.provenance["linear_solve_relative_residual"] < 2.0e-11
    assert not result.supports_feature_depth


def test_finite_density_closes_all_ion_collision_orders_and_broadens():
    density = 0.2 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    first = _model(density).solve()
    second = _model(density).solve()
    np.testing.assert_array_equal(
        first.distribution.velocity_sqrt_eV,
        second.distribution.velocity_sqrt_eV)
    np.testing.assert_array_equal(first.distribution.weight, second.distribution.weight)
    assert first.expected_collision_count_lower_bound > 0.0
    assert first.unresolved_probability < 2.0e-11
    assert (
        first.ion_arrival_probability + first.escaped_probability
        == pytest.approx(1.0, abs=2.0e-11)
    )
    assert first.expected_charge_exchange_count_lower_bound > 0.0
    assert first.resolved_fast_neutral_arrivals_per_source_ion > 0.0
    assert first.probability_ledger_relative_residual < 2.0e-11
    assert first.maximum_resolved_energy_ledger_relative_residual < 2.0e-10
    assert first.provenance["ion_collision_order_closed"] is True


def test_density_jvp_matches_centered_finite_difference():
    density = 0.04 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    direction = 0.3 * density
    result, tangent = _model(density).density_jvp(direction)
    step = 2.0e-5
    plus = _model(density + step * direction).solve()
    minus = _model(density - step * direction).solve()

    def finite(name):
        return (getattr(plus, name) - getattr(minus, name)) / (2.0 * step)

    assert tangent.ion_arrival_probability_tangent == pytest.approx(
        finite("ion_arrival_probability"), rel=2.0e-3, abs=2.0e-7)
    assert tangent.expected_collision_count_tangent == pytest.approx(
        finite("expected_collision_count_lower_bound"),
        rel=2.0e-3,
        abs=2.0e-7,
    )
    finite_energy = (
        plus.distribution.mean_energy_eV
        - minus.distribution.mean_energy_eV
    ) / (2.0 * step)
    assert tangent.mean_impact_energy_tangent_eV == pytest.approx(
        finite_energy, rel=3.0e-3, abs=2.0e-4)
    assert tangent.distribution_weight_tangent.shape == result.distribution.weight.shape


def test_moving_sheath_convergence_receipt_is_fail_closed():
    density = 0.05 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    coarse = _model(density, steps_per_period=24, steps_per_transit=24).solve()
    fine = _model(density, steps_per_period=48, steps_per_transit=48).solve()
    receipt = certify_moving_sheath_convergence(
        coarse,
        fine,
        mean_energy_relative_limit=0.2,
        rms_angle_relative_limit=0.2,
        collision_count_relative_limit=0.2,
        arrival_probability_relative_limit=0.1,
    )
    assert receipt.passed
    assert receipt.probability_ledger_residual < 2.0e-11


def test_boundary_export_retains_physical_open_gates():
    density = 0.03 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    result = _model(density).solve()
    boundary = result.to_boundary_state(ion_name="Ar+", ion_mass_amu=39.948)
    assert boundary.get("Ar+").flux_m2_s == pytest.approx(
        result.arriving_ion_flux_m2_s)
    assert boundary.provenance["supports_feature_depth"] is False
    assert boundary.provenance["fast_neutral_boundary_is_lower_bound"] is True
    assert result.provenance["supports_generator_power_inversion"] is False
    assert result.provenance["fast_neutral_transport_closed"] is False
