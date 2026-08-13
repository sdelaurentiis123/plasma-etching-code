import numpy as np
import pytest

from petch.reactor_global import (
    ArgonBornMayerPhelpsCollisionModel,
    DeterministicArgonCollisionalSheathTransfer,
    DeterministicCollisionalRFSheath,
    DiagnosticConditionedRFSheathTransfer,
)
from petch.sheath import CollisionlessWaveformSheath, PeriodicSheathVoltage


BOHR_RADIUS_M = 5.29177210903e-11
BOLTZMANN_J_K = 1.380649e-23


def _static_sheath(voltage_v=1000.0, thickness_m=0.01):
    waveform = PeriodicSheathVoltage.sinusoidal(
        dc_v=voltage_v,
        amplitude_v=0.0,
        frequency_hz=2.0e6,
        source="manufactured static collisional-sheath gate",
    )
    return CollisionlessWaveformSheath(
        waveform=waveform,
        Te_eV=4.0,
        ion_mass_amu=39.948,
        thickness_m=thickness_m,
    )


def _model(
    gas_number_density_m3,
    *,
    maximum_collision_order=1,
    phase_count=4,
    initial_thermal_radial_order=1,
    initial_thermal_azimuth_order=4,
    position_quadrature_order=5,
    hazard_quadrature_order=6,
    impact_quadrature_order=2,
    collision_azimuth_order=2,
):
    return DeterministicCollisionalRFSheath(
        sheath=_static_sheath(),
        collision_model=ArgonBornMayerPhelpsCollisionModel(),
        gas_number_density_m3=gas_number_density_m3,
        neutral_gas_temperature_K=500.0,
        source_ion_flux_m2_s=1.0e19,
        phase_count=phase_count,
        initial_thermal_radial_order=initial_thermal_radial_order,
        initial_thermal_azimuth_order=initial_thermal_azimuth_order,
        position_quadrature_order=position_quadrature_order,
        hazard_quadrature_order=hazard_quadrature_order,
        impact_quadrature_order=impact_quadrature_order,
        collision_azimuth_order=collision_azimuth_order,
        maximum_collision_order=maximum_collision_order,
        steps_per_period=64,
        steps_per_transit=128,
    )


def test_argon_born_mayer_source_inference_reproduces_independent_geometry():
    model = ArgonBornMayerPhelpsCollisionModel()
    angle = model.center_of_mass_scattering_angle_rad(
        1000.0, 5.4 * BOHR_RADIUS_M)
    assert np.rad2deg(angle) == pytest.approx(0.1, abs=2.0e-10)
    # Khrabrov & Kaganovich independently print R_cx ~= 8.09 a0.  It was not
    # used to infer the Born-Mayer potential parameters.
    assert model.collision_radius_m(1000.0) / BOHR_RADIUS_M == pytest.approx(
        8.09, rel=0.04
    )
    theta = np.asarray([
        model.center_of_mass_scattering_angle_rad(
            1000.0, fraction * model.collision_radius_m(1000.0))
        for fraction in (0.1, 0.3, 0.6, 0.9)
    ])
    assert np.all(np.diff(theta) < 0.0)


def test_charge_label_swap_reproduces_published_half_angle_kinematics():
    from petch.reactor_global.collisional_sheath import (
        _equal_mass_collision_velocities,
    )

    theta = np.deg2rad(4.0)
    projectile, target = _equal_mass_collision_velocities(
        np.array([0.0, 0.0, np.sqrt(1000.0)]), theta, 0.0)
    # The fast neutral after charge exchange is the projectile momentum
    # branch and therefore exits at theta/2 in the lab frame.
    neutral_angle = np.arctan2(
        np.linalg.norm(projectile[:2]), projectile[2])
    ion_angle = np.arctan2(np.linalg.norm(target[:2]), target[2])
    assert neutral_angle == pytest.approx(0.5 * theta, abs=2.0e-15)
    assert ion_angle == pytest.approx(
        0.5 * (np.pi - theta), abs=2.0e-15)
    assert np.dot(projectile, projectile) + np.dot(
        target, target) == pytest.approx(1000.0, abs=2.0e-12)


def test_collisionless_limit_recovers_finite_transit_sheath_and_thermal_core():
    result = _model(0.0).solve()
    assert result.ion_arrival_probability == pytest.approx(1.0, abs=2e-12)
    assert result.uncollided_arrival_probability == pytest.approx(
        1.0, abs=2e-12)
    assert result.unresolved_probability == pytest.approx(0.0, abs=1e-15)
    assert result.escaped_probability == pytest.approx(0.0, abs=1e-15)
    assert result.expected_collision_count_lower_bound == pytest.approx(0.0)
    assert result.distribution.mean_energy_eV == pytest.approx(
        result.collisionless_reference_mean_normal_energy_eV
        + BOLTZMANN_J_K * 500.0 / 1.602176634e-19,
        abs=0.04,
    )
    assert result.probability_ledger_relative_residual < 2.0e-12
    assert result.maximum_resolved_energy_ledger_relative_residual < 2.0e-10


def test_one_collision_expansion_is_deterministic_conservative_and_broadens():
    density = 1.33322 / (BOLTZMANN_J_K * 500.0)  # 10 mTorr at 500 K
    first = _model(density).solve()
    second = _model(density).solve()
    np.testing.assert_array_equal(
        first.distribution.velocity_sqrt_eV,
        second.distribution.velocity_sqrt_eV,
    )
    np.testing.assert_array_equal(
        first.distribution.weight, second.distribution.weight)
    assert first.mean_total_optical_depth > 1.0
    assert 0.0 < first.unresolved_probability < 1.0
    assert first.ion_arrival_probability < 1.0
    assert first.expected_charge_exchange_count_lower_bound > 0.0
    assert (
        first.expected_fast_neutral_birth_energy_lower_bound_eV_per_source_ion
        > 0.0
    )
    assert first.resolved_fast_neutral_arrivals_per_source_ion > 0.0
    assert first.resolved_fast_neutral_distribution is not None
    assert first.resolved_fast_neutral_flux_m2_s > 0.0
    assert (
        first.resolved_fast_neutral_arrivals_per_source_ion
        + first.unresolved_fast_neutral_collisions_per_source_ion
        + first.escaped_fast_neutrals_per_source_ion
        == pytest.approx(
            first.expected_fast_neutral_birth_count_lower_bound,
            abs=2.0e-12,
        )
    )
    assert first.fast_neutral_lineage_ledger_relative_residual < 2.0e-11
    assert np.sqrt(
        first.distribution.mean_squared_polar_angle_rad2
    ) > np.deg2rad(0.1)
    assert first.probability_ledger_relative_residual < 2.0e-11
    assert first.maximum_resolved_energy_ledger_relative_residual < 2.0e-10
    assert not first.supports_fast_neutral_wafer_flux
    assert not first.supports_feature_depth


def test_adding_a_collision_order_reduces_only_the_unresolved_mass():
    density = 0.35 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    order_one = _model(density, maximum_collision_order=1).solve()
    order_two = _model(density, maximum_collision_order=2).solve()
    assert order_two.unresolved_probability < order_one.unresolved_probability
    assert order_two.ion_arrival_probability > order_one.ion_arrival_probability
    assert order_two.probability_ledger_relative_residual < 2.0e-11


def test_neutral_density_jvp_matches_centered_finite_difference():
    density = 0.12 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    direction = 0.35 * density
    model = _model(density)
    result, tangent = model.density_jvp(direction)
    step = 2.0e-5
    plus = _model(density + step * direction).solve()
    minus = _model(density - step * direction).solve()

    def finite(name):
        return (getattr(plus, name) - getattr(minus, name)) / (2.0 * step)

    assert tangent.ion_arrival_probability_tangent == pytest.approx(
        finite("ion_arrival_probability"), rel=2.0e-5, abs=2.0e-9)
    assert tangent.unresolved_probability_tangent == pytest.approx(
        finite("unresolved_probability"), rel=2.0e-5, abs=2.0e-9)
    assert tangent.expected_collision_count_tangent == pytest.approx(
        finite("expected_collision_count_lower_bound"),
        rel=2.0e-5,
        abs=2.0e-9,
    )
    finite_mean = (
        plus.distribution.mean_energy_eV
        - minus.distribution.mean_energy_eV
    ) / (2.0 * step)
    assert tangent.mean_impact_energy_tangent_eV == pytest.approx(
        finite_mean, rel=3.0e-5, abs=2.0e-5)
    assert tangent.probability_ledger_tangent_residual < 2.0e-10
    assert tangent.distribution_weight_tangent.shape == result.distribution.weight.shape


def test_resolved_ions_export_through_common_boundary_without_promotion():
    density = 0.05 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    result = _model(density).solve()
    boundary = result.to_boundary_state(
        ion_name="Ar+", ion_mass_amu=39.948)
    ion = boundary.get("Ar+")
    neutral = boundary.get("Ar_fast_neutral")
    assert ion.flux_m2_s == pytest.approx(result.arriving_ion_flux_m2_s)
    assert ion.mean_energy_eV == pytest.approx(
        result.distribution.mean_energy_eV)
    assert neutral.flux_m2_s == pytest.approx(
        result.resolved_fast_neutral_flux_m2_s)
    assert neutral.charge_number == 0
    assert boundary.provenance["fast_neutral_boundary_is_lower_bound"] is True
    assert boundary.provenance["supports_feature_depth"] is False


def test_density_jvp_is_defined_at_the_collisionless_boundary():
    direction = 0.03 * 1.33322 / (BOLTZMANN_J_K * 500.0)
    _, tangent = _model(0.0).density_jvp(direction)
    step = 1.0e-5
    plus = _model(step * direction).solve()
    base = _model(0.0).solve()
    forward = (
        plus.uncollided_arrival_probability
        - base.uncollided_arrival_probability
    ) / step
    assert tangent.uncollided_arrival_probability_tangent == pytest.approx(
        forward, rel=2.0e-5, abs=2.0e-7)
    assert tangent.expected_collision_count_tangent > 0.0


def test_power_closed_argon_projection_lifts_without_changing_legacy_path():
    legacy = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Ar+": 39.948},
        electrode_area_m2=np.pi * 0.15 ** 2,
        plasma_potential_eV=18.0,
        frequency_hz=2.0e6,
        phase_count=16,
        steps_per_period=64,
        steps_per_transit=128,
        source="manufactured power-closed Ar gate",
    ).predict(
        positive_ion_flux_m2_s={"Ar+": 1.0e19},
        electron_temperature_eV=4.0,
        electron_density_m3=1.0e17,
        delivered_bias_power_W=20.0,
    )
    lift = DeterministicArgonCollisionalSheathTransfer(
        solver_kind="collision_order_reference",
        initial_thermal_radial_order=1,
        initial_thermal_azimuth_order=2,
        position_quadrature_order=3,
        hazard_quadrature_order=4,
        impact_quadrature_order=2,
        collision_azimuth_order=2,
        maximum_collision_order=1,
        steps_per_period=64,
        steps_per_transit=128,
    ).project(
        legacy,
        pressure_Pa=0.2 * 1.33322,
        gas_temperature_K=500.0,
    )
    assert lift.collisional.mean_total_optical_depth > 0.0
    assert lift.collisionless is legacy
    assert lift.collisionless.power_closure_relative_residual < 2.0e-8
    assert lift.collisionless_reference_relative_residual < 5.0e-4
    assert not lift.supports_feature_depth
    assert lift.to_boundary_state().get("Ar+").flux_m2_s > 0.0


def test_argon_projection_rejects_cross_chemistry_borrowing():
    wrong_species = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"CF3+": 69.0},
        electrode_area_m2=0.01,
        plasma_potential_eV=10.0,
        frequency_hz=1.0e6,
        phase_count=16,
        steps_per_period=64,
        steps_per_transit=64,
    ).project_from_bias_dc_component(
        positive_ion_flux_m2_s={"CF3+": 1.0e19},
        electron_temperature_eV=3.0,
        electron_density_m3=1.0e17,
        bias_dc_component_v=100.0,
    )
    with pytest.raises(ValueError, match="Ar\\+"):
        DeterministicArgonCollisionalSheathTransfer().project(
            wrong_species,
            pressure_Pa=1.0,
            gas_temperature_K=400.0,
        )
