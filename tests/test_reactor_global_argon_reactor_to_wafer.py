import numpy as np
import pytest

from petch.reactor_global import (
    ArgonGlobalCondition,
    ArgonReactorToWaferCondition,
    ArgonTransportState,
    CylindricalReactor,
    DeterministicArgonCollisionalSheathTransfer,
    DeterministicArgonReactorToWaferModel,
    FixedArgonTransportProvider,
    LeeLiebermanArgonGlobalModel,
    PASCAL_PER_MTORR,
    maxwellian_floating_sheath_potential_eV,
)


def _condition():
    return ArgonReactorToWaferCondition(
        global_condition=ArgonGlobalCondition(
            condition_id="manufactured-full-Ar-stack",
            absorbed_power_W=500.0,
            pressure_Pa=8.0 * PASCAL_PER_MTORR,
            gas_temperature_K=500.0,
            geometry=CylindricalReactor(radius_m=0.15, length_m=0.075),
            ion_wall_energy_factor_Te=5.0,
            ion_wall_energy_source="Lee--Lieberman published range member",
            ion_wall_energy_evidence="published_range_member",
            absorbed_power_source="manufactured absorbed bulk power",
            absorbed_power_evidence="assumed",
            absorbed_power_boundary_kind="manufactured_test",
        ),
        delivered_bias_power_W=20.0,
        bias_frequency_hz=2.0e6,
    )


def _model():
    transport = FixedArgonTransportProvider(ArgonTransportState(
        ion_mean_free_path_m=0.004,
        ambipolar_diffusion_m2_s=50.0,
        metastable_effective_diffusion_m2_s=1.0,
        source="manufactured Ar transport",
        evidence_kind="assumed",
    ))
    collisions = DeterministicArgonCollisionalSheathTransfer(
        initial_thermal_radial_order=1,
        initial_thermal_azimuth_order=2,
        position_quadrature_order=3,
        hazard_quadrature_order=4,
        impact_quadrature_order=2,
        collision_azimuth_order=2,
        maximum_collision_order=1,
        steps_per_period=64,
        steps_per_transit=64,
    )
    return DeterministicArgonReactorToWaferModel(
        global_model=LeeLiebermanArgonGlobalModel(transport),
        collisional_transfer=collisions,
        sheath_phase_count=16,
        sheath_steps_per_period=64,
        sheath_steps_per_transit=64,
    )


def test_maxwellian_floating_potential_has_expected_argon_scale():
    assert maxwellian_floating_sheath_potential_eV(4.0) == pytest.approx(
        18.718, rel=2.0e-4)


def test_absorbed_knobs_close_to_common_collisional_wafer_boundary():
    first = _model().solve(_condition())
    second = _model().solve(_condition())
    assert first.global_plasma.maximum_normalized_residual < 1.0e-8
    assert first.wafer.collisionless.power_closure_relative_residual < 2.0e-8
    assert first.wafer.collisional.mean_total_optical_depth > 0.0
    assert first.maximum_conservation_residual < 2.0e-8
    assert first.boundary.get("Ar+").flux_m2_s > 0.0
    np.testing.assert_array_equal(
        first.boundary.get("Ar+").velocity_sqrt_eV,
        second.boundary.get("Ar+").velocity_sqrt_eV,
    )
    assert not first.supports_equipment_prediction
    assert not first.supports_feature_depth
    assert first.provenance["generator_forward_power_inversion_closed"] is False
    assert first.provenance["ion_collision_order_closed"] is True
    assert first.provenance["fast_neutral_wafer_transport_closed"] is False
