import numpy as np
import pytest

from petch.reactor_global import (
    ArgonBornMayerPhelpsCollisionModel,
    ArgonCurrentDrivenReactorToWaferCondition,
    ArgonGlobalCondition,
    ArgonTransportState,
    CylindricalReactor,
    DeterministicCurrentDrivenArgonReactorToWaferModel,
    FixedArgonTransportProvider,
    LeeLiebermanArgonGlobalModel,
    PASCAL_PER_MTORR,
    PeriodicCurrentDensity,
)
from petch.sheath import bohm_speed


def _condition(*, current_evidence="assumed"):
    return ArgonCurrentDrivenReactorToWaferCondition(
        global_condition=ArgonGlobalCondition(
            condition_id="manufactured-current-driven-Ar-stack",
            absorbed_power_W=500.0,
            pressure_Pa=2.0 * PASCAL_PER_MTORR,
            gas_temperature_K=500.0,
            geometry=CylindricalReactor(radius_m=0.15, length_m=0.075),
            ion_wall_energy_factor_Te=5.0,
            ion_wall_energy_source="manufactured validated wall-energy model",
            ion_wall_energy_evidence="validated_model",
            absorbed_power_source="manufactured measured absorbed power",
            absorbed_power_evidence="measured",
            absorbed_power_boundary_kind="manufactured_test",
        ),
        sheath_current_density=PeriodicCurrentDensity(
            fundamental_frequency_hz=2.0e6,
            harmonic_number=np.array([1, 2]),
            # This amplitude gives a few-hundred-volt sheath for the
            # reactor-derived edge density.  A tiny manufactured current
            # creates a sub-volt, sub-micron sheath and pathological trapped
            # low-energy test ordinates that are not representative of an
            # etch bias.
            sine_A_m2=np.array([-40.0, 4.5]),
            cosine_A_m2=np.array([0.0, -2.5]),
            source="manufactured de-embedded electrode current",
            evidence_kind=current_evidence,
        ),
    )


def _model(*, transport_evidence="validated_model"):
    transport = FixedArgonTransportProvider(ArgonTransportState(
        ion_mean_free_path_m=0.004,
        ambipolar_diffusion_m2_s=50.0,
        metastable_effective_diffusion_m2_s=1.0,
        source="manufactured validated Ar transport",
        evidence_kind=transport_evidence,
    ))
    return DeterministicCurrentDrivenArgonReactorToWaferModel(
        global_model=LeeLiebermanArgonGlobalModel(transport),
        collision_model=ArgonBornMayerPhelpsCollisionModel(),
        sheath_phase_quadrature_count=1024,
        phase_node_count=4,
        position_node_count=3,
        total_energy_node_count=3,
        transverse_fraction_node_count=3,
        initial_thermal_radial_order=1,
        output_azimuth_order=2,
        impact_quadrature_order=1,
        collision_azimuth_order=2,
        collision_event_quadrature_order=2,
        steps_per_period=32,
        steps_per_transit=32,
        maximum_transit_periods=12.0,
    )


def test_global_bohm_flux_drives_moving_collisional_wafer_boundary():
    first = _model().solve(_condition())
    second = _model().solve(_condition())

    expected_density = first.global_plasma.axial_ion_flux_m2_s / bohm_speed(
        first.global_plasma.electron_temperature_eV, 39.948)
    assert first.sheath_edge_density_m3 == pytest.approx(expected_density)
    assert first.bohm_flux_seam_relative_residual < 2.0e-15
    assert first.maximum_conservation_residual < 2.0e-8
    assert first.wafer.source_ion_flux_m2_s == pytest.approx(
        first.global_plasma.axial_ion_flux_m2_s)
    assert first.boundary.get("Ar+").flux_m2_s == pytest.approx(
        first.wafer.arriving_ion_flux_m2_s)
    np.testing.assert_array_equal(
        first.boundary.get("Ar+").velocity_sqrt_eV,
        second.boundary.get("Ar+").velocity_sqrt_eV,
    )
    assert first.provenance["moving_electron_front_resolved"] is True
    assert first.provenance["RF_phase_is_kinetic_state"] is True
    assert first.provenance["ion_collision_order_closed"] is True
    assert first.provenance["fast_neutral_transport_closed"] is False
    assert not first.supports_feature_depth
    assert not first.supports_equipment_prediction


def test_prediction_evidence_is_fail_closed_at_each_reactor_seam():
    assumed_current = _model().solve(_condition(current_evidence="assumed"))
    assert assumed_current.global_plasma.supports_prediction
    assert not assumed_current.sheath.current.supports_predictive_boundary
    assert not assumed_current.supports_resolved_ion_boundary_prediction

    assumed_transport = _model(transport_evidence="assumed").solve(
        _condition(current_evidence="measured_sheath_current"))
    assert not assumed_transport.global_plasma.supports_prediction
    assert assumed_transport.sheath.current.supports_predictive_boundary
    assert not assumed_transport.supports_resolved_ion_boundary_prediction

    evidenced_seams = _model().solve(
        _condition(current_evidence="measured_sheath_current"))
    assert evidenced_seams.global_plasma.supports_prediction
    assert evidenced_seams.sheath.current.supports_predictive_boundary
    assert (
        evidenced_seams.wafer
        .below_born_mayer_support_collision_probability_lower_bound
        > 0.0
    )
    assert (
        evidenced_seams.provenance[
            "low_energy_ion_angular_scattering_closed"]
        is True
    )
    assert evidenced_seams.supports_resolved_ion_boundary_prediction


def test_current_driven_stack_rejects_generator_setpoint_substitution():
    with pytest.raises(ValueError, match="current-driven Ar reactor condition"):
        ArgonCurrentDrivenReactorToWaferCondition(
            global_condition=_condition().global_condition,
            sheath_current_density=20.0,
        )
