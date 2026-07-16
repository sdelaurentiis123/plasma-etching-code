import numpy as np
import pytest

from petch.charging_stationarity_3d import (
    PROFILE_STATIONARITY_CONTRACT_DRAFT,
    ProfileChargingStationarityBlock3D,
    ProfileChargingStationarityContract3D,
    assess_profile_charging_stationarity_3d,
)


def _contract(**updates):
    values = dict(
        potential_drift_tolerance_v=0.1,
        current_relative_l1_tolerance=0.08,
        transported_flux_relative_l1_tolerance=0.08,
        profile_velocity_relative_l1_tolerance=0.08,
        profile_increment_tolerance_m=0.02,
        minimum_independent_replicates=2,
        confidence_multiplier=2.0)
    values.update(updates)
    return ProfileChargingStationarityContract3D(**values)


def _block(start, end, *, epoch=10, current_shift=0.0, flux_shift=0.0,
           velocity_shift=0.0, error=0.0, hard=True):
    positive = np.array([2.0, 1.0])
    negative = np.array([1.0 + current_shift, 2.0 - current_shift])
    velocity = np.array([1.0, -0.5]) + velocity_shift
    duration = 0.01
    return ProfileChargingStationarityBlock3D(
        potential_start_v=np.array([start, start]),
        potential_end_v=np.array([end, end]),
        positive_face_current_density_a_m2=positive,
        negative_face_current_density_a_m2=negative,
        net_face_current_standard_error_a_m2=np.full(2, error),
        face_area_m2=np.array([1.0, 2.0]),
        species_face_flux_m2_s={"ion": np.array([4.0, 2.0]) + flux_shift},
        species_face_flux_standard_error_m2_s={"ion": np.full(2, error)},
        profile_velocity_m_s=velocity,
        profile_velocity_standard_error_m_s=np.full(2, error),
        profile_increment_m=velocity * duration,
        independent_replicates=2,
        scoring_sampling_epochs=(epoch, epoch + 1),
        duration_s=duration,
        exact_hard_visibility=hard)


def test_profile_stationarity_passes_stable_independent_blocks():
    first = _block(0.0, 0.05, epoch=10)
    second = _block(0.05, 0.08, epoch=20)

    result = assess_profile_charging_stationarity_3d(first, second, _contract())

    assert result.passed
    assert result.contract_revision == PROFILE_STATIONARITY_CONTRACT_DRAFT
    assert not result.diagnostics["experimental_claim_authorized"]
    assert result.potential_drift_upper_v == pytest.approx(0.03)
    assert result.current_relative_l1_upper == pytest.approx(0.0)
    assert result.transported_flux_relative_l1_upper == pytest.approx(0.0)


def test_profile_stationarity_uses_uncertainty_as_an_upper_bound():
    first = _block(0.0, 0.01, epoch=10, error=0.1)
    second = _block(0.01, 0.02, epoch=20, error=0.1)

    result = assess_profile_charging_stationarity_3d(first, second, _contract())

    assert not result.passed
    assert "independent kinetic-current change exceeds tolerance" in result.reasons
    assert "delivered species-flux change exceeds tolerance" in result.reasons
    assert "predicted profile-velocity change exceeds tolerance" in result.reasons


def test_profile_stationarity_refuses_reused_scoring_samples_and_nonconsecutive_blocks():
    first = _block(0.0, 0.05, epoch=10)
    with pytest.raises(ValueError, match="disjoint"):
        assess_profile_charging_stationarity_3d(
            first, _block(0.05, 0.06, epoch=11), _contract())
    with pytest.raises(ValueError, match="consecutive"):
        assess_profile_charging_stationarity_3d(
            first, _block(0.06, 0.07, epoch=20), _contract())


def test_profile_stationarity_reports_each_physical_gate_without_softening():
    first = _block(0.0, 0.05, epoch=10)
    second = _block(
        0.05, 0.30, epoch=20, current_shift=0.5, flux_shift=1.0,
        velocity_shift=0.5, hard=False)

    result = assess_profile_charging_stationarity_3d(
        first, second, _contract(profile_increment_tolerance_m=0.001))

    assert not result.passed
    assert len(result.reasons) == 6
    assert "exact hard-visibility operator was not used in both blocks" in result.reasons
    assert "second-block potential drift exceeds tolerance" in result.reasons
