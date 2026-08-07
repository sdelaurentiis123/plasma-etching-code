import numpy as np
import pytest

from petch.reactor_global import (
    ArgonGlobalCondition,
    ArgonTransportState,
    CylindricalReactor,
    DirectDriveRFPowerBoundary,
    FixedArgonTransportProvider,
    LeeLiebermanArgonGlobalModel,
    MatchedRFPowerBoundary,
    MeasuredAbsorbedPowerBoundary,
    time_average_real_power_W,
)


def _geometry():
    return CylindricalReactor(radius_m=0.10, length_m=0.10)


def _condition_from_power(power):
    return ArgonGlobalCondition.from_power_estimate(
        condition_id="power-boundary-test",
        power_estimate=power,
        pressure_Pa=1.0,
        gas_temperature_K=600.0,
        geometry=_geometry(),
        ion_wall_energy_factor_Te=5.0,
        ion_wall_energy_source="test measured wall closure",
        ion_wall_energy_evidence="measured",
    )


def _transport():
    return FixedArgonTransportProvider(ArgonTransportState(
        ion_mean_free_path_m=0.004,
        ambipolar_diffusion_m2_s=50.0,
        metastable_effective_diffusion_m2_s=1.0,
        source="test measured transport closure",
        evidence_kind="measured",
    ))


def test_direct_absorbed_power_measurement_preserves_uncertainty_and_evidence():
    estimate = MeasuredAbsorbedPowerBoundary(
        absorbed_power_W=80.0,
        absolute_uncertainty_W=4.0,
        source="calorimeter-1",
        method="plasma-on minus plasma-off water calorimetry",
    ).estimate()
    assert estimate.lower_W == 76.0
    assert estimate.upper_W == 84.0
    assert estimate.point_W == 80.0
    assert estimate.supports_prediction
    assert estimate.evidence_kind == "measured"

    solution = LeeLiebermanArgonGlobalModel(_transport()).solve(
        _condition_from_power(estimate))
    assert solution.absorbed_power_evidence == "measured"
    assert "calorimeter-1" in solution.absorbed_power_source
    assert solution.absorbed_power_boundary_kind == (
        "direct_absorbed_power_measurement")
    assert solution.supports_prediction


def test_matched_rf_power_does_not_become_absorbed_power_without_loss_closure():
    estimate = MatchedRFPowerBoundary(
        forward_power_W=100.0,
        reflected_power_W=4.0,
        hardware_loss_lower_W=10.0,
        hardware_loss_upper_W=30.0,
        measurement_source="directional coupler",
        loss_source="unresolved match plus coil loss interval",
    ).estimate()
    assert estimate.lower_W == 66.0
    assert estimate.upper_W == 86.0
    assert estimate.point_W is None
    assert not estimate.supports_prediction
    assert estimate.evidence_kind == "unresolved"
    with pytest.raises(ValueError, match="no point estimate"):
        _condition_from_power(estimate)


def test_measured_matched_rf_loss_chain_supports_prediction():
    estimate = MatchedRFPowerBoundary(
        forward_power_W=100.0,
        reflected_power_W=4.0,
        hardware_loss_lower_W=17.0,
        hardware_loss_upper_W=19.0,
        hardware_loss_point_W=18.0,
        measurement_source="calibrated directional coupler",
        loss_source="plasma-off measured match and coil dissipation",
        loss_evidence="measured",
    ).estimate()
    assert estimate.lower_W == 77.0
    assert estimate.upper_W == 79.0
    assert estimate.point_W == 78.0
    assert estimate.supports_prediction


def test_direct_drive_output_still_requires_downstream_hardware_loss():
    time = np.arange(4096, dtype=float) / 4096.0
    voltage = np.sqrt(2.0) * 100.0 * np.sin(2.0 * np.pi * time)
    current = np.sqrt(2.0) * 2.0 * np.sin(
        2.0 * np.pi * time - np.deg2rad(60.0))
    output_power = time_average_real_power_W(voltage, current)
    assert np.isclose(output_power, 100.0, rtol=1.0e-12)

    estimate = DirectDriveRFPowerBoundary(
        output_real_power_W=output_power,
        hardware_loss_lower_W=5.0,
        hardware_loss_upper_W=15.0,
        measurement_source="simultaneous output-node v(t), i(t)",
        loss_source="unresolved coil/window loss interval",
    ).estimate()
    assert np.isclose(estimate.lower_W, 85.0, rtol=1.0e-12)
    assert np.isclose(estimate.upper_W, 95.0, rtol=1.0e-12)
    assert not estimate.supports_prediction


@pytest.mark.parametrize(
    "boundary",
    [
        MatchedRFPowerBoundary(
            forward_power_W=100.0,
            reflected_power_W=100.0,
            hardware_loss_lower_W=0.0,
            hardware_loss_upper_W=1.0,
            measurement_source="bad reflection",
            loss_source="test",
        ),
        DirectDriveRFPowerBoundary(
            output_real_power_W=100.0,
            hardware_loss_lower_W=0.0,
            hardware_loss_upper_W=100.0,
            measurement_source="bad complete loss",
            loss_source="test",
        ),
    ],
)
def test_power_boundaries_fail_closed_on_nonphysical_inputs(boundary):
    with pytest.raises(ValueError):
        boundary.estimate()


def test_waveform_power_rejects_shape_and_sign_errors():
    with pytest.raises(ValueError, match="equal finite"):
        time_average_real_power_W([1.0, 2.0], [1.0])
    with pytest.raises(ValueError, match="positive"):
        time_average_real_power_W([1.0, -1.0], [-1.0, 1.0])
