import math

import numpy as np
import pytest

from scripts.audit_bosch_dynamic_wall_calibration import (
    BoschExactWallResponseTable,
    _dynamic_steps,
    _inputs,
    _wall_nodes,
)


def test_exact_response_table_preserves_nodes_and_refuses_extrapolation():
    keys = ("synthetic_01", "synthetic_02")
    nodes = _wall_nodes(9)
    multiplier = np.exp(nodes)
    point_scale = np.linspace(1.0, 2.0, 89)
    silicon = np.stack([
        np.stack([
            value * point_scale,
            2.0 * value * point_scale,
        ])
        for value in multiplier
    ])
    oxide = silicon / 20.0
    table = BoschExactWallResponseTable(
        experiment_keys=keys,
        log_wall_multiplier_nodes=nodes,
        silicon_depth_m=silicon,
        oxide_loss_m=oxide,
    )

    predictions = table.interpolate(np.array([1.0, 4.0]))
    assert np.allclose(predictions[0].silicon_depth_m, point_scale, rtol=1e-14)
    assert np.allclose(
        predictions[1].silicon_depth_m, 8.0 * point_scale, rtol=1e-14)
    assert np.allclose(
        predictions[1].oxide_loss_m, 0.4 * point_scale, rtol=1e-14)
    with pytest.raises(ValueError, match="invalid Bosch wall multipliers"):
        table.interpolate(np.array([0.24, 1.0]))


def test_dynamic_history_carries_unmeasured_processed_wafer_without_using_outcome():
    measurements, process_traces, lot_type_by_date = _inputs()
    coefficients = np.array([
        0.0, 0.0, 0.0, math.log(0.2), math.log(0.05), 0.2,
    ])
    steps, _static, _dynamic = _dynamic_steps(
        coefficients, process_traces, lot_type_by_date)
    measurement_keys = {item.experiment_key for item in measurements}

    assert len(process_traces) == 76
    assert len(measurements) == 75
    assert set(steps) - measurement_keys == {"2024-07-02_07"}
    assert steps["2024-07-02_08"].start_state == (
        steps["2024-07-02_07"].end_state)
    assert steps["2024-07-05_01"].start_state.occupancy == 0.0
    assert steps["2024-07-02_10"].end_state.occupancy > (
        steps["2024-07-02_01"].end_state.occupancy)
