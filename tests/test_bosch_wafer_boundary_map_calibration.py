import math

import numpy as np
import pytest

from scripts.audit_bosch_wafer_boundary_map_calibration import (
    BoschExactWallIonResponseTable,
    _load_preregistration,
    _nodes,
)


def _synthetic_table():
    wall = _nodes(13)
    ion = _nodes(13)
    wafer = np.arange(3, dtype=float)[None, None, :, None]
    point = np.arange(89, dtype=float)[None, None, None, :]
    wall_value = wall[:, None, None, None]
    ion_value = ion[None, :, None, None]
    silicon = 1.0e-6 * (
        40.0 + 2.0 * wall_value + 3.0 * ion_value
        + 0.1 * wafer + 0.001 * point)
    oxide = 1.0e-6 * (
        0.6 + 0.02 * wall_value + 0.03 * ion_value
        + 0.001 * wafer + 0.00001 * point)
    return BoschExactWallIonResponseTable(
        experiment_keys=("a", "b", "c"),
        log_wall_multiplier_nodes=wall,
        log_ion_factor_nodes=ion,
        silicon_depth_m=silicon,
        oxide_loss_m=oxide,
    )


def test_v8_preregistration_hash_firewall_is_closed():
    payload = _load_preregistration()

    assert payload["target_firewall"]["heldout_outcomes_read_at_freeze"] is False
    assert payload["calibration_exposure"]["v8_fit_started_before_this_freeze"] is False


def test_wall_ion_response_tensor_interpolates_independent_local_queries():
    table = _synthetic_table()
    wall_query = np.exp(np.array([-0.7, 0.0, 0.8]))
    conditioned = table.condition_on_wall(wall_query)
    log_ion_query = np.linspace(-1.0, 1.0, 3 * 89).reshape(3, 89)
    predictions, silicon_derivative, oxide_derivative = conditioned.evaluate(
        np.exp(log_ion_query))
    silicon = np.stack([item.silicon_depth_m for item in predictions])
    oxide = np.stack([item.oxide_loss_m for item in predictions])
    wafer = np.arange(3, dtype=float)[:, None]
    point = np.arange(89, dtype=float)[None]
    expected_silicon = 1.0e-6 * (
        40.0 + 2.0 * np.log(wall_query)[:, None]
        + 3.0 * log_ion_query + 0.1 * wafer + 0.001 * point)
    expected_oxide = 1.0e-6 * (
        0.6 + 0.02 * np.log(wall_query)[:, None]
        + 0.03 * log_ion_query + 0.001 * wafer + 0.00001 * point)

    assert silicon == pytest.approx(expected_silicon, rel=2.0e-13, abs=1e-18)
    assert oxide == pytest.approx(expected_oxide, rel=2.0e-13, abs=1e-18)
    assert silicon_derivative == pytest.approx(
        np.full((3, 89), 3.0e-6), rel=2.0e-13, abs=1e-18)
    assert oxide_derivative == pytest.approx(
        np.full((3, 89), 0.03e-6), rel=2.0e-13, abs=1e-18)


def test_wall_ion_response_rejects_extrapolation():
    table = _synthetic_table()
    with pytest.raises(ValueError, match="wall-multiplier"):
        table.condition_on_wall(np.array([0.2, 1.0, 1.0]))
    conditioned = table.condition_on_wall(np.ones(3))
    with pytest.raises(ValueError, match="ion-factor"):
        conditioned.evaluate(np.full((3, 89), math.exp(1.5)))
