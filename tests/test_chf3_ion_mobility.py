import math

import numpy as np
import pytest

from petch.reactor_global.chf3_ion_mobility import (
    load_basurto_2002_chf2_chf3_mobility_model,
)
from scripts.digitize_basurto_2002_chf2_chf3_mobility import (
    CSV_PATH,
    MANIFEST_PATH,
    csv_text,
    manifest_text,
)


def test_committed_basurto_digitization_replays_from_pixels():
    payload = csv_text()
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)


def test_measured_nodes_replay_and_interpolant_is_positive_c1():
    model = load_basurto_2002_chf2_chf3_mobility_model()
    np.testing.assert_allclose(
        model.reduced_mobility(model.reduced_field_Td),
        model.reduced_mobility_cm2_V_s,
        rtol=0.0,
        atol=2.0e-15,
    )
    probe = np.geomspace(*model.reduced_field_support_Td, 1001)
    assert np.all(model.reduced_mobility(probe) > 0.0)
    assert np.all(np.isfinite(
        model.reduced_mobility_derivative_cm2_V_s_per_Td(probe)))


def test_model_fails_closed_outside_digitized_support():
    model = load_basurto_2002_chf2_chf3_mobility_model()
    with pytest.raises(ValueError, match="outside digitized support"):
        model.reduced_mobility(40.0)
    with pytest.raises(ValueError, match="outside digitized support"):
        model.reduced_mobility(500.0)


def test_target_pressure_scale_is_collisional_and_not_a_fake_iead():
    model = load_basurto_2002_chf2_chf3_mobility_model()
    density = 30.0 * 0.133322368 / (1.380649e-23 * 293.15)
    low = model.evaluate(
        reduced_field_Td=100.0,
        total_neutral_density_m3=density,
    )
    high = model.evaluate(
        reduced_field_Td=400.0,
        total_neutral_density_m3=density,
    )

    assert 0.48 < low.reduced_mobility_cm2_V_s < 0.50
    assert 125.0 < low.drift_speed_m_s < 135.0
    assert 0.08e-3 < low.drift_relaxation_length_m < 0.13e-3
    assert 0.7e-3 < high.drift_relaxation_length_m < 1.5e-3
    assert high.drift_relaxation_length_m > low.drift_relaxation_length_m
    assert low.supports_measured_swarm_transport is True
    assert low.supports_elastic_differential_cross_section is False
    assert low.supports_target_sheath_iead is False
    assert low.supports_absolute_depth_prediction is False
    assert model.provenance["coefficient_selected_from_depth_target"] is None


def test_reduced_mobility_conversion_has_correct_density_scaling():
    model = load_basurto_2002_chf2_chf3_mobility_model()
    state_a = model.evaluate(
        reduced_field_Td=200.0,
        total_neutral_density_m3=1.0e21,
    )
    state_b = model.evaluate(
        reduced_field_Td=200.0,
        total_neutral_density_m3=2.0e21,
    )
    assert math.isclose(
        state_a.actual_mobility_m2_V_s,
        2.0 * state_b.actual_mobility_m2_V_s,
    )
    assert math.isclose(state_a.drift_speed_m_s, state_b.drift_speed_m_s)
