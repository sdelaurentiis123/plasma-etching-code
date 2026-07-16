import pytest

from petch.experimental_data import jeon_2022_condition_wall_duration_s


def test_jeon_pulse_exposure_basis_is_explicit_and_changes_only_wall_duration():
    assert jeon_2022_condition_wall_duration_s(600.0, 1.0, "unspecified") == 600.0
    assert jeon_2022_condition_wall_duration_s(600.0, 0.5, "wall_time") == 600.0
    assert jeon_2022_condition_wall_duration_s(600.0, 0.5, "rf_on_time") == 1200.0
    with pytest.raises(ValueError, match="pulse-exposure-basis"):
        jeon_2022_condition_wall_duration_s(600.0, 0.5, "unspecified")


@pytest.mark.parametrize(
    "duration,duty,basis",
    [(0.0, 0.5, "wall_time"), (1.0, 0.0, "wall_time"),
     (1.0, 1.1, "wall_time"), (1.0, 0.5, "guessed")])
def test_jeon_pulse_exposure_contract_rejects_invalid_inputs(duration, duty, basis):
    with pytest.raises(ValueError):
        jeon_2022_condition_wall_duration_s(duration, duty, basis)
