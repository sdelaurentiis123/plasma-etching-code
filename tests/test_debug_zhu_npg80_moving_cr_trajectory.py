import pytest

from scripts import audit_zhu_npg80_moving_cr_profiles as board
from scripts.debug_zhu_npg80_moving_cr_trajectory import _production_job


def test_debug_driver_selects_exact_v7_production_cell():
    job = _production_job(
        width_nm=320.0,
        scenario_name="ion_low_tail_0p0",
        selectivity=14.0,
    )
    spec = board._job_spec(job)

    assert spec["width_nm"] == 320.0
    assert spec["scenario"]["name"] == "ion_low_tail_0p0"
    assert spec["selectivity"] == 14.0
    assert spec["duration_s"] == 1200.0
    assert spec["mesh_spacing_nm"] == 10.0
    assert board._cache_path(spec).name == (
        "w320_s14.000_ion_low_tail_0p0_c9b10e259e2c5b39.json"
    )


def test_debug_driver_refuses_nonboard_parameter_choices():
    with pytest.raises(ValueError, match="not in board"):
        _production_job(
            width_nm=319.0,
            scenario_name="ion_low_tail_0p0",
            selectivity=14.0,
        )
    with pytest.raises(ValueError, match="not in"):
        _production_job(
            width_nm=320.0,
            scenario_name="not-a-scenario",
            selectivity=14.0,
        )
