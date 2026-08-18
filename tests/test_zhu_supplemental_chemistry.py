import numpy as np
import pytest

from petch.reactor_global import RateContext
from petch.reactor_global.zhu_supplemental_chemistry import (
    KOKKORIS_2009_TARGET_PRESSURE_PA,
    ZHU_RECIPE_PRESSURE_PA,
    build_zhu_supplemental_chemistry,
)


def test_supplemental_network_is_closed_and_does_not_duplicate_parent_rows():
    chemistry = build_zhu_supplemental_chemistry()
    names = {reaction.name for reaction in chemistry.network.reactions}
    assert len(chemistry.network.species) == 45
    assert len(chemistry.network.reactions) == 46
    assert chemistry.parent_sf6_rows_replaced == (
        "G1", "G2", "G3", "G8", "G9", "G10", "G17", "G18")
    assert not any(
        name == f"kokkoris_2009_{row}"
        or name.startswith(f"kokkoris_2009_{row}_")
        for row in chemistry.parent_sf6_rows_replaced
        for name in names
    )
    assert chemistry.pressure_specific_rows_excluded == ("G35", "G36", "G37")
    assert not any(
        f"kokkoris_2009_{row}" in name
        for row in chemistry.pressure_specific_rows_excluded
        for name in names
    )
    chemistry.network.assert_closed_conservation()


def test_two_kokkoris_eedf_assumptions_are_exposed_not_averaged():
    druyvesteyn = build_zhu_supplemental_chemistry(
        kokkoris_eedf_shape="druyvesteyn")
    maxwellian = build_zhu_supplemental_chemistry(
        kokkoris_eedf_shape="maxwellian")
    d_rates = {
        reaction.name.rsplit("_", 1)[0]:
        reaction.rate_coefficient.coefficient_si(RateContext(3.0))
        for reaction in druyvesteyn.network.reactions
        if reaction.name.startswith("kokkoris_2009_G")
        and reaction.name.endswith("_druyvesteyn")
    }
    m_rates = {
        reaction.name.rsplit("_", 1)[0]:
        reaction.rate_coefficient.coefficient_si(RateContext(3.0))
        for reaction in maxwellian.network.reactions
        if reaction.name.startswith("kokkoris_2009_G")
        and reaction.name.endswith("_maxwellian")
    }
    assert d_rates.keys() == m_rates.keys()
    ratios = np.asarray([d_rates[name] / m_rates[name] for name in d_rates])
    assert ratios.min() < 0.5
    assert ratios.max() > 1.05


def test_pressure_specific_neutral_falloff_is_not_silently_transferred():
    chemistry = build_zhu_supplemental_chemistry()
    assert ZHU_RECIPE_PRESSURE_PA == pytest.approx(3.999671052631579)
    assert ZHU_RECIPE_PRESSURE_PA / KOKKORIS_2009_TARGET_PRESSURE_PA == pytest.approx(1.9998355263157894)
    assert chemistry.supports_target_pressure_falloff is False


def test_declared_limits_remain_false_until_next_physics_layers_land():
    chemistry = build_zhu_supplemental_chemistry()
    assert chemistry.supports_measured_parent_eedf is True
    assert chemistry.supports_complete_daughter_eedf is False
    assert chemistry.supports_complete_oxygen_heavy_chemistry is False
    with pytest.raises(ValueError):
        build_zhu_supplemental_chemistry(kokkoris_eedf_shape="interpolated")
