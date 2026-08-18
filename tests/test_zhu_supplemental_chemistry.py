import numpy as np
import pytest

from petch.reactor_global import RateContext
from petch.reactor_global.zhu_supplemental_chemistry import (
    KOKKORIS_2009_TARGET_PRESSURE_PA,
    ZHU_RECIPE_PRESSURE_PA,
    build_zhu_supplemental_chemistry,
    zhu_reactor_species,
)


def test_supplemental_network_is_closed_and_does_not_duplicate_parent_rows():
    chemistry = build_zhu_supplemental_chemistry()
    names = {reaction.name for reaction in chemistry.network.reactions}
    assert len(chemistry.network.species) == 66
    assert len(chemistry.network.reactions) == 259
    assert chemistry.sandia_rows_selected == tuple(range(20, 39))
    assert chemistry.chf3_f_rate_branch == "voloshin_350K"
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
    assert chemistry.supports_oxygen_daughter_chemistry is True
    assert chemistry.supports_sf6_o2_titration_chemistry is True
    assert chemistry.supports_chf3_neutral_chain is True
    assert chemistry.supports_complete_cross_ion_recombination is False
    with pytest.raises(ValueError):
        build_zhu_supplemental_chemistry(kokkoris_eedf_shape="interpolated")


def test_conflicting_chf3_f_rates_are_explicit_branches_not_duplicates():
    voloshin = build_zhu_supplemental_chemistry(
        chf3_f_rate_branch="voloshin_350K")
    lim = build_zhu_supplemental_chemistry(chf3_f_rate_branch="lim_700K")
    v_reaction = next(
        reaction for reaction in voloshin.network.reactions
        if reaction.name == "voloshin_2007_R13")
    l_reaction = next(
        reaction for reaction in lim.network.reactions
        if reaction.name == "lim_2014_R26")
    context = RateContext(3.0, 350.0)
    assert (
        v_reaction.rate_coefficient.coefficient_si(context)
        / l_reaction.rate_coefficient.coefficient_si(context)
    ) == pytest.approx(1.82e-12 / 1.58e-13)
    with pytest.raises(ValueError):
        build_zhu_supplemental_chemistry(chf3_f_rate_branch="average")


def test_every_parent_negative_ion_has_a_volume_loss_path():
    chemistry = build_zhu_supplemental_chemistry()
    negative_names = {
        species.name for species in zhu_reactor_species()
        if species.role == "negative_ion"
    }
    consumed = {
        name
        for reaction in chemistry.network.reactions
        for name in reaction.reactants
        if name in negative_names
    }
    assert consumed == negative_names


def test_oxygen_titration_releases_fluorine_without_losing_atoms():
    chemistry = build_zhu_supplemental_chemistry()
    reaction = next(
        item for item in chemistry.network.reactions
        if item.name == "pateau_2014_R134_O")
    assert reaction.reactants == {"SF3": 1, "O": 1}
    assert reaction.products == {"SOF3": 1}
    fluorine_releasing = next(
        item for item in chemistry.network.reactions
        if item.name == "pateau_2014_R119_O")
    assert fluorine_releasing.products == {"SOF2": 1, "F": 1}
