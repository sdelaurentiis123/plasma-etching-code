import pytest

from petch.reactor_global import RateContext
from petch.reactor_global.lim_2014_chf3_oxygen_chemistry import (
    build_lim_2014_daughter_chemistry,
)
from petch.reactor_global.network import ReactionNetwork
from petch.reactor_global.zhu_supplemental_chemistry import zhu_reactor_species


def test_landed_subset_has_every_declared_table_row_once():
    chemistry = build_lim_2014_daughter_chemistry()
    assert len(chemistry.species) == 7
    assert len(chemistry.reactions) == 62
    assert chemistry.electron_rows_selected == (
        6, 10, 11, 12, 13, 14, 17, 20, 21, 23, 24, 25)
    assert chemistry.neutral_rows_selected == tuple(range(26, 76))
    names = {reaction.name for reaction in chemistry.reactions}
    assert "voloshin_2007_R13" in names
    assert "lim_2014_R26" not in names
    assert chemistry.supports_target_temperature_transfer is False


def test_full_zhu_species_basis_conserves_every_landed_reaction():
    chemistry = build_lim_2014_daughter_chemistry()
    network = ReactionNetwork(
        species=zhu_reactor_species(), reactions=chemistry.reactions)
    network.assert_closed_conservation()


def test_conflicting_temperature_branches_keep_printed_coefficients():
    context = RateContext(electron_temperature_eV=3.0, gas_temperature_K=350.0)
    voloshin = build_lim_2014_daughter_chemistry(
        chf3_f_rate_branch="voloshin_350K")
    lim = build_lim_2014_daughter_chemistry(chf3_f_rate_branch="lim_700K")
    v = voloshin.reactions[12].rate_coefficient.coefficient_si(context)
    l = lim.reactions[12].rate_coefficient.coefficient_si(context)
    assert v == pytest.approx(1.82e-18)
    assert l == pytest.approx(1.58e-19)
    assert v / l == pytest.approx(1.82e-12 / 1.58e-13)


def test_representative_neutral_and_electron_rows_match_source():
    chemistry = build_lim_2014_daughter_chemistry()
    by_name = {reaction.name: reaction for reaction in chemistry.reactions}
    assert by_name["lim_2014_R43"].reactants == {"CF3": 1, "H": 1}
    assert by_name["lim_2014_R43"].products == {"CF2": 1, "HF": 1}
    assert by_name["lim_2014_R47"].rate_coefficient.coefficient_si(
        RateContext(3.0, 350.0)) == pytest.approx(3.20e-17)
    assert by_name["lim_2014_R10"].electron_energy_loss_eV == 3.80
