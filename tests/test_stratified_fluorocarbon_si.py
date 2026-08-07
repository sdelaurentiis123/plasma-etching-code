import numpy as np
import pytest

from petch.stratified_fluorocarbon_si import (
    StratifiedSiEvents,
    StratifiedSiParameters,
    StratifiedSiState,
    advance_stratified_si,
    stratified_ion_energies,
)


def _seeded_state():
    return StratifiedSiState(
        film_c_atoms_m2=4.0e19,
        film_f_atoms_m2=6.0e19,
        film_cf_bonds_m2=5.0e19,
        film_cc_crosslinks_m2=2.0e19,
        transport_c_atoms_m2=5.0e19,
        transport_si_atoms_m2=4.0e19,
        transport_si_f_bonds_m2=8.0e19,
        reaction_front_f_bonds_m2=1.0e20,
        cumulative_removed_si_atoms_m2=1.0e19,
        cumulative_drawn_si_atoms_m2=5.0e19,
    )


def test_atom_and_bond_ledgers_close_for_all_internal_transfers():
    parameters = StratifiedSiParameters()
    result = advance_stratified_si(
        _seeded_state(),
        StratifiedSiEvents(
            deposited_film_c_atoms_m2=4.0e18,
            deposited_film_f_atoms_m2=8.0e18,
            formed_film_cf_bonds_m2=8.0e18,
            formed_film_cc_crosslinks_m2=2.0e18,
            adsorbed_atomic_f_to_front_m2=3.0e18,
            broken_film_cf_bonds_m2=5.0e18,
            broken_film_cc_crosslinks_m2=1.0e18,
            transferred_film_f_to_front_m2=4.0e18,
            transferred_film_c_to_transport_m2=2.0e18,
            promoted_bulk_si_to_transport_m2=3.0e18,
            promoted_front_f_bonds_to_transport_m2=6.0e18,
            recycled_transport_f_bonds_to_front_m2=2.0e18,
            removed_film_c_atoms_m2=1.0e18,
            removed_film_f_atoms_m2=2.0e18,
            exported_film_cf_bonds_m2=1.0e18,
            removed_transport_c_atoms_m2=1.0e18,
            removed_transport_si_atoms_m2=2.0e18,
            exported_transport_si_f_bonds_m2=4.0e18,
        ),
        parameters,
    )
    assert result.maximum_absolute_ledger_residual < 1.0e4
    assert np.allclose(
        result.state.cumulative_drawn_si_atoms_m2,
        result.state.transport_si_atoms_m2
        + result.state.cumulative_removed_si_atoms_m2,
    )


def test_newly_promoted_silicon_cannot_be_removed_in_the_same_step():
    with pytest.raises(
            ValueError, match="start-of-step inventory"):
        advance_stratified_si(
            StratifiedSiState.bare(),
            StratifiedSiEvents(
                promoted_bulk_si_to_transport_m2=1.0e18,
                removed_transport_si_atoms_m2=1.0e18,
            ),
        )


def test_f_transfer_requires_cf_scission_and_forms_the_sif_ledger():
    state = StratifiedSiState(
        film_c_atoms_m2=1.0e19,
        film_f_atoms_m2=1.0e19,
        film_cf_bonds_m2=1.0e19,
    )
    with pytest.raises(ValueError, match="requires C-F scission"):
        advance_stratified_si(
            state,
            StratifiedSiEvents(
                transferred_film_f_to_front_m2=1.0e18),
        )
    result = advance_stratified_si(
        state,
        StratifiedSiEvents(
            broken_film_cf_bonds_m2=1.0e18,
            transferred_film_f_to_front_m2=1.0e18,
        ),
    )
    assert result.state.film_f_atoms_m2 == pytest.approx(9.0e18)
    assert result.state.reaction_front_f_bonds_m2 == pytest.approx(1.0e18)
    assert result.bond_ledger_residual_bonds_m2["Si-F"] == 0.0


def test_transport_capacity_is_depth_density_derived_and_exceeds_one_ml():
    parameters = StratifiedSiParameters()
    one_si_100_monolayer_m2 = 6.78e18
    assert (
        parameters.transport_capacity_atoms_m2
        > 30.0 * one_si_100_monolayer_m2
    )
    state = StratifiedSiState(
        transport_c_atoms_m2=0.6
        * parameters.transport_capacity_atoms_m2,
        transport_si_atoms_m2=0.4
        * parameters.transport_capacity_atoms_m2,
        cumulative_drawn_si_atoms_m2=0.4
        * parameters.transport_capacity_atoms_m2,
    )
    state.validate(parameters)
    with pytest.raises(ValueError, match="capacity"):
        StratifiedSiState(
            transport_c_atoms_m2=0.7
            * parameters.transport_capacity_atoms_m2,
            transport_si_atoms_m2=0.4
            * parameters.transport_capacity_atoms_m2,
            cumulative_drawn_si_atoms_m2=0.4
            * parameters.transport_capacity_atoms_m2,
        ).validate(parameters)


def test_csda_transmission_has_finite_range_and_slant_path():
    parameters = StratifiedSiParameters()
    bare = stratified_ion_energies(
        StratifiedSiState.bare(), 200.0, 1.0, parameters)
    assert bare.at_reaction_front_energy_eV == pytest.approx(200.0)

    film = StratifiedSiState(
        film_c_atoms_m2=0.5
        * parameters.film_atom_density_m3 * 1.0e-9,
        film_f_atoms_m2=0.5
        * parameters.film_atom_density_m3 * 1.0e-9,
    )
    normal = stratified_ion_energies(film, 200.0, 1.0, parameters)
    slanted = stratified_ion_energies(film, 200.0, 0.5, parameters)
    assert 0.0 < normal.after_film_energy_eV < 200.0
    assert slanted.after_film_energy_eV == 0.0
    assert slanted.at_reaction_front_energy_eV == 0.0


def test_state_rejects_impossible_valence_and_dead_si_ledger():
    with pytest.raises(ValueError, match="C-F bond"):
        StratifiedSiState(
            film_c_atoms_m2=1.0,
            film_f_atoms_m2=10.0,
            film_cf_bonds_m2=5.0,
        ).validate(StratifiedSiParameters())
    with pytest.raises(ValueError, match="cumulative silicon"):
        StratifiedSiState(
            transport_si_atoms_m2=1.0,
            cumulative_drawn_si_atoms_m2=0.0,
        ).validate(StratifiedSiParameters())
