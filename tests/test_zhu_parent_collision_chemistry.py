import math
import os
from pathlib import Path

import numpy as np
import pytest

from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)
from petch.reactor_global.zhu_parent_collision_chemistry import (
    build_zhu_parent_collision_chemistry,
    zhu_parent_collision_species,
)


LOCAL_SOURCE_CANDIDATES = (
    Path(os.environ.get("PETCH_SONG_2026_O2_WORKBOOK", "")),
    Path("/private/tmp/o2_song_2026_supplement.xlsx"),
)


def _local_source() -> Path:
    for path in LOCAL_SOURCE_CANDIDATES:
        if str(path) and path.is_file():
            return path
    pytest.skip("licensed Song 2026 O2 source workbook not supplied")


def _grid(deck) -> ElectronEnergyGrid:
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV and process.energy_loss_eV < 120.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, .0001, .001, .01, .1, 1.0, 10.0, 40.0, 120.0),
        (40, 60, 90, 48, 60, 120, 120, 150),
        inserted_boundaries_eV=thresholds,
    )


def test_parent_product_species_have_explicit_atoms_and_charge():
    species = zhu_parent_collision_species()
    by_name = {item.name: item for item in species}
    assert len(species) == len(by_name) == 44
    assert by_name["SF4++"].charge_number == 2
    assert by_name["SF4++"].composition == {"S": 1, "F": 4}
    assert by_name["F2-"].charge_number == -1
    assert by_name["CHF2+"].composition == {"C": 1, "H": 1, "F": 2}
    assert by_name["O2+"].composition == {"O": 2}


def test_all_parent_collision_rows_have_one_conserved_heavy_mapping():
    replay = build_zhu_parent_collision_chemistry(_local_source())
    assert replay.mixed_deck.targets == ("CHF3", "O2", "SF6")
    assert len(replay.mixed_deck.processes) == 49
    assert len(replay.collision_chemistry.mappings) == 46
    assert replay.sf6_replay.supports_direct_attachment_products
    assert replay.sf6_neutral_dissociation_closure == "dominant_SF5_plus_F"
    assert replay.o2_positive_ionization_closure == "all_O2plus"
    assert replay.supports_complete_daughter_eedf is False
    assert replay.supports_feature_depth is False


def test_solved_eepf_emits_atom_and_charge_closed_parent_sources():
    replay = build_zhu_parent_collision_chemistry(_local_source())
    feed = {"CHF3": 55.0 / 61.0, "SF6": 5.0 / 61.0, "O2": 1.0 / 61.0}
    gas_density = 3.99967104 / (1.380649e-23 * 293.15)
    densities = {item.name: 0.0 for item in replay.species}
    densities.update({name: fraction * gas_density for name, fraction in feed.items()})
    densities["e"] = 1.0e16
    condition = TwoTermBoltzmannCondition(
        reduced_electric_field_Td=225.0,
        gas_temperature_K=293.15,
        target_mole_fractions=feed,
        growth_model="temporal_growth",
        angular_field_frequency_over_density_m3_s=(
            2.0 * math.pi * 13.56e6 / gas_density
        ),
    )
    solution = DeterministicTwoTermBoltzmannSolver(
        _grid(replay.mixed_deck), replay.mixed_deck
    ).solve(
        condition,
        relative_tolerance=3.0e-6,
        maximum_tail_population_fraction=2.0e-6,
    )
    state = replay.collision_chemistry.evaluate(
        solution, condition, densities, closure_relative_tolerance=2.0e-7)
    by_name = {item.name: item for item in replay.species}
    for element in ("C", "H", "F", "S", "O"):
        residual = sum(
            state.species_sources_m3_s[name]
            * species.composition.get(element, 0)
            for name, species in by_name.items()
        )
        scale = sum(
            abs(state.species_sources_m3_s[name])
            * species.composition.get(element, 0)
            for name, species in by_name.items()
        )
        assert abs(residual) <= 2.0e-15 * max(scale, 1.0)
    charge_residual = sum(
        state.species_sources_m3_s[name] * species.charge_number
        for name, species in by_name.items()
    )
    charge_scale = sum(
        abs(state.species_sources_m3_s[name] * species.charge_number)
        for name, species in by_name.items()
    )
    assert abs(charge_residual) <= 2.0e-15 * max(charge_scale, 1.0)
    assert state.relative_electron_growth_closure < 2.0e-7
    assert any(
        "SF4++" in name and rate > 0.0
        for name, rate in state.event_rates_m3_s.items()
    )
