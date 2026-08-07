from pathlib import Path

import numpy as np
import pytest

from petch.interaction_data import load_kounis_melas_2024_tables
from petch.si_cl_ale_depth import (
    DEEPMD_CELL_LENGTH_ANGSTROM,
    DEEPMD_SI_ATOMS_PER_MATERIAL_ML,
    DIAMOND_SI_LATTICE_CONSTANT_ANGSTROM,
    VellaHaoAleBoundary,
    deepmd_cell_area_cm2,
    predict_vella_hao_ale_depth,
    printed_rom_chlorine_creation_per_ar,
    silicon_atomic_density_cm3,
    steady_cycle_removal_material_ml,
)
from petch.surface_interaction_table import SurfaceInteractionDomainError


DATA = (
    Path(__file__).parents[1]
    / "data"
    / "surface_interactions"
    / "kounis_melas_2024"
)


def test_material_layer_conversion_uses_atoms_area_and_diamond_si_density():
    density = silicon_atomic_density_cm3()
    area = deepmd_cell_area_cm2()
    material_layer_nm = (
        DEEPMD_SI_ATOMS_PER_MATERIAL_ML / area / density * 1.0e7)

    assert DEEPMD_CELL_LENGTH_ANGSTROM == 6 * DIAMOND_SI_LATTICE_CONSTANT_ANGSTROM
    assert np.isclose(material_layer_nm, 0.13575, rtol=0.0, atol=1.0e-15)
    assert np.isclose(
        material_layer_nm,
        DIAMOND_SI_LATTICE_CONSTANT_ANGSTROM / 4.0 / 10.0,
        rtol=0.0,
        atol=1.0e-15,
    )


def test_steady_cycle_removal_is_the_last_three_atom_counted_increments():
    table = load_kounis_melas_2024_tables(DATA).ale_cycles

    assert np.isclose(
        steady_cycle_removal_material_ml(table, 60.0),
        0.5648148148148148,
    )
    assert np.isclose(
        steady_cycle_removal_material_ml(table, 80.0),
        0.7824074074074074,
    )
    assert np.isclose(
        steady_cycle_removal_material_ml(table, 100.0),
        0.9583333333333333,
    )


def test_no_fit_absolute_depth_board_closes_dimensional_si_atom_ledger():
    tables = load_kounis_melas_2024_tables(DATA)
    results = [
        predict_vella_hao_ale_depth(
            energy, tables.ale_cycles, tables.sputtering)
        for energy in (60.0, 80.0, 100.0)
    ]

    assert np.allclose(
        [item.total_depth_nm for item in results],
        [0.48421003142705055, 0.945257376794254, 1.4006484721614578],
        rtol=0.0,
        atol=1.0e-14,
    )
    assert all(
        abs(item.dimensional_atom_balance_residual_cm2) < 1.0
        for item in results
    )
    assert all(
        item.experimental_positive_ion_fluence_cm2
        > item.simulated_ar_fluence_cm2
        for item in results
    )
    assert np.all(np.asarray([
        item.source_reported_tail_depth_uncertainty_nm
        for item in results
    ]) > 0.0)


def test_absolute_depth_transfer_refuses_unreleased_energy_and_negative_tail():
    tables = load_kounis_melas_2024_tables(DATA)

    with pytest.raises(SurfaceInteractionDomainError, match="ion_energy"):
        predict_vella_hao_ale_depth(
            40.0, tables.ale_cycles, tables.sputtering)
    with pytest.raises(SurfaceInteractionDomainError, match="ion_energy"):
        predict_vella_hao_ale_depth(
            110.0, tables.ale_cycles, tables.sputtering)
    with pytest.raises(ValueError, match="smaller than"):
        predict_vella_hao_ale_depth(
            80.0,
            tables.ale_cycles,
            tables.sputtering,
            VellaHaoAleBoundary(
                positive_ion_flux_cm2_s=1.0e15,
                ion_bombardment_duration_s=1.0,
            ),
        )


def test_printed_transient_rom_sicl2_law_creates_chlorine_at_partial_coverage():
    yield_sicl2 = 0.02

    assert printed_rom_chlorine_creation_per_ar(0.0, 0.0, yield_sicl2) == 0.0
    assert printed_rom_chlorine_creation_per_ar(1.0, 1.0, yield_sicl2) == 0.0
    assert np.isclose(
        printed_rom_chlorine_creation_per_ar(0.5, 0.5, yield_sicl2),
        yield_sicl2,
    )
    with pytest.raises(ValueError, match="coverages"):
        printed_rom_chlorine_creation_per_ar(1.1, 0.0, yield_sicl2)
