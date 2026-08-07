import numpy as np
import pytest

from petch.guo_c4f8_sio2 import (
    GuoC4F8ArSiO2Mechanism,
    GuoIncidentComposition,
    GuoIonQuadrature,
    GuoSourceLawUnderspecified,
    GuoTmlState,
    formula_atoms,
    ion_enhanced_angular,
    nearest_neighbor_probability,
    physical_sputtering_angular,
    physical_sputtering_angular_literal,
)


def test_formula_parser_and_incident_boundary_refuse_hidden_species_mass():
    assert formula_atoms("C3F4+") == {"C": 3.0, "F": 4.0}
    assert formula_atoms("SiF3") == {"Si": 1.0, "F": 3.0}
    assert formula_atoms("Ar") == {}
    with pytest.raises(ValueError, match="exceed"):
        GuoIncidentComposition({}, {"CF+": 0.7, "Ar+": 0.4})
    with pytest.raises(ValueError, match="unsupported"):
        formula_atoms("C4F6*")


def test_guo_nearest_neighbor_equation_is_symmetric_and_counts_bonds():
    state = GuoTmlState(0.30, 0.40, 0.15, 0.15, 0.08)
    species = ("Si", "O", "C", "F", "V")

    for first in species:
        for second in species:
            assert nearest_neighbor_probability(
                state, first, second
            ) == pytest.approx(
                nearest_neighbor_probability(state, second, first))

    total_unordered_bonds = sum(
        nearest_neighbor_probability(state, first, second)
        for index, first in enumerate(species)
        for second in species[index:]
    )
    valence_density = (
        4.0 * state.si + 2.0 * state.o + 4.0 * state.c
        + state.f + state.vacancy
    )
    assert total_unordered_bonds == pytest.approx(
        0.5 * valence_density)


def test_guo_angular_laws_reproduce_source_limits_and_repair_table_typo():
    cosine = np.cos(np.deg2rad([0.0, 25.0, 65.0, 85.0, 90.0]))
    physical = physical_sputtering_angular(cosine)
    enhanced = ion_enhanced_angular(cosine)

    assert physical[0] == pytest.approx(1.02, abs=2.0e-12)
    assert physical[2] > physical[0]
    assert physical[-1] == 0.0
    assert np.all(physical >= 0.0)
    assert enhanced[0] == 1.0
    assert enhanced[1] == pytest.approx(1.0)
    assert np.all(enhanced >= 0.0)

    literal = physical_sputtering_angular_literal(cosine)
    assert literal[0] == pytest.approx(1.02, abs=2.0e-12)
    assert np.all(literal[1:] < 0.0)


def test_low_energy_generic_ion_incorporation_refuses_missing_threshold():
    with pytest.raises(GuoSourceLawUnderspecified, match="72 eV"):
        GuoC4F8ArSiO2Mechanism(
            GuoIncidentComposition({"CF2": 10.0}, {}),
            GuoIonQuadrature.monoenergetic(60.0),
        )


def test_translating_layer_reaches_atom_balanced_steady_state_without_depth_fit():
    composition = GuoIncidentComposition(
        {
            "C3F4": 9.5e16 / 1.2e16,
            "C2F3": 6.8e16 / 1.2e16,
            "CF": 4.4e16 / 1.2e16,
            "CF2": 9.4e16 / 1.2e16,
            "CF3": 8.4e15 / 1.2e16,
            "O": 7.7e16 / 1.2e16,
        },
        {},
    )
    mechanism = GuoC4F8ArSiO2Mechanism(
        composition,
        GuoIonQuadrature.monoenergetic(350.0),
    )
    result = mechanism.solve_steady_state(maximum_coordinate=2.0e4)

    assert sum(result.state.as_array()[:4]) == pytest.approx(1.0)
    assert result.state.vacancy >= 0.0
    assert abs(result.elemental_derivative_residual) < 1.0e-12
    assert result.steady_state_residual < 2.0e-8
    assert result.sio2_yield_per_ion > 0.0
    assert not result.source_extrapolation["beyond_source_fit_energy"]
    assert result.movement_atoms_per_ion == pytest.approx(
        sum(result.removed_atoms_per_ion.values())
        - sum(result.incoming_atoms_per_ion.values())
    )
