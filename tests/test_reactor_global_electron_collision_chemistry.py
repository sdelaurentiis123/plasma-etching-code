import pytest

from petch.reactor_global.chlorine import lee_lieberman_chlorine_species
from petch.reactor_global.electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionHeavyMapping,
)
from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)


def _processes():
    return (
        ElectronCollisionProcess(
            kind="ELASTIC",
            target="Cl2",
            product=None,
            electron_energy_eV=(0.0, 100.0),
            cross_section_m2=(2.0e-20, 2.0e-20),
            mass_ratio=7.68e-6,
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            target="Cl2",
            product="Cl2(v=1)",
            electron_energy_eV=(0.0, 2.0, 5.0, 100.0),
            cross_section_m2=(0.0, 0.0, 1.0e-20, 1.0e-20),
            energy_loss_eV=2.0,
        ),
        ElectronCollisionProcess(
            kind="ATTACHMENT",
            target="Cl2",
            product="Cl- + Cl",
            electron_energy_eV=(0.0, 0.2, 1.0, 100.0),
            cross_section_m2=(2.0e-20, 2.0e-20, 0.0, 0.0),
            energy_loss_eV=0.0,
        ),
        ElectronCollisionProcess(
            kind="IONIZATION",
            target="Cl2",
            product="Cl2+",
            electron_energy_eV=(0.0, 10.0, 15.0, 100.0),
            cross_section_m2=(0.0, 0.0, 2.0e-20, 2.0e-20),
            energy_loss_eV=10.0,
        ),
    )


def _deck():
    return ElectronCollisionDeck(
        processes=_processes(),
        payload_sha256="d" * 64,
        source_database="manufactured chemistry test",
        retrieved_at="2026-08-08",
        source_reference="test only",
    )


def _mappings():
    common = {"source": "manufactured test", "evidence_kind": "manufactured"}
    return (
        ElectronCollisionHeavyMapping(
            process_index=1,
            reaction_name="vibrational_excitation_untracked_return",
            heavy_reactants={"Cl2": 1},
            heavy_products={"Cl2": 1},
            **common,
        ),
        ElectronCollisionHeavyMapping(
            process_index=2,
            reaction_name="dissociative_attachment",
            heavy_reactants={"Cl2": 1},
            heavy_products={"Cl-": 1, "Cl": 1},
            **common,
        ),
        ElectronCollisionHeavyMapping(
            process_index=3,
            reaction_name="molecular_ionization",
            heavy_reactants={"Cl2": 1},
            heavy_products={"Cl2+": 1},
            **common,
        ),
    )


def test_collision_chemistry_closes_growth_atoms_charge_and_power():
    deck = _deck()
    chemistry = ElectronCollisionChemistry(
        deck, lee_lieberman_chlorine_species(), _mappings())
    condition = TwoTermBoltzmannCondition(
        reduced_electric_field_Td=100.0,
        gas_temperature_K=300.0,
        target_mole_fractions={"Cl2": 1.0},
        growth_model="temporal_growth",
    )
    solution = DeterministicTwoTermBoltzmannSolver(
        ElectronEnergyGrid.linear(80.0, 480), deck).solve(
            condition,
            relative_tolerance=1.0e-8,
            maximum_tail_population_fraction=1.0e-5,
        )
    densities = {
        "e": 2.0e16,
        "Cl2": 3.0e20,
        "Cl": 1.0e20,
        "Cl2+": 2.0e16,
        "Cl+": 1.0e15,
        "Cl-": 1.0e15,
    }
    state = chemistry.evaluate(solution, condition, densities)
    # The two ledgers deliberately recompute the same physical moment by
    # independent summation paths. Judge their closure relative to the
    # O(1e21--1e22) event-rate scale rather than with an absolute tolerance
    # smaller than a few floating-point ulps at that scale.
    assert state.relative_electron_growth_closure < 5.0e-15
    assert state.species_sources_m3_s["e"] > 0.0
    assert state.species_sources_m3_s["Cl2+"] > 0.0
    assert state.species_sources_m3_s["Cl-"] > 0.0
    chlorine_source = (
        2.0 * state.species_sources_m3_s["Cl2"]
        + state.species_sources_m3_s["Cl"]
        + 2.0 * state.species_sources_m3_s["Cl2+"]
        + state.species_sources_m3_s["Cl+"]
        + state.species_sources_m3_s["Cl-"]
    )
    charge_source = (
        -state.species_sources_m3_s["e"]
        + state.species_sources_m3_s["Cl2+"]
        + state.species_sources_m3_s["Cl+"]
        - state.species_sources_m3_s["Cl-"]
    )
    source_scale = max(abs(value) for value in state.species_sources_m3_s.values())
    assert abs(chlorine_source) / source_scale < 2.0e-15
    assert abs(charge_source) / source_scale < 2.0e-15
    assert state.collisional_field_power_gain_W_m3 > 0.0
    assert not state.supports_reactor_state_prediction
    assert not state.supports_feature_depth


def test_collision_chemistry_refuses_incomplete_or_nonconservative_mapping():
    deck = _deck()
    with pytest.raises(ValueError, match="every non-momentum"):
        ElectronCollisionChemistry(
            deck, lee_lieberman_chlorine_species(), _mappings()[:-1])
    broken = list(_mappings())
    broken[-1] = ElectronCollisionHeavyMapping(
        process_index=3,
        reaction_name="broken_ionization",
        heavy_reactants={"Cl2": 1},
        heavy_products={"Cl": 2},
        source="manufactured test",
        evidence_kind="manufactured",
    )
    with pytest.raises(ValueError, match="does not conserve charge"):
        ElectronCollisionChemistry(
            deck, lee_lieberman_chlorine_species(), tuple(broken))
