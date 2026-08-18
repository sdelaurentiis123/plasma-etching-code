import math
from pathlib import Path
import subprocess
import sys

from petch.reactor_global.c4f6_electron_collisions import (
    C4F6_MASS_AMU,
    load_lan_jeon_2014_c4f6_replay,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)


ROOT = Path(__file__).resolve().parents[1]


def _grid(deck):
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV and process.energy_loss_eV < 200.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, .001, .01, .1, 1.0, 10.0, 50.0, 100.0, 200.0),
        (8, 12, 24, 48, 96, 96, 80, 80),
        inserted_boundaries_eV=thresholds,
    )


def test_committed_lan_jeon_tables_and_figure7_board_replay():
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "extract_c4f6_electron_evidence.py"),
        "--check",
    ], check=True, cwd=ROOT)
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "digitize_lan_jeon_2014_figure7.py"),
        "--check",
    ], check=True, cwd=ROOT)


def test_c4f6_working_set_topology_units_and_evidence_boundary_are_exact():
    replay = load_lan_jeon_2014_c4f6_replay()
    assert replay.process_labels == (
        "Qm", "Qa", "Qv1", "Qv2", "Qex1", "Qex2", "Qex3", "Qex4",
        "Qdiss", "Qi",
    )
    assert len(replay.raw_deck.processes) == 10
    momentum = replay.raw_deck.processes[0]
    assert (momentum.electron_energy_eV[0], momentum.electron_energy_eV[-1]) == (
        0.0001, 362.0,
    )
    assert (momentum.cross_section_m2[0], momentum.cross_section_m2[-1]) == (
        83.9e-20, 1.1e-20,
    )
    assert momentum.mass_ratio == 5.48579909065e-4 / C4F6_MASS_AMU
    assert replay.derived_deck.processes[0].electron_energy_eV[0] == 0.0
    assert replay.derived_deck.processes[-1].electron_energy_eV[-1] == 200.0
    assert replay.raw_deck.processes[1].electron_number_change == -1
    assert replay.raw_deck.processes[-1].electron_number_change == 1
    assert replay.supports_resolved_primary_chemistry is False
    assert replay.supports_reactor_state_prediction is False
    assert replay.supports_feature_depth is False


def test_c4f6_local_flux_solution_is_deterministic_but_not_pt_drift_grade():
    replay = load_lan_jeon_2014_c4f6_replay()
    solver = DeterministicTwoTermBoltzmannSolver(
        _grid(replay.derived_deck), replay.derived_deck)
    solution = solver.solve(TwoTermBoltzmannCondition(
        reduced_electric_field_Td=298.048,
        gas_temperature_K=300.0,
        target_mole_fractions={"C4F6": 1.0},
        growth_model="temporal_growth",
        initial_electron_temperature_eV=.2,
    ), relative_tolerance=2.0e-6, maximum_tail_population_fraction=1.0e-6)
    drift = (
        solution.transport_moments
        .flux_reduced_mobility_m_inv_V_inv_s_inv
        * 298.048e-21
    )
    assert math.isclose(drift, 185406.47811890102, rel_tol=2.0e-11)
    assert math.isclose(
        solution.distribution.mean_energy_eV,
        3.2529321426768774,
        rel_tol=2.0e-11,
    )
    assert solution.supports_direct_swarm_grade is False


def test_committed_c4f6_swarm_audit_fails_closed_on_depth_authority():
    path = ROOT / "results" / "curated" / "c4f6_electron_swarm_v1" / "audit.json"
    import json
    audit = json.loads(path.read_text(encoding="utf-8"))
    transport = audit["transport_definition_diagnostic"]
    convergence = audit["numerical_convergence"]
    certification = audit["certification"]
    assert transport["measurement_equivalent_grade"] is False
    assert .07 < transport["mean_absolute_relative_residual"] < .09
    assert convergence["maximum_absolute_flux_drift_relative_change"] < .001
    assert certification["supports_use_as_c4f6_component_collision_input"] is True
    assert certification["supports_unique_krueger_reactor_state"] is False
    assert certification["supports_feature_depth"] is False
