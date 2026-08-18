import subprocess
import sys
from pathlib import Path

import numpy as np

from petch.reactor_global.chf3_electron_collisions import (
    CHF3_MASS_AMU,
    derive_nist_evaluated_chf3_replay,
    load_kushner_zhang_2000_chf3_replay,
    load_nist_1999_chf3_drift_curve,
    load_nist_1999_chf3_momentum_transfer,
    load_nist_1999_chf3_total_scattering,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)


ROOT = Path(__file__).resolve().parents[1]


def _grid(deck, cells_scale=1):
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV
        and process.energy_loss_eV < 120.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, 0.005, 0.01, 0.1, 1.0, 10.0, 40.0, 120.0),
        tuple(cells_scale * value for value in (8, 8, 36, 48, 96, 96, 120)),
        inserted_boundaries_eV=thresholds,
    )


def test_committed_source_tables_and_manifests_replay_without_source_files():
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "extract_chf3_electron_evidence.py"),
        "--check",
    ], check=True, cwd=ROOT)


def test_author_working_set_topology_units_and_electron_balance_are_exact():
    replay = load_kushner_zhang_2000_chf3_replay()
    assert len(replay.raw_deck.processes) == 20
    assert replay.process_labels[:4] == ("MOM", "VIB14", "VIB25", "VIB36")
    assert replay.process_labels[-2:] == ("ATT1", "ATT2")
    momentum = replay.raw_deck.processes[0]
    assert momentum.electron_energy_eV[0] == 0.01
    assert momentum.cross_section_m2[0] == 2.135e-17
    assert momentum.mass_ratio == (
        5.48579909065e-4 / CHF3_MASS_AMU
    )
    assert replay.derived_deck.processes[0].electron_energy_eV[0] == 0.0
    assert replay.derived_deck.processes[1].electron_energy_eV[-1] == 120.0
    assert [process.electron_number_change for process in replay.raw_deck.processes[-2:]] == [-1, 0]
    assert replay.supports_independent_branch_validation is False
    assert replay.supports_feature_depth is False


def test_nist_tables_are_exact_si_transcriptions_with_distinct_evidence():
    total = load_nist_1999_chf3_total_scattering()
    momentum = load_nist_1999_chf3_momentum_transfer()
    drift = load_nist_1999_chf3_drift_curve()
    assert total.electron_energy_eV[[0, -1]].tolist() == [0.005, 600.0]
    assert total.cross_section_m2[[0, -1]].tolist() == [
        3321.2e-20, 4.6e-20,
    ]
    assert momentum.electron_energy_eV.tolist() == [10, 15, 20, 25, 30]
    assert momentum.cross_section_m2.tolist() == [
        15.5e-20, 12.4e-20, 11.4e-20, 10.5e-20, 10.1e-20,
    ]
    assert drift.reduced_electric_field_Td[[0, -1]].tolist() == [0.4, 250.0]
    assert drift.drift_velocity_m_s[[0, -1]].tolist() == [220.0, 169000.0]
    assert drift.supports_independent_grade_of_working_set is False


def test_evaluated_backbone_reconstructs_total_scattering_once_below_9ev():
    replay = derive_nist_evaluated_chf3_replay()
    momentum = replay.derived_deck.processes[0]
    total = replay.total_scattering
    probe = total.electron_energy_eV[total.electron_energy_eV <= 9.0]
    assembled = np.interp(
        probe, momentum.electron_energy_eV, momentum.cross_section_m2)
    for process in replay.derived_deck.processes[1:]:
        assembled += np.interp(
            probe, process.electron_energy_eV, process.cross_section_m2)
    np.testing.assert_allclose(
        assembled,
        total.cross_section_m2[:probe.size],
        rtol=2.0e-15,
        atol=1.0e-30,
    )
    np.testing.assert_allclose(
        np.interp(
            replay.momentum_transfer.electron_energy_eV,
            momentum.electron_energy_eV,
            momentum.cross_section_m2,
        ),
        replay.momentum_transfer.cross_section_m2,
        rtol=0.0,
        atol=0.0,
    )


def test_high_energy_closure_is_an_explicit_sensitivity_only_above_30ev():
    constant = derive_nist_evaluated_chf3_replay(
        high_energy_closure="constant_join_ratio")
    tapered = derive_nist_evaluated_chf3_replay(
        high_energy_closure="linear_return_to_working_set_at_120eV")
    first = constant.derived_deck.processes[0]
    second = tapered.derived_deck.processes[0]
    probe = np.linspace(0.0, 30.0, 301)
    np.testing.assert_allclose(
        np.interp(probe, first.electron_energy_eV, first.cross_section_m2),
        np.interp(probe, second.electron_energy_eV, second.cross_section_m2),
        rtol=0.0,
        atol=0.0,
    )
    assert first.cross_section_m2[-1] > second.cross_section_m2[-1]


def test_evaluated_chf3_deck_closes_near_cancellation_and_60td_transport():
    replay = derive_nist_evaluated_chf3_replay()
    solver = DeterministicTwoTermBoltzmannSolver(
        _grid(replay.derived_deck), replay.derived_deck)
    low = solver.solve(TwoTermBoltzmannCondition(
        reduced_electric_field_Td=2.5,
        gas_temperature_K=298.0,
        target_mole_fractions={"CHF3": 1.0},
        growth_model="temporal_growth",
        initial_electron_temperature_eV=0.05,
    ), relative_tolerance=2.0e-6, maximum_tail_population_fraction=1.0e-6)
    high = solver.solve(TwoTermBoltzmannCondition(
        reduced_electric_field_Td=60.0,
        gas_temperature_K=298.0,
        target_mole_fractions={"CHF3": 1.0},
        growth_model="temporal_growth",
        initial_electron_temperature_eV=0.3,
    ), relative_tolerance=2.0e-6, maximum_tail_population_fraction=1.0e-6)
    assert abs(low.net_growth_rate_coefficient_m3_s) < 1.0e-19
    predicted_drift = (
        high.transport_moments.flux_reduced_mobility_m_inv_V_inv_s_inv
        * 60.0e-21
    )
    measured_drift = load_nist_1999_chf3_drift_curve().drift_velocity_m_s[
        np.where(load_nist_1999_chf3_drift_curve().reduced_electric_field_Td == 60.0)[0][0]
    ]
    assert abs(predicted_drift / measured_drift - 1.0) < 0.04
    assert high.supports_direct_swarm_grade is False
    assert replay.supports_reactor_state_prediction is False
