from pathlib import Path
import sys

import numpy as np


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import charging_task0a_ensemble_response as task0a  # noqa: E402
import charging_task0b_electron_boundary as task0b  # noqa: E402


def test_task0a_direction_set_is_fixed_normalized_and_spans_required_classes(tmp_path):
    rng = np.random.default_rng(3)
    jacobian = rng.normal(size=(7, 7))
    residual = np.arange(7, dtype=float)
    archive = tmp_path / "jacobian.npz"
    np.savez(archive, jacobian=jacobian, residual=residual)

    loaded_residual, directions = task0a.direction_set(archive)

    assert np.array_equal(loaded_residual, residual)
    assert len(directions) == 10
    assert [kind for _, kind, _ in directions].count("worst_coordinate") == 2
    assert [kind for _, kind, _ in directions].count("random") == 5
    assert [kind for _, kind, _ in directions].count("dominant") == 3
    assert np.allclose([np.linalg.norm(vector) for _, _, vector in directions], 1.0)


def test_task0b_barrier_is_exponential_only_while_electrons_are_repelled():
    voltage = np.array([-8.0, -4.0, 0.0, 4.0])

    factor = task0b.barrier_factor(voltage, temperature_eV=4.0)

    assert np.allclose(factor, [np.exp(-2.0), np.exp(-1.0), 1.0, 1.0])


def test_switch_decomposition_keeps_signed_condition_effect_in_synthetic_case():
    nonswitch = np.eye(2)
    injected_switch = np.diag([0.0, -0.99])
    total = nonswitch + injected_switch

    metrics = task0a.response_decomposition_metrics(total, injected_switch)

    assert np.isclose(metrics["condition"], 100.0)
    assert np.isclose(metrics["no_switch_condition"], 1.0)
    assert metrics["signed_log_condition_change"] > 0.0
    assert metrics["switch_component_energy_fraction"] > 0.0

    beneficial_switch = np.diag([0.0, 0.99])
    beneficial = task0a.response_decomposition_metrics(total + beneficial_switch, beneficial_switch)
    assert beneficial["signed_log_condition_change"] < 0.0
