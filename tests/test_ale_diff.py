"""Gate the differentiable ALE model (src/petch/ale_diff.py): autograd correctness, forward parity
with the numpy reference, and gradient-based inverse recovery of the ion energy."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import torch
from petch.ale_diff import epc_torch, dEPC_dE, invert_energy_for_epc
from petch.ale import run_ale


def test_forward_parity_with_numpy():
    """The fast torch transient+leak model matches the accurate numpy ale.py within ~5%."""
    for E in (17.5, 20.0, 25.0, 30.0):
        t = float(epc_torch(torch.tensor(E, dtype=torch.float64)).detach())
        n = run_ale(E)["epc"]
        assert abs(t - n) < 0.05 * max(n, 0.5), (E, t, n)


def test_autograd_matches_finite_difference():
    """Reverse-mode dEPC/dE equals the central finite difference (autograd is correct)."""
    for E in (17.5, 22.5, 27.5):
        _, g = dEPC_dE(E)
        h = 1e-3
        fd = (float(epc_torch(torch.tensor(E + h, dtype=torch.float64)).detach())
              - float(epc_torch(torch.tensor(E - h, dtype=torch.float64)).detach())) / (2 * h)
        assert abs(g - fd) < 1e-3, (E, g, fd)


def test_gradient_positive_in_sputter_regime():
    """EPC rises with energy above the window -> positive gradient there."""
    assert dEPC_dE(25.0)[1] > 0.1


def test_inverse_recovers_energy_self_consistent():
    """Gradient-based inversion recovers the energy that produced a target EPC (same model)."""
    for E_true in (18.0, 21.0, 26.0):
        target = float(epc_torch(torch.tensor(E_true, dtype=torch.float64)).detach())
        E_sol, _, _ = invert_energy_for_epc(target, E0=23.0)
        assert abs(E_sol - E_true) < 0.1, (E_true, E_sol)
