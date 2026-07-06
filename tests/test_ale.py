"""Gate the Si/Cl2/Ar+ ALE ROM (src/petch/ale.py) against Vella-Graves 2025. The qualitative gates
(the ALE window + self-limitation loss + synergy) are the physics; they must hold from the yields
alone with nothing tuned."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from petch.ale import run_ale, synergy, yields


def test_ale_window_flat_and_sublimonolayer():
    """Inside the 15-20 eV window EPC is flat and sub-monolayer (< 1.36 A); it stays under ~1 A."""
    epc = {E: run_ale(E)["epc"] for E in (15, 17.5, 20)}
    assert all(0.2 < v < 1.36 for v in epc.values()), epc      # sub-monolayer window
    assert epc[20] > epc[15]                                    # gentle rise across the window


def test_self_limitation_lost_above_window():
    """Above 20 eV the bare-Si sputter channel turns on and EPC rises steeply (loss of ALE)."""
    assert run_ale(22.5)["epc"] > 1.4 * run_ale(20)["epc"]
    assert run_ale(30)["epc"] > 3.0 * run_ale(20)["epc"]


def test_synergy_high_in_window_collapses_outside():
    """Kanarik synergy S -> ~1 inside the window (ideal ALE), collapses in the sputter regime."""
    assert synergy(17.5)["synergy"] > 0.9
    assert synergy(30)["synergy"] < 0.4


def test_yields_clamped_nonnegative():
    """Below threshold every yield is clamped to 0 (no negative sputter yields)."""
    y = yields(10.0)
    assert all(y[k] >= 0.0 for k in ("Y_Cl", "Y_SiCl", "Y_SiCl2", "Y_Si"))


def test_absolute_epc_near_paper_at_anchor():
    """Absolute EPC at the 17.5/20 eV anchors is within ~25% of the ROM figure (0.7 / 0.9 A)."""
    assert abs(run_ale(17.5)["epc"] - 0.7) < 0.25
    assert abs(run_ale(20.0)["epc"] - 0.9) < 0.25
