from hashlib import sha256
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "curated" / "zhu_npg80_open_reactor_v1"


def test_conserved_open_reactor_receipt_is_pinned_and_not_overpromoted():
    payload = json.loads((RESULT_DIR / "central.json").read_text())
    seed = RESULT_DIR / "hydrogen_closed_continuation.json"
    declared = payload["input"]["continuation_state"]
    assert declared["path"] == str(seed.relative_to(ROOT))
    assert declared["sha256"] == sha256(seed.read_bytes()).hexdigest()
    assert payload["input"]["feature_or_sem_target_used"] is False
    assert payload["numerics"]["maximum_normalized_residual"] < 2.0e-6
    assert payload["certification"]["conserved_open_reactor_equations_solved"]
    assert payload["certification"]["wafer_flux_prediction"] is False
    assert payload["certification"]["feature_depth_prediction"] is False


def test_receipt_exposes_daughter_collision_basis_shortfall():
    state = json.loads((RESULT_DIR / "central.json").read_text())["state"]
    fraction = state["electron_collision_basis_neutral_fraction"]
    assert fraction == pytest.approx(0.3755114622864711)
    assert state["implied_total_neutral_reduced_electric_field_Td"] == (
        pytest.approx(state["reduced_electric_field_Td"] * fraction)
    )
    assert state["total_axial_positive_ion_flux_m2_s"] == pytest.approx(
        2.2891234955809833e19)
    assert set(state["neutral_thermal_flux_m2_s"]) == {
        name
        for name, density in state["densities_m3"].items()
        if name not in state["axial_positive_ion_flux_m2_s"]
        and name not in {
            "e", "F-", "F2-", "O-", "SF2-", "SF3-", "SF4-", "SF5-", "SF6-"
        }
    }
