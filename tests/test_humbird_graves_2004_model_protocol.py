import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = (
    ROOT
    / "data"
    / "experimental"
    / "humbird_graves_2004"
    / "model_protocol.json"
)


def test_calibration_and_held_out_panels_are_disjoint_and_depth_is_forbidden():
    protocol = json.loads(PROTOCOL.read_text())
    calibration = set(protocol["calibration"]["source_panels"])
    held_out = set(protocol["held_out"]["source_panels"])
    assert calibration.isdisjoint(held_out)
    assert held_out == {
        "20_percent_F_surface_state",
        "20_percent_F_etch_yield",
    }
    assert "Krueger 825 nm endpoint" in (
        protocol["calibration"]["forbidden_observables"])


def test_protocol_requires_stratification_ledgers_and_nonzero_si_residence():
    protocol = json.loads(PROTOCOL.read_text())
    requirements = " ".join(protocol["model_form_requirements"])
    assert "Si-C transport-layer" in requirements
    assert "Si-F reaction-front" in requirements
    assert "finite CSDA residual-energy" in requirements
    assert protocol["held_out"]["gates"][
        "newly_promoted_silicon_may_not_leave_in_same_step"]
