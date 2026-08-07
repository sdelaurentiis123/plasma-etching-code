import json
from pathlib import Path

from scripts.audit_levinson_1997_feature_identifiability import build_audit


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated"
    / "levinson_1997_feature_identifiability" / "audit.json"
)


def test_committed_levinson_audit_is_exact_replay():
    assert json.loads(AUDIT.read_text(encoding="utf-8")) == build_audit()


def test_article_cannot_be_relabeled_as_absolute_depth_gate():
    audit = build_audit()
    depth = audit["absolute_depth_identifiability"]

    assert not depth["case_specific_ion_flux_published"]
    assert not depth["case_specific_exposure_time_published"]
    assert not depth["case_specific_ion_fluence_published"]
    assert not depth["apparatus_design_current_is_case_boundary"]
    assert not depth["target_depth_may_select_simulation_time"]
    assert not depth["article_identifies_absolute_feature_depth_prediction"]
    assert not audit["allowed_now"]["absolute_feature_depth_test"]


def test_controlled_beam_board_keeps_all_caption_cases():
    cases = build_audit()["figure11_cases"]
    assert [case["figure_panel"] for case in cases] == ["11a", "11b", "11c"]
    assert [case["ion_energy_eV"] for case in cases] == [100.0, 100.0, 100.0]
    assert [case["unshadowed_ion_to_neutral_flux_ratio"] for case in cases] == [
        0.004, 0.008, 0.008
    ]
    assert all(case["measured_depth_to_width_ratio"] > 0 for case in cases)
