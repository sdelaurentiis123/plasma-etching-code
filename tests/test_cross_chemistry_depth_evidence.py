from scripts.audit_cross_chemistry_depth_evidence import build_audit


def test_cross_chemistry_ledger_keeps_evidence_classes_separate():
    audit = build_audit()
    boards = {row["board"]: row for row in audit["boards"]}
    tinacba = boards["Tinacba 2021 SF5+ beam"]
    assert tinacba["point_count"] == 4
    assert tinacba["fit_to_compared_depth_or_yield"] is False
    assert tinacba["feature_profile_test"] is False

    vella = boards["Vella-Hao / Kounis-Melas Si ALE"]
    assert vella["point_count"] == 3
    assert vella["fit_to_compared_depth_or_yield"] is False
    assert vella["feature_profile_test"] is False


def test_no_held_out_feature_pass_is_claimed_before_prediction():
    audit = build_audit()
    boards = {row["board"]: row for row in audit["boards"]}
    yoshie = boards["Yoshie 2023 cyclic SF6/C4F8 trench"]
    assert yoshie["point_count"] == 49
    assert yoshie["value_blind_held_out"] is True
    assert yoshie["model_predictions_completed"] == 0
    assert yoshie["formal_pass_label"] is False

    summary = audit["summary"]
    assert summary[
        "formal_held_out_feature_or_profile_predictions_completed"] == 0
    assert summary["formal_held_out_feature_or_profile_passes"] == 0
    assert summary[
        "independent_non_krueger_chemistry_families_with_absolute_surface_depth_evidence"] == 2
