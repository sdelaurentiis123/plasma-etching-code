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

    levinson = boards["Levinson 1997 Ar+/Cl2/Si MIBE features"]
    assert levinson["point_count"] == 3
    assert levinson["reported_absolute_feature_depths"] is True
    assert levinson["absolute_depth_or_depth_per_dose"] is False
    assert levinson["feature_profile_test"] is False
    assert levinson["original_pixels_archived"] is False
    assert levinson["formal_pass_label"] is False

    woo = boards["Woo 2024 CF4/C4F6/He patterned SiO2"]
    assert woo["point_count"] == 10
    assert woo["absolute_depth_or_depth_per_dose"] is True
    assert woo["original_pixels_archived"] is True
    assert woo["model_predictions_completed"] == 0
    assert woo["source_internal_consistency_passed"] is False
    assert woo["formal_pass_label"] is False

    mahorowala = boards["Mahorowala 1998 Cl2/oxide-mask Si features"]
    assert mahorowala["source_run_count"] == 13
    assert mahorowala["point_count"] == 11
    assert mahorowala["derived_poly_si_depth_range_nm"] == [112.5, 459.375]
    assert mahorowala["source_profile_panel_count"] == 13
    assert mahorowala["source_feature_profiles_available"] is True
    assert mahorowala["original_resolution_visual_audit_passed"] is True
    assert mahorowala["feature_profile_test"] is False
    assert mahorowala["model_predictions_completed"] == 0
    assert mahorowala["formal_pass_label"] is False


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
    assert summary[
        "controlled_feature_boards_blocked_by_missing_time_or_fluence"] == 1
    assert summary["fixed_time_chlorine_absolute_rate_targets_ready"] == 11
    assert summary["fixed_time_chlorine_sem_profile_targets_ready"] == 13
    assert summary[
        "c4f6_same_reactor_absolute_patterned_rate_targets_ready"] == 10
