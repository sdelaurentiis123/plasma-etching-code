from scripts.audit_zhu_npg80_daughter_collision_gap import build_receipt


def test_daughter_collision_gap_replays_reported_parent_basis():
    receipt = build_receipt()
    board = receipt["power_board"]
    assert [row["absorbed_power_W"] for row in board] == [60, 90, 105, 120]
    assert board[0]["energy_resolved_two_term_eedf_fraction"] > 0.5
    assert board[-1]["energy_resolved_two_term_eedf_fraction"] < 0.13
    for row in board:
        assert abs(sum((
            row["energy_resolved_two_term_eedf_fraction"],
            row["electron_rate_regression_only_fraction"],
            row["electron_transparent_fraction"],
        )) - 1.0) < 2.0e-14


def test_hf_is_dominant_missing_120W_electron_target():
    receipt = build_receipt()
    state = receipt["accepted_120W_state"]
    ranked = state["ranked_missing_targets"]
    assert ranked[0]["species"] == "HF"
    assert ranked[0]["electron_kinetic_tier"] == (
        "electron_rate_regression_only"
    )
    assert ranked[0]["neutral_fraction"] > 0.59
    assert ranked[0]["reduced_electron_reaction_count"] == 1
    assert ranked[1]["species"] == "F2"
    assert ranked[1]["electron_kinetic_tier"] == (
        "electron_rate_regression_only"
    )


def test_ranked_milestones_quantify_minimum_collision_basis_work():
    receipt = build_receipt()
    milestones = receipt["accepted_120W_state"][
        "collision_basis_closure_milestones"
    ]
    assert milestones[0]["threshold"] == 0.5
    assert milestones[0]["minimum_ranked_daughter_targets"] == ["HF"]
    assert milestones[1]["threshold"] == 0.75
    assert milestones[1]["minimum_ranked_daughter_targets"] == ["HF", "F2"]
    assert receipt["sem_or_depth_target_used"] is False
    assert receipt["priority_verdict"]["supports_feature_depth_prediction"] is False
