from scripts.audit_zhu_npg80_absorbed_power_ensemble import build_receipt


def test_power_states_conserve_and_close_sheath_fixed_points():
    receipt = build_receipt()
    assert receipt["sem_or_depth_target_used"] is False
    assert len(receipt["state_board"]) == 4
    certification = receipt["certification"]
    assert certification["all_accepted_states_conserve_below_2e_6"] is True
    assert certification["all_sheath_fixed_points_below_0p01_V"] is True
    assert certification["supports_unique_profile_depth"] is False


def test_absorbed_power_response_and_depth_envelope_are_not_collapsed():
    receipt = build_receipt()
    response = receipt["power_response"]
    assert response["ion_flux_monotonically_increases_with_absorbed_power"]
    assert response["F_flux_monotonically_increases_with_absorbed_power"]
    assert response["central_ion_flux_120W_to_60W_ratio"] > 1.5
    low, high = response["required_surface_yield_envelope"]
    assert low < 0.9
    assert high > 1.6
    verdict = receipt["corrected_depth_verdict"]
    assert verdict["original_frozen_binary_call_preserved"] is True
    assert verdict["current_clearance_supported_over_full_power_envelope"] is False


def test_extended_field_coordinate_exposes_daughter_collision_gap():
    audit = build_receipt()["field_domain_audit"]
    assert audit["accepted_120W_represented_field_Td"] > 600.0
    assert audit["accepted_120W_implied_total_neutral_field_Td"] < 100.0
    assert audit["accepted_120W_collision_basis_fraction"] < 0.15
    assert build_receipt()["independent_feature_response"][
        "cw_layout_depth_ratio"
    ] > 1.45
