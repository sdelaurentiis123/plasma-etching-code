from scripts.audit_zhu_npg80_self_bias_global_ensemble import build_receipt


def test_all_self_bias_states_close_without_target_outcomes():
    receipt = build_receipt()
    assert receipt["sem_or_depth_target_used"] is False
    assert len(receipt["state_board"]) == 7
    assert max(
        row["fixed_point_residual_V"] for row in receipt["state_board"]
    ) < 0.01
    certification = receipt["certification"]
    assert certification["all_reactor_conservation_gates_passed"] is True
    assert certification["all_axisymmetric_global_residuals_below_1_percent"] is True
    assert certification["supports_unique_profile_depth"] is False


def test_fixed_power_voltage_response_and_history_integration_are_explicit():
    receipt = build_receipt()
    response = receipt["voltage_response"]
    assert response[
        "ion_flux_monotonically_decreases_with_bias_at_fixed_absorbed_power"
    ] is True
    assert response[
        "F_flux_monotonically_decreases_with_bias_at_fixed_absorbed_power"
    ] is True
    assert response["ion_flux_200V_to_400V_ratio"] > 1.1
    assert response["F_flux_200V_to_400V_ratio"] > 2.5
    history = receipt["exact_ngp80_conditioning_threshold_history"]
    assert history["endpoints_are_censor_thresholds"] is True
    assert history["simpson_node_bias_magnitude_V"] == [300.0, 250.0, 200.0]
