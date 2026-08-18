import math

from scripts.audit_zhu_npg80_tio2_depth_gate import build_receipt


def test_blind_depth_gate_is_target_free_and_atom_counted():
    receipt = build_receipt()

    assert receipt["frozen_before_specific_condition_sem"] is True
    assert receipt["sem_target_used"] is False
    assert receipt["measured_depth_target_used"] is False
    assert receipt["coefficient_selected_from_target"] is None
    assert receipt["blind_forecast"]["predicted_tio2_depth_nm"] == 700.0
    required = receipt["central_reactor_dose_clearance_gate"][
        "required_blanket_formula_units_per_positive_ion"
    ]
    assert math.isclose(required[0], 0.624481965032084)
    assert math.isclose(required[1], 0.797415432271738)


def test_depth_gate_keeps_rate_mask_and_machine_transfer_boundaries_visible():
    receipt = build_receipt()

    analog = receipt["independent_tio2_process_comparison"]
    assert analog["same_recipe_or_machine"] is False
    assert analog["twenty_minute_film_capped_depth_comparison_nm"] == [
        682.5,
        700.0,
    ]
    mask = receipt["cr_mask_survival_straddle"]
    assert mask["closest_feature_witness_supported_tio2_nm"] == 630.0
    assert mask["source_power_sweep_supported_tio2_nm"] > 700.0
    sf6 = receipt["independent_sf6_direction_evidence"]
    assert sf6["supports_direction_not_magnitude"] is True
    assert sf6["reported_tio2_rate_nm_min"]["SF6"] > (
        sf6["reported_tio2_rate_nm_min"]["CHF3"]
    )
    assert receipt["inputs"]["reactor_sensitivity_state"][
        "certified_as_local_wafer_flux"
    ] is False
