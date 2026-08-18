import math

from scripts.audit_zhu_npg80_tio2_pre_sem import (
    build_receipt,
    load_manifest,
)


def test_operator_corrections_override_screenshot_values():
    manifest = load_manifest()

    assert manifest["displayed_then_corrected"]["pressure_Torr"] == {
        "screenshot": 0.045,
        "operator_correction": 0.03,
    }
    assert manifest["displayed_then_corrected"]["SF6_sccm"] == {
        "screenshot": 2.0,
        "operator_correction": 5.0,
    }
    assert manifest["process"]["pressure_Torr"] == 0.03
    assert manifest["process"]["gases_sccm"]["SF6"] == 5.0
    assert manifest["experimental_program"]["research_group"] == (
        "Nanfang Yu Group, Columbia University"
    )
    assert manifest["stack"]["film_deposition_method"] == (
        "atomic layer deposition (ALD)"
    )
    assert manifest["stack"]["substrate_material"] == "fused silica"
    assert manifest["pattern"]["layout_id"] is None


def test_pre_sem_receipt_fixes_recipe_arithmetic_without_using_target():
    receipt = build_receipt(load_manifest())
    recipe = receipt["authoritative_recipe"]
    gates = receipt["necessary_conditions_for_full_clear"]

    assert receipt["sem_target_used"] is False
    assert receipt["measured_depth_target_used"] is False
    assert receipt["coefficient_selected_from_target"] is None
    assert recipe["etch_time_min"] == 20.0
    assert recipe["film_deposition_method"] == "atomic layer deposition (ALD)"
    assert recipe["substrate_material"] == "fused silica"
    assert recipe["total_process_flow_sccm"] == 61.0
    assert math.isclose(recipe["gas_flow_fractions"]["CHF3"], 55.0 / 61.0)
    assert math.isclose(recipe["gas_flow_fractions"]["SF6"], 5.0 / 61.0)
    assert math.isclose(recipe["gas_flow_fractions"]["O2"], 1.0 / 61.0)
    assert gates["minimum_effective_tio2_rate_nm_min"] == 35.0
    assert math.isclose(
        gates["minimum_zero_margin_tio2_to_cr_selectivity"],
        700.0 / 45.0,
    )


def test_adjacent_literature_scenario_exposes_mask_survival_failure():
    receipt = build_receipt(load_manifest())
    comparison = receipt["adjacent_literature_comparison"]

    assert comparison["rate_implied_clear_time_min"] == 17.5
    assert comparison["mask_supported_tio2_relief_nm"] == 630.0
    assert comparison["mask_required_for_700nm_clear_nm"] == 50.0
    assert comparison["mask_shortfall_for_full_clear_nm"] == 5.0
    assert comparison["mask_exhaustion_time_at_comparison_rate_min"] == 15.75
    assert comparison["time_after_comparison_mask_exhaustion_min"] == 4.25
    assert (
        receipt["identifiability_gates"]["supports_absolute_depth_prediction"]
        is False
    )


def test_same_group_nature_device_is_context_not_a_transferred_target():
    evidence = build_receipt(load_manifest())["adjacent_same_group_device_evidence"]

    assert evidence["source"] == "holman-2026-nature-metasurface-tweezers"
    assert evidence["reported_tio2_meta_atom_height_nm"] == 750.0
    assert evidence["reported_width_range_nm"] == [100.0, 190.0]
    assert evidence["reported_unit_cell_nm"] == 290.0
    assert evidence["used_as_condition_matched_target"] is False
