import math

from scripts.audit_zhu_npg80_tio2_pre_sem import (
    build_receipt,
    load_manifest,
)
from scripts.audit_zhu_npg80_cf3_collision_scale import (
    build_receipt as build_cf3_collision_receipt,
)
from scripts.audit_zhu_npg80_chf2_mobility_scale import (
    build_receipt as build_chf2_mobility_receipt,
)
from scripts.audit_zhu_npg80_tio2_analog_board import (
    build_receipt as build_tio2_analog_receipt,
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


def test_machine_family_bias_evidence_enables_sensitivity_not_fake_prediction():
    receipt = build_receipt(load_manifest())
    evidence = receipt["machine_family_self_bias_evidence"]
    gates = receipt["identifiability_gates"]

    assert evidence["available"] is True
    assert evidence["matched_chemistry_reduced_drive_anchor_V"] == 276.0
    assert evidence["exact_tool_conditioning_drift_is_censored"] is True
    assert gates["machine_family_self_bias_sensitivity_available"] is True
    assert gates["achieved_dc_self_bias_measured"] is False
    assert gates["supports_absolute_depth_prediction"] is False


def test_cf3_collision_scale_is_target_free_and_partial_by_construction():
    receipt = build_cf3_collision_receipt()
    sample = next(
        item for item in receipt["samples"]
        if item["laboratory_energy_eV"] == 200.0
    )
    gates = receipt["identifiability_gates"]

    assert receipt["sem_target_used"] is False
    assert receipt["measured_depth_target_used"] is False
    assert receipt["energy_convention"]["laboratory_to_relative_factor"] == (
        70.0 / 139.0
    )
    assert 0.18 < sample["reactive_optical_depth_per_1mm"]["central"] < 0.20
    assert gates["measured_cf3_chf3_reactive_kernel_available"] is True
    assert gates["elastic_or_momentum_transfer_kernel_available"] is False
    assert gates["supports_target_iead"] is False
    assert gates["supports_absolute_depth_prediction"] is False


def test_chf2_mobility_scale_is_measured_target_free_and_partial():
    receipt = build_chf2_mobility_receipt()
    low = next(
        item for item in receipt["samples"]
        if item["reduced_field_Td"] == 100.0
    )
    gates = receipt["identifiability_gates"]

    assert receipt["sem_target_used"] is False
    assert receipt["measured_depth_target_used"] is False
    assert 0.48 < low["reduced_mobility_cm2_V_s"] < 0.50
    assert 0.08e-3 < low["drift_relaxation_length_m"] < 0.13e-3
    assert gates["measured_chf2_chf3_swarm_mobility_available"] is True
    assert gates["elastic_differential_cross_section_available"] is False
    assert gates["supports_target_iead"] is False
    assert gates["supports_absolute_depth_prediction"] is False


def test_pre_sem_receipt_includes_measured_chf2_mobility_without_overclaim():
    receipt = build_receipt(load_manifest())
    evidence = receipt["measured_molecular_mobility_evidence"]
    gates = receipt["identifiability_gates"]

    assert evidence["available"] is True
    assert evidence["elastic_differential_cross_section"] is False
    assert gates["measured_chf2_chf3_swarm_mobility_available"] is True
    assert gates["complete_molecular_ion_transport_available"] is False
    assert gates["supports_absolute_depth_prediction"] is False


def test_exact_tio2_analog_board_brackets_clearance_without_transfer():
    receipt = build_tio2_analog_receipt()
    closest = receipt["closest_stack_witness"]
    sweep = receipt["source_power_sweep_interpolation"]
    depths = receipt["source_feature_depth_board"]
    gates = receipt["identifiability_gates"]

    assert receipt["sem_target_used"] is False
    assert closest["system_model"] == "Fluor Z401S"
    assert closest["dc_bias_V_signed"] == -950.0
    assert receipt["target_similarity"]["same_cr_mask_thickness"] is True
    assert 1.06 < receipt["target_similarity"][
        "source_to_target_reduced_drive_ratio"] < 1.07
    assert 52.0 < sweep["source_system_tio2_rate_nm_min"] < 53.0
    assert sweep["source_system_cr_residual_after_700nm_nm"] > 5.0
    assert depths["minimum_implied_rate_nm_min"] < 35.0
    assert depths["maximum_implied_rate_nm_min"] > 40.0
    assert sweep["supports_target_prediction"] is False
    assert gates["supports_absolute_target_depth_prediction"] is False


def test_pre_sem_receipt_links_exact_tio2_board_without_promoting_it():
    receipt = build_receipt(load_manifest())
    analogs = receipt["audited_tio2_process_analogs"]
    gates = receipt["identifiability_gates"]

    assert analogs["available"] is True
    assert analogs["closest_source_dc_bias_V_signed"] == -950.0
    assert analogs["transferred_as_target_coefficient"] is False
    assert gates["adjacent_tio2_process_response_board_available"] is True
    assert gates["tio2_surface_law_measured_or_validated_for_condition"] is False


def test_pre_sem_receipt_links_binary_clearance_call_without_profile_overclaim():
    receipt = build_receipt(load_manifest())
    forecast = receipt["blind_tio2_clearance_forecast"]
    gates = receipt["identifiability_gates"]

    assert forecast["available"] is True
    assert forecast["forecast_type"] == "film-capped binary clearance call"
    assert forecast["predicted_film_capped_tio2_depth_nm"] == 700.0
    assert forecast["promoted_to_absolute_feature_profile_prediction"] is False
    assert gates["target_free_binary_clearance_forecast_frozen"] is True
    assert gates["supports_absolute_depth_prediction"] is False
