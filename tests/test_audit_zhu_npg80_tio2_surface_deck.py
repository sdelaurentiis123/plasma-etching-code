from scripts.audit_zhu_npg80_tio2_surface_deck import build


def test_oxford_surface_deck_is_executable_but_fail_closed():
    audit = build()

    assert audit["target_sem_used"] is False
    assert audit["target_depth_used"] is False
    assert audit["executable_contract"]["reduced_sensitivity_execution_available"] is True
    assert audit["executable_contract"]["silent_sio2_coefficient_transfer_allowed"] is False
    assert audit["executable_contract"]["numerical_defaults_supplied"] is False
    assert audit["executable_contract"]["competitive_oxygen_state_implemented"] is True
    assert audit["supports_absolute_oxford_profile_prediction"] is False
    assert audit["supports_atomic_accuracy"] is False


def test_every_required_surface_coefficient_remains_explicit():
    audit = build()
    slots = audit["parameter_slots"]

    assert set(audit["target_parameters_not_identified"]) == set(slots)
    assert all(
        status["numerical_target_value_identified"] is False
        for status in slots.values()
    )
    assert "competitive_oxygen_blocking_or_cleanup_surface_state" not in (
        audit["unresolved_model_form"]
    )
    assert "chemistry_dependent_roughness_evolution" in audit["unresolved_model_form"]
    bare = slots["bare_sio2_yield"]
    assert bare["reference_formula_unit_yield_at_276_eV"] == 0.192143
    assert bare["reference_evidence_class"] == "digitized_semiempirical_fit_curve"
    assert bare["reference_transferable_as_target_coefficient"] is False
    assert bare["numerical_target_value_identified"] is False
