from scripts.audit_zhu_npg80_tio2_surface_topology import build


def test_surface_topology_audit_is_target_free_and_fail_closed():
    audit = build()

    assert audit["target_sem_used"] is False
    assert audit["target_depth_used"] is False
    assert audit["supports_absolute_oxford_profile_prediction"] is False
    assert audit["supports_atomic_accuracy"] is False
    assert audit["current_conditional_mechanism"][
        "physical_within_declared_scope"
    ] is True
    assert audit["current_conditional_mechanism"][
        "sufficient_for_target_surface_prediction"
    ] is False


def test_independent_boards_force_energy_passivation_mask_and_loading_physics():
    audit = build()
    required = audit["minimum_physical_model_contract"]
    missing = set(audit["missing_from_current_conditional_mechanism"])

    for name in (
        "ion_energy_dependent_surface_yield",
        "neutral_radical_surface_reactions",
        "fluorinated_tio2_surface_inventory",
        "passivation_inventory_with_physical_thickness",
        "ion_assisted_passivation_removal",
        "cr_mask_geometry_evolution",
        "pattern_dependent_neutral_and_ion_transport",
    ):
        assert required[name]["required"] is True
        assert name in missing
    assert audit["experimental_discriminants"]["spacing"][
        "sharp_threshold_identified"
    ] is False
    assert audit["deterministic_differentiable_state_design"][
        "monte_carlo_required"
    ] is False
