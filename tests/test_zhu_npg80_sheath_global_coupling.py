from scripts.audit_zhu_npg80_sheath_global_coupling import build_receipt


def test_sheath_global_fixed_point_is_target_free_and_conserved():
    receipt = build_receipt()
    assert receipt["target_outcome_used"] is False
    assert receipt["fixed_point"]["converged_below_0p01_V"] is True
    assert abs(receipt["wall_ledger"][
        "particle_closure_relative_residual"
    ]) < 2.0e-14
    assert receipt["certification"][
        "supports_unique_absolute_depth_prediction"
    ] is False
    assert receipt["axisymmetric_wafer_state"][
        "absolute_target_wafer_flux_supported"
    ] is False


def test_wall_resolution_strengthens_but_does_not_fit_clearance():
    receipt = build_receipt()
    effect = receipt["effect_of_corrected_wall_power"]
    assert effect["charged_wall_power_ratio"] < 1.0
    assert effect["electron_density_ratio"] > 1.0
    assert effect["positive_ion_flux_ratio"] > 1.0
    assert effect["neutral_F_flux_ratio"] > 1.0
    required = receipt["conditional_tio2_dose"][
        "required_blanket_formula_units_per_positive_ion"
    ]
    assert 0.9 < min(required) < max(required) < 1.3
    optic_required = receipt["conditional_tio2_dose"][
        "required_central_3mm_formula_units_per_positive_ion"
    ]
    assert max(optic_required) < max(required)
