import math

from scripts.audit_zhu_npg80_geometry_correction import build_receipt


def test_source_correct_geometry_is_target_free_and_conservation_gated():
    receipt = build_receipt()
    old = receipt["immutable_preregistered_v1"]
    new = receipt["source_geometry_v2"]
    assert old["forecast_changed_retroactively"] is False
    assert old["state"]["electrode_diameter_mm"] == 170.0
    assert new["central"]["electrode_diameter_mm"] == 240.0
    assert new["published_rate_alternate"]["electrode_diameter_mm"] == 240.0
    assert new["central"]["maximum_normalized_residual"] < 2.0e-6
    assert new["published_rate_alternate"]["maximum_normalized_residual"] < 2.0e-6


def test_geometry_dominates_total_ion_flux_over_published_rate_branch():
    receipt = build_receipt()
    geometry = receipt["geometry_effect_central_over_v1"]
    chemistry = receipt["published_chf3_f_rate_branch_effect_at_240mm"]
    assert math.isclose(geometry["positive_ion_flux_ratio"], 0.5464487948209091)
    assert math.isclose(
        chemistry["positive_ion_flux_ratio"],
        1.0118524175432844,
        rel_tol=1.0e-12,
    )
    assert geometry["neutral_F_flux_ratio"] < 0.4
    assert chemistry["neutral_F_flux_ratio"] > 1.3


def test_revised_dose_gate_fails_closed_on_profile_authority():
    receipt = build_receipt()
    required = receipt["revised_conditional_clearance_ledger"][
        "central_required_formula_units_per_positive_ion"
    ]
    assert math.isclose(required[0], 1.142800516628002)
    assert math.isclose(required[1], 1.4592683520019103)
    assert receipt["verdict"]["current_reactor_only_call"] == (
        "clearance remains unresolved"
    )
    assert receipt["verdict"]["spatial_profile_authorized"] is False
