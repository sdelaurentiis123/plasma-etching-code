from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated"
    / "zhu_npg80_profile_convergence_refinement_v1" / "audit.json"
)


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_refinement_is_target_free_and_checksum_bound():
    audit = _load(AUDIT)

    assert audit["target_sem_used"] is False
    assert audit["target_depth_used"] is False
    assert audit["supports_absolute_target_profile_prediction"] is False
    assert audit["supports_atomic_accuracy"] is False
    for receipt in audit["inputs"]["cases"].values():
        path = ROOT / receipt["path"]
        assert sha256(path.read_bytes()).hexdigest() == receipt["sha256"]
    for name in (
        "refinement_preregistration", "reference_case",
        "external_physical_shape_witness",
    ):
        receipt = audit["inputs"][name]
        path = ROOT / receipt["path"]
        assert sha256(path.read_bytes()).hexdigest() == receipt["sha256"]


def test_refinement_separates_physical_shape_from_missing_surface_physics():
    audit = _load(AUDIT)
    diagnostic = audit["conditional_shape_diagnostic"]

    assert diagnostic["surface_mechanism_growth_enabled"] is False
    assert diagnostic["surface_product_redeposition_enabled"] is False
    motion = diagnostic["material_motion_invariant"]
    assert motion["allowed_normal_velocity"] == "removal-only and nonnegative"
    assert motion["negative_removal_velocity_rejected_by_mechanism"] is True
    assert motion[
        "wider_lower_section_requires_slow_lower-wall_recession_not_growth"
    ] is True
    assert motion["pinned_materials"] == ["Cr mask", "fused-silica substrate"]
    assert "self-shadowing" in diagnostic["candidate_mechanism_if_present"]
    witness = diagnostic["external_physical_shape_witness"]
    assert witness["bibkey"] == "ji-2024-tio2-hierarchical"
    assert witness["same_feed_constituents"] is True
    assert witness["same_reactor_or_condition"] is False
    assert witness["coefficient_transfer_allowed"] is False
    assert "cannot validate" in diagnostic["interpretation_limit"]
    assert diagnostic["ultrafine_bottom_minus_top_cd_nm"] == (
        diagnostic["ultrafine_bottom_cd_nm"]
        - diagnostic["ultrafine_top_cd_nm"]
    )
    assert diagnostic["qualitative_footing_present"] == (
        diagnostic["ultrafine_bottom_minus_top_cd_nm"] > 0.0
    )
    if audit["bottom_cd_numerically_certified_for_sentinel"]:
        assert "numerically certified" in diagnostic["classification"]
    else:
        assert "numerically unresolved" in diagnostic["classification"]


def test_refinement_uses_the_frozen_numerical_ladder():
    audit = _load(AUDIT)
    cases = {
        name: _load(ROOT / receipt["path"])
        for name, receipt in audit["inputs"]["cases"].items()
    }

    assert cases["fine_dt4_cuda_replicate"]["case_specification"] == {
        "mesh_spacing_nm": 10.0,
        "maximum_step_s": 4.0,
        "required_execution_device": "cuda",
    }
    assert cases["fine_dt2"]["case_specification"] == {
        "mesh_spacing_nm": 10.0,
        "maximum_step_s": 2.0,
        "required_execution_device": "cuda",
    }
    assert cases["ultrafine_dt2"]["case_specification"] == {
        "mesh_spacing_nm": 5.0,
        "maximum_step_s": 2.0,
        "required_execution_device": "cuda",
    }
    for case in cases.values():
        assert case["target_sem_used"] is False
        assert case["target_depth_used"] is False
        assert case["profile"]["tio2_clearance_detected"] is False
        assert case["profile"]["validity"][
            "parameter_evidence_supports_prediction"] is False


def test_failed_bottom_gate_is_localized_without_reclassifying_it_as_passed():
    audit = _load(AUDIT)
    diagnostic = audit["post_result_localization_diagnostic"]

    assert audit["gate_results"]["ultrafine_grid_cd"] is False
    assert audit["all_numerical_gates_pass"] is False
    assert audit["bottom_cd_numerically_certified_for_sentinel"] is False
    assert diagnostic["changes_frozen_gate_result"] is False
    assert diagnostic["exploratory_not_preregistered"] is True
    assert diagnostic["sentinel_cleared_tio2"] is False
    assert diagnostic["remaining_tio2_below_etched_floor_nm"] > 400.0
    assert diagnostic["independent_pillar_bottom_exists"] is False
    assert diagnostic["body_maximum_10nm_to_5nm_width_change_nm"] < 3.0
    assert diagnostic["near_floor_maximum_10nm_to_5nm_width_change_nm"] > 10.0
    assert diagnostic["first_fraction_exceeding_frozen_5nm_cd_gate"] == 0.85
    assert diagnostic["body_profile_within_frozen_timestep_cd_tolerance"] is True
    assert diagnostic["nonphysical_growth_or_interface_drift_detected"] is False
    assert diagnostic["target_pillar_bottom_certified"] is False
