from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "zhu_npg80_conditional_profiles_v1"
    / "audit.json"
)


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_conditional_profile_board_is_blind_and_rate_factorized():
    audit = _load(AUDIT)

    assert audit["target_sem_used"] is False
    assert audit["target_depth_used"] is False
    assert audit["coefficient_selected_from_target"] is None
    assert "not an Oxford absolute-depth prediction" in audit[
        "surface_scale_status"
    ]
    assert audit["geometry_evolution_mode"] == (
        "single evolving TiO2 gas-solid union with pinned Cr and fused silica"
    )
    factorization = audit["trajectory_factorization"]
    assert factorization["governing_invariance"] == (
        "geometry depends on blanket_rate_times_time"
    )
    assert factorization["exact_for_declared_rate_normalized_law"] is True
    assert factorization["independent_trajectories"] == 28
    assert factorization["reported_profile_endpoints"] == 56
    resolution = audit["numerical_resolution_status"]
    assert resolution["production_mesh_spacing_nm"] == 20.0
    assert resolution["minimum_lateral_cells_for_cd_claim"] == 6.0
    assert resolution["production_mesh_convergence_certified"] is False


def test_conditional_profile_board_covers_every_frozen_width_and_scenario():
    profiles = _load(AUDIT)["profiles"]
    widths = {80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 320.0}
    scenarios = {
        "ion_low_tail_0p0",
        "ion_low_tail_0p65",
        "ion_high_tail_0p0",
        "ion_high_tail_0p65",
    }
    rates = {34.125, 43.46666666666667}

    assert len(profiles) == len(widths) * len(scenarios) * len(rates)
    board = {
        (
            row["width_nm"],
            row["transport_scenario"]["name"],
            row["blanket_rate_nm_min"],
        )
        for row in profiles
    }
    assert board == {
        (width, scenario, rate)
        for width in widths
        for scenario in scenarios
        for rate in rates
    }


def test_conditional_profile_trajectory_caches_are_checksum_and_spec_bound():
    audit = _load(AUDIT)
    receipts = audit["trajectory_cache_receipts"]

    assert len(receipts) == 28
    for receipt in receipts:
        path = ROOT / receipt["path"]
        assert path.is_file()
        assert sha256(path.read_bytes()).hexdigest() == receipt["sha256"]
        trajectory = _load(path)
        assert trajectory["schema"] == (
            "petch.zhu-npg80-conditional-profile-trajectory.v1"
        )
        assert trajectory["job_spec"]["model_revision"] == (
            "external-union-active-band-dose-factorization-v1"
        )
        assert len(trajectory["profiles"]) == 2


def test_conditional_profiles_conserve_and_keep_surface_scale_nonpredictive():
    profiles = _load(AUDIT)["profiles"]

    for row in profiles:
        assert row["dose_equivalence_exact_for_declared_surface_law"] is True
        assert row["maximum_transport_relative_particle_balance_error"] < 2e-12
        assert (
            row["maximum_state_remap_relative_conservation_residual"] < 1e-12
        )
        validity = row["validity"]
        assert validity["within_declared_scope"] is True
        assert validity["parameter_evidence_supports_prediction"] is False
        assert set(validity["nonpredictive_parameters"]) == {
            "blanket_removal_velocity_m_s",
            "bulk_material_unit_density_m3",
        }


def test_conditional_profile_clearance_and_geometry_claims_are_separated():
    profiles = _load(AUDIT)["profiles"]

    for row in profiles:
        profile = row["profile"]
        assert 0.0 <= profile["etched_depth_nm"] <= 700.0
        assert profile["vertical_relief_grid_resolved"] is True
        expected_lateral_resolution = row["width_nm"] not in {80.0, 320.0}
        assert profile["lateral_feature_grid_resolved"] is (
            expected_lateral_resolution
        )
        assert profile["cd_metrics_grid_resolved"] is expected_lateral_resolution
        for key in ("top_cd_nm", "middle_cd_nm", "bottom_cd_nm"):
            assert np.isfinite(profile[key])
            assert 0.0 < profile[key] <= 400.0
        assert np.isfinite(profile["bow_nm"])
        assert 0.0 < profile["sidewall_angle_from_wafer_deg"] <= 90.0
        maximum_asymmetry = max(
            abs(section["width_x_nm"] - section["width_y_nm"])
            for section in profile["cross_section"]
        )
        if profile["cd_metrics_grid_resolved"]:
            assert maximum_asymmetry < 3.0
        else:
            assert maximum_asymmetry < row["mesh_spacing_nm"]

        if row["tio2_clearance_detected"]:
            assert profile["etched_depth_nm"] == 700.0
            assert row["post_clearance_profile_identified"] is False
            assert row["profile_geometry_status"].startswith(
                "last_pre_clearance_geometry"
            )
            lower, upper = row["clearance_time_bracket_s"]
            assert 0.0 < lower < upper <= 1200.0
        else:
            assert profile["etched_depth_nm"] < 700.0
            assert row["post_clearance_profile_identified"] is True
            assert row["profile_geometry_status"] == "endpoint_geometry"
            assert row["clearance_time_bracket_s"] is None
            assert np.isclose(
                row["accepted_process_equivalent_duration_s"], 1200.0
            )
