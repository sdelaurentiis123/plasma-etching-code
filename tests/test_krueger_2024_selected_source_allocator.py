import importlib.util
from pathlib import Path

import numpy as np
import pytest

from petch.neutral_radiosity_3d import DiffuseFormFactors3D


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_selected_source_allocator",
    ROOT / "scripts" / "krueger_2024_selected_source_allocator.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _factors(rows, escape, rays=16):
    source = []
    target = []
    fraction = []
    for left, values in enumerate(rows):
        for right, value in values.items():
            source.append(left)
            target.append(right)
            fraction.append(value)
    return DiffuseFormFactors3D(
        len(rows), np.asarray(source), np.asarray(target),
        np.asarray(fraction), np.asarray(escape), rays)


def test_exact_radiosity_source_row_decomposition_closes_signed_change():
    coarse = _factors(
        ({1: 0.20}, {2: 0.30}, {0: 0.10}), (0.80, 0.70, 0.90), 8)
    fine = _factors(
        ({1: 0.30}, {2: 0.20}, {0: 0.15}), (0.70, 0.80, 0.85), 16)
    weights = np.asarray([
        [1.0, 0.0, 0.0],
        [0.0, 0.25, 0.75],
    ])
    result = MODULE.exact_radiosity_source_row_decomposition(
        coarse, fine, np.asarray([1.0, 2.0, 3.0]),
        np.asarray([4.0, 2.0, 1.0]), np.asarray([0.2, 0.4, 0.6]),
        weights)
    np.testing.assert_allclose(
        np.sum(result["contributions"], axis=1),
        weights @ (
            result["fine_incident_rate"] - result["coarse_incident_rate"]),
        rtol=2e-12, atol=2e-12)
    assert np.max(result["relative_closure_error"]) < 2e-12


def test_source_row_ranking_selects_concentrated_causal_row():
    authority = np.asarray([[8.0, 1.0, 0.5, 0.1], [4.0, 0.5, 0.2, 0.1]])
    replicates = np.repeat(authority[None, :, :], 8, axis=0)
    ranking = MODULE.rank_source_rows(
        authority, replicates, np.ones(2), target_concentration=0.75,
        maximum_selected_fraction=0.5)
    assert ranking["selected_source_faces"].tolist() == [0]
    assert ranking["target_reached_within_cap"]
    assert ranking["selected_score_fraction"] > 0.75


def test_source_row_ranking_reports_diffuse_cap_blocker():
    authority = np.ones((1, 20))
    replicates = np.repeat(authority[None, :, :], 8, axis=0)
    ranking = MODULE.rank_source_rows(
        authority, replicates, np.ones(1), target_concentration=0.9,
        maximum_selected_fraction=0.25)
    assert ranking["selected_face_count"] == 5
    assert not ranking["target_reached_within_cap"]
    assert ranking["selected_score_fraction"] == pytest.approx(0.25)


def test_stage_a_contract_refuses_a_passing_or_heldout_artifact(tmp_path):
    base = {
        "schema": MODULE.closure.SCHEMA,
        "stage": "stage_a",
        "status": "bounded_precision_hold",
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
        },
        "operator": MODULE.closure.OPERATOR,
        "sampling": {
            "ray_levels": MODULE.closure.RAY_LEVELS,
            "replicate_seeds": MODULE.closure.REPLICATE_SEEDS,
        },
        "stage_a": {
            "all_gates_pass": False,
            "gates": {
                "exact_nested_sampling_extension": True,
                "exact_replicate_count": True,
                "row_closure": True,
            },
        },
        "checkpoint": {"checkpoint_sha256": "abc"},
    }
    path = tmp_path / "audit.json"
    MODULE._write_json_atomic(path, base)
    assert MODULE._validate_stage_a_contract(base, path)[
        "checkpoint_sha256"] == "abc"
    passing = {**base, "status": MODULE.closure.STAGE_A_PASS_STATUS}
    with pytest.raises(ValueError, match="bounded Stage-A hold"):
        MODULE._validate_stage_a_contract(passing, path)
    heldout = {**base, "data_firewall": {
        **base["data_firewall"], "held_out_observations_loaded": True}}
    with pytest.raises(ValueError, match="firewall"):
        MODULE._validate_stage_a_contract(heldout, path)


def test_cli_controls_refuse_global_or_unbounded_escalation():
    args = MODULE.parse_args([])
    assert args.maximum_selected_fraction == 0.25
    assert args.maximum_wall_s == 120.0
    assert MODULE.BASE_RAY_LEVEL == 16
    assert MODULE.SELECTED_RAY_LEVEL == 32
    assert MODULE.HORIZON_DIVISOR == 1024
    with pytest.raises(SystemExit):
        MODULE.parse_args(["--maximum-selected-fraction", "0.5"])
    with pytest.raises(SystemExit):
        MODULE.parse_args(["--maximum-wall-s", "121"])


def test_projected_support_counts_one_periodic_y_cell_at_forty_nm():
    class Scheme:
        patch_key = np.asarray([
            [1, 2, 1, 0, 0, 0],  # y is tangential: expected 40 nm x 20 nm
            [1, 1, 1, 0, 0, 0],  # y is normal: expected 40 nm x 40 nm
        ])
        contribution_patch_index = np.asarray([0, 1])
        contribution_face_index = np.asarray([0, 1])
        contribution_area_m2 = np.asarray([8.0e-16, 16.0e-16])

    support = MODULE._dominant_axis_projected_support_fraction(
        Scheme(), np.asarray([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
        patch_scale_m=40.0e-9, periodic_y_extent_m=20.0e-9)
    np.testing.assert_allclose(support, [1.0, 1.0])
