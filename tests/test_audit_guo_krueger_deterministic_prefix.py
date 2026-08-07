import json
from pathlib import Path

import pytest

from scripts.audit_guo_krueger_deterministic_prefix import build_report


ROOT = Path(__file__).resolve().parents[1]
RESULTS = (
    ROOT
    / "results"
    / "curated"
    / "guo_krueger_deterministic_prefix"
)


def _audit(n_steps, depths, *, dx_um=0.01):
    times = [0.0, 0.0625, 0.125, 0.25, 0.5]
    history = []
    for index, (time, depth) in enumerate(zip(times, depths)):
        history.append(
            {
                "step": index,
                "physical_time_s": time,
                "metrics": {
                    "etch_depth_nm": depth,
                    "mask_opening_nm": 90.0 - 8.0 * time,
                    "asymmetry_cell_count": 0.0,
                    "mirrored_node_sign_mismatch_pair_count": 0.0,
                    "mirrored_material_label_mismatch_pair_count": 0.0,
                    "maximum_subcell_interface_asymmetry_cells": 0.003,
                },
                "maximum_material_ledger_residual_units_m2": 0.0,
                "maximum_neutral_radiosity_relative_balance_error": 1.0e-12,
                "max_velocity_m_s": 1.0e-8,
                "raw_maximum_face_velocity_m_s": 1.0e-8,
                "rejected_trials": [],
            }
        )
    return {
        "status": "complete",
        "config_hash": f"synthetic-{n_steps}-{dx_um}",
        "configuration": {
            "duration_s": 0.5,
            "n_steps": n_steps,
            "dx_um": dx_um,
            "maximum_accepted_steps": n_steps + 8,
            "radiosity_backend": "deterministic_extruded_2d",
        },
        "history": history,
        "topology_events": [],
        "extrusion_projection_max_deviation_mesh_units": 1.0e-15,
    }


def _write(directory: Path, audit: dict):
    directory.mkdir()
    (directory / "audit.json").write_text(
        json.dumps(audit), encoding="utf-8"
    )


def test_early_relative_depth_gate_is_not_hidden_by_endpoint_pass(tmp_path):
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    _write(coarse, _audit(8, [0.0, 0.80, 1.4, 3.0, 6.10]))
    _write(fine, _audit(16, [0.0, 0.65, 1.4, 3.0, 6.20]))

    report = build_report(coarse, fine, "time")

    assert report["gates"]["terminal_depth_abs_relative"]["passed"]
    assert not report["gates"][
        "maximum_matched_depth_abs_relative"
    ]["passed"]
    assert not report["all_gates_passed"]


def test_identical_physics_except_time_refinement_passes(tmp_path):
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    _write(coarse, _audit(8, [0.0, 0.65, 1.4, 3.0, 6.10]))
    _write(fine, _audit(16, [0.0, 0.65, 1.4, 3.0, 6.20]))

    report = build_report(coarse, fine, "time")

    assert report["all_gates_passed"]
    assert report["health"]["maximum_raw_to_resolved_speed_ratio"] == 1.0


def test_unauthorized_configuration_change_is_refused(tmp_path):
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    reference = _audit(8, [0.0, 0.65, 1.4, 3.0, 6.10])
    refined = _audit(16, [0.0, 0.65, 1.4, 3.0, 6.20])
    reference["configuration"]["seed"] = 241
    refined["configuration"]["seed"] = 242
    _write(coarse, reference)
    _write(fine, refined)

    with pytest.raises(ValueError, match="unauthorized configuration"):
        build_report(coarse, fine, "time")


def test_subcell_asymmetry_is_reported_but_not_mislabeled_as_cell_failure(
        tmp_path):
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    _write(coarse, _audit(8, [0.0, 0.65, 1.4, 3.0, 6.10]))
    _write(fine, _audit(16, [0.0, 0.65, 1.4, 3.0, 6.20]))

    report = build_report(coarse, fine, "time")

    assert report["all_gates_passed"]
    assert (
        report["health"]["maximum_subcell_interface_asymmetry_cells"]
        == 0.003
    )
    assert "maximum_subcell_interface_asymmetry_cells" not in report["gates"]


@pytest.mark.parametrize(
    "metric",
    (
        "asymmetry_cell_count",
        "mirrored_node_sign_mismatch_pair_count",
        "mirrored_material_label_mismatch_pair_count",
    ),
)
def test_each_resolved_symmetry_disagreement_fails_closed(tmp_path, metric):
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    reference = _audit(8, [0.0, 0.65, 1.4, 3.0, 6.10])
    refined = _audit(16, [0.0, 0.65, 1.4, 3.0, 6.20])
    refined["history"][-1]["metrics"][metric] = 1
    _write(coarse, reference)
    _write(fine, refined)

    report = build_report(coarse, fine, "time")

    gate = {
        "asymmetry_cell_count": "maximum_asymmetric_cell_count",
        "mirrored_node_sign_mismatch_pair_count": (
            "maximum_mirrored_node_sign_mismatch_pair_count"
        ),
        "mirrored_material_label_mismatch_pair_count": (
            "maximum_mirrored_material_label_mismatch_pair_count"
        ),
    }[metric]
    assert not report["gates"][gate]["passed"]
    assert not report["all_gates_passed"]


def test_committed_time_ladder_preserves_failures_before_final_pass():
    original = json.loads((RESULTS / "time_gate.json").read_text())
    refined = json.loads((RESULTS / "time_gate_refined.json").read_text())
    final = json.loads((RESULTS / "time_gate_final.json").read_text())

    assert not original["all_gates_passed"]
    assert not refined["all_gates_passed"]
    assert final["all_gates_passed"]
    assert original["gates"][
        "maximum_matched_depth_abs_relative"
    ]["value"] > refined["gates"][
        "maximum_matched_depth_abs_relative"
    ]["value"] > final["gates"][
        "maximum_matched_depth_abs_relative"
    ]["value"]
    assert final["reference"]["nominal_step_s"] == 0.015625
    assert final["refined"]["nominal_step_s"] == 0.0078125
