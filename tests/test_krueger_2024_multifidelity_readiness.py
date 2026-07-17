import importlib.util
from hashlib import sha256
import json
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).parents[1] / "scripts" /
    "krueger_2024_multifidelity_readiness.py")
_SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_multifidelity_readiness", _SCRIPT)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


def _write(path, payload):
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _content_sha(payload, field):
    clean = dict(payload)
    clean.pop(field, None)
    return sha256(json.dumps(
        clean, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _inputs(tmp_path):
    target_sha = _MODULE._sha(_MODULE.TARGETS)
    candidate = {
        "effective_mask_crosslinked_growth_fraction": 0.9,
        "oxide_etch_yield_scale": 0.55,
    }
    launch = {
        "schema": "petch.krueger-2024.r19-response-check-launch.v1",
        "protocol_id": "K24-PETCH-R1.9",
        "authority": False,
        "held_out_profile_data_read": False,
        "ten_nm_sequence_closes_after_this_evaluation": True,
        "candidate": candidate,
        "target_nm": {"mask_opening": 45.0, "etch_depth": 825.0},
        "base_only_input_sha256": {_MODULE.TARGETS.name: target_sha},
        "executable_source_sha256": {
            "src/petch/feature_step_3d.py": "old-feature",
            "src/petch/material_mechanism_3d.py": "same-material",
        },
    }
    launch_path = _write(tmp_path / "launch.json", launch)
    run = {
        "status": "complete",
        "configuration": {
            "boundary_case": "base", "duration_s": 60.0, "dx_um": 0.01,
            **candidate,
        },
    }
    run_path = _write(tmp_path / "run.json", run)
    evaluation = {
        "schema": "petch.krueger-2024.r19-response-check-evaluation.v1",
        "held_out_profile_data_read": False,
        "ten_nm_sequence_closed": True,
        "decision": "reject_response_model",
        "trajectory_contract": {"pass": True},
        "actual_nm": {"mask_opening": 45.1, "etch_depth": 853.2},
        "target_error_nm": {"mask_opening": 0.1, "etch_depth": 28.2},
        "inputs": {
            "launch_manifest": {"sha256": _MODULE._sha(launch_path)},
            "run_audit": {"sha256": _MODULE._sha(run_path)},
        },
    }
    evaluation["evaluation_sha256"] = _content_sha(evaluation, "evaluation_sha256")
    evaluation_path = _write(tmp_path / "evaluation.json", evaluation)
    multires = {
        "paired_10nm_vs_5nm": {"initial": {
            "etch_depth_nm_s": {"relative_to_fine": -0.001},
            "mask_opening_nm_s": {"relative_to_fine": -0.01},
            "maximum_feature_width_nm_s": {"relative_to_fine": 0.7},
        }},
    }
    multires_path = _write(tmp_path / "multires.json", multires)
    summary = {
        "schema": "petch.krueger_2024_cuda_profile_summary.v1",
        "held_out_profile_data_read": False,
        "calibration_performed": False,
        "inputs": {
            "paired_10nm_5nm_initial_audit": {
                "sha256": _MODULE._sha(multires_path)}},
        "operator_receipts": {"surface_state_remap_backend": "legacy_knn"},
    }
    summary["summary_sha256"] = _content_sha(summary, "summary_sha256")
    summary_path = _write(tmp_path / "summary.json", summary)
    return launch_path, evaluation_path, run_path, multires_path, summary_path


def test_readiness_blocks_stale_epoch_and_short_rate_substitution(tmp_path):
    inputs = _inputs(tmp_path)
    result = _MODULE.derive(
        *inputs, current_revision="abc123", current_sources={
            "src/petch/feature_step_3d.py": "new-feature",
            "src/petch/material_mechanism_3d.py": "same-material",
            "src/petch/surface_mesh_3d.py": "new-surface",
        })

    assert result["status"] == "blocked_before_parameter_proposal"
    assert result["held_out_profile_data_read"] is False
    codes = {item["code"] for item in result["blockers"]}
    assert "remap_backend_not_selected" in codes
    assert "r19_response_belongs_to_prior_operator_epoch" in codes
    assert "no_current_epoch_high_fidelity_endpoint_anchor" in codes
    assert "no_paired_endpoint_discrepancy" in codes
    assert result["current_center"]["parameters"][
        "effective_mask_crosslinked_growth_fraction"] == pytest.approx(0.9)
    assert result["readiness_sha256"] == _content_sha(result, "readiness_sha256")


def test_readiness_refuses_tampered_evaluation(tmp_path):
    inputs = list(_inputs(tmp_path))
    evaluation = json.loads(inputs[1].read_text(encoding="utf-8"))
    evaluation["actual_nm"]["etch_depth"] = 800.0
    inputs[1].write_text(json.dumps(evaluation), encoding="utf-8")

    with pytest.raises(ValueError, match="valid closed receipt"):
        _MODULE.derive(
            *inputs, current_revision="abc123",
            current_sources={"src/petch/feature_step_3d.py": "new"})
