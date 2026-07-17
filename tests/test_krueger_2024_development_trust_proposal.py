import importlib.util
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest


_SCRIPT = (
    Path(__file__).parents[1] / "scripts" /
    "krueger_2024_development_trust_proposal.py")
_SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_development_trust_proposal", _SCRIPT)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
derive = _MODULE.derive


def _content_sha(payload, field="proposal_sha256"):
    clean = dict(payload)
    clean.pop(field, None)
    return sha256(json.dumps(
        clean, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _write_inputs(tmp_path, *, endpoint_parameters=(0.4, 0.8), endpoint_response=(30.0, 900.0)):
    current_10 = np.asarray([0.3, 0.9])
    prior_update = np.asarray(endpoint_parameters) - current_10
    model = {
        "schema": "petch.krueger-2024.base-axisymmetric-secant.v1",
        "protocol_id": "K24-PETCH-R1.6",
        "protocol_sha256": "old-protocol-checksum",
        "base_target_table_sha256": _MODULE._sha(_MODULE.BASE_TARGETS),
        "held_out_profile_data_read": False,
        "current_axisymmetric_endpoint": {
            "fraction": float(current_10[0]),
            "oxide_etch_yield_scale": float(current_10[1]),
        },
        "derivation": {
            "parameter_order": list(_MODULE.PARAMETERS),
            "response_order": list(_MODULE.RESPONSES),
            "jacobian": [[100.0, 0.0], [0.0, 1000.0]],
            "condition_number": 10.0,
            "parameter_update": prior_update.tolist(),
        },
        "proposed_configuration": {
            _MODULE.PARAMETERS[0]: endpoint_parameters[0],
            _MODULE.PARAMETERS[1]: endpoint_parameters[1],
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
        },
    }
    model["proposal_sha256"] = _content_sha(model)
    model_path = tmp_path / "response.json"
    model_path.write_text(json.dumps(model), encoding="utf-8")

    endpoint = {
        "status": "complete",
        "config_hash": "a" * 64,
        "configuration": {
            "boundary_case": "base",
            "duration_s": 60.0,
            "dx_um": 0.005,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "topology_change_policy": "continue_gas_cavity",
            _MODULE.PARAMETERS[0]: endpoint_parameters[0],
            _MODULE.PARAMETERS[1]: endpoint_parameters[1],
        },
        "final_metrics": dict(zip(_MODULE.RESPONSES, endpoint_response)),
    }
    endpoint_path = tmp_path / "audit.json"
    endpoint_path.write_text(json.dumps(endpoint), encoding="utf-8")
    receipt = {
        "schema": "petch.krueger-2024.mixed-operator-topology-continuation.v1",
        "scientific_status": "development evidence only; synthetic test",
        "held_out_profile_data_read": False,
        "copied_file_sha256": {
            endpoint_path.name: sha256(endpoint_path.read_bytes()).hexdigest()},
        "final_metrics": dict(zip(_MODULE.RESPONSES, endpoint_response)),
    }
    receipt_path = tmp_path / "receipt.json"
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    return model_path, endpoint_path, receipt_path


def test_development_proposal_separates_full_and_direction_preserving_clipped_steps(tmp_path):
    model, endpoint, receipt = _write_inputs(tmp_path)

    result = derive(model, endpoint, receipt)

    # Full model step is [+0.15, -0.075].  Previous evaluated magnitudes are [0.1, 0.1],
    # so the first coordinate clips the entire direction by 2/3.
    assert result["derivation"]["full_newton_secant_step"] == pytest.approx({
        _MODULE.PARAMETERS[0]: 0.15,
        _MODULE.PARAMETERS[1]: -0.075,
    })
    assert result["trust_region"]["applied_direction_scale"] == pytest.approx(2.0 / 3.0)
    assert result["trust_region"]["clipped"] is True
    assert result["derivation"]["safeguarded_step"] == pytest.approx({
        _MODULE.PARAMETERS[0]: 0.1,
        _MODULE.PARAMETERS[1]: -0.05,
    })
    assert result["proposed_configuration"][_MODULE.PARAMETERS[0]] == pytest.approx(0.5)
    assert result["proposed_configuration"][_MODULE.PARAMETERS[1]] == pytest.approx(0.75)
    assert result["authority"] is False
    assert result["held_out_profile_data_read"] is False
    assert result["proposal_sha256"] == _content_sha(result)


def test_development_proposal_leaves_full_step_unclipped_inside_last_step_box(tmp_path):
    model, endpoint, receipt = _write_inputs(
        tmp_path, endpoint_response=(40.0, 850.0))

    result = derive(model, endpoint, receipt)

    assert result["trust_region"]["clipped"] is False
    assert result["trust_region"]["applied_direction_scale"] == pytest.approx(1.0)
    assert result["derivation"]["safeguarded_candidate"] == pytest.approx(
        result["derivation"]["full_unclipped_candidate"])
    assert result["derivation"]["full_predicted_response"] == pytest.approx({
        "mask_opening_nm": 45.0,
        "etch_depth_nm": 825.0,
    })


def test_development_proposal_refuses_tampered_response_model(tmp_path):
    model, endpoint, receipt = _write_inputs(tmp_path)
    payload = json.loads(model.read_text(encoding="utf-8"))
    payload["derivation"]["jacobian"][0][0] = 99.0
    model.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="proposal_sha256"):
        derive(model, endpoint, receipt)


def test_development_proposal_refuses_unbound_or_mismatched_endpoint(tmp_path):
    model, endpoint, receipt = _write_inputs(tmp_path)
    endpoint_payload = json.loads(endpoint.read_text(encoding="utf-8"))
    endpoint_payload["configuration"][_MODULE.PARAMETERS[0]] = 0.41
    endpoint.write_text(json.dumps(endpoint_payload), encoding="utf-8")

    with pytest.raises(ValueError, match="receipt does not bind"):
        derive(model, endpoint, receipt)


def test_development_proposal_refuses_authority_like_receipt(tmp_path):
    model, endpoint, receipt = _write_inputs(tmp_path)
    payload = json.loads(receipt.read_text(encoding="utf-8"))
    payload["scientific_status"] = "authoritative validation endpoint"
    receipt.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="development-only"):
        derive(model, endpoint, receipt)


def test_development_proposal_refuses_when_fine_cell_gate_already_passes(tmp_path):
    model, endpoint, receipt = _write_inputs(
        tmp_path, endpoint_response=(44.0, 829.0))

    with pytest.raises(ValueError, match="no proposal earned"):
        derive(model, endpoint, receipt)
