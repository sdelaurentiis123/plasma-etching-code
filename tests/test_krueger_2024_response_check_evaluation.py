import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "krueger_2024_response_check_evaluation.py"
SPEC = importlib.util.spec_from_file_location("krueger_response_evaluation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _launch():
    return {
        "schema": "petch.krueger-2024.r19-response-check-launch.v1",
        "protocol_id": "K24-PETCH-R1.9",
        "authority": False,
        "held_out_profile_data_read": False,
        "ten_nm_sequence_closes_after_this_evaluation": True,
        "candidate": {
            "effective_mask_crosslinked_growth_fraction": 0.9004722559883319,
            "oxide_etch_yield_scale": 0.5586489665864749,
        },
        "target_nm": {"mask_opening": 45.0, "etch_depth": 825.0},
        "committed_predictions_nm": {
            "same_operator_local_response": {
                "mask_opening": 44.1514, "etch_depth": 822.9295}},
        "committed_model_error_gate_nm": {
            "mask_opening": 2.153, "etch_depth": 12.905},
    }


def _audit(*, opening=45.2, depth=824.0, status="complete", ledger=0.0):
    return {
        "status": status,
        "config_hash": "manufactured",
        "configuration": {
            "boundary_case": "base",
            "duration_s": 60.0,
            "dx_um": 0.01,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "topology_change_policy": "continue_gas_cavity",
            "effective_mask_crosslinked_growth_fraction": 0.9004722559883319,
            "oxide_etch_yield_scale": 0.5586489665864749,
        },
        "history": [{
            "physical_time_s": 60.0,
            "metrics": {"mask_opening_nm": opening, "etch_depth_nm": depth},
            "maximum_material_ledger_residual_units_m2": ledger,
            "maximum_neutral_radiosity_relative_balance_error": 2e-12,
            "validity": {"within_declared_scope": True},
        }],
        "final_metrics": {"mask_opening_nm": opening, "etch_depth_nm": depth},
        "terminal_event": None,
    }


def _write(tmp_path, name, payload):
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_strong_response_pass_closes_ten_nm_sequence(tmp_path):
    result = MODULE.evaluate(
        _write(tmp_path, "launch.json", _launch()),
        _write(tmp_path, "audit.json", _audit()))

    assert result["decision"] == "strong_response_model_pass"
    assert result["ten_nm_sequence_closed"] is True
    assert result["authority"] is False
    assert result["held_out_profile_data_read"] is False
    assert result["response_gates"]["target_within_one_fine_cell"] is True


def test_numerical_contract_failure_cannot_pass_response(tmp_path):
    result = MODULE.evaluate(
        _write(tmp_path, "launch.json", _launch()),
        _write(tmp_path, "audit.json", _audit(ledger=1e-15)))

    assert result["decision"] == "reject_numerical_contract"
    assert result["ten_nm_sequence_closed"] is True


def test_changed_candidate_is_refused(tmp_path):
    launch = _launch()
    launch["candidate"]["oxide_etch_yield_scale"] = 0.6

    try:
        MODULE.evaluate(
            _write(tmp_path, "launch.json", launch),
            _write(tmp_path, "audit.json", _audit()))
    except ValueError as error:
        assert "sealed field" in str(error)
    else:
        raise AssertionError("changed candidate was not refused")
