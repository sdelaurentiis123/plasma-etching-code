import json
import importlib.util
from hashlib import sha256
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).parents[1] / "scripts" / "krueger_2024_calibration.py"
_SPEC = importlib.util.spec_from_file_location("krueger_2024_calibration", _SCRIPT)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
propose = _MODULE.propose

_CORRECTION_SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "krueger_2024_coupled_correction.py")
_CORRECTION_SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_coupled_correction", _CORRECTION_SCRIPT)
_CORRECTION_MODULE = importlib.util.module_from_spec(_CORRECTION_SPEC)
_CORRECTION_SPEC.loader.exec_module(_CORRECTION_MODULE)
derive_coupled_correction = _CORRECTION_MODULE.derive
derive_path_correction = _CORRECTION_MODULE.derive_path_resolved
derive_axisymmetric_secant = _CORRECTION_MODULE.derive_axisymmetric_secant

_GRID_SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "krueger_2024_grid_correction.py")
_GRID_SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_grid_correction", _GRID_SCRIPT)
_GRID_MODULE = importlib.util.module_from_spec(_GRID_SPEC)
_GRID_SPEC.loader.exec_module(_GRID_MODULE)
derive_grid_correction = _GRID_MODULE.derive


def _audit(path, fraction, *, opening=None, clogged=False):
    configuration = {
        "boundary_case": "base",
        "charging": "disabled_for_Krueger_2024_calibration_and_transfer",
        "duration_s": 60.0,
        "dx_um": 0.01,
        "geometry": {
            "substrate_top_um": 1.8,
            "domain_height_um": 2.8,
        },
    }
    if fraction:
        configuration["effective_mask_crosslinked_growth_fraction"] = fraction
    payload = {
        "configuration": configuration,
        "status": "terminal_feature_clogged" if clogged else "complete",
        "history": [{"physical_time_s": 60.0}],
        "final_metrics": {"mask_opening_nm": opening},
        "terminal_event": ({"kind": "feature_clogged"} if clogged else None),
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_krueger_calibration_proposal_accepts_closure_as_censored_zero(tmp_path):
    zero = _audit(tmp_path / "zero.json", 0.0, opening=12.0, clogged=True)
    one = _audit(tmp_path / "one.json", 1.0, opening=60.0)

    result = propose(zero, one)

    assert result["endpoints"][0]["mask_opening_nm"] == 0.0
    assert result["endpoints"][0]["endpoint_kind"] == "resolved_closure_event"
    assert result["interpolation"][
        "proposed_effective_mask_crosslinked_growth_fraction"] == 0.75
    assert result["held_out_profile_data_read"] is False


def test_krueger_calibration_proposal_refuses_unbracketed_endpoint(tmp_path):
    zero = _audit(tmp_path / "zero.json", 0.0, opening=10.0)
    one = _audit(tmp_path / "one.json", 1.0, opening=20.0)

    with pytest.raises(ValueError, match="do not bracket"):
        propose(zero, one)


def _complete_coupled_audit(
        path, fraction, yield_scale, opening, depth, *, axisymmetric=False, dx=0.01):
    configuration = {
        "boundary_case": "base",
        "duration_s": 60.0,
        "dx_um": dx,
        "effective_mask_crosslinked_growth_fraction": fraction,
        "oxide_etch_yield_scale": yield_scale,
    }
    if axisymmetric:
        configuration.update({
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
        })
    payload = {
        "config_hash": "a" * 64,
        "configuration": configuration,
        "status": "complete",
        "final_metrics": {
            "mask_opening_nm": opening,
            "etch_depth_nm": depth,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_krueger_coupled_correction_uses_only_base_response_and_certified_azimuth(
        tmp_path):
    endpoints = [
        _complete_coupled_audit(tmp_path / "f0.json", 0.0, 1.0, 13.0, 981.0),
        _complete_coupled_audit(tmp_path / "fm.json", 0.7, 1.0, 15.0, 1191.0),
        _complete_coupled_audit(tmp_path / "f1.json", 1.0, 1.0, 59.0, 1421.0),
    ]
    coupled = _complete_coupled_audit(
        tmp_path / "coupled.json", 0.90, 0.61, 40.0, 900.0)
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
        "variants": {
            "single_published_plane": {
                "material_profile_rate": {
                    "2": {"net_volume_rate_mesh_units3_s": -1.0}}},
            "axisymmetric_uniform_16": {
                "material_profile_rate": {
                    "2": {"net_volume_rate_mesh_units3_s": -1.06}}},
        },
    }), encoding="utf-8")

    result = derive_coupled_correction(endpoints, coupled, azimuth)

    proposed = result["proposed_configuration"]
    assert 0.90 < proposed["effective_mask_crosslinked_growth_fraction"] < 1.0
    assert 0.0 < proposed["oxide_etch_yield_scale"] < 0.61
    assert proposed["ion_azimuthal_order"] == 16
    assert result["held_out_profile_data_read"] is False


def test_krueger_coupled_correction_refuses_uncertified_azimuth(tmp_path):
    endpoints = [
        _complete_coupled_audit(tmp_path / "f0.json", 0.0, 1.0, 13.0, 981.0),
        _complete_coupled_audit(tmp_path / "fm.json", 0.7, 1.0, 15.0, 1191.0),
        _complete_coupled_audit(tmp_path / "f1.json", 1.0, 1.0, 59.0, 1421.0),
    ]
    coupled = _complete_coupled_audit(
        tmp_path / "coupled.json", 0.90, 0.61, 40.0, 900.0)
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": False,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")

    with pytest.raises(ValueError, match="has not certified"):
        derive_coupled_correction(endpoints, coupled, azimuth)


def test_krueger_path_correction_solves_bounded_two_by_two_base_secant(tmp_path):
    endpoints = [
        _complete_coupled_audit(tmp_path / "f0.json", 0.0, 1.0, 13.0, 981.0),
        _complete_coupled_audit(tmp_path / "fm.json", 0.7, 1.0, 15.0, 1191.0),
        _complete_coupled_audit(tmp_path / "f1.json", 1.0, 1.0, 59.0, 1421.0),
    ]
    plane = _complete_coupled_audit(
        tmp_path / "plane.json", 0.90, 0.61, 40.0, 900.0)
    axisymmetric = _complete_coupled_audit(
        tmp_path / "axisymmetric.json", 0.934, 0.519, 50.3, 750.2,
        axisymmetric=True)
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")

    result = derive_path_correction(endpoints, plane, axisymmetric, azimuth)

    proposed = result["proposed_configuration"]
    assert 0.90 < proposed["effective_mask_crosslinked_growth_fraction"] < 0.934
    assert 0.519 < proposed["oxide_etch_yield_scale"] < 0.61
    assert result["derivation"]["condition_number"] < 1e6
    assert result["held_out_profile_data_read"] is False


def test_krueger_final_axisymmetric_secant_uses_two_paths_and_no_heldout(tmp_path):
    endpoints = [
        _complete_coupled_audit(tmp_path / "f0.json", 0.0, 1.0, 13.0, 981.0),
        _complete_coupled_audit(tmp_path / "fm.json", 0.7, 1.0, 15.0, 1191.0),
        _complete_coupled_audit(tmp_path / "f1.json", 1.0, 1.0, 59.0, 1421.0),
    ]
    previous = _complete_coupled_audit(
        tmp_path / "previous.json", 0.934, 0.519, 50.3, 750.2,
        axisymmetric=True)
    current = _complete_coupled_audit(
        tmp_path / "current.json", 0.921, 0.575, 49.9, 853.3,
        axisymmetric=True)
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")

    result = derive_axisymmetric_secant(endpoints, previous, current, azimuth)

    proposed = result["proposed_configuration"]
    assert 0.85 < proposed["effective_mask_crosslinked_growth_fraction"] < 0.921
    assert 0.50 < proposed["oxide_etch_yield_scale"] < 0.575
    assert result["derivation"]["condition_number"] < 100.0
    assert result["held_out_profile_data_read"] is False


def test_krueger_final_axisymmetric_secant_refuses_plane_endpoint(tmp_path):
    endpoints = [
        _complete_coupled_audit(tmp_path / "f0.json", 0.0, 1.0, 13.0, 981.0),
        _complete_coupled_audit(tmp_path / "fm.json", 0.7, 1.0, 15.0, 1191.0),
        _complete_coupled_audit(tmp_path / "f1.json", 1.0, 1.0, 59.0, 1421.0),
    ]
    previous = _complete_coupled_audit(
        tmp_path / "previous.json", 0.934, 0.519, 50.3, 750.2)
    current = _complete_coupled_audit(
        tmp_path / "current.json", 0.921, 0.575, 49.9, 853.3,
        axisymmetric=True)
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")

    with pytest.raises(ValueError, match="certified 16-node"):
        derive_axisymmetric_secant(endpoints, previous, current, azimuth)


def _final_derivation(path, fraction=0.9, yield_scale=0.5):
    payload = {
        "schema": "petch.krueger-2024.base-axisymmetric-secant.v1",
        "protocol_sha256": _GRID_MODULE._sha(_GRID_MODULE.PROTOCOL),
        "held_out_profile_data_read": False,
        "proposed_configuration": {
            "effective_mask_crosslinked_growth_fraction": fraction,
            "oxide_etch_yield_scale": yield_scale,
        },
        "derivation": {"jacobian": [[100.0, 0.0], [0.0, 1000.0]]},
    }
    payload["proposal_sha256"] = sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_krueger_grid_correction_reuses_frozen_response_and_base_only(tmp_path):
    derivation = _final_derivation(tmp_path / "derivation.json")
    base10 = _complete_coupled_audit(
        tmp_path / "base10.json", 0.9, 0.5, 43.0, 838.0,
        axisymmetric=True)
    base5 = _complete_coupled_audit(
        tmp_path / "base5.json", 0.9, 0.5, 55.0, 850.0,
        axisymmetric=True, dx=0.005)

    result = derive_grid_correction(derivation, base10, base5)

    proposed = result["proposed_configuration"]
    assert proposed["effective_mask_crosslinked_growth_fraction"] == pytest.approx(0.8)
    assert proposed["oxide_etch_yield_scale"] == pytest.approx(0.475)
    assert result["held_out_profile_data_read"] is False


def test_krueger_grid_correction_refuses_when_refinement_already_passes(tmp_path):
    derivation = _final_derivation(tmp_path / "derivation.json")
    base10 = _complete_coupled_audit(
        tmp_path / "base10.json", 0.9, 0.5, 43.0, 838.0,
        axisymmetric=True)
    base5 = _complete_coupled_audit(
        tmp_path / "base5.json", 0.9, 0.5, 46.0, 828.0,
        axisymmetric=True, dx=0.005)

    with pytest.raises(ValueError, match="already passes"):
        derive_grid_correction(derivation, base10, base5)
