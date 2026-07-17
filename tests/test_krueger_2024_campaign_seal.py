import importlib.util
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).parents[1]


def _module(name):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FREEZE = _module("krueger_2024_freeze")
CAMPAIGN = _module("krueger_2024_transfer_campaign")
GRID = _module("krueger_2024_grid_correction")


def _base(path, dx):
    payload = {
        "status": "complete",
        "config_hash": "a" * 64,
        "configuration": {
            "boundary_case": "base",
            "duration_s": 60.0,
            "n_steps": 60,
            "dx_um": dx,
            "n_position": 16,
            "neutral_transverse_order": 5,
            "neutral_normal_order": 2,
            "neutral_direction_polar_order": 8,
            "neutral_direction_azimuthal_order": 16,
            "ion_energy_bin_eV": 250.0,
            "ion_angle_bin_deg": 0.25,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "ballistic_transport": "face_gather",
            "ballistic_face_quadrature_points": 3,
            "compressed_boundary_quadrature": True,
            "neutral_speed_quadrature": "analytic_speed_marginal",
            "neutral_tensor_velocity_quadrature_active": False,
            "radiosity_rays_per_face": 8,
            "radiosity_relative_tolerance": 1e-12,
            "radiosity_maximum_iterations": 2000,
            "radiosity_enabled": True,
            "seed": 241,
            "adaptive_profile_timestep": True,
            "minimum_step_duration_s": 0.001,
            "target_displacement_cells": 0.35,
            "maximum_displacement_cells": 0.75,
            "adaptive_shrink_factor": 0.5,
            "adaptive_growth_factor": 1.5,
            "adaptive_safety_factor": 0.9,
            "maximum_accepted_steps": 10000,
            "profile_reinitialization": "cr2",
            "topology_change_policy": "continue_gas_cavity",
            "surface_state_remap_backend": "common_refinement",
            "geometry": {
                "cell_width_um": 0.13,
                "cell_length_um": 0.02,
                "domain_height_um": 2.8,
                "opening_width_um": 0.09,
                "mask_thickness_um": 0.85,
                "substrate_top_um": 1.8,
                "initial_etched_depth_um": 0.0,
            },
            "effective_mask_crosslinked_growth_fraction": 0.92,
            "oxide_etch_yield_scale": 0.575,
        },
        "final_metrics": {
            "mask_opening_nm": 45.0,
            "etch_depth_nm": 825.0,
            "maximum_feature_width_nm": 90.0,
            "remaining_mask_thickness_nm": 850.0,
            "asymmetry_cell_count": 0.0,
        },
        "history": [{
            "maximum_material_ledger_residual_units_m2": 0.0,
            "maximum_neutral_radiosity_relative_balance_error": 1e-12,
            "rejected_trials": [],
            "unresolved_gas_cavity_volume_upper_bound_m3": 0.0,
        }],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _calibration(path):
    payload = {
        "schema": "petch.krueger-2024.base-axisymmetric-secant.v1",
        "protocol_sha256": FREEZE._sha(FREEZE.PROTOCOL),
        "held_out_profile_data_read": False,
        "proposed_configuration": {
            "effective_mask_crosslinked_growth_fraction": 0.92,
            "oxide_etch_yield_scale": 0.575,
        },
        "derivation": {
            "jacobian": [[-100.0, 0.0], [0.0, 1000.0]],
        },
    }
    payload["proposal_sha256"] = sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_krueger_grid_correction_is_r17_base_only_and_checksum_bound(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    fine = json.loads(base5.read_text(encoding="utf-8"))
    fine["final_metrics"].update({
        "mask_opening_nm": 43.0,
        "etch_depth_nm": 850.0,
    })
    base5.write_text(json.dumps(fine), encoding="utf-8")
    calibration = _calibration(tmp_path / "calibration.json")

    result = GRID.derive(calibration, base10, base5)

    assert result["protocol_id"] == "K24-PETCH-R1.9"
    assert result["held_out_profile_data_read"] is False
    assert result["derivation"]["fine_grid_target_error"] == [-2.0, 25.0]
    assert result["proposed_configuration"][
        "effective_mask_crosslinked_growth_fraction"] == pytest.approx(0.90)
    assert result["proposed_configuration"][
        "oxide_etch_yield_scale"] == pytest.approx(0.55)
    GRID._verify_embedded_sha256(result, "proposal_sha256")


def test_krueger_grid_correction_refuses_tampered_calibration(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    fine = json.loads(base5.read_text(encoding="utf-8"))
    fine["final_metrics"]["etch_depth_nm"] = 850.0
    base5.write_text(json.dumps(fine), encoding="utf-8")
    calibration = _calibration(tmp_path / "calibration.json")
    payload = json.loads(calibration.read_text(encoding="utf-8"))
    payload["derivation"]["jacobian"][1][1] = 900.0
    calibration.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="proposal_sha256"):
        GRID.derive(calibration, base10, base5)


def _supporting_evidence(tmp_path):
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")
    charging = tmp_path / "charging.json"
    charging.write_text(json.dumps({
        "schema": "petch.krueger-2024.charging-causality.v1",
        "held_out_profile_data_read": False,
        "conservation": {"maximum_charge_conservation_relative_error": 0.0},
        "paired_exact_hard_visibility_audit": {
            "charged_over_zero_floor_ion_flux": 1.0},
    }), encoding="utf-8")
    return azimuth, charging


def test_krueger_freeze_binds_refinement_and_causal_evidence_without_heldout(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    calibration = _calibration(tmp_path / "calibration.json")
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")
    charging = tmp_path / "charging.json"
    charging.write_text(json.dumps({
        "schema": "petch.krueger-2024.charging-causality.v1",
        "held_out_profile_data_read": False,
        "conservation": {"maximum_charge_conservation_relative_error": 1e-16},
        "paired_exact_hard_visibility_audit": {
            "charged_over_zero_floor_ion_flux": 0.995},
    }), encoding="utf-8")

    result = FREEZE.freeze(base10, base5, calibration, azimuth, charging)

    assert result["frozen_physics"]["ion_azimuthal_order"] == 16
    assert result["base_grid_difference"]["etch_depth_nm"] == 0.0
    assert result["held_out_profile_data_read"] is False
    assert "data/experimental/krueger_2024/transfer_observations.csv" not in (
        result["boundary_data_sha256"])
    assert result["source_sha256"]
    assert result["boundary_data_sha256"]
    assert len(result["reveal_sha256"]) == 64


def test_krueger_freeze_refuses_nonconservative_base_endpoint(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    calibration = _calibration(tmp_path / "calibration.json")
    payload = json.loads(base5.read_text(encoding="utf-8"))
    payload["history"][0]["maximum_material_ledger_residual_units_m2"] = 1e-30
    base5.write_text(json.dumps(payload), encoding="utf-8")
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")
    charging = tmp_path / "charging.json"
    charging.write_text(json.dumps({
        "schema": "petch.krueger-2024.charging-causality.v1",
        "held_out_profile_data_read": False,
        "conservation": {"maximum_charge_conservation_relative_error": 0.0},
        "paired_exact_hard_visibility_audit": {
            "charged_over_zero_floor_ion_flux": 1.0},
    }), encoding="utf-8")

    with pytest.raises(ValueError, match="exactly closed"):
        FREEZE.freeze(base10, base5, calibration, azimuth, charging)


def test_krueger_freeze_refuses_mixed_refinement_operator(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    calibration = _calibration(tmp_path / "calibration.json")
    payload = json.loads(base5.read_text(encoding="utf-8"))
    payload["configuration"]["radiosity_rays_per_face"] = 4
    base5.write_text(json.dumps(payload), encoding="utf-8")
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")
    charging = tmp_path / "charging.json"
    charging.write_text(json.dumps({
        "schema": "petch.krueger-2024.charging-causality.v1",
        "held_out_profile_data_read": False,
        "conservation": {"maximum_charge_conservation_relative_error": 0.0},
        "paired_exact_hard_visibility_audit": {
            "charged_over_zero_floor_ion_flux": 1.0},
    }), encoding="utf-8")

    with pytest.raises(ValueError, match="one numerical operator"):
        FREEZE.freeze(base10, base5, calibration, azimuth, charging)


def test_krueger_freeze_refuses_legacy_or_mixed_remap_operator(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    calibration = _calibration(tmp_path / "calibration.json")
    azimuth, charging = _supporting_evidence(tmp_path)
    payload = json.loads(base5.read_text(encoding="utf-8"))
    payload["configuration"]["surface_state_remap_backend"] = "legacy_knn"
    base5.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="R1.9 operator"):
        FREEZE.freeze(base10, base5, calibration, azimuth, charging)

    payload["configuration"]["surface_state_remap_backend"] = "indexed_knn"
    base5.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="one numerical operator"):
        FREEZE.freeze(base10, base5, calibration, azimuth, charging)


def test_krueger_freeze_accepts_only_monotone_fine_grid_timestep_refinement(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    payload = json.loads(base5.read_text(encoding="utf-8"))
    payload["configuration"]["minimum_step_duration_s"] = 0.00025
    payload["history"][0].update({
        "reassigned_unresolved_material_nodes": 1,
        "unresolved_material_volume_upper_bound_m3": (0.005e-6) ** 3,
    })
    base5.write_text(json.dumps(payload), encoding="utf-8")
    calibration = _calibration(tmp_path / "calibration.json")
    azimuth, charging = _supporting_evidence(tmp_path)

    result = FREEZE.freeze(base10, base5, calibration, azimuth, charging)

    assert result["protocol_id"] == "K24-PETCH-R1.9"
    assert result["authority_numerics"]["dx_um"] == 0.005
    assert result["authority_numerics"]["minimum_timestep_ratio_to_proposal"] == 0.25
    assert result["authority_numerics"]["topology_change_policy"] == "continue_gas_cavity"
    assert result["authority_numerics"][
        "surface_state_remap_backend"] == "common_refinement"
    assert result["base_endpoints"]["5nm"][
        "reassigned_unresolved_material_node_count"] == 1

    payload["configuration"]["minimum_step_duration_s"] = 0.002
    base5.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="monotonically refine"):
        FREEZE.freeze(base10, base5, calibration, azimuth, charging)


def test_krueger_freeze_refuses_underpriced_subcell_material_closure(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    payload = json.loads(base5.read_text(encoding="utf-8"))
    payload["history"][0].update({
        "reassigned_unresolved_material_nodes": 1,
        "unresolved_material_volume_upper_bound_m3": 2.0 * (0.005e-6) ** 3,
    })
    base5.write_text(json.dumps(payload), encoding="utf-8")
    calibration = _calibration(tmp_path / "calibration.json")
    azimuth, charging = _supporting_evidence(tmp_path)

    with pytest.raises(ValueError, match="subcell-material bound"):
        FREEZE.freeze(base10, base5, calibration, azimuth, charging)


def test_krueger_freeze_binds_the_single_fine_grid_correction(tmp_path):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    fine_payload = json.loads(base5.read_text(encoding="utf-8"))
    fine_payload["configuration"].update({
        "effective_mask_crosslinked_growth_fraction": 0.91,
        "oxide_etch_yield_scale": 0.56,
    })
    base5.write_text(json.dumps(fine_payload), encoding="utf-8")
    calibration = _calibration(tmp_path / "calibration.json")
    grid = tmp_path / "grid.json"
    grid_payload = {
        "schema": "petch.krueger-2024.base-grid-correction.v1",
        "protocol_sha256": FREEZE._sha(FREEZE.PROTOCOL),
        "held_out_profile_data_read": False,
        "calibration_derivation": {"sha256": FREEZE._sha(calibration)},
        "proposed_configuration": {
            "effective_mask_crosslinked_growth_fraction": 0.91,
            "oxide_etch_yield_scale": 0.56,
        },
    }
    grid_payload["proposal_sha256"] = sha256(json.dumps(
        grid_payload, sort_keys=True,
        separators=(",", ":")).encode("utf-8")).hexdigest()
    grid.write_text(json.dumps(grid_payload), encoding="utf-8")
    azimuth = tmp_path / "azimuth.json"
    azimuth.write_text(json.dumps({
        "axisymmetric_order_refinement_pass": True,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "held_out_profile_data_read": False,
    }), encoding="utf-8")
    charging = tmp_path / "charging.json"
    charging.write_text(json.dumps({
        "schema": "petch.krueger-2024.charging-causality.v1",
        "held_out_profile_data_read": False,
        "conservation": {"maximum_charge_conservation_relative_error": 0.0},
        "paired_exact_hard_visibility_audit": {
            "charged_over_zero_floor_ion_flux": 1.0},
    }), encoding="utf-8")

    result = FREEZE.freeze(
        base10, base5, calibration, azimuth, charging, grid)

    assert result["frozen_physics"][
        "effective_mask_crosslinked_growth_fraction"] == 0.91
    assert result["grid_correction"]["sha256"] == FREEZE._sha(grid)


def test_krueger_transfer_supervisor_has_complete_blind_case_matrix():
    assert len(CAMPAIGN._selected("oxygen")) == 4
    assert len(CAMPAIGN._selected("power")) == 4
    assert len(CAMPAIGN._selected("all")) == 8
    assert CAMPAIGN.SUCCESS_STATUSES == {"complete"}


def test_krueger_transfer_supervisor_verifies_frozen_code_and_boundary_data(
        tmp_path, monkeypatch):
    source = tmp_path / "engine.py"
    boundary = tmp_path / "boundary.csv"
    source.write_text("engine-v1\n", encoding="utf-8")
    boundary.write_text("flux,1\n", encoding="utf-8")
    monkeypatch.setattr(CAMPAIGN, "ROOT", tmp_path)
    freeze = {
        "source_sha256": {"engine.py": CAMPAIGN._sha(source)},
        "boundary_data_sha256": {"boundary.csv": CAMPAIGN._sha(boundary)},
    }

    CAMPAIGN._verify_frozen_files(freeze)
    source.write_text("engine-v2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="engine.py"):
        CAMPAIGN._verify_frozen_files(freeze)


def test_krueger_transfer_supervisor_refuses_mutated_reveal_before_subprocess(tmp_path):
    freeze = tmp_path / "freeze.json"
    freeze.write_text(json.dumps({
        "schema": "petch.krueger-2024.frozen-physics-reveal.v2",
        "held_out_profile_data_read": False,
        "reveal_sha256": "0" * 64,
    }), encoding="utf-8")
    args = SimpleNamespace(
        freeze=str(freeze), output_root=str(tmp_path / "runs"), case_set="all",
        transport_device="cpu", max_wall_s=1.0, maximum_resume_count=0)

    with pytest.raises(ValueError, match="checksum"):
        CAMPAIGN.run(args)


def test_krueger_transfer_replays_the_sealed_remap_backend(tmp_path, monkeypatch):
    base10 = _base(tmp_path / "base10.json", 0.01)
    base5 = _base(tmp_path / "base5.json", 0.005)
    calibration = _calibration(tmp_path / "calibration.json")
    azimuth, charging = _supporting_evidence(tmp_path)
    frozen = FREEZE.freeze(base10, base5, calibration, azimuth, charging)
    freeze_path = tmp_path / "freeze.json"
    freeze_path.write_text(json.dumps(frozen), encoding="utf-8")
    commands = []

    def fake_run(command, cwd, check):
        commands.append(command)
        output = Path(command[command.index("--output") + 1])
        output.mkdir(parents=True, exist_ok=True)
        (output / "audit.json").write_text(json.dumps({
            "status": "complete",
            "config_hash": "b" * 64,
            "final_metrics": {},
        }), encoding="utf-8")

    monkeypatch.setattr(CAMPAIGN, "_selected", lambda _: (("base", ()),))
    monkeypatch.setattr(CAMPAIGN.subprocess, "run", fake_run)
    args = SimpleNamespace(
        freeze=str(freeze_path), output_root=str(tmp_path / "runs"),
        case_set="all", transport_device="cuda:0", max_wall_s=1.0,
        maximum_resume_count=0)

    CAMPAIGN.run(args)

    command = commands[0]
    assert command[command.index("--transport-device") + 1] == "cuda:0"
    assert command[command.index("--surface-state-remap-backend") + 1] == (
        "common_refinement")
