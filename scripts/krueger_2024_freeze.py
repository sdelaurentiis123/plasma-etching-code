#!/usr/bin/env python3
"""Freeze the Krüger physics reveal before any held-out profile is opened.

The artifact binds the two calibrated physical scalars, the certified 3-D boundary closure, the
10/5 nm base endpoints, causal charging audit, protocol, and executable source.  It deliberately
does not import or read the transfer-observation table.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md"
BASE_TARGETS = ROOT / "data" / "experimental" / "krueger_2024" / "base_case_metrics.csv"
FROZEN_SOURCE = tuple(sorted((ROOT / "src" / "petch").glob("*.py"))) + (
    ROOT / "scripts" / "krueger_2024_coupled_correction.py",
    ROOT / "scripts" / "krueger_2024_grid_correction.py",
    ROOT / "scripts" / "krueger_2024_trench_pilot.py",
    ROOT / "scripts" / "krueger_2024_transfer_campaign.py",
)
FROZEN_BOUNDARY_DATA = tuple(
    ROOT / "data" / "experimental" / "krueger_2024" / name
    for name in (
        "base_case_boundary_fluxes.csv",
        "digitized_figure16_metadata.json",
        "digitized_figure16a_transfer_fluxes.csv",
        "digitized_figure16b_power_ieads.csv",
        "digitized_figure4_iead.csv",
        "digitized_figure4_iead_metadata.json",
    )
)
R19_FIXED_PAIR = {
    "effective_mask_crosslinked_growth_fraction": 0.9004722559883319,
    "oxide_etch_yield_scale": 0.5586489665864749,
}


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _verify_embedded_sha256(payload, field):
    canonical = dict(payload)
    claimed = canonical.pop(field, None)
    actual = sha256(json.dumps(
        canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if claimed != actual:
        raise ValueError(f"{field} does not match its artifact content")


def _same_parameter_pair(candidate, expected):
    names = tuple(R19_FIXED_PAIR)
    try:
        return all(float(candidate[name]) == float(expected[name]) for name in names)
    except (KeyError, TypeError, ValueError):
        return False


def _load_calibration_derivation(path, endpoint_pair):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    _verify_embedded_sha256(payload, "proposal_sha256")
    schema = payload.get("schema")
    proposed = payload.get("proposed_configuration", {})
    common_valid = (
        payload.get("protocol_sha256") == _sha(PROTOCOL)
        and payload.get("held_out_profile_data_read") is False
        and _same_parameter_pair(proposed, endpoint_pair))
    if schema == "petch.krueger-2024.base-axisymmetric-secant.v1":
        valid = common_valid
    elif schema == "petch.krueger-2024.development-trust-region-proposal.v1":
        targets = payload.get("calibration_targets", {})
        contract = payload.get("next_evaluation_contract", {})
        valid = (
            common_valid
            and payload.get("protocol_id") == "K24-PETCH-R1.9"
            and payload.get("authority") is False
            and targets.get("sha256") == _sha(BASE_TARGETS)
            and int(contract.get("count", -1)) == 1
            and proposed.get("topology_change_policy") == "continue_gas_cavity"
            and _same_parameter_pair(proposed, R19_FIXED_PAIR))
    else:
        valid = False
    if not valid:
        raise ValueError(
            "10 nm endpoint is not bound to its final base-only derivation")
    return path, payload


def _load_launch_manifest(path, *, expected_schema, authority_candidate,
                          endpoint_pair):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    _verify_embedded_sha256(payload, "launch_sha256")
    source = payload.get("source_epoch", {})
    if (payload.get("schema") != expected_schema
            or payload.get("protocol_id") != "K24-PETCH-R1.9"
            or payload.get("authority_candidate") is not authority_candidate
            or payload.get("held_out_profile_data_read") is not False
            or payload.get("calibration_performed_by_this_launch") is not False
            or source.get("git_dirty") is not False
            or not source.get("git_revision")
            or not source.get("archive_sha256")
            or not _same_parameter_pair(
                payload.get("fixed_parameters", {}), endpoint_pair)):
        raise ValueError(f"invalid checksum-bound base launch manifest: {path}")
    return path, payload


def _bind_launch_manifests(launch_10nm_path, launch_5nm_path, endpoint_pair):
    path_5, launch_5 = _load_launch_manifest(
        launch_5nm_path,
        expected_schema="petch.krueger-2024.base-authority-retry.v2",
        authority_candidate=True, endpoint_pair=endpoint_pair)
    path_10, launch_10 = _load_launch_manifest(
        launch_10nm_path,
        expected_schema="petch.krueger-2024.base-refinement-companion.v1",
        authority_candidate=False, endpoint_pair=endpoint_pair)
    source_5 = launch_5["source_epoch"]
    source_10 = launch_10["source_epoch"]
    source_keys = ("git_revision", "git_dirty", "archive_sha256")
    if any(source_5.get(key) != source_10.get(key) for key in source_keys):
        raise ValueError("10/5 nm launch manifests do not share one clean source epoch")
    if (launch_10.get("paired_authority_launch_sha256")
            != launch_5["launch_sha256"]):
        raise ValueError("10 nm launch is not paired to the 5 nm authority launch")
    authority_operator = launch_5.get("numerical_operator", {})
    companion_difference = launch_10.get("numerical_difference_from_authority", {})
    if (float(authority_operator.get("duration_s", np.nan)) != 60.0
            or int(authority_operator.get("n_steps", -1)) != 60
            or float(authority_operator.get("dx_um", np.nan)) != 0.005
            or authority_operator.get("surface_state_remap_backend")
            != "common_refinement"
            or float(companion_difference.get("dx_um", np.nan)) != 0.01
            or companion_difference.get("all_other_physical_and_numerical_controls")
            != "identical to the paired 5-nm launch"):
        raise ValueError("paired launch manifests do not declare the R1.9 refinement operator")

    executable = launch_5.get("executable_source_sha256", {})
    required = {
        "scripts/krueger_2024_trench_pilot.py",
        "src/petch/boundary_transport_3d.py",
        "src/petch/feature_step_3d.py",
        "src/petch/surface_common_refinement_3d.py",
    }
    if not required.issubset(executable):
        raise ValueError("5 nm launch manifest omits executable operator checksums")
    mismatch = {}
    for name, digest in executable.items():
        source_path = ROOT / name
        current = _sha(source_path) if source_path.is_file() else None
        if current != digest:
            mismatch[name] = {"launch": digest, "current": current}
    if mismatch:
        raise ValueError(
            f"current executable source does not match the base launch: {mismatch}")
    return {
        "10nm": {
            "path_name": path_10.name,
            "sha256": _sha(path_10),
            "launch_sha256": launch_10["launch_sha256"],
        },
        "5nm": {
            "path_name": path_5.name,
            "sha256": _sha(path_5),
            "launch_sha256": launch_5["launch_sha256"],
        },
        "shared_source_epoch": {
            key: source_5[key] for key in source_keys
        },
    }


def _load_complete(path, expected_dx):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    config = payload.get("configuration", {})
    if (payload.get("status") != "complete"
            or config.get("boundary_case") != "base"
            or not np.isclose(float(config.get("duration_s", np.nan)), 60.0)
            or not np.isclose(float(config.get("dx_um", np.nan)), expected_dx)
            or config.get("ion_azimuthal_closure") != "axisymmetric_uniform"
            or int(config.get("ion_azimuthal_order", -1)) != 16
            or config.get("ballistic_transport") != "face_gather"
            or int(config.get("ballistic_face_quadrature_points", -1)) != 3
            or config.get("surface_state_remap_backend") not in (
                "indexed_knn", "common_refinement")
            or not config.get("compressed_boundary_quadrature", False)
            or not config.get("radiosity_enabled", False)
            or not config.get("adaptive_profile_timestep", False)):
        raise ValueError(f"base endpoint does not satisfy the R1.9 operator: {path}")
    history = payload.get("history", ())
    if not history:
        raise ValueError(f"base endpoint has no conservation history: {path}")
    maximum_ledger_residual = max(
        float(item.get("maximum_material_ledger_residual_units_m2", 0.0))
        for item in history)
    maximum_radiosity_error = max(
        float(item.get("maximum_neutral_radiosity_relative_balance_error", 0.0))
        for item in history)
    if maximum_ledger_residual != 0.0:
        raise ValueError(f"base endpoint material ledger is not exactly closed: {path}")
    if maximum_radiosity_error > 1e-9:
        raise ValueError(f"base endpoint radiosity balance is unresolved: {path}")
    cell_volume_m3 = (float(expected_dx) * 1e-6) ** 3
    for item in history:
        count = int(item.get("reassigned_unresolved_material_nodes", 0))
        bound = float(item.get("unresolved_material_volume_upper_bound_m3", 0.0))
        if count < 0 or bound < 0.0 or bound > count * cell_volume_m3 * (1.0 + 1e-12):
            raise ValueError(f"base endpoint subcell-material bound is invalid: {path}")
    metrics = payload["final_metrics"]
    return payload, {
        "path_name": path.name,
        "sha256": _sha(path),
        "config_hash": payload["config_hash"],
        "dx_um": float(expected_dx),
        "mask_opening_nm": float(metrics["mask_opening_nm"]),
        "etch_depth_nm": float(metrics["etch_depth_nm"]),
        "maximum_feature_width_nm": float(metrics["maximum_feature_width_nm"]),
        "remaining_mask_thickness_nm": float(metrics["remaining_mask_thickness_nm"]),
        "asymmetry_cell_count": float(metrics["asymmetry_cell_count"]),
        "maximum_material_ledger_residual_units_m2": maximum_ledger_residual,
        "maximum_neutral_radiosity_relative_balance_error": maximum_radiosity_error,
        "adaptive_retry_count": int(sum(
            len(item.get("rejected_trials", ())) for item in history)),
        "reassigned_unresolved_material_node_count": int(sum(
            int(item.get("reassigned_unresolved_material_nodes", 0))
            for item in history)),
        "maximum_unresolved_material_volume_upper_bound_m3": float(max(
            float(item.get("unresolved_material_volume_upper_bound_m3", 0.0))
            for item in history)),
        "cumulative_unresolved_material_volume_upper_bound_m3": float(sum(
            float(item.get("unresolved_material_volume_upper_bound_m3", 0.0))
            for item in history)),
        "maximum_unresolved_cavity_volume_upper_bound_m3": float(max(
            float(item.get("unresolved_gas_cavity_volume_upper_bound_m3", 0.0))
            for item in history)),
    }


def freeze(base_10nm_path, base_5nm_path, calibration_path, azimuth_path,
           charging_path, grid_correction_path=None, *, launch_10nm_path,
           launch_5nm_path):
    base_10, summary_10 = _load_complete(base_10nm_path, 0.01)
    base_5, summary_5 = _load_complete(base_5nm_path, 0.005)
    config_10 = base_10["configuration"]
    config_5 = base_5["configuration"]
    if config_5.get("topology_change_policy") != "continue_gas_cavity":
        raise ValueError(
            "5 nm authority endpoint did not use the R1.9 gas-cavity continuation policy")
    parameter_names = (
        "effective_mask_crosslinked_growth_fraction", "oxide_etch_yield_scale")
    pair_10 = {name: float(config_10[name]) for name in parameter_names}
    pair_5 = {name: float(config_5[name]) for name in parameter_names}
    calibration_path, calibration = _load_calibration_derivation(
        calibration_path, pair_10)
    launch_manifests = _bind_launch_manifests(
        launch_10nm_path, launch_5nm_path, pair_10)
    operator_keys = (
        "duration_s", "n_steps", "n_position", "compressed_boundary_quadrature",
        "neutral_speed_quadrature", "neutral_tensor_velocity_quadrature_active",
        "neutral_transverse_order", "neutral_normal_order",
        "neutral_direction_polar_order", "neutral_direction_azimuthal_order",
        "ion_energy_bin_eV", "ion_angle_bin_deg", "ion_azimuthal_closure",
        "ion_azimuthal_order", "ballistic_transport",
        "ballistic_face_quadrature_points", "radiosity_rays_per_face",
        "radiosity_relative_tolerance", "radiosity_maximum_iterations",
        "radiosity_enabled", "seed", "adaptive_profile_timestep",
        "target_displacement_cells",
        "maximum_displacement_cells", "adaptive_shrink_factor",
        "adaptive_growth_factor", "adaptive_safety_factor",
        "maximum_accepted_steps", "profile_reinitialization",
        "surface_state_remap_backend",
    )
    inconsistent = {
        key: {"10nm": config_10.get(key), "5nm": config_5.get(key)}
        for key in operator_keys if config_10.get(key) != config_5.get(key)
    }
    physical_geometry_keys = (
        "cell_width_um", "cell_length_um", "domain_height_um", "opening_width_um",
        "mask_thickness_um", "substrate_top_um", "initial_etched_depth_um")
    inconsistent_geometry = {
        key: {
            "10nm": config_10.get("geometry", {}).get(key),
            "5nm": config_5.get("geometry", {}).get(key),
        }
        for key in physical_geometry_keys
        if config_10.get("geometry", {}).get(key)
        != config_5.get("geometry", {}).get(key)
    }
    if inconsistent or inconsistent_geometry:
        raise ValueError(
            "10/5 nm endpoints do not share one numerical operator: "
            f"settings={inconsistent}, geometry={inconsistent_geometry}")
    minimum_step_10 = float(config_10["minimum_step_duration_s"])
    minimum_step_5 = float(config_5["minimum_step_duration_s"])
    if (not np.isfinite(minimum_step_10) or not np.isfinite(minimum_step_5)
            or minimum_step_10 <= 0.0 or minimum_step_5 <= 0.0
            or minimum_step_5 > minimum_step_10):
        raise ValueError("5 nm minimum timestep must monotonically refine the 10 nm safety floor")
    changed = any(not np.isclose(pair_10[name], pair_5[name]) for name in parameter_names)
    grid_correction = None
    if changed:
        if grid_correction_path is None:
            raise ValueError("a changed 5 nm pair requires a checksum-bound grid correction")
        grid_correction_path = Path(grid_correction_path)
        grid_payload = json.loads(grid_correction_path.read_text(encoding="utf-8"))
        _verify_embedded_sha256(grid_payload, "proposal_sha256")
        grid_pair = grid_payload.get("proposed_configuration", {})
        if (grid_payload.get("schema")
                != "petch.krueger-2024.base-grid-correction.v1"
                or grid_payload.get("protocol_sha256") != _sha(PROTOCOL)
                or grid_payload.get("held_out_profile_data_read") is not False
                or grid_payload.get("calibration_derivation", {}).get("sha256")
                != _sha(calibration_path)
                or any(not np.isclose(float(grid_pair.get(name, np.nan)), pair_5[name])
                       for name in parameter_names)):
            raise ValueError("5 nm endpoint is not bound to its grid-correction derivation")
        grid_correction = {
            "path_name": grid_correction_path.name,
            "sha256": _sha(grid_correction_path),
            "proposal_sha256": grid_payload["proposal_sha256"],
        }
    elif grid_correction_path is not None:
        raise ValueError("unchanged 5 nm pair must not claim a grid correction")

    azimuth_path = Path(azimuth_path)
    azimuth = json.loads(azimuth_path.read_text(encoding="utf-8"))
    if (not azimuth.get("axisymmetric_order_refinement_pass")
            or int(azimuth.get("production_azimuthal_order", -1)) != 16
            or int(azimuth.get("reference_azimuthal_order", -1)) != 32
            or azimuth.get("held_out_profile_data_read") is not False):
        raise ValueError("R1.5 azimuth evidence is incomplete")
    charging_path = Path(charging_path)
    charging = json.loads(charging_path.read_text(encoding="utf-8"))
    if (charging.get("schema") != "petch.krueger-2024.charging-causality.v1"
            or charging.get("held_out_profile_data_read") is not False
            or charging["conservation"]["maximum_charge_conservation_relative_error"] > 1e-12):
        raise ValueError("charging causal evidence is incomplete or nonconservative")

    target = {"mask_opening_nm": 45.0, "etch_depth_nm": 825.0}
    error_5 = {
        name: summary_5[name] - value for name, value in target.items()}
    if (abs(error_5["mask_opening_nm"]) > 5.0
            or abs(error_5["etch_depth_nm"]) > 5.0):
        raise ValueError("5 nm calibration endpoint is outside the one-cell R1.7 freeze gate")

    payload = {
        "schema": "petch.krueger-2024.frozen-physics-reveal.v2",
        "protocol_id": "K24-PETCH-R1.9",
        "protocol_sha256": _sha(PROTOCOL),
        "scientific_status": (
            "frozen two-observable calibration reveal for held-out transfer; formal predictive "
            "validation remains conditional on boundary/mechanism evidence and uncertainty"),
        "base_target_table_sha256": _sha(BASE_TARGETS),
        "frozen_physics": {
            **pair_5,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "charging_policy": (
                "disabled for this high-energy Krueger transfer after paired causal null audit; "
                "the unified charging operator remains available for charge-sensitive cases"),
        },
        "proposal_numerics": {
            "role": "10 nm multi-fidelity proposal and discrepancy model; never held-out authority",
            "dx_um": 0.01,
            "minimum_step_duration_s": minimum_step_10,
            **{key: config_10[key] for key in operator_keys},
            "geometry": {
                key: config_10["geometry"][key] for key in physical_geometry_keys},
        },
        "authority_numerics": {
            "role": "uniform 5 nm calibration confirmation and held-out authority",
            "dx_um": 0.005,
            "minimum_step_duration_s": minimum_step_5,
            "minimum_timestep_ratio_to_proposal": (
                minimum_step_5 / minimum_step_10),
            **{key: config_5[key] for key in operator_keys},
            "geometry": {
                key: config_5["geometry"][key] for key in physical_geometry_keys},
            "subcell_material_policy": (
                "suppress only newly born per-material components with fewer than eight unique "
                "periodic nodes; retain prior ownership; resolved/existing splits still refuse"),
            "topology_change_policy": config_5["topology_change_policy"],
        },
        "base_endpoints": {"10nm": summary_10, "5nm": summary_5},
        "base_launch_manifests": launch_manifests,
        "base_calibration_derivation": {
            "schema": calibration["schema"],
            "path_name": calibration_path.name,
            "sha256": _sha(calibration_path),
            "proposal_sha256": calibration["proposal_sha256"],
        },
        "base_endpoint_difference": {
            name: summary_5[name] - summary_10[name]
            for name in (
                "mask_opening_nm", "etch_depth_nm", "maximum_feature_width_nm",
                "remaining_mask_thickness_nm", "asymmetry_cell_count")
        },
        "base_grid_difference": (
            {
                name: summary_5[name] - summary_10[name]
                for name in (
                    "mask_opening_nm", "etch_depth_nm", "maximum_feature_width_nm",
                    "remaining_mask_thickness_nm", "asymmetry_cell_count")
            }
            if not changed else None),
        "base_difference_interpretation": (
            "pure fixed-parameter 10-to-5-nm grid difference"
            if not changed else
            "combined grid-plus-single-correction endpoint difference; the first fixed-pair 5 nm "
            "endpoint is checksum-bound inside grid_correction"),
        "base_5nm_target_error": error_5,
        "grid_correction": grid_correction,
        "azimuth_audit": {"path_name": azimuth_path.name, "sha256": _sha(azimuth_path)},
        "charging_causal_audit": {
            "path_name": charging_path.name,
            "sha256": _sha(charging_path),
            "charged_over_zero_floor_ion_flux": float(
                charging["paired_exact_hard_visibility_audit"]
                ["charged_over_zero_floor_ion_flux"]),
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): _sha(path) for path in FROZEN_SOURCE
        },
        "boundary_data_sha256": {
            str(path.relative_to(ROOT)): _sha(path) for path in FROZEN_BOUNDARY_DATA
        },
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["reveal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-10nm", required=True)
    parser.add_argument("--base-5nm", required=True)
    parser.add_argument("--calibration-derivation", required=True)
    parser.add_argument("--azimuth-audit", required=True)
    parser.add_argument("--charging-audit", required=True)
    parser.add_argument("--launch-10nm", required=True)
    parser.add_argument("--launch-5nm", required=True)
    parser.add_argument("--grid-correction")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = freeze(
        args.base_10nm, args.base_5nm, args.calibration_derivation,
        args.azimuth_audit, args.charging_audit, args.grid_correction,
        launch_10nm_path=args.launch_10nm,
        launch_5nm_path=args.launch_5nm)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(payload["reveal_sha256"])


if __name__ == "__main__":
    main()
