#!/usr/bin/env python3
"""Produce the sealed physics reveal under protocol amendment R1.10 (10 nm authority).

Amendment R1.10 (preregistered 2026-07-20 with every held-out observation sealed) moves the
freeze fidelity to the uniform 10 nm operator and exercises the one permitted base-only
coupled correction.  This tool applies the same conservation, operator, launch-manifest, and
calibration-derivation gates as the R1.9 freeze, against a single 10 nm authority endpoint,
and additionally REQUIRES the +/-5 nm base tolerance on both declared observables before a
reveal is produced.  The emitted artifact keeps the frozen-physics-reveal.v2 schema so the
transfer campaign consumes it unchanged, with the amendment bound by protocol hash.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from krueger_2024_freeze import (  # noqa: F401
    BASE_TARGETS, FROZEN_BOUNDARY_DATA, FROZEN_SOURCE, PROTOCOL, ROOT,
    _load_complete, _sha, _verify_embedded_sha256)


PARAMETER_NAMES = (
    "effective_mask_crosslinked_growth_fraction", "oxide_etch_yield_scale")
BASE_TOLERANCE_NM = 5.0
TARGETS_NM = {"mask_opening_nm": 45.0, "etch_depth_nm": 825.0}


def _load_correction_derivation(path, endpoint_pair):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    _verify_embedded_sha256(payload, "proposal_sha256")
    proposed = payload.get("proposed_configuration", {})
    if (payload.get("schema") != "petch.krueger-2024.r4-secant.v1"
            or payload.get("protocol_sha256") != _sha(PROTOCOL)
            or payload.get("held_out_profile_data_read") is not False
            or any(not np.isclose(float(proposed.get(name, np.nan)),
                                  endpoint_pair[name])
                   for name in PARAMETER_NAMES)):
        raise ValueError(
            "10 nm endpoint is not bound to its R4 secant derivation")
    return path, payload


def _load_r110_launch_manifest(path, endpoint_pair):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    _verify_embedded_sha256(payload, "launch_sha256")
    source = payload.get("source_epoch", {})
    operator = payload.get("numerical_operator", {})
    if (payload.get("schema") != "petch.krueger-2024.r4-base-authority.v1"
            or payload.get("protocol_id") != "K24-PETCH-R4"
            or payload.get("authority_candidate") is not True
            or payload.get("held_out_profile_data_read") is not False
            or payload.get("calibration_performed_by_this_launch") is not False
            or source.get("git_dirty") is not False
            or not source.get("git_revision")
            or not source.get("archive_sha256")
            or any(not np.isclose(
                float(payload.get("fixed_parameters", {}).get(name, np.nan)),
                endpoint_pair[name]) for name in PARAMETER_NAMES)
            or float(operator.get("duration_s", np.nan)) != 60.0
            or int(operator.get("n_steps", -1)) != 60
            or float(operator.get("dx_um", np.nan)) != 0.01
            or operator.get("surface_state_remap_backend") != "common_refinement"):
        raise ValueError(f"invalid R1.10 authority launch manifest: {path}")
    executable = payload.get("executable_source_sha256", {})
    required = {
        "scripts/krueger_2024_trench_pilot.py",
        "src/petch/boundary_transport_3d.py",
        "src/petch/feature_step_3d.py",
        "src/petch/surface_common_refinement_3d.py",
    }
    if not required.issubset(executable):
        raise ValueError("R1.10 launch manifest omits executable operator checksums")
    mismatch = {
        name: {"launch": digest, "current": _sha(ROOT / name)
               if (ROOT / name).is_file() else None}
        for name, digest in executable.items()
        if ((ROOT / name).is_file() and _sha(ROOT / name) != digest)
        or not (ROOT / name).is_file()
    }
    if mismatch:
        raise ValueError(
            f"current executable source does not match the R1.10 launch: {mismatch}")
    return path, payload


def freeze_r110(base_10nm_path, calibration_path, azimuth_path, charging_path, *,
                launch_10nm_path):
    base_10, summary_10 = _load_complete(base_10nm_path, 0.01)
    config_10 = base_10["configuration"]
    if config_10.get("topology_change_policy") != "continue_gas_cavity":
        raise ValueError(
            "R1.10 authority endpoint did not use the gas-cavity continuation policy")
    pair_10 = {name: float(config_10[name]) for name in PARAMETER_NAMES}

    error = {
        name: summary_10[name] - TARGETS_NM[name]
        for name in ("mask_opening_nm", "etch_depth_nm")
    }
    outside = {
        name: value for name, value in error.items()
        if abs(value) > BASE_TOLERANCE_NM
    }
    if outside:
        raise ValueError(
            "corrected base endpoint is outside the declared +/-5 nm freeze "
            f"tolerance: {outside}; a second miss ends the campaign as a failed "
            "calibration")

    calibration_path, calibration = _load_correction_derivation(
        calibration_path, pair_10)
    launch_path, launch = _load_r110_launch_manifest(launch_10nm_path, pair_10)

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
            or charging["conservation"]["maximum_charge_conservation_relative_error"]
            > 1e-12):
        raise ValueError("charging causal audit is incomplete")

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
        "target_displacement_cells", "maximum_displacement_cells",
        "adaptive_shrink_factor", "adaptive_growth_factor",
        "adaptive_safety_factor", "maximum_accepted_steps",
        "profile_reinitialization", "surface_state_remap_backend",
    )
    physical_geometry_keys = (
        "cell_width_um", "cell_length_um", "domain_height_um", "opening_width_um",
        "mask_thickness_um", "substrate_top_um", "initial_etched_depth_um")

    payload = {
        "schema": "petch.krueger-2024.frozen-physics-reveal.v2",
        "protocol_id": "K24-PETCH-R4",
        "protocol_sha256": _sha(PROTOCOL),
        "protocol_amendment": {
            "id": "R4",
            "authority_fidelity_dx_um": 0.01,
            "five_nm_status": "post-hoc refinement confirmation pending",
            "parameter_note": (
                "the frozen pair is declared operator-relative at 10 nm fidelity; "
                "grid-relative effective-parameter compensation is acknowledged and "
                "tested post hoc by the 5 nm confirmation under the frozen pair"),
        },
        "scientific_status": (
            "frozen two-observable calibration reveal for held-out transfer at 10 nm "
            "numerical fidelity under amendment R1.10; formal predictive validation "
            "remains conditional on boundary/mechanism evidence and uncertainty"),
        "base_target_table_sha256": _sha(BASE_TARGETS),
        "frozen_physics": {
            **pair_10,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "charging_policy": (
                "disabled for this high-energy Krueger transfer after paired causal "
                "null audit; the unified charging operator remains available for "
                "charge-sensitive cases"),
        },
        "authority_numerics": {
            "role": ("uniform 10 nm calibration confirmation and held-out authority "
                     "under amendment R1.10"),
            "dx_um": 0.01,
            "minimum_step_duration_s": float(config_10["minimum_step_duration_s"]),
            **{key: config_10[key] for key in operator_keys},
            "geometry": {
                key: config_10["geometry"][key] for key in physical_geometry_keys},
            "topology_change_policy": config_10["topology_change_policy"],
        },
        "base_endpoints": {"10nm": summary_10},
        "base_target_error_nm": error,
        "base_launch_manifests": {
            "10nm_authority": {
                "path_name": launch_path.name,
                "sha256": _sha(launch_path),
                "launch_sha256": launch["launch_sha256"],
            },
            "shared_source_epoch": dict(launch["source_epoch"]),
        },
        "base_calibration_derivation": {
            "schema": calibration["schema"],
            "path_name": calibration_path.name,
            "sha256": _sha(calibration_path),
            "proposal_sha256": calibration["proposal_sha256"],
        },
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
    parser.add_argument("--calibration", required=True)
    parser.add_argument("--azimuth-audit", required=True)
    parser.add_argument("--charging-audit", required=True)
    parser.add_argument("--launch-10nm", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = freeze_r110(
        args.base_10nm, args.calibration, args.azimuth_audit, args.charging_audit,
        launch_10nm_path=args.launch_10nm)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print(output)
    print(payload["reveal_sha256"])


if __name__ == "__main__":
    main()
