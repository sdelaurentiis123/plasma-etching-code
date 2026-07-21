#!/usr/bin/env python3
"""Produce the sealed physics reveal under successor protocol K24-PETCH-R5.

R5 (preregistered 2026-07-21 with every held-out observation sealed) replaces the sampled
transport with the deterministic analytic-occlusion extruded operator and calibrates by an
exact-Jacobian Newton step (base + two probes).  This tool applies the same gate structure
as the R1.10 freeze -- +/-5 nm base tolerance on both declared observables, calibration-
derivation binding, launch-manifest executable checksums against current source, azimuth and
charging evidence -- against a single 10 nm authority endpoint produced by the deterministic
operator.  The base and probe runs informing the Jacobian may precede the frozen candidate's
epoch (provenance-only configuration enrichment); the gate applies to the candidate's own
endpoint on its own recorded epoch.  Emits frozen-physics-reveal.v3 for the R5 transfer
campaign.
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


def _load_jacobian_derivation(path, endpoint_pair):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    _verify_embedded_sha256(payload, "proposal_sha256")
    candidate = payload.get("candidate", {})
    if (payload.get("schema") != "petch.krueger-2024.r5-jacobian.v1"
            or payload.get("protocol_id") != "K24-PETCH-R5"
            or payload.get("protocol_sha256") != _sha(PROTOCOL)
            or payload.get("held_out_profile_data_read") is not False
            or not np.isclose(
                float(candidate.get("fraction", np.nan)),
                endpoint_pair["effective_mask_crosslinked_growth_fraction"])
            or not np.isclose(
                float(candidate.get("yield_scale", np.nan)),
                endpoint_pair["oxide_etch_yield_scale"])
            or float(payload.get("jacobian_condition", np.inf)) > 1.0e4
            or len(payload.get("inputs", ())) != 3):
        raise ValueError("10 nm endpoint is not bound to its R5 Jacobian derivation")
    return path, payload


def _load_r5_launch_manifest(path, endpoint_pair):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    _verify_embedded_sha256(payload, "launch_sha256")
    source = payload.get("source_epoch", {})
    operator = payload.get("numerical_operator", {})
    exchange = operator.get("deterministic_exchange", {})
    if (payload.get("schema") != "petch.krueger-2024.r5-base-authority.v1"
            or payload.get("protocol_id") != "K24-PETCH-R5"
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
            or float(operator.get("dx_um", np.nan)) != 0.01
            or operator.get("radiosity_backend") != "deterministic_extruded_2d"
            or exchange.get("exchange_method") != "analytic_occlusion"
            or float(exchange.get("exchange_geometry_tolerance", np.nan)) != 1.0e-9
            or int(exchange.get("maximum_refinement_level", -1)) != 24
            or operator.get("surface_state_remap_backend") != "common_refinement"):
        raise ValueError(f"invalid R5 authority launch manifest: {path}")
    executable = payload.get("executable_source_sha256", {})
    required = {
        "scripts/krueger_2024_trench_pilot.py",
        "src/petch/boundary_transport_3d.py",
        "src/petch/feature_step_3d.py",
        "src/petch/deterministic_exchange_2d.py",
        "src/petch/extruded_exchange_3d.py",
        "src/petch/surface_common_refinement_3d.py",
    }
    if not required.issubset(executable):
        raise ValueError("R5 launch manifest omits executable operator checksums")
    mismatch = {
        name: {"launch": digest, "current": _sha(ROOT / name)
               if (ROOT / name).is_file() else None}
        for name, digest in executable.items()
        if ((ROOT / name).is_file() and _sha(ROOT / name) != digest)
        or not (ROOT / name).is_file()
    }
    if mismatch:
        raise ValueError(
            f"current executable source does not match the R5 launch: {mismatch}")
    return path, payload


def freeze_r5(base_path, calibration_path, azimuth_path, charging_path, *,
              launch_path):
    base, summary = _load_complete(base_path, 0.01)
    config = base["configuration"]
    if config.get("topology_change_policy") != "continue_gas_cavity":
        raise ValueError(
            "R5 authority endpoint did not use the gas-cavity continuation policy")
    if config.get("radiosity_backend") != "deterministic_extruded_2d":
        raise ValueError("R5 authority endpoint did not use the deterministic backend")
    pair = {name: float(config[name]) for name in PARAMETER_NAMES}

    error = {
        name: summary[name] - TARGETS_NM[name]
        for name in ("mask_opening_nm", "etch_depth_nm")
    }
    outside = {
        name: value for name, value in error.items()
        if abs(value) > BASE_TOLERANCE_NM
    }
    if outside:
        raise ValueError(
            "candidate base endpoint is outside the declared +/-5 nm freeze "
            f"tolerance: {outside}; per R5 a second candidate miss ends the round "
            "as a failed calibration")

    calibration_path, calibration = _load_jacobian_derivation(calibration_path, pair)
    launch_manifest_path, launch = _load_r5_launch_manifest(launch_path, pair)

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
        "ballistic_face_quadrature_points", "radiosity_backend",
        "radiosity_relative_tolerance", "radiosity_maximum_iterations",
        "radiosity_enabled", "deterministic_exchange",
        "adaptive_profile_timestep",
        "target_displacement_cells", "maximum_displacement_cells",
        "adaptive_shrink_factor", "adaptive_growth_factor",
        "adaptive_safety_factor", "maximum_accepted_steps",
        "profile_reinitialization", "surface_state_remap_backend",
    )
    physical_geometry_keys = (
        "cell_width_um", "cell_length_um", "domain_height_um", "opening_width_um",
        "mask_thickness_um", "substrate_top_um", "initial_etched_depth_um")

    payload = {
        "schema": "petch.krueger-2024.frozen-physics-reveal.v3",
        "protocol_id": "K24-PETCH-R5",
        "protocol_sha256": _sha(PROTOCOL),
        "protocol_amendment": {
            "id": "R5",
            "authority_fidelity_dx_um": 0.01,
            "operator": "deterministic_analytic_occlusion_extruded_mean_field",
            "five_nm_status": "post-hoc refinement confirmation pending",
            "parameter_note": (
                "the frozen pair is declared operator-relative at 10 nm fidelity "
                "under the sampling-free deterministic exchange; grid-relative "
                "effective-parameter compensation is acknowledged and tested post "
                "hoc by the 5 nm confirmation under the frozen pair"),
        },
        "scientific_status": (
            "frozen two-observable calibration reveal for held-out transfer at 10 nm "
            "numerical fidelity under protocol R5; formal predictive validation "
            "remains conditional on boundary/mechanism evidence and uncertainty"),
        "base_target_table_sha256": _sha(BASE_TARGETS),
        "frozen_physics": {
            **pair,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
            "charging_policy": (
                "disabled for this high-energy Krueger transfer after paired causal "
                "null audit; the unified charging operator remains available for "
                "charge-sensitive cases"),
        },
        "authority_numerics": {
            "role": ("uniform 10 nm calibration confirmation and held-out authority "
                     "under protocol K24-PETCH-R5"),
            "dx_um": 0.01,
            "minimum_step_duration_s": float(config["minimum_step_duration_s"]),
            **{key: config[key] for key in operator_keys},
            "geometry": {
                key: config["geometry"][key] for key in physical_geometry_keys},
            "topology_change_policy": config["topology_change_policy"],
        },
        "base_endpoints": {"10nm": summary},
        "base_target_error_nm": error,
        "base_launch_manifests": {
            "10nm_authority": {
                "path_name": launch_manifest_path.name,
                "sha256": _sha(launch_manifest_path),
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
    parser.add_argument("--base", required=True)
    parser.add_argument("--calibration", required=True)
    parser.add_argument("--azimuth-audit", required=True)
    parser.add_argument("--charging-audit", required=True)
    parser.add_argument("--launch", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = freeze_r5(
        args.base, args.calibration, args.azimuth_audit, args.charging_audit,
        launch_path=args.launch)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print(output)
    print(payload["reveal_sha256"])


if __name__ == "__main__":
    main()
