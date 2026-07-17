#!/usr/bin/env python3
"""Derive the bounded R1.3--R1.6 base-only coupled corrections.

Only the two declared base calibration observables are consumed.  The mask fraction correction
uses the local derivative of the three already-completed mask-fraction endpoints.  The oxide-yield
correction uses its declared multiplicative-rate role, corrected for the fitted fraction response
and for the measured plane-to-axisymmetric substrate-rate ratio.  No transfer/held-out table is
opened by this program.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_TARGETS = ROOT / "data" / "experimental" / "krueger_2024" / "base_case_metrics.csv"
PROTOCOL = ROOT / "KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md"


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _read_complete_base(path):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    config = payload.get("configuration", {})
    if (payload.get("status") != "complete"
            or config.get("boundary_case") != "base"
            or not np.isclose(float(config.get("duration_s", np.nan)), 60.0)
            or not np.isclose(float(config.get("dx_um", np.nan)), 0.01)):
        raise ValueError(f"not a completed 60 s, 10 nm base endpoint: {path}")
    return payload, {
        "path_name": path.name,
        "sha256": _sha(path),
        "fraction": float(config.get("effective_mask_crosslinked_growth_fraction", 0.0)),
        "oxide_etch_yield_scale": float(config.get("oxide_etch_yield_scale", 1.0)),
        "mask_opening_nm": float(payload["final_metrics"]["mask_opening_nm"]),
        "etch_depth_nm": float(payload["final_metrics"]["etch_depth_nm"]),
    }


def _targets():
    with BASE_TARGETS.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {
        row["metric"]: float(row["value"])
        for row in rows
        if row["split"] == "calibration" and row["evidence_type"] == "experiment"
        and row["metric"] in {"mask_opening", "etch_depth"}
    }
    if set(selected) != {"mask_opening", "etch_depth"}:
        raise ValueError("the two declared base calibration targets are not uniquely available")
    return selected


def derive(endpoint_paths, coupled_path, azimuth_path):
    endpoints = [_read_complete_base(path)[1] for path in endpoint_paths]
    endpoints.sort(key=lambda item: item["fraction"])
    fraction = np.asarray([item["fraction"] for item in endpoints], dtype=float)
    if (len(endpoints) != 3 or np.unique(fraction).size != 3
            or not np.isclose(fraction[0], 0.0) or not np.isclose(fraction[-1], 1.0)
            or any(not np.isclose(item["oxide_etch_yield_scale"], 1.0)
                   for item in endpoints)):
        raise ValueError("fraction response requires three distinct unit-yield endpoints at 0/mid/1")

    _, coupled = _read_complete_base(coupled_path)
    azimuth_path = Path(azimuth_path)
    azimuth = json.loads(azimuth_path.read_text(encoding="utf-8"))
    if (not azimuth.get("axisymmetric_order_refinement_pass")
            or int(azimuth.get("production_azimuthal_order", -1)) != 16
            or int(azimuth.get("reference_azimuthal_order", -1)) != 32
            or azimuth.get("held_out_profile_data_read") is not False):
        raise ValueError("azimuth audit has not certified the R1.5 16-to-32 closure")

    targets = _targets()
    opening_fit = np.polyfit(
        fraction, [item["mask_opening_nm"] for item in endpoints], 2)
    depth_fit = np.polyfit(
        fraction, [item["etch_depth_nm"] for item in endpoints], 2)
    old_fraction = coupled["fraction"]
    opening_slope = float(2.0 * opening_fit[0] * old_fraction + opening_fit[1])
    if not np.isfinite(opening_slope) or opening_slope <= 0.0:
        raise ValueError("local base opening response is not positive")
    new_fraction = float(np.clip(
        old_fraction
        + (targets["mask_opening"] - coupled["mask_opening_nm"]) / opening_slope,
        0.0, 1.0))

    depth_fraction_factor = float(
        np.polyval(depth_fit, new_fraction) / np.polyval(depth_fit, old_fraction))
    variants = azimuth["variants"]
    plane_rate = abs(float(variants["single_published_plane"]
                           ["material_profile_rate"]["2"]
                           ["net_volume_rate_mesh_units3_s"]))
    axisymmetric_rate = abs(float(variants["axisymmetric_uniform_16"]
                                  ["material_profile_rate"]["2"]
                                  ["net_volume_rate_mesh_units3_s"]))
    if plane_rate <= 0.0 or axisymmetric_rate <= 0.0:
        raise ValueError("substrate endpoint rates must be nonzero")
    azimuth_rate_factor = axisymmetric_rate / plane_rate
    new_yield_scale = float(
        coupled["oxide_etch_yield_scale"] * targets["etch_depth"]
        / (coupled["etch_depth_nm"] * depth_fraction_factor * azimuth_rate_factor))
    if not np.isfinite(new_yield_scale) or new_yield_scale <= 0.0:
        raise ValueError("derived oxide yield scale is invalid")

    payload = {
        "schema": "petch.krueger-2024.base-coupled-correction.v1",
        "protocol_id": "K24-PETCH-R1.5",
        "protocol_sha256": _sha(PROTOCOL),
        "scientific_status": "single permitted base-only coupled correction; not validation",
        "base_target_table_sha256": _sha(BASE_TARGETS),
        "calibration_targets": {
            "mask_opening_nm": targets["mask_opening"],
            "etch_depth_nm": targets["etch_depth"],
        },
        "unit_yield_fraction_endpoints": endpoints,
        "coupled_endpoint": coupled,
        "azimuth_audit_sha256": _sha(azimuth_path),
        "derivation": {
            "opening_response": "quadratic-through-three-completed-unit-yield-base-endpoints",
            "opening_local_slope_nm_per_fraction": opening_slope,
            "depth_fraction_factor": depth_fraction_factor,
            "axisymmetric_to_plane_substrate_net_rate_factor": azimuth_rate_factor,
            "oxide_scale_rule": (
                "multiplicative energetic-SiO2 rate adapter, adjusted by the measured base "
                "fraction response and frozen-endpoint azimuth rate ratio"),
        },
        "proposed_configuration": {
            "effective_mask_crosslinked_growth_fraction": new_fraction,
            "oxide_etch_yield_scale": new_yield_scale,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
        },
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["proposal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def derive_path_resolved(
        endpoint_paths, plane_coupled_path, axisymmetric_coupled_path, azimuth_path):
    """Derive the bounded R1.6 secant update after one full axisymmetric path."""
    endpoints = [_read_complete_base(path)[1] for path in endpoint_paths]
    endpoints.sort(key=lambda item: item["fraction"])
    fraction = np.asarray([item["fraction"] for item in endpoints], dtype=float)
    if (len(endpoints) != 3 or np.unique(fraction).size != 3
            or not np.isclose(fraction[0], 0.0) or not np.isclose(fraction[-1], 1.0)
            or any(not np.isclose(item["oxide_etch_yield_scale"], 1.0)
                   for item in endpoints)):
        raise ValueError("path correction requires the same three unit-yield base endpoints")
    _, plane = _read_complete_base(plane_coupled_path)
    axisymmetric_payload, axisymmetric = _read_complete_base(axisymmetric_coupled_path)
    config = axisymmetric_payload["configuration"]
    if (config.get("ion_azimuthal_closure") != "axisymmetric_uniform"
            or int(config.get("ion_azimuthal_order", -1)) != 16):
        raise ValueError("path-resolved endpoint must use the certified 16-node closure")

    azimuth_path = Path(azimuth_path)
    azimuth = json.loads(azimuth_path.read_text(encoding="utf-8"))
    if (not azimuth.get("axisymmetric_order_refinement_pass")
            or int(azimuth.get("production_azimuthal_order", -1)) != 16
            or int(azimuth.get("reference_azimuthal_order", -1)) != 32
            or azimuth.get("held_out_profile_data_read") is not False):
        raise ValueError("azimuth audit has not certified the R1.5 16-to-32 closure")

    targets = _targets()
    opening_fit = np.polyfit(
        fraction, [item["mask_opening_nm"] for item in endpoints], 2)
    depth_fit = np.polyfit(
        fraction, [item["etch_depth_nm"] for item in endpoints], 2)
    old_fraction = axisymmetric["fraction"]
    old_yield = axisymmetric["oxide_etch_yield_scale"]
    opening_fraction_slope = float(
        2.0 * opening_fit[0] * old_fraction + opening_fit[1])
    depth_fraction_slope = float(
        (2.0 * depth_fit[0] * old_fraction + depth_fit[1]) * old_yield)

    delta_fraction = axisymmetric["fraction"] - plane["fraction"]
    delta_yield = (
        axisymmetric["oxide_etch_yield_scale"]
        - plane["oxide_etch_yield_scale"])
    midpoint_fraction = 0.5 * (axisymmetric["fraction"] + plane["fraction"])
    midpoint_opening_fraction_slope = float(
        2.0 * opening_fit[0] * midpoint_fraction + opening_fit[1])
    if abs(delta_yield) <= 1e-12:
        raise ValueError("path secant requires distinct oxide-yield scales")
    # This secant contains the finite plane-to-axisymmetric path transition.  It is used only to
    # estimate the sign/magnitude of mask/oxide coupling; the next full axisymmetric endpoint is
    # the arbiter.  The oxide diagonal uses its exact declared multiplicative role locally.
    opening_yield_slope = float(
        (axisymmetric["mask_opening_nm"] - plane["mask_opening_nm"]
         - midpoint_opening_fraction_slope * delta_fraction) / delta_yield)
    depth_yield_slope = float(
        axisymmetric["etch_depth_nm"] / old_yield)
    jacobian = np.asarray([
        [opening_fraction_slope, opening_yield_slope],
        [depth_fraction_slope, depth_yield_slope],
    ], dtype=float)
    response = np.asarray([
        targets["mask_opening"] - axisymmetric["mask_opening_nm"],
        targets["etch_depth"] - axisymmetric["etch_depth_nm"],
    ], dtype=float)
    if not np.all(np.isfinite(jacobian)) or np.linalg.cond(jacobian) > 1e6:
        raise ValueError("base-only path secant is singular or ill-conditioned")
    update = np.linalg.solve(jacobian, response)
    new_fraction = float(np.clip(old_fraction + update[0], 0.0, 1.0))
    new_yield = float(old_yield + update[1])
    if not np.isfinite(new_yield) or new_yield <= 0.0:
        raise ValueError("path-resolved oxide-yield proposal is invalid")

    payload = {
        "schema": "petch.krueger-2024.base-path-correction.v1",
        "protocol_id": "K24-PETCH-R1.6",
        "protocol_sha256": _sha(PROTOCOL),
        "scientific_status": "bounded base-only path-resolved correction; not validation",
        "base_target_table_sha256": _sha(BASE_TARGETS),
        "calibration_targets": {
            "mask_opening_nm": targets["mask_opening"],
            "etch_depth_nm": targets["etch_depth"],
        },
        "unit_yield_fraction_endpoints": endpoints,
        "plane_coupled_endpoint": plane,
        "axisymmetric_coupled_endpoint": axisymmetric,
        "azimuth_audit_sha256": _sha(azimuth_path),
        "derivation": {
            "method": "two-by-two-base-only-local-secant",
            "response_order": ["mask_opening_nm", "etch_depth_nm"],
            "parameter_order": [
                "effective_mask_crosslinked_growth_fraction",
                "oxide_etch_yield_scale",
            ],
            "jacobian": jacobian.tolist(),
            "condition_number": float(np.linalg.cond(jacobian)),
            "requested_response": response.tolist(),
            "parameter_update": update.tolist(),
            "opening_yield_secant_caveat": (
                "includes the finite plane-to-axisymmetric path transition and is used only as "
                "a proposal; the next complete 16-node axisymmetric base endpoint verifies it"),
            "depth_yield_rule": "local multiplicative energetic-SiO2 rate adapter",
        },
        "proposed_configuration": {
            "effective_mask_crosslinked_growth_fraction": new_fraction,
            "oxide_etch_yield_scale": new_yield,
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
        },
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["proposal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def derive_axisymmetric_secant(
        endpoint_paths, previous_axisymmetric_path, current_axisymmetric_path,
        azimuth_path):
    """Derive the final R1.6 update from two complete axisymmetric profile paths.

    The two complete axisymmetric paths supply one measured two-parameter secant.  A second,
    independent response direction is supplied by the already-completed unit-yield mask-fraction
    sweep.  Subtracting that fraction response from the measured path secant identifies the local
    oxide-yield response.  Only the two preregistered base observables enter this calculation.
    """
    endpoints = [_read_complete_base(path)[1] for path in endpoint_paths]
    endpoints.sort(key=lambda item: item["fraction"])
    fraction = np.asarray([item["fraction"] for item in endpoints], dtype=float)
    if (len(endpoints) != 3 or np.unique(fraction).size != 3
            or not np.isclose(fraction[0], 0.0) or not np.isclose(fraction[-1], 1.0)
            or any(not np.isclose(item["oxide_etch_yield_scale"], 1.0)
                   for item in endpoints)):
        raise ValueError("final secant requires the same three unit-yield base endpoints")

    previous_payload, previous = _read_complete_base(previous_axisymmetric_path)
    current_payload, current = _read_complete_base(current_axisymmetric_path)
    for payload in (previous_payload, current_payload):
        config = payload["configuration"]
        if (config.get("ion_azimuthal_closure") != "axisymmetric_uniform"
                or int(config.get("ion_azimuthal_order", -1)) != 16):
            raise ValueError("final secant endpoints must use the certified 16-node closure")

    azimuth_path = Path(azimuth_path)
    azimuth = json.loads(azimuth_path.read_text(encoding="utf-8"))
    if (not azimuth.get("axisymmetric_order_refinement_pass")
            or int(azimuth.get("production_azimuthal_order", -1)) != 16
            or int(azimuth.get("reference_azimuthal_order", -1)) != 32
            or azimuth.get("held_out_profile_data_read") is not False):
        raise ValueError("azimuth audit has not certified the R1.5 16-to-32 closure")

    targets = _targets()
    opening_fit = np.polyfit(
        fraction, [item["mask_opening_nm"] for item in endpoints], 2)
    depth_fit = np.polyfit(
        fraction, [item["etch_depth_nm"] for item in endpoints], 2)
    previous_parameters = np.asarray([
        previous["fraction"], previous["oxide_etch_yield_scale"]], dtype=float)
    current_parameters = np.asarray([
        current["fraction"], current["oxide_etch_yield_scale"]], dtype=float)
    parameter_secant = current_parameters - previous_parameters
    if abs(parameter_secant[1]) <= 1e-12:
        raise ValueError("axisymmetric secant requires distinct oxide-yield scales")

    midpoint_fraction = float(0.5 * (
        previous["fraction"] + current["fraction"]))
    midpoint_yield = float(0.5 * (
        previous["oxide_etch_yield_scale"]
        + current["oxide_etch_yield_scale"]))
    fraction_response = np.asarray([
        2.0 * opening_fit[0] * midpoint_fraction + opening_fit[1],
        (2.0 * depth_fit[0] * midpoint_fraction + depth_fit[1]) * midpoint_yield,
    ], dtype=float)
    observed_response = np.asarray([
        current["mask_opening_nm"] - previous["mask_opening_nm"],
        current["etch_depth_nm"] - previous["etch_depth_nm"],
    ], dtype=float)
    yield_response = (
        observed_response - fraction_response * parameter_secant[0]
    ) / parameter_secant[1]
    jacobian = np.column_stack([fraction_response, yield_response])
    condition_number = float(np.linalg.cond(jacobian))
    if not np.all(np.isfinite(jacobian)) or condition_number > 1e6:
        raise ValueError("final base-only axisymmetric secant is singular or ill-conditioned")

    requested_response = np.asarray([
        targets["mask_opening"] - current["mask_opening_nm"],
        targets["etch_depth"] - current["etch_depth_nm"],
    ], dtype=float)
    update = np.linalg.solve(jacobian, requested_response)
    proposed = current_parameters + update
    if (not np.all(np.isfinite(proposed)) or proposed[0] < 0.0 or proposed[0] > 1.0
            or proposed[1] <= 0.0):
        raise ValueError("final axisymmetric secant proposal lies outside physical bounds")

    payload = {
        "schema": "petch.krueger-2024.base-axisymmetric-secant.v1",
        "protocol_id": "K24-PETCH-R1.6",
        "protocol_sha256": _sha(PROTOCOL),
        "scientific_status": "final bounded 10 nm base-only correction; not validation",
        "base_target_table_sha256": _sha(BASE_TARGETS),
        "calibration_targets": {
            "mask_opening_nm": targets["mask_opening"],
            "etch_depth_nm": targets["etch_depth"],
        },
        "unit_yield_fraction_endpoints": endpoints,
        "previous_axisymmetric_endpoint": previous,
        "current_axisymmetric_endpoint": current,
        "azimuth_audit_sha256": _sha(azimuth_path),
        "derivation": {
            "method": "axisymmetric-path-secant-plus-independent-fraction-response",
            "response_order": ["mask_opening_nm", "etch_depth_nm"],
            "parameter_order": [
                "effective_mask_crosslinked_growth_fraction",
                "oxide_etch_yield_scale",
            ],
            "midpoint_fraction_response": fraction_response.tolist(),
            "measured_parameter_secant": parameter_secant.tolist(),
            "measured_observable_secant": observed_response.tolist(),
            "inferred_oxide_yield_response": yield_response.tolist(),
            "jacobian": jacobian.tolist(),
            "condition_number": condition_number,
            "requested_response": requested_response.tolist(),
            "parameter_update": update.tolist(),
        },
        "proposed_configuration": {
            "effective_mask_crosslinked_growth_fraction": float(proposed[0]),
            "oxide_etch_yield_scale": float(proposed[1]),
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
        },
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["proposal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint", action="append", required=True)
    parser.add_argument("--coupled")
    parser.add_argument("--axisymmetric-coupled")
    parser.add_argument("--previous-axisymmetric")
    parser.add_argument("--azimuth-audit", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.previous_axisymmetric is not None:
        if args.axisymmetric_coupled is None:
            parser.error("--previous-axisymmetric requires --axisymmetric-coupled")
        payload = derive_axisymmetric_secant(
            args.endpoint, args.previous_axisymmetric, args.axisymmetric_coupled,
            args.azimuth_audit)
    elif args.axisymmetric_coupled is None:
        if args.coupled is None:
            parser.error("--coupled is required for the initial coupled correction")
        payload = derive(args.endpoint, args.coupled, args.azimuth_audit)
    else:
        if args.coupled is None:
            parser.error("--coupled is required for the path-resolved correction")
        payload = derive_path_resolved(
            args.endpoint, args.coupled, args.axisymmetric_coupled,
            args.azimuth_audit)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps(payload["proposed_configuration"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
