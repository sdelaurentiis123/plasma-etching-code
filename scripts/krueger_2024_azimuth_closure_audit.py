#!/usr/bin/env python3
"""Bounded endpoint audit of the Krüger IEAD's missing 3-D azimuthal closure.

The publication supplies energy and polar angle at the wafer but no azimuthal distribution.  The
reactor and nominal trench have no preferred wafer-plane direction, so the declared 3-D closure is
uniform azimuth.  This audit freezes one evolved surface state and compares the historical
single-plane quadrature with uniform 8-, 16-, and 32-point azimuthal rings.  It changes no geometry and
reads no held-out profile.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from krueger_2024_endpoint_operator_audit import (
    GATES,
    _evaluate,
    _relative_operator_error,
    _sha256,
    _surface_flux_by_species,
)
from krueger_2024_trench_pilot import _load_checkpoint


def _write_json(path, payload):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def run(args):
    source = Path(args.source)
    audit_path = source / "audit.json"
    checkpoint_path = source / "checkpoint.npz"
    pilot = json.loads(audit_path.read_text(encoding="utf-8"))
    if pilot.get("status") != "complete":
        raise RuntimeError("azimuth audit requires a completed 60 s base endpoint")
    geometry, state, fingerprint, _ = _load_checkpoint(checkpoint_path)
    base_config = dict(pilot["configuration"])

    configs = {}
    plane = dict(base_config)
    plane.pop("ion_azimuthal_closure", None)
    plane.pop("ion_azimuthal_order", None)
    configs["single_published_plane"] = plane
    for order in (8, 16, 32):
        config = dict(base_config)
        config["ion_azimuthal_closure"] = "axisymmetric_uniform"
        config["ion_azimuthal_order"] = order
        configs[f"axisymmetric_uniform_{order}"] = config

    evaluated = {}
    variants = {}
    for name, config in configs.items():
        result, boundary, wall = _evaluate(
            geometry, state, fingerprint,
            boundary_mode="angular_8x16", ion_bins=(250.0, 0.25),
            face_points=3, pilot_config=config,
            radiosity_rays=int(args.radiosity_rays), seed=int(args.seed),
            ballistic_transport="face_gather",
            transport_device=str(args.transport_device))
        evaluated[name] = result
        variants[name] = {
            "wall_time_s": float(wall),
            "ion_quadrature_node_count": int(boundary.get("ions").weight.size),
            "ion_azimuthal_closure": boundary.provenance.get(
                "ion_azimuthal_closure"),
            "ion_azimuthal_order": boundary.provenance.get("ion_azimuthal_order"),
        }

    reference = evaluated["axisymmetric_uniform_32"]
    reference_index = np.asarray(reference.active_face_index, dtype=int)
    reference_centroid = np.asarray(reference.active_face_centroid, dtype=float)
    reference_area = np.asarray(reference.active_face_area, dtype=float)
    for name, result in evaluated.items():
        active_index = np.asarray(result.active_face_index, dtype=int)
        material = np.asarray(result.face_material_id, dtype=int)[active_index]
        velocity = np.asarray(result.face_velocity_mesh_units_s, dtype=float)[active_index]
        area = np.asarray(result.active_face_area, dtype=float)
        variants[name]["material_profile_rate"] = {
            str(int(material_id)): {
                "active_area_mesh_units2": float(np.sum(area[material == material_id])),
                "net_volume_rate_mesh_units3_s": float(np.sum(
                    velocity[material == material_id] * area[material == material_id])),
                "gross_volume_rate_mesh_units3_s": float(np.sum(
                    np.abs(velocity[material == material_id]) * area[material == material_id])),
                "area_weighted_mean_velocity_mesh_units_s": float(np.sum(
                    velocity[material == material_id] * area[material == material_id])
                    / np.sum(area[material == material_id])),
            }
            for material_id in np.unique(material)
        }
    comparisons = {}
    for label, candidate_name in (
            ("historical_plane_to_axisymmetric_32", "single_published_plane"),
            ("axisymmetric_8_to_32", "axisymmetric_uniform_8"),
            ("axisymmetric_16_to_32", "axisymmetric_uniform_16")):
        candidate = evaluated[candidate_name]
        if (not np.array_equal(candidate.active_face_index, reference_index)
                or not np.array_equal(candidate.active_face_centroid, reference_centroid)):
            raise RuntimeError("azimuth variants do not share one exact active-face mesh")
        comparison = {
            "instantaneous_profile_velocity": _relative_operator_error(
                np.asarray(candidate.face_velocity_mesh_units_s)[reference_index],
                np.asarray(reference.face_velocity_mesh_units_s)[reference_index],
                reference_area),
            "species_surface_flux": {},
        }
        candidate_flux = _surface_flux_by_species(candidate)
        reference_flux = _surface_flux_by_species(reference)
        if set(candidate_flux) != set(reference_flux):
            raise RuntimeError("azimuth variants expose different species")
        comparison["species_surface_flux"] = {
            name: _relative_operator_error(
                candidate_flux[name][reference_index],
                reference_flux[name][reference_index], reference_area)
            for name in sorted(reference_flux)
        }
        comparisons[label] = comparison

    refinement = comparisons["axisymmetric_16_to_32"][
        "instantaneous_profile_velocity"]
    gate_pass = all(refinement[name] <= limit for name, limit in GATES.items())
    payload = {
        "schema": "petch.krueger-2024.azimuth-closure-audit.v3",
        "scientific_status": (
            "bounded fixed-endpoint numerical/closure audit; not held-out validation"),
        "source_audit_sha256": _sha256(audit_path),
        "source_checkpoint_sha256": _sha256(checkpoint_path),
        "declared_physical_closure": (
            "uniform azimuth from reactor/trench rotational symmetry; publication supplies "
            "only energy and polar angle"),
        "variants": variants,
        "comparisons": comparisons,
        "production_azimuthal_order": 16,
        "reference_azimuthal_order": 32,
        "axisymmetric_order_refinement_gates": dict(GATES),
        "axisymmetric_order_refinement_pass": bool(gate_pass),
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["audit_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_json(output, payload)
    print(json.dumps({
        "axisymmetric_order_refinement_pass": gate_pass,
        "historical_plane_velocity_error": comparisons[
            "historical_plane_to_axisymmetric_32"]["instantaneous_profile_velocity"],
        "axisymmetric_8_to_32_velocity_error": comparisons[
            "axisymmetric_8_to_32"]["instantaneous_profile_velocity"],
        "axisymmetric_16_to_32_velocity_error": refinement,
    }, indent=2, sort_keys=True))
    return payload


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--radiosity-rays", type=int, default=8)
    parser.add_argument("--seed", type=int, default=5241)
    parser.add_argument("--transport-device", default="cpu")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
