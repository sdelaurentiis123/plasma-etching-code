#!/usr/bin/env python3
"""Geometry-only physical-patch support sensitivity on the sealed Krueger base mesh.

This audit performs no transport, chemistry, calibration, or profile motion.  It records how the
predeclared physical-patch mean-support rule partitions the already-sealed checkpoint at every
required patch scale.  The conservative integrated gate continues to contain every patch.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT / "src", ROOT / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import krueger_2024_replicated_form_factor_closure as closure  # noqa: E402
from petch.surface_patch_convergence_3d import (  # noqa: E402
    aggregate_surface_field_on_physical_patches_3d,
)


SCHEMA = "petch.krueger-2024.physical-patch-support-sensitivity.v1"
SOURCE_PATHS = (
    "scripts/krueger_2024_patch_support_sensitivity.py",
    "scripts/krueger_2024_replicated_form_factor_closure.py",
    "src/petch/surface_patch_convergence_3d.py",
)


def _array_digest(value):
    array = np.ascontiguousarray(np.asarray(value))
    digest = sha256()
    digest.update(array.dtype.str.encode("ascii") + b"\0")
    digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _mesh_identity(direct):
    receipt = {
        name: _array_digest(direct[name])
        for name in (
            "verts", "faces", "face_area_m2", "gas_normals", "face_material_id")
    }
    encoded = json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return {"arrays": receipt, "combined_sha256": sha256(encoded).hexdigest()}


def _support_table(direct, source):
    common = {
        "face_field": np.ones(len(direct["faces"]), dtype=float),
        "face_area_m2": direct["face_area_m2"],
        "verts": direct["verts"],
        "faces": direct["faces"],
        "face_gas_normals": direct["gas_normals"],
        "face_material_id": direct["face_material_id"],
        "mesh_length_unit_m": source["geometry"].mesh_length_unit_m,
        **closure._physical_patch_aggregation_kwargs(source),
    }
    output = []
    for scale in closure.PATCH_SCALES_M:
        receipt = aggregate_surface_field_on_physical_patches_3d(
            patch_scale_m=scale, **common)
        rows = []
        for threshold in closure.PATCH_SUPPORT_SENSITIVITY_THRESHOLDS:
            eligible = receipt.patch_projected_support_fraction >= threshold
            excluded = ~eligible
            rows.append({
                "minimum_mean_support_fraction": float(threshold),
                "eligible_mean_patch_count": int(np.count_nonzero(eligible)),
                "excluded_mean_patch_count": int(np.count_nonzero(excluded)),
                "excluded_surface_area_m2": float(np.sum(
                    receipt.patch_area_m2[excluded])),
                "excluded_surface_area_fraction": float(np.sum(
                    receipt.patch_area_m2[excluded]) / np.sum(receipt.patch_area_m2)),
                "excluded_projected_support_area_m2": float(np.sum(
                    receipt.patch_projected_support_area_m2[excluded])),
                "excluded_projected_support_fraction": float(np.sum(
                    receipt.patch_projected_support_area_m2[excluded])
                    / np.sum(receipt.patch_projected_support_area_m2)),
            })
        output.append({
            "patch_scale_m": float(scale),
            "patch_count": len(receipt.patch_key),
            "scheme_sha256": receipt.scheme_sha256,
            "total_surface_area_m2": float(np.sum(receipt.patch_area_m2)),
            "total_projected_support_area_m2": float(np.sum(
                receipt.patch_projected_support_area_m2)),
            "represented_nominal_projected_area_m2": {
                "minimum": float(np.min(receipt.patch_nominal_projected_area_m2)),
                "maximum": float(np.max(receipt.patch_nominal_projected_area_m2)),
                "total": float(np.sum(receipt.patch_nominal_projected_area_m2)),
            },
            "projected_support_fraction": {
                "minimum": float(np.min(receipt.patch_projected_support_fraction)),
                "median": float(np.median(receipt.patch_projected_support_fraction)),
                "maximum": float(np.max(receipt.patch_projected_support_fraction)),
            },
            "surface_area_inventory_relative_error": float(abs(
                np.sum(receipt.patch_area_m2) - np.sum(direct["face_area_m2"]))
                / np.sum(direct["face_area_m2"])),
            "thresholds": rows,
        })
    return output


def run(source_path, output_path):
    source = closure._load_sealed_base_source(source_path)
    direct = closure._direct_geometry_from_checkpoint(source)
    operator = closure._physical_patch_operator_contract(source)
    table = _support_table(direct, source)
    payload = {
        "schema": SCHEMA,
        "status": "complete",
        "scope": "sealed base checkpoint geometry only; no transport/chemistry/profile/heldout",
        "checkpoint": {
            "path_name": Path(source["checkpoint_path"]).name,
            "sha256": closure._sha256(source["checkpoint_path"]),
            "source_audit_sha256": closure._sha256(source["audit_path"]),
        },
        "source_manifest": closure._hash_manifest(SOURCE_PATHS),
        "mesh_identity": _mesh_identity(direct),
        "physical_patch_operator": operator,
        "primary_threshold_binding": {
            "minimum_mean_support_fraction": float(
                closure.DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
            "selection": "predeclared geometry-independent round fraction",
            "selected_before_any_transport_or_chemistry_score": True,
            "old_selected_source_allocator_threshold_0p10_is_executable_input": False,
            "integrated_all_patch_gate_unchanged": True,
        },
        "sensitivity_thresholds": [
            float(value) for value in closure.PATCH_SUPPORT_SENSITIVITY_THRESHOLDS],
        "patch_scales": table,
    }
    output_path = Path(output_path)
    closure._write_json_atomic(output_path, payload)
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path,
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "remote_artifacts"))
    parser.add_argument(
        "--output", type=Path,
        default=(ROOT / "results" / "curated"
                 / "krueger_2024_physical_patch_support_sensitivity_2026-07-17.json"))
    return parser.parse_args(argv)


if __name__ == "__main__":
    _args = parse_args()
    _payload = run(_args.source, _args.output)
    print(json.dumps({
        "status": _payload["status"],
        "output": str(_args.output),
        "mesh_sha256": _payload["mesh_identity"]["combined_sha256"],
    }, indent=2, sort_keys=True))
