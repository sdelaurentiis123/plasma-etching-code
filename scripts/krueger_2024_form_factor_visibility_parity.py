#!/usr/bin/env python3
"""Bounded full-event Warp/float64 visibility parity on the base Krueger checkpoint.

This audit traces one small finite-area/cosine-direction Sobol population through both the fast Warp
classifier and the cell-by-cell float64 reference.  It does not evaluate chemistry, move geometry,
load held-out transfer observations, or select a production ray count.  Exact per-ray event/target
agreement is the predeclared promotion gate; a timeout or mismatch keeps replay hardening in place.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT / "src", ROOT / "scripts"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from krueger_2024_frozen_checkpoint_2x2 import (  # noqa: E402
    EvaluationDeadlineExceeded,
    _hard_deadline,
    _load_source,
    _sha256,
    _write_json_atomic,
)
from petch.boundary_transport_3d import (  # noqa: E402
    _diffuse_form_factor_ray_samples_3d,
    _trace_diffuse_form_factor_events_warp_3d,
    trace_diffuse_form_factor_events_float64_3d,
    trace_diffuse_form_factor_events_warp_cellwise_3d,
)
from petch.feature_step_3d import _face_material_ids, _surface_gas_normals  # noqa: E402
from petch.threed import extract_mesh_3d  # noqa: E402


SCHEMA = "petch.krueger-2024.form-factor-visibility-parity.v1"


def _array_sha256(*arrays):
    digest = sha256()
    for index, supplied in enumerate(arrays):
        array = np.ascontiguousarray(np.asarray(supplied))
        digest.update(f"{index}:{array.dtype.str}:{array.shape}\n".encode("ascii"))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def compare_visibility_events(
        source_face, fast_face, reference, face_area_m2, *,
        fast_termination=None, fast_wrap_count=None):
    """Compare event identity and conservative row distributions without a loose tolerance."""
    source = np.asarray(source_face, dtype=int)
    fast = np.asarray(fast_face, dtype=int)
    area = np.asarray(face_area_m2, dtype=float)
    termination = (
        np.where(fast >= 0, 1, 2).astype(np.int8)
        if fast_termination is None else np.asarray(fast_termination, dtype=np.int8))
    fast_wrap = (
        np.zeros(len(fast), dtype=int)
        if fast_wrap_count is None else np.asarray(fast_wrap_count, dtype=int))
    if (source.ndim != 1 or fast.shape != source.shape
            or termination.shape != source.shape or fast_wrap.shape != source.shape
            or reference.hit_face.shape != source.shape
            or area.ndim != 1 or len(area) == 0 or np.any(area <= 0.0)
            or np.any(source < 0) or np.any(source >= len(area))
            or np.any(~np.isin(termination, (1, 2, 3, 4)))
            or np.any(np.isin(termination, (1, 4)) != (fast >= 0))):
        raise ValueError("incompatible visibility event arrays")
    reference_hit = reference.termination == 1
    fast_hit = termination == 1
    classification_mismatch = termination != reference.termination
    target_mismatch = fast_hit & reference_hit & (fast != reference.hit_face)
    event_mismatch = classification_mismatch | target_mismatch
    counts = np.bincount(source, minlength=len(area))
    if np.any(counts <= 0) or np.any(counts != counts[0]):
        raise ValueError("parity audit requires an equal ray count on every source face")
    rays_per_face = int(counts[0])
    row_total_variation = np.zeros(len(area))
    escape_category = len(area)

    def event_category(face, event):
        category = np.asarray(face, dtype=int).copy()
        category[event == 2] = escape_category
        category[event == 3] = escape_category + 1
        category[event == 4] = escape_category + 2
        return category

    for face_index in range(len(area)):
        selected = source == face_index
        fast_category = event_category(fast[selected], termination[selected])
        reference_category = event_category(
            reference.hit_face[selected], reference.termination[selected])
        fast_count = np.bincount(fast_category, minlength=len(area) + 3)
        reference_count = np.bincount(reference_category, minlength=len(area) + 3)
        row_total_variation[face_index] = (
            0.5 * np.sum(np.abs(fast_count - reference_count)) / rays_per_face)

    wrap_exhaustion = reference.termination == 3
    solid_facing = reference.termination == 4
    return {
        "ray_count": len(source),
        "face_count": len(area),
        "rays_per_face": rays_per_face,
        "reference_hit_count": int(np.sum(reference_hit)),
        "reference_open_escape_count": int(np.sum(reference.termination == 2)),
        "reference_wrap_exhaustion_count": int(np.sum(wrap_exhaustion)),
        "reference_solid_facing_intersection_count": int(np.sum(solid_facing)),
        "fast_hit_count": int(np.sum(fast_hit)),
        "fast_open_escape_count": int(np.sum(termination == 2)),
        "fast_wrap_exhaustion_count": int(np.sum(termination == 3)),
        "fast_solid_facing_intersection_count": int(np.sum(termination == 4)),
        "classification_mismatch_count": int(np.sum(classification_mismatch)),
        "target_face_mismatch_count": int(np.sum(target_mismatch)),
        "any_event_mismatch_count": int(np.sum(event_mismatch)),
        "false_fast_escape_count": int(np.sum((termination == 2) & reference_hit)),
        "false_fast_hit_count": int(np.sum(fast_hit & (~reference_hit))),
        "maximum_source_row_total_variation": float(np.max(row_total_variation)),
        "area_weighted_mean_source_row_total_variation": float(
            np.sum(area * row_total_variation) / np.sum(area)),
        "maximum_reference_wrap_count": int(np.max(reference.wrap_count, initial=0)),
        "maximum_fast_wrap_count": int(np.max(fast_wrap, initial=0)),
        "fast_event_sha256": _array_sha256(fast, termination, fast_wrap),
        "reference_event_sha256": _array_sha256(
            reference.hit_face, reference.termination, reference.wrap_count),
        "all_gates_pass": bool(
            not np.any(event_mismatch) and not np.any(wrap_exhaustion)
            and not np.any(solid_facing)),
    }


def mismatch_geometry_diagnostics(
        source_face, fast_face, reference, faces, centroids, normals,
        material_id, origin, direction, dx, *, fast_termination=None):
    """Classify whether target disagreements are local mesh ownership or remote events."""
    source = np.asarray(source_face, dtype=int)
    fast = np.asarray(fast_face, dtype=int)
    triangles = np.asarray(faces, dtype=int)
    centroid = np.asarray(centroids, dtype=float)
    normal = np.asarray(normals, dtype=float)
    material = np.asarray(material_id, dtype=int)
    ray_direction = np.asarray(direction, dtype=float)
    ray_origin = np.asarray(origin, dtype=float)
    termination = (
        np.where(fast >= 0, 1, 2).astype(np.int8)
        if fast_termination is None else np.asarray(fast_termination, dtype=np.int8))
    both = (termination == 1) & (reference.termination == 1)
    target_mismatch = np.flatnonzero(both & (fast != reference.hit_face))
    shared_vertex = np.zeros(len(target_mismatch), dtype=int)
    centroid_distance_dx = np.zeros(len(target_mismatch))
    normal_alignment = np.zeros(len(target_mismatch))
    same_material = np.zeros(len(target_mismatch), dtype=bool)
    for local, ray_index in enumerate(target_mismatch):
        left = int(fast[ray_index])
        right = int(reference.hit_face[ray_index])
        shared_vertex[local] = len(set(triangles[left]) & set(triangles[right]))
        centroid_distance_dx[local] = np.linalg.norm(
            centroid[left] - centroid[right]) / float(dx)
        normal_alignment[local] = np.dot(normal[left], normal[right])
        same_material[local] = material[left] == material[right]
    fast_hit = np.isin(termination, (1, 4))
    fast_cosine = np.full(len(fast), np.nan)
    fast_cosine[fast_hit] = -np.einsum(
        "ij,ij->i", ray_direction[fast_hit], normal[fast[fast_hit]])
    refusal = np.flatnonzero(np.isin(reference.termination, (3, 4)))

    def quantiles(value):
        if len(value) == 0:
            return {"q05": 0.0, "q50": 0.0, "q95": 0.0, "maximum": 0.0}
        return {
            "q05": float(np.quantile(value, 0.05)),
            "q50": float(np.quantile(value, 0.50)),
            "q95": float(np.quantile(value, 0.95)),
            "maximum": float(np.max(value)),
        }

    selected = np.concatenate((target_mismatch[:32], refusal[:32]))
    selected = np.unique(selected)
    records = []
    for ray_index in selected:
        records.append({
            "ray_index": int(ray_index),
            "source_face": int(source[ray_index]),
            "fast_face": int(fast[ray_index]),
            "reference_face": int(reference.hit_face[ray_index]),
            "reference_termination": int(reference.termination[ray_index]),
            "reference_wrap_count": int(reference.wrap_count[ray_index]),
            "reference_hit_cosine": float(reference.hit_cosine[ray_index]),
            "fast_gas_normal_cosine": (
                None if not np.isfinite(fast_cosine[ray_index])
                else float(fast_cosine[ray_index])),
            "direction": ray_direction[ray_index].tolist(),
            "origin": ray_origin[ray_index].tolist(),
            "reference_hit_position": reference.hit_position[ray_index].tolist(),
        })
    return {
        "target_mismatch_count": len(target_mismatch),
        "shared_edge_count": int(np.sum(shared_vertex >= 2)),
        "shared_vertex_only_count": int(np.sum(shared_vertex == 1)),
        "no_shared_vertex_count": int(np.sum(shared_vertex == 0)),
        "same_material_count": int(np.sum(same_material)),
        "centroid_distance_dx": quantiles(centroid_distance_dx),
        "target_normal_alignment": quantiles(normal_alignment),
        "fast_gas_normal_invalid_count": int(np.sum(fast_cosine[fast_hit] < -2e-6)),
        "first_mismatch_or_refusal_records": records,
    }


def run(args):
    started = perf_counter()
    destination = Path(args.output)
    source = _load_source("r19", args.r19_source)
    vertices, faces, centroids, areas = extract_mesh_3d(
        source["geometry"].phi, source["geometry"].dx)
    normals = np.asarray(_surface_gas_normals(
        vertices, faces, centroids, source["geometry"]), dtype=float)
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    material_id = _face_material_ids(centroids, source["geometry"])
    domain = (np.asarray(source["geometry"].phi.shape) - 1) * source["geometry"].dx
    rays_per_face = int(args.rays_per_face)
    if rays_per_face <= 0 or rays_per_face & (rays_per_face - 1):
        raise ValueError("rays_per_face must be a positive power of two")
    source_face, origin, direction, launch_diagnostics = (
        _diffuse_form_factor_ray_samples_3d(
        vertices, faces, centroids, normals, rays_per_face=rays_per_face,
        seed=int(args.seed), ray_offset=float(args.ray_offset_dx) * source["geometry"].dx,
        source_sampling="triangle_area", return_launch_diagnostics=True))
    payload = {
        "schema": SCHEMA,
        "status": "running",
        "scientific_scope": (
            "bounded full-event visibility parity only; no chemistry, geometry motion, "
            "calibration, or held-out transfer observation"),
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
        },
        "checkpoint": {
            "audit_sha256": _sha256(source["audit_path"]),
            "checkpoint_sha256": _sha256(source["checkpoint_path"]),
            "physical_time_s": source["checkpoint_metadata"].get("physical_time_s"),
        },
        "configuration": {
            "seed": int(args.seed),
            "rays_per_face": rays_per_face,
            "sampling_dimension": 4,
            "source_sampling": "triangle_area",
            "periodic_lateral": True,
            "ray_offset_dx": float(args.ray_offset_dx),
            "maximum_wraps": int(args.maximum_wraps),
            "maximum_wall_s": float(args.maximum_wall_s),
            "device": "cpu",
            "launch_diagnostics": launch_diagnostics,
        },
        "mesh": {
            "face_count": len(faces),
            "vertex_count": len(vertices),
            "mesh_sha256": _array_sha256(
                vertices, faces, centroids, areas, normals),
        },
        "comparison": None,
    }
    _write_json_atomic(destination, payload)
    try:
        with _hard_deadline(float(args.maximum_wall_s)):
            legacy_fast = _trace_diffuse_form_factor_events_warp_3d(
                vertices, faces, origin, direction, domain, True, "cpu")
            candidate = trace_diffuse_form_factor_events_warp_cellwise_3d(
                origin, direction, vertices, faces, normals,
                domain_size=domain, periodic_lateral=True,
                maximum_wraps=int(args.maximum_wraps), device="cpu")
            reference = trace_diffuse_form_factor_events_float64_3d(
                origin, direction, vertices, faces, normals,
                domain_size=domain, periodic_lateral=True,
                maximum_wraps=int(args.maximum_wraps))
        physical_area = areas * source["geometry"].mesh_length_unit_m ** 2
        legacy_comparison = compare_visibility_events(
            source_face, legacy_fast, reference, physical_area)
        legacy_comparison["mismatch_geometry"] = mismatch_geometry_diagnostics(
            source_face, legacy_fast, reference, faces, centroids, normals,
            material_id, origin, direction, source["geometry"].dx)
        candidate_comparison = compare_visibility_events(
            source_face, candidate.hit_face, reference, physical_area,
            fast_termination=candidate.termination,
            fast_wrap_count=candidate.wrap_count)
        candidate_comparison["mismatch_geometry"] = mismatch_geometry_diagnostics(
            source_face, candidate.hit_face, reference, faces, centroids, normals,
            material_id, origin, direction, source["geometry"].dx,
            fast_termination=candidate.termination)
        payload["comparison"] = {
            "legacy_apply_bc": legacy_comparison,
            "cellwise_candidate": candidate_comparison,
            "promotion_gate_pass": candidate_comparison["all_gates_pass"],
        }
        payload["status"] = (
            "pass" if candidate_comparison["all_gates_pass"] else "parity_failure")
    except EvaluationDeadlineExceeded as error:
        payload["status"] = "bounded_timeout"
        payload["refusal"] = {
            "type": type(error).__name__, "reason": str(error),
            "physics_conclusion_permitted": False,
        }
    except (RuntimeError, ValueError) as error:
        payload["status"] = "authority_refusal"
        payload["refusal"] = {
            "type": type(error).__name__, "reason": str(error),
            "physics_conclusion_permitted": False,
        }
    payload["total_wall_time_s"] = float(perf_counter() - started)
    _write_json_atomic(destination, payload)
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--r19-source",
        default=ROOT / "results" / "krueger_2024_r19_response_check" / "remote_artifacts")
    parser.add_argument(
        "--output",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "form_factor_visibility_parity" / "audit.json"))
    parser.add_argument("--seed", type=int, default=10241)
    parser.add_argument("--rays-per-face", type=int, default=8)
    parser.add_argument("--ray-offset-dx", type=float, default=1e-3)
    parser.add_argument("--maximum-wraps", type=int, default=1024)
    parser.add_argument("--maximum-wall-s", type=float, default=180.0)
    return parser.parse_args(argv)


if __name__ == "__main__":
    result = run(parse_args())
    comparisons = result.get("comparison") or {}
    print(json.dumps({
        "status": result["status"],
        "total_wall_time_s": result["total_wall_time_s"],
        "promotion_gate_pass": comparisons.get("promotion_gate_pass"),
        "cellwise_candidate": {
            key: comparisons.get("cellwise_candidate", {}).get(key)
            for key in (
                "ray_count", "any_event_mismatch_count",
                "classification_mismatch_count", "target_face_mismatch_count",
                "reference_wrap_exhaustion_count",
                "reference_solid_facing_intersection_count",
                "maximum_reference_wrap_count",
                "area_weighted_mean_source_row_total_variation",
                "all_gates_pass")
        },
        "legacy_apply_bc": {
            key: comparisons.get("legacy_apply_bc", {}).get(key)
            for key in (
                "any_event_mismatch_count", "classification_mismatch_count",
                "target_face_mismatch_count",
                "area_weighted_mean_source_row_total_variation")
        },
    }, indent=2, sort_keys=True))
