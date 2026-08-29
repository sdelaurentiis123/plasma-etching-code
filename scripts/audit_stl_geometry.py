#!/usr/bin/env python3
"""Audit and safely repair a supplied STL without assigning undocumented units."""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.stl_import import (
    StlMesh,
    diagnose_mesh,
    drop_degenerate_faces,
    read_stl,
    write_stl,
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _diagnostics_payload(mesh):
    item = diagnose_mesh(mesh)
    return {
        "n_vertices": item.n_vertices,
        "n_faces": item.n_faces,
        "n_degenerate_faces": item.n_degenerate_faces,
        "n_boundary_edges": item.n_boundary_edges,
        "n_nonmanifold_edges": item.n_nonmanifold_edges,
        "consistently_oriented": item.consistently_oriented,
        "signed_volume_file_units_cubed": item.signed_volume,
        "is_watertight": item.is_watertight,
        "outward_oriented": item.outward_oriented,
        "failure_reason": item.failure_reason(),
    }


def _submesh(mesh, face_indices):
    faces = np.asarray(mesh.faces, dtype=int)[np.asarray(face_indices, dtype=int)]
    used, inverse = np.unique(faces.reshape(-1), return_inverse=True)
    return StlMesh(
        np.asarray(mesh.vertices)[used], inverse.reshape(-1, 3),
        None if mesh.file_normals is None
        else np.asarray(mesh.file_normals)[np.asarray(face_indices, dtype=int)],
    )


def _components(mesh):
    edges = defaultdict(list)
    for face_index, face in enumerate(np.asarray(mesh.faces, dtype=int)):
        for first, second in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edges[tuple(sorted((int(first), int(second))))].append(face_index)
    adjacency = [set() for _ in mesh.faces]
    for incident in edges.values():
        for first in incident:
            adjacency[first].update(second for second in incident if second != first)
    unseen = set(range(len(mesh.faces)))
    output = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        unseen.remove(seed)
        indices = []
        while stack:
            current = stack.pop()
            indices.append(current)
            neighbors = adjacency[current] & unseen
            stack.extend(sorted(neighbors, reverse=True))
            unseen.difference_update(neighbors)
        output.append(_submesh(mesh, indices))
    return sorted(output, key=lambda item: abs(item.signed_volume), reverse=True)


def _planar_extrusion_payload(*, extent, components, candidates):
    payload = {
        "candidate_axes": list(candidates),
        "unique_candidate": len(candidates) == 1,
        "candidate_axis": candidates[0] if len(candidates) == 1 else None,
        "extrusion_thickness_file_units": None,
        "projected_solid_area_file_units_squared": None,
        "projected_bounding_box_area_file_units_squared": None,
        "projected_fill_fraction": None,
        "component_projected_area_file_units_squared": None,
    }
    if len(candidates) != 1:
        return payload
    axis = int(candidates[0])
    thickness = float(extent[axis])
    planar_axes = [index for index in range(3) if index != axis]
    bounding_area = float(np.prod(np.asarray(extent)[planar_axes]))
    component_areas = [
        abs(float(item["signed_volume_file_units_cubed"])) / thickness
        for item in components
    ]
    projected_area = float(sum(component_areas))
    payload.update({
        "extrusion_thickness_file_units": thickness,
        "projected_solid_area_file_units_squared": projected_area,
        "projected_bounding_box_area_file_units_squared": bounding_area,
        "projected_fill_fraction": projected_area / bounding_area,
        "component_projected_area_file_units_squared": component_areas,
    })
    return payload


def build(source):
    source = Path(source)
    raw = read_stl(source)
    repaired, repair = drop_degenerate_faces(raw)
    lower, upper = (np.asarray(item, dtype=float) for item in repaired.bounds)
    extent = upper - lower
    triangles = np.asarray(repaired.triangles, dtype=float)
    edges = np.concatenate((
        triangles[:, 1] - triangles[:, 0],
        triangles[:, 2] - triangles[:, 1],
        triangles[:, 0] - triangles[:, 2],
    ))
    lengths = np.linalg.norm(edges, axis=1)
    lengths = lengths[lengths > 0.0]
    components = []
    for index, component in enumerate(_components(repaired)):
        component_lower, component_upper = (
            np.asarray(item, dtype=float) for item in component.bounds)
        components.append({
            "component": index,
            "face_count": len(component.faces),
            "vertex_count": len(component.vertices),
            "bounds_file_units": [
                component_lower.tolist(), component_upper.tolist(),
            ],
            "extent_file_units": (component_upper - component_lower).tolist(),
            "signed_volume_file_units_cubed": component.signed_volume,
            "diagnostics": _diagnostics_payload(component),
        })
    unique_per_axis = [
        len(np.unique(np.asarray(repaired.vertices)[:, axis]))
        for axis in range(3)
    ]
    extrusion_candidates = []
    for axis, unique_count in enumerate(unique_per_axis):
        if unique_count != 2:
            continue
        if all(
            np.isclose(item["extent_file_units"][axis], extent[axis])
            for item in components
        ):
            extrusion_candidates.append(axis)
    planar = _planar_extrusion_payload(
        extent=extent, components=components, candidates=extrusion_candidates)
    planar["unique_coordinate_count_by_axis"] = unique_per_axis
    return repaired, {
        "schema": "petch.stl-geometry-audit.v2",
        "source": {
            "filename": source.name,
            "size_bytes": source.stat().st_size,
            "sha256": _hash(source),
        },
        "units": {
            "declared_by_stl": None,
            "physical_scale_identified": False,
            "warning": (
                "STL stores coordinates but no units; no physical depth, CD, "
                "or etch time may be inferred until the sender identifies scale"
            ),
        },
        "raw_diagnostics": _diagnostics_payload(raw),
        "repair": repair,
        "clean_diagnostics": _diagnostics_payload(repaired),
        "bounds_file_units": [lower.tolist(), upper.tolist()],
        "extent_file_units": extent.tolist(),
        "surface_area_file_units_squared": float(np.sum(repaired.face_areas)),
        "volume_file_units_cubed": repaired.signed_volume,
        "connected_component_count": len(components),
        "components": components,
        "edge_length_file_units": {
            "minimum": float(np.min(lengths)),
            "p01": float(np.percentile(lengths, 1)),
            "median": float(np.median(lengths)),
            "p99": float(np.percentile(lengths, 99)),
            "maximum": float(np.max(lengths)),
        },
        "planar_extrusion": planar,
        "simulation_readiness": {
            "geometry_topology_ready": (
                diagnose_mesh(repaired).failure_reason() is None),
            "physical_scale_ready": False,
            "mask_or_void_polarity_ready": False,
            "material_stack_ready": False,
            "reactor_boundary_ready": False,
            "surface_mechanism_ready": False,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--repaired-output", type=Path, required=True)
    args = parser.parse_args()
    repaired, audit = build(args.source)
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.repaired_output.parent.mkdir(parents=True, exist_ok=True)
    write_stl(args.repaired_output, repaired, binary=True, name="petch_repaired")
    audit["repaired"] = {
        "filename": args.repaired_output.name,
        "sha256": _hash(args.repaired_output),
        "size_bytes": args.repaired_output.stat().st_size,
    }
    args.audit_output.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.audit_output)


if __name__ == "__main__":
    main()
