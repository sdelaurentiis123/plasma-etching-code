#!/usr/bin/env python3
"""Recover the physical scale of the clock-gate STL from its source GDS cell."""
from __future__ import annotations

import argparse
from functools import reduce
from hashlib import sha256
import itertools
import json
from math import gcd
from pathlib import Path
import re
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MAIN_REPO = ROOT.parents[1]
if not (MAIN_REPO / "src" / "petch").is_dir():
    MAIN_REPO = MAIN_REPO / "plasma-etching-code"
sys.path.insert(0, str(MAIN_REPO / "src"))

from petch.gds_import import read_gds  # noqa: E402
from petch.stl_import import drop_degenerate_faces, read_stl  # noqa: E402


STL = ROOT / "raw" / "CLKGATE X1 25x Reduced 3D.STL"
OUTPUT = ROOT / "results" / "source_scale_audit.json"
HEADER = re.compile(
    r"stlbn (?P<cell>.+?)_25x_reduced_L(?P<layer>[0-9]+)_3D",
    re.IGNORECASE,
)
REFERENCE_URL = (
    "https://raw.githubusercontent.com/mflowgen/freepdk-45nm/"
    "master/stdcells.gds"
)
LAYER_MAP_URL = (
    "https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts/"
    "blob/master/flow/platforms/nangate45/FreePDK45.lyt"
)
NANOSCRIBE_UNITS_URL = (
    "https://www.nanoscribe.com/en/contact-support/support/cad-model-creation/"
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _unique_rows(values, decimals):
    return np.unique(np.round(np.asarray(values, dtype=float), decimals), axis=0)


def _directed_nearest(first, second):
    distance = np.linalg.norm(
        np.asarray(first)[:, None, :] - np.asarray(second)[None, :, :], axis=2)
    return np.min(distance, axis=1)


def _minimum_point_to_segments(points, polygon):
    points = np.asarray(points, dtype=float)
    start = np.asarray(polygon, dtype=float)
    end = np.roll(start, -1, axis=0)
    direction = end - start
    length_squared = np.sum(direction ** 2, axis=1)
    if np.any(length_squared == 0.0):
        raise ValueError("source polygon contains a zero-length edge")
    offset = points[:, None, :] - start[None, :, :]
    fraction = np.sum(
        offset * direction[None, :, :], axis=2) / length_squared[None, :]
    fraction = np.clip(fraction, 0.0, 1.0)
    nearest = start[None, :, :] + fraction[:, :, None] * direction[None, :, :]
    return float(np.min(np.linalg.norm(points[:, None, :] - nearest, axis=2)))


def _minimum_polygon_gap(polygons):
    distances = []
    for first_index, first in enumerate(polygons):
        for second in polygons[first_index + 1:]:
            distances.append(min(
                _minimum_point_to_segments(first, second),
                _minimum_point_to_segments(second, first),
            ))
    return min(value for value in distances if value > 0.0)


def _fit_similarity(source_um, target_file_units):
    """Fit every axis permutation/reflection with one isotropic scale."""
    candidates = []
    target_extent = np.ptp(target_file_units, axis=0)
    for permutation in itertools.permutations((0, 1)):
        permuted = source_um[:, permutation]
        source_extent = np.ptp(permuted, axis=0)
        axis_scale = target_extent / source_extent
        scale = float(np.mean(axis_scale))
        for signs in itertools.product((-1.0, 1.0), repeat=2):
            oriented = permuted * np.asarray(signs)
            offset = (
                np.min(target_file_units, axis=0)
                - np.min(oriented * scale, axis=0)
            )
            predicted = oriented * scale + offset
            forward = _directed_nearest(predicted, target_file_units)
            reverse = _directed_nearest(target_file_units, predicted)
            candidates.append({
                "axis_permutation": list(permutation),
                "axis_sign": list(signs),
                "file_units_per_source_um": scale,
                "axis_scale_file_units_per_source_um": axis_scale.tolist(),
                "translation_file_units": offset.tolist(),
                "maximum_bidirectional_vertex_residual_file_units": float(max(
                    np.max(forward), np.max(reverse))),
                "rms_source_to_stl_vertex_residual_file_units": float(
                    np.sqrt(np.mean(forward ** 2))),
            })
    return min(
        candidates,
        key=lambda item: item[
            "maximum_bidirectional_vertex_residual_file_units"],
    )


def build(reference_gds):
    header = STL.read_bytes()[:80].decode("ascii", errors="replace").rstrip()
    match = HEADER.fullmatch(header)
    if match is None:
        raise ValueError(f"unrecognized STL provenance header: {header!r}")
    cell_name = match.group("cell")
    layer = int(match.group("layer"))

    library = read_gds(reference_gds)
    if cell_name not in library.cells:
        raise ValueError(f"source GDS does not contain {cell_name!r}")
    polygons = [
        polygon for polygon in library.cells[cell_name].polygons
        if polygon.layer == layer and polygon.datatype == 0
    ]
    if not polygons:
        raise ValueError(f"source GDS has no layer {layer}/0 polygons")
    source_polygons_um = [
        polygon.vertices_db * library.database_unit_m * 1.0e6
        for polygon in polygons
    ]
    source_db = np.concatenate([polygon.vertices_db for polygon in polygons])
    source_um = np.concatenate(source_polygons_um)
    source_unique_um = _unique_rows(source_um, 9)

    mesh, repair = drop_degenerate_faces(read_stl(STL))
    vertices = np.asarray(mesh.vertices, dtype=float)
    unique_counts = [len(np.unique(vertices[:, axis])) for axis in range(3)]
    extrusion_axes = [
        axis for axis, count in enumerate(unique_counts) if count == 2]
    if len(extrusion_axes) != 1:
        raise ValueError(f"STL does not have one extrusion axis: {unique_counts}")
    extrusion_axis = extrusion_axes[0]
    planar_axes = [axis for axis in range(3) if axis != extrusion_axis]
    target_unique = _unique_rows(vertices[:, planar_axes], 9)
    fit = _fit_similarity(source_unique_um, target_unique)

    # The source geometry sits on the Nangate 5 nm manufacturing grid.  Use
    # integer database coordinates so the result is exact rather than inferred
    # from binary floating-point differences.
    coordinate = np.unique(source_db.reshape(-1))
    differences = np.diff(coordinate)
    grid_db = reduce(gcd, (int(value) for value in differences if value), 0)
    source_grid_um = grid_db * library.database_unit_m * 1.0e6
    edge_length_um = np.concatenate([
        np.linalg.norm(np.roll(vertices, -1, axis=0) - vertices, axis=1)
        for vertices in source_polygons_um
    ])
    minimum_edge_um = float(np.min(edge_length_um[edge_length_um > 0.0]))
    minimum_gap_um = float(_minimum_polygon_gap(source_polygons_um))

    file_unit_um = 1000.0
    physical_enlargement = (
        fit["file_units_per_source_um"] * file_unit_um)
    bounds = np.asarray(mesh.bounds, dtype=float)
    extent = bounds[1] - bounds[0]
    physical_extent_um = extent * file_unit_um
    source_extent_um = np.ptp(source_unique_um, axis=0)
    maximum_residual_um = (
        fit["maximum_bidirectional_vertex_residual_file_units"]
        * file_unit_um)

    if not np.isclose(physical_enlargement, 40.0, rtol=0.0, atol=1e-4):
        raise RuntimeError("source-to-STL enlargement is not the expected exact 40x")
    if len(source_unique_um) != len(target_unique) or maximum_residual_um > 1e-3:
        raise RuntimeError("source GDS and STL vertices do not match at sub-nm tolerance")

    return {
        "schema": "petch.partner.clockgate-source-scale-audit.v1",
        "source_stl": {
            "path": str(STL.relative_to(ROOT)),
            "sha256": _hash(STL),
            "binary_header": header,
            "repaired_degenerate_face_count": int(repair["removed_face_count"]),
        },
        "reference_gds": {
            "download_url": REFERENCE_URL,
            "local_sha256": _hash(reference_gds),
            "library_name": library.name,
            "cell_name": cell_name,
            "layer": layer,
            "datatype": 0,
            "polygon_count": len(polygons),
            "unique_vertex_count": len(source_unique_um),
            "database_unit_m": library.database_unit_m,
            "source_footprint_extent_um": source_extent_um.tolist(),
            "source_grid_um": float(source_grid_um),
            "layer_semantics": "FreePDK45 Metal 1",
            "layer_map_source": LAYER_MAP_URL,
        },
        "exact_geometry_match": {
            **fit,
            "stl_unique_planar_vertex_count": len(target_unique),
            "maximum_bidirectional_vertex_residual_um": maximum_residual_um,
            "physical_enlargement_from_source_gds": physical_enlargement,
            "match_grade": "exact within binary-STL float32 roundoff",
        },
        "recovered_units": {
            "stl_file_unit": "millimeter",
            "micrometers_per_stl_file_unit": file_unit_um,
            "evidence_class": "strongly inferred by exact public-source match and print workflow",
            "declared_in_stl_standard": False,
            "nanoscribe_import_convention_source": NANOSCRIBE_UNITS_URL,
            "nanoscribe_import_note": (
                "DeScribe interprets raw STL units as micrometers, so this "
                "millimeter-authored STL requires 1000x import scaling"
            ),
            "filename_25x_interpretation": (
                "Nanoscribe 25x objective or process preset; not geometry scale"
            ),
            "filename_25x_cannot_be_geometry_scale": (
                "the geometry itself is exactly 40x the source Metal-1 polygons"
            ),
        },
        "physical_geometry": {
            "extrusion_axis": extrusion_axis,
            "planar_axes": planar_axes,
            "extent_um_in_stl_axis_order": physical_extent_um.tolist(),
            "mask_height_um": float(physical_extent_um[extrusion_axis]),
            "footprint_extent_um": physical_extent_um[planar_axes].tolist(),
            "source_grid_after_enlargement_um": float(
                source_grid_um * physical_enlargement),
            "minimum_source_polygon_gap_um": minimum_gap_um,
            "minimum_physical_opening_um": float(
                minimum_gap_um * physical_enlargement),
            "minimum_source_rectilinear_edge_um": minimum_edge_um,
            "minimum_physical_rectilinear_edge_um": float(
                minimum_edge_um * physical_enlargement),
            "minimum_opening_mask_aspect_ratio": float(
                physical_extent_um[extrusion_axis]
                / (minimum_gap_um * physical_enlargement)),
            "solid_component_count": 10,
        },
        "process_interpretation": {
            "stl_solids_are": "40x enlarged FreePDK45 Metal-1 polygons",
            "working_mask_polarity": "printed solids protect silicon",
            "working_target": "raised silicon replicas of the Metal-1 wiring",
            "alternate_polarity": (
                "treat solids as openings only if the desired silicon result is trenches"
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-gds", type=Path, required=True)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = json.dumps(
        build(args.reference_gds), indent=2, sort_keys=True) + "\n"
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    elif not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("source-scale audit is stale")
    else:
        print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
