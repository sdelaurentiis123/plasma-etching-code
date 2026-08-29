#!/usr/bin/env python3
"""Checksum and decode Freddie Zhu's supplied Oxford mask-layout GDSII file."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.gds_import import GdsArrayReference, GdsReference, read_gds


DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
SOURCE = DATA / "layout" / "gds_sem_image_450nm.gds"
OUTPUT = DATA / "layout" / "gds_geometry_audit.json"
SQUARE_NAME = re.compile(r"M\d+_S(?P<width>[0-9]+)$")


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _bounds(polygons):
    points = np.concatenate([item.vertices_db for item in polygons], axis=0)
    lower = points.min(axis=0)
    upper = points.max(axis=0)
    return np.asarray([lower[0], lower[1], upper[0], upper[1]], dtype=float)


def _transformed_bounds(reference: GdsReference, polygons):
    points = np.concatenate(
        [reference.transform(item.vertices_db) for item in polygons], axis=0)
    lower = points.min(axis=0)
    upper = points.max(axis=0)
    return np.asarray([lower[0], lower[1], upper[0], upper[1]], dtype=float)


def _array_bounds(array: GdsArrayReference, polygons):
    corners = (
        array.origin_db,
        array.origin_db + (array.columns - 1) * array.column_vector_db,
        array.origin_db + (array.rows - 1) * array.row_vector_db,
        array.origin_db
        + (array.columns - 1) * array.column_vector_db
        + (array.rows - 1) * array.row_vector_db,
    )
    return np.asarray([
        _transformed_bounds(array.as_reference(origin), polygons)
        for origin in corners
    ])


def _top_bounds(library, top_cell):
    cell = library.cells[top_cell]
    bounds = []
    if cell.polygons:
        bounds.append(_bounds(cell.polygons))
    for reference in cell.references:
        child = library.cells[reference.cell_name]
        if child.references or child.arrays:
            raise ValueError("nested transformed GDS bounds are not implemented")
        bounds.append(_transformed_bounds(reference, child.polygons))
    for array in cell.arrays:
        child = library.cells[array.cell_name]
        if child.references or child.arrays:
            raise ValueError("nested transformed GDS bounds are not implemented")
        bounds.extend(_array_bounds(array, child.polygons))
    values = np.asarray(bounds)
    return np.asarray([
        np.min(values[:, 0]), np.min(values[:, 1]),
        np.max(values[:, 2]), np.max(values[:, 3]),
    ])


def build():
    library = read_gds(SOURCE)
    if not library.geometry_complete:
        raise RuntimeError("supplied GDS contains unsupported geometric elements")
    if library.top_cells != ("main",):
        raise RuntimeError(f"expected one 'main' top cell, got {library.top_cells}")
    scale_nm = library.database_unit_m * 1.0e9
    top = library.cells["main"]
    square_cells = []
    special_cells = []
    for name, cell in sorted(library.cells.items()):
        if name == "main":
            continue
        bounds = _bounds(cell.polygons) * scale_nm
        polygon_payload = [
            {
                "layer": int(polygon.layer),
                "datatype": int(polygon.datatype),
                "vertices_nm": (
                    np.asarray(polygon.vertices_db, dtype=float) * scale_nm
                ).tolist(),
            }
            for polygon in cell.polygons
        ]
        match = SQUARE_NAME.fullmatch(name)
        item = {
            "cell_name": name,
            "bounds_nm": bounds.tolist(),
            "width_x_nm": float(bounds[2] - bounds[0]),
            "width_y_nm": float(bounds[3] - bounds[1]),
            "polygons": polygon_payload,
        }
        if match:
            declared_width = float(match.group("width"))
            if not np.allclose(
                [item["width_x_nm"], item["width_y_nm"]],
                declared_width,
                rtol=0.0,
                atol=1.0e-12,
            ):
                raise RuntimeError(f"cell name/geometry mismatch for {name}")
            item["declared_width_nm"] = declared_width
            square_cells.append(item)
        else:
            special_cells.append(item)

    steps = {
        tuple(np.round(np.concatenate((
            array.column_vector_db, array.row_vector_db,
        )) * scale_nm, 12))
        for array in top.arrays
    }
    if steps != {(350.0, 0.0, 0.0, 350.0)}:
        raise RuntimeError(f"GDS array pitch is not unique: {sorted(steps)}")
    expanded = library.expanded_reference_counts("main")
    uniform_arrays = [
        {
            "cell_name": item.cell_name,
            "columns": int(item.columns),
            "rows": int(item.rows),
            "instance_count": int(item.instance_count),
            "origin_nm": (item.origin_db * scale_nm).tolist(),
        }
        for item in top.arrays
        if item.columns == 285 and item.rows == 285
    ]
    bounds_nm = _top_bounds(library, "main") * scale_nm
    layers = Counter(
        polygon.layer
        for cell in library.cells.values()
        for polygon in cell.polygons
    )
    return {
        "schema": "petch.zhu-npg80-gds-geometry-audit.v1",
        "condition_id": "zhu-2026-npg80-tio2-chf3-sf6-o2-20min",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": _hash(SOURCE),
            "size_bytes": SOURCE.stat().st_size,
            "received_filename": "GDS SEM Image 450nm.gds",
            "target_sem_used": False,
        },
        "gdsii": {
            "library_name": library.name,
            "database_unit_m": library.database_unit_m,
            "database_unit_nm": scale_nm,
            "user_unit_per_database_unit": (
                library.user_unit_per_database_unit),
            "top_cell": "main",
            "cell_count": len(library.cells),
            "geometry_complete_for_boundary_sref_aref": (
                library.geometry_complete),
            "top_bounds_nm": bounds_nm.tolist(),
            "top_extent_um": [
                float((bounds_nm[2] - bounds_nm[0]) * 1.0e-3),
                float((bounds_nm[3] - bounds_nm[1]) * 1.0e-3),
            ],
            "polygon_layer_counts": {
                str(key): int(value) for key, value in sorted(layers.items())
            },
            "top_direct_polygon_count": len(top.polygons),
            "top_sref_count": len(top.references),
            "top_aref_count": len(top.arrays),
            "expanded_leaf_instance_count": int(sum(expanded.values())),
            "expanded_leaf_counts": expanded,
            "unique_array_pitch_nm": {
                "x": 350.0,
                "y": 350.0,
            },
            "uniform_285_by_285_arrays": uniform_arrays,
        },
        "exact_mask_primitives": {
            "square_cells": square_cells,
            "special_cells": special_cells,
        },
        "simulation_board": {
            "pitch_nm": 350.0,
            "square_widths_nm": sorted(
                item["declared_width_nm"] for item in square_cells),
            "special_cells": sorted(item["cell_name"] for item in special_cells),
            "replaces_same_group_geometry_prior": True,
            "does_not_change_frozen_recipe": True,
            "does_not_identify_surface_coefficients": True,
        },
        "interpretation_limits": {
            "gds_has_no_material_semantics": True,
            "mask_polarity_confirmed_by_operator": False,
            "gds_filename_450nm_is_not_array_pitch": True,
            "internal_array_pitch_nm": 350.0,
            "library_name_indicates_700nm_zep_dose_test": True,
            "permitted": (
                "exact geometry/density board conditional on polygons being Cr islands"
            ),
            "not_permitted": (
                "claim that the GDS fixes self-bias, TiO2/Cr rates, or the final SEM"
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build())
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
        return
    if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("Zhu GDS geometry audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
