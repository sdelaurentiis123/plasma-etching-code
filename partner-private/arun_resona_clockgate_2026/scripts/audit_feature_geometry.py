#!/usr/bin/env python3
"""Prove the recovered clock-gate polygons enter the layered 3-D engine."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MAIN_REPO = ROOT.parents[1]
if not (MAIN_REPO / "src" / "petch").is_dir():
    MAIN_REPO = MAIN_REPO / "plasma-etching-code"
sys.path.insert(0, str(MAIN_REPO / "src"))

from petch.feature_step_3d import make_footprint_mask_geometry_3d  # noqa: E402
from petch.gds_import import read_gds  # noqa: E402
from petch.mask_footprints import polygon_union_footprint_levelset  # noqa: E402


OUTPUT = ROOT / "results" / "feature_geometry_audit.json"
CELL = "CLKGATE_X1"
LAYER = 11
DATATYPE = 0
ENLARGEMENT = 40.0
PRODUCTION_DX_UM = 0.2
PILOT_DX_UM = 0.4


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _area(vertices):
    following = np.roll(vertices, -1, axis=0)
    return 0.5 * abs(float(np.sum(
        vertices[:, 0] * following[:, 1]
        - following[:, 0] * vertices[:, 1]
    )))


def _source_polygons_um(reference_gds):
    library = read_gds(reference_gds)
    polygons = [
        polygon.vertices_db * library.database_unit_m * 1.0e6
        for polygon in library.cells[CELL].polygons
        if polygon.layer == LAYER and polygon.datatype == DATATYPE
    ]
    if len(polygons) != 10:
        raise ValueError("reference GDS no longer has ten Metal-1 polygons")
    combined = np.concatenate(polygons)
    lower = np.min(combined, axis=0)
    physical = [(vertices - lower) * ENLARGEMENT for vertices in polygons]
    extent = (np.max(combined, axis=0) - lower) * ENLARGEMENT
    return physical, extent


def build(reference_gds):
    polygons, extent = _source_polygons_um(reference_gds)
    width, length = map(float, extent)
    exact_area = sum(_area(vertices) for vertices in polygons)
    exact_fraction = exact_area / (width * length)

    production_footprint = polygon_union_footprint_levelset(
        cell_width=width,
        cell_length=length,
        dx=PRODUCTION_DX_UM,
        polygons=polygons,
    )
    pilot_footprint = polygon_union_footprint_levelset(
        cell_width=width,
        cell_length=length,
        dx=PILOT_DX_UM,
        polygons=polygons,
    )
    geometry = make_footprint_mask_geometry_3d(
        cell_width=width,
        cell_length=length,
        domain_height=38.4,
        dx=PILOT_DX_UM,
        mask_footprint_levelset_xy=pilot_footprint,
        film_thickness=6.0,
        mask_thickness=30.0,
        base_top=1.2,
        film_material_id=1,
        mask_material_id=2,
        base_material_id=3,
    )
    production_shape = [
        production_footprint.shape[0],
        production_footprint.shape[1],
        int(round(38.4 / PRODUCTION_DX_UM)) + 1,
    ]
    production_voxels = int(np.prod(production_shape))

    return {
        "schema": "petch.partner.clockgate-feature-geometry-audit.v1",
        "reference_gds": {
            "sha256": _hash(reference_gds),
            "cell": CELL,
            "layer": LAYER,
            "datatype": DATATYPE,
        },
        "physical_geometry": {
            "footprint_extent_um": [width, length],
            "mask_height_um": 30.0,
            "polygon_count": len(polygons),
            "exact_mask_area_um2": exact_area,
            "exact_mask_area_fraction": exact_fraction,
            "working_polarity": "positive footprint protects silicon",
        },
        "production_footprint": {
            "dx_um": PRODUCTION_DX_UM,
            "shape_with_periodic_endpoints": list(production_footprint.shape),
            "minimum_opening_cells": 2.6 / PRODUCTION_DX_UM,
            "periodic_x": bool(np.array_equal(
                production_footprint[0, :], production_footprint[-1, :])),
            "periodic_y": bool(np.array_equal(
                production_footprint[:, 0], production_footprint[:, -1])),
            "finite": bool(np.all(np.isfinite(production_footprint))),
            "has_mask_and_opening": bool(
                np.any(production_footprint > 0.0)
                and np.any(production_footprint < 0.0)),
            "full_3d_shape_at_38p4_um": production_shape,
            "full_3d_voxel_count": production_voxels,
            "float64_bytes_per_scalar_field": production_voxels * 8,
        },
        "pilot_layered_geometry": {
            "dx_um": PILOT_DX_UM,
            "shape": list(geometry.phi.shape),
            "voxel_count": int(geometry.phi.size),
            "material_ids": sorted(map(int, geometry.material_levelsets)),
            "finite": bool(
                np.all(np.isfinite(geometry.phi))
                and all(np.all(np.isfinite(field))
                        for field in geometry.material_levelsets.values())),
            "contains_gas": bool(np.any(geometry.material_id == 0)),
            "contains_film_mask_and_base": bool(all(
                np.any(geometry.material_id == material_id)
                for material_id in (1, 2, 3))),
        },
        "readiness": {
            "exact_layout_to_layered_3d_engine_ready": True,
            "target_sem_used": False,
            "absolute_surface_response_certified_for_target_tool": False,
            "next_step": "run the evidence-bounded SF6/O2 boundary/surface ensemble on the exact layered geometry",
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
    rendered = json.dumps(build(args.reference_gds), indent=2, sort_keys=True) + "\n"
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    elif not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("feature-geometry audit is stale")
    else:
        print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
