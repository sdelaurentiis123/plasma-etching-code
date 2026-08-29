#!/usr/bin/env python3
"""Certify the four non-square primitive cells in Freddie's supplied GDS."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.feature_step_3d import make_footprint_mask_geometry_3d
from petch.mask_footprints import (
    centered_cross_footprint_levelset,
    centered_inverse_square_hole_footprint_levelset,
    centered_rectangle_footprint_levelset,
)


DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
GDS_AUDIT = DATA / "layout" / "gds_geometry_audit.json"
PREREGISTRATION = DATA / "gds_square_profile_preregistration.json"
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_gds_special_geometry_v1"
    / "audit.json"
)
PITCH_NM = 350.0
DX_NM = 10.0


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _footprint(name):
    common = {"pitch": PITCH_NM * 1e-3, "dx": DX_NM * 1e-3}
    if name == "RECT_250x105":
        return centered_rectangle_footprint_levelset(
            cell_width=common["pitch"],
            cell_length=common["pitch"],
            dx=common["dx"],
            rectangle_width=0.250,
            rectangle_length=0.105,
        ), 250.0 * 105.0 / PITCH_NM ** 2
    if name == "CROSS_250x105":
        return centered_cross_footprint_levelset(
            **common, outer_width=0.250, arm_width=0.105,
        ), (2.0 * 250.0 * 105.0 - 105.0 ** 2) / PITCH_NM ** 2
    if name == "INVHOLE_105":
        return centered_inverse_square_hole_footprint_levelset(
            **common, opening_width=0.105,
        ), 1.0 - 105.0 ** 2 / PITCH_NM ** 2
    if name == "INVHOLE_250":
        return centered_inverse_square_hole_footprint_levelset(
            **common, opening_width=0.250,
        ), 1.0 - 250.0 ** 2 / PITCH_NM ** 2
    raise ValueError(f"unsupported special GDS primitive: {name}")


def build():
    gds = _load(GDS_AUDIT)
    preregistration = _load(PREREGISTRATION)
    names = gds["simulation_board"]["special_cells"]
    rows = []
    for name in names:
        footprint, exact_fraction = _footprint(name)
        geometry = make_footprint_mask_geometry_3d(
            cell_width=PITCH_NM * 1e-3,
            cell_length=PITCH_NM * 1e-3,
            domain_height=1.05,
            dx=DX_NM * 1e-3,
            mask_footprint_levelset_xy=footprint,
            film_thickness=0.700,
            mask_thickness=0.045,
            base_top=0.100,
            film_material_id=1,
            mask_material_id=2,
            base_material_id=3,
        )
        center = tuple(size // 2 for size in footprint.shape)
        rows.append({
            "cell_name": name,
            "analytic_mask_area_fraction": float(exact_fraction),
            "nodal_positive_mask_fraction": float(np.mean(footprint >= 0.0)),
            "center_is_mask": bool(footprint[center] >= 0.0),
            "footprint_periodic_x": bool(np.allclose(
                footprint[0, :], footprint[-1, :])),
            "footprint_periodic_y": bool(np.allclose(
                footprint[:, 0], footprint[:, -1])),
            "geometry_shape": list(geometry.phi.shape),
            "material_ids": sorted(map(int, geometry.material_levelsets)),
            "finite_levelsets": bool(
                np.all(np.isfinite(geometry.phi))
                and all(np.all(np.isfinite(field))
                        for field in geometry.material_levelsets.values())),
        })
    return {
        "schema": "petch.zhu-npg80-gds-special-geometry-audit.v1",
        "condition_id": preregistration["condition_id"],
        "gds_sha256": gds["source"]["sha256"],
        "gds_audit_sha256": _hash(GDS_AUDIT),
        "mesh_spacing_nm": DX_NM,
        "pitch_nm": PITCH_NM,
        "mask_polarity_assumption": preregistration[
            "exact_layout_geometry"]["mask_polarity_assumption"],
        "mask_polarity_confirmed_by_operator": False,
        "target_sem_used": False,
        "profile_evolution_geometry_ready": True,
        "profile_evolution_surface_response_still_conditional": True,
        "rows": rows,
    }


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build())
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
        return
    if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("exact-GDS special geometry audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
