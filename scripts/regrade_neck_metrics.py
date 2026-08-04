"""Regrade archived feature checkpoints with the community-standard CD triple.

Campaign 5 established that ``mask_opening_nm`` minimises over the mask band
alone, so a pinch at the mask top and a neck part-way down the mask report the
same scalar.  Krueger's ``w_m`` and the SEM both neck 200-270 nm below the mask
top, so a run can miss on aperture SIZE, on neck LOCATION, or on both, and the
legacy scalar cannot tell them apart.  This tool replays archived
``checkpoint.npz`` interfaces through the same measurement the pilot now emits.

Usage:
    python scripts/regrade_neck_metrics.py RUN_DIR [RUN_DIR ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from petch.feature_geometry_state_3d import FeatureGeometry3D  # noqa: E402

from krueger_2024_trench_pilot import measure_krueger_metrics  # noqa: E402


def geometry_from_checkpoint(path):
    """Rebuild the dense geometry state stored beside a curated audit."""
    data = np.load(path)
    metadata = json.loads(str(data["metadata_json"]))
    levelsets = {
        int(key.rsplit("_", 1)[1]): np.asarray(data[key], dtype=float)
        for key in data.files if key.startswith("material_levelset_")}
    geometry = FeatureGeometry3D(
        phi=np.asarray(data["phi"], dtype=float),
        material_id=np.asarray(data["material_id"], dtype=int),
        dx=float(metadata["dx"]),
        mesh_length_unit_m=float(metadata["mesh_length_unit_m"]),
        mesh_origin_m=tuple(np.asarray(data["mesh_origin_m"], dtype=float)),
        material_levelsets=levelsets or None)
    return geometry, metadata


def regrade(run_dir, *, substrate_top_um=1.8):
    run_dir = Path(run_dir)
    geometry, metadata = geometry_from_checkpoint(run_dir / "checkpoint.npz")
    metrics = measure_krueger_metrics(
        geometry, substrate_top_um=float(substrate_top_um),
        aperture_profile_points=None)
    mask_top = float(metrics["mask_top_z_um"])
    profile = metrics["aperture_profile"]
    mask_band = [item for item in profile if item["z_um"] >= substrate_top_um]
    target_band = [
        item for item in profile
        if 200.0 <= (mask_top - item["z_um"]) * 1e3 <= 270.0]
    mask_neck = min(mask_band, key=lambda item: item["width_nm"]) if mask_band else None
    return {
        "run": run_dir.name,
        "step": metadata.get("step"),
        "physical_time_s": metadata.get("physical_time_s"),
        "dx_um": metadata.get("dx"),
        **{name: metrics[name] for name in (
            "top_cd_nm", "neck_cd_nm", "neck_z_um",
            "neck_depth_from_mask_top_nm", "mask_opening_nm",
            "etch_depth_nm", "remaining_mask_thickness_nm", "mask_top_z_um")},
        # Krueger's w_m lives in the mask; the global neck may sit in the etched
        # trench, so both are reported rather than conflated.
        "mask_band_neck_nm": (
            None if mask_neck is None else float(mask_neck["width_nm"])),
        "mask_band_neck_depth_nm": (
            None if mask_neck is None
            else (mask_top - float(mask_neck["z_um"])) * 1e3),
        # Aperture where the SEM (200 nm) and MCFPM (271 nm) actually neck.
        "target_band_min_nm": (
            None if not target_band
            else min(item["width_nm"] for item in target_band)),
        "target_band_max_nm": (
            None if not target_band
            else max(item["width_nm"] for item in target_band)),
        "aperture_profile": metrics["aperture_profile"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dirs", nargs="+", type=Path)
    parser.add_argument("--substrate-top-um", type=float, default=1.8)
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args(argv)

    rows = [regrade(item, substrate_top_um=args.substrate_top_um)
            for item in args.run_dirs]
    header = (f"{'run':24s} {'topCD':>7s} {'neckCD':>7s} {'neck z':>7s} "
              f"{'maskCD':>7s} {'mask z':>7s} {'200-270nm':>12s} {'depth':>7s}")
    print(header)
    print("-" * len(header))
    for row in rows:
        band = (
            "n/a" if row["target_band_min_nm"] is None
            else f"{row['target_band_min_nm']:.1f}-{row['target_band_max_nm']:.1f}")
        print(f"{row['run']:24s} {row['top_cd_nm']:7.1f} {row['neck_cd_nm']:7.1f} "
              f"{row['neck_depth_from_mask_top_nm']:7.0f} "
              f"{row['mask_band_neck_nm']:7.1f} {row['mask_band_neck_depth_nm']:7.0f} "
              f"{band:>12s} {row['etch_depth_nm']:7.1f}")
    if args.json_out is not None:
        args.json_out.write_text(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
