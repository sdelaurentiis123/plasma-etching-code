#!/usr/bin/env python3
"""Render presentation PNG frames from krueger_2024_presentation_frames npz dumps.

Usage: render_presentation_frames.py <run_dir> [--out <dir>] [--max-frames N] [--scale S]
Reads <run_dir>/frames/frame_*.npz and <run_dir>/audit.json (for per-step metrics),
writes styled cross-section PNGs to <out> (default <run_dir>/rendered).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GAS = "#0b1220"
OXIDE = "#1f8a8c"
MASK = "#b5348f"
EDGE = "#f5f0e8"


def render_frame(npz_path, metrics, out_path, dpi=160, tiles=3, z_crop_nm=(400.0, 2750.0)):
    data = np.load(npz_path)
    step, t_s, dx_um, _unit = data["scalars"]
    phi = data["phi"]           # (nx, nz), >0 = solid
    nx, nz = phi.shape

    mask_ids = [k for k in data.files if k.startswith("material_")]
    solids = {k: np.asarray(data[k]) for k in mask_ids}
    # The mask material is the one whose solid support reaches the domain top.
    def top_extent(field):
        rows = np.any(field > 0.0, axis=0)
        return np.max(np.nonzero(rows)[0]) if rows.any() else -1
    mask_key = max(solids, key=lambda k: top_extent(solids[k]))

    # Tile the periodic cell laterally so the frame reads as a trench array.
    def tiled(field):
        return np.concatenate([field[:-1]] * tiles + [field[-1:]], axis=0)
    phi_t = tiled(phi)
    solids_t = {k: tiled(v) for k, v in solids.items()}
    x_nm = np.arange(phi_t.shape[0]) * dx_um * 1e3
    z_nm = np.arange(nz) * dx_um * 1e3

    fig, ax = plt.subplots(figsize=(3.6, 7.4), dpi=dpi)
    fig.patch.set_facecolor(GAS)
    ax.set_facecolor(GAS)
    for key in (k for k in solids_t if k != mask_key):
        ax.contourf(x_nm, z_nm, solids_t[key].T, levels=[0.0, np.inf], colors=[OXIDE])
    ax.contourf(x_nm, z_nm, solids_t[mask_key].T, levels=[0.0, np.inf], colors=[MASK])
    ax.contour(x_nm, z_nm, phi_t.T, levels=[0.0], colors=[EDGE], linewidths=0.9)

    ax.set_xlim(x_nm[0], x_nm[-1])
    ax.set_ylim(*z_crop_nm)
    ax.set_aspect("equal")
    ax.tick_params(colors="#8fa3bf", labelsize=7)
    for spine in ax.spines.values():
        spine.set_color("#31405c")
    ax.set_xlabel("x (nm)", color="#8fa3bf", fontsize=8)
    ax.set_ylabel("z (nm)", color="#8fa3bf", fontsize=8)

    label = f"t = {t_s:5.1f} s"
    if metrics is not None:
        label += (f"   depth {metrics['etch_depth_nm']:4.0f} nm"
                  f"   opening {metrics['mask_opening_nm']:3.0f} nm")
    fig.suptitle(label, color=EDGE, fontsize=9.5, family="monospace", y=0.995)
    fig.tight_layout(pad=0.4, rect=(0, 0, 1, 0.975))
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--max-frames", type=int, default=160)
    args = parser.parse_args()

    frames = sorted((args.run_dir / "frames").glob("frame_*.npz"))
    if not frames:
        raise SystemExit("no frames found")
    out = args.out or (args.run_dir / "rendered")
    out.mkdir(parents=True, exist_ok=True)

    metrics_by_step = {}
    audit_path = args.run_dir / "audit.json"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text())
        for item in audit.get("history", ()):
            metrics_by_step[int(item["step"])] = item["metrics"]

    if len(frames) > args.max_frames:
        keep = np.unique(np.linspace(0, len(frames) - 1, args.max_frames).astype(int))
        frames = [frames[i] for i in keep]
    for index, path in enumerate(frames):
        step = int(path.stem.split("_")[1])
        render_frame(path, metrics_by_step.get(step), out / f"r_{index:04d}.png")
    print(f"rendered {len(frames)} frames -> {out}")


if __name__ == "__main__":
    main()
