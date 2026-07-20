#!/usr/bin/env python3
"""Run the standard Krueger trench pilot while saving a mid-plane cross-section frame
after every accepted step, for presentation animations.

This wraps the campaign pilot's checkpoint writer without modifying it: the physics,
operator, and receipts are exactly the pilot's. Frames land in <output>/frames/ as
compressed npz (float32 phi and per-material level-set slices at the y mid-plane).
All pilot CLI arguments apply unchanged.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import krueger_2024_trench_pilot as pilot  # noqa: E402

_FRAMES_DIR = None
_ORIGINAL_CHECKPOINT = pilot._checkpoint


def _frame_checkpoint(path, geometry, state, fingerprint, step, physical_time_s,
                      next_step_duration_s):
    j = geometry.phi.shape[1] // 2
    arrays = {
        "phi": np.asarray(geometry.phi[:, j, :], dtype=np.float32),
        "scalars": np.asarray(
            [float(step), float(physical_time_s), float(geometry.dx),
             float(geometry.mesh_length_unit_m)], dtype=np.float64),
    }
    for material_id, field in dict(geometry.material_levelsets or {}).items():
        arrays[f"material_{int(material_id)}"] = np.asarray(
            field[:, j, :], dtype=np.float32)
    out = _FRAMES_DIR / f"frame_{int(step):05d}.npz"
    with out.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    return _ORIGINAL_CHECKPOINT(
        path, geometry, state, fingerprint, step, physical_time_s,
        next_step_duration_s)


def main():
    global _FRAMES_DIR
    args = pilot.parse_args()
    _FRAMES_DIR = Path(args.output) / "frames"
    _FRAMES_DIR.mkdir(parents=True, exist_ok=True)
    pilot._checkpoint = _frame_checkpoint
    pilot.run(args)


if __name__ == "__main__":
    main()
