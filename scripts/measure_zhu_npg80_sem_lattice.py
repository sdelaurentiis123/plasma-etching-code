#!/usr/bin/env python3
"""Lattice-based measurement of Freddie's 0817 SEM dose series.

Per plan-view frame (10 deg stage tilt, InLens):
  * pitch from the 2-D power spectrum (robust to lobed pillar shapes);
  * site occupancy: every lattice site is classified standing / lost from the
    mean brightness in a disc around the site -> survival fraction;
  * CD from the lattice-averaged unit cell: the mean pillar image is
    thresholded at half of (top - floor) and the equivalent-square width and
    the axis-projected width are reported in nm.
Tilted (15 deg) frames are catalogued only.  Writes lattice_measurements.json
and a per-dose summary; no model quantity enters.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
SEM = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80" / "sem_0817"
FOOTER = 683


def fft_pitch(img):
    f = np.abs(np.fft.fftshift(np.fft.fft2(img - img.mean())))
    h, w = f.shape
    cy, cx = h // 2, w // 2
    f[cy - 3:cy + 4, cx - 3:cx + 4] = 0
    # search radial ring for the strongest peak away from DC, restricted to
    # spatial periods between 12 and 400 px
    yy, xx = np.mgrid[:h, :w]
    r = np.hypot(yy - cy, xx - cx)
    period = np.where(r > 0, np.sqrt((w / np.maximum(xx - cx, 1e-9)) ** 2 * 0 + 1), 1)
    # use axis-aligned peaks: lattice is close to axis aligned in these frames
    prof_x = f[cy - 2:cy + 3, :].sum(0)
    prof_y = f[:, cx - 2:cx + 3].sum(1)
    def best(prof, n):
        k = np.arange(len(prof)) - len(prof) // 2
        mask = (np.abs(k) >= n / 400) & (np.abs(k) <= n / 12)
        idx = np.argmax(np.where(mask, prof, 0))
        return abs(k[idx])
    kx, ky = best(prof_x, w), best(prof_y, h)
    px_x = w / kx if kx else None
    px_y = h / ky if ky else None
    return px_x, px_y, kx, ky


def lattice_sites(img, px_x, px_y):
    """Estimate lattice phase by maximizing summed brightness on the grid."""
    h, w = img.shape
    sm = ndimage.gaussian_filter(img, 2)
    best = None
    for ox in np.linspace(0, px_x, 12, endpoint=False):
        for oy in np.linspace(0, px_y, 12, endpoint=False):
            xs = np.arange(ox, w, px_x); ys = np.arange(oy, h, px_y)
            X, Y = np.meshgrid(xs, ys)
            v = ndimage.map_coordinates(sm, [Y.ravel(), X.ravel()], order=1, mode="nearest")
            s = v.mean()
            if best is None or s > best[0]:
                best = (s, ox, oy)
    _, ox, oy = best
    xs = np.arange(ox, w, px_x); ys = np.arange(oy, h, px_y)
    return xs, ys


def measure(path: Path, pixel_nm: float):
    img = np.asarray(Image.open(path).convert("L"), float)[:FOOTER]
    px_x, px_y, kx, ky = fft_pitch(img)
    if not px_x or not px_y or px_x < 14 or px_y < 14:
        return {"error": "no lattice"}
    xs, ys = lattice_sites(img, px_x, px_y)
    # keep interior sites
    m = int(min(px_x, px_y) * 0.55)
    xs = xs[(xs > m) & (xs < img.shape[1] - m)]
    ys = ys[(ys > m) & (ys < img.shape[0] - m)]
    if len(xs) < 2 or len(ys) < 2:
        return {"error": "too few sites"}
    rad = int(min(px_x, px_y) * 0.28)
    yy, xx = np.mgrid[-rad:rad + 1, -rad:rad + 1]
    disc = (xx ** 2 + yy ** 2) <= rad ** 2
    site_mean = []
    cells = []
    half = int(min(px_x, px_y) // 2) - 1
    for y in ys:
        for x in xs:
            yi, xi = int(round(y)), int(round(x))
            patch = img[yi - rad:yi + rad + 1, xi - rad:xi + rad + 1]
            if patch.shape != disc.shape:
                continue
            site_mean.append(patch[disc].mean())
            cell = img[yi - half:yi + half + 1, xi - half:xi + half + 1]
            if cell.shape == (2 * half + 1, 2 * half + 1):
                cells.append(cell)
    site_mean = np.array(site_mean)
    floor = np.percentile(img, 20)
    top = np.percentile(img, 97)
    standing = site_mean > floor + 0.45 * (top - floor)
    # survival only over the array's own extent: lattice rows/cols that hold
    # at least one standing pillar bound the block; empty field outside the
    # written block must not count as loss.  Islands inside do count.
    grid = standing.reshape(len(ys), len(xs)) if standing.size == len(ys) * len(xs) else None
    if grid is not None and grid.any():
        rows = np.where(grid.any(1))[0]; cols = np.where(grid.any(0))[0]
        inner = grid[rows.min():rows.max() + 1, cols.min():cols.max() + 1]
        survival = float(inner.mean())
        block_sites = int(inner.size)
    else:
        survival = float(standing.mean()); block_sites = int(standing.size)
    cd = None
    if standing.sum() >= 3:
        stand_cells = np.array([c for c, s in zip(cells, standing) if s])
        mean_cell = stand_cells.mean(0)
        lvl = floor + 0.5 * (np.percentile(mean_cell, 95) - floor)
        pm = mean_cell > lvl
        # keep the central component only
        lab, n = ndimage.label(pm)
        if n:
            c = lab[half, half]
            if c == 0:
                sizes = ndimage.sum(pm, lab, range(1, n + 1))
                c = int(np.argmax(sizes)) + 1
            pm = lab == c
            area = pm.sum()
            cols = pm.any(0).sum(); rows = pm.any(1).sum()
            cd = {
                "equivalent_square_nm": float(np.sqrt(area) * pixel_nm),
                "x_extent_nm": float(cols * pixel_nm),
                "y_extent_nm": float(rows * pixel_nm / np.cos(np.radians(10.0))),
            }
    return {
        "pitch_x_nm": float(px_x * pixel_nm),
        "pitch_y_nm": float(px_y * pixel_nm / np.cos(np.radians(10.0))),
        "sites": int(len(site_mean)),
        "block_sites": block_sites,
        "survival_fraction": survival,
        "cd": cd,
    }


def main():
    manifest = json.loads((SEM / "manifest.json").read_text())
    out = []
    for fr in manifest["frames"]:
        meta = fr["meta"]
        px = meta["pixel_size_nm"]
        tilted = fr["series_tilt_label_deg"] > 0 or (meta["stage_tilt_deg"] or 0) > 12
        rec = {"file": fr["file"], "dose": fr["dose_index"], "pixel_nm": px,
               "tilt_deg": meta["stage_tilt_deg"], "plan_view": not tilted}
        if not tilted and px and 1.5 <= px <= 30:
            rec.update(measure(SEM / "originals" / fr["file"], px))
        out.append(rec)
    (SEM / "lattice_measurements.json").write_text(json.dumps(out, indent=1))
    # per-dose summary from frames with a valid lattice and >= 20 sites
    by = {}
    for r in out:
        if r.get("cd") is None or r.get("sites", 0) < 20:
            continue
        by.setdefault(r["dose"], []).append(r)
    summary = {}
    for d, rs in sorted(by.items()):
        cds = [r["cd"]["equivalent_square_nm"] for r in rs]
        xs = [r["cd"]["x_extent_nm"] for r in rs]
        pit = [r["pitch_x_nm"] for r in rs]
        surv = [r["survival_fraction"] for r in rs]
        summary[d] = {
            "frames": len(rs),
            "cd_equivalent_square_nm": [float(np.median(cds)), float(np.min(cds)), float(np.max(cds))],
            "cd_x_extent_nm": [float(np.median(xs)), float(np.min(xs)), float(np.max(xs))],
            "pitch_nm_median": float(np.median(pit)),
            "survival_fraction": [float(np.median(surv)), float(np.min(surv)), float(np.max(surv))],
        }
    (SEM / "dose_summary.json").write_text(json.dumps(summary, indent=1))
    for d, s in summary.items():
        print(f"dose{d}: n={s['frames']} CD {s['cd_equivalent_square_nm'][0]:.0f} nm "
              f"[{s['cd_equivalent_square_nm'][1]:.0f}-{s['cd_equivalent_square_nm'][2]:.0f}] "
              f"x-extent {s['cd_x_extent_nm'][0]:.0f} pitch {s['pitch_nm_median']:.0f} "
              f"survival {s['survival_fraction'][0]:.2f} [{s['survival_fraction'][1]:.2f}-{s['survival_fraction'][2]:.2f}]")


if __name__ == "__main__":
    main()
