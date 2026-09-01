#!/usr/bin/env python3
"""Pillar height from tilted SEM geometry (Freddie's 0817 set).

Stage tilt theta about the image X axis (Zeiss "Tilt = In Y"): a vertical
sidewall of height h projects to h*sin(theta) along image Y; the pillar top of
width a projects to a*cos(theta).  Along X nothing foreshortens.  So per
standing pillar:  Y_extent = a*cos(theta) + h*sin(theta),  X_extent = a
  =>  h = (Y_extent - X_extent*cos(theta)) / sin(theta).
Lattice-averaged unit cells (half-height threshold) supply X/Y extents; the
same is done for stub sites in loss zones (remaining height).  Also writes an
annotated frame with the two lines drawn on real pillars.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
SEM = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80" / "sem_0817"
FOOTER = 683

import importlib.util
spec = importlib.util.spec_from_file_location("lat", ROOT / "scripts" / "measure_zhu_npg80_sem_lattice.py")
lat = importlib.util.module_from_spec(spec); spec.loader.exec_module(lat)


def cells_and_extents(img, px_x, px_y, pixel_nm, theta_deg):
    xs, ys = lat.lattice_sites(img, px_x, px_y)
    m = int(min(px_x, px_y) * 0.55)
    xs = xs[(xs > m) & (xs < img.shape[1] - m)]
    ys = ys[(ys > m) & (ys < img.shape[0] - m)]
    half = int(min(px_x, px_y) // 2) - 1
    floor = np.percentile(img, 20); top = np.percentile(img, 97)
    rad = int(min(px_x, px_y) * 0.28)
    yy, xx = np.mgrid[-rad:rad + 1, -rad:rad + 1]; disc = (xx**2 + yy**2) <= rad**2
    rows = []
    for y in ys:
        for x in xs:
            yi, xi = int(round(y)), int(round(x))
            patch = img[yi - rad:yi + rad + 1, xi - rad:xi + rad + 1]
            cell = img[yi - half:yi + half + 1, xi - half:xi + half + 1]
            if patch.shape != disc.shape or cell.shape != (2*half+1, 2*half+1):
                continue
            bright = patch[disc].mean() > floor + 0.45 * (top - floor)
            lvl = floor + 0.5 * (np.percentile(cell, 95) - floor)
            pm = cell > lvl
            lab, n = ndimage.label(pm)
            if not n:
                continue
            c = lab[half, half]
            if c == 0:
                sizes = ndimage.sum(pm, lab, range(1, n + 1)); c = int(np.argmax(sizes)) + 1
            pm = lab == c
            if pm.sum() < 30:
                continue
            xe = pm.any(0).sum() * pixel_nm; ye = pm.any(1).sum() * pixel_nm
            th = np.radians(theta_deg)
            h = (ye - xe * np.cos(th)) / np.sin(th)
            rows.append({"x": xi, "y": yi, "standing": bool(bright), "x_extent_nm": xe,
                         "y_extent_nm": ye, "height_nm": h,
                         "row0": yi - half + int(np.where(pm.any(1))[0].min()),
                         "row1": yi - half + int(np.where(pm.any(1))[0].max()),
                         "col0": xi - half + int(np.where(pm.any(0))[0].min()),
                         "col1": xi - half + int(np.where(pm.any(0))[0].max())})
    return rows


def main():
    manifest = json.loads((SEM / "manifest.json").read_text())
    results = []
    for fr in manifest["frames"]:
        px = fr["meta"]["pixel_size_nm"]; th = fr["meta"]["stage_tilt_deg"] or 0
        if not px or px > 7.5 or th < 5:
            continue
        img = np.asarray(Image.open(SEM / "originals" / fr["file"]).convert("L"), float)[:FOOTER]
        px_x, px_y, _, _ = lat.fft_pitch(img)
        if not px_x or not px_y or px_x < 14 or px_y < 14:
            continue
        rows = cells_and_extents(img, px_x, px_y, px, th)
        stand = [r for r in rows if r["standing"]]
        stubs = [r for r in rows if not r["standing"] and r["x_extent_nm"] > 20]
        if len(stand) >= 5:
            hs = np.array([r["height_nm"] for r in stand])
            results.append({"file": fr["file"], "dose": fr["dose_index"], "tilt_deg": th,
                            "pixel_nm": px, "n_standing": len(stand),
                            "height_nm_median": float(np.median(hs)),
                            "height_nm_p25": float(np.percentile(hs, 25)),
                            "height_nm_p75": float(np.percentile(hs, 75)),
                            "x_extent_nm_median": float(np.median([r["x_extent_nm"] for r in stand])),
                            "y_extent_nm_median": float(np.median([r["y_extent_nm"] for r in stand])),
                            "n_stubs": len(stubs),
                            "stub_height_nm_median": float(np.median([r["height_nm"] for r in stubs])) if len(stubs) >= 3 else None})
    (SEM / "height_measurements.json").write_text(json.dumps(results, indent=1))
    for r in results:
        stub = f" stubs n={r['n_stubs']} h~{r['stub_height_nm_median']:.0f}" if r["stub_height_nm_median"] else ""
        print(f"{r['file']:<18} T{r['tilt_deg']:>4.0f} {r['pixel_nm']:.2f}nm/px n={r['n_standing']:>3} "
              f"X {r['x_extent_nm_median']:.0f} Y {r['y_extent_nm_median']:.0f} -> h {r['height_nm_median']:.0f} nm "
              f"[{r['height_nm_p25']:.0f}-{r['height_nm_p75']:.0f}]{stub}")
    # weighted summary: 15-degree frames carry 3.9x less error amplification than 10-degree
    t15 = [r for r in results if r["tilt_deg"] >= 12]; t10 = [r for r in results if r["tilt_deg"] < 12]
    def summ(rs):
        hs = np.array([r["height_nm_median"] for r in rs]); return (float(np.median(hs)), float(np.percentile(hs, 25)), float(np.percentile(hs, 75)), len(rs))
    if t15: print("15-deg frames: h median %.0f nm [%.0f-%.0f], n=%d" % summ(t15))
    if t10: print("10-deg frames: h median %.0f nm [%.0f-%.0f], n=%d" % summ(t10))

    # annotated frame: draw the two lines on real pillars
    fname = "dose9_15d_02.tif"
    fr = [f for f in manifest["frames"] if f["file"] == fname][0]
    px = fr["meta"]["pixel_size_nm"]; th = fr["meta"]["stage_tilt_deg"]
    img = np.asarray(Image.open(SEM / "originals" / fname).convert("L"), float)[:FOOTER]
    px_x, px_y, _, _ = lat.fft_pitch(img)
    rows = cells_and_extents(img, px_x, px_y, px, th)
    canvas = Image.open(SEM / "originals" / fname).convert("RGB")
    d = ImageDraw.Draw(canvas)
    for r in rows:
        col = (45, 212, 191) if r["standing"] else (245, 176, 74)
        d.rectangle([r["col0"], r["row0"], r["col1"], r["row1"]], outline=col, width=1)
        d.line([r["col1"] + 3, r["row0"], r["col1"] + 3, r["row1"]], fill=(226, 106, 208), width=2)  # Y extent
        d.line([r["col0"], r["row1"] + 3, r["col1"], r["row1"] + 3], fill=(110, 168, 254), width=2)  # X extent
    stand = [r for r in rows if r["standing"]]
    hs = np.median([r["height_nm"] for r in stand]) if stand else float("nan")
    d.rectangle([8, 8, 560, 60], fill=(5, 10, 20))
    d.text((14, 12), f"{fname}  tilt {th:.0f} deg  {px:.2f} nm/px", fill=(232, 236, 244))
    d.text((14, 28), "pink = Y extent (a cos t + h sin t), blue = X extent (a)", fill=(232, 236, 244))
    d.text((14, 44), f"standing pillars: median h = {hs:.0f} nm", fill=(45, 212, 191))
    canvas.save(SEM / "height_annotated_dose9_15d_02.png")
    print("annotated ->", SEM / "height_annotated_dose9_15d_02.png")


if __name__ == "__main__":
    main()
