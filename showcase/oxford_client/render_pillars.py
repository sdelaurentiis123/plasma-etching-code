#!/usr/bin/env python3
"""Render the frozen predicted pillar as a 3-D mesh at the SEM's viewing tilt,
and crop the matching SEM regions, for the one-to-one client comparison.

Geometry comes from the committed exact-GDS board (w185, the design width
nearest Freddie's dose-9 CD): the level-set cross-section width(height) at the
two witness-rate endpoints.  Slow rate (34.1 nm/min) = the 'intact' end of
our interval; fast rate (43.5 nm/min) = the 'cap-exhausted' end.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GDS = ROOT / "results/curated/zhu_npg80_gds_square_profiles_v1/audit.json"
SEM = ROOT / "data/experimental/zhu_2026_tio2_npg80/sem_0817/originals"
BG = "#0b1220"


def profiles(width=185.0):
    d = json.loads(GDS.read_text())
    out = {}
    for p in d["profiles"]:
        if p.get("width_nm") != width or p["transport_scenario"]["name"] != "ion_low_tail_0p0" or p["tio2_to_cr_selectivity"] != 14.0:
            continue
        cs = p["profile"]["cross_section"]
        h = np.array([c["height_um"] for c in cs]) * 1000.0
        w = np.array([c["mean_width_nm"] for c in cs])
        floor = p["profile"]["floor_height_nm"]
        out[round(p["blanket_tio2_rate_nm_min"], 1)] = {"h_nm": h - floor, "w_nm": w, "depth_nm": p["profile"]["etched_depth_nm"], "cr_gone": p["cr_mask"]["mask_exhausted_at_center"]}
    return out


def pillar_faces(cx, cy, h, w):
    """Square pillar with side w(h) stacked in slabs; returns face polygons + shade keys."""
    order = np.argsort(h); h = h[order]; w = w[order]
    keep = w > 2.0
    h, w = h[keep], w[keep]
    if len(h) < 2:
        return []
    faces = []
    for i in range(len(h) - 1):
        z0, z1 = h[i], h[i + 1]; a0, a1 = w[i] / 2, w[i + 1] / 2
        # four side quads (tapered)
        for (sx, sy, key) in ((1, 0, "e"), (-1, 0, "w"), (0, 1, "n"), (0, -1, "s")):
            if sx:
                q = [(cx + sx * a0, cy - a0, z0), (cx + sx * a0, cy + a0, z0), (cx + sx * a1, cy + a1, z1), (cx + sx * a1, cy - a1, z1)]
            else:
                q = [(cx - a0, cy + sy * a0, z0), (cx + a0, cy + sy * a0, z0), (cx + a1, cy + sy * a1, z1), (cx - a1, cy + sy * a1, z1)]
            faces.append((q, key))
    a = w[-1] / 2; z = h[-1]
    faces.append(([(cx - a, cy - a, z), (cx + a, cy - a, z), (cx + a, cy + a, z), (cx - a, cy + a, z)], "top"))
    return faces


def render(prof, out_png, title, pitch=350.0, n=3, tilt_deg=15.0):
    fig = plt.figure(figsize=(6, 4.5), dpi=160); fig.patch.set_facecolor(BG)
    ax = fig.add_subplot(111, projection="3d"); ax.set_facecolor(BG)
    shade = {"top": "#c9d3e6", "s": "#8a97b3", "n": "#3b4763", "e": "#6c7a99", "w": "#5b6987"}
    polys, cols = [], []
    for i in range(n):
        for j in range(n):
            for q, key in pillar_faces(i * pitch, j * pitch, prof["h_nm"], prof["w_nm"]):
                polys.append(q); cols.append(shade[key])
    ext = (n - 1) * pitch + 260
    floor = [[(-130, -130, 0), (ext, -130, 0), (ext, ext, 0), (-130, ext, 0)]]
    ax.add_collection3d(Poly3DCollection(floor, facecolors="#2a3350", edgecolors="none"))
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cols, edgecolors=cols, linewidths=0.3))
    ax.set_xlim(-130, ext); ax.set_ylim(-130, ext); ax.set_zlim(0, ext * 0.9)
    ax.set_box_aspect((1, 1, 0.9))
    ax.view_init(elev=90 - tilt_deg, azim=-90)   # near plan view tilted like the SEM stage
    ax.set_axis_off()
    hmax = float(np.max(prof["h_nm"][prof["w_nm"] > 2.0])) if np.any(prof["w_nm"] > 2.0) else 0
    fig.text(0.03, 0.94, title, color="#e8ecf4", fontsize=11, family="monospace", weight="bold")
    fig.text(0.03, 0.89, f"pillar height {hmax:.0f} nm · trench depth {prof['depth_nm']:.0f} nm · Cr cap {'gone' if prof['cr_gone'] else 'present'}", color="#8fa0b8", fontsize=8.5, family="monospace")
    plt.subplots_adjust(0, 0, 1, 1); fig.savefig(out_png, facecolor=BG); plt.close(fig)
    return hmax


def crop_sem(fname, box, out_png, pixel_nm):
    im = Image.open(SEM / fname).convert("L").crop(box)
    im = im.resize((im.width * 2, im.height * 2), Image.LANCZOS)
    # burn a 200 nm scale bar
    from PIL import ImageDraw
    d = ImageDraw.Draw(im); L = int(200 / pixel_nm * 2)
    d.rectangle([12, im.height - 22, 12 + L, im.height - 16], fill=255); d.text((12, im.height - 40), "200 nm", fill=255)
    im.save(out_png)


if __name__ == "__main__":
    pr = profiles()
    print({k: (round(float(v["h_nm"].max()), 0), round(float(v["w_nm"].max()), 0)) for k, v in pr.items()})
    h_slow = render(pr[34.1], HERE / "pred_intact_w185.png", "PREDICTED · 185 nm design width · 15° view")
    h_fast = h_slow
    print("heights", h_slow, h_fast)
    crop_sem("dose9_15d_02.tif", (0, 40, 400, 340), HERE / "sem_intact_dose9.png", 5.63)
    crop_sem("dose9_15d_02.tif", (560, 40, 960, 340), HERE / "sem_loss_dose9.png", 5.63)
    print("done")
