#!/usr/bin/env python3
"""Pillar height from floor occlusion at two stage tilts.

Along the tilt direction (image Y) each standing pillar hides the floor behind
it over h*sin(theta).  With lattice pitch p and top extent a along Y:
    floor_gap(theta) = (p - a) * cos(theta) - h * sin(theta)
The period-averaged Y profile of an intact region gives floor_gap as the
fraction of the period at floor level.  Frames at 10 deg and 15 deg on the
same dose give two equations for (a, h); alternatively a is taken from the
untilted X extent (tip-to-tip) and h follows from one tilt.
"""
from __future__ import annotations

import json, sys
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SEM = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80" / "sem_0817"
FOOTER = 683

# intact-region column windows (x0,x1) chosen from the whole-set review so the
# profile crosses standing pillars only; (None) = full width
WINDOWS = {
    "dose9_15d_02.tif": (10, 330), "dose9_15d_07.tif": (10, 260), "dose9_15d_08.tif": (10, 260),
    "dose7_15d_02.tif": (600, 1020),
    "dose9_01.tif": None, "dose9_02.tif": None, "dose9_04.tif": None, "dose9_08.tif": None,
    "dose7_06.tif": None, "dose7_07.tif": None, "dose7_08.tif": None,
    "dose4_02.tif": None, "dose4_04.tif": None, "dose3_04.tif": None, "dose3_09.tif": None,
    "dose6_05.tif": None, "dose5_04.tif": None, "dose8_10.tif": None, "dose1_07.tif": None,
}


def period_profile(img, x0, x1):
    col = img[:, x0:x1]
    # Y period from spectrum of column-mean profile
    prof = col.mean(1); prof = prof - prof.mean()
    f = np.abs(np.fft.rfft(prof)); k = np.arange(len(f))
    mask = (k >= len(prof) / 500) & (k <= len(prof) / 14)
    kp = int(np.argmax(np.where(mask, f, 0))); period = len(prof) / kp
    # fold every column's profile at the period with the column's own best phase
    n = int(round(period)); acc = np.zeros(n); cnt = 0
    for j in range(col.shape[1]):
        p = col[:, j]
        m = (len(p) // n) * n
        folded = p[:m].reshape(-1, n).mean(0)
        acc += folded; cnt += 1
    mean = acc / cnt
    # rotate so the floor minimum sits at index 0
    mean = np.roll(mean, -int(np.argmin(mean)))
    return period, mean


def floor_gap_fraction(mean):
    lo, hi = mean.min(), np.percentile(mean, 95)
    thr = lo + 0.18 * (hi - lo)
    # contiguous floor run around index 0 (already the minimum)
    n = len(mean); below = mean < thr
    run = 0
    for i in range(n):
        if below[i]: run += 1
        else: break
    for i in range(n - 1, 0, -1):
        if below[i]: run += 1
        else: break
    return run / n, thr


def main():
    manifest = {f["file"]: f for f in json.loads((SEM / "manifest.json").read_text())["frames"]}
    out = []
    for fname, win in WINDOWS.items():
        fr = manifest[fname]; px = fr["meta"]["pixel_size_nm"]; th = fr["meta"]["stage_tilt_deg"]
        img = np.asarray(Image.open(SEM / "originals" / fname).convert("L"), float)[:FOOTER]
        x0, x1 = (0, img.shape[1]) if win is None else win
        period, mean = period_profile(img, x0, x1)
        frac, thr = floor_gap_fraction(mean)
        p_proj_nm = period * px                      # projected pitch along Y
        p_true_nm = p_proj_nm / np.cos(np.radians(th))
        gap_nm = frac * p_proj_nm
        out.append({"file": fname, "dose": fr["dose_index"], "tilt_deg": th, "pixel_nm": px,
                    "period_px": period, "pitch_true_nm": p_true_nm, "floor_gap_nm": gap_nm,
                    "occluded_extent_nm": p_proj_nm - gap_nm, "profile": mean.tolist()})
        print(f"{fname:<18} T{th:>4.0f} pitch {p_true_nm:5.0f} nm  floor gap {gap_nm:5.0f} nm  occluded (a cos t + h sin t) {p_proj_nm-gap_nm:5.0f} nm")
    # solve per dose with both tilts
    print()
    for dose in sorted({r["dose"] for r in out}):
        rs = [r for r in out if r["dose"] == dose]
        t10 = [r for r in rs if r["tilt_deg"] < 12]; t15 = [r for r in rs if r["tilt_deg"] >= 12]
        if t10 and t15:
            E10 = np.median([r["occluded_extent_nm"] for r in t10]); E15 = np.median([r["occluded_extent_nm"] for r in t15])
            c10, s10 = np.cos(np.radians(10)), np.sin(np.radians(10)); c15, s15 = np.cos(np.radians(15)), np.sin(np.radians(15))
            A = np.array([[c10, s10], [c15, s15]]); a, h = np.linalg.solve(A, [E10, E15])
            print(f"dose{dose}: two-tilt solve -> top extent a = {a:.0f} nm, height h = {h:.0f} nm   (E10={E10:.0f}, E15={E15:.0f})")
        for r in rs:
            # single-tilt with a = tip-to-tip from plan measurement (dose summary)
            pass
    summ = json.loads((SEM / "dose_summary.json").read_text())
    print()
    for r in out:
        a = summ.get(str(r["dose"]), {}).get("cd_x_extent_nm", [None])[0]
        if a:
            th = np.radians(r["tilt_deg"]); h = (r["occluded_extent_nm"] - a * np.cos(th)) / np.sin(th)
            r["height_from_tip_to_tip_nm"] = h
            print(f"{r['file']:<18} a(tip-to-tip)={a:.0f}  -> h = {h:.0f} nm")
    (SEM / "occlusion_height.json").write_text(json.dumps(out, indent=1))
    # figure: mean period profiles for dose9 at both tilts with the floor run marked
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), dpi=110)
    for ax, key in zip(axs, ["dose9_04.tif", "dose9_15d_02.tif"]):
        r = [x for x in out if x["file"] == key][0]; m = np.array(r["profile"]); y = np.arange(len(m)) * r["pixel_nm"]
        ax.plot(y, m, lw=1.5); frac, thr = floor_gap_fraction(m); ax.axhline(thr, color="orange", ls="--", lw=1)
        ax.axvspan(0, r["floor_gap_nm"], color="orange", alpha=.15)
        ax.set_title(f"{key} tilt {r['tilt_deg']:.0f}°: one period along Y; floor gap {r['floor_gap_nm']:.0f} nm")
        ax.set_xlabel("Y within period (nm, projected)"); ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(SEM / "occlusion_profiles_dose9.png"); print("figure saved")


if __name__ == "__main__":
    main()
