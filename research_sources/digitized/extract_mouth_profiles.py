"""Digitize the mask-region aperture profile from Krüger 2024 Fig. 7.

Fig. 7(a) is an actual MCFPM output (converged single-feature parameter set,
60 s); Fig. 7(b) is the base-case experimental SEM used to set the target
metrics. Both are cropped out of the published PDF at 600 dpi.

Crops (reproduce with pdftoppm on tmp/pdfs/krueger-2024.pdf, page 8):
  sim: -r 600 -x 920  -y 2650 -W 800 -H 2250
  sem: -r 600 -x 1740 -y 2650 -W 320 -H 2250

Vertical calibration: the AC mask is 850 nm thick (paper Sec. V, thesis
p. 181). Horizontal calibration: the aperture at the mask/oxide interface is
w_t = 90 nm (Table IV target, matched to <1 % by the converged run).

Outputs two CSVs of aperture-vs-depth-into-mask.
"""
from __future__ import annotations

import csv
import pathlib

import numpy as np
from PIL import Image

HERE = pathlib.Path(__file__).resolve().parent


def simulated_profile() -> list[tuple[float, float]]:
    im = np.array(Image.open(HERE / "krueger_fig7a_mcfpm_crop.png").convert("RGB")).astype(int)
    r, g, b = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    magenta = (r > 200) & (b > 200) & (g < 170)          # AC mask
    blue = (b - np.maximum(r, g)) > 40                    # fluorocarbon polymer
    cyan = (np.minimum(g, b) - r) > 40                    # SiO2
    material = magenta | blue | cyan
    cx = 327                                              # feature centre column
    mask_top, mask_bottom = 60, 1150                      # rows bounding the AC mask
    rows: list[tuple[int, int]] = []
    for y in range(mask_top + 10, mask_bottom - 10):
        idx = np.where(material[y])[0]
        left, right = idx[idx < cx], idx[idx > cx]
        if left.size == 0 or right.size == 0:
            continue
        rows.append((y, int(right.min() - left.max() - 1)))
    base = float(np.median([n for y, n in rows if y > mask_bottom - 80]))
    px_per_nm = (mask_bottom - mask_top) / 850.0
    return [((y - mask_top) / px_per_nm, n * 90.0 / base) for y, n in rows]


def experimental_profile() -> list[tuple[float, float]]:
    im = np.array(Image.open(HERE / "krueger_fig7b_sem_crop.png").convert("L")).astype(float)
    h, w = im.shape
    mask_top, mask_bottom = 85, 1500                      # from the row-brightness step
    nm_per_px = 850.0 / (mask_bottom - mask_top)
    out: list[tuple[float, float]] = []
    for y in range(mask_top + 25, mask_bottom):
        row = im[max(0, y - 6):y + 7, :].mean(axis=0)
        c = 175                                           # feature centre column
        li = int(np.argmax(row[60:c - 5]) + 60)           # bright left sidewall
        ri = int(np.argmax(row[c + 5:305]) + c + 5)       # bright right sidewall
        interior = row[li:ri + 1]
        half = (interior.min() + min(row[li], row[ri])) / 2.0
        cen = int(np.argmin(interior)) + li
        left = cen
        while left > li and row[left] < half:
            left -= 1
        right = cen
        while right < ri and row[right] < half:
            right += 1
        out.append(((y - mask_top) * nm_per_px, (right - left) * nm_per_px))
    return out


def _write(name: str, rows: list[tuple[float, float]]) -> None:
    with open(HERE / name, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["depth_into_mask_nm", "aperture_nm"])
        for depth, aperture in rows:
            writer.writerow([f"{depth:.2f}", f"{aperture:.2f}"])


if __name__ == "__main__":
    _write("krueger_fig7a_simulated_aperture.csv", simulated_profile())
    _write("krueger_fig7b_experimental_aperture.csv", experimental_profile())
    print("wrote CSVs to", HERE)
