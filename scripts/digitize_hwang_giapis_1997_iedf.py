#!/usr/bin/env python3
"""Reproduce the Fig. 4(a) IEDF digitization used by the Nozawa replay.

The extraction follows the continuous black curve in the 240-dpi rendering.
Calibration points are stated in source-image pixels and the output retains the
unnormalized curve height so a reviewer can audit the normalized bin masses.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "refs" / "HG_jvstb97.pdf"
PAGE_IMAGE = ROOT / "tmp" / "pdfs" / "hwang_giapis_1997" / "source_pages" / "page-004.png"
OUTPUT = ROOT / "data" / "experimental" / "hwang_giapis_1997"

# Source rendering calibration, not display coordinates.  The energy tick
# centers at 0 and 80 eV define x(E); the zero-IEDF baseline defines y=0.
X_AT_0_EV = 1301.6
X_AT_80_EV = 1731.0
Y_AT_ZERO = 709.5
CURVE_Y_MIN = 275
CURVE_Y_MAX_EXCLUSIVE = 708
DARK_THRESHOLD = 180
MAXIMUM_COLUMN_STEP_PX = 35.0


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def render_page():
    PAGE_IMAGE.parent.mkdir(parents=True, exist_ok=True)
    if not PAGE_IMAGE.exists():
        subprocess.run([
            "pdftoppm", "-f", "4", "-l", "4", "-png", "-r", "240",
            "-singlefile", str(PDF), str(PAGE_IMAGE.with_suffix("")),
        ], check=True)


def follow_curve(image):
    x_min = int(np.floor(X_AT_0_EV))
    x_max = int(np.ceil(X_AT_80_EV)) + 1
    previous_y = 706.0
    curve_y = []
    for x in range(x_min, x_max):
        candidates = np.flatnonzero(
            image[CURVE_Y_MIN:CURVE_Y_MAX_EXCLUSIVE, x] < DARK_THRESHOLD)
        candidates = candidates + CURVE_Y_MIN
        if candidates.size:
            distance = np.abs(candidates - previous_y)
            nearest = float(distance.min())
            if nearest <= MAXIMUM_COLUMN_STEP_PX:
                previous_y = float(np.median(candidates[distance <= nearest + 2.0]))
        curve_y.append(previous_y)
    return np.arange(x_min, x_max, dtype=float), np.asarray(curve_y)


def main():
    render_page()
    image = np.asarray(Image.open(PAGE_IMAGE).convert("L"))
    if image.shape != (2640, 2040):
        raise RuntimeError(f"unexpected 240-dpi page shape: {image.shape}")
    source_x, source_y = follow_curve(image)
    lower = np.arange(0.0, 80.0, 1.0)
    upper = lower + 1.0
    center = 0.5 * (lower + upper)
    sample_x = X_AT_0_EV + center * (X_AT_80_EV - X_AT_0_EV) / 80.0
    height = Y_AT_ZERO - np.interp(sample_x, source_x, source_y)
    height = np.convolve(
        np.pad(height, (1, 1), mode="edge"), np.ones(3) / 3.0, mode="valid")
    height = np.maximum(height, 0.0)
    probability = height / height.sum()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT / "fig4a_ion_energy_distribution.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "normal_energy_lower_eV", "normal_energy_upper_eV",
            "probability_mass", "digitized_curve_height_px",
            "source_pdf_sha256", "source_pdf_page", "source_figure",
        ])
        for values in zip(lower, upper, probability, height):
            writer.writerow([
                f"{values[0]:.1f}", f"{values[1]:.1f}",
                f"{values[2]:.12g}", f"{values[3]:.9g}",
                sha256(PDF), 4, "Fig. 4(a)",
            ])
    manifest = {
        "manifest_id": "HWANG-GIAPIS-1997-FIG4A-IEDF-R1",
        "source": {
            "citation": "Hwang & Giapis, JVST B 15, 70 (1997)",
            "doi": "10.1116/1.589258",
            "pdf_sha256": sha256(PDF),
            "pdf_page": 4,
            "print_page": 73,
            "figure": "Fig. 4(a)",
            "render_dpi": 240,
            "rendered_page_sha256": sha256(PAGE_IMAGE),
        },
        "pixel_calibration": {
            "x_at_0_eV": X_AT_0_EV,
            "x_at_80_eV": X_AT_80_EV,
            "y_at_zero_iedf": Y_AT_ZERO,
            "curve_y_search": [CURVE_Y_MIN, CURVE_Y_MAX_EXCLUSIVE],
            "dark_threshold": DARK_THRESHOLD,
            "maximum_column_step_px": MAXIMUM_COLUMN_STEP_PX,
            "energy_bin_width_eV": 1.0,
            "smoothing": "three-bin centered moving average with endpoint replication",
        },
        "digitization_uncertainty": {
            "energy_eV": 0.5,
            "curve_vertical_px": 2.0,
            "note": "curve mass is normalized after extraction; original paper reports no numerical table",
        },
        "derived_checks": {
            "probability_mass_sum": float(probability.sum()),
            "mean_normal_energy_eV": float(np.dot(center, probability)),
            "low_energy_mass_below_39_eV": float(probability[center < 39.0].sum()),
            "low_to_high_peak_height_ratio": float(
                height[(center >= 3.0) & (center <= 12.0)].max()
                / height[(center >= 62.0) & (center <= 74.0)].max()),
        },
        "output": {
            "path": str(csv_path.relative_to(ROOT)),
            "sha256": sha256(csv_path),
        },
    }
    (OUTPUT / "digitization_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest["derived_checks"], indent=2))


if __name__ == "__main__":
    main()
