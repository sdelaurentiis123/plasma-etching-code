#!/usr/bin/env python3
"""Reproduce the Fig. 4(b) electron-energy digitization used by the 2-D replay.

The source paper sampled the electron energy distribution and the Fig. 5(b)
angular distribution independently at the sheath lower boundary.  This script
extracts only the Fig. 4(b) energy curve.  The angular law remains the paper's
explicit ``cos(theta)**0.6`` fit and is represented analytically in the engine.
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
PAGE_IMAGE = (
    ROOT / "tmp" / "pdfs" / "hwang_giapis_1997"
    / "source_pages" / "page-004.png")
OUTPUT = ROOT / "data" / "experimental" / "hwang_giapis_1997"

# Source-render pixel calibration at 240 dpi.  The curve becomes thinner than
# the paper/rendering resolution above roughly 12 eV; heights below 15 pixels
# are therefore declared unresolved and set to zero instead of following text
# or plot-frame pixels.
X_AT_0_EV = 1247.0
X_AT_20_EV = 1783.0
Y_AT_ZERO = 1416.5
CURVE_X_MIN = 1254
CURVE_X_MAX_EXCLUSIVE = 1784
CURVE_Y_MIN = 900
CURVE_Y_MAX_EXCLUSIVE = 1408
DARK_THRESHOLD = 170
MAXIMUM_COLUMN_STEP_PX = 35.0
MINIMUM_RESOLVED_HEIGHT_PX = 15.0


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
    previous_y = 1200.0
    curve_y = []
    for x in range(CURVE_X_MIN, CURVE_X_MAX_EXCLUSIVE):
        candidates = np.flatnonzero(
            image[CURVE_Y_MIN:CURVE_Y_MAX_EXCLUSIVE, x] < DARK_THRESHOLD)
        candidates = candidates + CURVE_Y_MIN
        if candidates.size:
            distance = np.abs(candidates - previous_y)
            nearest = float(distance.min())
            if nearest <= MAXIMUM_COLUMN_STEP_PX:
                previous_y = float(np.median(
                    candidates[distance <= nearest + 2.0]))
        curve_y.append(previous_y)
    return (
        np.arange(CURVE_X_MIN, CURVE_X_MAX_EXCLUSIVE, dtype=float),
        np.asarray(curve_y),
    )


def main():
    render_page()
    image = np.asarray(Image.open(PAGE_IMAGE).convert("L"))
    if image.shape != (2640, 2040):
        raise RuntimeError(f"unexpected 240-dpi page shape: {image.shape}")
    source_x, source_y = follow_curve(image)
    lower = np.arange(0.0, 20.0, 0.25)
    upper = lower + 0.25
    center = 0.5 * (lower + upper)
    sample_x = (
        X_AT_0_EV
        + center * (X_AT_20_EV - X_AT_0_EV) / 20.0)
    height = Y_AT_ZERO - np.interp(sample_x, source_x, source_y)
    height = np.convolve(
        np.pad(height, (2, 2), mode="edge"),
        np.ones(5) / 5.0,
        mode="valid")
    height = np.where(
        height >= MINIMUM_RESOLVED_HEIGHT_PX, height, 0.0)
    if height.sum() <= 0.0:
        raise RuntimeError("digitized EEDF has no resolved probability mass")
    probability = height / height.sum()

    OUTPUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT / "fig4b_electron_energy_distribution.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "energy_lower_eV", "energy_upper_eV", "probability_mass",
            "digitized_curve_height_px", "source_pdf_sha256",
            "source_pdf_page", "source_figure",
        ])
        for values in zip(lower, upper, probability, height):
            writer.writerow([
                f"{values[0]:.2f}", f"{values[1]:.2f}",
                f"{values[2]:.12g}", f"{values[3]:.9g}",
                sha256(PDF), 4, "Fig. 4(b)",
            ])

    manifest = {
        "manifest_id": "HWANG-GIAPIS-1997-FIG4B-EEDF-R1",
        "source": {
            "citation": "Hwang & Giapis, JVST B 15, 70 (1997)",
            "doi": "10.1116/1.589258",
            "pdf_sha256": sha256(PDF),
            "pdf_page": 4,
            "print_page": 73,
            "figure": "Fig. 4(b)",
            "render_dpi": 240,
            "rendered_page_sha256": sha256(PAGE_IMAGE),
        },
        "pixel_calibration": {
            "x_at_0_eV": X_AT_0_EV,
            "x_at_20_eV": X_AT_20_EV,
            "y_at_zero_eedf": Y_AT_ZERO,
            "curve_x_search": [CURVE_X_MIN, CURVE_X_MAX_EXCLUSIVE],
            "curve_y_search": [CURVE_Y_MIN, CURVE_Y_MAX_EXCLUSIVE],
            "dark_threshold": DARK_THRESHOLD,
            "maximum_column_step_px": MAXIMUM_COLUMN_STEP_PX,
            "minimum_resolved_height_px": MINIMUM_RESOLVED_HEIGHT_PX,
            "energy_bin_width_eV": 0.25,
            "smoothing": (
                "five-bin centered moving average with endpoint replication"),
            "tail_policy": (
                "curve heights below the declared pixel-resolution threshold "
                "are zero; no tail law is invented"),
        },
        "digitization_uncertainty": {
            "energy_eV": 0.125,
            "curve_vertical_px": 3.0,
            "note": (
                "curve mass is normalized after extraction; the original "
                "paper reports no numerical EEDF table"),
        },
        "derived_checks": {
            "probability_mass_sum": float(probability.sum()),
            "mean_energy_eV": float(np.dot(center, probability)),
            "probability_mass_below_5_eV": float(
                probability[center < 5.0].sum()),
            "highest_resolved_energy_eV": float(
                upper[np.flatnonzero(height)[-1]]),
        },
        "angular_companion": {
            "source_figure": "Fig. 5(b)",
            "declared_fit": "cos(theta)**0.6",
            "digitized": False,
            "reason": "the paper provides the analytic fit explicitly",
        },
        "output": {
            "path": str(csv_path.relative_to(ROOT)),
            "sha256": sha256(csv_path),
        },
    }
    manifest_path = OUTPUT / "fig4b_digitization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest["derived_checks"], indent=2))


if __name__ == "__main__":
    main()
