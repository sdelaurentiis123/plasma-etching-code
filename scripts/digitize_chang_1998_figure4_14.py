#!/usr/bin/env python3
"""Checksum-backed PIL replay of Chang Figure 4.14 oxide markers.

The source pixels are not redistributed.  This script renders the locally
versioned thesis at 600 dpi, checks the exact render, locates the dense cores
of the seven filled-square SiO2 markers, and maps their centers through the
printed axes.  Polysilicon triangles and both drawn guide curves are ignored.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PDF = ROOT / "research_sources" / "chang_thesis.pdf"
DEFAULT_OUTPUT = ROOT / "data" / "experimental" / "chang_1998_figure4_14"
SOURCE_PDF_SHA256 = (
    "ef5c511a7fd8e6d6a0bf721874e57ba18194184bf3db90457188ed42c4bd3b4b"
)
RENDER_SHA256 = (
    "766bcd772ffbc8fcd9dc9b37be3fbe5bbc8f9182ebcad2c13b4415761bed78bc"
)
PDF_PAGE = 99
PRINT_PAGE = 99
RENDER_DPI = 600
RENDER_SIZE = (5100, 6600)

# Full-range axis intersections on the original render.
X_ZERO_PX = 1540.0
X_NINETY_PX = 3264.0
Y_ZERO_PX = 4160.0
Y_POINT_TWO_PX = 2538.0

# Visual seeds checked at original resolution.  The reported values are not
# taken from these coordinates: each center is recovered from dense source
# pixels in a +/-60 px neighborhood.
NOMINAL_ANGLES_DEG = np.asarray((0, 15, 30, 45, 60, 75, 90), dtype=float)
MARKER_SEEDS_XY_PX = (
    (1540, 3561), (1822, 3562), (2108, 3545), (2394, 3424),
    (2682, 3014), (2969, 3175), (3259, 4153),
)
DARK_THRESHOLD = 90
DENSE_WINDOW_PX = 11
DENSE_WINDOW_MINIMUM_DARK_PIXELS = 115


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _render(pdf: Path, directory: Path) -> Path:
    prefix = directory / "chang_figure4_14_page099"
    subprocess.run(
        [
            "pdftoppm", "-f", str(PDF_PAGE), "-l", str(PDF_PAGE),
            "-r", str(RENDER_DPI), "-png", "-singlefile", str(pdf),
            str(prefix),
        ],
        check=True,
    )
    return prefix.with_suffix(".png")


def _dense_marker_center(gray: np.ndarray, seed_xy):
    seed_x, seed_y = seed_xy
    radius = 60
    region = (
        gray[seed_y - radius:seed_y + radius + 1,
             seed_x - radius:seed_x + radius + 1]
        < DARK_THRESHOLD
    )
    windows = np.lib.stride_tricks.sliding_window_view(
        region, (DENSE_WINDOW_PX, DENSE_WINDOW_PX)
    )
    density = np.sum(windows, axis=(-1, -2))
    local_y, local_x = np.nonzero(
        density >= DENSE_WINDOW_MINIMUM_DARK_PIXELS
    )
    if local_x.size < 400:
        raise ValueError(f"oxide marker pixel audit failed at seed {seed_xy}")
    offset = DENSE_WINDOW_PX // 2
    full_x = local_x + seed_x - radius + offset
    full_y = local_y + seed_y - radius + offset
    return {
        "center_x_px": float(np.median(full_x)),
        "center_y_px": float(np.median(full_y)),
        "dense_pixel_count": int(full_x.size),
        "dense_x_span_px": [int(np.min(full_x)), int(np.max(full_x))],
        "dense_y_span_px": [int(np.min(full_y)), int(np.max(full_y))],
    }


def _extract(render: Path):
    if _hash(render) != RENDER_SHA256:
        raise ValueError("600-dpi Figure 4.14 page render checksum changed")
    image = Image.open(render).convert("L")
    if image.size != RENDER_SIZE:
        raise ValueError(f"unexpected rendered page size: {image.size}")
    gray = np.asarray(image)
    markers = [
        _dense_marker_center(gray, seed) for seed in MARKER_SEEDS_XY_PX
    ]
    x_scale = 90.0 / (X_NINETY_PX - X_ZERO_PX)
    y_scale = 0.2 / (Y_ZERO_PX - Y_POINT_TWO_PX)
    extracted = []
    for nominal, marker in zip(NOMINAL_ANGLES_DEG, markers):
        pixel_angle = (marker["center_x_px"] - X_ZERO_PX) * x_scale
        if abs(pixel_angle - nominal) > 0.6:
            raise ValueError(
                f"marker x-position does not reconcile with {nominal:g} deg"
            )
        yield_value = (Y_ZERO_PX - marker["center_y_px"]) * y_scale
        # The 90-degree square is clipped by and centered on the zero axis;
        # dense-pixel centering alone gives the upper half-pixel inventory.
        axis_clipped = bool(nominal == 90.0)
        if axis_clipped:
            yield_value = 0.0
        extracted.append({
            "angle_deg": float(nominal),
            "pixel_recovered_angle_deg": float(pixel_angle),
            "yield_sio2_formula_per_ar_ion": float(max(yield_value, 0.0)),
            "axis_clipped_zero_marker": axis_clipped,
            **marker,
        })
    return extracted


def _write(output: Path, extracted):
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "sio2_angular_yield_100eV_r90.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "angle_deg", "pixel_recovered_angle_deg",
            "yield_sio2_formula_per_ar_ion", "axis_clipped_zero_marker",
            "center_x_px", "center_y_px", "dense_pixel_count",
        ))
        writer.writeheader()
        for row in extracted:
            writer.writerow({name: row[name] for name in writer.fieldnames})

    manifest = {
        "schema": "petch.chang-1998-figure4.14-sio2-angular-yield.v1",
        "source": {
            "citation": (
                "J. P. Chang, Study of plasma-surface kinetics and simulation "
                "of feature profile evolution in chlorine etching of patterned "
                "polysilicon, MIT PhD thesis (1998), Figure 4.14"
            ),
            "pdf_path": "research_sources/chang_thesis.pdf",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": PDF_PAGE,
            "print_page": PRINT_PAGE,
            "render_dpi": RENDER_DPI,
            "render_size_px": list(RENDER_SIZE),
            "render_sha256": RENDER_SHA256,
        },
        "experiment": {
            "projectile": "Ar+",
            "target": "SiO2",
            "ion_energy_eV": 100.0,
            "atomic_cl_to_ar_ion_flux_ratio": 90.0,
            "yield_unit": "removed SiO2 formula unit per incident Ar+",
        },
        "pixel_calibration": {
            "angle_zero_x_px": X_ZERO_PX,
            "angle_ninety_x_px": X_NINETY_PX,
            "yield_zero_y_px": Y_ZERO_PX,
            "yield_point_two_y_px": Y_POINT_TWO_PX,
            "dark_threshold": DARK_THRESHOLD,
            "dense_window_px": DENSE_WINDOW_PX,
            "dense_window_minimum_dark_pixels": (
                DENSE_WINDOW_MINIMUM_DARK_PIXELS),
        },
        "extraction": {
            "method": (
                "original-resolution visual seeds followed by PIL/NumPy "
                "dense-square-core localization; guide curves ignored"
            ),
            "markers": extracted,
            "ninety_degree_handling": (
                "filled square is clipped by the zero axis; the plotted zero "
                "and axis intersection are authoritative"
            ),
            "feature_depth_used": False,
            "reactor_state_used": False,
        },
        "scope": {
            "valid_use": (
                "incidence-angle interpolation for 100 eV Ar+/Cl etching of "
                "SiO2 at atomic-Cl/Ar+=90"
            ),
            "invalid_use": (
                "Cl+ or Cl2+ projectile identity, energy scaling, fluorocarbon "
                "oxide etching, or target-reactor calibration"
            ),
        },
        "output": {
            "csv_path": (
                "data/experimental/chang_1998_figure4_14/"
                "sio2_angular_yield_100eV_r90.csv"
            ),
            "csv_sha256": _hash(csv_path),
        },
    }
    manifest_path = output / "digitization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return csv_path, manifest_path


def _overlay(render: Path, path: Path, extracted):
    image = Image.open(render).convert("RGB")
    draw = ImageDraw.Draw(image)
    draw.line(
        [(X_ZERO_PX, Y_ZERO_PX), (X_NINETY_PX, Y_ZERO_PX)],
        fill=(0, 180, 0), width=8,
    )
    draw.line(
        [(X_ZERO_PX, Y_ZERO_PX), (X_ZERO_PX, Y_POINT_TWO_PX)],
        fill=(0, 180, 0), width=8,
    )
    for row in extracted:
        x = row["center_x_px"]
        y = row["center_y_px"]
        draw.ellipse((x - 18, y - 18, x + 18, y + 18), outline=(220, 0, 0), width=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.crop((1200, 2300, 3550, 4450)).save(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if _hash(arguments.pdf) != SOURCE_PDF_SHA256:
        raise ValueError("Chang thesis PDF checksum changed")
    with tempfile.TemporaryDirectory(prefix="petch-chang-fig414-") as directory:
        render = _render(arguments.pdf, Path(directory))
        extracted = _extract(render)
        csv_path, manifest_path = _write(arguments.output, extracted)
        if arguments.overlay is not None:
            _overlay(render, arguments.overlay, extracted)
    if arguments.check:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["output"]["csv_sha256"] != _hash(csv_path):
            raise ValueError("written Figure 4.14 CSV checksum mismatch")
    print(json.dumps({
        "csv": str(csv_path),
        "manifest": str(manifest_path),
        "rows": extracted,
    }, indent=2))


if __name__ == "__main__":
    main()
