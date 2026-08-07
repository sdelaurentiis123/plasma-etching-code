#!/usr/bin/env python3
"""Digitize the pure-Ar curve in Lee--Lieberman (1994) Figure 3.

The source plot is a scanned log--log figure.  Axis centers and coarse curve
seeds were selected by full-resolution visual inspection.  PIL then recenters
each seed on the dark-pixel stroke, making the final pixel-to-data conversion
deterministic and auditable.  The copyrighted source raster and QA overlay are
not committed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "research_sources" / "digitized"
CSV_PATH = OUTPUT_DIRECTORY / "lee_lieberman_1994_figure3_argon.csv"
MANIFEST_PATH = (
    OUTPUT_DIRECTORY / "lee_lieberman_1994_figure3_argon_manifest.json")

SOURCE_RENDER_SHA256 = (
    "e8f26ee36f1fde182ff4c1f6b7a74ac6f8b4f2791d009b26b978bcc2b2713fae")
SOURCE_RENDER_SIZE = (1990, 1527)
VISUALLY_INSPECTED_OVERLAY_SHA256 = (
    "ea4be63f42e71933683eb340d2ad8fe0f1254132563228dcf92fe885b2d9ce2c")

# Bottom/left log-axis calibration.  Interior decade ticks at x =
# 720, 1051, and 1382 px give 331 px/decade.  The 1 and 10 eV y ticks
# are centered at y = 1268 and 487 px, respectively.
X_AT_0P1_MTORR = 389.0
X_PIXELS_PER_DECADE = 331.0
Y_AT_1_EV = 1268.0
Y_PIXELS_PER_DECADE = 781.0

# Pressure samples and visually selected coarse Ar-stroke centers.  The
# recentering routine below uses the raster, not these integer seeds, for the
# final y coordinate.  Above about 70 mTorr the Ar and O2 strokes overlap.
CURVE_SEEDS = (
    (0.4, 613, False),
    (0.5, 644, False),
    (0.7, 686, False),
    (1.0, 725, False),
    (1.5, 767, False),
    (2.0, 793, False),
    (3.0, 828, False),
    (5.0, 868, False),
    (7.0, 892, False),
    (10.0, 917, False),
    (15.0, 944, False),
    (20.0, 963, False),
    (30.0, 988, False),
    (50.0, 1020, False),
    (70.0, 1037, True),
    (100.0, 1062, True),
    (150.0, 1087, True),
    (200.0, 1104, True),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def pressure_to_x(pressure_mTorr: float) -> float:
    return float(
        X_AT_0P1_MTORR
        + X_PIXELS_PER_DECADE * (np.log10(pressure_mTorr) + 1.0))


def y_to_temperature_eV(y_pixel: float) -> float:
    return float(10.0 ** (
        (Y_AT_1_EV - y_pixel) / Y_PIXELS_PER_DECADE))


def _dark_group_near(
        gray: np.ndarray, *, x_pixel: int, y_seed: int,
        half_window_y: int = 18, threshold: int = 160,
        ) -> tuple[float, float]:
    """Return darkness-weighted center and half-height near one curve seed."""
    column_slice = gray[:, max(0, x_pixel - 2):x_pixel + 3]
    darkness_count = np.sum(column_slice < threshold, axis=1)
    lower = max(0, y_seed - half_window_y)
    upper = min(gray.shape[0], y_seed + half_window_y + 1)
    indices = np.flatnonzero(darkness_count[lower:upper] > 0) + lower
    if indices.size == 0:
        raise RuntimeError(
            f"no dark stroke near pixel ({x_pixel}, {y_seed})")
    groups = np.split(indices, np.flatnonzero(np.diff(indices) > 1) + 1)
    group = max(
        groups,
        key=lambda values: int(np.sum(darkness_count[values])),
    )
    weights = darkness_count[group].astype(float)
    center = float(np.average(group, weights=weights))
    half_height = float(0.5 * (group[-1] - group[0] + 1))
    return center, half_height


def digitize(image: Image.Image) -> list[dict[str, object]]:
    gray = np.asarray(image.convert("L"))
    rows: list[dict[str, object]] = []
    for pressure, seed_y, overlap in CURVE_SEEDS:
        x_float = pressure_to_x(pressure)
        center_y, half_height_y = _dark_group_near(
            gray, x_pixel=int(round(x_float)), y_seed=seed_y)
        rows.append({
            "pressure_mTorr": pressure,
            "electron_temperature_eV": y_to_temperature_eV(center_y),
            "pixel_x": x_float,
            "pixel_y": center_y,
            "stroke_half_height_pixels": half_height_y,
            "argon_oxygen_strokes_overlap": overlap,
        })
    return rows


def _write_csv(rows: list[dict[str, object]]) -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_manifest(
        source: Path, rows: list[dict[str, object]],
        overlay_sha256: str | None) -> None:
    manifest = {
        "source": {
            "bibkey": "lee-lieberman-1994-global",
            "figure": "Figure 3, pure-Ar solid heavy curve",
            "source_pdf_sha256":
                "f2049e7041984d658d23688e8e8112a8d8e8a524172a8d2e335be8fde7fc2e23",
            "render_sha256": SOURCE_RENDER_SHA256,
            "render_size_pixels": list(SOURCE_RENDER_SIZE),
            "render_resolution_dpi": 180,
            "source_raster_committed": False,
        },
        "conditions": {
            "absorbed_power_W": 1000.0,
            "flow_sccm": 35.0,
            "surface_recombination_coefficient": 0.0,
            "reactor_radius_m": 0.1525,
            "reactor_length_m": 0.075,
            "gas_temperature_K": 600.0,
        },
        "axis_calibration": {
            "x_scale": "log10 pressure in mTorr",
            "x_at_0.1_mTorr_pixels": X_AT_0P1_MTORR,
            "x_pixels_per_decade": X_PIXELS_PER_DECADE,
            "x_interior_major_ticks_pixels": [720.0, 1051.0, 1382.0],
            "y_scale": "log10 electron temperature in eV",
            "y_at_1_eV_pixels": Y_AT_1_EV,
            "y_at_10_eV_pixels": 487.0,
            "y_pixels_per_decade": Y_PIXELS_PER_DECADE,
        },
        "method": (
            "Full-resolution visual curve identification and axis calibration; "
            "PIL five-column dark-stroke recentering at preregistered pressure "
            "samples; explicit overlap flag where Ar and O2 strokes merge."
        ),
        "uncertainty": {
            "isolated_stroke_relative_temperature_percent": 2.0,
            "overlapped_stroke_relative_temperature_percent": 3.5,
            "basis": (
                "stroke half-height plus a two-pixel axis-center allowance; "
                "overlap points carry an additional centerline ambiguity"
            ),
        },
        "row_count": len(rows),
        "outputs": {
            "csv": str(CSV_PATH.relative_to(ROOT)),
            "overlay_committed": False,
            "overlay_sha256": overlay_sha256,
        },
        "visual_inspection": {
            "status": (
                "passed"
                if overlay_sha256 == VISUALLY_INSPECTED_OVERLAY_SHA256
                else "pending"
            ),
            "inspection_date": (
                "2026-08-07"
                if overlay_sha256 == VISUALLY_INSPECTED_OVERLAY_SHA256
                else None
            ),
            "inspection_resolution_pixels": (
                list(SOURCE_RENDER_SIZE)
                if overlay_sha256 == VISUALLY_INSPECTED_OVERLAY_SHA256
                else None
            ),
            "finding": (
                "All markers centered on the heavy Ar stroke; orange markers "
                "correctly flag the Ar/O2 overlap region."
                if overlay_sha256 == VISUALLY_INSPECTED_OVERLAY_SHA256
                else "Inspect the full-resolution overlay before a pass."
            ),
        },
        "rights": (
            "Only derived coordinates and calibration metadata are committed; "
            "the source raster and annotated overlay remain local."
        ),
        "source_render_input": "local non-redistributed source render",
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _draw_overlay(
        image: Image.Image, rows: list[dict[str, object]], output: Path) -> str:
    overlay = image.convert("RGB")
    draw = ImageDraw.Draw(overlay)
    for x in (389.0, 720.0, 1051.0, 1382.0):
        draw.line((x, 240, x, 1280), fill=(0, 180, 255), width=1)
    for y in (487.0, 1268.0):
        draw.line((380, y, 1730, y), fill=(0, 180, 255), width=1)
    for row in rows:
        x = float(row["pixel_x"])
        y = float(row["pixel_y"])
        color = (
            (255, 140, 0)
            if bool(row["argon_oxygen_strokes_overlap"])
            else (255, 0, 180)
        )
        draw.ellipse((x - 7, y - 7, x + 7, y + 7), outline=color, width=2)
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)
    return _sha256(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", type=Path, required=True)
    parser.add_argument("--overlay", type=Path)
    arguments = parser.parse_args()

    if _sha256(arguments.render) != SOURCE_RENDER_SHA256:
        raise RuntimeError("source render hash does not match audited Figure 3")
    image = Image.open(arguments.render)
    if image.size != SOURCE_RENDER_SIZE:
        raise RuntimeError(
            f"unexpected source render size {image.size!r}")
    rows = digitize(image)
    _write_csv(rows)
    overlay_sha256 = None
    if arguments.overlay is not None:
        overlay_sha256 = _draw_overlay(image, rows, arguments.overlay)
    _write_manifest(arguments.render, rows, overlay_sha256)


if __name__ == "__main__":
    main()
