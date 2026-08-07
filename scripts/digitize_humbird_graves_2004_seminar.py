#!/usr/bin/env python3
"""Replay the Humbird--Graves 2004 seminar surface-state digitization.

The source is a primary-author UC Berkeley seminar, not the peer-reviewed
journal figure.  The raster images are not redistributed.  This script:

1. verifies the two source-image hashes and dimensions;
2. verifies every selected trace center against full-resolution RGB pixels;
3. maps the frozen pixel centers through explicit linear axis calibrations;
4. reproduces the committed CSV and manifest; and
5. optionally draws full-resolution QA overlays.

The selected points intentionally avoid automatic curve fitting.  They are
fixed, sparse validation observations with conservative pixel allowances.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "experimental" / "humbird_graves_2004"
CSV_PATH = OUTPUT_DIR / "seminar_surface_state_curves.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"

SEMINAR_PAGE = (
    "https://www.slideserve.com/moeshe/"
    "molecular-dynamics-simulations-of-plasma-surface-interactions-and-etching"
)


@dataclass(frozen=True)
class Source:
    key: str
    title: str
    url: str
    sha256: str
    width_px: int = 1023
    height_px: int = 777


@dataclass(frozen=True)
class Axis:
    key: str
    x_zero_px: float
    x_full_px: float
    x_full_value: float
    y_zero_px: float
    y_full_px: float
    y_full_value: float

    def x_value(self, pixel: float) -> float:
        return (
            (pixel - self.x_zero_px)
            * self.x_full_value
            / (self.x_full_px - self.x_zero_px)
        )

    def y_value(self, pixel: float) -> float:
        return (
            (self.y_zero_px - pixel)
            * self.y_full_value
            / (self.y_zero_px - self.y_full_px)
        )


@dataclass(frozen=True)
class Point:
    source: str
    panel: str
    axis: str
    boundary: str
    quantity: str
    trace: str
    x_px: int
    y_px: int
    uncertainty: float


SOURCES = {
    "slide20": Source(
        key="slide20",
        title="Thermal CF2 / Ar+: 9/1 (Si impacts)",
        url=(
            "https://image1.slideserve.com/3215203/"
            "thermal-cf-2-ar-9-1-si-impacts-l.jpg"
        ),
        sha256=(
            "20ee4c5c2b9b0f9e522343073d10cbca800a33806e694c93a11607b2985ace1c"
        ),
    ),
    "slide21": Source(
        key="slide21",
        title="Thermal F & CF2 / 200 eV Ar+",
        url=(
            "https://image1.slideserve.com/3215203/"
            "thermal-f-cf-2-200-ev-ar-l.jpg"
        ),
        sha256=(
            "ff8d8610a8785c37046c86b69702314725be72a784635a5e57cdd75b1a425bd5"
        ),
    ),
}


AXES = {
    "s20_20ev_state": Axis(
        "s20_20ev_state", 195.0, 486.0, 1500.0, 413.0, 212.0, 10.0),
    "s20_200ev_state": Axis(
        "s20_200ev_state", 195.0, 483.0, 1200.0, 705.0, 499.0, 20.0),
    "s20_200ev_yield": Axis(
        "s20_200ev_yield", 620.0, 916.0, 1200.0, 706.0, 483.0, 0.2),
    "s21_10f_state": Axis(
        "s21_10f_state", 200.0, 499.0, 1500.0, 411.0, 204.0, 50.0),
    "s21_10f_yield": Axis(
        "s21_10f_yield", 603.0, 911.0, 1500.0, 412.0, 191.0, 0.4),
    "s21_20f_state": Axis(
        "s21_20f_state", 200.0, 499.0, 1500.0, 689.0, 482.0, 80.0),
    "s21_20f_yield": Axis(
        "s21_20f_yield", 603.0, 911.0, 1500.0, 704.0, 451.0, 0.5),
}


def _points_for_trace(source, panel, axis, boundary, quantity, trace, pairs,
                      uncertainty):
    return [
        Point(
            source, panel, axis, boundary, quantity, trace, x, y, uncertainty)
        for x, y in pairs
    ]


# Centers were selected on the original 1023 x 777 rasters.  The colored
# traces use direct RGB-component isolation.  Black centers were selected
# after full-resolution inspection; sparse sampling avoids labels and axes.
POINTS = (
    _points_for_trace(
        "slide20", "20_eV_surface_state", "s20_20ev_state",
        "CF2/Ar+=9/1; normal 20 eV Ar+; Si", "surface_F_uptake_ML",
        "magenta",
        [(234, 318), (273, 288), (311, 283), (350, 269),
         (389, 253), (428, 247)],
        0.20,
    )
    + _points_for_trace(
        "slide20", "20_eV_surface_state", "s20_20ev_state",
        "CF2/Ar+=9/1; normal 20 eV Ar+; Si", "surface_C_uptake_ML",
        "blue",
        [(234, 366), (273, 348), (311, 339), (350, 327),
         (389, 317), (428, 310)],
        0.20,
    )
    + _points_for_trace(
        "slide20", "200_eV_surface_state", "s20_200ev_state",
        "CF2/Ar+=9/1; normal 200 eV Ar+; Si", "surface_F_uptake_ML",
        "magenta",
        [(243, 648), (291, 636), (339, 636), (387, 635), (435, 630)],
        0.40,
    )
    + _points_for_trace(
        "slide20", "200_eV_surface_state", "s20_200ev_state",
        "CF2/Ar+=9/1; normal 200 eV Ar+; Si", "surface_C_uptake_ML",
        "blue",
        [(243, 635), (291, 595), (339, 566), (387, 548), (435, 532)],
        0.40,
    )
    + _points_for_trace(
        "slide20", "200_eV_surface_state", "s20_200ev_state",
        "CF2/Ar+=9/1; normal 200 eV Ar+; Si", "cumulative_Si_etch_ML",
        "black",
        [(243, 649), (291, 609), (339, 582), (387, 564), (435, 552)],
        0.50,
    )
    + _points_for_trace(
        "slide20", "200_eV_etch_yield", "s20_200ev_yield",
        "CF2/Ar+=9/1; normal 200 eV Ar+; Si", "Si_etch_yield_per_ion",
        "black",
        [(669, 524), (719, 575), (768, 618), (817, 644), (867, 662)],
        0.006,
    )
    + _points_for_trace(
        "slide21", "10_percent_F_surface_state", "s21_10f_state",
        "CF2/F/Ar+=8/1/1; normal 200 eV Ar+; Si", "surface_F_uptake_ML",
        "magenta",
        [(240, 387), (280, 379), (320, 381), (359, 384),
         (399, 384), (439, 382), (479, 380)],
        0.75,
    )
    + _points_for_trace(
        "slide21", "10_percent_F_surface_state", "s21_10f_state",
        "CF2/F/Ar+=8/1/1; normal 200 eV Ar+; Si", "surface_C_uptake_ML",
        "blue",
        [(240, 390), (280, 376), (320, 372), (359, 367),
         (399, 365), (439, 362), (479, 362)],
        0.75,
    )
    + _points_for_trace(
        "slide21", "10_percent_F_surface_state", "s21_10f_state",
        "CF2/F/Ar+=8/1/1; normal 200 eV Ar+; Si", "cumulative_Si_etch_ML",
        "black",
        [(240, 365), (280, 332), (320, 307), (359, 285),
         (399, 263), (439, 245), (479, 230)],
        1.00,
    )
    + _points_for_trace(
        "slide21", "10_percent_F_etch_yield", "s21_10f_yield",
        "CF2/F/Ar+=8/1/1; normal 200 eV Ar+; Si",
        "Si_etch_yield_per_ion", "black",
        [(644, 260), (685, 281), (726, 317), (767, 334),
         (808, 336), (849, 347), (891, 353)],
        0.010,
    )
    + _points_for_trace(
        "slide21", "20_percent_F_surface_state", "s21_20f_state",
        "CF2/F/Ar+=7/2/1; normal 200 eV Ar+; Si", "surface_F_uptake_ML",
        "magenta",
        [(240, 674), (280, 673), (320, 671), (359, 670),
         (399, 672), (439, 673), (479, 674)],
        1.20,
    )
    + _points_for_trace(
        "slide21", "20_percent_F_surface_state", "s21_20f_state",
        "CF2/F/Ar+=7/2/1; normal 200 eV Ar+; Si", "surface_C_uptake_ML",
        "blue",
        # C and F overlap at several intermediate columns; retain only
        # locations where a blue center remains independently visible.
        [(240, 677), (359, 668), (399, 670), (439, 669), (479, 670)],
        1.20,
    )
    + _points_for_trace(
        "slide21", "20_percent_F_surface_state", "s21_20f_state",
        "CF2/F/Ar+=7/2/1; normal 200 eV Ar+; Si", "cumulative_Si_etch_ML",
        "black",
        [(240, 645), (280, 610), (320, 583), (359, 557),
         (399, 533), (439, 511), (479, 486)],
        1.50,
    )
    + _points_for_trace(
        "slide21", "20_percent_F_etch_yield", "s21_20f_yield",
        "CF2/F/Ar+=7/2/1; normal 200 eV Ar+; Si",
        "Si_etch_yield_per_ion", "black",
        [(644, 505), (685, 564), (726, 594), (767, 601),
         (808, 612), (849, 599), (891, 594)],
        0.015,
    )
)


CSV_FIELDS = [
    "source_slide",
    "source_panel",
    "boundary_condition",
    "quantity",
    "cf2_fluence_1e15_cm2",
    "digitized_value",
    "digitization_uncertainty",
    "marker_center_x_px",
    "marker_center_y_px",
    "trace_color",
    "evidence_grade",
    "source_image_sha256",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _trace_mask(rgb: np.ndarray, trace: str) -> np.ndarray:
    red = rgb[:, :, 0].astype(int)
    green = rgb[:, :, 1].astype(int)
    blue = rgb[:, :, 2].astype(int)
    if trace == "magenta":
        return (
            (red > 145)
            & (blue > 95)
            & (green < 155)
            & ((red - green) > 35)
            & ((blue - green) > 20)
        )
    if trace == "blue":
        return (
            (blue > 55)
            & (red < 130)
            & (green < 150)
            & ((blue - red) > 20)
        )
    if trace == "black":
        return np.max(rgb, axis=2) < 115
    raise ValueError(f"unknown trace color {trace}")


def _load_and_verify(paths: dict[str, Path]):
    images = {}
    for key, source in SOURCES.items():
        path = paths[key]
        actual_hash = _sha256(path)
        if actual_hash != source.sha256:
            raise RuntimeError(
                f"{key} SHA-256 mismatch: {actual_hash} != {source.sha256}")
        image = Image.open(path).convert("RGB")
        if image.size != (source.width_px, source.height_px):
            raise RuntimeError(
                f"{key} dimensions {image.size} do not match "
                f"{source.width_px} x {source.height_px}")
        images[key] = image
    return images


def _verify_points(images):
    for point in POINTS:
        rgb = np.asarray(images[point.source])
        mask = _trace_mask(rgb, point.trace)
        radius = 4
        x0 = max(0, point.x_px - radius)
        x1 = min(mask.shape[1], point.x_px + radius + 1)
        y0 = max(0, point.y_px - radius)
        y1 = min(mask.shape[0], point.y_px + radius + 1)
        if not np.any(mask[y0:y1, x0:x1]):
            raise RuntimeError(
                f"no {point.trace} trace pixel near "
                f"{point.source}:{point.x_px},{point.y_px}")


def _rows():
    rows = []
    for point in POINTS:
        axis = AXES[point.axis]
        rows.append({
            "source_slide": point.source,
            "source_panel": point.panel,
            "boundary_condition": point.boundary,
            "quantity": point.quantity,
            "cf2_fluence_1e15_cm2": f"{axis.x_value(point.x_px):.3f}",
            "digitized_value": f"{axis.y_value(point.y_px):.6f}",
            "digitization_uncertainty": f"{point.uncertainty:.6f}",
            "marker_center_x_px": str(point.x_px),
            "marker_center_y_px": str(point.y_px),
            "trace_color": point.trace,
            "evidence_grade": "primary_author_seminar_not_peer_reviewed",
            "source_image_sha256": SOURCES[point.source].sha256,
        })
    return rows


def _csv_text(rows) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _manifest():
    return {
        "schema_version": 1,
        "source": {
            "title": (
                "Molecular Dynamics Simulations of Plasma-Surface "
                "Interactions and Etching"
            ),
            "authors": ["David Humbird", "David Graves"],
            "affiliation": "University of California, Berkeley",
            "event": "FLCC Research Seminar",
            "date": "2004-02-23",
            "landing_page": SEMINAR_PAGE,
            "evidence_grade": "primary_author_seminar_not_peer_reviewed",
            "relationship_to_peer_reviewed_sources": [
                "doi:10.1063/1.1644338",
                "doi:10.1063/1.1736321",
                "doi:10.1063/1.1769602",
            ],
        },
        "rasters": {
            key: {
                "title": item.title,
                "url": item.url,
                "sha256": item.sha256,
                "width_px": item.width_px,
                "height_px": item.height_px,
            }
            for key, item in SOURCES.items()
        },
        "digitization": {
            "method": (
                "PIL/NumPy full-resolution RGB-component isolation; "
                "sparse manually reconciled trace centers; explicit linear "
                "axis transforms; source pixels verified within 4 px"
            ),
            "curve_fit_used": False,
            "model_values_visible_or_used": False,
            "uncertainty_semantics": (
                "conservative digitization allowance only; the seminar does "
                "not publish MD ensemble or potential uncertainty"
            ),
            "copyrighted_rasters_redistributed": False,
        },
        "axis_calibrations": {
            key: {
                "x_zero_px": axis.x_zero_px,
                "x_full_px": axis.x_full_px,
                "x_full_value": axis.x_full_value,
                "y_zero_px": axis.y_zero_px,
                "y_full_px": axis.y_full_px,
                "y_full_value": axis.y_full_value,
            }
            for key, axis in AXES.items()
        },
        "observation_count": len(POINTS),
        "limitations": [
            "The plotted calculations are classical molecular dynamics, not exact quantum dynamics.",
            "The seminar raster is primary-author evidence but is not the archival journal figure.",
            "The source does not publish statistical error bars for these traces.",
            "Monolayer normalization and the empirical interatomic potential are inherited from the source.",
            "The curves constrain surface-state topology and transient response; they are not reactor flux measurements.",
        ],
    }


def _json_text(value) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def _check_or_write(path: Path, expected: str, check: bool):
    if check:
        if not path.exists():
            raise RuntimeError(f"missing committed artifact: {path}")
        actual = path.read_text()
        if actual != expected:
            raise RuntimeError(f"committed artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(expected)


def _draw_overlays(images, destination: Path):
    destination.mkdir(parents=True, exist_ok=True)
    colors = {"magenta": (255, 0, 255), "blue": (0, 180, 255), "black": (255, 0, 0)}
    for key, source_image in images.items():
        image = source_image.copy()
        draw = ImageDraw.Draw(image)
        for point in POINTS:
            if point.source != key:
                continue
            radius = 5
            color = colors[point.trace]
            draw.ellipse(
                (point.x_px - radius, point.y_px - radius,
                 point.x_px + radius, point.y_px + radius),
                outline=color,
                width=2,
            )
        image.save(destination / f"{key}_digitization_overlay.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slide-20", type=Path, required=True)
    parser.add_argument("--slide-21", type=Path, required=True)
    parser.add_argument("--overlay-dir", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    images = _load_and_verify({
        "slide20": args.slide_20,
        "slide21": args.slide_21,
    })
    _verify_points(images)
    _check_or_write(CSV_PATH, _csv_text(_rows()), args.check)
    _check_or_write(MANIFEST_PATH, _json_text(_manifest()), args.check)
    if args.overlay_dir is not None:
        _draw_overlays(images, args.overlay_dir)


if __name__ == "__main__":
    main()
