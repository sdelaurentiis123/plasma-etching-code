#!/usr/bin/env python3
"""Reproduce Ji et al. Figure-3(b1-b3) TiO2 morphology responses.

The source is CC BY, but source pixels are not committed here. The numerical
table replays full-page marker centers on a checksum-pinned 600 dpi Poppler
render. It is a same-gas mechanism-response board, not an Oxford coefficient.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = (
    ROOT / "data" / "experimental" / "ji_2024_tio2_hierarchical"
)
CSV_PATH = OUTPUT_DIRECTORY / "figure3_rf_morphology_response.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure3_digitization_manifest.json"
SOURCE_PDF_SHA256 = (
    "92a57c600e93e4113e574ef286aae4a248d787c3445f8d21e2ff151c73647e81"
)
RENDER_SHA256 = (
    "976db637feed609a182301266654e675a42d1422f6e7c08cff17ea3744d6c148"
)
RENDER_SIZE = (4961, 7016)
RENDER_DPI = 600


@dataclass(frozen=True)
class Series:
    panel: str
    observable: str
    unit: str
    y_ticks_px: tuple[float, ...]
    y_tick_values: tuple[float, ...]
    marker_centers_px: tuple[tuple[float, float], ...]
    digitization_bound: float


RF_POWER_W = (90.0, 120.0, 150.0, 180.0, 210.0)
SERIES = (
    Series(
        "3b1",
        "upper_triangle_height",
        "nm",
        (1961.5, 1863.5, 1765.5, 1667.5, 1569.0),
        (150.0, 200.0, 250.0, 300.0, 350.0),
        (
            (3944.0, 1974.0),
            (4074.5, 1850.0),
            (4205.0, 1726.5),
            (4335.5, 1573.5),
            (4466.0, 1596.5),
        ),
        2.0,
    ),
    Series(
        "3b2",
        "upper_tip_corner_radius",
        "nm",
        (2537.5, 2448.5, 2359.0, 2270.0, 2180.5),
        (20.0, 40.0, 60.0, 80.0, 100.0),
        (
            (3939.5, 2181.0),
            (4070.0, 2248.0),
            (4200.0, 2560.0),
            (4330.5, 2583.0),
            (4461.0, 2583.0),
        ),
        1.0,
    ),
    Series(
        "3b3",
        "interfeature_gap",
        "nm",
        (3149.0, 3060.5, 2971.0, 2882.0, 2793.0),
        (20.0, 40.0, 60.0, 80.0, 100.0),
        (
            (3938.0, 2811.0),
            (4069.5, 2882.0),
            (4201.5, 2927.0),
            (4333.0, 3051.0),
            (4465.0, 3158.0),
        ),
        1.0,
    ),
)
FIELDNAMES = (
    "panel",
    "rf_power_W",
    "observable",
    "value",
    "unit",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "digitization_absolute_bound",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _axis_coefficients(series: Series) -> tuple[float, float]:
    design = np.column_stack((
        np.ones(len(series.y_ticks_px)),
        np.asarray(series.y_ticks_px, dtype=float),
    ))
    intercept, slope = np.linalg.lstsq(
        design,
        np.asarray(series.y_tick_values, dtype=float),
        rcond=None,
    )[0]
    return float(intercept), float(slope)


def value_at_y(series: Series, y_px: float) -> float:
    intercept, slope = _axis_coefficients(series)
    return intercept + slope * float(y_px)


def rows() -> list[dict[str, str]]:
    output = []
    for series in SERIES:
        for power, (x_px, y_px) in zip(RF_POWER_W, series.marker_centers_px):
            output.append({
                "panel": series.panel,
                "rf_power_W": f"{power:.1f}",
                "observable": series.observable,
                "value": f"{value_at_y(series, y_px):.9g}",
                "unit": series.unit,
                "marker_center_x_full_page_px": f"{x_px:.1f}",
                "marker_center_y_full_page_px": f"{y_px:.1f}",
                "digitization_absolute_bound": (
                    f"{series.digitization_bound:.1f}"
                ),
            })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _values(observable: str) -> list[float]:
    return [
        float(row["value"])
        for row in rows()
        if row["observable"] == observable
    ]


def manifest(csv_digest: str) -> dict[str, object]:
    calibrations = {}
    for series in SERIES:
        intercept, slope = _axis_coefficients(series)
        calibrations[series.panel] = {
            "observable": series.observable,
            "intercept": intercept,
            "value_per_vertical_pixel": slope,
            "tick_replay": [
                value_at_y(series, value) for value in series.y_ticks_px
            ],
            "tick_values": list(series.y_tick_values),
        }
    height = _values("upper_triangle_height")
    radius = _values("upper_tip_corner_radius")
    gap = _values("interfeature_gap")
    return {
        "manifest_id": "JI-2024-TIO2-FIG3-RF-MORPHOLOGY-R1",
        "source": {
            "citation": (
                "X. Ji et al., Micromachines 15, 1160 (2024)"
            ),
            "doi": "10.3390/mi15091160",
            "license": "CC BY 4.0",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 6,
            "figure": "Figure 3(b1-b3)",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "source_pixels_committed": False,
        },
        "experiment": {
            "material": "800 nm electron-beam-deposited TiO2",
            "mask": "60 nm Cr",
            "geometry": "200 nm nominal CD, 300 nm pitch grating",
            "source_power_W": 350.0,
            "rf_power_W": list(RF_POWER_W),
            "feed_sccm": {"SF6": 40.0, "CHF3": 10.0, "O2": 5.0},
            "pressure_mTorr": 10.0,
            "temperature_C": 40.0,
            "etch_time_s": None,
        },
        "pixel_calibration": calibrations,
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy red-marker localization; "
                "linear ordinate calibration; original-resolution visual audit"
            ),
            "point_count": len(rows()),
            "height_digitization_absolute_bound_nm": 2.0,
            "radius_and_gap_digitization_absolute_bound_nm": 1.0,
            "source_measurement_uncertainty_reported": False,
            "visual_audit_status": "passed_original_resolution",
        },
        "derived_checks": {
            "upper_height_peaks_at_180W": (
                int(np.argmax(height)) == RF_POWER_W.index(180.0)
            ),
            "corner_radius_nonincreasing": all(
                left >= right for left, right in zip(radius, radius[1:])
            ),
            "gap_strictly_decreases": all(
                left > right for left, right in zip(gap, gap[1:])
            ),
            "gap_change_90W_to_210W_nm": gap[-1] - gap[0],
        },
        "physics_use": {
            "valid": (
                "same-gas TiO2/Cr response validation for a surface mechanism "
                "with mask erosion and passivation/growth"
            ),
            "removal_only_model_can_reproduce_strict_gap_narrowing": False,
            "required_model_topology": [
                "evolving Cr mask",
                "fluorinated or passivating surface inventory",
                "positive deposited or retained surface-volume channel",
                "ion-assisted removal of TiO2 and passivating material",
            ],
            "not_valid": [
                "an Oxford NPG80 TiO2 coefficient",
                "a Freddie target profile",
                "an etch-time or absolute-depth calibration",
            ],
        },
        "output": {
            "path": str(CSV_PATH.relative_to(ROOT)),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Ji Figure-3 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Ji Figure-3 manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Ji source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Ji 600-dpi page-6 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Ji render size: {image.size}")
    rgb = np.asarray(image)
    for series in SERIES:
        for x_px, y_px in series.marker_centers_px:
            x = int(round(x_px))
            y = int(round(y_px))
            patch = rgb[y - 25:y + 26, x - 25:x + 26]
            red = (
                (patch[:, :, 0] > 180)
                & (patch[:, :, 1] < 120)
                & (patch[:, :, 2] < 120)
            )
            if int(red.sum()) < 250:
                raise RuntimeError(
                    f"insufficient red-marker support at {series.panel}: "
                    f"{int(red.sum())}"
                )
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {"3b1": "#377eb8", "3b2": "#4daf4a", "3b3": "#984ea3"}
    for series in SERIES:
        for x_px, y_px in series.marker_centers_px:
            radius = 28
            draw.ellipse(
                (
                    x_px - radius,
                    y_px - radius,
                    x_px + radius,
                    y_px + radius,
                ),
                outline=colors[series.panel],
                width=7,
            )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render", type=Path)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.source_pdf:
        verify_source_pdf(args.source_pdf)
    image = verify_render(args.render) if args.render else None
    if args.overlay:
        if image is None:
            parser.error("--overlay requires --render")
        draw_overlay(image, args.overlay)
    payload = csv_text()
    if args.write:
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(payload, encoding="utf-8")
        MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")
    if args.check:
        verify_committed_files()
    if not any((args.write, args.check, args.source_pdf, args.render)):
        parser.error("select --write, --check, --source-pdf, or --render")


if __name__ == "__main__":
    main()
