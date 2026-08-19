#!/usr/bin/env python3
"""Reproduce Benck Figure 14(a)'s 5 sccm CF2/CF feed board.

The AIP source pixels are not redistributed.  The committed table is rebuilt
from full-page coordinates on a checksum-pinned 600 dpi Poppler render.  The
5 sccm filled-circle series is retained because it aligns with the published
5 sccm ion feed board used by the light-ion inverse; the separate 10 sccm
open-circle sensitivity is not silently mixed into that condition.
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
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "benck_2003_c4f6"
CSV_PATH = OUTPUT_DIRECTORY / "figure14a_cf2_cf_feed_ratio.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure14a_digitization_manifest.json"
SOURCE_PDF_SHA256 = (
    "96ef064afb0f804d3d853adcc275a63cd52507b0ea7e39f921d6cf6209e86d08"
)
RENDER_SHA256 = (
    "cba82d310e6da875d397839fb489405900338fb2f691bce6227f49d2805bddee"
)
RENDER_SIZE = (5192, 6745)
RENDER_DPI = 600

# Full-page major ticks on both sides of the slightly rotated scan.  A 2D
# affine calibration removes the page tilt instead of treating it as signal.
LEFT_X_PX = 925.5
RIGHT_X_PX = 2429.5
LEFT_Y_TICKS_PX = (1729.5, 1545.5, 1361.5, 1177.5, 993.5, 813.5)
RIGHT_Y_TICKS_PX = (1719.5, 1537.5, 1353.5, 1169.5, 985.5, 805.5)
CF2_CF_AT_TICKS = (0.0, 4.0, 8.0, 12.0, 16.0, 20.0)


@dataclass(frozen=True)
class PixelPoint:
    c4f6_percent: int
    x_px: float
    y_px: float


PIXEL_POINTS = (
    PixelPoint(25, 1266.8, 1228.1),
    PixelPoint(50, 1609.3, 1095.2),
    PixelPoint(75, 1950.4, 1020.5),
    PixelPoint(100, 2289.2, 942.0),
)
FIELDNAMES = (
    "c4f6_feed_percent",
    "ar_feed_percent",
    "flow_sccm",
    "pressure_Pa",
    "power_W",
    "neutral_CF2_to_CF_density_ratio",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "marker",
    "digitization_absolute_bound",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _axis_coefficients() -> tuple[float, float, float]:
    x = np.asarray(
        [LEFT_X_PX] * 6 + [RIGHT_X_PX] * 6,
        dtype=float,
    )
    y = np.asarray(LEFT_Y_TICKS_PX + RIGHT_Y_TICKS_PX, dtype=float)
    value = np.asarray(CF2_CF_AT_TICKS * 2, dtype=float)
    design = np.column_stack((np.ones_like(x), y, x))
    intercept, y_slope, x_slope = np.linalg.lstsq(
        design, value, rcond=None
    )[0]
    return float(intercept), float(y_slope), float(x_slope)


def ratio_at_pixel(x_px: float, y_px: float) -> float:
    intercept, y_slope, x_slope = _axis_coefficients()
    return intercept + y_slope * float(y_px) + x_slope * float(x_px)


def rows() -> list[dict[str, str]]:
    return [
        {
            "c4f6_feed_percent": str(point.c4f6_percent),
            "ar_feed_percent": str(100 - point.c4f6_percent),
            "flow_sccm": "5.0",
            "pressure_Pa": "1.33",
            "power_W": "200.0",
            "neutral_CF2_to_CF_density_ratio": (
                f"{ratio_at_pixel(point.x_px, point.y_px):.9g}"
            ),
            "marker_center_x_full_page_px": f"{point.x_px:.1f}",
            "marker_center_y_full_page_px": f"{point.y_px:.1f}",
            "marker": "filled_circle",
            # 12 pixels is 0.261 ratio unit.  This bounds marker-center
            # transcription but is kept separate from the plotted fit bars,
            # which the paper says are not true statistical uncertainty.
            "digitization_absolute_bound": "0.27",
        }
        for point in PIXEL_POINTS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_digest: str) -> dict[str, object]:
    intercept, y_slope, x_slope = _axis_coefficients()
    return {
        "manifest_id": "BENCK-2003-C4F6-FIG14A-CF2-CF-R1",
        "source": {
            "citation": (
                "E. C. Benck, A. Goyette, and Y. Wang, Journal of "
                "Applied Physics 94, 1382-1389 (2003)"
            ),
            "doi": "10.1063/1.1586978",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 8,
            "print_page": 1389,
            "figure": "Figure 14(a)",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "redistribution": (
                "source pixels are not committed; the scan carries an AIP "
                "redistribution notice"
            ),
        },
        "experiment": {
            "reactor": "modified inductively coupled GEC reference cell",
            "pressure_Pa": 1.33,
            "power_W": 200.0,
            "flow_sccm": 5.0,
            "series": "filled circles",
            "observable": "line-of-sight neutral CF2/CF density ratio",
            "paper_caution": (
                "the diagnostic overweights cooler outside-plasma path and "
                "surface-produced CF2; the in-plasma ratio is probably lower"
            ),
        },
        "pixel_calibration": {
            "full_page_left_x_px": LEFT_X_PX,
            "full_page_right_x_px": RIGHT_X_PX,
            "full_page_left_y_major_tick_px": list(LEFT_Y_TICKS_PX),
            "full_page_right_y_major_tick_px": list(RIGHT_Y_TICKS_PX),
            "CF2_CF_at_major_ticks": list(CF2_CF_AT_TICKS),
            "ratio_intercept": intercept,
            "ratio_per_vertical_pixel": y_slope,
            "ratio_per_horizontal_pixel_scan_tilt": x_slope,
            "left_major_tick_replay": [
                ratio_at_pixel(LEFT_X_PX, y) for y in LEFT_Y_TICKS_PX
            ],
            "right_major_tick_replay": [
                ratio_at_pixel(RIGHT_X_PX, y) for y in RIGHT_Y_TICKS_PX
            ],
            "axis_fit": (
                "least-squares 2D affine fit to twelve linear-axis ticks; "
                "horizontal term removes scan rotation"
            ),
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy dark-pixel axis and "
                "filled-marker localization; original-resolution visual audit"
            ),
            "retained_series": "5 sccm filled circles",
            "excluded_series": (
                "10 sccm open circles retained as qualitative flow "
                "sensitivity and not mixed with the 5 sccm ion board"
            ),
            "vertical_pixel_allowance": 12.0,
            "digitization_absolute_bound": 0.27,
            "plotted_error_bars_are_statistical_uncertainty": False,
            "visual_audit_status": "passed_original_resolution",
        },
        "derived_checks": {
            "point_count": len(PIXEL_POINTS),
            "ratio_increases_from_25_to_100_percent": True,
            "source_text_approximate_endpoint_range": "10 to 16",
        },
        "claim_boundary": {
            "valid": (
                "condition-resolved neutral-ratio constraint for the same "
                "Benck 5 sccm feed board"
            ),
            "not_valid": [
                "an absolute CF or CF2 density",
                "a stable C4F6 parent density",
                "a volume-local in-plasma ratio without line-of-sight bias",
                "a Krueger boundary or feature-depth coefficient",
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
        raise RuntimeError("committed Benck Figure-14a CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Benck Figure-14a manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Benck source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Benck 600-dpi page-8 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Benck render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    for point in PIXEL_POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark = int(np.sum(gray[y - 30:y + 31, x - 30:x + 31] < 96))
        if dark < 1200:
            raise RuntimeError(
                f"insufficient filled-marker support at "
                f"{point.c4f6_percent}%: {dark}"
            )
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for point in PIXEL_POINTS:
        radius = 34
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline="#e41a1c",
            width=6,
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
