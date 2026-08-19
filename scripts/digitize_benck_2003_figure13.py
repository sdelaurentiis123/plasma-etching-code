#!/usr/bin/env python3
"""Reproduce Benck Figure 13 plasma-potential oscillation boards.

The AIP source pixels are not redistributed. The committed table is rebuilt
from full-page coordinates on a checksum-pinned 600 dpi Poppler render. Both
panels are retained: feed/flow response at 1.33 Pa and pressure response for
pure and 50% C4F6 at 5 sccm.
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
CSV_PATH = OUTPUT_DIRECTORY / "figure13_plasma_potential_oscillation.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure13_digitization_manifest.json"
SOURCE_PDF_SHA256 = (
    "96ef064afb0f804d3d853adcc275a63cd52507b0ea7e39f921d6cf6209e86d08"
)
RENDER_SHA256 = (
    "77d5990061addab5a38fc24dea4fd501b11dfaaac6219c1400140dfe23c62e51"
)
RENDER_SIZE = (5181, 6745)
RENDER_DPI = 600

# The scanned page is slightly rotated. Each panel therefore uses a 2D affine
# calibration through seven major ticks on both vertical axes. Coordinates
# are full-page pixels; values are 0, 4, ..., 24 V from bottom to top.
PANEL_A_LEFT_X_PX = 2901.0
PANEL_A_RIGHT_X_PX = 4495.0
PANEL_A_LEFT_Y_TICKS_PX = (
    1716.4, 1551.5, 1386.6, 1223.1, 1058.7, 896.8, 727.7,
)
PANEL_A_RIGHT_Y_TICKS_PX = (
    1698.0, 1534.0, 1370.0, 1206.0, 1042.0, 878.0, 714.0,
)
PANEL_B_LEFT_X_PX = 2901.0
PANEL_B_RIGHT_X_PX = 4495.0
PANEL_B_LEFT_Y_TICKS_PX = (
    3077.2, 2914.6, 2751.7, 2587.6, 2423.0, 2261.0, 2094.7,
)
PANEL_B_RIGHT_Y_TICKS_PX = (
    3062.0, 2902.0, 2738.0, 2574.0, 2410.0, 2246.0, 2082.0,
)
VPP_AT_TICKS = (0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0)


@dataclass(frozen=True)
class PixelPoint:
    panel: str
    series: str
    c4f6_feed_percent: int
    flow_sccm: float
    pressure_Pa: float
    x_px: float
    y_px: float
    marker: str
    overlap_note: str = ""


PIXEL_POINTS = (
    PixelPoint("13a", "5_sccm", 25, 5.0, 1.33, 3303.3, 1438.5, "filled_circle"),
    PixelPoint(
        "13a", "10_sccm", 25, 10.0, 1.33, 3303.3, 1400.0,
        "open_circle", "partially overlaps the 5 sccm marker",
    ),
    PixelPoint("13a", "5_sccm", 50, 5.0, 1.33, 3651.6, 1262.1, "filled_circle"),
    PixelPoint("13a", "10_sccm", 50, 10.0, 1.33, 3651.8, 1142.7, "open_circle"),
    PixelPoint("13a", "5_sccm", 75, 5.0, 1.33, 4004.3, 1149.0, "filled_circle"),
    PixelPoint("13a", "5_sccm", 100, 5.0, 1.33, 4358.9, 886.7, "filled_circle"),
    PixelPoint("13b", "100_percent_C4F6", 100, 5.0, 0.67, 3291.5, 2472.3, "filled_circle"),
    PixelPoint("13b", "50_percent_C4F6", 50, 5.0, 0.67, 3292.6, 2945.7, "open_circle"),
    PixelPoint("13b", "100_percent_C4F6", 100, 5.0, 1.33, 3634.0, 2253.6, "filled_circle"),
    PixelPoint("13b", "50_percent_C4F6", 50, 5.0, 1.33, 3632.9, 2627.8, "open_circle"),
    PixelPoint(
        "13b", "100_percent_C4F6", 100, 5.0, 2.66, 4322.7, 2352.0,
        "filled_circle", "partially overlaps the 50 percent marker",
    ),
    PixelPoint(
        "13b", "50_percent_C4F6", 50, 5.0, 2.66, 4322.7, 2398.0,
        "open_circle", "partially overlaps the 100 percent marker",
    ),
)
FIELDNAMES = (
    "panel",
    "series",
    "c4f6_feed_percent",
    "ar_feed_percent",
    "flow_sccm",
    "pressure_Pa",
    "power_W",
    "plasma_potential_peak_to_peak_V",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "marker",
    "digitization_absolute_bound_V",
    "overlap_note",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _axis_coefficients(panel: str) -> tuple[float, float, float]:
    if panel == "13a":
        left_x = PANEL_A_LEFT_X_PX
        right_x = PANEL_A_RIGHT_X_PX
        left_y = PANEL_A_LEFT_Y_TICKS_PX
        right_y = PANEL_A_RIGHT_Y_TICKS_PX
    elif panel == "13b":
        left_x = PANEL_B_LEFT_X_PX
        right_x = PANEL_B_RIGHT_X_PX
        left_y = PANEL_B_LEFT_Y_TICKS_PX
        right_y = PANEL_B_RIGHT_Y_TICKS_PX
    else:
        raise ValueError(f"unknown Figure-13 panel {panel!r}")
    x = np.asarray([left_x] * 7 + [right_x] * 7, dtype=float)
    y = np.asarray(left_y + right_y, dtype=float)
    value = np.asarray(VPP_AT_TICKS * 2, dtype=float)
    design = np.column_stack((np.ones_like(x), y, x))
    intercept, y_slope, x_slope = np.linalg.lstsq(
        design, value, rcond=None
    )[0]
    return float(intercept), float(y_slope), float(x_slope)


def vpp_at_pixel(panel: str, x_px: float, y_px: float) -> float:
    intercept, y_slope, x_slope = _axis_coefficients(panel)
    return intercept + y_slope * float(y_px) + x_slope * float(x_px)


def rows() -> list[dict[str, str]]:
    output = []
    for point in PIXEL_POINTS:
        output.append({
            "panel": point.panel,
            "series": point.series,
            "c4f6_feed_percent": str(point.c4f6_feed_percent),
            "ar_feed_percent": str(100 - point.c4f6_feed_percent),
            "flow_sccm": f"{point.flow_sccm:.1f}",
            "pressure_Pa": f"{point.pressure_Pa:.2f}",
            "power_W": "200.0",
            "plasma_potential_peak_to_peak_V": (
                f"{vpp_at_pixel(point.panel, point.x_px, point.y_px):.9g}"
            ),
            "marker_center_x_full_page_px": f"{point.x_px:.1f}",
            "marker_center_y_full_page_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            # Twenty pixels is about 0.49 V. The rounded 0.6 V allowance also
            # covers manual separation of the two overlapping marker pairs.
            "digitization_absolute_bound_V": "0.6",
            "overlap_note": point.overlap_note,
        })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _tick_replay(panel: str, x: float, y_ticks) -> list[float]:
    return [vpp_at_pixel(panel, x, y) for y in y_ticks]


def manifest(csv_digest: str) -> dict[str, object]:
    a_intercept, a_y_slope, a_x_slope = _axis_coefficients("13a")
    b_intercept, b_y_slope, b_x_slope = _axis_coefficients("13b")
    values = [float(row["plasma_potential_peak_to_peak_V"]) for row in rows()]
    return {
        "manifest_id": "BENCK-2003-C4F6-FIG13-VPP-R1",
        "source": {
            "citation": (
                "E. C. Benck, A. Goyette, and Y. Wang, Journal of "
                "Applied Physics 94, 1382-1389 (2003)"
            ),
            "doi": "10.1063/1.1586978",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 7,
            "print_page": 1388,
            "figure": "Figure 13(a,b)",
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
            "power_W": 200.0,
            "observable": (
                "plasma-potential peak-to-peak oscillation inferred from the "
                "effective widths of measured ion-energy distributions"
            ),
            "panel_13a": (
                "C4F6 feed fraction and 5/10 sccm flow at 1.33 Pa"
            ),
            "panel_13b": (
                "0.67/1.33/2.66 Pa pressure at 5 sccm for pure and 50% C4F6"
            ),
        },
        "pixel_calibration": {
            "axis_fit": (
                "separate least-squares 2D affine fits to fourteen linear-axis "
                "ticks per panel; horizontal term removes scan rotation"
            ),
            "panel_13a": {
                "intercept": a_intercept,
                "volts_per_vertical_pixel": a_y_slope,
                "volts_per_horizontal_pixel_scan_tilt": a_x_slope,
                "left_tick_replay_V": _tick_replay(
                    "13a", PANEL_A_LEFT_X_PX, PANEL_A_LEFT_Y_TICKS_PX
                ),
                "right_tick_replay_V": _tick_replay(
                    "13a", PANEL_A_RIGHT_X_PX, PANEL_A_RIGHT_Y_TICKS_PX
                ),
            },
            "panel_13b": {
                "intercept": b_intercept,
                "volts_per_vertical_pixel": b_y_slope,
                "volts_per_horizontal_pixel_scan_tilt": b_x_slope,
                "left_tick_replay_V": _tick_replay(
                    "13b", PANEL_B_LEFT_X_PX, PANEL_B_LEFT_Y_TICKS_PX
                ),
                "right_tick_replay_V": _tick_replay(
                    "13b", PANEL_B_RIGHT_X_PX, PANEL_B_RIGHT_Y_TICKS_PX
                ),
            },
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy dark-pixel axis and marker "
                "localization; original-resolution visual audit"
            ),
            "point_count": len(PIXEL_POINTS),
            "overlapping_marker_pairs": 2,
            "digitization_absolute_bound_V": 0.6,
            "visual_audit_status": "passed_original_resolution",
        },
        "derived_checks": {
            "minimum_V": min(values),
            "maximum_V": max(values),
            "five_sccm_feed_series_increases": all(
                earlier < later
                for earlier, later in zip(
                    [values[index] for index in (0, 2, 4)],
                    [values[index] for index in (2, 4, 5)],
                )
            ),
            "pure_C4F6_pressure_response_is_nonmonotonic": True,
            "fifty_percent_C4F6_pressure_response_increases": True,
        },
        "claim_boundary": {
            "valid": (
                "independent voltage-response validation board for a C4F6/Ar "
                "reactor and sheath model at the Benck apparatus"
            ),
            "not_valid": [
                "a direct voltage waveform measurement",
                "a target-machine self-bias",
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
        raise RuntimeError("committed Benck Figure-13 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Benck Figure-13 manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Benck source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Benck 600-dpi page-7 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Benck render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    for point in PIXEL_POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark = int(np.sum(gray[y - 45:y + 46, x - 45:x + 46] < 190))
        if dark < 180:
            raise RuntimeError(
                f"insufficient marker support at {point.panel}/{point.series}: {dark}"
            )
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for point in PIXEL_POINTS:
        radius = 42
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline="#e41a1c" if point.marker == "filled_circle" else "#377eb8",
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
