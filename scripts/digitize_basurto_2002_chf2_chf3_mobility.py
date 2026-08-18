#!/usr/bin/env python3
"""Reproduce Basurto 2002 Figure-1 CHF2+ in CHF3 digitization."""
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
OUTPUT_DIR = ROOT / "data" / "experimental" / "basurto_2002_chf3"
CSV_PATH = OUTPUT_DIR / "figure1_chf2_chf3_reduced_mobility.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "a959267e28687a1497ddb96669bdc096d096446a3c08765e1a5a5c3ecdd53b48"
)
RENDER_SHA256 = (
    "44f37a67b73d60bdcf5a7dbc510d7ce6b656f6b77ea7b83f0d11972014f2dc74"
)
RENDER_SIZE = (5100, 6600)
RENDER_DPI = 600

# Coordinates below use a (2200, 200) full-page crop origin. Both axes are
# logarithmic. Tick centers were localized from the bilateral plot ticks and
# thick border centerlines; the three labeled major ticks overdetermine each
# calibration.
CROP_ORIGIN = (2200.0, 200.0)
X_TICKS_CROP_PX = np.array([951.0, 1394.5, 2028.0])
REDUCED_FIELD_AT_X_TICKS_TD = np.array([20.0, 100.0, 1000.0])
Y_TICKS_CROP_PX = np.array([246.0, 733.0, 1010.0])
REDUCED_MOBILITY_AT_Y_TICKS_CM2_V_S = np.array([5.0, 1.0, 0.4])


@dataclass(frozen=True)
class PixelPoint:
    x_crop_px: float
    y_crop_px: float

    @property
    def x_full_page_px(self) -> float:
        return self.x_crop_px + CROP_ORIGIN[0]

    @property
    def y_full_page_px(self) -> float:
        return self.y_crop_px + CROP_ORIGIN[1]


# Open-circle centers for CHF2+ in CHF3 only. Hough candidates were followed
# by an original-resolution visual overlay audit. The overlapping symbols and
# one genuine low mobility point near 230 Td are preserved, not smoothed.
PIXEL_POINTS = tuple(PixelPoint(*point) for point in (
    (1175.5, 922.5), (1202.5, 926.5), (1229.5, 920.5),
    (1252.5, 921.5), (1275.5, 939.5), (1296.5, 936.5),
    (1316.5, 928.5), (1333.5, 942.5), (1364.5, 941.5),
    (1394.5, 950.5), (1419.5, 953.5), (1445.5, 954.5),
    (1466.5, 937.5), (1505.5, 939.5), (1555.5, 920.5),
    (1585.5, 914.5), (1622.5, 936.5), (1647.5, 895.5),
    (1666.5, 887.5), (1696.5, 850.5), (1739.5, 826.5),
    (1775.5, 785.5), (1808.5, 776.5),
))

FIELDNAMES = (
    "reduced_field_Td",
    "reduced_mobility_cm2_V_s",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "digitization_reduced_field_relative_bound",
    "digitization_reduced_mobility_relative_bound",
    "source_measurement_relative_uncertainty_lower",
    "source_measurement_relative_uncertainty_upper",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _log_axis_coefficients(
    pixels: np.ndarray, values: np.ndarray,
) -> tuple[float, float]:
    matrix = np.column_stack((np.ones(pixels.size), np.log10(values)))
    intercept, slope = np.linalg.lstsq(matrix, pixels, rcond=None)[0]
    return float(intercept), float(slope)


def reduced_field_at_crop_pixel(x_px: float) -> float:
    intercept, slope = _log_axis_coefficients(
        X_TICKS_CROP_PX, REDUCED_FIELD_AT_X_TICKS_TD)
    return 10.0 ** ((float(x_px) - intercept) / slope)


def reduced_mobility_at_crop_pixel(y_px: float) -> float:
    intercept, slope = _log_axis_coefficients(
        Y_TICKS_CROP_PX, REDUCED_MOBILITY_AT_Y_TICKS_CM2_V_S)
    return 10.0 ** ((float(y_px) - intercept) / slope)


def rows() -> list[dict[str, str]]:
    _, x_slope = _log_axis_coefficients(
        X_TICKS_CROP_PX, REDUCED_FIELD_AT_X_TICKS_TD)
    _, y_slope = _log_axis_coefficients(
        Y_TICKS_CROP_PX, REDUCED_MOBILITY_AT_Y_TICKS_CM2_V_S)
    x_bound = 10.0 ** (10.0 / abs(x_slope)) - 1.0
    y_bound = 10.0 ** (10.0 / abs(y_slope)) - 1.0
    return [{
        "reduced_field_Td": f"{reduced_field_at_crop_pixel(p.x_crop_px):.8g}",
        "reduced_mobility_cm2_V_s": (
            f"{reduced_mobility_at_crop_pixel(p.y_crop_px):.8g}"),
        "marker_center_x_full_page_px": f"{p.x_full_page_px:.1f}",
        "marker_center_y_full_page_px": f"{p.y_full_page_px:.1f}",
        "digitization_reduced_field_relative_bound": f"{x_bound:.8g}",
        "digitization_reduced_mobility_relative_bound": f"{y_bound:.8g}",
        "source_measurement_relative_uncertainty_lower": "0.02",
        "source_measurement_relative_uncertainty_upper": "0.04",
    } for p in PIXEL_POINTS]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_digest: str) -> dict:
    x_intercept, x_slope = _log_axis_coefficients(
        X_TICKS_CROP_PX, REDUCED_FIELD_AT_X_TICKS_TD)
    y_intercept, y_slope = _log_axis_coefficients(
        Y_TICKS_CROP_PX, REDUCED_MOBILITY_AT_Y_TICKS_CM2_V_S)
    return {
        "manifest_id": "BASURTO-2002-CHF2-CHF3-FIG1-R1",
        "source": {
            "citation": (
                "E. Basurto and J. de Urquijo, Journal of Applied Physics "
                "91, 36-39 (2002)"
            ),
            "title": "Mobility of CF3+ in CF4, CHF2+ in CHF3, and C+ in Ar",
            "doi": "10.1063/1.1421034",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 2,
            "print_page": 37,
            "figure": "Figure 1, CHF2+ in CHF3 open-circle series",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "redistribution": "source pixels are not committed",
        },
        "pixel_calibration": {
            "crop_origin_full_page_px": list(CROP_ORIGIN),
            "x_tick_centers_crop_px": X_TICKS_CROP_PX.tolist(),
            "reduced_field_at_x_ticks_Td": (
                REDUCED_FIELD_AT_X_TICKS_TD.tolist()),
            "x_axis": "log10",
            "x_pixel_intercept": x_intercept,
            "x_pixels_per_log10_reduced_field": x_slope,
            "y_tick_centers_crop_px": Y_TICKS_CROP_PX.tolist(),
            "reduced_mobility_at_y_ticks_cm2_V_s": (
                REDUCED_MOBILITY_AT_Y_TICKS_CM2_V_S.tolist()),
            "y_axis": "log10",
            "y_pixel_intercept": y_intercept,
            "y_pixels_per_log10_reduced_mobility": y_slope,
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy tick localization; "
                "OpenCV Hough candidates followed by original-resolution "
                "visual overlay audit"
            ),
            "marker": "open circles, CHF2+ in CHF3 only",
            "point_count": len(PIXEL_POINTS),
            "center_placement_allowance_px": 10.0,
            "digitization_reduced_field_relative_bound": (
                10.0 ** (10.0 / abs(x_slope)) - 1.0),
            "digitization_reduced_mobility_relative_bound": (
                10.0 ** (10.0 / abs(y_slope)) - 1.0),
            "source_measurement_relative_uncertainty_range": [0.02, 0.04],
            "source_measurement_uncertainty_kept_separate": True,
            "visual_audit_status": "passed_original_resolution",
        },
        "claim_boundary": {
            "valid": (
                "measured CHF2+ reduced mobility in CHF3 over the plotted "
                "reduced-field support"
            ),
            "not_valid": [
                "CF3+ mobility in CHF3",
                "an elastic differential cross section",
                "a target-reactor ion species fraction",
                "a target sheath IEAD",
                "a TiO2 etch-depth fit",
            ],
        },
        "output": {
            "path": str(CSV_PATH.relative_to(ROOT)),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    return json.dumps(manifest(sha256(
        csv_payload.encode("utf-8")).hexdigest()), indent=2) + "\n"


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Basurto Figure-1 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Basurto digitization manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Basurto source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Basurto 600-dpi page render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Basurto render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    for point in PIXEL_POINTS:
        x = int(round(point.x_full_page_px))
        y = int(round(point.y_full_page_px))
        dark = int(np.sum(gray[y - 25:y + 26, x - 25:x + 26] < 120))
        if dark < 120:
            raise RuntimeError(f"insufficient marker support at ({x}, {y})")
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for index, point in enumerate(PIXEL_POINTS):
        x = point.x_full_page_px
        y = point.y_full_page_px
        draw.ellipse((x - 25, y - 25, x + 25, y + 25),
                     outline="#e41a1c", width=6)
        draw.text((x + 28, y - 28), str(index), fill="#e41a1c")
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render", type=Path)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    image = None
    if args.source_pdf is not None:
        verify_source_pdf(args.source_pdf)
    if args.render is not None:
        image = verify_render(args.render)
    if args.overlay is not None:
        if image is None:
            raise SystemExit("--overlay requires --render")
        draw_overlay(image, args.overlay)
    if args.check:
        verify_committed_files()
        return
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


if __name__ == "__main__":
    main()
