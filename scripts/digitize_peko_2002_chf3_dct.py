#!/usr/bin/env python3
"""Reproduce Peko 2002 Figure-2 CF3+ + CHF3 DCT digitization.

The AIP/NIST source pixels are not redistributed.  The committed table is
reconstructed from marker centers on a 400-dpi Poppler render.  When supplied
locally, the source PDF and render are checksum-verified and a visual overlay
can be emitted for original-resolution audit.
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
OUTPUT_DIR = ROOT / "data" / "experimental" / "peko_2002_chf3"
CSV_PATH = OUTPUT_DIR / "figure2_cf3_chf3_dct_sum.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "c0ab4fac39a611e364efc29c12632590d37e319cd8274382ae13be5a2d33d99c"
)
RENDER_SHA256 = (
    "b4b1cc4422a697fdae471b5c86b6b23316d7a13fb9077a5b567899e6b0c9ac07"
)
RENDER_SIZE = (3400, 4446)
RENDER_DPI = 400

# Full-page pixel coordinates.  The x axis is logarithmic in collision
# energy; the y axis is linear in cross section.
X_TICKS_PX = np.array([2012.5, 2260.5, 2404.5, 2507.0, 2587.5, 2838.0])
ENERGY_AT_X_TICKS_EV = np.array([20.0, 40.0, 60.0, 80.0, 100.0, 200.0])
Y_TICKS_PX = np.array([
    530.5, 658.5, 786.5, 915.0, 1043.5, 1173.5, 1301.5, 1429.5,
])
CROSS_SECTION_AT_Y_TICKS_A2 = np.array([
    8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0,
])


@dataclass(frozen=True)
class PixelPoint:
    x_px: float
    y_px: float


# Full-page coordinates: crop origin (1650, 250) plus original-resolution
# star centers.  A PIL overlay was inspected at original resolution.
PIXEL_POINTS = tuple(
    PixelPoint(x + 1650.0, y + 250.0)
    for x, y in (
        (371, 970), (405, 944), (447, 916), (476, 863), (512, 808),
        (570, 738), (611, 679), (656, 613), (698, 546), (771, 494),
        (837, 455), (895, 422), (940, 388), (979, 370), (1023, 352),
        (1055, 328), (1077, 310), (1101, 295), (1131, 284),
        (1153, 266), (1179, 278),
    )
)

FIELDNAMES = (
    "collision_energy_eV",
    "dct_sum_cross_section_A2",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "digitization_energy_relative_bound",
    "digitization_cross_section_A2_bound",
    "source_measurement_relative_uncertainty",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _x_coefficients() -> tuple[float, float]:
    matrix = np.column_stack((
        np.ones(X_TICKS_PX.size), np.log10(ENERGY_AT_X_TICKS_EV)))
    intercept, slope = np.linalg.lstsq(
        matrix, X_TICKS_PX, rcond=None)[0]
    return float(intercept), float(slope)


def _y_coefficients() -> tuple[float, float]:
    matrix = np.column_stack((np.ones(Y_TICKS_PX.size), Y_TICKS_PX))
    intercept, slope = np.linalg.lstsq(
        matrix, CROSS_SECTION_AT_Y_TICKS_A2, rcond=None)[0]
    return float(intercept), float(slope)


def energy_at_pixel(x_px: float) -> float:
    intercept, slope = _x_coefficients()
    return 10.0 ** ((float(x_px) - intercept) / slope)


def cross_section_at_pixel(y_px: float) -> float:
    intercept, slope = _y_coefficients()
    return intercept + slope * float(y_px)


def rows() -> list[dict[str, str]]:
    _, x_slope = _x_coefficients()
    _, y_slope = _y_coefficients()
    energy_relative_bound = 10.0 ** (6.0 / x_slope) - 1.0
    cross_section_bound = abs(6.0 * y_slope)
    return [
        {
            "collision_energy_eV": f"{energy_at_pixel(point.x_px):.8g}",
            "dct_sum_cross_section_A2": (
                f"{cross_section_at_pixel(point.y_px):.8g}"
            ),
            "marker_center_x_full_page_px": f"{point.x_px:.1f}",
            "marker_center_y_full_page_px": f"{point.y_px:.1f}",
            "digitization_energy_relative_bound": (
                f"{energy_relative_bound:.8g}"
            ),
            "digitization_cross_section_A2_bound": (
                f"{cross_section_bound:.8g}"
            ),
            "source_measurement_relative_uncertainty": "0.25",
        }
        for point in PIXEL_POINTS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_digest: str) -> dict:
    x_intercept, x_slope = _x_coefficients()
    y_intercept, y_slope = _y_coefficients()
    return {
        "manifest_id": "PEKO-2002-CHF3-FIG2-DCT-R1",
        "source": {
            "citation": (
                "B. L. Peko, R. L. Champion, M. V. V. S. Rao, and J. K. "
                "Olthoff, Journal of Applied Physics 92, 1657-1662 (2002)"
            ),
            "title": "Measured cross sections and ion energies for a CHF3 discharge",
            "doi": "10.1063/1.1491276",
            "official_pdf": (
                "https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=9195"
            ),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 3,
            "print_page": 1659,
            "figure": "Figure 2, summed DCT series",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "redistribution": "source pixels are not committed",
        },
        "pixel_calibration": {
            "x_tick_centers_full_page_px": X_TICKS_PX.tolist(),
            "collision_energy_at_x_ticks_eV": ENERGY_AT_X_TICKS_EV.tolist(),
            "x_axis": "log10",
            "x_pixel_intercept": x_intercept,
            "x_pixels_per_log10_energy": x_slope,
            "y_tick_centers_full_page_px": Y_TICKS_PX.tolist(),
            "cross_section_at_y_ticks_A2": (
                CROSS_SECTION_AT_Y_TICKS_A2.tolist()
            ),
            "y_axis": "linear",
            "cross_section_intercept_A2": y_intercept,
            "cross_section_A2_per_vertical_pixel": y_slope,
        },
        "digitization": {
            "method": (
                "400-dpi Poppler render; PIL/NumPy dark-pixel tick "
                "localization and density-peak marker detection; "
                "original-resolution PIL overlay visually audited"
            ),
            "marker": "asterisk, summed DCT curve only",
            "point_count": len(PIXEL_POINTS),
            "center_placement_allowance_px": 6.0,
            "digitization_energy_relative_bound": (
                10.0 ** (6.0 / x_slope) - 1.0
            ),
            "digitization_cross_section_A2_bound": abs(6.0 * y_slope),
            "source_measurement_relative_uncertainty": 0.25,
            "source_measurement_uncertainty_kept_separate": True,
            "visual_audit_status": "passed_original_resolution",
        },
        "claim_boundary": {
            "valid": (
                "summed dissociative-charge-transfer cross section for "
                "CF3+ + CHF3 over the plotted collision-energy support"
            ),
            "not_valid": [
                "an elastic or total momentum-transfer cross section",
                "an angular scattering law",
                "SF5+ transport in SF6",
                "a species-resolved target-reactor ion flux",
                "a TiO2 etch-depth fit",
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
        raise RuntimeError("committed Peko Figure-2 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Peko digitization manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Peko source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Peko 400-dpi page-3 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Peko render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    for point in PIXEL_POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark = int(np.sum(gray[y - 18:y + 19, x - 18:x + 19] < 120))
        if dark < 65:
            raise RuntimeError(f"insufficient marker support at ({x}, {y})")
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for index, point in enumerate(PIXEL_POINTS):
        radius = 18
        draw.ellipse(
            (point.x_px - radius, point.y_px - radius,
             point.x_px + radius, point.y_px + radius),
            outline="#e41a1c", width=5,
        )
        draw.text((point.x_px + 20, point.y_px - 20), str(index), fill="#e41a1c")
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render", type=Path)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()

    payload = csv_text()
    if args.check:
        verify_committed_files()
    else:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(payload, encoding="utf-8")
        MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")
    if args.source_pdf is not None:
        verify_source_pdf(args.source_pdf)
    if args.render is not None:
        image = verify_render(args.render)
        if args.overlay is not None:
            draw_overlay(image, args.overlay)
    elif args.overlay is not None:
        raise SystemExit("--overlay requires --render")


if __name__ == "__main__":
    main()
