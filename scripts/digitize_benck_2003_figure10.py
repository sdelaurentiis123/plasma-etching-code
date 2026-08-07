#!/usr/bin/env python3
"""Reproduce Benck Figure 10 and cross-check it against Figure 9.

Figure 10 sweeps pressure in a 50/50 C4F6/Ar ICP.  The 1.33 Pa column is the
same physical condition as the 50% column of Figure 9, so the two independently
drawn panels provide an unusually useful pixel-level digitization cross-check.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.digitize_benck_2003_figure9 import (
    CSV_PATH as FIGURE9_CSV_PATH,
    RENDER_SHA256,
    RENDER_SIZE,
    SOURCE_PDF_SHA256,
    verify_render as verify_page_render,
    verify_source_pdf,
)


OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "benck_2003_c4f6"
CSV_PATH = OUTPUT_DIRECTORY / "figure10_pressure_mass_resolved_ion_current.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure10_digitization_manifest.json"

Y_TICKS_PX = (1065.5, 1637.5, 2213.5)
LOG10_CURRENT_AT_TICKS = (-1.0, -2.0, -3.0)
X_AT_PRESSURE_PA = {
    0.67: 3333.0,
    1.33: 3513.0,
    2.66: 3878.0,
}


@dataclass(frozen=True)
class PixelPoint:
    pressure_Pa: float
    species: str
    y_px: float
    marker: str

    @property
    def x_px(self) -> float:
        return X_AT_PRESSURE_PA[self.pressure_Pa]


PIXEL_POINTS = (
    PixelPoint(0.67, "total_positive_ion_current", 852.6, "filled_circle"),
    PixelPoint(1.33, "total_positive_ion_current", 873.0, "filled_circle"),
    PixelPoint(2.66, "total_positive_ion_current", 925.8, "filled_circle"),
    PixelPoint(0.67, "Ar+", 999.0, "open_circle"),
    PixelPoint(1.33, "Ar+", 999.0, "open_circle"),
    PixelPoint(2.66, "Ar+", 1121.0, "open_circle"),
    PixelPoint(0.67, "CF+", 1349.5, "filled_square"),
    PixelPoint(1.33, "CF+", 1447.5, "filled_square"),
    PixelPoint(2.66, "CF+", 1553.1, "filled_square"),
    PixelPoint(0.67, "CF2+", 1492.7, "filled_down_triangle"),
    PixelPoint(1.33, "CF2+", 1601.3, "filled_down_triangle"),
    PixelPoint(2.66, "CF2+", 1716.6, "filled_down_triangle"),
    PixelPoint(0.67, "CF3+", 1783.7, "filled_up_triangle"),
    PixelPoint(1.33, "CF3+", 1823.6, "filled_up_triangle"),
    PixelPoint(2.66, "CF3+", 1789.0, "filled_up_triangle"),
)

FIELDNAMES = (
    "pressure_Pa",
    "pressure_mTorr",
    "c4f6_feed_percent",
    "ar_feed_percent",
    "species",
    "ion_current_density_mA_cm2",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "marker",
    "digitization_relative_bound",
    "source_transmission_relative_uncertainty",
)


def _log_axis_coefficients() -> tuple[float, float]:
    y = np.asarray(Y_TICKS_PX, dtype=float)
    value = np.asarray(LOG10_CURRENT_AT_TICKS, dtype=float)
    slope = float(
        np.sum((y - np.mean(y)) * (value - np.mean(value)))
        / np.sum((y - np.mean(y)) ** 2)
    )
    intercept = float(np.mean(value) - slope * np.mean(y))
    return intercept, slope


def current_at_pixel(y_px: float) -> float:
    intercept, slope = _log_axis_coefficients()
    return 10.0 ** (intercept + slope * float(y_px))


def rows() -> list[dict[str, str]]:
    return [
        {
            "pressure_Pa": f"{point.pressure_Pa:.2f}",
            "pressure_mTorr": f"{point.pressure_Pa / 0.133322368:.3f}",
            "c4f6_feed_percent": "50",
            "ar_feed_percent": "50",
            "species": point.species,
            "ion_current_density_mA_cm2": (
                f"{current_at_pixel(point.y_px):.9g}"
            ),
            "marker_center_x_full_page_px": f"{point.x_px:.1f}",
            "marker_center_y_full_page_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            "digitization_relative_bound": "0.102",
            "source_transmission_relative_uncertainty": "0.20",
        }
        for point in PIXEL_POINTS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _figure9_50_percent() -> dict[str, float]:
    with FIGURE9_CSV_PATH.open(newline="", encoding="utf-8") as stream:
        return {
            row["species"]: float(row["ion_current_density_mA_cm2"])
            for row in csv.DictReader(stream)
            if row["c4f6_feed_percent"] == "50"
        }


def cross_figure_relative_differences() -> dict[str, float]:
    figure9 = _figure9_50_percent()
    figure10 = {
        row["species"]: float(row["ion_current_density_mA_cm2"])
        for row in rows()
        if row["pressure_Pa"] == "1.33"
    }
    return {
        species: figure10[species] / figure9[species] - 1.0
        for species in sorted(figure9)
    }


def manifest(csv_digest: str) -> dict[str, object]:
    intercept, slope = _log_axis_coefficients()
    differences = cross_figure_relative_differences()
    return {
        "manifest_id": "BENCK-2003-C4F6-FIG10-PRESSURE-ION-CURRENT-R1",
        "source": {
            "citation": (
                "E. C. Benck, A. Goyette, and Y. Wang, Journal of "
                "Applied Physics 94, 1382-1389 (2003)"
            ),
            "doi": "10.1063/1.1586978",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 6,
            "print_page": 1387,
            "figure": "Figure 10",
            "render_dpi": 600,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "redistribution": "source pixels are not committed",
        },
        "experiment": {
            "c4f6_feed_percent": 50,
            "ar_feed_percent": 50,
            "power_W": 200.0,
            "pressure_setpoints_Pa": [0.67, 1.33, 2.66],
            "surface_bias": "grounded water-cooled lower electrode",
            "oxygen_added": False,
            "corrected_transmission_relative_uncertainty": 0.20,
        },
        "pixel_calibration": {
            "full_page_x_px_by_pressure_Pa": {
                f"{key:.2f}": value for key, value in X_AT_PRESSURE_PA.items()
            },
            "full_page_y_major_tick_px": list(Y_TICKS_PX),
            "log10_current_at_major_ticks": list(LOG10_CURRENT_AT_TICKS),
            "log10_current_intercept": intercept,
            "log10_current_per_vertical_pixel": slope,
            "major_tick_replay_mA_cm2": [
                current_at_pixel(y_px) for y_px in Y_TICKS_PX
            ],
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy axis localization; "
                "full-resolution glyph centers; original-resolution overlay"
            ),
            "series_policy": (
                "same five unambiguous total/Ar+/CF+/CF2+/CF3+ series as "
                "Figure 9; exclude the overlapping remaining glyphs"
            ),
            "vertical_pixel_allowance": 24.0,
            "digitization_relative_bound": 0.102,
            "maximum_exact_bound_from_24_pixels": (
                10.0 ** (abs(slope) * 24.0) - 1.0
            ),
            "visual_audit_status": "passed_original_resolution",
        },
        "cross_figure_reconciliation": {
            "condition": "50% C4F6, 50% Ar, 200 W, 1.33 Pa",
            "figure9_column": "50_percent_C4F6",
            "figure10_column": "1.33_Pa",
            "signed_relative_differences": differences,
            "maximum_absolute_relative_difference": max(
                abs(value) for value in differences.values()
            ),
            "passed_within_individual_digitization_bound": all(
                abs(value) <= 0.102 for value in differences.values()
            ),
            "interpretation": (
                "the same condition was independently drawn in two panels; "
                "agreement is a digitization cross-check, not a second "
                "physical experiment"
            ),
        },
        "claim_boundary": {
            "valid": (
                "held-out pressure-response target for a 50/50 C4F6/Ar "
                "reactor provider at this apparatus and power"
            ),
            "not_valid": [
                "a Krueger reactor boundary",
                "a biased-sheath species-resolved IEAD",
                "a stable-neutral C4F6 flux",
                "a surface-yield or feature-depth calibration",
            ],
        },
        "output": {
            "path": (
                "data/experimental/benck_2003_c4f6/"
                "figure10_pressure_mass_resolved_ion_current.csv"
            ),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Benck Figure-10 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Benck Figure-10 manifest is stale")


def verify_render(path: Path) -> Image.Image:
    image = verify_page_render(path)
    gray = np.asarray(image.convert("L"))
    for point in PIXEL_POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark = int(np.sum(gray[y - 24:y + 25, x - 24:x + 25] < 112))
        if dark < 80:
            raise RuntimeError(
                f"insufficient Figure-10 marker support at {point.species} "
                f"{point.pressure_Pa} Pa: {dark}"
            )
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {
        "total_positive_ion_current": "#e41a1c",
        "Ar+": "#377eb8",
        "CF+": "#4daf4a",
        "CF2+": "#ff7f00",
        "CF3+": "#984ea3",
    }
    for point in PIXEL_POINTS:
        color = colors[point.species]
        radius = 26
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline=color,
            width=5,
        )
        draw.line(
            (point.x_px - 16, point.y_px, point.x_px + 16, point.y_px),
            fill=color,
            width=3,
        )
        draw.line(
            (point.x_px, point.y_px - 16, point.x_px, point.y_px + 16),
            fill=color,
            width=3,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render", type=Path)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()

    verify_committed_files()
    if args.source_pdf is not None:
        verify_source_pdf(args.source_pdf)
    image = verify_render(args.render) if args.render is not None else None
    if args.overlay is not None:
        if image is None:
            raise ValueError("--overlay requires --render")
        draw_overlay(image, args.overlay)
    print(
        json.dumps(
            {
                "status": "verified",
                "point_count": len(PIXEL_POINTS),
                "maximum_cross_figure_relative_difference": max(
                    abs(value)
                    for value in cross_figure_relative_differences().values()
                ),
                "source_pdf_verified": args.source_pdf is not None,
                "render_verified": image is not None,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
