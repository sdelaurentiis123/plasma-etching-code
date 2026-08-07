#!/usr/bin/env python3
"""Reproduce the Benck 2003 C4F6/Ar mass-resolved ion-current board.

The source is an AIP scan that is not redistributed.  This script reproduces
the committed numerical table from full-page 600-dpi marker coordinates and,
when given the locally held render, verifies its checksum and draws a visual
overlay.  Only the five reactor-boundary series with unambiguous glyph
continuity are transcribed: total, Ar+, CF+, CF2+, and CF3+.

These measurements grade a C4F6 reactor provider.  They are not a boundary
condition for Krueger's different, biased, oxygen-containing reactor.
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
CSV_PATH = OUTPUT_DIRECTORY / "figure9_mass_resolved_ion_current.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "96ef064afb0f804d3d853adcc275a63cd52507b0ea7e39f921d6cf6209e86d08"
)
RENDER_SHA256 = (
    "f423f6992e049ba40f47a540eb56450e9c080dddb4b03ad45ee8883b109a71fd"
)
RENDER_SIZE = (5209, 6760)
RENDER_DPI = 600

# Full-page major-tick centers on the logarithmic ordinate.
Y_TICKS_PX = (4182.5, 4753.5, 5333.5)
LOG10_CURRENT_AT_TICKS = (-1.0, -2.0, -3.0)

# Full-page x-tick centers for the four published feed setpoints.
X_AT_C4F6_PERCENT = {
    25: 1233.5,
    50: 1461.5,
    75: 1689.5,
    100: 1915.7,
}


@dataclass(frozen=True)
class PixelPoint:
    c4f6_percent: int
    species: str
    y_px: float
    marker: str

    @property
    def x_px(self) -> float:
        return X_AT_C4F6_PERCENT[self.c4f6_percent]


PIXEL_POINTS = (
    PixelPoint(25, "total_positive_ion_current", 3938.0, "filled_circle"),
    PixelPoint(50, "total_positive_ion_current", 3988.0, "filled_circle"),
    PixelPoint(75, "total_positive_ion_current", 4006.2, "filled_circle"),
    PixelPoint(100, "total_positive_ion_current", 4070.4, "filled_circle"),
    PixelPoint(25, "Ar+", 3999.5, "open_circle"),
    PixelPoint(50, "Ar+", 4115.6, "open_circle"),
    PixelPoint(75, "Ar+", 4271.1, "open_circle"),
    PixelPoint(25, "CF+", 4687.3, "filled_square"),
    PixelPoint(50, "CF+", 4571.7, "filled_square"),
    PixelPoint(75, "CF+", 4485.7, "filled_square"),
    PixelPoint(100, "CF+", 4469.5, "filled_square"),
    PixelPoint(25, "CF2+", 4878.7, "filled_down_triangle"),
    PixelPoint(50, "CF2+", 4724.8, "filled_down_triangle"),
    PixelPoint(75, "CF2+", 4637.3, "filled_down_triangle"),
    PixelPoint(100, "CF2+", 4616.4, "filled_down_triangle"),
    PixelPoint(25, "CF3+", 5114.2, "filled_up_triangle"),
    PixelPoint(50, "CF3+", 4950.4, "filled_up_triangle"),
    PixelPoint(75, "CF3+", 4862.3, "filled_up_triangle"),
    PixelPoint(100, "CF3+", 4783.0, "filled_up_triangle"),
)

FIELDNAMES = (
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


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


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
            "c4f6_feed_percent": str(point.c4f6_percent),
            "ar_feed_percent": str(100 - point.c4f6_percent),
            "species": point.species,
            "ion_current_density_mA_cm2": (
                f"{current_at_pixel(point.y_px):.9g}"
            ),
            "marker_center_x_full_page_px": f"{point.x_px:.1f}",
            "marker_center_y_full_page_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            # Twenty-four vertical pixels span a full ambiguous marker radius
            # at the overlapping 25% and 100% columns.  On this log axis that
            # is a 10.08% multiplicative allowance.
            "digitization_relative_bound": "0.101",
            # Benck estimates corrected ion transmission as uniform to 20%.
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


def manifest(csv_digest: str) -> dict[str, object]:
    intercept, slope = _log_axis_coefficients()
    tick_replay = [
        current_at_pixel(y_px) for y_px in Y_TICKS_PX
    ]
    return {
        "manifest_id": "BENCK-2003-C4F6-FIG9-ION-CURRENT-R1",
        "source": {
            "citation": (
                "E. C. Benck, A. Goyette, and Y. Wang, Journal of "
                "Applied Physics 94, 1382-1389 (2003)"
            ),
            "title": (
                "Ion energy distribution and optical measurements in "
                "high-density, inductively coupled C4F6 discharges"
            ),
            "doi": "10.1063/1.1586978",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 6,
            "print_page": 1387,
            "figure": "Figure 9",
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
            "pressure_mTorr": 10.0,
            "power_W": 200.0,
            "surface_bias": "grounded water-cooled lower electrode",
            "oxygen_added": False,
            "sampling": (
                "10 um grounded side orifice; Faraday cup at the same "
                "radius and height normalizes total current"
            ),
            "energy_resolution_eV": 1.0,
            "estimated_energy_scale_uncertainty_eV": 1.0,
            "corrected_transmission_relative_uncertainty": 0.20,
            "ordinate": (
                "singly charged positive-ion current density as printed, "
                "mA/cm2; no particle-flux conversion is applied"
            ),
        },
        "pixel_calibration": {
            "full_page_x_px_by_c4f6_percent": {
                str(key): value for key, value in X_AT_C4F6_PERCENT.items()
            },
            "full_page_y_major_tick_px": list(Y_TICKS_PX),
            "log10_current_at_major_ticks": list(LOG10_CURRENT_AT_TICKS),
            "log10_current_intercept": intercept,
            "log10_current_per_vertical_pixel": slope,
            "major_tick_replay_mA_cm2": tick_replay,
            "axis_fit": (
                "least-squares affine fit of log10(current) to the three "
                "full-page 10^-1, 10^-2, and 10^-3 major-tick centers"
            ),
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy dark-pixel axis and "
                "x-tick localization; full-resolution glyph-center "
                "transcription; original-resolution visual overlay"
            ),
            "series_policy": (
                "transcribe only total, Ar+, CF+, CF2+, and CF3+, whose "
                "glyph identity and curve continuity are unambiguous; "
                "exclude overlapping C/F/Cx and SiFx markers rather than "
                "guessing them"
            ),
            "vertical_pixel_allowance": 24.0,
            "digitization_relative_bound": 0.101,
            "maximum_exact_bound_from_24_pixels": (
                10.0 ** (abs(slope) * 24.0) - 1.0
            ),
            "source_measurement_uncertainty_kept_separate": True,
            "visual_audit_status": "passed_original_resolution",
        },
        "derived_checks": {
            "point_count": len(PIXEL_POINTS),
            "series_point_count": {
                name: sum(point.species == name for point in PIXEL_POINTS)
                for name in (
                    "total_positive_ion_current", "Ar+", "CF+", "CF2+", "CF3+"
                )
            },
            "ar_100_percent_c4f6_marker": "not_plotted",
            "total_current_decreases": True,
            "cf_and_cf2_level_between_75_and_100_percent": True,
            "cf3_continues_to_increase": True,
        },
        "claim_boundary": {
            "valid": (
                "quantitative held-out reactor-provider target for total and "
                "mass-resolved positive-ion current in this C4F6/Ar ICP"
            ),
            "not_valid": [
                "a boundary normalization for the Krueger reactor",
                "a species mixture for a high-power biased CCP",
                "a stable-C4F6 neutral flux",
                "a C4F6/SiO2 surface-yield law",
                "a feature-depth fit or prediction",
            ],
        },
        "output": {
            "path": (
                "data/experimental/benck_2003_c4f6/"
                "figure9_mass_resolved_ion_current.csv"
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
        raise RuntimeError("committed Benck Figure-9 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Benck Figure-9 manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Benck source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Benck 600-dpi page-6 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Benck render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    axis_checks = {
        "left": np.mean(np.min(gray[3830:5640, 1000:1018], axis=1) < 96),
        "right": np.mean(np.min(gray[3830:5640, 2460:2478], axis=1) < 96),
        "top": np.mean(np.min(gray[3828:3845, 1000:2478], axis=0) < 96),
        "bottom": np.mean(np.min(gray[5625:5642, 1000:2478], axis=0) < 96),
    }
    failed = {
        key: float(value)
        for key, value in axis_checks.items()
        if value < 0.65
    }
    if failed:
        raise RuntimeError(f"Benck Figure-9 axis verification failed: {failed}")
    for point in PIXEL_POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark = int(np.sum(gray[y - 24:y + 25, x - 24:x + 25] < 112))
        if dark < 80:
            raise RuntimeError(
                f"insufficient marker support at {point.species} "
                f"{point.c4f6_percent}%: {dark}"
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
                "source_pdf_verified": args.source_pdf is not None,
                "render_verified": image is not None,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
