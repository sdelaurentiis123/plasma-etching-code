#!/usr/bin/env python3
"""Reproduce the unambiguous pure-C4F6 markers in Lan--Jeon Figure 7.

The source PDF is not redistributed.  The committed table is derived from a
checksum-pinned 600-dpi full-page render using log-axis calibration and
full-resolution marker centers.  Only the 18 filled circles whose identity
and continuity remain unambiguous after the C4F6/Ar curves merge are kept.

Because these drift measurements were used to regress the collision set,
they are a source-replay test, not independent validation of its branches.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "lan_jeon_2014_c4f6"
CSV_PATH = OUTPUT_DIRECTORY / "figure7_pure_c4f6_drift.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure7_digitization_manifest.json"
SOURCE_PDF_SHA256 = (
    "82d672f5b60611894a4584aa503e4c26d66aec3c66792550de46fb660f77aeb6"
)
RENDER_SHA256 = (
    "6cf27de7bb64e8cf9a34969b7ed9d44f5508bac0c6705eec3ef34f52fb540033"
)
RENDER_SIZE = (4961, 7017)
RENDER_DPI = 600

# Full-page major-tick centers located from the four axis spines/tick runs.
X_TICKS_PX = (576.0, 978.5, 1381.5, 1783.5, 2186.0)
LOG10_FIELD_AT_TICKS = (-1.0, 0.0, 1.0, 2.0, 3.0)
Y_TICKS_PX = (2084.5, 1644.5, 1204.0)
LOG10_PRINTED_W_AT_TICKS = (-1.0, 0.0, 1.0)


@dataclass(frozen=True)
class PixelPoint:
    x_px: float
    y_px: float


# Crop origin was (450, 700). Candidate glyphs were found with a Hough circle
# transform, then every retained center was inspected at original resolution.
# The four lower-field candidates are excluded because mixture glyphs overlap.
PIXEL_POINTS = tuple(PixelPoint(x + 450.0, y + 700.0) for x, y in (
    (1367.5, 508.5),
    (1392.5, 483.5),
    (1415.5, 457.5),
    (1435.5, 438.5),
    (1451.5, 407.5),
    (1478.5, 390.5),
    (1497.5, 377.5),
    (1524.5, 363.5),
    (1552.5, 340.5),
    (1576.5, 325.5),
    (1596.5, 313.5),
    (1615.5, 299.5),
    (1645.5, 276.5),
    (1673.5, 257.5),
    (1696.5, 239.5),
    (1716.5, 226.5),
    (1735.5, 212.5),
    (1767.5, 178.5),
))
FIELDNAMES = (
    "reduced_electric_field_Td",
    "drift_velocity_m_s",
    "marker_center_x_full_page_px",
    "marker_center_y_full_page_px",
    "marker",
    "field_digitization_relative_bound",
    "drift_digitization_relative_bound",
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _axis_coefficients(pixels, log_values) -> tuple[float, float]:
    pixel = np.asarray(pixels, dtype=float)
    value = np.asarray(log_values, dtype=float)
    slope = float(
        np.sum((pixel - np.mean(pixel)) * (value - np.mean(value)))
        / np.sum((pixel - np.mean(pixel)) ** 2)
    )
    intercept = float(np.mean(value) - slope * np.mean(pixel))
    return intercept, slope


def field_at_pixel(x_px: float) -> float:
    intercept, slope = _axis_coefficients(X_TICKS_PX, LOG10_FIELD_AT_TICKS)
    return 10.0 ** (intercept + slope * float(x_px))


def drift_at_pixel(y_px: float) -> float:
    intercept, slope = _axis_coefficients(
        Y_TICKS_PX, LOG10_PRINTED_W_AT_TICKS)
    # Figure ordinate is W / (10^6 cm/s), so one printed unit is 10^4 m/s.
    return 1.0e4 * 10.0 ** (intercept + slope * float(y_px))


def _relative_pixel_bound(slope: float, allowance_px: float = 5.0) -> float:
    return 10.0 ** (abs(slope) * allowance_px) - 1.0


def rows() -> list[dict[str, str]]:
    _, x_slope = _axis_coefficients(X_TICKS_PX, LOG10_FIELD_AT_TICKS)
    _, y_slope = _axis_coefficients(
        Y_TICKS_PX, LOG10_PRINTED_W_AT_TICKS)
    field_bound = _relative_pixel_bound(x_slope)
    drift_bound = _relative_pixel_bound(y_slope)
    return [{
        "reduced_electric_field_Td": f"{field_at_pixel(point.x_px):.9g}",
        "drift_velocity_m_s": f"{drift_at_pixel(point.y_px):.9g}",
        "marker_center_x_full_page_px": f"{point.x_px:.1f}",
        "marker_center_y_full_page_px": f"{point.y_px:.1f}",
        "marker": "filled_circle",
        "field_digitization_relative_bound": f"{field_bound:.9g}",
        "drift_digitization_relative_bound": f"{drift_bound:.9g}",
    } for point in PIXEL_POINTS]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_digest: str) -> dict[str, object]:
    x_intercept, x_slope = _axis_coefficients(
        X_TICKS_PX, LOG10_FIELD_AT_TICKS)
    y_intercept, y_slope = _axis_coefficients(
        Y_TICKS_PX, LOG10_PRINTED_W_AT_TICKS)
    return {
        "manifest_id": "LAN-JEON-2014-C4F6-FIG7-PURE-DRIFT-R1",
        "source": {
            "citation": (
                "P.-T. Lan and B.-H. Jeon, Journal of the Korean Physical "
                "Society 64, 1320-1326 (2014)"
            ),
            "doi": "10.3938/jkps.64.1320",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 6,
            "print_page": 1325,
            "figure": "Figure 7",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "source_measurement": "Goyette et al. reference [17]",
            "redistribution": "source pixels are not committed",
        },
        "pixel_calibration": {
            "full_page_x_major_tick_px": list(X_TICKS_PX),
            "log10_field_Td_at_major_ticks": list(LOG10_FIELD_AT_TICKS),
            "field_log10_intercept": x_intercept,
            "field_log10_per_horizontal_pixel": x_slope,
            "full_page_y_major_tick_px": list(Y_TICKS_PX),
            "log10_printed_W_at_major_ticks": list(
                LOG10_PRINTED_W_AT_TICKS),
            "drift_log10_intercept": y_intercept,
            "drift_log10_per_vertical_pixel": y_slope,
            "ordinate_conversion": "1 printed unit = 1e6 cm/s = 1e4 m/s",
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy axis localization; "
                "Hough-assisted candidate localization; original-resolution "
                "visual acceptance and full-page marker centers"
            ),
            "series_policy": (
                "retain only the 18 unambiguous filled-circle pure-C4F6 "
                "markers; exclude four lower-field candidates overlapped by "
                "C4F6/Ar mixture glyphs"
            ),
            "marker_center_allowance_px": 5.0,
            "field_digitization_relative_bound": _relative_pixel_bound(
                x_slope),
            "drift_digitization_relative_bound": _relative_pixel_bound(
                y_slope),
            "source_measurement_uncertainty": "not reported in this paper",
            "visual_audit_status": "passed_original_resolution",
        },
        "evidence_boundary": {
            "point_count": len(PIXEL_POINTS),
            "field_range_Td": [
                field_at_pixel(PIXEL_POINTS[0].x_px),
                field_at_pixel(PIXEL_POINTS[-1].x_px),
            ],
            "supports_source_replay": True,
            "supports_independent_validation_of_regressed_set": False,
            "supports_resolved_product_branching": False,
            "supports_reactor_state_prediction": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
            "note": (
                "The same drift data constrained the printed effective "
                "cross-section set; agreement catches transcription/solver "
                "errors but is not a held-out physical validation."
            ),
        },
        "output": {
            "path": (
                "data/experimental/lan_jeon_2014_c4f6/"
                "figure7_pure_c4f6_drift.csv"
            ),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = sha256(csv_payload.encode()).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def _verify_render(path: Path) -> Image.Image:
    if _sha(path) != RENDER_SHA256:
        raise RuntimeError("Lan--Jeon Figure-7 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Lan--Jeon render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    axis_checks = {
        "left": np.mean(np.min(gray[825:2090, 570:585], axis=1) < 96),
        "right": np.mean(np.min(gray[825:2090, 2300:2315], axis=1) < 96),
        "top": np.mean(np.min(gray[825:840, 570:2315], axis=0) < 96),
        "bottom": np.mean(np.min(gray[2078:2092, 570:2315], axis=0) < 96),
    }
    failed = {key: float(value) for key, value in axis_checks.items() if value < .65}
    if failed:
        raise RuntimeError(f"Lan--Jeon Figure-7 axis verification failed: {failed}")
    for index, point in enumerate(PIXEL_POINTS):
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark = int(np.sum(gray[y - 13:y + 14, x - 13:x + 14] < 112))
        if dark < 110:
            raise RuntimeError(
                f"insufficient filled-circle support at marker {index}: {dark}")
    return image


def _draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for point in PIXEL_POINTS:
        radius = 20
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline="#e41a1c",
            width=4,
        )
        draw.line(
            (point.x_px - 12, point.y_px, point.x_px + 12, point.y_px),
            fill="#377eb8",
            width=3,
        )
        draw.line(
            (point.x_px, point.y_px - 12, point.x_px, point.y_px + 12),
            fill="#377eb8",
            width=3,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def _write() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def _check() -> None:
    payload = csv_text()
    if not CSV_PATH.exists() or CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Lan--Jeon Figure-7 CSV is stale")
    if (
        not MANIFEST_PATH.exists()
        or MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload)
    ):
        raise RuntimeError("committed Lan--Jeon Figure-7 manifest is stale")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render", type=Path)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()
    if args.write:
        _write()
    if args.check:
        _check()
    if args.source_pdf is not None and _sha(args.source_pdf) != SOURCE_PDF_SHA256:
        raise RuntimeError("Lan--Jeon source PDF checksum changed")
    image = _verify_render(args.render) if args.render is not None else None
    if args.overlay is not None:
        if image is None:
            raise ValueError("--overlay requires --render")
        _draw_overlay(image, args.overlay)
    if not args.write and not args.check and args.source_pdf is None and image is None:
        parser.error("select --write, --check, --source-pdf, or --render")
    print(json.dumps({
        "status": "verified",
        "point_count": len(PIXEL_POINTS),
        "source_pdf_verified": args.source_pdf is not None,
        "render_verified": image is not None,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
