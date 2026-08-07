#!/usr/bin/env python3
"""Replay the Karahashi 2007 Figure-9 product-timing digitization.

Figure 9 reports arbitrary-unit time-of-flight traces for products desorbed
from SiO2 by a pulsed, mass-selected CF3+ beam.  Only the time coordinate of
the highest measured marker in each panel is digitized.  The ordinate is not
used: it has no numerical labels and the source does not establish a common
gain across the three product panels.
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
SOURCE_PDF = (
    ROOT / "research_sources" / "karahashi_2007_hyomen_kagaku_28_60.pdf"
)
DEFAULT_RENDER = (
    ROOT / "tmp" / "pdfs" / "karahashi_figures" / "pages5_6_600dpi.png"
)
OUTPUT_DIR = ROOT / "data" / "experimental" / "karahashi_2007"
CSV_PATH = OUTPUT_DIR / "figure9_cf3_product_peak_times.csv"
MANIFEST_PATH = OUTPUT_DIR / "figure9_product_timing_manifest.json"

SOURCE_PDF_SHA256 = (
    "093b18b91b0a6d910fc414779ee8320b7a046ac4cad38ef5de0b7f2dd25a2d79"
)
RENDER_SHA256 = (
    "021f39694507753be117f19fe463b596402a4c701b0507e8fb99d2fa2b25ce33"
)
RENDER_SIZE = (4961, 7017)
CROP_BOUNDS = (2800, 2250, 4400, 6050)

# Coordinates in the exact crop above.  All panels share the same horizontal
# scale.  The vertical axes independently locate the three plot interiors.
X_AT_0_MS = 211.0
X_AT_1_MS = 1240.0
PANEL_AXES_Y = {
    "SiF": (297.0, 1186.0),
    "SiF2": (1533.0, 2423.0),
    "SiF4": (2738.0, 3627.0),
}
PIXEL_CENTER_BOUND = 6.0
TIME_DIGITIZATION_BOUND_MS = (
    PIXEL_CENTER_BOUND / (X_AT_1_MS - X_AT_0_MS)
)

# This value is filled only after the generated overlay has been inspected at
# original resolution.  The overlay is a QA artifact, not source evidence.
VISUALLY_INSPECTED_OVERLAY_SHA256 = (
    "4c4ad44675186aeb1a717018ea23e9c6c9fa8c7e49d97e7132532c068253bede"
)


@dataclass(frozen=True)
class PeakMarker:
    product: str
    x_px: float
    y_px: float


PEAK_MARKERS = (
    PeakMarker("SiF", 457.0, 412.0),
    PeakMarker("SiF2", 519.0, 1648.0),
    PeakMarker("SiF4", 767.0, 2867.0),
)

FIELDNAMES = (
    "incident_species",
    "target",
    "energy_eV",
    "incidence_angle_deg",
    "beam_pulse_fwhm_us",
    "desorbed_product",
    "peak_measured_sample_time_ms",
    "marker_center_x_crop_px",
    "marker_center_y_crop_px",
    "time_digitization_bound_ms",
    "ordinate_status",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _time_at_pixel(x_px: float) -> float:
    return (x_px - X_AT_0_MS) / (X_AT_1_MS - X_AT_0_MS)


def rows() -> list[dict[str, str]]:
    return [
        {
            "incident_species": "CF3+",
            "target": "SiO2",
            "energy_eV": "1000",
            "incidence_angle_deg": "30",
            "beam_pulse_fwhm_us": "100",
            "desorbed_product": marker.product,
            "peak_measured_sample_time_ms": (
                f"{_time_at_pixel(marker.x_px):.4f}"
            ),
            "marker_center_x_crop_px": f"{marker.x_px:.1f}",
            "marker_center_y_crop_px": f"{marker.y_px:.1f}",
            "time_digitization_bound_ms": (
                f"{TIME_DIGITIZATION_BOUND_MS:.3f}"
            ),
            "ordinate_status": (
                "arbitrary_units_not_digitized_no_cross_panel_comparison"
            ),
        }
        for marker in PEAK_MARKERS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=FIELDNAMES, lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    peak_times = {
        row["desorbed_product"]: float(row["peak_measured_sample_time_ms"])
        for row in rows()
    }
    return {
        "manifest_id": "KARAHASHI-2007-FIG9-CF3-PRODUCT-TIMING-R1",
        "source": {
            "citation": (
                "K. Karahashi, Hyomen Kagaku 28, 60-66 (2007), Figure 9"
            ),
            "local_pdf": (
                "research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf"
            ),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 5,
            "print_page": 64,
            "figure": "Figure 9",
            "render_dpi": 600,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "figure_crop_bounds_px": list(CROP_BOUNDS),
        },
        "experiment_scope": {
            "target": "SiO2",
            "incident_ion": "CF3+",
            "beam_selection": "mass selected",
            "neutral_radical_flux": "none",
            "ion_energy_eV": 1000,
            "incidence_angle_deg": 30,
            "beam_pulse_fwhm_us": 100,
            "detected_products": ["SiF", "SiF2", "SiF4"],
            "quantity": (
                "time-of-flight detector intensity in arbitrary units"
            ),
        },
        "pixel_calibration": {
            "crop_bounds_full_render_px": list(CROP_BOUNDS),
            "x_at_0_ms": X_AT_0_MS,
            "x_at_1_ms": X_AT_1_MS,
            "panel_axis_y_bounds_crop_px": {
                product: list(bounds)
                for product, bounds in PANEL_AXES_Y.items()
            },
            "transform": (
                "time_ms = (x_px - 211.0) / (1240.0 - 211.0)"
            ),
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; exact Pillow crop; dark-pixel "
                "axis localization; plus-marker template scan; "
                "original-resolution visual overlay"
            ),
            "transcribed_quantity": (
                "time coordinate of the highest measured plus marker in "
                "each product panel"
            ),
            "pixel_center_bound": PIXEL_CENTER_BOUND,
            "time_digitization_bound_ms": TIME_DIGITIZATION_BOUND_MS,
            "ordinate_policy": (
                "do not digitize intensity because the ordinate is "
                "unlabelled arbitrary units and common inter-panel gain is "
                "not established"
            ),
            "full_resolution_visual_inspection": {
                "date": "2026-08-06",
                "passed": (
                    VISUALLY_INSPECTED_OVERLAY_SHA256
                    != "PENDING_ORIGINAL_RESOLUTION_VISUAL_INSPECTION"
                ),
                "overlay_sha256": VISUALLY_INSPECTED_OVERLAY_SHA256,
                "assertion": (
                    "all three crosshairs are centered on the highest "
                    "measured plus marker in their respective panels"
                ),
            },
        },
        "derived_checks": {
            "point_count": len(PEAK_MARKERS),
            "peak_measured_sample_time_ms": peak_times,
            "species_time_order": ["SiF", "SiF2", "SiF4"],
            "sif4_peak_is_near_half_millisecond": (
                0.45 <= peak_times["SiF4"] <= 0.60
            ),
        },
        "source_text_cross_checks": {
            "sif_path": "predominantly prompt collision-cascade ejection",
            "sif2_path": (
                "substantial thermally activated desorption of a "
                "collision-generated precursor"
            ),
            "sif4_path": (
                "precursor diffusion followed by desorption; the source "
                "describes a roughly 0.5 ms / hundreds-of-microseconds lag"
            ),
            "delay_depends_on_incident_energy": True,
        },
        "claim_boundary": {
            "valid": [
                (
                    "species ordering and approximate measured-peak timing "
                    "for this one beam condition"
                ),
                (
                    "a required distinction between prompt collision "
                    "products and delayed precursor transport"
                ),
            ],
            "not_valid": [
                "an escape probability",
                "a diffusion coefficient or diffusion length",
                "a first-order residence time",
                "a prompt-versus-delayed branching fraction",
                "an absolute intensity or cross-panel product yield",
                (
                    "an instrument-response-deconvolved surface clock; the "
                    "trace includes the 100 us beam pulse, flight time, "
                    "detector geometry, and product velocity distribution"
                ),
                "a C4F6-plasma or Krueger-reactor boundary parameter",
                "a target-depth calibration",
            ],
            "production_escape_parameter_use": False,
        },
        "output": {
            "path": (
                "data/experimental/karahashi_2007/"
                "figure9_cf3_product_peak_times.csv"
            ),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    csv_sha = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(
        manifest(csv_sha), indent=2, ensure_ascii=False
    ) + "\n"


def _assert_source() -> None:
    if _sha256(SOURCE_PDF) != SOURCE_PDF_SHA256:
        raise RuntimeError("Karahashi source PDF checksum does not match")


def verify_render(render_path: Path) -> Image.Image:
    if _sha256(render_path) != RENDER_SHA256:
        raise RuntimeError("600-dpi Figure 9 page render checksum does not match")
    image = Image.open(render_path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected render size {image.size}; expected {RENDER_SIZE}"
        )
    crop = image.crop(CROP_BOUNDS)
    gray = np.asarray(crop.convert("L"))
    checks = {
        "left_axis": float(np.mean(gray[296:3629, 209:214] < 96)),
        "right_axis": float(np.mean(gray[296:3629, 1238:1243] < 96)),
        "panel_a_top": float(np.mean(gray[295:300, 209:1243] < 96)),
        "panel_a_bottom": float(np.mean(gray[1184:1189, 209:1243] < 96)),
        "panel_b_top": float(np.mean(gray[1531:1536, 209:1243] < 96)),
        "panel_b_bottom": float(np.mean(gray[2421:2426, 209:1243] < 96)),
        "panel_c_top": float(np.mean(gray[2736:2741, 209:1243] < 96)),
        "panel_c_bottom": float(np.mean(gray[3625:3630, 209:1243] < 96)),
    }
    # The shared vertical axes are interrupted by panel gaps, and the raster
    # line occupies only two or three rows/columns inside each five-pixel
    # verification band.
    limits = {
        "left_axis": 0.40,
        "right_axis": 0.40,
        "panel_a_top": 0.38,
        "panel_a_bottom": 0.38,
        "panel_b_top": 0.38,
        "panel_b_bottom": 0.38,
        "panel_c_top": 0.38,
        "panel_c_bottom": 0.38,
    }
    failed = {
        name: value
        for name, value in checks.items()
        if value < limits[name]
    }
    if failed:
        raise RuntimeError(f"axis dark-pixel verification failed: {failed}")
    return crop


def verify_committed_files() -> None:
    expected_csv = csv_text()
    expected_manifest = manifest_text(expected_csv)
    if CSV_PATH.read_text(encoding="utf-8") != expected_csv:
        raise RuntimeError(f"{CSV_PATH.relative_to(ROOT)} is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != expected_manifest:
        raise RuntimeError(f"{MANIFEST_PATH.relative_to(ROOT)} is stale")


def write_committed_files() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def draw_overlay(crop: Image.Image, output: Path) -> None:
    overlay = crop.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {"SiF": "#e41a1c", "SiF2": "#377eb8", "SiF4": "#4daf4a"}
    for marker in PEAK_MARKERS:
        radius = 18
        color = colors[marker.product]
        draw.ellipse(
            (
                marker.x_px - radius,
                marker.y_px - radius,
                marker.x_px + radius,
                marker.y_px + radius,
            ),
            outline=color,
            width=4,
        )
        draw.line(
            (
                marker.x_px - 2 * radius,
                marker.y_px,
                marker.x_px + 2 * radius,
                marker.y_px,
            ),
            fill=color,
            width=3,
        )
        draw.line(
            (
                marker.x_px,
                marker.y_px - 2 * radius,
                marker.x_px,
                marker.y_px + 2 * radius,
            ),
            fill=color,
            width=3,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", type=Path, default=DEFAULT_RENDER)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    _assert_source()
    crop = verify_render(args.render)
    if args.write:
        write_committed_files()
    else:
        verify_committed_files()
    if args.overlay is not None:
        draw_overlay(crop, args.overlay)
    print(json.dumps(manifest(hashlib.sha256(
        csv_text().encode("utf-8")
    ).hexdigest()), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
