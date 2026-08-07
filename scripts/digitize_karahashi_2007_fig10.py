#!/usr/bin/env python3
"""Reproduce Karahashi 2007 Figure 10 CF3+ product fractions.

Figure 10 reports the energy-dependent composition of desorbed
silicon-fluoride products under mass-selected CF3+ bombardment of SiO2.
The ordinate is SiFx / sum(SiFx), so the result constrains product branching,
not absolute etch yield or the probability that a buried product escapes.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDF = (
    ROOT / "research_sources" / "karahashi_2007_hyomen_kagaku_28_60.pdf")
DEFAULT_RENDER = (
    ROOT / "tmp" / "pdfs" / "karahashi_figures" / "page6_600dpi.png")
OUTPUT_DIR = ROOT / "data" / "experimental" / "karahashi_2007"
CSV_PATH = OUTPUT_DIR / "figure10_cf3_product_fractions.csv"
MANIFEST_PATH = OUTPUT_DIR / "figure10_product_digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "093b18b91b0a6d910fc414779ee8320b7a046ac4cad38ef5de0b7f2dd25a2d79")
RENDER_SHA256 = (
    "e9e49375c53aa95e819302f434b68068f784219d118e2cf2fa3bb78c9c19a352")
VISUALLY_INSPECTED_OVERLAY_SHA256 = (
    "71b556bf7f71afe18c0680ae678ae5f4ce8de082e62bedca778354d964c239fe")
RENDER_SIZE = (4961, 7017)
CROP_BOUNDS = (650, 500, 2200, 1900)

# Axis centers in the exact crop. The chart extends to the unlabeled 2500 eV
# top tick; measured markers use the printed 500, 1000, and 2000 eV setpoints.
X_AT_0_EV = 355.5
X_AT_2500_EV = 1413.5
Y_AT_0_PERCENT = 1058.0
Y_AT_100_PERCENT = 144.0


@dataclass(frozen=True)
class PixelPoint:
    product: str
    energy_eV: int
    x_px: float
    y_px: float
    marker: str


PIXEL_POINTS = (
    PixelPoint("SiF", 500, 567.5, 860.0, "up_triangle"),
    PixelPoint("SiF2", 500, 567.5, 528.5, "circle"),
    PixelPoint("SiF4", 500, 567.5, 876.0, "square"),
    PixelPoint("SiF", 1000, 779.0, 756.0, "up_triangle"),
    PixelPoint("SiF2", 1000, 779.0, 629.0, "circle"),
    PixelPoint("SiF4", 1000, 779.0, 876.0, "square"),
    PixelPoint("SiF", 2000, 1202.0, 630.0, "up_triangle"),
    PixelPoint("SiF2", 2000, 1202.0, 646.5, "circle"),
    PixelPoint("SiF4", 2000, 1202.0, 993.5, "square"),
)

FIELDNAMES = (
    "incident_species",
    "target",
    "energy_eV",
    "ion_incidence_angle_deg",
    "ion_incidence_angle_status",
    "desorbed_product",
    "product_fraction_percent_of_detected_sifx",
    "marker_center_x_crop_px",
    "marker_center_y_crop_px",
    "marker",
    "digitization_uncertainty_percentage_points",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _energy_at_pixel(x_px: float) -> float:
    return 2500.0 * (x_px - X_AT_0_EV) / (
        X_AT_2500_EV - X_AT_0_EV)


def _percent_at_pixel(y_px: float) -> float:
    return 100.0 * (Y_AT_0_PERCENT - y_px) / (
        Y_AT_0_PERCENT - Y_AT_100_PERCENT)


def rows() -> list[dict[str, str]]:
    return [
        {
            "incident_species": "CF3+",
            "target": "SiO2",
            "energy_eV": str(point.energy_eV),
            "ion_incidence_angle_deg": "",
            "ion_incidence_angle_status": "unreported_in_source",
            "desorbed_product": point.product,
            "product_fraction_percent_of_detected_sifx": (
                f"{_percent_at_pixel(point.y_px):.3f}"),
            "marker_center_x_crop_px": f"{point.x_px:.1f}",
            "marker_center_y_crop_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            # 6.5 pixels covers line/marker overlap at 500 and 2000 eV.
            "digitization_uncertainty_percentage_points": "0.72",
        }
        for point in PIXEL_POINTS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    by_energy: dict[int, float] = {}
    for row in rows():
        energy = int(row["energy_eV"])
        by_energy[energy] = by_energy.get(energy, 0.0) + float(
            row["product_fraction_percent_of_detected_sifx"])
    x_offsets = [
        abs(_energy_at_pixel(point.x_px) - point.energy_eV)
        for point in PIXEL_POINTS
    ]
    return {
        "manifest_id": "KARAHASHI-2007-FIG10-CF3-PRODUCTS-R2",
        "source": {
            "citation": (
                "K. Karahashi, Hyomen Kagaku 28, 60-66 (2007), Figure 10"),
            "local_pdf": (
                "research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf"),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 6,
            "print_page": 65,
            "figure": "Figure 10",
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
            "energies_eV": [500, 1000, 2000],
            "ion_incidence_angle_deg": None,
            "ion_incidence_angle_status": "unreported_in_source",
            "detected_products": ["SiF", "SiF2", "SiF4"],
            "quantity": "SiFx / sum(SiFx), percent",
            "absolute_product_yield_reported": False,
        },
        "pixel_calibration": {
            "crop_bounds_full_render_px": list(CROP_BOUNDS),
            "x_at_0_eV": X_AT_0_EV,
            "x_at_2500_eV": X_AT_2500_EV,
            "y_at_zero_percent": Y_AT_0_PERCENT,
            "y_at_100_percent": Y_AT_100_PERCENT,
            "transform": {
                "energy_eV": (
                    "2500 * (x_px - 355.5) / (1413.5 - 355.5)"),
                "percent": (
                    "100 * (1058.0 - y_px) / (1058.0 - 144.0)"),
            },
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; exact PIL crop; dark-pixel axis "
                "localization; marker-center and connected-contour "
                "transcription; full-resolution visual overlay"
            ),
            "uncertainty_percentage_points": 0.72,
            "uncertainty_basis": (
                "6.5-pixel allowance covering line/marker overlap, rounded "
                "upward through the linear ordinate transform"
            ),
            "nominal_energy_policy": (
                "use printed beam-energy setpoints; retain marker x pixels "
                "as an independent placement check"
            ),
            "full_resolution_visual_inspection": {
                "date": "2026-08-06",
                "passed": True,
                "overlay_sha256": VISUALLY_INSPECTED_OVERLAY_SHA256,
                "assertion": (
                    "all nine color-coded crosshairs and rings are centered "
                    "on their source markers, including the overlapping "
                    "500 eV and 2000 eV pairs"
                ),
            },
        },
        "derived_checks": {
            "point_count": len(PIXEL_POINTS),
            "maximum_marker_x_energy_offset_eV": round(max(x_offsets), 3),
            "fraction_sum_percent_by_energy": {
                str(energy): round(total, 3)
                for energy, total in sorted(by_energy.items())
            },
            "dominant_product_at_or_below_1000_eV": "SiF2",
            "sif_fraction_increases_with_energy": True,
            "sif4_fraction_decreases_with_energy": True,
        },
        "source_text_cross_checks": {
            "figure10_ion_incidence_angle": "unreported",
            "main_product": "SiF2",
            "below_or_equal_1000_eV_dominant_path": (
                "thermally activated desorption of a collision-cascade-"
                "formed precursor"),
            "sif_path": "prompt physical collision-cascade ejection",
            "sif4_delay_order_ms": 0.5,
            "co_detected_as_oxygen_containing_product": True,
        },
        "claim_boundary": {
            "valid": [
                (
                    "energy-resolved branching among detected SiF, SiF2, "
                    "and SiF4 products for CF3+ on SiO2"),
                (
                    "a product-identity constraint for a condition-matched "
                    "reactive-ion event model"),
            ],
            "not_valid": [
                "an absolute etch yield",
                "a product escape probability or diffusion length",
                "a prompt-versus-delayed numerical partition",
                "an accounting of CO in the SiFx-normalized ordinate",
                "a universal branch law for CF+, CF2+, or plasma mixtures",
                (
                    "a condition-matched join to the normal-incidence "
                    "Figure-4 etch yields"
                ),
                "a target-depth calibration",
            ],
            "production_escape_parameter_use": False,
        },
        "output": {
            "path": (
                "data/experimental/karahashi_2007/"
                "figure10_cf3_product_fractions.csv"),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    csv_sha = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(
        manifest(csv_sha), indent=2, ensure_ascii=False) + "\n"


def _assert_source() -> None:
    if _sha256(SOURCE_PDF) != SOURCE_PDF_SHA256:
        raise RuntimeError("Karahashi source PDF checksum does not match")


def verify_render(render_path: Path) -> Image.Image:
    if _sha256(render_path) != RENDER_SHA256:
        raise RuntimeError("600-dpi Figure 10 page render checksum does not match")
    image = Image.open(render_path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected render size {image.size}; expected {RENDER_SIZE}")
    crop = image.crop(CROP_BOUNDS)
    gray = np.asarray(crop.convert("L"))
    checks = {
        "left_axis": float(np.mean(gray[142:1061, 354:358] < 96)),
        "right_axis": float(np.mean(gray[142:1061, 1412:1416] < 96)),
        "top_axis": float(np.mean(gray[142:147, 354:1416] < 96)),
        "bottom_axis": float(np.mean(gray[1056:1061, 354:1416] < 96)),
    }
    failed = {name: value for name, value in checks.items() if value < 0.70}
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
    for point in PIXEL_POINTS:
        radius = 12
        color = colors[point.product]
        draw.line(
            (point.x_px - radius, point.y_px,
             point.x_px + radius, point.y_px),
            fill=color,
            width=3,
        )
        draw.line(
            (point.x_px, point.y_px - radius,
             point.x_px, point.y_px + radius),
            fill=color,
            width=3,
        )
        draw.ellipse(
            (point.x_px - 18, point.y_px - 18,
             point.x_px + 18, point.y_px + 18),
            outline=color,
            width=3,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", type=Path, default=DEFAULT_RENDER)
    parser.add_argument("--overlay", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    _assert_source()
    crop = verify_render(arguments.render)
    if arguments.write:
        write_committed_files()
    elif arguments.check or CSV_PATH.exists() or MANIFEST_PATH.exists():
        verify_committed_files()
    if arguments.overlay:
        draw_overlay(crop, arguments.overlay)
    print(manifest_text(csv_text()), end="")


if __name__ == "__main__":
    main()
