#!/usr/bin/env python3
"""Audit and reproduce the Karahashi 2007 Figure 4 beam-yield digitization.

The source figure reports SiO2 etching yield for mass-selected F+, CF+, CF2+,
and CF3+ beams with no incident neutral-radical flux.  This script keeps the
source pixels, axis transform, marker centers, and plotted error-bar caps in
one replayable record.  It can:

* regenerate the committed CSV and manifest from the pixel record;
* verify those files and the source checksums without changing them; and
* draw a high-resolution QA overlay for visual comparison with the source.

The nominal energy values are the plotted experimental setpoints.  The pixel
x coordinates are retained as an independent transcription check; the energy
column is not inferred from small marker-placement offsets.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDF = ROOT / "research_sources" / "karahashi_2007_hyomen_kagaku_28_60.pdf"
DEFAULT_RENDER = ROOT / "tmp" / "pdfs" / "karahashi_figures" / "page4_600dpi.png"
OUTPUT_DIR = ROOT / "data" / "experimental" / "karahashi_2007"
CSV_PATH = OUTPUT_DIR / "figure4_reactive_ion_yields.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"

SOURCE_PDF_SHA256 = "093b18b91b0a6d910fc414779ee8320b7a046ac4cad38ef5de0b7f2dd25a2d79"
RENDER_SHA256 = "075cd97e5d7b5238f028f546d5c63cfc0c96cd6ff33a769b640f130d73b5fe1d"
RENDER_SIZE = (4961, 7017)

# Axis centers in the full 600-dpi page render.  Each line is 3--5 pixels
# thick; these centers were obtained from the repeated dark-pixel maxima.
X_AT_0_EV = 973.0
X_AT_2500_EV = 2080.0
Y_AT_0_YIELD = 1524.0
Y_AT_2P5_YIELD = 650.0

# The source plots error bars but does not identify them as a statistical
# confidence interval in the accessible text.  Keep them as "plotted" bounds.
# Marker centers and cap centers are in full-render pixels.
@dataclass(frozen=True)
class PixelPoint:
    species: str
    energy_eV: int
    x_px: float
    y_px: float
    upper_cap_y_px: float
    lower_cap_y_px: float
    marker: str


PIXEL_POINTS = (
    PixelPoint("F+", 500, 1194.0, 1444.5, 1428.0, 1461.5, "down_triangle"),
    PixelPoint("F+", 1000, 1419.0, 1411.0, 1393.0, 1430.0, "down_triangle"),
    PixelPoint("F+", 1500, 1640.0, 1408.5, 1387.0, 1430.0, "down_triangle"),
    PixelPoint("F+", 2000, 1862.5, 1371.5, 1346.5, 1396.5, "down_triangle"),
    PixelPoint("CF+", 750, 1304.5, 1349.0, 1326.0, 1373.0, "up_triangle"),
    PixelPoint("CF+", 1000, 1415.0, 1288.0, 1264.0, 1330.0, "up_triangle"),
    PixelPoint("CF+", 1250, 1529.5, 1288.0, 1260.0, 1323.0, "up_triangle"),
    PixelPoint("CF+", 1500, 1639.0, 1284.5, 1236.0, 1333.0, "up_triangle"),
    PixelPoint("CF+", 2000, 1861.5, 1261.0, 1226.0, 1297.0, "up_triangle"),
    PixelPoint("CF2+", 500, 1195.0, 1335.5, 1311.5, 1360.0, "square"),
    PixelPoint("CF2+", 750, 1305.5, 1185.5, 1158.0, 1233.0, "square"),
    PixelPoint("CF2+", 1000, 1416.0, 1106.0, 1059.5, 1152.0, "square"),
    PixelPoint("CF2+", 1500, 1638.5, 1088.5, 1038.0, 1140.0, "square"),
    PixelPoint("CF2+", 2000, 1858.0, 1076.5, 1022.0, 1130.5, "square"),
    PixelPoint("CF3+", 250, 1083.0, 1342.0, 1319.0, 1367.0, "circle"),
    PixelPoint("CF3+", 500, 1194.5, 1222.0, 1185.0, 1259.0, "circle"),
    PixelPoint("CF3+", 750, 1305.0, 1110.5, 1061.5, 1158.0, "circle"),
    PixelPoint("CF3+", 1000, 1416.0, 1010.0, 958.0, 1059.5, "circle"),
    PixelPoint("CF3+", 1250, 1527.5, 970.0, 908.5, 1029.0, "circle"),
    PixelPoint("CF3+", 1500, 1636.5, 869.0, 800.0, 931.0, "circle"),
    PixelPoint("CF3+", 2000, 1860.0, 910.5, 846.0, 973.0, "circle"),
)

FIELDNAMES = (
    "species",
    "energy_eV",
    "yield_sio2_per_ion",
    "plotted_lower_yield",
    "plotted_upper_yield",
    "marker_center_x_px",
    "marker_center_y_px",
    "upper_cap_y_px",
    "lower_cap_y_px",
    "marker",
    "digitization_yield_uncertainty",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _yield_at_pixel(y_px: float) -> float:
    return 2.5 * (Y_AT_0_YIELD - y_px) / (Y_AT_0_YIELD - Y_AT_2P5_YIELD)


def _energy_at_pixel(x_px: float) -> float:
    return 2500.0 * (x_px - X_AT_0_EV) / (X_AT_2500_EV - X_AT_0_EV)


def rows() -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    for point in PIXEL_POINTS:
        center = _yield_at_pixel(point.y_px)
        lower = _yield_at_pixel(point.lower_cap_y_px)
        upper = _yield_at_pixel(point.upper_cap_y_px)
        result.append(
            {
                "species": point.species,
                "energy_eV": str(point.energy_eV),
                "yield_sio2_per_ion": f"{center:.4f}",
                "plotted_lower_yield": f"{lower:.4f}",
                "plotted_upper_yield": f"{upper:.4f}",
                "marker_center_x_px": f"{point.x_px:.1f}",
                "marker_center_y_px": f"{point.y_px:.1f}",
                "upper_cap_y_px": f"{point.upper_cap_y_px:.1f}",
                "lower_cap_y_px": f"{point.lower_cap_y_px:.1f}",
                "marker": point.marker,
                # 3.5 pixels covers line thickness plus center-placement
                # ambiguity.  Rounded upward from 3.5 * 2.5 / 874.
                "digitization_yield_uncertainty": "0.011",
            }
        )
    return result


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    x_offsets = [
        abs(_energy_at_pixel(point.x_px) - point.energy_eV)
        for point in PIXEL_POINTS
    ]
    return {
        "manifest_id": "KARAHASHI-2007-FIG4-REACTIVE-ION-YIELDS-R1",
        "source": {
            "citation": (
                "K. Karahashi, Surface Science Society of Japan 28, 60-65 (2007); "
                "open review of the mass-selected beam experiments reported in "
                "JVST A 22, 1166 (2004), DOI 10.1116/1.1761119"
            ),
            "local_pdf": "research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 4,
            "print_page": 63,
            "figure": "Figure 4",
            "render_dpi": 600,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
        },
        "experiment_scope": {
            "target": "SiO2",
            "incident_ions": ["F+", "CF+", "CF2+", "CF3+"],
            "neutral_radical_flux": "none",
            "beam_selection": "mass selected",
            "reported_energy_domain_eV": [250, 2000],
            "quantity": "etched SiO2 units per incident ion",
        },
        "pixel_calibration": {
            "x_at_0_eV": X_AT_0_EV,
            "x_at_2500_eV": X_AT_2500_EV,
            "y_at_zero_yield": Y_AT_0_YIELD,
            "y_at_2p5_yield": Y_AT_2P5_YIELD,
            "transform": {
                "energy_eV": "2500 * (x_px - 973.0) / (2080.0 - 973.0)",
                "yield": "2.5 * (1524.0 - y_px) / (1524.0 - 650.0)",
            },
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL dark-pixel axis localization; "
                "marker-center and plotted-cap transcription; connected-contour "
                "and full-resolution visual-overlay reconciliation"
            ),
            "nominal_energy_policy": (
                "use plotted experimental setpoints; retain marker x pixels only "
                "as an independent placement check"
            ),
            "yield_digitization_uncertainty": 0.011,
            "uncertainty_basis": (
                "3.5-pixel vertical center-placement allowance, rounded upward "
                "after the linear y-axis transform"
            ),
            "plotted_bounds_policy": (
                "transcribe the source error-bar caps but do not call them a "
                "confidence interval because the accessible source text does not "
                "define their statistical meaning"
            ),
        },
        "derived_checks": {
            "point_count": len(PIXEL_POINTS),
            "maximum_marker_x_energy_offset_eV": round(max(x_offsets), 3),
            "cf3_yield_at_1000_eV": float(
                next(
                    row["yield_sio2_per_ion"]
                    for row in rows()
                    if row["species"] == "CF3+" and row["energy_eV"] == "1000"
                )
            ),
            "cf3_maximum_digitized_yield": max(
                float(row["yield_sio2_per_ion"])
                for row in rows()
                if row["species"] == "CF3+"
            ),
            "cf3_maximum_energy_eV": 1500,
        },
        "claim_boundary": {
            "valid": (
                "direct radical-free single-species beam benchmark over the "
                "digitized species/energy support"
            ),
            "not_valid": [
                "a universal yield ceiling",
                "a fluorocarbon-plasma total-yield ceiling",
                "an ion-mixture law for the Krueger reactor",
                "an extrapolation above 2000 eV",
            ],
        },
        "output": {
            "path": "data/experimental/karahashi_2007/figure4_reactive_ion_yields.csv",
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    csv_sha = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(csv_sha), indent=2, ensure_ascii=False) + "\n"


def _assert_source() -> None:
    if _sha256(SOURCE_PDF) != SOURCE_PDF_SHA256:
        raise RuntimeError("Karahashi source PDF checksum does not match the audited source")


def _axis_dark_fraction(gray: np.ndarray, *, horizontal: bool, coordinate: int) -> float:
    if horizontal:
        values = gray[coordinate, int(X_AT_0_EV) : int(X_AT_2500_EV) + 1]
    else:
        values = gray[int(Y_AT_2P5_YIELD) : int(Y_AT_0_YIELD) + 1, coordinate]
    return float(np.mean(values < 96))


def verify_render(render_path: Path) -> Image.Image:
    if _sha256(render_path) != RENDER_SHA256:
        raise RuntimeError("600-dpi Figure 4 page render checksum does not match")
    image = Image.open(render_path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected render size {image.size}; expected {RENDER_SIZE}")
    gray = np.asarray(image.convert("L"))
    checks = {
        "left axis": _axis_dark_fraction(gray, horizontal=False, coordinate=973),
        "right axis": _axis_dark_fraction(gray, horizontal=False, coordinate=2080),
        "top axis": _axis_dark_fraction(gray, horizontal=True, coordinate=650),
        "bottom axis": _axis_dark_fraction(gray, horizontal=True, coordinate=1524),
    }
    failed = {name: value for name, value in checks.items() if value < 0.75}
    if failed:
        raise RuntimeError(f"axis dark-pixel verification failed: {failed}")
    return image


def verify_committed_files() -> None:
    expected_csv = csv_text()
    expected_manifest = manifest_text(expected_csv)
    if CSV_PATH.read_text(encoding="utf-8") != expected_csv:
        raise RuntimeError(f"{CSV_PATH.relative_to(ROOT)} is not reproduced by this script")
    if MANIFEST_PATH.read_text(encoding="utf-8") != expected_manifest:
        raise RuntimeError(f"{MANIFEST_PATH.relative_to(ROOT)} is not reproduced by this script")


def write_committed_files() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def _draw_cross(draw: ImageDraw.ImageDraw, x: float, y: float, color: str) -> None:
    radius = 12
    draw.line((x - radius, y, x + radius, y), fill=color, width=3)
    draw.line((x, y - radius, x, y + radius), fill=color, width=3)
    draw.ellipse((x - 18, y - 18, x + 18, y + 18), outline=color, width=3)


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {"F+": "#e41a1c", "CF+": "#377eb8", "CF2+": "#4daf4a", "CF3+": "#984ea3"}
    for point in PIXEL_POINTS:
        color = colors[point.species]
        _draw_cross(draw, point.x_px, point.y_px, color)
        draw.line(
            (point.x_px + 25, point.upper_cap_y_px, point.x_px + 25, point.lower_cap_y_px),
            fill=color,
            width=3,
        )
        draw.line(
            (point.x_px + 18, point.upper_cap_y_px, point.x_px + 32, point.upper_cap_y_px),
            fill=color,
            width=3,
        )
        draw.line(
            (point.x_px + 18, point.lower_cap_y_px, point.x_px + 32, point.lower_cap_y_px),
            fill=color,
            width=3,
        )
    crop = overlay.crop((900, 600, 2150, 1580))
    output.parent.mkdir(parents=True, exist_ok=True)
    crop.save(output)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--render",
        type=Path,
        default=DEFAULT_RENDER,
        help="600-dpi PNG of source PDF page 4",
    )
    parser.add_argument(
        "--overlay",
        type=Path,
        help="optional output path for a visually inspectable QA crop",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="regenerate the committed CSV and manifest (default: verify only)",
    )
    return parser


def main() -> None:
    args = _parser().parse_args()
    _assert_source()
    image = verify_render(args.render)
    if args.write:
        write_committed_files()
    verify_committed_files()
    if args.overlay is not None:
        draw_overlay(image, args.overlay)
    print(
        json.dumps(
            {
                "status": "verified",
                "points": len(PIXEL_POINTS),
                "cf3_1000_eV": manifest(
                    hashlib.sha256(csv_text().encode()).hexdigest()
                )["derived_checks"]["cf3_yield_at_1000_eV"],
                "cf3_max": manifest(
                    hashlib.sha256(csv_text().encode()).hexdigest()
                )["derived_checks"]["cf3_maximum_digitized_yield"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
