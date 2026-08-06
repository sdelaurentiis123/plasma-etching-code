#!/usr/bin/env python3
"""Audit and reproduce Takada et al. (2005) Figure 3 SiO2 yields.

Figure 3 reports radical-free co-incidence of a normal Ar+ beam and either a
45-degree C5F8 molecular beam or a CF2 radical beam.  The x axis is logarithmic.
This script retains the full-page 600-dpi marker pixels, replays the axis
transform, verifies the archived source/render checksums, and produces a visual
overlay.  It is deliberately an analog surface-physics record: C5F8 is not
C4F6, and these data are not a boundary condition for the Krueger reactor.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDF = ROOT / "research_sources" / "takada_2005_radical_free_etching.pdf"
DEFAULT_RENDER = (
    ROOT / "tmp" / "pdfs" / "takada_2005" / "page3_600dpi.png")
OUTPUT_DIR = ROOT / "data" / "experimental" / "takada_2005"
CSV_PATH = OUTPUT_DIR / "figure3_sio2_coincidence_yields.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "9034f445f575b85c0b9c95d79b81699c42eafa47a347a7e73a6de73cbe222e25")
RENDER_SHA256 = (
    "dc1b45e63764c6fa562a3b6b1f7d652977eb96f07650668125372543cf014d83")
RENDER_SIZE = (4912, 7132)

# Full-page pixel centers of the scanned Figure-3 axes.  Repeated dark-pixel
# maxima put the vertical borders at x=844 and x=2084 and the 0.0 and 1.4
# major ticks at y=1740.5 and y=989.0, respectively.
X_AT_R_0P1 = 844.0
X_AT_R_10 = 2084.0
Y_AT_0_YIELD = 1740.5
Y_AT_1P4_YIELD = 989.0


@dataclass(frozen=True)
class PixelPoint:
    co_incident_species: str
    flux_ratio: float
    x_px: float
    y_px: float
    marker: str


# SiO2 series only.  The source also plots Si yields; they are not consumed by
# the oxide-depth question and are excluded rather than partially transcribed.
PIXEL_POINTS = (
    PixelPoint("C5F8", 0.25, 1088.0, 1381.0, "filled_circle"),
    PixelPoint("C5F8", 0.50, 1272.0, 1195.0, "filled_circle"),
    PixelPoint("C5F8", 1.00, 1460.0, 1098.0, "filled_circle"),
    PixelPoint("C5F8", 2.50, 1708.0, 1172.0, "filled_circle"),
    PixelPoint("C5F8", 10.0, 2084.0, 1314.0, "filled_circle"),
    PixelPoint("CF2", 0.25, 1088.0, 1498.0, "filled_square"),
    PixelPoint("CF2", 1.00, 1460.0, 1395.0, "filled_square"),
    PixelPoint("CF2", 2.50, 1708.0, 1472.0, "filled_square"),
)

FIELDNAMES = (
    "coincident_species",
    "flux_ratio_to_ar_ion",
    "ar_ion_energy_eV",
    "yield_sio2_per_ar_ion",
    "marker_center_x_px",
    "marker_center_y_px",
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
    return 1.4 * (
        Y_AT_0_YIELD - y_px) / (Y_AT_0_YIELD - Y_AT_1P4_YIELD)


def _ratio_at_pixel(x_px: float) -> float:
    log10_ratio = -1.0 + 2.0 * (
        x_px - X_AT_R_0P1) / (X_AT_R_10 - X_AT_R_0P1)
    return 10.0 ** log10_ratio


def rows() -> list[dict[str, str]]:
    return [
        {
            "coincident_species": point.co_incident_species,
            "flux_ratio_to_ar_ion": f"{point.flux_ratio:g}",
            "ar_ion_energy_eV": "400",
            "yield_sio2_per_ar_ion": f"{_yield_at_pixel(point.y_px):.4f}",
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            # Eight vertical pixels is wider than the scanned marker edge and
            # axis-center ambiguity: 8 * 1.4 / 751.5 = 0.0149.
            "digitization_yield_uncertainty": "0.015",
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
    ratio_offsets = [
        abs(math.log10(_ratio_at_pixel(point.x_px) / point.flux_ratio))
        for point in PIXEL_POINTS
    ]
    c5f8_rows = [
        row for row in rows() if row["coincident_species"] == "C5F8"]
    return {
        "manifest_id": "TAKADA-2005-FIG3-SIO2-COINCIDENCE-YIELDS-R1",
        "source": {
            "citation": (
                "N. Takada, H. Toyoda, and H. Sugai, Transactions of the "
                "Materials Research Society of Japan 30[1], 319-322 (2005), "
                "Evidence of Radical-free Etching of SiO2 by Fluorocarbon "
                "Molecule under Ion Bombardment"
            ),
            "related_journal_article": (
                "N. Takada et al., Journal of Applied Physics 97, 013534 "
                "(2005), DOI 10.1063/1.1829400"
            ),
            "url": (
                "https://www.mrs-j.org/pub/tmrsj/vol30_no1/"
                "vol30_no1_319.pdf"
            ),
            "local_pdf": (
                "research_sources/takada_2005_radical_free_etching.pdf"
            ),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 3,
            "print_page": 321,
            "figure": "Figure 3",
            "render_dpi": 600,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
        },
        "experiment_scope": {
            "target": "thermally grown 20 nm SiO2 on Si",
            "ion": "mass-selected Ar+",
            "ion_energy_eV": 400,
            "ion_incidence": "parallel to surface normal",
            "coincident_species": {
                "C5F8": "stable molecule, 45-degree beam",
                "CF2": "radical from thermal HFPO decomposition",
            },
            "gas_phase_radicals": "absent for the C5F8 series",
            "quantity": "removed SiO2 units per incident Ar+ ion",
        },
        "pixel_calibration": {
            "x_at_flux_ratio_0p1": X_AT_R_0P1,
            "x_at_flux_ratio_10": X_AT_R_10,
            "y_at_zero_yield": Y_AT_0_YIELD,
            "y_at_1p4_yield": Y_AT_1P4_YIELD,
            "x_scale": "log10",
            "transform": {
                "flux_ratio": (
                    "10 ** (-1 + 2 * (x_px - 844.0) / (2084.0 - 844.0))"
                ),
                "yield": (
                    "1.4 * (1740.5 - y_px) / (1740.5 - 989.0)"
                ),
            },
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy dark-pixel axis "
                "localization; full-resolution marker-center transcription; "
                "visual-overlay reconciliation"
            ),
            "nominal_flux_ratio_policy": (
                "use plotted experimental setpoints; retain marker x pixels "
                "as an independent log-axis placement check"
            ),
            "yield_digitization_uncertainty": 0.015,
            "uncertainty_basis": (
                "eight-pixel vertical allowance, wider than the scanned "
                "marker-edge and axis-center ambiguity"
            ),
        },
        "text_cross_checks": {
            "c5f8_yield_at_ratio_1_rounded": 1.2,
            "c5f8_yield_at_ratio_0p25_rounded": 0.67,
            "c5f8_yield_at_900eV_ratio_1_not_in_figure3": 2.5,
        },
        "derived_checks": {
            "point_count": len(PIXEL_POINTS),
            "maximum_abs_log10_ratio_pixel_offset": round(
                max(ratio_offsets), 5),
            "c5f8_maximum_400eV_yield": max(
                float(row["yield_sio2_per_ar_ion"])
                for row in c5f8_rows
            ),
            "c5f8_maximum_400eV_flux_ratio": 1.0,
        },
        "claim_boundary": {
            "valid": (
                "direct radical-free evidence that a stable fluorocarbon "
                "molecule can enhance ion-assisted SiO2 removal and can also "
                "drive a non-monotone deposition/etch balance"
            ),
            "not_valid": [
                "a C4F6 surface-yield law",
                "a boundary flux for the Krueger reactor",
                "a high-aspect-ratio transport measurement",
                "an extrapolation in ion energy or molecular identity",
            ],
        },
        "output": {
            "path": (
                "data/experimental/takada_2005/"
                "figure3_sio2_coincidence_yields.csv"
            ),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    csv_sha = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(csv_sha), indent=2) + "\n"


def _assert_source() -> None:
    if _sha256(SOURCE_PDF) != SOURCE_PDF_SHA256:
        raise RuntimeError("Takada source PDF checksum does not match")


def verify_render(render_path: Path) -> Image.Image:
    if _sha256(render_path) != RENDER_SHA256:
        raise RuntimeError("600-dpi Takada page-3 render checksum does not match")
    image = Image.open(render_path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected render size {image.size}; expected {RENDER_SIZE}")
    gray = np.asarray(image.convert("L"))
    checks = {
        "left axis": np.mean(
            np.min(gray[930:1851, 841:848], axis=1) < 96),
        "right axis": np.mean(
            np.min(gray[930:1851, 2081:2088], axis=1) < 96),
        "top axis": np.mean(
            np.min(gray[932:939, 844:2085], axis=0) < 96),
        "bottom axis": np.mean(
            np.min(gray[1844:1851, 844:2085], axis=0) < 96),
    }
    failed = {
        name: float(value) for name, value in checks.items() if value < 0.70}
    if failed:
        raise RuntimeError(f"axis dark-pixel verification failed: {failed}")
    return image


def write_committed_files() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Takada CSV is not reproduced by script")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Takada manifest is not reproduced by script")


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {"C5F8": "#e41a1c", "CF2": "#377eb8"}
    for point in PIXEL_POINTS:
        color = colors[point.co_incident_species]
        radius = 18
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline=color,
            width=4,
        )
        draw.line(
            (point.x_px - 12, point.y_px, point.x_px + 12, point.y_px),
            fill=color,
            width=3,
        )
        draw.line(
            (point.x_px, point.y_px - 12, point.x_px, point.y_px + 12),
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
    image = verify_render(args.render)
    if args.write:
        write_committed_files()
    verify_committed_files()
    if args.overlay is not None:
        draw_overlay(image, args.overlay)

    c5f8 = {
        float(row["flux_ratio_to_ar_ion"]): float(
            row["yield_sio2_per_ar_ion"])
        for row in rows()
        if row["coincident_species"] == "C5F8"
    }
    print(json.dumps({
        "status": "verified",
        "points": len(PIXEL_POINTS),
        "c5f8_ratio_0p25": c5f8[0.25],
        "c5f8_ratio_1": c5f8[1.0],
        "c5f8_ratio_10": c5f8[10.0],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
