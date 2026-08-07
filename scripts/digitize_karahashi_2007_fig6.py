#!/usr/bin/env python3
"""Audit and reproduce Karahashi 2007 Figure 6 angular yields.

Figure 6 reports mass-selected CF+, CF2+, CF3+, and Ar+ bombardment of
SiO2 versus incidence angle. The source does not report the ion energy.
Arts et al. (Applied Physics Reviews 2021) reproduce this exact panel and
explicitly label the kinetic energy unknown. This script therefore preserves
the useful marker geometry while enforcing a non-production claim boundary.

The normal-incidence points independently match Karahashi Figure 4 at
1000 eV. That is recorded as a cross-figure inference, not silently promoted
to a reported condition.
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
    ROOT / "research_sources" / "karahashi_2007_hyomen_kagaku_28_60.pdf")
DEFAULT_RENDER = (
    ROOT / "tmp" / "pdfs" / "karahashi_figures" / "page4_600dpi.png")
OUTPUT_DIR = ROOT / "data" / "experimental" / "karahashi_2007"
CSV_PATH = OUTPUT_DIR / "figure6_angular_yields_unknown_energy.csv"
MANIFEST_PATH = OUTPUT_DIR / "figure6_angular_digitization_manifest.json"
FIGURE4_CSV_PATH = OUTPUT_DIR / "figure4_reactive_ion_yields.csv"

SOURCE_PDF_SHA256 = (
    "093b18b91b0a6d910fc414779ee8320b7a046ac4cad38ef5de0b7f2dd25a2d79")
RENDER_SHA256 = (
    "075cd97e5d7b5238f028f546d5c63cfc0c96cd6ff33a769b640f130d73b5fe1d")
RENDER_SIZE = (4961, 7017)
CROP_BOUNDS = (2700, 550, 4750, 1650)

# Axis centers in the exact crop above. The 0, 30, 60, and 90 degree tick
# centers give the x transform. The 0, 1, and 2 yield tick/axis centers give
# the y transform. Text glyphs in this PDF are baseline-shifted relative to
# the ticks, so calibration uses the ruled tick pixels, not glyph centers.
X_AT_0_DEG = 348.5
X_AT_90_DEG = 1422.0
Y_AT_0_YIELD = 993.5
Y_AT_2_YIELD = 206.0


@dataclass(frozen=True)
class PixelPoint:
    species: str
    angle_deg: int
    x_px: float
    y_px: float
    marker: str


PIXEL_POINTS = (
    PixelPoint("CF+", 0, 348.5, 729.5, "diamond"),
    PixelPoint("CF+", 30, 712.0, 662.5, "diamond"),
    PixelPoint("CF+", 45, 894.5, 547.0, "diamond"),
    PixelPoint("CF+", 60, 1059.0, 399.0, "diamond"),
    PixelPoint("CF+", 75, 1241.0, 582.5, "diamond"),
    PixelPoint("CF2+", 0, 349.0, 533.0, "up_triangle"),
    PixelPoint("CF2+", 30, 711.5, 512.0, "up_triangle"),
    PixelPoint("CF2+", 45, 894.0, 479.5, "up_triangle"),
    PixelPoint("CF2+", 60, 1059.0, 383.0, "up_triangle"),
    PixelPoint("CF2+", 75, 1241.0, 629.5, "up_triangle"),
    PixelPoint("CF3+", 0, 347.5, 414.5, "circle"),
    PixelPoint("CF3+", 30, 712.0, 349.0, "circle"),
    PixelPoint("CF3+", 45, 878.0, 366.0, "circle"),
    PixelPoint("CF3+", 60, 1059.0, 249.5, "circle"),
    PixelPoint("CF3+", 75, 1240.5, 546.5, "circle"),
    PixelPoint("Ar+", 0, 348.0, 793.0, "square"),
    PixelPoint("Ar+", 30, 712.0, 680.0, "square"),
    PixelPoint("Ar+", 50, 941.5, 594.0, "square"),
    PixelPoint("Ar+", 60, 1059.0, 264.0, "square"),
    PixelPoint("Ar+", 75, 1239.0, 249.5, "square"),
)

FIELDNAMES = (
    "species",
    "angle_deg_from_normal",
    "yield_sio2_per_ion",
    "marker_center_x_crop_px",
    "marker_center_y_crop_px",
    "marker",
    "digitization_yield_uncertainty",
    "ion_energy_eV",
    "ion_energy_status",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _yield_at_pixel(y_px: float) -> float:
    return 2.0 * (Y_AT_0_YIELD - y_px) / (
        Y_AT_0_YIELD - Y_AT_2_YIELD)


def _angle_at_pixel(x_px: float) -> float:
    return 90.0 * (x_px - X_AT_0_DEG) / (
        X_AT_90_DEG - X_AT_0_DEG)


def rows() -> list[dict[str, str]]:
    return [
        {
            "species": point.species,
            "angle_deg_from_normal": str(point.angle_deg),
            "yield_sio2_per_ion": f"{_yield_at_pixel(point.y_px):.4f}",
            "marker_center_x_crop_px": f"{point.x_px:.1f}",
            "marker_center_y_crop_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            # 5.5 pixels allows for line overlap at the crowded 60 degree
            # markers and rounds upward through the linear y transform.
            "digitization_yield_uncertainty": "0.015",
            "ion_energy_eV": "",
            "ion_energy_status": "not_reported_by_source",
        }
        for point in PIXEL_POINTS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _figure4_values() -> dict[tuple[str, int], float]:
    with FIGURE4_CSV_PATH.open(newline="", encoding="utf-8") as handle:
        return {
            (row["species"], int(row["energy_eV"])): float(
                row["yield_sio2_per_ion"])
            for row in csv.DictReader(handle)
        }


def _normal_values() -> dict[str, float]:
    return {
        row["species"]: float(row["yield_sio2_per_ion"])
        for row in rows()
        if row["angle_deg_from_normal"] == "0"
        and row["species"] != "Ar+"
    }


def _candidate_energy_rmse() -> dict[str, float]:
    normal = _normal_values()
    figure4 = _figure4_values()
    candidates = (750, 1000, 1500, 2000)
    result = {}
    for energy in candidates:
        errors = [
            normal[species] - figure4[(species, energy)]
            for species in ("CF+", "CF2+", "CF3+")
        ]
        result[str(energy)] = float(
            np.sqrt(np.mean(np.square(errors))))
    return result


def _ratio(species: str, angle: int) -> float:
    by_key = {
        (row["species"], int(row["angle_deg_from_normal"])): float(
            row["yield_sio2_per_ion"])
        for row in rows()
    }
    return by_key[(species, angle)] / by_key[(species, 0)]


def manifest(csv_sha256: str) -> dict[str, object]:
    x_offsets = [
        abs(_angle_at_pixel(point.x_px) - point.angle_deg)
        for point in PIXEL_POINTS
    ]
    candidate_rmse = _candidate_energy_rmse()
    return {
        "manifest_id": "KARAHASHI-2007-FIG6-ANGULAR-YIELDS-R1",
        "source": {
            "citation": (
                "K. Karahashi, Hyomen Kagaku 28, 60-66 (2007), Figure 6; "
                "data attributed there to K.-i. Yanai et al., Journal of "
                "Applied Physics 97, 053302 (2005), DOI 10.1063/1.1854726"
            ),
            "local_pdf": (
                "research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf"),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 4,
            "print_page": 63,
            "figure": "Figure 6",
            "caption_reports_energy": False,
            "body_reports_figure6_energy": False,
            "render_dpi": 600,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "figure_crop_bounds_px": list(CROP_BOUNDS),
        },
        "independent_condition_audit": {
            "source": (
                "J. Arts et al., Applied Physics Reviews 8, 041305 (2021), "
                "DOI 10.1063/5.0058904, Figure 24(d)"
            ),
            "finding": (
                "The review reproduces the Karahashi angular panel and "
                "explicitly labels its kinetic energy unknown."
            ),
            "reported_energy_eV": None,
        },
        "experiment_scope": {
            "target": "SiO2",
            "incident_ions": ["CF+", "CF2+", "CF3+", "Ar+"],
            "neutral_radical_flux": "none",
            "beam_selection": "mass selected",
            "incidence_angles_deg": {
                "CFx+": [0, 30, 45, 60, 75],
                "Ar+": [0, 30, 50, 60, 75],
            },
            "quantity": "etched SiO2 units per incident ion",
            "ion_energy_eV": None,
            "ion_energy_status": "not_reported_by_source",
        },
        "pixel_calibration": {
            "crop_bounds_full_render_px": list(CROP_BOUNDS),
            "x_at_0_deg": X_AT_0_DEG,
            "x_at_90_deg": X_AT_90_DEG,
            "y_at_zero_yield": Y_AT_0_YIELD,
            "y_at_2_yield": Y_AT_2_YIELD,
            "transform": {
                "angle_deg": (
                    "90 * (x_px - 348.5) / (1422.0 - 348.5)"),
                "yield": (
                    "2 * (993.5 - y_px) / (993.5 - 206.0)"),
            },
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; exact PIL crop; dark-pixel axis "
                "localization; marker-center transcription; Hough-circle and "
                "connected-contour checks; full-resolution visual overlay"
            ),
            "yield_digitization_uncertainty": 0.015,
            "uncertainty_basis": (
                "5.5-pixel allowance including overlap at the crowded "
                "60-degree markers, rounded upward through the y transform"
            ),
            "nominal_angle_policy": (
                "use printed angular setpoints; retain marker x pixels as "
                "independent placement check"
            ),
        },
        "derived_checks": {
            "point_count": len(PIXEL_POINTS),
            "maximum_marker_x_angle_offset_deg": round(max(x_offsets), 3),
            "yield_ratio_60_to_0": {
                species: round(_ratio(species, 60), 4)
                for species in ("CF+", "CF2+", "CF3+", "Ar+")
            },
            "normal_incidence_cross_figure_candidate_rmse": candidate_rmse,
            "best_cross_figure_candidate_energy_eV": 1000,
            "best_candidate_status": (
                "strong_inference_only_not_source_reported"),
        },
        "claim_boundary": {
            "valid": [
                (
                    "qualitative species ordering of angular anisotropy at "
                    "the source's unreported energy condition"),
                (
                    "digitized 60-to-0 ratios with the missing-energy warning "
                    "carried alongside every row"),
                (
                    "cross-figure evidence that 1000 eV is the most likely "
                    "condition")
            ],
            "not_valid": [
                "an energy-resolved production angular law",
                "an assumption that the source reported 1000 eV",
                "extrapolation to Krueger's 471-4821 eV IEAD",
                "a fit target for absolute feature depth",
            ],
            "production_surface_model_use": False,
        },
        "output": {
            "path": (
                "data/experimental/karahashi_2007/"
                "figure6_angular_yields_unknown_energy.csv"),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    csv_sha = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(csv_sha), indent=2, ensure_ascii=False) + "\n"


def _assert_source() -> None:
    if _sha256(SOURCE_PDF) != SOURCE_PDF_SHA256:
        raise RuntimeError("Karahashi source PDF checksum does not match")


def verify_render(render_path: Path) -> Image.Image:
    if _sha256(render_path) != RENDER_SHA256:
        raise RuntimeError("600-dpi Figure 6 page render checksum does not match")
    image = Image.open(render_path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected render size {image.size}; expected {RENDER_SIZE}")
    crop = image.crop(CROP_BOUNDS)
    gray = np.asarray(crop.convert("L"))
    checks = {
        "left_axis": float(np.mean(
            gray[84:997, 347:351] < 96)),
        "right_axis": float(np.mean(
            gray[84:997, 1420:1424] < 96)),
        "top_axis": float(np.mean(
            gray[84:89, 348:1423] < 96)),
        "bottom_axis": float(np.mean(
            gray[991:997, 348:1423] < 96)),
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


def _draw_cross(
        draw: ImageDraw.ImageDraw, x: float, y: float, color: str) -> None:
    radius = 11
    draw.line((x - radius, y, x + radius, y), fill=color, width=3)
    draw.line((x, y - radius, x, y + radius), fill=color, width=3)
    draw.ellipse(
        (x - 17, y - 17, x + 17, y + 17), outline=color, width=3)


def draw_overlay(crop: Image.Image, output: Path) -> None:
    overlay = crop.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {
        "CF+": "#377eb8",
        "CF2+": "#4daf4a",
        "CF3+": "#e41a1c",
        "Ar+": "#984ea3",
    }
    for point in PIXEL_POINTS:
        _draw_cross(
            draw, point.x_px, point.y_px, colors[point.species])
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main():
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
