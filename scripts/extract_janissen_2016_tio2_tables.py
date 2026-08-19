#!/usr/bin/env python3
"""Replay the visually audited Janissen/Ha TiO2 supplementary tables."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "experimental" / "janissen_2016_tio2"
TABLE_S31_CSV = OUTPUT_DIR / "table_s3_1_chf3_rie_optimization.csv"
TABLE_S32_CSV = OUTPUT_DIR / "table_s3_2_tio2_nanofabrication.csv"
TABLE_S33_CSV = OUTPUT_DIR / "table_s3_3_feature_dimensions.csv"
FIGURE_S33_CSV = OUTPUT_DIR / "figure_s3_3_oxygen_profile_response.csv"
MANIFEST_PATH = OUTPUT_DIR / "extraction_manifest.json"

SOURCE_PDF_SHA256 = (
    "d4e8afedd1349a91b14ee10f589e9cae01aab48eb31d2be44433cd26e6fce912"
)
RENDER_SIZE = (2055, 2884)
RENDER_DPI = 300
PAGE_54_RENDER_SHA256 = (
    "83e144798df436a8bf6b411b5ce599f9044c4079afd659dba384ea66fa4ed4a2"
)
PAGE_55_RENDER_SHA256 = (
    "d79e4a59189b58befa426d13e67b8e8a354c5178b010439ede32f670c0e03e4b"
)
PAGE_56_RENDER_SHA256 = (
    "787c95b14b45b7c544487504eda0f80adc8ec99659ee39363d034f58385498f5"
)
TABLE_BOXES = {
    "pdf_page_54_table_s3_1": (280, 400, 2015, 1020),
    "pdf_page_55_table_s3_2": (280, 325, 2015, 2360),
    "pdf_page_56_table_s3_3": (280, 280, 2015, 890),
}

# Independent 400-dpi replay of thesis print page 48 (PDF page 59).  The TU
# Delft portal regenerated container metadata on the 2026-08-19 download, so
# both that exact PDF and the rendered page are pinned here rather than
# replacing the earlier table-source checksum.
FIGURE_S33_SOURCE_PDF_SHA256 = (
    "097ca9738beea5ae04438a50f8bcb3b24012d6be398036df99e17bbfe17b3ed1"
)
FIGURE_S33_RENDER_SHA256 = (
    "0e50197123bec731e9a4e4d7c82bedfc2320cc280f4339024ded5ae4f4aae0ee"
)

S31_FIELDS = (
    "sample", "batch", "mask_height_nm", "mask_diameter_nm", "CHF3_sccm",
    "O2_sccm", "Ar_sccm", "rf_power_W", "pressure_ubar",
    "tio2_etch_rate_nm_min", "cr_etch_rate_nm_min",
    "tio2_to_cr_selectivity",
)
S31_ROWS = (
    ("Ra1", "Ra", 100, 535, 50, 5, 30, 200, 50, 68, 3.7, 18),
    ("Ra2", "Ra", 100, 535, 50, 5, 30, 165, 50, 58, 3.2, 18),
    ("Ra3", "Ra", 100, 535, 50, 5, 30, 100, 50, 30, 1.7, 18),
    ("Ra4", "Ra", 100, 535, 50, 5, 30, 200, 10, 30, 8.3, 4),
    ("Rb1", "Rb", 100, 345, 50, 5, 30, 200, 50, 77, 4.1, 19),
    ("Rb2", "Rb", 100, 345, 50, 5, 0, 200, 50, 73, 3.5, 21),
)

S32_FIELDS = (
    "figure", "system_id", "system_model", "mask_height_nm",
    "mask_diameter_nm", "CHF3_sccm", "O2_sccm", "SF6_sccm", "CH4_sccm",
    "Ar_sccm", "He_sccm", "icp_power_W", "rf_power_W", "pressure_ubar",
    "dc_bias_V_signed", "etch_time_min", "sample_holder_temperature_C",
    "chamber_temperature_C", "optimized",
)
S32_ROWS = (
    ("3.2a", "F1", "Fluor Z401S", 45, 175, 50, 0.5, "", "", "", "",
     "", 200, 50, -950, 11, "", "", True),
    ("3.2b", "F2", "Fluor Z401S", 100, 255, 50, 8, "", "", "", "",
     "", 200, 50, -1100, 15, "", "", True),
    ("3.2c", "F1", "Fluor Z401S", 100, 535, 50, 5, "", "", 30, "",
     "", 165, 50, -855, 15, "", "", False),
    ("3.2d", "F3", "AMS100 I-speeder", 100, 535, "", "", 15, 30, 50,
     "", 2500, 300, 30, -50, 10, 0, 200, False),
    ("3.2e", "F4", "Plasmalab System 100", 130, 185, "", "", 50, "",
     "", 100, 800, 200, 50, -475, 12, "", 25, False),
    ("3.2f", "F4", "Plasmalab System 100", 120, 205, "", "", 20, "",
     "", 100, 300, 300, 100, -835, 4, "", 25, True),
    ("3.3", "F2", "Fluor Z401S", 60, 220, 50, 4, "", "", "", "",
     "", 200, 50, -1100, 15, "", "", True),
    ("S3.4", "F2", "Fluor Z401S", 30, 190, 50, 4, "", "", "", "",
     "", 200, 50, -1100, 8, "", "", False),
)

S33_FIELDS = (
    "figure", "etch_time_min", "average_top_diameter_nm",
    "average_bottom_diameter_nm", "average_height_nm", "average_top_roundness",
    "average_bottom_roundness", "average_volume_um3", "height_local_rsd_percent",
    "height_global_rsd_percent",
)
S33_ROWS = (
    ("3.3", 15, 151, 215, 652, 0.97, 0.98, 0.017, 0.6, 1.4),
    ("S3.4", 8, 149, 187, 273, 0.98, 0.99, 0.006, 1.3, 3.1),
)

FIGURE_S33_FIELDS = (
    "figure_panel", "O2_sccm", "etch_time_min", "profile_class",
    "top_diameter_nm", "bottom_diameter_nm", "waist_diameter_nm",
    "upper_height_nm", "lower_height_nm", "total_height_nm",
    "reported_sidewall_angle_deg", "upper_sidewall_angle_deg",
    "lower_sidewall_angle_deg",
    "tio2_rate_from_rounded_height_nm_min", "cr_rate_digitized_nm_min",
    "cr_rate_digitization_uncertainty_nm_min",
    "selectivity_digitized", "selectivity_digitization_uncertainty",
    "selectivity_from_height_and_cr_rate",
)


def _figure_s33_row(
        panel, oxygen, profile, top, bottom, waist, upper_height,
        lower_height, total_height, sidewall_angle, upper_angle, lower_angle,
        cr_rate, plotted_selectivity):
    tio2_rate = float(total_height) / 15.0
    return (
        panel, oxygen, 15.0, profile, top, bottom, waist, upper_height,
        lower_height, total_height, sidewall_angle, upper_angle, lower_angle,
        tio2_rate, cr_rate, 0.1, plotted_selectivity, 0.5,
        tio2_rate / cr_rate,
    )


FIGURE_S33_ROWS = (
    _figure_s33_row(
        "S3.3a", 0.0, "positive_sidewall", 285, 395, "", "", "",
        500, 84, "", "", 2.7, 12.0),
    _figure_s33_row(
        "S3.3b", 0.5, "vertical_sidewall", 275, 275, "", "", "",
        545, 90, "", "", 2.7, 13.5),
    _figure_s33_row(
        "S3.3c", 1.0, "negative_sidewall", 305, 260, "", "", "",
        520, -88, "", "", 2.5, 14.0),
    _figure_s33_row(
        "S3.3d", 5.0, "symmetric_hourglass", 260, 260, 175, 555, 555,
        1110, "", -86, 86, 3.4, 22.0),
    _figure_s33_row(
        "S3.3e", 10.0, "asymmetric_hourglass", 245, 370, 135, 385,
        1085, 1470, "", -82, 84, 3.5, 28.0),
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _csv_text(fields: tuple[str, ...], rows: tuple[tuple, ...]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(fields)
    writer.writerows(rows)
    return stream.getvalue()


def outputs() -> dict[Path, str]:
    return {
        TABLE_S31_CSV: _csv_text(S31_FIELDS, S31_ROWS),
        TABLE_S32_CSV: _csv_text(S32_FIELDS, S32_ROWS),
        TABLE_S33_CSV: _csv_text(S33_FIELDS, S33_ROWS),
        FIGURE_S33_CSV: _csv_text(FIGURE_S33_FIELDS, FIGURE_S33_ROWS),
    }


def manifest(payloads: dict[Path, str]) -> dict:
    return {
        "manifest_id": "HA-2016-TIO2-SUPPLEMENTARY-DATA-R2",
        "source": {
            "citation": (
                "S. Ha et al., Nanoscale 8, 10739-10748 (2016); "
                "supplement reproduced in S. Ha, TU Delft thesis (2018)"
            ),
            "doi": "10.1039/C6NR00898D",
            "thesis_url": (
                "https://pure.tudelft.nl/ws/portalfiles/portal/47058891/"
                "PhD_Thesis_S.Ha.pdf"
            ),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "source_pixel_redistribution": "not committed",
        },
        "visual_audit": {
            "render_dpi": RENDER_DPI,
            "render_size_px": list(RENDER_SIZE),
            "pdf_page_54_sha256": PAGE_54_RENDER_SHA256,
            "pdf_page_55_sha256": PAGE_55_RENDER_SHA256,
            "pdf_page_56_sha256": PAGE_56_RENDER_SHA256,
            "table_boxes_full_page_px": {
                key: list(value) for key, value in TABLE_BOXES.items()
            },
            "status": "passed_original_resolution",
            "method": (
                "300-dpi Poppler renders; original-resolution visual reading; "
                "PIL crop-box overlay; values transcribed from printed tables"
            ),
            "figure_s3_3": {
                "source_pdf_sha256": FIGURE_S33_SOURCE_PDF_SHA256,
                "source_pdf_page": 59,
                "source_print_page": 48,
                "render_dpi": 400,
                "render_size_px": [2739, 3845],
                "render_sha256": FIGURE_S33_RENDER_SHA256,
                "caption_values": (
                    "verbatim dimensions and profile classes; TiO2 rate "
                    "calculated as rounded caption height divided by the "
                    "printed 15 minute etch time"
                ),
                "plot_values": (
                    "Cr rate and selectivity read to the plotted precision; "
                    "uncertainty bounds cover marker thickness and axis reading"
                ),
                "axis_calibration_pixels_in_2300x700_plot_crop": {
                    "TiO2_rate": {
                        "x_ticks_0_5_10": [216, 409, 602],
                        "y_ticks_0_25_50_75_100": [382, 303, 225, 146, 67],
                    },
                    "Cr_rate": {
                        "x_ticks_0_5_10": [874, 1066, 1260],
                        "y_ticks_0_1_2_3_4_5": [382, 323, 264, 206, 147, 89],
                    },
                    "selectivity": {
                        "x_ticks_0_5_10": [1531, 1725, 1918],
                        "y_ticks_0_10_20_30": [382, 285, 187, 90],
                    },
                },
                "status": "passed_original_resolution_and_ratio_crosscheck",
            },
        },
        "outputs": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(text.encode("utf-8")).hexdigest(),
            }
            for path, text in payloads.items()
        ],
        "claim_boundary": {
            "valid": [
                "printed process settings and rates in Tables S3.1-S3.3",
                "measured feature dimensions and printed RSD values",
                "Figure S3.3 caption dimensions and profile classes",
                "Figure S3.3 Cr-rate and selectivity values within digitization uncertainty",
                "equipment-specific DC self-bias witnesses",
            ],
            "not_valid": [
                "transfer of Fluor Z401S voltage to Oxford NPG80",
                "transfer from single-crystal rutile to ALD TiO2",
                "a CHF3/SF6/O2 surface coefficient",
                "an uncertainty-free absolute-depth prediction",
                "transfer of the single-crystal oxygen sweep to ALD TiO2 or an SF6-containing feed",
            ],
        },
    }


def manifest_text(payloads: dict[Path, str]) -> str:
    return json.dumps(manifest(payloads), indent=2) + "\n"


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Janissen/Ha source PDF checksum changed")


def verify_render(path: Path, expected_sha: str) -> Image.Image:
    if _sha256(path) != expected_sha:
        raise RuntimeError("Janissen/Ha page render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Janissen/Ha render size: {image.size}")
    return image


def draw_overlay(image: Image.Image, box: tuple[int, int, int, int], path: Path):
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    draw.rectangle(box, outline="#e41a1c", width=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(path)


def verify_committed_files() -> None:
    payloads = outputs()
    for path, expected in payloads.items():
        if path.read_text(encoding="utf-8") != expected:
            raise RuntimeError(f"committed table is stale: {path}")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payloads):
        raise RuntimeError("committed Janissen/Ha manifest is stale")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--page-54-render", type=Path)
    parser.add_argument("--page-55-render", type=Path)
    parser.add_argument("--page-56-render", type=Path)
    parser.add_argument("--overlay-dir", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.source_pdf is not None:
        verify_source_pdf(args.source_pdf)
    images = {}
    if args.page_54_render is not None:
        images[54] = verify_render(args.page_54_render, PAGE_54_RENDER_SHA256)
    if args.page_55_render is not None:
        images[55] = verify_render(args.page_55_render, PAGE_55_RENDER_SHA256)
    if args.page_56_render is not None:
        images[56] = verify_render(args.page_56_render, PAGE_56_RENDER_SHA256)
    if args.overlay_dir is not None:
        if set(images) != {54, 55, 56}:
            raise SystemExit("--overlay-dir requires all three page renders")
        draw_overlay(images[54], TABLE_BOXES["pdf_page_54_table_s3_1"],
                     args.overlay_dir / "page54-table-s3-1-overlay.png")
        draw_overlay(images[55], TABLE_BOXES["pdf_page_55_table_s3_2"],
                     args.overlay_dir / "page55-table-s3-2-overlay.png")
        draw_overlay(images[56], TABLE_BOXES["pdf_page_56_table_s3_3"],
                     args.overlay_dir / "page56-table-s3-3-overlay.png")
    if args.check:
        verify_committed_files()
        return
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payloads = outputs()
    for path, content in payloads.items():
        path.write_text(content, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payloads), encoding="utf-8")


if __name__ == "__main__":
    main()
