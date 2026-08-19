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
    }


def manifest(payloads: dict[Path, str]) -> dict:
    return {
        "manifest_id": "JANISSEN-2016-TIO2-SUPPLEMENTARY-TABLES-R1",
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
                "equipment-specific DC self-bias witnesses",
            ],
            "not_valid": [
                "transfer of Fluor Z401S voltage to Oxford NPG80",
                "transfer from single-crystal rutile to ALD TiO2",
                "a CHF3/SF6/O2 surface coefficient",
                "an uncertainty-free absolute-depth prediction",
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
