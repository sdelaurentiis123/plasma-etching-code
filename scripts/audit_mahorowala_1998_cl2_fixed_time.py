#!/usr/bin/env python3
"""Reproduce Mahorowala's fixed-time Cl2/oxide-mask rate board.

The 1998 MIT thesis is viewable but may not be redistributed.  This script
therefore checksum-binds an authorized local copy and its Poppler renders,
reproduces the committed Table-2.2 transcription, and draws source-audit
overlays without committing source pixels.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = (
    ROOT / "data" / "experimental" / "mahorowala_1998_cl2"
)
CSV_PATH = OUTPUT_DIRECTORY / "table2_2_oxide_mask_fixed_time.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "audit_manifest.json"

SOURCE_PDF_SHA256 = (
    "e561dfe9780c0e27439b2b7288788c18333b8228152b6cc40baf72ba2edf4b6a"
)
TABLE_RENDER_SHA256 = (
    "af753b4aabb5fbe8de7f82c0feedac1273290b67635154e007105aa4c8b20367"
)
PROFILE_RENDER_SHA256 = (
    "3b00616b38ffe691a8c8c91d3bf6139bebebad9e923956bb7c88141dde71b156"
)
RENDER_SIZE = (5100, 6600)

FIELDNAMES = (
    "run",
    "inductive_power_W",
    "rf_bias_power_W",
    "cl2_flow_sccm",
    "pressure_mTorr",
    "etch_time_s",
    "poly_si_etch_rate_A_min",
    "oxide_etch_rate_A_min",
    "selectivity",
    "derived_poly_si_removed_nm",
    "derived_oxide_removed_nm",
    "source_profile_panel",
    "quantitative_status",
)

SOURCE_ROWS = (
    (1, 400, 80, 100, 3150, 125, "25.20"),
    (2, 550, 80, 175, 3675, 250, "14.70"),
    (3, 250, 80, 25, 2000, 25, "80.00"),
    (4, 400, 20, 175, 1400, 225, "6.22"),
    (5, 550, 20, 100, 1300, 100, "13.00"),
    (6, 550, 80, 25, 2375, 200, "11.88"),
    (7, 250, 80, 175, 2125, 200, "10.63"),
    (8, 550, 140, 100, None, None, ""),
    (9, 400, 140, 25, 2900, 475, "6.11"),
    (10, 400, 20, 25, 900, 225, "4.00"),
    (11, 250, 20, 100, 1075, 25, "43.00"),
    (12, 400, 140, 175, None, None, ""),
    (13, 250, 140, 100, 2900, 300, "9.67"),
)

# Original-pixel Table-2.2 grid support on the 600-dpi full-page render.
TABLE_X_GRID_PX = (563, 873, 1362, 1773, 2220, 2630, 3192, 3643, 4246)
TABLE_TOP_PX = 4445
TABLE_BOTTOM_PX = 6003


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows() -> list[dict[str, str]]:
    output = []
    for run, top, bias, flow, poly_rate, oxide_rate, selectivity in SOURCE_ROWS:
        usable = poly_rate is not None
        output.append(
            {
                "run": str(run),
                "inductive_power_W": str(top),
                "rf_bias_power_W": str(bias),
                "cl2_flow_sccm": str(flow),
                "pressure_mTorr": "10",
                "etch_time_s": "75",
                "poly_si_etch_rate_A_min": "" if not usable else str(poly_rate),
                "oxide_etch_rate_A_min": "" if not usable else str(oxide_rate),
                "selectivity": selectivity,
                "derived_poly_si_removed_nm": (
                    "" if not usable else f"{0.125 * poly_rate:.3f}"
                ),
                "derived_oxide_removed_nm": (
                    "" if not usable else f"{0.125 * oxide_rate:.3f}"
                ),
                "source_profile_panel": f"figure2.4_panel{run}",
                "quantitative_status": (
                    "usable" if usable else "overetched_no_rate"
                ),
            }
        )
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_digest: str) -> dict[str, object]:
    usable = [row for row in rows() if row["quantitative_status"] == "usable"]
    depths = [float(row["derived_poly_si_removed_nm"]) for row in usable]
    return {
        "schema": "petch.mahorowala-1998-cl2-fixed-time.v1",
        "audit_id": "MAHOROWALA-1998-CL2-TABLE2.2-R1",
        "source": {
            "citation": (
                "A. P. Mahorowala, Feature profile evolution during the "
                "high density plasma etching of polysilicon, MIT PhD "
                "thesis (1998)"
            ),
            "handle": "http://hdl.handle.net/1721.1/50514",
            "dspace_item_uuid": "9af24424-813c-401d-9ca1-33cb863d86c3",
            "dspace_bitstream_uuid": "c04d2da3-cbc2-4898-97ed-c83c58293e05",
            "bitstream_name": "42415621-MIT.pdf",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_pages": 200,
            "rights": (
                "MIT thesis pixels are viewed locally and not redistributed"
            ),
        },
        "source_locations": {
            "table": {
                "name": "Table 2.2",
                "pdf_page": 39,
                "print_page": 39,
                "render_dpi": 600,
                "render_sha256": TABLE_RENDER_SHA256,
                "render_size_px": list(RENDER_SIZE),
            },
            "profile_montage": {
                "name": "Figure 2.4",
                "pdf_page": 42,
                "print_page": 42,
                "render_dpi": 600,
                "render_sha256": PROFILE_RENDER_SHA256,
                "render_size_px": list(RENDER_SIZE),
                "panel_count": 13,
            },
        },
        "experiment": {
            "reactor": "Lam TCP 9400SE inductively coupled plasma",
            "chemistry": "pure Cl2",
            "wafer_temperature_C": 60.0,
            "pressure_mTorr": 10.0,
            "etch_time_s": 75.0,
            "initial_poly_si_thickness_nm": 500.0,
            "oxide_mask_thickness_nm": 200.0,
            "line_width_nm": 250.0,
            "figure2_4_spacing_nm": 310.0,
            "sample_classes": (
                "alternate lines/spaces; approximately 60% exposed polysilicon"
            ),
        },
        "transcription": {
            "method": (
                "source Table 2.2 transcription checked on a checksum-bound "
                "600-dpi Poppler render; Figure 2.4 inspected at original "
                "render resolution"
            ),
            "source_run_count": 13,
            "usable_quantitative_run_count": len(usable),
            "overetched_no_rate_runs": [8, 12],
            "reported_measurement_uncertainty": None,
            "derived_depth_formula": (
                "rate_A_per_min * 75 s / 60 s_per_min * 0.1 nm_per_A"
            ),
            "derived_poly_si_depth_range_nm": [min(depths), max(depths)],
            "visual_audit_status": "passed_original_resolution",
        },
        "identifiability": {
            "absolute_time_published": True,
            "absolute_rate_published": True,
            "species_resolved_wafer_flux_published": False,
            "measured_iead_published": False,
            "measured_iad_published": False,
            "center_condition_estimates_in_thesis": (
                "2 mA/cm2 ion current and about 100-120 eV; estimates, not "
                "per-run measurements"
            ),
            "valid_use": (
                "absolute fixed-time rate/depth and SEM-profile board for a "
                "reactor provider or facility-conditioned transfer"
            ),
            "invalid_use": (
                "first-principles knobs-to-depth validation before the "
                "species/IEAD boundary is independently identified"
            ),
            "feature_profile_formal_pass_granted": False,
        },
        "output": {
            "path": (
                "data/experimental/mahorowala_1998_cl2/"
                "table2_2_oxide_mask_fixed_time.csv"
            ),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    return json.dumps(
        manifest(sha256(csv_payload.encode("utf-8")).hexdigest()),
        indent=2,
    ) + "\n"


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Mahorowala Table-2.2 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Mahorowala manifest is stale")


def _verify_image(path: Path, expected_digest: str) -> Image.Image:
    if _sha256(path) != expected_digest:
        raise RuntimeError(f"render checksum changed: {path}")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"render size changed: {image.size}")
    return image


def verify_table_render(path: Path) -> Image.Image:
    image = _verify_image(path, TABLE_RENDER_SHA256)
    gray = np.asarray(image.convert("L"))
    for x_px in TABLE_X_GRID_PX:
        support = int(
            np.sum(
                gray[
                    TABLE_TOP_PX - 4 : TABLE_BOTTOM_PX + 5,
                    x_px - 4 : x_px + 5,
                ]
                < 112
            )
        )
        if support < 5000:
            raise RuntimeError(
                f"insufficient Table-2.2 grid support at x={x_px}: {support}"
            )
    return image


def draw_table_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    draw.rectangle(
        (
            TABLE_X_GRID_PX[0],
            TABLE_TOP_PX,
            TABLE_X_GRID_PX[-1],
            TABLE_BOTTOM_PX,
        ),
        outline="#e41a1c",
        width=8,
    )
    for x_px in TABLE_X_GRID_PX:
        draw.line(
            (x_px, TABLE_TOP_PX, x_px, TABLE_BOTTOM_PX),
            fill="#377eb8",
            width=3,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--table-render", type=Path)
    parser.add_argument("--profile-render", type=Path)
    parser.add_argument("--table-overlay", type=Path)
    args = parser.parse_args()

    verify_committed_files()
    if args.source_pdf is not None and _sha256(args.source_pdf) != SOURCE_PDF_SHA256:
        raise RuntimeError("Mahorowala source PDF checksum changed")
    table_image = (
        verify_table_render(args.table_render)
        if args.table_render is not None
        else None
    )
    if args.profile_render is not None:
        _verify_image(args.profile_render, PROFILE_RENDER_SHA256)
    if args.table_overlay is not None:
        if table_image is None:
            raise ValueError("--table-overlay requires --table-render")
        draw_table_overlay(table_image, args.table_overlay)
    print(
        json.dumps(
            {
                "status": "verified",
                "source_run_count": len(SOURCE_ROWS),
                "usable_quantitative_run_count": sum(
                    row["quantitative_status"] == "usable" for row in rows()
                ),
                "source_pdf_verified": args.source_pdf is not None,
                "table_render_verified": table_image is not None,
                "profile_render_verified": args.profile_render is not None,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
