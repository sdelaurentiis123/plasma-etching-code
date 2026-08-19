#!/usr/bin/env python3
"""Transcribe and verify Ji et al. Figure-4 spacing-response annotations."""
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
    ROOT / "data" / "experimental" / "ji_2024_tio2_hierarchical"
)
CSV_PATH = OUTPUT_DIRECTORY / "figure4_spacing_morphology_response.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure4_digitization_manifest.json"
SOURCE_PDF_SHA256 = (
    "92a57c600e93e4113e574ef286aae4a248d787c3445f8d21e2ff151c73647e81"
)
RENDER_SHA256 = (
    "9c6348025b0214c5d43f482e14d48f9b1fa445aec9a4bc44ce7517365582cbdc"
)
RENDER_SIZE = (4961, 7016)
RENDER_DPI = 600

# Values are printed in yellow on the five source SEM panels and independently
# represented by Figure 4(f)'s grouped bars. They are transcribed at the
# source's printed precision; no experimental error bars are shown.
PRINTED = (
    # panel, designed gap, h1, h2, L, theta
    ("4a", 750.0, 414.0, 255.0, 261.0, 31.6),
    ("4b", 530.0, 399.0, 258.0, 264.0, 31.0),
    ("4c", 350.0, 399.0, 261.0, 261.0, 32.0),
    ("4d", 100.0, 277.0, 411.0, 276.0, 38.0),
    ("4e", 70.0, 267.0, 386.0, 264.0, 40.0),
)
ANNOTATION_REGIONS = {
    "4a": (1750, 1330, 2260, 1870),
    "4b": (2740, 1330, 3240, 1870),
    "4c": (3720, 1330, 4230, 1870),
    "4d": (1750, 2000, 2260, 2540),
    "4e": (1750, 2680, 2260, 3240),
}
FIELDNAMES = (
    "panel",
    "designed_gap_nm",
    "upper_triangle_height_h1_nm",
    "lower_rectangle_height_h2_nm",
    "feature_width_L_nm",
    "upper_triangle_angle_theta_deg",
    "value_origin",
    "source_measurement_uncertainty_reported",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rows() -> list[dict[str, str]]:
    return [
        {
            "panel": panel,
            "designed_gap_nm": f"{gap:.1f}",
            "upper_triangle_height_h1_nm": f"{h1:.1f}",
            "lower_rectangle_height_h2_nm": f"{h2:.1f}",
            "feature_width_L_nm": f"{width:.1f}",
            "upper_triangle_angle_theta_deg": f"{angle:.1f}",
            "value_origin": "printed_SEM_annotation_cross_checked_to_bar_chart",
            "source_measurement_uncertainty_reported": "false",
        }
        for panel, gap, h1, h2, width, angle in PRINTED
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_digest: str) -> dict[str, object]:
    gap = np.asarray([row[1] for row in PRINTED], dtype=float)
    h1 = np.asarray([row[2] for row in PRINTED], dtype=float)
    h2 = np.asarray([row[3] for row in PRINTED], dtype=float)
    width = np.asarray([row[4] for row in PRINTED], dtype=float)
    angle = np.asarray([row[5] for row in PRINTED], dtype=float)
    # The prose draws the regime boundary at 100 nm, but the plotted 100-nm
    # datum has already shifted strongly away from the 350--750 nm cluster.
    # Preserve that source-level ambiguity instead of manufacturing a sharp
    # threshold from only five samples.
    high_gap = gap > 100.0
    loading_gap = gap <= 100.0
    return {
        "manifest_id": "JI-2024-TIO2-FIG4-SPACING-MORPHOLOGY-R1",
        "source": {
            "citation": "X. Ji et al., Micromachines 15, 1160 (2024)",
            "doi": "10.3390/mi15091160",
            "license": "CC BY 4.0",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 7,
            "figure": "Figure 4(a-f)",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "source_pixels_committed": False,
        },
        "experiment": {
            "material": "800 nm electron-beam-deposited TiO2",
            "mask": "60 nm Cr",
            "nominal_feature_cd_nm": 200.0,
            "designed_gap_nm": gap.tolist(),
            "source_power_W": 350.0,
            "rf_power_W": 120.0,
            "feed_sccm": {"SF6": 40.0, "CHF3": 10.0, "O2": 5.0},
            "pressure_mTorr": 10.0,
            "temperature_C": 40.0,
            "etch_time_s": None,
        },
        "transcription": {
            "method": (
                "yellow SEM annotations transcribed at printed precision and "
                "cross-checked against the grouped bar chart on the same page"
            ),
            "point_count": len(PRINTED) * 4,
            "annotation_regions_full_page_px": {
                key: list(value) for key, value in ANNOTATION_REGIONS.items()
            },
            "source_measurement_uncertainty_reported": False,
            "visual_audit_status": "passed_original_resolution",
        },
        "derived_checks": {
            "source_text_boundary_claim": (
                "gaps_greater_than_100nm_stable_and_gaps_smaller_than_100nm_change"
            ),
            "digitized_100nm_point_already_shifted": True,
            "strict_threshold_identified": False,
            "empirical_transition_bracket_nm": [100.0, 350.0],
            "high_gap_h1_range_nm": [
                float(np.min(h1[high_gap])), float(np.max(h1[high_gap]))
            ],
            "high_gap_h2_range_nm": [
                float(np.min(h2[high_gap])), float(np.max(h2[high_gap]))
            ],
            "all_widths_within_15nm": float(np.ptp(width)) <= 15.0,
            "loading_gap_h1_below_every_high_gap_h1": bool(
                np.max(h1[loading_gap]) < np.min(h1[high_gap])
            ),
            "loading_gap_angle_above_every_high_gap_angle": bool(
                np.min(angle[loading_gap]) > np.max(angle[high_gap])
            ),
        },
        "freddie_geometry_implication": {
            "conditional_board_pitch_nm": 400.0,
            "conditional_board_width_nm": [
                80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 320.0
            ],
            "corresponding_gap_nm": [
                320.0, 280.0, 240.0, 200.0, 160.0, 120.0, 80.0
            ],
            "member_at_or_below_changed_100nm_point_width_nm": [320.0],
            "member_in_unsampled_transition_interval_width_nm": [280.0],
            "threshold_transfer_allowed": False,
            "use": (
                "predeclare enhanced pattern-loading sensitivity for the widest/"
                "smallest-gap board member and uncertainty for the 120 nm-gap member; "
                "do not impose a sharp Ji threshold on Oxford"
            ),
        },
        "claim_boundary": {
            "valid": (
                "same-gas TiO2/Cr pattern-loading response validation and a "
                "topology requirement for passivation-aware feature transport"
            ),
            "not_valid": [
                "an Oxford pattern-loading threshold",
                "a Freddie target geometry measurement",
                "an absolute etch-depth or time calibration",
            ],
        },
        "output": {
            "path": str(CSV_PATH.relative_to(ROOT)),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Ji Figure-4 CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Ji Figure-4 manifest is stale")


def verify_source_pdf(path: Path) -> None:
    if _sha256(path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Ji source PDF checksum changed")


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Ji 600-dpi page-7 render checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected Ji render size: {image.size}")
    rgb = np.asarray(image)
    for panel, (x0, y0, x1, y1) in ANNOTATION_REGIONS.items():
        patch = rgb[y0:y1, x0:x1]
        yellow = (
            (patch[:, :, 0] > 150)
            & (patch[:, :, 1] > 130)
            & (patch[:, :, 2] < 100)
        )
        if int(yellow.sum()) < 1000:
            raise RuntimeError(
                f"insufficient printed-annotation support at {panel}: "
                f"{int(yellow.sum())}"
            )
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for panel, region in ANNOTATION_REGIONS.items():
        draw.rectangle(region, outline="#377eb8", width=8)
        draw.text((region[0] + 8, region[1] + 8), panel, fill="#377eb8")
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render", type=Path)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.source_pdf:
        verify_source_pdf(args.source_pdf)
    image = verify_render(args.render) if args.render else None
    if args.overlay:
        if image is None:
            parser.error("--overlay requires --render")
        draw_overlay(image, args.overlay)
    payload = csv_text()
    if args.write:
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(payload, encoding="utf-8")
        MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")
    if args.check:
        verify_committed_files()
    if not any((args.write, args.check, args.source_pdf, args.render)):
        parser.error("select --write, --check, --source-pdf, or --render")


if __name__ == "__main__":
    main()
