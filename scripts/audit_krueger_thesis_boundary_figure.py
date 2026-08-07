#!/usr/bin/env python3
"""Audit the boundary information visible in Krueger thesis Figure 6.17.

This is deliberately a label/completeness audit, not a curve digitization.
The scientific result is the absence of two required boundary variables, so
inventing values from nearby curves would defeat the purpose of the audit.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "results" / "curated" / "krueger_thesis_boundary_visual_audit"
    / "audit.json")

PDF_SHA256 = (
    "4929018bfeb1cfdbbbc5ba5aecfd6a6010a38ee883dc0f1a574c090bd132c89a")
RENDER_SHA256 = (
    "c63e53b63481ff2a23482d5c74a3d5be395916100d9912a12601e1c88d4dc919")
RENDER_SIZE = (3400, 4400)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def verify_source_files(source_pdf: Path, rendered_page: Path) -> None:
    """Fail closed unless the exact visually inspected source is supplied."""
    if _sha256(source_pdf) != PDF_SHA256:
        raise ValueError("Krueger thesis PDF checksum mismatch")
    if _sha256(rendered_page) != RENDER_SHA256:
        raise ValueError("Krueger Figure 6.17 render checksum mismatch")
    with Image.open(rendered_page) as image:
        if image.size != RENDER_SIZE:
            raise ValueError(
                "Krueger Figure 6.17 render dimensions mismatch: "
                f"{image.size} != {RENDER_SIZE}")
        image.verify()


def build_audit() -> dict:
    """Return the frozen human-vision/PIL completeness audit."""
    return {
        "audit_id": "KRUEGER-2024-THESIS-FIG6.17-BOUNDARY-VISUAL-R1",
        "source": {
            "title": (
                "Modeling and Optimization of High Aspect Ratio Plasma "
                "Etching"),
            "author": "Florian Krueger",
            "year": 2024,
            "doi": "10.7302/23106",
            "figure": "6.17",
            "printed_page": 204,
            "pdf_page_1_based": 228,
            "pdf_sha256": PDF_SHA256,
            "render": {
                "dpi": 400,
                "pixel_width": RENDER_SIZE[0],
                "pixel_height": RENDER_SIZE[1],
                "sha256": RENDER_SHA256,
            },
        },
        "method": {
            "type": "full_resolution_human_vision_plus_PIL_integrity",
            "vision_task": (
                "read every curve label and inspect the figure for a stable "
                "C4F6 parent curve or species-resolved positive-ion curves"),
            "pil_task": (
                "open the exact 400 dpi page render, verify image integrity "
                "and dimensions, and bind it to the source by SHA-256"),
            "numerical_curve_digitization_performed": False,
            "reason_no_curve_digitization": (
                "the audit asks whether required boundary variables are "
                "published, not for values of the curves that are present"),
        },
        "figure_condition": {
            "low_frequency_power_kW": 6.0,
            "high_frequency_power_kW": 2.5,
            "o2_to_c4f6_feed_ratios": [0.5, 1.0, 1.5, 2.5],
            "scope": (
                "oxygen-flow transfer study; this is not the 8 kW base-case "
                "boundary"),
        },
        "visual_transcription": {
            "curve_labels": [
                "CF2", "C3F4", "O", "C2F3", "CF", "Ions", "CF3", "CO", "C",
            ],
            "stable_c4f6_parent_curve_present": False,
            "positive_ion_species_resolved": False,
            "positive_ion_label": "Ions",
        },
        "claim_boundary": {
            "proves": [
                (
                    "Figure 6.17 does not publish stable C4F6 wafer flux for "
                    "the plotted transfer conditions"),
                (
                    "Figure 6.17 reports positive ions only as one aggregate "
                    "curve"),
            ],
            "does_not_prove": [
                "stable C4F6 wafer flux is zero",
                "the aggregate ion flux is species independent",
                "the plotted transfer-condition fluxes equal the base case",
                "any numerical C4F6 or species-resolved ion boundary",
            ],
            "depth_consequence": (
                "the thesis figure does not supply the stable-parent flux or "
                "positive-ion composition needed to close an independently "
                "predictive Krueger absolute-depth calculation"),
        },
    }


def audit_text() -> str:
    return json.dumps(build_audit(), indent=2) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--rendered-page", type=Path)
    args = parser.parse_args()
    if (args.source_pdf is None) != (args.rendered_page is None):
        parser.error("--source-pdf and --rendered-page must be supplied together")
    if args.source_pdf is not None:
        verify_source_files(args.source_pdf, args.rendered_page)
    payload = audit_text()
    if OUTPUT.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Figure 6.17 boundary audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
