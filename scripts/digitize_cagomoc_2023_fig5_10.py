#!/usr/bin/env python3
"""Replay Cagomoc 2023 dissertation Figure 5.10 digitization.

The figure is a classical-MD mechanism board, not an experiment.  It reports
steady-state Si removal from flat SiO2 under alternating 2000 eV CF3+ and CF3
radical injections.  The plotted radical/ion ratios are source setpoints; x
pixels are retained only as a placement check.
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
OUTPUT_DIR = ROOT / "data" / "surface_interactions" / "cagomoc_2023"
CSV_PATH = OUTPUT_DIR / "figure5_10_cf3_radical_ion_yield.csv"
MANIFEST_PATH = OUTPUT_DIR / "source_manifest.json"
DEFAULT_RENDER = (
    ROOT / "tmp" / "pdfs" / "cagomoc_pages" / "pdf121_fig5_10.png"
)

SOURCE_URL = (
    "https://ir.library.osaka-u.ac.jp/repo/ouka/all/91922/"
    "33239_Dissertation.pdf"
)
SOURCE_PDF_SHA256 = (
    "f1eb74b1b42bc12cc89fe60426b327c873fd75f68c2012c2f04e3c54adcade9f"
)
EXTRACTED_TEXT_SHA256 = (
    "904fc9293090f6cf666e2d9ce5723a74473aaf6a7c730ff842738052327fcc01"
)
RENDER_SHA256 = (
    "e8888593f82d107a15b329e91acf20af1bb60e59ac3fb39aebc4d11f52adacf0"
)
RENDER_SIZE = (3311, 4682)

# Full-page 400-dpi coordinates.  The right-hand y ticks independently locate
# 0, 1, 2, 3 and 4 at 1422, 1209.5, 997, 784.5 and 572 pixels.  The x
# calibration uses the centers of the 0 and 300 setpoint markers; intermediate
# marker offsets are retained and bounded below.
X_AT_RATIO_0 = 927.5
X_AT_RATIO_300 = 2573.0
Y_AT_YIELD_0 = 1422.0
Y_AT_YIELD_4 = 572.0
PIXEL_CENTER_BOUND = 5.5
YIELD_DIGITIZATION_BOUND = (
    PIXEL_CENTER_BOUND * 4.0 / (Y_AT_YIELD_0 - Y_AT_YIELD_4)
)


@dataclass(frozen=True)
class PixelPoint:
    radical_to_ion_ratio: int
    x_px: float
    y_px: float
    center_policy: str


POINTS = (
    PixelPoint(0, 927.7, 1037.5, "blue_disk_distance_core"),
    PixelPoint(25, 1066.8, 583.7, "blue_disk_distance_core"),
    PixelPoint(50, 1203.1, 610.1, "blue_disk_distance_core"),
    PixelPoint(100, 1478.5, 940.5, "blue_disk_distance_core"),
    # The dashed y=0 guide occludes both marker centers.  The dissertation text
    # independently identifies etch stop, so the axis center is the honest
    # value and the ambiguity remains in the pixel uncertainty.
    PixelPoint(200, 2027.0, 1422.0, "zero_axis_center_text_confirmed"),
    PixelPoint(300, 2573.0, 1422.0, "zero_axis_center_text_confirmed"),
)

FIELDNAMES = (
    "radical_to_ion_flux_ratio",
    "si_removal_yield_per_cf3_ion",
    "marker_center_x_px",
    "marker_center_y_px",
    "center_policy",
    "digitization_yield_bound",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _yield_at_pixel(y_px: float) -> float:
    return (
        4.0
        * (Y_AT_YIELD_0 - y_px)
        / (Y_AT_YIELD_0 - Y_AT_YIELD_4)
    )


def _ratio_at_pixel(x_px: float) -> float:
    return (
        300.0
        * (x_px - X_AT_RATIO_0)
        / (X_AT_RATIO_300 - X_AT_RATIO_0)
    )


def rows() -> list[dict[str, str]]:
    result = []
    for point in POINTS:
        result.append(
            {
                "radical_to_ion_flux_ratio": str(
                    point.radical_to_ion_ratio
                ),
                "si_removal_yield_per_cf3_ion": (
                    f"{max(0.0, _yield_at_pixel(point.y_px)):.4f}"
                ),
                "marker_center_x_px": f"{point.x_px:.1f}",
                "marker_center_y_px": f"{point.y_px:.1f}",
                "center_policy": point.center_policy,
                "digitization_yield_bound": (
                    f"{YIELD_DIGITIZATION_BOUND:.3f}"
                ),
            }
        )
    return result


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    x_offsets = [
        abs(
            _ratio_at_pixel(point.x_px)
            - point.radical_to_ion_ratio
        )
        for point in POINTS
    ]
    return {
        "manifest_id": "CAGOMOC-2023-FIG5.10-CF3-RADICAL-ION-R1",
        "source": {
            "citation": (
                "C. M. D. Cagomoc, Molecular Dynamics Simulation of SiO2 "
                "and SiN Etching for 3D NAND Memory Device Applications, "
                "Osaka University doctoral dissertation (2023)"
            ),
            "doi": "10.18910/91922",
            "handle": "https://hdl.handle.net/11094/91922",
            "official_pdf": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "local_extracted_text": (
                "research_sources/thesis_extracts/"
                "cagomoc_2023_dissertation.txt"
            ),
            "extracted_text_sha256": EXTRACTED_TEXT_SHA256,
            "pdf_page": 121,
            "printed_page": 106,
            "figure": "Figure 5.10",
            "render_dpi": 400,
            "render_sha256": RENDER_SHA256,
            "render_size_px": list(RENDER_SIZE),
            "pdf_redistributed": False,
        },
        "method_scope": {
            "evidence_class": (
                "classical molecular dynamics mechanism constraint; "
                "not an experimental validation target"
            ),
            "target": "flat SiO2",
            "ion": "CF3+ represented as a fast neutral at surface impact",
            "ion_energy_eV": 2000,
            "ion_incidence": "normal",
            "radical": "CF3",
            "radical_energy_eV": 0.5,
            "radical_incidence": (
                "normal for computational convenience, not an isotropic "
                "300 K gas distribution"
            ),
            "radical_injection_batch": "25 or 50 at a time",
            "reported_ratios": [0, 25, 50, 100, 200, 300],
            "reported_quantity": (
                "steady-state removed Si atoms per CF3+ ion injection"
            ),
        },
        "pixel_calibration": {
            "x_at_ratio_0": X_AT_RATIO_0,
            "x_at_ratio_300": X_AT_RATIO_300,
            "y_at_zero_yield": Y_AT_YIELD_0,
            "y_at_four_yield": Y_AT_YIELD_4,
            "right_axis_tick_centers_y_px": [
                1422.0,
                1209.5,
                997.0,
                784.5,
                572.0,
            ],
            "transform": {
                "ratio": (
                    "300 * (x_px - 927.5) / (2573.0 - 927.5)"
                ),
                "yield": (
                    "4 * (1422.0 - y_px) / (1422.0 - 572.0)"
                ),
            },
            "maximum_setpoint_placement_offset_ratio_units": max(
                x_offsets
            ),
        },
        "digitization": {
            "method": (
                "400-dpi Poppler render; Pillow RGB/dimension check; "
                "dark-pixel axis/tick localization; saturated-blue "
                "distance-core marker localization; original-resolution "
                "visual reconciliation"
            ),
            "marker_center_bound_px": PIXEL_CENTER_BOUND,
            "yield_bound": YIELD_DIGITIZATION_BOUND,
            "zero_point_policy": (
                "200:1 and 300:1 centers are occluded by the dashed zero "
                "guide; use its calibrated center because the caption and "
                "body independently state etch stop"
            ),
            "nominal_x_policy": (
                "source ratio setpoints are authoritative; x pixels only "
                "audit placement"
            ),
            "csv": CSV_PATH.name,
            "csv_sha256": csv_sha256,
        },
        "mechanism_constraints": [
            (
                "The response is non-monotone: an intermediate CF3 supply "
                "nearly doubles Si removal, while excess CF3 causes FC "
                "accumulation and etch stop."
            ),
            (
                "The same source reports few-nanometer C/F/O/Si mixing "
                "layers whose thickness varies with energy and material."
            ),
            (
                "Nanohole yields can be lower than flat yields because "
                "products redeposit on sidewalls and ions intercept tapered "
                "surfaces."
            ),
        ],
        "claim_boundaries": [
            (
                "Do not use this curve to identify Krueger's C4F6 reactor "
                "flux, ion mixture, IEAD, or 60 s depth."
            ),
            (
                "The source is classical MD, not a direct beam measurement "
                "or DFT-trained potential."
            ),
            (
                "The simulation removes non-covalently bonded products at "
                "finite cycle boundaries and acknowledges that slow "
                "in-hole redeposition may be underestimated."
            ),
            (
                "The 0.5 eV normal radical injection is a computational "
                "surrogate; Appendix D compares deposited layers with "
                "0.026 eV radicals but does not validate angular transfer "
                "into a feature."
            ),
        ],
        "other_visually_audited_pages": [
            {
                "pdf_page": 117,
                "figure": "5.7",
                "render_dpi": 300,
                "render_sha256": (
                    "7a9c1f39d237bfbaf458f27f7205499af3541c93ac1c0dcf"
                    "0319148483e92e09"
                ),
            },
            {
                "pdf_page": 118,
                "figure": "5.8",
                "render_dpi": 300,
                "render_sha256": (
                    "d1c6b57e3724a0267c53e5ea2e4b7a700c2e004e3d5f125"
                    "4008ffa11fa4a198a"
                ),
            },
            {
                "pdf_page": 120,
                "figure": "5.9",
                "render_dpi": 300,
                "render_sha256": (
                    "979c03673229b393d234d71e9d26fd08e3c5b70e2a6d986"
                    "c3bc45ced0e6e4fa2"
                ),
            },
        ],
    }


def _verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError(f"render checksum mismatch: {path}")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"render dimensions {image.size} != {RENDER_SIZE}"
        )
    array = np.asarray(image)
    dark = array.mean(axis=2) < 80
    for y_px in (572, 784, 997, 1209, 1422):
        if int(dark[y_px - 2 : y_px + 3, 2620:2640].sum()) < 40:
            raise RuntimeError(f"right-axis tick absent near y={y_px}")
    blue = (
        (array[:, :, 2] > 150)
        & (array[:, :, 0] < 100)
        & (array[:, :, 1] < 180)
    )
    for point in POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        if int(blue[y - 18 : y + 19, x - 18 : x + 19].sum()) < 150:
            raise RuntimeError(
                f"blue marker support absent near ratio "
                f"{point.radical_to_ion_ratio}"
            )
    return image


def _overlay(image: Image.Image, path: Path) -> None:
    canvas = image.copy()
    draw = ImageDraw.Draw(canvas)
    draw.line(
        [(X_AT_RATIO_0, Y_AT_YIELD_0), (X_AT_RATIO_300, Y_AT_YIELD_0)],
        fill=(255, 0, 255),
        width=3,
    )
    draw.line(
        [(X_AT_RATIO_0, Y_AT_YIELD_0), (X_AT_RATIO_0, Y_AT_YIELD_4)],
        fill=(255, 0, 255),
        width=3,
    )
    for point in POINTS:
        radius = 18
        draw.ellipse(
            [
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ],
            outline=(255, 0, 255),
            width=3,
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", type=Path, default=DEFAULT_RENDER)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    image = _verify_render(args.render)
    payload = csv_text()
    report = json.dumps(
        manifest(hashlib.sha256(payload.encode("utf-8")).hexdigest()),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.check:
        if CSV_PATH.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"stale digitization: {CSV_PATH}")
        if MANIFEST_PATH.read_text(encoding="utf-8") != report:
            raise RuntimeError(f"stale manifest: {MANIFEST_PATH}")
    else:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(payload, encoding="utf-8")
        MANIFEST_PATH.write_text(report, encoding="utf-8")
    if args.overlay:
        _overlay(image, args.overlay)
    print(
        f"verified {len(POINTS)} Figure 5.10 points; "
        f"yield bound ±{YIELD_DIGITIZATION_BOUND:.3f}"
    )


if __name__ == "__main__":
    main()
