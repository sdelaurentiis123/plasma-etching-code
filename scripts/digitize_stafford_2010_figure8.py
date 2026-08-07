#!/usr/bin/env python3
"""Reproduce Stafford et al. (2010) Figure 8 chlorine wall data.

Figure 8 collects spinning-wall measurements of the Cl-atom recombination
probability on plasma-conditioned stainless steel and anodized aluminum.  The
marker shape identifies pressure; fill identifies material.  Individual ICP
powers are not encoded, so this digitization retains only the published
100--600 W range and must not manufacture a per-point power value.

The resulting table is an evidence envelope for a state-dependent wall model,
not a scalar calibration target and not a Lam-reactor boundary.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "stafford_2010"
CSV_PATH = OUTPUT_DIRECTORY / "figure8_chlorine_wall_recombination.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_URL = (
    "https://iupac.org/publications/pac/pdf/2010/pdf/8206x1301.pdf")
SOURCE_PDF_SHA256 = (
    "2f74388576d435d9b3ae3843d5fa14f6a941ef61124b406f3ee7a7496e464b08")
RENDER_SHA256 = (
    "296bc9dfb1426272449d62694ba773bb57cbb692da8f78a1805a839cfb04e447")
RENDER_SIZE = (4500, 6450)
PDF_PAGE = 10
PRINT_PAGE = 1310
RENDER_DPI = 600

# Full-page vector-axis centers.  The x axis is linear from 0 to 0.9.
# The y axis is logarithmic; the 0.1 and 0.01 major ticks are 649 pixels
# apart.  Values below use the centers of the repeated vector strokes.
X_AT_ZERO = 1563.5
X_AT_POINT_NINE = 3175.0
Y_AT_POINT_ONE = 905.5
Y_AT_POINT_ZERO_ONE = 1554.5
HORIZONTAL_PIXEL_ALLOWANCE = 5.0
VERTICAL_PIXEL_ALLOWANCE = 8.0

# The full-page points were audited in a crop beginning at (1500, 650).
_CROP_X = 1500.0
_CROP_Y = 650.0


@dataclass(frozen=True)
class PixelPoint:
    material: str
    pressure_mTorr: float
    x_px: float
    y_px: float
    marker: str


def _point(
        material: str,
        pressure_mTorr: float,
        crop_x_px: float,
        crop_y_px: float,
        marker: str,
) -> PixelPoint:
    return PixelPoint(
        material=material,
        pressure_mTorr=pressure_mTorr,
        x_px=_CROP_X + crop_x_px,
        y_px=_CROP_Y + crop_y_px,
        marker=marker,
    )


# All 39 markers were identified visually on the original 600-dpi render.
# PIL/NumPy connected components supplied isolated centers.  Hough-circle
# replay separated the circle clusters crossed by the dashed eye guides.
# The two nearly overlapping AA squares at x ~= 0.78 were reconciled on a
# nearest-neighbor crop; the larger vertical allowance covers that overlap.
PIXEL_POINTS = (
    # 20 mTorr: inverted triangles.
    _point("anodized_aluminum", 20.0, 252.6, 935.3, "open_down_triangle"),
    _point("anodized_aluminum", 20.0, 282.6, 905.2, "open_down_triangle"),
    _point("anodized_aluminum", 20.0, 395.5, 905.3, "open_down_triangle"),
    _point("anodized_aluminum", 20.0, 450.4, 879.1, "open_down_triangle"),
    _point("stainless_steel", 20.0, 252.8, 1050.4, "filled_down_triangle"),
    _point("stainless_steel", 20.0, 282.8, 1084.2, "filled_down_triangle"),
    _point("stainless_steel", 20.0, 395.4, 1119.0, "filled_down_triangle"),
    _point("stainless_steel", 20.0, 450.4, 1164.1, "filled_down_triangle"),
    # 10 mTorr: upright triangles.
    _point("anodized_aluminum", 10.0, 494.2, 810.8, "open_up_triangle"),
    _point("anodized_aluminum", 10.0, 542.2, 790.4, "open_up_triangle"),
    _point("anodized_aluminum", 10.0, 689.6, 772.3, "open_up_triangle"),
    _point("anodized_aluminum", 10.0, 707.5, 742.5, "open_up_triangle"),
    _point("stainless_steel", 10.0, 493.4, 981.9, "filled_up_triangle"),
    _point("stainless_steel", 10.0, 547.0, 960.7, "filled_up_triangle"),
    _point("stainless_steel", 10.0, 689.8, 1001.0, "filled_up_triangle"),
    _point("stainless_steel", 10.0, 707.8, 967.3, "filled_up_triangle"),
    # 5 mTorr: circles.
    _point("anodized_aluminum", 5.0, 618.5, 772.5, "open_circle"),
    _point("anodized_aluminum", 5.0, 618.5, 810.5, "open_circle"),
    _point("anodized_aluminum", 5.0, 689.5, 638.5, "open_circle"),
    _point("anodized_aluminum", 5.0, 689.5, 709.5, "open_circle"),
    _point("anodized_aluminum", 5.0, 851.5, 615.5, "open_circle"),
    _point("anodized_aluminum", 5.0, 832.5, 629.5, "open_circle"),
    _point("anodized_aluminum", 5.0, 833.5, 659.5, "open_circle"),
    _point("stainless_steel", 5.0, 618.7, 992.8, "filled_circle"),
    _point("stainless_steel", 5.0, 689.7, 901.0, "filled_circle"),
    _point("stainless_steel", 5.0, 850.3, 798.9, "filled_circle"),
    _point("stainless_steel", 5.0, 833.5, 867.5, "filled_circle"),
    # 1.25 mTorr: squares.
    _point("anodized_aluminum", 1.25, 1100.3, 274.1, "open_square"),
    _point("anodized_aluminum", 1.25, 1100.2, 309.0, "open_square"),
    _point("anodized_aluminum", 1.25, 1296.4, 338.0, "open_square"),
    _point("anodized_aluminum", 1.25, 1296.2, 374.0, "open_square"),
    _point("anodized_aluminum", 1.25, 1242.9, 511.3, "open_square"),
    _point("anodized_aluminum", 1.25, 1242.9, 541.6, "open_square"),
    _point("anodized_aluminum", 1.25, 1459.5, 475.0, "open_square"),
    _point("anodized_aluminum", 1.25, 1459.5, 493.5, "open_square"),
    _point("stainless_steel", 1.25, 1101.0, 580.6, "filled_square"),
    _point("stainless_steel", 1.25, 1244.0, 762.5, "filled_square"),
    _point("stainless_steel", 1.25, 1298.8, 691.6, "filled_square"),
    _point("stainless_steel", 1.25, 1458.5, 738.0, "filled_square"),
)

FIELDNAMES = (
    "material",
    "pressure_mTorr",
    "reported_icp_power_range_W",
    "cl_to_cl2_density_ratio",
    "cl_recombination_probability",
    "marker_center_x_px",
    "marker_center_y_px",
    "marker",
    "observable_basis",
    "digitization_ratio_uncertainty",
    "digitization_log10_gamma_uncertainty",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def ratio_at_pixel(x_px: float) -> float:
    return 0.9 * (x_px - X_AT_ZERO) / (X_AT_POINT_NINE - X_AT_ZERO)


def gamma_at_pixel(y_px: float) -> float:
    pixels_per_decade = Y_AT_POINT_ZERO_ONE - Y_AT_POINT_ONE
    log10_gamma = -1.0 + (Y_AT_POINT_ONE - y_px) / pixels_per_decade
    return 10.0 ** log10_gamma


def rows() -> list[dict[str, str]]:
    ratio_uncertainty = (
        0.9 * HORIZONTAL_PIXEL_ALLOWANCE
        / (X_AT_POINT_NINE - X_AT_ZERO)
    )
    log_uncertainty = (
        VERTICAL_PIXEL_ALLOWANCE
        / (Y_AT_POINT_ZERO_ONE - Y_AT_POINT_ONE)
    )
    output = []
    for point in PIXEL_POINTS:
        output.append({
            "material": point.material,
            "pressure_mTorr": f"{point.pressure_mTorr:g}",
            "reported_icp_power_range_W": "100-600",
            "cl_to_cl2_density_ratio": f"{ratio_at_pixel(point.x_px):.6f}",
            "cl_recombination_probability": (
                f"{gamma_at_pixel(point.y_px):.7f}"
            ),
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            "observable_basis": "spinning_wall_LH_recombination",
            "digitization_ratio_uncertainty": f"{ratio_uncertainty:.6f}",
            "digitization_log10_gamma_uncertainty": (
                f"{log_uncertainty:.6f}"
            ),
        })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    decoded = rows()
    by_material = {
        material: [
            float(row["cl_recombination_probability"])
            for row in decoded
            if row["material"] == material
        ]
        for material in ("stainless_steel", "anodized_aluminum")
    }
    ratios = [
        float(row["cl_to_cl2_density_ratio"]) for row in decoded]
    return {
        "manifest_id": "STAFFORD-2010-FIG8-CL-WALL-RECOMBINATION-R1",
        "source": {
            "citation": (
                "L. Stafford, J. Guha, R. Khare, S. Mattei, O. Boudreault, "
                "B. Clain, and V. M. Donnelly, Pure and Applied Chemistry "
                "82, 1301-1315 (2010), Experimental and modeling study of "
                "O and Cl atoms surface recombination reactions in O2 and "
                "Cl2 plasmas"
            ),
            "doi": "10.1351/PAC-CON-09-11-02",
            "official_iupac_url": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": PDF_PAGE,
            "print_page": PRINT_PAGE,
            "figure": "Figure 8",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
            "figure_data_lineage": (
                "review figure reproduced from the authors' refs. 31 and 34"
            ),
        },
        "experiment": {
            "plasma": "inductively coupled Cl2",
            "pressure_range_mTorr": [1.25, 20.0],
            "icp_power_range_W": [100.0, 600.0],
            "surfaces": [
                "plasma-conditioned stainless steel",
                "plasma-conditioned anodized aluminum",
            ],
            "measurement": (
                "spinning-wall Langmuir-Hinshelwood recombination; net Cl2 "
                "desorption extrapolated to zero reaction time and divided "
                "by Cl-atom impingement flux"
            ),
            "individual_power_per_marker": "not encoded in Figure 8",
        },
        "pixel_calibration": {
            "x_at_0_ratio": X_AT_ZERO,
            "x_at_0.9_ratio": X_AT_POINT_NINE,
            "y_at_0.1_gamma": Y_AT_POINT_ONE,
            "y_at_0.01_gamma": Y_AT_POINT_ZERO_ONE,
            "x_scale": "linear",
            "y_scale": "log10",
        },
        "digitization": {
            "marker_count": len(PIXEL_POINTS),
            "method": (
                "official IUPAC PDF; 600-dpi Poppler full-page render; "
                "original-resolution human-vision reconciliation; PIL/NumPy "
                "dark-axis and marker replay; connected-component centers; "
                "Hough-circle separation where eye-guide dashes overlap"
            ),
            "horizontal_pixel_allowance": HORIZONTAL_PIXEL_ALLOWANCE,
            "vertical_pixel_allowance": VERTICAL_PIXEL_ALLOWANCE,
            "ratio_uncertainty": (
                0.9 * HORIZONTAL_PIXEL_ALLOWANCE
                / (X_AT_POINT_NINE - X_AT_ZERO)
            ),
            "log10_gamma_uncertainty": (
                VERTICAL_PIXEL_ALLOWANCE
                / (Y_AT_POINT_ZERO_ONE - Y_AT_POINT_ONE)
            ),
            "measurement_uncertainty": (
                "not assigned: Figure 8 has no error bars; digitization "
                "uncertainty is not experimental uncertainty"
            ),
        },
        "source_internal_checks": {
            "ratio_range": [min(ratios), max(ratios)],
            "stainless_steel_gamma_range": [
                min(by_material["stainless_steel"]),
                max(by_material["stainless_steel"]),
            ],
            "anodized_aluminum_gamma_range": [
                min(by_material["anodized_aluminum"]),
                max(by_material["anodized_aluminum"]),
            ],
            "paper_text_stainless_steel_range": [0.004, 0.03],
            "paper_text_anodized_aluminum_relative_statement": (
                "about a factor of 2 larger than stainless steel"
            ),
        },
        "claim_boundary": {
            "valid": [
                "conditioned-surface chlorine recombination evidence envelope",
                "dependence on Cl/Cl2 density ratio and surface condition",
                "wall-model validation target within the published domain",
            ],
            "not_valid": [
                "one constant gamma_Cl independent of plasma and wall state",
                "extrapolation outside the measured pressure, power, and ratio",
                "an individual power value for any Figure-8 marker",
                "a proprietary Lam reactor boundary",
                "a feature-depth calibration",
            ],
            "model_warning": (
                "Figure 9 of the paper shows that no one constant gamma_Cl "
                "adequately fits the full pressure-dependent dissociation set"
            ),
        },
        "output": {
            "path": (
                "data/experimental/stafford_2010/"
                "figure8_chlorine_wall_recombination.csv"
            ),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def _render_page(pdf_path: Path) -> Image.Image:
    if _sha256(pdf_path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Stafford source PDF checksum does not match")
    with tempfile.TemporaryDirectory(
            prefix="petch-stafford-vision-") as directory:
        prefix = Path(directory) / "page10_600dpi"
        subprocess.run(
            [
                "pdftoppm",
                "-f", str(PDF_PAGE),
                "-l", str(PDF_PAGE),
                "-r", str(RENDER_DPI),
                "-png",
                "-singlefile",
                str(pdf_path),
                str(prefix),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        render_path = prefix.with_suffix(".png")
        if _sha256(render_path) != RENDER_SHA256:
            raise RuntimeError(
                "Stafford 600-dpi page render checksum does not match")
        image = Image.open(render_path).convert("RGB")
        image.load()
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected render size {image.size}; expected {RENDER_SIZE}")
    return image


def _vision_audit(image: Image.Image) -> None:
    gray = np.asarray(image.convert("L"))
    dark = gray < 120
    checks = {
        "left axis": np.mean(
            np.min(gray[710:2012, 1559:1569], axis=1) < 150),
        "bottom axis": np.mean(
            np.min(gray[2005:2015, 1560:3178], axis=0) < 150),
        "top axis": np.mean(
            np.min(gray[708:718, 1560:3178], axis=0) < 150),
        "right axis": np.mean(
            np.min(gray[710:2012, 3170:3180], axis=1) < 150),
    }
    failed = {
        name: float(value)
        for name, value in checks.items()
        if value < 0.70
    }
    if failed:
        raise RuntimeError(f"axis dark-pixel verification failed: {failed}")

    for point in PIXEL_POINTS:
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        half_width = 18
        count = int(np.count_nonzero(
            dark[
                y - half_width:y + half_width + 1,
                x - half_width:x + half_width + 1,
            ]
        ))
        minimum = 120 if point.material == "anodized_aluminum" else 450
        if count < minimum:
            raise RuntimeError(
                f"failed to recover {point.material} {point.marker} marker "
                f"near ({point.x_px:.1f}, {point.y_px:.1f}): {count} pixels"
            )


def _write_overlay(image: Image.Image, path: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for index, point in enumerate(PIXEL_POINTS, start=1):
        radius = 7
        box = (
            point.x_px - radius,
            point.y_px - radius,
            point.x_px + radius,
            point.y_px + radius,
        )
        draw.ellipse(box, outline=(255, 0, 0), width=3)
        draw.text(
            (point.x_px + 9, point.y_px - 9),
            str(index),
            fill=(255, 0, 0),
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(path)


def write_committed_files() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def check_committed_files() -> None:
    expected_csv = csv_text()
    expected_manifest = manifest_text(expected_csv)
    if CSV_PATH.read_text(encoding="utf-8") != expected_csv:
        raise RuntimeError(f"{CSV_PATH} is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != expected_manifest:
        raise RuntimeError(f"{MANIFEST_PATH} is stale")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--overlay", type=Path)
    arguments = parser.parse_args()
    if arguments.write:
        write_committed_files()
    if arguments.check:
        check_committed_files()
    if arguments.source_pdf is not None:
        image = _render_page(arguments.source_pdf)
        _vision_audit(image)
        if arguments.overlay is not None:
            _write_overlay(image, arguments.overlay)
    elif arguments.overlay is not None:
        parser.error("--overlay requires --source-pdf")
    if not (
        arguments.write
        or arguments.check
        or arguments.source_pdf is not None
    ):
        parser.error("choose --write, --check, or --source-pdf")


if __name__ == "__main__":
    main()
