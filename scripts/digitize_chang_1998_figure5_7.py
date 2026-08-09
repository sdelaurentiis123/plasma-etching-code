#!/usr/bin/env python3
"""Reproduce the Chang Figure 5.7 Cl2+ yield-line digitization.

The source prints the fitted Cl2+ slope (0.22 Si/ion/sqrt(eV)) but not its
intercept.  This audit renders the checksum-bound thesis page at 600 dpi,
verifies the affine plot calibration in the original pixels, and recovers the
printed fit-line intercept with a narrow dark-pixel corridor.  No feature etch
rate or reactor result enters the extraction.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PDF = ROOT / "research_sources" / "chang_thesis.pdf"
DEFAULT_OUTPUT = ROOT / "data" / "experimental" / "chang_1998_figure5_7"
SOURCE_PDF_SHA256 = (
    "ef5c511a7fd8e6d6a0bf721874e57ba18194184bf3db90457188ed42c4bd3b4b"
)
RENDER_SHA256 = (
    "34c7bf8df4fcc77ec6f87fd682639b54e7a2538271da49b397cf9bb34ed9460a"
)
PDF_PAGE = 113
PRINT_PAGE = 113
RENDER_DPI = 600
RENDER_SIZE = (5100, 6600)

# Axis intersections located on the original-resolution render.  Each pair is
# separated by the full printed axis range, which minimizes calibration noise.
X_ZERO_PX = 1684.0
X_FORTY_PX = 3522.0
Y_ZERO_PX = 4247.0
Y_FOUR_PX = 2571.0

# The long Hough segment on the Cl2+ printed regression line seeds a two-pixel
# dark-pixel corridor.  The final line is refit from all source pixels in that
# corridor; these endpoints are not converted directly into physics values.
SEED_LINE = (1926.0, 4231.0, 2677.0, 2719.0)
FIT_X_RANGE_PX = (1920, 2765)
FIT_Y_RANGE_PX = (2550, 4240)
CORRIDOR_HALF_WIDTH_PX = 2.0
PRINTED_SLOPE = 0.22
PRINTED_SLOPE_ROUNDING_HALF_WIDTH = 0.005
VALID_ENERGY_DOMAIN_EV = (26.0, 625.0)


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _render(pdf: Path, directory: Path) -> Path:
    prefix = directory / "chang_figure5_7_page113"
    subprocess.run(
        [
            "pdftoppm", "-f", str(PDF_PAGE), "-l", str(PDF_PAGE),
            "-r", str(RENDER_DPI), "-png", "-singlefile", str(pdf),
            str(prefix),
        ],
        check=True,
    )
    return prefix.with_suffix(".png")


def _dark_fraction_near_axis(gray: np.ndarray, *, horizontal: bool, pixel: int):
    if horizontal:
        region = gray[pixel - 2:pixel + 3, int(X_ZERO_PX):int(X_FORTY_PX) + 1]
    else:
        region = gray[int(Y_FOUR_PX):int(Y_ZERO_PX) + 1, pixel - 2:pixel + 3]
    return float(np.mean(region < 160))


def _extract(render: Path):
    if _hash(render) != RENDER_SHA256:
        raise ValueError("600-dpi Figure 5.7 page render checksum changed")
    image = Image.open(render).convert("L")
    if image.size != RENDER_SIZE:
        raise ValueError(f"unexpected rendered page size: {image.size}")
    gray = np.asarray(image)
    axis_darkness = {
        "x_axis_y0": _dark_fraction_near_axis(
            gray, horizontal=True, pixel=int(Y_ZERO_PX)),
        "x_axis_y4": _dark_fraction_near_axis(
            gray, horizontal=True, pixel=int(Y_FOUR_PX)),
        "y_axis_x0": _dark_fraction_near_axis(
            gray, horizontal=False, pixel=int(X_ZERO_PX)),
        "y_axis_x40": _dark_fraction_near_axis(
            gray, horizontal=False, pixel=int(X_FORTY_PX)),
    }
    if min(axis_darkness.values()) < 0.10:
        raise ValueError(f"axis pixel audit failed: {axis_darkness}")

    x1, y1, x2, y2 = SEED_LINE
    seed_slope = (y2 - y1) / (x2 - x1)
    seed_intercept = y1 - seed_slope * x1
    pixel_y, pixel_x = np.nonzero(gray < 120)
    selected = (
        (pixel_x >= FIT_X_RANGE_PX[0])
        & (pixel_x <= FIT_X_RANGE_PX[1])
        & (pixel_y >= FIT_Y_RANGE_PX[0])
        & (pixel_y <= FIT_Y_RANGE_PX[1])
        & (
            np.abs(pixel_y - (seed_slope * pixel_x + seed_intercept))
            <= CORRIDOR_HALF_WIDTH_PX
        )
    )
    count = int(np.count_nonzero(selected))
    if count < 1000:
        raise ValueError(f"too few fit-line pixels survived: {count}")
    pixel_slope, pixel_intercept = np.polyfit(
        pixel_x[selected], pixel_y[selected], 1)

    x_pixels_per_sqrt_eV = (X_FORTY_PX - X_ZERO_PX) / 40.0
    y_pixels_per_yield = (Y_ZERO_PX - Y_FOUR_PX) / 4.0
    digitized_slope = (
        -float(pixel_slope) * x_pixels_per_sqrt_eV / y_pixels_per_yield)
    x_intercept_px = (Y_ZERO_PX - float(pixel_intercept)) / float(pixel_slope)
    threshold_sqrt_eV = (
        (x_intercept_px - X_ZERO_PX) / x_pixels_per_sqrt_eV)
    threshold_eV = threshold_sqrt_eV ** 2

    if abs(digitized_slope - PRINTED_SLOPE) > PRINTED_SLOPE_ROUNDING_HALF_WIDTH:
        raise ValueError(
            "pixel-recovered slope does not reconcile with printed slope: "
            f"{digitized_slope} versus {PRINTED_SLOPE}")
    return {
        "axis_dark_fraction": axis_darkness,
        "fit_pixel_count": count,
        "pixel_line_slope": float(pixel_slope),
        "pixel_line_intercept": float(pixel_intercept),
        "digitized_slope_si_per_ion_per_sqrt_eV": digitized_slope,
        "threshold_sqrt_eV": threshold_sqrt_eV,
        "threshold_energy_eV": threshold_eV,
        "x_intercept_px": x_intercept_px,
    }


def _write_outputs(output: Path, extracted: dict):
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "cl2plus_poly_si_yield_fit.csv"
    fields = (
        "source_curve", "yield_law", "printed_slope_si_per_ion_per_sqrt_eV",
        "digitized_slope_si_per_ion_per_sqrt_eV", "threshold_sqrt_eV",
        "threshold_energy_eV", "minimum_valid_energy_eV",
        "maximum_valid_energy_eV", "source_surface_state",
        "supports_chlorinated_poly_si_prediction",
    )
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerow({
            "source_curve": "Chang thesis Figure 5.7 Cl2+ (Balooch)",
            "yield_law": "Y=0.22*max(sqrt(E_eV)-sqrt(Eth_eV),0)",
            "printed_slope_si_per_ion_per_sqrt_eV": PRINTED_SLOPE,
            "digitized_slope_si_per_ion_per_sqrt_eV": (
                extracted["digitized_slope_si_per_ion_per_sqrt_eV"]),
            "threshold_sqrt_eV": extracted["threshold_sqrt_eV"],
            "threshold_energy_eV": extracted["threshold_energy_eV"],
            "minimum_valid_energy_eV": VALID_ENERGY_DOMAIN_EV[0],
            "maximum_valid_energy_eV": VALID_ENERGY_DOMAIN_EV[1],
            "source_surface_state": (
                "Balooch approximately 1e-4 Torr chlorine background; "
                "Chang interprets surface as highly chlorinated"),
            "supports_chlorinated_poly_si_prediction": True,
        })

    manifest = {
        "schema": "petch.chang-1998-figure5.7-cl2plus-yield.v1",
        "source": {
            "citation": (
                "J. P. Chang, Study of plasma-surface kinetics and simulation "
                "of feature profile evolution in chlorine etching of patterned "
                "polysilicon, MIT PhD thesis (1998), Figure 5.7"
            ),
            "underlying_cl2plus_data": "Balooch et al. (1996), as plotted by Chang",
            "pdf_path": "research_sources/chang_thesis.pdf",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": PDF_PAGE,
            "print_page": PRINT_PAGE,
            "render_dpi": RENDER_DPI,
            "render_size_px": list(RENDER_SIZE),
            "render_sha256": RENDER_SHA256,
        },
        "pixel_calibration": {
            "sqrt_energy_zero_x_px": X_ZERO_PX,
            "sqrt_energy_40_x_px": X_FORTY_PX,
            "yield_zero_y_px": Y_ZERO_PX,
            "yield_four_y_px": Y_FOUR_PX,
            "axis_dark_fraction": extracted["axis_dark_fraction"],
        },
        "fit_line_extraction": {
            "method": (
                "original-resolution visual identification; long-line Hough seed; "
                "PIL/NumPy dark-pixel ordinary-least-squares replay inside a "
                "two-pixel corridor"
            ),
            "seed_line_xyxy_px": list(SEED_LINE),
            "fit_x_range_px": list(FIT_X_RANGE_PX),
            "fit_y_range_px": list(FIT_Y_RANGE_PX),
            "corridor_half_width_px": CORRIDOR_HALF_WIDTH_PX,
            "fit_pixel_count": extracted["fit_pixel_count"],
            "pixel_line_slope": extracted["pixel_line_slope"],
            "pixel_line_intercept": extracted["pixel_line_intercept"],
            "x_intercept_px": extracted["x_intercept_px"],
        },
        "reconciled_yield_law": {
            "law": "Y_Si_per_Cl2plus=0.22*max(sqrt(E_eV)-sqrt(25.999_eV),0)",
            "printed_slope_si_per_ion_per_sqrt_eV": PRINTED_SLOPE,
            "printed_slope_rounding_half_width": PRINTED_SLOPE_ROUNDING_HALF_WIDTH,
            "pixel_recovered_slope_si_per_ion_per_sqrt_eV": (
                extracted["digitized_slope_si_per_ion_per_sqrt_eV"]),
            "pixel_recovered_threshold_sqrt_eV": extracted["threshold_sqrt_eV"],
            "pixel_recovered_threshold_energy_eV": extracted["threshold_energy_eV"],
            "declared_threshold_digitization_bound_eV": 2.0,
            "valid_energy_domain_eV": list(VALID_ENERGY_DOMAIN_EV),
        },
        "scope": {
            "surface_state": (
                "highly chlorinated poly-Si; source text attributes the elevated "
                "Cl2+ curve to approximately 1e-4 Torr chlorine background"
            ),
            "valid_use": (
                "species-resolved Cl2+ removal yield on a highly chlorinated "
                "poly-Si surface inside the plotted energy support"
            ),
            "invalid_use": (
                "bare-Si sputtering, oxygenated/passivated surfaces, incidence-angle "
                "response, or extrapolation outside the plotted energy support"
            ),
            "feature_depth_used": False,
            "reactor_state_used": False,
        },
        "output": {
            "csv_path": (
                "data/experimental/chang_1998_figure5_7/"
                "cl2plus_poly_si_yield_fit.csv"),
            "csv_sha256": _hash(csv_path),
        },
    }
    manifest_path = output / "digitization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return csv_path, manifest_path


def _write_overlay(render: Path, path: Path, extracted: dict):
    image = Image.open(render).convert("RGB")
    draw = ImageDraw.Draw(image)
    draw.line(
        [(X_ZERO_PX, Y_ZERO_PX), (X_FORTY_PX, Y_ZERO_PX)],
        fill=(0, 180, 0), width=8,
    )
    draw.line(
        [(X_ZERO_PX, Y_ZERO_PX), (X_ZERO_PX, Y_FOUR_PX)],
        fill=(0, 180, 0), width=8,
    )
    m = extracted["pixel_line_slope"]
    b = extracted["pixel_line_intercept"]
    xlo, xhi = FIT_X_RANGE_PX
    draw.line(
        [(xlo, m * xlo + b), (xhi, m * xhi + b)],
        fill=(220, 0, 0), width=10,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    image.crop((1100, 2350, 3650, 4700)).save(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if _hash(args.source_pdf) != SOURCE_PDF_SHA256:
        raise ValueError("Chang thesis checksum changed")
    with tempfile.TemporaryDirectory(prefix="petch-chang57-") as temporary:
        render = _render(args.source_pdf, Path(temporary))
        extracted = _extract(render)
        if args.overlay is not None:
            _write_overlay(render, args.overlay, extracted)
    if args.check:
        expected_csv = args.output_dir / "cl2plus_poly_si_yield_fit.csv"
        expected_manifest = args.output_dir / "digitization_manifest.json"
        before = (
            _hash(expected_csv) if expected_csv.exists() else None,
            _hash(expected_manifest) if expected_manifest.exists() else None,
        )
        with tempfile.TemporaryDirectory(prefix="petch-chang57-output-") as temporary:
            csv_path, manifest_path = _write_outputs(Path(temporary), extracted)
            after = (_hash(csv_path), _hash(manifest_path))
        if before != after:
            raise ValueError(f"committed digitization outputs changed: {before} != {after}")
    else:
        _write_outputs(args.output_dir, extracted)
    print(json.dumps(extracted, sort_keys=True))


if __name__ == "__main__":
    main()
