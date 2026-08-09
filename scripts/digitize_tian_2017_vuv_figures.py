#!/usr/bin/env python3
"""Digitize Tian (2017) Figures 5.10(a,c) and 5.12(c).

The figures are *equipment-model outputs*, not measurements.  They provide an
independent, wavelength-resolved regression target for the deterministic
Ar/Cl2 reactor/radiation stack: substrate flux, photon/ion ratio, spectral
composition, and line trapping.  They cannot validate absolute depth.

Pixel coordinates below were isolated from a 300-dpi Poppler render and are
replayed against checksum-pinned full-page images with PIL/NumPy.  No curve
interpolation is presented as source data; only visible square markers are
retained.
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
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
CSV_PATH = OUTPUT_DIRECTORY / "digitized_figures_5_10_5_12.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_URL = (
    "https://cpseg.eecs.umich.edu/pub/theses/tian_peng_phd_thesis.pdf")
SOURCE_PDF_SHA256 = (
    "4d260ab9e85240bd051ccb3ba32cd047b4ac1ddb6c309dbbca4093822e37790b")
RENDER_DPI = 300
RENDER_SIZE = (2550, 3300)
RENDER_SHA256 = {
    171: "c0fc7cc6f06aef081ce46f81116bf06963d8886d20754cbad85dbb0aba0da63c",
    173: "21f084f3cd7c309db86dc656be68c2ff991481aac56f990b80847c15a7f3c863",
}


@dataclass(frozen=True)
class PixelPoint:
    figure: str
    panel: str
    series: str
    cl2_fraction_percent: float
    x_px: float
    y_px: float
    marker_state: str = "isolated"


def _points(
    figure: str,
    panel: str,
    series: str,
    fractions: tuple[float, ...],
    xs: tuple[float, ...],
    ys: tuple[float, ...],
    *,
    overlaps: tuple[float, ...] = (),
) -> tuple[PixelPoint, ...]:
    if not (len(fractions) == len(xs) == len(ys)):
        raise ValueError("pixel-coordinate arrays must have equal lengths")
    return tuple(
        PixelPoint(
            figure,
            panel,
            series,
            fraction,
            x,
            y,
            "overlapping_markers" if fraction in overlaps else "isolated",
        )
        for fraction, x, y in zip(fractions, xs, ys)
    )


FRACTIONS = (5., 10., 15., 20., 30., 40., 50., 60., 80., 90., 95.)
FIG510_AC_X = (
    1037., 1073.5, 1109.5, 1146., 1219., 1291.5,
    1364.5, 1437.5, 1583.5, 1656.5, 1693.5,
)
FIG510_C_X = (
    1032.4, 1063.5, 1095., 1126.5, 1189.5, 1252.,
    1315.4, 1378.5, 1504.5, 1567.5, 1598.9,
)
FIG512_C_X = (
    1018., 1054.5, 1090.5, 1127., 1200., 1273.,
    1346., 1419., 1565., 1638., 1674.,
)

PIXEL_POINTS = (
    *_points(
        "5.10", "a", "total_photon_flux", FRACTIONS, FIG510_AC_X,
        (299.5, 363.5, 413.5, 457.5, 506.5, 549.5,
         576.5, 601.5, 647.5, 683., 688.5)),
    *_points(
        "5.10", "a", "Ar_104.8_nm_photon_flux",
        FRACTIONS[:-1], FIG510_AC_X[:-1],
        (394., 455., 502., 542.3, 582.9, 619.6,
         628.5, 666.5, 730.1, 811.5),
        overlaps=(90.,)),
    *_points(
        "5.10", "a", "Ar_106.7_nm_photon_flux", FRACTIONS, FIG510_AC_X,
        (322.9, 388.8, 445., 494., 553.1, 604.7,
         666.6, 680.4, 766., 811.5, 841.),
        overlaps=(90.,)),
    *_points(
        "5.10", "a", "Cl_139_nm_photon_flux", FRACTIONS, FIG510_AC_X,
        (652., 627., 628.5, 632.6, 658.5, 680.,
         683.4, 701.5, 699.6, 699.7, 702.9)),
    *_points(
        "5.10", "c", "total_photon_to_total_ion_flux_ratio_beta",
        FRACTIONS, FIG510_C_X,
        (1889.5, 2107.5, 2210.5, 2273.6, 2326.4, 2356.7,
         2389.5, 2399.5, 2421.5, 2433., 2433.4)),
    *_points(
        "5.10", "c", "Cl_139_nm_fraction_of_three_line_spectrum",
        FRACTIONS[1:], FIG510_C_X[1:],
        (2441., 2423.5, 2399., 2376., 2356.7,
         2305.1, 2292., 2134., 1936.5, 1903.5)),
    *_points(
        "5.12", "c", "Ar_106.7_nm_trapping_factor", FRACTIONS,
        FIG512_C_X,
        (1834., 1905.9, 1976., 2044., 2128., 2210.,
         2228., 2303., 2385., 2446., 2465.)),
    *_points(
        "5.12", "c", "Ar_104.8_nm_trapping_factor", FRACTIONS,
        FIG512_C_X,
        (2201.4, 2204.5, 2212.5, 2220.5, 2233.5, 2249.5,
         2253.5, 2274.5, 2312.5, 2369.5, 2418.5)),
    *_points(
        "5.12", "c", "Cl_139_nm_trapping_factor", FRACTIONS,
        FIG512_C_X,
        (2431., 2395., 2370.9, 2358.9, 2360., 2362.,
         2365., 2371., 2359., 2346.5, 2346.5)),
)

FIELDNAMES = (
    "figure", "panel", "series", "cl2_fraction_percent", "value",
    "units", "marker_center_x_px", "marker_center_y_px", "marker_state",
    "digitization_vertical_uncertainty_px", "evidence_type",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _value_and_units(point: PixelPoint) -> tuple[float, str]:
    if point.figure == "5.10" and point.panel == "a":
        # Printed ticks: y=250 -> 10^2 and y=916 -> 10^-2 in units of
        # 1e14 cm^-2 s^-1.
        log10_plot_value = 2.0 - 4.0 * (point.y_px - 250.0) / (916.0 - 250.0)
        return 1.0e14 * 10.0 ** log10_plot_value, "cm^-2 s^-1"
    if point.figure == "5.10" and point.panel == "c":
        # Left beta axis is 0--0.5; the right spectral-fraction axis is 0--1.
        span = 2458.0 - 1793.0
        normalized = (2458.0 - point.y_px) / span
        if point.series.endswith("ratio_beta"):
            return 0.5 * normalized, "dimensionless"
        return normalized, "dimensionless"
    if point.figure == "5.12" and point.panel == "c":
        # Printed ticks are 0 at y=2480 and 100 at y=2359.  The plotted Cl
        # 139-nm series is explicitly multiplied by ten in the source.
        plotted = 100.0 * (2480.0 - point.y_px) / (2480.0 - 2359.0)
        if point.series == "Cl_139_nm_trapping_factor":
            plotted /= 10.0
        return plotted, "dimensionless"
    raise ValueError(f"unknown figure/panel {point.figure}{point.panel}")


def rows() -> list[dict[str, str]]:
    output = []
    for point in PIXEL_POINTS:
        value, units = _value_and_units(point)
        output.append({
            "figure": point.figure,
            "panel": point.panel,
            "series": point.series,
            "cl2_fraction_percent": f"{point.cl2_fraction_percent:g}",
            "value": f"{value:.9g}",
            "units": units,
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "marker_state": point.marker_state,
            "digitization_vertical_uncertainty_px": "4",
            "evidence_type": "source_equipment_model_digitized",
        })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _series_values(name: str) -> list[float]:
    return [
        float(row["value"]) for row in rows() if row["series"] == name
    ]


def manifest(csv_sha256: str) -> dict[str, object]:
    total = _series_values("total_photon_flux")
    beta = _series_values("total_photon_to_total_ion_flux_ratio_beta")
    trap_1067 = _series_values("Ar_106.7_nm_trapping_factor")
    trap_1048 = _series_values("Ar_104.8_nm_trapping_factor")
    trap_139 = _series_values("Cl_139_nm_trapping_factor")
    return {
        "manifest_id": "TIAN-2017-VUV-FIG5.10-FIG5.12-R1",
        "source": {
            "citation": (
                "P. Tian, Controlling Photon and Ion Fluxes in Low Pressure "
                "Low Temperature Plasmas, PhD dissertation, University of "
                "Michigan (2017; published 2018)"),
            "official_url": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_pages": [171, 173],
            "print_pages": [153, 155],
            "figures": ["5.10(a,c)", "5.12(c)"],
            "render_dpi": RENDER_DPI,
            "render_size_px": list(RENDER_SIZE),
            "render_sha256_by_pdf_page": {
                str(page): digest for page, digest in RENDER_SHA256.items()
            },
        },
        "experiment_context": {
            "evidence_type": "source_equipment_model_digitized",
            "not_measurement": True,
            "reactor": "HPEM Ar/Cl2 ICP, 22.5 cm diameter, 12 cm height",
            "pressure_mTorr": 20.0,
            "flow_sccm": 200.0,
            "frequency_MHz": 10.0,
            "continuous_wave_power_W": 150.0,
            "mixture_sweep": "Ar/Cl2 = 95/5 through 5/95",
        },
        "pixel_calibration": {
            "figure_5_10_a": {
                "x_at_0_percent": 1000.0,
                "x_at_100_percent": 1730.0,
                "y_at_1e16_cm-2_s-1": 250.0,
                "y_at_1e12_cm-2_s-1": 916.0,
                "y_scale": "log10",
            },
            "figure_5_10_c": {
                "x_at_0_percent": 1000.0,
                "x_at_100_percent": 1631.0,
                "left_y_top_beta": 0.5,
                "right_y_top_fraction": 1.0,
                "y_top_px": 1793.0,
                "y_bottom_px": 2458.0,
            },
            "figure_5_12_c": {
                "x_at_0_percent": 981.0,
                "x_at_100_percent": 1711.0,
                "y_at_0": 2480.0,
                "y_at_100": 2359.0,
                "cl_139_printed_multiplier": 10.0,
            },
        },
        "digitization": {
            "marker_count": len(PIXEL_POINTS),
            "method": (
                "300-dpi Poppler full-page render; full-resolution visual "
                "inspection; PIL/NumPy dark-axis and color/marker-component "
                "isolation; printed-axis replay; optional QA overlay"),
            "vertical_pixel_uncertainty": 4.0,
            "measurement_uncertainty": (
                "not applicable: these curves are source equipment-model "
                "outputs; the pixel allowance is digitization uncertainty"),
        },
        "source_internal_checks": {
            "total_photon_flux_decreases_across_mixture_sweep": all(
                left > right for left, right in zip(total, total[1:])),
            "beta_decreases_across_mixture_sweep": all(
                left > right for left, right in zip(beta, beta[1:])),
            "ar_106_7_trapping_factor_endpoint_ratio": (
                trap_1067[0] / trap_1067[-1]),
            "ar_104_8_trapping_factor_endpoint_ratio": (
                trap_1048[0] / trap_1048[-1]),
            "cl_139_trapping_factor_min_max": [min(trap_139), max(trap_139)],
        },
        "claim_boundary": {
            "valid": [
                "regression target for wavelength-resolved reactor output",
                "regression target for deterministic line escape/trapping",
                "mixture-trend comparison for VUV/ion ratio",
            ],
            "not_valid": [
                "experimental validation",
                "Lam chamber calibration",
                "surface photo-etch yield",
                "feature-depth calibration or validation",
                "atomic branching data independent of the HPEM mechanism",
            ],
        },
        "output": {
            "path": (
                "data/experimental/tian_2017_vuv_figures/"
                "digitized_figures_5_10_5_12.csv"),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def write_committed_files() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Tian digitization CSV is not reproduced")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Tian digitization manifest is not reproduced")


def _render_pages(pdf_path: Path) -> dict[int, Image.Image]:
    if _sha256(pdf_path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Tian source PDF checksum does not match")
    output: dict[int, Image.Image] = {}
    with tempfile.TemporaryDirectory(prefix="petch-tian-vuv-vision-") as directory:
        for page in RENDER_SHA256:
            prefix = Path(directory) / f"page{page}_300dpi"
            subprocess.run(
                [
                    "pdftoppm", "-f", str(page), "-l", str(page),
                    "-r", str(RENDER_DPI), "-png", "-singlefile",
                    str(pdf_path), str(prefix),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            path = prefix.with_suffix(".png")
            if _sha256(path) != RENDER_SHA256[page]:
                raise RuntimeError(f"Tian page {page} render checksum mismatch")
            image = Image.open(path).convert("RGB")
            image.load()
            if image.size != RENDER_SIZE:
                raise RuntimeError(f"unexpected Tian render size {image.size}")
            output[page] = image
    return output


def _vision_audit(images: dict[int, Image.Image]) -> None:
    axes = {
        171: [
            (1000, 1730, 250, 916),
            (1000, 1730, 1022, 1687),
            (1000, 1631, 1793, 2458),
        ],
        173: [(981, 1711, 1814, 2480)],
    }
    for page, boxes in axes.items():
        gray = np.asarray(images[page].convert("L"))
        for x0, x1, y0, y1 in boxes:
            checks = (
                np.mean(np.min(gray[y0:y1 + 1, x0 - 2:x0 + 3], axis=1) < 150),
                np.mean(np.min(gray[y0:y1 + 1, x1 - 2:x1 + 3], axis=1) < 150),
                np.mean(np.min(gray[y0 - 2:y0 + 3, x0:x1 + 1], axis=0) < 150),
                np.mean(np.min(gray[y1 - 2:y1 + 3, x0:x1 + 1], axis=0) < 150),
            )
            if min(checks) < 0.90:
                raise RuntimeError(
                    f"Tian page {page} axis audit failed for {x0,x1,y0,y1}: "
                    f"{checks}")

    for point in PIXEL_POINTS:
        page = 171 if point.figure == "5.10" else 173
        gray = np.asarray(images[page].convert("L"))
        x = int(round(point.x_px))
        y = int(round(point.y_px))
        dark_pixels = int(np.count_nonzero(
            gray[y - 10:y + 11, x - 10:x + 11] < 170))
        minimum = 90 if point.marker_state == "overlapping_markers" else 120
        if dark_pixels < minimum:
            raise RuntimeError(
                f"failed to recover {point.figure}{point.panel} "
                f"{point.series} marker at {point.cl2_fraction_percent:g}%: "
                f"{dark_pixels} dark pixels")


def draw_overlays(images: dict[int, Image.Image], directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    colors = {
        "total_photon_flux": "#ff00ff",
        "Ar_104.8_nm_photon_flux": "#00ffff",
        "Ar_106.7_nm_photon_flux": "#ff7f00",
        "Cl_139_nm_photon_flux": "#7fff00",
        "total_photon_to_total_ion_flux_ratio_beta": "#ff00ff",
        "Cl_139_nm_fraction_of_three_line_spectrum": "#00ffff",
        "Ar_106.7_nm_trapping_factor": "#ff7f00",
        "Ar_104.8_nm_trapping_factor": "#ff00ff",
        "Cl_139_nm_trapping_factor": "#00ffff",
    }
    for page, image in images.items():
        overlay = image.copy()
        draw = ImageDraw.Draw(overlay)
        for point in PIXEL_POINTS:
            point_page = 171 if point.figure == "5.10" else 173
            if point_page != page:
                continue
            radius = 13
            draw.ellipse(
                (point.x_px - radius, point.y_px - radius,
                 point.x_px + radius, point.y_px + radius),
                outline=colors[point.series], width=3)
        overlay.save(directory / f"tian_page{page}_qa.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--overlay-dir", type=Path)
    args = parser.parse_args()

    if args.write:
        write_committed_files()
    verify_committed_files()
    if args.source_pdf is not None:
        images = _render_pages(args.source_pdf)
        _vision_audit(images)
        if args.overlay_dir is not None:
            draw_overlays(images, args.overlay_dir)
    print(json.dumps(manifest(hashlib.sha256(csv_text().encode()).hexdigest()), indent=2))


if __name__ == "__main__":
    main()
