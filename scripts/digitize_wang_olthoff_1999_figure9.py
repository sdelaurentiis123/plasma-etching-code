#!/usr/bin/env python3
"""Reproduce Wang and Olthoff (1999) Figure 9 chlorine ion fluxes.

The source reports absolute total current through a sampling orifice and
mass-resolved IED intensities scaled to that total.  This script preserves
that distinction, checks the official NIST PDF and a 600-dpi page render,
replays the logarithmic axes, and emits an optional full-page QA overlay.

The experiment is an independent reactor/sheath validation target.  It is not
an electron-impact branching measurement and must not be used to choose a
Cl2 ionization branching fraction.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = (
    ROOT / "data" / "experimental" / "wang_olthoff_1999")
CSV_PATH = OUTPUT_DIRECTORY / "figure9_chlorine_ion_flux.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_URL = "https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=1506"
SOURCE_PDF_SHA256 = (
    "17702f0ffb904cca42760867c693c382b174c68f02ef7b5064ec95933adea460")
RENDER_SHA256 = (
    "44a1953b90353711dfb6798b357ecb90eb0fbc10abb1c49d725a21a39a54b12a")
RENDER_SIZE = (5096, 6716)
PDF_PAGE = 6
RENDER_DPI = 600

# Full-page axes.  Repeated four-pixel vector strokes put both panels at the
# same x positions.  Major log-y ticks are separated by 404 pixels in panel
# (a) and 200 pixels in panel (b).
X_AT_0_PA = 3073.5
X_AT_10_PA = 4437.5
PANEL_A_Y_AT_10 = 681.5
PANEL_A_Y_AT_1 = 1085.5
PANEL_B_Y_AT_10 = 1981.5
PANEL_B_Y_AT_1 = 2181.5


@dataclass(frozen=True)
class PixelPoint:
    panel: str
    feed: str
    pressure_Pa: float
    species: str
    x_px: float
    y_px: float
    marker: str
    observable_basis: str


# Marker centers were isolated at the full 600-dpi resolution.  Nominal
# pressures come from the printed experimental setpoints; x pixels are retained
# as an independent placement check.  The top-panel Cl+ and total markers
# overlap increasingly with pressure, so their distinct square/circle centers
# were reconciled against a 4x nearest-neighbor crop before being frozen.
PIXEL_POINTS = (
    PixelPoint(
        "a", "100% Cl2", 2.7, "total_positive_ion", 3442.0, 733.0,
        "open_circle", "direct_total_current"),
    PixelPoint(
        "a", "100% Cl2", 4.0, "total_positive_ion", 3620.0, 674.0,
        "open_circle", "direct_total_current"),
    PixelPoint(
        "a", "100% Cl2", 6.7, "total_positive_ion", 3987.0, 588.0,
        "open_circle", "direct_total_current"),
    PixelPoint(
        "a", "100% Cl2", 2.7, "Cl+", 3442.0, 765.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "a", "100% Cl2", 4.0, "Cl+", 3620.0, 690.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "a", "100% Cl2", 6.7, "Cl+", 3987.0, 596.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "a", "100% Cl2", 2.7, "Cl2+", 3444.0, 1076.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "a", "100% Cl2", 4.0, "Cl2+", 3621.0, 1141.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "a", "100% Cl2", 6.7, "Cl2+", 3987.0, 1289.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 0.67, "total_positive_ion", 3170.0,
        1954.0, "open_circle", "direct_total_current"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 1.3, "total_positive_ion", 3257.0,
        1944.0, "open_circle", "direct_total_current"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 2.7, "total_positive_ion", 3446.0,
        1920.0, "open_circle", "direct_total_current"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 4.0, "total_positive_ion", 3624.0,
        1894.0, "open_circle", "direct_total_current"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 6.7, "total_positive_ion", 3988.0,
        1859.0, "open_circle", "direct_total_current"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 0.67, "Cl+", 3170.0, 2014.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 1.3, "Cl+", 3257.0, 2007.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 2.7, "Cl+", 3447.0, 1979.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 4.0, "Cl+", 3622.0, 1953.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 6.7, "Cl+", 3990.0, 1918.0,
        "filled_square", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 0.67, "Cl2+", 3171.0, 2281.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 1.3, "Cl2+", 3257.0, 2265.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 2.7, "Cl2+", 3449.0, 2261.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 4.0, "Cl2+", 3624.0, 2298.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
    PixelPoint(
        "b", "20% Cl2 + 80% Ar", 6.7, "Cl2+", 3991.0, 2369.0,
        "filled_circle", "mass_resolved_ied_scaled_to_total"),
)

FIELDNAMES = (
    "panel",
    "feed",
    "pressure_Pa",
    "listed_net_power_to_matching_network_W",
    "estimated_plasma_dissipated_power_W",
    "species",
    "ion_flux_mA_cm2",
    "marker_center_x_px",
    "marker_center_y_px",
    "marker",
    "observable_basis",
    "digitization_log10_flux_uncertainty",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _pressure_at_pixel(x_px: float) -> float:
    return 10.0 * (x_px - X_AT_0_PA) / (X_AT_10_PA - X_AT_0_PA)


def _flux_at_pixel(panel: str, y_px: float) -> float:
    if panel == "a":
        y_at_10 = PANEL_A_Y_AT_10
        pixels_per_decade = PANEL_A_Y_AT_1 - PANEL_A_Y_AT_10
    elif panel == "b":
        y_at_10 = PANEL_B_Y_AT_10
        pixels_per_decade = PANEL_B_Y_AT_1 - PANEL_B_Y_AT_10
    else:
        raise ValueError(f"unknown panel {panel!r}")
    log10_flux = 1.0 + (y_at_10 - y_px) / pixels_per_decade
    return 10.0 ** log10_flux


def rows() -> list[dict[str, str]]:
    output = []
    for point in PIXEL_POINTS:
        log_uncertainty = 8.0 / (
            PANEL_A_Y_AT_1 - PANEL_A_Y_AT_10
            if point.panel == "a"
            else PANEL_B_Y_AT_1 - PANEL_B_Y_AT_10
        )
        output.append({
            "panel": point.panel,
            "feed": point.feed,
            "pressure_Pa": f"{point.pressure_Pa:g}",
            "listed_net_power_to_matching_network_W": "300",
            "estimated_plasma_dissipated_power_W": "240",
            "species": point.species,
            "ion_flux_mA_cm2": f"{_flux_at_pixel(point.panel, point.y_px):.5f}",
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "marker": point.marker,
            "observable_basis": point.observable_basis,
            "digitization_log10_flux_uncertainty": f"{log_uncertainty:.5f}",
        })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _series(panel: str, species: str) -> list[float]:
    return [
        float(row["ion_flux_mA_cm2"])
        for row in rows()
        if row["panel"] == panel and row["species"] == species
    ]


def manifest(csv_sha256: str) -> dict[str, object]:
    x_offsets = [
        abs(_pressure_at_pixel(point.x_px) - point.pressure_Pa)
        for point in PIXEL_POINTS
    ]
    total_a = _series("a", "total_positive_ion")
    atomic_a = _series("a", "Cl+")
    molecular_a = _series("a", "Cl2+")
    return {
        "manifest_id": "WANG-OLTHOFF-1999-FIG9-CHLORINE-ION-FLUX-R1",
        "source": {
            "citation": (
                "Y. Wang and J. K. Olthoff, Journal of Applied Physics 85, "
                "6358-6365 (1999), Ion Energy Distributions in Inductively "
                "Coupled Radio-Frequency Discharges in Argon, Nitrogen, "
                "Oxygen, Chlorine, and Their Mixtures"
            ),
            "doi": "10.1063/1.370138",
            "official_nist_url": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": PDF_PAGE,
            "print_page": 6363,
            "figure": "Figure 9",
            "render_dpi": RENDER_DPI,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
        },
        "experiment": {
            "reactor": "inductively coupled GEC rf reference cell",
            "coil": "five-turn planar coil behind quartz window",
            "frequency_MHz": 13.56,
            "flow_sccm": 10.0,
            "window_to_grounded_electrode_gap_cm": 4.13,
            "listed_power_node": "net power to the coil matching network",
            "listed_power_W": 300.0,
            "plasma_dissipated_fraction": "approximately 80%",
            "estimated_plasma_dissipated_power_W": 240.0,
            "measurement": {
                "total": (
                    "total sampled ion current divided by the 10 um orifice "
                    "area"
                ),
                "species": (
                    "mass-resolved IED relative intensities scaled to the "
                    "measured total"
                ),
            },
        },
        "pixel_calibration": {
            "x_at_0_Pa": X_AT_0_PA,
            "x_at_10_Pa": X_AT_10_PA,
            "panel_a_y_at_10_mA_cm2": PANEL_A_Y_AT_10,
            "panel_a_y_at_1_mA_cm2": PANEL_A_Y_AT_1,
            "panel_b_y_at_10_mA_cm2": PANEL_B_Y_AT_10,
            "panel_b_y_at_1_mA_cm2": PANEL_B_Y_AT_1,
            "x_scale": "linear",
            "y_scale": "log10",
        },
        "digitization": {
            "method": (
                "official NIST PDF; 600-dpi Poppler full-page render; "
                "full-resolution visual/PIL marker identity audit; dark-axis "
                "and marker-pixel replay; nearest-neighbor overlap inspection"
            ),
            "vertical_pixel_allowance": 8.0,
            "panel_a_log10_flux_uncertainty": round(
                8.0 / (PANEL_A_Y_AT_1 - PANEL_A_Y_AT_10), 8),
            "panel_b_log10_flux_uncertainty": round(
                8.0 / (PANEL_B_Y_AT_1 - PANEL_B_Y_AT_10), 8),
            "maximum_pressure_setpoint_pixel_offset_Pa": round(
                max(x_offsets), 5),
            "measurement_uncertainty": (
                "not assigned: Figure 9 publishes no chlorine error bars; "
                "the pixel allowance is digitization uncertainty only"
            ),
        },
        "source_internal_checks": {
            "panel_a_species_sum_to_total_fraction": [
                (atomic + molecular) / total
                for total, atomic, molecular
                in zip(total_a, atomic_a, molecular_a)
            ],
            "panel_a_Cl_plus_fraction": [
                atomic / total
                for total, atomic in zip(total_a, atomic_a)
            ],
            "panel_a_total_flux_monotone_in_pressure": all(
                left < right for left, right in zip(total_a, total_a[1:])
            ),
            "panel_a_Cl2_plus_flux_monotone_decreasing": all(
                left > right
                for left, right in zip(molecular_a, molecular_a[1:])
            ),
        },
        "claim_boundary": {
            "valid": [
                "absolute total positive-ion flux versus pressure",
                "reactor-output Cl+ and Cl2+ composition at the grounded plane",
                "held-out validation of a complete reactor/sheath model",
            ],
            "not_valid": [
                "an electron-impact Cl2+ versus Cl+ branching cross section",
                "a Lam proprietary reactor boundary",
                "a capacitively coupled dielectric-etch boundary",
                "a feature-depth calibration",
            ],
            "surface_condition_warning": (
                "the paper attributes its pure-Cl2 species disagreement with "
                "Woodworth et al. to reactor surface conditions"
            ),
        },
        "output": {
            "path": (
                "data/experimental/wang_olthoff_1999/"
                "figure9_chlorine_ion_flux.csv"
            ),
            "sha256": csv_sha256,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def _render_page(pdf_path: Path) -> Image.Image:
    if _sha256(pdf_path) != SOURCE_PDF_SHA256:
        raise RuntimeError("Wang-Olthoff source PDF checksum does not match")
    with tempfile.TemporaryDirectory(
            prefix="petch-wang-olthoff-vision-") as directory:
        prefix = Path(directory) / "page6_600dpi"
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
                "Wang-Olthoff 600-dpi page render checksum does not match")
        image = Image.open(render_path).convert("RGB")
        image.load()
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected render size {image.size}; expected {RENDER_SIZE}")
    return image


def _vision_audit(image: Image.Image) -> None:
    gray = np.asarray(image.convert("L"))
    dark = gray < 110
    checks = {
        "panel a left axis": np.mean(
            np.min(gray[470:1370, 3070:3079], axis=1) < 150),
        "panel a bottom axis": np.mean(
            np.min(gray[1360:1372, 3070:4440], axis=0) < 150),
        "panel b left axis": np.mean(
            np.min(gray[1615:2520, 3070:3080], axis=1) < 150),
        "panel b bottom axis": np.mean(
            np.min(gray[2508:2525, 3070:4440], axis=0) < 150),
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
        minimum = 130 if point.marker == "open_circle" else 450
        if count < minimum:
            raise RuntimeError(
                f"failed to recover {point.panel} {point.species} marker "
                f"at {point.pressure_Pa:g} Pa: {count} dark pixels")


def write_committed_files() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    payload = csv_text()
    CSV_PATH.write_text(payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_text(payload), encoding="utf-8")


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Wang-Olthoff CSV is not reproduced")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Wang-Olthoff manifest is not reproduced")


def draw_overlay(image: Image.Image, path: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {
        "total_positive_ion": "#e41a1c",
        "Cl+": "#377eb8",
        "Cl2+": "#4daf4a",
    }
    for point in PIXEL_POINTS:
        color = colors[point.species]
        radius = 22
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline=color,
            width=4,
        )
        draw.line(
            (
                point.x_px - 14, point.y_px,
                point.x_px + 14, point.y_px,
            ),
            fill=color,
            width=3,
        )
        draw.line(
            (
                point.x_px, point.y_px - 14,
                point.x_px, point.y_px + 14,
            ),
            fill=color,
            width=3,
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--overlay", type=Path)
    arguments = parser.parse_args()

    image = None
    if arguments.pdf is not None:
        image = _render_page(arguments.pdf)
        _vision_audit(image)
    if arguments.write:
        write_committed_files()
    verify_committed_files()
    if arguments.overlay is not None:
        if image is None:
            raise ValueError("--overlay requires --pdf")
        draw_overlay(image, arguments.overlay)


if __name__ == "__main__":
    main()
