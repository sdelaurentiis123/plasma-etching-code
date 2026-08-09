#!/usr/bin/env python3
"""Digitize Wise et al. (1996) Figure 3 direct GEC-ICP markers.

Figure 3 overlays the source simulation curves with Langmuir-probe marker
measurements of radial electron density, electron temperature, and plasma
potential at 180 W and 20 mTorr chlorine.  Only the visible square markers are
retained here.  Pixel centers were isolated on a checksum-pinned 300-dpi full
page render and replayed with PIL; no curve pixels are represented as data.
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

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "wise_1996_gec_icp"
CSV_PATH = OUTPUT_DIRECTORY / "figure3_radial_measurements.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "figure3_digitization_manifest.json"
SOURCE_URL = (
    "https://www.chee.uh.edu/sites/chbe/files/faculty/economou/"
    "apl_mpres_1996.pdf"
)
SOURCE_PDF_SHA256 = (
    "25bfedd2ebf45bcb62b237ea0a5e27d1e495850673d73f3ded4c8b8321999e65"
)
RENDER_PAGE = 3
RENDER_DPI = 300
RENDER_SIZE = (2550, 3300)
RENDER_SHA256 = (
    "067184f05194b8667387f67de51ca89821838d2d82539f58d48b1e7ed1927d6a"
)
RADIAL_DISTANCE_M = (0.0, 0.0127, 0.0254, 0.0381, 0.0508, 0.0635, 0.0762)


@dataclass(frozen=True)
class Marker:
    observable: str
    units: str
    scale: float
    x_px: float
    y_px: float
    marker_state: str = "isolated_square"


MARKERS = (
    # Top panel: plot y=181 -> 16 and y=509 -> 0 in printed 1e16 m^-3.
    *(Marker("electron_density", "m^-3", 1.0e16, x, y, state)
      for x, y, state in zip(
          (563., 605., 647., 688., 731., 772., 814.),
          (220., 247.5, 322., 366., 398.5, 443., 468.),
          ("axis_overlapping_square",) + ("isolated_square",) * 6)),
    # Middle panel: y=687 -> 4 eV and y=1014 -> 0 eV.
    *(Marker("electron_temperature", "eV", 1.0, x, y, state)
      for x, y, state in zip(
          (563., 600.5, 642.5, 684.5, 726.5, 768., 809.5),
          (701., 699., 696.5, 708., 731.5, 750.5, 770.5),
          ("axis_overlapping_square",) + ("isolated_square",) * 6)),
    # Bottom panel: y=1148 -> 20 V and y=1477 -> 0 V.
    *(Marker("plasma_potential", "V", 1.0, x, y, state)
      for x, y, state in zip(
          (559., 598., 641.5, 681.5, 723.5, 765.5, 807.5),
          (1199., 1206., 1222., 1233., 1258., 1281., 1304.),
          ("axis_overlapping_square",) + ("isolated_square",) * 6)),
)

PANEL_CALIBRATION = {
    "electron_density": {
        "y_top_px": 181.0, "y_bottom_px": 509.0,
        "top_value": 16.0, "bottom_value": 0.0,
    },
    "electron_temperature": {
        "y_top_px": 687.0, "y_bottom_px": 1014.0,
        "top_value": 4.0, "bottom_value": 0.0,
    },
    "plasma_potential": {
        "y_top_px": 1148.0, "y_bottom_px": 1477.0,
        "top_value": 20.0, "bottom_value": 0.0,
    },
}

FIELDNAMES = (
    "observable", "radial_distance_m", "value", "units",
    "marker_center_x_px", "marker_center_y_px", "marker_state",
    "digitization_vertical_uncertainty_px", "source_uncertainty",
    "evidence_type",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _marker_value(marker: Marker) -> float:
    calibration = PANEL_CALIBRATION[marker.observable]
    fraction = (
        (marker.y_px - calibration["y_top_px"])
        / (calibration["y_bottom_px"] - calibration["y_top_px"])
    )
    plotted = (
        calibration["top_value"]
        + fraction
        * (calibration["bottom_value"] - calibration["top_value"])
    )
    return marker.scale * plotted


def rows():
    output = []
    for panel_index, observable in enumerate(PANEL_CALIBRATION):
        panel_markers = [
            marker for marker in MARKERS if marker.observable == observable]
        if len(panel_markers) != len(RADIAL_DISTANCE_M):
            raise RuntimeError("Wise Figure 3 marker count mismatch")
        for radius, marker in zip(RADIAL_DISTANCE_M, panel_markers):
            output.append({
                "observable": observable,
                "radial_distance_m": f"{radius:.4f}",
                "value": f"{_marker_value(marker):.9g}",
                "units": marker.units,
                "marker_center_x_px": f"{marker.x_px:.1f}",
                "marker_center_y_px": f"{marker.y_px:.1f}",
                "marker_state": marker.marker_state,
                "digitization_vertical_uncertainty_px": "3",
                "source_uncertainty": "not_printed",
                "evidence_type": "direct_reactor_measurement_digitized",
            })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str):
    by_observable = {
        observable: [
            _marker_value(marker) for marker in MARKERS
            if marker.observable == observable
        ]
        for observable in PANEL_CALIBRATION
    }
    return {
        "manifest_id": "WISE-1996-FIG3-GEC-ICP-RADIAL-R1",
        "source": {
            "citation": (
                "R. S. Wise, D. P. Lymberopoulos, and D. J. Economou, "
                "Rapid two-dimensional self-consistent simulation of "
                "inductively coupled plasma and comparison with experimental "
                "data, Applied Physics Letters 68, 2499-2501 (1996)"
            ),
            "official_url": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 3,
            "print_page": 2501,
            "figure": "3",
            "render_dpi": RENDER_DPI,
            "render_size_px": list(RENDER_SIZE),
            "render_sha256": RENDER_SHA256,
        },
        "experiment_context": {
            "reactor": "GEC reference-cell chlorine ICP",
            "plasma_power_W": 180.0,
            "pressure_mTorr": 20.0,
            "flow_sccm": 20.0,
            "frequency_MHz": 13.56,
            "wafer_present": False,
            "substrate_bias": False,
            "diagnostic": "radially translated Langmuir probe",
            "marker_source_reference": 15,
        },
        "pixel_calibration_300dpi": {
            "x_positions": {
                "radial_distance_m": list(RADIAL_DISTANCE_M),
                "note": (
                    "visible marker train uses 0.5-inch (0.0127 m) increments; "
                    "axis spans approximately 0--0.14 m"
                ),
            },
            **PANEL_CALIBRATION,
        },
        "digitization": {
            "marker_count": len(MARKERS),
            "method": (
                "300-dpi checksum-pinned Poppler full-page render; "
                "full-resolution visual inspection; PIL grayscale component "
                "isolation and 11x11 filled-square density scan; printed-axis "
                "replay; PIL QA overlay"
            ),
            "vertical_pixel_uncertainty": 3.0,
            "axis_overlapping_marker_count": sum(
                marker.marker_state == "axis_overlapping_square"
                for marker in MARKERS),
            "source_measurement_uncertainty": (
                "not printed in Figure 3; the paper separately warns that "
                "Langmuir-probe electron-temperature errors can be substantial"
            ),
        },
        "source_internal_checks": {
            "electron_density_decreases_after_axis": all(
                left > right for left, right in zip(
                    by_observable["electron_density"],
                    by_observable["electron_density"][1:],
                )
            ),
            "plasma_potential_decreases_radially": all(
                left > right for left, right in zip(
                    by_observable["plasma_potential"],
                    by_observable["plasma_potential"][1:],
                )
            ),
            "electron_temperature_peak_eV": max(
                by_observable["electron_temperature"]),
        },
        "claim_boundary": {
            "valid": [
                "independent radial reactor-state grade at the printed GEC condition",
                "shape/width grade for a chlorine charged spatial model",
                "electron-energy and Boltzmann-potential localization",
            ],
            "not_valid": [
                "Lam reactor calibration",
                "wafer flux or IEAD validation",
                "feature-depth calibration or validation",
                "absolute uncertainty-weighted chi-square because error bars are absent",
            ],
        },
        "output": {
            "path": str(CSV_PATH.relative_to(ROOT)),
            "sha256": csv_sha256,
        },
    }


def _render(pdf: Path, directory: Path) -> Path:
    prefix = directory / "wise_1996_page3"
    subprocess.run([
        "pdftoppm", "-f", str(RENDER_PAGE), "-l", str(RENDER_PAGE),
        "-png", "-r", str(RENDER_DPI), "-singlefile", str(pdf), str(prefix),
    ], check=True)
    image = prefix.with_suffix(".png")
    if Image.open(image).size != RENDER_SIZE or _sha256(image) != RENDER_SHA256:
        raise ValueError("Wise render checksum/size mismatch")
    return image


def write(*, pdf: Path | None = None, overlay: Path | None = None):
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    text = csv_text()
    CSV_PATH.write_text(text, encoding="utf-8")
    digest = hashlib.sha256(text.encode()).hexdigest()
    MANIFEST_PATH.write_text(
        json.dumps(manifest(digest), indent=2) + "\n", encoding="utf-8")
    if pdf is not None:
        if _sha256(pdf) != SOURCE_PDF_SHA256:
            raise ValueError("Wise source PDF checksum mismatch")
        with tempfile.TemporaryDirectory() as temporary:
            image_path = _render(pdf, Path(temporary))
            image = Image.open(image_path).convert("RGB")
            draw = ImageDraw.Draw(image)
            colors = {
                "electron_density": "red",
                "electron_temperature": "lime",
                "plasma_potential": "deepskyblue",
            }
            for marker in MARKERS:
                radius = 10
                draw.ellipse((
                    marker.x_px - radius, marker.y_px - radius,
                    marker.x_px + radius, marker.y_px + radius,
                ), outline=colors[marker.observable], width=3)
            target = (
                overlay if overlay is not None
                else Path("/private/tmp/wise_1996_figure3_qa.png"))
            target.parent.mkdir(parents=True, exist_ok=True)
            image.save(target)
            print(target)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path)
    parser.add_argument("--overlay", type=Path)
    arguments = parser.parse_args()
    write(pdf=arguments.pdf, overlay=arguments.overlay)


if __name__ == "__main__":
    main()
