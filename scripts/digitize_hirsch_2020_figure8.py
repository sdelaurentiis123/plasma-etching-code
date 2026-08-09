#!/usr/bin/env python3
"""Checksum-backed PIL extraction of Hirsch 2020 Figure-8 PAE curve.

The red curve is an author-assumed anti-synergy response constrained by the
paper's limiting arguments; it is not direct photon-flux data.  This replay
keeps that evidence class explicit.  Source pixels are not redistributed.
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
DEFAULT_OUTPUT = ROOT / "data" / "experimental" / "hirsch_2020_pae_iae"
SOURCE_PDF_SHA256 = (
    "36a5644eb3314e5baa0e999c3023e0f7ec547cf77527c055ac2386c649ef9b62"
)
RENDER_SHA256 = (
    "7ab90ee667b6715b73736644973a22c9799567e4398bcb1f554b804504a95020"
)
PDF_PAGE = 9
RENDER_DPI = 220
RENDER_SIZE = (1870, 2475)

# Axis centers established on the original render.  The top and bottom black
# strokes occupy y=1362--1363 and y=1732--1733; x=302 is the vertical axis and
# x=745 the right boundary.
X_ZERO_PX = 302.0
X_HUNDRED_PX = 745.0
Y_ZERO_PX = 1733.0
Y_ONE_PX = 1363.0
SAMPLED_DUTY_PERCENT = tuple(range(5, 91, 5))
X_HALF_WINDOW_PX = 3
RED_MINIMUM = 150
RED_DOMINANCE = 1.5


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _render(pdf: Path, directory: Path) -> Path:
    prefix = directory / "hirsch_2020_page009"
    subprocess.run(
        [
            "pdftoppm", "-f", str(PDF_PAGE), "-l", str(PDF_PAGE),
            "-r", str(RENDER_DPI), "-png", "-singlefile", str(pdf),
            str(prefix),
        ],
        check=True,
    )
    return prefix.with_suffix(".png")


def _vertical_clusters(values: np.ndarray) -> list[np.ndarray]:
    unique = np.unique(values)
    groups = []
    start = 0
    for index in range(1, unique.size):
        if unique[index] > unique[index - 1] + 2:
            groups.append(unique[start:index])
            start = index
    groups.append(unique[start:])
    return groups


def _extract(render: Path) -> list[dict[str, float | int | str]]:
    if _hash(render) != RENDER_SHA256:
        raise ValueError("220-dpi Hirsch Figure-8 render checksum changed")
    image = Image.open(render).convert("RGB")
    if image.size != RENDER_SIZE:
        raise ValueError(f"unexpected rendered page size: {image.size}")
    pixels = np.asarray(image)
    red = (
        (pixels[:, :, 0] > RED_MINIMUM)
        & (pixels[:, :, 0] > RED_DOMINANCE * pixels[:, :, 1])
        & (pixels[:, :, 0] > RED_DOMINANCE * pixels[:, :, 2])
    )
    x_scale = (X_HUNDRED_PX - X_ZERO_PX) / 100.0
    y_scale = Y_ZERO_PX - Y_ONE_PX
    extracted: list[dict[str, float | int | str]] = [
        {
            "dc_bias_duty_cycle_percent": 0.0,
            "relative_pae_yield": 1.0,
            "center_x_px": X_ZERO_PX,
            "center_y_px": Y_ONE_PX,
            "line_y_min_px": Y_ONE_PX,
            "line_y_max_px": Y_ONE_PX,
            "evidence": "textual limiting constraint",
        }
    ]
    for duty in SAMPLED_DUTY_PERCENT:
        target_x = X_ZERO_PX + x_scale * duty
        center_x = int(round(target_x))
        x0 = center_x - X_HALF_WINDOW_PX
        x1 = center_x + X_HALF_WINDOW_PX + 1
        local_y, local_x = np.where(red[int(Y_ONE_PX):int(Y_ZERO_PX) + 1, x0:x1])
        if local_y.size == 0:
            raise ValueError(f"no red curve pixels recovered at {duty}% duty")
        source_y = local_y + int(Y_ONE_PX)
        # The red annotation lies below the curve near 5--10%.  The curve is
        # the uppermost connected vertical cluster at every sampled abscissa.
        cluster = min(_vertical_clusters(source_y), key=lambda item: np.mean(item))
        mask = (source_y >= cluster.min()) & (source_y <= cluster.max())
        center_y = float(np.median(source_y[mask]))
        relative_yield = float(np.clip(
            (Y_ZERO_PX - center_y) / y_scale, 0.0, 1.0))
        extracted.append({
            "dc_bias_duty_cycle_percent": float(duty),
            "relative_pae_yield": relative_yield,
            "center_x_px": float(target_x),
            "center_y_px": center_y,
            "line_y_min_px": int(cluster.min()),
            "line_y_max_px": int(cluster.max()),
            "evidence": "PIL digitization of author-assumed red curve",
        })
    return extracted


def _write(output: Path, extracted) -> tuple[Path, Path]:
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "figure8_relative_pae_yield.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "dc_bias_duty_cycle_percent", "relative_pae_yield",
            "center_x_px", "center_y_px", "line_y_min_px",
            "line_y_max_px", "evidence",
        ))
        writer.writeheader()
        writer.writerows(extracted)
    manifest = {
        "schema": "petch.hirsch-2020-figure8-pae-anti-synergy.v1",
        "source": {
            "citation": (
                "E. W. Hirsch, W. H. Du, and D. J. Economou, Evidence for "
                "anti-synergism between ion-assisted etching and in-plasma "
                "photoassisted etching of silicon in a high-density chlorine "
                "plasma, JVST A 38, 023009 (2020), DOI 10.1116/1.5138189"
            ),
            "official_url": (
                "https://www.chee.uh.edu/sites/chbe/files/faculty/economou/"
                "pae-iae-anti-synergism.pdf"
            ),
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": PDF_PAGE,
            "journal_page": "023009-9",
            "render_dpi": RENDER_DPI,
            "render_size_px": list(RENDER_SIZE),
            "render_sha256": RENDER_SHA256,
        },
        "experiment": {
            "chemistry": "10% Cl2 / 90% Ar",
            "pressure_mTorr": 60.0,
            "flow_case": "high flow, 250 sccm",
            "ion_flux_cm2_s": 1.4e17,
            "bias": "-70 V pulsed DC at 10 kHz",
            "surface": "silicon",
        },
        "pixel_calibration": {
            "duty_zero_x_px": X_ZERO_PX,
            "duty_hundred_x_px": X_HUNDRED_PX,
            "relative_yield_zero_y_px": Y_ZERO_PX,
            "relative_yield_one_y_px": Y_ONE_PX,
            "red_minimum": RED_MINIMUM,
            "red_dominance": RED_DOMINANCE,
            "x_half_window_px": X_HALF_WINDOW_PX,
        },
        "extraction": {
            "method": (
                "original-resolution PIL/NumPy red-pixel isolation; uppermost "
                "vertical cluster rejects the nearby red annotation"
            ),
            "points": extracted,
            "feature_depth_used": False,
            "reactor_parameter_fit_used": False,
        },
        "evidence_class": {
            "curve": (
                "author-assumed smooth relative PAE-yield response constrained "
                "to unity at high Cl/ion ratio and nearly zero at low ratio"
            ),
            "not": (
                "a direct spectral photon-flux measurement or a universal RF "
                "suppression law"
            ),
        },
        "scope": {
            "valid_use": (
                "pulsed-DC high-flow sensitivity for the existence and scale "
                "of PAE/IAE anti-synergy"
            ),
            "invalid_use": (
                "direct transfer to 13.56 MHz RF, pure Cl2, another reactor, "
                "or absolute depth without a spectral photon boundary"
            ),
        },
        "output": {
            "csv_path": (
                "data/experimental/hirsch_2020_pae_iae/"
                "figure8_relative_pae_yield.csv"
            ),
            "csv_sha256": _hash(csv_path),
        },
    }
    manifest_path = output / "digitization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return csv_path, manifest_path


def _overlay(render: Path, output: Path, extracted) -> None:
    image = Image.open(render).convert("RGB")
    draw = ImageDraw.Draw(image)
    draw.line(
        [(X_ZERO_PX, Y_ZERO_PX), (X_HUNDRED_PX, Y_ZERO_PX)],
        fill=(0, 180, 0), width=3,
    )
    draw.line(
        [(X_ZERO_PX, Y_ZERO_PX), (X_ZERO_PX, Y_ONE_PX)],
        fill=(0, 180, 0), width=3,
    )
    for row in extracted:
        x = row["center_x_px"]
        y = row["center_y_px"]
        draw.ellipse((x - 5, y - 5, x + 5, y + 5), outline=(0, 0, 255), width=2)
    output.parent.mkdir(parents=True, exist_ok=True)
    image.crop((220, 1300, 790, 1790)).save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overlay", type=Path)
    arguments = parser.parse_args()
    if _hash(arguments.pdf) != SOURCE_PDF_SHA256:
        raise ValueError("Hirsch 2020 source PDF checksum changed")
    with tempfile.TemporaryDirectory(prefix="petch-hirsch-fig8-") as directory:
        render = _render(arguments.pdf, Path(directory))
        extracted = _extract(render)
        csv_path, manifest_path = _write(arguments.output, extracted)
        if arguments.overlay is not None:
            _overlay(render, arguments.overlay, extracted)
    print(csv_path)
    print(manifest_path)


if __name__ == "__main__":
    main()
