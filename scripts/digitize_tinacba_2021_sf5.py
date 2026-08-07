#!/usr/bin/env python3
"""Replay Tinacba et al. 2021 Figure 8 SF5+ beam-yield digitization.

The source compares independently calculated molecular-dynamics yields with
mass/energy-selected SF5+ beam measurements on Si and SiO2.  Energies are the
plotted setpoints; marker x pixels are retained only as a placement check.
The experimental yield was obtained by the authors from profilometer depth,
Faraday-cup ion dose, and material number density.
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
OUTPUT_DIR = ROOT / "data" / "experimental" / "tinacba_2021"
CSV_PATH = OUTPUT_DIR / "figure8_sf5_md_experiment.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"
DEFAULT_RENDER = ROOT / "tmp" / "sources" / "tinacba_2021" / "page8_300dpi.png"

SOURCE_URL = (
    "https://researchmap.jp/satoshi_hamaguchi/published_papers/"
    "37820817/attachment_file.pdf"
)
SOURCE_PDF_SHA256 = (
    "c0be3b475aa17b396c1f788baee14ba37b9026b264bb870dc2553055f27b31ad"
)
RENDER_SHA256 = (
    "9a0867caa7991dbcf8495baf9ed5a47a4d140fe8055b578777664edcf0800a02"
)
RENDER_SIZE = (2550, 3375)

Y_AT_ZERO = 2484.0
Y_AT_SEVEN = 1808.5
YIELD_DIGITIZATION_BOUND = 0.04
LEFT_X_AT_ZERO = 306.5
LEFT_X_AT_2200 = 1233.5
RIGHT_X_AT_ZERO = 1394.5
RIGHT_X_AT_2200 = 2320.0


@dataclass(frozen=True)
class PixelPoint:
    material: str
    series: str
    energy_eV: int
    x_px: float
    y_px: float
    marker: str


POINTS = (
    PixelPoint("Si", "sf5_md", 150, 370.0, 2289.0, "open_diamond"),
    PixelPoint("Si", "sf5_md", 300, 433.0, 2257.0, "open_diamond"),
    PixelPoint("Si", "sf5_md", 500, 517.0, 2209.0, "open_diamond"),
    PixelPoint("Si", "sf5_md", 1000, 728.0, 2074.0, "open_diamond"),
    PixelPoint("Si", "sf5_md", 1500, 939.0, 1922.0, "open_diamond"),
    PixelPoint("Si", "sf5_md", 2000, 1149.0, 1864.0, "open_diamond"),
    PixelPoint("Si", "sf5_experiment", 150, 370.5, 2314.5, "filled_circle"),
    PixelPoint("Si", "sf5_experiment", 2000, 1149.5, 1888.5, "filled_circle"),
    PixelPoint("SiO2", "sf5_md", 150, 1457.5, 2405.0, "open_diamond"),
    PixelPoint("SiO2", "sf5_md", 300, 1521.0, 2353.0, "open_diamond"),
    PixelPoint("SiO2", "sf5_md", 500, 1605.0, 2320.0, "open_diamond"),
    PixelPoint("SiO2", "sf5_md", 1000, 1815.0, 2277.0, "open_diamond"),
    PixelPoint("SiO2", "sf5_md", 1500, 2026.0, 2225.0, "open_diamond"),
    PixelPoint("SiO2", "sf5_md", 2000, 2236.0, 2161.0, "open_diamond"),
    # At 150 eV the filled experimental circle and open MD diamond overlap.
    PixelPoint("SiO2", "sf5_experiment", 150, 1457.5, 2405.0, "filled_circle"),
    PixelPoint("SiO2", "sf5_experiment", 2000, 2235.5, 2174.5, "filled_circle"),
)

FIELDNAMES = (
    "material",
    "series",
    "energy_eV",
    "si_removal_yield_per_sf5_ion",
    "marker_center_x_px",
    "marker_center_y_px",
    "marker",
    "digitization_yield_bound",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _yield_at_pixel(y_px: float) -> float:
    return 7.0 * (Y_AT_ZERO - y_px) / (Y_AT_ZERO - Y_AT_SEVEN)


def _energy_at_pixel(material: str, x_px: float) -> float:
    if material == "Si":
        left, right = LEFT_X_AT_ZERO, LEFT_X_AT_2200
    else:
        left, right = RIGHT_X_AT_ZERO, RIGHT_X_AT_2200
    return 2200.0 * (x_px - left) / (right - left)


def rows() -> list[dict[str, str]]:
    out = []
    for point in POINTS:
        out.append(
            {
                "material": point.material,
                "series": point.series,
                "energy_eV": str(point.energy_eV),
                "si_removal_yield_per_sf5_ion": (
                    f"{_yield_at_pixel(point.y_px):.4f}"
                ),
                "marker_center_x_px": f"{point.x_px:.1f}",
                "marker_center_y_px": f"{point.y_px:.1f}",
                "marker": point.marker,
                "digitization_yield_bound": (
                    f"{YIELD_DIGITIZATION_BOUND:.2f}"
                ),
            }
        )
    return out


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    x_errors = [
        abs(_energy_at_pixel(point.material, point.x_px) - point.energy_eV)
        for point in POINTS
    ]
    return {
        "manifest_id": "TINACBA-2021-FIG8-SF5-MD-EXPERIMENT-R1",
        "source": {
            "citation": (
                "E. J. C. Tinacba, T. Ito, K. Karahashi, M. Isobe, and "
                "S. Hamaguchi, JVST B 39, 043203 (2021)"
            ),
            "doi": "10.1116/6.0001230",
            "retrieval_url": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 8,
            "print_page": "043203-7",
            "figure": "Figure 8",
            "render_dpi": 300,
            "render_sha256": RENDER_SHA256,
            "render_size_px": list(RENDER_SIZE),
        },
        "experiment": {
            "projectile": "mass-selected SF5+",
            "incidence": "normal",
            "energy_measurement": "energy-mass analyzer",
            "dose_measurement": "Faraday cup at sample position",
            "depth_measurement": "Dektak3ST contact profilometer",
            "radical_beam": "none",
            "reported_energy_range_eV": [150, 2000],
            "yield_equation": "Y = d*N/D",
        },
        "atomistic_provider": {
            "method": "modified Stillinger-Weber molecular dynamics",
            "surface_temperature_K": 300,
            "impacts_per_trajectory": 4000,
            "sulfur_model": (
                "DFT-informed S-F carrier model; S-S, S-Si, and S-O chemical "
                "reactions intentionally suppressed"
            ),
            "fit_to_beam_depth_or_yield": False,
        },
        "pixel_calibration": {
            "y_at_zero_yield": Y_AT_ZERO,
            "y_at_seven_yield": Y_AT_SEVEN,
            "left_x_at_zero_eV": LEFT_X_AT_ZERO,
            "left_x_at_2200_eV": LEFT_X_AT_2200,
            "right_x_at_zero_eV": RIGHT_X_AT_ZERO,
            "right_x_at_2200_eV": RIGHT_X_AT_2200,
        },
        "digitization": {
            "method": (
                "300-dpi Poppler render; Pillow dark-axis verification; "
                "full-resolution marker-center transcription and QA overlay"
            ),
            "point_count": len(POINTS),
            "yield_bound": YIELD_DIGITIZATION_BOUND,
            "yield_bound_basis": (
                "3.9 source pixels on the linear 0-7 axis; covers line width, "
                "anti-aliasing, and center placement"
            ),
            "maximum_setpoint_x_offset_eV": round(max(x_errors), 3),
            "energy_policy": (
                "use plotted nominal setpoints; retain x pixels as a check"
            ),
        },
        "claim_boundary": {
            "supports": (
                "retrospective independent validation of the source's "
                "DFT-informed atomistic provider at mass-selected beam points"
            ),
            "does_not_support": [
                "an SF6 reactor boundary",
                "an SFx+ mixture law",
                "off-normal incidence",
                "neutral-F synergy",
                "sulfur-surface chemistry",
                "a Krueger C4F6 depth fit",
            ],
        },
        "output": {
            "path": str(CSV_PATH.relative_to(ROOT)),
            "sha256": csv_sha256,
        },
    }


def manifest_text(payload: str) -> str:
    return json.dumps(
        manifest(hashlib.sha256(payload.encode()).hexdigest()),
        indent=2,
    ) + "\n"


def verify_render(path: Path) -> Image.Image:
    if _sha256(path) != RENDER_SHA256:
        raise RuntimeError("Figure 8 render checksum mismatch")
    image = Image.open(path).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(f"unexpected render size: {image.size}")
    gray = np.asarray(image.convert("L"))
    checks = {
        "top": np.mean(gray[1808, 306:2321] < 96),
        "bottom": np.mean(gray[2484, 306:2321] < 96),
        "left_si": np.mean(gray[1808:2485, 307] < 96),
        "right_sio2": np.mean(gray[1808:2485, 2320] < 96),
    }
    if checks["top"] < 0.8 or checks["bottom"] < 0.8:
        raise RuntimeError(f"horizontal axis verification failed: {checks}")
    if checks["left_si"] < 0.8 or checks["right_sio2"] < 0.8:
        raise RuntimeError(f"vertical axis verification failed: {checks}")
    return image


def verify_committed() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Tinacba CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed Tinacba manifest is stale")


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {"sf5_md": "#e41a1c", "sf5_experiment": "#377eb8"}
    for point in POINTS:
        color = colors[point.series]
        radius = 14
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline=color,
            width=3,
        )
        draw.line(
            (
                point.x_px - radius - 5,
                point.y_px,
                point.x_px + radius + 5,
                point.y_px,
            ),
            fill=color,
            width=2,
        )
        draw.line(
            (
                point.x_px,
                point.y_px - radius - 5,
                point.x_px,
                point.y_px + radius + 5,
            ),
            fill=color,
            width=2,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.crop((270, 1760, 2350, 2525)).save(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render", type=Path, default=DEFAULT_RENDER)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()
    image = verify_render(args.render)
    verify_committed()
    if args.overlay:
        draw_overlay(image, args.overlay)
    print(
        json.dumps(
            {
                "status": "verified",
                "points": len(POINTS),
                "maximum_setpoint_x_offset_eV": manifest("")["digitization"][
                    "maximum_setpoint_x_offset_eV"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
