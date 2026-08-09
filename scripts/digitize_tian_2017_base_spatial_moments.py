#!/usr/bin/env python3
"""Reproduce coarse spatial moments from Tian 2017 Figures 5.2--5.3.

This is a deliberately low-rank, source-model digitization for deterministic
radiation model discovery.  It retains only contour-supported hot-core,
warm-bulk, and cold-shell moments; it is not an experimental field and is not
presented as a pixel-perfect reconstruction of the HPEM mesh.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
    / "base_case_spatial_moments.json"
)
SOURCE_PDF_SHA256 = (
    "4d260ab9e85240bd051ccb3ba32cd047b4ac1ddb6c309dbbca4093822e37790b")
SOURCE_URL = "https://cpseg.eecs.umich.edu/pub/theses/tian_peng_phd_thesis.pdf"
RENDER_DPI = 300
RENDER_SIZE = [2550, 3300]
RENDER_SHA256 = {
    163: "2b20db4ab6176869c5fa866b54e30c9daad4bbe8cf4f5542e498fdf5fc225e12",
    164: "7513344258a23940fbc180e0ca7df2a83e57496e787374ecb8b4c78d11a24a2d",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def payload() -> dict[str, object]:
    return {
        "manifest_id": "TIAN-2017-FIG5.2-FIG5.3-SPATIAL-MOMENTS-R1",
        "source": {
            "citation": (
                "P. Tian, Controlling Photon and Ion Fluxes in Low Pressure "
                "Low Temperature Plasmas, PhD dissertation, University of "
                "Michigan (2017; published 2018)"),
            "official_url": SOURCE_URL,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_pages": [163, 164],
            "print_pages": [145, 146],
            "figures": ["5.2(b)", "5.3(b), middle row"],
            "render_dpi": RENDER_DPI,
            "render_size_px": RENDER_SIZE,
            "render_sha256_by_pdf_page": {
                str(page): digest for page, digest in RENDER_SHA256.items()
            },
        },
        "condition": {
            "Ar_fraction": 0.8,
            "Cl2_fraction": 0.2,
            "pressure_mTorr": 20.0,
            "flow_sccm": 200.0,
            "frequency_MHz": 10.0,
            "continuous_wave_power_W": 150.0,
            "reactor_radius_cm": 11.25,
            "plasma_height_cm": 12.0,
        },
        "pixel_calibration_300dpi": {
            "figure_5_2_b_gas_temperature": {
                "plot_left_px": 1318,
                "plot_right_px": 1705,
                "plot_top_12cm_px": 563,
                "plot_bottom_0cm_px": 1001,
            },
            "figure_5_3_b_Ar_1s4_density": {
                "plot_left_px": 1320,
                "plot_right_px": 1711,
                "plot_top_12cm_px": 1139,
                "plot_bottom_0cm_px": 1555,
            },
            "visual_QA": (
                "full-resolution pages inspected; zone edges replayed as "
                "PIL overlays against both panels"),
        },
        "axisymmetric_zone_field": {
            "radial_edges_cm": [0.0, 7.5, 10.0, 11.25],
            "axial_edges_cm": [0.0, 2.5, 6.0, 10.5, 12.0],
            "cell_zone_index_radial_by_axial": [
                [2, 1, 0, 1],
                [2, 1, 1, 1],
                [2, 2, 2, 2],
            ],
            "zones": [
                {
                    "zone_index": 0,
                    "name": "hot_emitting_core",
                    "gas_temperature_K": 800.0,
                    "Ar_1s4_emitter_density_cm3": 4.7e9,
                    "visible_contour_support": (
                        "Tgas labels 795 and 830 K; Ar(1s4) maximum label "
                        "4.7e9 cm^-3"),
                },
                {
                    "zone_index": 1,
                    "name": "warm_bulk",
                    "gas_temperature_K": 650.0,
                    "Ar_1s4_emitter_density_cm3": 1.2e9,
                    "visible_contour_support": (
                        "Tgas bounded by printed 590 and 720 K contours; "
                        "Ar(1s4) bounded by 0.8 and 1.6e9 cm^-3 labels"),
                },
                {
                    "zone_index": 2,
                    "name": "cold_boundary_shell",
                    "gas_temperature_K": 400.0,
                    "Ar_1s4_emitter_density_cm3": 2.0e8,
                    "visible_contour_support": (
                        "outside printed 455 K contour; Ar(1s4) outer label "
                        "0.2e9 cm^-3"),
                },
            ],
        },
        "digitization": {
            "evidence_type": "source_equipment_model_digitized",
            "not_measurement": True,
            "method": (
                "300-dpi checksum-pinned Poppler render; full-resolution "
                "visual inspection; printed contour anchors; three-zone "
                "axisymmetric moment compression; PIL boundary replay"),
            "model_reduction": True,
            "fitted_to_trapping_factor": False,
            "uncertainty": (
                "zone averages and boundaries are bracketed by visible "
                "contours; no source mesh values or uncertainty were printed"),
        },
        "claim_boundary": {
            "valid": [
                "base-case spatial-moment sensitivity",
                "deterministic zonal radiation model-discovery input",
                "cold-absorber/hot-emitter mechanism check",
            ],
            "not_valid": [
                "experimental field validation",
                "pixel-perfect reconstruction of the HPEM mesh",
                "Ar(1s2) emitter field independent of Ar(1s4)",
                "mixture-sweep spatial closure",
                "feature-depth calibration",
            ],
        },
    }


def text() -> str:
    return json.dumps(payload(), indent=2) + "\n"


def _render(pdf: Path) -> dict[int, Image.Image]:
    if _sha256(pdf) != SOURCE_PDF_SHA256:
        raise RuntimeError("Tian source PDF checksum does not match")
    images = {}
    with tempfile.TemporaryDirectory(prefix="petch-tian-spatial-") as folder:
        for page, expected_hash in RENDER_SHA256.items():
            prefix = Path(folder) / f"page{page}"
            subprocess.run([
                "pdftoppm", "-f", str(page), "-l", str(page),
                "-r", str(RENDER_DPI), "-png", "-singlefile",
                str(pdf), str(prefix),
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            image_path = prefix.with_suffix(".png")
            if _sha256(image_path) != expected_hash:
                raise RuntimeError(f"Tian page {page} render checksum mismatch")
            images[page] = Image.open(image_path).convert("RGB").copy()
    return images


def _coordinate(value_cm: float, low_px: int, high_px: int, maximum_cm: float) -> int:
    return round(low_px + value_cm / maximum_cm * (high_px - low_px))


def qa_overlays(pdf: Path, output_directory: Path) -> None:
    images = _render(pdf)
    output_directory.mkdir(parents=True, exist_ok=True)
    panels = {
        163: (1318, 1705, 1001, 563),
        164: (1320, 1711, 1555, 1139),
    }
    radial_edges = (7.5, 10.0)
    axial_edges = (2.5, 6.0, 10.5)
    for page, image in images.items():
        draw = ImageDraw.Draw(image)
        left, right, bottom, top = panels[page]
        draw.rectangle((left, top, right, bottom), outline=(255, 255, 255), width=3)
        for radius in radial_edges:
            x = _coordinate(radius, left, right, 11.25)
            draw.line((x, top, x, bottom), fill=(255, 255, 255), width=4)
        for axial in axial_edges:
            y = _coordinate(axial, bottom, top, 12.0)
            draw.line((left, y, right, y), fill=(255, 255, 255), width=4)
        image.save(output_directory / f"tian_page_{page}_spatial_zone_QA.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--pdf", type=Path)
    parser.add_argument("--qa-directory", type=Path)
    args = parser.parse_args()
    if args.verify and OUTPUT.read_text(encoding="utf-8") != text():
        raise RuntimeError("committed Tian spatial moments are not reproduced")
    if args.pdf or args.qa_directory:
        if not (args.pdf and args.qa_directory):
            raise ValueError("--pdf and --qa-directory are required together")
        qa_overlays(args.pdf, args.qa_directory)
    if not args.verify:
        print(text(), end="")


if __name__ == "__main__":
    main()
