#!/usr/bin/env python3
"""Replay the Depla--Van Bever TiO2 Ar-sputter reference curve.

Figure 1 reports *total atoms emitted per incident Ar ion*.  The red TiO2
curve is a semi-empirical fit, not a direct low-energy beam measurement.  This
replay digitizes only the 125--400 eV part of that fitted curve and also emits
the stoichiometric formula-unit conversion (three atoms per TiO2 unit).

The distinction matters: this is a useful bare-oxide physical-sputter scale,
but it is not the reactive CHF3/SF6/O2 removal law needed by the Oxford run.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
import math
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "experimental" / "depla_2024_tio2_sputter"
CSV_PATH = OUTPUT_DIR / "figure1_tio2_ar_sputter_fit.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "af5e2530985073d01c401dd0107797080b65fd24208c7d18805a9f0f628ab5e3"
)
FIGURE1_RENDER_SHA256 = (
    "4d402bb4a0ebaa2b6491cbf0dd5d4e658379e1a38e7d5bdce33b6359710c5c86"
)
FIGURE2_RENDER_SHA256 = (
    "c96ed3c82b9f9c847c8c5cab81a43500ce51f43e2ba4cc6265587c6bb92909da"
)
FIGURE1_CROP_SHA256 = (
    "7bae96dd92a4725f27d1736a19e92876b9490ef4fabd912b764c7b699963548b"
)
RENDER_SIZE = (3400, 4400)
CROP_SIZE = (1200, 1450)
CROP_BOX_FULL_RENDER = (2050, 2350, 3250, 3800)

# Log-log axes in the 400-dpi crop.  A 280-pixel displacement is one decade.
X_AT_100_EV = 19.0
X_AT_1_KEV = 299.0
X_AT_10_KEV = 579.0
Y_AT_0P1_ATOM_PER_ION = 585.0
Y_AT_1_ATOM_PER_ION = 305.0
Y_AT_10_ATOM_PER_ION = 25.0

# Cubic least-squares trace of the pure-red TiO2 fitted curve in crop pixels.
# Support was restricted to x=44..183 pixels; marker-contaminated columns were
# rejected by a median-absolute-residual filter.  This polynomial is only a
# smooth pixel interpolation, not a sputtering model.
CURVE_PIXEL_POLYNOMIAL = (
    -5.35989911e-06,
    5.34245709e-03,
    -2.33135012,
    6.11181418e02,
)
FIT_SUPPORT_X_PX = (44, 183)
FIT_KEPT_COLUMNS = 93
FIT_PIXEL_MAD = 0.789
ENERGIES_EV = (125, 150, 175, 200, 225, 250, 276, 300, 350, 400)

FIELDS = (
    "argon_ion_energy_eV",
    "curve_x_crop_px",
    "curve_y_crop_px",
    "fitted_total_atom_yield_per_ion",
    "stoichiometric_tio2_formula_units_per_ion",
    "digitization_relative_uncertainty",
    "evidence_class",
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _x_from_energy(energy_eV: float) -> float:
    return X_AT_100_EV + 280.0 * math.log10(energy_eV / 100.0)


def _y_from_x(x_px: float) -> float:
    a3, a2, a1, a0 = CURVE_PIXEL_POLYNOMIAL
    return ((a3 * x_px + a2) * x_px + a1) * x_px + a0


def _yield_from_y(y_px: float) -> float:
    return 0.1 * 10.0 ** ((Y_AT_0P1_ATOM_PER_ION - y_px) / 280.0)


def rows() -> list[dict[str, str]]:
    result = []
    for energy in ENERGIES_EV:
        x_px = _x_from_energy(energy)
        y_px = _y_from_x(x_px)
        atom_yield = _yield_from_y(y_px)
        result.append({
            "argon_ion_energy_eV": str(energy),
            "curve_x_crop_px": f"{x_px:.6f}",
            "curve_y_crop_px": f"{y_px:.6f}",
            "fitted_total_atom_yield_per_ion": f"{atom_yield:.6f}",
            "stoichiometric_tio2_formula_units_per_ion": (
                f"{atom_yield / 3.0:.6f}"
            ),
            # Three pixels includes line thickness, curve isolation, and axis
            # center placement: 10**(3/280)-1 = 2.50%, rounded upward.
            "digitization_relative_uncertainty": "0.03",
            "evidence_class": "digitized_semiempirical_fit_curve",
        })
    return result


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha: str) -> dict[str, object]:
    by_energy = {int(row["argon_ion_energy_eV"]): row for row in rows()}
    return {
        "manifest_id": "DEPLA-2024-TIO2-AR-SPUTTER-FIT-R1",
        "source": {
            "citation": (
                "D. Depla and J. Van Bever, Calculation of oxide sputter "
                "yields, Vacuum 222 (2024) 112994"
            ),
            "doi": "10.1016/j.vacuum.2024.112994",
            "primary_record": (
                "https://biblio.ugent.be/publication/01HN0EVQ4K9X7ZKEMNFYB2R2RP"
            ),
            "source_pdf_sha256": SOURCE_PDF_SHA256,
            "source_pdf_committed": False,
        },
        "visual_audit": {
            "render_dpi": 400,
            "render_size_px": list(RENDER_SIZE),
            "figure1_pdf_page": 10,
            "figure1_render_sha256": FIGURE1_RENDER_SHA256,
            "figure1_crop_box_full_render_px": list(CROP_BOX_FULL_RENDER),
            "figure1_crop_size_px": list(CROP_SIZE),
            "figure1_crop_sha256": FIGURE1_CROP_SHA256,
            "figure2_pdf_page": 15,
            "figure2_render_sha256": FIGURE2_RENDER_SHA256,
            "status": "passed_original_resolution",
        },
        "quantity": {
            "projectile": "Ar+",
            "target": "stoichiometric TiO2 reference",
            "source_vertical_axis": "total target atoms emitted per incident ion",
            "formula_unit_conversion": "divide total atom yield by 3",
            "incidence_angle": (
                "source compilation normalizes literature oblique data using "
                "its stated <=35 degree inverse-cosine approximation"
            ),
        },
        "source_evidence": {
            "fit_model": "Seah-Nunney semi-empirical model",
            "fit_free_parameter": "oxygen surface binding energy",
            "tio2_oxygen_surface_binding_energy_eV": 3.52,
            "tio2_oxygen_surface_binding_energy_95pct_eV": 1.14,
            "low_energy_input": (
                "Crist 2023 reports etch rates from 0.2--3 keV without ion "
                "current density; Depla and Van Bever infer current density "
                "from metal references before converting to yield"
            ),
            "high_energy_input": "Bach 1974 crystalline-target IBE data",
        },
        "pixel_replay": {
            "axis_type": "log10-log10",
            "x_calibration": {
                "100_eV_px": X_AT_100_EV,
                "1000_eV_px": X_AT_1_KEV,
                "10000_eV_px": X_AT_10_KEV,
            },
            "y_calibration": {
                "0.1_atom_per_ion_px": Y_AT_0P1_ATOM_PER_ION,
                "1_atom_per_ion_px": Y_AT_1_ATOM_PER_ION,
                "10_atom_per_ion_px": Y_AT_10_ATOM_PER_ION,
            },
            "curve_rgb": [255, 0, 0],
            "curve_pixel_polynomial_high_to_low": list(CURVE_PIXEL_POLYNOMIAL),
            "fit_support_x_px": list(FIT_SUPPORT_X_PX),
            "kept_columns": FIT_KEPT_COLUMNS,
            "pixel_residual_mad": FIT_PIXEL_MAD,
            "digitization_relative_uncertainty": 0.03,
        },
        "derived_checks": {
            "tio2_formula_units_per_ion_at_200_eV": float(
                by_energy[200]["stoichiometric_tio2_formula_units_per_ion"]
            ),
            "tio2_formula_units_per_ion_at_276_eV": float(
                by_energy[276]["stoichiometric_tio2_formula_units_per_ion"]
            ),
            "tio2_formula_units_per_ion_at_400_eV": float(
                by_energy[400]["stoichiometric_tio2_formula_units_per_ion"]
            ),
        },
        "output": {
            "path": str(CSV_PATH.relative_to(ROOT)),
            "sha256": csv_sha,
        },
        "claim_boundary": {
            "valid": [
                "bare stoichiometric TiO2 Ar-sputter reference scale",
                "shape of the source's semi-empirical TiO2 fit over 125--400 eV",
                "evidence that pure physical sputtering is much smaller than the Oxford effective-removal requirement",
            ],
            "not_valid": [
                "direct low-energy TiO2 beam-yield measurement",
                "reactive CHF3/SF6/O2 TiO2 yield coefficient",
                "CFx+, HF+, SFx+, O+, or mixed-ion sputter yield",
                "ALD-film-specific yield or Oxford NPG80 target coefficient",
                "lower or upper bound on ion-assisted chemical removal",
            ],
        },
    }


def manifest_text(payload: str) -> str:
    return json.dumps(
        manifest(sha256(payload.encode("utf-8")).hexdigest()), indent=2
    ) + "\n"


def verify_pixels(source_pdf: Path, figure1: Path, crop: Path, figure2: Path) -> None:
    expected = (
        (source_pdf, SOURCE_PDF_SHA256, None),
        (figure1, FIGURE1_RENDER_SHA256, RENDER_SIZE),
        (crop, FIGURE1_CROP_SHA256, CROP_SIZE),
        (figure2, FIGURE2_RENDER_SHA256, RENDER_SIZE),
    )
    for path, digest, size in expected:
        if _sha(path) != digest:
            raise RuntimeError(f"source checksum changed: {path}")
        if size is not None and Image.open(path).size != size:
            raise RuntimeError(f"source image size changed: {path}")
    with Image.open(figure1).convert("RGB") as full, Image.open(crop).convert("RGB") as cut:
        if ImageChops.difference(full.crop(CROP_BOX_FULL_RENDER), cut).getbbox():
            raise RuntimeError("Figure 1 crop no longer matches pinned full render")


def draw_overlay(crop: Path, output: Path) -> None:
    image = Image.open(crop).convert("RGB")
    draw = ImageDraw.Draw(image)
    for row in rows():
        x_px = float(row["curve_x_crop_px"])
        y_px = float(row["curve_y_crop_px"])
        draw.ellipse((x_px - 4, y_px - 4, x_px + 4, y_px + 4), outline="blue", width=2)
    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--figure1-render", type=Path)
    parser.add_argument("--figure1-crop", type=Path)
    parser.add_argument("--figure2-render", type=Path)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")

    pixel_args = (
        args.source_pdf, args.figure1_render, args.figure1_crop, args.figure2_render
    )
    if any(value is not None for value in pixel_args):
        if any(value is None for value in pixel_args):
            parser.error("all four source/render paths are required together")
        verify_pixels(*pixel_args)
    if args.overlay:
        if args.figure1_crop is None:
            parser.error("--overlay requires the verified Figure 1 crop")
        draw_overlay(args.figure1_crop, args.overlay)

    payload = csv_text()
    receipt = manifest_text(payload)
    if args.write:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(payload, encoding="utf-8")
        MANIFEST_PATH.write_text(receipt, encoding="utf-8")
        print(CSV_PATH.relative_to(ROOT))
        print(MANIFEST_PATH.relative_to(ROOT))
        return
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise SystemExit("Depla TiO2 digitized curve is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != receipt:
        raise SystemExit("Depla TiO2 digitization manifest is stale")
    print("PASS Depla 2024 TiO2 Ar-sputter fit replay")


if __name__ == "__main__":
    main()
