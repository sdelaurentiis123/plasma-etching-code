#!/usr/bin/env python3
"""Transcribe the pure-Cl2 swarm tables from Gonzalez-Magana (2018).

The source publishes native numeric tables, so this is a visual table
transcription rather than curve digitization. The exact publisher PDF and the
three 300-dpi table-page renders are hash gated. Only the derived CSV and audit
manifest are committed; the copyrighted PDF and renders remain local.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDF_SHA256 = (
    "fc42a2fede6094acc8344780a790ac6e0dae0d55cce2dcad6b8f5cae2a6dba5d"
)
SOURCE_RENDER_SHA256 = {
    6: "3e9d7c6c26b1d2551ee4e08091c8d71511d5ab79ecbfb5248c92ab8104a46ae8",
    7: "d37c469236097421437f8e8f3ad42e83c4de0282f96726847c482288bbe614ff",
    8: "4cace83fa2e94ad78dd0f9f4ccd076b1ebe76331a06a513db4ffaae7ce801eed",
}
OUTPUT_PATH = (
    ROOT / "research_sources" / "digitized"
    / "gonzalez_magana_2018_pure_cl2_swarm.csv"
)
PACKAGE_OUTPUT_PATH = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "gonzalez_magana_2018_pure_cl2_swarm.csv"
)
MANIFEST_PATH = (
    ROOT / "research_sources" / "digitized"
    / "gonzalez_magana_2018_pure_cl2_swarm_manifest.json"
)

# (E/N in Td, value in the table's printed unit).
DRIFT_VELOCITY_ROWS = (
    (100, "8.93"),
    (120, "9.97"),
    (140, "11.24"),
    (160, "11.46"),
    (180, "12.55"),
    (200, "13.63"),
    (220, "14.38"),
    (240, "15.12"),
    (260, "16.27"),
    (270, "16.60"),
    (280, "17.10"),
    (290, "17.66"),
    (300, "18.09"),
    (310, "18.59"),
    (320, "19.24"),
    (330, "19.82"),
    (340, "20.52"),
    (360, "21.72"),
    (380, "22.61"),
    (400, "23.66"),
    (420, "24.63"),
    (440, "25.71"),
    (460, "26.83"),
)
EFFECTIVE_IONIZATION_ROWS = (
    (100, "-2.540"),
    (120, "-2.380"),
    (140, "-2.260"),
    (160, "-2.180"),
    (180, "-1.940"),
    (200, "-1.560"),
    (220, "-1.260"),
    (240, "-0.847"),
    (260, "-0.452"),
    (270, "-0.189"),
    (280, "0.002"),
    (290, "0.237"),
    (300, "0.491"),
    (310, "0.890"),
    (320, "1.096"),
    (340, "1.618"),
    (360, "2.220"),
    (380, "2.789"),
    (400, "3.490"),
    (420, "4.000"),
    (440, "4.610"),
)
LONGITUDINAL_DIFFUSION_ROWS = (
    (240, "1.38"),
    (270, "1.59"),
    (300, "1.96"),
    (310, "1.87"),
    (320, "2.05"),
    (330, "1.98"),
    (340, "2.07"),
    (360, "2.26"),
)

FIELDNAMES = (
    "observation_id",
    "observable",
    "reduced_field_Td",
    "value_si",
    "si_unit",
    "published_value",
    "published_unit",
    "relative_uncertainty_min",
    "relative_uncertainty_max",
    "source_table",
    "source_pdf_page",
    "source_print_page",
    "measurement_method",
    "gas_temperature_K_min",
    "gas_temperature_K_max",
    "pressure_Torr_min",
    "pressure_Torr_max",
    "supports_cross_section_validation",
    "supports_reactor_state_prediction",
    "supports_wafer_flux",
    "supports_feature_depth",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _scientific(value: Decimal) -> str:
    return format(value, ".6E").replace("E", "e")


def _common_row() -> dict[str, str]:
    return {
        "measurement_method": "pulsed_townsend_transient",
        "gas_temperature_K_min": "293",
        "gas_temperature_K_max": "298",
        "pressure_Torr_min": "1.1",
        "pressure_Torr_max": "6.1",
        "supports_cross_section_validation": "true",
        "supports_reactor_state_prediction": "false",
        "supports_wafer_flux": "false",
        "supports_feature_depth": "false",
    }


def rows() -> tuple[dict[str, str], ...]:
    output: list[dict[str, str]] = []
    for reduced_field, published in DRIFT_VELOCITY_ROWS:
        output.append({
            "observation_id": f"gonzalez_2018_cl2_W_{reduced_field:03d}Td",
            "observable": "electron_drift_velocity",
            "reduced_field_Td": str(reduced_field),
            "value_si": _scientific(Decimal(published) * Decimal("1e4")),
            "si_unit": "m s^-1",
            "published_value": published,
            "published_unit": "10^6 cm s^-1",
            "relative_uncertainty_min": "0.02",
            "relative_uncertainty_max": "0.02",
            "source_table": "Table A1",
            "source_pdf_page": "6" if reduced_field <= 240 else "7",
            "source_print_page": "5" if reduced_field <= 240 else "6",
            **_common_row(),
        })
    for reduced_field, published in EFFECTIVE_IONIZATION_ROWS:
        output.append({
            "observation_id": (
                f"gonzalez_2018_cl2_alpha_minus_eta_{reduced_field:03d}Td"
            ),
            "observable": "effective_ionization_coefficient",
            "reduced_field_Td": str(reduced_field),
            "value_si": _scientific(Decimal(published) * Decimal("1e-21")),
            "si_unit": "m^2",
            "published_value": published,
            "published_unit": "10^-17 cm^2",
            "relative_uncertainty_min": "0.05",
            "relative_uncertainty_max": "0.09",
            "source_table": "Table A2",
            "source_pdf_page": "7",
            "source_print_page": "6",
            **_common_row(),
        })
    for reduced_field, published in LONGITUDINAL_DIFFUSION_ROWS:
        output.append({
            "observation_id": (
                f"gonzalez_2018_cl2_NDL_{reduced_field:03d}Td"
            ),
            "observable": "density_normalized_longitudinal_diffusion",
            "reduced_field_Td": str(reduced_field),
            "value_si": _scientific(Decimal(published) * Decimal("1e24")),
            "si_unit": "m^-1 s^-1",
            "published_value": published,
            "published_unit": "10^22 cm^-1 s^-1",
            "relative_uncertainty_min": "0.10",
            "relative_uncertainty_max": "0.15",
            "source_table": "Table A3",
            "source_pdf_page": "8",
            "source_print_page": "7",
            **_common_row(),
        })
    return tuple(output)


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    return {
        "source": {
            "bibkey": "gonzalez-magana-de-urquijo-2018-cl2",
            "doi": "10.1088/1361-6595/aac95e",
            "publisher_pdf_sha256": SOURCE_PDF_SHA256,
            "publisher_pdf_committed": False,
        },
        "visual_audit": {
            "status": "passed",
            "inspection_date": "2026-08-08",
            "render_resolution_dpi": 300,
            "pages": [
                {
                    "pdf_page": page,
                    "print_page": page - 1,
                    "render_sha256": digest,
                }
                for page, digest in SOURCE_RENDER_SHA256.items()
            ],
            "finding": (
                "All 52 pure-Cl2 numeric cells were read from Tables A1-A3; "
                "signs, decimal precision, headings, and units were checked "
                "against the 300-dpi page renders."
            ),
        },
        "transcription": {
            "method": "manual_native_table_transcription_with_pixel_audit",
            "row_count": len(rows()),
            "observable_counts": {
                "electron_drift_velocity": len(DRIFT_VELOCITY_ROWS),
                "effective_ionization_coefficient": len(
                    EFFECTIVE_IONIZATION_ROWS),
                "density_normalized_longitudinal_diffusion": len(
                    LONGITUDINAL_DIFFUSION_ROWS),
            },
            "csv_sha256": csv_sha256,
            "si_conversions": {
                "10^6 cm s^-1_to_m s^-1": "multiply by 1e4",
                "10^-17 cm^2_to_m^2": "multiply by 1e-21",
                "10^22 cm^-1 s^-1_to_m^-1 s^-1": "multiply by 1e24",
            },
        },
        "measurement_uncertainty": {
            "electron_drift_velocity_relative": [0.02, 0.02],
            "effective_ionization_coefficient_relative": [0.05, 0.09],
            "density_normalized_longitudinal_diffusion_relative": [0.10, 0.15],
            "semantics": (
                "paper-wide typical ranges; no point-specific sigma is "
                "invented"
            ),
        },
        "source_conflicts_preserved": {
            "pure_cl2_range": (
                "body says 100-420 Td; Table A1 prints W through 460 Td"
            ),
            "resolution": (
                "all printed table markers are retained and the conflict is "
                "not silently repaired"
            ),
        },
        "use_boundary": {
            "supports_electron_cross_section_validation": True,
            "supports_reactor_state_prediction": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
            "reason": (
                "swarm transport constrains the electron-collision deck; it "
                "does not measure a plasma-reactor state or wafer boundary"
            ),
        },
        "rights": (
            "Only derived table values and audit metadata are committed; the "
            "publisher PDF and rendered pages remain local."
        ),
    }


def _verify_renders(source_pdf: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="gonzalez-2018-render-") as name:
        prefix = Path(name) / "page"
        subprocess.run(
            [
                "pdftoppm", "-f", "6", "-l", "8", "-png", "-r", "300",
                str(source_pdf), str(prefix),
            ],
            check=True,
        )
        for page, expected in SOURCE_RENDER_SHA256.items():
            rendered = prefix.with_name(f"page-{page}.png")
            if _sha256(rendered) != expected:
                raise RuntimeError(
                    f"source page {page} render does not match visual audit"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument(
        "--package-output", type=Path, default=PACKAGE_OUTPUT_PATH)
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    arguments = parser.parse_args()

    if _sha256(arguments.pdf) != SOURCE_PDF_SHA256:
        raise RuntimeError("publisher PDF hash does not match audited source")
    _verify_renders(arguments.pdf)

    payload = csv_text()
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    for path in (arguments.output, arguments.package_output):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(payload, encoding="utf-8")
    arguments.manifest.parent.mkdir(parents=True, exist_ok=True)
    arguments.manifest.write_text(
        json.dumps(manifest(digest), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
