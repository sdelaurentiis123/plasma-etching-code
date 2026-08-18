#!/usr/bin/env python3
"""Extract the author-tabulated CHF3 collision set and NIST swarm table.

The Kushner--Zhang workbook is a public author-supplied source artifact.  It
is not committed here; this script checksum-gates it and converts its cm^2
tables to an auditable long-form SI table.  The NIST drift values are a
visually checked transcription of Table 6, not OCR output.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KZ_DIR = ROOT / "data" / "experimental" / "kushner_zhang_2000_chf3"
NIST_DIR = (
    ROOT / "data" / "experimental"
    / "christophorou_olthoff_1999_chf3"
)
KZ_CSV = KZ_DIR / "cross_sections.csv"
KZ_MANIFEST = KZ_DIR / "extraction_manifest.json"
NIST_CSV = NIST_DIR / "table6_drift_velocity.csv"
NIST_TABLE4_CSV = NIST_DIR / "table4_total_scattering.csv"
NIST_TABLE5_CSV = NIST_DIR / "table5_momentum_transfer.csv"
NIST_MANIFEST = NIST_DIR / "extraction_manifest.json"

KZ_XLS_SHA256 = (
    "5835f36cb7d0e01c7f2c9c8e8206beccf51f52420f58b6ee4b0c8072f5750c87"
)
KZ_LEGEND_PDF_SHA256 = (
    "47cce6d46b67a791d694273ae591397e3ba340884a49d4ddea7e50e8399bd6e8"
)
KZ_PAPER_PDF_SHA256 = (
    "5a6ea61d002f0c2d6385e77869f33ef5cc4968bbf5d41a94760efa59037f61a2"
)
NIST_PDF_SHA256 = (
    "2699b108bf21eb292f6d49a948d3d867ee49341fc831701085b0495e9b1a9ffa"
)
NIST_PAGE11_RENDER_SHA256 = (
    "3179fabd1990d9e85e9088095ef73df69872835a8430dfbe73323c0838142ca6"
)
NIST_PAGE8_RENDER_SHA256 = (
    "5ee4792cbd2c913fb1a73c2d13cd43e29b5ce7c4c7590b4fbcdcf566238c4b1c"
)

# label, energy column, sigma column, kind, product, energy loss, delta-e.
# Thresholds and products follow the author's companion chf3_xsec.pdf.
PROCESS_DEFINITIONS = (
    ("MOM", 0, 1, "MOMENTUM", "CHF3", None, 0),
    ("VIB14", 3, 4, "EXCITATION", "CHF3(v1,v4)", 0.37, 0),
    ("VIB25", 3, 5, "EXCITATION", "CHF3(v2,v5)", 0.18, 0),
    ("VIB36", 3, 6, "EXCITATION", "CHF3(v3,v6)", 0.13, 0),
    ("NEU1", 8, 9, "EXCITATION", "CF3 + H", 11.0, 0),
    ("NEU2", 8, 10, "EXCITATION", "CHF2 + F", 13.0, 0),
    ("NEU3", 8, 11, "EXCITATION", "CF2 + H + F", 23.6, 0),
    ("NEU4", 8, 12, "EXCITATION", "CHF + F + F", 35.0, 0),
    ("NEU5", 8, 13, "EXCITATION", "CF + H + F + F", 19.5, 0),
    ("NEU6", 8, 14, "EXCITATION", "CF + H + F2", 19.5, 0),
    ("NEU1_ADD_ON", 8, 15, "EXCITATION", "CF3 + H", 11.0, 0),
    ("ION1", 17, 18, "IONIZATION", "CF3+ + H", 15.2, 1),
    ("ION2", 17, 19, "IONIZATION", "CHF2+ + F", 16.8, 1),
    ("ION3", 17, 20, "IONIZATION", "CF2+ + HF", 17.6, 1),
    ("ION4", 17, 21, "IONIZATION", "CHF+ + F + F", 19.8, 1),
    ("ION5", 17, 22, "IONIZATION", "CF+ + HF + F", 20.9, 1),
    ("ION6", 17, 23, "IONIZATION", "CH+ + F2 + F", 33.5, 1),
    ("ION7", 17, 24, "IONIZATION", "F+ + CHF2", 37.0, 1),
    ("ATT1", 26, 27, "ATTACHMENT", "CHF2 + F-", 0.0, -1),
    # ATT2 is ion-pair production: the incident electron remains.  It is an
    # energy-loss channel, not net electron attachment.
    ("ATT2", 26, 28, "EXCITATION", "CHF2+ + F-", 11.5, 0),
)

NIST_TABLE6_ROWS = (
    (0.40, 0.022), (0.45, 0.024), (0.50, 0.026), (0.60, 0.030),
    (0.70, 0.034), (0.80, 0.038), (0.90, 0.042), (1.00, 0.046),
    (1.50, 0.065), (2.00, 0.085), (2.50, 0.105), (3.00, 0.125),
    (4.00, 0.166), (5.00, 0.208), (6.00, 0.252), (7.00, 0.296),
    (8.00, 0.342), (9.00, 0.390), (10.0, 0.440), (15.0, 0.720),
    (20.0, 1.09), (25.0, 1.52), (30.0, 2.02), (40.0, 3.38),
    (50.0, 4.98), (60.0, 6.10), (70.0, 6.92), (80.0, 7.49),
    (90.0, 8.12), (100.0, 8.66), (150.0, 11.6), (200.0, 14.3),
    (250.0, 16.9),
)

NIST_TABLE4_ROWS = (
    (0.005, 3321.2), (0.006, 2767.6), (0.007, 2372.3),
    (0.008, 2075.8), (0.009, 1845.1), (0.010, 1660.6),
    (0.020, 830.3), (0.030, 553.5), (0.040, 415.2), (0.050, 332.1),
    (0.060, 276.8), (0.070, 237.2), (0.080, 207.6), (0.090, 184.5),
    (0.100, 166.1), (0.200, 82.9), (0.300, 55.4), (0.400, 41.5),
    (0.500, 35.0), (0.600, 31.3), (0.700, 29.2), (0.800, 27.7),
    (0.900, 26.5), (1.00, 25.6), (1.50, 23.5), (2.00, 22.1),
    (3.00, 20.3), (4.00, 19.9), (5.00, 20.5), (6.00, 21.3),
    (7.00, 21.9), (8.00, 21.8), (9.00, 21.4), (10.0, 20.6),
    (15.0, 18.7), (20.0, 19.0), (30.0, 18.6), (40.0, 17.7),
    (50.0, 16.7), (60.0, 15.9), (70.0, 15.0), (80.0, 14.2),
    (90.0, 13.4), (100.0, 12.7), (200.0, 8.8), (300.0, 7.0),
    (400.0, 5.9), (500.0, 5.2), (600.0, 4.6),
)

NIST_TABLE5_ROWS = (
    (10.0, 15.5), (15.0, 12.4), (20.0, 11.4), (25.0, 10.5),
    (30.0, 10.1),
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _csv_text(fields, rows) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(fields)
    writer.writerows(rows)
    return stream.getvalue()


def _verify(path: Path, expected: str, label: str) -> None:
    if _sha(path) != expected:
        raise RuntimeError(f"{label} checksum changed")


def extract_kz_xls(path: Path) -> str:
    _verify(path, KZ_XLS_SHA256, "Kushner--Zhang workbook")
    try:
        import xlrd
    except ImportError as exc:
        raise RuntimeError("xlrd is required only to replay the source XLS") from exc
    sheet = xlrd.open_workbook(path).sheet_by_index(0)
    if sheet.nrows != 208 or sheet.ncols != 29:
        raise RuntimeError("Kushner--Zhang workbook shape changed")
    rows = []
    for label, energy_col, sigma_col, kind, product, loss, delta_e in (
        PROCESS_DEFINITIONS
    ):
        count = 0
        for row_index in range(7, sheet.nrows):
            energy = sheet.cell_value(row_index, energy_col)
            sigma = sheet.cell_value(row_index, sigma_col)
            if isinstance(energy, (int, float)) and isinstance(
                sigma, (int, float)
            ):
                rows.append((
                    label, kind, "CHF3", product,
                    "" if loss is None else format(loss, ".12g"), delta_e,
                    format(float(energy), ".12g"),
                    format(float(sigma) * 1.0e-4, ".12g"),
                    energy_col, sigma_col, "swarm_regressed_working_set",
                ))
                count += 1
        if count < 2:
            raise RuntimeError(f"missing source rows for {label}")
    return _csv_text((
        "process_label", "kind", "target", "product", "energy_loss_eV",
        "electron_number_change", "energy_eV", "cross_section_m2",
        "source_energy_column_zero_based",
        "source_cross_section_column_zero_based", "evidence_class",
    ), rows)


def nist_table6_text() -> str:
    return _csv_text((
        "reduced_electric_field_Td", "drift_velocity_m_s",
        "source_printed_velocity_1e6_cm_s", "gas_temperature_K",
        "evidence_class",
    ), tuple(
        (format(field, ".12g"), format(velocity * 1.0e4, ".12g"),
         format(velocity, ".12g"), "298", "recommended_measured_swarm_fit")
        for field, velocity in NIST_TABLE6_ROWS
    ))


def nist_cross_section_text(rows, evidence_class: str) -> str:
    return _csv_text((
        "electron_energy_eV", "cross_section_m2",
        "source_printed_cross_section_1e_minus_20_m2", "evidence_class",
    ), tuple(
        (format(energy, ".12g"), format(value * 1.0e-20, ".12g"),
         format(value, ".12g"), evidence_class)
        for energy, value in rows
    ))


def kz_manifest_text(csv_payload: str) -> str:
    process_rows = {}
    for label, *_ in PROCESS_DEFINITIONS:
        process_rows[label] = sum(
            1 for row in csv_payload.splitlines()[1:]
            if row.split(",", 1)[0] == label
        )
    data = {
        "manifest_id": "KUSHNER-ZHANG-2000-CHF3-XSEC-R1",
        "source": {
            "citation": "M. J. Kushner and D. Zhang, J. Appl. Phys. 88, 3231 (2000)",
            "doi": "10.1063/1.1289187",
            "author_data_url": "https://cpseg.eecs.umich.edu/pub/data/chf3_xsec.xls",
            "author_legend_url": "https://cpseg.eecs.umich.edu/pub/data/chf3_xsec.pdf",
            "xls_sha256": KZ_XLS_SHA256,
            "legend_pdf_sha256": KZ_LEGEND_PDF_SHA256,
            "paper_pdf_sha256": KZ_PAPER_PDF_SHA256,
            "source_artifacts_committed": False,
        },
        "conversion": "author cm^2 values multiplied by 1e-4 exactly",
        "process_row_counts": process_rows,
        "output": {
            "path": str(KZ_CSV.relative_to(ROOT)),
            "sha256": sha256(csv_payload.encode()).hexdigest(),
        },
        "claim_boundary": {
            "valid": [
                "exact replay of the author-tabulated working set",
                "deterministic two-term source-reproduction calculations",
                "species-resolved electron-impact source terms within support",
            ],
            "not_valid": [
                "independent validation against the swarm data used to regress the set",
                "direct measurement of neutral-dissociation branches",
                "a unique target-reactor plasma state or wafer flux",
                "an atomic-accuracy feature-depth claim",
            ],
            "notable_source_limits": [
                "neutral dissociation was scaled/regressed to swarm data",
                "CF neutral branching was arbitrarily divided equally",
                "attachment is weak and uncertain",
                "published total-ionization measurements disagree by about 2x",
            ],
        },
    }
    return json.dumps(data, indent=2) + "\n"


def nist_manifest_text(
    drift_payload: str,
    total_scattering_payload: str,
    momentum_payload: str,
) -> str:
    data = {
        "manifest_id": "CHRISTOPHOROU-OLTHOFF-1999-CHF3-TABLE6-R1",
        "source": {
            "citation": "L. G. Christophorou and J. K. Olthoff, J. Phys. Chem. Ref. Data 28, 967 (1999)",
            "doi": "10.1063/1.556050",
            "nist_url": "https://www.nist.gov/publications/electron-interactions-plasma-processing-gases-update-cf4-chf3-c2f6-and-c3f8",
            "pdf_sha256": NIST_PDF_SHA256,
        },
        "visual_audit": {
            "render_dpi": 400,
            "render_size_px": [3400, 4406],
            "renders": {
                "pdf_page_8_tables_4_5_sha256": NIST_PAGE8_RENDER_SHA256,
                "pdf_page_11_table_6_sha256": NIST_PAGE11_RENDER_SHA256,
            },
            "status": "passed_original_resolution",
        },
        "conversions": [
            "printed 10^6 cm/s values multiplied by 1e4 to m/s",
            "printed 1e-20 m^2 cross sections multiplied by 1e-20",
        ],
        "outputs": [
            {"path": str(NIST_TABLE4_CSV.relative_to(ROOT)),
             "sha256": sha256(total_scattering_payload.encode()).hexdigest(),
             "row_count": len(NIST_TABLE4_ROWS)},
            {"path": str(NIST_TABLE5_CSV.relative_to(ROOT)),
             "sha256": sha256(momentum_payload.encode()).hexdigest(),
             "row_count": len(NIST_TABLE5_ROWS)},
            {"path": str(NIST_CSV.relative_to(ROOT)),
             "sha256": sha256(drift_payload.encode()).hexdigest(),
             "row_count": len(NIST_TABLE6_ROWS)},
        ],
        "claim_boundary": {
            "valid": [
                "recommended total electron-scattering curve",
                "suggested calculated momentum-transfer values from 10-30 eV",
                "recommended 298 K pure-CHF3 drift-velocity curve",
            ],
            "not_valid": [
                "a statistically independent grade of a swarm-regressed deck",
                "target-reactor power coupling, plasma density, or wafer flux",
            ],
        },
    }
    return json.dumps(data, indent=2) + "\n"


def _payloads(kz_csv: str | None = None) -> dict[Path, str]:
    if kz_csv is None:
        kz_csv = KZ_CSV.read_text(encoding="utf-8")
    nist_csv = nist_table6_text()
    total_scattering = nist_cross_section_text(
        NIST_TABLE4_ROWS, "recommended_total_scattering")
    momentum = nist_cross_section_text(
        NIST_TABLE5_ROWS, "suggested_calculated_momentum_transfer")
    return {
        KZ_CSV: kz_csv,
        KZ_MANIFEST: kz_manifest_text(kz_csv),
        NIST_TABLE4_CSV: total_scattering,
        NIST_TABLE5_CSV: momentum,
        NIST_CSV: nist_csv,
        NIST_MANIFEST: nist_manifest_text(
            nist_csv, total_scattering, momentum),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-xls", type=Path)
    parser.add_argument("--legend-pdf", type=Path)
    parser.add_argument("--paper-pdf", type=Path)
    parser.add_argument("--nist-pdf", type=Path)
    parser.add_argument("--nist-page-render", type=Path)
    parser.add_argument("--nist-page-8-render", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.legend_pdf is not None:
        _verify(args.legend_pdf, KZ_LEGEND_PDF_SHA256, "CHF3 legend PDF")
    if args.paper_pdf is not None:
        _verify(args.paper_pdf, KZ_PAPER_PDF_SHA256, "CHF3 paper PDF")
    if args.nist_pdf is not None:
        _verify(args.nist_pdf, NIST_PDF_SHA256, "NIST CHF3 PDF")
    if args.nist_page_render is not None:
        _verify(args.nist_page_render, NIST_PAGE11_RENDER_SHA256,
                "NIST Table 6 render")
    if args.nist_page_8_render is not None:
        _verify(args.nist_page_8_render, NIST_PAGE8_RENDER_SHA256,
                "NIST Tables 4/5 render")
    kz_csv = extract_kz_xls(args.source_xls) if args.source_xls else None
    payloads = _payloads(kz_csv)
    if args.check:
        for path, expected in payloads.items():
            if path.read_text(encoding="utf-8") != expected:
                raise RuntimeError(f"committed CHF3 evidence is stale: {path}")
        return
    if kz_csv is None:
        raise SystemExit("initial extraction requires --source-xls")
    for path, payload in payloads.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
