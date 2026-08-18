#!/usr/bin/env python3
"""Materialize visually checked NIST SF6 electron-interaction tables.

The numerical rows below are a manual transcription of Tables 9, 14, 15,
17, 20, 28, and 36 in Christophorou and Olthoff, J. Phys. Chem. Ref. Data
29, 267 (2000), doi:10.1063/1.1288407.  The source PDF is not redistributed;
when supplied, it and the 400 dpi table-page renders are checksum gated.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "data" / "experimental"
    / "christophorou_olthoff_2000_sf6"
)
PDF_SHA256 = (
    "b13fac820570646e6b72a59e573fb85ba99814f72cbb0962193b969036a54ad0"
)
TABLE_PAGE_RENDER_SHA256 = {
    "table9_pdf_page17": "503fea7ec6fa6834104fbcb6df00531592c276ce968d3f435b69a3ff344b6183",
    "table14_pdf_page23": "6f7cebb10d612883b06adca52788a656e66c5824f47b9be08e418b650c3f27bf",
    "table15_pdf_page26": "eb9641df2b3f9736618eef1e70044421fc560613edb8dabcc95193ba3acef5a6",
    "table17_pdf_page28": "a2f83b2c43162d155f22c2e6781d5f286f7e94512cd2759824f5e1750d30290a",
    "table20_pdf_page33": "82c9da2506adb36c1cfbfc247c02b8affed8bef92edc9f799505f62f40ca5838",
    "table28_pdf_page44": "7b9fb4d0adff153fd41427c1daaeff7ff59da696cb6625aa0dec94f5361e0d64",
    "table36_pdf_page52": "1ac7b23d2a85d27aeb29bb30fe2a6a542840b5f9143bf0481ed51fdcf8f27ffb",
    "table37_pdf_page53": "e81d3e81942928d01e0e3c7a80074a1f9369b4a94e2f2a51902ba565a192c97b",
}

TABLE9 = (
    (.035, 379.8), (.040, 344.7), (.045, 315.3), (.050, 290.1),
    (.060, 249.4), (.070, 217.6), (.080, 192.7), (.090, 173.1),
    (.10, 157.4), (.15, 109.5), (.20, 84.3), (.25, 68.7),
    (.30, 58.1), (.35, 50.4), (.40, 44.6), (.45, 40.3),
    (.50, 37.2), (.60, 32.8), (.70, 29.4), (.80, 26.8),
    (.90, 24.9), (1.0, 23.4), (1.5, 22.2), (2.0, 22.8),
    (2.5, 23.4), (3.0, 23.1), (3.5, 22.7), (4.0, 22.4),
    (4.5, 22.6), (5.0, 23.7), (6.0, 28.0), (7.0, 30.7),
    (8.0, 29.2), (9.0, 27.6), (10.0, 27.3), (11.0, 29.6),
    (12.0, 32.9), (13.0, 28.9), (14.0, 26.6), (15.0, 26.4),
    (18.0, 26.7), (20.0, 27.6), (22.0, 28.7), (25.0, 29.4),
    (30.0, 29.4), (35.0, 29.3), (40.0, 29.3), (45.0, 29.3),
    (50.0, 29.3), (60.0, 28.8), (70.0, 28.2), (80.0, 27.2),
    (90.0, 26.5), (100.0, 25.8), (150.0, 22.7), (200.0, 20.2),
    (250.0, 18.0), (300.0, 16.3), (350.0, 15.1), (400.0, 14.1),
    (450.0, 13.1), (500.0, 12.4), (600.0, 11.2), (700.0, 10.2),
    (800.0, 9.41), (900.0, 8.73), (1000.0, 8.15), (1500.0, 6.07),
    (2000.0, 4.80), (2500.0, 3.94), (3000.0, 3.36),
    (3500.0, 2.95), (4000.0, 2.64),
)

TABLE14 = (
    (2.75, 16.0), (3.0, 15.4), (3.5, 14.5), (4.0, 14.0),
    (4.5, 13.9), (5.0, 14.1), (6.0, 15.1), (7.0, 15.5),
    (8.0, 14.8), (9.0, 14.4), (10.0, 15.1), (11.0, 16.7),
    (12.0, 17.6), (13.0, 17.1), (14.0, 15.8), (15.0, 14.9),
    (16.0, 14.5), (17.0, 14.7), (18.0, 15.0), (19.0, 15.4),
    (20.0, 15.7), (22.0, 15.7), (25.0, 15.0), (27.0, 14.3),
    (30.0, 13.2), (35.0, 11.5), (40.0, 10.3), (45.0, 9.37),
    (50.0, 8.65), (60.0, 7.69), (70.0, 7.06), (75.0, 6.74),
    (80.0, 6.46), (90.0, 6.03), (100.0, 5.70), (125.0, 4.92),
    (150.0, 4.16), (200.0, 2.98), (250.0, 2.23), (300.0, 1.76),
    (350.0, 1.47), (400.0, 1.28), (450.0, 1.13), (500.0, 1.02),
    (600.0, .82), (700.0, .66),
)

TABLE15 = (
    (.09, 1.9), (.10, 7.0), (.12, 21.3), (.15, 30.6), (.17, 30.6),
    (.20, 34.9), (.22, 35.5), (.25, 33.6), (.28, 30.9), (.30, 29.4),
    (.35, 26.8), (.40, 25.4), (.45, 23.5), (.50, 21.6), (.60, 19.7),
    (.70, 18.1), (.80, 16.4), (.90, 15.1), (1.0, 13.9), (1.2, 12.6),
    (1.5, 10.8), (2.0, 8.0), (2.5, 5.6), (3.0, 3.8), (3.5, 2.8),
    (4.0, 2.3), (4.5, 2.0), (5.0, 2.4), (6.0, 4.4), (7.0, 6.5),
    (8.0, 4.6), (9.0, 3.1), (10.0, 2.5), (11.0, 3.5), (12.0, 6.3),
    (13.0, 2.4), (14.0, .5),
)

TABLE17 = (
    (16.5, .020), (17.0, .035), (17.5, .055), (18.0, .084),
    (19.0, .155), (20.0, .240), (22.0, .457), (25.0, 1.04),
    (30.0, 1.93), (35.0, 2.87), (40.0, 3.47), (45.0, 3.79),
    (50.0, 4.35), (60.0, 5.09), (70.0, 5.65), (75.0, 5.77),
    (80.0, 5.95), (90.0, 6.28), (100.0, 6.53), (125.0, 6.87),
    (150.0, 6.97), (200.0, 6.83), (250.0, 6.48), (300.0, 6.04),
    (350.0, 5.60), (400.0, 5.16), (450.0, 4.75), (500.0, 4.36),
    (550.0, 4.00), (600.0, 3.65),
)

TABLE20 = (
    (15.0, .8), (16.0, 1.2), (17.0, 1.6), (18.0, 1.8),
    (19.0, 2.2), (20.0, 2.7), (22.0, 3.5), (25.0, 3.7),
    (30.0, 3.1), (35.0, 2.4), (40.0, 2.3), (45.0, 2.7),
    (50.0, 2.8), (60.0, 2.4), (70.0, 2.1), (75.0, 1.7),
    (80.0, 1.5), (90.0, 1.1), (100.0, .87), (125.0, .53),
    (150.0, .23), (200.0, .07),
)

TABLE28 = (
    (.0001, 7617), (.0002, 5283), (.0003, 4284), (.0004, 3692),
    (.0005, 3280), (.0006, 2968), (.0007, 2724), (.0008, 2529),
    (.0009, 2369), (.001, 2237), (.002, 1511), (.003, 1202),
    (.004, 993), (.005, 859), (.006, 760), (.007, 683), (.008, 621),
    (.009, 569), (.010, 526), (.015, 383), (.020, 304), (.025, 257),
    (.030, 221), (.035, 190), (.040, 171), (.045, 149), (.050, 132),
    (.060, 109), (.070, 92.7), (.080, 82.9), (.090, 74.3), (.10, 51.4),
    (.12, 32.9), (.14, 20.2), (.15, 16.7), (.16, 13.1), (.18, 8.72),
    (.20, 6.01), (.22, 4.69), (.25, 4.38), (.28, 4.40), (.30, 4.40),
    (.35, 4.12), (.40, 3.46), (.45, 2.75), (.50, 2.15), (.60, 1.25),
    (.70, .722), (.80, .416), (.90, .245), (1.0, .147), (1.2, .060),
    (1.5, .020), (2.0, .0043), (2.25, .0020), (2.5, .0019),
    (2.75, .0017), (3.0, .0010), (3.5, .0018), (4.0, .0092),
    (4.5, .0290), (5.0, .0514), (5.5, .0493), (6.0, .0317),
    (6.5, .0162), (7.0, .0088), (7.5, .0066), (8.0, .0099),
    (8.5, .0143), (9.0, .0159), (9.5, .0144), (10.0, .0120),
    (10.5, .0142), (11.0, .0227), (11.5, .0252), (12.0, .0206),
    (12.5, .0128), (13.0, .0066), (13.5, .0041), (14.0, .0035),
    (14.5, .0031), (15.0, .0030),
)

# Square brackets are deduced, plain values at 275--1000 Td are recommended,
# and parentheses above 1000 Td are suggested.  10^-17 V cm2 equals one Td;
# 10^6 cm/s equals 10^4 m/s.
TABLE36 = (
    (0, 0.0, "deduced"), (25, 4.1, "deduced"),
    (50, 6.8, "deduced"), (100, 10.2, "deduced"),
    (150, 12.1, "deduced"), (200, 13.5, "deduced"),
    (250, 15.6, "deduced"), (275, 17.0, "recommended"),
    (300, 18.3, "recommended"), (350, 20.5, "recommended"),
    (400, 22.6, "recommended"), (450, 24.6, "recommended"),
    (500, 26.4, "recommended"), (550, 28.2, "recommended"),
    (600, 30.0, "recommended"), (650, 31.7, "recommended"),
    (700, 33.4, "recommended"), (750, 35.1, "recommended"),
    (800, 36.8, "recommended"), (850, 38.5, "recommended"),
    (900, 40.1, "recommended"), (950, 41.8, "recommended"),
    (1000, 43.4, "recommended"), (1500, 58.3, "suggested"),
    (2000, 71.7, "suggested"), (2500, 83.9, "suggested"),
    (3000, 95.4, "suggested"), (3500, 106.3, "suggested"),
    (4000, 116.8, "suggested"),
)

TABLE35 = (
    (200, -55.3), (250, -32.8), (300, -16.1), (350, -2.43),
    (400, 10.9), (450, 25.8), (500, 39.3), (550, 51.9),
    (600, 63.8), (650, 75.2), (700, 87.0), (750, 98.8),
    (800, 110), (850, 122), (900, 132), (950, 143),
    (1000, 154), (1250, 204), (1500, 250), (2000, 338),
    (2500, 413), (3000, 478), (3500, 531), (4000, 578),
)

# Table 37: 10^-12 cm3/s = 10^-18 m3/s.
TABLE37 = (
    (100, 887), (150, 702), (200, 632), (250, 626),
    (300, 637), (350, 648), (400, 440), (450, 647),
    (500, 647), (550, 637), (600, 621), (650, 590),
)

TABLES = {
    "table9_total_scattering.csv": (
        "recommended_total_scattering", TABLE9),
    "table14_momentum_transfer.csv": (
        "suggested_elastic_momentum_transfer", TABLE14),
    "table15_vibrational_excitation.csv": (
        "deduced_total_vibrational_excitation", TABLE15),
    "table17_total_ionization.csv": (
        "recommended_total_ionization", TABLE17),
    "table20_total_neutral_dissociation.csv": (
        "deduced_total_neutral_dissociation_requires_confirmation", TABLE20),
    "table28_total_attachment.csv": (
        "suggested_room_temperature_total_attachment", TABLE28),
}


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _csv_text(fields, rows) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(fields)
    writer.writerows(rows)
    return stream.getvalue()


def _payloads() -> dict[str, str]:
    payloads = {}
    for name, (evidence, rows) in TABLES.items():
        payloads[name] = _csv_text(
            ("electron_energy_eV", "cross_section_m2", "evidence_class"),
            ((energy, sigma * 1.0e-20, evidence) for energy, sigma in rows),
        )
    payloads["table36_drift_velocity.csv"] = _csv_text(
        (
            "reduced_electric_field_Td", "drift_velocity_m_s",
            "recommendation_class", "gas_temperature_K",
        ),
        (
            (field, velocity * 1.0e4, grade, "293-300")
            for field, velocity, grade in TABLE36
        ),
    )
    payloads["table35_effective_ionization.csv"] = _csv_text(
        (
            "reduced_electric_field_Td",
            "effective_ionization_coefficient_m2",
            "evidence_class",
        ),
        (
            (field, coefficient * 1.0e-22, "recommended_swarm_fit")
            for field, coefficient in TABLE35
        ),
    )
    payloads["table37_attachment_rate.csv"] = _csv_text(
        (
            "reduced_electric_field_Td",
            "attachment_rate_coefficient_m3_s",
            "evidence_class",
        ),
        (
            (field, rate * 1.0e-18, "assessed_from_eta_over_n_times_drift")
            for field, rate in TABLE37
        ),
    )
    return payloads


def _manifest(payloads: dict[str, str]) -> str:
    value = {
        "schema": "petch.christophorou_olthoff_2000_sf6.v1",
        "source": {
            "title": "Electron Interactions with SF6",
            "authors": ["L. G. Christophorou", "J. K. Olthoff"],
            "journal": "J. Phys. Chem. Ref. Data 29, 267-330 (2000)",
            "doi": "10.1063/1.1288407",
            "nist_url": "https://srd.nist.gov/jpcrdreprint/1.1288407.pdf",
            "pdf_sha256": PDF_SHA256,
            "source_pdf_committed": False,
        },
        "transcription": {
            "method": "manual_double_column_transcription_visually_checked_at_400_dpi",
            "si_conversion": "cross sections 1e-20 m2; drift 1e6 cm/s to 1e4 m/s",
            "render_sha256": TABLE_PAGE_RENDER_SHA256,
            "render_pixels_committed": False,
        },
        "tables": {
            name: {
                "sha256": sha256(text.encode()).hexdigest(),
                "row_count": len(text.splitlines()) - 1,
            }
            for name, text in payloads.items()
        },
        "use_boundaries": {
            "total_scattering_is_not_elastic_momentum_transfer": True,
            "table15_is_deduced_not_directly_measured": True,
            "table20_is_approximate_and_requires_confirmation": True,
            "table28_is_room_temperature": True,
            "table36_direct_recommended_interval_Td": [275, 1000],
            "table35_is_recommended_effective_townsend_fit": True,
            "table37_is_derived_product_not_independent_attachment_data": True,
            "supports_collision_input": True,
            "supports_unique_reactor_state": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
        },
    }
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--render-directory", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.source_pdf is not None and _sha(args.source_pdf) != PDF_SHA256:
        raise RuntimeError("NIST SF6 source PDF checksum changed")
    if args.render_directory is not None:
        for label, expected in TABLE_PAGE_RENDER_SHA256.items():
            page = label.split("page", 1)[1]
            path = args.render_directory / f"page_{page}.png"
            if _sha(path) != expected:
                raise RuntimeError(f"NIST SF6 render checksum changed: {path}")

    payloads = _payloads()
    payloads["extraction_manifest.json"] = _manifest(payloads)
    if args.check:
        for name, expected in payloads.items():
            path = OUTPUT / name
            if path.read_text(encoding="utf-8") != expected:
                raise RuntimeError(f"committed SF6 transcription changed: {path}")
        print("NIST SF6 committed transcriptions match generator")
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, text in payloads.items():
        (OUTPUT / name).write_text(text, encoding="utf-8")
    print(json.dumps({
        name: sha256(text.encode()).hexdigest()
        for name, text in payloads.items()
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
