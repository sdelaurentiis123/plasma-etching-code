#!/usr/bin/env python3
"""Materialize the NIST-evaluated CFx ionization tables.

The rows are a manual transcription of Tables 31--33 in Christophorou,
Olthoff, and Rao, J. Phys. Chem. Ref. Data 25, 1341 (1996),
doi:10.1063/1.555986.  The official NIST scan and 240-dpi renders of journal
pages 1383--1384 were visually checked at original resolution.  When local
copies are supplied, both source and renders are checksum gated.
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
    / "christophorou_olthoff_rao_1996_cfx"
)
PDF_SHA256 = (
    "381e368840e28c84bb03eb4e684691e5c330c8c4f592bc666ae13bc292a39244"
)
RENDER_SHA256 = {
    "page-43.png": "c82b0fe08649cc53365ce560f5c9c802a943f837f27b1d2627dd021cc285634f",
    "page-44.png": "40e36245eaaa74bdb3aab785458e43f567f8dd9f8f553815453cf3d1643f62fe",
}

ENERGY = (
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28,
    30, 32, 34, 36, 38, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100,
    120, 140, 160, 180, 200,
)

# Table 31, units 1e-20 m2. None is a blank source cell, not measured zero.
TABLE31 = {
    "CF3": (
        .015, .029, .041, .060, .099, .111, .145, .157, .167, .194,
        .204, .270, .303, .315, .320, .325, .329, .335, .338, .346,
        .350, .358, .360, .372, .374, .380, .376, .368, .365, .350,
        .342, .333, .318, .306, .292,
    ),
    "CF2": (
        .05, .09, .15, .18, .26, .35, .39, .42, .47, .55,
        .64, .69, .73, .78, .82, .87, .89, .91, .93, .96,
        .98, .99, 1.01, 1.03, 1.03, 1.05, 1.03, .99, .96, .91,
        .86, .78, .67, .58, .49,
    ),
    "CF": (
        None, None, .03, .07, .13, .18, .23, .28, .33, .40,
        .45, .55, .63, .70, .76, .81, .86, .91, .95, .99,
        1.01, 1.08, 1.15, 1.18, 1.23, 1.25, 1.25, 1.26, 1.25,
        1.23, 1.14, 1.04, .90, .79, .67,
    ),
}

# Table 32. The CF3+ parent column duplicates Table 31 at lower precision and
# is therefore kept as a source check, not a second curve.
TABLE32_CF2 = {
    18: .06, 19: .12, 20: .17, 22: .25, 24: .31, 26: .34,
    28: .40, 30: .49, 32: .53, 34: .56, 36: .59, 38: .61,
    40: .63, 45: .65, 50: .67, 55: .71, 60: .72, 65: .74,
    70: .76, 80: .79, 90: .78, 100: .78, 120: .78, 140: .77,
    160: .76, 180: .74, 200: .73,
}
TABLE32_CF = {
    22: .04, 24: .10, 26: .15, 28: .20, 30: .26, 32: .31,
    34: .34, 36: .36, 38: .37, 40: .40, 45: .45, 50: .53,
    55: .58, 60: .62, 65: .65, 70: .68, 80: .70, 90: .72,
    100: .73, 120: .75, 140: .77, 160: .76, 180: .74, 200: .72,
}

# Table 33. The CF2+ parent column duplicates Table 31 at lower precision.
TABLE33_CF = {
    15: .04, 16: .09, 17: .13, 18: .18, 19: .20, 20: .23,
    22: .31, 24: .36, 26: .40, 28: .43, 30: .48, 32: .62,
    34: .74, 36: .80, 38: .84, 40: .88, 45: .97, 50: 1.02,
    55: 1.08, 60: 1.11, 65: 1.16, 70: 1.19, 80: 1.22,
    90: 1.25, 100: 1.28, 120: 1.24, 140: 1.18, 160: 1.12,
    180: 1.05, 200: .93,
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
    parent_rows = []
    uncertainty = {"CF": .15, "CF2": .16, "CF3": .18}
    for target, values in TABLE31.items():
        if len(values) != len(ENERGY):
            raise RuntimeError(f"Table 31 length changed for {target}")
        parent_rows.extend(
            (
                target, f"{target}+", energy, sigma * 1.0e-20,
                uncertainty[target], "measured_parent_ionization",
            )
            for energy, sigma in zip(ENERGY, values)
            if sigma is not None
        )
    payloads = {
        "table31_parent_ionization.csv": _csv_text(
            (
                "target_neutral", "product_ion", "electron_energy_eV",
                "cross_section_m2", "relative_uncertainty", "evidence_class",
            ),
            parent_rows,
        ),
        "table32_cf3_dissociative_ionization.csv": _csv_text(
            (
                "target_neutral", "product_ion", "neutral_coproduct",
                "electron_energy_eV", "cross_section_m2",
                "relative_uncertainty", "evidence_class",
            ),
            (
                *(
                    ("CF3", "CF2+", "F", energy, sigma * 1.0e-20, .20,
                     "measured_dissociative_ionization")
                    for energy, sigma in TABLE32_CF2.items()
                ),
                *(
                    ("CF3", "CF+", "2F", energy, sigma * 1.0e-20, .20,
                     "measured_dissociative_ionization")
                    for energy, sigma in TABLE32_CF.items()
                ),
                ("CF3", "F+", "CF2", 70, .35e-20, .30,
                 "single_energy_anchor_not_curve"),
            ),
        ),
        "table33_cf2_dissociative_ionization.csv": _csv_text(
            (
                "target_neutral", "product_ion", "neutral_coproduct",
                "electron_energy_eV", "cross_section_m2",
                "relative_uncertainty", "evidence_class",
            ),
            (
                *(
                    ("CF2", "CF+", "F", energy, sigma * 1.0e-20, .16,
                     "measured_net_CFplus_including_two_onsets")
                    for energy, sigma in TABLE33_CF.items()
                ),
                ("CF2", "F+", "CF", 70, .60e-20, .30,
                 "single_energy_anchor_mixed_single_and_ion_pair"),
            ),
        ),
    }
    return payloads


def _manifest(payloads: dict[str, str]) -> str:
    value = {
        "schema": "petch.christophorou_olthoff_rao_1996_cfx.v1",
        "source": {
            "title": "Electron Interactions with CF4",
            "authors": [
                "L. G. Christophorou", "J. K. Olthoff", "M. V. V. S. Rao",
            ],
            "journal": "J. Phys. Chem. Ref. Data 25, 1341-1388 (1996)",
            "doi": "10.1063/1.555986",
            "nist_url": "https://www.nist.gov/system/files/documents/srd/jpcrd512.pdf",
            "pdf_sha256": PDF_SHA256,
            "source_pdf_committed": False,
        },
        "transcription": {
            "method": "manual_table_transcription_visually_checked_at_240_dpi",
            "source_journal_pages": [1383, 1384],
            "si_conversion": "printed cross sections 1e-20 m2",
            "render_sha256": RENDER_SHA256,
            "render_pixels_committed": False,
        },
        "tables": {
            name: {
                "sha256": sha256(text.encode()).hexdigest(),
                "row_count": len(text.splitlines()) - 1,
            }
            for name, text in payloads.items()
        },
        "source_checks": {
            "table31_relative_uncertainty": {"CF": .15, "CF2": .16, "CF3": .18},
            "table32_CF2plus_and_CFplus_relative_uncertainty": .20,
            "table33_CFplus_relative_uncertainty": .16,
            "table32_Fplus_70eV_relative_uncertainty": .30,
            "table33_Fplus_70eV_relative_uncertainty": .30,
        },
        "use_boundaries": {
            "supports_secondary_CFx_electron_ionization": True,
            "supports_secondary_neutral_dissociation": False,
            "table32_Fplus_is_one_energy_anchor": True,
            "table33_Fplus_is_one_energy_mixed_channel_anchor": True,
            "table33_CFplus_combines_single_and_ion_pair_onsets": True,
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
        raise RuntimeError("NIST CF4 source PDF checksum changed")
    if args.render_directory is not None:
        for name, expected in RENDER_SHA256.items():
            path = args.render_directory / name
            if _sha(path) != expected:
                raise RuntimeError(f"NIST CFx render checksum changed: {path}")

    payloads = _payloads()
    payloads["extraction_manifest.json"] = _manifest(payloads)
    if args.check:
        for name, expected in payloads.items():
            path = OUTPUT / name
            if path.read_text(encoding="utf-8") != expected:
                raise RuntimeError(f"committed CFx transcription changed: {path}")
        print("NIST CFx committed transcriptions match generator")
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
