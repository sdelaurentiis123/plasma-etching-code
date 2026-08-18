#!/usr/bin/env python3
"""Materialize the visually checked Lan--Jeon C4F6 collision tables.

Tables 1 and 2 are complete numerical source tables in Lan and Jeon,
J. Korean Phys. Soc. 64, 1320--1326 (2014), doi:10.3938/jkps.64.1320.
The source PDF is not redistributed. When supplied, its bytes and the 600 dpi
renders used for visual transcription are checksum gated.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "experimental" / "lan_jeon_2014_c4f6"
PDF_SHA256 = "82d672f5b60611894a4584aa503e4c26d66aec3c66792550de46fb660f77aeb6"
RENDER_SHA256 = {
    "table1_pdf_page3_600dpi": "ccf3f0d9da9699d7074cd36b229fc6225405836599c9dc618085dc7ee9cc3bbc",
    "table2_pdf_page5_600dpi": "59fce2e88c65d14c0fc165267a17197d79657cb7f7ee702d8dc4e4857231abb0",
    "figure7_pdf_page6_600dpi": "6cf27de7bb64e8cf9a34969b7ed9d44f5508bac0c6705eec3ef34f52fb540033",
}

# Source unit: 1e-16 cm2 = 1e-20 m2.
TABLE1_MOMENTUM = (
    (.0001, 83.9), (.13, 83.5), (.2, 80.0), (.25, 77.3),
    (.32, 72.8), (.37, 71.0), (.42, 66.7), (.5, 65.5),
    (.53, 63.9), (.55, 63.8), (.57, 59.8), (.577, 56.0),
    (.65, 2.3), (.7, .3), (.73, .095), (.74, .0915),
    (.75, .09), (.86, .09), (1.08, .09), (1.12, .092),
    (1.14, .1), (1.44, 20.1), (1.55, 24.1), (1.73, 23.8),
    (1.95, 21.8), (2.09, 21.1), (2.5, 19.1), (2.7, 17.25),
    (2.9, 15.6), (3.33, 12.89), (5.3, 7.3), (6.0, 7.6),
    (7.0, 8.0), (8.84, 8.6), (14.5, 8.4), (46.9, 4.3),
    (362.0, 1.1),
)

TABLE2 = {
    "Qa": (
        (.0001, 90.0), (.04, 10.0), (.08, 3.38), (.1, 2.94),
        (.14, 3.08), (.18, 3.15), (.22, 2.66), (.25, 2.10),
        (.35, .77), (.46, .34), (.5, .26), (.6, .19), (.65, .17),
        (.75, .17), (.9, .18), (1.0, .2), (1.1, .19), (1.25, .1),
        (1.4, .011), (1.45, .004), (1.6, .0034), (1.7, .003),
        (2.0, .002), (2.28, .0017), (2.55, .0014), (2.75, .0017),
        (3.0, .001), (3.65, .0007), (4.15, .0006), (5.1, .00044),
        (5.7, .0003), (6.0, .00003), (6.5, .00004), (7.0, 0.0),
    ),
    "Qv1": (
        (.25, 0.0), (.35, .16), (.46, 10.50), (.5, 26.6),
        (.6, 28.16), (.65, 26.6), (.75, 17.9), (.9, 11.55),
        (1.0, 8.2), (1.1, 6.35), (1.25, 4.4), (1.4, 3.3),
        (1.45, 3.0), (1.6, 2.35), (1.7, 2.06), (2.0, .7),
        (2.28, .22), (2.55, .07), (2.75, .037), (3.0, .0125),
        (3.65, .0025), (4.15, .00075), (5.1, .00036), (5.7, .0001),
        (6.0, .00008), (6.5, .00005), (7.0, .00003),
    ),
    "Qv2": (
        (.35, 0.0), (.46, .01), (.5, .015), (.6, .08), (.65, .234),
        (.75, .92), (.9, 8.2), (1.0, 18.25), (1.1, 16.9),
        (1.25, 15.0), (1.4, 7.35), (1.45, 6.7), (1.6, 5.5),
        (1.7, 4.4), (2.0, 2.35), (2.28, 2.0), (2.55, 2.06),
        (2.75, 4.23), (3.0, 7.76), (3.65, 10.13), (4.15, 11.12),
        (5.1, 12.5), (5.7, 12.75), (6.0, 12.49), (6.5, 12.25),
        (7.0, 12.0), (7.7, 10.6), (10.0, 3.6), (12.0, .9),
        (14.0, .28), (17.0, .05), (19.0, .02), (22.0, .006),
        (35.0, .00029), (60.0, .00001),
    ),
    "Qex1": (
        (2.0, 0.0), (3.0, .00008), (3.65, .0025), (4.15, .025),
        (5.1, .86), (5.7, 2.64), (6.0, 2.9), (6.5, 2.82),
        (7.0, 2.74), (7.7, 2.4), (10.0, .1), (12.0, .011),
        (14.0, .0016),
    ),
    "Qex2": (
        (3.0, 0.0), (6.0, .00151), (7.0, .026), (7.7, .2),
        (10.0, 1.59), (12.0, 2.14), (14.0, 2.14), (17.0, 2.0),
        (19.0, 1.8), (22.0, 1.4), (35.0, .446), (60.0, .096),
        (100.0, .015), (110.0, .01),
    ),
    "Qex3": (
        (6.5, 0.0), (7.0, .0035), (7.7, .0038), (10.0, 1.09),
        (12.0, 10.32), (14.0, 12.5), (17.0, 11.8), (19.0, 10.35),
        (22.0, 6.3), (35.0, 1.4), (60.0, .23), (100.0, .05),
    ),
    "Qex4": (
        (7.7, .00001), (10.0, .6), (12.0, 1.52), (14.0, 2.8),
        (17.0, 3.1), (19.0, 3.0), (22.0, 2.6), (35.0, 1.3),
        (60.0, .45), (100.0, .17), (110.0, .14), (200.0, .04),
    ),
    "Qdiss": (
        (10.0, .0007), (12.0, .004), (14.0, .01), (17.0, .02),
        (19.0, .02), (22.0, .037), (35.0, .08), (60.0, .13),
        (100.0, .17),
    ),
    "Qi": (
        (12.0, 0.0), (14.0, .6), (17.0, 1.3), (19.0, 1.8),
        (22.0, 2.6), (35.0, 5.6), (60.0, 9.0), (100.0, 10.7),
        (110.0, 10.5), (200.0, 9.4),
    ),
}

PROCESS_METADATA = {
    "Qa": ("ATTACHMENT", 0.0, "measured_attachment_imported_unchanged"),
    "Qv1": ("EXCITATION", .25, "swarm_regressed_effective_vibration"),
    "Qv2": ("EXCITATION", .35, "cC4F8_analog_then_swarm_regressed"),
    "Qex1": ("EXCITATION", 2.0, "swarm_regressed_effective_excitation"),
    "Qex2": ("EXCITATION", 3.0, "swarm_regressed_effective_excitation"),
    "Qex3": ("EXCITATION", 6.5, "swarm_regressed_effective_excitation"),
    "Qex4": ("EXCITATION", 7.7, "cC4F8_analog_then_swarm_regressed"),
    "Qdiss": ("EXCITATION", 10.0, "cC4F8_analog_imported_unchanged"),
    "Qi": ("IONIZATION", 12.0, "measured_total_ionization_imported_unchanged"),
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
    payloads = {
        "table1_momentum_transfer.csv": _csv_text(
            ("electron_energy_eV", "cross_section_m2", "evidence_class"),
            (
                (energy, sigma * 1.0e-20, "swarm_regressed_from_measured_drift")
                for energy, sigma in TABLE1_MOMENTUM
            ),
        )
    }
    rows = []
    for label, points in TABLE2.items():
        kind, threshold, evidence = PROCESS_METADATA[label]
        rows.extend(
            (label, kind, threshold, energy, sigma * 1.0e-20, evidence)
            for energy, sigma in points
        )
    payloads["table2_inelastic.csv"] = _csv_text(
        (
            "process_label", "kind", "energy_loss_eV", "electron_energy_eV",
            "cross_section_m2", "evidence_class",
        ),
        rows,
    )
    manifest = {
        "schema": "petch.lan-jeon-2014-c4f6-electron-evidence.v1",
        "citation": {
            "authors": "Phan-Thi Lan and Byong-Hoon Jeon",
            "title": "Determination of the Electron Collision Cross-section Set for the C4F6 Molecule by Using an Electron Swarm Study",
            "journal": "Journal of the Korean Physical Society 64, 1320-1326 (2014)",
            "doi": "10.3938/jkps.64.1320",
        },
        "source_pdf_sha256": PDF_SHA256,
        "visual_audit_render_sha256": RENDER_SHA256,
        "source_units": "1e-16 cm2",
        "si_conversion": "1e-16 cm2 = 1e-20 m2",
        "table1_point_count": len(TABLE1_MOMENTUM),
        "table2_point_count_by_process": {
            label: len(points) for label, points in TABLE2.items()
        },
        "evidence_boundary": {
            "swarm_fit_domain_Td": [0.35, 1200.0],
            "reported_typical_drift_reproduction": "within +/-5%",
            "reported_exception_0p5pct_20_35Td_percent": 21.0,
            "reported_exception_10pct_120_230Td_percent": -16.0,
            "figure8_present_in_source_pdf": False,
            "note": (
                "The text references Figure 8, but the archived seven-page "
                "publisher PDF contains Figures 1-7 only. Figure 7 remains a "
                "combined measured/source-replay board."
            ),
            "supports_independent_validation_of_regressed_rows": False,
            "supports_reactor_state_prediction": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
        },
    }
    manifest["generated_payload_sha256"] = {
        name: sha256(text.encode()).hexdigest() for name, text in payloads.items()
    }
    payloads["manifest.json"] = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    return payloads


def _verify_optional(path: Path | None, expected: str, label: str) -> None:
    if path is not None and _sha(path) != expected:
        raise RuntimeError(f"{label} checksum changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--table1-render", type=Path)
    parser.add_argument("--table2-render", type=Path)
    parser.add_argument("--figure7-render", type=Path)
    args = parser.parse_args()
    _verify_optional(args.source_pdf, PDF_SHA256, "source PDF")
    _verify_optional(
        args.table1_render, RENDER_SHA256["table1_pdf_page3_600dpi"],
        "Table 1 render",
    )
    _verify_optional(
        args.table2_render, RENDER_SHA256["table2_pdf_page5_600dpi"],
        "Table 2 render",
    )
    _verify_optional(
        args.figure7_render, RENDER_SHA256["figure7_pdf_page6_600dpi"],
        "Figure 7 render",
    )
    payloads = _payloads()
    if args.check:
        for name, text in payloads.items():
            path = OUTPUT / name
            if not path.exists() or path.read_text(encoding="utf-8") != text:
                raise SystemExit(f"committed Lan--Jeon payload is stale: {path}")
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, text in payloads.items():
        (OUTPUT / name).write_text(text, encoding="utf-8")
    print(f"wrote {len(payloads)} Lan--Jeon evidence payloads to {OUTPUT}")


if __name__ == "__main__":
    main()
