#!/usr/bin/env python3
"""Audit NIST Cl2-attachment energy support on the Lam Te board.

The audit evaluates only the printed 0.05--11.8 eV cross-section support. It
does not extrapolate either tail, tune a coefficient, infer absorbed power,
or use an etched depth. The 1e-6 kernel tolerance is inherited from the
pre-existing evaluated-cross-section provider rather than selected from the
resulting board.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from statistics import median

from petch.reactor_global import (
    MalyshevMeasuredElectronTemperatureProvider,
    RateContext,
    nist_cl2_dissociative_attachment_cross_section_support,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "results" / "curated" / "reactor_global_chlorine"
CSV_PATH = OUTPUT_DIRECTORY / "malyshev_1998_attachment_energy_support.csv"
AUDIT_PATH = OUTPUT_DIRECTORY / "malyshev_1998_attachment_energy_support.json"
REPORT_PATH = OUTPUT_DIRECTORY.joinpath(
    "MALYSHEV_1998_ATTACHMENT_ENERGY_SUPPORT.md")
SOURCE_TABLE = (
    ROOT / "research_sources" / "digitized"
    / "christophorou_olthoff_1999_table16_cl2_attachment.csv"
)
KERNEL_TOLERANCE = 1.0e-6

FIELDNAMES = (
    "window_to_wafer_gap_cm",
    "pressure_mTorr",
    "tcp_source_power_W",
    "electron_temperature_eV",
    "attachment_rate_support_m3_s",
    "attachment_incident_energy_support_eV_m3_s",
    "support_conditioned_mean_incident_energy_eV",
    "rate_kernel_missing_below_0p05_fraction",
    "rate_kernel_missing_above_11p8_fraction",
    "incident_energy_kernel_missing_below_0p05_fraction",
    "incident_energy_kernel_missing_above_11p8_fraction",
    "rate_support_complete_at_inherited_tolerance",
    "incident_energy_support_complete_at_inherited_tolerance",
    "supports_absorbed_power_solve",
    "supports_wafer_flux",
    "supports_feature_depth",
)


def build_outputs() -> tuple[str, dict[str, object], str]:
    temperature_provider = (
        MalyshevMeasuredElectronTemperatureProvider.from_package_data())
    attachment = (
        nist_cl2_dissociative_attachment_cross_section_support())
    rows: list[dict[str, str]] = []
    numeric_rows: list[dict[str, float | bool]] = []

    for marker in temperature_provider.markers:
        temperature = marker.electron_temperature_eV
        context = RateContext(temperature)
        rate = attachment.tabulated_rate_coefficient_si(context)
        energy = (
            attachment.tabulated_incident_energy_moment_eV_m3_s(context))
        rate_low, rate_high = attachment.rate_kernel_missing_fractions(
            temperature)
        energy_low, energy_high = (
            attachment.incident_energy_kernel_missing_fractions(temperature))
        rate_complete = max(rate_low, rate_high) <= KERNEL_TOLERANCE
        energy_complete = max(energy_low, energy_high) <= KERNEL_TOLERANCE
        numeric = {
            "temperature": temperature,
            "rate": rate,
            "energy": energy,
            "mean_energy": energy / rate,
            "rate_low": rate_low,
            "rate_high": rate_high,
            "energy_low": energy_low,
            "energy_high": energy_high,
            "rate_complete": rate_complete,
            "energy_complete": energy_complete,
        }
        numeric_rows.append(numeric)
        rows.append({
            "window_to_wafer_gap_cm": (
                f"{marker.window_to_wafer_gap_cm:g}"),
            "pressure_mTorr": f"{marker.pressure_mTorr:g}",
            "tcp_source_power_W": f"{marker.tcp_source_power_W:.3f}",
            "electron_temperature_eV": f"{temperature:.8g}",
            "attachment_rate_support_m3_s": f"{rate:.12e}",
            "attachment_incident_energy_support_eV_m3_s": (
                f"{energy:.12e}"),
            "support_conditioned_mean_incident_energy_eV": (
                f"{energy / rate:.12e}"),
            "rate_kernel_missing_below_0p05_fraction": f"{rate_low:.12e}",
            "rate_kernel_missing_above_11p8_fraction": f"{rate_high:.12e}",
            "incident_energy_kernel_missing_below_0p05_fraction": (
                f"{energy_low:.12e}"),
            "incident_energy_kernel_missing_above_11p8_fraction": (
                f"{energy_high:.12e}"),
            "rate_support_complete_at_inherited_tolerance": (
                str(rate_complete).lower()),
            "incident_energy_support_complete_at_inherited_tolerance": (
                str(energy_complete).lower()),
            "supports_absorbed_power_solve": "false",
            "supports_wafer_flux": "false",
            "supports_feature_depth": "false",
        })

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    csv_payload = stream.getvalue()

    def bounds(key: str) -> dict[str, float]:
        values = [float(row[key]) for row in numeric_rows]
        return {
            "minimum": min(values),
            "median": median(values),
            "maximum": max(values),
        }

    rate_complete_count = sum(
        bool(row["rate_complete"]) for row in numeric_rows)
    energy_complete_count = sum(
        bool(row["energy_complete"]) for row in numeric_rows)
    audit: dict[str, object] = {
        "audit_id": "MALYSHEV-1998-CL2-ATTACHMENT-ENERGY-SUPPORT-R1",
        "claim_class": (
            "measured-Te support diagnostic; not an absorbed-power, reactor-"
            "state, wafer-flux, or feature-depth prediction"
        ),
        "coefficient_selection_target": None,
        "reactor_fit_target": None,
        "feature_depth_target": None,
        "source_table": {
            "identity": "Christophorou--Olthoff 1999 Table 16",
            "energy_support_eV": [
                attachment.electron_energy_eV[0],
                attachment.electron_energy_eV[-1],
            ],
            "row_count": len(attachment.electron_energy_eV),
            "csv_sha256": hashlib.sha256(
                SOURCE_TABLE.read_bytes()).hexdigest(),
            "scalar_relative_uncertainty": None,
        },
        "kernel_tolerance": {
            "value": KERNEL_TOLERANCE,
            "basis": (
                "inherited from the pre-existing evaluated Maxwellian cross-"
                "section provider before this board was evaluated"
            ),
            "interpretation": (
                "constant-cross-section EEDF-kernel exposure diagnostic; not "
                "a bound on the unknown cross-section-weighted tail error"
            ),
        },
        "marker_accounting": {
            "lam_measured_Te_markers": len(numeric_rows),
            "particle_rate_support_complete": rate_complete_count,
            "incident_energy_support_complete": energy_complete_count,
        },
        "electron_temperature_eV": bounds("temperature"),
        "support_only_attachment_rate_m3_s": bounds("rate"),
        "support_only_incident_energy_moment_eV_m3_s": bounds("energy"),
        "support_conditioned_mean_incident_energy_eV": bounds("mean_energy"),
        "rate_kernel_missing_below_table": bounds("rate_low"),
        "rate_kernel_missing_above_table": bounds("rate_high"),
        "incident_energy_kernel_missing_below_table": bounds("energy_low"),
        "incident_energy_kernel_missing_above_table": bounds("energy_high"),
        "electron_power_closure": {
            "status": "blocked",
            "demonstrated_blocker": (
                "Table 16 does not close either Maxwellian tail on any of "
                "the 62 measured Lam electron-temperature markers at the "
                "inherited tolerance"
            ),
            "other_unclosed_channels": [
                "Cl2 elastic momentum-transfer energy exchange",
                "Cl2 vibrational excitation",
                "Cl2 non-dissociative electronic excitation",
                "electron-impact detachment from Cl-",
                "ion-pair energy partition",
                "species-resolved molecular ionization branch",
            ],
            "next_required_evidence": (
                "primary attachment cross sections or energy-weighted rates "
                "covering the Lam EEDF, plus the remaining channel moments"
            ),
        },
        "csv_sha256": hashlib.sha256(
            csv_payload.encode("utf-8")).hexdigest(),
    }

    temperatures = audit["electron_temperature_eV"]
    rate_high = audit["rate_kernel_missing_above_table"]
    energy_high = audit["incident_energy_kernel_missing_above_table"]
    mean_energy = audit["support_conditioned_mean_incident_energy_eV"]
    report = f"""# Lam chlorine attachment-energy support audit

**Claim class: measured-Te support diagnostic; not a power or depth result**

The 500-dpi-audited NIST Table 16 supplies 42 Cl2 dissociative-attachment
cross sections over `0.05--11.8 eV`. All 62 Lam Figure-3 electron-temperature
markers (`{temperatures['minimum']:.4f}--{temperatures['maximum']:.4f} eV`)
were evaluated without extrapolating that table.

At the inherited `1e-6` support tolerance, particle-rate support is complete
for `{rate_complete_count}/62` markers and incident-energy support is complete
for `{energy_complete_count}/62`. The missing constant-cross-section kernel
above 11.8 eV spans `{rate_high['minimum']:.6g}--{rate_high['maximum']:.6g}`
for `<sigma v>` and
`{energy_high['minimum']:.6g}--{energy_high['maximum']:.6g}`
for `<sigma v E>`. These are EEDF-kernel exposure diagnostics, not bounds on
the unknown cross-section-weighted error.

On printed support alone, the collision-conditioned incident energy spans
`{mean_energy['minimum']:.6g}--{mean_energy['maximum']:.6g} eV`. That number is
reported as a partial moment, not substituted for a complete attachment
energy loss.

## Verdict

Table 16 advances the chlorine ledger because it separates the attachment
particle moment from the electron-removal energy moment. It does **not** close
the Lam electron-power balance: both tails remain exposed on every measured
temperature marker, and elastic, vibrational, non-dissociative electronic,
detachment, ion-pair, and molecular branching channels remain open.

No coefficient was selected against a reactor observable or feature depth.
No absorbed-power, wafer-flux, or etched-depth prediction is supported by this
audit.
"""
    return csv_payload, audit, report


def main() -> None:
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    csv_payload, audit, report = build_outputs()
    json_payload = json.dumps(audit, indent=2) + "\n"
    if args.check:
        if CSV_PATH.read_text(encoding="utf-8") != csv_payload:
            raise SystemExit("attachment-energy CSV is stale")
        if AUDIT_PATH.read_text(encoding="utf-8") != json_payload:
            raise SystemExit("attachment-energy JSON is stale")
        if REPORT_PATH.read_text(encoding="utf-8") != report:
            raise SystemExit("attachment-energy report is stale")
        return
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    CSV_PATH.write_text(csv_payload, encoding="utf-8")
    AUDIT_PATH.write_text(json_payload, encoding="utf-8")
    REPORT_PATH.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
