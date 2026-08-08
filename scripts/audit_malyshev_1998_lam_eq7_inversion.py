#!/usr/bin/env python3
"""Freeze the measured-state Malyshev Eq.-7 wall-return inversion board."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import median

from petch.reactor_global import (
    MalyshevMeasuredChlorineDissociationProvider,
    MalyshevMeasuredElectronTemperatureProvider,
    malyshev_1998_eq7_wall_return_inversion,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "results" / "curated" / "reactor_global_chlorine"
CSV_PATH = OUTPUT_DIRECTORY / "malyshev_1998_lam_eq7_inversion.csv"
AUDIT_PATH = OUTPUT_DIRECTORY / "malyshev_1998_lam_eq7_inversion.json"

FIELDNAMES = (
    "source_figure",
    "window_to_wafer_gap_cm",
    "pressure_mTorr",
    "tcp_source_power_W",
    "relative_cl2_density_percent",
    "cl2_dissociation_percent",
    "electron_temperature_eV",
    "electron_temperature_method",
    "volume_average_electron_density_m3",
    "electron_density_method",
    "hamilton_neutral_dissociation_rate_m3_s",
    "lee_dissociative_attachment_rate_m3_s",
    "attachment_fraction_of_cl2_destruction_rate",
    "electron_driven_cl2_destruction_frequency_s_inv",
    "required_wall_return_frequency_s_inv",
    "cl_to_cl2_number_density_ratio",
    "reported_cl2_uncertainty_lower_frequency_s_inv",
    "reported_cl2_uncertainty_upper_frequency_s_inv",
    "reproduced_relative_cl2_density_percent",
    "supports_wall_probability_inference",
    "supports_wafer_flux",
    "supports_feature_depth",
)


def build_outputs() -> tuple[str, dict[str, object]]:
    provider = (
        MalyshevMeasuredChlorineDissociationProvider.from_package_data())
    rows = []
    excluded = {
        "diagnostic_flow_check": 0,
        "nonphysical_or_zero_derived_dissociation": 0,
        "electron_temperature_support_missing": 0,
        "electron_density_support_missing": 0,
    }
    for marker in provider.markers:
        if marker.validation_role == "diagnostic_flow_check":
            excluded["diagnostic_flow_check"] += 1
            continue
        if not marker.supports_eq7_inversion:
            excluded["nonphysical_or_zero_derived_dissociation"] += 1
            continue
        try:
            state = malyshev_1998_eq7_wall_return_inversion(marker)
        except ValueError as error:
            message = str(error)
            if "Figure-3" in message:
                excluded["electron_temperature_support_missing"] += 1
            elif "Figure-11" in message:
                excluded["electron_density_support_missing"] += 1
            else:
                raise
            continue
        total_rate = (
            state.hamilton_neutral_dissociation_rate_m3_s
            + state.lee_dissociative_attachment_rate_m3_s
        )
        upper = state.reported_cl2_uncertainty_upper_frequency_s_inv
        rows.append({
            "source_figure": marker.source_figure,
            "window_to_wafer_gap_cm": f"{marker.window_to_wafer_gap_cm:g}",
            "pressure_mTorr": f"{marker.pressure_mTorr:g}",
            "tcp_source_power_W": f"{marker.tcp_source_power_W:.3f}",
            "relative_cl2_density_percent": (
                f"{marker.relative_cl2_density_percent:.4f}"),
            "cl2_dissociation_percent": (
                f"{marker.cl2_dissociation_percent:.4f}"),
            "electron_temperature_eV": (
                f"{state.electron_temperature_state.electron_temperature.value:.8g}"
            ),
            "electron_temperature_method": (
                state.electron_temperature_state.method),
            "volume_average_electron_density_m3": (
                f"{state.electron_density_state.volume_average_electron_density.value:.8e}"
            ),
            "electron_density_method": state.electron_density_state.method,
            "hamilton_neutral_dissociation_rate_m3_s": (
                f"{state.hamilton_neutral_dissociation_rate_m3_s:.8e}"),
            "lee_dissociative_attachment_rate_m3_s": (
                f"{state.lee_dissociative_attachment_rate_m3_s:.8e}"),
            "attachment_fraction_of_cl2_destruction_rate": (
                f"{state.lee_dissociative_attachment_rate_m3_s / total_rate:.8f}"
            ),
            "electron_driven_cl2_destruction_frequency_s_inv": (
                f"{state.electron_driven_cl2_destruction_frequency_s_inv:.8e}"
            ),
            "required_wall_return_frequency_s_inv": (
                f"{state.required_wall_return_frequency_s_inv:.8e}"),
            "cl_to_cl2_number_density_ratio": (
                f"{state.cl_to_cl2_number_density_ratio:.8e}"),
            "reported_cl2_uncertainty_lower_frequency_s_inv": (
                f"{state.reported_cl2_uncertainty_lower_frequency_s_inv:.8e}"
            ),
            "reported_cl2_uncertainty_upper_frequency_s_inv": (
                "" if upper is None else f"{upper:.8e}"),
            "reproduced_relative_cl2_density_percent": (
                f"{state.reproduced_relative_cl2_density_percent:.8f}"),
            "supports_wall_probability_inference": "false",
            "supports_wafer_flux": "false",
            "supports_feature_depth": "false",
        })

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    csv_payload = stream.getvalue()

    frequencies = [
        float(row["required_wall_return_frequency_s_inv"])
        for row in rows
    ]
    finite_envelope = [
        value for row, value in zip(rows, frequencies)
        if row["reported_cl2_uncertainty_upper_frequency_s_inv"]
    ]
    method_pairs: dict[str, int] = {}
    for row in rows:
        key = (
            f"Te:{row['electron_temperature_method']}|"
            f"ne:{row['electron_density_method']}"
        )
        method_pairs[key] = method_pairs.get(key, 0) + 1

    temperature_provider = (
        MalyshevMeasuredElectronTemperatureProvider.from_package_data())
    ten_mTorr_large_gap_maximum_Te_eV = max(
        marker.electron_temperature_eV
        for marker in temperature_provider.markers
        if marker.window_to_wafer_gap_cm == 11.0
        and marker.pressure_mTorr == 10.0
    )
    printed_kdis_cm3_s = (
        4.52e-8
        * math.exp(-7.40 / ten_mTorr_large_gap_maximum_Te_eV)
    )
    printed_kd_upper_from_attachment_ratio_cm3_s = (
        printed_kdis_cm3_s * (1.0 + 1.0 / 7.0)
    )
    footnote_kd_cm3_s = 7.0e-9
    audit = {
        "audit_id": "MALYSHEV-1998-LAM-EQ7-MEASURED-STATE-INVERSION-R1",
        "claim_class": (
            "measured-state diagnostic inversion; not a wall fit and not "
            "reactor or feature-depth prediction"
        ),
        "coefficient_selection_target": None,
        "feature_depth_target": None,
        "source_equation": (
            "Malyshev et al. Eq. 7: relative Cl2 = "
            "1/(1 + kd*ne/(2*kr))"
        ),
        "electron_destruction_evidence": {
            "neutral_dissociation": (
                "Hamilton 2018 eight-state Maxwellian rate; no scalar "
                "physical uncertainty published"
            ),
            "dissociative_attachment": (
                "Lee-Lieberman retained published-compilation rate"
            ),
        },
        "marker_accounting": {
            "audited_marker_total": len(provider.markers),
            "successful_inversions": len(rows),
            "excluded": excluded,
            "electron_state_method_pairs": method_pairs,
        },
        "required_wall_return_frequency_s_inv": {
            "minimum_all": min(frequencies),
            "median_all": median(frequencies),
            "maximum_all": max(frequencies),
            "finite_reported_cl2_upper_envelope_count": len(finite_envelope),
            "minimum_finite_envelope_subset": min(finite_envelope),
            "median_finite_envelope_subset": median(finite_envelope),
            "maximum_finite_envelope_subset": max(finite_envelope),
        },
        "source_internal_consistency": {
            "status": "printed_nominal_values_do_not_share_one_condition",
            "figure3_maximum_Te_eV_at_11cm_10mTorr": (
                ten_mTorr_large_gap_maximum_Te_eV),
            "printed_kdis_at_that_Te_cm3_s": printed_kdis_cm3_s,
            "printed_kd_upper_using_source_attachment_1_over_7_cm3_s": (
                printed_kd_upper_from_attachment_ratio_cm3_s),
            "footnote14_kd_cm3_s": footnote_kd_cm3_s,
            "footnote_to_printed_upper_ratio": (
                footnote_kd_cm3_s
                / printed_kd_upper_from_attachment_ratio_cm3_s
            ),
            "footnote14_ne_cm3": 1.0e11,
            "footnote14_destruction_frequency_s_inv": 700.0,
            "printed_upper_destruction_frequency_s_inv": (
                printed_kd_upper_from_attachment_ratio_cm3_s * 1.0e11
            ),
            "interpretation": (
                "Footnote 14 is algebraically self-consistent as 7e-9 times "
                "1e11 = 700/s, but its kd is over six times the maximum from "
                "the article's printed kdis law across the measured 11 cm, "
                "10 mTorr Te board even after adding the source-stated "
                "maximum 1/7 attachment contribution. The nominal footnote "
                "numbers are quarantined from calibration."
            ),
        },
        "uncertainty_boundary": (
            "The source's +/-25% absolute-density accuracy statement is "
            "propagated as a non-statistical envelope. Te and volume-average "
            "ne lack complete measurement uncertainties, Hamilton publishes "
            "no scalar physical uncertainty, and attachment is compiled; "
            "therefore no formal pass/fail or predictive interval is claimed."
        ),
        "inversion_boundary": (
            "The output is kr required by the paper's fast-reaction Eq. 7. "
            "Mapping kr to gamma requires independently validated neutral "
            "diffusion and wall-state physics. Volume-average ne cannot be "
            "used as a local sheath density or wafer flux."
        ),
        "csv_sha256": hashlib.sha256(
            csv_payload.encode("utf-8")).hexdigest(),
    }
    return csv_payload, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    csv_payload, audit = build_outputs()
    if args.write:
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(csv_payload, encoding="utf-8")
        AUDIT_PATH.write_text(
            json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    else:
        print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
