#!/usr/bin/env python3
"""Bound Cl- electron-detachment leverage on the atomic-Cl source replay.

This is a frozen-state timescale audit, not a coupled reactor solve. It uses
the Kemaneci 2014 Table 4 Maxwellian fits only inside their 0.5--10 eV domain
and compares them with the already active Lee/Kemaneci mutual-neutralization
loss. No feature observable participates.
"""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"
INPUT = OUTPUT / "malyshev_1998_eedf_hamilton_atomic_cl_source_replay.json"
JSON_OUTPUT = OUTPUT / "malyshev_1998_detachment_importance.json"
REPORT_OUTPUT = OUTPUT / "MALYSHEV_1998_DETACHMENT_IMPORTANCE.md"
MUTUAL_NEUTRALIZATION_M3_S = 5.0e-14


def _single_detachment_m3_s(temperature_eV: float) -> float:
    return 9.02e-15 * temperature_eV ** 0.92 * math.exp(
        -4.88 / temperature_eV)


def _double_detachment_m3_s(temperature_eV: float) -> float:
    return 3.62e-15 * temperature_eV ** 0.72 * math.exp(
        -25.38 / temperature_eV)


def audit() -> dict[str, object]:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = []
    for source_row in source["rows"]:
        temperature = float(source_row["mean_energy_temperature_proxy_eV"])
        if not 0.5 <= temperature <= 10.0:
            raise RuntimeError("Kemaneci detachment fit is out of domain")
        electron_density = float(source_row["electron_density_m3"])
        electronegativity = float(source_row["electronegativity"])
        single = _single_detachment_m3_s(temperature)
        double = _double_detachment_m3_s(temperature)
        electron_detachment_frequency = electron_density * (single + double)
        # Quasineutrality gives n(Cl2+) + n(Cl+) = ne + n(Cl-).
        neutralization_frequency = (
            MUTUAL_NEUTRALIZATION_M3_S
            * electron_density * (1.0 + electronegativity)
        )
        rows.append({
            "absorbed_fraction_sensitivity": source_row[
                "absorbed_fraction_sensitivity"],
            "source_power_W": source_row["source_power_W"],
            "mean_energy_temperature_proxy_eV": temperature,
            "electron_density_m3": electron_density,
            "electronegativity": electronegativity,
            "single_detachment_rate_coefficient_m3_s": single,
            "double_detachment_rate_coefficient_m3_s": double,
            "electron_detachment_loss_frequency_s_inv": (
                electron_detachment_frequency),
            "mutual_neutralization_loss_frequency_s_inv": (
                neutralization_frequency),
            "detachment_to_neutralization_loss_ratio": (
                electron_detachment_frequency / neutralization_frequency),
        })
    maximum_ratio = max(
        row["detachment_to_neutralization_loss_ratio"] for row in rows)
    return {
        "schema": "petch.malyshev_1998_detachment_importance.v1",
        "claim_class": "frozen-state source-replay timescale sensitivity",
        "source_replay": str(INPUT.relative_to(ROOT)),
        "rate_source": (
            "Kemaneci et al. 2014 Table 4 reactions 23--24; rates trace "
            "Fritioff et al. 2003 electron-impact detachment measurements"
        ),
        "single_detachment_formula": (
            "9.02e-15 Te^0.92 exp(-4.88/Te) m3/s"),
        "double_detachment_formula": (
            "3.62e-15 Te^0.72 exp(-25.38/Te) m3/s"),
        "fit_temperature_domain_eV": [0.5, 10.0],
        "temperature_input_boundary": (
            "2/3 of non-Maxwellian EEPF mean energy; Maxwellian-fit "
            "sensitivity, not an exact EEPF collision moment"),
        "mutual_neutralization_rate_coefficient_m3_s": (
            MUTUAL_NEUTRALIZATION_M3_S),
        "maximum_detachment_to_neutralization_loss_ratio": maximum_ratio,
        "verdict": (
            "low frozen-state leverage: detachment is below 4% of the "
            "already active mutual-neutralization loss on this board"),
        "coupled_response_is_formal_bound": False,
        "feature_depth_used": False,
        "supports_reactor_state_prediction": False,
        "supports_wafer_flux": False,
        "supports_feature_depth": False,
        "rows": rows,
    }


def main() -> None:
    result = audit()
    JSON_OUTPUT.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    table = []
    for row in result["rows"]:
        table.append(
            "| {absorbed_fraction_sensitivity:.2f} | {source_power_W:.0f} | "
            "{mean_energy_temperature_proxy_eV:.3f} | {electronegativity:.3f} | "
            "{detachment_to_neutralization_loss_ratio:.4f} |".format(**row)
        )
    REPORT_OUTPUT.write_text(
        "# Malyshev source-replay detachment importance\n\n"
        "## Verdict\n\n"
        "The Kemaneci/Fritioff electron-detachment fits contribute only "
        f"{100.0 * result['maximum_detachment_to_neutralization_loss_ratio']:.1f}% "
        "of the already active mutual-neutralization loss at worst. Electron "
        "detachment is therefore a low-leverage explanation for the current "
        "17--24 percentage-point Cl2 composition-proxy miss. This is a "
        "frozen-state timescale diagnostic, not a formal bound on the fully "
        "coupled response.\n\n"
        "| absorbed fraction | source W | 2/3 mean-E eV | electronegativity | "
        "detachment / neutralization |\n"
        "|---:|---:|---:|---:|---:|\n"
        + "\n".join(table)
        + "\n\nThe rates are Maxwellian fits evaluated at the EEPF "
        "mean-energy-equivalent temperature. No feature depth or reactor "
        "diagnostic selected a coefficient.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
