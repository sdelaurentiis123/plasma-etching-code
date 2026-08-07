#!/usr/bin/env python3
"""Invert Mahoney density only to diagnose the missing absorbed-power input.

This is target-informed diagnosis, not validation and not a fitted production
boundary.  Its method is frozen in the adjacent preregistration before use.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from petch.reactor_global import (
    ArgonGlobalCondition,
    CylindricalReactor,
    LeeLiebermanArgonGlobalModel,
    LeeLiebermanArgonTransportProvider,
    PASCAL_PER_MTORR,
)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = (
    ROOT / "research_sources" / "digitized"
    / "mahoney_1994_table1_argon_100W.csv")
OUTPUT_DIRECTORY = ROOT / "results" / "curated" / "reactor_global_argon"
RESULTS_CSV = OUTPUT_DIRECTORY / "mahoney_1994_power_inversion.csv"
SUMMARY_JSON = OUTPUT_DIRECTORY / "mahoney_1994_power_inversion.json"
REPORT_MD = OUTPUT_DIRECTORY / "MAHONEY_1994_POWER_INVERSION.md"

GAS_TEMPERATURE_ENDPOINTS_K = (300.0, 600.0)
WALL_ENERGY_FACTORS_TE = (5.0, 8.0)
DENSITY_MULTIPLIERS = (2.0, 5.0)
NET_RF_POWER_W = 100.0
ROOT_BRACKET_W = (1.0, 500.0)
HOPWOOD_CONTEXT_FRACTION = (0.70, 0.90)


def load_reference() -> list[dict[str, object]]:
    with REFERENCE.open(newline="", encoding="utf-8") as stream:
        rows = []
        for row in csv.DictReader(stream):
            rows.append({
                "condition_id": row["condition_id"],
                "pressure_mTorr": float(row["pressure_mTorr"]),
                "pump_mode": row["pump_mode"],
                "measured_peak_electron_density_m3":
                    float(row["peak_electron_density_1e10_cm3"]) * 1.0e16,
            })
    if len(rows) != 5:
        raise RuntimeError("Mahoney Table I reference does not match freeze")
    return rows


def _condition(
        target: dict[str, object], *, absorbed_power_W: float,
        gas_temperature_K: float, wall_factor: float,
        ) -> ArgonGlobalCondition:
    return ArgonGlobalCondition(
        condition_id=str(target["condition_id"]),
        absorbed_power_W=absorbed_power_W,
        pressure_Pa=(
            float(target["pressure_mTorr"]) * PASCAL_PER_MTORR),
        gas_temperature_K=gas_temperature_K,
        geometry=CylindricalReactor(radius_m=0.114, length_m=0.137),
        ion_wall_energy_factor_Te=wall_factor,
        ion_wall_energy_source=(
            "lee-lieberman-1994-global stated 5--8 Te range"),
        ion_wall_energy_evidence="published_range_member",
        absorbed_power_source=(
            "target-inverted from mahoney-1994-planar-icp density; "
            "diagnostic only"),
        absorbed_power_evidence="sensitivity",
        absorbed_power_boundary_kind="target_inverted_diagnostic",
    )


def _inferred_power_W(
        model: LeeLiebermanArgonGlobalModel, target: dict[str, object],
        *, density_multiplier: float, gas_temperature_K: float,
        wall_factor: float) -> float:
    target_density = (
        density_multiplier
        * float(target["measured_peak_electron_density_m3"]))

    def log_density_error(log_power: float) -> float:
        power = float(np.exp(log_power))
        solution = model.solve(_condition(
            target,
            absorbed_power_W=power,
            gas_temperature_K=gas_temperature_K,
            wall_factor=wall_factor,
        ),
            initial_electron_density_m3=1.0e15,
            initial_metastable_density_m3=1.0e13,
        )
        return float(np.log(solution.electron_density_m3 / target_density))

    lower_log = float(np.log(ROOT_BRACKET_W[0]))
    upper_log = float(np.log(ROOT_BRACKET_W[1]))
    lower_error = log_density_error(lower_log)
    upper_error = log_density_error(upper_log)
    if lower_error >= 0.0 or upper_error <= 0.0:
        raise RuntimeError(
            "absorbed-power inversion did not bracket the target density")
    return float(np.exp(brentq(
        log_density_error,
        lower_log,
        upper_log,
        xtol=1.0e-12,
        rtol=1.0e-12,
    )))


def run_diagnostic() -> tuple[list[dict[str, object]], dict[str, object]]:
    reference = load_reference()
    model = LeeLiebermanArgonGlobalModel(
        LeeLiebermanArgonTransportProvider())
    rows: list[dict[str, object]] = []
    members: list[dict[str, object]] = []

    for gas_temperature_K in GAS_TEMPERATURE_ENDPOINTS_K:
        for wall_factor in WALL_ENERGY_FACTORS_TE:
            member_rows = []
            for target in reference:
                lower_power = _inferred_power_W(
                    model,
                    target,
                    density_multiplier=DENSITY_MULTIPLIERS[0],
                    gas_temperature_K=gas_temperature_K,
                    wall_factor=wall_factor,
                )
                upper_power = _inferred_power_W(
                    model,
                    target,
                    density_multiplier=DENSITY_MULTIPLIERS[1],
                    gas_temperature_K=gas_temperature_K,
                    wall_factor=wall_factor,
                )
                row = {
                    "condition_id": target["condition_id"],
                    "pump_mode": target["pump_mode"],
                    "pressure_mTorr": target["pressure_mTorr"],
                    "gas_temperature_K": gas_temperature_K,
                    "wall_energy_factor_Te": wall_factor,
                    "absorbed_power_for_2x_density_W": lower_power,
                    "absorbed_power_for_5x_density_W": upper_power,
                    "transfer_fraction_for_2x_density":
                        lower_power / NET_RF_POWER_W,
                    "transfer_fraction_for_5x_density":
                        upper_power / NET_RF_POWER_W,
                    "upper_endpoint_exceeds_net_power":
                        upper_power > NET_RF_POWER_W,
                }
                rows.append(row)
                member_rows.append(row)

            intersection_lower = max(
                float(row["transfer_fraction_for_2x_density"])
                for row in member_rows)
            intersection_upper = min(
                min(
                    float(row["transfer_fraction_for_5x_density"]),
                    1.0,
                )
                for row in member_rows
            )
            intersection_exists = intersection_lower <= intersection_upper
            hopwood_overlap_lower = max(
                intersection_lower, HOPWOOD_CONTEXT_FRACTION[0])
            hopwood_overlap_upper = min(
                intersection_upper, HOPWOOD_CONTEXT_FRACTION[1])
            member = {
                "gas_temperature_K": gas_temperature_K,
                "wall_energy_factor_Te": wall_factor,
                "constant_transfer_fraction_intersection": (
                    [intersection_lower, intersection_upper]
                    if intersection_exists else None
                ),
                "constant_transfer_fraction_exists": intersection_exists,
                "overlaps_hopwood_70_90_percent_context": bool(
                    intersection_exists
                    and hopwood_overlap_lower <= hopwood_overlap_upper),
            }
            members.append(member)

    summary = {
        "diagnostic": (
            "Mahoney 1994 target-inverted absorbed-power intersection"),
        "claim_class": (
            "target-informed diagnostic; not validation and not a "
            "production power boundary"),
        "coefficient_selection_target": (
            "Mahoney peak electron density times source-stated 2--5 "
            "diagnostic interval"),
        "net_rf_power_W": NET_RF_POWER_W,
        "density_multiplier_interval": list(DENSITY_MULTIPLIERS),
        "root_bracket_W": list(ROOT_BRACKET_W),
        "hopwood_external_context_fraction":
            list(HOPWOOD_CONTEXT_FRACTION),
        "members": members,
    }
    return rows, summary


def write_results(
        rows: list[dict[str, object]], summary: dict[str, object]) -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    SUMMARY_JSON.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = []
    for member in summary["members"]:
        interval = member["constant_transfer_fraction_intersection"]
        interval_text = (
            f"{100.0 * interval[0]:.1f}–{100.0 * interval[1]:.1f}%"
            if interval is not None else "none"
        )
        lines.append(
            f"| {member['gas_temperature_K']:.0f} | "
            f"{member['wall_energy_factor_Te']:.0f} | "
            f"{interval_text} | "
            f"{'yes' if member['overlaps_hopwood_70_90_percent_context'] else 'no'} |"
        )
    report = f"""# Mahoney 1994 absorbed-power diagnostic inversion

**Claim class: target-informed diagnosis, not validation**

The unchanged argon model was inverted only to find absorbed-power intervals
that place its center density at two to five times Mahoney's measured peak
electron density. That multiplier is the source's stated electron-versus-ion
diagnostic interval.

| gas T (K) | ion wall energy (Te) | one constant transfer fraction satisfying all 5 rows | overlaps Hopwood 70–90% context |
|---:|---:|---:|:---:|
{chr(10).join(lines)}

The intersection uses all five rows and is clipped at the physical 100 W net
RF ceiling. Hopwood's 70–90% range comes from a different planar ICP and is
shown only as external context.

These fractions were selected from the Mahoney density and therefore may not
be used as a predictive power provider or to relabel the frozen independent
board. The diagnostic asks whether a future, independently measured
hardware-loss fraction could be represented as one constant under each
declared gas-temperature/wall-energy sensitivity corner.
"""
    REPORT_MD.write_text(report, encoding="utf-8")


def main() -> None:
    rows, summary = run_diagnostic()
    write_results(rows, summary)


if __name__ == "__main__":
    main()
