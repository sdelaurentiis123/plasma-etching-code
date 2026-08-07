#!/usr/bin/env python3
"""Grade the unchanged argon global model on Mahoney et al. 1994.

The condition, sensitivities, and limits are frozen in
``results/curated/reactor_global_argon/PREREGISTRATION.md``.  In particular,
the paper reports net rf generator power rather than absorbed plasma power.
This board therefore treats 100 W as an all-net-power-absorbed upper-bound
scenario; it does not infer or fit a transfer efficiency.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

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
RESULTS_CSV = OUTPUT_DIRECTORY / "mahoney_1994_upper_bound.csv"
GRADE_JSON = OUTPUT_DIRECTORY / "mahoney_1994_upper_bound_grade.json"
REPORT_MD = OUTPUT_DIRECTORY / "MAHONEY_1994_INDEPENDENT.md"

GAS_TEMPERATURE_ENDPOINTS_K = (300.0, 600.0)
WALL_ENERGY_FACTORS_TE = (5.0, 8.0)
TEMPERATURE_MAPE_LIMIT_PERCENT = 30.0
TEMPERATURE_MAX_APE_LIMIT_PERCENT = 50.0
DENSITY_DIAGNOSTIC_RATIO_RANGE = (1.0, 5.0)
DENSITY_SHAPE_LOG_RMSE_LIMIT = float(np.log(2.0))
RESIDUAL_LIMIT = 1.0e-8
SHAPE_SEQUENCE = (
    "mahoney-10mT-cryo-100W",
    "mahoney-20mT-mechanical-100W",
    "mahoney-50mT-mechanical-100W",
    "mahoney-100mT-mechanical-100W",
)


def load_reference() -> list[dict[str, object]]:
    with REFERENCE.open(newline="", encoding="utf-8") as stream:
        rows = []
        for row in csv.DictReader(stream):
            rows.append({
                "condition_id": row["condition_id"],
                "pressure_mTorr": float(row["pressure_mTorr"]),
                "pump_mode": row["pump_mode"],
                "net_rf_power_W": float(row["net_rf_power_W"]),
                "measured_peak_electron_density_m3":
                    float(row["peak_electron_density_1e10_cm3"]) * 1.0e16,
                "measured_bulk_electron_temperature_eV":
                    float(row["bulk_electron_temperature_eV"]),
            })
    if (
        len(rows) != 5
        or {str(row["condition_id"]) for row in rows}
        != {
            "mahoney-10mT-cryo-100W",
            "mahoney-20mT-cryo-100W",
            "mahoney-20mT-mechanical-100W",
            "mahoney-50mT-mechanical-100W",
            "mahoney-100mT-mechanical-100W",
        }
        or any(float(row["net_rf_power_W"]) != 100.0 for row in rows)
    ):
        raise RuntimeError("Mahoney Table I reference does not match freeze")
    return rows


def _shape_rows(
        member_rows: list[dict[str, object]],
        ) -> list[dict[str, object]]:
    by_id = {
        str(row["condition_id"]): row
        for row in member_rows
    }
    return [by_id[condition_id] for condition_id in SHAPE_SEQUENCE]


def run_board() -> tuple[list[dict[str, object]], dict[str, object]]:
    reference = load_reference()
    provider = LeeLiebermanArgonTransportProvider()
    model = LeeLiebermanArgonGlobalModel(provider)
    output_rows: list[dict[str, object]] = []
    members: list[dict[str, object]] = []

    for gas_temperature_K in GAS_TEMPERATURE_ENDPOINTS_K:
        for wall_factor in WALL_ENERGY_FACTORS_TE:
            previous = {
                "electron_temperature_eV": 3.0,
                "electron_density_m3": 1.0e17,
                "metastable_density_m3": 1.0e15,
            }
            member_rows: list[dict[str, object]] = []
            for target in reference:
                condition = ArgonGlobalCondition(
                    condition_id=str(target["condition_id"]),
                    absorbed_power_W=float(target["net_rf_power_W"]),
                    pressure_Pa=(
                        float(target["pressure_mTorr"])
                        * PASCAL_PER_MTORR
                    ),
                    gas_temperature_K=gas_temperature_K,
                    geometry=CylindricalReactor(
                        radius_m=0.114,
                        length_m=0.137,
                    ),
                    ion_wall_energy_factor_Te=wall_factor,
                    ion_wall_energy_source=(
                        "lee-lieberman-1994-global stated 5--8 Te range"
                    ),
                    ion_wall_energy_evidence="published_range_member",
                    absorbed_power_source=(
                        "mahoney-1994-planar-icp net incident-minus-reflected "
                        "RF power, run as all-absorbed upper bound"
                    ),
                    absorbed_power_evidence="sensitivity",
                    absorbed_power_boundary_kind=(
                        "matched_rf_net_power_upper_bound"
                    ),
                )
                solution = model.solve(
                    condition,
                    initial_electron_temperature_eV=float(
                        previous["electron_temperature_eV"]),
                    initial_electron_density_m3=float(
                        previous["electron_density_m3"]),
                    initial_metastable_density_m3=float(
                        previous["metastable_density_m3"]),
                )
                previous = {
                    "electron_temperature_eV":
                        solution.electron_temperature_eV,
                    "electron_density_m3": solution.electron_density_m3,
                    "metastable_density_m3": solution.metastable_density_m3,
                }
                measured_temperature = float(
                    target["measured_bulk_electron_temperature_eV"])
                measured_density = float(
                    target["measured_peak_electron_density_m3"])
                temperature_error_percent = 100.0 * (
                    solution.electron_temperature_eV - measured_temperature
                ) / measured_temperature
                density_ratio = solution.electron_density_m3 / measured_density
                row = {
                    "condition_id": target["condition_id"],
                    "pump_mode": target["pump_mode"],
                    "pressure_mTorr": target["pressure_mTorr"],
                    "net_rf_power_W": target["net_rf_power_W"],
                    "assumed_absorbed_power_upper_bound_W":
                        target["net_rf_power_W"],
                    "gas_temperature_K": gas_temperature_K,
                    "wall_energy_factor_Te": wall_factor,
                    "measured_bulk_electron_temperature_eV":
                        measured_temperature,
                    "model_electron_temperature_eV":
                        solution.electron_temperature_eV,
                    "temperature_error_percent": temperature_error_percent,
                    "measured_peak_electron_density_m3": measured_density,
                    "model_center_electron_density_m3":
                        solution.electron_density_m3,
                    "model_to_measured_density_ratio": density_ratio,
                    "model_metastable_density_m3":
                        solution.metastable_density_m3,
                    "model_axial_ion_flux_m2_s":
                        solution.axial_ion_flux_m2_s,
                    "maximum_normalized_residual":
                        solution.maximum_normalized_residual,
                    "transport_evidence": solution.transport.evidence_kind,
                    "absorbed_power_source":
                        solution.absorbed_power_source,
                    "absorbed_power_evidence":
                        solution.absorbed_power_evidence,
                    "absorbed_power_boundary_kind":
                        solution.absorbed_power_boundary_kind,
                    "supports_prediction": solution.supports_prediction,
                }
                member_rows.append(row)
                output_rows.append(row)

            temperature_absolute_errors = np.abs(np.asarray([
                float(row["temperature_error_percent"])
                for row in member_rows
            ]))
            ordered = _shape_rows(member_rows)
            measured_density = np.asarray([
                float(row["measured_peak_electron_density_m3"])
                for row in ordered
            ])
            modeled_density = np.asarray([
                float(row["model_center_electron_density_m3"])
                for row in ordered
            ])
            measured_shape = measured_density / measured_density[0]
            modeled_shape = modeled_density / modeled_density[0]
            shape_log_rmse = float(np.sqrt(np.mean(
                np.log(modeled_shape / measured_shape) ** 2
            )))
            modeled_temperature = np.asarray([
                float(row["model_electron_temperature_eV"])
                for row in ordered
            ])
            density_ratios = np.asarray([
                float(row["model_to_measured_density_ratio"])
                for row in member_rows
            ])
            maximum_residual = max(
                float(row["maximum_normalized_residual"])
                for row in member_rows
            )
            member = {
                "gas_temperature_K": gas_temperature_K,
                "wall_energy_factor_Te": wall_factor,
                "temperature_mean_absolute_percent_error":
                    float(np.mean(temperature_absolute_errors)),
                "temperature_maximum_absolute_percent_error":
                    float(np.max(temperature_absolute_errors)),
                "temperature_strictly_decreases_with_pressure":
                    bool(np.all(np.diff(modeled_temperature) < 0.0)),
                "density_nondecreasing_with_pressure":
                    bool(np.all(np.diff(modeled_density) >= 0.0)),
                "minimum_model_to_measured_density_ratio":
                    float(np.min(density_ratios)),
                "maximum_model_to_measured_density_ratio":
                    float(np.max(density_ratios)),
                "all_density_ratios_inside_source_diagnostic_interval":
                    bool(np.all(
                        (density_ratios
                         >= DENSITY_DIAGNOSTIC_RATIO_RANGE[0])
                        & (density_ratios
                           <= DENSITY_DIAGNOSTIC_RATIO_RANGE[1])
                    )),
                "normalized_density_shape_log_rmse": shape_log_rmse,
                "maximum_normalized_residual": maximum_residual,
                "all_positive_finite": bool(all(
                    np.isfinite(float(row[key])) and float(row[key]) > 0.0
                    for row in member_rows
                    for key in (
                        "model_electron_temperature_eV",
                        "model_center_electron_density_m3",
                        "model_metastable_density_m3",
                        "model_axial_ion_flux_m2_s",
                    )
                )),
            }
            member["passed"] = bool(
                member["temperature_mean_absolute_percent_error"]
                <= TEMPERATURE_MAPE_LIMIT_PERCENT
                and member["temperature_maximum_absolute_percent_error"]
                <= TEMPERATURE_MAX_APE_LIMIT_PERCENT
                and member["temperature_strictly_decreases_with_pressure"]
                and member["density_nondecreasing_with_pressure"]
                and member[
                    "all_density_ratios_inside_source_diagnostic_interval"
                ]
                and member["normalized_density_shape_log_rmse"]
                <= DENSITY_SHAPE_LOG_RMSE_LIMIT
                and member["maximum_normalized_residual"] <= RESIDUAL_LIMIT
                and member["all_positive_finite"]
            )
            members.append(member)

    grade = {
        "gate": (
            "Mahoney et al. 1994 independent pure-Ar ICP plasma-state board"
        ),
        "claim_class": (
            "independent condition-specific validation under a declared "
            "all-net-rf-power-absorbed upper-bound scenario"
        ),
        "coefficient_selection_target": None,
        "reference_rows": len(reference),
        "gas_temperature_sensitivity_endpoints_K":
            list(GAS_TEMPERATURE_ENDPOINTS_K),
        "wall_energy_range_members_Te": list(WALL_ENERGY_FACTORS_TE),
        "limits": {
            "temperature_mean_absolute_percent_error":
                TEMPERATURE_MAPE_LIMIT_PERCENT,
            "temperature_maximum_absolute_percent_error":
                TEMPERATURE_MAX_APE_LIMIT_PERCENT,
            "model_to_measured_density_ratio":
                list(DENSITY_DIAGNOSTIC_RATIO_RANGE),
            "normalized_density_shape_log_rmse":
                DENSITY_SHAPE_LOG_RMSE_LIMIT,
            "maximum_normalized_residual": RESIDUAL_LIMIT,
            "temperature_strict_decrease_required": True,
            "density_nondecrease_required": True,
            "positive_finite_outputs_required": True,
        },
        "members": members,
        "passed": all(bool(member["passed"]) for member in members),
        "power_boundary": (
            "Mahoney reports net generator power, not absorbed plasma power. "
            "No transfer fraction is inferred here; 100 W is an upper-bound "
            "scenario and cannot establish absolute density prediction."
        ),
        "measurement_boundary": (
            "Mahoney states the electron-density diagnostic can read 2--5x "
            "below ion-density determinations. The board freezes [1,5] for "
            "model center density / measured peak electron density."
        ),
        "predictive_boundary": (
            "A pass validates the unchanged argon closure only within the "
            "declared temperature, wall-energy, diagnostic, and power-input "
            "brackets. It does not validate absorbed-power transfer, C4F6 "
            "chemistry, or Krueger depth."
        ),
    }
    return output_rows, grade


def write_results(
        rows: list[dict[str, object]], grade: dict[str, object]) -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    GRADE_JSON.write_text(
        json.dumps(grade, indent=2) + "\n", encoding="utf-8")

    member_lines = []
    for member in grade["members"]:
        member_lines.append(
            f"| {member['gas_temperature_K']:.0f} | "
            f"{member['wall_energy_factor_Te']:.0f} | "
            f"{member['temperature_mean_absolute_percent_error']:.2f}% | "
            f"{member['temperature_maximum_absolute_percent_error']:.2f}% | "
            f"{member['minimum_model_to_measured_density_ratio']:.2f}–"
            f"{member['maximum_model_to_measured_density_ratio']:.2f} | "
            f"{member['normalized_density_shape_log_rmse']:.3f} | "
            f"{member['maximum_normalized_residual']:.2e} | "
            f"{'PASS' if member['passed'] else 'FAIL'} |"
        )
    report = f"""# Mahoney 1994 independent argon ICP board

**Verdict: {'PASS' if grade['passed'] else 'FAIL'}**

The unchanged five-reaction argon closure was run on the chamber geometry and
all five 100 W operating rows transcribed from Mahoney et al. Table I. No
Mahoney temperature, density, flux, or profile selected a coefficient.

| gas T (K) | ion wall energy (Te) | Te MAPE | Te max APE | model/measured density | density-shape log RMSE | max residual | verdict |
|---:|---:|---:|---:|---:|---:|---:|:---:|
{chr(10).join(member_lines)}

Frozen limits were Te MAPE <= {TEMPERATURE_MAPE_LIMIT_PERCENT:.0f}%, Te maximum
APE <= {TEMPERATURE_MAX_APE_LIMIT_PERCENT:.0f}%, model/measured peak-density
ratio inside [{DENSITY_DIAGNOSTIC_RATIO_RANGE[0]:.0f},
{DENSITY_DIAGNOSTIC_RATIO_RANGE[1]:.0f}], normalized density-shape log RMSE
<= ln(2), the source trends, positive finite state variables, and normalized
balance residual <= {RESIDUAL_LIMIT:.0e}. Both neutral-temperature sensitivity
endpoints and both published ion-wall-energy endpoints had to pass.

## Claim boundary

Mahoney reports net generator power after reflected-power subtraction, not
calorimetric absorbed plasma power, and explicitly notes coil/matching-network
losses. This board therefore runs 100 W as an all-net-power-absorbed
upper-bound scenario. It cannot by itself validate absolute plasma density or
a net-to-absorbed transfer closure.

The paper also states that its electron-density diagnostic can read two to five
times below ion-density determinations. That source-declared interval is the
absolute-density comparison boundary; it was frozen before the run.

This is an independent plasma-state test of the argon chemistry/transport
closure. It does not validate C4F6 chemistry, a sheath IEAD, or Krueger depth.
"""
    REPORT_MD.write_text(report, encoding="utf-8")


def main() -> None:
    rows, grade = run_board()
    write_results(rows, grade)
    if not grade["passed"]:
        raise SystemExit("Mahoney independent argon board failed frozen gate")


if __name__ == "__main__":
    main()
