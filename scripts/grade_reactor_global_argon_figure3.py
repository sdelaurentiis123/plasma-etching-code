#!/usr/bin/env python3
"""Grade the no-fit argon global model against Lee--Lieberman Figure 3."""
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
DIGITIZED = (
    ROOT / "research_sources" / "digitized"
    / "lee_lieberman_1994_figure3_argon.csv")
OUTPUT_DIRECTORY = ROOT / "results" / "curated" / "reactor_global_argon"
RESULTS_CSV = OUTPUT_DIRECTORY / "figure3_reproduction.csv"
GRADE_JSON = OUTPUT_DIRECTORY / "figure3_grade.json"
REPORT_MD = OUTPUT_DIRECTORY / "FIGURE3_REPRODUCTION.md"

WALL_ENERGY_FACTORS = (5.0, 8.0)
MEAN_ABSOLUTE_PERCENT_ERROR_LIMIT = 10.0
MAXIMUM_ABSOLUTE_PERCENT_ERROR_LIMIT = 20.0
RESIDUAL_LIMIT = 1.0e-8


def load_reference() -> list[dict[str, float | bool]]:
    with DIGITIZED.open(newline="", encoding="utf-8") as stream:
        rows = []
        for row in csv.DictReader(stream):
            rows.append({
                "pressure_mTorr": float(row["pressure_mTorr"]),
                "reference_electron_temperature_eV":
                    float(row["electron_temperature_eV"]),
                "reference_strokes_overlap":
                    row["argon_oxygen_strokes_overlap"] == "True",
            })
    if len(rows) < 3:
        raise RuntimeError("Figure 3 reference is incomplete")
    return rows


def run_reproduction() -> tuple[list[dict[str, object]], dict[str, object]]:
    reference = load_reference()
    provider = LeeLiebermanArgonTransportProvider()
    output_rows: list[dict[str, object]] = []
    member_grades: list[dict[str, object]] = []

    for wall_factor in WALL_ENERGY_FACTORS:
        previous = {
            "electron_temperature_eV": 3.0,
            "electron_density_m3": 1.0e17,
            "metastable_density_m3": 1.0e15,
        }
        member_rows: list[dict[str, object]] = []
        for target in reference:
            pressure_mTorr = float(target["pressure_mTorr"])
            condition = ArgonGlobalCondition(
                condition_id=(
                    f"lee-lieberman-figure3-{pressure_mTorr:g}-mTorr"
                    f"-wall-{wall_factor:g}Te"
                ),
                absorbed_power_W=1000.0,
                pressure_Pa=pressure_mTorr * PASCAL_PER_MTORR,
                gas_temperature_K=600.0,
                geometry=CylindricalReactor(
                    radius_m=0.1525, length_m=0.075),
                ion_wall_energy_factor_Te=wall_factor,
                ion_wall_energy_source=(
                    "lee-lieberman-1994-global stated 5--8 Te range"),
                ion_wall_energy_evidence="published_range_member",
            )
            solution = LeeLiebermanArgonGlobalModel(provider).solve(
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
            reference_temperature = float(
                target["reference_electron_temperature_eV"])
            percent_error = 100.0 * (
                solution.electron_temperature_eV - reference_temperature
            ) / reference_temperature
            row = {
                "wall_energy_factor_Te": wall_factor,
                "pressure_mTorr": pressure_mTorr,
                "reference_electron_temperature_eV": reference_temperature,
                "reference_strokes_overlap":
                    bool(target["reference_strokes_overlap"]),
                "model_electron_temperature_eV":
                    solution.electron_temperature_eV,
                "percent_error": percent_error,
                "electron_density_m3": solution.electron_density_m3,
                "metastable_density_m3": solution.metastable_density_m3,
                "axial_ion_flux_m2_s": solution.axial_ion_flux_m2_s,
                "ion_mean_free_path_m":
                    solution.transport.ion_mean_free_path_m,
                "ambipolar_diffusion_m2_s":
                    solution.transport.ambipolar_diffusion_m2_s,
                "metastable_effective_diffusion_m2_s":
                    solution.transport.metastable_effective_diffusion_m2_s,
                "maximum_normalized_residual":
                    solution.maximum_normalized_residual,
                "supports_prediction": solution.supports_prediction,
            }
            member_rows.append(row)
            output_rows.append(row)

        temperatures = np.asarray([
            float(row["model_electron_temperature_eV"])
            for row in member_rows
        ])
        absolute_errors = np.abs(np.asarray([
            float(row["percent_error"]) for row in member_rows
        ]))
        maximum_residual = max(
            float(row["maximum_normalized_residual"])
            for row in member_rows)
        member_grade = {
            "wall_energy_factor_Te": wall_factor,
            "mean_absolute_percent_error": float(np.mean(absolute_errors)),
            "maximum_absolute_percent_error": float(np.max(absolute_errors)),
            "monotonically_decreasing_with_pressure":
                bool(np.all(np.diff(temperatures) < 0.0)),
            "maximum_normalized_residual": maximum_residual,
            "all_positive_finite": bool(all(
                np.isfinite(float(row[key])) and float(row[key]) > 0.0
                for row in member_rows
                for key in (
                    "model_electron_temperature_eV",
                    "electron_density_m3",
                    "metastable_density_m3",
                    "axial_ion_flux_m2_s",
                )
            )),
        }
        member_grade["passed"] = bool(
            member_grade["mean_absolute_percent_error"]
            <= MEAN_ABSOLUTE_PERCENT_ERROR_LIMIT
            and member_grade["maximum_absolute_percent_error"]
            <= MAXIMUM_ABSOLUTE_PERCENT_ERROR_LIMIT
            and member_grade["monotonically_decreasing_with_pressure"]
            and member_grade["maximum_normalized_residual"] <= RESIDUAL_LIMIT
            and member_grade["all_positive_finite"]
        )
        member_grades.append(member_grade)

    grade = {
        "gate": "Lee--Lieberman 1994 Figure 3 pure-Ar reproduction",
        "claim_class": "published-model reproduction, not independent validation",
        "coefficient_selection_target": None,
        "reference_points": len(reference),
        "wall_energy_range_members_Te": list(WALL_ENERGY_FACTORS),
        "limits": {
            "mean_absolute_percent_error":
                MEAN_ABSOLUTE_PERCENT_ERROR_LIMIT,
            "maximum_absolute_percent_error":
                MAXIMUM_ABSOLUTE_PERCENT_ERROR_LIMIT,
            "maximum_normalized_residual": RESIDUAL_LIMIT,
            "monotonic_decrease_required": True,
            "positive_finite_outputs_required": True,
        },
        "members": member_grades,
        "worst_member_mean_absolute_percent_error": max(
            float(member["mean_absolute_percent_error"])
            for member in member_grades),
        "worst_member_maximum_absolute_percent_error": max(
            float(member["maximum_absolute_percent_error"])
            for member in member_grades),
        "passed": all(bool(member["passed"]) for member in member_grades),
        "predictive_boundary": (
            "Transport is evidence-classed published_model because the NIST "
            "Ar self-diffusion correlation is extrapolated from 418 K to the "
            "source model's 600 K. Independent reactor data remain required."
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
            f"| {member['wall_energy_factor_Te']:.1f} | "
            f"{member['mean_absolute_percent_error']:.3f}% | "
            f"{member['maximum_absolute_percent_error']:.3f}% | "
            f"{'yes' if member['monotonically_decreasing_with_pressure'] else 'no'} | "
            f"{member['maximum_normalized_residual']:.3e} | "
            f"{'PASS' if member['passed'] else 'FAIL'} |"
        )
    report = f"""# Lee--Lieberman Figure 3 argon reproduction

**Verdict: {'PASS' if grade['passed'] else 'FAIL'}**

This is a no-fit reproduction of the published global-model curve, not an
independent reactor validation. The 18 reference points were digitized before
running this grade. Neither a target electron temperature, reactor density,
ion flux, etch rate, nor feature depth selected any coefficient.

| ion wall energy member | MAPE | maximum APE | monotonic | maximum balance residual | verdict |
|---:|---:|---:|:---:|---:|:---:|
{chr(10).join(member_lines)}

Frozen limits were MAPE <= {MEAN_ABSOLUTE_PERCENT_ERROR_LIMIT:.0f}%, maximum
APE <= {MAXIMUM_ABSOLUTE_PERCENT_ERROR_LIMIT:.0f}%, normalized particle/power
residual <= {RESIDUAL_LIMIT:.0e}, strictly decreasing electron temperature
with pressure, and positive finite densities/fluxes. Both endpoints of the
source's published 5--8 Te ion-wall-energy range had to pass.

## Claim boundary

The transport closure is still labeled `published_model`: Phelps' Ar+-Ar law
is source-backed and energy-averaged, but the NIST Ar-in-Ar self-diffusion
correlation is extrapolated from its stated 418 K ceiling to Lee and
Lieberman's 600 K condition. Independent, condition-specific reactor
measurements are the next gate. No Krueger depth result is fitted or changed
by this board.
"""
    REPORT_MD.write_text(report, encoding="utf-8")


def main() -> None:
    rows, grade = run_reproduction()
    write_results(rows, grade)
    if not grade["passed"]:
        raise SystemExit("Figure 3 reproduction failed the frozen gate")


if __name__ == "__main__":
    main()
