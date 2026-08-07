#!/usr/bin/env python3
"""Build a fast Maxwellian Cl2 state-rate table from Hamilton's exact data.

The source OPJ was exported with liborigin 3.0.4.  This script integrates the
full 50,000-point state-resolved cross sections analytically, on the authors'
own Figure-5 temperature grid, and checks their independently supplied total
Maxwellian rate.  No reactor observable is used.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from petch.reactor_global import (
    ElectronMaxwellianCrossSectionRateCoefficient,
    RateContext,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRECTORY = ROOT / "research_sources" / "digitized"
STATE_CROSS_SECTIONS = (
    SOURCE_DIRECTORY / "hamilton_2018_cl2_state_cross_sections.csv")
REFERENCE_RATES = (
    SOURCE_DIRECTORY / "hamilton_2018_cl2_reference_rate_coefficients.csv")
PACKAGE_DIRECTORY = ROOT / "src" / "petch" / "reactor_global" / "data"
OUTPUT_RATES = (
    PACKAGE_DIRECTORY / "hamilton_2018_cl2_state_maxwellian_rates.csv")
AUDIT_DIRECTORY = (
    ROOT / "results" / "curated" / "reactor_global_chlorine")
AUDIT_JSON = AUDIT_DIRECTORY / "hamilton_2018_dissociation_rate_audit.json"
AUDIT_MD = AUDIT_DIRECTORY / "HAMILTON_2018_DISSOCIATION_RATE_AUDIT.md"

TEMPERATURE_DOMAIN_EV = (0.3, 5.0)
NUMERICAL_REPRODUCTION_RELATIVE_LIMIT = 0.01
SOURCE_CROSS_SECTION_TAIL_LIMIT = 1.0e-10
STATE_DEFINITIONS = (
    ("a_3Pi_u", "a_3Pi_u_m2", 3.252),
    ("A_1Pi_u", "A_1Pi_u_m2", 4.348),
    ("b_3Pi_g", "b_3Pi_g_m2", 6.498),
    ("B_1Pi_g", "B_1Pi_g_m2", 7.537),
    ("C_1Delta_g", "C_1Delta_g_m2", 7.790),
    ("c_3Sigma_g_minus", "c_3Sigma_g_minus_m2", 7.257),
    ("D_1Sigma_g_plus", "D_1Sigma_g_plus_m2", 8.228),
    ("e_3Sigma_u_plus", "e_3Sigma_u_plus_m2", 9.219),
)


def load_source_cross_sections():
    with STATE_CROSS_SECTIONS.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    energies = tuple(float(row["electron_energy_eV"]) for row in rows)
    providers = {}
    for state, column, threshold_eV in STATE_DEFINITIONS:
        providers[state] = ElectronMaxwellianCrossSectionRateCoefficient(
            electron_energy_eV=energies,
            cross_section_m2=tuple(float(row[column]) for row in rows),
            threshold_eV=threshold_eV,
            relative_uncertainty=None,
            source=(
                "hamilton-2018-cl2-dissociation exact OPJ state cross "
                f"section for {state}"
            ),
            evidence_kind="semi_empirical",
            maximum_kernel_tail_fraction=(
                SOURCE_CROSS_SECTION_TAIL_LIMIT),
        )
    return energies, providers


def load_reference_rates() -> list[dict[str, float]]:
    with REFERENCE_RATES.open(newline="", encoding="utf-8") as stream:
        rows = []
        for row in csv.DictReader(stream):
            temperature = float(row["effective_electron_temperature_eV"])
            if TEMPERATURE_DOMAIN_EV[0] <= temperature <= (
                    TEMPERATURE_DOMAIN_EV[1]):
                rows.append({
                    "electron_temperature_eV": temperature,
                    "hamilton_total_rate_m3_s":
                        float(row["hamilton_x1_m3_s"]),
                })
    if len(rows) < 200:
        raise RuntimeError("Hamilton reference rate domain is incomplete")
    return rows


def build_rate_table():
    _, providers = load_source_cross_sections()
    reference = load_reference_rates()
    reference_by_temperature = {
        row["electron_temperature_eV"]: row["hamilton_total_rate_m3_s"]
        for row in reference
    }
    output_rows = []
    temperatures = sorted({
        TEMPERATURE_DOMAIN_EV[0],
        TEMPERATURE_DOMAIN_EV[1],
        *(row["electron_temperature_eV"] for row in reference),
    })
    for temperature in temperatures:
        context = RateContext(temperature)
        state_rates = {
            state: provider.coefficient_si(context)
            for state, provider in providers.items()
        }
        summed_rate = sum(state_rates.values())
        reference_rate = reference_by_temperature.get(temperature)
        output_rows.append({
            "electron_temperature_eV": temperature,
            **{
                f"{state}_m3_s": state_rates[state]
                for state, _, _ in STATE_DEFINITIONS
            },
            "summed_state_rate_m3_s": summed_rate,
            "hamilton_reference_total_rate_m3_s": reference_rate,
            "relative_reproduction_error": (
                summed_rate / reference_rate - 1.0
                if reference_rate is not None else None),
        })
    return output_rows


def grade_rate_table(rows):
    reference_rows = [
        row for row in rows
        if row["relative_reproduction_error"] is not None]
    errors = np.abs(np.asarray([
        row["relative_reproduction_error"] for row in reference_rows]))
    maximum_index = int(np.argmax(errors))
    maximum_row = reference_rows[maximum_index]
    return {
        "gate": (
            "analytic integration of Hamilton state cross sections versus "
            "the authors' supplied Maxwellian x=1 total rate"),
        "claim_class": (
            "numerical source reproduction, not experimental validation"),
        "coefficient_selection_target": None,
        "temperature_domain_eV": list(TEMPERATURE_DOMAIN_EV),
        "state_count": len(STATE_DEFINITIONS),
        "temperature_points": len(reference_rows),
        "runtime_table_points": len(rows),
        "relative_error_limit": (
            NUMERICAL_REPRODUCTION_RELATIVE_LIMIT),
        "maximum_absolute_relative_error": float(errors[maximum_index]),
        "mean_absolute_relative_error": float(np.mean(errors)),
        "worst_temperature_eV":
            maximum_row["electron_temperature_eV"],
        "passed": bool(
            np.max(errors) <= NUMERICAL_REPRODUCTION_RELATIVE_LIMIT),
        "source_boundary": (
            "Hamilton cross sections are fixed-nuclei R-matrix calculations "
            "for Cl2(v=0), extended above the ionization potential by "
            "transition-specific scaling; no scalar physical uncertainty "
            "was published"),
    }


def write_results(rows, grade):
    PACKAGE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    package_fields = [
        key for key in rows[0]
        if key not in {
            "hamilton_reference_total_rate_m3_s",
            "relative_reproduction_error",
        }
    ]
    with OUTPUT_RATES.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=package_fields, lineterminator="\n",
            extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    AUDIT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    AUDIT_JSON.write_text(
        json.dumps(grade, indent=2) + "\n", encoding="utf-8")
    report = f"""# Hamilton 2018 Cl2 dissociation-rate audit

**Verdict: {'PASS' if grade['passed'] else 'FAIL'}**

The analytic Maxwellian integrator consumed all 50,000 source points for each
of eight dissociative excited states and was compared against the authors'
independently supplied Figure-5 `x=1` total-rate array at
{grade['temperature_points']} temperatures from
{TEMPERATURE_DOMAIN_EV[0]:.1f} to {TEMPERATURE_DOMAIN_EV[1]:.1f} eV.

- maximum absolute relative error:
  `{100.0 * grade['maximum_absolute_relative_error']:.4f}%`
- mean absolute relative error:
  `{100.0 * grade['mean_absolute_relative_error']:.4f}%`
- engineering reproduction limit: `1.0000%`
- worst temperature:
  `{grade['worst_temperature_eV']:.6f} eV`

This is numerical source reproduction, not experimental validation and not a
rate fit. No reactor density, ion flux, etch rate, or feature depth selected a
coefficient.

## Physics gained

Each state retains its Table-2 vertical excitation energy, so particle
production and the electron-energy sink can be summed consistently instead
of treating Lee's `3.824 eV` Arrhenius exponent as a physical event energy.
All eight retained states dissociate to two ground-state Cl atoms.

## Boundary retained

The source uses fixed-nuclei R-matrix calculations for Cl2(v=0), then
transition-specific high-energy scaling. Hamilton et al. explicitly say
Cosby's vibrationally distributed experiment is not directly comparable and
publish no scalar uncertainty. The resulting rates are therefore
`semi_empirical`, Maxwellian-only, and fail outside the paper's stated
industrial 0.3--5 eV domain.
"""
    AUDIT_MD.write_text(report, encoding="utf-8")


def main():
    rows = build_rate_table()
    grade = grade_rate_table(rows)
    write_results(rows, grade)
    if not grade["passed"]:
        raise SystemExit("Hamilton rate reproduction failed")


if __name__ == "__main__":
    main()
