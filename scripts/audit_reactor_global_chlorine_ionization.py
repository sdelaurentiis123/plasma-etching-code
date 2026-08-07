#!/usr/bin/env python3
"""Compare legacy chlorine ionization fits with evaluated cross sections.

This audit does not tune either representation and defines no retrospective
pass band.  It quantifies what changes when the evaluated NIST/Hayes evidence
is used, and records the molecular-ion branching information that remains
missing.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from petch.reactor_global import (
    RateContext,
    build_lee_lieberman_chlorine_particle_network,
    nist_hayes_atomic_chlorine_ionization_rate,
    nist_molecular_chlorine_total_ionization_rate,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = (
    ROOT / "results" / "curated" / "reactor_global_chlorine")
RESULTS_CSV = OUTPUT_DIRECTORY / "evaluated_ionization_comparison.csv"
SUMMARY_JSON = OUTPUT_DIRECTORY / "evaluated_ionization_comparison.json"
REPORT_MD = OUTPUT_DIRECTORY / "EVALUATED_IONIZATION_COMPARISON.md"
TEMPERATURES_EV = (2.0, 3.0, 4.0, 5.0)


def run_audit() -> tuple[list[dict[str, float]], dict[str, object]]:
    network = build_lee_lieberman_chlorine_particle_network()
    reactions = {reaction.name: reaction for reaction in network.reactions}
    evaluated_atomic = nist_hayes_atomic_chlorine_ionization_rate()
    evaluated_molecular = nist_molecular_chlorine_total_ionization_rate()

    rows: list[dict[str, float]] = []
    for temperature_eV in TEMPERATURES_EV:
        context = RateContext(temperature_eV)
        lee_atomic = reactions[
            "e_Cl_ionization"].rate_coefficient.coefficient_si(context)
        nist_atomic = evaluated_atomic.coefficient_si(context)
        lee_molecular = sum(
            reactions[name].rate_coefficient.coefficient_si(context)
            for name in (
                "e_Cl2_nondissociative_ionization",
                "e_Cl2_dissociative_ionization",
            )
        )
        nist_molecular = evaluated_molecular.coefficient_si(context)
        rows.append({
            "electron_temperature_eV": temperature_eV,
            "nist_hayes_atomic_cl_rate_m3_s": nist_atomic,
            "lennon_atomic_cl_rate_m3_s": lee_atomic,
            "nist_to_lennon_atomic_ratio": nist_atomic / lee_atomic,
            "nist_cl2_total_ionization_rate_m3_s": nist_molecular,
            "lee_cl2_sum_ionization_rate_m3_s": lee_molecular,
            "nist_to_lee_cl2_total_ratio": nist_molecular / lee_molecular,
        })

    atomic_ratios = [
        row["nist_to_lennon_atomic_ratio"] for row in rows]
    molecular_ratios = [
        row["nist_to_lee_cl2_total_ratio"] for row in rows]
    summary: dict[str, object] = {
        "audit": (
            "evaluated chlorine ionization versus Lee--Lieberman rate fits"),
        "claim_class": (
            "no-fit evidence comparison; no retrospective pass threshold"),
        "coefficient_selection_target": None,
        "electron_temperature_eV": list(TEMPERATURES_EV),
        "lower_temperature_exclusion": (
            "Lennon Eq. 6 is valid only for Te > 1.296 eV; no extrapolated "
            "atomic-fit comparison is reported below 2 eV"),
        "atomic_cl_nist_to_lennon_ratio_range": [
            min(atomic_ratios), max(atomic_ratios)],
        "molecular_cl2_nist_to_lee_total_ratio_range": [
            min(molecular_ratios), max(molecular_ratios)],
        "molecular_total_ionization_evidence": (
            "NIST evaluated average of two incompatible measurement sets"),
        "molecular_species_branching": "unresolved",
        "branching_missing_quantity": (
            "electron-impact Cl2+ versus Cl+ production fraction"),
        "predictive_use": {
            "atomic_cl_positive_ion_source": "ready within quoted +/-14%",
            "molecular_total_positive_ion_source": (
                "ready as aggregate with no defensible scalar uncertainty"),
            "species_resolved_molecular_positive_ion_source": "not ready",
        },
    }
    return rows, summary


def write_results(
        rows: list[dict[str, float]], summary: dict[str, object]) -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    SUMMARY_JSON.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    table_rows = []
    for row in rows:
        table_rows.append(
            f"| {row['electron_temperature_eV']:.1f} | "
            f"{row['nist_to_lennon_atomic_ratio']:.4f} | "
            f"{row['nist_to_lee_cl2_total_ratio']:.4f} |"
        )
    atomic_range = summary["atomic_cl_nist_to_lennon_ratio_range"]
    molecular_range = summary[
        "molecular_cl2_nist_to_lee_total_ratio_range"]
    report = f"""# Evaluated chlorine ionization comparison

**Claim class: no-fit evidence comparison; no retrospective pass threshold**

| Maxwellian Te (eV) | NIST/Hayes Cl / Lennon | NIST Cl2 total / Lee sum |
|---:|---:|---:|
{chr(10).join(table_rows)}

Across 2--5 eV, the evaluated atomic-Cl rate is
`{atomic_range[0]:.4f}--{atomic_range[1]:.4f}` times Lennon's analytic fit.
The evaluated molecular total-ionization rate is
`{molecular_range[0]:.4f}--{molecular_range[1]:.4f}` times the sum of Lee's
non-dissociative and dissociative ionization fits.

These ratios are diagnostics, not a fitted agreement gate. The atomic
measurement carries a quoted ±14% absolute scale uncertainty. NIST's
molecular total is an evaluated average of two measurement sets whose
difference exceeds their combined quoted uncertainties, so assigning it one
smaller uncertainty would be false precision.

## Remaining failure

NIST states that partial electron-impact ionization measurements do not exist,
so the relative production of Cl2+ and Cl+ is unknown. Total positive-ion
production is now evidence-backed; species-resolved ion flux, sheath mass
transport, and feature delivery are not closed by this table. No feature
depth or reactor state selected any coefficient in this audit.
"""
    REPORT_MD.write_text(report, encoding="utf-8")


def main() -> None:
    rows, summary = run_audit()
    write_results(rows, summary)


if __name__ == "__main__":
    main()
