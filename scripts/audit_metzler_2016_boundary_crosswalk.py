#!/usr/bin/env python3
"""Crosswalk Metzler's absolute-depth and per-ion figures.

Figures 6.6 and 6.9(b) repeat the same nominal 5 A C4F8/Ar ALE conditions:
25 or 30 eV and 20, 40, or 60 s of ion exposure.  The former reports depth
per cycle; the latter reports substrate units removed per incident Ar ion.
Their ratio therefore contains the *author-normalized total ion fluence*:

    Phi_i = depth * substrate number density / (substrate units / ion).

This is a useful recovery of a total-fluence boundary that was previously
missed.  It is deliberately not called an independent current measurement:
Figure 6.9 is itself a source-derived normalization, the two figures may
represent different experimental repeats, SiO2 number density is declared
rather than reported for the PECVD film, and neither an IEDF nor an IADF is
published.  The audit propagates digitization intervals and retains every
plotted replicate so those limitations remain visible.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import fmean, stdev


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "metzler_2016"
DEPTH_SOURCE = DATA / "figures6_5_6_6_cyclic_depth.csv"
YIELD_SOURCE = DATA / "figure6_9_cycle_averaged_yield.csv"
OUTPUT = (
    ROOT / "results" / "curated" / "metzler_2016_boundary_crosswalk"
    / "audit.json"
)

ELEMENTARY_CHARGE_C = 1.602176634e-19
SILICON_LATTICE_CONSTANT_M = 5.431e-10
NUMBER_DENSITY_M3 = {
    # The petch oxide convention: formula units, consistent with the Figure
    # 6.9 ordinate "SiO2/Ar+".  This is a declared density assumption, not a
    # density measurement of Metzler's PECVD stack.
    "SiO2": 2.2e28,
    # Diamond cubic silicon has eight atoms in the conventional cubic cell.
    "Si": 8.0 / SILICON_LATTICE_CONSTANT_M ** 3,
}
PDF_SHA256 = (
    "ea5701d0bcf67b56403625253f7bb619c0e3e3a5e0a9cfd2ff6f1e435fa90f62"
)


def _rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _key(row):
    return (
        int(float(row["energy_eV"])),
        row["substrate"],
        float(row["etch_step_s"]),
        int(row["replicate"]),
    )


def _depth_index():
    rows = [
        row for row in _rows(DEPTH_SOURCE)
        if row["panel"] == "6.6"
    ]
    return {_key(row): row for row in rows}


def _interval(depth_A, depth_uncertainty_A, response, response_uncertainty,
              duration_s, number_density_m3):
    """Return conservative independent-placement bounds on flux."""
    low = (
        max(depth_A - depth_uncertainty_A, 0.0) * 1e-10
        * number_density_m3
        / (response + response_uncertainty)
        / duration_s
    )
    high = (
        (depth_A + depth_uncertainty_A) * 1e-10
        * number_density_m3
        / max(response - response_uncertainty, 1e-300)
        / duration_s
    )
    return low, high


def _point_rows():
    depth = _depth_index()
    rows = []
    for source in _rows(YIELD_SOURCE):
        if source["panel"] != "6.9b":
            continue
        key = _key(source)
        depth_row = depth.get(key)
        if depth_row is None:
            raise ValueError(f"no Figure 6.6 crosswalk row for {key}")
        material = source["substrate"]
        number_density = NUMBER_DENSITY_M3[material]
        duration_s = float(source["etch_step_s"])
        depth_A = float(depth_row["etch_depth_A_per_cycle"])
        depth_uncertainty_A = float(
            depth_row["digitization_uncertainty_depth_A"]
        )
        response = float(source["substrate_units_per_incident_ion"])
        response_uncertainty = float(source["digitization_uncertainty_ratio"])
        fluence = depth_A * 1e-10 * number_density / response
        flux = fluence / duration_s
        low, high = _interval(
            depth_A,
            depth_uncertainty_A,
            response,
            response_uncertainty,
            duration_s,
            number_density,
        )
        rows.append({
            "energy_eV": int(float(source["energy_eV"])),
            "substrate": material,
            "etch_step_s": duration_s,
            "replicate": int(source["replicate"]),
            "depth_observation_id": depth_row["observation_id"],
            "yield_observation_id": source["observation_id"],
            "etch_depth_A_per_cycle": depth_A,
            "substrate_units_per_incident_ion": response,
            "declared_substrate_number_density_m3": number_density,
            "author_normalized_ion_fluence_m2": fluence,
            "author_normalized_average_ion_flux_m2_s": flux,
            "digitization_only_flux_interval_m2_s": [low, high],
            "equivalent_singly_charged_current_density_mA_cm2":
                flux * ELEMENTARY_CHARGE_C * 0.1,
        })
    return rows


def _group_metrics(rows):
    groups = {}
    for row in rows:
        key = (row["substrate"], row["energy_eV"])
        groups.setdefault(key, []).append(row)
    output = []
    for (substrate, energy), members in sorted(groups.items()):
        values = [
            row["author_normalized_average_ion_flux_m2_s"]
            for row in members
        ]
        mean = fmean(values)
        output.append({
            "substrate": substrate,
            "energy_eV": energy,
            "points": len(values),
            "mean_author_normalized_flux_m2_s": mean,
            "coefficient_of_variation": (
                stdev(values) / mean if len(values) > 1 else 0.0
            ),
            "range_over_mean": (max(values) - min(values)) / mean,
        })
    return output


def _material_metrics(rows):
    output = {}
    for substrate in ("SiO2", "Si"):
        # The duplicate 30 eV/40 s Si marker is a plotted replicate, not an
        # extra process duration.  Retain it in point rows but average the two
        # replicates before forming the material-level flux.
        conditions = {}
        for row in rows:
            if row["substrate"] != substrate:
                continue
            key = (row["energy_eV"], row["etch_step_s"])
            conditions.setdefault(key, []).append(
                row["author_normalized_average_ion_flux_m2_s"]
            )
        condition_values = [
            fmean(values) for _, values in sorted(conditions.items())
        ]
        mean = fmean(condition_values)
        output[substrate] = {
            "independent_nominal_conditions": len(condition_values),
            "mean_author_normalized_flux_m2_s": mean,
            "coefficient_of_variation_across_conditions":
                stdev(condition_values) / mean,
            "equivalent_singly_charged_current_density_mA_cm2":
                mean * ELEMENTARY_CHARGE_C * 0.1,
        }
    return output


def build_report():
    points = _point_rows()
    group_metrics = _group_metrics(points)
    material_metrics = _material_metrics(points)
    sio2_groups = [
        row for row in group_metrics if row["substrate"] == "SiO2"
    ]
    return {
        "audit_id": "METZLER-2016-BOUNDARY-CROSSWALK-R1",
        "operation": (
            "algebraic crosswalk of two already published observables; "
            "no surface parameter, flux, density, energy distribution, or "
            "etch depth was optimized"
        ),
        "source": {
            "citation": (
                "D. Metzler, High Precision Plasma Etch for Pattern "
                "Transfer: Towards Fluorocarbon Based Atomic Layer Etching, "
                "PhD thesis, University of Maryland (2016), Figures 6.6 "
                "and 6.9(b)"
            ),
            "source_pdf_sha256": PDF_SHA256,
            "depth_path": str(DEPTH_SOURCE.relative_to(ROOT)),
            "yield_path": str(YIELD_SOURCE.relative_to(ROOT)),
            "full_resolution_visual_receipt": {
                "poppler_resolution_dpi": 240,
                "figure_6_6_pdf_page": 142,
                "figure_6_6_png_sha256": (
                    "52845fe4c47ca9ba8f270af1a0de4fc71a04cc489849a51044cf1c5cc4ebc5fa"
                ),
                "figure_6_9_pdf_page": 147,
                "figure_6_9_png_sha256": (
                    "86e665a6d1fc417ac9ae6ff8b7b900e36f69871a5b84d0467756b7cf2d5824e2"
                ),
                "condition_mapping_visually_confirmed": (
                    "both panels state 5 A FC deposition, 25/30 eV, and "
                    "20-60 s etch-step length"
                ),
            },
        },
        "equation": (
            "author_normalized_ion_fluence = "
            "etch_depth * substrate_number_density / "
            "substrate_units_per_incident_ion"
        ),
        "density_conventions": {
            "SiO2_formula_units_m3": NUMBER_DENSITY_M3["SiO2"],
            "Si_atoms_m3": NUMBER_DENSITY_M3["Si"],
            "Si_derivation": (
                "8 atoms / (5.431e-10 m)^3 for diamond-cubic silicon"
            ),
            "SiO2_limitation": (
                "2.2e28 formula units/m3 is the declared petch oxide "
                "density; Metzler did not report a density measurement for "
                "the PECVD stack"
            ),
        },
        "points": points,
        "group_metrics": group_metrics,
        "material_metrics": material_metrics,
        "gates": {
            "six_SiO2_conditions_recover_duration_linear_fluence": {
                "threshold_maximum_within_energy_range_over_mean": 0.02,
                "observed": max(
                    row["range_over_mean"] for row in sio2_groups
                ),
                "passes": max(
                    row["range_over_mean"] for row in sio2_groups
                ) < 0.02,
            },
        },
        "claim_boundary": {
            "established": [
                (
                    "the paired source observables recover a density-"
                    "conditioned, author-normalized total ion fluence for "
                    "each repeated 5 A condition"
                ),
                (
                    "the six SiO2 crosswalks scale linearly with 20/40/60 s "
                    "to within two percent inside each nominal energy"
                ),
            ],
            "not_established": [
                "an independently measured ion current or flux",
                "a universal flux shared by the separately processed Si and SiO2 samples",
                "ion species fractions",
                "an ion-energy or ion-angle distribution",
                "an absolute deposited-film C/F inventory independent of the source normalization",
                "permission to reuse a target yield as a validation prediction",
            ],
            "valid_use": (
                "boundary reconstruction or calibration evidence with the "
                "source-derived label; a prediction gate must use a different "
                "surface/depth observable than the ratio used here"
            ),
        },
    }


def canonical_payload(report):
    return json.dumps(report, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = canonical_payload(build_report())
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(payload, encoding="utf-8")
        print(f"wrote {OUTPUT.relative_to(ROOT)}")
    elif args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"stale Metzler boundary crosswalk: {OUTPUT}")
        print(f"verified {OUTPUT.relative_to(ROOT)}")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
