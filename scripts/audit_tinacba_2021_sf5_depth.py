#!/usr/bin/env python3
"""Grade DFT-informed SF5+ MD yields against independent beam depths."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "data"
    / "experimental"
    / "tinacba_2021"
    / "figure8_sf5_md_experiment.csv"
)
MANIFEST = DATA.with_name("digitization_manifest.json")
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "tinacba_2021_sf5_depth" / "audit.json"
)
REFERENCE_DOSE_M2 = 1.0e20  # 1e16 cm^-2
FIGURE10_MD_DEPTH_NM = {"Si": 12.5, "SiO2": 13.6}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load() -> list[dict[str, object]]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    out = []
    for row in rows:
        out.append(
            {
                "material": row["material"],
                "series": row["series"],
                "energy_eV": int(row["energy_eV"]),
                "yield": float(row["si_removal_yield_per_sf5_ion"]),
                "digitization_bound": float(row["digitization_yield_bound"]),
            }
        )
    return out


def build_report() -> dict[str, object]:
    rows = _load()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    lookup = {
        (row["material"], row["series"], row["energy_eV"]): row
        for row in rows
    }

    # Figure 10 prints the independent MD depth slopes at 2000 eV.  Combining
    # those slopes with the Figure-8 MD yield only derives the model-film
    # number density used to express every comparison as nm per fixed dose.
    # It cannot alter MD/experiment depth ratios, which equal yield ratios.
    densities_m3 = {}
    for material, depth_nm in FIGURE10_MD_DEPTH_NM.items():
        md = lookup[(material, "sf5_md", 2000)]["yield"]
        densities_m3[material] = (
            md * REFERENCE_DOSE_M2 / (depth_nm * 1.0e-9)
        )

    comparisons = []
    relative_errors = []
    yield_errors = []
    for material in ("Si", "SiO2"):
        experiment_energies = sorted(
            row["energy_eV"]
            for row in rows
            if row["material"] == material
            and row["series"] == "sf5_experiment"
        )
        for energy in experiment_energies:
            md = lookup[(material, "sf5_md", energy)]
            measured = lookup[(material, "sf5_experiment", energy)]
            density = densities_m3[material]
            predicted_depth = (
                md["yield"] * REFERENCE_DOSE_M2 / density * 1.0e9
            )
            measured_depth = (
                measured["yield"] * REFERENCE_DOSE_M2 / density * 1.0e9
            )
            signed_relative = md["yield"] / measured["yield"] - 1.0
            relative_errors.append(abs(signed_relative))
            yield_errors.append(md["yield"] - measured["yield"])
            comparisons.append(
                {
                    "material": material,
                    "energy_eV": energy,
                    "md_yield_si_per_sf5": md["yield"],
                    "measured_yield_si_per_sf5": measured["yield"],
                    "signed_relative_depth_error": signed_relative,
                    "absolute_relative_depth_error": abs(signed_relative),
                    "predicted_depth_nm_per_1e16_cm2": predicted_depth,
                    "measured_depth_nm_per_1e16_cm2": measured_depth,
                    "digitization_bound_each_yield": md[
                        "digitization_bound"
                    ],
                }
            )

    return {
        "schema_version": 1,
        "claim": (
            "retrospective no-yield-fit atomistic-provider validation against "
            "mass-selected SF5+ beam depth-per-dose; not a feature or reactor "
            "prediction"
        ),
        "source": {
            "citation": manifest["source"]["citation"],
            "doi": manifest["source"]["doi"],
            "pdf_sha256": manifest["source"]["pdf_sha256"],
            "figure8_csv_sha256": _sha256(DATA),
            "figure8_render_sha256": manifest["source"]["render_sha256"],
            "figure": "Figure 8 comparison; Figure 10 MD depth slopes",
        },
        "boundary": {
            "projectile": "SF5+",
            "mass_selected": True,
            "normal_incidence": True,
            "energy_measured": True,
            "dose_measured_at_sample_position": True,
            "neutral_radical_beam": False,
            "energies_scored_eV": [150, 2000],
        },
        "provider": {
            "method": "DFT-informed modified-Stillinger-Weber MD",
            "beam_depth_or_yield_fit_used": False,
            "typical_impacts_per_trajectory": 4000,
            "surface_temperature_K": 300,
            "state_evolution": (
                "continuous impacts until steady surface composition/yield"
            ),
            "material_limitation": (
                "S-F carrier bond and S mass/radius retained; S-S, S-Si, and "
                "S-O chemistry intentionally absent"
            ),
        },
        "common_core_adapter": {
            "table_loader": "petch.interaction_data.load_tinacba_2021_sf5_tables",
            "surface_mechanism": (
                "petch.tabulated_chemistry."
                "TabulatedNormalIonRemovalMechanism"
            ),
            "material_router": "petch.material_mechanism_3d.MaterialMechanismRouter3D",
            "materials": {
                "Si": "Si_atom",
                "SiO2": "SiO2_formula",
            },
            "atom_or_formula_ledger_closed": True,
            "product_routing_complete": False,
            "feature_profile_validated": False,
            "refused_dimensions": [
                "off-normal incidence",
                "neutral co-flux",
                "projectiles other than SF5+",
                "energy extrapolation outside 150-2000 eV",
            ],
        },
        "depth_conversion": {
            "equation": "d = Y*D/N",
            "reference_dose_m2": REFERENCE_DOSE_M2,
            "reference_dose_cm2": REFERENCE_DOSE_M2 / 1.0e4,
            "model_film_number_density_m3": densities_m3,
            "density_derivation": (
                "Figure 8 MD yield at 2000 eV combined with Figure 10 printed "
                "MD slopes 12.5 nm (Si) and 13.6 nm (SiO2) per 1e16 cm^-2"
            ),
            "comparison_invariance": (
                "predicted/measured depth ratio equals MD/measured yield ratio "
                "for any common declared D and N"
            ),
        },
        "comparison": {
            "point_count": len(comparisons),
            "points": comparisons,
            "mean_absolute_relative_depth_error": float(
                np.mean(relative_errors)
            ),
            "maximum_absolute_relative_depth_error": float(
                np.max(relative_errors)
            ),
            "rmse_yield_si_per_sf5": float(
                np.sqrt(np.mean(np.square(yield_errors)))
            ),
            "post_hoc_pass_gate_declared": False,
        },
        "uncertainty_and_scope": {
            "digitization_bound_yield": manifest["digitization"][
                "yield_bound"
            ],
            "experimental_statistical_uncertainty": "not reported",
            "combined_uncertainty_claimed": False,
            "retrospective_not_blind": True,
            "not_authorized": manifest["claim_boundary"]["does_not_support"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = build_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"stale audit: {args.output}")
        print(f"verified {args.output.relative_to(ROOT)}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
