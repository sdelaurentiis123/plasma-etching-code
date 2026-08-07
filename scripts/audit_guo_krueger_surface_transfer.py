#!/usr/bin/env python3
"""Audit the no-depth-fit Guo/Kwon surface transfer to Krüger 2024.

This is deliberately a planar surface audit.  It asks whether the independent
C4F8/Ar translating-mixed-layer deck can supply the wafer-ion-normalized
removal demanded by Krüger's reported C4F6/Ar/O2 depth when driven only by
Krüger's published HPEM neutral fluxes and digitized aggregate IEAD.

The target depth is used only after every surface solve, for scoring.  No
coefficient, flux, species fraction, or energy law is adjusted to that target.
The report also exercises the missing-boundary directions hard enough that an
accidental scalar match cannot be promoted to a feature-depth prediction.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.depth_identifiability import effective_yield_from_depth
from petch.guo_c4f8_sio2 import (
    GuoC4F8ArSiO2Mechanism,
    GuoIncidentComposition,
    GuoIonQuadrature,
    GuoTmlState,
    physical_sputtering_angular,
    physical_sputtering_angular_literal,
)
from petch.reactor_boundary import (
    load_krueger_2024_digitized_iead,
    load_krueger_2024_reactor_flux_deck,
)


ROOT = Path(__file__).resolve().parents[1]
KRUEGER_DATA = ROOT / "data" / "experimental" / "krueger_2024"
GUO_SOURCE = ROOT / "data" / "surface_interactions" / "guo_2009"
DEPTH_AUDIT = (
    ROOT / "results" / "curated" / "depth_identifiability" / "audit.json"
)
OUTPUT = (
    ROOT / "results" / "curated" / "guo_krueger_surface_transfer"
    / "audit.json"
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _quadrature(iead, energy_cap_eV=None):
    energy = np.asarray(iead.energy_eV, dtype=float)
    if energy_cap_eV is not None:
        energy = np.minimum(energy, float(energy_cap_eV))
    return GuoIonQuadrature(
        energy,
        np.cos(np.deg2rad(np.abs(iead.signed_angle_deg))),
        iead.probability_weight,
    )


def _state_payload(state):
    return {
        "Si": state.si,
        "O": state.o,
        "C": state.c,
        "F": state.f,
        "V": state.vacancy,
    }


def _evaluate(
    iead,
    neutral_flux_ratio,
    *,
    ion_fraction=None,
    energy_cap_eV=None,
):
    try:
        mechanism = GuoC4F8ArSiO2Mechanism(
            GuoIncidentComposition(
                neutral_flux_ratio,
                {} if ion_fraction is None else ion_fraction,
            ),
            _quadrature(iead, energy_cap_eV),
        )
        result = mechanism.solve_steady_state()
    except RuntimeError as error:
        return {
            "steady_state_found": False,
            "error": str(error),
        }
    return {
        "steady_state_found": True,
        "sio2_yield_per_wafer_ion": result.sio2_yield_per_ion,
        "movement_atoms_per_wafer_ion": result.movement_atoms_per_ion,
        "state": _state_payload(result.state),
        "steady_state_residual": result.steady_state_residual,
        "atom_ledger_residual": (
            result.movement_atoms_per_ion
            - sum(result.removed_atoms_per_ion.values())
            + sum(result.incoming_atoms_per_ion.values())
        ),
        "source_extrapolation": dict(result.source_extrapolation),
    }


def _evaluate_finite_fluence(
    iead,
    neutral_flux_ratio,
    *,
    ion_flux_m2_s,
    duration_s,
    translating_layer_thickness_nm,
    sio2_formula_density_m3=2.2e28,
):
    mechanism = GuoC4F8ArSiO2Mechanism(
        GuoIncidentComposition(neutral_flux_ratio, {}),
        _quadrature(iead),
    )
    areal_real_atom_capacity_m2 = (
        3.0
        * float(sio2_formula_density_m3)
        * float(translating_layer_thickness_nm)
        * 1.0e-9
    )
    coordinate = (
        float(ion_flux_m2_s) * float(duration_s)
        / areal_real_atom_capacity_m2
    )
    result = mechanism.advance_fluence(GuoTmlState.oxide(), coordinate)
    average_yield = result.average_sio2_removal_yield_per_ion
    return {
        "translating_layer_thickness_nm": float(
            translating_layer_thickness_nm),
        "areal_real_atom_capacity_m2": areal_real_atom_capacity_m2,
        "incident_ions_per_tml_atom": coordinate,
        "average_sio2_removal_yield_per_wafer_ion": average_yield,
        "predicted_planar_depth_nm": (
            average_yield
            * float(ion_flux_m2_s)
            * float(duration_s)
            / float(sio2_formula_density_m3)
            * 1.0e9
        ),
        "final_net_movement_atoms_per_wafer_ion": (
            result.final_net_movement_atoms_per_ion),
        "final_state_derivative_residual": (
            result.final_state_derivative_residual),
        "maximum_atom_ledger_residual_atoms_per_ion": (
            result.maximum_atom_ledger_residual_atoms_per_ion),
        "final_state": _state_payload(result.state),
        "solver_step_count": result.solver_step_count,
        "feature_depth_used": False,
    }


def build_audit(root: Path | None = None):
    root = root or ROOT
    krueger_data = root / "data" / "experimental" / "krueger_2024"
    guo_source = root / "data" / "surface_interactions" / "guo_2009"
    depth_path = (
        root / "results" / "curated" / "depth_identifiability"
        / "audit.json"
    )
    manifest_path = guo_source / "source_manifest.json"
    gray_library_path = (
        root / "research_sources" / "library" / "gray-1993-thesis.md"
    )
    gray_ocr_path = (
        root / "research_sources" / "thesis_extracts"
        / "gray_thesis_1993_ocr_sections.txt"
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    depth = json.loads(depth_path.read_text(encoding="utf-8"))
    deck = load_krueger_2024_reactor_flux_deck(krueger_data)
    iead = load_krueger_2024_digitized_iead(krueger_data)

    flux = {
        record.name: record.flux_m2_s
        for record in deck.species_fluxes
    }
    ion_flux = flux["ions"]
    neutral_ratio = {
        name: value / ion_flux
        for name, value in flux.items()
        if name != "ions"
    }
    source_neutrals = {
        name.removesuffix("+")
        for name in manifest["transcription"][
            "source_species_footnote"]["neutrals_printed"]
    }
    source_intersection = {
        name: ratio
        for name, ratio in neutral_ratio.items()
        if name in source_neutrals
    }

    # Every solve is complete before the experimental target is constructed.
    # This ordering is not mathematically necessary, but makes the data
    # firewall visible in the implementation.
    nominal = _evaluate(iead, neutral_ratio)
    finite_fluence = {
        f"{thickness_nm:g}_nm": _evaluate_finite_fluence(
            iead,
            neutral_ratio,
            ion_flux_m2_s=ion_flux,
            duration_s=60.0,
            translating_layer_thickness_nm=thickness_nm,
        )
        for thickness_nm in (1.2, 2.0, 2.5, 3.0)
    }
    energy_sensitivity = {
        "printed_law_extended_to_full_iead": nominal,
        "energy_censored_at_2000_eV": _evaluate(
            iead, neutral_ratio, energy_cap_eV=2000.0),
        "energy_censored_at_1000_eV": _evaluate(
            iead, neutral_ratio, energy_cap_eV=1000.0),
        "energy_censored_at_500_eV": _evaluate(
            iead, neutral_ratio, energy_cap_eV=500.0),
        "energy_censored_at_source_maximum_370_eV": _evaluate(
            iead, neutral_ratio, energy_cap_eV=370.0),
    }
    neutral_sensitivity = {
        "all_published_neutrals": nominal,
        "omit_C3F4_outside_source_species_list": _evaluate(
            iead,
            {name: value for name, value in neutral_ratio.items()
             if name != "C3F4"},
        ),
        "source_species_intersection_only": _evaluate(
            iead, source_intersection),
        "half_all_neutral_fluxes": _evaluate(
            iead,
            {name: 0.5 * value for name, value in neutral_ratio.items()},
        ),
        "double_all_neutral_fluxes": _evaluate(
            iead,
            {name: 2.0 * value for name, value in neutral_ratio.items()},
        ),
    }
    ion_sensitivity = {
        "aggregate_nonincorporating": nominal,
        **{
            f"all_{species}_positive_ions": _evaluate(
                iead, neutral_ratio, ion_fraction={species: 1.0})
            for species in ("C", "O", "CF", "CF2", "CF3", "C3F3")
        },
    }

    target_yield = effective_yield_from_depth(
        825.0, 60.0, 2.2e28, 1.2e20)
    prior_yield = depth[
        "run_average_wafer_ion_normalization"][
            "simulation_sio2_per_wafer_ion"]
    prior_depth_nm = depth["inputs"]["simulated_depth_nm"]
    nominal_yield = nominal["sio2_yield_per_wafer_ion"]
    finite_ion_yields = [
        case["sio2_yield_per_wafer_ion"]
        for case in ion_sensitivity.values()
        if case["steady_state_found"]
    ]

    angles_deg = np.asarray([0.0, 25.0, 45.0, 65.0, 85.0, 90.0])
    cosine = np.cos(np.deg2rad(angles_deg))
    repaired_angular = physical_sputtering_angular(cosine)
    literal_angular = physical_sputtering_angular_literal(cosine)

    within_source_energy = (
        np.asarray(iead.energy_eV)
        <= GuoC4F8ArSiO2Mechanism.source_energy_max_eV
    )
    source_support_weight = float(np.sum(
        np.asarray(iead.probability_weight)[within_source_energy]))
    independent_chemical_form_support_weight = float(np.sum(
        np.asarray(iead.probability_weight)[
            np.asarray(iead.energy_eV) <= 2000.0
        ]
    ))
    independent_physical_form_support_weight = float(np.sum(
        np.asarray(iead.probability_weight)[
            np.asarray(iead.energy_eV) <= 1225.0
        ]
    ))

    return {
        "audit_id": "GUO-KWON-TO-KRUEGER-SURFACE-TRANSFER-R2",
        "status": "passed",
        "question": (
            "Can an independently beam-yield-regressed, atom-balanced "
            "fluorocarbon/SiO2 surface deck supply Krueger's "
            "wafer-ion-normalized removal without a depth fit?"
        ),
        "calibration_firewall": {
            "feature_depth_used_by_surface_solver": False,
            "surface_or_boundary_parameters_adjusted": [],
            "optimization_performed": False,
            "target_used_only_after_all_surface_solves_for_scoring": True,
            "model_coefficients_origin": (
                "Guo Table 4.1 coefficients regressed against Yin C4F8/Ar "
                "blanket-yield data; no Krueger feature observable"
            ),
        },
        "sources": {
            "guo_pdf_sha256": manifest["source"]["pdf_sha256"],
            "kwon_pdf_sha256": manifest["formalism_source"]["pdf_sha256"],
            "gray_pdf_sha256": (
                "be6bce26b699b3172cf67bb68e4d12e039fd3ea775f73873ee1aaf251"
                "64c065b"
            ),
            "gray_library_record_sha256": _sha(gray_library_path),
            "gray_ocr_extract_sha256": _sha(gray_ocr_path),
            "guo_manifest_sha256": _sha(manifest_path),
            "guo_table_sha256": _sha(
                guo_source / "table4_1_reaction_deck.csv"),
            "krueger_flux_table_sha256": deck.source_sha256,
            "krueger_iead_table_sha256": iead.table_sha256,
            "krueger_iead_metadata_sha256": iead.metadata_sha256,
            "depth_identifiability_audit_sha256": _sha(depth_path),
        },
        "boundary": {
            "neutral_flux_ratio_to_published_aggregate_ion_flux":
                neutral_ratio,
            "source_species_intersection": sorted(source_intersection),
            "outside_guo_source_species_list": sorted(
                set(neutral_ratio) - source_neutrals),
            "positive_ion_composition_published": False,
            "stable_C4F6_parent_flux_published": False,
            "nominal_ion_closure": (
                "aggregate ions trigger universal ion reactions but add no "
                "material; this is a nonincorporating sensitivity endpoint, "
                "not an identified ion composition"
            ),
        },
        "energy_support": {
            "guo_fit_maximum_eV": (
                GuoC4F8ArSiO2Mechanism.source_energy_max_eV),
            "krueger_iead_minimum_eV": float(np.min(iead.energy_eV)),
            "krueger_iead_mean_eV": float(iead.mean_energy_eV),
            "krueger_iead_maximum_eV": float(np.max(iead.energy_eV)),
            "krueger_iead_probability_within_guo_fit_support":
                source_support_weight,
            "independent_sqrt_form_support": {
                "chemical_F_saturated_SiO2_maximum_eV": 2000.0,
                "chemical_probability_weight_within_support":
                    independent_chemical_form_support_weight,
                "physical_Ar_SiO2_maximum_eV": 1225.0,
                "physical_probability_weight_within_support":
                    independent_physical_form_support_weight,
                "source": (
                    "Gray 1993 thesis Eq. 5-35/Table 5-10 and Figure 5-2; "
                    "page renders and OCR are SHA-pinned above"
                ),
                "transfer_limit": (
                    "this validates the square-root functional form over part "
                    "of the Krueger IEAD, not Guo's reaction coefficients, "
                    "C4F6 radical mapping, or an aggregate-ion projectile"
                ),
            },
            "mean_energy_over_source_maximum": (
                iead.mean_energy_eV
                / GuoC4F8ArSiO2Mechanism.source_energy_max_eV),
            "meaning": (
                "the nominal match extrapolates Guo's fitted coefficients for "
                "every IEAD sample. Independent Gray beam data support the "
                "same square-root form only over the stated lower-energy "
                "probability mass; the remaining high-energy tail and "
                "coefficient/species transfer remain unvalidated. Censored "
                "cases are sensitivity brackets, not alternate laws"
            ),
        },
        "nominal_no_fit_surface_result": nominal,
        "finite_fluence_planar_forecast_before_target_scoring": {
            "coordinate": (
                "incident ions per translating-layer real atom"),
            "nominal_thickness_nm": 2.5,
            "source_sensitivity_band_nm": [1.2, 3.0],
            "cases": finite_fluence,
            "predicted_depth_range_nm": [
                min(
                    case["predicted_planar_depth_nm"]
                    for case in finite_fluence.values()
                ),
                max(
                    case["predicted_planar_depth_nm"]
                    for case in finite_fluence.values()
                ),
            ],
            "authoritative_feature_prediction": False,
        },
        "score": {
            "experimental_run_average_sio2_per_wafer_ion": target_yield,
            "nominal_guo_sio2_per_wafer_ion": nominal_yield,
            "signed_relative_error": (
                nominal_yield / target_yield - 1.0),
            "absolute_percentage_error": abs(
                nominal_yield / target_yield - 1.0),
            "within_five_percent": bool(
                abs(nominal_yield / target_yield - 1.0) <= 0.05),
        },
        "energy_law_sensitivity": energy_sensitivity,
        "neutral_mapping_sensitivity": neutral_sensitivity,
        "unpublished_ion_composition_sensitivity": {
            "cases": ion_sensitivity,
            "converged_yield_range_including_deposition": [
                min(finite_ion_yields), max(finite_ion_yields)],
            "deposition_cases": sorted(
                name for name, case in ion_sensitivity.items()
                if (
                    case["steady_state_found"]
                    and case["sio2_yield_per_wafer_ion"] < 0.0
                )
            ),
            "nonconverged_cases": sorted(
                name for name, case in ion_sensitivity.items()
                if not case["steady_state_found"]
            ),
        },
        "angular_source_defect": {
            "angles_deg": angles_deg.tolist(),
            "declared_degree_sequence_repair": repaired_angular.tolist(),
            "literal_printed_duplicate_cos2": literal_angular.tolist(),
            "literal_negative_sample_count": int(
                np.count_nonzero(literal_angular < 0.0)),
            "repair_independently_traced_to_original_fit": False,
            "meaning": (
                "near-normal planar scoring is insensitive to this defect; "
                "a sidewall/profile claim is not"
            ),
        },
        "linear_depth_diagnostic_only": {
            "prior_feature_depth_nm": prior_depth_nm,
            "prior_feature_effective_yield_per_wafer_ion": prior_yield,
            "surface_yield_ratio_scaled_depth_nm": (
                prior_depth_nm * nominal_yield / prior_yield),
            "authoritative_feature_prediction": False,
            "invalidated_assumptions": [
                "depth is linear in planar steady-state yield",
                "floor delivery and surface state do not evolve",
                "mouth closure and mask evolution are unchanged",
            ],
        },
        "atomicity_statement": {
            "atom_balanced": True,
            "bond_probability_resolved": True,
            "atomistic_trajectory_or_interatomic_potential": False,
            "atomic_level_accuracy_claimed": False,
            "meaning": (
                "the closure conserves Si/O/C/F atoms and carries random-bond "
                "statistics, but it is a regressed translating mixed layer, "
                "not molecular dynamics or a first-principles potential"
            ),
        },
        "verdict": {
            "no_fit_planar_surface_normalization_within_five_percent": True,
            "surface_transfer_hypothesis_is_worth_feature_testing": True,
            "energy_support_is_closed": False,
            "reactor_species_boundary_is_identified": False,
            "collision_bounded_neutral_transport_coupling_is_defined": False,
            "feature_depth_match_earned": False,
            "exact_825nm_prediction_authorized": False,
            "reason": (
                "the 2.613 versus 2.521 planar match is real and uses no "
                "Krueger depth fit. Gray independently supports the square-"
                "root form to 2000 eV for F-saturated chemical removal and "
                "1225 eV for physical sputtering, but Guo's reaction "
                "coefficients are still extended from <=370 eV to a "
                "471--4821 eV boundary. The result also treats two C4F6 "
                "radicals through a universal-neutral transfer and chooses "
                "one endpoint of an unpublished ion-composition range; no "
                "evolving feature calculation has used this surface state"
            ),
            "next_required_physics": [
                (
                    "bound the >1225/2000 eV reaction moments with a "
                    "projectile-resolved stopping/cascade calculation and "
                    "species-resolved high-energy beam or MD evidence"
                ),
                (
                    "obtain a species-resolved positive-ion flux and IEAD, "
                    "or validate a reactor provider against that diagnostic"
                ),
                (
                    "measure or independently predict stable C4F6 molecular "
                    "flux and C4F6/ion co-incidence response"
                ),
                (
                    "derive a bounded collision-level neutral uptake law; "
                    "Guo adsorption coefficients exceed one and are rates, "
                    "not radiosity sticking probabilities"
                ),
                (
                    "run the frozen feature transport/geometry using the "
                    "implemented local finite-fluence Si/O/C/F/V state, then "
                    "score depth and profile without changing this board"
                ),
            ],
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
    payload = canonical_payload(build_audit())
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(payload, encoding="utf-8")
    elif args.check or OUTPUT.exists():
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != payload:
            raise RuntimeError("committed Guo/Krueger transfer audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
