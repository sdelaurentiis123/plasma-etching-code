#!/usr/bin/env python3
"""Bound what Krueger's aggregate IEAD can identify about ion composition.

The published observable is one *combined* positive-ion IEAD at each power.
Without species-resolved kernels its forward operator has identical columns:
every composition on the simplex reproduces the same published aggregate.
This audit records that structural rank result, then propagates measured
C4F6/Ar composition bands from Benck et al. as explicitly cross-reactor
sensitivities through the no-depth-fit Guo/Kwon surface model.

No ion fraction is selected using Krueger's 825 nm endpoint.  The endpoint is
constructed only after all surface solves and is used solely for scoring.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.depth_identifiability import effective_yield_from_depth
from petch.guo_c4f8_sio2 import (
    GuoC4F8ArSiO2Mechanism,
    GuoIncidentComposition,
    GuoIonQuadrature,
)
from petch.reactor_boundary import (
    load_krueger_2024_digitized_iead,
    load_krueger_2024_reactor_flux_deck,
    load_krueger_2024_transfer_boundary_data,
)


ROOT = Path(__file__).resolve().parents[1]
KRUEGER = ROOT / "data" / "experimental" / "krueger_2024"
BENCK = (
    ROOT / "data" / "experimental" / "benck_2003_c4f6"
    / "figure9_mass_resolved_ion_current.csv"
)
BENCK_MANIFEST = BENCK.with_name("digitization_manifest.json")
KRUEGER_THESIS_TEXT = (
    ROOT / "research_sources" / "thesis_extracts"
    / "krueger_thesis_2024.txt"
)
OUTPUT = (
    ROOT / "results" / "curated" / "krueger_ion_mixture_inverse"
    / "audit.json"
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _quadrature(iead) -> GuoIonQuadrature:
    return GuoIonQuadrature(
        iead.energy_eV,
        np.cos(np.deg2rad(np.abs(iead.signed_angle_deg))),
        iead.probability_weight,
    )


def _benck_composition(feed_percent: int) -> dict[str, object]:
    with BENCK.open(newline="", encoding="utf-8") as stream:
        selected = {
            row["species"]: float(row["ion_current_density_mA_cm2"])
            for row in csv.DictReader(stream)
            if int(row["c4f6_feed_percent"]) == int(feed_percent)
        }
    required = {"total_positive_ion_current", "Ar+", "CF+", "CF2+", "CF3+"}
    if set(selected) != required:
        raise ValueError(f"incomplete Benck composition at {feed_percent}% C4F6")
    total = selected["total_positive_ion_current"]
    fraction = {
        name: selected[name] / total
        for name in ("Ar+", "CF+", "CF2+", "CF3+")
    }
    unresolved = 1.0 - sum(fraction.values())
    if unresolved < 0.0:
        raise ValueError("Benck resolved ion currents exceed total current")
    return {
        "c4f6_feed_percent": int(feed_percent),
        "ar_feed_percent": 100 - int(feed_percent),
        "resolved_fraction_of_total_positive_current": fraction,
        "unresolved_positive_ion_fraction": unresolved,
        "resolved_plus_unresolved_closure": sum(fraction.values()) + unresolved,
    }


def _surface_case(iead, neutral_ratio, ion_fraction) -> dict[str, object]:
    mechanism = GuoC4F8ArSiO2Mechanism(
        GuoIncidentComposition(neutral_ratio, ion_fraction),
        _quadrature(iead),
    )
    solved = mechanism.solve_steady_state()
    return {
        "ion_fraction": dict(ion_fraction),
        "inert_remainder_fraction": (
            1.0 - sum(float(value) for value in ion_fraction.values())
        ),
        "sio2_yield_per_wafer_ion": solved.sio2_yield_per_ion,
        "movement_atoms_per_wafer_ion": solved.movement_atoms_per_ion,
        "steady_state_residual": solved.steady_state_residual,
        "atom_ledger_residual": (
            solved.movement_atoms_per_ion
            - sum(solved.removed_atoms_per_ion.values())
            + sum(solved.incoming_atoms_per_ion.values())
        ),
        "surface_state": dict(solved.state.fractions()),
        "feature_depth_used": False,
    }


def build_audit() -> dict[str, object]:
    base_iead = load_krueger_2024_digitized_iead(KRUEGER)
    transfer = load_krueger_2024_transfer_boundary_data(KRUEGER)
    deck = load_krueger_2024_reactor_flux_deck(KRUEGER)
    metadata = json.loads(
        (KRUEGER / "digitized_figure16_metadata.json").read_text(
            encoding="utf-8"))
    flux = {item.name: item.flux_m2_s for item in deck.species_fluxes}
    neutral_ratio = {
        name: value / flux["ions"]
        for name, value in flux.items() if name != "ions"
    }

    candidate_species = (
        "Ar+", "O+", "CF+", "CF2+", "CF3+",
        "C2F3+", "C2F4+", "C2F5+", "C3F5+", "C3F6+", "C3F7+",
    )
    # The source publishes only F=sum_s w_s K_s and no K_s.  The only design
    # matrix authorized by the published data therefore repeats F in every
    # column.  Its composition-sensitive contrast has exact rank zero.
    published_operator = np.ones((4, len(candidate_species)), dtype=float)
    aggregate_rank = int(np.linalg.matrix_rank(published_operator))
    simplex_contrast = published_operator[:, 1:] - published_operator[:, [0]]
    contrast_rank = int(np.linalg.matrix_rank(simplex_contrast))

    benck_boards = {
        str(feed): _benck_composition(feed) for feed in (50, 75)
    }
    surface_solves: dict[str, dict[str, dict[str, object]]] = {}
    for feed, board in benck_boards.items():
        resolved = board["resolved_fraction_of_total_positive_current"]
        known_reactive = {
            species.removesuffix("+"): resolved[species]
            for species in ("CF+", "CF2+", "CF3+")
        }
        unresolved = float(board["unresolved_positive_ion_fraction"])
        cases = {}
        for closure in ("inert", "CF", "CF2", "CF3", "C3F3"):
            composition = dict(known_reactive)
            if closure != "inert":
                composition[closure] = composition.get(closure, 0.0) + unresolved
            cases[f"unresolved_as_{closure}"] = _surface_case(
                base_iead, neutral_ratio, composition)
        surface_solves[feed] = cases

    # Data firewall: construct the Krueger endpoint only after every mixture
    # and surface state above has been fixed and solved.
    target_yield = effective_yield_from_depth(
        825.0, 60.0, 2.2e28, flux["ions"])
    for cases in surface_solves.values():
        for case in cases.values():
            case["score_only_after_solve"] = {
                "target_sio2_per_wafer_ion": target_yield,
                "signed_relative_error": (
                    case["sio2_yield_per_wafer_ion"] / target_yield - 1.0
                ),
            }

    energy_bin_eV = float(
        metadata["iead_digitization"]["joint_bins"]["energy_eV"])
    same_author_mass_shift_eV = 60.0
    return {
        "audit_id": "KRUEGER-ION-MIXTURE-INVERSE-R1",
        "status": "structurally_nonidentifiable_from_published_iead",
        "question": (
            "Can Krueger's combined positive-ion IEAD and power sweep identify "
            "the species mixture without using feature depth?"
        ),
        "calibration_firewall": {
            "feature_depth_used_to_choose_mixture": False,
            "optimization_against_825_nm_performed": False,
            "target_constructed_only_after_all_surface_solves": True,
        },
        "published_inverse_problem": {
            "candidate_positive_ions_from_thesis_surface_deck": list(
                candidate_species),
            "published_channels_per_power": 1,
            "power_nodes_kw": sorted(
                transfer.iead_by_low_frequency_power_kw),
            "species_resolved_iead_kernels_published": False,
            "aggregate_operator_rank": aggregate_rank,
            "composition_contrast_rank": contrast_rank,
            "composition_nullity": len(candidate_species) - 1 - contrast_rank,
            "meaning": (
                "Every candidate column is the same combined IEAD under the "
                "published evidence. The power sweep changes that aggregate "
                "with sheath power but does not label any species column."
            ),
        },
        "figure_resolution_gate": {
            "declared_energy_compression_bin_eV": energy_bin_eV,
            "same_author_different_condition_Oplus_to_CF3plus_mean_shift_eV": (
                same_author_mass_shift_eV),
            "mass_shift_per_digitization_bin": (
                same_author_mass_shift_eV / energy_bin_eV),
            "same_author_condition_is_not_krueger_base_process": True,
            "verdict": (
                "The apparent approximately 250 eV ladder in the committed "
                "digitization is the binning operator, not a species comb."
            ),
        },
        "cross_reactor_measured_composition_sensitivities": {
            "source": (
                "Benck, Goyette and Wang 2003, C4F6/Ar ICP, 10 mTorr, "
                "200 W, no O2, grounded sampling surface"
            ),
            "transfer_to_krueger_is_predictive": False,
            "boards": benck_boards,
            "surface_propagation": surface_solves,
        },
        "scientific_verdict": {
            "unique_krueger_mixture_recovered": False,
            "published_plot_can_support_exact_species_inversion": False,
            "code_can_propagate_any_declared_mixture": True,
            "next_non_target_fit_test": (
                "Run the Benck-conditioned mixture bands through the 0.25 s "
                "feature prefix, then the full 60 s profile ensemble. Those "
                "are cross-reactor sensitivities, not a Krueger boundary claim."
            ),
            "clean_closure": (
                "Obtain the original species-resolved HPEM PCMCM output or "
                "mass-resolved wafer ion measurements."
            ),
        },
        "sources": {
            "krueger_base_iead_sha256": base_iead.table_sha256,
            "krueger_figure16_iead_sha256": transfer.iead_table_sha256,
            "krueger_figure16_metadata_sha256": transfer.metadata_sha256,
            "krueger_thesis_text_sha256": _sha(KRUEGER_THESIS_TEXT),
            "benck_figure9_csv_sha256": _sha(BENCK),
            "benck_manifest_sha256": _sha(BENCK_MANIFEST),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    audit = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit["scientific_verdict"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
