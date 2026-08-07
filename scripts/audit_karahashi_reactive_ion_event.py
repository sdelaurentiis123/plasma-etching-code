#!/usr/bin/env python3
"""Audit what Karahashi's CF3+ beam data do and do not close.

The audit keeps three evidence layers separate:

* source-backed molecular inventory and mass-partitioned fragment energy;
* measured normal-incidence SiO2 removal yield from Figure 4;
* conditional SiFx product fractions from Figure 10, whose incidence angle
  is unreported.

Figure 4 and Figure 10 are joined only as an explicitly non-production
stoichiometric projection.  No depth, escape length, ion mixture, or surface
reservoir parameter is fitted.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from petch.experimental_data import (
    KARAHASHI_2007_FIGURE4_SHA256,
    KARAHASHI_2007_FIGURE10_SHA256,
    load_karahashi_2007_cf3_product_fractions,
    load_karahashi_2007_reactive_ion_yields,
)
from petch.reactive_ion_beam import Karahashi2007ReactiveIonYieldTable
from petch.reactive_ion_event import (
    Karahashi2007CF3ProductBranchTable,
    Karahashi2007ReactiveIonEventKernel,
    mass_partitioned_projectile_fragments,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "karahashi_2007"
OUTPUT = (
    ROOT / "results" / "curated" / "karahashi_reactive_ion_event"
    / "audit.json"
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def build_audit(root: Path = ROOT) -> dict[str, object]:
    data = root / "data" / "experimental" / "karahashi_2007"
    yield_path = data / "figure4_reactive_ion_yields.csv"
    product_path = data / "figure10_cf3_product_fractions.csv"
    yield_manifest = data / "digitization_manifest.json"
    product_manifest = data / "figure10_product_digitization_manifest.json"
    source_pdf = (
        root / "research_sources"
        / "karahashi_2007_hyomen_kagaku_28_60.pdf"
    )

    yield_table = Karahashi2007ReactiveIonYieldTable.from_observations(
        load_karahashi_2007_reactive_ion_yields(yield_path),
        source_table_sha256=KARAHASHI_2007_FIGURE4_SHA256,
    )
    product_table = Karahashi2007CF3ProductBranchTable(
        load_karahashi_2007_cf3_product_fractions(product_path),
        source_table_sha256=KARAHASHI_2007_FIGURE10_SHA256,
    )
    kernel = Karahashi2007ReactiveIonEventKernel(
        yield_table, product_table)

    energy_board = []
    for energy_eV in (500.0, 1000.0, 2000.0):
        fragments = mass_partitioned_projectile_fragments(
            "CF3+", energy_eV)
        fragment_rows = {
            fragment.element: {
                "multiplicity": fragment.multiplicity,
                "energy_per_fragment_eV": float(
                    fragment.energy_per_fragment_eV),
                "total_element_energy_eV": float(
                    fragment.multiplicity
                    * fragment.energy_per_fragment_eV),
            }
            for fragment in fragments
        }
        reconstructed_energy = sum(
            row["total_element_energy_eV"]
            for row in fragment_rows.values()
        )
        outcome = kernel.evaluate_conditional_cf3(
            energy_eV,
            acknowledge_unresolved_incidence_angle=True,
        )
        minimum_other_f_sink = max(
            outcome.unresolved_f_balance_lower, 0.0)
        minimum_surface_f_source = max(
            -outcome.unresolved_f_balance_upper, 0.0)
        energy_board.append({
            "energy_eV": energy_eV,
            "mass_partitioned_fragments": fragment_rows,
            "fragment_energy_reconstruction_eV": reconstructed_energy,
            "fragment_energy_residual_eV": (
                reconstructed_energy - energy_eV),
            "normal_incidence_figure4_yield_sio2_per_ion": (
                outcome.removed_sio2_formula_per_ion),
            "figure10_fraction_sum_before_normalization": (
                1.0 / outcome.branch_normalization_factor),
            "figure10_normalized_conditional_fractions": dict(
                outcome.normalized_product_fraction),
            "hypothetical_sifx_particles_per_ion": dict(
                outcome.conditional_sifx_particles_per_ion),
            "conditional_sifx_f_demand_per_ion": (
                outcome.required_f_atoms_per_ion),
            "conditional_sifx_f_demand_interval": [
                outcome.required_f_atoms_lower,
                outcome.required_f_atoms_upper,
            ],
            "incident_f_atoms_per_cf3_ion": 3.0,
            "unresolved_f_balance_incident_minus_conditional_sifx": (
                outcome.unresolved_f_balance_per_ion),
            "unresolved_f_balance_interval": [
                outcome.unresolved_f_balance_lower,
                outcome.unresolved_f_balance_upper,
            ],
            "minimum_f_sink_outside_sifx_within_digitization_interval": (
                minimum_other_f_sink),
            "minimum_resident_surface_f_source_within_digitization_interval": (
                minimum_surface_f_source),
            "condition_match_status": outcome.condition_match_status,
            "production_eligible": outcome.production_eligible,
        })

    maximum_energy_residual = max(
        abs(row["fragment_energy_residual_eV"])
        for row in energy_board)
    return {
        "audit_id": "KARAHASHI-REACTIVE-ION-EVENT-R1",
        "status": "completed_nonproduction_constraint",
        "question": (
            "Do Karahashi's mass-selected CF3+ data uniquely close a "
            "production reactive-ion event for feature-scale depth prediction?"
        ),
        "answer": False,
        "calibration_firewall": {
            "target_depth_used": False,
            "reactor_boundary_used": False,
            "fitted_parameters": [],
            "interpolated_product_branches": False,
            "figure4_figure10_join": (
                "explicitly_hypothetical_because_product_incidence_angle_"
                "is_unreported"
            ),
        },
        "source_receipts": {
            "figure4_data_path": str(yield_path.relative_to(root)),
            "figure4_data_sha256": _sha256(yield_path),
            "figure4_loader_checksum": KARAHASHI_2007_FIGURE4_SHA256,
            "figure4_manifest_path": str(yield_manifest.relative_to(root)),
            "figure4_manifest_sha256": _sha256(yield_manifest),
            "figure10_data_path": str(product_path.relative_to(root)),
            "figure10_data_sha256": _sha256(product_path),
            "figure10_loader_checksum": KARAHASHI_2007_FIGURE10_SHA256,
            "figure10_manifest_path": str(product_manifest.relative_to(root)),
            "figure10_manifest_sha256": _sha256(product_manifest),
            "source_pdf_path": str(source_pdf.relative_to(root)),
            "source_pdf_sha256": _sha256(source_pdf),
        },
        "evidence_partition": {
            "source_backed": [
                "mass-selected CF3+ ion identity",
                "CF3+ fragments into constituent atoms at impact",
                "kinetic energy is distributed among fragments by mass",
                "Figure-4 normal-incidence total SiO2 removal yields",
                "Figure-10 conditional fractions among detected SiFx products",
                "SiF is primarily prompt collision-cascade product",
                "SiF2 and SiF4 include delayed precursor-mediated desorption",
                "CO is observed as an oxygen-containing product",
            ],
            "derived_without_fit": [
                "atomic formula inventory C1F3",
                "mass-weighted fragment energies",
                "conditional SiFx fluorine-demand intervals",
            ],
            "not_measured_or_not_closed": [
                "Figure-10 ion incidence angle",
                "absolute product-resolved yields",
                "resident surface C/F inventory and exchange per impact",
                "oxygen coproduct branching and absolute CO yield",
                "prompt-versus-delayed numerical branching",
                "precursor diffusion, trapping, escape, and reincorporation",
                "film-state-dependent stopping and reaction probability",
                "reactor ion mixture and species-resolved IEAD",
            ],
        },
        "energy_board": energy_board,
        "numerical_gates": {
            "maximum_fragment_energy_residual_eV": maximum_energy_residual,
            "fragment_energy_tolerance_eV": 5.0e-13,
            "fragment_energy_conservation_pass": bool(
                maximum_energy_residual <= 5.0e-13),
            "all_conditional_rows_marked_nonproduction": all(
                not row["production_eligible"] for row in energy_board),
            "default_kernel_refuses_unresolved_condition_join": True,
        },
        "model_boundary": {
            "event_kernel_ready_for_feature_solver": False,
            "reason": (
                "A scalar total yield and a condition-unmatched product "
                "fraction do not determine state-resolved atom transfer. "
                "The surface must carry finite C/F/O/Si and delayed-precursor "
                "inventories before this evidence can be consumed without "
                "double counting the steady beam-conditioned film."
            ),
            "csda_fragment_transport_status": (
                "available_as_an_analytic_ZBL_Lindhard_straight_path_"
                "approximation_but_not_claimed_as_atomistic_MD"
            ),
            "production_integration_action": "refused_pending_independent_closure",
        },
        "next_discriminating_measurements": [
            "condition-matched ion incidence angle for the Figure-10 branches",
            "absolute species-resolved product yields versus energy and angle",
            "dose-resolved surface C/F/O/Si composition or film thickness",
            "time-resolved product branching sufficient to identify precursor "
            "creation and loss rates",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = json.dumps(
        build_audit(), indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"{arguments.output} is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
