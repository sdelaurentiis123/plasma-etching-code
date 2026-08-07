#!/usr/bin/env python3
"""Test the Guo/Kwon surface closure on independent reactive-ion beams.

Karahashi's mass-selected, radical-free beams remove the reactor and feature
transport ambiguities: each row supplies one ion identity, one normal-incidence
energy, and one measured SiO2 yield.  This script evaluates the frozen Guo
reaction deck at all 21 digitized points.  No coefficient or incident
composition is fitted to the Karahashi observations.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.experimental_data import (
    KARAHASHI_2007_FIGURE4_SHA256,
    load_karahashi_2007_reactive_ion_yields,
)
from petch.guo_c4f8_sio2 import (
    GuoC4F8ArSiO2Mechanism,
    GuoIncidentComposition,
    GuoIonQuadrature,
)


ROOT = Path(__file__).resolve().parents[1]
KARAHASHI_DIRECTORY = (
    ROOT / "data" / "experimental" / "karahashi_2007")
GUO_DIRECTORY = ROOT / "data" / "surface_interactions" / "guo_2009"
OUTPUT = (
    ROOT / "results" / "curated" / "guo_karahashi_transfer" / "audit.json")


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _metrics(rows):
    observed = np.asarray([row["observed_yield_sio2_per_ion"] for row in rows])
    predicted = np.asarray([row["predicted_yield_sio2_per_ion"] for row in rows])
    error = predicted - observed
    relative = error / observed
    return {
        "point_count": len(rows),
        "mean_absolute_error_sio2_per_ion": float(np.mean(np.abs(error))),
        "root_mean_square_error_sio2_per_ion": float(
            np.sqrt(np.mean(error ** 2))),
        "mean_absolute_relative_error": float(np.mean(np.abs(relative))),
        "maximum_absolute_relative_error": float(np.max(np.abs(relative))),
        "within_plotted_source_interval_count": sum(
            row["prediction_within_plotted_source_interval"] for row in rows),
        "correct_net_etch_sign_count": sum(
            row["predicted_yield_sio2_per_ion"] > 0.0 for row in rows),
    }


def _state_payload(result):
    return {
        "Si": result.state.si,
        "O": result.state.o,
        "C": result.state.c,
        "F": result.state.f,
        "V": result.state.vacancy,
    }


def build_audit(root: Path | None = None):
    root = root or ROOT
    karahashi_directory = (
        root / "data" / "experimental" / "karahashi_2007")
    guo_directory = root / "data" / "surface_interactions" / "guo_2009"
    data_path = karahashi_directory / "figure4_reactive_ion_yields.csv"
    digitization_manifest_path = (
        karahashi_directory / "digitization_manifest.json")
    guo_manifest_path = guo_directory / "source_manifest.json"
    guo_deck_path = guo_directory / "table4_1_reaction_deck.csv"

    observations = load_karahashi_2007_reactive_ion_yields(data_path)
    point_rows = []
    for observation in observations:
        mechanism = GuoC4F8ArSiO2Mechanism(
            GuoIncidentComposition(
                neutral_flux_ratio={},
                ion_fraction={observation.species: 1.0},
            ),
            GuoIonQuadrature.monoenergetic(observation.energy_eV),
        )
        result = mechanism.solve_steady_state_complementarity_bdf(
            maximum_coordinate=30.0)
        prediction = result.sio2_yield_per_ion
        error = prediction - observation.yield_sio2_per_ion
        atom_ledger_residual = (
            result.movement_atoms_per_ion
            - sum(result.removed_atoms_per_ion.values())
            + sum(result.incoming_atoms_per_ion.values())
        )
        point_rows.append({
            "species": observation.species,
            "energy_eV": observation.energy_eV,
            "observed_yield_sio2_per_ion": (
                observation.yield_sio2_per_ion),
            "plotted_lower_yield": observation.plotted_lower_yield,
            "plotted_upper_yield": observation.plotted_upper_yield,
            "digitization_yield_uncertainty": (
                observation.digitization_yield_uncertainty),
            "predicted_yield_sio2_per_ion": prediction,
            "signed_error_sio2_per_ion": error,
            "signed_relative_error": (
                error / observation.yield_sio2_per_ion),
            "prediction_within_plotted_source_interval": bool(
                observation.plotted_lower_yield
                <= prediction
                <= observation.plotted_upper_yield
            ),
            "correct_net_etch_sign": bool(prediction > 0.0),
            "state": _state_payload(result),
            "movement_atoms_per_ion": result.movement_atoms_per_ion,
            "incoming_atoms_per_ion": dict(result.incoming_atoms_per_ion),
            "removed_atoms_per_ion": dict(result.removed_atoms_per_ion),
            "atom_ledger_residual": atom_ledger_residual,
            "steady_state_residual": result.steady_state_residual,
            "steady_state_solver": (
                result.source_extrapolation["steady_state_solver"]),
            "beyond_guo_source_fit_energy": (
                result.source_extrapolation["beyond_source_fit_energy"]),
        })

    species_metrics = {
        species: _metrics([
            row for row in point_rows if row["species"] == species])
        for species in ("F+", "CF+", "CF2+", "CF3+")
    }
    source_supported = [
        row for row in point_rows
        if not row["beyond_guo_source_fit_energy"]
    ]
    extrapolated = [
        row for row in point_rows
        if row["beyond_guo_source_fit_energy"]
    ]
    maximum_atom_residual = max(
        abs(row["atom_ledger_residual"]) for row in point_rows)
    maximum_steady_residual = max(
        row["steady_state_residual"] for row in point_rows)

    return {
        "audit_id": "GUO-KARAHASHI-REACTIVE-ION-TRANSFER-R1",
        "status": "completed_mixed_transfer",
        "question": (
            "Does the frozen atom-balanced Guo/Kwon C4F8-Ar/SiO2 closure "
            "predict independent mass-selected reactive-ion beam yields?"
        ),
        "calibration_firewall": {
            "karahashi_yields_used_by_surface_solver": False,
            "fitted_parameters": [],
            "incident_neutral_flux": "exactly_zero_as_in_beam_experiment",
            "incident_ion_fraction": (
                "exactly_one_for_the_measured_mass_selected_species"),
            "ion_energy": "measured_experimental_setpoint",
            "incidence": "normal",
            "comparison_operation": (
                "evaluate_frozen_Guo_deck_then_score_against_observation"),
        },
        "source_receipts": {
            "karahashi_data_path": str(data_path.relative_to(root)),
            "karahashi_data_sha256": _sha256(data_path),
            "loader_pinned_karahashi_data_sha256": (
                KARAHASHI_2007_FIGURE4_SHA256),
            "karahashi_digitization_manifest_path": str(
                digitization_manifest_path.relative_to(root)),
            "karahashi_digitization_manifest_sha256": _sha256(
                digitization_manifest_path),
            "guo_source_manifest_path": str(
                guo_manifest_path.relative_to(root)),
            "guo_source_manifest_sha256": _sha256(guo_manifest_path),
            "guo_reaction_deck_path": str(guo_deck_path.relative_to(root)),
            "guo_reaction_deck_sha256": _sha256(guo_deck_path),
        },
        "experiment": {
            "target": "SiO2",
            "species": ["F+", "CF+", "CF2+", "CF3+"],
            "mass_selected": True,
            "neutral_radical_flux": "none",
            "gas_phase_reactions": "none",
            "incidence": "normal",
            "energy_support_eV": [250.0, 2000.0],
            "measurement_uncertainty_semantics": (
                "source_plots_error_bars_but_accessible_text_does_not_define_"
                "their_statistical_meaning"),
        },
        "model_evidence_boundary": {
            "guo_source_fit_maximum_energy_eV": (
                GuoC4F8ArSiO2Mechanism.source_energy_max_eV),
            "points_inside_guo_source_energy_support": len(source_supported),
            "points_beyond_guo_source_energy_support": len(extrapolated),
            "surface_evidence_ceiling": "L1_yield_regressed",
            "atomistic_trajectory_or_interatomic_potential": False,
            "first_principles_claimed": False,
        },
        "all_point_metrics": _metrics(point_rows),
        "per_species_metrics": species_metrics,
        "point_board": point_rows,
        "numerical_and_atomic_gates": {
            "maximum_steady_state_residual": maximum_steady_residual,
            "steady_state_residual_tolerance": 2.0e-8,
            "all_steady_states_pass": bool(
                maximum_steady_residual <= 2.0e-8),
            "maximum_absolute_atom_ledger_residual": maximum_atom_residual,
            "atom_ledger_tolerance": 1.0e-12,
            "all_atom_ledgers_pass": bool(
                maximum_atom_residual <= 1.0e-12),
        },
        "verdict": {
            "species_resolved_transfer_validated": False,
            "reason": (
                "The frozen closure does not reproduce the full independent "
                "species ladder. It overpredicts F+ and CF+ removal, and "
                "predicts net deposition for the measured 250 and 500 eV "
                "CF3+ etch points. Some CF2+/CF3+ high-energy points are "
                "closer, but selecting those after inspection is not a "
                "validation gate."
            ),
            "krueger_ion_identity_envelope_implication": (
                "The all-CF2+ and all-CF3+ Krueger depth-scale bracket is a "
                "boundary sensitivity, not an independently validated "
                "reactive-ion closure. Its apparent agreement cannot identify "
                "the reactor ion mixture."
            ),
            "physics_defects_exposed": [
                (
                    "generic ion incorporation lacks a species- and energy-"
                    "resolved threshold/branch closure"),
                (
                    "the square-root reaction moments do not transfer "
                    "quantitatively across F+, CF+, CF2+, and CF3+"),
                (
                    "most beam points are above the <=370 eV Guo/Yin fit "
                    "support, so a high-energy stopping/fragmentation law is "
                    "still missing"),
            ],
            "next_required_evidence": [
                (
                    "fit no Karahashi point in isolation; derive or import a "
                    "species-resolved reactive-ion stopping, implantation, "
                    "fragmentation, and chemical-removal closure"),
                (
                    "calibrate only on a declared subset and preregister "
                    "held-out species/energy blocks before claiming transfer"),
                (
                    "retain the measured tabulation as the highest-evidence "
                    "closure inside its exact species/energy/angle support"),
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
    arguments = parser.parse_args()
    payload = canonical_payload(build_audit())
    if arguments.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(payload, encoding="utf-8")
    elif arguments.check or OUTPUT.exists():
        if not OUTPUT.exists() or OUTPUT.read_text(
                encoding="utf-8") != payload:
            raise RuntimeError("committed Guo/Karahashi audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
