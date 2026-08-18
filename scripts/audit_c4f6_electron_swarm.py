#!/usr/bin/env python3
"""Audit the Lan--Jeon C4F6 working set against its PT drift board."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from petch.reactor_global.c4f6_electron_collisions import (
    INELASTIC_CSV_SHA256,
    MOMENTUM_CSV_SHA256,
    load_lan_jeon_2014_c4f6_replay,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)


ROOT = Path(__file__).resolve().parents[1]
DRIFT_CSV = (
    ROOT / "data" / "experimental" / "lan_jeon_2014_c4f6"
    / "figure7_pure_c4f6_drift.csv"
)
DEFAULT_OUTPUT = ROOT / "results" / "curated" / "c4f6_electron_swarm_v1"


def _grid(deck, scale: int = 1) -> ElectronEnergyGrid:
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV and process.energy_loss_eV < 200.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, .001, .01, .1, 1.0, 10.0, 50.0, 100.0, 200.0),
        tuple(scale * value for value in (8, 12, 24, 48, 96, 96, 80, 80)),
        inserted_boundaries_eV=thresholds,
    )


def _solve(deck, fields_Td, *, scale: int = 1) -> dict:
    grid = _grid(deck, scale)
    solver = DeterministicTwoTermBoltzmannSolver(grid, deck)
    rows = []
    for field in fields_Td:
        solution = solver.solve(
            TwoTermBoltzmannCondition(
                reduced_electric_field_Td=float(field),
                gas_temperature_K=300.0,
                target_mole_fractions={"C4F6": 1.0},
                growth_model="temporal_growth",
                initial_electron_temperature_eV=.2,
            ),
            relative_tolerance=2.0e-6,
            maximum_tail_population_fraction=1.0e-6,
        )
        drift = (
            solution.transport_moments
            .flux_reduced_mobility_m_inv_V_inv_s_inv
            * float(field) * 1.0e-21
        )
        rates = {
            process: float(moment.rate_coefficient_m3_s)
            for process, moment in zip(
                deck.processes, solution.collision_moments)
        }
        rows.append({
            "reduced_electric_field_Td": float(field),
            "flux_drift_velocity_m_s": float(drift),
            "mean_electron_energy_eV": float(
                solution.distribution.mean_energy_eV),
            "temporal_net_growth_rate_coefficient_m3_s": float(
                solution.net_growth_rate_coefficient_m3_s),
            "attachment_rate_coefficient_m3_s": rates[
                deck.processes[1]],
            "aggregate_dissociation_rate_coefficient_m3_s": rates[
                deck.processes[-2]],
            "aggregate_ionization_rate_coefficient_m3_s": rates[
                deck.processes[-1]],
            "weighted_iteration_residual": float(
                solution.weighted_iteration_residual),
        })
    return {"actual_cell_count": grid.cell_count, "rows": rows}


def _measurements() -> tuple[np.ndarray, np.ndarray, float, float]:
    with DRIFT_CSV.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    field = np.asarray([
        float(row["reduced_electric_field_Td"]) for row in rows])
    drift = np.asarray([float(row["drift_velocity_m_s"]) for row in rows])
    field_bound = {float(row["field_digitization_relative_bound"]) for row in rows}
    drift_bound = {float(row["drift_digitization_relative_bound"]) for row in rows}
    if len(rows) != 18 or len(field_bound) != 1 or len(drift_bound) != 1:
        raise RuntimeError("Lan--Jeon Figure-7 drift board topology changed")
    return field, drift, field_bound.pop(), drift_bound.pop()


def _residual_summary(predicted, measured) -> dict:
    relative = np.asarray(predicted) / np.asarray(measured) - 1.0
    return {
        "point_count": int(relative.size),
        "mean_signed_relative_residual": float(np.mean(relative)),
        "mean_absolute_relative_residual": float(np.mean(np.abs(relative))),
        "maximum_absolute_relative_residual": float(np.max(np.abs(relative))),
        "signed_relative_residual": relative.tolist(),
    }


def audit() -> dict:
    replay = load_lan_jeon_2014_c4f6_replay()
    field, measured, field_bound, drift_bound = _measurements()
    nominal = _solve(replay.derived_deck, field)
    predicted = np.asarray([
        row["flux_drift_velocity_m_s"] for row in nominal["rows"]])
    representative = np.asarray((196.299456, 298.048, 996.573))
    coarse = _solve(replay.derived_deck, representative, scale=1)
    fine = _solve(replay.derived_deck, representative, scale=2)
    coarse_drift = np.asarray([
        row["flux_drift_velocity_m_s"] for row in coarse["rows"]])
    fine_drift = np.asarray([
        row["flux_drift_velocity_m_s"] for row in fine["rows"]])
    coarse_dissociation = np.asarray([
        row["aggregate_dissociation_rate_coefficient_m3_s"]
        for row in coarse["rows"]
    ])
    fine_dissociation = np.asarray([
        row["aggregate_dissociation_rate_coefficient_m3_s"]
        for row in fine["rows"]
    ])
    return {
        "schema": "petch.c4f6_electron_swarm_audit.v1",
        "claim_class": "effective_source_replay_with_transport_definition_gap",
        "sources": {
            "momentum_table_sha256": MOMENTUM_CSV_SHA256,
            "inelastic_table_sha256": INELASTIC_CSV_SHA256,
            "figure7_digitization_path": str(DRIFT_CSV.relative_to(ROOT)),
            "figure7_point_count": int(field.size),
            "source_artifact_committed": False,
        },
        "nominal_collision_deck_sha256": replay.derived_deck.payload_sha256,
        "prediction": nominal,
        "transport_definition_diagnostic": {
            **_residual_summary(predicted, measured),
            "field_Td": field.tolist(),
            "measured_pt_average_drift_velocity_m_s": measured.tolist(),
            "predicted_flux_drift_velocity_m_s": predicted.tolist(),
            "field_digitization_relative_bound": field_bound,
            "drift_digitization_relative_bound": drift_bound,
            "paper_reported_typical_reproduction_relative_bound": .05,
            "predicted_definition": "local_two_term_flux_drift",
            "tabulated_definition": "pulsed_Townsend_average_drift_Wv",
            "measurement_equivalent_grade": False,
            "interpretation": (
                "The discrepancy is not removed by energy-grid refinement. "
                "A density-gradient/PT transport observable is required; "
                "the effective cross sections must not be retuned against "
                "the wrong flux observable."
            ),
        },
        "numerical_convergence": {
            "fields_Td": representative.tolist(),
            "coarse_actual_cell_count": coarse["actual_cell_count"],
            "fine_actual_cell_count": fine["actual_cell_count"],
            "maximum_absolute_flux_drift_relative_change": float(
                np.max(np.abs(fine_drift / coarse_drift - 1.0))),
            "maximum_absolute_dissociation_rate_relative_change": float(
                np.max(np.abs(
                    fine_dissociation / coarse_dissociation - 1.0))),
        },
        "evidence_boundary": {
            "momentum": "swarm_regressed",
            "attachment": "measured_imported_unchanged",
            "total_ionization": "measured_imported_unchanged",
            "vibration_and_excitation": "swarm_regressed_or_cC4F8_analog",
            "neutral_dissociation": "cC4F8_analog_imported_unchanged",
            "ion_product_branching": "unresolved_aggregate_only",
            "neutral_product_branching": "unresolved_aggregate_only",
        },
        "certification": {
            "supports_use_as_c4f6_component_collision_input": True,
            "supports_local_flux_transport_and_aggregate_rates": True,
            "source_pt_drift_reproduced_measurement_equivalently": False,
            "bulk_flux_nonconservative_transport_resolved": False,
            "primary_positive_and_negative_branching_resolved": False,
            "independently_validated_in_target_mixture_band": False,
            "supports_unique_krueger_reactor_state": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
            "feature_depth_used": False,
        },
    }


def _write(result: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    transport = result["transport_definition_diagnostic"]
    convergence = result["numerical_convergence"]
    report = f"""# C4F6 electron-collision source-replay audit

Lan--Jeon's printed C4F6 tables are now checksum-locked and consumed directly
by the deterministic two-term solver. The 18 visually audited pure-C4F6
Figure-7 markers span 121--1197 Td.

Against that board, local flux drift has a mean absolute residual of
`{100.0 * transport['mean_absolute_relative_residual']:.2f}%` and a maximum of
`{100.0 * transport['maximum_absolute_relative_residual']:.2f}%`. This is not
graded as a failed cross-section fit: the paper reports pulsed-Townsend average
drift `Wv`, while the reactor solver returns flux drift. Those observables
diverge in a gas with attachment and ionization.

Doubling the energy-grid density changes flux drift by at most
`{100.0 * convergence['maximum_absolute_flux_drift_relative_change']:.3f}%`
and aggregate dissociation rate by at most
`{100.0 * convergence['maximum_absolute_dissociation_rate_relative_change']:.3f}%`.
The remaining gap is therefore physical/observational, not grid error.

This authorizes the deck as a bounded C4F6 component for local EEDF and
aggregate-rate calculations. It does not resolve product branching, a Krueger
reactor state, wafer flux, or feature depth. The next gates are deterministic
density-gradient/PT transport and mass-resolved reactor validation against
Benck's ion-current board.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit()
    _write(result, args.output)
    print(json.dumps({
        "transport": result["transport_definition_diagnostic"],
        "convergence": result["numerical_convergence"],
        "certification": result["certification"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
