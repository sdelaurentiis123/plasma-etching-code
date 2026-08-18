#!/usr/bin/env python3
"""Audit the NIST-evaluated SF6 aggregate electron deck.

The audit deliberately separates source replay from measurement-equivalent
validation.  The local two-term solver returns flux transport and temporal
growth, while the pure-SF6 review tabulates transit/Townsend observables in a
strongly nonconservative gas.  Those are useful diagnostics but not silently
treated as the same observable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)
from petch.reactor_global.sf6_electron_collisions import (
    NIST_SF6_TABLE14_CSV_SHA256,
    NIST_SF6_TABLE15_CSV_SHA256,
    NIST_SF6_TABLE17_CSV_SHA256,
    NIST_SF6_TABLE20_CSV_SHA256,
    NIST_SF6_TABLE28_CSV_SHA256,
    NIST_SF6_TABLE35_CSV_SHA256,
    NIST_SF6_TABLE36_CSV_SHA256,
    NIST_SF6_TABLE37_CSV_SHA256,
    NIST_SF6_TABLE9_CSV_SHA256,
    derive_nist_evaluated_sf6_replay,
    load_nist_2000_sf6_attachment_rate_curve,
    load_nist_2000_sf6_drift_curve,
    load_nist_2000_sf6_effective_ionization_curve,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "curated" / "sf6_electron_swarm_v1"
DIRECT_DRIFT_INTERVAL_TD = (275.0, 1000.0)
SOURCE_CRITICAL_FIELD_TD = 359.3
SOURCE_CRITICAL_FIELD_UNCERTAINTY_TD = 3.0


def _grid(deck, *, low_energy_scale: int = 5) -> ElectronEnergyGrid:
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV and process.energy_loss_eV < 200.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, .0001, .001, .01, .1, 1.0, 10.0, 40.0, 120.0, 200.0),
        (
            8 * low_energy_scale,
            12 * low_energy_scale,
            18 * low_energy_scale,
            36 * low_energy_scale,
            48, 96, 96, 120, 160,
        ),
        inserted_boundaries_eV=thresholds,
    )


def _solve(deck, fields_Td, *, low_energy_scale: int = 5) -> dict:
    grid = _grid(deck, low_energy_scale=low_energy_scale)
    solver = DeterministicTwoTermBoltzmannSolver(grid, deck)
    rows = []
    for field in fields_Td:
        solution = solver.solve(
            TwoTermBoltzmannCondition(
                reduced_electric_field_Td=float(field),
                gas_temperature_K=296.0,
                target_mole_fractions={"SF6": 1.0},
                growth_model="temporal_growth",
                initial_electron_temperature_eV=float(
                    np.clip(.015 * float(field), 2.0, 12.0)
                ),
            ),
            relative_tolerance=3.0e-6,
            maximum_tail_population_fraction=2.0e-6,
        )
        drift = (
            solution.transport_moments
            .flux_reduced_mobility_m_inv_V_inv_s_inv
            * float(field) * 1.0e-21
        )
        attachment = next(
            moment.rate_coefficient_m3_s
            for moment in solution.collision_moments
            if moment.process_kind == "ATTACHMENT"
        )
        ionization = next(
            moment.rate_coefficient_m3_s
            for moment in solution.collision_moments
            if moment.process_kind == "IONIZATION"
        )
        rows.append({
            "reduced_electric_field_Td": float(field),
            "flux_drift_velocity_m_s": float(drift),
            "temporal_net_growth_rate_coefficient_m3_s": float(
                solution.net_growth_rate_coefficient_m3_s
            ),
            "flux_growth_over_drift_m2": float(
                solution.net_growth_rate_coefficient_m3_s / drift
            ),
            "attachment_rate_coefficient_m3_s": float(attachment),
            "ionization_rate_coefficient_m3_s": float(ionization),
            "mean_electron_energy_eV": float(
                solution.distribution.mean_energy_eV
            ),
            "weighted_iteration_residual": float(
                solution.weighted_iteration_residual
            ),
        })
    return {
        "actual_cell_count": grid.cell_count,
        "rows": rows,
    }


def _map(rows: list[dict], key: str) -> dict[float, float]:
    return {
        row["reduced_electric_field_Td"]: row[key] for row in rows
    }


def _residual_summary(predicted, measured) -> dict:
    relative = np.asarray(predicted) / np.asarray(measured) - 1.0
    return {
        "point_count": int(relative.size),
        "median_absolute_relative_residual": float(
            np.median(np.abs(relative))
        ),
        "maximum_absolute_relative_residual": float(
            np.max(np.abs(relative))
        ),
        "signed_relative_residual": relative.tolist(),
    }


def _zero_crossing(rows: list[dict]) -> float:
    field = np.asarray([
        row["reduced_electric_field_Td"] for row in rows
    ])
    growth = np.asarray([
        row["temporal_net_growth_rate_coefficient_m3_s"] for row in rows
    ])
    for index in range(field.size - 1):
        if growth[index] <= 0.0 < growth[index + 1]:
            fraction = -growth[index] / (growth[index + 1] - growth[index])
            return float(field[index] + fraction * (field[index + 1] - field[index]))
    raise RuntimeError("SF6 temporal growth did not cross zero")


def audit() -> dict:
    drift = load_nist_2000_sf6_drift_curve()
    effective = load_nist_2000_sf6_effective_ionization_curve()
    attachment = load_nist_2000_sf6_attachment_rate_curve()
    nominal = derive_nist_evaluated_sf6_replay()
    fields = np.unique(np.concatenate((
        drift.reduced_electric_field_Td[drift.recommended_mask],
        effective.reduced_electric_field_Td[
            effective.reduced_electric_field_Td <= 1000.0
        ],
        attachment.reduced_electric_field_Td,
    )))
    prediction = _solve(nominal.derived_deck, fields)
    rows = prediction["rows"]

    drift_map = _map(rows, "flux_drift_velocity_m_s")
    drift_field = drift.reduced_electric_field_Td[drift.recommended_mask]
    drift_measured = drift.drift_velocity_m_s[drift.recommended_mask]
    drift_predicted = [drift_map[float(item)] for item in drift_field]

    attach_map = _map(rows, "attachment_rate_coefficient_m3_s")
    attach_predicted = [
        attach_map[float(item)]
        for item in attachment.reduced_electric_field_Td
    ]

    effective_map = _map(rows, "flux_growth_over_drift_m2")
    effective_mask = effective.reduced_electric_field_Td <= 1000.0
    effective_field = effective.reduced_electric_field_Td[effective_mask]
    effective_measured = effective.effective_ionization_coefficient_m2[
        effective_mask
    ]
    effective_predicted = [
        effective_map[float(item)] for item in effective_field
    ]

    representative = np.asarray((350.0, 650.0, 1000.0))
    coarse = _solve(nominal.derived_deck, representative)
    fine = _solve(
        nominal.derived_deck,
        representative,
        low_energy_scale=6,
    )
    coarse_drift = np.asarray([
        row["flux_drift_velocity_m_s"] for row in coarse["rows"]
    ])
    fine_drift = np.asarray([
        row["flux_drift_velocity_m_s"] for row in fine["rows"]
    ])
    coarse_attach = np.asarray([
        row["attachment_rate_coefficient_m3_s"] for row in coarse["rows"]
    ])
    fine_attach = np.asarray([
        row["attachment_rate_coefficient_m3_s"] for row in fine["rows"]
    ])

    alternate = derive_nist_evaluated_sf6_replay(
        vibrational_energy_loss_eV=.117
    )
    alternate_prediction = _solve(alternate.derived_deck, representative)
    alternate_drift = np.asarray([
        row["flux_drift_velocity_m_s"]
        for row in alternate_prediction["rows"]
    ])

    return {
        "schema": "petch.sf6_electron_swarm_audit.v1",
        "claim_class": (
            "evaluated_source_replay_with_nonconservative_transport_"
            "definition_diagnostics"
        ),
        "sources": {
            "table9_sha256": NIST_SF6_TABLE9_CSV_SHA256,
            "table14_sha256": NIST_SF6_TABLE14_CSV_SHA256,
            "table15_sha256": NIST_SF6_TABLE15_CSV_SHA256,
            "table17_sha256": NIST_SF6_TABLE17_CSV_SHA256,
            "table20_sha256": NIST_SF6_TABLE20_CSV_SHA256,
            "table28_sha256": NIST_SF6_TABLE28_CSV_SHA256,
            "table35_sha256": NIST_SF6_TABLE35_CSV_SHA256,
            "table36_sha256": NIST_SF6_TABLE36_CSV_SHA256,
            "table37_sha256": NIST_SF6_TABLE37_CSV_SHA256,
            "source_artifacts_committed": False,
        },
        "nominal_collision_deck_sha256": (
            nominal.derived_deck.payload_sha256
        ),
        "prediction": prediction,
        "source_replay": {
            "attachment_rate": {
                **_residual_summary(
                    attach_predicted,
                    attachment.attachment_rate_coefficient_m3_s,
                ),
                "evidence_dependency": (
                    "Table 37 is assessed from eta/N times drift and is not "
                    "independent of the reviewed attachment evidence"
                ),
            },
        },
        "transport_definition_diagnostics": {
            "drift": {
                **_residual_summary(drift_predicted, drift_measured),
                "field_interval_Td": list(DIRECT_DRIFT_INTERVAL_TD),
                "predicted_definition": "local_two_term_flux_drift",
                "tabulated_definition": (
                    "reviewed swarm drift; flux-versus-bulk not identified"
                ),
                "measurement_equivalent_grade": False,
            },
            "effective_ionization": {
                "field_Td": effective_field.tolist(),
                "predicted_flux_growth_over_drift_m2": effective_predicted,
                "tabulated_effective_townsend_m2": effective_measured.tolist(),
                "predicted_definition": "temporal_growth_rate_over_flux_drift",
                "tabulated_definition": "steady_state_spatial_Townsend_alpha_minus_eta",
                "measurement_equivalent_grade": False,
            },
            "critical_field": {
                "predicted_temporal_zero_Td": _zero_crossing(rows),
                "source_spatial_townsend_zero_Td": SOURCE_CRITICAL_FIELD_TD,
                "source_uncertainty_Td": SOURCE_CRITICAL_FIELD_UNCERTAINTY_TD,
                "measurement_equivalent_grade": False,
            },
        },
        "closure_sensitivity": {
            "fields_Td": representative.tolist(),
            "nominal_vibrational_loss_eV": .095,
            "alternate_vibrational_loss_eV": .117,
            "flux_drift_relative_change": (
                alternate_drift / coarse_drift - 1.0
            ).tolist(),
        },
        "numerical_convergence": {
            "fields_Td": representative.tolist(),
            "coarse_actual_cell_count": coarse["actual_cell_count"],
            "fine_actual_cell_count": fine["actual_cell_count"],
            "maximum_absolute_flux_drift_relative_change": float(
                np.max(np.abs(fine_drift / coarse_drift - 1.0))
            ),
            "maximum_absolute_attachment_rate_relative_change": float(
                np.max(np.abs(fine_attach / coarse_attach - 1.0))
            ),
        },
        "certification": {
            "supports_use_as_sf6_component_collision_input": True,
            "primary_positive_and_negative_branching_resolved": False,
            "bulk_flux_nonconservative_transport_resolved": False,
            "independently_validated_in_target_mixture_band": False,
            "supports_unique_target_reactor_state": False,
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
    attach = result["source_replay"]["attachment_rate"]
    drift = result["transport_definition_diagnostics"]["drift"]
    critical = result["transport_definition_diagnostics"]["critical_field"]
    report = f"""# SF6 evaluated electron-collision audit

The exact NIST recommended/suggested/deduced tables are now replayable in SI
units. The aggregate deck avoids total-scattering double counting and carries
explicit vibration, neutral-dissociation, ionization, and attachment closures.

The source-derived attachment-rate replay has a median absolute residual of
`{100.0 * attach['median_absolute_relative_residual']:.2f}%` (maximum
`{100.0 * attach['maximum_absolute_relative_residual']:.2f}%`). This is a
source-consistency check, not independent validation.

The raw flux-drift comparison has a median residual of
`{100.0 * drift['median_absolute_relative_residual']:.2f}%`. It is deliberately
not graded as measurement-equivalent: in strongly attaching SF6 the solver's
local flux drift/temporal growth and the source's swarm drift/spatial Townsend
observables diverge. The corresponding critical fields are
`{critical['predicted_temporal_zero_Td']:.1f} Td` temporal versus
`{critical['source_spatial_townsend_zero_Td']:.1f} +/- {critical['source_uncertainty_Td']:.1f} Td`
spatial. Closing that distinction requires a spatial-growth/bulk-transport
solver, not parameter fitting.

Representative-grid maximum changes are
`{100.0 * result['numerical_convergence']['maximum_absolute_flux_drift_relative_change']:.3f}%`
for flux drift and
`{100.0 * result['numerical_convergence']['maximum_absolute_attachment_rate_relative_change']:.3f}%`
for attachment rate.

This authorizes SF6 as a bounded component of the mixed-gas EEDF. It does not
by itself authorize a unique Oxford reactor state, wafer flux, or feature depth.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit()
    _write(result, args.output)
    print(json.dumps({
        "attachment": result["source_replay"]["attachment_rate"],
        "drift": result["transport_definition_diagnostics"]["drift"],
        "critical_field": result[
            "transport_definition_diagnostics"
        ]["critical_field"],
        "convergence": result["numerical_convergence"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
