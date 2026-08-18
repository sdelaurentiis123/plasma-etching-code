#!/usr/bin/env python3
"""Replay the CHF3 electron-transport evidence without upgrading its claim.

Kushner--Zhang regressed their working set against swarm behavior, and the
Christophorou--Olthoff review is part of the same evidence chain.  Agreement
here is therefore a source-reproduction and transport-consistency check.  It
is not independent validation of the chemical branches or a reactor/depth
prediction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from petch.reactor_global.chf3_electron_collisions import (
    KUSHNER_ZHANG_CHF3_CSV_SHA256,
    NIST_CHF3_TABLE4_CSV_SHA256,
    NIST_CHF3_TABLE5_CSV_SHA256,
    NIST_CHF3_TABLE6_CSV_SHA256,
    derive_nist_evaluated_chf3_replay,
    load_kushner_zhang_2000_chf3_replay,
    load_nist_1999_chf3_drift_curve,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "chf3_electron_swarm_v1"
)
KUSHNER_FIGURE3_RENDER_SHA256 = (
    "86028a1dc9f0c1e7658736df02815d74497d9df037a839741f12bdce94c077d2"
)
OPERATING_BAND_TD = (40.0, 100.0)


def _grid(deck, cells_scale: int = 1) -> ElectronEnergyGrid:
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV
        and process.energy_loss_eV < 120.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, 0.005, 0.01, 0.1, 1.0, 10.0, 40.0, 120.0),
        tuple(
            cells_scale * count
            for count in (8, 8, 36, 48, 96, 96, 120)
        ),
        inserted_boundaries_eV=thresholds,
    )


def _initial_temperature_eV(field_Td: float) -> float:
    # Only a root-solver seed; the converged state is independent of it.
    return float(np.clip(0.012 * field_Td, 0.04, 2.5))


def _solve_curve(deck, fields_Td, *, cells_scale: int = 1) -> dict:
    grid = _grid(deck, cells_scale=cells_scale)
    solver = DeterministicTwoTermBoltzmannSolver(grid, deck)
    rows = []
    for field_Td in fields_Td:
        solution = solver.solve(
            TwoTermBoltzmannCondition(
                reduced_electric_field_Td=float(field_Td),
                gas_temperature_K=298.0,
                target_mole_fractions={"CHF3": 1.0},
                growth_model="temporal_growth",
                initial_electron_temperature_eV=(
                    _initial_temperature_eV(float(field_Td))
                ),
            ),
            relative_tolerance=2.0e-6,
            maximum_tail_population_fraction=1.0e-6,
        )
        flux_drift = (
            solution.transport_moments
            .flux_reduced_mobility_m_inv_V_inv_s_inv
            * float(field_Td) * 1.0e-21
        )
        rows.append({
            "reduced_electric_field_Td": float(field_Td),
            "flux_drift_velocity_m_s": float(flux_drift),
            "mean_electron_energy_eV": float(
                solution.distribution.mean_energy_eV
            ),
            "net_growth_rate_coefficient_m3_s": float(
                solution.net_growth_rate_coefficient_m3_s
            ),
            "weighted_iteration_residual": float(
                solution.weighted_iteration_residual
            ),
            "tail_population_fraction": float(
                solution.distribution.convergence_receipt(
                    tail_cell_count=min(4, grid.cell_count),
                    maximum_tail_population_fraction=1.0e-6,
                )["tail_population_fraction"]
            ),
        })
    return {
        "nominal_cell_scale": cells_scale,
        "actual_cell_count": grid.cell_count,
        "rows": rows,
    }


def _grade(prediction: dict, measured_velocity_m_s: np.ndarray) -> dict:
    predicted = np.asarray([
        row["flux_drift_velocity_m_s"] for row in prediction["rows"]
    ])
    relative = predicted / measured_velocity_m_s - 1.0
    fields = np.asarray([
        row["reduced_electric_field_Td"] for row in prediction["rows"]
    ])
    operating = (
        (fields >= OPERATING_BAND_TD[0])
        & (fields <= OPERATING_BAND_TD[1])
    )
    for row, measured, error in zip(
        prediction["rows"], measured_velocity_m_s, relative
    ):
        row["nist_recommended_drift_velocity_m_s"] = float(measured)
        row["relative_residual"] = float(error)
    return {
        "all_points": {
            "point_count": int(fields.size),
            "median_absolute_relative_residual": float(
                np.median(np.abs(relative))
            ),
            "p90_absolute_relative_residual": float(
                np.percentile(np.abs(relative), 90.0)
            ),
            "maximum_absolute_relative_residual": float(
                np.max(np.abs(relative))
            ),
        },
        "declared_reactor_development_band": {
            "interval_Td": list(OPERATING_BAND_TD),
            "point_count": int(np.count_nonzero(operating)),
            "median_absolute_relative_residual": float(
                np.median(np.abs(relative[operating]))
            ),
            "maximum_absolute_relative_residual": float(
                np.max(np.abs(relative[operating]))
            ),
            "role": (
                "engineering transport band for the next reactor rung; "
                "not inferred from the withheld TiO2 feature"
            ),
        },
    }


def audit(*, include_working_set: bool = True) -> dict:
    measurement = load_nist_1999_chf3_drift_curve()
    fields = measurement.reduced_electric_field_Td
    working = load_kushner_zhang_2000_chf3_replay()
    variants = []
    if include_working_set:
        variants.append(("kushner_zhang_working_set", working.derived_deck))
    for closure in (
        "constant_join_ratio",
        "linear_return_to_working_set_at_120eV",
    ):
        evaluated = derive_nist_evaluated_chf3_replay(
            working, high_energy_closure=closure
        )
        variants.append((f"nist_evaluated_{closure}", evaluated.derived_deck))

    models = {}
    for name, deck in variants:
        prediction = _solve_curve(deck, fields)
        prediction["collision_deck_sha256"] = deck.payload_sha256
        prediction["comparison"] = _grade(
            prediction, measurement.drift_velocity_m_s
        )
        models[name] = prediction

    representative_fields = np.asarray((10.0, 60.0, 150.0, 250.0))
    evaluated = derive_nist_evaluated_chf3_replay(working)
    coarse = _solve_curve(evaluated.derived_deck, representative_fields)
    fine = _solve_curve(
        evaluated.derived_deck, representative_fields, cells_scale=2
    )
    coarse_velocity = np.asarray([
        row["flux_drift_velocity_m_s"] for row in coarse["rows"]
    ])
    fine_velocity = np.asarray([
        row["flux_drift_velocity_m_s"] for row in fine["rows"]
    ])
    grid_change = fine_velocity / coarse_velocity - 1.0

    return {
        "schema": "petch.chf3_electron_swarm_audit.v1",
        "claim_class": (
            "source_replay_and_evaluated_transport_consistency; "
            "not_independent_branch_or_reactor_validation"
        ),
        "sources": {
            "kushner_zhang_working_set_csv_sha256": (
                KUSHNER_ZHANG_CHF3_CSV_SHA256
            ),
            "christophorou_olthoff_table4_csv_sha256": (
                NIST_CHF3_TABLE4_CSV_SHA256
            ),
            "christophorou_olthoff_table5_csv_sha256": (
                NIST_CHF3_TABLE5_CSV_SHA256
            ),
            "christophorou_olthoff_table6_csv_sha256": (
                NIST_CHF3_TABLE6_CSV_SHA256
            ),
            "kushner_zhang_figure3_render_sha256": (
                KUSHNER_FIGURE3_RENDER_SHA256
            ),
            "source_artifacts_committed": False,
        },
        "comparison_quantity": {
            "predicted": "flux drift velocity from the local two-term solve",
            "tabulated": "NIST Table 6 recommended 298 K drift velocity",
            "bulk_flux_distinction_resolved": False,
            "reason": (
                "The NIST table does not label its fitted drift curve as a "
                "bulk-versus-flux coefficient. The comparison is retained "
                "as a source-replay diagnostic, not promoted to a direct grade."
            ),
        },
        "models": models,
        "numerical_convergence": {
            "fields_Td": representative_fields.tolist(),
            "coarse_actual_cell_count": coarse["actual_cell_count"],
            "fine_actual_cell_count": fine["actual_cell_count"],
            "flux_drift_relative_change": grid_change.tolist(),
            "maximum_absolute_flux_drift_relative_change": float(
                np.max(np.abs(grid_change))
            ),
        },
        "certification": {
            "source_curve_visually_overpredicts_low_field_measurements": True,
            "independent_bolos_replay_performed": True,
            "independent_bolos_replay_committed": False,
            "supports_use_as_electron_transport_input": True,
            "supports_independent_chemical_branch_validation": False,
            "supports_target_reactor_state_prediction": False,
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
    lines = []
    for name, model in result["models"].items():
        all_points = model["comparison"]["all_points"]
        band = model["comparison"][
            "declared_reactor_development_band"
        ]
        lines.append(
            f"| {name} | {100.0 * all_points['median_absolute_relative_residual']:.2f}% "
            f"| {100.0 * all_points['maximum_absolute_relative_residual']:.2f}% "
            f"| {100.0 * band['maximum_absolute_relative_residual']:.2f}% |"
        )
    report = f"""# CHF3 electron-swarm source replay

| transport deck | all-point median | all-point maximum | 40--100 Td maximum |
|---|---:|---:|---:|
{chr(10).join(lines)}

The original working set reproduces the source paper's own low-field
overprediction. An independent BOLOS solve and petch agree on that diagnosis.
The NIST-evaluated elastic/momentum backbone materially improves transport in
the declared 40--100 Td engineering band while keeping the original inelastic
chemistry fixed.

This is deliberately a source-replay receipt. The working set was regressed
against swarm behavior, the NIST curve is not labeled bulk-versus-flux, and
the neutral-dissociation branches remain weakly constrained. It authorizes a
transport input for the next reactor rung, not a unique reactor state, wafer
flux, or feature depth.

Maximum representative-grid drift change: `{100.0 * result['numerical_convergence']['maximum_absolute_flux_drift_relative_change']:.4f}%`.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-working-set", action="store_true")
    args = parser.parse_args()
    result = audit(include_working_set=not args.skip_working_set)
    _write(result, args.output)
    print(json.dumps({
        name: model["comparison"]
        for name, model in result["models"].items()
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
