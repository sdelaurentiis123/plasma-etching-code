#!/usr/bin/env python3
"""Audit the measured-feed electron-kinetic rung for the blind NPG80 case.

This closes recipe gas fractions and the 13.56 MHz/30 mTorr reduced frequency,
but deliberately scans bulk reduced field: forward generator power is not a
bulk E/N measurement.  No TiO2 outcome enters this calculation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from petch.reactor_global.chf3_electron_collisions import (
    derive_nist_evaluated_chf3_replay,
    load_kushner_zhang_2000_chf3_replay,
)
from petch.reactor_global.electron_collision_mixture import (
    compose_electron_collision_decks,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
)
from petch.reactor_global.o2_electron_collisions import (
    SONG_2026_O2_ARTICLE_DOI,
    SONG_2026_O2_DATASET_DOI,
    SONG_2026_O2_WORKBOOK_SHA256,
    load_song_2026_o2_replay,
)
from petch.reactor_global.sf6_electron_collisions import (
    derive_nist_evaluated_sf6_replay,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_feed_electron_kinetics_v1"
)
PRESSURE_PA = 3.99967104
GAS_TEMPERATURE_K = 293.15
RF_FREQUENCY_HZ = 13.56e6
BOLTZMANN_J_K = 1.380649e-23
FEED_SCCM = {"CHF3": 55.0, "SF6": 5.0, "O2": 1.0}
FEED_FRACTIONS = {
    species: flow / sum(FEED_SCCM.values())
    for species, flow in FEED_SCCM.items()
}
CANDIDATE_BULK_FIELDS_TD = (40.0, 60.0, 80.0, 100.0)


def _grid(deck, *, scale: int = 5) -> ElectronEnergyGrid:
    if scale not in {4, 5, 6}:
        raise ValueError("electron grid scale must be 4, 5, or 6")
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV and process.energy_loss_eV < 120.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, .0001, .001, .01, .1, 1.0, 10.0, 40.0, 120.0),
        (
            8 * scale,
            12 * scale,
            18 * scale,
            36 * scale,
            48 + 12 * (scale - 4),
            96 + 24 * (scale - 4),
            96 + 24 * (scale - 4),
            120 + 30 * (scale - 4),
        ),
        inserted_boundaries_eV=thresholds,
    )


def _collision_rows(solution) -> list[dict]:
    rows = []
    for moment in solution.collision_moments:
        if moment.process_kind not in {"IONIZATION", "ATTACHMENT"}:
            continue
        fraction = FEED_FRACTIONS[moment.target]
        electron_change = 1 if moment.process_kind == "IONIZATION" else -1
        rows.append({
            "target": moment.target,
            "kind": moment.process_kind,
            "product": moment.product,
            "rate_coefficient_m3_s": float(moment.rate_coefficient_m3_s),
            "feed_fraction_weighted_rate_coefficient_m3_s": float(
                fraction * moment.rate_coefficient_m3_s
            ),
            "electron_source_coefficient_m3_s": float(
                electron_change * fraction * moment.rate_coefficient_m3_s
            ),
        })
    return rows


def _solve(
    deck,
    fields_Td,
    *,
    angular_frequency_over_density_m3_s: float,
    grid_scale: int = 5,
) -> dict:
    grid = _grid(deck, scale=grid_scale)
    solver = DeterministicTwoTermBoltzmannSolver(grid, deck)
    rows = []
    for field in fields_Td:
        solution = solver.solve(
            TwoTermBoltzmannCondition(
                reduced_electric_field_Td=float(field),
                gas_temperature_K=GAS_TEMPERATURE_K,
                target_mole_fractions=FEED_FRACTIONS,
                growth_model="temporal_growth",
                initial_electron_temperature_eV=float(
                    np.clip(.012 * float(field), .3, 2.5)
                ),
                angular_field_frequency_over_density_m3_s=(
                    angular_frequency_over_density_m3_s
                ),
            ),
            relative_tolerance=3.0e-6,
            maximum_tail_population_fraction=2.0e-6,
        )
        collision_rows = _collision_rows(solution)
        reconstructed_growth = sum(
            item["electron_source_coefficient_m3_s"]
            for item in collision_rows
        )
        flux_drift = (
            solution.transport_moments
            .flux_reduced_mobility_m_inv_V_inv_s_inv
            * float(field) * 1.0e-21
        )
        rows.append({
            "reduced_electric_field_Td": float(field),
            "mean_electron_energy_eV": float(
                solution.distribution.mean_energy_eV
            ),
            "flux_drift_velocity_m_s": float(flux_drift),
            "net_growth_rate_coefficient_m3_s": float(
                solution.net_growth_rate_coefficient_m3_s
            ),
            "reconstructed_net_growth_rate_coefficient_m3_s": float(
                reconstructed_growth
            ),
            "net_growth_reconstruction_error_m3_s": float(
                reconstructed_growth
                - solution.net_growth_rate_coefficient_m3_s
            ),
            "reduced_field_power_gain_eV_m3_s": float(
                solution.transport_moments.reduced_field_power_gain_eV_m3_s
            ),
            "weighted_iteration_residual": float(
                solution.weighted_iteration_residual
            ),
            "iteration_count": solution.iteration_count,
            "growth_root_evaluations": solution.growth_root_evaluations,
            "electron_creation_and_attachment": collision_rows,
        })
    return {
        "grid_scale": grid_scale,
        "actual_cell_count": grid.cell_count,
        "angular_frequency_over_density_m3_s": (
            angular_frequency_over_density_m3_s
        ),
        "rows": rows,
    }


def _primary_decks(
    workbook: Path,
    *,
    chf3_closure: str = "constant_join_ratio",
    sf6_vibrational_loss_eV: float = .095,
    o2_tail: str = "inverse_energy",
    o2_dissociation_onset: str = "zero_until_first_tabulated_energy",
):
    chf3 = derive_nist_evaluated_chf3_replay(
        load_kushner_zhang_2000_chf3_replay(),
        high_energy_closure=chf3_closure,
    )
    sf6 = derive_nist_evaluated_sf6_replay(
        vibrational_energy_loss_eV=sf6_vibrational_loss_eV,
    )
    o2 = load_song_2026_o2_replay(
        workbook,
        high_energy_tail_closure=o2_tail,
        dissociation_onset_closure=o2_dissociation_onset,
    )
    mixed = compose_electron_collision_decks(
        (chf3.derived_deck, sf6.derived_deck, o2.derived_deck),
        retrieved_at="2026-08-18",
        mixture_name="Zhu NPG80 55 CHF3 / 5 SF6 / 1 O2 feed",
    )
    return chf3, sf6, o2, mixed


def _scalar_rows(result: dict) -> list[dict]:
    keys = (
        "reduced_electric_field_Td",
        "mean_electron_energy_eV",
        "flux_drift_velocity_m_s",
        "net_growth_rate_coefficient_m3_s",
        "reduced_field_power_gain_eV_m3_s",
    )
    return [{key: row[key] for key in keys} for row in result["rows"]]


def _relative_change(candidate: dict, nominal: dict) -> dict:
    candidate_rows = _scalar_rows(candidate)
    nominal_rows = _scalar_rows(nominal)
    if [row["reduced_electric_field_Td"] for row in candidate_rows] != [
        row["reduced_electric_field_Td"] for row in nominal_rows
    ]:
        raise RuntimeError("electron sensitivity field grids differ")
    quantities = (
        "mean_electron_energy_eV",
        "flux_drift_velocity_m_s",
        "net_growth_rate_coefficient_m3_s",
        "reduced_field_power_gain_eV_m3_s",
    )
    changes = {}
    for quantity in quantities:
        base = np.asarray([row[quantity] for row in nominal_rows])
        other = np.asarray([row[quantity] for row in candidate_rows])
        relative = other / base - 1.0
        changes[quantity] = {
            "relative_change": relative.tolist(),
            "maximum_absolute_relative_change": float(
                np.max(np.abs(relative))
            ),
        }
    return changes


def audit(workbook: Path) -> dict:
    number_density = PRESSURE_PA / (BOLTZMANN_J_K * GAS_TEMPERATURE_K)
    angular_frequency = 2.0 * math.pi * RF_FREQUENCY_HZ
    reduced_frequency = angular_frequency / number_density
    chf3, sf6, o2, nominal_deck = _primary_decks(workbook)
    nominal_rf = _solve(
        nominal_deck,
        CANDIDATE_BULK_FIELDS_TD,
        angular_frequency_over_density_m3_s=reduced_frequency,
    )
    representative = (60.0, 100.0)
    nominal_rf_representative = {
        **nominal_rf,
        "rows": [
            row for row in nominal_rf["rows"]
            if row["reduced_electric_field_Td"] in representative
        ],
    }
    dc = _solve(
        nominal_deck,
        representative,
        angular_frequency_over_density_m3_s=0.0,
    )
    sensitivities = {}
    variants = (
        (
            "chf3_high_energy_tail",
            {"chf3_closure": "linear_return_to_working_set_at_120eV"},
        ),
        (
            "sf6_vibrational_loss",
            {"sf6_vibrational_loss_eV": .117},
        ),
        (
            "o2_high_energy_tail",
            {"o2_tail": "linear_to_zero_at_30eV"},
        ),
        (
            "o2_dissociation_onset",
            {"o2_dissociation_onset": "linear_from_physical_threshold"},
        ),
    )
    for name, arguments in variants:
        *_, variant_deck = _primary_decks(workbook, **arguments)
        prediction = _solve(
            variant_deck,
            representative,
            angular_frequency_over_density_m3_s=reduced_frequency,
        )
        sensitivities[name] = {
            "variant_arguments": arguments,
            "prediction": _scalar_rows(prediction),
            "relative_to_nominal": _relative_change(
                prediction, nominal_rf_representative
            ),
        }
    coarse_failure = None
    try:
        coarse = _solve(
            nominal_deck,
            representative,
            angular_frequency_over_density_m3_s=reduced_frequency,
            grid_scale=4,
        )
    except (FloatingPointError, RuntimeError, ValueError) as exc:
        coarse = None
        coarse_failure = f"{type(exc).__name__}: {exc}"
    fine = _solve(
        nominal_deck,
        representative,
        angular_frequency_over_density_m3_s=reduced_frequency,
        grid_scale=6,
    )
    return {
        "schema": "petch.zhu_npg80_feed_electron_kinetics_audit.v1",
        "condition_id": "zhu-2026-npg80-tio2-chf3-sf6-o2-20min",
        "claim_class": (
            "measured_recipe_feed_and_frequency_conditioned_local_electron_"
            "kinetics; bulk_field_scanned_not_inferred"
        ),
        "source": {
            "song_2026_o2_workbook_sha256": SONG_2026_O2_WORKBOOK_SHA256,
            "song_2026_o2_article_doi": SONG_2026_O2_ARTICLE_DOI,
            "song_2026_o2_dataset_doi": SONG_2026_O2_DATASET_DOI,
            "song_2026_o2_license": "CC BY-NC 4.0",
            "song_2026_o2_source_artifact_committed": False,
            "chf3_derived_deck_sha256": chf3.derived_deck.payload_sha256,
            "sf6_derived_deck_sha256": sf6.derived_deck.payload_sha256,
            "o2_derived_deck_sha256": o2.derived_deck.payload_sha256,
            "mixture_deck_sha256": nominal_deck.payload_sha256,
        },
        "recipe_condition": {
            "feed_sccm": FEED_SCCM,
            "feed_mole_fractions": FEED_FRACTIONS,
            "pressure_Pa": PRESSURE_PA,
            "gas_temperature_K_assumption": GAS_TEMPERATURE_K,
            "gas_number_density_m3": number_density,
            "rf_frequency_Hz": RF_FREQUENCY_HZ,
            "angular_frequency_rad_s": angular_frequency,
            "angular_frequency_over_density_m3_s": reduced_frequency,
            "candidate_bulk_reduced_field_Td": list(
                CANDIDATE_BULK_FIELDS_TD
            ),
            "candidate_field_status": (
                "development scan; not a recipe measurement and not inferred "
                "from the withheld SEM"
            ),
        },
        "nominal_finite_frequency_prediction": nominal_rf,
        "dc_vs_rf_diagnostic": {
            "dc_prediction": _scalar_rows(dc),
            "rf_prediction": _scalar_rows(nominal_rf_representative),
            "relative_dc_to_rf": _relative_change(
                dc, nominal_rf_representative
            ),
            "interpretation": (
                "DC is a diagnostic only; the recipe is 13.56 MHz and the "
                "finite-frequency operator is the nominal model."
            ),
        },
        "cross_section_closure_sensitivities": sensitivities,
        "numerical_convergence": {
            "coarse_cell_count": (
                _grid(nominal_deck, scale=4).cell_count
            ),
            "nominal_cell_count": nominal_rf["actual_cell_count"],
            "fine_cell_count": fine["actual_cell_count"],
            "coarse_failure": coarse_failure,
            "coarse_relative_to_nominal": (
                None
                if coarse is None
                else _relative_change(coarse, nominal_rf_representative)
            ),
            "fine_relative_to_nominal": _relative_change(
                fine, nominal_rf_representative
            ),
            "interpretation": (
                "The coarse grid is deliberately retained as a rejection "
                "test because the SF6 attachment resonance creates a "
                "low-energy boundary layer."
            ),
        },
        "physics_findings": {
            "feed_is_electronegative_in_scanned_local_state": all(
                row["net_growth_rate_coefficient_m3_s"] < 0.0
                for row in nominal_rf["rows"]
            ),
            "forward_power_does_not_close_bulk_field": True,
            "electron_density_not_closed_by_eedf_alone": True,
            "requires_global_particle_and_power_balance_next": True,
            "feature_depth_used": False,
        },
        "certification": {
            "supports_species_resolved_feed_eedf_input": True,
            "supports_unique_bulk_field": False,
            "supports_unique_electron_density": False,
            "supports_unique_wafer_flux": False,
            "supports_feature_depth": False,
        },
    }


def _write(result: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _check(output: Path) -> None:
    path = output / "audit.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        payload.get("schema")
        != "petch.zhu_npg80_feed_electron_kinetics_audit.v1"
        or payload.get("source", {}).get("song_2026_o2_workbook_sha256")
        != SONG_2026_O2_WORKBOOK_SHA256
        or payload.get("source", {}).get(
            "song_2026_o2_source_artifact_committed"
        ) is not False
        or payload.get("physics_findings", {}).get("feature_depth_used")
        is not False
    ):
        raise RuntimeError("committed NPG80 electron-kinetic audit is stale")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-workbook", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        _check(args.output)
        return
    if args.source_workbook is None:
        parser.error("--source-workbook is required unless --check is used")
    _write(audit(args.source_workbook), args.output)


if __name__ == "__main__":
    main()
