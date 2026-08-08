#!/usr/bin/env python3
"""Train one Lam power fraction at 300 W and forecast 500 W out of sample."""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from petch.reactor_global import (
    DeterministicTwoTermBoltzmannSolver,
    EEDFChlorineAbsorbedPowerModel,
    ElectronEnergyGrid,
    HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256,
    ReactorDiagnosticConditionedPowerFraction,
    ReactionNetwork,
    build_lee_lieberman_chlorine_particle_network,
    load_legacy_siglo_hamilton_comsol_chlorine_replay,
)
try:
    from scripts.audit_malyshev_1998_eedf_source_replay import (
        _condition,
        _providers,
        _reference,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...py`` execution.
    from audit_malyshev_1998_eedf_source_replay import (
        _condition,
        _providers,
        _reference,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"
JSON_NAME = "malyshev_1998_power_fraction_transfer.json"
CSV_NAME = "malyshev_1998_power_fraction_transfer.csv"
REPORT_NAME = "MALYSHEV_1998_POWER_FRACTION_TRANSFER.md"
TRAINING_POWER_W = 300.0
HELD_OUT_POWER_W = 500.0
FRACTION_BRACKET = (0.30, 0.50)


def _build_model(collision_deck: Path, atomic_cl_momentum: Path):
    replay = load_legacy_siglo_hamilton_comsol_chlorine_replay(
        collision_deck,
        atomic_cl_momentum,
        maximum_energy_eV=200.0,
    )
    thresholds = tuple(sorted({
        float(process.energy_loss_eV)
        for process in replay.derived_deck.processes
        if process.energy_loss_eV is not None
        and 0.0 < float(process.energy_loss_eV) < 200.0
    }))
    grid = ElectronEnergyGrid.piecewise_linear(
        (0.0, 0.5, 5.0, 20.0, 80.0, 200.0),
        (80, 120, 100, 60, 40),
        inserted_boundaries_eV=thresholds,
    )
    lee = build_lee_lieberman_chlorine_particle_network()
    heavy = ReactionNetwork(species=lee.species, reactions=lee.reactions[6:8])
    model = EEDFChlorineAbsorbedPowerModel(
        DeterministicTwoTermBoltzmannSolver(grid, replay.derived_deck),
        replay.collision_chemistry,
        heavy,
        electron_electron_coulomb_model="isotropic_classical_debye",
    )
    return replay, grid, model


def _solve(
    model,
    *,
    source_power_W: float,
    absorbed_fraction: float,
    seed=None,
    conditioned_boundary: ReactorDiagnosticConditionedPowerFraction | None = None,
):
    charged, neutral, wall_energy = _providers()
    condition = _condition(source_power_W, absorbed_fraction)
    if conditioned_boundary is not None:
        condition = replace(
            condition,
            absorbed_power=conditioned_boundary.estimate(
                source_power_W,
                source=(
                    "Malyshev 1998 forward TCP power into matching network; "
                    "inline power meters"
                ),
            ),
        )
    seed_arguments = {} if seed is None else {
        "initial_densities_m3": seed.densities_m3,
        "initial_exhaust_loss_frequency_s_inv": (
            seed.exhaust_loss_frequency_s_inv),
        "initial_reduced_electric_field_Td": seed.reduced_electric_field_Td,
    }
    return model.solve(
        condition,
        charged_transport_provider=charged,
        neutral_wall_transport_provider=neutral,
        wall_energy_provider=wall_energy,
        residual_tolerance=2.0e-7,
        maximum_evaluations=1200,
        maximum_tail_population_fraction=1.0e-6,
        **seed_arguments,
    )


def _row(solution, reference, *, role: str, absorbed_fraction: float):
    densities = solution.densities_m3
    cl2_proxy = 100.0 * densities["Cl2"] / (
        densities["Cl2"] + densities["Cl"])
    density_reference = float(reference["electron_density_m3"])
    temperature_reference = float(reference["electron_temperature_eV"])
    temperature_proxy = 2.0 / 3.0 * solution.mean_electron_energy_eV
    cl2_reference = reference["relative_cl2_density_percent"]
    total_flux = sum(solution.axial_positive_ion_flux_m2_s.values())
    electron = solution.electron_solution
    return {
        "validation_role": role,
        "source_power_W": (
            TRAINING_POWER_W
            if role == "calibration_training"
            else HELD_OUT_POWER_W),
        "absorbed_fraction": absorbed_fraction,
        "absorbed_power_W": absorbed_fraction * (
            TRAINING_POWER_W
            if role == "calibration_training"
            else HELD_OUT_POWER_W),
        "reduced_electric_field_Td": solution.reduced_electric_field_Td,
        "electron_density_m3": densities["e"],
        "measured_volume_average_electron_density_m3": density_reference,
        "electron_density_percent_error": 100.0 * (
            densities["e"] / density_reference - 1.0),
        "mean_electron_energy_eV": solution.mean_electron_energy_eV,
        "mean_energy_temperature_proxy_eV": temperature_proxy,
        "measured_oes_electron_temperature_eV": temperature_reference,
        "temperature_proxy_percent_error": 100.0 * (
            temperature_proxy / temperature_reference - 1.0),
        "modeled_relative_cl2_density_percent_proxy": cl2_proxy,
        "measured_relative_cl2_density_percent": cl2_reference,
        "relative_cl2_proxy_error_percentage_point": (
            None if cl2_reference is None else cl2_proxy - cl2_reference),
        "total_positive_ion_axial_flux_m2_s": total_flux,
        "cl2plus_axial_flux_m2_s": (
            solution.axial_positive_ion_flux_m2_s["Cl2+"]),
        "clplus_axial_flux_m2_s": (
            solution.axial_positive_ion_flux_m2_s["Cl+"]),
        "clplus_ion_fraction": (
            solution.axial_positive_ion_flux_m2_s["Cl+"] / total_flux),
        "maximum_normalized_residual": solution.maximum_normalized_residual,
        "solver_evaluations": solution.solver_evaluations,
        "coulomb_logarithm": electron.coulomb_logarithm,
        "coulomb_nonlinear_iterations": electron.nonlinear_iteration_count,
        "electron_growth_root_evaluations": (
            electron.growth_root_evaluations),
        "supports_absorbed_power_measurement": False,
        "supports_wafer_flux": False,
        "supports_feature_depth": False,
    }


def audit(collision_deck: Path, atomic_cl_momentum: Path):
    replay, grid, model = _build_model(collision_deck, atomic_cl_momentum)
    training_reference = _reference(TRAINING_POWER_W)
    target_density = float(training_reference["electron_density_m3"])
    solutions = {}
    trace = []

    def residual(fraction: float) -> float:
        fraction = float(fraction)
        if fraction in solutions:
            solution = solutions[fraction]
        else:
            seed = None
            if solutions:
                nearest = min(solutions, key=lambda value: abs(value - fraction))
                seed = solutions[nearest]
            solution = _solve(
                model,
                source_power_W=TRAINING_POWER_W,
                absorbed_fraction=fraction,
                seed=seed,
            )
            solutions[fraction] = solution
            trace.append({
                "absorbed_fraction": fraction,
                "absorbed_power_W": fraction * TRAINING_POWER_W,
                "electron_density_m3": solution.densities_m3["e"],
                "log_density_residual": float(np.log(
                    solution.densities_m3["e"] / target_density)),
                "solver_evaluations": solution.solver_evaluations,
                "maximum_normalized_residual": (
                    solution.maximum_normalized_residual),
            })
        return float(np.log(solution.densities_m3["e"] / target_density))

    lower_residual = residual(FRACTION_BRACKET[0])
    upper_residual = residual(FRACTION_BRACKET[1])
    if lower_residual * upper_residual >= 0.0:
        raise RuntimeError("training density does not bracket a power fraction")
    fitted_fraction, root = brentq(
        residual,
        *FRACTION_BRACKET,
        xtol=1.0e-7,
        rtol=1.0e-12,
        maxiter=20,
        full_output=True,
        disp=True,
    )
    fitted_fraction = float(fitted_fraction)
    training_solution = solutions.get(fitted_fraction)
    if training_solution is None:
        nearest = min(solutions, key=lambda value: abs(value - fitted_fraction))
        training_solution = _solve(
            model,
            source_power_W=TRAINING_POWER_W,
            absorbed_fraction=fitted_fraction,
            seed=solutions[nearest],
        )
        solutions[fitted_fraction] = training_solution
        trace.append({
            "absorbed_fraction": fitted_fraction,
            "absorbed_power_W": fitted_fraction * TRAINING_POWER_W,
            "electron_density_m3": training_solution.densities_m3["e"],
            "log_density_residual": float(np.log(
                training_solution.densities_m3["e"] / target_density)),
            "solver_evaluations": training_solution.solver_evaluations,
            "maximum_normalized_residual": (
                training_solution.maximum_normalized_residual),
        })

    boundary = ReactorDiagnosticConditionedPowerFraction(
        absorbed_fraction=fitted_fraction,
        calibration_condition_id=(
            "malyshev-gap6.5cm-pressure2mTorr-source300W"),
        calibration_observable="volume-average electron density",
        calibration_source=(
            "Malyshev et al. 1998 Figure 11; measured along 1.35-cm line and "
            "converted by the authors' reported spatial assumption"
        ),
    )
    held_out_solution = _solve(
        model,
        source_power_W=HELD_OUT_POWER_W,
        absorbed_fraction=fitted_fraction,
        seed=training_solution,
        conditioned_boundary=boundary,
    )
    rows = [
        _row(
            training_solution,
            training_reference,
            role="calibration_training",
            absorbed_fraction=fitted_fraction,
        ),
        _row(
            held_out_solution,
            _reference(HELD_OUT_POWER_W),
            role="held_out_reactor_diagnostic_forecast",
            absorbed_fraction=fitted_fraction,
        ),
    ]
    return {
        "schema": "petch.malyshev_1998_power_fraction_transfer.v1",
        "claim_class": (
            "reactor-diagnostic-conditioned equipment transfer; not an "
            "absorbed-power measurement or knobs-to-wafer validation"),
        "model_variant": (
            "legacy_siglo_hamilton_plus_comsol_nist_atomic_cl_plus_"
            "isotropic_ee"),
        "collision_identity": {
            "raw_collision_payload_sha256": (
                replay.molecular_replay.raw_payload_sha256),
            "derived_collision_deck_sha256": replay.derived_deck.payload_sha256,
            "atomic_momentum_payload_sha256": (
                replay.atomic_momentum_payload_sha256),
            "hamilton_state_cross_sections_sha256": (
                HAMILTON_2018_CL2_STATE_CROSS_SECTIONS_SHA256),
            "raw_collision_bytes_committed": False,
        },
        "energy_grid": {
            "family": "threshold_aligned_piecewise_linear_v1",
            "actual_cell_count": grid.cell_count,
            "maximum_energy_eV": grid.boundaries_eV[-1],
        },
        "calibration": {
            "model": "constant_source_to_plasma_absorbed_fraction",
            "training_power_W": TRAINING_POWER_W,
            "training_observable": "volume_average_electron_density_m3",
            "fraction_search_bracket": list(FRACTION_BRACKET),
            "fitted_absorbed_fraction": fitted_fraction,
            "fitted_absorbed_power_W": fitted_fraction * TRAINING_POWER_W,
            "root_converged": bool(root.converged),
            "root_iterations": int(root.iterations),
            "training_trace": trace,
            "temperature_used_for_selection": False,
            "dissociation_used_for_selection": False,
            "held_out_500W_used_for_selection": False,
            "feature_depth_used_for_selection": False,
        },
        "transfer": {
            "held_out_power_W": HELD_OUT_POWER_W,
            "absorbed_fraction_frozen_from_training": fitted_fraction,
            "formal_pass_threshold": None,
            "reason_no_formal_threshold": (
                "Malyshev reports no electron-density uncertainty in this "
                "article and the diagnostic does not directly measure "
                "absorbed plasma power"),
        },
        "evidence_boundary": {
            "power": (
                "forward TCP power was measured into the matching network; "
                "the fitted fraction is an effective reactor-model closure, "
                "not measured coil/window/plasma absorption"),
            "temperature": (
                "2/3 mean electron energy is a proxy, not the OES forward "
                "observable"),
            "relative_cl2": (
                "model nCl2/(nCl2+nCl) is not Malyshev's rare-gas-reduced "
                "absolute Cl2 observable"),
            "flux": (
                "global axial positive-ion flux is not a local wafer flux "
                "and has no IED/IAD"),
        },
        "supports_absorbed_power_measurement": False,
        "supports_reactor_state_prediction": False,
        "supports_wafer_flux": False,
        "supports_feature_depth": False,
        "rows": rows,
    }


def _write(result, output_directory: Path):
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / JSON_NAME).write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with (output_directory / CSV_NAME).open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(result["rows"][0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(result["rows"])
    training, held_out = result["rows"]
    report = f"""# Malyshev 1998 constant power-fraction transfer

## Verdict

One effective source-to-plasma fraction, `{result['calibration']['fitted_absorbed_fraction']:.6f}`,
was inferred from the 300 W volume-average electron density only. With that
fraction frozen, the untouched 500 W density error is
`{held_out['electron_density_percent_error']:+.2f}%`. This closes the power-scaling
residual numerically, but it is **not a direct absorbed-power validation**:
Malyshev reports forward TCP power into the matching network and no density
uncertainty in this article.

| role | source W | absorbed W | ne error | energy-proxy error | Cl2-proxy error | axial ion flux m-2 s-1 |
|---|---:|---:|---:|---:|---:|---:|
| training | {training['source_power_W']:.0f} | {training['absorbed_power_W']:.2f} | {training['electron_density_percent_error']:+.2f}% | {training['temperature_proxy_percent_error']:+.2f}% | n/a | {training['total_positive_ion_axial_flux_m2_s']:.3e} |
| held out | {held_out['source_power_W']:.0f} | {held_out['absorbed_power_W']:.2f} | {held_out['electron_density_percent_error']:+.2f}% | {held_out['temperature_proxy_percent_error']:+.2f}% | {held_out['relative_cl2_proxy_error_percentage_point']:+.2f} pp | {held_out['total_positive_ion_axial_flux_m2_s']:.3e} |

## Interpretation boundary

- Neither temperature, dissociation, the held-out 500 W condition, nor any
  feature depth selected the fraction.
- The held-out density transfer is descriptive because the source provides no
  uncertainty for that observable here.
- The remaining temperature and Cl2 residuals show that a power fraction alone
  does not close the reactor physics.
- The axial ion flux is still volume-model output, not a wafer boundary; no
  species-resolved IED/IAD exists yet.
- Every prediction/depth support flag remains false.
"""
    (output_directory / REPORT_NAME).write_text(report, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("collision_deck", type=Path)
    parser.add_argument("atomic_cl_momentum", type=Path)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    arguments = parser.parse_args()
    result = audit(arguments.collision_deck, arguments.atomic_cl_momentum)
    if not arguments.no_write:
        _write(result, arguments.output_directory)
    if arguments.summary_only:
        print(
            f"fraction={result['calibration']['fitted_absorbed_fraction']:.6f}"
        )
        for row in result["rows"]:
            print(
                f"{row['validation_role']} {row['source_power_W']:.0f} W: "
                f"ne={row['electron_density_percent_error']:+.2f}% "
                f"energy={row['temperature_proxy_percent_error']:+.2f}% "
                f"Cl2={row['relative_cl2_proxy_error_percentage_point']}"
            )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
