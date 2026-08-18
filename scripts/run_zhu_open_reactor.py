#!/usr/bin/env python3
"""Run one target-free NPG80 open-reactor sensitivity state.

The supplied recipe fixes flow, pressure, frequency, and forward-power demand.
It does not fix absorbed power or the geometric/wall/sheath closures exposed
by this command.  The defaults are a central development state, not a
certified machine identification and not a fit to the withheld SEM.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping

from petch.reactor_global import (
    CylindricalReactor,
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    ZhuOpenReactorCondition,
    ZhuOpenReactorModel,
    build_zhu_parent_collision_chemistry,
    build_zhu_supplemental_chemistry,
    standard_volume_flow_molecules_s,
)


FEED_SCCM = {"CHF3": 55.0, "SF6": 5.0, "O2": 1.0}


def _grid(deck) -> ElectronEnergyGrid:
    thresholds = tuple(sorted({
        process.energy_loss_eV
        for process in deck.processes
        if process.energy_loss_eV and process.energy_loss_eV < 120.0
    }))
    return ElectronEnergyGrid.piecewise_linear(
        (0.0, .0001, .001, .01, .1, 1.0, 10.0, 40.0, 120.0),
        (40, 60, 90, 48, 60, 120, 120, 150),
        inserted_boundaries_eV=thresholds,
    )


def _feed_molecules_s() -> dict[str, float]:
    return {
        species: standard_volume_flow_molecules_s(
            flow,
            standard_temperature_K=273.15,
            standard_pressure_Pa=101325.0,
        )
        for species, flow in FEED_SCCM.items()
    }


def _continuation_state(
    path: Path,
    *,
    model: ZhuOpenReactorModel,
    condition: ZhuOpenReactorCondition,
) -> tuple[dict[str, float], float, float]:
    """Lift a prior solved state onto an expanded conserved species basis."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    state = payload.get("state")
    if not isinstance(state, Mapping):
        raise ValueError("continuation JSON has no state mapping")
    previous = state.get("densities_m3")
    if not isinstance(previous, Mapping):
        raise ValueError("continuation JSON has no state.densities_m3 mapping")
    usable = {
        str(name): float(value)
        for name, value in previous.items()
        if isinstance(value, (int, float))
        and math.isfinite(float(value))
        and float(value) > 0.0
    }
    target_density = condition.target_neutral_density_m3
    floor = max(1.0e4, 1.0e-16 * target_density)
    densities = {
        name: usable.get(name, floor) for name in model.species_order
    }
    new_neutrals = tuple(
        name for name in model.neutral_names if name not in usable)
    old_neutrals = tuple(
        name for name in model.neutral_names if name in usable)
    if new_neutrals:
        # Log-density derivatives disappear if newly added daughter species
        # begin at the numerical floor.  Reserve a small, pressure-neutral
        # seed inventory and rescale the converged old neutral composition.
        reserve = min(0.10, 0.005 * len(new_neutrals))
        old_total = sum(densities[name] for name in old_neutrals)
        if old_total <= 0.0:
            raise ValueError("continuation state has no usable neutral density")
        for name in old_neutrals:
            densities[name] *= (1.0 - reserve) * target_density / old_total
        for name in new_neutrals:
            densities[name] = reserve * target_density / len(new_neutrals)
    else:
        total = sum(densities[name] for name in model.neutral_names)
        for name in model.neutral_names:
            densities[name] *= target_density / total
    exhaust = float(state.get("exhaust_loss_frequency_s_inv", float("nan")))
    field = float(state.get("reduced_electric_field_Td", float("nan")))
    if not math.isfinite(exhaust) or exhaust <= 0.0:
        raise ValueError("continuation state has no valid exhaust frequency")
    if not math.isfinite(field) or field <= 0.0:
        raise ValueError("continuation state has no valid reduced field")
    return densities, exhaust, field


def run(args) -> dict:
    parent = build_zhu_parent_collision_chemistry(args.source_workbook)
    supplemental = build_zhu_supplemental_chemistry(
        kokkoris_eedf_shape=args.kokkoris_eedf_shape,
        chf3_f_rate_branch=args.chf3_f_rate_branch,
    )
    solver = DeterministicTwoTermBoltzmannSolver(
        _grid(parent.mixed_deck), parent.mixed_deck)
    geometry = CylindricalReactor(
        radius_m=.5e-3 * args.electrode_diameter_mm,
        length_m=1.0e-3 * args.plasma_height_mm,
    )
    condition = ZhuOpenReactorCondition(
        condition_id="zhu-2026-npg80-tio2-chf3-sf6-o2-20min",
        geometry=geometry,
        neutral_control_volume_m3=geometry.volume_m3,
        pressure_Pa=3.99967104,
        gas_temperature_K=args.gas_temperature_K,
        feed_molecules_s=_feed_molecules_s(),
        absorbed_power_W=args.absorbed_power_W,
        source_frequency_hz=13.56e6,
        reduced_field_bounds_Td=(205.0, 600.0),
        ion_temperature_eV=args.ion_temperature_eV,
        ion_momentum_mean_free_path_m=1.0e-6 * args.ion_mfp_um,
        mean_positive_ion_wall_energy_eV=args.mean_wall_ion_energy_eV,
        neutral_reduced_diffusivity_m_inv_s=(
            args.neutral_reduced_diffusivity_m_inv_s),
        neutral_wall_probabilities={
            "F": args.f_wall_probability,
            "H": args.h_wall_probability,
            "O": args.o_wall_probability,
            "O(1d)": args.excited_o_wall_probability,
            "C": args.c_wall_probability,
            "CF3": args.cf3_wall_probability,
            "CF2": args.cf2_wall_probability,
            "CF": args.cf_wall_probability,
        },
        source=(
            "Zhu operator recipe manifest; sccm standard state explicitly "
            "taken as 273.15 K and 101325 Pa"
        ),
        absorbed_power_source=(
            "declared target-free sensitivity; generator setpoint is 150 W"
        ),
        machine_closure_source=(
            "Oxford-80 family geometry/voltage and target-pressure mobility "
            "sensitivity; not target-tool diagnostics"
        ),
    )
    model = ZhuOpenReactorModel(solver, parent, supplemental)
    initial_densities = None
    initial_exhaust = None
    continuation_field = None
    if args.initial_state_json is not None:
        initial_densities, initial_exhaust, continuation_field = (
            _continuation_state(
                args.initial_state_json, model=model, condition=condition))
    initial_field = (
        args.initial_field_Td
        if args.initial_field_Td is not None
        else (continuation_field if continuation_field is not None else 240.0)
    )
    solution = model.solve(
        condition,
        initial_densities_m3=initial_densities,
        initial_exhaust_loss_frequency_s_inv=initial_exhaust,
        initial_reduced_electric_field_Td=initial_field,
        maximum_evaluations=args.maximum_evaluations,
        residual_tolerance=args.residual_tolerance,
    )
    densities = dict(solution.densities_m3)
    positive_total = sum(
        densities[species.name] * species.charge_number
        for species in model.species if species.role == "positive_ion")
    negative_total = sum(
        densities[species.name] * -species.charge_number
        for species in model.species if species.role == "negative_ion")
    return {
        "schema": "petch.zhu-open-reactor-state.v1",
        "condition_id": solution.condition_id,
        "input": {
            "feed_sccm": FEED_SCCM,
            "pressure_Pa": condition.pressure_Pa,
            "gas_temperature_K": condition.gas_temperature_K,
            "forward_power_setpoint_W": 150.0,
            "absorbed_power_sensitivity_W": condition.absorbed_power_W,
            "electrode_diameter_mm": args.electrode_diameter_mm,
            "plasma_height_mm": args.plasma_height_mm,
            "ion_mfp_um": args.ion_mfp_um,
            "mean_all_wall_ion_energy_eV": args.mean_wall_ion_energy_eV,
            "kokkoris_eedf_shape": args.kokkoris_eedf_shape,
            "chf3_f_rate_branch": args.chf3_f_rate_branch,
            "continuation_state_json": (
                None if args.initial_state_json is None
                else str(args.initial_state_json)
            ),
            "feature_or_sem_target_used": False,
        },
        "state": {
            "reduced_electric_field_Td": solution.reduced_electric_field_Td,
            "mean_electron_energy_eV": solution.mean_electron_energy_eV,
            "electron_density_m3": densities["e"],
            "electronegativity": negative_total / densities["e"],
            "positive_charge_density_m3": positive_total,
            "exhaust_loss_frequency_s_inv": (
                solution.exhaust_loss_frequency_s_inv),
            "densities_m3": densities,
            "axial_positive_ion_flux_m2_s": dict(
                solution.axial_positive_ion_flux_m2_s),
            "neutral_thermal_flux_m2_s": dict(
                solution.neutral_thermal_flux_m2_s),
            "total_axial_positive_ion_flux_m2_s": (
                solution.total_axial_positive_ion_flux_m2_s),
            "electron_collision_basis_neutral_fraction": (
                solution.electron_collision_basis_neutral_fraction),
        },
        "power_density_W_m3": {
            "absorbed": solution.absorbed_power_density_W_m3,
            "parent_collision": (
                solution.parent_collision_power_density_W_m3),
            "supplemental_collision": (
                solution.supplemental_collision_power_density_W_m3),
            "charged_wall": solution.charged_wall_power_density_W_m3,
        },
        "numerics": {
            "solver_evaluations": solution.solver_evaluations,
            "maximum_normalized_residual": (
                solution.maximum_normalized_residual),
            "normalized_residuals": dict(solution.normalized_residuals),
        },
        "certification": {
            "conserved_open_reactor_equations_solved": True,
            "unique_machine_state": False,
            "wafer_flux_prediction": False,
            "feature_depth_prediction": False,
            "reason": (
                "absorbed power, active geometry, wall state, all-wall sheath "
                "energy, daughter-gas electron transport, and several "
                "cross-family reactions remain machine sensitivity inputs"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-workbook", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--initial-state-json", type=Path)
    parser.add_argument("--absorbed-power-W", type=float, default=90.0)
    parser.add_argument("--electrode-diameter-mm", type=float, default=170.0)
    parser.add_argument("--plasma-height-mm", type=float, default=30.0)
    parser.add_argument("--gas-temperature-K", type=float, default=350.0)
    parser.add_argument("--ion-temperature-eV", type=float, default=.03)
    parser.add_argument("--ion-mfp-um", type=float, default=100.0)
    parser.add_argument("--mean-wall-ion-energy-eV", type=float, default=250.0)
    parser.add_argument(
        "--neutral-reduced-diffusivity-m-inv-s", type=float, default=6.0e20)
    parser.add_argument("--f-wall-probability", type=float, default=.05)
    parser.add_argument("--h-wall-probability", type=float, default=.05)
    parser.add_argument("--o-wall-probability", type=float, default=.10)
    parser.add_argument("--excited-o-wall-probability", type=float, default=1.0)
    parser.add_argument("--c-wall-probability", type=float, default=1.0)
    parser.add_argument("--cf3-wall-probability", type=float, default=.05)
    parser.add_argument("--cf2-wall-probability", type=float, default=.10)
    parser.add_argument("--cf-wall-probability", type=float, default=.10)
    parser.add_argument(
        "--kokkoris-eedf-shape",
        choices=("druyvesteyn", "maxwellian"),
        default="druyvesteyn",
    )
    parser.add_argument(
        "--chf3-f-rate-branch",
        choices=("voloshin_350K", "lim_700K"),
        default="voloshin_350K",
    )
    parser.add_argument("--initial-field-Td", type=float)
    parser.add_argument("--maximum-evaluations", type=int, default=1000)
    parser.add_argument("--residual-tolerance", type=float, default=2.0e-6)
    args = parser.parse_args()
    result = run(args)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()
