#!/usr/bin/env python3
"""Run the preregistered Jeon calibration condition through the unified 3-D engine.

This is an explicitly non-predictive baseline until measured/tabulated IEDF, radical composition,
surface parameters, mask interaction, and initial geometry evidence replace the named closures.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from petch.experimental_boundary import (
    Jeon2022BoundaryClosure, build_jeon_2022_boundary_state,
)
from petch.experimental_data import (
    build_jeon_2022_dimensionless_targets,
    load_jeon_2022_electron_bias_controls,
    load_jeon_2022_plasma_controls,
    load_jeon_2022_trench_depths,
)
from petch.feature_step_3d import make_rectangular_trench_geometry_3d, solve_feature_3d
from petch.surface_kinetics import (
    EnergeticYield, ParameterEvidence, ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
)


ROOT = Path(__file__).parents[1]
DATA = ROOT / "data" / "experimental" / "jeon_2022"


def _evidence(source, evidence_type, *, supports=False, note=""):
    return ParameterEvidence(
        source, evidence_type, note=note,
        supports_prediction_within_declared_domain=supports)


def baseline_mechanism(complex_probability, deposition_probability):
    """Literature-order baseline with every unsupported number labeled as a closure."""
    evidence = {
        "site_density_m2": _evidence(
            "Kaler et al., J. Phys. D 50, 234001 (2017), CF2 uptake model S=1e15 cm^-2",
            "source_model_assumption"),
        "bulk_formula_density_m3": _evidence(
            "derived from fused-SiO2 density 2200 kg/m3 and molar mass 60.0843 g/mol",
            "material_constant_derived", supports=True),
        "polymer_monolayer_density_m2": _evidence(
            "Kaler et al. 2017 CF2 uptake site density 1e15 cm^-2",
            "source_model_assumption"),
        "complex_formation_probability": _evidence(
            "aggregate-radical calibration closure; Jeong species-resolved sticking unmeasured",
            "calibration_closure"),
        "polymer_deposition_probability_on_substrate": _evidence(
            "aggregate-radical calibration closure; Kaler CF2-on-SiO2 beam value was 0.19",
            "calibration_closure"),
        "polymer_deposition_probability_on_polymer": _evidence(
            "Kaler et al. 2017 fitted sFC/sox ratio 0.015/0.19",
            "source_ratio_applied_to_closure"),
        "oxygen_polymer_etch_probability": _evidence(
            "zero incident O channel in the reduced Ar/C4F8 baseline", "declared_absent_channel"),
        "bare_sio2_yield": _evidence(
            "order baseline constrained by Kaler 2017 65 eV sputter threshold",
            "nonpredictive_literature_order"),
        "complex_sio2_yield": _evidence(
            "order baseline: Takada et al. JAP 97, 013534 CF2/Ar+ yield near 900 eV",
            "nonpredictive_cross_condition"),
        "polymer_sputter_yield": _evidence(
            "order baseline; no matched Jeong polymer sputter table", "calibration_closure"),
    }
    return ReducedSiO2FluorocarbonMechanism(ReducedSiO2FluorocarbonParameters(
        site_density_m2=1e19,
        bulk_formula_density_m3=2.205e28,
        polymer_monolayer_density_m2=1e19,
        complex_formation_probability={"FC_total": complex_probability},
        polymer_deposition_probability_on_substrate={"FC_total": deposition_probability},
        polymer_deposition_probability_on_polymer={
            "FC_total": deposition_probability * 0.015 / 0.19},
        oxygen_species="O", oxygen_polymer_etch_probability=0.0,
        bare_sio2_yield=EnergeticYield(
            0.20, 65.0, 900.0, energy_exponent=0.5,
            angular_model="kress_1999", angular_parameter=9.3),
        complex_sio2_yield=EnergeticYield(
            1.60, 65.0, 900.0, energy_exponent=0.5,
            angular_model="chang_sawin_1997"),
        polymer_sputter_yield=EnergeticYield(
            0.50, 20.0, 900.0, energy_exponent=0.5,
            angular_model="kress_1999", angular_parameter=9.3),
        evidence=evidence))


def _floor_height(phi, dx):
    line = phi[phi.shape[0] // 2, phi.shape[1] // 2]
    transition = np.flatnonzero((line[:-1] >= 0.0) & (line[1:] < 0.0))
    if transition.size != 1:
        raise RuntimeError(f"center line has {transition.size} solid-to-gas transitions")
    index = int(transition[0])
    return dx * (index + line[index] / (line[index] - line[index + 1]))


def run_width(width_nm, args, plasma_control, electron_bias):
    width_um = width_nm * 1e-3
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=args.pitch_um, cell_length=args.cell_length_um,
        domain_height=args.domain_height_um, dx=args.dx_um, opening_width=width_um,
        mask_thickness=args.mask_thickness_um, substrate_top=args.substrate_top_um,
        etched_depth=args.initial_depth_um)
    source_z = (geometry.phi.shape[2] - 1) * geometry.dx
    closure = Jeon2022BoundaryClosure(
        ion_name="Ar+", ion_mass_amu=39.948,
        ion_normal_energy_eV=[electron_bias.self_bias_magnitude_v],
        ion_normal_energy_weight=[1.0], ion_tangential_temperature_eV=0.026,
        neutral_flux_fraction={"FC_total": 1.0}, neutral_mass_amu={"FC_total": 50.005},
        neutral_temperature_K=300.0,
        provenance={
            "model": "self_bias_monoenergy_and_aggregate_radical_baseline",
            "warning": "not a measured IEDF or species composition"},
        supports_prediction_within_declared_domain=False)
    boundary = build_jeon_2022_boundary_state(
        plasma_control, electron_bias, closure,
        reference_plane_m=source_z * geometry.mesh_length_unit_m,
        n_transverse_ion=3, n_transverse_neutral=3, n_normal_neutral=4)
    mechanism = baseline_mechanism(args.complex_probability, args.deposition_probability)
    initial_floor = _floor_height(geometry.phi, geometry.dx)
    started = perf_counter()
    result = solve_feature_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "FC_total": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=args.duration_s,
        n_steps=args.steps, source_bounds=(0.0, args.pitch_um, 0.0, args.cell_length_um),
        source_z=source_z, n_position=args.source_positions, seed=args.seed,
        cfl_number=0.3, reinitialize=True, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": args.form_factor_rays, "seed": args.seed + 1000,
            "periodic_lateral": True,
            "domain_size": np.asarray(geometry.phi.shape) * geometry.dx,
            "nonetchable_reaction_probability_by_material": {
                2: {"FC_total": args.mask_reaction_probability}},
        },
        ballistic_transport=args.ballistic_transport,
        ballistic_face_quadrature_points=args.ballistic_face_quadrature_points)
    wall = perf_counter() - started
    floor = _floor_height(result.geometry.phi, geometry.dx)
    initial_depth_nm = (args.substrate_top_um - initial_floor) * 1000.0
    final_depth_nm = (args.substrate_top_um - floor) * 1000.0
    return {
        "width_nm": width_nm,
        "initialized_depth_nm": initial_depth_nm,
        "total_depth_nm": final_depth_nm,
        "etched_increment_nm": final_depth_nm - initial_depth_nm,
        "wall_time_s": wall,
        "steps": args.steps,
        "within_declared_scope": result.validity.within_declared_scope,
        "parameter_evidence_supports_prediction": (
            result.validity.parameter_evidence_supports_prediction),
        "nonpredictive_parameters": list(result.validity.nonpredictive_parameters),
        "known_limitations": list(result.validity.known_limitations),
        "maximum_neutral_balance_error": max(
            step.diagnostics["neutral_radiosity"]["FC_total"]["relative_balance_error"]
            for step in result.steps),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--widths-nm", type=float, nargs="+", default=[200.0])
    parser.add_argument("--duration-s", type=float, default=1000.0)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--pitch-um", type=float, default=0.50)
    parser.add_argument("--cell-length-um", type=float, default=0.10)
    parser.add_argument("--mask-thickness-um", type=float, default=0.70)
    parser.add_argument("--substrate-top-um", type=float, default=1.40)
    parser.add_argument("--domain-height-um", type=float, default=2.35)
    parser.add_argument("--initial-depth-um", type=float, default=0.03)
    parser.add_argument("--complex-probability", type=float, default=1e-3)
    parser.add_argument("--deposition-probability", type=float, default=5e-4)
    parser.add_argument("--mask-reaction-probability", type=float, default=1e-3)
    parser.add_argument("--source-positions", type=int, default=32)
    parser.add_argument("--form-factor-rays", type=int, default=32)
    parser.add_argument(
        "--ballistic-transport", choices=("forward", "face_gather"),
        default="face_gather")
    parser.add_argument(
        "--ballistic-face-quadrature-points", type=int, choices=(1, 3), default=3)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    plasma = load_jeon_2022_plasma_controls(DATA / "digitized_plasma_controls.csv")
    electron = load_jeon_2022_electron_bias_controls(
        DATA / "digitized_electron_bias_controls.csv")
    depths = load_jeon_2022_trench_depths(DATA / "digitized_trench_depths.csv")
    plasma_control = next(item for item in plasma
                          if item.condition_family == "gas_fraction_cw"
                          and item.c4f8_fraction == 0.2)
    electron_bias = next(item for item in electron
                         if item.condition_family == "gas_fraction_cw"
                         and item.c4f8_fraction == 0.2)
    runs = [run_width(width, args, plasma_control, electron_bias)
            for width in args.widths_nm]
    calibration_targets = {
        item.trench_width_nm: item
        for item in build_jeon_2022_dimensionless_targets(depths)
        if item.split == "calibration"
        and item.observable == "width_shape_depth_over_200nm"
    }
    width_shape = None
    by_width = {item["width_nm"]: item for item in runs}
    if 200.0 in by_width and by_width[200.0]["etched_increment_nm"] > 0.0:
        reference = by_width[200.0]["etched_increment_nm"]
        scored = []
        for width in sorted(by_width):
            if width not in calibration_targets:
                continue
            prediction = float(by_width[width]["etched_increment_nm"] / reference)
            target = calibration_targets[width]
            scored.append({
                "width_nm": width,
                "prediction": prediction,
                "target": float(target.value),
                "digitization_interval": [
                    float(target.digitization_lower), float(target.digitization_upper)],
                "within_digitization_interval": bool(
                    target.digitization_lower <= prediction <= target.digitization_upper),
                "log_residual": float(np.log(prediction / target.value)),
            })
        width_shape = {
            "normalization": "simulated_increment_over_same_run_200nm_increment",
            "absolute_duration_not_scored": True,
            "points": scored,
            "log_rmse": float(np.sqrt(np.mean([
                item["log_residual"] ** 2 for item in scored]))) if scored else None,
        }
    output = {
        "campaign": "jeon_2022_unified_nonpredictive_baseline",
        "closures": {
            "iedf": "self_bias_monoenergy_nonpredictive",
            "neutral_composition": "aggregate_FC_total_nonpredictive",
            "pitch_um": args.pitch_um,
            "periodic_cell_length_um": args.cell_length_um,
            "mask_thickness_um": args.mask_thickness_um,
            "mask_reaction_probability": args.mask_reaction_probability,
            "initial_bare_sidewall_depth_um": args.initial_depth_um,
        },
        "numerics": {
            "dx_um": args.dx_um,
            "duration_s": args.duration_s,
            "steps": args.steps,
            "source_positions": args.source_positions,
            "form_factor_rays": args.form_factor_rays,
            "ballistic_transport": args.ballistic_transport,
            "ballistic_face_quadrature_points": args.ballistic_face_quadrature_points,
            "seed": args.seed,
        },
        "runs": runs,
        "calibration_width_shape": width_shape,
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
