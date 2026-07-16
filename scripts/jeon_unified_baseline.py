#!/usr/bin/env python3
"""Run the preregistered Jeon calibration condition through the unified 3-D engine.

This is an explicitly non-predictive baseline until measured/tabulated IEDF, radical composition,
surface parameters, mask interaction, and initial geometry evidence replace the named closures.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from petch.boundary_state import (
    IonEnergyTransverseMaxwellianDensity, PlasmaBoundaryState, SpeciesBoundaryState,
    qmc_boundary_proposal,
)
from petch.experimental_boundary import (
    Jeon2022BoundaryClosure, build_jeon_2022_boundary_state,
)
from petch.experimental_data import (
    build_jeon_2022_dimensionless_targets,
    jeon_2022_condition_wall_duration_s,
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


def _write_json_atomic(path: Path, payload):
    """Persist resumable campaign evidence without exposing a half-written JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _canonical_input_deck(args):
    """Return the complete replay deck and a content hash independent of output location."""
    deck = {}
    for key, value in sorted(vars(args).items()):
        if key in {"output", "quiet"}:
            continue
        if isinstance(value, Path):
            value = str(value)
        deck[key] = value
    encoded = json.dumps(
        deck, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return deck, hashlib.sha256(encoded).hexdigest()


def _evidence(source, evidence_type, *, supports=False, note=""):
    return ParameterEvidence(
        source, evidence_type, note=note,
        supports_prediction_within_declared_domain=supports)


def _wall_time_boundary(boundary, condition_family, pulse_off_ms):
    """Convert Jeon's on-state Bohm closure to the paper's wall-time pulse measure."""
    if not str(condition_family).startswith("pulse_off_"):
        duty = 1.0
    else:
        # Both pulse series declare a fixed 1 ms pulse-on time.  Their plotted neutral/ion ratio
        # already includes the ion-off interval, while the Bohm flux reconstructed from electron
        # density is the pulse-on value.  Scaling the complete ratio-bearing boundary by duty gives
        # ion_on*duty and neutral=ratio*ion_on*duty without changing the reported ratio.
        duty = 1.0 / (1.0 + float(pulse_off_ms))
    if duty == 1.0:
        return boundary, duty
    species = tuple(SpeciesBoundaryState(
        name=item.name,
        charge_number=item.charge_number,
        mass_amu=item.mass_amu,
        flux_m2_s=item.flux_m2_s * duty,
        velocity_sqrt_eV=item.velocity_sqrt_eV,
        weight=item.weight,
        phase_rad=item.phase_rad,
        position_m=item.position_m,
        density_model=item.density_model,
        provenance={**dict(item.provenance), "wall_time_duty_factor": duty},
    ) for item in boundary.species)
    return PlasmaBoundaryState(
        species,
        boundary.reference_plane_m,
        provenance={
            **dict(boundary.provenance),
            "pulse_on_ms": 1.0,
            "pulse_off_ms": float(pulse_off_ms),
            "wall_time_duty_factor": duty,
        },
    ), duty


def _qmc_ion_boundary(
        boundary, *, ion_name, normal_energy_eV, energy_halfwidth_eV,
        tangential_temperature_eV, log2_samples, seed):
    """Replace a tensor ion rule by an equal-weight scrambled-Sobol flux rule.

    Hard mask visibility is discontinuous in launch angle, for which low-order tensor
    Gauss-Hermite rules retain a spuriously large exactly-vertical atom.  This proposal samples the
    same continuous transverse Maxwellian and a declared one-bin normal-energy closure.  Independent
    scrambles and sample-level refinement quantify the remaining integration error.
    """
    ion = boundary.get(ion_name)
    lower = max(float(normal_energy_eV) - float(energy_halfwidth_eV), 0.0)
    upper = float(normal_energy_eV) + float(energy_halfwidth_eV)
    if not upper > lower:
        raise ValueError("QMC ion-energy interval must have positive width")
    density = IonEnergyTransverseMaxwellianDensity(
        np.asarray([lower, upper]), np.asarray([1.0]),
        float(tangential_temperature_eV))
    template = SpeciesBoundaryState(
        name=ion.name, charge_number=ion.charge_number, mass_amu=ion.mass_amu,
        flux_m2_s=ion.flux_m2_s,
        velocity_sqrt_eV=[[0.0, 0.0, np.sqrt(0.5 * (lower + upper))]],
        weight=[1.0], density_model=density,
        provenance=dict(ion.provenance))
    sampled = qmc_boundary_proposal(
        template, int(log2_samples), seed=int(seed), name=ion.name)
    qmc_ion = SpeciesBoundaryState(
        name=sampled.name, charge_number=sampled.charge_number,
        mass_amu=sampled.mass_amu, flux_m2_s=ion.flux_m2_s,
        velocity_sqrt_eV=sampled.velocity_sqrt_eV, weight=sampled.weight,
        density_model=density,
        provenance={
            **dict(ion.provenance),
            "numerical_rule": "scrambled_sobol_continuous_iedf_iadf",
            "normal_energy_interval_eV": [lower, upper],
            "log2_samples": int(log2_samples), "seed": int(seed),
        })
    species = tuple(qmc_ion if item.name == ion_name else item for item in boundary.species)
    return PlasmaBoundaryState(
        species, boundary.reference_plane_m,
        provenance={
            **dict(boundary.provenance),
            "ion_numerical_rule": "scrambled_sobol_continuous_iedf_iadf",
            "ion_qmc_log2_samples": int(log2_samples),
            "ion_qmc_seed": int(seed),
        })


def baseline_mechanism(
        complex_probability, deposition_probability, complex_removal_reaction_order=1, *,
        bare_reference_yield=0.20, complex_reference_yield=1.60,
        polymer_reference_yield=0.50):
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
        "complex_removal_reaction_order": _evidence(
            "W. Guo, MIT PhD thesis (2009), Sec. 4.3; nearest-neighbour "
            "COF2 event proportional to C-O times (C-F)^2",
            "source_model_topology",
            note="bounded campaign choice: linear order 1 or neighbour-bond order 2"),
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
            bare_reference_yield, 65.0, 900.0, energy_exponent=0.5,
            angular_model="kress_1999", angular_parameter=9.3),
        complex_sio2_yield=EnergeticYield(
            complex_reference_yield, 65.0, 900.0, energy_exponent=0.5,
            angular_model="chang_sawin_1997"),
        polymer_sputter_yield=EnergeticYield(
            polymer_reference_yield, 20.0, 900.0, energy_exponent=0.5,
            angular_model="kress_1999", angular_parameter=9.3),
        complex_removal_reaction_order=complex_removal_reaction_order,
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
    if args.ion_iad_component_sigma_deg is None:
        ion_tangential_temperature_eV = 0.026
        ion_iad_component_sigma_deg = float(np.rad2deg(np.sqrt(
            ion_tangential_temperature_eV
            / (2.0 * electron_bias.self_bias_magnitude_v))))
        ion_iad_source = "room-temperature tangential baseline"
    else:
        ion_iad_component_sigma_deg = float(args.ion_iad_component_sigma_deg)
        if (not np.isfinite(ion_iad_component_sigma_deg)
                or ion_iad_component_sigma_deg <= 0.0):
            raise ValueError("ion IAD component sigma must be positive and finite")
        sigma_rad = float(np.deg2rad(ion_iad_component_sigma_deg))
        ion_tangential_temperature_eV = float(
            2.0 * electron_bias.self_bias_magnitude_v * sigma_rad * sigma_rad)
        ion_iad_source = "calibration closure; Jeon IADF unmeasured"
    closure = Jeon2022BoundaryClosure(
        ion_name="Ar+", ion_mass_amu=39.948,
        ion_normal_energy_eV=[electron_bias.self_bias_magnitude_v],
        ion_normal_energy_weight=[1.0],
        ion_tangential_temperature_eV=ion_tangential_temperature_eV,
        neutral_flux_fraction={"FC_total": 1.0}, neutral_mass_amu={"FC_total": 50.005},
        neutral_temperature_K=300.0,
        provenance={
            "model": "self_bias_monoenergy_and_aggregate_radical_baseline",
            "warning": "not a measured IEDF, IADF, or species composition",
            "ion_iad_component_sigma_deg": ion_iad_component_sigma_deg,
            "ion_iad_source": ion_iad_source},
        supports_prediction_within_declared_domain=False)
    boundary = build_jeon_2022_boundary_state(
        plasma_control, electron_bias, closure,
        reference_plane_m=source_z * geometry.mesh_length_unit_m,
        n_transverse_ion=args.ion_transverse_quadrature,
        n_transverse_neutral=3, n_normal_neutral=4)
    if args.ion_quadrature == "qmc":
        boundary = _qmc_ion_boundary(
            boundary, ion_name="Ar+",
            normal_energy_eV=electron_bias.self_bias_magnitude_v,
            energy_halfwidth_eV=args.ion_energy_halfwidth_ev,
            tangential_temperature_eV=ion_tangential_temperature_eV,
            log2_samples=args.ion_qmc_log2,
            seed=(args.seed + 2000 if args.ion_qmc_seed is None else args.ion_qmc_seed))
    boundary, wall_time_duty = _wall_time_boundary(
        boundary, args.condition_family, args.pulse_off_ms)
    condition_duration_s = jeon_2022_condition_wall_duration_s(
        args.duration_s, wall_time_duty, args.pulse_exposure_basis)
    mechanism = baseline_mechanism(
        args.complex_probability, args.deposition_probability,
        args.complex_removal_reaction_order,
        bare_reference_yield=args.bare_reference_yield,
        complex_reference_yield=args.complex_reference_yield,
        polymer_reference_yield=args.polymer_reference_yield)
    initial_floor = _floor_height(geometry.phi, geometry.dx)
    started = perf_counter()
    result = solve_feature_3d(
        geometry, boundary,
        {"Ar+": "energetic_bombardment", "FC_total": "neutral_reactant"},
        mechanism, etchable_material_ids=(1,), duration_s=condition_duration_s,
        n_steps=args.steps, source_bounds=(0.0, args.pitch_um, 0.0, args.cell_length_um),
        source_z=source_z, n_position=args.source_positions, seed=args.seed,
        cfl_number=0.3, reinitialize=args.reinitialize, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": args.form_factor_rays, "seed": args.seed + 1000,
            "periodic_lateral": True,
            "domain_size": (np.asarray(geometry.phi.shape) - 1) * geometry.dx,
            "nonetchable_reaction_probability_by_material": {
                2: {"FC_total": args.mask_reaction_probability}},
        },
        ballistic_transport=args.ballistic_transport,
        ballistic_face_quadrature_points=args.ballistic_face_quadrature_points,
        reinitialization_method=args.reinitialization_method)
    wall = perf_counter() - started
    floor = _floor_height(result.geometry.phi, geometry.dx)
    initial_depth_nm = (args.substrate_top_um - initial_floor) * 1000.0
    final_depth_nm = (args.substrate_top_um - floor) * 1000.0
    final_centroid = result.steps[-1].next_active_face_centroid
    final_area = result.steps[-1].next_active_face_area
    final_state = result.surface_state
    floor_band = final_centroid[:, 2] <= float(np.min(final_centroid[:, 2]) + geometry.dx)

    def weighted_mean(value, selected=None):
        value = np.asarray(value, dtype=float)
        use = np.ones(value.shape, dtype=bool) if selected is None else np.asarray(selected, dtype=bool)
        return float(np.dot(value[use], final_area[use]) / final_area[use].sum())

    state_summary = {
        "area_weighted_complex_fraction": weighted_mean(final_state.complex_fraction),
        "area_weighted_polymer_units_m2": weighted_mean(final_state.polymer_units_m2),
        "floor_complex_fraction": weighted_mean(final_state.complex_fraction, floor_band),
        "floor_polymer_units_m2": weighted_mean(final_state.polymer_units_m2, floor_band),
        "maximum_polymer_units_m2": float(np.max(final_state.polymer_units_m2)),
    }
    history = None
    if args.history_every:
        history = []
        for step_index, step_result in enumerate(result.steps, start=1):
            if step_index % args.history_every == 0 or step_index == len(result.steps):
                step_floor = _floor_height(step_result.geometry.phi, geometry.dx)
                active_velocity = step_result.face_velocity_mesh_units_s[
                    step_result.active_face_index]
                active_centroid = step_result.active_face_centroid
                active_area = step_result.active_face_area
                active_floor = active_centroid[:, 2] <= (
                    float(np.min(active_centroid[:, 2])) + geometry.dx)
                history.append({
                    "step": step_index,
                    "time_s": condition_duration_s * step_index / args.steps,
                    "etched_increment_nm": (
                        initial_floor - step_floor) * 1000.0,
                    "ion_hit_probability": step_result.transport.hit_probability["Ar+"],
                    "floor_area_weighted_velocity_um_s": float(np.dot(
                        active_velocity[active_floor], active_area[active_floor])
                        / active_area[active_floor].sum()),
                    "maximum_face_velocity_um_s": float(np.max(active_velocity)),
                    "centerline_extended_velocity_um_s": step_result.diagnostics.get(
                        "centerline_extended_velocity_mesh_units_s"),
                    "centerline_interface_fraction": step_result.diagnostics.get(
                        "centerline_interface_fraction"),
                    "centerline_advected_interface_fraction": step_result.diagnostics.get(
                        "centerline_advected_interface_fraction"),
                    "centerline_reinitialized_interface_fraction": step_result.diagnostics.get(
                        "centerline_reinitialized_interface_fraction"),
                })
    return {
        "width_nm": width_nm,
        "initialized_depth_nm": initial_depth_nm,
        "total_depth_nm": final_depth_nm,
        "etched_increment_nm": final_depth_nm - initial_depth_nm,
        "wall_time_s": wall,
        "steps": args.steps,
        "condition_wall_duration_s": condition_duration_s,
        "reference_duration_s": args.duration_s,
        "pulse_exposure_basis": args.pulse_exposure_basis,
        "within_declared_scope": result.validity.within_declared_scope,
        "parameter_evidence_supports_prediction": (
            result.validity.parameter_evidence_supports_prediction),
        "nonpredictive_parameters": list(result.validity.nonpredictive_parameters),
        "known_limitations": list(result.validity.known_limitations),
        "wall_time_duty_factor": wall_time_duty,
        "maximum_neutral_balance_error": max(
            step.diagnostics["neutral_radiosity"]["FC_total"]["relative_balance_error"]
            for step in result.steps),
        "surface_state": state_summary,
        "depth_history": history,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--widths-nm", type=float, nargs="+", default=[200.0])
    parser.add_argument(
        "--condition-family",
        choices=("gas_fraction_cw", "pulse_off_20pct", "pulse_off_80pct"),
        default="gas_fraction_cw")
    parser.add_argument("--c4f8-fraction", type=float, default=0.2)
    parser.add_argument("--pulse-off-ms", type=float, default=0.0)
    parser.add_argument("--duration-s", type=float, default=1000.0)
    parser.add_argument(
        "--pulse-exposure-basis",
        choices=("unspecified", "wall_time", "rf_on_time"), default="unspecified",
        help=("meaning of --duration-s for a pulsed condition; Jeon did not report whether "
              "wall time or cumulative RF-on time was held fixed"))
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--pitch-um", type=float, default=0.50)
    parser.add_argument("--cell-length-um", type=float, default=0.10)
    parser.add_argument("--mask-thickness-um", type=float, default=0.70)
    parser.add_argument("--substrate-top-um", type=float, default=1.40)
    parser.add_argument("--domain-height-um", type=float, default=2.35)
    parser.add_argument("--initial-depth-um", type=float, default=0.0)
    parser.add_argument("--complex-probability", type=float, default=1e-3)
    parser.add_argument("--deposition-probability", type=float, default=5e-4)
    parser.add_argument(
        "--complex-removal-reaction-order", type=int, choices=(1, 2), default=1,
        help=("complex-coverage activation order: 1 is the legacy reduced law; "
              "2 is the nearest-neighbour mixing-layer law"))
    parser.add_argument(
        "--bare-reference-yield", type=float, default=0.20,
        help="bare-SiO2 energetic yield at 900 eV (development closure; bounded to [0, 5])")
    parser.add_argument(
        "--complex-reference-yield", type=float, default=1.60,
        help=("FC-complex-assisted SiO2 energetic yield at 900 eV "
              "(development closure; bounded to [0, 5])"))
    parser.add_argument(
        "--polymer-reference-yield", type=float, default=0.50,
        help="polymer sputter yield at 900 eV (development closure; bounded to [0, 5])")
    parser.add_argument("--mask-reaction-probability", type=float, default=1e-3)
    parser.add_argument(
        "--ion-iad-component-sigma-deg", type=float,
        help=("standard deviation of either transverse ion-angle component; "
              "default preserves the 0.026 eV tangential baseline"))
    parser.add_argument(
        "--ion-transverse-quadrature", type=int, default=3,
        help=("Gauss-Hermite order per transverse ion-velocity component; "
              "must be convergence-tested for broad IADFs and narrow openings"))
    parser.add_argument(
        "--ion-quadrature", choices=("gauss_hermite", "qmc"),
        default="gauss_hermite")
    parser.add_argument("--ion-qmc-log2", type=int, default=12)
    parser.add_argument("--ion-qmc-seed", type=int)
    parser.add_argument(
        "--ion-energy-halfwidth-ev", type=float, default=20.0,
        help="one-bin IEDF half-width used only by QMC; default is Jeon's bias digitization bound")
    parser.add_argument("--source-positions", type=int, default=32)
    parser.add_argument("--form-factor-rays", type=int, default=32)
    parser.add_argument(
        "--ballistic-transport", choices=("forward", "face_gather"),
        default="face_gather")
    parser.add_argument(
        "--ballistic-face-quadrature-points", type=int, choices=(1, 3), default=3)
    parser.add_argument(
        "--reinitialize", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--reinitialization-method", choices=("skfmm", "fsm", "cr2"), default="cr2")
    parser.add_argument("--history-every", type=int, default=0)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--output", type=Path,
        help="optional JSON artifact path")
    parser.add_argument(
        "--quiet", action="store_true",
        help="suppress the final stdout JSON (progress lines are still emitted)")
    parser.add_argument(
        "--allow-unscored-widths", action="store_true",
        help=("permit diagnostic simulated openings that do not match a registered Jeon width; "
              "such runs can never populate calibration_width_shape"))
    args = parser.parse_args()
    if (args.ion_transverse_quadrature <= 0 or args.ion_qmc_log2 < 0
            or not np.isfinite(args.ion_energy_halfwidth_ev)
            or args.ion_energy_halfwidth_ev <= 0.0):
        raise ValueError("invalid ion quadrature controls")
    reference_yields = (
        args.bare_reference_yield, args.complex_reference_yield,
        args.polymer_reference_yield)
    if any(not np.isfinite(value) or value < 0.0 or value > 5.0
           for value in reference_yields):
        raise ValueError("development reference yields must lie in [0, 5]")
    plasma = load_jeon_2022_plasma_controls(DATA / "digitized_plasma_controls.csv")
    electron = load_jeon_2022_electron_bias_controls(
        DATA / "digitized_electron_bias_controls.csv")
    depths = load_jeon_2022_trench_depths(DATA / "digitized_trench_depths.csv")
    matching_plasma = [
        item for item in plasma
        if item.condition_family == args.condition_family
        and np.isclose(item.c4f8_fraction, args.c4f8_fraction, rtol=0.0, atol=1e-12)
        and np.isclose(item.pulse_off_ms, args.pulse_off_ms, rtol=0.0, atol=1e-12)]
    matching_electron = [
        item for item in electron
        if item.condition_family == args.condition_family
        and np.isclose(item.c4f8_fraction, args.c4f8_fraction, rtol=0.0, atol=1e-12)
        and np.isclose(item.pulse_off_ms, args.pulse_off_ms, rtol=0.0, atol=1e-12)]
    if len(matching_plasma) != 1 or len(matching_electron) != 1:
        raise ValueError(
            "requested Jeon condition must match exactly one plasma and electron/bias record")
    width_targets = {
        item.trench_width_nm: item
        for item in build_jeon_2022_dimensionless_targets(depths)
        if item.observable == "width_shape_depth_over_200nm"
        and item.condition_family == args.condition_family
        and np.isclose(item.c4f8_fraction, args.c4f8_fraction, rtol=0.0, atol=1e-12)
        and np.isclose(item.pulse_off_ms, args.pulse_off_ms, rtol=0.0, atol=1e-12)
    }
    unscored_widths = sorted(
        float(width) for width in args.widths_nm if float(width) not in width_targets)
    if unscored_widths and not args.allow_unscored_widths:
        raise ValueError(
            "requested simulated widths do not match the registered Jeon geometry: "
            f"{unscored_widths}; use --allow-unscored-widths only for labeled diagnostics")
    input_deck, input_deck_sha256 = _canonical_input_deck(args)
    plasma_control = matching_plasma[0]
    electron_bias = matching_electron[0]
    runs = []
    for run_index, width in enumerate(args.widths_nm, start=1):
        run = run_width(width, args, plasma_control, electron_bias)
        runs.append(run)
        print(
            f"completed width {width:g} nm ({run_index}/{len(args.widths_nm)}): "
            f"depth={run['etched_increment_nm']:.6g} nm, "
            f"wall={run['wall_time_s']:.2f} s",
            flush=True)
        if args.output is not None:
            _write_json_atomic(args.output, {
                "schema_version": "jeon_2022_unified_baseline_v2",
                "campaign": "jeon_2022_unified_nonpredictive_baseline",
                "status": "incomplete",
                "input_deck": input_deck,
                "input_deck_sha256": input_deck_sha256,
                "completed_widths": run_index,
                "requested_widths": len(args.widths_nm),
                "runs": runs,
            })
    width_shape = None
    by_width = {item["width_nm"]: item for item in runs}
    if 200.0 in by_width and by_width[200.0]["etched_increment_nm"] > 0.0:
        reference = by_width[200.0]["etched_increment_nm"]
        scored = []
        for width in sorted(by_width):
            if width not in width_targets:
                continue
            prediction = float(by_width[width]["etched_increment_nm"] / reference)
            target = width_targets[width]
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
            "split": sorted({item.split for item in width_targets.values()}),
            "requested_widths_nm": sorted(float(width) for width in by_width),
            "registered_widths_not_run_nm": sorted(
                float(width) for width in set(width_targets) - set(by_width)),
            "unscored_simulated_widths_nm": unscored_widths,
            "complete_registered_width_coverage": bool(
                not unscored_widths and set(by_width) == set(width_targets)),
            "points": scored,
            "log_rmse": float(np.sqrt(np.mean([
                item["log_residual"] ** 2 for item in scored]))) if scored else None,
        }
    output = {
        "schema_version": "jeon_2022_unified_baseline_v2",
        "campaign": "jeon_2022_unified_nonpredictive_baseline",
        "input_deck": input_deck,
        "input_deck_sha256": input_deck_sha256,
        "condition": {
            "condition_family": args.condition_family,
            "c4f8_fraction": args.c4f8_fraction,
            "pulse_off_ms": args.pulse_off_ms,
            "pulse_on_ms": (1.0 if args.condition_family.startswith("pulse_off_") else None),
            "pulse_exposure_basis": args.pulse_exposure_basis,
            "source_exposure_protocol": "not_reported",
        },
        "closures": {
            "iedf": "self_bias_monoenergy_nonpredictive",
            "neutral_composition": "aggregate_FC_total_nonpredictive",
            "surface_mechanism": "reduced_si_o2_fluorocarbon",
            "surface_reaction_probabilities": {
                "complex_formation_on_sio2": {
                    "value": args.complex_probability,
                    "bounds": [0.0, 1.0],
                    "source": "Jeon 20% C4F8 CW calibration-only development closure",
                },
                "polymer_deposition_on_sio2": {
                    "value": args.deposition_probability,
                    "bounds": [0.0, 1.0],
                    "source": "Jeon 20% C4F8 CW calibration-only development closure",
                },
                "polymer_deposition_on_polymer": {
                    "value": args.deposition_probability * 0.015 / 0.19,
                    "bounds": [0.0, 1.0],
                    "source": (
                        "substrate probability scaled by the declared 0.015/0.19 "
                        "polymer/substrate sticking-ratio closure"),
                },
            },
            "complex_removal_reaction_order": args.complex_removal_reaction_order,
            "complex_removal_reaction_order_bounds": [1, 2],
            "complex_removal_reaction_order_source": (
                "W. Guo, MIT PhD thesis (2009), Sec. 4.3"),
            "energetic_reference_yields_at_900eV": {
                "bare_sio2": {
                    "value": args.bare_reference_yield,
                    "bounds": [0.0, 5.0],
                    "source": (
                        "development closure; 65 eV threshold and angular order from "
                        "Kaler et al., J. Phys. D 50, 234001 (2017)"),
                },
                "fc_complex_sio2": {
                    "value": args.complex_reference_yield,
                    "bounds": [0.0, 5.0],
                    "source": (
                        "development closure; energy scale and order constrained by "
                        "Takada et al., J. Appl. Phys. 97, 013534 (2005)"),
                },
                "polymer": {
                    "value": args.polymer_reference_yield,
                    "bounds": [0.0, 5.0],
                    "source": (
                        "development closure; no Jeong-matched polymer sputter table"),
                },
                "bounds_semantics": (
                    "declared conservative development-screen bounds, not source uncertainty"),
            },
            "pitch_um": args.pitch_um,
            "periodic_cell_length_um": args.cell_length_um,
            "mask_thickness_um": args.mask_thickness_um,
            "mask_reaction_probability": args.mask_reaction_probability,
            "ion_iad_component_sigma_deg": (
                float(np.rad2deg(np.sqrt(
                    0.026 / (2.0 * electron_bias.self_bias_magnitude_v))))
                if args.ion_iad_component_sigma_deg is None
                else args.ion_iad_component_sigma_deg),
            "ion_transverse_quadrature": args.ion_transverse_quadrature,
            "ion_quadrature": args.ion_quadrature,
            "ion_qmc_log2": args.ion_qmc_log2,
            "ion_qmc_seed": (
                args.seed + 2000 if args.ion_qmc_seed is None else args.ion_qmc_seed),
            "ion_energy_halfwidth_eV": args.ion_energy_halfwidth_ev,
            "initial_bare_sidewall_depth_um": args.initial_depth_um,
        },
        "numerics": {
            "dx_um": args.dx_um,
            "duration_s": args.duration_s,
            "duration_s_semantics": (
                "cumulative_rf_on_time" if args.pulse_exposure_basis == "rf_on_time"
                else "wall_time" if args.pulse_exposure_basis == "wall_time"
                else "continuous_wave_identity_or_unspecified"),
            "steps": args.steps,
            "source_positions": args.source_positions,
            "form_factor_rays": args.form_factor_rays,
            "ballistic_transport": args.ballistic_transport,
            "ballistic_face_quadrature_points": args.ballistic_face_quadrature_points,
            "reinitialize": args.reinitialize,
            "reinitialization_method": args.reinitialization_method,
            "seed": args.seed,
        },
        "runs": runs,
        "width_shape_score": width_shape,
        "calibration_width_shape": (
            width_shape if width_shape is not None
            and width_shape["split"] == ["calibration"]
            and width_shape["complete_registered_width_coverage"] else None),
    }
    if args.output is not None:
        _write_json_atomic(args.output, output)
    if not args.quiet:
        print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
