#!/usr/bin/env python3
"""Bounded moving-profile development gate for the de Boer SF6/O2 operator.

This is not an experimental validation.  Every de Boer Figure-9 marker has already been exposed,
and the paper does not declare enough boundary inputs to run an absolute predictive case.  The
purpose of this script is narrower and mechanical: run the source-correct common Belen operator,
neutral radiosity, and certified grazing-ion response on an actually moving trench, then compare
that trajectory with the frozen-geometry rate-curve counterfactual.

Two nested sampling levels are run with identical physics.  The output records conservative
response ledgers, surface fixed-point residuals, profile depth at every step, and the initial/final
level sets needed for visual inspection.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from deboer_feature3d import (  # noqa: E402
    build_common_belen_si_mechanism,
    ion_species,
    thermal_neutral,
)
from petch.boundary_state import PlasmaBoundaryState  # noqa: E402
from petch.charged_surface_response_3d import (  # noqa: E402
    GrazingSpecularIonReflection3D,
)
from petch.feature_step_3d import make_rectangular_trench_geometry_3d  # noqa: E402
from petch.physical_api import PhysicalProcess  # noqa: E402


DEFAULT_STATIC_CURVE = (
    ROOT / "results" / "deboer_2002_belen_reflection_pilot" / "audit.json")


def _sha256(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _git_revision():
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _atomic_json(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _floor_depth_um(geometry, substrate_top_um):
    center = (geometry.phi.shape[0] // 2, geometry.phi.shape[1] // 2)
    line = np.asarray(geometry.phi[center], dtype=float)
    crossing = np.flatnonzero((line[:-1] >= 0.0) & (line[1:] < 0.0))
    if crossing.size != 1:
        raise RuntimeError(
            "moving profile does not have one resolved silicon-floor centerline crossing")
    lower = int(crossing[0])
    fraction = float(line[lower] / (line[lower] - line[lower + 1]))
    floor_z = (lower + fraction) * geometry.dx
    return float(substrate_top_um - floor_z)


def _static_depth_um(curve, *, opening_um, initial_depth_um, duration_s):
    aspect_ratio = np.asarray(curve["aspect_ratio"], dtype=float)
    rate = np.asarray(curve["raw_rate_m_s"], dtype=float)
    if (aspect_ratio[0] != 0.0 or np.any(np.diff(aspect_ratio) <= 0.0)
            or np.any(rate <= 0.0) or initial_depth_um / opening_um >= aspect_ratio[-1]):
        raise ValueError("static rate curve does not support the moving-profile initial state")
    log_rate = PchipInterpolator(aspect_ratio, np.log(rate), extrapolate=False)

    def rhs(_time, depth):
        ar = float(depth[0]) / float(opening_um)
        if ar >= aspect_ratio[-1]:
            raise RuntimeError("static counterfactual left the audited aspect-ratio support")
        return [1e6 * float(np.exp(log_rate(ar)))]

    solution = solve_ivp(
        rhs, (0.0, float(duration_s)), [float(initial_depth_um)],
        rtol=1e-9, atol=1e-12, max_step=max(float(duration_s) / 100.0, 1e-6))
    if not solution.success:
        raise RuntimeError("static rate-curve integration failed: " + solution.message)
    return float(solution.y[0, -1])


def _curve_from_audit(path):
    payload = json.loads(Path(path).read_text())
    selected = float(payload["selected_s_F"])
    candidate = next(
        item for item in payload["candidate_results"]
        if float(item["s_F"]) == selected)
    curve = candidate["curve"]
    if not curve.get("ion_reflection", False):
        raise ValueError("static counterfactual must use the reflected-ion common operator")
    return selected, curve


def _case(*, opening_um, dx_um, initial_aspect_ratio, velocity_log2,
          n_position, radiosity_rays, seed, duration_s, n_steps, s_f):
    initial_depth_um = float(opening_um) * float(initial_aspect_ratio)
    substrate_top_um = initial_depth_um + max(4.0 * dx_um, 0.05)
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=2.0 * opening_um,
        cell_length=max(6.0 * dx_um, 0.06),
        domain_height=substrate_top_um + 0.05 + max(6.0 * dx_um, 0.06),
        dx=dx_um, opening_width=opening_um, mask_thickness=0.05,
        substrate_top=substrate_top_um, etched_depth=initial_depth_um)
    domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(domain[2])
    reference_plane_m = source_z * geometry.mesh_length_unit_m
    boundary = PlasmaBoundaryState((
        ion_species(
            2e19, reference_plane_m, energy_eV=40.0, iad_sigma_deg=3.0,
            log2=velocity_log2, seed=seed + 1),
        thermal_neutral(
            "F", 19.0, 2e20, reference_plane_m,
            log2=velocity_log2, seed=seed + 2),
        thermal_neutral(
            "O", 16.0, 4e19, reference_plane_m,
            log2=velocity_log2, seed=seed + 3),
    ), reference_plane_m=reference_plane_m,
       provenance={"case": "deboer_2002_evolving_profile_development"})
    response = GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
        1, ion_species_name="ion")
    process = PhysicalProcess(
        geometry, boundary,
        {"ion": "energetic_bombardment", "F": "neutral_reactant",
         "O": "neutral_reactant"},
        build_common_belen_si_mechanism(s_F=s_f), (1,),
        duration_s, n_steps,
        (0.0, float(domain[0]), 0.0, float(domain[1])), source_z,
        solver_options={
            "n_position": n_position,
            "seed": seed,
            "cfl_number": 0.3,
            "reinitialize": True,
            "reinitialization_method": "cr2",
            "transport_device": "cpu",
            "neutral_radiosity_options": {
                "rays_per_face": radiosity_rays,
                "seed": seed + 100,
                "periodic_lateral": True,
                "domain_size": domain,
                "nonetchable_reaction_probability_by_material": {
                    2: {"F": 1e-3, "O": 1e-3}},
            },
            "neutral_surface_fixed_point_tolerance": 1e-3,
            "neutral_surface_fixed_point_max_iterations": 20,
            "charged_surface_response": response,
            "charged_surface_response_options": {
                "fixed_dt": 0.005,
                "max_steps": 256,
                "trajectory_adaptive_horizon": True,
                "trajectory_emergency_max_steps": 4096,
                "max_bounces": 64,
                "relative_tail_tolerance": 1e-8,
                "adaptive_bounce_extension": True,
                "emergency_max_bounces": 512,
                "periodic_lateral": True,
            },
        })
    started = perf_counter()
    initial_phi = np.asarray(geometry.phi).copy()
    result = process.run()
    wall_time_s = perf_counter() - started
    depth = [initial_depth_um]
    for step in result.steps:
        depth.append(_floor_depth_um(step.geometry, substrate_top_um))
    step_audit = [{
        "step": index + 1,
        "time_s": (index + 1) * duration_s / n_steps,
        "floor_depth_um": depth[index + 1],
        "surface_fixed_point_iterations": step.diagnostics[
            "neutral_surface_fixed_point_iterations"],
        "surface_fixed_point_residual": step.diagnostics[
            "neutral_surface_fixed_point_residual"],
        "response_bounces": step.diagnostics["charged_surface_response_bounces"],
        "response_reimpact_events": step.diagnostics[
            "charged_surface_response_reimpact_events"],
        "response_relative_charge_error": step.diagnostics[
            "charged_surface_response_relative_charge_error"],
        "response_maximum_energy_error": step.diagnostics[
            "charged_surface_response_maximum_energy_error"],
        "response_tail_l1_error_bound": step.diagnostics[
            "charged_surface_response_tail_l1_error_bound"],
        "state_remap_maximum_relative_conservation_residual": max(
            item["max_relative_conservation_residual"]
            for item in step.state_remap_diagnostics["materials"].values()),
    } for index, step in enumerate(result.steps)]
    return result, initial_phi, {
        "velocity_log2_samples": int(velocity_log2),
        "n_position": int(n_position),
        "radiosity_rays_per_face": int(radiosity_rays),
        "wall_time_s": wall_time_s,
        "floor_depth_um": depth,
        "final_floor_depth_um": depth[-1],
        "profile_increment_um": depth[-1] - depth[0],
        "step_audit": step_audit,
        "manifest": dict(result.run_manifest),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results" / "deboer_2002_evolving_profile_development")
    parser.add_argument("--static-curve", type=Path, default=DEFAULT_STATIC_CURVE)
    parser.add_argument("--opening-um", type=float, default=0.10)
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--initial-aspect-ratio", type=float, default=2.0)
    parser.add_argument("--duration-s", type=float, default=150.0)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--seed", type=int, default=41)
    args = parser.parse_args()
    if (args.opening_um <= 0.0 or args.dx_um <= 0.0
            or args.initial_aspect_ratio < 0.0 or args.duration_s <= 0.0
            or args.steps <= 0):
        raise ValueError("invalid moving-profile development controls")

    selected_s_f, curve = _curve_from_audit(args.static_curve)
    initial_depth_um = args.opening_um * args.initial_aspect_ratio
    static_depth_um = _static_depth_um(
        curve, opening_um=args.opening_um, initial_depth_um=initial_depth_um,
        duration_s=args.duration_s)
    levels = (
        {"name": "base", "velocity_log2": 9, "n_position": 16, "radiosity_rays": 32},
        {"name": "refined", "velocity_log2": 10, "n_position": 32, "radiosity_rays": 64},
        {"name": "fine", "velocity_log2": 11, "n_position": 64, "radiosity_rays": 128},
    )
    outputs = {}; fine_result = None; fine_initial_phi = None
    for level in levels:
        result, initial_phi, audit = _case(
            opening_um=args.opening_um, dx_um=args.dx_um,
            initial_aspect_ratio=args.initial_aspect_ratio,
            velocity_log2=level["velocity_log2"], n_position=level["n_position"],
            radiosity_rays=level["radiosity_rays"], seed=args.seed,
            duration_s=args.duration_s, n_steps=args.steps, s_f=selected_s_f)
        outputs[level["name"]] = audit
        fine_result = result; fine_initial_phi = initial_phi
        print(
            f"{level['name']}: depth {initial_depth_um:.6f} -> "
            f"{audit['final_floor_depth_um']:.6f} um in {audit['wall_time_s']:.2f} s",
            flush=True)

    base = outputs["base"]; refined = outputs["refined"]; fine = outputs["fine"]
    coarse_sampling_delta_um = (
        refined["final_floor_depth_um"] - base["final_floor_depth_um"])
    sampling_delta_um = (
        fine["final_floor_depth_um"] - refined["final_floor_depth_um"])
    evolving_minus_static_um = fine["final_floor_depth_um"] - static_depth_um
    maximum_charge_error = max(
        step["response_relative_charge_error"]
        for level in outputs.values() for step in level["step_audit"])
    maximum_energy_error = max(
        step["response_maximum_energy_error"]
        for level in outputs.values() for step in level["step_audit"])
    maximum_tail_bound = max(
        step["response_tail_l1_error_bound"]
        for level in outputs.values() for step in level["step_audit"])
    maximum_remap_error = max(
        step["state_remap_maximum_relative_conservation_residual"]
        for level in outputs.values() for step in level["step_audit"])
    report = {
        "schema": "deboer-2002-evolving-profile-development-v1",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "static_curve_path": str(args.static_curve.resolve().relative_to(ROOT)),
        "static_curve_sha256": _sha256(args.static_curve),
        "data_status": (
            "DEVELOPMENT_ONLY: moving-profile mechanism isolation; no experimental validation"),
        "configuration": {
            "opening_um": args.opening_um,
            "dx_um": args.dx_um,
            "initial_aspect_ratio": args.initial_aspect_ratio,
            "initial_depth_um": initial_depth_um,
            "duration_s": args.duration_s,
            "steps": args.steps,
            "seed": args.seed,
            "selected_s_F_from_single_legal_calibration": selected_s_f,
            "ion_energy_eV": 40.0,
            "iad_component_sigma_deg": 3.0,
            "ion_reflection": dict(
                GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
                    1, ion_species_name="ion").provenance),
            "photon_channel": (
                "absent: Figure 9 declares neither photon flux/spectrum nor a photon-assisted "
                "silicon yield/state law"),
        },
        "static_rate_counterfactual_final_depth_um": static_depth_um,
        "levels": outputs,
        "sampling_refinement_delta_um": sampling_delta_um,
        "coarse_sampling_refinement_delta_um": coarse_sampling_delta_um,
        "sampling_delta_contraction_ratio": (
            abs(sampling_delta_um) / max(abs(coarse_sampling_delta_um), 1e-300)),
        "evolving_minus_static_um": evolving_minus_static_um,
        "moving_vs_static_resolved_at_latest_sampling_delta": bool(
            abs(evolving_minus_static_um) > 2.0 * abs(sampling_delta_um)),
        "maximum_response_relative_charge_error": maximum_charge_error,
        "maximum_response_energy_error": maximum_energy_error,
        "maximum_response_tail_l1_error_bound": maximum_tail_bound,
        "maximum_state_remap_relative_conservation_error": maximum_remap_error,
        "conservation_gate_passed": bool(
            maximum_charge_error <= 1e-10 and maximum_energy_error <= 1e-10
            and maximum_tail_bound <= 2e-8 and maximum_remap_error <= 1e-8),
        "earned_conclusion": (
            "quantifies moving-profile/history departure from the identical static common operator; "
            "does not score or validate Figure-9 experimental markers"),
    }
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    _atomic_json(output / "audit.json", report)
    np.savez_compressed(
        output / "profiles.npz",
        initial_phi=np.asarray(fine_initial_phi),
        final_phi=np.asarray(fine_result.geometry.phi),
        material_id=np.asarray(fine_result.geometry.material_id),
        dx_um=np.asarray(args.dx_um),
        mesh_length_unit_m=np.asarray(fine_result.geometry.mesh_length_unit_m))
    print(json.dumps({
        "static_final_depth_um": static_depth_um,
        "fine_evolving_final_depth_um": fine["final_floor_depth_um"],
        "evolving_minus_static_um": evolving_minus_static_um,
        "sampling_refinement_delta_um": sampling_delta_um,
        "conservation_gate_passed": report["conservation_gate_passed"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
