#!/usr/bin/env python3
"""Bounded charging-on/off diagnostic for the measured-depth Jeon trench geometry.

This is a mechanism diagnostic, not an experimental-validation result.  It reuses the production
physical-time charging operator on one fixed geometry and pairs the zero-charge and charged endpoint
audits with the same unused scramble epoch.  No root solver, profile calibration, or held-out depth
observation enters the charge trajectory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from jeon_unified_baseline import DATA, baseline_mechanism
from petch.boundary_state import (
    IonEnergyTransverseMaxwellianDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    maxwellian_electron_boundary_state,
)
from petch.charging_coevolution_3d import integrate_surface_charging_to_saturation_3d
from petch.charging_poisson_3d import NodalPoissonSystem3D
from petch.experimental_boundary import Jeon2022BoundaryClosure, build_jeon_2022_boundary_state
from petch.experimental_data import (
    load_jeon_2022_electron_bias_controls,
    load_jeon_2022_plasma_controls,
)
from petch.feature_step_3d import (
    _face_material_ids,
    _surface_gas_normals,
    make_rectangular_trench_geometry_3d,
)
from petch.threed import extract_mesh_3d


def _cell_center_average(field):
    value = np.asarray(field, dtype=float)
    return sum(
        value[i:i + value.shape[0] - 1,
              j:j + value.shape[1] - 1,
              k:k + value.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0


def _poisson_system(geometry, mask_epsilon_r):
    """Q1 dielectric feature cell with the plasma/source plane held at zero potential."""
    fields = {
        material: _cell_center_average(levelset)
        for material, levelset in geometry.material_levelsets.items()
    }
    materials = np.asarray(sorted(fields), dtype=int)
    stack = np.stack([fields[int(material)] for material in materials])
    owner = materials[np.argmax(stack, axis=0)]
    solid = np.max(stack, axis=0) >= 0.0
    epsilon_r = np.ones(owner.shape)
    epsilon_r[solid & (owner == 1)] = 3.9
    epsilon_r[solid & (owner == 2)] = float(mask_epsilon_r)
    fixed = np.zeros(geometry.phi.shape, dtype=bool)
    fixed[:, :, -1] = True
    return NodalPoissonSystem3D(
        epsilon_r, geometry.dx * geometry.mesh_length_unit_m, fixed,
        periodic_axes=(0, 1))


def _charged_boundary(source_z_m):
    plasma = next(
        item for item in load_jeon_2022_plasma_controls(
            DATA / "digitized_plasma_controls.csv")
        if item.condition_family == "gas_fraction_cw" and item.c4f8_fraction == 0.2)
    electron_control = next(
        item for item in load_jeon_2022_electron_bias_controls(
            DATA / "digitized_electron_bias_controls.csv")
        if item.condition_family == "gas_fraction_cw" and item.c4f8_fraction == 0.2)
    closure = Jeon2022BoundaryClosure(
        ion_name="Ar+", ion_mass_amu=39.948,
        ion_normal_energy_eV=[electron_control.self_bias_magnitude_v],
        ion_normal_energy_weight=[1.0], ion_tangential_temperature_eV=0.026,
        neutral_flux_fraction={"FC_total": 1.0},
        neutral_mass_amu={"FC_total": 50.005}, neutral_temperature_K=300.0,
        provenance={
            "model": "Jeon diagnostic self-bias energy closure",
            "warning": "self-bias is not a measured IEDF"},
        supports_prediction_within_declared_domain=False)
    full = build_jeon_2022_boundary_state(
        plasma, electron_control, closure, reference_plane_m=source_z_m,
        n_transverse_ion=3, n_transverse_neutral=3, n_normal_neutral=4)
    discrete_ion = full.get("Ar+")
    # The Jeon adapter intentionally accepts tabulated IEDFs, including a one-point closure.
    # Randomized-QMC forward transport additionally needs a continuous density.  For this bounded
    # diagnostic, interpret the measured self-bias uncertainty as a one-bin ion-energy density;
    # this changes neither the engine nor the central energy and avoids inventing an RF-sheath IEDF.
    ion_energy_half_width_eV = max(
        float(electron_control.self_bias_digitization_uncertainty_v), 1.0)
    ion_density = IonEnergyTransverseMaxwellianDensity(
        normal_energy_edges_eV=[
            max(0.0, electron_control.self_bias_magnitude_v - ion_energy_half_width_eV),
            electron_control.self_bias_magnitude_v + ion_energy_half_width_eV,
        ],
        probability_mass=[1.0],
        tangential_temperature_eV=closure.ion_tangential_temperature_eV,
    )
    ion = SpeciesBoundaryState(
        name=discrete_ion.name,
        charge_number=discrete_ion.charge_number,
        mass_amu=discrete_ion.mass_amu,
        flux_m2_s=discrete_ion.flux_m2_s,
        velocity_sqrt_eV=discrete_ion.velocity_sqrt_eV,
        weight=discrete_ion.weight,
        density_model=ion_density,
        provenance={
            **dict(discrete_ion.provenance),
            "continuous_qmc_closure": "uniform ion energy over self-bias digitization uncertainty",
            "ion_energy_half_width_eV": ion_energy_half_width_eV,
        },
    )
    electron = maxwellian_electron_boundary_state(
        3.0, ion.flux_m2_s, n_transverse=5, n_normal=8,
        reference_plane_m=source_z_m).get("electron")
    return PlasmaBoundaryState(
        (ion, electron), source_z_m,
        provenance={
            "source": "Jeon 2022 diagnostics",
            "electron_temperature_eV": 3.0,
            "equal_ion_electron_particle_flux": True,
            "self_bias_is_not_iedf": True,
            "purpose": "bounded charging mechanism diagnostic"})


def _ion_flux_by_face(transport, face_count):
    value = np.zeros(int(face_count))
    for population in transport.surface_fluxes.energetic_fluxes:
        if population.name == "Ar+":
            value += np.asarray(population.flux_m2_s, dtype=float)
    return value


def _floor_flux(transport, centroids, normals, material, areas, floor_z, dx):
    selected = (
        (material == 1)
        & (normals[:, 2] > 0.5)
        & (centroids[:, 2] <= float(floor_z) + 1.5 * float(dx)))
    if not np.any(selected):
        raise RuntimeError("measured-depth geometry has no resolved SiO2 floor faces")
    flux = _ion_flux_by_face(transport, len(areas))
    return float(np.dot(flux[selected], areas[selected]) / np.sum(areas[selected]))


def _atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--width-nm", type=float, default=60.0)
    parser.add_argument("--measured-depth-nm", type=float, default=295.8)
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--timestep-ns", type=float, default=125.0)
    parser.add_argument("--phase-space-log2-samples", type=int, default=8)
    parser.add_argument("--n-position", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1702)
    parser.add_argument("--mask-epsilon-r", type=float, default=4.0)
    parser.add_argument("--audit-epoch", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    substrate_top = 2.4
    mask_thickness = 0.7
    domain_height = 3.3
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=0.1, domain_height=domain_height,
        dx=args.dx_um, opening_width=args.width_nm * 1e-3,
        mask_thickness=mask_thickness, substrate_top=substrate_top,
        etched_depth=args.measured_depth_nm * 1e-3)
    source_z = (geometry.phi.shape[2] - 1) * geometry.dx
    boundary = _charged_boundary(source_z * geometry.mesh_length_unit_m)
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    material = _face_material_ids(centroids, geometry)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    poisson = _poisson_system(geometry, args.mask_epsilon_r)
    initial_sigma = np.zeros(len(faces))
    floor_z = substrate_top - args.measured_depth_nm * 1e-3

    common = dict(
        poisson_system=poisson, boundary=boundary, verts=verts, faces=faces, areas=areas,
        face_centroids=centroids, face_gas_normals=normals, face_material_id=material,
        source_bounds=(0.0, 0.5, 0.0, 0.1), source_z=source_z,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        patch_scales_m=(max(2.0 * args.dx_um, 0.04) * 1e-6, 0.2e-6),
        potential_rate_tolerance_v_s=1.0e3,
        timestep_s=args.timestep_ns * 1e-9,
        current_balance_tolerance=0.08, timestep_policy="fixed",
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        mesh_origin_m=geometry.mesh_origin_m,
        n_position=args.n_position, seed=args.seed,
        trajectory_fixed_dt=0.005, trajectory_max_steps=1024,
        trajectory_adaptive_horizon=True, trajectory_emergency_max_steps=8192,
        phase_space_log2_samples=args.phase_space_log2_samples,
        periodic_lateral=True, transport_estimator="forward",
        transport_device="cpu", stop_on_saturation=False,
        scramble_mode="fresh", compatible_q1_charge_state=True)

    started = perf_counter()

    def progress(*, potential_v, history_item, accepted_steps, physical_time_s, **_):
        evaluated_step = int(accepted_steps) + 1
        if (evaluated_step == 1 or evaluated_step % 5 == 0
                or evaluated_step == args.steps + 1):
            print(
                f"charging evaluation {evaluated_step}/{args.steps + 1}: "
                f"accepted={accepted_steps}, t={physical_time_s * 1e6:.3f} us, "
                f"max|V|={np.max(np.abs(potential_v)):.4g} V, "
                f"node_rms={history_item['rms_relative_current_imbalance_node']:.4g}",
                flush=True)

    march = integrate_surface_charging_to_saturation_3d(
        initial_sigma_c_per_m2=initial_sigma, maximum_steps=args.steps,
        initial_sampling_epoch=0, progress_callback=progress, **common)
    zero_audit = integrate_surface_charging_to_saturation_3d(
        initial_sigma_c_per_m2=initial_sigma, maximum_steps=0,
        initial_sampling_epoch=args.audit_epoch, **common)
    charged_audit = integrate_surface_charging_to_saturation_3d(
        initial_sigma_c_per_m2=march.sigma_c_per_m2, maximum_steps=0,
        initial_sampling_epoch=args.audit_epoch, **common)

    zero_floor = _floor_flux(
        zero_audit.final_step.transport, centroids, normals, material, areas,
        floor_z, geometry.dx)
    charged_floor = _floor_flux(
        charged_audit.final_step.transport, centroids, normals, material, areas,
        floor_z, geometry.dx)
    output = {
        "schema": "petch.jeon-2022.charging-diagnostic.v1",
        "scientific_status": "development mechanism diagnostic; not validation",
        "geometry": {
            "width_nm": args.width_nm,
            "measured_depth_nm": args.measured_depth_nm,
            "oxide_thickness_um": substrate_top,
            "mask_thickness_um": mask_thickness,
            "dx_um": args.dx_um,
            "resolved_opening_cells": args.width_nm * 1e-3 / args.dx_um,
            "mask_epsilon_r_closure": args.mask_epsilon_r,
        },
        "integration": {
            "steps": args.steps,
            "timestep_ns": args.timestep_ns,
            "physical_time_us": march.physical_time_s * 1e6,
            "scramble_mode": "fresh",
            "phase_space_log2_samples": args.phase_space_log2_samples,
            "n_position": args.n_position,
            "paired_audit_epoch": args.audit_epoch,
            "wall_time_s": perf_counter() - started,
        },
        "paired_exact_operator_audit": {
            "zero_charge_floor_ion_flux_m2_s": zero_floor,
            "charged_floor_ion_flux_m2_s": charged_floor,
            "charged_over_zero_floor_ion_flux": charged_floor / zero_floor,
            "maximum_absolute_potential_v": float(np.max(np.abs(charged_audit.potential_v))),
            "retained_node_rms_relative_current_imbalance": charged_audit.diagnostics[
                "retained_node_rms_relative_current_imbalance"],
            "retained_node_max_relative_current_imbalance": charged_audit.diagnostics[
                "retained_node_max_relative_current_imbalance"],
            "converged": charged_audit.converged,
        },
        "conservation": {
            "maximum_charge_conservation_relative_error": max(
                item["charge_conservation_relative_error"] for item in march.history),
            "maximum_face_to_node_update_relative_error": max(
                item["face_to_node_update_relative_error"] for item in march.history),
            "maximum_lineage_replay_fraction": max(
                item["transport_lineage_replay_fraction"] for item in march.history),
        },
    }
    _atomic_json(args.output, output)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
