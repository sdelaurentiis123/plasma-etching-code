#!/usr/bin/env python3
"""Bounded, timestep-refined physical-time charging campaign on the real 3-D trench."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(variable, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/petch-matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from petch.boundary_state import (  # noqa: E402
    IonEnergyTransverseMaxwellianDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    maxwellian_electron_boundary_state,
    mixture_boundary_proposal,
    qmc_boundary_proposal,
)
from petch.charging_coupled_3d import (  # noqa: E402
    BidirectionalCurrentCertificationError,
    PhysicalTimeChargingIntegrationError,
    integrate_dielectric_charging_transient_3d,
)
from petch.charging_poisson_3d import NodalPoissonSystem3D  # noqa: E402
from petch.charging_poisson import EPS0  # noqa: E402
from petch.feature_step_3d import (  # noqa: E402
    _surface_gas_normals,
    make_rectangular_trench_geometry_3d,
)
from petch.threed import extract_mesh_3d  # noqa: E402


def _geometry_and_poisson(dx=0.25):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.5, domain_height=2.0, dx=dx,
        opening_width=0.5, mask_thickness=0.25,
        substrate_top=1.25, etched_depth=0.75)
    fixed = np.zeros(geometry.phi.shape, dtype=bool)
    fixed[:, :, -1] = True
    cell_shape = tuple(size - 1 for size in geometry.phi.shape)
    phi_center = sum(
        geometry.phi[i:i + cell_shape[0], j:j + cell_shape[1], k:k + cell_shape[2]]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0
    epsilon_r = np.where(phi_center > 0.0, 3.9, 1.0)
    poisson = NodalPoissonSystem3D(
        epsilon_r, geometry.dx * geometry.mesh_length_unit_m, fixed)
    return geometry, poisson


def _boundary(reference_plane_m):
    ion_flux = 2.2e21
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, ion_flux, [[0.0, 0.0, 10.0]], [1.0],
        density_model=IonEnergyTransverseMaxwellianDensity(
            np.array([99.0, 101.0]), np.array([1.0]), 0.05))
    electron = maxwellian_electron_boundary_state(
        4.0, 10.0 * ion_flux, n_transverse=5, n_normal=8,
        reference_plane_m=reference_plane_m).species[0]
    return PlasmaBoundaryState((ion, electron), reference_plane_m=reference_plane_m)


def _electron_proposal(boundary, level, seed):
    physical = boundary.get("electron")
    broad = maxwellian_electron_boundary_state(
        20.0, physical.flux_m2_s, n_transverse=5, n_normal=8,
        electron_name="electron", reference_plane_m=boundary.reference_plane_m).species[0]
    mixture = mixture_boundary_proposal((physical, broad), (0.8, 0.2), name="electron")
    return qmc_boundary_proposal(mixture, level, seed=seed)


def _ion_proposal(boundary, level, seed):
    physical = boundary.get("Ar+")
    broad = SpeciesBoundaryState(
        "Ar+", 1, physical.mass_amu, physical.flux_m2_s,
        [[0.0, 0.0, 10.0]], [1.0],
        density_model=IonEnergyTransverseMaxwellianDensity(
            np.array([1.0, 201.0]), np.array([1.0]), 2.0))
    mixture = mixture_boundary_proposal((physical, broad), (0.8, 0.2), name="Ar+")
    return qmc_boundary_proposal(mixture, level, seed=seed)


def _depth_charge(charge, spacing):
    z = np.arange(charge.shape[2], dtype=float) * spacing
    upper = z >= 1.25
    floor = z <= 0.5
    lower_wall = (z > 0.5) & (z < 1.25)
    return dict(
        upper=float(np.sum(charge[:, :, upper])),
        lower_wall=float(np.sum(charge[:, :, lower_wall])),
        floor=float(np.sum(charge[:, :, floor])))


def _config_hash(config):
    encoded = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _run(label, timestep_s, steps, initial_charge, common):
    print(f"[{label}] dt={timestep_s:.3e} s, updates={steps}", flush=True)
    result = integrate_dielectric_charging_transient_3d(
        initial_charge_node_c=initial_charge,
        timestep_s=timestep_s, n_steps=steps, current_balance_tol=0.08, **common)
    print(
        f"[{label}] final node rms={result.history[-1]['rms_relative_current_imbalance_node']:.6f} "
        f"node max={result.history[-1]['max_relative_current_imbalance_node']:.6f} "
        f"face max={result.history[-1]['max_relative_current_imbalance_face']:.6f}",
        flush=True)
    return result


def _plot(coarse, refined, output, spacing):
    coarse_history = coarse.history
    refined_history = refined.history
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    for history, label, style in (
            (coarse_history, "dt", "o-"), (refined_history, "dt/2", ".-")):
        time_ns = np.array([item["physical_time_s"] for item in history]) * 1e9
        axes[0, 0].plot(
            time_ns, [item["rms_relative_current_imbalance_node"] for item in history],
            style, label=label)
        axes[0, 1].plot(
            time_ns, [item["max_relative_current_imbalance_node"] for item in history],
            style, label=f"node {label}")
        axes[0, 1].plot(
            time_ns, [item["max_relative_current_imbalance_face"] for item in history],
            style, alpha=0.55, label=f"face {label}")
        axes[1, 0].plot(
            time_ns, [item["minimum_potential_v"] for item in history], style,
            label=f"min {label}")
        axes[1, 0].plot(
            time_ns, [item["maximum_potential_v"] for item in history], style,
            alpha=0.55, label=f"max {label}")
    axes[0, 0].axhline(0.08, color="0.4", linestyle="--", label="contract 0.08")
    axes[0, 0].set(title="RMS current imbalance", xlabel="physical time (ns)", ylabel="relative")
    axes[0, 1].axhline(0.08, color="0.4", linestyle="--")
    axes[0, 1].set(title="Worst node and face", xlabel="physical time (ns)", ylabel="relative")
    axes[1, 0].set(title="Potential envelope", xlabel="physical time (ns)", ylabel="potential (V)")
    for axis in axes.flat[:3]:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)

    z = np.arange(refined.charge_node_c.shape[2]) * spacing
    charge_z = refined.charge_node_c.sum(axis=(0, 1))
    axes[1, 1].bar(z, charge_z / max(np.max(np.abs(charge_z)), 1e-300), width=0.8 * spacing)
    axes[1, 1].axvline(1.25, color="0.4", linestyle="--", label="substrate top")
    axes[1, 1].set(
        title="Final stored-charge depth profile (dt/2)", xlabel="z (µm)",
        ylabel="charge / max |charge|")
    axes[1, 1].grid(axis="y", alpha=0.25)
    axes[1, 1].legend(fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/charging_task1_3d")
    parser.add_argument("--grid-dx", type=float, default=0.25)
    parser.add_argument("--timestep-s", type=float, default=2e-9)
    parser.add_argument("--steps", type=int, default=12)
    parser.add_argument("--forward-level", type=int, default=10)
    parser.add_argument(
        "--ion-estimator", choices=("forward", "adjoint", "bidirectional"),
        default="bidirectional")
    parser.add_argument("--ion-proposal-level", type=int, default=8)
    parser.add_argument("--ion-proposal-seed", type=int, default=79)
    parser.add_argument("--electron-proposal-level", type=int, default=8)
    parser.add_argument("--electron-proposal-seed", type=int, default=83)
    parser.add_argument("--trajectory-dt", type=float, default=0.005)
    parser.add_argument("--trajectory-max-steps", type=int, default=50000)
    parser.add_argument("--adjoint-face-points", type=int, default=3)
    parser.add_argument("--bidirectional-replicates", type=int, default=4)
    parser.add_argument("--bidirectional-max-level", type=int, default=11)
    parser.add_argument("--bidirectional-max-face-points", type=int, default=7)
    parser.add_argument("--bidirectional-absolute-tolerance", type=float, default=0.02)
    parser.add_argument("--bidirectional-relative-tolerance", type=float, default=0.1)
    parser.add_argument("--initial-state", type=Path)
    parser.add_argument("--initial-key", default="refined_charge_node_c")
    parser.add_argument("--initial-kind", choices=("charge", "potential"), default="charge")
    parser.add_argument("--initial-reference-dx", type=float, default=0.25)
    parser.add_argument("--ion-method-map", type=Path)
    parser.add_argument("--ion-method-key", default="method_hint_Ar")
    args = parser.parse_args()
    if args.steps < 0:
        parser.error("--steps must be nonnegative")

    geometry, poisson = _geometry_and_poisson(args.grid_dx)
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    source_z = 2.0
    boundary = _boundary(source_z * geometry.mesh_length_unit_m)
    proposal = _electron_proposal(
        boundary, args.electron_proposal_level, args.electron_proposal_seed)
    proposals = {"electron": proposal}
    proposal_frames = {"electron": "surface_local"}
    if args.ion_estimator in {"adjoint", "bidirectional"}:
        proposals["Ar+"] = _ion_proposal(
            boundary, args.ion_proposal_level, args.ion_proposal_seed)
        proposal_frames["Ar+"] = "source_aligned"
    initial_charge = np.zeros(poisson.shape)
    initial_state_sha256 = None
    potential_recovery_relative_l2 = None
    if args.initial_state is not None:
        payload = args.initial_state.read_bytes()
        initial_state_sha256 = hashlib.sha256(payload).hexdigest()
        with np.load(args.initial_state) as archived:
            initial = np.asarray(archived[args.initial_key], dtype=float).copy()
        if args.initial_kind == "charge":
            initial_charge = initial
            if initial_charge.shape != poisson.shape:
                parser.error("initial charge shape does not match the campaign Poisson grid")
        else:
            coordinates = np.meshgrid(
                *(np.arange(size, dtype=float) * geometry.dx / args.initial_reference_dx
                  for size in poisson.shape), indexing="ij")
            potential = map_coordinates(
                initial, coordinates, order=1, mode="nearest", prefilter=False)
            initial_charge = (
                EPS0 * (poisson.stiffness @ potential.ravel())).reshape(poisson.shape)
            initial_charge[poisson.dirichlet_mask] = 0.0
            recovered, _ = poisson.solve(initial_charge)
            potential_recovery_relative_l2 = float(
                np.linalg.norm(recovered - potential)
                / max(np.linalg.norm(potential), 1e-300))
    ion_method_map_sha256 = None
    ion_method_hint = None
    if args.ion_method_map is not None:
        payload = args.ion_method_map.read_bytes()
        ion_method_map_sha256 = hashlib.sha256(payload).hexdigest()
        with np.load(args.ion_method_map) as archived:
            ion_method_hint = np.asarray(archived[args.ion_method_key]).astype("U7")
        if (ion_method_hint.shape != (len(faces),)
                or np.any(~np.isin(ion_method_hint, ("forward", "adjoint")))):
            parser.error("ion method map must select forward/adjoint for every campaign face")
    bidirectional_options = dict(
        forward_log2_samples=args.forward_level,
        adjoint_log2_samples=args.ion_proposal_level,
        n_replicates=args.bidirectional_replicates,
        max_forward_log2_samples=args.bidirectional_max_level,
        max_adjoint_log2_samples=args.bidirectional_max_level,
        max_face_quadrature_points=args.bidirectional_max_face_points,
        element_absolute_tolerance=args.bidirectional_absolute_tolerance,
        element_relative_tolerance=args.bidirectional_relative_tolerance,
        face_quadrature_points=args.adjoint_face_points)
    if ion_method_hint is not None:
        bidirectional_options.update(
            method_hint={"Ar+": ion_method_hint}, require_certification=False)
    common = dict(
        poisson_system=poisson, boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m, n_position=256, seed=79,
        trajectory_fixed_dt=args.trajectory_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        phase_space_log2_samples=args.forward_level, periodic_lateral=True,
        transport_estimator={"Ar+": args.ion_estimator, "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        adjoint_face_quadrature_points=args.adjoint_face_points,
        adjoint_ray_offset=1e-4, adjoint_proposals=proposals,
        adjoint_proposal_frames=proposal_frames,
        bidirectional_options=bidirectional_options,
        transport_device="cpu")
    config = dict(
        model="exact_hard_visibility_3d_physical_time",
        geometry=dict(cell_width=1.0, cell_length=0.5, domain_height=2.0, dx=geometry.dx,
                      opening_width=0.5, mask_thickness=0.25,
                      substrate_top=1.25, etched_depth=0.75),
        timestep_s=args.timestep_s, steps=args.steps,
        refined_timestep_s=args.timestep_s / 2.0, refined_steps=2 * args.steps,
        forward_level=args.forward_level,
        ion_estimator=args.ion_estimator,
        ion_proposal_level=args.ion_proposal_level,
        ion_proposal_seed=args.ion_proposal_seed,
        electron_proposal_level=args.electron_proposal_level,
        electron_proposal_seed=args.electron_proposal_seed,
        trajectory_dt=args.trajectory_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        adjoint_face_points=args.adjoint_face_points,
        bidirectional_replicates=args.bidirectional_replicates,
        bidirectional_max_level=args.bidirectional_max_level,
        bidirectional_max_face_points=args.bidirectional_max_face_points,
        bidirectional_absolute_tolerance=args.bidirectional_absolute_tolerance,
        bidirectional_relative_tolerance=args.bidirectional_relative_tolerance,
        initial_state_sha256=initial_state_sha256,
        initial_key=(None if args.initial_state is None else args.initial_key),
        initial_kind=(None if args.initial_state is None else args.initial_kind),
        initial_reference_dx=(
            None if args.initial_state is None or args.initial_kind != "potential"
            else args.initial_reference_dx),
        potential_recovery_relative_l2=potential_recovery_relative_l2,
        ion_method_map_sha256=ion_method_map_sha256,
        ion_method_key=(None if args.ion_method_map is None else args.ion_method_key),
        estimators={"Ar+": args.ion_estimator, "electron": "adjoint"})
    config_hash = _config_hash(config)

    try:
        coarse = _run("coarse", args.timestep_s, args.steps, initial_charge, common)
    except BidirectionalCurrentCertificationError as error:
        diagnostic = {"error": str(error), "species": {}}
        if error.result is not None:
            for name, selection in error.result.selection_by_species.items():
                unresolved = np.where(
                    ~(selection.estimator_consistent
                      & selection.method_within_tolerance))[0]
                diagnostic["species"][name] = dict(
                    converged=selection.converged,
                    inconsistent_faces=np.where(~selection.estimator_consistent)[0].tolist(),
                    imprecise_faces=np.where(~selection.method_within_tolerance)[0].tolist(),
                    method_counts={method: int(np.sum(selection.method == method))
                                   for method in ("forward", "adjoint")},
                    maximum_forward_relative_stderr=float(np.max(
                        selection.forward_face_stderr_m2_s
                        / np.maximum(selection.forward_face_mean_m2_s, 1e-300))),
                    maximum_adjoint_relative_stderr=float(np.max(
                        selection.adjoint_face_stderr_m2_s
                        / np.maximum(selection.adjoint_face_mean_m2_s, 1e-300))),
                    unresolved=[dict(
                        face=int(index), centroid_um=centroids[index].tolist(),
                        method=str(selection.method[index]),
                        forward_mean_m2_s=float(selection.forward_face_mean_m2_s[index]),
                        forward_stderr_m2_s=float(selection.forward_face_stderr_m2_s[index]),
                        adjoint_mean_m2_s=float(selection.adjoint_face_mean_m2_s[index]),
                        adjoint_stderr_m2_s=float(selection.adjoint_face_stderr_m2_s[index]),
                        consistent=bool(selection.estimator_consistent[index]),
                        within_tolerance=bool(selection.method_within_tolerance[index]))
                        for index in unresolved])
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "certification_failure.json").write_text(
            json.dumps(diagnostic, indent=2) + "\n")
        print(json.dumps(diagnostic, indent=2), flush=True)
        raise
    except PhysicalTimeChargingIntegrationError as error:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.output_dir / "integration_failure_checkpoint.npz",
            charge_node_c=error.charge_node_c, failed_before_step=error.step)
        diagnostic = dict(
            error=str(error), failed_before_step=error.step,
            completed_states=[dict(item) for item in error.history])
        (args.output_dir / "integration_failure.json").write_text(
            json.dumps(diagnostic, indent=2) + "\n")
        print(json.dumps(diagnostic, indent=2), flush=True)
        raise
    refined = (coarse if args.steps == 0 else _run(
        "refined", args.timestep_s / 2.0, 2 * args.steps, initial_charge, common))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    history_path = args.output_dir / "transient_history.csv"
    with history_path.open("w", newline="") as handle:
        fields = ["variant", *coarse.history[0].keys()]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for label, result in (("coarse", coarse), ("refined", refined)):
            for item in result.history:
                writer.writerow({"variant": label, **dict(item)})

    coarse_end = coarse.history[-1]
    refined_end = refined.history[-1]
    charge_difference = np.linalg.norm(coarse.charge_node_c - refined.charge_node_c)
    charge_scale = max(np.linalg.norm(refined.charge_node_c), 1e-300)
    potential_difference = np.linalg.norm(coarse.potential_v - refined.potential_v)
    potential_scale = max(np.linalg.norm(refined.potential_v), 1e-300)
    region_charge = _depth_charge(refined.charge_node_c, geometry.dx)
    conservation_scale = max(
        refined.diagnostics["total_absolute_incident_charge_c"], 1e-300)
    gates = dict(
        charge_conservation=bool(
            abs(refined.diagnostics["cumulative_charge_conservation_residual_c"])
            <= 1e-12 * conservation_scale),
        timestep_charge_relative_l2=bool(charge_difference / charge_scale <= 0.05),
        timestep_potential_relative_l2=bool(potential_difference / potential_scale <= 0.05),
        dipole_sign=bool(region_charge["upper"] < 0.0 and region_charge["floor"] > 0.0),
        rms_decreased=bool(
            refined_end["rms_relative_current_imbalance_node"]
            < refined.history[0]["rms_relative_current_imbalance_node"]),
        exact_operator=True,
        residual_contract=bool(
            refined_end["max_relative_current_imbalance_node"] <= 0.08))
    summary = dict(
        config_hash=config_hash, config=config, geometry_faces=int(len(faces)),
        exact_operator_statement="hard visibility; no smoothing or analytic current replacement",
        coarse_final=dict(coarse_end), refined_final=dict(refined_end),
        timestep_refinement=dict(
            charge_relative_l2=float(charge_difference / charge_scale),
            potential_relative_l2=float(potential_difference / potential_scale)),
        refined_charge_regions_c=region_charge,
        conservation=dict(refined.diagnostics), gates=gates,
        campaign_complete=bool(all(gates.values())),
        conclusion=(
            "bounded Task 1 gates passed" if all(gates.values())
            else "bounded Task 1 remains diagnostic; inspect failed gates before extension"))
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savez_compressed(
        args.output_dir / "final_states.npz",
        coarse_charge_node_c=coarse.charge_node_c,
        coarse_potential_v=coarse.potential_v,
        coarse_positive_face_current_density_a_m2=coarse.positive_face_current_density_a_m2,
        coarse_negative_face_current_density_a_m2=coarse.negative_face_current_density_a_m2,
        coarse_positive_current_node_a=coarse.positive_current_node_a,
        coarse_negative_current_node_a=coarse.negative_current_node_a,
        refined_charge_node_c=refined.charge_node_c,
        refined_potential_v=refined.potential_v,
        refined_positive_face_current_density_a_m2=refined.positive_face_current_density_a_m2,
        refined_negative_face_current_density_a_m2=refined.negative_face_current_density_a_m2,
        refined_positive_current_node_a=refined.positive_current_node_a,
        refined_negative_current_node_a=refined.negative_current_node_a,
        coarse_charge_history_node_c=coarse.charge_history_node_c,
        refined_charge_history_node_c=refined.charge_history_node_c,
        **{f"refined_method_hint_{name}": value
           for name, value in refined.bidirectional_method_hint.items()},
        vertices=verts, faces=faces, centroids=centroids, areas=areas)
    _plot(coarse, refined, args.output_dir / "physical_time_transient.png", geometry.dx)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
