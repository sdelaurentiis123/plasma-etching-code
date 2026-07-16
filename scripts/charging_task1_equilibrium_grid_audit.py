#!/usr/bin/env python3
"""Exact-operator charging balance audit on raw faces, Q1 nodes, and physical patches.

The refined grid evaluates the same trilinearly interpolated voltage state as the archived coarse
endpoint. It therefore measures operator/discretization sensitivity before attempting a costly
refined-grid transient; it does not claim that the interpolated state is a refined-grid equilibrium.
"""
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
    current_balance_metrics_3d,
    integrate_dielectric_charging_transient_3d,
)
from petch.charging_poisson import EPS0  # noqa: E402
from petch.charging_poisson_3d import NodalPoissonSystem3D  # noqa: E402
from petch.feature_step_3d import (  # noqa: E402
    _surface_gas_normals,
    make_rectangular_trench_geometry_3d,
)
from petch.threed import extract_mesh_3d  # noqa: E402


GEOMETRY = dict(
    cell_width=1.0, cell_length=0.5, domain_height=2.0,
    opening_width=0.5, mask_thickness=0.25,
    substrate_top=1.25, etched_depth=0.75)


def _hash_json(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _geometry_and_poisson(dx):
    geometry = make_rectangular_trench_geometry_3d(dx=dx, **GEOMETRY)
    fixed = np.zeros(geometry.phi.shape, dtype=bool)
    fixed[:, :, -1] = True
    cell_shape = tuple(size - 1 for size in geometry.phi.shape)
    phi_center = sum(
        geometry.phi[i:i + cell_shape[0], j:j + cell_shape[1], k:k + cell_shape[2]]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0
    epsilon_r = np.where(phi_center > 0.0, 3.9, 1.0)
    poisson = NodalPoissonSystem3D(
        epsilon_r, geometry.dx * geometry.mesh_length_unit_m, fixed,
        periodic_axes=(0, 1))
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    return geometry, poisson, verts, faces, centroids, areas, normals


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


def _interpolate_reference_potential(reference, reference_dx, shape, dx):
    coordinates = np.meshgrid(
        *(np.arange(size, dtype=float) * dx / reference_dx for size in shape),
        indexing="ij")
    return map_coordinates(
        np.asarray(reference, dtype=float), coordinates, order=1,
        mode="nearest", prefilter=False)


def _charge_reproducing_potential(poisson, potential):
    flat_potential = np.asarray(potential, dtype=float).ravel()
    charge = (EPS0 * (poisson.stiffness @ flat_potential)).reshape(poisson.shape)
    charge[poisson.dirichlet_mask] = 0.0
    recovered, _ = poisson.solve(charge)
    scale = max(float(np.linalg.norm(potential)), 1e-300)
    return charge, float(np.linalg.norm(recovered - potential) / scale)


def _patch_groups(centroids, normals, width):
    """Return fixed-physical-width surface patches, pooled across the invariant y direction."""
    centroids = np.asarray(centroids, dtype=float)
    normals = np.asarray(normals, dtype=float)
    if normals.shape != centroids.shape:
        raise ValueError("face centroids and normals must match")
    # Classify by the nearest analytic trench surface. Marching-cubes corner triangles are diagonal,
    # so dominant-normal classification changes ownership merely because a tie breaks differently
    # after refinement. Nearest declared surface gives the same physical partition on both grids.
    surface_distance = np.column_stack((
        np.abs(centroids[:, 2] - 0.5),
        np.abs(centroids[:, 0] - 0.25),
        np.abs(centroids[:, 0] - 0.75),
        np.abs(centroids[:, 2] - 1.5),
    ))
    region_names = np.asarray(("floor", "left-wall", "right-wall", "top"), dtype=object)
    region = region_names[np.argmin(surface_distance, axis=1)]
    coordinate = np.empty(len(centroids))
    horizontal = np.isin(region, ("floor", "top"))
    coordinate[horizontal] = centroids[horizontal, 0]
    side = ~horizontal
    region[side & (centroids[:, 0] < 0.5)] = "left-wall"
    region[side & (centroids[:, 0] >= 0.5)] = "right-wall"
    coordinate[side] = centroids[side, 2]
    bin_index = np.floor((coordinate + 1e-12) / float(width)).astype(int)
    keys = [(str(name), int(index)) for name, index in zip(region, bin_index)]
    unique = sorted(set(keys))
    lookup = {key: index for index, key in enumerate(unique)}
    labels = np.asarray([lookup[key] for key in keys], dtype=int)
    names = [f"{name}:{index * width:.3f}-{(index + 1) * width:.3f}um"
             for name, index in unique]
    return labels, names


def _metric_row(dx, basis, metrics):
    return dict(
        dx_um=float(dx), basis=basis, active_count=metrics.active_count,
        rms_relative_imbalance=metrics.rms_relative_imbalance,
        maximum_relative_imbalance=metrics.maximum_relative_imbalance,
        throughput_weighted_rms_relative_imbalance=(
            metrics.throughput_weighted_rms_relative_imbalance),
        global_relative_imbalance=metrics.global_relative_imbalance)


def _pilot_method_map(
        args, dx, reference_potential, reference_dx, raw_dir, run_tag):
    geometry, poisson, verts, faces, centroids, areas, normals = _geometry_and_poisson(dx)
    potential = _interpolate_reference_potential(
        reference_potential, reference_dx, poisson.shape, dx)
    charge, _ = _charge_reproducing_potential(poisson, potential)
    path = raw_dir / f"method_map_dx_{dx:.6f}_{run_tag}.npz"
    if path.exists():
        return np.load(path)["method_hint_Ar"]
    boundary = _boundary(2.0 * geometry.mesh_length_unit_m)
    proposals = {
        "Ar+": _ion_proposal(boundary, args.ion_proposal_level, args.pilot_ion_seed),
        "electron": _electron_proposal(
            boundary, args.electron_proposal_level, args.pilot_electron_seed),
    }
    print(f"dx={dx:.3f} um certifying independent ion method map", flush=True)
    try:
        result = integrate_dielectric_charging_transient_3d(
            poisson_system=poisson, initial_charge_node_c=charge,
            boundary=boundary, verts=verts, faces=faces, areas=areas,
            source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=2.0,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            mesh_length_unit_m=geometry.mesh_length_unit_m,
            n_position=256, seed=args.pilot_forward_seed,
            trajectory_fixed_dt=args.trajectory_dt,
            trajectory_max_steps=args.trajectory_max_steps,
            phase_space_log2_samples=args.pilot_base_level,
            periodic_lateral=True,
            transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
            face_centroids=centroids, face_gas_normals=normals,
            adjoint_face_quadrature_points=args.face_points,
            adjoint_ray_offset=1e-4, adjoint_proposals=proposals,
            adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
            bidirectional_options=dict(
                forward_log2_samples=args.pilot_base_level,
                adjoint_log2_samples=args.pilot_base_level,
                n_replicates=4,
                max_forward_log2_samples=args.pilot_max_level,
                max_adjoint_log2_samples=args.pilot_max_level,
                max_face_quadrature_points=args.pilot_max_face_points,
                element_absolute_tolerance=args.pilot_absolute_tolerance,
                element_relative_tolerance=args.pilot_relative_tolerance,
                face_quadrature_points=args.face_points),
            transport_device="cpu", timestep_s=1e-9, n_steps=0,
            current_balance_tol=0.08)
    except BidirectionalCurrentCertificationError as error:
        diagnostic = dict(error=str(error), dx_um=dx, faces=len(faces))
        if error.result is not None:
            selection = error.result.selection_by_species["Ar+"]
            diagnostic.update(
                inconsistent_faces=np.where(~selection.estimator_consistent)[0].tolist(),
                imprecise_faces=np.where(~selection.method_within_tolerance)[0].tolist())
        (raw_dir / f"method_map_failure_dx_{dx:.6f}_{run_tag}.json").write_text(
            json.dumps(diagnostic, indent=2) + "\n")
        raise
    method = np.asarray(result.bidirectional_method_hint["Ar+"])
    sampling = result.bidirectional_sampling_provenance["Ar+"]
    np.savez_compressed(
        path, method_hint_Ar=method,
        forward_log2_samples_Ar=sampling.forward_log2_samples,
        adjoint_log2_samples_by_face_Ar=sampling.adjoint_log2_samples_by_face,
        face_quadrature_points_by_face_Ar=sampling.face_quadrature_points_by_face,
        replicate_seeds_Ar=sampling.replicate_seeds)
    return method


def _run_replicate(
        args, dx, replicate, reference_potential, reference_dx, raw_dir,
        method_hint, run_tag):
    geometry, poisson, verts, faces, centroids, areas, normals = _geometry_and_poisson(dx)
    potential = _interpolate_reference_potential(
        reference_potential, reference_dx, poisson.shape, dx)
    charge, recovery_error = _charge_reproducing_potential(poisson, potential)
    raw_path = raw_dir / f"dx_{dx:.6f}_replicate_{replicate:02d}_{run_tag}.npz"
    if raw_path.exists():
        archived = np.load(raw_path)
        return dict(
            positive_face=archived["positive_face_current_density_a_m2"],
            negative_face=archived["negative_face_current_density_a_m2"],
            positive_node=archived["positive_current_node_a"],
            negative_node=archived["negative_current_node_a"],
            centroids=archived["centroids"], areas=archived["areas"],
            normals=archived["normals"], potential=archived["potential_v"],
            recovery_error=float(archived["potential_recovery_relative_l2"]),
            face_count=int(len(archived["areas"])))
    boundary = _boundary(2.0 * geometry.mesh_length_unit_m)
    proposals = {
        "Ar+": _ion_proposal(
            boundary, args.ion_proposal_level, args.ion_seed_base + replicate),
        "electron": _electron_proposal(
            boundary, args.electron_proposal_level, args.electron_seed_base + replicate),
    }
    print(
        f"dx={dx:.3f} um replicate={replicate + 1}/{args.replicates} "
        f"faces={len(faces)}", flush=True)
    result = integrate_dielectric_charging_transient_3d(
        poisson_system=poisson, initial_charge_node_c=charge,
        boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=2.0,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        n_position=256, seed=79, trajectory_fixed_dt=args.trajectory_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        periodic_lateral=True,
        phase_space_log2_samples=args.forward_level,
        transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        adjoint_face_quadrature_points=args.face_points,
        adjoint_ray_offset=1e-4, adjoint_proposals=proposals,
        adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
        bidirectional_options=dict(
            forward_log2_samples=args.forward_level,
            adjoint_log2_samples=args.ion_proposal_level,
            n_replicates=4,
            method_hint={"Ar+": np.asarray(method_hint)},
            require_certification=False,
            face_quadrature_points=args.face_points),
        transport_device="cpu", timestep_s=1e-9, n_steps=0,
        current_balance_tol=0.08)
    np.savez_compressed(
        raw_path,
        positive_face_current_density_a_m2=result.positive_face_current_density_a_m2,
        negative_face_current_density_a_m2=result.negative_face_current_density_a_m2,
        positive_current_node_a=result.positive_current_node_a,
        negative_current_node_a=result.negative_current_node_a,
        centroids=centroids, areas=areas, normals=normals,
        potential_v=result.potential_v,
        potential_recovery_relative_l2=np.asarray(recovery_error))
    return dict(
        positive_face=result.positive_face_current_density_a_m2,
        negative_face=result.negative_face_current_density_a_m2,
        positive_node=result.positive_current_node_a,
        negative_node=result.negative_current_node_a,
        centroids=centroids, areas=areas, normals=normals,
        potential=result.potential_v, recovery_error=recovery_error,
        face_count=int(len(faces)))


def _summarize_resolution(dx, replicates, patch_widths):
    positive_face_density = np.stack([item["positive_face"] for item in replicates])
    negative_face_density = np.stack([item["negative_face"] for item in replicates])
    positive_node = np.stack([item["positive_node"] for item in replicates])
    negative_node = np.stack([item["negative_node"] for item in replicates])
    areas_m2 = np.asarray(replicates[0]["areas"]) * 1e-12
    positive_face_a = positive_face_density * areas_m2[None, :]
    negative_face_a = negative_face_density * areas_m2[None, :]
    mean_positive_face = np.mean(positive_face_a, axis=0)
    mean_negative_face = np.mean(negative_face_a, axis=0)
    mean_positive_node = np.mean(positive_node, axis=0)
    mean_negative_node = np.mean(negative_node, axis=0)
    metric_rows = [
        _metric_row(dx, "raw-face", current_balance_metrics_3d(
            mean_positive_face, mean_negative_face)),
        _metric_row(dx, "compatible-node", current_balance_metrics_3d(
            mean_positive_node, mean_negative_node)),
    ]
    patch_rows = []
    patch_by_width = {}
    for width in patch_widths:
        group, names = _patch_groups(
            replicates[0]["centroids"], replicates[0]["normals"], width)
        metrics = current_balance_metrics_3d(
            mean_positive_face, mean_negative_face, group=group)
        metric_rows.append(_metric_row(dx, f"patch-{width:g}um", metrics))
        replicate_net = []
        for positive, negative in zip(positive_face_a, negative_face_a):
            sample = current_balance_metrics_3d(positive, negative, group=group)
            replicate_net.append(sample.positive_current_a - sample.negative_current_a)
        replicate_net = np.asarray(replicate_net)
        stderr = (np.std(replicate_net, axis=0, ddof=1) / np.sqrt(len(replicates))
                  if len(replicates) > 1 else np.full(len(names), np.nan))
        throughput = metrics.positive_current_a + metrics.negative_current_a
        confidence_half_width = np.divide(
            2.0 * stderr, throughput, out=np.full_like(stderr, np.nan),
            where=throughput > 0.0)
        rows = []
        for index, name in enumerate(names):
            row = dict(
                dx_um=float(dx), patch_width_um=float(width), patch=name,
                positive_current_a=float(metrics.positive_current_a[index]),
                negative_current_a=float(metrics.negative_current_a[index]),
                signed_relative_imbalance=float(metrics.signed_relative_imbalance[index]),
                net_current_stderr_a=float(stderr[index]),
                relative_two_sigma_half_width=float(confidence_half_width[index]),
                resolved_outside_contract=bool(
                    abs(metrics.signed_relative_imbalance[index])
                    - confidence_half_width[index] > 0.08),
            )
            rows.append(row); patch_rows.append(row)
        patch_by_width[width] = rows
    return dict(
        dx=float(dx), face_count=replicates[0]["face_count"],
        potential_recovery_relative_l2=max(item["recovery_error"] for item in replicates),
        metric_rows=metric_rows, patch_rows=patch_rows,
        patch_by_width=patch_by_width,
        potential=replicates[0]["potential"])


def _grid_comparison(resolutions, width):
    coarse, refined = resolutions[0], resolutions[-1]
    coarse_rows = {item["patch"]: item for item in coarse["patch_by_width"][width]}
    refined_rows = {item["patch"]: item for item in refined["patch_by_width"][width]}
    common = sorted(set(coarse_rows) & set(refined_rows))
    difference = np.asarray([
        refined_rows[name]["signed_relative_imbalance"]
        - coarse_rows[name]["signed_relative_imbalance"] for name in common])
    persistent = [name for name in common
                  if coarse_rows[name]["resolved_outside_contract"]
                  and refined_rows[name]["resolved_outside_contract"]
                  and np.sign(coarse_rows[name]["signed_relative_imbalance"])
                  == np.sign(refined_rows[name]["signed_relative_imbalance"])]
    return dict(
        patch_width_um=float(width), common_patch_count=len(common),
        signed_relative_rms_difference=float(np.sqrt(np.mean(difference ** 2))),
        signed_relative_maximum_difference=float(np.max(np.abs(difference))),
        coarse_maximum=max(abs(item["signed_relative_imbalance"])
                           for item in coarse_rows.values()),
        refined_maximum=max(abs(item["signed_relative_imbalance"])
                            for item in refined_rows.values()),
        refined_resolved_outside_contract_count=sum(
            item["resolved_outside_contract"] for item in refined_rows.values()),
        persistent_resolved_same_sign_count=len(persistent),
        persistent_resolved_same_sign_patches=persistent)


def _plot(resolutions, metric_rows, patch_width, output):
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.5), constrained_layout=True)
    bases = ["raw-face", "compatible-node", "patch-0.25um", "patch-0.5um", "patch-1um"]
    x = np.arange(len(bases)); width = 0.34
    for offset, resolution in zip((-0.5, 0.5), resolutions):
        lookup = {row["basis"]: row for row in resolution["metric_rows"]}
        axes[0, 0].bar(
            x + offset * width,
            [lookup[name]["maximum_relative_imbalance"] for name in bases],
            width=width, label=f"dx={resolution['dx']:.3f} um")
        axes[0, 1].bar(
            x + offset * width,
            [lookup[name]["rms_relative_imbalance"] for name in bases],
            width=width, label=f"dx={resolution['dx']:.3f} um")
    for axis, title in ((axes[0, 0], "Worst local/patch imbalance"),
                        (axes[0, 1], "RMS local/patch imbalance")):
        axis.axhline(0.08, color="0.25", linestyle="--", label="contract 0.08")
        axis.set_xticks(x, ["faces", "nodes", "0.25 um", "0.5 um", "1.0 um"], rotation=20)
        axis.set(ylabel="relative imbalance", title=title)
        axis.grid(axis="y", alpha=0.25); axis.legend(fontsize=8)

    all_names = sorted(set().union(*(
        {item["patch"] for item in resolution["patch_by_width"][patch_width]}
        for resolution in resolutions)))
    px = np.arange(len(all_names))
    for resolution, marker in zip(resolutions, ("o", "s")):
        lookup = {item["patch"]: item
                  for item in resolution["patch_by_width"][patch_width]}
        value = [lookup[name]["signed_relative_imbalance"] if name in lookup else np.nan
                 for name in all_names]
        axes[1, 0].plot(px, value, marker=marker, label=f"dx={resolution['dx']:.3f} um")
    axes[1, 0].axhline(0.08, color="0.25", linestyle="--")
    axes[1, 0].axhline(-0.08, color="0.25", linestyle="--")
    axes[1, 0].set_xticks(px, all_names, rotation=65, ha="right", fontsize=7)
    axes[1, 0].set(
        title=f"Signed balance in fixed {patch_width:g} um patches",
        ylabel="(ion - electron) / total")
    axes[1, 0].grid(axis="y", alpha=0.25); axes[1, 0].legend(fontsize=8)

    for resolution in resolutions:
        middle_y = resolution["potential"].shape[1] // 2
        image = resolution["potential"][:, middle_y, :].T
        axes[1, 1].plot(
            np.linspace(0.0, 1.0, image.shape[1]), image[np.argmin(
                np.abs(np.linspace(0.0, 2.0, image.shape[0]) - 1.0))],
            label=f"z=1 um, dx={resolution['dx']:.3f} um")
    axes[1, 1].set(
        title="Same interpolated voltage state", xlabel="x (um)", ylabel="potential (V)")
    axes[1, 1].grid(alpha=0.25); axes[1, 1].legend(fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--checkpoint", type=Path,
        default=ROOT / "results/charging_task1_3d_cont6_strict16/final_states.npz")
    parser.add_argument("--checkpoint-key", default="refined_potential_v")
    parser.add_argument("--reference-dx", type=float, default=0.25)
    parser.add_argument("--grid-spacings", default="0.25,0.125")
    parser.add_argument("--patch-widths", default="0.25,0.5,1.0")
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--forward-level", type=int, default=9)
    parser.add_argument("--ion-proposal-level", type=int, default=9)
    parser.add_argument("--electron-proposal-level", type=int, default=9)
    parser.add_argument("--ion-seed-base", type=int, default=1701)
    parser.add_argument("--electron-seed-base", type=int, default=2701)
    parser.add_argument("--pilot-forward-seed", type=int, default=79)
    parser.add_argument("--pilot-ion-seed", type=int, default=83)
    parser.add_argument("--pilot-electron-seed", type=int, default=89)
    parser.add_argument("--pilot-base-level", type=int, default=8)
    parser.add_argument("--pilot-max-level", type=int, default=13)
    parser.add_argument("--pilot-max-face-points", type=int, default=15)
    parser.add_argument("--pilot-absolute-tolerance", type=float, default=0.02)
    parser.add_argument("--pilot-relative-tolerance", type=float, default=0.1)
    parser.add_argument("--trajectory-dt", type=float, default=0.005)
    parser.add_argument("--trajectory-max-steps", type=int, default=50000)
    parser.add_argument("--face-points", type=int, choices=(1, 3, 7), default=3)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "results/charging_task1_equilibrium_grid_audit")
    args = parser.parse_args()
    spacings = tuple(float(item) for item in args.grid_spacings.split(","))
    patch_widths = tuple(float(item) for item in args.patch_widths.split(","))
    if (len(spacings) < 2 or any(value <= 0.0 for value in spacings)
            or any(value <= 0.0 for value in patch_widths)
            or args.replicates <= 0):
        parser.error("invalid grid, patch, or replicate controls")
    checkpoint_bytes = args.checkpoint.read_bytes()
    checkpoint_sha256 = hashlib.sha256(checkpoint_bytes).hexdigest()
    with np.load(args.checkpoint) as archived:
        reference_potential = np.asarray(archived[args.checkpoint_key], dtype=float).copy()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    config = dict(
        model="exact_hard_visibility_equilibrium_grid_audit",
        checkpoint_sha256=checkpoint_sha256, checkpoint_name=args.checkpoint.name,
        checkpoint_key=args.checkpoint_key, reference_dx_um=args.reference_dx,
        grid_spacings_um=spacings, patch_widths_um=patch_widths,
        replicates=args.replicates,
        forward_level=args.forward_level,
        ion_proposal_level=args.ion_proposal_level,
        electron_proposal_level=args.electron_proposal_level,
        ion_seeds=[args.ion_seed_base + index for index in range(args.replicates)],
        electron_seeds=[args.electron_seed_base + index for index in range(args.replicates)],
        pilot=dict(
            forward_seed=args.pilot_forward_seed, ion_seed=args.pilot_ion_seed,
            electron_seed=args.pilot_electron_seed, base_level=args.pilot_base_level,
            max_level=args.pilot_max_level, max_face_points=args.pilot_max_face_points,
            absolute_tolerance=args.pilot_absolute_tolerance,
            relative_tolerance=args.pilot_relative_tolerance),
        trajectory_dt=args.trajectory_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        face_points=args.face_points, geometry=GEOMETRY,
        estimator_map={"Ar+": "frozen_bidirectional_pilot", "electron": "adjoint"})
    config_hash = _hash_json(config)
    run_tag = config_hash[:12]
    resolutions = []
    for dx in spacings:
        method_hint = _pilot_method_map(
            args, dx, reference_potential, args.reference_dx, raw_dir, run_tag)
        samples = [_run_replicate(
            args, dx, replicate, reference_potential, args.reference_dx, raw_dir,
            method_hint, run_tag)
            for replicate in range(args.replicates)]
        resolution = _summarize_resolution(dx, samples, patch_widths)
        resolution["ion_method_counts"] = {
            name: int(np.sum(method_hint == name)) for name in ("forward", "adjoint")}
        resolutions.append(resolution)
    metric_rows = [row for result in resolutions for row in result["metric_rows"]]
    patch_rows = [row for result in resolutions for row in result["patch_rows"]]
    comparisons = [_grid_comparison(resolutions, width) for width in patch_widths]
    physically_coarsened = [
        item for item in comparisons if item["patch_width_um"] >= 2.0 * spacings[0]]
    if any(item["persistent_resolved_same_sign_count"] > 0
           for item in physically_coarsened):
        conclusion = (
            "resolved above-contract wall imbalance persists with the same sign on fixed physical "
            "patches spanning at least two coarse cells; at this mapped state the floor is not "
            "solely a raw-face or nodal discretization artifact")
    elif all(item["refined_maximum"] <= 0.08 for item in physically_coarsened):
        conclusion = (
            "fixed-size refined patches satisfy the contract while raw elements do not; "
            "promote a discretization-scale contract review with further refinement")
    else:
        conclusion = (
            "fixed-size patch balance is not grid-stable at available precision; "
            "no equilibrium/discretization attribution yet")
    summary = dict(
        config_hash=config_hash, config=config,
        exact_operator_statement="hard visibility; full kinetic currents; no smoothing",
        state_statement=(
            "same archived voltage field, trilinearly interpolated; this is an operator refinement "
            "audit, not a refined-grid steady-state claim"),
        resolution=[dict(
            dx_um=item["dx"], face_count=item["face_count"],
            ion_method_counts=item["ion_method_counts"],
            potential_recovery_relative_l2=item["potential_recovery_relative_l2"],
            metrics=item["metric_rows"]) for item in resolutions],
        fixed_patch_grid_comparisons=comparisons,
        conclusion=conclusion,
        contract_changed=False)
    encoded_config = json.dumps(config, indent=2) + "\n"
    encoded_summary = json.dumps(summary, indent=2) + "\n"
    for name in ("config.json", f"config_{run_tag}.json"):
        (args.output_dir / name).write_text(encoded_config)
    for name in ("summary.json", f"summary_{run_tag}.json"):
        (args.output_dir / name).write_text(encoded_summary)
    for name in ("balance_metrics.csv", f"balance_metrics_{run_tag}.csv"):
        with (args.output_dir / name).open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(metric_rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(metric_rows)
    for name in ("patch_balance.csv", f"patch_balance_{run_tag}.csv"):
        with (args.output_dir / name).open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(patch_rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(patch_rows)
    for name in ("equilibrium_grid_audit.png", f"equilibrium_grid_audit_{run_tag}.png"):
        _plot(resolutions, metric_rows, patch_widths[0], args.output_dir / name)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
