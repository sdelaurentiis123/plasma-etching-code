#!/usr/bin/env python3
"""Bounded, restartable charging audit on a Krüger 2024 profile checkpoint.

This is a causal development test, not a held-out validation. It uses the published aggregate ion
flux/IEAD, closes the global electron current exactly as MCFPM does, and evolves local charge with
the common fresh-scramble physical-time/Q1-Poisson engine. The paired final audit asks whether the
charged field reduces ion delivery to the feature floor relative to the identical zero-charge
geometry. It never reads a held-out profile and it never invokes a frozen-map root solver.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from petch.boundary_state import PlasmaBoundaryState
from petch.boundary_transport_3d import BoundaryTransport3DResult
from petch.charging_coevolution_3d import integrate_surface_charging_to_saturation_3d
from petch.feature_step_3d import (
    FeatureGeometry3D, _face_material_ids, _surface_gas_normals)
from petch.krueger_replay_3d import make_krueger_2024_poisson_system_3d
from petch.reactor_boundary import (
    append_global_current_balance_maxwellian_electrons,
    build_krueger_2024_development_boundary,
)
from petch.threed import extract_mesh_3d


DATA = ROOT / "data" / "experimental" / "krueger_2024"
ELECTRON_TEMPERATURE_BOUNDS_EV = (3.4, 3.8)
MASK_RELATIVE_PERMITTIVITY_BOUNDS = (2.9, 3.8)


def _sha256_path(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _jsonable(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return str(value)


def _load_profile_geometry(path):
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"]))
        layers = {
            int(name.rsplit("_", 1)[1]): np.asarray(archive[name])
            for name in archive.files if name.startswith("material_levelset_")}
        geometry = FeatureGeometry3D(
            archive["phi"], archive["material_id"], metadata["dx"],
            metadata["mesh_length_unit_m"], tuple(archive["mesh_origin_m"]),
            material_levelsets=layers)
    return geometry, metadata


def _charged_boundary(reference_plane_m, electron_temperature_eV):
    full = build_krueger_2024_development_boundary(
        DATA, reference_plane_m=reference_plane_m,
        neutral_direction_polar_order=8,
        neutral_direction_azimuthal_order=16,
        ion_energy_bin_eV=250.0,
        ion_angle_bin_deg=0.25,
        ion_azimuthal_closure="axisymmetric_uniform")
    ion = full.get("ions")
    charged_heavy = PlasmaBoundaryState(
        species=(ion,), reference_plane_m=full.reference_plane_m,
        provenance={
            **dict(full.provenance),
            "selected_for_fixed_geometry_charge_audit": ("ions",),
            "neutral_species_omitted_because_they_carry_no_current": True,
        })
    return append_global_current_balance_maxwellian_electrons(
        charged_heavy,
        electron_temperature_eV=electron_temperature_eV,
        temperature_source=(
            "Krueger PhD thesis (2024), Fig. 6.3: HPEM bulk Te=3.4-3.8 eV; "
            "Sec. 2.2.2 global net-neutral charged flux; Wang and Kushner, "
            "JAP 107, 023309 (2010), thermal Maxwellian/Lambertian closure"),
        temperature_evidence_kind="published_HPEM_output_and_MCFPM_closure",
        n_transverse=5,
        n_normal=8)


def _ion_flux_by_face(transport: BoundaryTransport3DResult, face_count):
    value = np.zeros(int(face_count), dtype=float)
    for population in transport.surface_fluxes.energetic_fluxes:
        if population.name == "ions":
            value += np.asarray(population.flux_m2_s, dtype=float)
    return value


def _floor_flux(transport, centroids, normals, material, areas, dx):
    candidates = (material == 1) & (normals[:, 2] > 0.5)
    if not np.any(candidates):
        raise RuntimeError("Krueger checkpoint has no resolved upward-facing SiO2 surface")
    floor_z = float(np.min(centroids[candidates, 2]))
    selected = candidates & (centroids[:, 2] <= floor_z + 1.5 * float(dx))
    flux = _ion_flux_by_face(transport, len(areas))
    return (
        float(np.dot(flux[selected], areas[selected]) / np.sum(areas[selected])),
        floor_z,
        int(np.count_nonzero(selected)),
    )


def _atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _patch_balance_record(item):
    """Persist patch metrics without duplicating the full face-to-patch map in every audit."""
    group = np.asarray(item.group, dtype=np.int64)
    result = {
        key: value for key, value in item.__dict__.items() if key != "group"
    }
    result.update(
        group_count=int(np.unique(group).size),
        face_count=int(group.size),
        group_map_sha256=sha256(group.astype("<i8", copy=False).tobytes()).hexdigest(),
    )
    return result


def _load_charge_checkpoint(path, *, profile_sha256, face_count):
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"]))
        sigma = np.asarray(archive["sigma_c_per_m2"], dtype=float)
    if (
        metadata.get("profile_checkpoint_sha256") != profile_sha256
        or sigma.shape != (int(face_count),)
        or np.any(~np.isfinite(sigma))
    ):
        raise ValueError("charge checkpoint does not match the profile surface")
    return sigma, metadata


def _atomic_charge_checkpoint(path, sigma, metadata):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(
            stream,
            sigma_c_per_m2=np.asarray(sigma, dtype=float),
            metadata_json=np.asarray(json.dumps(_jsonable(metadata), sort_keys=True)))
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--charge-checkpoint", type=Path)
    parser.add_argument("--resume-charge", action="store_true")
    parser.add_argument("--steps", type=int, default=32)
    parser.add_argument("--timestep-us", type=float, default=1.0)
    parser.add_argument("--phase-space-log2-samples", type=int, default=9)
    parser.add_argument("--n-position", type=int, default=32)
    parser.add_argument("--seed", type=int, default=4241)
    parser.add_argument("--audit-epoch-offset", type=int, default=100000)
    parser.add_argument("--electron-temperature-eV", type=float, default=3.6)
    parser.add_argument("--mask-relative-permittivity", type=float, default=3.3)
    parser.add_argument("--trajectory-fixed-dt", type=float, default=0.001)
    parser.add_argument("--trajectory-max-steps", type=int, default=4096)
    parser.add_argument("--trajectory-emergency-max-steps", type=int, default=65536)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument(
        "--compatible-q1-charge-state", action="store_true",
        help=(
            "Canonicalize P0 face charge into the minimum-norm Q1-visible representative. "
            "Off by default: the raw conservative face ledger produces the identical Q1 load "
            "without the optional dense rank-revealing diagnostic."))
    args = parser.parse_args()

    if (
        int(args.steps) != args.steps or args.steps < 0
        or not np.isfinite(args.timestep_us) or args.timestep_us <= 0.0
        or not ELECTRON_TEMPERATURE_BOUNDS_EV[0]
        <= args.electron_temperature_eV
        <= ELECTRON_TEMPERATURE_BOUNDS_EV[1]
        or not MASK_RELATIVE_PERMITTIVITY_BOUNDS[0]
        <= args.mask_relative_permittivity
        <= MASK_RELATIVE_PERMITTIVITY_BOUNDS[1]
    ):
        raise ValueError("charging controls fall outside their declared source bounds")

    profile_sha = _sha256_path(args.profile_checkpoint)
    geometry, profile_metadata = _load_profile_geometry(args.profile_checkpoint)
    x, y, z = geometry.coordinate_arrays
    source_z = float(z[-1])
    source_bounds = (float(x[0]), float(x[-1]), float(y[0]), float(y[-1]))
    boundary = _charged_boundary(
        geometry.mesh_origin_m[2] + source_z * geometry.mesh_length_unit_m,
        args.electron_temperature_eV)
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    material = _face_material_ids(centroids, geometry)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    poisson = make_krueger_2024_poisson_system_3d(
        geometry, mask_relative_permittivity=args.mask_relative_permittivity)

    charge_path = (
        args.output.with_name(args.output.stem + "_charge_state.npz")
        if args.charge_checkpoint is None else args.charge_checkpoint)
    initial_sigma = np.zeros(len(faces), dtype=float)
    initial_epoch = 0
    prior_physical_time_s = 0.0
    prior_steps = 0
    if args.resume_charge:
        initial_sigma, prior = _load_charge_checkpoint(
            charge_path, profile_sha256=profile_sha, face_count=len(faces))
        expected_restart_controls = {
            "electron_temperature_eV": args.electron_temperature_eV,
            "mask_relative_permittivity": args.mask_relative_permittivity,
            "compatible_q1_charge_state": bool(args.compatible_q1_charge_state),
        }
        mismatched = {
            key: (prior.get(key), expected)
            for key, expected in expected_restart_controls.items()
            if prior.get(key) != expected
        }
        if mismatched:
            raise ValueError(
                f"charge checkpoint controls do not match the resumed operator: {mismatched}")
        initial_epoch = int(prior["resume_sampling_epoch"])
        prior_physical_time_s = float(prior["cumulative_physical_time_s"])
        prior_steps = int(prior["cumulative_accepted_steps"])

    common = dict(
        poisson_system=poisson,
        boundary=boundary,
        verts=verts,
        faces=faces,
        areas=areas,
        face_centroids=centroids,
        face_gas_normals=normals,
        face_material_id=material,
        source_bounds=source_bounds,
        source_z=source_z,
        potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=geometry.dx,
        patch_scales_m=(2.0 * geometry.dx * geometry.mesh_length_unit_m, 90.0e-9),
        potential_rate_tolerance_v_s=1.0e3,
        timestep_s=args.timestep_us * 1.0e-6,
        current_balance_tolerance=0.08,
        timestep_policy="fixed",
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        mesh_origin_m=geometry.mesh_origin_m,
        n_position=args.n_position,
        seed=args.seed,
        trajectory_fixed_dt=args.trajectory_fixed_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        trajectory_adaptive_horizon=True,
        trajectory_emergency_max_steps=args.trajectory_emergency_max_steps,
        phase_space_log2_samples=args.phase_space_log2_samples,
        periodic_lateral=True,
        transport_estimator="forward",
        transport_device=args.transport_device,
        stop_on_saturation=False,
        scramble_mode="fresh",
        compatible_q1_charge_state=bool(args.compatible_q1_charge_state),
    )

    started = perf_counter()

    def progress(*, potential_v, history_item, accepted_steps, physical_time_s, **_):
        evaluated = int(accepted_steps) + 1
        if evaluated == 1 or evaluated % 4 == 0 or evaluated == args.steps + 1:
            print(
                f"charge eval {evaluated}/{args.steps + 1}: "
                f"accepted={accepted_steps}, segment_t={physical_time_s * 1e6:.3f} us, "
                f"max|V|={np.max(np.abs(potential_v)):.5g} V, "
                f"node_rms={history_item['rms_relative_current_imbalance_node']:.5g}",
                flush=True)

    march = integrate_surface_charging_to_saturation_3d(
        initial_sigma_c_per_m2=initial_sigma,
        maximum_steps=args.steps,
        initial_sampling_epoch=initial_epoch,
        progress_callback=progress,
        **common)
    audit_epoch = int(
        march.diagnostics["resume_sampling_epoch"] + args.audit_epoch_offset)
    zero_audit = integrate_surface_charging_to_saturation_3d(
        initial_sigma_c_per_m2=np.zeros(len(faces)),
        maximum_steps=0,
        initial_sampling_epoch=audit_epoch,
        **common)
    charged_audit = integrate_surface_charging_to_saturation_3d(
        initial_sigma_c_per_m2=march.sigma_c_per_m2,
        maximum_steps=0,
        initial_sampling_epoch=audit_epoch,
        **common)

    zero_floor, floor_z, floor_face_count = _floor_flux(
        zero_audit.final_step.transport, centroids, normals, material, areas, geometry.dx)
    charged_floor, charged_floor_z, charged_floor_face_count = _floor_flux(
        charged_audit.final_step.transport, centroids, normals, material, areas, geometry.dx)
    if floor_z != charged_floor_z or floor_face_count != charged_floor_face_count:
        raise RuntimeError("paired floor audit changed its geometric measurement set")

    cumulative_time = prior_physical_time_s + march.physical_time_s
    cumulative_steps = prior_steps + march.accepted_steps
    checkpoint_metadata = {
        "schema": "petch.krueger-2024.fixed-geometry-charge-checkpoint.v1",
        "profile_checkpoint_sha256": profile_sha,
        "resume_sampling_epoch": int(march.diagnostics["resume_sampling_epoch"]),
        "cumulative_physical_time_s": cumulative_time,
        "cumulative_accepted_steps": cumulative_steps,
        "electron_temperature_eV": args.electron_temperature_eV,
        "mask_relative_permittivity": args.mask_relative_permittivity,
        "phase_space_log2_samples": args.phase_space_log2_samples,
        "seed": args.seed,
        "compatible_q1_charge_state": bool(args.compatible_q1_charge_state),
    }
    _atomic_charge_checkpoint(
        charge_path, march.sigma_c_per_m2, checkpoint_metadata)

    output = {
        "schema": "petch.krueger-2024.charging-causality.v1",
        "scientific_status": (
            "bounded fixed-geometry causal development audit; not calibration reveal or "
            "held-out experimental validation"),
        "profile_checkpoint": {
            "path_name": args.profile_checkpoint.name,
            "sha256": profile_sha,
            "profile_metadata": profile_metadata,
        },
        "source_physics": {
            "global_current_closure": (
                "electron flux balances the charge-weighted positive-ion flux globally"),
            "local_current_closure": "none; resolved by kinetic trajectories and Q1 Poisson",
            "electron_temperature_eV": args.electron_temperature_eV,
            "electron_temperature_bounds_eV": ELECTRON_TEMPERATURE_BOUNDS_EV,
            "electron_temperature_source": (
                "Krueger 2024 thesis Fig. 6.3, HPEM bulk Te=3.4-3.8 eV"),
            "electron_distribution": "flux half-Maxwellian with Lambertian angular marginal",
            "electron_flux_method_source": (
                "Krueger 2024 thesis Sec. 2.2.2; Wang and Kushner, JAP 107, 023309 (2010)"),
            "unpublished_source_boundary": (
                "the chapter-6 wafer EEAD/high-energy-electron fraction is not published; "
                "thermal-only electrons are a declared sensitivity endpoint"),
            "volume_boltzmann_electron_term": False,
            "secondary_electron_emission": "disabled; not declared for the chapter-6 replay",
            "surface_conductivity": "zero-mobility insulating endpoint; source value unpublished",
        },
        "electrostatics": {
            "sio2_relative_permittivity": 3.9,
            "mask_relative_permittivity": args.mask_relative_permittivity,
            "mask_relative_permittivity_bounds": MASK_RELATIVE_PERMITTIVITY_BOUNDS,
            "mask_permittivity_source": (
                "a-C:H literature range, e.g. DOI 10.1143/JJAP.42.259; central value is "
                "a declared development closure, not a Krueger measurement"),
            "lateral_boundary": "periodic",
            "top_boundary": "natural_zero_normal_displacement",
            "bottom_boundary": "grounded",
        },
        "integration": {
            "segment_steps": args.steps,
            "segment_timestep_us": args.timestep_us,
            "segment_physical_time_us": march.physical_time_s * 1e6,
            "cumulative_physical_time_us": cumulative_time * 1e6,
            "cumulative_accepted_steps": cumulative_steps,
            "scramble_mode": "fresh",
            "phase_space_log2_samples": args.phase_space_log2_samples,
            "n_position": args.n_position,
            "trajectory_fixed_dt": args.trajectory_fixed_dt,
            "trajectory_max_steps": args.trajectory_max_steps,
            "trajectory_emergency_max_steps": args.trajectory_emergency_max_steps,
            "transport_device": args.transport_device,
            "paired_unused_audit_epoch": audit_epoch,
            "compatible_q1_charge_state": bool(args.compatible_q1_charge_state),
            "surface_charge_state_authority": (
                "minimum_norm_Q1_visible_face_representative"
                if args.compatible_q1_charge_state
                else "raw_conservative_P0_face_deposition_ledger"),
            "wall_time_s": perf_counter() - started,
            "restartable_charge_checkpoint_name": charge_path.name,
        },
        "paired_exact_hard_visibility_audit": {
            "floor_z_um": floor_z,
            "floor_face_count": floor_face_count,
            "zero_charge_floor_ion_flux_m2_s": zero_floor,
            "charged_floor_ion_flux_m2_s": charged_floor,
            "charged_over_zero_floor_ion_flux": charged_floor / zero_floor,
            "maximum_absolute_potential_v": float(
                np.max(np.abs(charged_audit.potential_v))),
            "minimum_potential_v": float(np.min(charged_audit.potential_v)),
            "maximum_potential_v": float(np.max(charged_audit.potential_v)),
            "retained_node_rms_relative_current_imbalance": charged_audit.diagnostics[
                "retained_node_rms_relative_current_imbalance"],
            "retained_node_max_relative_current_imbalance": charged_audit.diagnostics[
                "retained_node_max_relative_current_imbalance"],
            "patch_balance": [
                _patch_balance_record(item) for item in charged_audit.patch_balance],
            "signed_r2_converged": charged_audit.converged,
        },
        "conservation": {
            "maximum_charge_conservation_relative_error": max(
                (item.get("charge_conservation_relative_error", 0.0)
                 for item in march.history),
                default=0.0),
            "maximum_face_to_node_update_relative_error": max(
                (item.get("face_to_node_update_relative_error", 0.0)
                 for item in march.history),
                default=0.0),
            "maximum_lineage_replay_fraction": max(
                (item.get("transport_lineage_replay_fraction", 0.0)
                 for item in march.history),
                default=0.0),
        },
        "history": march.history,
        "held_out_profile_data_read": False,
    }
    _atomic_json(args.output, output)
    print(json.dumps(_jsonable({
        "cumulative_physical_time_us": cumulative_time * 1e6,
        **output["paired_exact_hard_visibility_audit"],
    }), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
