#!/usr/bin/env python3
"""Bounded, calibration-only Krueger short multiresolution audit.

This audit locates spatial-discretization error without opening any held-out
profile observations and without authorizing a long refinement campaign.  It
compares the same published *base* boundary, chemistry parameters, random seed
epochs, and physical timestep schedule at two representative states:

* the initial open feature, advanced for a short physical interval;
* an archived late 5 nm checkpoint, first restricted to an aligned grid with a
  material-local conservative surface-state remap, then advanced briefly.

The finest 5 nm case owns the timestep schedule.  Coarser cases replay that
schedule exactly.  Every case has a process-level wall timeout and an accepted
step ceiling.  A nominal 20 nm case is intentionally refused for the published
130 x 20 nm periodic cell: an isotropic 20 nm nodal grid cannot represent both
periods, so running it would change the physical domain rather than refine the
same problem.

Outputs are development evidence only.  They may diagnose rate, transport,
geometry, or remap sensitivity, but they do not calibrate parameters and do not
read Krueger held-out outcomes.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import krueger_2024_trench_pilot as pilot
from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.feature_step_3d import (
    FeatureGeometry3D,
    SurfaceTopologyChangeError,
    _face_material_ids,
    _periodic_physical_volume_topology_signature,
    _surface_mesh_fingerprint,
    advance_feature_step_3d,
    conservative_remap_surface_state,
    make_rectangular_trench_geometry_3d,
)
from petch.reactor_boundary import build_krueger_2024_development_boundary
from petch.threed import extract_mesh_3d


DATA = ROOT / "data" / "experimental" / "krueger_2024"
DEFAULT_SOURCE = (
    ROOT / "results" / "krueger_2024_base_calibration_r17"
    / "fine_fixed_pair_development" / "checkpoint.npz"
)
CELL_EXTENT_UM = np.asarray((0.13, 0.02, 2.8), dtype=float)
ETCHABLE = (1, 2)
CALIBRATION = {
    "effective_mask_crosslinked_growth_fraction": 0.8934059741411972,
    "oxide_etch_yield_scale": 0.5667632723491973,
}
OPERATOR = {
    "boundary_case": "base_calibration_only",
    "n_position": 16,
    "neutral_speed_quadrature": "analytic_speed_marginal",
    "neutral_tensor_velocity_quadrature_active": False,
    "neutral_direction_polar_order": 8,
    "neutral_direction_azimuthal_order": 16,
    "ion_energy_bin_eV": 250.0,
    "ion_angle_bin_deg": 0.25,
    "ion_azimuthal_order": 16,
    "ballistic_transport": "face_gather",
    "ballistic_face_quadrature_points": 3,
    "radiosity_rays_per_face": 8,
    "radiosity_relative_tolerance": 1e-12,
    "radiosity_maximum_iterations": 2000,
    "profile_reinitialization": "cr2",
    "profile_periodic_lateral": True,
    "charging": "disabled_for_Krueger_2024_calibration_and_transfer",
}
REMAP_BACKENDS = (
    "legacy_knn", "indexed_knn", "partitioned_overlap", "common_refinement")


def operator_contract(surface_state_remap_backend):
    backend = str(surface_state_remap_backend)
    if backend not in REMAP_BACKENDS:
        raise ValueError("unknown surface-state remap backend")
    return dict(OPERATOR, surface_state_remap_backend=backend)


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


def _write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _sha256(path):
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_source_provenance(checkpoint, metadata):
    """Verify that a late checkpoint came from the declared base calibration operator."""
    checkpoint = Path(checkpoint)
    audit_path = checkpoint.with_name("audit.json")
    if not audit_path.is_file():
        raise ValueError("late checkpoint requires its sibling audit.json provenance")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    configuration = dict(audit.get("configuration", {}))
    required = {
        "boundary_case": "base",
        "oxygen_to_fluorocarbon_ratio": None,
        "low_frequency_power_kw": None,
        "effective_mask_crosslinked_growth_fraction": (
            CALIBRATION["effective_mask_crosslinked_growth_fraction"]
        ),
        "oxide_etch_yield_scale": CALIBRATION["oxide_etch_yield_scale"],
        "dx_um": float(metadata["dx"]),
        "seed": 241,
    }
    mismatch = {
        name: {"required": value, "observed": configuration.get(name)}
        for name, value in required.items()
        if configuration.get(name) != value
    }
    history = audit.get("history", ())
    if not history:
        mismatch["history"] = {"required": "nonempty", "observed": "missing"}
    else:
        audited_time = float(history[-1]["physical_time_s"])
        if not np.isclose(
                audited_time, float(metadata["physical_time_s"]), rtol=0.0, atol=1e-12):
            mismatch["physical_time_s"] = {
                "required": float(metadata["physical_time_s"]),
                "observed": audited_time,
            }
    if mismatch:
        raise ValueError(
            "late checkpoint is not the declared base-calibration source: "
            + json.dumps(mismatch, sort_keys=True)
        )
    return {
        "audit_name": audit_path.name,
        "audit_sha256": _sha256(audit_path),
        "audit_config_hash": audit.get("config_hash"),
        "audit_status": audit.get("status"),
        "verified_configuration": required,
    }


def grid_contract(dx_um, extents_um=CELL_EXTENT_UM):
    """Return whether one isotropic nodal spacing preserves every domain endpoint."""
    dx = float(dx_um)
    extents = np.asarray(extents_um, dtype=float)
    if not np.isfinite(dx) or dx <= 0.0 or extents.shape != (3,):
        raise ValueError("grid contract requires a positive spacing and three extents")
    requested_intervals = extents / dx
    nearest_intervals = np.rint(requested_intervals).astype(int)
    # FeatureGeometry3D requires at least three nodes, so the constructor realizes at least two
    # intervals even when one physical period is only one requested dx wide.
    realized_intervals = np.maximum(nearest_intervals, 2)
    endpoint_error = np.abs(realized_intervals * dx - extents)
    minimum_intervals = nearest_intervals >= 2
    compatible = bool(
        np.all(endpoint_error <= 128.0 * np.finfo(float).eps * np.maximum(extents, 1.0))
        and np.all(minimum_intervals)
    )
    realized = realized_intervals * dx
    reasons = []
    if np.any(endpoint_error > 128.0 * np.finfo(float).eps * np.maximum(extents, 1.0)):
        reasons.append("spacing does not divide every physical periodic/domain extent")
    if np.any(~minimum_intervals):
        reasons.append("at least two intervals per axis are required by the 3-D engine")
    return {
        "dx_nm": dx * 1000.0,
        "compatible": compatible,
        "requested_extent_um": extents,
        "nearest_requested_interval_count": nearest_intervals,
        "realized_interval_count": realized_intervals,
        "realized_extent_um": realized,
        "endpoint_error_nm": endpoint_error * 1000.0,
        "reason": "; ".join(reasons) if reasons else None,
    }


def build_plan(source_checkpoint, *, levels_nm=(20.0, 10.0, 5.0),
               initial_duration_s=0.5, late_duration_s=0.1, seed=241,
               surface_state_remap_backend="legacy_knn"):
    source = Path(source_checkpoint)
    levels = []
    for value in levels_nm:
        contract = grid_contract(float(value) / 1000.0)
        levels.append({
            **contract,
            "initial_case": "eligible" if contract["compatible"] else "blocked",
            "late_case": "eligible" if contract["compatible"] else "blocked",
            "scientific_use": (
                "paired spatial-refinement evidence"
                if contract["compatible"] else
                "none; running this spacing would change the physical domain"
            ),
        })
    source_record = {
        "path_name": source.name,
        "exists": source.is_file(),
        "sha256": _sha256(source) if source.is_file() else None,
    }
    if source.is_file():
        source_geometry, source_state, source_fingerprint, metadata = (
            pilot._load_checkpoint(source)
        )
        source_record.update({
            "dx_nm": float(metadata["dx"]) * 1000.0,
            "physical_time_s": float(metadata["physical_time_s"]),
            "step": int(metadata["step"]),
            "mesh_fingerprint": str(metadata["fingerprint"]),
        })
        source_record["provenance"] = _validate_source_provenance(source, metadata)
        for level in levels:
            if not level["compatible"]:
                level["late_restriction_preflight"] = {
                    "status": "blocked_incompatible_grid",
                    "reason": level["reason"],
                }
                continue
            try:
                _, _, _, restriction = restrict_checkpoint_state(
                    source_geometry,
                    source_state,
                    source_fingerprint,
                    float(level["dx_nm"]) / 1000.0,
                )
            except ValueError as error:
                level["late_case"] = "blocked"
                level["scientific_use"] = (
                    "initial-state refinement only; late restriction is not topology-equivalent"
                )
                level["late_restriction_preflight"] = {
                    "status": "blocked_non_equivalent_checkpoint",
                    "reason": str(error),
                }
            else:
                level["late_restriction_preflight"] = {
                    "status": "eligible",
                    "method": restriction["method"],
                    "topology": restriction["topology"],
                }
    return {
        "campaign": "krueger_2024_bounded_multiresolution_audit",
        "scientific_status": "calibration-only numerical diagnosis; no held-out outcomes read",
        "source_checkpoint": source_record,
        "levels": levels,
        "state_windows": {
            "initial": {"start_s": 0.0, "duration_s": float(initial_duration_s)},
            "late": {
                "start_s": source_record.get("physical_time_s"),
                "duration_s": float(late_duration_s),
            },
        },
        "pairing_contract": {
            "boundary": operator_contract(surface_state_remap_backend),
            "calibration_parameters": CALIBRATION,
            "seed": int(seed),
            "schedule_owner": "5 nm",
            "coarser_schedule": "exact replay of accepted 5 nm physical timesteps",
            "final_scoring": "same frozen operator summary at both endpoints",
        },
        "required_outputs": (
            "depth/opening/width increments and rates",
            "area-weighted velocity by material",
            "area-integrated neutral and energetic incident rates",
            "surface-state intensive means and conservative integrals",
            "material-exchange and radiosity ledgers",
            "restriction/remap provenance and conservation",
        ),
    }


def _active_mesh(geometry):
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    material = _face_material_ids(centroids, geometry)
    active = np.flatnonzero(np.isin(material, ETCHABLE))
    if not active.size:
        raise ValueError("restricted geometry contains no active material surface")
    return verts, faces, centroids, areas, material, active


def _state_summary(state, area, material, mesh_length_unit_m):
    fields = dict(state.conservative_surface_fields())
    modes = dict(state.surface_field_remap_modes())
    scale = float(mesh_length_unit_m) ** 2
    output = {}
    for material_id in sorted(set(np.asarray(material, dtype=int))):
        selected = np.asarray(material) == material_id
        local_area = np.asarray(area, dtype=float)[selected]
        record = {}
        for name, values in fields.items():
            local = np.asarray(values, dtype=float)[selected]
            if modes[name] == "conservative":
                record[name] = {
                    "mode": "conservative",
                    "area_integral_physical": float(np.dot(local, local_area) * scale),
                }
            else:
                record[name] = {
                    "mode": "intensive",
                    "area_weighted_mean": float(np.dot(local, local_area) / local_area.sum()),
                    "minimum": float(np.min(local)),
                    "maximum": float(np.max(local)),
                }
        output[str(int(material_id))] = record
    return output


def restrict_checkpoint_state(geometry, state, fingerprint, target_dx_um):
    """Restrict an aligned uniform checkpoint and conservatively remap face state."""
    target_dx = float(target_dx_um)
    if target_dx < geometry.dx:
        raise ValueError("checkpoint restriction cannot invent a finer state")
    factor_float = target_dx / geometry.dx
    factor = int(round(factor_float))
    if factor <= 0 or not np.isclose(factor_float, factor, rtol=0.0, atol=1e-12):
        raise ValueError("target spacing must be an integer multiple of source spacing")
    extents = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    contract = grid_contract(target_dx, extents)
    if not contract["compatible"]:
        raise ValueError(
            f"target {target_dx * 1000:g} nm cannot preserve checkpoint extents: "
            f"{contract['reason']}"
        )
    if factor == 1:
        return geometry, state, fingerprint, {
            "method": "identity",
            "source_dx_nm": geometry.dx * 1000.0,
            "target_dx_nm": target_dx * 1000.0,
            "topology": _periodic_physical_volume_topology_signature(geometry, ETCHABLE),
        }

    index = tuple(np.arange(0, size, factor, dtype=int) for size in geometry.phi.shape)
    if any(axis[-1] != size - 1 for axis, size in zip(index, geometry.phi.shape)):
        raise ValueError("aligned restriction did not retain every physical endpoint")
    selector = np.ix_(*index)
    if geometry.material_levelsets is None:
        raise ValueError("late checkpoint lacks material-local level sets")
    layers = {
        int(material_id): np.asarray(field)[selector]
        for material_id, field in geometry.material_levelsets.items()
    }
    material_ids = np.asarray(sorted(layers), dtype=int)
    layer_stack = np.stack([layers[int(material_id)] for material_id in material_ids])
    union = np.max(layer_stack, axis=0)
    winner = material_ids[np.argmax(layer_stack, axis=0)]
    owner = np.where(union >= 0.0, winner, 0)
    restricted = FeatureGeometry3D(
        union,
        owner,
        target_dx,
        geometry.mesh_length_unit_m,
        geometry.mesh_origin_m,
        material_levelsets=layers,
    )

    old_topology = _periodic_physical_volume_topology_signature(geometry, ETCHABLE)
    new_topology = _periodic_physical_volume_topology_signature(restricted, ETCHABLE)
    if old_topology != new_topology:
        raise ValueError(
            "aligned grid restriction changes physical topology from "
            f"{old_topology} to {new_topology}; late coarse burst is not faithful"
        )
    old = _active_mesh(geometry)
    new = _active_mesh(restricted)
    old_fp = _surface_mesh_fingerprint(old[0], old[1], old[5], old[4], geometry)
    if str(fingerprint) != old_fp:
        raise ValueError("source checkpoint fingerprint does not match reconstructed surface")
    periodic_lengths = tuple(extents[:2]) + (None,)
    remapped, diagnostics = conservative_remap_surface_state(
        state,
        old[2][old[5]], old[3][old[5]], old[4][old[5]],
        new[2][new[5]], new[3][new[5]], new[4][new[5]],
        dx=target_dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        maximum_distance=2.0 * target_dx,
        old_triangles=old[0][old[1][old[5]]],
        periodic_lengths=periodic_lengths,
    )
    new_fp = _surface_mesh_fingerprint(new[0], new[1], new[5], new[4], restricted)
    return restricted, remapped, new_fp, {
        "method": "aligned_nodal_material_levelset_restriction_then_conservative_face_remap",
        "source_dx_nm": geometry.dx * 1000.0,
        "target_dx_nm": target_dx * 1000.0,
        "restriction_factor": factor,
        "source_shape": geometry.phi.shape,
        "target_shape": restricted.phi.shape,
        "physical_extents_um": extents,
        "topology": old_topology,
        "surface_state_remap": diagnostics,
    }


def _operator(geometry, seed, device):
    domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(domain[2])
    boundary = build_krueger_2024_development_boundary(
        DATA,
        n_transverse_neutral=5,
        n_normal_neutral=8,
        reference_plane_m=source_z * geometry.mesh_length_unit_m,
        neutral_direction_polar_order=8,
        neutral_direction_azimuthal_order=16,
        ion_energy_bin_eV=250.0,
        ion_angle_bin_deg=0.25,
        ion_azimuthal_closure="axisymmetric_uniform",
        ion_azimuthal_order=16,
    )
    mechanism = build_krueger_2024_material_router_3d(**CALIBRATION)
    role = {
        species.name: (
            "energetic_bombardment" if species.charge_number != 0 else "neutral_reactant"
        )
        for species in boundary.species
    }
    radiosity = {
        "rays_per_face": 8,
        "seed": int(seed) + 10000,
        "periodic_lateral": True,
        "domain_size": domain,
        "relative_tolerance": 1e-12,
        "maximum_iterations": 2000,
    }
    return boundary, mechanism, role, radiosity, domain, source_z


def _advance(geometry, state, fingerprint, *, duration_s, seed, device,
             topology_policy, remap_backend="legacy_knn"):
    boundary, mechanism, role, radiosity, domain, source_z = _operator(
        geometry, seed, device
    )
    return advance_feature_step_3d(
        geometry,
        boundary,
        role,
        mechanism,
        etchable_material_ids=ETCHABLE,
        duration_s=float(duration_s),
        source_bounds=(0.0, domain[0], 0.0, domain[1]),
        source_z=source_z,
        surface_state=state,
        surface_state_mesh_fingerprint=fingerprint,
        n_position=16,
        seed=int(seed),
        cfl_number=0.25,
        reinitialize=True,
        reinitialization_method="cr2",
        profile_periodic_lateral=True,
        topology_change_policy=str(topology_policy),
        surface_state_remap_backend=str(remap_backend),
        transport_device=str(device),
        neutral_radiosity_options=radiosity,
        ballistic_transport="face_gather",
        ballistic_face_quadrature_points=3,
    )


def _operator_summary(result):
    area = np.asarray(result.active_face_area, dtype=float)
    active = np.asarray(result.active_face_index, dtype=int)
    material = np.asarray(result.face_material_id, dtype=int)[active]
    velocity = np.asarray(result.face_velocity_mesh_units_s, dtype=float)[active]
    physical_area = area * result.geometry.mesh_length_unit_m ** 2
    fluxes = result.transport.surface_fluxes
    neutral = {}
    for name, values in fluxes.neutral_flux_m2_s.items():
        value = np.asarray(values, dtype=float)
        neutral[name] = {
            "area_weighted_mean_m2_s": float(np.dot(value, area) / area.sum()),
            "incident_rate_s": float(np.dot(value, physical_area)),
        }
    energetic = {}
    for population in fluxes.energetic_fluxes:
        value = np.asarray(population.flux_m2_s, dtype=float)
        energetic[population.name] = {
            "area_weighted_mean_m2_s": float(np.dot(value, area) / area.sum()),
            "incident_rate_s": float(np.dot(value, physical_area)),
        }
    by_material = {}
    for material_id in sorted(set(material)):
        selected = material == material_id
        local_area = area[selected]
        local_velocity = velocity[selected] * result.geometry.mesh_length_unit_m
        by_material[str(int(material_id))] = {
            "surface_area_m2": float(local_area.sum() * result.geometry.mesh_length_unit_m ** 2),
            "area_weighted_signed_velocity_m_s": float(
                np.dot(local_velocity, local_area) / local_area.sum()
            ),
            "area_weighted_absolute_velocity_m_s": float(
                np.dot(np.abs(local_velocity), local_area) / local_area.sum()
            ),
            "maximum_absolute_velocity_m_s": float(np.max(np.abs(local_velocity))),
        }
    radiosity = result.diagnostics["neutral_radiosity"]
    return {
        "active_face_count": int(area.size),
        "maximum_velocity_m_s": float(result.diagnostics["max_velocity_m_s"]),
        "raw_maximum_face_velocity_m_s": float(
            result.diagnostics["raw_maximum_face_velocity_m_s"]
        ),
        "velocity_by_material": by_material,
        "neutral_flux": neutral,
        "energetic_flux": energetic,
        "hit_probability": dict(result.transport.hit_probability),
        "state": _state_summary(
            result.surface.state, area, material, result.geometry.mesh_length_unit_m
        ),
        "maximum_material_ledger_residual_units_m2": pilot._maximum_ledger_residual(
            result.surface.material_exchange
        ),
        "maximum_radiosity_relative_balance_error": max(
            (float(item["relative_balance_error"]) for item in radiosity.values()),
            default=0.0,
        ),
    }


def _complete_analytic_initial_metrics(metrics):
    """Fill observables undefined before the substrate has moved one resolved row.

    At t=0 the rectangular CSG declares the feature and mask opening to be the same 90 nm
    interval.  This is geometry, not an experimental target or fitted value.
    """
    output = dict(metrics)
    opening = float(output["mask_opening_nm"])
    completed = []
    for name in ("top_feature_width_nm", "maximum_feature_width_nm"):
        if not np.isfinite(output[name]):
            output[name] = opening
            completed.append(name)
    return output, {
        "method": "analytic_rectangular_CSG_at_exact_zero_depth",
        "completed_fields": tuple(completed),
        "source_field": "mask_opening_nm",
    }


def _worker(args):
    output = Path(args.case_output)
    output.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    dx = float(args.worker_dx_nm) / 1000.0
    contract = grid_contract(dx)
    if not contract["compatible"]:
        payload = {
            "status": "blocked_incompatible_grid",
            "grid_contract": contract,
            "message": "case was not executed because it would change the physical domain",
        }
        _write_json(output / "audit.json", payload)
        return payload

    if args.worker_state == "initial":
        geometry = make_rectangular_trench_geometry_3d(
            cell_width=0.13,
            cell_length=0.02,
            domain_height=2.8,
            dx=dx,
            opening_width=0.09,
            mask_thickness=0.85,
            substrate_top=1.8,
            etched_depth=0.0,
        )
        state = None
        fingerprint = None
        absolute_start = 0.0
        source_step = 0
        restriction = {"method": "analytic_initial_geometry"}
        burst_duration = float(args.initial_duration_s)
    else:
        source = Path(args.source_checkpoint)
        geometry, state, fingerprint, metadata = pilot._load_checkpoint(source)
        source_provenance = _validate_source_provenance(source, metadata)
        absolute_start = float(metadata["physical_time_s"])
        source_step = int(metadata["step"])
        try:
            geometry, state, fingerprint, restriction = restrict_checkpoint_state(
                geometry, state, fingerprint, dx
            )
        except ValueError as error:
            payload = {
                "status": "blocked_non_equivalent_checkpoint",
                "state": "late",
                "dx_nm": float(args.worker_dx_nm),
                "grid_contract": contract,
                "source_checkpoint": {
                    "name": source.name,
                    "sha256": _sha256(source),
                    "physical_time_s": absolute_start,
                    "step": source_step,
                    "provenance": source_provenance,
                },
                "message": str(error),
                "scientific_action": (
                    "do not fabricate a coarse late state; use an earlier common-topology "
                    "checkpoint or compare instantaneous operators after AMR transfer exists"
                ),
            }
            _write_json(output / "audit.json", payload)
            return payload
        burst_duration = float(args.late_duration_s)

    start_metrics = pilot.measure_krueger_metrics(geometry, substrate_top_um=1.8)
    initial_metric_completion = None
    if args.worker_state == "initial":
        start_metrics, initial_metric_completion = _complete_analytic_initial_metrics(
            start_metrics
        )
    frozen = _advance(
        geometry,
        state,
        fingerprint,
        duration_s=0.0,
        seed=int(args.seed) + source_step,
        device=args.device,
        topology_policy=args.topology_policy,
        remap_backend=args.surface_state_remap_backend,
    )
    frozen_summary = _operator_summary(frozen)
    if state is None:
        # The zero-motion operator initializes the declared material state.  Use that exact state for
        # the burst while retaining the original geometry.
        state = frozen.next_surface_state
        fingerprint = frozen.next_surface_state_mesh_fingerprint
        geometry = frozen.geometry

    prescribed = None
    if args.worker_schedule:
        schedule_payload = json.loads(Path(args.worker_schedule).read_text(encoding="utf-8"))
        prescribed = [float(value) for value in schedule_payload["accepted_step_duration_s"]]
        if not np.isclose(sum(prescribed), burst_duration, rtol=0.0, atol=1e-12):
            raise ValueError("reference timestep schedule does not sum to requested burst duration")

    elapsed_physical = 0.0
    accepted = []
    records = []
    next_dt = min(float(args.maximum_step_s), burst_duration)
    status = "running"
    topology_event = None
    while elapsed_physical < burst_duration - 1e-15:
        if len(accepted) >= int(args.maximum_accepted_steps):
            status = "accepted_step_budget"
            break
        if perf_counter() - started >= float(args.max_wall_s):
            status = "wall_budget_checkpoint"
            break
        if prescribed is not None:
            if len(accepted) >= len(prescribed):
                raise RuntimeError("prescribed schedule ended before the burst")
            dt = prescribed[len(accepted)]
        else:
            dt = min(next_dt, burst_duration - elapsed_physical)
        rejected = []
        while True:
            try:
                result = _advance(
                    geometry,
                    state,
                    fingerprint,
                    duration_s=dt,
                    seed=int(args.seed) + source_step + len(accepted),
                    device=args.device,
                    topology_policy=args.topology_policy,
                    remap_backend=args.surface_state_remap_backend,
                )
            except SurfaceTopologyChangeError as error:
                topology_event = {
                    "event_kind": error.event_kind,
                    "old_topology": error.old_topology,
                    "new_topology": error.new_topology,
                    "message": str(error),
                }
                status = "topology_event"
                result = None
                break
            except (ValueError, RuntimeError) as error:
                retryable = str(error).startswith((
                    "surface topology changed under ",
                    "surface remap distance ",
                    "material surface appeared or disappeared",
                    "surface contraction exceeds bounded coverage capacity",
                ))
                if prescribed is not None or not retryable or dt <= float(args.minimum_step_s):
                    raise
                rejected.append({"duration_s": dt, "reason": str(error)})
                dt = max(float(args.minimum_step_s), 0.5 * dt)
                continue
            displacement = float(result.diagnostics["max_displacement_mesh_units"])
            limit = float(args.maximum_displacement_cells) * geometry.dx
            if displacement > limit:
                if prescribed is not None:
                    raise RuntimeError(
                        "coarse replay rejects the fine-owned timestep schedule on displacement"
                    )
                if dt <= float(args.minimum_step_s):
                    raise RuntimeError("displacement remains unresolved at minimum timestep")
                rejected.append({
                    "duration_s": dt,
                    "reason": f"displacement {displacement:g} exceeds {limit:g}",
                })
                dt = max(float(args.minimum_step_s), 0.5 * dt)
                continue
            break
        if result is None:
            break
        geometry = result.geometry
        state = result.next_surface_state
        fingerprint = result.next_surface_state_mesh_fingerprint
        elapsed_physical += dt
        accepted.append(dt)
        records.append({
            "step": len(accepted),
            "absolute_physical_time_s": absolute_start + elapsed_physical,
            "burst_time_s": elapsed_physical,
            "accepted_step_duration_s": dt,
            "rejected_trials": rejected,
            "metrics": pilot.measure_krueger_metrics(geometry, substrate_top_um=1.8),
            "operator": _operator_summary(result),
        })
        displacement = float(result.diagnostics["max_displacement_mesh_units"])
        if prescribed is None:
            target = float(args.target_displacement_cells) * geometry.dx
            factor = (
                float(args.adaptive_growth_factor)
                if displacement == 0.0 else
                float(np.clip(0.9 * target / displacement, 0.5, args.adaptive_growth_factor))
            )
            next_dt = float(np.clip(
                dt * factor, float(args.minimum_step_s), float(args.maximum_step_s)
            ))
        pilot._checkpoint(
            output / "checkpoint.npz",
            geometry,
            state,
            fingerprint,
            source_step + len(accepted),
            absolute_start + elapsed_physical,
            next_dt,
        )
    if elapsed_physical >= burst_duration - 1e-15:
        status = "complete"

    end_metrics = pilot.measure_krueger_metrics(geometry, substrate_top_um=1.8)
    frozen_end_summary = None
    if status == "complete":
        frozen_end = _advance(
            geometry,
            state,
            fingerprint,
            duration_s=0.0,
            seed=int(args.seed) + source_step + len(accepted),
            device=args.device,
            topology_policy=args.topology_policy,
            remap_backend=args.surface_state_remap_backend,
        )
        frozen_end_summary = _operator_summary(frozen_end)
    increments = {
        name: float(end_metrics[name]) - float(start_metrics[name])
        for name in (
            "etch_depth_nm", "mask_opening_nm", "top_feature_width_nm",
            "maximum_feature_width_nm", "remaining_mask_thickness_nm",
        )
    }
    rates = {
        name.replace("_nm", "_nm_s"): value / elapsed_physical
        for name, value in increments.items()
        if elapsed_physical > 0.0
    }
    payload = {
        "status": status,
        "state": args.worker_state,
        "dx_nm": float(args.worker_dx_nm),
        "grid_contract": contract,
        "operator": operator_contract(args.surface_state_remap_backend),
        "calibration_parameters": CALIBRATION,
        "seed_epoch_start": int(args.seed) + source_step,
        "source_checkpoint": (
            None if args.worker_state == "initial" else {
                "name": Path(args.source_checkpoint).name,
                "sha256": _sha256(args.source_checkpoint),
                "physical_time_s": absolute_start,
                "step": source_step,
                "provenance": source_provenance,
            }
        ),
        "restriction": restriction,
        "initial_metric_completion": initial_metric_completion,
        "schedule": {
            "owner": (
                "self_adaptive_5nm_reference"
                if prescribed is None and np.isclose(float(args.worker_dx_nm), 5.0)
                else "unpaired_self_adaptive_debug"
                if prescribed is None else "replayed_5nm"
            ),
            "accepted_step_duration_s": accepted,
            "sum_s": float(sum(accepted)),
        },
        "start_metrics": start_metrics,
        "end_metrics": end_metrics,
        "increments": increments,
        "rates": rates,
        "frozen_start_operator": frozen_summary,
        "frozen_end_operator": frozen_end_summary,
        "step_records": records,
        "topology_event": topology_event,
        "wall_time_s": perf_counter() - started,
    }
    _write_json(output / "audit.json", payload)
    _write_json(
        output / "schedule.json",
        {"accepted_step_duration_s": accepted, "sum_s": float(sum(accepted))},
    )
    return payload


def _worker_command(args, state, dx_nm, case_output, schedule=None):
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker-state", state,
        "--worker-dx-nm", f"{float(dx_nm):g}",
        "--case-output", str(case_output),
        "--source-checkpoint", str(args.source_checkpoint),
        "--initial-duration-s", f"{float(args.initial_duration_s):.17g}",
        "--late-duration-s", f"{float(args.late_duration_s):.17g}",
        "--maximum-step-s", f"{float(args.maximum_step_s):.17g}",
        "--minimum-step-s", f"{float(args.minimum_step_s):.17g}",
        "--target-displacement-cells", f"{float(args.target_displacement_cells):.17g}",
        "--maximum-displacement-cells", f"{float(args.maximum_displacement_cells):.17g}",
        "--adaptive-growth-factor", f"{float(args.adaptive_growth_factor):.17g}",
        "--maximum-accepted-steps", str(int(args.maximum_accepted_steps)),
        "--max-wall-s", f"{float(args.max_wall_s):.17g}",
        "--seed", str(int(args.seed)),
        "--device", str(args.device),
        "--topology-policy", str(args.topology_policy),
        "--surface-state-remap-backend", str(args.surface_state_remap_backend),
    ]
    if schedule is not None:
        command.extend(("--worker-schedule", str(schedule)))
    return command


def _plot(output, records):
    if not records:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 3, figsize=(11.0, 6.4), constrained_layout=True)
    measures = (
        ("etch_depth_nm_s", "Depth rate (nm/s)"),
        ("mask_opening_nm_s", "Opening rate (nm/s)"),
        ("maximum_feature_width_nm_s", "Max-width rate (nm/s)"),
    )
    for row, state in enumerate(("initial", "late")):
        selected = sorted(
            (item for item in records if item.get("state") == state and item.get("status") == "complete"),
            key=lambda item: item["dx_nm"], reverse=True,
        )
        for axis, (field, label) in zip(axes[row], measures):
            if selected:
                x = [item["dx_nm"] for item in selected]
                y = [item["rates"][field] for item in selected]
                axis.plot(x, y, "o-", color="#1f77b4")
                for xi, yi in zip(x, y):
                    axis.annotate(f"{yi:.3g}", (xi, yi), xytext=(0, 6),
                                  textcoords="offset points", ha="center", fontsize=8)
            axis.set_title(f"{state.capitalize()}: {label}")
            axis.set_xlabel("Uniform spacing (nm)")
            axis.set_ylabel(label)
            axis.grid(alpha=0.25)
            axis.invert_xaxis()
    figure.savefig(Path(output) / "resolution_diagnostic.png", dpi=180)
    plt.close(figure)


def _aggregate(output, plan):
    output = Path(output)
    records = []
    for state in ("initial", "late"):
        for dx_nm in (5, 10):
            path = output / f"{state}_{dx_nm}nm" / "audit.json"
            if path.is_file():
                records.append(json.loads(path.read_text(encoding="utf-8")))
    paired = {}
    for state in ("initial", "late"):
        by_dx = {item["dx_nm"]: item for item in records if item.get("state") == state}
        if 5.0 in by_dx and 10.0 in by_dx:
            fine = by_dx[5.0]
            coarse = by_dx[10.0]
            fields = set(fine.get("rates", ())) & set(coarse.get("rates", ()))
            paired[state] = {
                field: {
                    "fine_5nm": fine["rates"][field],
                    "coarse_10nm": coarse["rates"][field],
                    "coarse_minus_fine": coarse["rates"][field] - fine["rates"][field],
                    "relative_to_fine": (
                        (coarse["rates"][field] - fine["rates"][field])
                        / max(abs(fine["rates"][field]), 1e-30)
                    ),
                }
                for field in sorted(fields)
            }
    payload = {
        "plan": plan,
        "cases": records,
        "paired_10nm_vs_5nm": paired,
        "claim": (
            "diagnostic only; AMR must reproduce the 5 nm short-burst operator before it may "
            "replace the uniform reference"
        ),
    }
    _write_json(output / "audit.json", payload)
    _plot(output, records)
    return payload


def _worker_environment(device):
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    # ``transport_device`` controls Warp trajectory kernels, while the level-set
    # kernels select their device at module import from PETCH_DEVICE.  A worker
    # subprocess must bind both to the one declared campaign device; otherwise
    # a nominal CUDA run quietly traces on GPU and redistances on CPU.
    environment["PETCH_DEVICE"] = str(device)
    return environment


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    plan = build_plan(
        args.source_checkpoint,
        initial_duration_s=args.initial_duration_s,
        late_duration_s=args.late_duration_s,
        seed=args.seed,
        surface_state_remap_backend=args.surface_state_remap_backend,
    )
    _write_json(output / "plan.json", plan)
    if not args.execute:
        return {"status": "planned", "plan": plan}
    if not Path(args.source_checkpoint).is_file() and args.phase in ("late", "all"):
        raise ValueError("late audit requires the archived 5 nm checkpoint")

    states = ("initial", "late") if args.phase == "all" else (args.phase,)
    environment = _worker_environment(args.device)
    execution = []
    for state in states:
        planned_levels = {item["dx_nm"]: item for item in plan["levels"]}
        if (state == "late" and planned_levels[10.0]["late_case"] != "eligible"
                and not args.allow_unpaired_fine):
            execution.append({
                "state": "late",
                "status": "skipped_no_faithful_pair",
                "reason": planned_levels[10.0]["late_restriction_preflight"]["reason"],
                "action": (
                    "supply an earlier 5 nm checkpoint whose 10 nm restriction preserves topology, "
                    "or explicitly request --allow-unpaired-fine for a 5 nm-only diagnostic"
                ),
            })
            continue
        schedule = None
        execution_levels = (
            (5.0,) if state == "late" and planned_levels[10.0]["late_case"] != "eligible"
            else (5.0, 10.0)
        )
        for dx_nm in execution_levels:
            case_output = output / f"{state}_{int(dx_nm)}nm"
            command = _worker_command(args, state, dx_nm, case_output, schedule=schedule)
            try:
                completed = subprocess.run(
                    command,
                    cwd=ROOT,
                    env=environment,
                    check=False,
                    timeout=float(args.max_wall_s) + 30.0,
                    text=True,
                    capture_output=True,
                )
            except subprocess.TimeoutExpired as error:
                execution.append({
                    "state": state,
                    "dx_nm": dx_nm,
                    "status": "hard_wall_timeout",
                    "timeout_s": float(args.max_wall_s) + 30.0,
                    "stdout_tail": (error.stdout or "")[-4000:],
                    "stderr_tail": (error.stderr or "")[-4000:],
                })
                break
            execution.append({
                "state": state,
                "dx_nm": dx_nm,
                "returncode": completed.returncode,
                "stdout_tail": completed.stdout[-4000:],
                "stderr_tail": completed.stderr[-4000:],
            })
            audit = case_output / "audit.json"
            if completed.returncode != 0 or not audit.is_file():
                break
            record = json.loads(audit.read_text(encoding="utf-8"))
            if record.get("status") != "complete":
                break
            if dx_nm == 5.0:
                schedule = case_output / "schedule.json"
    _write_json(output / "execution.json", execution)
    return _aggregate(output, plan)


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default=ROOT / "results" / "krueger_2024_multiresolution_audit",
    )
    parser.add_argument("--source-checkpoint", default=DEFAULT_SOURCE)
    parser.add_argument("--phase", choices=("initial", "late", "all"), default="all")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--initial-duration-s", type=float, default=0.5)
    parser.add_argument("--late-duration-s", type=float, default=0.1)
    parser.add_argument("--maximum-step-s", type=float, default=0.025)
    parser.add_argument("--minimum-step-s", type=float, default=0.00025)
    parser.add_argument("--target-displacement-cells", type=float, default=0.35)
    parser.add_argument("--maximum-displacement-cells", type=float, default=0.75)
    parser.add_argument("--adaptive-growth-factor", type=float, default=1.5)
    parser.add_argument("--maximum-accepted-steps", type=int, default=200)
    parser.add_argument("--max-wall-s", type=float, default=900.0)
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument("--device", default="cpu")
    parser.add_argument(
        "--allow-unpaired-fine",
        action="store_true",
        help="run a 5 nm late diagnostic even when no topology-equivalent 10 nm pair exists",
    )
    parser.add_argument(
        "--topology-policy",
        choices=("refuse", "continue_gas_cavity"),
        default="refuse",
    )
    parser.add_argument(
        "--surface-state-remap-backend",
        choices=REMAP_BACKENDS,
        default="legacy_knn",
        help="explicit state-transfer operator; recorded in every plan and worker receipt",
    )
    parser.add_argument("--worker-state", choices=("initial", "late"), help=argparse.SUPPRESS)
    parser.add_argument("--worker-dx-nm", type=float, help=argparse.SUPPRESS)
    parser.add_argument("--worker-schedule", help=argparse.SUPPRESS)
    parser.add_argument("--case-output", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    controls = np.asarray((
        args.initial_duration_s,
        args.late_duration_s,
        args.maximum_step_s,
        args.minimum_step_s,
        args.target_displacement_cells,
        args.maximum_displacement_cells,
        args.adaptive_growth_factor,
        args.max_wall_s,
    ), dtype=float)
    if (np.any(~np.isfinite(controls)) or np.any(controls <= 0.0)
            or args.minimum_step_s > args.maximum_step_s
            or args.target_displacement_cells > args.maximum_displacement_cells
            or int(args.maximum_accepted_steps) <= 0):
        parser.error("invalid positive bounded-run controls")
    if bool(args.worker_state) != bool(args.case_output) or (
            args.worker_state and args.worker_dx_nm is None):
        parser.error("internal worker arguments must be supplied together")
    return args


if __name__ == "__main__":
    arguments = parse_args()
    result = _worker(arguments) if arguments.worker_state else run(arguments)
    print(json.dumps({
        "status": result.get("status", "aggregated"),
        "output": str(arguments.case_output or arguments.output),
    }, sort_keys=True))
