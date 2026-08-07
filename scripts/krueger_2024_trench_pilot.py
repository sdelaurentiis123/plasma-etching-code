#!/usr/bin/env python3
"""Bounded Krüger 2024 calibration/transfer replay through the unified 3-D engine.

This is a development replay, not a predictive validation claim. It uses the published HPEM
wafer fluxes, the digitized joint positive-ion IEAD, and the reduced Appendix-B oxide/mask
mechanisms. The source does not publish the species-resolved energetic mixture or the missing
three-dimensional ion azimuth, so those closures remain explicit in the output.

The script is intentionally operationally bounded:

* ``--benchmark-only`` measures one complete zero-motion transport/chemistry update;
* the evolving run checkpoints after every profile step and honors ``--max-wall-s``;
* every checkpoint contains the geometry and conservative material-surface state;
* ``--topology-change-policy continue_gas_cavity`` accepts only resolved periodic cavity
  enclosure/opening and records every event; the default still refuses every topology change;
* no charging solve, parameter fitting, or target-dependent velocity scaling occurs.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np
from scipy.ndimage import label

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.feature_step_3d import (
    FeatureGeometry3D,
    SurfaceTopologyChangeError,
    advance_feature_step_3d,
    make_rectangular_trench_geometry_3d,
)
from petch.material_mechanism_3d import MaterialSurfaceState3D
from petch.reactor_boundary import (
    build_krueger_2024_development_boundary,
    build_krueger_2024_transfer_boundary,
)


DATA = ROOT / "data" / "experimental" / "krueger_2024"
TARGET = {
    "mask_opening_nm": 45.0,
    "top_feature_width_nm": 90.0,
    "maximum_feature_width_nm": 90.0,
    "etch_depth_nm": 825.0,
    "remaining_mask_thickness_nm": 850.0,
    "asymmetry_cell_count": 0.0,
}


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


def _linear_zero(a, fa, b, fb):
    denominator = float(fa) - float(fb)
    if denominator == 0.0:
        return 0.5 * (float(a) + float(b))
    return float(a) + (float(b) - float(a)) * float(fa) / denominator


def _gas_interval(field, x, center_index):
    values = np.asarray(field, dtype=float)
    if values[center_index] >= 0.0:
        return None
    left = int(center_index)
    right = int(center_index)
    while left > 0 and values[left - 1] < 0.0:
        left -= 1
    while right + 1 < len(values) and values[right + 1] < 0.0:
        right += 1
    if left == 0 or right + 1 == len(values):
        return None
    x_left = _linear_zero(
        x[left - 1], values[left - 1], x[left], values[left])
    x_right = _linear_zero(
        x[right], values[right], x[right + 1], values[right + 1])
    return x_left, x_right


def _periodic_component_roots(field):
    """Label a duplicate-endpoint x/y-periodic Boolean volume."""
    occupied = np.asarray(field, dtype=bool)
    component, count = label(occupied)
    parent = np.arange(int(count) + 1)

    def find(index):
        index = int(index)
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    def union(left, right):
        left = find(left)
        right = find(right)
        if left and right and left != right:
            parent[right] = left

    for axis in (0, 1):
        first = np.take(component, 0, axis=axis)
        last = np.take(component, -1, axis=axis)
        selected = (first > 0) & (last > 0)
        for left, right in zip(first[selected], last[selected]):
            union(left, right)
    roots = np.zeros(int(count) + 1, dtype=int)
    for index in range(1, int(count) + 1):
        roots[index] = find(index)
    return roots[component]


def _mask_throat_connectivity(
        geometry, *, substrate_top_um, opening_center_um, opening_width_um):
    """Return whether exterior gas reaches the lowest resolved mask plane.

    The Krueger ``w_m`` observable concerns the mask throat, not whether a
    shallow etched cavity already owns a full gas node below the substrate.
    Requiring the latter falsely closes a real 6.8 nm etch on a 10 nm grid.
    The first grid plane strictly inside the mask is the stable lower witness;
    the upper witness is exterior gas at the top domain boundary.
    """
    # Geometry stores duplicate periodic x/y endpoints. Drop them and merge the
    # remaining opposite faces exactly as the production topology contract does.
    gas = np.asarray(geometry.phi[:-1, :-1, :], dtype=float) <= 0.0
    components = _periodic_component_roots(gas)
    upper = {
        int(value) for value in np.unique(components[:, :, -1])
        if value > 0
    }
    x, _, z = geometry.coordinate_arrays
    x = x[:-1]
    exit_candidates = np.flatnonzero(
        z > float(substrate_top_um) + 0.25 * geometry.dx)
    if not exit_candidates.size:
        raise ValueError("no resolved mask-interior exit plane")
    exit_index = int(exit_candidates[0])
    x_roi = np.flatnonzero(
        np.abs(x - float(opening_center_um))
        <= 0.5 * float(opening_width_um) + geometry.dx)
    if not x_roi.size:
        raise ValueError("declared opening ROI has no resolved lateral node")
    exit_roots = {
        int(value)
        for value in np.unique(components[x_roi, :, exit_index])
        if value > 0
    }
    return {
        "exterior_to_mask_exit_open": bool(upper & exit_roots),
        "mask_exit_z_um": float(z[exit_index]),
        "shared_component_count": int(len(upper & exit_roots)),
    }


def measure_krueger_metrics(
        geometry, *, substrate_top_um, opening_center_um=0.065,
        opening_width_um=0.09, aperture_profile_points=64):
    """Measure the six Fig. 7 scalars from the resolved center cross-section."""
    unit_to_nm = geometry.mesh_length_unit_m * 1.0e9
    x, y, z = geometry.coordinate_arrays
    center_index = int(np.argmin(np.abs(x - float(opening_center_um))))
    y_rows = range(1, len(y) - 1) if len(y) > 2 else range(len(y))

    floor_z = []
    for j in y_rows:
        values = geometry.phi[center_index, j]
        crossing = np.flatnonzero((values[:-1] >= 0.0) & (values[1:] < 0.0))
        for k in crossing:
            z_cross = _linear_zero(z[k], values[k], z[k + 1], values[k + 1])
            if z_cross <= substrate_top_um + 0.5 * geometry.dx:
                floor_z.append(z_cross)
    floor = float(np.median(floor_z)) if floor_z else np.nan
    etch_depth = max(0.0, float(substrate_top_um) - floor) if np.isfinite(floor) else np.nan

    mask_layer = (
        None if geometry.material_levelsets is None
        else np.asarray(geometry.material_levelsets[2], dtype=float))
    if mask_layer is None:
        raise ValueError("Krueger mask opening requires mask material level set 2")
    mask_top_z = []
    if mask_layer is not None:
        outside_opening = np.flatnonzero(
            np.abs(x - float(opening_center_um)) >= 0.055)
        for i in outside_opening:
            for j in y_rows:
                values = mask_layer[i, j]
                crossing = np.flatnonzero(
                    (values[:-1] >= 0.0) & (values[1:] < 0.0))
                for k in crossing:
                    top = _linear_zero(
                        z[k], values[k], z[k + 1], values[k + 1])
                    if top >= substrate_top_um:
                        mask_top_z.append(top)
    mask_top = float(np.median(mask_top_z)) if mask_top_z else np.nan
    mask_height = (
        max(0.0, mask_top - float(substrate_top_um))
        if np.isfinite(mask_top) else np.nan)

    mask_widths = []
    feature_widths = []
    feature_profile = []
    top_widths = []
    left_offsets = []
    right_offsets = []
    for k, z_value in enumerate(z):
        if (np.isfinite(mask_top)
                and substrate_top_um + 0.25 * geometry.dx
                <= z_value <= mask_top - 0.25 * geometry.dx):
            mask_intervals = [
                _gas_interval(mask_layer[:, j, k], x, center_index)
                for j in y_rows]
            mask_intervals = [
                item for item in mask_intervals if item is not None]
            if mask_intervals:
                per_y_width = np.asarray([
                    item[1] - item[0] for item in mask_intervals], dtype=float)
                mask_widths.append({
                    "width": float(np.mean(per_y_width)),
                    "z": float(z_value),
                    "cross_y_span": float(np.ptp(per_y_width)),
                    "per_y_widths": per_y_width,
                })
        intervals = []
        for j in y_rows:
            interval = _gas_interval(
                geometry.phi[:, j, k], x, center_index)
            if interval is not None:
                intervals.append(interval)
        if not intervals:
            continue
        left = float(np.mean([item[0] for item in intervals]))
        right = float(np.mean([item[1] for item in intervals]))
        width = right - left
        if (np.isfinite(floor)
                and floor + 0.25 * geometry.dx
                <= z_value <= substrate_top_um - 0.25 * geometry.dx):
            feature_widths.append(width)
            feature_profile.append({"width": width, "z": float(z_value)})
            left_offsets.append(float(opening_center_um) - left)
            right_offsets.append(right - float(opening_center_um))
            if z_value >= substrate_top_um - 1.5 * geometry.dx:
                top_widths.append(width)

    # Before the floor moves by one resolved row, the gas interval immediately above the
    # substrate is the physically declared top opening and is the only meaningful width.
    if not feature_widths:
        k = int(np.argmin(np.abs(z - (substrate_top_um + 0.5 * geometry.dx))))
        intervals = [
            _gas_interval(geometry.phi[:, j, k], x, center_index)
            for j in y_rows]
        intervals = [item for item in intervals if item is not None]
        if intervals:
            feature_widths = [
                float(np.mean([item[1] - item[0] for item in intervals]))]
            top_widths = list(feature_widths)
            feature_profile = [
                {"width": feature_widths[0], "z": float(z[k])}]

    asymmetry_cells = (
        float(np.max(np.abs(np.asarray(left_offsets) - np.asarray(right_offsets)))
              / geometry.dx)
        if left_offsets else 0.0)
    minimum_mask = (
        min(mask_widths, key=lambda item: item["width"])
        if mask_widths else None)
    connectivity = _mask_throat_connectivity(
        geometry, substrate_top_um=float(substrate_top_um),
        opening_center_um=float(opening_center_um),
        opening_width_um=float(opening_width_um))
    pocket_width_nm = (
        float(minimum_mask["width"]) * unit_to_nm
        if minimum_mask is not None else np.nan)
    paper_opening_nm = (
        pocket_width_nm
        if connectivity["exterior_to_mask_exit_open"] else 0.0)

    # Community-standard CD triple (Top CD / Necking CD / neck depth). The
    # legacy mask_opening_nm minimises over the MASK band only, so a pinch at
    # the mask top and a mid-mask neck collapse into one number; the frozen-
    # geometry probe (RESULTS_MOUTH_EQUILIBRIUM_PROBE_2026-08-02) measured
    # closure 30-50x faster at the mask top than at the 200-250 nm band where
    # the SEM and MCFPM neck, so size and location must be reported apart.
    aperture_profile = [
        {"z_um": float(item["z"]),
         "width_nm": float(item["width"]) * unit_to_nm}
        for item in (list(mask_widths) + list(feature_profile))]
    aperture_profile.sort(key=lambda item: item["z_um"], reverse=True)
    top_cd_nm = np.nan
    if mask_widths:
        top_cd_nm = (
            float(max(mask_widths, key=lambda item: item["z"])["width"])
            * unit_to_nm)
    elif aperture_profile:
        top_cd_nm = float(aperture_profile[0]["width_nm"])
    neck_cd_nm = neck_z_um = neck_depth_nm = np.nan
    if aperture_profile:
        neck = min(aperture_profile, key=lambda item: item["width_nm"])
        neck_cd_nm = float(neck["width_nm"])
        neck_z_um = float(neck["z_um"])
        if np.isfinite(mask_top):
            neck_depth_nm = (float(mask_top) - neck_z_um) * unit_to_nm
    # The audit keeps a coarse trace; diagnostics pass None for every plane.
    stride = (
        1 if aperture_profile_points is None
        else max(1, len(aperture_profile) // int(aperture_profile_points)))
    return {
        "mask_opening_nm": paper_opening_nm,
        "mask_pocket_width_nm": pocket_width_nm,
        "mask_opening_connected_to_exterior": bool(
            connectivity["exterior_to_mask_exit_open"]),
        "mask_opening_throat_z_um": (
            float(minimum_mask["z"]) if minimum_mask is not None else np.nan),
        "mask_opening_cross_y_span_nm": (
            float(minimum_mask["cross_y_span"]) * unit_to_nm
            if minimum_mask is not None else np.nan),
        "mask_opening_per_y_widths_nm": (
            (minimum_mask["per_y_widths"] * unit_to_nm).tolist()
            if minimum_mask is not None else []),
        "mask_exit_z_um": connectivity["mask_exit_z_um"],
        "top_feature_width_nm": (
            float(np.mean(top_widths)) * unit_to_nm if top_widths else np.nan),
        "maximum_feature_width_nm": (
            float(np.max(feature_widths)) * unit_to_nm if feature_widths else np.nan),
        "top_cd_nm": top_cd_nm,
        "neck_cd_nm": neck_cd_nm,
        "neck_z_um": neck_z_um,
        "neck_depth_from_mask_top_nm": neck_depth_nm,
        "aperture_profile": aperture_profile[::stride],
        "etch_depth_nm": float(etch_depth) * unit_to_nm,
        "remaining_mask_thickness_nm": float(mask_height) * unit_to_nm,
        "asymmetry_cell_count": asymmetry_cells,
        "floor_z_um": floor,
        "mask_top_z_um": mask_top,
    }


def _maximum_ledger_residual(exchange):
    values = []
    for name in exchange.removed_units_m2:
        residual = np.asarray(exchange.residual_units_m2(name), dtype=float)
        values.append(float(np.max(np.abs(residual))) if residual.size else 0.0)
    return max(values, default=0.0)


def _checkpoint(
        path, geometry, state, fingerprint, step, physical_time_s,
        next_step_duration_s):
    arrays = {
        "phi": np.asarray(geometry.phi),
        "material_id": np.asarray(geometry.material_id),
        "mesh_origin_m": np.asarray(geometry.mesh_origin_m),
    }
    for material_id, field in dict(geometry.material_levelsets or {}).items():
        arrays[f"material_levelset_{material_id}"] = np.asarray(field)
    for name, field in state.fields.items():
        arrays[f"state_{name}"] = np.asarray(field)
    metadata = {
        "step": int(step),
        "physical_time_s": float(physical_time_s),
        "next_step_duration_s": float(next_step_duration_s),
        "dx": float(geometry.dx),
        "mesh_length_unit_m": float(geometry.mesh_length_unit_m),
        "fingerprint": str(fingerprint),
        "state_upper_bounds": dict(state.upper_bounds),
        "state_remap_modes": dict(state.remap_modes),
    }
    arrays["metadata_json"] = np.asarray(json.dumps(metadata, sort_keys=True))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)


def _load_checkpoint(path):
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"]))
        layers = {
            int(name.rsplit("_", 1)[1]): archive[name]
            for name in archive.files if name.startswith("material_levelset_")}
        fields = {
            name[len("state_"):]: archive[name]
            for name in archive.files if name.startswith("state_")}
        geometry = FeatureGeometry3D(
            archive["phi"], archive["material_id"], metadata["dx"],
            metadata["mesh_length_unit_m"], tuple(archive["mesh_origin_m"]),
            material_levelsets=layers)
        state = MaterialSurfaceState3D(
            fields, metadata["state_upper_bounds"], metadata["state_remap_modes"])
    return geometry, state, metadata["fingerprint"], metadata


def _write_json(path, payload):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _plot(
        output_directory, initial_geometry, final_geometry, history, *,
        substrate_top_um, mask_thickness_um):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    row = final_geometry.phi.shape[1] // 2
    x = np.arange(final_geometry.phi.shape[0]) * final_geometry.dx * 1000.0
    z = np.arange(final_geometry.phi.shape[2]) * final_geometry.dx * 1000.0
    figure, axis = plt.subplots(figsize=(5.5, 8.0), constrained_layout=True)
    axis.pcolormesh(
        x, z, final_geometry.material_id[:, row, :].T,
        cmap=ListedColormap(["white", "#19a7ae", "#d84fbd"]),
        vmin=0, vmax=2, shading="nearest")
    axis.contour(
        x, z, initial_geometry.phi[:, row, :].T, levels=[0.0],
        colors="#777777", linestyles="--", linewidths=1.0)
    axis.contour(
        x, z, final_geometry.phi[:, row, :].T, levels=[0.0],
        colors="black", linewidths=1.3)
    axis.set_aspect("equal")
    axis.set_xlabel("Lateral position (nm)")
    axis.set_ylabel("Height (nm)")
    axis.set_title("Krüger campaign: dashed initial, black final")
    axis.set_ylim(
        max(0.0, (float(substrate_top_um) - 1.6) * 1000.0),
        (float(substrate_top_um) + float(mask_thickness_um) + 0.1) * 1000.0)
    figure.savefig(output_directory / "profile.png", dpi=180)
    plt.close(figure)

    time = np.asarray([row["physical_time_s"] for row in history], dtype=float)
    figure, axes = plt.subplots(3, 1, figsize=(7.0, 7.5), sharex=True,
                                constrained_layout=True)
    fields = (
        ("etch_depth_nm", "Etch depth (nm)"),
        ("mask_opening_nm", "Minimum mask opening (nm)"),
        ("maximum_feature_width_nm", "Maximum feature width (nm)"),
    )
    for axis, (field, label) in zip(axes, fields):
        axis.plot(time, [row["metrics"][field] for row in history], "o-", ms=3)
        axis.axhline(TARGET[field], color="#d95f02", ls="--", label="experiment")
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
    axes[0].legend()
    axes[-1].set_xlabel("Etch time (s)")
    figure.savefig(output_directory / "metric_trajectory.png", dpi=180)
    plt.close(figure)


def _git_state():
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "status", "--porcelain"], cwd=ROOT, check=True,
            capture_output=True, text=True).stdout.strip())
        return {"revision": revision, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"revision": None, "dirty": None}


def _configuration(args):
    configuration = {
        "boundary_case": str(args.boundary_case),
        "oxygen_to_fluorocarbon_ratio": (
            None if args.oxygen_ratio is None else float(args.oxygen_ratio)),
        "low_frequency_power_kw": (
            None if args.low_frequency_power_kw is None
            else float(args.low_frequency_power_kw)),
        "duration_s": float(args.duration_s),
        "n_steps": int(args.n_steps),
        "dx_um": float(args.dx_um),
        "n_position": int(args.n_position),
        "compressed_boundary_quadrature": bool(args.compress_boundary_quadrature),
        "neutral_speed_quadrature": "analytic_speed_marginal",
        "neutral_tensor_velocity_quadrature_active": False,
        "neutral_transverse_order": int(args.neutral_transverse_order),
        "neutral_normal_order": int(args.neutral_normal_order),
        "neutral_direction_polar_order": (
            int(args.neutral_direction_polar_order)
            if args.compress_boundary_quadrature else 12),
        "neutral_direction_azimuthal_order": (
            int(args.neutral_direction_azimuthal_order)
            if args.compress_boundary_quadrature else 24),
        "ion_energy_bin_eV": (
            float(args.ion_energy_bin_eV)
            if args.compress_boundary_quadrature else None),
        "ion_angle_bin_deg": (
            float(args.ion_angle_bin_deg)
            if args.compress_boundary_quadrature else None),
        "ion_azimuthal_closure": "axisymmetric_uniform",
        "ion_azimuthal_order": int(args.ion_azimuthal_order),
        "ballistic_transport": str(args.ballistic_transport),
        "transport_device": str(args.transport_device),
        "ballistic_face_quadrature_points": int(args.face_quadrature_points),
        "radiosity_rays_per_face": int(args.radiosity_rays),
        "radiosity_relative_tolerance": float(args.radiosity_tolerance),
        "radiosity_maximum_iterations": int(args.radiosity_max_iterations),
        "radiosity_enabled": not args.no_radiosity,
        "seed": int(args.seed),
        "radiosity_backend": str(args.radiosity_backend),
        "adaptive_profile_timestep": bool(args.adaptive_profile_timestep),
        "minimum_step_duration_s": float(args.minimum_step_s),
        "target_displacement_cells": float(args.target_displacement_cells),
        "maximum_displacement_cells": float(args.maximum_displacement_cells),
        "adaptive_shrink_factor": float(args.adaptive_shrink_factor),
        "adaptive_growth_factor": float(args.adaptive_growth_factor),
        "adaptive_safety_factor": float(args.adaptive_safety_factor),
        "maximum_accepted_steps": int(args.maximum_accepted_steps),
        "topology_change_policy": str(args.topology_change_policy),
        "surface_state_remap_backend": str(args.surface_state_remap_backend),
        "geometry": {
            "cell_width_um": 0.13,
            "cell_length_um": 0.02,
            "domain_height_um": float(args.domain_height_um),
            "opening_width_um": 0.09,
            "mask_thickness_um": 0.85,
            "substrate_top_um": float(args.substrate_top_um),
            "initial_etched_depth_um": 0.0,
        },
        "transport": (
            f"collisionless_{args.ballistic_transport}_joint_IEAD"),
        "charging": "disabled_for_Krueger_2024_calibration_and_transfer",
        "profile_reinitialization": "cr2",
        "surface_model": str(args.surface_model),
        "grazing_ion_reflection": str(args.grazing_ion_reflection),
        "mixed_layer_volatilization_yield": float(
            args.mixed_layer_volatilization_yield),
        "yield_energy_model": str(args.yield_energy_model),
        "deposition_layer_depth_nm": float(args.deposition_layer_depth_nm),
        "oxygen_half_saturation_flux_m2_s": float(
            args.oxygen_half_saturation_flux_m2_s),
    }
    # Preserve byte-for-byte compatibility with the already-running all-fresh checkpoint.  The
    # calibration closure enters the fingerprint only when promoted away from its historical zero.
    if float(args.effective_mask_crosslinked_growth_fraction) != 0.0:
        configuration["effective_mask_crosslinked_growth_fraction"] = float(
            args.effective_mask_crosslinked_growth_fraction)
    if float(args.oxide_etch_yield_scale) != 1.0:
        configuration["oxide_etch_yield_scale"] = float(
            args.oxide_etch_yield_scale)
    if str(args.guo_aggregate_ion_formula) != "unresolved":
        configuration["guo_aggregate_ion_formula"] = str(
            args.guo_aggregate_ion_formula)
    if str(args.surface_model) == "guo_tml":
        configuration["guo_translating_layer_thickness_nm"] = float(
            args.guo_translating_layer_thickness_nm)
    if args.radiosity_backend == "deterministic_extruded_2d":
        configuration["deterministic_exchange"] = {
            "exchange_method": str(args.exchange_method),
            "exchange_relative_tolerance": float(args.exchange_relative_tolerance),
            "exchange_geometry_tolerance": float(args.exchange_geometry_tolerance),
            "maximum_refinement_level": 24,
            "extrusion_projection_guard_cells": 0.02,
        }
    return configuration


def _monotone_resume_refinement(previous, current):
    """Allow only provenance-logged tightening of numerical safety limits on resume."""
    previous = json.loads(json.dumps(previous))
    current = json.loads(json.dumps(current))
    changes = {}
    horizon_keys_present = (
        "duration_s" in previous,
        "n_steps" in previous,
        "duration_s" in current,
        "n_steps" in current,
    )
    if any(horizon_keys_present) and not all(horizon_keys_present):
        return False, {}
    if all(horizon_keys_present):
        old_duration = float(previous.pop("duration_s"))
        new_duration = float(current.pop("duration_s"))
        old_steps = int(previous.pop("n_steps"))
        new_steps = int(current.pop("n_steps"))
        if (
            old_duration <= 0.0
            or old_steps <= 0
            or new_duration < old_duration
            or new_steps < old_steps
            or not np.isclose(
                new_duration / new_steps,
                old_duration / old_steps,
                rtol=0.0,
                atol=8.0 * np.finfo(float).eps
                * max(new_duration / new_steps, old_duration / old_steps),
            )
        ):
            return False, {}
        if new_duration != old_duration or new_steps != old_steps:
            changes["physical_time_horizon"] = {
                "old_duration_s": old_duration,
                "new_duration_s": new_duration,
                "old_nominal_steps": old_steps,
                "new_nominal_steps": new_steps,
                "nominal_step_duration_s": old_duration / old_steps,
                "classification": "monotone_execution_horizon_extension",
            }
    old_remap_declared = "surface_state_remap_backend" in previous
    new_remap_declared = "surface_state_remap_backend" in current
    old_remap_backend = str(previous.pop(
        "surface_state_remap_backend", "legacy_knn"))
    new_remap_backend = str(current.pop(
        "surface_state_remap_backend", "legacy_knn"))
    if old_remap_backend != new_remap_backend:
        # A remap backend is part of the physical/numerical evolution operator.  Never splice
        # two such operators into one trajectory, even if both conserve their own ledgers.
        return False, {}
    if not old_remap_declared and new_remap_declared:
        changes["surface_state_remap_backend_declaration"] = {
            "old": "implicit legacy_knn",
            "new": "explicit legacy_knn",
            "classification": "provenance_only_operator_declaration",
        }
    elif old_remap_declared != new_remap_declared:
        return False, {}
    old_topology_policy = str(previous.pop(
        "topology_change_policy", "refuse"))
    new_topology_policy = str(current.pop(
        "topology_change_policy", "refuse"))
    if old_topology_policy != new_topology_policy:
        if not (
                old_topology_policy == "refuse"
                and new_topology_policy == "continue_gas_cavity"):
            return False, {}
        changes["topology_change_policy"] = {
            "old": old_topology_policy,
            "new": new_topology_policy,
            "classification": "explicit_gas_cavity_continuation_scope",
        }
    old_minimum = float(previous.pop("minimum_step_duration_s"))
    new_minimum = float(current.pop("minimum_step_duration_s"))
    if new_minimum > old_minimum:
        return False, {}
    if new_minimum != old_minimum:
        changes["minimum_step_duration_s"] = {
            "old": old_minimum,
            "new": new_minimum,
            "classification": "monotone_numerical_refinement",
        }
    old_maximum = int(previous.pop("maximum_accepted_steps"))
    new_maximum = int(current.pop("maximum_accepted_steps"))
    if new_maximum < old_maximum:
        return False, {}
    if new_maximum != old_maximum:
        changes["maximum_accepted_steps"] = {
            "old": old_maximum,
            "new": new_maximum,
            "classification": "monotone_numerical_budget_extension",
        }
    return previous == current, changes


def _resume_topology_events(existing):
    """Retain a refused event as provenance when its exact pre-event checkpoint resumes."""
    events = list(existing.get("topology_events", ()))
    prior = existing.get("terminal_event")
    if (prior is not None
            and not any(
                item.get("source") == "prior_terminal_refusal"
                and item.get("physical_time_lower_s")
                == prior.get("physical_time_lower_s")
                for item in events)):
        events.append(dict(
            prior, accepted=False, source="prior_terminal_refusal"))
    return events


def _accepted_topology_event(event, *, physical_time_s, step_duration_s, step):
    """Attach one accepted geometry event to its physical-time interval."""
    if event is None:
        return None
    return dict(
        event,
        physical_time_lower_s=float(physical_time_s - step_duration_s),
        physical_time_upper_s=float(physical_time_s),
        accepted_step=int(step),
        source="accepted_feature_step")


def _visibility_history_summary(neutral_radiosity_diagnostics):
    """Collapse one shared form-factor receipt without multiplying it by species count."""
    diagnostics = tuple(neutral_radiosity_diagnostics.values())

    def maximum(name):
        return max((int(item.get(name, 0)) for item in diagnostics), default=0)

    def maximum_float(name):
        return max((float(item.get(name, 0.0)) for item in diagnostics), default=0.0)

    return {
        "maximum_visibility_float64_evaluated_count": maximum(
            "visibility_float64_evaluated_count"),
        "maximum_visibility_recovered_hit_count": maximum(
            "visibility_recovered_hit_count"),
        "maximum_visibility_source_relaunch_count": maximum(
            "visibility_source_relaunch_count"),
        "maximum_visibility_source_relaunch_distance_um": maximum_float(
            "visibility_maximum_source_relaunch_distance"),
        "maximum_visibility_overlap_skip_count": maximum(
            "visibility_overlap_skip_count"),
        "maximum_visibility_overlap_skip_depth_um": maximum_float(
            "visibility_maximum_overlap_skip_depth"),
        "maximum_visibility_unclassified_ray_count": maximum(
            "visibility_unclassified_ray_count"),
        "maximum_visibility_derived_horizon_extension_count": maximum(
            "visibility_derived_horizon_extension_count"),
        "maximum_visibility_wrap_count": maximum(
            "visibility_maximum_wrap_count"),
        "maximum_visibility_final_horizon_wraps": maximum(
            "visibility_final_maximum_wraps"),
        "maximum_visibility_source_support_face_count": maximum(
            "visibility_source_support_face_count"),
        "maximum_visibility_source_support_area_fraction": maximum_float(
            "visibility_source_support_area_fraction"),
        "maximum_visibility_source_support_distance_um": maximum_float(
            "visibility_maximum_source_support_distance"),
    }


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    geometry_control = np.asarray(
        [args.substrate_top_um, args.domain_height_um, args.dx_um], dtype=float)
    if (np.any(~np.isfinite(geometry_control)) or np.any(geometry_control <= 0.0)
            or not np.isfinite(args.oxide_etch_yield_scale)
            or args.oxide_etch_yield_scale <= 0.0
            or float(args.domain_height_um) - float(args.substrate_top_um)
            < 0.85 + 2.0 * float(args.dx_um)):
        raise ValueError(
            "Krüger domain must contain positive substrate depth, the 0.85 um mask, "
            "and resolved vacuum headspace")
    config = _configuration(args)
    config_hash = sha256(
        json.dumps(config, sort_keys=True).encode("utf-8")).hexdigest()

    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.13, cell_length=0.02,
        domain_height=float(args.domain_height_um),
        dx=args.dx_um, opening_width=0.09, mask_thickness=0.85,
        substrate_top=float(args.substrate_top_um), etched_depth=0.0)
    initial_geometry = geometry
    state = None
    fingerprint = None
    start_step = 0
    physical_time = 0.0
    nominal_step_duration = float(args.duration_s) / int(args.n_steps)
    next_step_duration = nominal_step_duration
    history = [{
        "step": 0,
        "physical_time_s": 0.0,
        "metrics": measure_krueger_metrics(
            geometry, substrate_top_um=float(args.substrate_top_um)),
        "step_wall_s": 0.0,
    }]
    checkpoint_path = output / "checkpoint.npz"
    audit_path = output / "audit.json"
    if args.resume:
        if not checkpoint_path.exists() or not audit_path.exists():
            raise ValueError("--resume requires checkpoint.npz and audit.json")
        geometry, state, fingerprint, metadata = _load_checkpoint(checkpoint_path)
        existing = json.loads(audit_path.read_text(encoding="utf-8"))
        resume_transitions = list(existing.get("resume_configuration_transitions", ()))
        if existing["config_hash"] != config_hash:
            compatible, changes = _monotone_resume_refinement(
                existing["configuration"], config)
            if not compatible or not changes:
                raise ValueError("resume configuration does not match the checkpoint")
            resume_transitions.append({
                "old_config_hash": existing["config_hash"],
                "new_config_hash": config_hash,
                "changes": changes,
                "reason": (
                    "resume from the exact accepted pre-event checkpoint under the explicit "
                    "gas-cavity continuation scope; boundary, chemistry, geometry state, and "
                    "surface inventories are unchanged"
                    if "topology_change_policy" in changes else
                    "resume from the exact accepted checkpoint with the historical implicit "
                    "legacy remap operator made explicit; operator and state are unchanged"
                    if "surface_state_remap_backend_declaration" in changes else
                    "resume from the exact accepted checkpoint after a numerical safety limit "
                    "was reached; physical operator and state are unchanged"),
            })
        start_step = int(metadata["step"])
        physical_time = float(metadata["physical_time_s"])
        next_step_duration = float(metadata.get(
            "next_step_duration_s", nominal_step_duration))
        history = list(existing["history"])
        topology_events = _resume_topology_events(existing)
    else:
        resume_transitions = []
        topology_events = []

    realized_domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(realized_domain[2])
    common_boundary = {
        "reference_plane_m": source_z * geometry.mesh_length_unit_m,
        "neutral_direction_polar_order": (
            int(args.neutral_direction_polar_order)
            if args.compress_boundary_quadrature else 12),
        "neutral_direction_azimuthal_order": (
            int(args.neutral_direction_azimuthal_order)
            if args.compress_boundary_quadrature else 24),
        "ion_energy_bin_eV": (
            float(args.ion_energy_bin_eV)
            if args.compress_boundary_quadrature else None),
        "ion_angle_bin_deg": (
            float(args.ion_angle_bin_deg)
            if args.compress_boundary_quadrature else None),
        "ion_azimuthal_closure": "axisymmetric_uniform",
        "ion_azimuthal_order": int(args.ion_azimuthal_order),
        # EXPLICIT TARGET FIT when != 1.0 (RESULTS_DEPTH_IDENTIFIABILITY):
        # scales only the aggregate ion flux. It is not a blanket calibration
        # or a measurement-derived correction.
        "ion_flux_normalization": float(args.ion_flux_normalization),
    }
    if args.boundary_case == "base":
        if args.oxygen_ratio is not None or args.low_frequency_power_kw is not None:
            raise ValueError("the base boundary does not accept transfer controls")
        boundary = build_krueger_2024_development_boundary(
            DATA, n_transverse_neutral=5, n_normal_neutral=8,
            **common_boundary)
    elif args.boundary_case == "oxygen_ratio":
        if args.oxygen_ratio is None:
            raise ValueError("oxygen-ratio transfer requires --oxygen-ratio")
        if (args.low_frequency_power_kw is not None
                and float(args.low_frequency_power_kw) != 6.0):
            raise ValueError("oxygen-ratio transfer is published only at 6 kW")
        boundary = build_krueger_2024_transfer_boundary(
            DATA, low_frequency_power_kw=6.0,
            oxygen_to_fluorocarbon_ratio=float(args.oxygen_ratio),
            **common_boundary)
    elif args.boundary_case == "power_sweep":
        if args.low_frequency_power_kw is None:
            raise ValueError("power transfer requires --low-frequency-power-kw")
        if args.oxygen_ratio is not None:
            raise ValueError("power transfer does not accept an oxygen-ratio override")
        boundary = build_krueger_2024_transfer_boundary(
            DATA, low_frequency_power_kw=float(args.low_frequency_power_kw),
            **common_boundary)
    else:  # argparse constrains this, but keep direct-call behavior fail-closed.
        raise ValueError("unknown Krüger boundary case")
    if str(args.surface_model) == "mixed_layer":
        mechanism = build_krueger_2024_material_router_3d(
            surface_model="mixed_layer",
            mixed_layer_volatilization_yield=float(
                args.mixed_layer_volatilization_yield))
    elif str(args.surface_model) == "guo_tml":
        mechanism = build_krueger_2024_material_router_3d(
            surface_model="guo_tml",
            guo_aggregate_ion_formula=(
                None
                if str(args.guo_aggregate_ion_formula) == "unresolved"
                else str(args.guo_aggregate_ion_formula)
            ),
            guo_translating_layer_thickness_nm=float(
                args.guo_translating_layer_thickness_nm),
            effective_mask_crosslinked_growth_fraction=float(
                args.effective_mask_crosslinked_growth_fraction),
            yield_energy_model=str(args.yield_energy_model),
            deposition_layer_depth_nm=float(args.deposition_layer_depth_nm),
            oxygen_half_saturation_flux_m2_s=(
                None if float(args.oxygen_half_saturation_flux_m2_s) <= 0.0
                else float(args.oxygen_half_saturation_flux_m2_s)))
    else:
        mechanism = build_krueger_2024_material_router_3d(
            effective_mask_crosslinked_growth_fraction=float(
                args.effective_mask_crosslinked_growth_fraction),
            oxide_etch_yield_scale=float(args.oxide_etch_yield_scale),
            yield_energy_model=str(args.yield_energy_model),
            deposition_layer_depth_nm=float(args.deposition_layer_depth_nm),
            oxygen_half_saturation_flux_m2_s=(
                None if float(args.oxygen_half_saturation_flux_m2_s) <= 0.0
                else float(args.oxygen_half_saturation_flux_m2_s)))
    role = {
        species.name: (
            "energetic_bombardment"
            if species.charge_number != 0 else "neutral_reactant")
        for species in boundary.species}
    radiosity = None
    if not args.no_radiosity:
        if args.radiosity_backend == "deterministic_extruded_2d":
            radiosity = {
                "form_factor_backend": "deterministic_extruded_2d",
                "periodic_lateral": True,
                "domain_size": realized_domain,
                "relative_tolerance": float(args.radiosity_tolerance),
                "maximum_iterations": int(args.radiosity_max_iterations),
                "deterministic_extruded_options": {
                    "exchange_method": str(args.exchange_method),
                    "exchange_relative_tolerance": float(
                        args.exchange_relative_tolerance),
                    # Headroom so a rare grazing fallback pair cannot exhaust the shared
                    # refinement budget at depth; declared and fingerprinted.
                    "maximum_refinement_level": 24,
                    # Section meshes are float64-exact after the extrusion projection;
                    # the float32-era default (128*eps32*scale ~ 1.5e-5) lets per-pair
                    # visibility slivers accumulate past the row-closure gate at depth
                    # (3.4e-4 measured at step 79). Declared and fingerprinted.
                    "geometry_tolerance": float(args.exchange_geometry_tolerance),
                },
            }
        else:
            radiosity = {
                "rays_per_face": int(args.radiosity_rays),
                "seed": int(args.seed) + 10000,
                "periodic_lateral": True,
                "domain_size": realized_domain,
                "relative_tolerance": float(args.radiosity_tolerance),
                "maximum_iterations": int(args.radiosity_max_iterations),
            }

    benchmark_only = bool(args.benchmark_only)
    accepted_step = int(start_step)
    run_started = perf_counter()
    status = "running"
    terminal_event = None
    maximum_projection_deviation = 0.0
    while (
            accepted_step < start_step + 1 if benchmark_only
            else physical_time < float(args.duration_s)):
        if accepted_step >= int(args.maximum_accepted_steps):
            raise RuntimeError(
                "pilot exhausted maximum accepted profile steps at "
                f"t={physical_time:.8g}/{float(args.duration_s):.8g} s")
        step_duration = (
            0.0 if benchmark_only else min(
                next_step_duration, float(args.duration_s) - physical_time))
        step_started = perf_counter()
        rejected_trials = []
        while True:
            try:
                result = advance_feature_step_3d(
                    geometry, boundary, role, mechanism,
                    etchable_material_ids=(1, 2), duration_s=step_duration,
                    source_bounds=(0.0, realized_domain[0], 0.0, realized_domain[1]),
                    source_z=source_z, surface_state=state,
                    surface_state_mesh_fingerprint=fingerprint,
                    n_position=int(args.n_position),
                    seed=int(args.seed) + accepted_step,
                    cfl_number=0.25, reinitialize=True, reinitialization_method="cr2",
                    profile_periodic_lateral=True,
                    transport_device=str(args.transport_device),
                    neutral_radiosity_options=radiosity,
                    ballistic_transport=args.ballistic_transport,
                    grazing_ion_reflection=(
                        {} if str(args.grazing_ion_reflection) == "literature_v1"
                        else None),
                    ballistic_face_quadrature_points=int(args.face_quadrature_points),
                    topology_change_policy=str(args.topology_change_policy),
                    surface_state_remap_backend=str(
                        args.surface_state_remap_backend))
            except (ValueError, RuntimeError) as error:
                message = str(error)
                retryable = (
                    message.startswith("surface topology changed under ")
                    or message.startswith("surface remap distance ")
                    or message.startswith("material surface appeared or disappeared")
                    or message.startswith(
                        "surface contraction exceeds bounded coverage capacity"))
                proposed = step_duration * float(args.adaptive_shrink_factor)
                at_minimum = step_duration <= float(args.minimum_step_s)
                if (isinstance(error, SurfaceTopologyChangeError)
                        and error.event_kind == "gas_cavity_enclosed"
                        and at_minimum):
                    terminal_event = {
                        "kind": "feature_clogged",
                        "geometry_event_kind": error.event_kind,
                        "topology_method": error.method,
                        "old_topology": error.old_topology,
                        "new_topology": error.new_topology,
                        "physical_time_lower_s": float(physical_time),
                        "physical_time_upper_s": float(
                            physical_time + step_duration),
                        "resolution_step_s": float(step_duration),
                        "classification": (
                            "resolved terminal process outcome; candidate closing step is not "
                            "remapped or accepted"),
                        "pre_event_metrics": measure_krueger_metrics(
                            geometry,
                            substrate_top_um=float(args.substrate_top_um)),
                        "message": message,
                    }
                    break
                if (benchmark_only or not args.adaptive_profile_timestep
                        or not retryable or at_minimum):
                    raise
                rejected_trials.append({
                    "duration_s": float(step_duration),
                    "reason": message,
                    "classification": "inline_recovery_retry",
                })
                step_duration = max(float(args.minimum_step_s), proposed)
                continue
            if args.adaptive_profile_timestep and not benchmark_only:
                displacement = float(
                    result.diagnostics["max_displacement_mesh_units"])
                limit = float(args.maximum_displacement_cells) * geometry.dx
                if displacement > limit:
                    proposed = (
                        step_duration * float(args.adaptive_safety_factor) * limit
                        / displacement)
                    if step_duration <= float(args.minimum_step_s):
                        raise RuntimeError(
                            "pilot coupling displacement remains unresolved at "
                            f"minimum dt={float(args.minimum_step_s):.8g} s")
                    rejected_trials.append({
                        "duration_s": float(step_duration),
                        "reason": (
                            f"coupling displacement {displacement:.8g} exceeds "
                            f"{limit:.8g} mesh units"),
                        "classification": "inline_recovery_retry",
                    })
                    step_duration = max(
                        float(args.minimum_step_s),
                        min(
                            step_duration * float(args.adaptive_shrink_factor),
                            proposed))
                    continue
            break
        if terminal_event is not None:
            status = "terminal_feature_clogged"
            break
        step_wall = perf_counter() - step_started
        geometry = result.geometry
        state = result.next_surface_state
        fingerprint = result.next_surface_state_mesh_fingerprint
        step_projection = result.state_remap_diagnostics.get(
            "extrusion_projection_deviation_mesh_units")
        if step_projection is not None:
            maximum_projection_deviation = max(
                maximum_projection_deviation, float(step_projection))
        if not benchmark_only:
            physical_time += step_duration
        accepted_step += 1
        if args.adaptive_profile_timestep and not benchmark_only:
            displacement = float(result.diagnostics["max_displacement_mesh_units"])
            target = float(args.target_displacement_cells) * geometry.dx
            if displacement > 0.0:
                factor = (
                    float(args.adaptive_safety_factor) * target / displacement)
                factor = float(np.clip(
                    factor, float(args.adaptive_shrink_factor),
                    float(args.adaptive_growth_factor)))
            else:
                factor = float(args.adaptive_growth_factor)
            next_step_duration = float(np.clip(
                step_duration * factor, float(args.minimum_step_s),
                nominal_step_duration))
        else:
            next_step_duration = nominal_step_duration
        radiosity_balance = max(
            (float(item["relative_balance_error"])
             for item in result.diagnostics["neutral_radiosity"].values()),
            default=0.0)
        radiosity_rays = max(
            (int(item["form_factor_rays_per_face"])
             for item in result.diagnostics["neutral_radiosity"].values()),
            default=0)
        radiosity_refinements = max(
            (int(item["form_factor_refinement_count"])
             for item in result.diagnostics["neutral_radiosity"].values()),
            default=0)
        visibility_summary = _visibility_history_summary(
            result.diagnostics["neutral_radiosity"])
        history.append({
            "step": int(accepted_step),
            "physical_time_s": float(physical_time),
            "accepted_step_duration_s": float(step_duration),
            "rejected_trials": rejected_trials,
            "metrics": measure_krueger_metrics(
                geometry, substrate_top_um=float(args.substrate_top_um)),
            "step_wall_s": float(step_wall),
            "max_velocity_m_s": float(result.diagnostics["max_velocity_m_s"]),
            "raw_maximum_face_velocity_m_s": float(
                result.diagnostics["raw_maximum_face_velocity_m_s"]),
            "cfl_substeps": int(result.diagnostics["cfl_substeps"]),
            "interior_gas_nucleation_suppressed_cells": int(
                result.diagnostics.get(
                    "interior_gas_nucleation_suppressed_cells", 0)),
            "removed_unresolved_solid_cells": int(
                result.diagnostics["removed_unresolved_solid_cells"]),
            "reassigned_unresolved_material_nodes": int(
                result.diagnostics["reassigned_unresolved_material_nodes"]),
            "unresolved_material_volume_upper_bound_m3": float(
                result.diagnostics[
                    "unresolved_material_volume_upper_bound_m3"]),
            "filled_unresolved_gas_cavity_cells": int(
                result.diagnostics["filled_unresolved_gas_cavity_cells"]),
            "unresolved_gas_cavity_volume_upper_bound_m3": float(
                result.diagnostics[
                    "unresolved_gas_cavity_volume_upper_bound_m3"]),
            "maximum_material_ledger_residual_units_m2": (
                _maximum_ledger_residual(result.surface.material_exchange)),
            "maximum_neutral_radiosity_relative_balance_error": radiosity_balance,
            "maximum_neutral_radiosity_rays_per_face": radiosity_rays,
            "maximum_neutral_radiosity_refinement_count": radiosity_refinements,
            **visibility_summary,
            "neutral_radiosity_solver_methods": {
                name: item["solver_method"]
                for name, item in result.diagnostics["neutral_radiosity"].items()},
            "maximum_neutral_radiosity_inactive_face_count": max(
                (int(item["inactive_face_count"])
                 for item in result.diagnostics["neutral_radiosity"].values()),
                default=0),
            "hit_probability": dict(result.transport.hit_probability),
            "topology_event": result.diagnostics["topology_event"],
            "validity": {
                "within_declared_scope": result.validity.within_declared_scope,
                "parameter_evidence_supports_prediction": (
                    result.validity.parameter_evidence_supports_prediction),
                "reasons": result.validity.reasons,
                "nonpredictive_parameters": result.validity.nonpredictive_parameters,
            },
        })
        if result.diagnostics["topology_event"] is not None:
            topology_events.append(_accepted_topology_event(
                result.diagnostics["topology_event"],
                physical_time_s=physical_time,
                step_duration_s=step_duration,
                step=accepted_step))
        _checkpoint(
            checkpoint_path, geometry, state, fingerprint, accepted_step,
            physical_time, next_step_duration)
        elapsed = perf_counter() - run_started
        status = "benchmark_complete" if benchmark_only else "running"
        payload = {
            "status": status,
            "config_hash": config_hash,
            "configuration": config,
            "target": TARGET,
            "history": history,
            "topology_events": topology_events,
            "resume_configuration_transitions": resume_transitions,
            "boundary_provenance": dict(boundary.provenance),
            "mechanism_provenance": mechanism.provenance,
            "git": _git_state(),
            "wall_time_s": elapsed,
            "extrusion_projection_max_deviation_mesh_units": maximum_projection_deviation,
            "scientific_status": (
                "calibrated development transfer only; the aggregate energetic ion mixture "
                "is unresolved, 3-D ion azimuth uses the R1.5 axisymmetric closure certified "
                "at 16 versus 32 nodes, and reduced chemistry omits explicit species-resolved "
                "crosslinking/redeposition channels"),
        }
        _write_json(audit_path, payload)
        print(
            f"step={accepted_step} t={physical_time:.3f}s dt={step_duration:.4g}s "
            f"depth={history[-1]['metrics']['etch_depth_nm']:.3f}nm "
            f"mask_opening={history[-1]['metrics']['mask_opening_nm']:.3f}nm "
            f"wall={step_wall:.2f}s", flush=True)
        if not benchmark_only and elapsed >= float(args.max_wall_s):
            status = "wall_budget_checkpoint"
            break
    if not benchmark_only and physical_time >= float(args.duration_s):
        status = "complete"

    if audit_path.exists():
        payload = json.loads(audit_path.read_text(encoding="utf-8"))
    else:
        payload = {
            "config_hash": config_hash,
            "configuration": config,
            "target": TARGET,
            "history": history,
            "topology_events": topology_events,
            "resume_configuration_transitions": resume_transitions,
            "boundary_provenance": dict(boundary.provenance),
            "mechanism_provenance": mechanism.provenance,
            "git": _git_state(),
        }
    payload["status"] = status
    payload["wall_time_s"] = perf_counter() - run_started
    payload["final_metrics"] = history[-1]["metrics"]
    payload["target_error"] = {
        name: history[-1]["metrics"][name] - target
        for name, target in TARGET.items()}
    payload["terminal_event"] = terminal_event
    payload["topology_events"] = topology_events
    _write_json(audit_path, payload)
    _plot(
        output, initial_geometry, geometry, history,
        substrate_top_um=float(args.substrate_top_um), mask_thickness_um=0.85)
    return payload


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", default=ROOT / "results" / "krueger_2024_trench_pilot")
    parser.add_argument("--duration-s", type=float, default=60.0)
    parser.add_argument("--n-steps", type=int, default=60)
    parser.add_argument("--dx-um", type=float, default=0.01)
    parser.add_argument(
        "--substrate-top-um", type=float, default=1.8,
        help="inert substrate depth below the initial interface; translated with the mask")
    parser.add_argument(
        "--domain-height-um", type=float, default=2.8,
        help="domain height; defaults retain 0.15 um vacuum above the translated mask")
    parser.add_argument("--n-position", type=int, default=16)
    parser.add_argument(
        "--boundary-case", choices=("base", "oxygen_ratio", "power_sweep"),
        default="base",
        help="base is calibration; oxygen_ratio and power_sweep are held-out transfer inputs")
    parser.add_argument("--oxygen-ratio", type=float)
    parser.add_argument("--low-frequency-power-kw", type=float)
    parser.add_argument(
        "--compress-boundary-quadrature", action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "march with the endpoint-certified angular-neutral and joint-IEAD reduction; "
            "disable for the 12x24 angular-refinement and exact digitized-IEAD reference"),
    )
    # Retained for explicit legacy tensor studies only.  Production reduction analytically
    # marginalizes neutral speed and resolves p(mu,phi)=2mu/(2pi) directly; 8x16 passed the frozen
    # endpoint operator gate against 12x24 refinement on the completed base trench.
    parser.add_argument("--neutral-transverse-order", type=int, default=5)
    parser.add_argument("--neutral-normal-order", type=int, default=2)
    parser.add_argument("--neutral-direction-polar-order", type=int, default=8)
    parser.add_argument("--neutral-direction-azimuthal-order", type=int, default=16)
    parser.add_argument("--ion-energy-bin-eV", type=float, default=250.0)
    parser.add_argument("--ion-angle-bin-deg", type=float, default=0.25)
    parser.add_argument(
        "--ion-azimuthal-order", type=int, default=16,
        help=(
            "uniform azimuthal quadrature order for the explicit axisymmetric lift of the "
            "published polar IEAD"))
    parser.add_argument(
        "--effective-mask-crosslinked-growth-fraction", type=float, default=0.0,
        help=(
            "fixed [0,1] base-SEM calibration blend between Appendix-B fresh and "
            "crosslinked mask-film radical attachment; held-out runs must reuse it"))
    parser.add_argument(
        "--oxide-etch-yield-scale", type=float, default=1.0,
        help=(
            "positive base-depth calibration multiplier on the published bare/complex SiO2 "
            "yield amplitudes; held-out runs must reuse it"))
    parser.add_argument(
        "--mixed-layer-volatilization-yield", type=float, default=1.0,
        help=(
            "mixed-layer absolute rate constant k_v (substrate removals per "
            "ion at the reference deposited energy); the single "
            "base-condition-anchored constant of the mixed-layer chemistry"))
    parser.add_argument(
        "--ion-flux-normalization", type=float, default=1.0,
        help=("DECLARED CALIBRATION of the aggregate ion flux (default 1.0 = "
              "published HPEM value untouched). Scales ions only; beam-measured "
              "yields and all neutral fluxes are unaffected."))
    parser.add_argument(
        "--grazing-ion-reflection", default="off",
        choices=("off", "literature_v1"),
        help=(
            "split grazing-incidence ion weight into single-bounce specular "
            "hot neutrals (P=0.95(1-cos^3), 0.90 energy retention; "
            "Helmer/Graves-bounded)"))
    parser.add_argument(
        "--surface-model", default="reduced",
        choices=("reduced", "mixed_layer", "guo_tml"),
        help=(
            "surface chemistry family: the reduced coverage/site-balance "
            "mechanisms, the element-resolved two-reservoir mixed-layer "
            "chemistry, or the source-fixed Guo/Kwon translating-layer "
            "transfer audit (no oxide yield/depth knob)"))
    parser.add_argument(
        "--guo-aggregate-ion-formula",
        default="unresolved",
        choices=("unresolved", "C", "CF", "CF2", "CF3", "C3F3"),
        help=(
            "explicit all-one-species sensitivity for Krueger's unpublished "
            "aggregate positive-ion row; unresolved is the nominal "
            "non-incorporating endpoint and no choice is a fitted mixture"
        ),
    )
    parser.add_argument(
        "--guo-translating-layer-thickness-nm",
        type=float,
        default=2.5,
        help=(
            "Guo finite-fluence TML capacity scale. The 2.5 nm default is "
            "source-fixed from Guo's profile discretization and is not fit "
            "to Krueger depth; 1.2-3.0 nm is the declared source sensitivity "
            "band."
        ),
    )
    parser.add_argument(
        "--yield-energy-model", default="threshold_power",
        choices=("threshold_power", "deposited_in_layer"),
        help=(
            "ion yield energy law: the published threshold power law, or the derived "
            "ZBL deposited-in-layer (Sigmund) form anchored at the reference energy"))
    parser.add_argument("--deposition-layer-depth-nm", type=float, default=1.5)
    parser.add_argument(
        "--oxygen-half-saturation-flux-m2-s", type=float, default=0.0,
        help="Langmuir half-saturation of the adsorbed-O channel; <=0 disables")
    parser.add_argument(
        "--ballistic-transport", choices=("forward", "face_gather"),
        default="face_gather")
    parser.add_argument(
        "--transport-device", default="cpu",
        help="Warp transport device (for example cpu or cuda:0); recorded in the run hash")
    parser.add_argument("--face-quadrature-points", type=int, default=3)
    parser.add_argument("--radiosity-rays", type=int, default=8)
    parser.add_argument(
        "--radiosity-backend", default="scrambled_qmc_3d",
        choices=("scrambled_qmc_3d", "deterministic_extruded_2d"))
    parser.add_argument("--radiosity-tolerance", type=float, default=1e-12)
    parser.add_argument(
        "--exchange-method", default="analytic_occlusion",
        choices=("analytic_occlusion", "adaptive_refinement"),
        help="deterministic extruded exchange construction: exact projective-interval "
             "occlusion with certified outer quadrature (per-pair fallback to adaptive "
             "refinement), or refinement everywhere (cross-check mode); recorded in the "
             "operator fingerprint")
    parser.add_argument(
        "--exchange-geometry-tolerance", type=float, default=1.0e-9,
        help="geometric predicate tolerance (mesh units) for the deterministic extruded "
             "exchange; float64 section meshes support 1e-9, and looser values leak "
             "per-pair visibility slivers that accumulate against the row-closure gate")
    parser.add_argument(
        "--exchange-relative-tolerance", type=float, default=1.0e-5,
        help="declared per-pair shadow-refinement budget for adaptive-refinement and "
             "taut-string fallback pairs; recorded in the operator fingerprint and the "
             "per-pair estimated_absolute_error receipt")
    parser.add_argument("--radiosity-max-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument("--max-wall-s", type=float, default=900.0)
    parser.add_argument(
        "--adaptive-profile-timestep", action=argparse.BooleanOptionalAction,
        default=True)
    parser.add_argument("--minimum-step-s", type=float, default=1e-5)
    parser.add_argument("--target-displacement-cells", type=float, default=0.35)
    parser.add_argument("--maximum-displacement-cells", type=float, default=0.75)
    parser.add_argument("--adaptive-shrink-factor", type=float, default=0.5)
    parser.add_argument("--adaptive-growth-factor", type=float, default=1.5)
    parser.add_argument("--adaptive-safety-factor", type=float, default=0.9)
    parser.add_argument("--maximum-accepted-steps", type=int, default=10000)
    parser.add_argument(
        "--topology-change-policy",
        choices=("refuse", "continue_gas_cavity"), default="refuse",
        help=(
            "refuse every resolved topology change, or explicitly continue only periodic "
            "gas-cavity enclosure/opening with conservative surface-state remap"))
    parser.add_argument(
        "--surface-state-remap-backend",
        choices=(
            "legacy_knn", "indexed_knn", "partitioned_overlap",
            "common_refinement"),
        default="legacy_knn",
        help=(
            "declared conservative surface-state transfer operator; it is fingerprinted and "
            "cannot change when resuming a trajectory"))
    parser.add_argument("--benchmark-only", action="store_true")
    parser.add_argument("--no-radiosity", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
