#!/usr/bin/env python3
"""Evolve blind, conditional TiO2 square-pillar profiles for the Oxford run.

The reactor boundary and geometry board are frozen before the target SEM.  The
surface scale in this rung comes only from the explicitly cross-machine
Janissen rate witnesses, so every output is a conditional profile envelope and
not a target-machine absolute-depth prediction.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.feature_step_3d import (
    FeatureGeometry3D,
    SurfaceTopologyChangeError,
    advance_feature_step_3d,
    make_square_pillar_mask_geometry_3d,
)
from petch.iadf_two_component import (
    build_two_component_boundary,
    kim_2025_reference_iadf,
)
from petch.rate_normalized_removal import (
    RateNormalizedRemovalMechanism,
    RateNormalizedRemovalParameters,
)
from petch.surface_kinetics import ParameterEvidence
from petch.tio2_square_pillar import tio2_formula_unit_density_m3


DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
PREREGISTRATION = DATA / "square_pillar_blind_preregistration.json"
ANALOG_BOARD = DATA / "janissen_tio2_analog_board.json"
REACTOR_DOSE = (
    ROOT / "results" / "curated" / "zhu_npg80_daughter_wafer_dose_v1"
    / "audit.json"
)
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_conditional_profiles_v1"
    / "audit.json"
)
CACHE_DIR = OUTPUT.parent / "trajectories"
MODEL_REVISION = "external-union-active-band-dose-factorization-v1"

FILM_MATERIAL = 1
MASK_MATERIAL = 2
BASE_MATERIAL = 3
BASE_TOP_UM = 0.1
FILM_THICKNESS_UM = 0.7
MASK_THICKNESS_UM = 0.045
FILM_TOP_UM = BASE_TOP_UM + FILM_THICKNESS_UM
DOMAIN_HEIGHT_UM = 1.05
SOURCE_Z_UM = 1.0
MESH_SPACING_NM = 20.0


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _render(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _trajectory_cache_spec(job, preregistration):
    width, scenario, target_rates, duration, dx = job
    return {
        "model_revision": MODEL_REVISION,
        "preregistration_sha256": _hash(PREREGISTRATION),
        "width_nm": float(width),
        "scenario": dict(scenario),
        "target_rates_nm_min": [float(value) for value in target_rates],
        "duration_s": float(duration),
        "mesh_spacing_nm": float(dx),
    }


def _trajectory_cache_path(spec):
    digest = sha256(_render(spec).encode("utf-8")).hexdigest()[:16]
    width = int(round(spec["width_nm"]))
    scenario = spec["scenario"]["name"]
    return CACHE_DIR / f"w{width:03d}_{scenario}_{digest}.json"


def _load_cached_trajectory(job, preregistration):
    spec = _trajectory_cache_spec(job, preregistration)
    path = _trajectory_cache_path(spec)
    if not path.exists():
        return None, path, spec
    payload = _load(path)
    if payload.get("job_spec") != spec:
        raise RuntimeError(f"trajectory cache specification mismatch: {path}")
    return payload["profiles"], path, spec


def _maximum_named_numeric(value, name):
    found = []
    if isinstance(value, dict) or hasattr(value, "items"):
        for key, item in dict(value).items():
            if key == name and np.isscalar(item):
                found.append(abs(float(item)))
            else:
                found.append(_maximum_named_numeric(item, name))
    elif isinstance(value, (tuple, list)):
        found.extend(_maximum_named_numeric(item, name) for item in value)
    return max(found, default=0.0)


def _linear_zero(x0, x1, y0, y1):
    if y1 == y0:
        return 0.5 * (x0 + x1)
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def _surface_height_um(column, z):
    """Highest gas/film zero crossing in one vertical signed-distance column."""
    column = np.asarray(column, dtype=float)
    crossings = []
    for index in range(len(z) - 1):
        if column[index] == 0.0:
            crossings.append(float(z[index]))
        if column[index] * column[index + 1] < 0.0:
            crossings.append(float(_linear_zero(
                z[index], z[index + 1], column[index], column[index + 1]
            )))
    if not crossings:
        return np.nan
    return max(crossings)


def _center_width_um(line, coordinate, center):
    """Width of the positive film interval containing the periodic-cell center."""
    line = np.asarray(line, dtype=float)
    inside = line >= 0.0
    middle = int(np.argmin(np.abs(coordinate - center)))
    if not inside[middle]:
        return 0.0
    left = middle
    while left > 0 and inside[left - 1]:
        left -= 1
    right = middle
    while right + 1 < len(line) and inside[right + 1]:
        right += 1
    left_edge = float(coordinate[left])
    right_edge = float(coordinate[right])
    if left > 0:
        left_edge = float(_linear_zero(
            coordinate[left - 1], coordinate[left], line[left - 1], line[left]
        ))
    if right + 1 < len(line):
        right_edge = float(_linear_zero(
            coordinate[right], coordinate[right + 1], line[right], line[right + 1]
        ))
    return max(0.0, right_edge - left_edge)


def _line_at_height(field, height, dx, *, axis, fixed_index):
    scaled = float(np.clip(height / dx, 0.0, field.shape[2] - 1.0))
    lower = int(np.floor(scaled))
    upper = min(lower + 1, field.shape[2] - 1)
    fraction = scaled - lower
    if axis == 0:
        first = field[:, fixed_index, lower]
        second = field[:, fixed_index, upper]
    else:
        first = field[fixed_index, :, lower]
        second = field[fixed_index, :, upper]
    return (1.0 - fraction) * first + fraction * second


def _profile_metrics(geometry, *, pitch_um, nominal_width_um):
    film = np.asarray(
        geometry.phi if geometry.material_levelsets is None
        else geometry.material_levelsets[FILM_MATERIAL],
        dtype=float,
    )
    nx, ny, nz = film.shape
    x = np.arange(nx) * geometry.dx
    y = np.arange(ny) * geometry.dx
    z = np.arange(nz) * geometry.dx
    center = 0.5 * pitch_um
    X, Y = np.meshgrid(x, y, indexing="ij")
    clearance = max(1.5 * geometry.dx, 0.1 * nominal_width_um)
    open_region = (
        (np.abs(X - center) > 0.5 * nominal_width_um + clearance)
        | (np.abs(Y - center) > 0.5 * nominal_width_um + clearance)
    )
    heights = np.asarray([
        _surface_height_um(film[ix, iy, :], z)
        for ix, iy in zip(*np.where(open_region))
    ])
    heights = heights[np.isfinite(heights)]
    if heights.size == 0:
        raise RuntimeError("no exposed-film surface heights resolved")
    floor_height = float(np.median(heights))
    floor_spread = float(np.percentile(heights, 95) - np.percentile(heights, 5))
    depth = float(np.clip(FILM_TOP_UM - floor_height, 0.0, FILM_THICKNESS_UM))
    center_y = int(np.argmin(np.abs(y - center)))
    center_x = int(np.argmin(np.abs(x - center)))
    # Sample only within the actual relief.  Sampling over ``max(depth, dx)``
    # would reach below the etched floor for a sub-cell smoke depth and mistake
    # the intact blanket for a full-cell pillar.  The result is still marked
    # non-authoritative until at least two vertical cells of relief exist.
    resolved_relief = depth >= 2.0 * geometry.dx
    fractions = np.linspace(0.10, 0.90, 17)
    cross_section = []
    for fraction in fractions:
        height = FILM_TOP_UM - fraction * depth
        width_x = _center_width_um(
            _line_at_height(
                film, height, geometry.dx, axis=0, fixed_index=center_y
            ),
            x,
            center,
        )
        width_y = _center_width_um(
            _line_at_height(
                film, height, geometry.dx, axis=1, fixed_index=center_x
            ),
            y,
            center,
        )
        cross_section.append({
            "relief_fraction_from_top": float(fraction),
            "height_um": float(height),
            "width_x_nm": float(width_x * 1.0e3),
            "width_y_nm": float(width_y * 1.0e3),
            "mean_width_nm": float(0.5 * (width_x + width_y) * 1.0e3),
        })
    widths = np.asarray([row["mean_width_nm"] for row in cross_section])
    top = float(np.mean(widths[:3]))
    middle_index = len(widths) // 2
    middle = float(np.mean(widths[middle_index - 1:middle_index + 2]))
    bottom = float(np.mean(widths[-3:]))
    lateral_change = 0.5 * abs(top - bottom)
    sidewall_angle = 90.0 if lateral_change == 0.0 else float(
        np.degrees(np.arctan2(depth * 1.0e3, lateral_change))
    )
    bow = float(middle - 0.5 * (top + bottom))
    return {
        "etched_depth_nm": depth * 1.0e3,
        "floor_height_nm": floor_height * 1.0e3,
        "floor_height_p95_minus_p05_nm": floor_spread * 1.0e3,
        "top_cd_nm": top,
        "middle_cd_nm": middle,
        "bottom_cd_nm": bottom,
        "sidewall_angle_from_wafer_deg": sidewall_angle,
        "bow_nm": bow,
        "cd_metrics_grid_resolved": bool(resolved_relief),
        "minimum_relief_for_cd_claim_nm": 2.0 * geometry.dx * 1.0e3,
        "cross_section": cross_section,
    }


def _scenario_inputs(preregistration, reactor):
    powered_drop = [
        float(row["powered_electrode_sheath_drop_V"])
        for row in reactor["power_board"]
    ]
    energy = {
        "low": 0.5 * min(powered_drop),
        "high": max(powered_drop),
    }
    tails = preregistration["deterministic_feature_transport"][
        "tail_fraction_sensitivity"
    ]
    return [
        {
            "name": f"ion_{energy_label}_tail_{str(tail).replace('.', 'p')}",
            "impact_energy_label": energy_label,
            "impact_energy_eV": float(energy_eV),
            "tail_fraction": float(tail),
        }
        for energy_label, energy_eV in energy.items()
        for tail in tails
    ]


def _boundary(scenario, preregistration):
    transport = preregistration["deterministic_feature_transport"]
    return build_two_component_boundary(
        kim_2025_reference_iadf(tail_fraction=scenario["tail_fraction"]),
        1.0,
        scenario["impact_energy_eV"],
        name=scenario["name"],
        n_polar=int(transport["polar_quadrature_nodes"]),
        azimuthal_order=int(transport["azimuthal_quadrature_nodes"]),
        reference_plane_m=SOURCE_Z_UM * 1.0e-6,
        extra_provenance={
            "evidence_class": "cross_machine_transport_sensitivity",
            "target_iead_measured": False,
        },
    )


def _mechanism(scenario_name, rate_nm_min, density_kg_m3):
    return RateNormalizedRemovalMechanism(RateNormalizedRemovalParameters(
        material_name="ALD TiO2",
        material_inventory_name="tio2_formula_units",
        projectile_species=(scenario_name,),
        reference_projectile_flux_m2_s=1.0,
        blanket_removal_velocity_m_s=float(rate_nm_min) * 1.0e-9 / 60.0,
        bulk_material_unit_density_m3=tio2_formula_unit_density_m3(density_kg_m3),
        evidence={
            "reference_projectile_flux_m2_s": ParameterEvidence(
                "unit-normalized conserved feature-boundary scenario",
                "calculated",
                supports_prediction_within_declared_domain=True,
            ),
            "blanket_removal_velocity_m_s": ParameterEvidence(
                "Janissen 2016 TiO2/CHF3/O2 cross-machine feature-rate witness",
                "cross_machine_process_analog",
                note="not measured on the Oxford NPG80 target run",
                supports_prediction_within_declared_domain=False,
            ),
            "bulk_material_unit_density_m3": ParameterEvidence(
                "measured cross-process ALD TiO2 density range",
                "measured_cross_process",
                supports_prediction_within_declared_domain=False,
            ),
        },
    ))


def _initial_profile_geometry(*, width_nm, dx_nm, preregistration):
    geometry_board = preregistration["inferred_geometry_board"]
    pitch_nm = float(geometry_board["pitch_nm"])
    layered_geometry = make_square_pillar_mask_geometry_3d(
        pitch=pitch_nm * 1.0e-3,
        domain_height=DOMAIN_HEIGHT_UM,
        dx=dx_nm * 1.0e-3,
        pillar_width=width_nm * 1.0e-3,
        film_thickness=FILM_THICKNESS_UM,
        mask_thickness=MASK_THICKNESS_UM,
        base_top=BASE_TOP_UM,
        mesh_length_unit_m=1.0e-6,
        film_material_id=FILM_MATERIAL,
        mask_material_id=MASK_MATERIAL,
        base_material_id=BASE_MATERIAL,
    )
    # Only TiO2 evolves in this conditional rung; Cr and fused silica are
    # explicitly pinned.  Evolve the external gas/solid union and retain the
    # material labels, rather than independently advecting TiO2's buried
    # solid/solid level-set contacts.  The latter would incorrectly apply a
    # sidewall recession speed to the TiO2/Cr contact within the extension
    # band.  Multi-material level sets remain required when more than one
    # material has a moving law.
    geometry = FeatureGeometry3D(
        layered_geometry.phi,
        layered_geometry.material_id,
        layered_geometry.dx,
        layered_geometry.mesh_length_unit_m,
        layered_geometry.mesh_origin_m,
        material_levelsets=None,
    )
    return geometry, pitch_nm


def _profile_snapshot(
        geometry, *, width_nm, dx_nm, pitch_nm, scenario, target_rate_nm_min,
        reference_rate_nm_min, requested_duration_s, reference_elapsed_s,
        accepted_steps, maximum_balance, maximum_remap, validity,
        clearance_reference_bracket_s=None):
    metrics = _profile_metrics(
        geometry,
        pitch_um=pitch_nm * 1.0e-3,
        nominal_width_um=width_nm * 1.0e-3,
    )
    cleared = clearance_reference_bracket_s is not None
    clearance_process_bracket = None
    if cleared:
        metrics["etched_depth_nm"] = FILM_THICKNESS_UM * 1.0e3
        clearance_process_bracket = [
            float(value) * float(reference_rate_nm_min)
            / float(target_rate_nm_min)
            for value in clearance_reference_bracket_s
        ]
    return {
        "width_nm": float(width_nm),
        "mesh_spacing_nm": float(dx_nm),
        "duration_s": float(requested_duration_s),
        "blanket_rate_nm_min": float(target_rate_nm_min),
        "transport_scenario": dict(scenario),
        "dose_equivalent_reference_rate_nm_min": float(reference_rate_nm_min),
        "dose_equivalent_reference_time_s": float(reference_elapsed_s),
        "dose_equivalence_exact_for_declared_surface_law": True,
        "accepted_profile_steps": int(accepted_steps),
        "accepted_process_equivalent_duration_s": float(
            reference_elapsed_s * reference_rate_nm_min / target_rate_nm_min
        ),
        "tio2_clearance_detected": bool(cleared),
        "clearance_time_bracket_s": clearance_process_bracket,
        "post_clearance_profile_identified": not cleared,
        "profile_geometry_status": (
            "last_pre_clearance_geometry; only depth=film thickness is identified"
            if cleared else "endpoint_geometry"
        ),
        "maximum_transport_relative_particle_balance_error": float(
            maximum_balance
        ),
        "maximum_state_remap_relative_conservation_residual": float(
            maximum_remap
        ),
        "validity": {
            "within_declared_scope": (
                True if validity is None else validity.within_declared_scope
            ),
            "parameter_evidence_supports_prediction": (
                False if validity is None
                else validity.parameter_evidence_supports_prediction
            ),
            "nonpredictive_parameters": (
                [] if validity is None else list(validity.nonpredictive_parameters)
            ),
            "known_limitations": (
                [] if validity is None else list(validity.known_limitations)
            ),
        },
        "profile": metrics,
    }


def _run_profile_dose_trajectory(
        *, width_nm, scenario, target_rates_nm_min, duration_s, dx_nm,
        preregistration):
    """Solve one geometry path and sample all rate-equivalent dose endpoints.

    The declared conditional law is linear in its transferred blanket rate and
    has no time-dependent state.  Therefore geometry depends on the product of
    rate and time.  Running the maximum rate once and stopping exactly at each
    lower-rate dose is the same continuous equation, not a surrogate or an
    interpolation between geometries.
    """
    rates = sorted(set(float(value) for value in target_rates_nm_min))
    if not rates or any(value <= 0.0 for value in rates):
        raise ValueError("target blanket rates must be positive")
    reference_rate = max(rates)
    targets = [
        {
            "rate_nm_min": rate,
            "reference_time_s": float(duration_s) * rate / reference_rate,
        }
        for rate in rates
    ]
    geometry, pitch_nm = _initial_profile_geometry(
        width_nm=width_nm,
        dx_nm=dx_nm,
        preregistration=preregistration,
    )
    boundary = _boundary(scenario, preregistration)
    density = float(np.mean(
        preregistration["surface_response_axes"]["ald_tio2_density_kg_m3"]
    ))
    mechanism = _mechanism(scenario["name"], reference_rate, density)
    state = None
    fingerprint = None
    elapsed = 0.0
    maximum_step_s = 8.0
    accepted_steps = 0
    maximum_balance = 0.0
    maximum_remap = 0.0
    validity = None
    profiles = []
    target_index = 0
    while target_index < len(targets):
        target = targets[target_index]
        remaining = float(target["reference_time_s"]) - elapsed
        if remaining <= 64.0 * np.finfo(float).eps * max(1.0, elapsed):
            profiles.append(_profile_snapshot(
                geometry,
                width_nm=width_nm,
                dx_nm=dx_nm,
                pitch_nm=pitch_nm,
                scenario=scenario,
                target_rate_nm_min=target["rate_nm_min"],
                reference_rate_nm_min=reference_rate,
                requested_duration_s=duration_s,
                reference_elapsed_s=elapsed,
                accepted_steps=accepted_steps,
                maximum_balance=maximum_balance,
                maximum_remap=maximum_remap,
                validity=validity,
            ))
            target_index += 1
            continue
        step_duration = min(maximum_step_s, remaining)
        try:
            step = advance_feature_step_3d(
                geometry,
                boundary,
                {scenario["name"]: "energetic_bombardment"},
                mechanism,
                etchable_material_ids=(FILM_MATERIAL,),
                duration_s=step_duration,
                source_bounds=(
                    0.0,
                    pitch_nm * 1.0e-3,
                    0.0,
                    pitch_nm * 1.0e-3,
                ),
                source_z=SOURCE_Z_UM,
                surface_state=state,
                surface_state_mesh_fingerprint=fingerprint,
                ballistic_transport="face_gather",
                ballistic_periodic_lateral=True,
                ballistic_face_quadrature_points=int(
                    preregistration["deterministic_feature_transport"][
                        "triangle_quadrature_points"
                    ]
                ),
                profile_periodic_lateral=True,
                topology_change_policy="continue_gas_cavity",
                surface_state_remap_backend="indexed_knn",
                reinitialization_method="cr2",
                transport_device="cpu",
            )
        except SurfaceTopologyChangeError as error:
            if error.event_kind != "domain_gas_breakthrough":
                raise
            clearance_bracket = [float(elapsed), float(elapsed + step_duration)]
            for unresolved in targets[target_index:]:
                profiles.append(_profile_snapshot(
                    geometry,
                    width_nm=width_nm,
                    dx_nm=dx_nm,
                    pitch_nm=pitch_nm,
                    scenario=scenario,
                    target_rate_nm_min=unresolved["rate_nm_min"],
                    reference_rate_nm_min=reference_rate,
                    requested_duration_s=duration_s,
                    reference_elapsed_s=elapsed,
                    accepted_steps=accepted_steps,
                    maximum_balance=maximum_balance,
                    maximum_remap=maximum_remap,
                    validity=validity,
                    clearance_reference_bracket_s=clearance_bracket,
                ))
            break
        geometry = step.geometry
        state = step.next_surface_state
        fingerprint = step.next_surface_state_mesh_fingerprint
        elapsed += step_duration
        accepted_steps += 1
        maximum_balance = max(
            maximum_balance,
            max(
                abs(
                    float(step.transport.hit_probability[name])
                    + float(step.transport.escape_probability[name])
                    + float(step.transport.truncation_probability[name])
                    - 1.0
                )
                for name in step.transport.hit_probability
            ),
        )
        maximum_remap = max(
            maximum_remap,
            _maximum_named_numeric(
                step.state_remap_diagnostics,
                "max_relative_conservation_residual",
            ),
        )
        validity = step.validity
    return profiles


def _run_profile_trajectory(*, width_nm, scenario, rate_nm_min, duration_s,
                            dx_nm, preregistration):
    """Compatibility wrapper for one conditional rate endpoint."""
    return _run_profile_dose_trajectory(
        width_nm=width_nm,
        scenario=scenario,
        target_rates_nm_min=(rate_nm_min,),
        duration_s=duration_s,
        dx_nm=dx_nm,
        preregistration=preregistration,
    )[0]


def _build(*, smoke=False):
    preregistration = _load(PREREGISTRATION)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_DOSE)
    scenarios = _scenario_inputs(preregistration, reactor)
    rates = [
        float(analog["source_feature_depth_board"]["minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"]["maximum_implied_rate_nm_min"]),
    ]
    if smoke:
        jobs = [(200.0, scenarios[0], (rates[0],), 60.0, MESH_SPACING_NM)]
    else:
        widths = tuple(float(value) for value in preregistration[
            "inferred_geometry_board"]["width_nm"])
        # Every width and every preregistered IADF condition receives both
        # cross-machine rate endpoints.  Rate homogeneity reduces the 56
        # reported profiles to 28 independent deterministic trajectories.
        jobs = [
            (width, scenario, tuple(rates), 1200.0, MESH_SPACING_NM)
            for width in widths
            for scenario in scenarios
        ]
    def execute(job):
        width, scenario, target_rates, duration, dx = job
        return _run_profile_dose_trajectory(
            width_nm=width,
            scenario=scenario,
            target_rates_nm_min=target_rates,
            duration_s=duration,
            dx_nm=dx,
            preregistration=preregistration,
        )

    # Smoke stays single-process for fast error traces. Production jobs are
    # independent deterministic periodic cells and therefore embarrassingly
    # parallel. ``map`` preserves the preregistered job order.
    if smoke:
        profile_groups = [execute(jobs[0])]
        cache_receipts = []
    else:
        profile_groups = [None] * len(jobs)
        cache_paths = [None] * len(jobs)
        missing = []
        for index, job in enumerate(jobs):
            cached, path, spec = _load_cached_trajectory(job, preregistration)
            cache_paths[index] = path
            if cached is None:
                missing.append((index, job, preregistration, spec, path))
            else:
                profile_groups[index] = cached
        worker_count = min(4, len(missing))
        if missing:
            with ProcessPoolExecutor(max_workers=worker_count) as pool:
                computed = pool.map(_execute_profile_job, missing)
                for (index, *_), group in zip(missing, computed):
                    profile_groups[index] = group
        if any(group is None for group in profile_groups):
            raise RuntimeError("conditional profile trajectory board is incomplete")
        cache_receipts = [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": _hash(path),
            }
            for path in cache_paths
        ]
    profiles = [profile for group in profile_groups for profile in group]
    return {
        "schema": "petch.zhu-npg80-conditional-profile-board.v1",
        "condition_id": preregistration["condition_id"],
        "smoke_only": bool(smoke),
        "target_sem_used": False,
        "target_depth_used": False,
        "coefficient_selected_from_target": None,
        "surface_scale_status": (
            "cross-machine conditional analog; not an Oxford absolute-depth prediction"
        ),
        "trajectory_factorization": {
            "governing_invariance": "geometry depends on blanket_rate_times_time",
            "exact_for_declared_rate_normalized_law": True,
            "independent_trajectories": len(jobs),
            "reported_profile_endpoints": len(profiles),
        },
        "geometry_evolution_mode": (
            "single evolving TiO2 gas-solid union with pinned Cr and fused silica"
        ),
        "mask_handling": (
            "45 nm Cr geometry is pinned during profile evolution; the separate blind board "
            "reports selectivity-conditioned exhaustion and forbids post-exhaustion claims"
        ),
        "inputs": {
            "preregistration_sha256": _hash(PREREGISTRATION),
            "analog_board_sha256": _hash(ANALOG_BOARD),
            "reactor_dose_sha256": _hash(REACTOR_DOSE),
        },
        "trajectory_cache_receipts": cache_receipts,
        "profiles": profiles,
    }


def _execute_profile_job(payload):
    index, job, preregistration, spec, path = payload
    width, scenario, target_rates, duration, dx = job
    profiles = _run_profile_dose_trajectory(
        width_nm=width,
        scenario=scenario,
        target_rates_nm_min=target_rates,
        duration_s=duration,
        dx_nm=dx,
        preregistration=preregistration,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_render({
        "schema": "petch.zhu-npg80-conditional-profile-trajectory.v1",
        "job_index": int(index),
        "job_spec": spec,
        "profiles": profiles,
    }), encoding="utf-8")
    return profiles


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = _build(smoke=args.smoke)
    rendered = _render(payload)
    if args.smoke:
        print(rendered, end="")
        return
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("conditional profile audit is stale")
        print(f"PASS {OUTPUT.relative_to(ROOT)}")
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
