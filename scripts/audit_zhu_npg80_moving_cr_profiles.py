#!/usr/bin/env python3
"""Evolve the blind Oxford TiO2 board with a moving Cr mask.

This is still a cross-machine conditional board.  It upgrades the former
pinned-mask geometry by routing TiO2 and Cr through independent dimensional
removal laws on the same deterministic transport mesh.  The TiO2 rate and
TiO2:Cr selectivity are the preregistered Janissen witness intervals; neither
is fitted to Freddie's withheld SEM.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
import multiprocessing
import os
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.feature_step_3d import (
    SurfaceTopologyChangeError,
    advance_feature_step_3d,
    make_square_pillar_mask_geometry_3d,
)
from petch.material_mechanism_3d import MaterialMechanismRouter3D
from petch.rate_normalized_removal import (
    RateNormalizedRemovalMechanism,
    RateNormalizedRemovalParameters,
)
from petch.surface_kinetics import ParameterEvidence
from petch.tio2_ion_dose import AVOGADRO_PER_MOL, tio2_formula_unit_density_m3

from scripts.audit_zhu_npg80_conditional_profiles import (
    BASE_MATERIAL,
    BASE_TOP_UM,
    DOMAIN_HEIGHT_UM,
    FILM_MATERIAL,
    FILM_THICKNESS_UM,
    FILM_TOP_UM,
    MASK_MATERIAL,
    MASK_THICKNESS_UM,
    SOURCE_Z_UM,
    _boundary,
    _hash,
    _load,
    _maximum_named_numeric,
    _profile_metrics,
    _scenario_inputs,
)


DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
PREREGISTRATION = DATA / "square_pillar_blind_preregistration.json"
ANALOG_BOARD = DATA / "janissen_tio2_analog_board.json"
REACTOR_DOSE = (
    ROOT / "results" / "curated" / "zhu_npg80_daughter_wafer_dose_v1"
    / "audit.json"
)
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_moving_cr_profiles_v1"
    / "audit.json"
)
CACHE_DIR = OUTPUT.parent / "trajectories"
# v5: the common engine continues through certified final Cr extinction.  A
# mask layer is retired only when it owns no material, supports no resolved
# volume cell, and has at most roundoff-scale positive level-set residue.
# Surviving TiO2 state is conservatively remapped, freshly exposed TiO2 is
# initialized by exact common refinement, and retired Cr inventory is receipted.
#
# v4: periodic strip projection is uniform in flux density, so incident rate
# is distributed in proportion to triangle area.  v3 divided equal event rate
# among strip triangles; an unequal-area marching-cubes sliver could therefore
# amplify its local removal velocity by millions and request millions of CFL
# substeps.  The corrected projection changes the transport operator, so no v3
# trajectory cache is eligible for reuse.
#
# v3: trajectories end as a declared physical event when the Cr mask thins
# below one vertical cell anywhere over the footprint interior.  v2 checked
# only the centre column; the mask thins non-uniformly, so edge columns went
# geometrically degenerate first and the marching-cubes topology certifier
# correctly refused the mesh ("unmatched interior edges").  The certifier is
# untouched; the trajectory now stops before the mesh can degenerate.
MODEL_REVISION = "two-material-moving-tio2-cr-dose-factorization-v5"
PRODUCTION_MESH_SPACING_NM = 10.0
CHROMIUM_MOLAR_MASS_KG_MOL = 51.9961e-3
CHROMIUM_REFERENCE_DENSITY_KG_M3 = 7190.0


def chromium_atom_density_m3(
    density_kg_m3: float = CHROMIUM_REFERENCE_DENSITY_KG_M3,
) -> float:
    density = float(density_kg_m3)
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("chromium density must be positive and finite")
    return density / CHROMIUM_MOLAR_MASS_KG_MOL * AVOGADRO_PER_MOL


def _render(payload) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _geometry(*, width_nm, dx_nm, preregistration):
    pitch_nm = float(preregistration["inferred_geometry_board"]["pitch_nm"])
    geometry = make_square_pillar_mask_geometry_3d(
        pitch=pitch_nm * 1.0e-3,
        domain_height=DOMAIN_HEIGHT_UM,
        dx=dx_nm * 1.0e-3,
        pillar_width=float(width_nm) * 1.0e-3,
        film_thickness=FILM_THICKNESS_UM,
        mask_thickness=MASK_THICKNESS_UM,
        base_top=BASE_TOP_UM,
        mesh_length_unit_m=1.0e-6,
        film_material_id=FILM_MATERIAL,
        mask_material_id=MASK_MATERIAL,
        base_material_id=BASE_MATERIAL,
    )
    return geometry, pitch_nm


def _rate_law(
        *, material_name, inventory, species, rate_nm_min, number_density_m3,
        rate_source, density_source):
    return RateNormalizedRemovalMechanism(RateNormalizedRemovalParameters(
        material_name=material_name,
        material_inventory_name=inventory,
        projectile_species=(species,),
        reference_projectile_flux_m2_s=1.0,
        blanket_removal_velocity_m_s=float(rate_nm_min) * 1.0e-9 / 60.0,
        bulk_material_unit_density_m3=float(number_density_m3),
        evidence={
            "reference_projectile_flux_m2_s": ParameterEvidence(
                "unit-normalized deterministic feature-boundary scenario",
                "calculated",
                supports_prediction_within_declared_domain=True,
            ),
            "blanket_removal_velocity_m_s": ParameterEvidence(
                rate_source,
                "cross_machine_process_analog",
                supports_prediction_within_declared_domain=False,
            ),
            "bulk_material_unit_density_m3": ParameterEvidence(
                density_source,
                "reference_material_constant_or_cross_process_range",
                supports_prediction_within_declared_domain=False,
            ),
        },
    ))


def _router(*, scenario_name, tio2_rate_nm_min, selectivity, density_kg_m3):
    cr_rate = float(tio2_rate_nm_min) / float(selectivity)
    tio2 = _rate_law(
        material_name="ALD TiO2",
        inventory="TiO2_formula_unit",
        species=scenario_name,
        rate_nm_min=tio2_rate_nm_min,
        number_density_m3=tio2_formula_unit_density_m3(density_kg_m3),
        rate_source="Janissen 2016 TiO2 feature-rate witness",
        density_source="measured cross-process ALD TiO2 density range",
    )
    chromium = _rate_law(
        material_name="Cr hard mask",
        inventory="Cr_atom",
        species=scenario_name,
        rate_nm_min=cr_rate,
        number_density_m3=chromium_atom_density_m3(),
        rate_source=(
            "Janissen 2016 cross-machine TiO2:Cr selectivity divided into "
            "the paired TiO2 rate"
        ),
        density_source="7190 kg/m3 bulk-Cr reference density",
    )
    return MaterialMechanismRouter3D(
        {FILM_MATERIAL: tio2, MASK_MATERIAL: chromium},
        provenance={
            FILM_MATERIAL: {
                "source": "janissen-2016-tio2-rie",
                "claim": "conditional TiO2 rate",
            },
            MASK_MATERIAL: {
                "source": "janissen-2016-tio2-rie",
                "claim": "conditional paired Cr rate from selectivity",
                "nguyen_2021_topology_warning": (
                    "nguyen-2021-cr-sf6-o2: "
                    "real Cr removal in F/O plasma resolves CrOx, CrFx, neutral "
                    "conversion, and ion-assisted inhibitor removal"
                ),
            },
        },
    )


def _zero_crossings(line, coordinate):
    line = np.asarray(line, dtype=float)
    crossing = []
    for index in range(len(line) - 1):
        if line[index] == 0.0:
            crossing.append(float(coordinate[index]))
        if line[index] * line[index + 1] < 0.0:
            fraction = -line[index] / (line[index + 1] - line[index])
            crossing.append(float(
                coordinate[index] + fraction * (coordinate[index + 1] - coordinate[index])
            ))
    return sorted(crossing)


def _mask_metrics(geometry, *, pitch_nm):
    if (geometry.material_levelsets is None
            or MASK_MATERIAL not in geometry.material_levelsets):
        return {
            "center_remaining_thickness_nm": 0.0,
            "center_top_height_nm": None,
            "center_bottom_height_nm": None,
            "vertical_cells_remaining": 0.0,
            "mask_exhausted_at_center": True,
            "mask_below_vertical_resolution_at_center": True,
            "material_layer_retired": True,
        }
    mask = np.asarray(geometry.material_levelsets[MASK_MATERIAL], dtype=float)
    coordinate = np.arange(mask.shape[0]) * geometry.dx
    z = np.arange(mask.shape[2]) * geometry.dx
    center = 0.5 * float(pitch_nm) * 1.0e-3
    middle = int(np.argmin(np.abs(coordinate - center)))
    crossing = _zero_crossings(mask[middle, middle, :], z)
    if len(crossing) < 2 or not np.any(mask[middle, middle, :] >= 0.0):
        return {
            "center_remaining_thickness_nm": 0.0,
            "center_top_height_nm": None,
            "center_bottom_height_nm": None,
            "vertical_cells_remaining": 0.0,
            "mask_exhausted_at_center": True,
            "mask_below_vertical_resolution_at_center": True,
            "material_layer_retired": False,
        }
    bottom, top = crossing[-2], crossing[-1]
    thickness = max(0.0, top - bottom)
    return {
        "center_remaining_thickness_nm": float(thickness * 1.0e3),
        "center_top_height_nm": float(top * 1.0e3),
        "center_bottom_height_nm": float(bottom * 1.0e3),
        "vertical_cells_remaining": float(thickness / geometry.dx),
        "mask_exhausted_at_center": False,
        "mask_below_vertical_resolution_at_center": bool(
            thickness < geometry.dx
        ),
        "material_layer_retired": False,
    }


def _mask_interior_minimum_thickness_nm(geometry, *, pitch_nm, width_nm):
    """Minimum remaining Cr thickness over the original mask footprint
    interior (eroded by one cell so corner rounding is excluded).

    Sub-cell thickness anywhere inside the footprint means the next level-set
    motion can produce a non-manifold extracted surface; the trajectory must
    end as a declared mask-exhaustion event before that happens.
    """
    if (geometry.material_levelsets is None
            or MASK_MATERIAL not in geometry.material_levelsets):
        return 0.0
    mask = np.asarray(geometry.material_levelsets[MASK_MATERIAL], dtype=float)
    dx = float(geometry.dx)
    coordinate = np.arange(mask.shape[0]) * dx
    center = 0.5 * float(pitch_nm) * 1.0e-3
    half = 0.5 * float(width_nm) * 1.0e-3 - dx
    if half <= 0.0:
        selected = np.zeros(mask.shape[0], dtype=bool)
        selected[int(np.argmin(np.abs(coordinate - center)))] = True
    else:
        selected = np.abs(coordinate - center) <= half + 1.0e-12
    columns = mask[np.ix_(selected, selected)]
    positive = columns >= 0.0
    a = columns[:, :, :-1]
    b = columns[:, :, 1:]
    both = positive[:, :, :-1] & positive[:, :, 1:]
    change = positive[:, :, :-1] ^ positive[:, :, 1:]
    denominator = np.where(change, np.abs(a - b), 1.0)
    segment = np.where(
        both, 1.0,
        np.where(
            change,
            np.maximum(a, b) / np.maximum(denominator, 1.0e-300),
            0.0,
        ),
    )
    return float(np.min(segment.sum(axis=2)) * dx * 1.0e3)


def _snapshot(
        geometry, *, width_nm, dx_nm, pitch_nm, scenario, rate_nm_min,
        selectivity, requested_duration_s, reference_rate_nm_min,
        reference_elapsed_s, accepted_steps, maximum_balance, maximum_remap,
        validity, terminal_reason, material_lifecycle_events):
    return {
        "width_nm": float(width_nm),
        "mesh_spacing_nm": float(dx_nm),
        "duration_s": float(requested_duration_s),
        "blanket_tio2_rate_nm_min": float(rate_nm_min),
        "tio2_to_cr_selectivity": float(selectivity),
        "blanket_cr_rate_nm_min": float(rate_nm_min / selectivity),
        "transport_scenario": dict(scenario),
        "dose_equivalent_reference_rate_nm_min": float(reference_rate_nm_min),
        "dose_equivalent_reference_time_s": float(reference_elapsed_s),
        "accepted_process_equivalent_duration_s": float(
            reference_elapsed_s * reference_rate_nm_min / rate_nm_min
        ),
        "accepted_profile_steps": int(accepted_steps),
        "terminal_reason": terminal_reason,
        "material_lifecycle_events": tuple(material_lifecycle_events),
        "profile": _profile_metrics(
            geometry,
            pitch_um=float(pitch_nm) * 1.0e-3,
            nominal_width_um=float(width_nm) * 1.0e-3,
        ),
        "cr_mask": _mask_metrics(geometry, pitch_nm=pitch_nm),
        "cr_mask_interior_minimum_thickness_nm": (
            _mask_interior_minimum_thickness_nm(
                geometry, pitch_nm=pitch_nm, width_nm=width_nm)),
        "maximum_transport_relative_particle_balance_error": float(maximum_balance),
        "maximum_state_remap_relative_conservation_residual": float(maximum_remap),
        "parameter_evidence_supports_prediction": (
            False if validity is None
            else bool(validity.parameter_evidence_supports_prediction)
        ),
    }


def _run_trajectory(
        *, width_nm, scenario, rates_nm_min, selectivity, duration_s, dx_nm,
        preregistration, maximum_step_s=None, transport_device="cpu"):
    rates = sorted(set(float(value) for value in rates_nm_min))
    reference_rate = max(rates)
    targets = [{
        "rate": rate,
        "reference_time": float(duration_s) * rate / reference_rate,
    } for rate in rates]
    geometry, pitch_nm = _geometry(
        width_nm=width_nm, dx_nm=dx_nm, preregistration=preregistration
    )
    boundary = _boundary(scenario, preregistration)
    density = float(np.mean(
        preregistration["surface_response_axes"]["ald_tio2_density_kg_m3"]
    ))
    router = _router(
        scenario_name=scenario["name"],
        tio2_rate_nm_min=reference_rate,
        selectivity=selectivity,
        density_kg_m3=density,
    )
    if maximum_step_s is None:
        maximum_step_s = 0.35 * float(dx_nm) / (reference_rate / 60.0)
    state = None
    fingerprint = None
    elapsed = 0.0
    accepted_steps = 0
    maximum_balance = 0.0
    maximum_remap = 0.0
    validity = None
    material_lifecycle_events = []
    profiles = []
    target_index = 0
    terminal_reason = "requested_duration"
    while target_index < len(targets):
        target = targets[target_index]
        remaining = target["reference_time"] - elapsed
        if remaining <= 64.0 * np.finfo(float).eps * max(1.0, elapsed):
            profiles.append(_snapshot(
                geometry, width_nm=width_nm, dx_nm=dx_nm, pitch_nm=pitch_nm,
                scenario=scenario, rate_nm_min=target["rate"],
                selectivity=selectivity, requested_duration_s=duration_s,
                reference_rate_nm_min=reference_rate,
                reference_elapsed_s=elapsed, accepted_steps=accepted_steps,
                maximum_balance=maximum_balance, maximum_remap=maximum_remap,
                validity=validity, terminal_reason=terminal_reason,
                material_lifecycle_events=material_lifecycle_events,
            ))
            target_index += 1
            continue
        interior_minimum_nm = _mask_interior_minimum_thickness_nm(
            geometry, pitch_nm=pitch_nm, width_nm=width_nm)
        step_duration = min(float(maximum_step_s), remaining)
        try:
            step = advance_feature_step_3d(
                geometry,
                boundary,
                {scenario["name"]: "energetic_bombardment"},
                router,
                etchable_material_ids=(FILM_MATERIAL, MASK_MATERIAL),
                duration_s=step_duration,
                source_bounds=(0.0, pitch_nm * 1.0e-3, 0.0, pitch_nm * 1.0e-3),
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
                topology_change_policy=(
                    "continue_gas_cavity_and_material_extinction"),
                surface_state_remap_backend="indexed_knn",
                reinitialization_method="cr2",
                transport_device=transport_device,
            )
        except RuntimeError as error:
            raise RuntimeError(
                "moving-mask step failed at reference elapsed "
                f"{elapsed:.3f} s after {accepted_steps} accepted steps; "
                f"cr interior minimum {interior_minimum_nm:.3f} nm; "
                f"cr centre metrics {_mask_metrics(geometry, pitch_nm=pitch_nm)}"
            ) from error
        except SurfaceTopologyChangeError as error:
            # A feature-floor breakthrough is an interpretable physical end
            # point.  Solid or material-component changes can instead expose a
            # discretization/remap defect, so never relabel them as successful
            # process completion here.
            if error.event_kind != "domain_gas_breakthrough":
                raise
            terminal_reason = error.event_kind
            for unresolved in targets[target_index:]:
                profiles.append(_snapshot(
                    geometry, width_nm=width_nm, dx_nm=dx_nm,
                    pitch_nm=pitch_nm, scenario=scenario,
                    rate_nm_min=unresolved["rate"], selectivity=selectivity,
                    requested_duration_s=duration_s,
                    reference_rate_nm_min=reference_rate,
                    reference_elapsed_s=elapsed,
                    accepted_steps=accepted_steps,
                    maximum_balance=maximum_balance,
                    maximum_remap=maximum_remap,
                    validity=validity, terminal_reason=terminal_reason,
                    material_lifecycle_events=material_lifecycle_events,
                ))
            break
        topology_event = step.diagnostics["topology_event"]
        if topology_event is not None:
            material_lifecycle_events.append({
                "kind": str(topology_event["kind"]),
                "reference_elapsed_s": float(elapsed + step_duration),
                "retired_material_ids": tuple(
                    topology_event.get("retired_material_ids", ())),
                "material_extinction_geometry": dict(
                    topology_event.get("material_extinction_geometry", {})),
                "state_retirement": dict(
                    step.state_remap_diagnostics.get("retired_materials", {})),
                "effective_surface_state_remap_backend": str(
                    step.state_remap_diagnostics[
                        "surface_state_remap_backend"]),
            })
        geometry = step.geometry
        state = step.next_surface_state
        fingerprint = step.next_surface_state_mesh_fingerprint
        elapsed += step_duration
        accepted_steps += 1
        validity = step.validity
        maximum_balance = max(maximum_balance, max(
            abs(float(step.transport.hit_probability[name])
                + float(step.transport.escape_probability[name])
                + float(step.transport.truncation_probability[name]) - 1.0)
            for name in step.transport.hit_probability
        ))
        maximum_remap = max(maximum_remap, _maximum_named_numeric(
            step.state_remap_diagnostics,
            "max_relative_conservation_residual",
        ))
    return profiles


def _job_spec(job):
    width, scenario, rates, selectivity, duration, dx = job
    return {
        "model_revision": MODEL_REVISION,
        "preregistration_sha256": _hash(PREREGISTRATION),
        "width_nm": float(width),
        "scenario": scenario,
        "rates_nm_min": list(rates),
        "selectivity": float(selectivity),
        "duration_s": float(duration),
        "mesh_spacing_nm": float(dx),
    }


def _cache_path(spec):
    digest = sha256(_render(spec).encode("utf-8")).hexdigest()[:16]
    return CACHE_DIR / (
        f"w{int(round(spec['width_nm'])):03d}_s{spec['selectivity']:.3f}_"
        f"{spec['scenario']['name']}_{digest}.json"
    )


def _execute(payload):
    job, transport_device = payload
    preregistration = _load(PREREGISTRATION)
    width, scenario, rates, selectivity, duration, dx = job
    spec = _job_spec(job)
    try:
        profiles = _run_trajectory(
            width_nm=width, scenario=scenario, rates_nm_min=rates,
            selectivity=selectivity, duration_s=duration, dx_nm=dx,
            preregistration=preregistration,
            transport_device=transport_device,
        )
    except Exception as error:
        raise RuntimeError(
            "moving-mask trajectory failed for "
            f"job_spec={json.dumps(spec, sort_keys=True)}"
        ) from error
    path = _cache_path(spec)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_render({
        "job_spec": spec,
        "execution": {"transport_device": str(transport_device)},
        "profiles": profiles,
    }), encoding="utf-8")
    return profiles, path, str(transport_device)


def _process_pool_options(transport_device):
    """Return a fork-safe executor configuration for the selected device."""
    if str(transport_device).startswith("cuda"):
        return {"mp_context": multiprocessing.get_context("spawn")}
    return {}


def build(*, smoke=False, transport_device="cpu", workers=None):
    preregistration = _load(PREREGISTRATION)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_DOSE)
    scenarios = _scenario_inputs(preregistration, reactor)
    rates = (
        float(analog["source_feature_depth_board"]["minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"]["maximum_implied_rate_nm_min"]),
    )
    selectivities = tuple(float(value) for value in
                          preregistration["surface_response_axes"][
                              "tio2_to_cr_selectivity"])
    if smoke:
        jobs = [(200.0, scenarios[0], (rates[0],), selectivities[0], 12.0, 20.0)]
    else:
        widths = tuple(float(value) for value in
                       preregistration["inferred_geometry_board"]["width_nm"])
        jobs = [
            (width, scenario, rates, selectivity, 1200.0, PRODUCTION_MESH_SPACING_NM)
            for width in widths for scenario in scenarios
            for selectivity in selectivities
        ]
    groups = []
    receipts = []
    missing = []
    for job in jobs:
        spec = _job_spec(job)
        path = _cache_path(spec)
        if path.exists():
            cached = _load(path)
            if cached.get("job_spec") != spec:
                raise RuntimeError(f"moving-mask cache mismatch: {path}")
            groups.append(cached["profiles"])
            receipts.append({
                "path": str(path.relative_to(ROOT)),
                "sha256": _hash(path),
                "transport_device": cached.get("execution", {}).get(
                    "transport_device", "unrecorded"
                ),
            })
        else:
            groups.append(None)
            receipts.append(None)
            missing.append((len(groups) - 1, job))
    if missing:
        if smoke:
            computed = [_execute((missing[0][1], transport_device))]
        else:
            if workers is None:
                workers = 1 if str(transport_device).startswith("cuda") else 4
            workers = int(workers)
            if workers <= 0:
                raise ValueError("workers must be a positive integer")
            payloads = [
                (item[1], str(transport_device)) for item in missing
            ]
            if workers == 1:
                computed = [_execute(payload) for payload in payloads]
            else:
                # CUDA primary contexts are not fork-safe.  A forked worker
                # inherits Warp's initialized runtime but cannot acquire the
                # parent's device context.  Spawn gives every deterministic
                # trajectory worker its own clean CUDA context.
                pool_options = _process_pool_options(transport_device)
                with ProcessPoolExecutor(
                    max_workers=min(workers, len(missing)),
                    **pool_options,
                ) as pool:
                    computed = list(pool.map(_execute, payloads))
        for (index, _), (profiles, path, device) in zip(missing, computed):
            groups[index] = profiles
            receipts[index] = {
                "path": str(path.relative_to(ROOT)),
                "sha256": _hash(path),
                "transport_device": str(device),
            }
    profiles = [profile for group in groups for profile in group]
    return {
        "schema": "petch.zhu-npg80-moving-cr-profile-board.v1",
        "condition_id": preregistration["condition_id"],
        "smoke_only": bool(smoke),
        "target_sem_used": False,
        "target_depth_used": False,
        "coefficient_selected_from_target": None,
        "mesh_spacing_nm": 20.0 if smoke else PRODUCTION_MESH_SPACING_NM,
        "moving_materials": ["ALD TiO2", "Cr hard mask"],
        "pinned_materials": ["fused-silica substrate"],
        "execution": {
            "trajectory_transport_devices": sorted(set(
                receipt["transport_device"] for receipt in receipts
            )),
            "execution_device_not_part_of_physics_spec": True,
            "cross_device_numerical_parity_required_before_combining": True,
        },
        "conditional_axes": {
            "tio2_rate_nm_min": list(rates),
            "tio2_to_cr_selectivity": list(selectivities),
            "transport_scenarios": scenarios,
        },
        "model_boundary": {
            "rate_normalized_not_microscopic": True,
            "cr_topology_source": "nguyen-2021-cr-sf6-o2",
            "nguyen_2021_cr_oxide_fluoride_state_unresolved": True,
            "supports_absolute_oxford_profile_prediction": False,
        },
        "trajectory_receipts": receipts,
        "profiles": profiles,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument(
        "--transport-device",
        default=os.environ.get("PETCH_TRANSPORT_DEVICE", "cpu"),
        help="deterministic Warp transport device, for example cpu or cuda:0",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=(
            int(os.environ["PETCH_PROFILE_WORKERS"])
            if "PETCH_PROFILE_WORKERS" in os.environ else None
        ),
        help="independent trajectory workers (default 1 on CUDA, 4 on CPU)",
    )
    args = parser.parse_args()
    if args.smoke:
        print(_render(build(
            smoke=True,
            transport_device=args.transport_device,
            workers=args.workers,
        )))
        return
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build(
        smoke=False,
        transport_device=args.transport_device,
        workers=args.workers,
    ))
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
        return
    if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("moving-Cr profile audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
