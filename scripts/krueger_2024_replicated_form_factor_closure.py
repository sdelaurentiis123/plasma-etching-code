#!/usr/bin/env python3
"""Bounded replicated form-factor closure on the sealed Krueger base checkpoint.

This audit answers one numerical question before another profile is evolved: does the
finite-area, certified-hard-visibility diffuse-neutral operator have enough randomized-QMC
precision at 8 and 16 rays per face for the shortest already-declared chemistry horizon?

The default is *Stage A only*.  It captures direct transport once, constructs eight independent
scrambled-Sobol form-factor replicates at each nested ray level, and scores the resulting
start-state chemistry response on integrated material observables and 20/40 nm physical patches.
It never reads a held-out Krueger observation and it never advances the profile.

Stage B is an explicitly authorized, resumable 22-job confirmation using the unified cached
surface-radiosity integrator.  It cannot be entered by a Stage-A run and requires a completed,
passing Stage-A artifact plus ``--authorize-stage-b``.  Every expensive operation has a hard
wall-time cap; failure at 16 rays writes ``bounded_precision_hold`` rather than silently escalating
to 32 rays or opening held-out data.
"""
from __future__ import annotations

import argparse
from collections.abc import Mapping
from dataclasses import asdict, is_dataclass
from hashlib import sha256
import json
import multiprocessing as mp
from pathlib import Path
import sys
from time import perf_counter
from unittest.mock import patch
from zipfile import BadZipFile

import numpy as np
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT / "src", ROOT / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from krueger_2024_endpoint_operator_audit import _evaluate  # noqa: E402
from krueger_2024_frozen_checkpoint_2x2 import (  # noqa: E402
    EvaluationDeadlineExceeded,
    _hard_deadline,
    _hash_manifest,
    _inputs_unchanged,
    _operator_config,
    _sha256,
    _snapshot_inputs,
    _write_json_atomic,
)
from krueger_2024_frozen_surface_chemistry import (  # noqa: E402
    BASE_INPUT_PATHS,
    PARAMETERS,
    surface_flux_sha256,
)
from krueger_2024_trench_pilot import (  # noqa: E402
    _load_checkpoint,
    measure_krueger_metrics,
)
from petch.amorphous_carbon_mask import (  # noqa: E402
    build_krueger_2024_material_router_3d,
)
from petch.boundary_transport_3d import (  # noqa: E402
    estimate_diffuse_form_factors_3d,
)
from petch.diffuse_form_factor_control_3d import (  # noqa: E402
    ReplicatedDiffuseFormFactors3D,
)
import petch.feature_step_3d as feature_step_module  # noqa: E402
from petch.feature_step_3d import (  # noqa: E402
    _face_material_ids,
    _surface_gas_normals,
)
from petch.surface_patch_convergence_3d import (  # noqa: E402
    DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION,
    aggregate_surface_field_on_physical_patches_3d,
    score_surface_field_refinement_at_physical_scales_3d,
    score_replicated_surface_field_at_physical_scales_3d,
)
from petch.surface_radiosity_coupling_3d import (  # noqa: E402
    SurfaceRadiosityOperatorCache3D,
    integrate_surface_radiosity_chemistry_3d,
)
from petch.neutral_radiosity_3d import DiffuseFormFactors3D  # noqa: E402
from petch.surface_kinetics import (  # noqa: E402
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)
from petch.threed import extract_mesh_3d  # noqa: E402


SCHEMA = "petch.krueger-2024.replicated-form-factor-closure.v3"
DIRECT_TRANSPORT_ARTIFACT_SCHEMA = "petch.direct-transport-artifact.v1"
FORM_FACTOR_ARTIFACT_SCHEMA = "petch.replicated-form-factor-artifact.v2"
STAGE_A_PASS_STATUS = "stage_a_level32_pass_stage_b_held"
REPLICATE_SEEDS = (104729, 130363, 155921, 181081, 205759, 230761, 256019, 280801)
RAY_LEVELS = (8, 16, 32)
DIAGNOSTIC_LEVEL_PAIR = (8, 16)
AUTHORITATIVE_LEVEL_PAIR = (16, 32)
FINAL_STAGE_A_LEVEL = 32
PATCH_SCALES_M = (20.0e-9, 40.0e-9)
PATCH_PERIODIC_AXES = (False, True, False)
PATCH_SUPPORT_SENSITIVITY_THRESHOLDS = (0.05, 0.075, 0.10, 0.25, 0.50, 0.75)
FROZEN_RESPONSE_HORIZON_FRACTION = 1.0 / 1024.0
OPERATOR = {
    "boundary_case": "base",
    "boundary_mode": "angular_8x16",
    "ion_energy_bin_eV": 250.0,
    "ion_angle_bin_deg": 0.25,
    "ion_azimuthal_closure": "axisymmetric_uniform",
    "ion_azimuthal_order": 16,
    "ballistic_transport": "face_gather",
    "ballistic_face_quadrature_points": 3,
    "n_position": 16,
    "transport_device": "cpu",
    "source_sampling": "triangle_area",
    "visibility_mode": "cellwise_certified",
    "periodic_lateral": True,
    "ray_offset_dx": 1.0e-3,
    "maximum_visibility_wraps": 1024,
    "maximum_exact_replay_wraps": 4096,
    "radiosity_relative_tolerance": 1.0e-12,
    "radiosity_maximum_iterations": 2000,
    "nested_rqmc_ray_levels": [8, 16, 32],
    "diagnostic_level_pair": [8, 16],
    "authoritative_level_pair": [16, 32],
    "final_stage_a_level": 32,
    "frozen_response_horizon_fraction_of_next_step": 1.0 / 1024.0,
}
GATES = {
    "minimum_replicate_count": 8,
    "maximum_form_factor_row_closure_error": 5.0e-13,
    "maximum_radiosity_relative_balance_error": 5.0e-12,
    "integrated_relative_tolerance": 0.01,
    "integrated_absolute_displacement_dx": 5.0e-4,
    "patch_relative_tolerance": 0.05,
    "patch_absolute_normalized": 2.5e-3,
    "maximum_embedded_relative_error": 0.01,
    "maximum_local_reaction_probability_change": 0.01,
    "tight_maximum_embedded_relative_error": 0.005,
    "tight_maximum_local_reaction_probability_change": 0.005,
    "minimum_chemistry_substep_s": 1.0e-4,
    "maximum_gross_displacement_dx": 0.05,
}
DEFAULT_BUDGETS = {
    "maximum_direct_transport_wall_s": 120.0,
    "maximum_form_factor_replicate_wall_s": 60.0,
    "maximum_endpoint_job_wall_s": 90.0,
    "maximum_stage_a_total_wall_s": 300.0,
    "maximum_total_wall_s": 300.0,
    "maximum_process_count": 4,
}
SOURCE_PATHS = (
    "scripts/krueger_2024_replicated_form_factor_closure.py",
    "scripts/krueger_2024_endpoint_operator_audit.py",
    "scripts/krueger_2024_frozen_checkpoint_2x2.py",
    "scripts/krueger_2024_frozen_surface_chemistry.py",
    "scripts/krueger_2024_trench_pilot.py",
) + tuple(sorted(
    str(path.relative_to(ROOT))
    for path in (ROOT / "src" / "petch").glob("*.py")
    if not path.name.startswith("._")
))
_EXCHANGE_ATTRIBUTES = {
    "removed": "removed_units_m2",
    "outgoing": "outgoing_units_m2",
    "unresolved": "unresolved_units_m2",
    "deposited": "deposited_units_m2",
}


def _jsonable(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return str(value)


def _canonical_sha256(payload):
    encoded = json.dumps(
        _jsonable(payload), sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode("utf-8")
    return sha256(encoded).hexdigest()


def _validate_budget(name, value, ceiling):
    value = float(value)
    if not np.isfinite(value) or value <= 0.0 or value > float(ceiling):
        raise ValueError(f"{name} must lie in (0, {float(ceiling):g}]")
    return value


def _inspect_sealed_base_source(source):
    """Read only base-run metadata and refuse a transfer case before checkpoint I/O."""
    source = Path(source)
    audit_path = source / "audit.json"
    checkpoint_path = source / "checkpoint.npz"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("status") != "complete":
        raise RuntimeError("Krueger source is not a completed run")
    config = audit.get("configuration", {})
    if config.get("boundary_case") != "base":
        raise ValueError(
            "replicated closure accepts only the sealed Krueger base boundary; "
            "checkpoint was not opened")
    return {
        "source": source,
        "audit_path": audit_path,
        "checkpoint_path": checkpoint_path,
        "audit": audit,
        "config": config,
    }


def _load_sealed_base_source(source):
    metadata = _inspect_sealed_base_source(source)
    geometry, state, fingerprint, checkpoint_metadata = _load_checkpoint(
        metadata["checkpoint_path"])
    geometry_config = metadata["config"].get("geometry", {})
    metrics = measure_krueger_metrics(
        geometry,
        substrate_top_um=float(geometry_config.get("substrate_top_um", 1.8)),
        opening_center_um=0.5 * float(geometry_config.get("cell_width_um", 0.13)),
        opening_width_um=float(geometry_config.get("opening_width_um", 0.09)),
    )
    return {
        **metadata,
        "geometry": geometry,
        "state": state,
        "fingerprint": fingerprint,
        "checkpoint_metadata": checkpoint_metadata,
        "metrics": metrics,
    }


def _validate_parameter_provenance(r17_source, r19_source):
    r17 = _inspect_sealed_base_source(r17_source)
    for label, config in (("r17", r17["config"]), ("r19", r19_source["config"])):
        observed = {name: float(config[name]) for name in PARAMETERS[label]}
        if observed != PARAMETERS[label]:
            raise ValueError(f"{label.upper()} parameters disagree with sealed provenance")
    return {
        "r17_audit_sha256": _sha256(r17["audit_path"]),
        "r19_audit_sha256": _sha256(r19_source["audit_path"]),
        "parameter_pairs": PARAMETERS,
    }


def _production_config(reference, parameters):
    config = _operator_config(reference, parameters)
    config.update(
        neutral_direction_polar_order=8,
        neutral_direction_azimuthal_order=16,
        ion_energy_bin_eV=OPERATOR["ion_energy_bin_eV"],
        ion_angle_bin_deg=OPERATOR["ion_angle_bin_deg"],
        ion_azimuthal_closure=OPERATOR["ion_azimuthal_closure"],
        ion_azimuthal_order=OPERATOR["ion_azimuthal_order"],
        ballistic_transport=OPERATOR["ballistic_transport"],
        ballistic_face_quadrature_points=OPERATOR[
            "ballistic_face_quadrature_points"],
        n_position=OPERATOR["n_position"],
        radiosity_enabled=True,
        radiosity_rays_per_face=RAY_LEVELS[0],
        radiosity_source_sampling=OPERATOR["source_sampling"],
        radiosity_relative_tolerance=OPERATOR[
            "radiosity_relative_tolerance"],
        radiosity_maximum_iterations=OPERATOR[
            "radiosity_maximum_iterations"],
        transport_device=OPERATOR["transport_device"],
        boundary_case="base",
        duration_s=0.0,
    )
    return config


def _capture_direct_transport_once(source, config, *, seed):
    """Capture the exact pre-radiosity transport at the production wrapper boundary."""
    captured_factor = []
    captured_base = []
    production_apply = feature_step_module._apply_diffuse_neutral_transport

    def capture_estimator(*args, **kwargs):
        # The desired object is the transport *before* radiosity.  An exact all-escape factor is a
        # neutral identity for this capture and avoids paying for an unrelated production trace.
        face_count = len(np.asarray(args[1]))
        rays = int(kwargs.get("rays_per_face", RAY_LEVELS[0]))
        captured_factor.append({"rays_per_face": rays, "trace_elided": True})
        return DiffuseFormFactors3D(
            face_count, np.asarray([], dtype=int), np.asarray([], dtype=int),
            np.asarray([], dtype=float), np.ones(face_count), rays)

    def capture_apply(transport, *args, **kwargs):
        captured_base.append(transport)
        return production_apply(transport, *args, **kwargs)

    with patch.object(
            feature_step_module, "estimate_diffuse_form_factors_3d",
            side_effect=capture_estimator), patch.object(
                feature_step_module, "_apply_diffuse_neutral_transport",
                side_effect=capture_apply):
        result, boundary, reported_wall = _evaluate(
            source["geometry"], source["state"], source["fingerprint"],
            boundary_mode=OPERATOR["boundary_mode"],
            ion_bins=(OPERATOR["ion_energy_bin_eV"], OPERATOR["ion_angle_bin_deg"]),
            face_points=OPERATOR["ballistic_face_quadrature_points"],
            pilot_config=config, radiosity_rays=RAY_LEVELS[0], seed=int(seed),
            ballistic_transport=OPERATOR["ballistic_transport"],
            n_position=OPERATOR["n_position"],
            transport_device=OPERATOR["transport_device"])
    if len(captured_base) != 1 or not captured_factor:
        raise RuntimeError("production evaluation did not expose one direct-radiosity boundary")
    verts, faces, centroids, areas = extract_mesh_3d(
        source["geometry"].phi, source["geometry"].dx)
    active = np.asarray(result.active_face_index, dtype=int)
    if (not np.array_equal(active, np.arange(len(faces)))
            or not np.array_equal(np.asarray(result.active_face_area), areas)):
        raise RuntimeError("replicated closure requires every checkpoint face to be active")
    normals = _surface_gas_normals(
        verts, faces, centroids, source["geometry"])
    material = _face_material_ids(centroids, source["geometry"])
    physical_area = areas * source["geometry"].mesh_length_unit_m ** 2
    roles = {
        species.name: (
            "energetic_bombardment" if species.charge_number != 0
            else "neutral_reactant")
        for species in boundary.species
    }
    return {
        "production_result": result,
        "boundary": boundary,
        "reported_wall_s": float(reported_wall),
        "direct_surface_fluxes": captured_base[0].surface_fluxes,
        "verts": verts,
        "faces": faces,
        "centroids": centroids,
        "areas_mesh": areas,
        "face_area_m2": physical_area,
        "gas_normals": normals,
        "face_material_id": material,
        "active_face_index": active,
        "species_role": roles,
        "production_form_factor_call_count": len(captured_factor),
        "production_form_factor_trace_elided": True,
    }


def _array_receipt(value):
    """Return a byte-level receipt for one non-object NumPy payload array."""
    array = np.asarray(value)
    if array.dtype.hasobject:
        raise TypeError("direct-transport artifacts refuse object arrays")
    contiguous = np.ascontiguousarray(array)
    return {
        "shape": list(array.shape),
        "dtype": array.dtype.str,
        "sha256": sha256(contiguous.view(np.uint8).tobytes()).hexdigest(),
    }


def _direct_transport_cache_identity(source, config, *, seed):
    """Bind a reusable direct solve to every input that can change its answer."""
    source_manifest = _hash_manifest(SOURCE_PATHS)
    base_input_manifest = _hash_manifest(BASE_INPUT_PATHS)
    identity = {
        "schema": DIRECT_TRANSPORT_ARTIFACT_SCHEMA,
        "checkpoint_sha256": _sha256(source["checkpoint_path"]),
        "source_audit_sha256": _sha256(source["audit_path"]),
        "source_manifest_sha256": _canonical_sha256(source_manifest),
        "base_input_manifest_sha256": _canonical_sha256(base_input_manifest),
        "operator": _jsonable(OPERATOR),
        "production_config": _jsonable(config),
        "transport_seed": int(seed),
    }
    return _jsonable(identity)


def _direct_geometry_from_checkpoint(source):
    """Reconstruct the cheap mesh contract used to certify a cached direct solve."""
    verts, faces, centroids, areas = extract_mesh_3d(
        source["geometry"].phi, source["geometry"].dx)
    normals = _surface_gas_normals(
        verts, faces, centroids, source["geometry"])
    material = _face_material_ids(centroids, source["geometry"])
    return {
        "verts": np.asarray(verts),
        "faces": np.asarray(faces),
        "centroids": np.asarray(centroids),
        "areas_mesh": np.asarray(areas),
        "face_area_m2": (
            np.asarray(areas)
            * float(source["geometry"].mesh_length_unit_m) ** 2),
        "gas_normals": np.asarray(normals),
        "face_material_id": np.asarray(material),
        "active_face_index": np.arange(len(faces), dtype=int),
    }


def _surface_flux_artifact_payload(fluxes):
    """Flatten both supported energetic representations without losing events."""
    arrays = {}
    neutral = []
    for index, (name, value) in enumerate(sorted(fluxes.neutral_flux_m2_s.items())):
        key = f"neutral_{index:03d}_flux_m2_s"
        arrays[key] = np.asarray(value)
        neutral.append({"name": str(name), "flux_key": key})
    energetic = []
    for index, item in enumerate(fluxes.energetic_fluxes):
        prefix = f"energetic_{index:03d}"
        if isinstance(item, EnergeticFlux):
            fields = {
                "flux_m2_s": item.flux_m2_s,
                "energy_eV": item.energy_eV,
                "cosine_incidence": item.cosine_incidence,
                "weight": item.weight,
            }
            kind = "EnergeticFlux"
            record = {"name": item.name, "kind": kind}
        elif isinstance(item, FaceResolvedEnergeticFlux):
            fields = {
                "event_face": item.event_face,
                "event_flux_m2_s": item.event_flux_m2_s,
                "event_energy_eV": item.event_energy_eV,
                "event_cosine_incidence": item.event_cosine_incidence,
            }
            if item.event_position is not None:
                fields["event_position"] = item.event_position
            if item.event_incident_direction is not None:
                fields["event_incident_direction"] = item.event_incident_direction
            kind = "FaceResolvedEnergeticFlux"
            record = {
                "name": item.name,
                "kind": kind,
                "face_count": int(item.face_count),
            }
        else:  # SurfaceFluxes already rejects this; retain a local artifact boundary.
            raise TypeError(f"unsupported energetic flux artifact type: {type(item)!r}")
        record["fields"] = {}
        for field, value in fields.items():
            key = f"{prefix}_{field}"
            arrays[key] = np.asarray(value)
            record["fields"][field] = key
        energetic.append(record)
    return arrays, {"neutral": neutral, "energetic": energetic}


def _surface_fluxes_from_artifact(arrays, layout):
    neutral = {
        item["name"]: np.asarray(arrays[item["flux_key"]]).copy()
        for item in layout.get("neutral", ())
    }
    energetic = []
    for item in layout.get("energetic", ()):
        fields = {
            name: np.asarray(arrays[key]).copy()
            for name, key in item["fields"].items()
        }
        if item["kind"] == "EnergeticFlux":
            energetic.append(EnergeticFlux(name=item["name"], **fields))
        elif item["kind"] == "FaceResolvedEnergeticFlux":
            energetic.append(FaceResolvedEnergeticFlux(
                name=item["name"], face_count=int(item["face_count"]), **fields))
        else:
            raise ValueError("unknown energetic flux type in direct-transport artifact")
    return SurfaceFluxes(neutral, tuple(energetic))


def _validate_direct_surface_fluxes(fluxes, *, face_count, species_role):
    """Reject structurally valid flux objects that belong to another surface."""
    neutral_names = set(fluxes.neutral_flux_m2_s)
    energetic_names = {item.name for item in fluxes.energetic_fluxes}
    if len(energetic_names) != len(fluxes.energetic_fluxes):
        raise ValueError("direct-transport artifact repeats an energetic species")
    expected_roles = {
        **{name: "neutral_reactant" for name in neutral_names},
        **{name: "energetic_bombardment" for name in energetic_names},
    }
    if dict(species_role) != expected_roles:
        raise ValueError("direct-transport species-role contract changed")
    for name, value in fluxes.neutral_flux_m2_s.items():
        if np.asarray(value).shape not in ((), (int(face_count),)):
            raise ValueError(f"cached neutral flux {name!r} has another face shape")
    for item in fluxes.energetic_fluxes:
        if isinstance(item, FaceResolvedEnergeticFlux):
            if item.face_count != int(face_count):
                raise ValueError("cached face-resolved energetic flux uses another mesh")
        elif np.asarray(item.flux_m2_s).shape not in ((), (int(face_count),)):
            raise ValueError("cached energetic flux has another face shape")


def _write_direct_transport_artifact(path, direct, *, identity):
    """Atomically persist the expensive pre-radiosity transport and its mesh contract."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    geometry_names = (
        "verts", "faces", "centroids", "areas_mesh", "face_area_m2",
        "gas_normals", "face_material_id", "active_face_index",
    )
    arrays = {name: np.asarray(direct[name]) for name in geometry_names}
    flux_arrays, flux_layout = _surface_flux_artifact_payload(
        direct["direct_surface_fluxes"])
    arrays.update(flux_arrays)
    array_manifest = {name: _array_receipt(value) for name, value in arrays.items()}
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)
    metadata = {
        "schema": DIRECT_TRANSPORT_ARTIFACT_SCHEMA,
        "identity": _jsonable(identity),
        "identity_sha256": _canonical_sha256(identity),
        "npz_sha256": _sha256(path),
        "array_manifest": array_manifest,
        "geometry_array_names": list(geometry_names),
        "surface_flux_layout": flux_layout,
        "surface_flux_sha256": surface_flux_sha256(
            direct["direct_surface_fluxes"]),
        "species_role": dict(sorted(direct["species_role"].items())),
        "reported_evaluator_wall_s": float(direct["reported_wall_s"]),
        "production_form_factor_call_count": int(
            direct["production_form_factor_call_count"]),
        "production_form_factor_trace_elided": bool(
            direct["production_form_factor_trace_elided"]),
    }
    _write_json_atomic(path.with_suffix(".json"), metadata)
    return metadata


def _read_direct_transport_artifact(path, source, *, identity):
    """Restore only a byte-certified artifact matching the current checkpoint mesh."""
    path = Path(path)
    metadata = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    if (metadata.get("schema") != DIRECT_TRANSPORT_ARTIFACT_SCHEMA
            or metadata.get("identity") != _jsonable(identity)
            or metadata.get("identity_sha256") != _canonical_sha256(identity)
            or metadata.get("npz_sha256") != _sha256(path)):
        raise ValueError("direct-transport artifact identity or byte digest changed")
    try:
        with np.load(path, allow_pickle=False) as data:
            names = set(data.files)
            expected_names = set(metadata.get("array_manifest", {}))
            if names != expected_names:
                raise ValueError("direct-transport artifact array set changed")
            arrays = {name: np.asarray(data[name]).copy() for name in names}
    except (BadZipFile, KeyError, OSError) as error:
        raise ValueError("direct-transport artifact payload is unreadable") from error
    observed_manifest = {
        name: _array_receipt(value) for name, value in arrays.items()
    }
    if observed_manifest != metadata.get("array_manifest"):
        raise ValueError("direct-transport artifact array receipt changed")

    expected_geometry = _direct_geometry_from_checkpoint(source)
    geometry_names = tuple(metadata.get("geometry_array_names", ()))
    if set(geometry_names) != set(expected_geometry):
        raise ValueError("direct-transport artifact geometry schema changed")
    for name, expected in expected_geometry.items():
        if not np.array_equal(arrays[name], np.asarray(expected)):
            raise ValueError(f"direct-transport artifact {name} changed")
    face_count = len(expected_geometry["faces"])
    areas_mesh = arrays["areas_mesh"]
    face_area_m2 = arrays["face_area_m2"]
    if (areas_mesh.shape != (face_count,) or face_area_m2.shape != (face_count,)
            or np.any(areas_mesh <= 0.0) or np.any(face_area_m2 <= 0.0)):
        raise ValueError("direct-transport artifact has invalid face areas")
    expected_physical_area = (
        areas_mesh * float(source["geometry"].mesh_length_unit_m) ** 2)
    if not np.array_equal(face_area_m2, expected_physical_area):
        raise ValueError("direct-transport artifact physical face areas changed")

    fluxes = _surface_fluxes_from_artifact(
        arrays, metadata.get("surface_flux_layout", {}))
    roles = metadata.get("species_role", {})
    _validate_direct_surface_fluxes(
        fluxes, face_count=face_count, species_role=roles)
    if surface_flux_sha256(fluxes) != metadata.get("surface_flux_sha256"):
        raise ValueError("direct-transport surface-flux digest changed")
    return {
        **expected_geometry,
        "reported_wall_s": float(metadata["reported_evaluator_wall_s"]),
        "direct_surface_fluxes": fluxes,
        "species_role": roles,
        "production_form_factor_call_count": int(
            metadata["production_form_factor_call_count"]),
        "production_form_factor_trace_elided": bool(
            metadata["production_form_factor_trace_elided"]),
    }, metadata


def _row_closure_error(factors):
    transferred = np.bincount(
        factors.source_face, weights=factors.transfer_fraction,
        minlength=factors.face_count)
    return float(np.max(np.abs(
        transferred + factors.escape_fraction - 1.0), initial=0.0))


def _integer_ray_counts(factors):
    """Recover the exact sampled hit/escape counts represented by one factor operator."""
    rays = int(factors.rays_per_face)
    transfer_float = np.asarray(factors.transfer_fraction) * rays
    escape_float = np.asarray(factors.escape_fraction) * rays
    transfer = np.rint(transfer_float).astype(np.int64)
    escape = np.rint(escape_float).astype(np.int64)
    error = max(
        float(np.max(np.abs(transfer_float - transfer), initial=0.0)),
        float(np.max(np.abs(escape_float - escape), initial=0.0)))
    if error > 5.0e-12:
        raise RuntimeError("form-factor fractions do not reconstruct integer ray counts")
    pair = (
        np.asarray(factors.source_face, dtype=np.int64) * factors.face_count
        + np.asarray(factors.target_face, dtype=np.int64))
    return {int(key): int(value) for key, value in zip(pair, transfer)}, escape, error


def _nested_sampling_extension_diagnostic(coarse, fine):
    """Prove that every fine scramble is the exact Sobol extension of its coarse twin."""
    if (coarse.replicate_seeds != fine.replicate_seeds
            or fine.rays_per_replicate != 2 * coarse.rays_per_replicate
            or len(coarse.replicate_form_factors)
            != len(fine.replicate_form_factors)):
        raise ValueError("nested ensembles do not share seeds or a doubled ray level")
    records = []
    for seed, left, right in zip(
            coarse.replicate_seeds, coarse.replicate_form_factors,
            fine.replicate_form_factors):
        left_count, left_escape, left_error = _integer_ray_counts(left)
        right_count, right_escape, right_error = _integer_ray_counts(right)
        keys = set(left_count) | set(right_count)
        hit_extension = np.asarray([
            right_count.get(key, 0) - left_count.get(key, 0)
            for key in keys], dtype=np.int64)
        escape_extension = right_escape - left_escape
        negative = int(np.sum(hit_extension < 0) + np.sum(escape_extension < 0))
        expected_added = int(
            left.face_count * (right.rays_per_face - left.rays_per_face))
        observed_added = int(np.sum(hit_extension) + np.sum(escape_extension))
        records.append({
            "seed": int(seed),
            "coarse_rays_per_face": int(left.rays_per_face),
            "fine_rays_per_face": int(right.rays_per_face),
            "maximum_integer_reconstruction_error": max(left_error, right_error),
            "negative_extension_count": negative,
            "expected_added_ray_count": expected_added,
            "observed_added_ray_count": observed_added,
            "pass": bool(negative == 0 and observed_added == expected_added),
        })
    return {"replicates": records, "all_gates_pass": all(item["pass"] for item in records)}


def _reciprocity_summary(diagnostic):
    return {
        "unordered_pair_count": diagnostic.unordered_pair_count,
        "one_sided_pair_count": diagnostic.one_sided_pair_count,
        "absolute_l1_rate_area_m2": diagnostic.absolute_l1_rate_area_m2,
        "relative_l1_error": diagnostic.relative_l1_error,
        "absolute_linf_rate_area_m2": diagnostic.absolute_linf_rate_area_m2,
        "relative_linf_error": diagnostic.relative_linf_error,
    }


def _estimate_replicated_level(
        direct, source, *, rays_per_face, seeds, maximum_replicate_wall_s,
        total_deadline_s, started):
    """Build one replicate at a time so every native trace has its own hard cap."""
    receipts = []
    geometry = source["geometry"]
    options = dict(
        rays_per_face=int(rays_per_face),
        source_sampling=OPERATOR["source_sampling"],
        visibility_mode=OPERATOR["visibility_mode"],
        maximum_visibility_wraps=OPERATOR["maximum_visibility_wraps"],
        maximum_visibility_replay_wraps=OPERATOR[
            "maximum_exact_replay_wraps"],
        periodic_lateral=OPERATOR["periodic_lateral"],
        domain_size=(np.asarray(geometry.phi.shape) - 1) * geometry.dx,
        ray_offset=OPERATOR["ray_offset_dx"] * geometry.dx,
        device=OPERATOR["transport_device"],
        return_visibility_receipt=True,
    )
    for seed in seeds:
        remaining = float(total_deadline_s) - (perf_counter() - started)
        with _hard_deadline(min(float(maximum_replicate_wall_s), remaining)):
            receipt = estimate_diffuse_form_factors_3d(
                direct["verts"], direct["faces"], direct["centroids"],
                direct["gas_normals"], seed=int(seed), **options)
        if receipt.visibility_mode != OPERATOR["visibility_mode"]:
            raise RuntimeError("form-factor replicate changed visibility authority")
        receipts.append(receipt)
    identity = {
        "schema": SCHEMA,
        "checkpoint_sha256": _sha256(source["checkpoint_path"]),
        "operator": OPERATOR,
        "rays_per_replicate": int(rays_per_face),
        "replicate_seeds": list(int(value) for value in seeds),
        "construction_call_count": len(receipts),
    }
    ensemble = ReplicatedDiffuseFormFactors3D(
        tuple(item.form_factors for item in receipts), tuple(seeds),
        direct["face_area_m2"], source_sampling=OPERATOR["source_sampling"],
        construction_identity=identity)
    row_errors = tuple(_row_closure_error(item.form_factors) for item in receipts)
    row_errors += (_row_closure_error(ensemble.mean_form_factors),)
    return ensemble, receipts, max(row_errors)


def _write_form_factor_ensemble(path, ensemble, *, checkpoint_sha256, ray_level):
    """Persist sampled factors so Stage B never pays Stage A's trace cost again."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {"face_area_m2": np.asarray(ensemble.face_area_m2)}
    for index, factors in enumerate(ensemble.replicate_form_factors):
        for name in (
                "source_face", "target_face", "transfer_fraction",
                "escape_fraction"):
            arrays[f"replicate_{index:02d}_{name}"] = np.asarray(
                getattr(factors, name))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)
    metadata = {
        "schema": FORM_FACTOR_ARTIFACT_SCHEMA,
        "checkpoint_sha256": str(checkpoint_sha256),
        "ray_level": int(ray_level),
        "replicate_seeds": ensemble.replicate_seeds,
        "replicate_count": len(ensemble.replicate_form_factors),
        "source_sampling": ensemble.source_sampling,
        "construction_identity": _jsonable(ensemble.construction_identity),
        "ensemble_sha256": ensemble.sha256,
        "npz_sha256": _sha256(path),
    }
    _write_json_atomic(path.with_suffix(".json"), metadata)
    return metadata


def _read_form_factor_ensemble(
        path, *, checkpoint_sha256, ray_level, replicate_seeds,
        expected_face_area_m2):
    """Load only a byte-certified Stage-A ensemble with the exact current mesh contract."""
    path = Path(path)
    metadata_path = path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected_seeds = tuple(int(value) for value in replicate_seeds)
    if (metadata.get("schema") != FORM_FACTOR_ARTIFACT_SCHEMA
            or metadata.get("checkpoint_sha256") != str(checkpoint_sha256)
            or int(metadata.get("ray_level", -1)) != int(ray_level)
            or tuple(metadata.get("replicate_seeds", ())) != expected_seeds
            or metadata.get("replicate_count") != len(expected_seeds)
            or metadata.get("source_sampling") != OPERATOR["source_sampling"]
            or metadata.get("npz_sha256") != _sha256(path)):
        raise ValueError("Stage-A form-factor artifact identity does not match Stage B")
    identity = metadata.get("construction_identity", {})
    if (identity.get("schema") != SCHEMA
            or identity.get("checkpoint_sha256") != str(checkpoint_sha256)
            or identity.get("operator") != OPERATOR
            or int(identity.get("rays_per_replicate", -1)) != int(ray_level)
            or tuple(identity.get("replicate_seeds", ())) != expected_seeds):
        raise ValueError("Stage-A form-factor construction identity changed")
    factors = []
    with np.load(path, allow_pickle=False) as data:
        area = np.asarray(data["face_area_m2"]).copy()
        if not np.array_equal(area, np.asarray(expected_face_area_m2)):
            raise ValueError("Stage-A form-factor artifact uses another face-area mesh")
        for index in range(len(expected_seeds)):
            factors.append(DiffuseFormFactors3D(
                len(area),
                np.asarray(data[f"replicate_{index:02d}_source_face"]),
                np.asarray(data[f"replicate_{index:02d}_target_face"]),
                np.asarray(data[f"replicate_{index:02d}_transfer_fraction"]),
                np.asarray(data[f"replicate_{index:02d}_escape_fraction"]),
                int(ray_level)))
    ensemble = ReplicatedDiffuseFormFactors3D(
        tuple(factors), expected_seeds, area,
        source_sampling=metadata["source_sampling"],
        construction_identity=metadata["construction_identity"])
    if ensemble.sha256 != metadata.get("ensemble_sha256"):
        raise ValueError("Stage-A form-factor ensemble digest does not reproduce")
    return ensemble, metadata


def _maximum_ledger_residual(exchange):
    maximum = 0.0
    for name, removed in exchange.removed_units_m2.items():
        source, outgoing, unresolved = np.broadcast_arrays(
            np.asarray(removed, dtype=float),
            np.asarray(exchange.outgoing_units_m2.get(name, 0.0), dtype=float),
            np.asarray(exchange.unresolved_units_m2.get(name, 0.0), dtype=float))
        residual = source - outgoing - unresolved
        maximum = max(maximum, float(np.max(np.abs(residual), initial=0.0)))
    return maximum


def _response_from_material_step(initial_state, step, horizon_s, balance_error):
    if step.product_populations:
        raise RuntimeError("replicated closure refuses emitted surface products")
    if not step.validity.within_declared_scope:
        raise RuntimeError("material mechanism left declared scope")
    ledger = _maximum_ledger_residual(step.material_exchange)
    if ledger != 0.0:
        raise RuntimeError("material mechanism ledger is not exactly closed")
    fields = {
        f"state_increment/{name}": (
            np.asarray(step.state.fields[name], dtype=float)
            - np.asarray(initial_state.fields[name], dtype=float))
        for name in sorted(initial_state.fields)
    }
    for kind, attribute in _EXCHANGE_ATTRIBUTES.items():
        for name, value in sorted(getattr(step.material_exchange, attribute).items()):
            fields[f"exchange/{kind}/{name}"] = np.asarray(value, dtype=float)
    fields["surface/integrated_recession_m"] = (
        np.asarray(step.etch_velocity_m_s, dtype=float) * float(horizon_s))
    fields["surface/integrated_growth_m"] = (
        np.asarray(step.normal_growth_velocity_m_s, dtype=float) * float(horizon_s))
    maximum_gross_displacement = float(np.max(
        fields["surface/integrated_recession_m"]
        + fields["surface/integrated_growth_m"], initial=0.0))
    return {
        "fields": fields,
        "maximum_radiosity_relative_balance_error": float(balance_error),
        "maximum_material_ledger_residual_units_m2": float(ledger),
        "maximum_gross_displacement_m": maximum_gross_displacement,
    }


def _instantaneous_response(
        direct, source, factors, mechanism, *, horizon_s, identity):
    cache = SurfaceRadiosityOperatorCache3D(
        direct["direct_surface_fluxes"], factors, direct["face_area_m2"],
        direct["active_face_index"], direct["face_material_id"],
        direct["species_role"], identity,
        relative_tolerance=OPERATOR["radiosity_relative_tolerance"],
        maximum_iterations=OPERATOR["radiosity_maximum_iterations"])
    solved = cache.solve(
        source["state"], mechanism, operator_identity_payload=identity)
    step = mechanism.advance_by_material(
        source["state"], solved.surface_fluxes, float(horizon_s),
        direct["face_material_id"])
    response = _response_from_material_step(
        source["state"], step, horizon_s,
        solved.maximum_relative_balance_error)
    for name, incident in sorted(solved.surface_fluxes.neutral_flux_m2_s.items()):
        incident = np.asarray(incident, dtype=float)
        response["fields"][f"radiosity/incident/{name}"] = incident
        response["fields"][f"radiosity/reacted/{name}"] = (
            incident * np.asarray(solved.reaction_probability[name], dtype=float))
    response["radiosity_species_diagnostics"] = {
        name: {
            "source_rate_s": value.source_rate_s,
            "reacted_rate_s": value.reacted_rate_s,
            "escaped_rate_s": value.escaped_rate_s,
            "relative_balance_error": value.relative_balance_error,
        }
        for name, value in sorted(solved.species_diagnostics.items())
    }
    return response


def _integrated_observables(response, direct):
    area = direct["face_area_m2"]
    material = direct["face_material_id"]
    recession = response["fields"]["surface/integrated_recession_m"]
    growth = response["fields"]["surface/integrated_growth_m"]

    def mean(values, material_id):
        selected = material == material_id
        selected_area = float(np.sum(area[selected]))
        if selected_area <= 0.0:
            raise ValueError(
                f"checkpoint lacks required material id {int(material_id)}")
        return float(np.sum(area[selected] * values[selected]) / selected_area)

    return {
        "sio2_mean_recession_m": mean(recession, 1),
        "mask_mean_recession_m": mean(recession, 2),
        "mask_mean_growth_m": mean(growth, 2),
        "mask_mean_net_recession_m": mean(recession - growth, 2),
    }


def _scalar_replicate_score(
        authority, replicates, *, absolute_tolerance, relative_tolerance,
        confidence_level=0.95):
    values = np.asarray(replicates, dtype=float)
    if values.ndim != 1 or len(values) < 4 or np.any(~np.isfinite(values)):
        raise ValueError("scalar score requires at least four finite replicates")
    authority = float(authority)
    mean = float(np.mean(values))
    critical = float(student_t.ppf(
        0.5 + 0.5 * float(confidence_level), len(values) - 1))
    half = critical * float(np.std(values, ddof=1)) / np.sqrt(len(values))
    bias = abs(authority - mean)
    scale = (
        float(absolute_tolerance)
        + float(relative_tolerance) * max(abs(authority), abs(mean)))
    combined = (half + bias) / scale
    return {
        "authority": authority,
        "replicate_mean": mean,
        "confidence_half_width": half,
        "authority_to_replicate_mean_bias": bias,
        "combined_mixed_normalized": combined,
        "pass": bool(combined <= 1.0),
    }


def _nested_scalar_score(
        coarse, fine, *, absolute_tolerance, relative_tolerance):
    coarse = float(coarse)
    fine = float(fine)
    scale = (
        float(absolute_tolerance)
        + float(relative_tolerance) * max(abs(coarse), abs(fine)))
    normalized = abs(fine - coarse) / scale
    return {
        "coarse": coarse,
        "fine": fine,
        "difference": fine - coarse,
        "mixed_normalized": normalized,
        "pass": bool(normalized <= 1.0),
    }


def _response_field_names(authority, replicates):
    names = set(authority["fields"])
    if any(set(item["fields"]) != names for item in replicates):
        raise RuntimeError("replicate responses changed conservative field contract")
    return tuple(sorted(names))


def _physical_patch_operator_contract(source):
    """Bind physical-patch support to the represented checkpoint domain.

    Krueger is represented by one periodic cell in the homogeneous y direction.  A 40 nm patch
    therefore has only 20 nm of represented y footprint on the production checkpoint; treating the
    unavailable repeated half as empty would misclassify every otherwise full 40 nm patch.
    """
    geometry = source["geometry"]
    shape = np.asarray(geometry.phi.shape, dtype=np.int64)
    origin = np.asarray(geometry.mesh_origin_m, dtype=float)
    unit = float(geometry.mesh_length_unit_m)
    dx = float(geometry.dx)
    if (shape.shape != (3,) or np.any(shape < 2) or origin.shape != (3,)
            or np.any(~np.isfinite(origin)) or not np.isfinite(unit) or unit <= 0.0
            or not np.isfinite(dx) or dx <= 0.0):
        raise ValueError("checkpoint lacks a valid physical-patch domain contract")
    domain_lengths = (shape - 1) * dx * unit
    periodic_lengths = np.where(
        np.asarray(PATCH_PERIODIC_AXES, dtype=bool), domain_lengths, 0.0)
    if np.any(periodic_lengths[np.asarray(PATCH_PERIODIC_AXES)] <= 0.0):
        raise ValueError("checkpoint periodic physical-patch domain has zero extent")
    return {
        "schema": "petch.physical-patch-support-operator.v1",
        "patch_scales_m": [float(value) for value in PATCH_SCALES_M],
        "projected_support_measure": "dominant-axis projected triangle overlap",
        "represented_nominal_footprint": (
            "tangential Cartesian patch footprint intersected with the represented "
            "periodic fundamental domain"),
        "integrated_gate_area_scale": "represented nominal projected footprint; every patch",
        "mean_gate_domain": "support-eligible patches only; all-patch mean retained diagnostic",
        "selected_source_allocator_artifact_role": (
            "historical diagnostic only; recomputation required under this operator"),
        "minimum_mean_support_fraction": float(
            DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
        "support_sensitivity_thresholds": [
            float(value) for value in PATCH_SUPPORT_SENSITIVITY_THRESHOLDS],
        "mesh_origin_m": origin.tolist(),
        "checkpoint_domain_lengths_m": domain_lengths.tolist(),
        "periodic_axes_xyz": [bool(value) for value in PATCH_PERIODIC_AXES],
        "periodic_domain_origin_m": origin.tolist(),
        "periodic_domain_lengths_m": periodic_lengths.tolist(),
    }


def _physical_patch_aggregation_kwargs(source):
    contract = _physical_patch_operator_contract(source)
    return {
        "mesh_origin_m": contract["mesh_origin_m"],
        "patch_origin_m": contract["mesh_origin_m"],
        "periodic_domain_origin_m": contract["periodic_domain_origin_m"],
        "periodic_domain_lengths_m": contract["periodic_domain_lengths_m"],
    }


def _support_inventory_summary(
        surface_area, projected_area, nominal_area, support_fraction,
        eligible_mask, *, threshold):
    surface_area = np.asarray(surface_area, dtype=float)
    projected_area = np.asarray(projected_area, dtype=float)
    nominal_area = np.asarray(nominal_area, dtype=float)
    support_fraction = np.asarray(support_fraction, dtype=float)
    eligible = np.asarray(eligible_mask, dtype=bool)
    if (surface_area.ndim != 1 or projected_area.shape != surface_area.shape
            or nominal_area.shape != surface_area.shape
            or support_fraction.shape != surface_area.shape
            or eligible.shape != surface_area.shape or not np.any(eligible)):
        raise ValueError("physical-patch mean gate has no support-eligible patch")
    excluded = ~eligible
    surface_total = float(np.sum(surface_area))
    projected_total = float(np.sum(projected_area))
    return {
        "minimum_mean_support_fraction": float(threshold),
        "eligible_mean_patch_count": int(np.count_nonzero(eligible)),
        "excluded_mean_patch_count": int(np.count_nonzero(excluded)),
        "excluded_mean_surface_area_m2": float(np.sum(
            surface_area[excluded])),
        "excluded_mean_surface_area_fraction": float(np.sum(
            surface_area[excluded]) / surface_total),
        "excluded_mean_projected_support_area_m2": float(np.sum(
            projected_area[excluded])),
        "excluded_mean_projected_support_fraction": float(np.sum(
            projected_area[excluded]) / projected_total),
        "represented_nominal_projected_area_m2": {
            "minimum": float(np.min(nominal_area)),
            "maximum": float(np.max(nominal_area)),
            "total": float(np.sum(nominal_area)),
        },
        "projected_support_fraction": {
            "minimum_all_patches": float(np.min(
                support_fraction)),
            "minimum_eligible": float(np.min(
                support_fraction[eligible])),
            "maximum": float(np.max(support_fraction)),
        },
    }


def _support_summary(authority, eligible_mask, *, threshold=None):
    threshold = (
        DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION if threshold is None
        else float(threshold))
    return _support_inventory_summary(
        authority.patch_area_m2,
        authority.patch_projected_support_area_m2,
        authority.patch_nominal_projected_area_m2,
        authority.patch_projected_support_fraction,
        eligible_mask, threshold=threshold)


def _replicated_patch_threshold_sensitivity(score):
    authority = score.authority
    mean_scale = (
        score.absolute_tolerance
        + score.relative_tolerance * np.maximum(
            np.abs(authority.mean_field),
            np.abs(score.replicate_mean_field)))
    mean_normalized = (
        score.mean_confidence_half_width
        + np.abs(authority.mean_field - score.replicate_mean_field)
    ) / mean_scale
    output = []
    for threshold in PATCH_SUPPORT_SENSITIVITY_THRESHOLDS:
        eligible = authority.patch_projected_support_fraction >= threshold
        if not np.any(eligible):
            output.append({
                "minimum_mean_support_fraction": float(threshold),
                "eligible_mean_patch_count": 0,
                "excluded_mean_patch_count": len(authority.patch_key),
                "gate_defined": False,
                "pass": False,
            })
            continue
        maximum_mean = float(np.max(mean_normalized[eligible], initial=0.0))
        output.append({
            **_support_summary(authority, eligible, threshold=threshold),
            "gate_defined": True,
            "maximum_support_eligible_mean_combined_mixed_normalized": maximum_mean,
            "maximum_integrated_combined_mixed_normalized": (
                score.maximum_integrated_combined_mixed_normalized),
            "pass": bool(
                maximum_mean <= 1.0
                and score.maximum_integrated_combined_mixed_normalized <= 1.0),
        })
    primary_pass = score.all_mixed_tolerances_pass
    return {
        "thresholds": output,
        "primary_threshold": float(DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
        "primary_gate_conclusion": bool(primary_pass),
        "gate_conclusion_stable_over_predeclared_thresholds": bool(all(
            item["pass"] == primary_pass for item in output)),
    }


def _patch_score_summary(score):
    return {
        "patch_scale_m": score.authority.patch_scale_m,
        "patch_count": len(score.authority.patch_key),
        "normalization_scale": None,
        **_support_summary(score.authority, score.mean_eligible_patch_mask),
        "maximum_integrated_combined_mixed_normalized": (
            score.maximum_integrated_combined_mixed_normalized),
        "maximum_support_eligible_mean_combined_mixed_normalized": (
            score.maximum_mean_combined_mixed_normalized),
        "maximum_all_patch_mean_combined_mixed_normalized_diagnostic": (
            score.maximum_mean_combined_all_patches_mixed_normalized),
        "support_threshold_sensitivity": (
            _replicated_patch_threshold_sensitivity(score)),
        "maximum_confidence_mixed_normalized": max(
            score.maximum_integrated_confidence_mixed_normalized,
            score.maximum_mean_confidence_mixed_normalized),
        "maximum_authority_bias_mixed_normalized": max(
            score.maximum_integrated_authority_bias_mixed_normalized,
            score.maximum_mean_authority_bias_mixed_normalized),
        "maximum_combined_mixed_normalized": max(
            score.maximum_integrated_combined_mixed_normalized,
            score.maximum_mean_combined_mixed_normalized),
        "pass": score.all_mixed_tolerances_pass,
        "scheme_sha256": score.authority.scheme_sha256,
    }


def _worst_face_uncertainty(authority, replicates, *, absolute, relative):
    authority = np.asarray(authority, dtype=float)
    values = np.asarray(replicates, dtype=float)
    critical = float(student_t.ppf(0.975, len(values) - 1))
    mean = np.mean(values, axis=0)
    half = critical * np.std(values, axis=0, ddof=1) / np.sqrt(len(values))
    bias = np.abs(authority - mean)
    denominator = absolute + relative * np.maximum(np.abs(authority), np.abs(mean))
    combined = (half + bias) / denominator
    index = int(np.argmax(combined))
    return {
        "face_index": index,
        "authority": float(authority[index]),
        "replicate_mean": float(mean[index]),
        "confidence_half_width": float(half[index]),
        "authority_bias": float(bias[index]),
        "combined_mixed_normalized": float(combined[index]),
        "gating": False,
    }


def _natural_response_field_scales(field_names, direct, source, *, dx_m):
    """Return predeclared physical scales, never scales inferred from sampled outcomes."""
    router = build_krueger_2024_material_router_3d()
    oxide = router.mechanisms[1].parameters
    mask = router.mechanisms[2].parameters
    inventory_scale = {
        "SiO2_formula_unit": oxide.bulk_formula_density_m3 * dx_m,
        "fluorocarbon_film_unit": oxide.polymer_monolayer_density_m2,
        "amorphous_carbon_atom": mask.bulk_carbon_atom_density_m3 * dx_m,
    }
    state_scale = {
        "m1__polymer_units_m2": oxide.polymer_monolayer_density_m2,
        "m1__removed_formula_units_m2": oxide.bulk_formula_density_m3 * dx_m,
        "m2__polymer_units_m2": mask.polymer_monolayer_density_m2,
        "m2__removed_carbon_atoms_m2": mask.bulk_carbon_atom_density_m3 * dx_m,
    }
    output = {}
    for name in field_names:
        if name.startswith("surface/"):
            scale = dx_m
        elif name.startswith("state_increment/"):
            state_name = name.split("/", 1)[1]
            upper = source["state"].upper_bounds.get(state_name)
            if upper is not None:
                scale = float(upper)
            elif state_name in state_scale:
                scale = float(state_scale[state_name])
            else:
                raise RuntimeError(
                    f"no predeclared physical scale for state field {state_name!r}")
        elif name.startswith("exchange/"):
            inventory_name = name.rsplit("/", 1)[1]
            if inventory_name not in inventory_scale:
                raise RuntimeError(
                    f"no predeclared physical scale for exchange {inventory_name!r}")
            scale = float(inventory_scale[inventory_name])
        elif name.startswith("radiosity/"):
            species = name.rsplit("/", 1)[1]
            if species not in direct["direct_surface_fluxes"].neutral_flux_m2_s:
                raise RuntimeError(f"no direct-flux scale for neutral {species!r}")
            scale = float(np.max(np.abs(
                direct["direct_surface_fluxes"].neutral_flux_m2_s[species]),
                initial=0.0))
        else:
            raise RuntimeError(f"no physical normalization contract for {name!r}")
        if not np.isfinite(scale) or scale <= 0.0:
            raise RuntimeError(f"invalid physical normalization scale for {name!r}")
        output[name] = scale
    return output


def _score_response_level(authority, replicates, direct, source, *, dx_m):
    fields = _response_field_names(authority, replicates)
    natural_scale = _natural_response_field_scales(
        fields, direct, source, dx_m=dx_m)
    absolute_displacement = GATES["integrated_absolute_displacement_dx"] * dx_m
    authority_integrated = _integrated_observables(authority, direct)
    replicate_integrated = [
        _integrated_observables(item, direct) for item in replicates]
    integrated = {
        name: _scalar_replicate_score(
            value, [item[name] for item in replicate_integrated],
            absolute_tolerance=absolute_displacement,
            relative_tolerance=GATES["integrated_relative_tolerance"])
        for name, value in authority_integrated.items()
    }
    area = direct["face_area_m2"]
    total_area = float(np.sum(area))
    for name in fields:
        if not name.startswith("radiosity/"):
            continue
        authority_mean = float(np.sum(
            area * np.asarray(authority["fields"][name])) / total_area)
        replicate_mean = [
            float(np.sum(area * np.asarray(item["fields"][name])) / total_area)
            for item in replicates]
        integrated[f"area_mean/{name}"] = _scalar_replicate_score(
            authority_mean, replicate_mean,
            absolute_tolerance=(
                GATES["integrated_absolute_displacement_dx"]
                * natural_scale[name]),
            relative_tolerance=GATES["integrated_relative_tolerance"])
    species_ledger = {}
    if "radiosity_species_diagnostics" in authority:
        species_names = set(authority["radiosity_species_diagnostics"])
        if any(set(item.get("radiosity_species_diagnostics", ())) != species_names
               for item in replicates):
            raise RuntimeError("replicate radiosity species ledgers differ")
        for species in sorted(species_names):
            source_rate = authority["radiosity_species_diagnostics"][species][
                "source_rate_s"]
            species_ledger[species] = {}
            for quantity in ("source_rate_s", "reacted_rate_s", "escaped_rate_s"):
                species_ledger[species][quantity] = _scalar_replicate_score(
                    authority["radiosity_species_diagnostics"][species][quantity],
                    [item["radiosity_species_diagnostics"][species][quantity]
                     for item in replicates],
                    absolute_tolerance=max(
                        GATES["integrated_absolute_displacement_dx"]
                        * abs(source_rate), np.finfo(float).tiny),
                    relative_tolerance=GATES["integrated_relative_tolerance"])
    patch = {}
    worst_face = {}
    for name in fields:
        authority_field = np.asarray(authority["fields"][name], dtype=float)
        replicate_field = np.stack([
            np.asarray(item["fields"][name], dtype=float) for item in replicates])
        normalization = natural_scale[name]
        scores = score_replicated_surface_field_at_physical_scales_3d(
            authority_field / normalization, replicate_field / normalization,
            direct["face_area_m2"], direct["verts"], direct["faces"],
            direct["gas_normals"], direct["face_material_id"], PATCH_SCALES_M,
            absolute_tolerance=GATES["patch_absolute_normalized"],
            relative_tolerance=GATES["patch_relative_tolerance"],
            mesh_length_unit_m=source["geometry"].mesh_length_unit_m,
            minimum_mean_support_fraction=(
                DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
            **_physical_patch_aggregation_kwargs(source))
        summaries = []
        for score in scores:
            summary = _patch_score_summary(score)
            summary["normalization_scale"] = normalization
            summaries.append(summary)
        patch[name] = summaries
        worst_face[name] = _worst_face_uncertainty(
            authority_field / normalization, replicate_field / normalization,
            absolute=GATES["patch_absolute_normalized"],
            relative=GATES["patch_relative_tolerance"])
    maximum_balance = max(
        [authority["maximum_radiosity_relative_balance_error"]]
        + [item["maximum_radiosity_relative_balance_error"] for item in replicates])
    exact_ledger = all(
        item["maximum_material_ledger_residual_units_m2"] == 0.0
        for item in (authority, *replicates))
    maximum_gross_displacement_dx = max(
        item["maximum_gross_displacement_m"] / dx_m
        for item in (authority, *replicates))
    integration_records = [
        item for item in (authority, *replicates)
        if "maximum_embedded_relative_error" in item]
    if integration_records and len(integration_records) != 1 + len(replicates):
        raise RuntimeError("responses mix instantaneous and integrated numerical contracts")
    integration_contract = bool(not integration_records or all(
        item["maximum_embedded_relative_error"]
        <= GATES["maximum_embedded_relative_error"]
        and item["maximum_local_reaction_probability_change"]
        <= GATES["maximum_local_reaction_probability_change"]
        and item["rejected_trial_exchange_contribution_is_zero"]
        for item in integration_records))
    gates = {
        "integrated_uncertainty": all(item["pass"] for item in integrated.values()),
        "radiosity_species_ledger_uncertainty": all(
            item["pass"] for species in species_ledger.values()
            for item in species.values()),
        "physical_patch_uncertainty": all(
            score["pass"] for scores in patch.values() for score in scores),
        "radiosity_balance": bool(
            maximum_balance <= GATES["maximum_radiosity_relative_balance_error"]),
        "exact_material_ledger": bool(exact_ledger),
        "gross_displacement_within_frozen_geometry_limit": bool(
            maximum_gross_displacement_dx
            <= GATES["maximum_gross_displacement_dx"]),
        "integration_numerical_contract": integration_contract,
    }
    return {
        "integrated_observables": integrated,
        "radiosity_species_ledgers": species_ledger,
        "physical_patch_fields": patch,
        "worst_face_diagnostics": worst_face,
        "maximum_radiosity_relative_balance_error": maximum_balance,
        "maximum_gross_displacement_dx": maximum_gross_displacement_dx,
        "gates": gates,
        "all_gates_pass": bool(all(gates.values())),
    }


def _nested_integrated_scores(coarse, fine, direct, *, dx_m):
    first = _integrated_observables(coarse, direct)
    second = _integrated_observables(fine, direct)
    absolute = GATES["integrated_absolute_displacement_dx"] * dx_m
    return {
        name: _nested_scalar_score(
            first[name], second[name], absolute_tolerance=absolute,
            relative_tolerance=GATES["integrated_relative_tolerance"])
        for name in sorted(first)
    }


def _refinement_patch_summary(item, *, normalization_scale):
    primary_pass = bool(
        item.integrated_mixed_normalized_linf <= 1.0
        and item.mean_mixed_normalized_linf <= 1.0
        and item.maximum_patch_area_relative_error <= 2.0e-10)
    sensitivity = []
    for threshold in PATCH_SUPPORT_SENSITIVITY_THRESHOLDS:
        eligible = item.common_patch_projected_support_fraction >= threshold
        if not np.any(eligible):
            sensitivity.append({
                "minimum_mean_support_fraction": float(threshold),
                "eligible_mean_patch_count": 0,
                "excluded_mean_patch_count": item.patch_count,
                "gate_defined": False,
                "pass": False,
            })
            continue
        maximum_mean = float(np.max(
            item.mean_mixed_normalized_by_patch[eligible], initial=0.0))
        sensitivity.append({
            **_support_inventory_summary(
                item.common_patch_surface_area_m2,
                item.common_patch_projected_support_area_m2,
                item.patch_nominal_projected_area_m2,
                item.common_patch_projected_support_fraction,
                eligible, threshold=threshold),
            "gate_defined": True,
            "mean_mixed_normalized_linf": maximum_mean,
            "integrated_mixed_normalized_linf": (
                item.integrated_mixed_normalized_linf),
            "pass": bool(
                maximum_mean <= 1.0
                and item.integrated_mixed_normalized_linf <= 1.0
                and item.maximum_patch_area_relative_error <= 2.0e-10),
        })
    return {
        "patch_scale_m": item.patch_scale_m,
        "patch_count": item.patch_count,
        "normalization_scale": normalization_scale,
        "integrated_mixed_normalized_linf": (
            item.integrated_mixed_normalized_linf),
        "mean_mixed_normalized_linf": item.mean_mixed_normalized_linf,
        "mean_all_patch_mixed_normalized_linf_diagnostic": (
            item.mean_all_patch_mixed_normalized_linf),
        "minimum_mean_support_fraction": item.minimum_mean_support_fraction,
        "eligible_mean_patch_count": item.eligible_mean_patch_count,
        "excluded_mean_patch_count": item.excluded_mean_patch_count,
        "excluded_mean_surface_area_fraction": (
            item.excluded_mean_surface_area_fraction),
        "excluded_mean_projected_support_fraction": (
            item.excluded_mean_projected_support_fraction),
        "represented_nominal_projected_area_m2": {
            "minimum": float(np.min(item.patch_nominal_projected_area_m2)),
            "maximum": float(np.max(item.patch_nominal_projected_area_m2)),
            "total": float(np.sum(item.patch_nominal_projected_area_m2)),
        },
        "maximum_patch_area_relative_error": (
            item.maximum_patch_area_relative_error),
        "support_threshold_sensitivity": {
            "thresholds": sensitivity,
            "primary_threshold": float(DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
            "primary_gate_conclusion": primary_pass,
            "gate_conclusion_stable_over_predeclared_thresholds": bool(all(
                record["pass"] == primary_pass for record in sensitivity)),
        },
        "pass": primary_pass,
    }


def _nested_patch_scores(coarse, fine, direct, source, *, dx_m):
    """Compare the 8/16 or nominal/tight authority fields at fixed physical scales."""
    names = _response_field_names(coarse, (fine,))
    natural_scale = _natural_response_field_scales(
        names, direct, source, dx_m=dx_m)
    output = {}
    for name in names:
        coarse_field = np.asarray(coarse["fields"][name], dtype=float)
        fine_field = np.asarray(fine["fields"][name], dtype=float)
        normalization = natural_scale[name]
        common = {
            "face_area_m2": direct["face_area_m2"],
            "verts": direct["verts"],
            "faces": direct["faces"],
            "face_gas_normals": direct["gas_normals"],
            "face_material_id": direct["face_material_id"],
            "mesh_length_unit_m": source["geometry"].mesh_length_unit_m,
            **_physical_patch_aggregation_kwargs(source),
        }
        scores = score_surface_field_refinement_at_physical_scales_3d(
            {"face_field": coarse_field / normalization, **common},
            {"face_field": fine_field / normalization, **common},
            PATCH_SCALES_M,
            absolute_tolerance=GATES["patch_absolute_normalized"],
            relative_tolerance=GATES["patch_relative_tolerance"],
            minimum_mean_support_fraction=(
                DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION))
        output[name] = tuple(
            _refinement_patch_summary(item, normalization_scale=normalization)
            for item in scores)
    return output


def _compare_response_refinement(coarse, fine, direct, source, *, dx_m):
    integrated = _nested_integrated_scores(coarse, fine, direct, dx_m=dx_m)
    patch = _nested_patch_scores(coarse, fine, direct, source, dx_m=dx_m)
    gates = {
        "integrated": all(item["pass"] for item in integrated.values()),
        "physical_patches": all(
            item["pass"] for values in patch.values() for item in values),
    }
    return {
        "integrated_observables": integrated,
        "physical_patch_fields": patch,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def _paired_scalar_refinement_score(
        coarse_authority, fine_authority, coarse_replicates, fine_replicates,
        *, absolute_tolerance, relative_tolerance):
    coarse_values = np.asarray(coarse_replicates, dtype=float)
    fine_values = np.asarray(fine_replicates, dtype=float)
    if coarse_values.shape != fine_values.shape or coarse_values.ndim != 1:
        raise ValueError("paired scalar refinement samples must match")
    paired = fine_values - coarse_values
    critical = float(student_t.ppf(0.975, len(paired) - 1))
    mean = float(np.mean(paired))
    half = critical * float(np.std(paired, ddof=1)) / np.sqrt(len(paired))
    authority_difference = float(fine_authority) - float(coarse_authority)
    scale = (
        float(absolute_tolerance)
        + float(relative_tolerance)
        * max(abs(float(coarse_authority)), abs(float(fine_authority))))
    combined = (abs(authority_difference) + half) / scale
    return {
        "authority_difference": authority_difference,
        "paired_replicate_mean_difference": mean,
        "paired_95pct_confidence_half_width": half,
        "authority_to_paired_mean_bias": abs(authority_difference - mean),
        "combined_mixed_normalized": combined,
        "pass": bool(combined <= 1.0),
    }


def _extended_integrated_observables(response, direct):
    output = dict(_integrated_observables(response, direct))
    area = direct["face_area_m2"]
    total_area = float(np.sum(area))
    for name, value in response["fields"].items():
        if name.startswith("radiosity/"):
            output[f"area_mean/{name}"] = float(
                np.sum(area * np.asarray(value, dtype=float)) / total_area)
    return output


def _paired_patch_refinement_score(
        coarse_authority, fine_authority, coarse_replicates, fine_replicates,
        direct, source, *, normalization, patch_scale_m):
    common = dict(
        face_area_m2=direct["face_area_m2"], verts=direct["verts"],
        faces=direct["faces"], face_gas_normals=direct["gas_normals"],
        face_material_id=direct["face_material_id"],
        patch_scale_m=patch_scale_m,
        mesh_length_unit_m=source["geometry"].mesh_length_unit_m,
        **_physical_patch_aggregation_kwargs(source))
    coarse_patch = aggregate_surface_field_on_physical_patches_3d(
        np.asarray(coarse_authority) / normalization, **common)
    fine_patch = aggregate_surface_field_on_physical_patches_3d(
        np.asarray(fine_authority) / normalization, **common)
    if (not np.array_equal(coarse_patch.patch_key, fine_patch.patch_key)
            or coarse_patch.scheme_sha256 != fine_patch.scheme_sha256
            or not np.allclose(
                coarse_patch.patch_area_m2, fine_patch.patch_area_m2,
                rtol=2e-10, atol=1e-30)
            or not np.allclose(
                coarse_patch.patch_projected_support_area_m2,
                fine_patch.patch_projected_support_area_m2,
                rtol=2e-10, atol=1e-30)
            or not np.allclose(
                coarse_patch.patch_nominal_projected_area_m2,
                fine_patch.patch_nominal_projected_area_m2,
                rtol=0.0, atol=0.0)):
        raise RuntimeError("paired refinement changed the physical-patch scheme")
    eligible = np.minimum(
        coarse_patch.patch_projected_support_fraction,
        fine_patch.patch_projected_support_fraction,
    ) >= DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION
    if not np.any(eligible):
        raise ValueError("paired physical-patch mean gate has no support-eligible patch")
    paired_face = (
        np.asarray(fine_replicates, dtype=float)
        - np.asarray(coarse_replicates, dtype=float)) / normalization
    paired_integrated = np.empty((len(paired_face), len(coarse_patch.patch_key)))
    for index, field in enumerate(paired_face):
        paired_integrated[index] = np.bincount(
            coarse_patch.contribution_patch_index,
            weights=(coarse_patch.contribution_area_m2
                     * field[coarse_patch.contribution_face_index]),
            minlength=len(coarse_patch.patch_key))
    paired_mean = paired_integrated / coarse_patch.patch_area_m2[None, :]
    critical = float(student_t.ppf(0.975, len(paired_face) - 1))
    integrated_half = critical * np.std(
        paired_integrated, axis=0, ddof=1) / np.sqrt(len(paired_face))
    mean_half = critical * np.std(
        paired_mean, axis=0, ddof=1) / np.sqrt(len(paired_face))
    authority_integrated = np.abs(
        fine_patch.integrated_field_area - coarse_patch.integrated_field_area)
    authority_mean = np.abs(fine_patch.mean_field - coarse_patch.mean_field)
    integrated_scale = (
        GATES["patch_absolute_normalized"]
        * coarse_patch.patch_nominal_projected_area_m2
        + GATES["patch_relative_tolerance"] * np.maximum(
            np.abs(coarse_patch.integrated_field_area),
            np.abs(fine_patch.integrated_field_area)))
    mean_scale = (
        GATES["patch_absolute_normalized"]
        + GATES["patch_relative_tolerance"] * np.maximum(
            np.abs(coarse_patch.mean_field), np.abs(fine_patch.mean_field)))
    maximum_integrated = float(np.max(
        (authority_integrated + integrated_half) / integrated_scale, initial=0.0))
    normalized_mean = (authority_mean + mean_half) / mean_scale
    maximum_mean = float(np.max(normalized_mean[eligible], initial=0.0))
    maximum_all_patch_mean = float(np.max(normalized_mean, initial=0.0))
    primary_pass = bool(maximum_integrated <= 1.0 and maximum_mean <= 1.0)
    sensitivity = []
    common_support_fraction = np.minimum(
        coarse_patch.patch_projected_support_fraction,
        fine_patch.patch_projected_support_fraction)
    for threshold in PATCH_SUPPORT_SENSITIVITY_THRESHOLDS:
        local_eligible = common_support_fraction >= threshold
        if not np.any(local_eligible):
            sensitivity.append({
                "minimum_mean_support_fraction": float(threshold),
                "eligible_mean_patch_count": 0,
                "excluded_mean_patch_count": len(coarse_patch.patch_key),
                "gate_defined": False,
                "pass": False,
            })
            continue
        local_maximum_mean = float(np.max(
            normalized_mean[local_eligible], initial=0.0))
        sensitivity.append({
            **_support_summary(
                coarse_patch, local_eligible, threshold=threshold),
            "gate_defined": True,
            "maximum_support_eligible_mean_authority_plus_paired_ci_mixed_normalized": (
                local_maximum_mean),
            "maximum_integrated_authority_plus_paired_ci_mixed_normalized": (
                maximum_integrated),
            "pass": bool(
                maximum_integrated <= 1.0 and local_maximum_mean <= 1.0),
        })
    return {
        "patch_scale_m": float(patch_scale_m),
        "patch_count": len(coarse_patch.patch_key),
        "normalization_scale": float(normalization),
        **_support_summary(coarse_patch, eligible),
        "maximum_integrated_authority_plus_paired_ci_mixed_normalized": (
            maximum_integrated),
        "maximum_support_eligible_mean_authority_plus_paired_ci_mixed_normalized": (
            maximum_mean),
        "maximum_mean_authority_plus_paired_ci_mixed_normalized": maximum_mean,
        "maximum_all_patch_mean_authority_plus_paired_ci_mixed_normalized_diagnostic": (
            maximum_all_patch_mean),
        "support_threshold_sensitivity": {
            "thresholds": sensitivity,
            "primary_threshold": float(DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
            "primary_gate_conclusion": primary_pass,
            "gate_conclusion_stable_over_predeclared_thresholds": bool(all(
                record["pass"] == primary_pass for record in sensitivity)),
        },
        "pass": primary_pass,
    }


def _paired_response_refinement(
        coarse_entry, fine_entry, direct, source, *, dx_m):
    coarse_authority = coarse_entry["authority"]
    fine_authority = fine_entry["authority"]
    coarse_replicates = coarse_entry["replicates"]
    fine_replicates = fine_entry["replicates"]
    if len(coarse_replicates) != len(fine_replicates):
        raise RuntimeError("nested response levels have different replicate counts")
    coarse_integrated = _extended_integrated_observables(coarse_authority, direct)
    fine_integrated = _extended_integrated_observables(fine_authority, direct)
    coarse_replicate_integrated = [
        _extended_integrated_observables(item, direct) for item in coarse_replicates]
    fine_replicate_integrated = [
        _extended_integrated_observables(item, direct) for item in fine_replicates]
    names = _response_field_names(coarse_authority, (fine_authority,))
    natural_scale = _natural_response_field_scales(
        names, direct, source, dx_m=dx_m)
    integrated = {}
    for name in sorted(coarse_integrated):
        if name.startswith("area_mean/radiosity/"):
            field_name = name.removeprefix("area_mean/")
            absolute = (
                GATES["integrated_absolute_displacement_dx"]
                * natural_scale[field_name])
        else:
            absolute = GATES["integrated_absolute_displacement_dx"] * dx_m
        integrated[name] = _paired_scalar_refinement_score(
            coarse_integrated[name], fine_integrated[name],
            [item[name] for item in coarse_replicate_integrated],
            [item[name] for item in fine_replicate_integrated],
            absolute_tolerance=absolute,
            relative_tolerance=GATES["integrated_relative_tolerance"])
    patch = {}
    for name in names:
        patch[name] = tuple(
            _paired_patch_refinement_score(
                coarse_authority["fields"][name], fine_authority["fields"][name],
                [item["fields"][name] for item in coarse_replicates],
                [item["fields"][name] for item in fine_replicates],
                direct, source, normalization=natural_scale[name],
                patch_scale_m=scale)
            for scale in PATCH_SCALES_M)
    gates = {
        "integrated_paired_refinement": all(
            item["pass"] for item in integrated.values()),
        "physical_patch_paired_refinement": all(
            item["pass"] for values in patch.values() for item in values),
    }
    return {
        "integrated_observables": integrated,
        "physical_patch_fields": patch,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def _paired_direction_score(level_responses, direct):
    output = {}
    for level in RAY_LEVELS:
        entry = level_responses[level]
        authority = {
            label: _integrated_observables(entry[label]["authority"], direct)[
                "sio2_mean_recession_m"]
            for label in ("r17", "r19")
        }
        replicates = []
        for index in range(len(REPLICATE_SEEDS)):
            r17 = _integrated_observables(
                entry["r17"]["replicates"][index], direct)[
                    "sio2_mean_recession_m"]
            r19 = _integrated_observables(
                entry["r19"]["replicates"][index], direct)[
                    "sio2_mean_recession_m"]
            replicates.append(r19 - r17)
        difference = authority["r19"] - authority["r17"]
        values = np.asarray(replicates)
        critical = float(student_t.ppf(0.975, len(values) - 1))
        mean = float(np.mean(values))
        half = critical * float(np.std(values, ddof=1)) / np.sqrt(len(values))
        output[level] = {
            "authority_r19_minus_r17_m": difference,
            "replicate_mean_r19_minus_r17_m": mean,
            "replicate_95pct_confidence_half_width_m": half,
            "direction": "r19_higher" if difference > 0.0 else (
                "r19_lower" if difference < 0.0 else "equal"),
            "replicate_mean_direction": "r19_higher" if mean > 0.0 else (
                "r19_lower" if mean < 0.0 else "equal"),
            "confidence_excludes_zero": bool(abs(mean) > half),
        }
    diagnostic_directions_agree = (
        output[DIAGNOSTIC_LEVEL_PAIR[0]]["direction"]
        == output[DIAGNOSTIC_LEVEL_PAIR[1]]["direction"] != "equal")
    authoritative_directions_agree = (
        output[AUTHORITATIVE_LEVEL_PAIR[0]]["direction"]
        == output[AUTHORITATIVE_LEVEL_PAIR[1]]["direction"] != "equal")
    gates = {
        "authoritative_16_to_32_direction_agreement": authoritative_directions_agree,
        "level32_paired_confidence_excludes_zero": output[FINAL_STAGE_A_LEVEL][
            "confidence_excludes_zero"],
        "level32_authority_and_replicate_mean_sign_agree": bool(
            output[FINAL_STAGE_A_LEVEL]["direction"]
            == output[FINAL_STAGE_A_LEVEL]["replicate_mean_direction"] != "equal"),
    }
    return {
        "levels": output,
        "diagnostic_8_to_16_direction_agreement": diagnostic_directions_agree,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def _build_level_responses(ensemble, direct, source, *, level, horizon_s):
    output = {}
    for label in ("r17", "r19"):
        mechanism = build_krueger_2024_material_router_3d(**PARAMETERS[label])
        common = {
            "schema": SCHEMA,
            "checkpoint_sha256": _sha256(source["checkpoint_path"]),
            "direct_surface_flux_sha256": surface_flux_sha256(
                direct["direct_surface_fluxes"]),
            "ensemble_sha256": ensemble.sha256,
            "ray_level": int(level),
            "parameter_label": label,
            "frozen_response_horizon_s": float(horizon_s),
            "frozen_response_horizon_fraction_of_next_step": (
                FROZEN_RESPONSE_HORIZON_FRACTION),
            "operator_epoch": "triangle-area-cellwise-certified-rqmc-8-16-32-v2",
        }
        authority_identity = {**common, "factor_member": "mean"}
        authority = _instantaneous_response(
            direct, source, ensemble.mean_form_factors, mechanism,
            horizon_s=horizon_s, identity=authority_identity)
        replicates = []
        for seed, factors in zip(
                ensemble.replicate_seeds, ensemble.replicate_form_factors):
            identity = {**common, "factor_member": f"seed:{seed}"}
            replicates.append(_instantaneous_response(
                direct, source, factors, mechanism,
                horizon_s=horizon_s, identity=identity))
        output[label] = {"authority": authority, "replicates": replicates}
    return output


def _stage_b_job_specs(seeds=REPLICATE_SEEDS):
    """Refuse Stage B during the v3 evidence-selection campaign."""
    del seeds
    raise RuntimeError(
        "Stage B is structurally held in the v3 8/16/32 Stage-A campaign")


def _response_from_integration(initial_state, result):
    fields = {
        f"state_increment/{name}": (
            np.asarray(result.state.fields[name], dtype=float)
            - np.asarray(initial_state.fields[name], dtype=float))
        for name in sorted(initial_state.fields)
    }
    for kind, attribute in _EXCHANGE_ATTRIBUTES.items():
        for name, value in sorted(getattr(result.material_exchange, attribute).items()):
            fields[f"exchange/{kind}/{name}"] = np.asarray(value, dtype=float)
    fields["surface/integrated_recession_m"] = np.asarray(
        result.integrated_recession_m, dtype=float)
    fields["surface/integrated_growth_m"] = np.asarray(
        result.integrated_growth_m, dtype=float)
    return {
        "fields": fields,
        "maximum_radiosity_relative_balance_error": float(
            result.diagnostics.maximum_radiosity_relative_balance_error),
        "maximum_material_ledger_residual_units_m2": float(
            _maximum_ledger_residual(result.material_exchange)),
        "maximum_gross_displacement_m": float(np.max(
            fields["surface/integrated_recession_m"]
            + fields["surface/integrated_growth_m"], initial=0.0)),
        "maximum_embedded_relative_error": float(
            result.diagnostics.maximum_accepted_embedded_relative_error),
        "maximum_local_reaction_probability_change": float(
            result.diagnostics.maximum_accepted_path_probability_change),
        "rejected_trial_exchange_contribution_is_zero": bool(
            result.diagnostics.rejected_trial_exchange_contribution_is_zero),
        "accepted_trial_count": int(result.diagnostics.accepted_trial_count),
        "rejected_trial_count": int(result.diagnostics.rejected_trial_count),
        "radiosity_solve_count": int(result.diagnostics.radiosity_solve_count),
        "minimum_accepted_chemistry_substep_s": float(
            result.diagnostics.minimum_accepted_chemistry_substep_s),
        "maximum_accepted_chemistry_substep_s": float(
            result.diagnostics.maximum_accepted_chemistry_substep_s),
    }


def _write_response_job(path, job_identity, response):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    field_names = tuple(sorted(response["fields"]))
    arrays = {
        f"field_{index:03d}": np.asarray(response["fields"][name])
        for index, name in enumerate(field_names)
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)
    metadata = {
        "job_identity": job_identity,
        "job_identity_sha256": _canonical_sha256(job_identity),
        "field_names": field_names,
        "scalars": {key: value for key, value in response.items() if key != "fields"},
        "npz_sha256": _sha256(path),
    }
    _write_json_atomic(path.with_suffix(".json"), metadata)


def _read_response_job(path, job_identity):
    path = Path(path)
    metadata_path = path.with_suffix(".json")
    if not path.exists() or not metadata_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (metadata.get("job_identity_sha256") != _canonical_sha256(job_identity)
            or metadata.get("npz_sha256") != _sha256(path)):
        return None
    with np.load(path, allow_pickle=False) as data:
        fields = {
            name: np.asarray(data[f"field_{index:03d}"]).copy()
            for index, name in enumerate(metadata["field_names"])
        }
    return {"fields": fields, **metadata["scalars"]}


def _integration_job_contract(response, *, tight, dx_m):
    embedded_limit = (
        GATES["tight_maximum_embedded_relative_error"]
        if tight else GATES["maximum_embedded_relative_error"])
    probability_limit = (
        GATES["tight_maximum_local_reaction_probability_change"]
        if tight else GATES["maximum_local_reaction_probability_change"])
    gates = {
        "radiosity_balance": bool(
            response["maximum_radiosity_relative_balance_error"]
            <= GATES["maximum_radiosity_relative_balance_error"]),
        "exact_material_ledger": bool(
            response["maximum_material_ledger_residual_units_m2"] == 0.0),
        "embedded_error": bool(
            response["maximum_embedded_relative_error"] <= embedded_limit),
        "local_probability_safety": bool(
            response["maximum_local_reaction_probability_change"]
            <= probability_limit),
        "rejected_trials_contribute_zero": bool(
            response["rejected_trial_exchange_contribution_is_zero"]),
        "gross_displacement": bool(
            response["maximum_gross_displacement_m"] / dx_m
            <= GATES["maximum_gross_displacement_dx"]),
        "accepted_step_contract": bool(
            response["accepted_trial_count"] > 0
            and response["radiosity_solve_count"]
            >= 2 * response["accepted_trial_count"]
            and response["minimum_accepted_chemistry_substep_s"]
            >= GATES["minimum_chemistry_substep_s"]),
    }
    return {"gates": gates, "all_gates_pass": all(gates.values())}


def _factor_for_job(job, ensembles):
    ensemble = ensembles[int(job["ray_level"])]
    if job["factor_member"] == "mean":
        return ensemble.mean_form_factors, ensemble.sha256
    seed = int(job["factor_member"].split(":", 1)[1])
    index = ensemble.replicate_seeds.index(seed)
    return ensemble.replicate_form_factors[index], ensemble.sha256


def _run_stage_b(
        args, payload, direct, source, ensembles, *, horizon_s, started):
    raise RuntimeError(
        "Stage B is structurally held in the v3 8/16/32 Stage-A campaign")
    job_directory = Path(args.output).parent / "stage_b_jobs"
    responses = {}
    for job in _stage_b_job_specs(args.replicate_seeds):
        factors, ensemble_sha = _factor_for_job(job, ensembles)
        identity = {
            "schema": SCHEMA,
            "stage": "stage_b",
            "checkpoint_sha256": _sha256(source["checkpoint_path"]),
            "direct_surface_flux_sha256": surface_flux_sha256(
                direct["direct_surface_fluxes"]),
            "ensemble_sha256": ensemble_sha,
            "horizon_s": float(horizon_s),
            "parameters": PARAMETERS[job["parameter_label"]],
            "operator": OPERATOR,
            "gates": GATES,
            "physical_patch_operator": payload["physical_patch_operator"],
            "source_epoch_sha256": _canonical_sha256(
                payload["provenance"]["source"]),
            "integration_controls": {
                "maximum_embedded_relative_error": (
                    GATES["tight_maximum_embedded_relative_error"]
                    if job["tight"] else GATES["maximum_embedded_relative_error"]),
                "maximum_local_reaction_probability_change": (
                    GATES["tight_maximum_local_reaction_probability_change"]
                    if job["tight"]
                    else GATES["maximum_local_reaction_probability_change"]),
                "minimum_chemistry_substep_s": GATES[
                    "minimum_chemistry_substep_s"],
            },
            **job,
        }
        artifact = job_directory / f"{job['name']}.npz"
        response = _read_response_job(artifact, identity)
        if response is None:
            remaining = float(args.maximum_total_wall_s) - (perf_counter() - started)
            with _hard_deadline(min(float(args.maximum_endpoint_job_wall_s), remaining)):
                mechanism = build_krueger_2024_material_router_3d(
                    **PARAMETERS[job["parameter_label"]])
                cache = SurfaceRadiosityOperatorCache3D(
                    direct["direct_surface_fluxes"], factors,
                    direct["face_area_m2"], direct["active_face_index"],
                    direct["face_material_id"], direct["species_role"], identity,
                    relative_tolerance=OPERATOR["radiosity_relative_tolerance"],
                    maximum_iterations=OPERATOR["radiosity_maximum_iterations"])
                result = integrate_surface_radiosity_chemistry_3d(
                    cache, mechanism, source["state"], horizon_s,
                    operator_identity_payload=identity,
                    maximum_embedded_relative_error=(
                        GATES["tight_maximum_embedded_relative_error"]
                        if job["tight"] else GATES["maximum_embedded_relative_error"]),
                    maximum_local_reaction_probability_change=(
                        GATES["tight_maximum_local_reaction_probability_change"]
                        if job["tight"]
                        else GATES["maximum_local_reaction_probability_change"]),
                    minimum_chemistry_substep_s=GATES[
                        "minimum_chemistry_substep_s"])
                response = _response_from_integration(source["state"], result)
            _write_response_job(artifact, identity, response)
        responses[job["name"]] = response
        payload["stage_b"]["completed_jobs"] = sorted(responses)
        _write_json_atomic(args.output, payload)
    # The 16-ray replicate ensemble is the uncertainty authority.  The 8-ray means and tight
    # means are retained as independent refinement/schedule checks in the final receipt.
    remaining = float(args.maximum_total_wall_s) - (perf_counter() - started)
    with _hard_deadline(remaining):
        level16 = {}
        nested = {}
        tight = {}
        dx_m = source["geometry"].dx * source["geometry"].mesh_length_unit_m
        for label in ("r17", "r19"):
            mean16 = responses[f"{label}_level16_mean_nominal"]
            replicates = [
                responses[f"{label}_level16_seed{seed}_nominal"]
                for seed in args.replicate_seeds]
            level16[label] = _score_response_level(
                mean16, replicates, direct, source, dx_m=dx_m)
            nested[label] = _compare_response_refinement(
                responses[f"{label}_level8_mean_nominal"], mean16,
                direct, source, dx_m=dx_m)
            tight[label] = _compare_response_refinement(
                mean16, responses[f"{label}_level16_mean_tight"],
                direct, source, dx_m=dx_m)
        job_contracts = {
            job["name"]: _integration_job_contract(
                responses[job["name"]], tight=job["tight"], dx_m=dx_m)
            for job in _stage_b_job_specs(args.replicate_seeds)
        }
        gates = {
            "level16_replicated_uncertainty": all(
                item["all_gates_pass"] for item in level16.values()),
            "nested_8_to_16_all_fields": all(
                item["all_gates_pass"] for item in nested.values()),
            "schedule_tightening_all_fields": all(
                item["all_gates_pass"] for item in tight.values()),
            "all_22_jobs_complete": len(responses) == 22,
            "all_job_numerical_contracts": all(
                item["all_gates_pass"] for item in job_contracts.values()),
        }
        payload["stage_b"].update({
            "level16_scores": level16,
            "nested_scores": nested,
            "schedule_tightening_scores": tight,
            "job_numerical_contracts": job_contracts,
            "gates": gates,
            "all_gates_pass": all(gates.values()),
        })
    payload["status"] = "pass" if all(gates.values()) else "bounded_precision_hold"
    return payload


def _primary_patch_support_contract_pass(*parameter_collections):
    records = []
    for collection in parameter_collections:
        for parameter in collection.values():
            for values in parameter.get("physical_patch_fields", {}).values():
                records.extend(values)
    return bool(records and all(
        float(record.get("minimum_mean_support_fraction", -1.0))
        == float(DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION)
        and int(record.get("eligible_mean_patch_count", 0)) > 0
        and "support_threshold_sensitivity" in record
        for record in records))


def _stage_a_gate(
        level_summary, authoritative_nested, authoritative_paired_nested,
        paired_direction, row_closure, nested_sampling, *, claimed_feature_extent_m):
    level32_parameters = level_summary[FINAL_STAGE_A_LEVEL]["parameter_scores"]
    level32_maximum_gross_displacement_dx = max(
        parameter["maximum_gross_displacement_dx"]
        for parameter in level32_parameters.values())
    gates = {
        "exact_replicate_count": all(
            item["replicate_count"] == GATES["minimum_replicate_count"]
            for item in level_summary.values()),
        "row_closure": all(
            value <= GATES["maximum_form_factor_row_closure_error"]
            for value in row_closure.values()),
        "level32_replicated_uncertainty": all(
            parameter["all_gates_pass"]
            for parameter in level32_parameters.values()),
        "authoritative_nested_16_to_32_all_fields": all(
            parameter["all_gates_pass"]
            for parameter in authoritative_nested.values()),
        "authoritative_paired_nested_16_to_32_all_fields": all(
            parameter["all_gates_pass"]
            for parameter in authoritative_paired_nested.values()),
        "level32_paired_direction": paired_direction["all_gates_pass"],
        "exact_nested_sampling_extension_16_to_32": nested_sampling[
            "16_to_32"]["all_gates_pass"],
        "level32_direct_gross_motion_at_frozen_horizon": bool(
            level32_maximum_gross_displacement_dx
            <= GATES["maximum_gross_displacement_dx"]),
        "primary_0p10_patch_support_contract": (
            _primary_patch_support_contract_pass(
                level32_parameters, authoritative_nested,
                authoritative_paired_nested)),
        "physical_patch_scale_resolves_claimed_feature": bool(
            min(PATCH_SCALES_M) <= float(claimed_feature_extent_m)),
    }
    return gates


def _stage_a_patch_records(stage_a):
    records = []
    for level_key, level in stage_a.get("levels", {}).items():
        gating = int(level_key) == FINAL_STAGE_A_LEVEL
        for parameter in level.get("parameter_scores", {}).values():
            for values in parameter.get("physical_patch_fields", {}).values():
                records.extend((record, gating) for record in values)
    for collection_name, gating in (
            ("diagnostic_nested_8_to_16_all_fields", False),
            ("diagnostic_paired_nested_8_to_16_all_fields", False),
            ("authoritative_nested_16_to_32_all_fields", True),
            ("authoritative_paired_nested_16_to_32_all_fields", True)):
        for parameter in stage_a.get(collection_name, {}).values():
            for values in parameter.get("physical_patch_fields", {}).values():
                records.extend((record, gating) for record in values)
    return records


def _stage_a_patch_receipts_complete(stage_a):
    records = _stage_a_patch_records(stage_a)
    required = {
        "patch_scale_m", "patch_count", "minimum_mean_support_fraction",
        "eligible_mean_patch_count", "excluded_mean_patch_count",
        "excluded_mean_surface_area_fraction",
        "excluded_mean_projected_support_fraction",
        "represented_nominal_projected_area_m2",
        "support_threshold_sensitivity", "pass",
    }
    if not records:
        return False
    expected_thresholds = tuple(float(value) for value in (
        PATCH_SUPPORT_SENSITIVITY_THRESHOLDS))
    for record, gating in records:
        if not isinstance(record, Mapping) or not required.issubset(record):
            return False
        if gating and record["pass"] is not True:
            return False
        if (float(record["minimum_mean_support_fraction"])
                != float(DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION)):
            return False
        count = int(record["patch_count"])
        eligible = int(record["eligible_mean_patch_count"])
        excluded = int(record["excluded_mean_patch_count"])
        if count <= 0 or eligible <= 0 or excluded < 0 or eligible + excluded != count:
            return False
        if any(
                not 0.0 <= float(record[name]) <= 1.0
                for name in (
                    "excluded_mean_surface_area_fraction",
                    "excluded_mean_projected_support_fraction")):
            return False
        nominal = record["represented_nominal_projected_area_m2"]
        if (not isinstance(nominal, Mapping)
                or any(float(nominal.get(name, 0.0)) <= 0.0
                       for name in ("minimum", "maximum", "total"))):
            return False
        sensitivity = record["support_threshold_sensitivity"]
        local = sensitivity.get("thresholds", []) if isinstance(
            sensitivity, Mapping) else []
        observed_thresholds = tuple(
            float(item.get("minimum_mean_support_fraction", -1.0))
            for item in local if isinstance(item, Mapping))
        if (observed_thresholds != expected_thresholds
                or float(sensitivity.get("primary_threshold", -1.0))
                != float(DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION)
                or bool(sensitivity.get("primary_gate_conclusion"))
                != bool(record["pass"])
                or "gate_conclusion_stable_over_predeclared_thresholds"
                not in sensitivity):
            return False
    return True


def _stage_a_v3_receipts_complete(stage_a):
    """Reject pre-v3 or partially populated artifacts even if labels were edited."""
    if not isinstance(stage_a, Mapping):
        return False
    levels = stage_a.get("levels", {})
    try:
        observed_levels = {int(level) for level in levels}
    except (TypeError, ValueError):
        return False
    if observed_levels != set(RAY_LEVELS):
        return False
    if stage_a.get("diagnostic_8_to_16_gating") is not False:
        return False
    for name in (
            "diagnostic_nested_8_to_16_all_fields",
            "diagnostic_paired_nested_8_to_16_all_fields",
            "authoritative_nested_16_to_32_all_fields",
            "authoritative_paired_nested_16_to_32_all_fields"):
        collection = stage_a.get(name)
        if not isinstance(collection, Mapping) or set(collection) != {"r17", "r19"}:
            return False
    nested = stage_a.get("nested_sampling_extension", {})
    if (not isinstance(nested, Mapping)
            or set(nested) != {"8_to_16", "16_to_32"}
            or not isinstance(nested["16_to_32"], Mapping)
            or nested["16_to_32"].get("all_gates_pass") is not True):
        return False
    candidate = stage_a.get("level32_final_candidate", {})
    if (not isinstance(candidate, Mapping)
            or candidate.get("ray_level") != FINAL_STAGE_A_LEVEL
            or candidate.get("directly_evaluated_not_linearly_projected") is not True
            or candidate.get("pass") is not True):
        return False
    gates = stage_a.get("gates", {})
    if (not isinstance(gates, Mapping) or not gates
            or not all(value is True for value in gates.values())):
        return False
    return True


def _load_stage_a_authorization(path, source, args):
    path = Path(path)
    audit = json.loads(path.read_text(encoding="utf-8"))
    if audit.get("schema") != SCHEMA or audit.get("status") != STAGE_A_PASS_STATUS:
        raise ValueError("Stage B requires a completed passing Stage-A artifact")
    if audit["checkpoint"]["checkpoint_sha256"] != _sha256(source["checkpoint_path"]):
        raise ValueError("Stage-A artifact uses another checkpoint")
    if tuple(audit["sampling"]["replicate_seeds"]) != tuple(args.replicate_seeds):
        raise ValueError("Stage-A artifact uses another replicate seed contract")
    if tuple(audit["sampling"]["ray_levels"]) != RAY_LEVELS:
        raise ValueError("Stage-A artifact uses another nested ray-level contract")
    if int(audit["direct_transport"]["transport_seed"]) != int(args.transport_seed):
        raise ValueError("Stage-A artifact uses another direct-transport seed")
    if (audit.get("operator") != OPERATOR or audit.get("gates") != GATES
            or audit.get("provenance", {}).get("source")
            != _hash_manifest(SOURCE_PATHS)
            or audit.get("provenance", {}).get("base_inputs")
            != _hash_manifest(BASE_INPUT_PATHS)):
        raise ValueError("Stage-A artifact uses another operator, gate, or source epoch")
    if audit.get("physical_patch_operator") != _physical_patch_operator_contract(source):
        raise ValueError("Stage-A artifact uses another physical-patch operator")
    firewall = audit.get("data_firewall", {})
    if (firewall.get("boundary_case") != "base"
            or firewall.get("held_out_observations_loaded") is not False
            or firewall.get("held_out_transfer_boundary_constructed") is not False
            or audit.get("stage_a", {}).get("all_gates_pass") is not True
            or not _stage_a_v3_receipts_complete(audit.get("stage_a", {}))
            or not _stage_a_patch_receipts_complete(audit.get("stage_a", {}))):
        raise ValueError("Stage-A artifact does not carry a passing sealed-base contract")
    expected_horizon = (
        FROZEN_RESPONSE_HORIZON_FRACTION
        * float(source["checkpoint_metadata"]["next_step_duration_s"]))
    candidate_horizon = audit["stage_a"]["level32_final_candidate"].get(
        "frozen_response_horizon_s")
    if (not np.isclose(
            float(audit["checkpoint"]["frozen_response_horizon_s"]),
            expected_horizon, rtol=0.0, atol=0.0)
            or candidate_horizon is None
            or not np.isclose(
                float(candidate_horizon), expected_horizon, rtol=0.0, atol=0.0)):
        raise ValueError("Stage-A artifact uses another chemistry horizon")
    return {"path": str(path), "sha256": _sha256(path), "audit": audit}


def run(args):
    if bool(getattr(args, "plan_only", False)):
        raise ValueError("plan-only receipt cannot enter the execution path")
    if (getattr(args, "stage", None) != "stage_a"
            or bool(getattr(args, "authorize_stage_b", False))
            or getattr(args, "stage_a_audit", None) is not None):
        raise ValueError(
            "v3 permits Stage A only; Stage B is structurally held regardless of result")
    started = perf_counter()
    source = _load_sealed_base_source(args.r19_source)
    physical_patch_operator = _physical_patch_operator_contract(source)
    parameters = _validate_parameter_provenance(args.r17_source, source)
    if source["checkpoint_metadata"].get("physical_time_s") != 60.0:
        raise ValueError("replicated closure requires the completed 60 s base checkpoint")
    next_step_s = float(source["checkpoint_metadata"]["next_step_duration_s"])
    horizon_s = FROZEN_RESPONSE_HORIZON_FRACTION * next_step_s
    snapshot = _snapshot_inputs(source["geometry"], source["state"])
    payload = {
        "schema": SCHEMA,
        "status": "running",
        "stage": args.stage,
        "scientific_scope": (
            "base-checkpoint form-factor sampling closure only; no profile motion, fitting, "
            "held-out observation, or transfer-boundary construction"),
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
            "non_base_refused_before_checkpoint_read": True,
        },
        "operator": OPERATOR,
        "sampling": {
            "ray_levels": RAY_LEVELS,
            "diagnostic_level_pair": DIAGNOSTIC_LEVEL_PAIR,
            "authoritative_level_pair": AUTHORITATIVE_LEVEL_PAIR,
            "final_stage_a_level": FINAL_STAGE_A_LEVEL,
            "replicate_seeds": args.replicate_seeds,
            "replicate_count": len(args.replicate_seeds),
            "nested_common_scramble_across_levels": True,
            "independent_scrambles_within_level": True,
            "source_sampling": OPERATOR["source_sampling"],
            "hard_visibility": True,
        },
        "execution_budget": {
            "maximum_direct_transport_wall_s": args.maximum_direct_transport_wall_s,
            "maximum_form_factor_replicate_wall_s": (
                args.maximum_form_factor_replicate_wall_s),
            "maximum_endpoint_job_wall_s": args.maximum_endpoint_job_wall_s,
            "maximum_total_wall_s": args.maximum_total_wall_s,
            "maximum_process_count": args.maximum_process_count,
            "actual_process_count": 1,
            "automatic_ray_escalation_above_32": False,
        },
        "gates": GATES,
        "physical_patch_operator": physical_patch_operator,
        "checkpoint": {
            "audit_sha256": _sha256(source["audit_path"]),
            "checkpoint_sha256": _sha256(source["checkpoint_path"]),
            "metadata": _jsonable(source["checkpoint_metadata"]),
            "metrics": _jsonable(source["metrics"]),
            "next_profile_step_s": next_step_s,
            "frozen_response_horizon_fraction": FROZEN_RESPONSE_HORIZON_FRACTION,
            "frozen_response_horizon_s": horizon_s,
        },
        "parameter_provenance": parameters,
        "provenance": {
            "source": _hash_manifest(SOURCE_PATHS),
            "base_inputs": _hash_manifest(BASE_INPUT_PATHS),
        },
        "stage_a": {},
        "stage_b": {
            "authorized": False,
            "campaign_status": "structurally_held_pending_separate_review",
            "automatic_entry": False,
            "job_count": 0,
            "completed_jobs": [],
        },
    }
    _write_json_atomic(args.output, payload)
    config = _production_config(source["config"], PARAMETERS["r19"])
    direct_cache_path = (
        Path(args.direct_transport_cache)
        if args.direct_transport_cache is not None
        else Path(args.output).parent / "direct_transport_cache.npz")
    direct_identity = _direct_transport_cache_identity(
        source, config, seed=args.transport_seed)
    try:
        direct_cache_hit = False
        direct_cache_miss_reason = None
        try:
            direct, direct_cache_metadata = _read_direct_transport_artifact(
                direct_cache_path, source, identity=direct_identity)
            direct_cache_hit = True
        except (FileNotFoundError, ValueError, json.JSONDecodeError) as error:
            direct_cache_miss_reason = type(error).__name__
            remaining = args.maximum_total_wall_s - (perf_counter() - started)
            with _hard_deadline(min(args.maximum_direct_transport_wall_s, remaining)):
                direct = _capture_direct_transport_once(
                    source, config, seed=args.transport_seed)
            direct_cache_metadata = _write_direct_transport_artifact(
                direct_cache_path, direct, identity=direct_identity)
        if not _inputs_unchanged(snapshot, source["geometry"], source["state"]):
            raise RuntimeError("direct transport mutated the sealed checkpoint")
        payload["direct_transport"] = {
            "transport_seed": int(args.transport_seed),
            "surface_flux_sha256": surface_flux_sha256(
                direct["direct_surface_fluxes"]),
            "reported_evaluator_wall_s": direct["reported_wall_s"],
            "production_form_factor_call_count": direct[
                "production_form_factor_call_count"],
            "production_form_factor_trace_elided": direct[
                "production_form_factor_trace_elided"],
            "input_checkpoint_unchanged": True,
            "cache": {
                "path_name": direct_cache_path.name,
                "cache_hit": direct_cache_hit,
                "cache_miss_reason": direct_cache_miss_reason,
                "artifact_npz_sha256": direct_cache_metadata["npz_sha256"],
                "artifact_metadata_sha256": _sha256(
                    direct_cache_path.with_suffix(".json")),
                "identity_sha256": direct_cache_metadata["identity_sha256"],
                "validated_against_current_checkpoint_mesh": True,
            },
        }
        ensembles = {}
        receipts = {}
        row_closure = {}
        factor_artifacts = {}
        nested_sampling = {}
        artifact_directory = Path(args.output).parent / "form_factor_ensembles"
        for level_index, level in enumerate(RAY_LEVELS):
            ensemble, local_receipts, closure = _estimate_replicated_level(
                direct, source, rays_per_face=level,
                seeds=args.replicate_seeds,
                maximum_replicate_wall_s=args.maximum_form_factor_replicate_wall_s,
                total_deadline_s=args.maximum_total_wall_s, started=started)
            if closure > GATES["maximum_form_factor_row_closure_error"]:
                raise RuntimeError(
                    f"level{level} form-factor row closure failed before downstream scoring")
            ensembles[level] = ensemble
            receipts[level] = local_receipts
            row_closure[level] = closure
            if level_index > 0:
                previous = RAY_LEVELS[level_index - 1]
                label = f"{previous}_to_{level}"
                diagnostic = _nested_sampling_extension_diagnostic(
                    ensembles[previous], ensemble)
                nested_sampling[label] = diagnostic
                if not diagnostic["all_gates_pass"]:
                    raise RuntimeError(
                        f"exact nested Sobol extension failed for {previous}->{level}")
            artifact_path = artifact_directory / f"level{level}.npz"
            metadata = _write_form_factor_ensemble(
                artifact_path, ensemble,
                checkpoint_sha256=_sha256(source["checkpoint_path"]),
                ray_level=level)
            factor_artifacts[level] = {
                "npz_relative_to_audit": str(
                    artifact_path.relative_to(Path(args.output).parent)),
                "npz_sha256": metadata["npz_sha256"],
                "metadata_sha256": _sha256(artifact_path.with_suffix(".json")),
            }
        if not _inputs_unchanged(snapshot, source["geometry"], source["state"]):
            raise RuntimeError("form-factor construction mutated the sealed checkpoint")
    except EvaluationDeadlineExceeded as error:
        payload["status"] = "bounded_timeout"
        payload["timeout"] = {
            "reason": str(error),
            "physics_conclusion_permitted": False,
        }
        payload["total_wall_time_s"] = perf_counter() - started
        _write_json_atomic(args.output, payload)
        return payload
    except (RuntimeError, ValueError) as error:
        payload["status"] = "authority_refusal"
        payload["refusal"] = {
            "type": type(error).__name__,
            "reason": str(error),
            "physics_conclusion_permitted": False,
        }
        payload["total_wall_time_s"] = perf_counter() - started
        _write_json_atomic(args.output, payload)
        return payload

    payload["form_factor_ensembles"] = {
        str(level): {
            "sha256": ensembles[level].sha256,
            "replicate_count": len(ensembles[level].replicate_form_factors),
            "rays_per_replicate": ensembles[level].rays_per_replicate,
            "total_rays_per_face": ensembles[level].total_rays_per_face,
            "maximum_row_closure_error": row_closure[level],
            "mean_raw_reciprocity": _reciprocity_summary(
                ensembles[level].mean_reciprocity),
            "replicate_raw_reciprocity": [
                _reciprocity_summary(item)
                for item in ensembles[level].replicate_reciprocity],
            "maximum_visibility_wrap_count": (
                max(item.maximum_wrap_count for item in receipts[level])
                if receipts[level] else None),
            "visibility_float64_evaluated_count": (
                sum(item.float64_evaluated_count for item in receipts[level])
                if receipts[level] else None),
            "launch_inset_count": (
                sum(item.launch_inset_count for item in receipts[level])
                if receipts[level] else None),
            "centroid_limit_count": (
                sum(item.centroid_limit_count for item in receipts[level])
                if receipts[level] else None),
            "artifact": factor_artifacts[level],
        }
        for level in RAY_LEVELS
    }
    try:
        remaining = args.maximum_total_wall_s - (perf_counter() - started)
        with _hard_deadline(remaining):
            level_responses = {}
            level_summary = {}
            dx_m = (
                source["geometry"].dx
                * source["geometry"].mesh_length_unit_m)
            for level in RAY_LEVELS:
                responses = _build_level_responses(
                    ensembles[level], direct, source, level=level,
                    horizon_s=horizon_s)
                level_responses[level] = responses
                scores = {
                    label: _score_response_level(
                        responses[label]["authority"],
                        responses[label]["replicates"], direct, source,
                        dx_m=dx_m)
                    for label in ("r17", "r19")
                }
                level_summary[level] = {
                    "replicate_count": len(args.replicate_seeds),
                    "parameter_scores": scores,
                }
            diagnostic_nested = {
                label: _compare_response_refinement(
                    level_responses[DIAGNOSTIC_LEVEL_PAIR[0]][label]["authority"],
                    level_responses[DIAGNOSTIC_LEVEL_PAIR[1]][label]["authority"],
                    direct, source, dx_m=dx_m)
                for label in ("r17", "r19")
            }
            diagnostic_paired_nested = {
                label: _paired_response_refinement(
                    level_responses[DIAGNOSTIC_LEVEL_PAIR[0]][label],
                    level_responses[DIAGNOSTIC_LEVEL_PAIR[1]][label],
                    direct, source, dx_m=dx_m)
                for label in ("r17", "r19")
            }
            authoritative_nested = {
                label: _compare_response_refinement(
                    level_responses[AUTHORITATIVE_LEVEL_PAIR[0]][label]["authority"],
                    level_responses[AUTHORITATIVE_LEVEL_PAIR[1]][label]["authority"],
                    direct, source, dx_m=dx_m)
                for label in ("r17", "r19")
            }
            authoritative_paired_nested = {
                label: _paired_response_refinement(
                    level_responses[AUTHORITATIVE_LEVEL_PAIR[0]][label],
                    level_responses[AUTHORITATIVE_LEVEL_PAIR[1]][label],
                    direct, source, dx_m=dx_m)
                for label in ("r17", "r19")
            }
            paired = _paired_direction_score(level_responses, direct)
            level32_motion_by_parameter = {
                label: float(level_summary[FINAL_STAGE_A_LEVEL][
                    "parameter_scores"][label]["maximum_gross_displacement_dx"])
                for label in ("r17", "r19")
            }
            maximum_level32_motion = max(level32_motion_by_parameter.values())
            level32_final_candidate = {
                "ray_level": FINAL_STAGE_A_LEVEL,
                "frozen_response_horizon_s": horizon_s,
                "horizon_fraction_of_next_step": FROZEN_RESPONSE_HORIZON_FRACTION,
                "directly_evaluated_not_linearly_projected": True,
                "parameter_maximum_gross_displacement_dx": (
                    level32_motion_by_parameter),
                "maximum_gross_displacement_dx": maximum_level32_motion,
                "maximum_gross_displacement_dx_limit": GATES[
                    "maximum_gross_displacement_dx"],
                "pass": bool(
                    maximum_level32_motion
                    <= GATES["maximum_gross_displacement_dx"]),
            }
            claimed_feature_extent_m = (
                float(source["metrics"]["mask_opening_nm"]) * 1.0e-9)
            if (not np.isfinite(claimed_feature_extent_m)
                    or claimed_feature_extent_m <= 0.0):
                raise RuntimeError(
                    "checkpoint has no positive feature extent for patch-scale gate")
            stage_a_gates = _stage_a_gate(
                level_summary, authoritative_nested,
                authoritative_paired_nested, paired, row_closure,
                nested_sampling,
                claimed_feature_extent_m=claimed_feature_extent_m)
            payload["stage_a"] = {
                "levels": level_summary,
                "diagnostic_8_to_16_gating": False,
                "diagnostic_nested_8_to_16_all_fields": diagnostic_nested,
                "diagnostic_paired_nested_8_to_16_all_fields": (
                    diagnostic_paired_nested),
                "authoritative_nested_16_to_32_all_fields": authoritative_nested,
                "authoritative_paired_nested_16_to_32_all_fields": (
                    authoritative_paired_nested),
                "paired_r19_minus_r17": paired,
                "nested_sampling_extension": nested_sampling,
                "level32_final_candidate": level32_final_candidate,
                "physical_patch_scale_contract": {
                    "operator": physical_patch_operator,
                    "claimed_feature_extent_m": claimed_feature_extent_m,
                    "reported_patch_scales_m": PATCH_SCALES_M,
                    "at_least_one_scale_no_larger_than_feature": bool(
                        min(PATCH_SCALES_M) <= claimed_feature_extent_m),
                },
                "gates": stage_a_gates,
                "all_gates_pass": all(stage_a_gates.values()),
                "failure_action": (
                    "bounded_precision_hold; no automatic level above 32 and no Stage B; "
                    "return for reviewed evidence-selection decision"),
            }
            payload["status"] = (
                STAGE_A_PASS_STATUS if all(stage_a_gates.values())
                else "bounded_precision_hold")
    except EvaluationDeadlineExceeded as error:
        payload["status"] = "bounded_timeout"
        payload["timeout"] = {
            "reason": str(error), "physics_conclusion_permitted": False}
    except (RuntimeError, ValueError) as error:
        payload["status"] = "authority_refusal"
        payload["refusal"] = {
            "type": type(error).__name__,
            "reason": str(error),
            "physics_conclusion_permitted": False,
        }
    payload["total_wall_time_s"] = perf_counter() - started
    payload["input_checkpoint_unchanged"] = _inputs_unchanged(
        snapshot, source["geometry"], source["state"])
    if not payload["input_checkpoint_unchanged"]:
        raise RuntimeError("replicated closure mutated the sealed checkpoint")
    _write_json_atomic(args.output, payload)
    print(json.dumps({
        "status": payload["status"],
        "stage": args.stage,
        "output": str(args.output),
        "total_wall_time_s": payload["total_wall_time_s"],
    }, indent=2, sort_keys=True))
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plan-only", action="store_true",
        help="print the sealed Stage-A plan without checkpoint or transport work")
    parser.add_argument(
        "--r17-source",
        default=(ROOT / "results" / "krueger_2024_base_calibration_r17"
                 / "axisym_candidate"))
    parser.add_argument(
        "--r19-source",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "remote_artifacts"))
    parser.add_argument(
        "--output",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "replicated_form_factor_closure" / "audit.json"))
    parser.add_argument("--stage", choices=("stage_a", "stage_b"), default="stage_a")
    parser.add_argument("--authorize-stage-b", action="store_true")
    parser.add_argument("--stage-a-audit")
    parser.add_argument("--transport-seed", type=int, default=241)
    parser.add_argument(
        "--direct-transport-cache",
        help=("persistent hash-bound pre-radiosity transport artifact; defaults "
              "to direct_transport_cache.npz beside --output"))
    parser.add_argument(
        "--replicate-seeds", type=int, nargs="+", default=REPLICATE_SEEDS)
    parser.add_argument(
        "--maximum-direct-transport-wall-s", type=float,
        default=DEFAULT_BUDGETS["maximum_direct_transport_wall_s"])
    parser.add_argument(
        "--maximum-form-factor-replicate-wall-s", type=float,
        default=DEFAULT_BUDGETS["maximum_form_factor_replicate_wall_s"])
    parser.add_argument(
        "--maximum-endpoint-job-wall-s", type=float,
        default=DEFAULT_BUDGETS["maximum_endpoint_job_wall_s"])
    parser.add_argument(
        "--maximum-total-wall-s", type=float,
        default=None)
    parser.add_argument(
        "--maximum-process-count", type=int,
        default=DEFAULT_BUDGETS["maximum_process_count"])
    args = parser.parse_args(argv)
    args.replicate_seeds = tuple(int(value) for value in args.replicate_seeds)
    if (len(args.replicate_seeds) != GATES["minimum_replicate_count"]
            or len(set(args.replicate_seeds)) != len(args.replicate_seeds)):
        parser.error("exactly eight distinct replicate seeds are required")
    try:
        args.maximum_direct_transport_wall_s = _validate_budget(
            "maximum_direct_transport_wall_s",
            args.maximum_direct_transport_wall_s,
            DEFAULT_BUDGETS["maximum_direct_transport_wall_s"])
        args.maximum_form_factor_replicate_wall_s = _validate_budget(
            "maximum_form_factor_replicate_wall_s",
            args.maximum_form_factor_replicate_wall_s,
            DEFAULT_BUDGETS["maximum_form_factor_replicate_wall_s"])
        args.maximum_endpoint_job_wall_s = _validate_budget(
            "maximum_endpoint_job_wall_s", args.maximum_endpoint_job_wall_s,
            DEFAULT_BUDGETS["maximum_endpoint_job_wall_s"])
        total_ceiling = DEFAULT_BUDGETS["maximum_stage_a_total_wall_s"]
        if args.maximum_total_wall_s is None:
            args.maximum_total_wall_s = total_ceiling
        args.maximum_total_wall_s = _validate_budget(
            "maximum_total_wall_s", args.maximum_total_wall_s, total_ceiling)
    except ValueError as error:
        parser.error(str(error))
    if (args.maximum_process_count < 1
            or args.maximum_process_count > DEFAULT_BUDGETS["maximum_process_count"]):
        parser.error("maximum_process_count must lie in [1, 4]")
    if args.stage != "stage_a":
        parser.error(
            "Stage B is structurally held in this v3 campaign regardless of Stage-A result")
    if args.authorize_stage_b or args.stage_a_audit is not None:
        parser.error("Stage-A runs cannot authorize or consume Stage B")
    return args


def _execution_plan_receipt(args):
    """Return a no-science receipt for reviewing the campaign before launch."""
    return {
        "schema": SCHEMA,
        "status": "plan_only_no_science_execution",
        "checkpoint_or_transport_execution_started": False,
        "stage": "stage_a",
        "stage_b": {
            "authorized": False,
            "campaign_status": "structurally_held_pending_separate_review",
        },
        "sampling": {
            "ray_levels": list(RAY_LEVELS),
            "diagnostic_level_pair": list(DIAGNOSTIC_LEVEL_PAIR),
            "diagnostic_pair_is_gating": False,
            "authoritative_level_pair": list(AUTHORITATIVE_LEVEL_PAIR),
            "final_stage_a_level": FINAL_STAGE_A_LEVEL,
            "replicate_seeds": list(args.replicate_seeds),
            "replicate_count": len(args.replicate_seeds),
            "exact_nested_sobol_extension": True,
            "source_sampling": "triangle_area",
            "visibility": "hard_cellwise_certified",
        },
        "response_horizon": {
            "definition": "checkpoint.next_step_duration_s / 1024",
            "fraction_of_next_step": FROZEN_RESPONSE_HORIZON_FRACTION,
            "shared_by_all_levels_and_replicates": True,
            "level32_direct_gross_motion_limit_dx": GATES[
                "maximum_gross_displacement_dx"],
        },
        "patch_gate": {
            "primary_minimum_mean_support_fraction": (
                DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
            "sensitivity_thresholds": list(
                PATCH_SUPPORT_SENSITIVITY_THRESHOLDS),
            "patch_scales_m": list(PATCH_SCALES_M),
        },
        "budgets": {
            "maximum_direct_transport_wall_s": (
                args.maximum_direct_transport_wall_s),
            "maximum_form_factor_replicate_wall_s": (
                args.maximum_form_factor_replicate_wall_s),
            "maximum_endpoint_job_wall_s": args.maximum_endpoint_job_wall_s,
            "maximum_total_wall_s": args.maximum_total_wall_s,
            "maximum_process_count": args.maximum_process_count,
        },
        "prerequisite_order": [
            "sealed base checkpoint and provenance",
            "hash-bound direct transport",
            "level 8 row closure",
            "level 16 row closure and exact 8-to-16 nested extension",
            "level 32 row closure and exact 16-to-32 nested extension",
            "level responses at one frozen horizon",
            "16-to-32 authority plus paired-scramble confidence",
            "level 32 final uncertainty and direct gross-motion gate",
            "bounded hold; no Stage B",
        ],
        "failure_policy": (
            "stop at the first failed construction prerequisite; never expand above "
            "32 rays and never enter Stage B automatically"),
    }


def _supervised_worker(args):
    run(args)


def _supervised_cli(args):
    """Put every native Warp/Numba call behind a parent-enforced process deadline."""
    if getattr(args, "plan_only", False):
        print(json.dumps(_execution_plan_receipt(args), indent=2, sort_keys=True))
        return 0
    process = mp.get_context("spawn").Process(
        target=_supervised_worker, args=(args,), daemon=False)
    process.start()
    process.join(float(args.maximum_total_wall_s))
    if process.is_alive():
        process.terminate()
        process.join(5.0)
        if process.is_alive():
            process.kill()
            process.join(5.0)
        payload = {
            "schema": SCHEMA,
            "status": "bounded_timeout",
            "stage": args.stage,
            "timeout": {
                "reason": "parent-enforced native-process wall deadline",
                "maximum_total_wall_s": float(args.maximum_total_wall_s),
                "physics_conclusion_permitted": False,
            },
        }
        _write_json_atomic(args.output, payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 2
    if process.exitcode != 0:
        raise RuntimeError(
            f"replicated closure worker exited with status {process.exitcode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_supervised_cli(parse_args()))
