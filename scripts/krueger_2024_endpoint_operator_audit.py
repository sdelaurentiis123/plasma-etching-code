#!/usr/bin/env python3
"""Certify the accelerated Krüger march against a refined endpoint operator.

The moving development replay uses a provenance-priced boundary quadrature and one
ballistic face point.  This audit freezes the final geometry and material state, then
evaluates independent zero-duration operators on that identical state:

* the historical tensor-velocity quadratures,
* analytically speed-marginalized neutral angular quadratures under refinement, and
* exact versus compressed digitized ion IEAD quadratures.

Zero duration makes all three evaluations exact geometry no-ops.  The comparison is
The 12x24 angular rule, exact digitized IEAD, and three face points form the reference.  The
comparison is therefore about incident flux and the instantaneous profile-velocity operator only;
it cannot be contaminated by remapping or a different evolved geometry.
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
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from krueger_2024_trench_pilot import _load_checkpoint, _maximum_ledger_residual
from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.feature_step_3d import advance_feature_step_3d
from petch.reactor_boundary import build_krueger_2024_development_boundary


DATA = ROOT / "data" / "experimental" / "krueger_2024"
GATES = {
    # Predeclared before inspecting the completed 60 s endpoint.  Net error is normalized by
    # the reference gross rate, so an accidental signed cancellation cannot make the gate blow up.
    "normalized_net_rate_error": 0.01,
    "area_weighted_l1_relative_error": 0.02,
    "area_weighted_rms_relative_error": 0.02,
    "maximum_relative_error": 0.10,
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


def _sha256(path):
    digest = sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_operator_error(candidate, reference, area):
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    area = np.asarray(area, dtype=float)
    if candidate.shape != reference.shape or candidate.shape != area.shape:
        raise ValueError("operator arrays must share one face mesh")
    difference = candidate - reference
    gross = float(np.sum(area * np.abs(reference)))
    squared = float(np.sum(area * reference ** 2))
    maximum = float(np.max(np.abs(reference))) if reference.size else 0.0
    tiny = np.finfo(float).tiny
    return {
        "candidate_net_rate": float(np.sum(area * candidate)),
        "reference_net_rate": float(np.sum(area * reference)),
        "reference_gross_rate": gross,
        "normalized_net_rate_error": float(
            abs(np.sum(area * difference)) / max(gross, tiny)),
        "area_weighted_l1_relative_error": float(
            np.sum(area * np.abs(difference)) / max(gross, tiny)),
        "area_weighted_rms_relative_error": float(
            np.sqrt(np.sum(area * difference ** 2) / max(squared, tiny))),
        "maximum_relative_error": float(
            np.max(np.abs(difference)) / max(maximum, tiny))
            if difference.size else 0.0,
    }


def _surface_flux_by_species(result):
    output = {
        name: np.asarray(value, dtype=float)
        for name, value in result.transport.surface_fluxes.neutral_flux_m2_s.items()
    }
    for population in result.transport.surface_fluxes.energetic_fluxes:
        if population.name in output:
            raise ValueError(f"duplicate surface-flux species {population.name!r}")
        output[population.name] = np.asarray(population.flux_m2_s, dtype=float)
    return output


def _evaluate(
        geometry, state, fingerprint, *, boundary_mode, ion_bins, face_points,
        pilot_config,
        radiosity_rays, seed, ballistic_transport=None, n_position=None,
        transport_device="cpu"):
    realized_domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(realized_domain[2])
    boundary_options = {}
    if boundary_mode == "legacy_compressed_tensor":
        boundary_options.update(
            n_transverse_neutral=int(pilot_config["neutral_transverse_order"]),
            n_normal_neutral=int(pilot_config["neutral_normal_order"]),
        )
    elif boundary_mode == "full_tensor":
        boundary_options.update(n_transverse_neutral=5, n_normal_neutral=8)
    elif boundary_mode.startswith("angular_"):
        polar, azimuthal = (
            int(value) for value in boundary_mode.removeprefix("angular_").split("x"))
        boundary_options.update(
            neutral_direction_polar_order=polar,
            neutral_direction_azimuthal_order=azimuthal,
        )
    else:
        raise ValueError(f"unknown endpoint boundary mode {boundary_mode!r}")
    if ion_bins is not None:
        boundary_options.update(
            ion_energy_bin_eV=float(ion_bins[0]),
            ion_angle_bin_deg=float(ion_bins[1]),
        )
    if pilot_config.get("ion_azimuthal_closure") is not None:
        boundary_options.update(
            ion_azimuthal_closure=str(
                pilot_config["ion_azimuthal_closure"]),
            ion_azimuthal_order=int(pilot_config["ion_azimuthal_order"]),
        )
    boundary = build_krueger_2024_development_boundary(
        DATA, reference_plane_m=source_z * geometry.mesh_length_unit_m,
        **boundary_options)
    mechanism = build_krueger_2024_material_router_3d(
        effective_mask_crosslinked_growth_fraction=float(
            pilot_config.get("effective_mask_crosslinked_growth_fraction", 0.0)),
        oxide_etch_yield_scale=float(
            pilot_config.get("oxide_etch_yield_scale", 1.0)))
    role = {
        species.name: (
            "energetic_bombardment"
            if species.charge_number != 0 else "neutral_reactant")
        for species in boundary.species
    }
    radiosity = None
    if bool(pilot_config["radiosity_enabled"]):
        radiosity = {
            "rays_per_face": int(radiosity_rays),
            "seed": int(seed) + 10000,
            "periodic_lateral": True,
            "domain_size": realized_domain,
            "relative_tolerance": float(pilot_config["radiosity_relative_tolerance"]),
            "maximum_iterations": int(pilot_config["radiosity_maximum_iterations"]),
        }
    started = perf_counter()
    result = advance_feature_step_3d(
        geometry, boundary, role, mechanism,
        etchable_material_ids=(1, 2), duration_s=0.0,
        source_bounds=(0.0, realized_domain[0], 0.0, realized_domain[1]),
        source_z=source_z, surface_state=state,
        surface_state_mesh_fingerprint=fingerprint,
        n_position=(
            int(pilot_config["n_position"])
            if n_position is None else int(n_position)), seed=int(seed),
        cfl_number=0.25, reinitialize=True, reinitialization_method="cr2",
        profile_periodic_lateral=True, transport_device=str(transport_device),
        neutral_radiosity_options=radiosity,
        ballistic_transport=(
            str(pilot_config["ballistic_transport"])
            if ballistic_transport is None else str(ballistic_transport)),
        ballistic_face_quadrature_points=int(face_points),
    )
    if not np.array_equal(result.geometry.phi, geometry.phi):
        raise RuntimeError("zero-duration endpoint audit changed the authoritative geometry")
    return result, boundary, perf_counter() - started


def run(args):
    source = Path(args.source)
    checkpoint = source / "checkpoint.npz"
    pilot_audit_path = source / "audit.json"
    pilot_audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if pilot_audit.get("status") != "complete":
        raise RuntimeError("endpoint certification requires a completed pilot")
    geometry, state, fingerprint, checkpoint_metadata = _load_checkpoint(checkpoint)
    pilot_config = pilot_audit["configuration"]

    # The 12x24 direct angular rule is the refinement reference for the speed-marginalized
    # neutral operator.  The old full tensor rule remains as an independent formulation check;
    # neither is silently assumed exact.  Full/compressed-ion twins separate IEAD reduction error.
    variants = (
        ("legacy_compressed_q3", "legacy_compressed_tensor", (500.0, 0.5), 3),
        ("tensor_full_q3", "full_tensor", None, 3),
        ("tensor_neutral_ion_compressed_q3", "full_tensor", (500.0, 0.5), 3),
        ("angular_4x8_full_ion_q3", "angular_4x8", None, 3),
        ("angular_6x12_full_ion_q3", "angular_6x12", None, 3),
        ("angular_8x16_full_ion_q3", "angular_8x16", None, 3),
        ("angular_12x24_full_ion_q3", "angular_12x24", None, 3),
        ("angular_6x12_compressed_ion_q3", "angular_6x12", (500.0, 0.5), 3),
        ("angular_8x16_compressed_ion_q3", "angular_8x16", (500.0, 0.5), 3),
        ("angular_8x16_ion_250x025_q3", "angular_8x16", (250.0, 0.25), 3),
        ("angular_8x16_ion_200x020_q3", "angular_8x16", (200.0, 0.2), 3),
    )
    evaluated = {}
    summary = {}
    for name, boundary_mode, ion_bins, face_points in variants:
        result, boundary, wall = _evaluate(
            geometry, state, fingerprint, boundary_mode=boundary_mode,
            ion_bins=ion_bins,
            face_points=face_points, pilot_config=pilot_config,
            radiosity_rays=int(args.radiosity_rays), seed=int(args.seed),
            ballistic_transport="face_gather", n_position=None)
        evaluated[name] = result
        summary[name] = {
            "boundary_mode": boundary_mode,
            "ion_quadrature_bins": None if ion_bins is None else {
                "energy_eV": float(ion_bins[0]),
                "angle_deg": float(ion_bins[1]),
            },
            "neutral_quadrature": dict(boundary.provenance["neutral_quadrature"]),
            "ballistic_transport": "face_gather",
            "ballistic_face_quadrature_points": face_points,
            "seed": int(args.seed),
            "boundary_quadrature_node_count": int(
                boundary.provenance["total_boundary_quadrature_nodes"]),
            "wall_time_s": float(wall),
            "maximum_material_ledger_residual_units_m2": float(
                _maximum_ledger_residual(result.surface.material_exchange)),
            "maximum_radiosity_balance_error": max(
                (float(item["relative_balance_error"])
                 for item in result.diagnostics["neutral_radiosity"].values()),
                default=0.0),
        }

    reference = evaluated["angular_12x24_full_ion_q3"]
    reference_index = np.asarray(reference.active_face_index, dtype=int)
    reference_centroid = np.asarray(reference.active_face_centroid, dtype=float)
    reference_area = np.asarray(reference.active_face_area, dtype=float)
    comparisons = {}
    pairs = (
        ("legacy_to_angular_reference", "legacy_compressed_q3",
         "angular_12x24_full_ion_q3"),
        ("tensor_to_angular_reference", "tensor_full_q3",
         "angular_12x24_full_ion_q3"),
        ("ion_compression_on_tensor_neutrals", "tensor_neutral_ion_compressed_q3",
         "tensor_full_q3"),
        ("angular_4x8_refinement", "angular_4x8_full_ion_q3",
         "angular_12x24_full_ion_q3"),
        ("angular_6x12_refinement", "angular_6x12_full_ion_q3",
         "angular_12x24_full_ion_q3"),
        ("angular_8x16_refinement", "angular_8x16_full_ion_q3",
         "angular_12x24_full_ion_q3"),
        ("ion_compression_on_angular_6x12", "angular_6x12_compressed_ion_q3",
         "angular_6x12_full_ion_q3"),
        ("ion_compression_on_angular_8x16", "angular_8x16_compressed_ion_q3",
         "angular_8x16_full_ion_q3"),
        ("ion_250x025_on_angular_8x16", "angular_8x16_ion_250x025_q3",
         "angular_8x16_full_ion_q3"),
        ("ion_200x020_on_angular_8x16", "angular_8x16_ion_200x020_q3",
         "angular_8x16_full_ion_q3"),
        ("candidate_to_reference", "angular_8x16_ion_250x025_q3",
         "angular_12x24_full_ion_q3"),
    )
    for label, candidate_name, reference_name in pairs:
        candidate = evaluated[candidate_name]
        pair_reference = evaluated[reference_name]
        if (not np.array_equal(candidate.active_face_index, reference_index)
                or not np.array_equal(pair_reference.active_face_index, reference_index)
                or not np.allclose(
                    candidate.active_face_centroid, reference_centroid,
                    rtol=0.0, atol=0.0)
                or not np.allclose(
                    pair_reference.active_face_centroid, reference_centroid,
                    rtol=0.0, atol=0.0)):
            raise RuntimeError("endpoint variants do not share an identical active-face mesh")
        candidate_velocity = np.asarray(
            candidate.face_velocity_mesh_units_s, dtype=float)[reference_index]
        pair_reference_velocity = np.asarray(
            pair_reference.face_velocity_mesh_units_s, dtype=float)[reference_index]
        velocity = _relative_operator_error(
            candidate_velocity, pair_reference_velocity, reference_area)
        candidate_flux = _surface_flux_by_species(candidate)
        reference_flux = _surface_flux_by_species(pair_reference)
        if set(candidate_flux) != set(reference_flux):
            raise RuntimeError("endpoint variants expose different flux species")
        flux = {
            name: _relative_operator_error(
                candidate_flux[name][reference_index],
                reference_flux[name][reference_index], reference_area)
            for name in sorted(candidate_flux)
        }
        comparisons[label] = {
            "candidate": candidate_name,
            "reference": reference_name,
            "profile_velocity": velocity,
            "incident_flux_by_species": flux,
        }

    march = comparisons["candidate_to_reference"]["profile_velocity"]
    gate_results = {
        name: bool(march[name] <= threshold)
        for name, threshold in GATES.items()
    }
    payload = {
        "status": "pass" if all(gate_results.values()) else "fail",
        "scientific_scope": (
            "instantaneous exact-geometry endpoint certification; this does not replace "
            "the moving-profile timestep/grid/sample refinement gates"),
        "source": {
            "pilot_audit": pilot_audit_path.name,
            "pilot_audit_sha256": _sha256(pilot_audit_path),
            "checkpoint": checkpoint.name,
            "checkpoint_sha256": _sha256(checkpoint),
            "checkpoint_metadata": checkpoint_metadata,
        },
        "predeclared_profile_velocity_gates": GATES,
        "profile_velocity_gate_results": gate_results,
        "variants": summary,
        "comparisons": comparisons,
    }
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "march_to_reference": march,
        "output": str(destination),
    }, indent=2, sort_keys=True))
    return payload


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source", default=ROOT / "results" / "krueger_2024_base_60s_adaptive")
    parser.add_argument(
        "--output",
        default=(ROOT / "results" / "krueger_2024_base_60s_adaptive"
                 / "endpoint_operator_audit.json"))
    parser.add_argument("--radiosity-rays", type=int, default=32)
    parser.add_argument("--seed", type=int, default=241)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
