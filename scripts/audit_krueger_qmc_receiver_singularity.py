#!/usr/bin/env python3
"""Reproduce the frozen Krüger QMC receiver-area singularity diagnosis.

This is a post-mortem diagnostic, not a convergence or prediction run. It
evaluates one zero-motion step on the exact accepted 10.249480 s checkpoint
with the production eight-ray diffuse-neutral operator and with diffuse
neutral transport disabled. The paired control changes no geometry, boundary,
surface state, energetic transport, or chemistry.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.feature_step_3d import (
    _periodic_lateral_surface_images,
    advance_feature_step_3d,
)
from petch.reactor_boundary import build_krueger_2024_development_boundary

import krueger_2024_trench_pilot as pilot


DEFAULT_INPUT = Path("/private/tmp/krueger_guo_transient_dt125_dx10")


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _scalar_at_face(value, face: int) -> float:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return float(array)
    return float(array[int(face)])


def _energetic_fluxes_at_face(surface_fluxes, face: int) -> dict[str, float]:
    return {
        population.name: _scalar_at_face(population.flux_m2_s, face)
        for population in surface_fluxes.energetic_fluxes
    }


def _neutral_fluxes_at_face(surface_fluxes, face: int) -> dict[str, float]:
    return {
        name: _scalar_at_face(value, face)
        for name, value in sorted(surface_fluxes.neutral_flux_m2_s.items())
    }


def _maximum_ledger_residual(exchange) -> float:
    return max(
        (
            float(np.max(np.abs(exchange.residual_units_m2(name))))
            for name in exchange.removed_units_m2
        ),
        default=0.0,
    )


def _build_physics(geometry, configuration: dict):
    domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(domain[2])
    boundary = build_krueger_2024_development_boundary(
        pilot.DATA,
        n_transverse_neutral=5,
        n_normal_neutral=8,
        reference_plane_m=source_z * geometry.mesh_length_unit_m,
        neutral_direction_polar_order=int(
            configuration["neutral_direction_polar_order"]
        ),
        neutral_direction_azimuthal_order=int(
            configuration["neutral_direction_azimuthal_order"]
        ),
        ion_energy_bin_eV=float(configuration["ion_energy_bin_eV"]),
        ion_angle_bin_deg=float(configuration["ion_angle_bin_deg"]),
        ion_azimuthal_closure="axisymmetric_uniform",
        ion_azimuthal_order=int(configuration["ion_azimuthal_order"]),
        ion_flux_normalization=1.0,
    )
    mechanism = build_krueger_2024_material_router_3d(
        surface_model="guo_tml",
        guo_aggregate_ion_formula=None,
        guo_translating_layer_thickness_nm=float(
            configuration["guo_translating_layer_thickness_nm"]
        ),
        effective_mask_crosslinked_growth_fraction=0.0,
        yield_energy_model=str(configuration["yield_energy_model"]),
        deposition_layer_depth_nm=float(
            configuration["deposition_layer_depth_nm"]
        ),
        oxygen_half_saturation_flux_m2_s=None,
    )
    role = {
        species.name: (
            "energetic_bombardment"
            if species.charge_number != 0
            else "neutral_reactant"
        )
        for species in boundary.species
    }
    return boundary, mechanism, role, domain, source_z


def _evaluate(
    geometry,
    state,
    fingerprint,
    *,
    boundary,
    mechanism,
    role,
    domain,
    source_z,
    configuration,
    checkpoint_step,
    radiosity_enabled,
):
    radiosity = None
    if radiosity_enabled:
        radiosity = {
            "rays_per_face": int(configuration["radiosity_rays_per_face"]),
            "seed": int(configuration["seed"]) + 10000,
            "periodic_lateral": True,
            "domain_size": domain,
            "relative_tolerance": float(
                configuration["radiosity_relative_tolerance"]
            ),
            "maximum_iterations": int(
                configuration["radiosity_maximum_iterations"]
            ),
        }
    return advance_feature_step_3d(
        geometry,
        boundary,
        role,
        mechanism,
        etchable_material_ids=(1, 2),
        duration_s=0.0,
        source_bounds=(0.0, domain[0], 0.0, domain[1]),
        source_z=source_z,
        surface_state=state,
        surface_state_mesh_fingerprint=fingerprint,
        n_position=int(configuration["n_position"]),
        seed=int(configuration["seed"]) + int(checkpoint_step),
        cfl_number=0.25,
        reinitialize=True,
        reinitialization_method="cr2",
        profile_periodic_lateral=True,
        transport_device=str(configuration["transport_device"]),
        neutral_radiosity_options=radiosity,
        ballistic_transport=str(configuration["ballistic_transport"]),
        ballistic_face_quadrature_points=int(
            configuration["ballistic_face_quadrature_points"]
        ),
        topology_change_policy=str(configuration["topology_change_policy"]),
        surface_state_remap_backend=str(
            configuration["surface_state_remap_backend"]
        ),
    )


def _narrow_band_assignments(result, geometry, local_face: int) -> dict:
    velocity = np.asarray(result.face_velocity_mesh_units_s, dtype=float)[
        result.active_face_index
    ]
    domain = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    _, centroids = _periodic_lateral_surface_images(
        velocity, result.active_face_centroid, domain
    )
    x, y, z = geometry.coordinate_arrays
    selected = np.abs(geometry.phi) < 4.0 * geometry.dx
    i, j, k = np.where(selected)
    points = np.stack((x[i], y[j], z[k]), axis=1)
    nearest = cKDTree(centroids).query(points)[1] % len(velocity)
    assigned = np.flatnonzero(nearest == int(local_face))
    return {
        "assigned_narrow_band_node_count": int(assigned.size),
        "assigned_narrow_band_nodes_mesh_units": [
            [float(value) for value in points[index]]
            for index in assigned
        ],
    }


def _face_receipt(result, geometry, global_face: int) -> dict:
    active = np.asarray(result.active_face_index, dtype=int)
    positions = np.flatnonzero(active == int(global_face))
    if positions.size != 1:
        raise RuntimeError("diagnostic face is not uniquely active")
    local = int(positions[0])
    area = float(result.active_face_area[local])
    etch = np.asarray(result.surface.etch_velocity_m_s, dtype=float)
    growth = np.asarray(
        result.surface.normal_growth_velocity_m_s, dtype=float
    )
    face_velocity_m_s = (
        float(result.face_velocity_mesh_units_s[global_face])
        * geometry.mesh_length_unit_m
    )
    return {
        "global_face_index": int(global_face),
        "active_local_face_index": local,
        "material_id": int(result.face_material_id[global_face]),
        "centroid_mesh_units": [
            float(value) for value in result.active_face_centroid[local]
        ],
        "area_mesh_units2": area,
        "area_m2": area * geometry.mesh_length_unit_m**2,
        "area_over_active_median": (
            area / float(np.median(result.active_face_area))
        ),
        "signed_face_velocity_m_s": face_velocity_m_s,
        "etch_velocity_m_s": float(etch[local]),
        "growth_velocity_m_s": float(growth[local]),
        "neutral_flux_m2_s": _neutral_fluxes_at_face(
            result.transport.surface_fluxes, global_face
        ),
        "energetic_flux_m2_s": _energetic_fluxes_at_face(
            result.transport.surface_fluxes, global_face
        ),
        **_narrow_band_assignments(result, geometry, local),
    }


def _run_receipt(input_directory: Path) -> dict:
    checkpoint_path = input_directory / "checkpoint.npz"
    audit_path = input_directory / "audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    configuration = audit["configuration"]
    if (
        configuration["surface_model"] != "guo_tml"
        or configuration["radiosity_backend"] != "scrambled_qmc_3d"
        or int(configuration["radiosity_rays_per_face"]) != 8
        or float(configuration["dx_um"]) != 0.01
    ):
        raise ValueError("input is not the sealed eight-ray 10 nm Guo checkpoint")
    geometry, state, fingerprint, metadata = pilot._load_checkpoint(
        checkpoint_path
    )
    if int(metadata["step"]) != 85 or not math.isclose(
        float(metadata["physical_time_s"]),
        10.249480050516949,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("input checkpoint identity/time is not the diagnosed state")
    boundary, mechanism, role, domain, source_z = _build_physics(
        geometry, configuration
    )
    with_radiosity = _evaluate(
        geometry,
        state,
        fingerprint,
        boundary=boundary,
        mechanism=mechanism,
        role=role,
        domain=domain,
        source_z=source_z,
        configuration=configuration,
        checkpoint_step=int(metadata["step"]),
        radiosity_enabled=True,
    )
    without_radiosity = _evaluate(
        geometry,
        state,
        fingerprint,
        boundary=boundary,
        mechanism=mechanism,
        role=role,
        domain=domain,
        source_z=source_z,
        configuration=configuration,
        checkpoint_step=int(metadata["step"]),
        radiosity_enabled=False,
    )
    velocity = np.abs(
        np.asarray(with_radiosity.face_velocity_mesh_units_s, dtype=float)
    )
    global_face = int(np.argmax(velocity))
    with_face = _face_receipt(with_radiosity, geometry, global_face)
    without_face = _face_receipt(without_radiosity, geometry, global_face)
    max_without = float(
        np.max(
            np.abs(
                np.asarray(
                    without_radiosity.face_velocity_mesh_units_s, dtype=float
                )
            )
        )
        * geometry.mesh_length_unit_m
    )
    max_with = float(velocity[global_face] * geometry.mesh_length_unit_m)
    neutral_at_face = with_face["neutral_flux_m2_s"]
    neutral_control = without_face["neutral_flux_m2_s"]
    gates = {
        "same_geometry_face_count": {
            "passed": (
                with_radiosity.diagnostics["face_count"]
                == without_radiosity.diagnostics["face_count"]
            )
        },
        "material_ledgers_exact": {
            "value": max(
                _maximum_ledger_residual(
                    with_radiosity.surface.material_exchange
                ),
                _maximum_ledger_residual(
                    without_radiosity.surface.material_exchange
                ),
            ),
            "limit": 0.0,
        },
        "control_target_neutral_flux_zero": {
            "value": max(neutral_control.values(), default=0.0),
            "limit": 0.0,
        },
        "qmc_target_neutral_flux_positive": {
            "value": max(neutral_at_face.values(), default=0.0),
            "minimum": 0.0,
        },
        "qmc_speed_spike_over_100x_control": {
            "value": max_with / max_without,
            "minimum": 100.0,
        },
        "target_area_below_1e_4_active_median": {
            "value": with_face["area_over_active_median"],
            "limit": 1.0e-4,
        },
    }
    for gate in gates.values():
        if "passed" in gate:
            continue
        if "limit" in gate:
            gate["passed"] = gate["value"] <= gate["limit"]
        else:
            gate["passed"] = gate["value"] > gate["minimum"]
    return {
        "schema": "petch.krueger.qmc-receiver-singularity.v1",
        "claim": (
            "Frozen paired control isolates the moving-profile speed spike to "
            "the eight-ray diffuse-neutral receiver tally on a vanishing-area "
            "mask triangle; it is not evidence for changing chemistry."
        ),
        "input": {
            "directory": str(input_directory),
            "checkpoint_sha256": _sha256(checkpoint_path),
            "audit_sha256": _sha256(audit_path),
            "config_hash": audit["config_hash"],
            "checkpoint_step": int(metadata["step"]),
            "physical_time_s": float(metadata["physical_time_s"]),
        },
        "operator": {
            "radiosity_backend": "scrambled_qmc_3d",
            "rays_per_source_face": 8,
            "receiver_estimator": (
                "categorical source-ray hit weight divided by receiver area"
            ),
            "control": (
                "identical boundary, energetic transport, surface state, "
                "chemistry, and geometry with diffuse neutral exchange disabled"
            ),
        },
        "with_radiosity": {
            "maximum_face_speed_m_s": max_with,
            "target_face": with_face,
        },
        "without_radiosity": {
            "maximum_face_speed_m_s": max_without,
            "same_target_face": without_face,
        },
        "speed_ratio_qmc_over_control": max_with / max_without,
        "gates": gates,
        "all_diagnostic_gates_passed": all(
            gate["passed"] for gate in gates.values()
        ),
        "decision": {
            "qmc_moving_profile_authority": False,
            "chemistry_change_authorized": False,
            "surface_smoothing_authorized": False,
            "krueger_authority_candidate": "deterministic_extruded_2d",
            "general_3d_reentry_requirement": (
                "replicated, physical-patch-local exchange/velocity convergence "
                "or an error-controlled hierarchical deterministic estimator"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-directory", type=Path, default=DEFAULT_INPUT
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = _run_receipt(args.input_directory)
    payload = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
