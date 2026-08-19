#!/usr/bin/env python3
"""Freeze the target-free Oxford TiO2 square-pillar prediction board.

The exact recipe and conserved reactor ensemble are already committed.  This
audit adds deterministic periodic 3-D feature transport on ideal vertical-wall
snapshots, then propagates explicitly labeled TiO2/Cr response sensitivities.
It never reads a target SEM, target depth, or coefficient selected from either.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.boundary_state import PlasmaBoundaryState
from petch.boundary_transport_3d import gather_boundary_state_ballistic_3d
from petch.feature_geometry_backend_3d import UniformFeatureGeometryBackend3D
from petch.feature_step_3d import (
    _surface_gas_normals,
    make_square_pillar_mask_geometry_3d,
)
from petch.iadf_two_component import (
    build_two_component_boundary,
    kim_2025_reference_iadf,
)
from petch.tio2_square_pillar import (
    integrate_square_pillar_depth,
    integrate_square_pillar_depth_from_blanket_rate,
)


DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
PREREGISTRATION = DATA / "square_pillar_blind_preregistration.json"
GEOMETRY_EVIDENCE = DATA / "device_geometry_evidence.json"
REACTOR_DOSE = (
    ROOT / "results" / "curated" / "zhu_npg80_daughter_wafer_dose_v1"
    / "audit.json"
)
ANALOG_BOARD = DATA / "janissen_tio2_analog_board.json"
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_square_pillar_blind_v1"
    / "audit.json"
)

FILM_MATERIAL = 1
MASK_MATERIAL = 2
BASE_MATERIAL = 3
BASE_TOP_UM = 0.1
FILM_THICKNESS_UM = 0.7
MASK_THICKNESS_UM = 0.045
FILM_TOP_UM = BASE_TOP_UM + FILM_THICKNESS_UM
MASK_TOP_UM = FILM_TOP_UM + MASK_THICKNESS_UM
DOMAIN_HEIGHT_UM = 1.05
SOURCE_Z_UM = 1.0


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _render(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _transport_scenarios(preregistration: dict, reactor: dict):
    transport = preregistration["deterministic_feature_transport"]
    powered_drop = [
        float(row["powered_electrode_sheath_drop_V"])
        for row in reactor["power_board"]
    ]
    energy = {
        "low": 0.5 * min(powered_drop),
        "high": max(powered_drop),
    }
    scenarios = []
    species = []
    for energy_label, energy_eV in energy.items():
        for tail_fraction in transport["tail_fraction_sensitivity"]:
            tail_label = str(tail_fraction).replace(".", "p")
            name = f"ion_{energy_label}_tail_{tail_label}"
            boundary = build_two_component_boundary(
                kim_2025_reference_iadf(tail_fraction=float(tail_fraction)),
                1.0,
                energy_eV,
                name=name,
                n_polar=int(transport["polar_quadrature_nodes"]),
                azimuthal_order=int(transport["azimuthal_quadrature_nodes"]),
                reference_plane_m=SOURCE_Z_UM * 1.0e-6,
                extra_provenance={
                    "evidence_class": "cross_machine_transport_sensitivity",
                    "target_iead_measured": False,
                    "impact_energy_label": energy_label,
                },
            )
            species.append(boundary.species[0])
            scenarios.append({
                "name": name,
                "impact_energy_label": energy_label,
                "impact_energy_eV": energy_eV,
                "tail_fraction": float(tail_fraction),
            })
    return scenarios, PlasmaBoundaryState(
        tuple(species),
        reference_plane_m=SOURCE_Z_UM * 1.0e-6,
        provenance={
            "source": "blind Oxford feature-transport sensitivity ensemble",
            "target_iead_measured": False,
        },
    )


def _region_integral(population, selected, areas):
    return float(np.dot(
        np.asarray(population.flux_m2_s)[selected],
        np.asarray(areas)[selected],
    ))


def _snapshot(
    *,
    pitch_nm: float,
    width_nm: float,
    depth_nm: float,
    dx_nm: float,
    boundary: PlasmaBoundaryState,
    face_quadrature_points: int,
):
    pitch_um = pitch_nm * 1.0e-3
    width_um = width_nm * 1.0e-3
    depth_um = depth_nm * 1.0e-3
    dx_um = dx_nm * 1.0e-3
    if depth_nm == 0.0:
        source_area = pitch_um ** 2
        floor_area = source_area - width_um ** 2
        return {
            "width_nm": width_nm,
            "etched_depth_nm": depth_nm,
            "mesh_spacing_nm": dx_nm,
            "surface_triangle_count": 0,
            "flat_limit_evaluated_analytically": True,
            "floor_projected_area_um2": floor_area,
            "mask_projected_area_um2": width_um ** 2,
            "floor_projected_area_relative_residual": 0.0,
            "mask_projected_area_relative_residual": 0.0,
            "maximum_incident_measure_relative_residual": 0.0,
            "scenario": {
                species.name: {
                    "floor_local_transmission": 1.0,
                    "mask_top_local_transmission": 1.0,
                    "floor_fraction_of_source_plane_dose": (
                        floor_area / source_area),
                    "film_sidewall_fraction_of_source_plane_dose": 0.0,
                    "total_incident_measure_relative_residual": 0.0,
                }
                for species in boundary.species
            },
        }
    geometry = make_square_pillar_mask_geometry_3d(
        pitch=pitch_um,
        domain_height=DOMAIN_HEIGHT_UM,
        dx=dx_um,
        pillar_width=width_um,
        film_thickness=FILM_THICKNESS_UM,
        mask_thickness=MASK_THICKNESS_UM,
        base_top=BASE_TOP_UM,
        etched_depth=depth_um,
        mesh_length_unit_m=1.0e-6,
        film_material_id=FILM_MATERIAL,
        mask_material_id=MASK_MATERIAL,
        base_material_id=BASE_MATERIAL,
    )
    surface = UniformFeatureGeometryBackend3D(
        geometry, periodic_axes=(0, 1)
    ).extract_surface()
    normals = _surface_gas_normals(
        surface.vertices_mesh,
        surface.faces,
        surface.centroids_mesh,
        geometry,
    )
    roles = {item.name: "energetic_bombardment" for item in boundary.species}
    result = gather_boundary_state_ballistic_3d(
        boundary,
        roles,
        surface.vertices_mesh,
        surface.faces,
        surface.areas_mesh2,
        surface.centroids_mesh,
        normals,
        source_bounds=(0.0, pitch_um, 0.0, pitch_um),
        source_z=SOURCE_Z_UM,
        mesh_length_unit_m=1.0e-6,
        face_quadrature_points=face_quadrature_points,
        periodic_lateral=True,
        domain_size=(pitch_um, pitch_um, DOMAIN_HEIGHT_UM),
        device="cpu",
    )

    floor_z = FILM_TOP_UM - depth_um
    tolerance = 1.75 * dx_um
    floor = (
        (surface.face_material_id == FILM_MATERIAL)
        & (np.abs(surface.centroids_mesh[:, 2] - floor_z) <= tolerance)
        & (normals[:, 2] > 0.05)
    )
    mask_top = (
        (surface.face_material_id == MASK_MATERIAL)
        & (surface.centroids_mesh[:, 2] > MASK_TOP_UM - dx_um)
        & (surface.centroids_mesh[:, 2] <= MASK_TOP_UM + dx_um)
        & (normals[:, 2] > 0.05)
    )
    film_sidewall = (
        (surface.face_material_id == FILM_MATERIAL)
        & (surface.centroids_mesh[:, 2] > floor_z + tolerance)
        & (surface.centroids_mesh[:, 2] < FILM_TOP_UM - 0.5 * dx_um)
        & (np.abs(normals[:, 2]) < 0.8)
    )
    floor_projected_area = float(np.dot(
        surface.areas_mesh2[floor], normals[floor, 2]))
    mask_projected_area = float(np.dot(
        surface.areas_mesh2[mask_top], normals[mask_top, 2]))
    expected_floor_area = pitch_um ** 2 - width_um ** 2
    expected_mask_area = width_um ** 2
    if floor_projected_area <= 0.0 or mask_projected_area <= 0.0:
        raise RuntimeError("square-pillar transport regions are unresolved")
    floor_area_residual = (
        floor_projected_area / expected_floor_area - 1.0
    )
    mask_area_residual = mask_projected_area / expected_mask_area - 1.0
    if max(abs(floor_area_residual), abs(mask_area_residual)) > 0.05:
        raise RuntimeError(
            "square-pillar projected area failed the 5% mesh gate: "
            f"width={width_nm:g} nm depth={depth_nm:g} nm dx={dx_nm:g} nm "
            f"floor_residual={floor_area_residual:.6g} "
            f"mask_residual={mask_area_residual:.6g}")

    source_area = pitch_um ** 2
    populations = {
        population.name: population
        for population in result.surface_fluxes.energetic_fluxes
    }
    by_scenario = {}
    maximum_conservation_residual = 0.0
    for species in boundary.species:
        population = populations[species.name]
        floor_dose = _region_integral(
            population, floor, surface.areas_mesh2)
        mask_dose = _region_integral(
            population, mask_top, surface.areas_mesh2)
        sidewall_dose = _region_integral(
            population, film_sidewall, surface.areas_mesh2)
        total_dose = float(np.dot(
            population.flux_m2_s, surface.areas_mesh2))
        expected_total = species.flux_m2_s * source_area
        conservation = total_dose / expected_total - 1.0
        maximum_conservation_residual = max(
            maximum_conservation_residual, abs(conservation))
        by_scenario[species.name] = {
            "floor_local_transmission": (
                floor_dose / (species.flux_m2_s * floor_projected_area)
            ),
            "mask_top_local_transmission": (
                mask_dose / (species.flux_m2_s * mask_projected_area)
            ),
            "floor_fraction_of_source_plane_dose": (
                floor_dose / expected_total
            ),
            "film_sidewall_fraction_of_source_plane_dose": (
                sidewall_dose / expected_total
            ),
            "total_incident_measure_relative_residual": conservation,
        }
    if maximum_conservation_residual > 2.0e-12:
        raise RuntimeError("periodic feature transport failed conservation")
    return {
        "width_nm": width_nm,
        "etched_depth_nm": depth_nm,
        "mesh_spacing_nm": dx_nm,
        "surface_triangle_count": int(len(surface.faces)),
        "flat_limit_evaluated_analytically": False,
        "floor_projected_area_um2": floor_projected_area,
        "mask_projected_area_um2": mask_projected_area,
        "floor_projected_area_relative_residual": floor_area_residual,
        "mask_projected_area_relative_residual": mask_area_residual,
        "maximum_incident_measure_relative_residual": (
            maximum_conservation_residual
        ),
        "scenario": by_scenario,
    }


def _transport_board(preregistration, boundary):
    geometry = preregistration["inferred_geometry_board"]
    options = preregistration["deterministic_feature_transport"]
    rows = []
    for width_nm in geometry["width_nm"]:
        for depth_nm in options["etched_depth_nodes_nm"]:
            rows.append(_snapshot(
                pitch_nm=float(geometry["pitch_nm"]),
                width_nm=float(width_nm),
                depth_nm=float(depth_nm),
                dx_nm=float(options["mesh_spacing_nm"]),
                boundary=boundary,
                face_quadrature_points=int(
                    options["triangle_quadrature_points"]),
            ))
    return rows


def _convergence_board(preregistration, boundary, coarse_rows):
    geometry = preregistration["inferred_geometry_board"]
    options = preregistration["deterministic_feature_transport"]
    sentinels = (
        (float(min(geometry["width_nm"])), 560.0),
        (float(max(geometry["width_nm"])), 560.0),
    )
    coarse = {
        (row["width_nm"], row["etched_depth_nm"]): row
        for row in coarse_rows
    }
    rows = []
    maximum = 0.0
    for width_nm, depth_nm in sentinels:
        fine = _snapshot(
            pitch_nm=float(geometry["pitch_nm"]),
            width_nm=width_nm,
            depth_nm=depth_nm,
            dx_nm=float(options["convergence_sentinel_spacing_nm"]),
            boundary=boundary,
            face_quadrature_points=int(options["triangle_quadrature_points"]),
        )
        comparison = []
        for name, fine_value in fine["scenario"].items():
            coarse_value = coarse[(width_nm, depth_nm)]["scenario"][name]
            relative = (
                coarse_value["floor_local_transmission"]
                / fine_value["floor_local_transmission"] - 1.0
            )
            maximum = max(maximum, abs(relative))
            comparison.append({
                "scenario": name,
                "coarse_floor_local_transmission": coarse_value[
                    "floor_local_transmission"],
                "fine_floor_local_transmission": fine_value[
                    "floor_local_transmission"],
                "coarse_to_fine_relative_change": relative,
            })
        rows.append({
            "width_nm": width_nm,
            "etched_depth_nm": depth_nm,
            "coarse_spacing_nm": options["mesh_spacing_nm"],
            "fine_spacing_nm": options["convergence_sentinel_spacing_nm"],
            "scenario": comparison,
        })
    return {
        "sentinels": rows,
        "maximum_floor_transmission_relative_change": maximum,
        "passed_5_percent": maximum < 0.05,
    }


def _trajectory_board(preregistration, reactor, scenarios, transport_rows):
    geometry = preregistration["inferred_geometry_board"]
    surface = preregistration["surface_response_axes"]
    depth_nodes = preregistration[
        "deterministic_feature_transport"]["etched_depth_nodes_nm"]
    lookup = {
        (row["width_nm"], row["etched_depth_nm"]): row
        for row in transport_rows
    }
    rows = []
    for width_nm in geometry["width_nm"]:
        for scenario in scenarios:
            floor = [
                lookup[(width_nm, depth)]["scenario"][scenario["name"]]
                ["floor_local_transmission"]
                for depth in depth_nodes
            ]
            mask = [
                lookup[(width_nm, depth)]["scenario"][scenario["name"]]
                ["mask_top_local_transmission"]
                for depth in depth_nodes
            ]
            for power in reactor["power_board"]:
                for density in surface["ald_tio2_density_kg_m3"]:
                    for removal_yield in surface[
                        "effective_tio2_formula_units_per_incident_positive_ion"
                    ]:
                        for selectivity in surface["tio2_to_cr_selectivity"]:
                            result = integrate_square_pillar_depth(
                                depth_nodes_nm=depth_nodes,
                                floor_transmission=floor,
                                mask_transmission=mask,
                                film_thickness_nm=700.0,
                                mask_thickness_nm=45.0,
                                positive_ion_flux_m2_s=power[
                                    "central_3mm_positive_ion_flux_m2_s"],
                                duration_s=1200.0,
                                mass_density_kg_m3=density,
                                formula_units_per_incident_ion=removal_yield,
                                tio2_to_cr_selectivity=selectivity,
                                integration_step_s=1.0,
                            )
                            rows.append({
                                "width_nm": width_nm,
                                "transport_scenario": scenario["name"],
                                "absorbed_power_sensitivity_W": power[
                                    "absorbed_power_sensitivity_W"],
                                "central_positive_ion_flux_m2_s": power[
                                    "central_3mm_positive_ion_flux_m2_s"],
                                "mass_density_kg_m3": density,
                                "effective_tio2_formula_units_per_ion": (
                                    removal_yield),
                                "tio2_to_cr_selectivity": selectivity,
                                **asdict(result),
                            })
    return rows


def _analog_trajectory_board(
    preregistration, analog, scenarios, transport_rows
):
    """Apply published adjacent rates without relabeling them as target yields."""
    geometry = preregistration["inferred_geometry_board"]
    depth_nodes = preregistration[
        "deterministic_feature_transport"]["etched_depth_nodes_nm"]
    lookup = {
        (row["width_nm"], row["etched_depth_nm"]): row
        for row in transport_rows
    }
    rates_nm_min = sorted({
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"]),
        float(analog["closest_stack_witness"]["implied_feature_rate_nm_min"]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"]),
    })
    selectivities = preregistration[
        "surface_response_axes"]["tio2_to_cr_selectivity"]
    rows = []
    for width_nm in geometry["width_nm"]:
        for scenario in scenarios:
            floor = [
                lookup[(width_nm, depth)]["scenario"][scenario["name"]]
                ["floor_local_transmission"]
                for depth in depth_nodes
            ]
            mask = [
                lookup[(width_nm, depth)]["scenario"][scenario["name"]]
                ["mask_top_local_transmission"]
                for depth in depth_nodes
            ]
            for rate_nm_min in rates_nm_min:
                for selectivity in selectivities:
                    result = integrate_square_pillar_depth_from_blanket_rate(
                        depth_nodes_nm=depth_nodes,
                        floor_transmission=floor,
                        mask_transmission=mask,
                        film_thickness_nm=700.0,
                        mask_thickness_nm=45.0,
                        blanket_tio2_rate_nm_s=rate_nm_min / 60.0,
                        duration_s=1200.0,
                        tio2_to_cr_selectivity=selectivity,
                        integration_step_s=1.0,
                    )
                    rows.append({
                        "width_nm": width_nm,
                        "transport_scenario": scenario["name"],
                        "cross_machine_tio2_rate_nm_min": rate_nm_min,
                        "tio2_to_cr_selectivity": selectivity,
                        "same_machine": False,
                        "same_chemistry": False,
                        "same_tio2_material_state": False,
                        "transferred_as_target_coefficient": False,
                        **asdict(result),
                    })
    return {
        "source": analog["source"],
        "evidence_class": "published_cross_machine_process_analog",
        "rates_nm_min": rates_nm_min,
        "selectivities": list(selectivities),
        "warning": (
            "This slice answers a conditional question: what would the blind "
            "Oxford transport operator predict if the adjacent Janissen TiO2 "
            "rate/selectivity transferred? It is not Freddie's surface law."
        ),
        "trajectory_board": rows,
        "width_summary": _summaries(geometry["width_nm"], rows),
    }


def _summaries(widths, trajectories):
    summaries = []
    for width in widths:
        rows = [row for row in trajectories if row["width_nm"] == width]
        exhaustion = [
            row["mask_exhaustion_time_s"] for row in rows
            if row["mask_exhaustion_time_s"] is not None
        ]
        summaries.append({
            "width_nm": width,
            "sensitivity_grid_size": len(rows),
            "mask_pinned_depth_envelope_nm": [
                min(row["mask_pinned_depth_nm"] for row in rows),
                max(row["mask_pinned_depth_nm"] for row in rows),
            ],
            "controlled_depth_envelope_nm": [
                min(row["controlled_depth_nm"] for row in rows),
                max(row["controlled_depth_nm"] for row in rows),
            ],
            "mask_survival_grid_fraction_not_probability": (
                sum(row["mask_survives_duration"] for row in rows) / len(rows)
            ),
            "mask_pinned_full_clear_grid_fraction_not_probability": (
                sum(row["mask_pinned_depth_nm"] >= 699.999 for row in rows)
                / len(rows)
            ),
            "mask_exhaustion_time_envelope_s": (
                None if not exhaustion else [min(exhaustion), max(exhaustion)]
            ),
        })
    return summaries


def build_receipt():
    preregistration = _load(PREREGISTRATION)
    geometry = _load(GEOMETRY_EVIDENCE)
    reactor = _load(REACTOR_DOSE)
    analog = _load(ANALOG_BOARD)
    if (
        preregistration["target_sem_used"]
        or preregistration["target_depth_used"]
        or not preregistration["frozen_before_specific_condition_sem"]
        or geometry["sem_target_used"]
        or geometry["measured_profile_target_used"]
        or reactor["target_outcome_used"]
        or reactor["measured_depth_target_used"]
        or analog["sem_target_used"]
        or analog["measured_depth_target_used"]
    ):
        raise ValueError("held-out target entered the blind square-pillar board")
    condition = preregistration["condition_id"]
    if any(item["condition_id"] != condition for item in (
        geometry, reactor, analog
    )):
        raise ValueError("square-pillar inputs belong to different conditions")

    scenarios, boundary = _transport_scenarios(preregistration, reactor)
    transport = _transport_board(preregistration, boundary)
    convergence = _convergence_board(preregistration, boundary, transport)
    trajectories = _trajectory_board(
        preregistration, reactor, scenarios, transport)
    analog_trajectories = _analog_trajectory_board(
        preregistration, analog, scenarios, transport)
    summaries = _summaries(
        preregistration["inferred_geometry_board"]["width_nm"],
        trajectories,
    )
    maximum_conservation = max(
        row["maximum_incident_measure_relative_residual"]
        for row in transport
    )
    return {
        "schema": "petch.zhu-npg80-square-pillar-blind-board.v1",
        "condition_id": condition,
        "target_sem_used": False,
        "target_depth_used": False,
        "coefficient_selected_from_target": None,
        "inputs": {
            "preregistration": {
                "path": str(PREREGISTRATION.relative_to(ROOT)),
                "sha256": _hash(PREREGISTRATION),
            },
            "geometry_evidence": {
                "path": str(GEOMETRY_EVIDENCE.relative_to(ROOT)),
                "sha256": _hash(GEOMETRY_EVIDENCE),
            },
            "reactor_wafer_dose": {
                "path": str(REACTOR_DOSE.relative_to(ROOT)),
                "sha256": _hash(REACTOR_DOSE),
            },
            "tio2_cr_analog_board": {
                "path": str(ANALOG_BOARD.relative_to(ROOT)),
                "sha256": _hash(ANALOG_BOARD),
            },
        },
        "transport_scenarios": scenarios,
        "transport_snapshot_board": transport,
        "transport_convergence": convergence,
        "trajectory_sensitivity_board": trajectories,
        "width_summary": summaries,
        "published_cross_machine_analog_slice": analog_trajectories,
        "certification": {
            "exact_recipe_and_stack_checksum_bound": True,
            "square_geometry_is_inferred_not_target_measured": True,
            "periodic_3d_transport_is_deterministic": True,
            "incident_measure_conserved_below_2e_12": (
                maximum_conservation < 2.0e-12),
            "maximum_incident_measure_relative_residual": (
                maximum_conservation),
            "transport_convergence_passed_5_percent": convergence[
                "passed_5_percent"],
            "surface_response_is_labeled_sensitivity_not_target_law": True,
            "target_tio2_surface_law_validated": False,
            "target_cr_surface_law_validated": False,
            "target_iead_measured": False,
            "supports_blind_physics_bounded_depth_board": True,
            "supports_unique_target_depth": False,
            "supports_unique_target_sem": False,
            "sidewall_profile_claimed": False,
            "target_outcome_used": False,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = _render(build_receipt())
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed blind square-pillar board is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
