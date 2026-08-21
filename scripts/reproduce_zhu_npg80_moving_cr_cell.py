#!/usr/bin/env python3
"""Reproduce one frozen Oxford moving-Cr blind-board trajectory.

The production audit intentionally executes a content-addressed board in
parallel.  This narrow driver makes a failed cell independently reproducible
without changing the job specification, consuming other missing cells, or
writing a board cache.  Its JSON output is diagnostic only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_zhu_npg80_moving_cr_profiles import (
    ANALOG_BOARD,
    FILM_MATERIAL,
    MASK_MATERIAL,
    PREREGISTRATION,
    REACTOR_DOSE,
    SOURCE_Z_UM,
    _load,
    _run_trajectory,
    _scenario_inputs,
    _boundary,
    _router,
)
from petch.feature_geometry_state_3d import FeatureGeometry3D
from petch.feature_step_3d import advance_feature_step_3d
from petch.material_mechanism_3d import MaterialSurfaceState3D


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width-nm", type=float, required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--selectivity", type=float, required=True)
    parser.add_argument("--duration-s", type=float, default=1200.0)
    parser.add_argument("--dx-nm", type=float, default=10.0)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--checkpoint-after-s", type=float, default=0.0)
    parser.add_argument("--resume-checkpoint", type=Path)
    args = parser.parse_args()

    preregistration = _load(PREREGISTRATION)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_DOSE)
    scenarios = _scenario_inputs(preregistration, reactor)
    matching = tuple(
        scenario for scenario in scenarios
        if str(scenario["name"]) == str(args.scenario)
    )
    if len(matching) != 1:
        available = tuple(str(scenario["name"]) for scenario in scenarios)
        raise ValueError(
            f"scenario {args.scenario!r} is not unique; available={available}"
        )
    rates = (
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"]),
    )
    density = float(np.mean(
        preregistration["surface_response_axes"]["ald_tio2_density_kg_m3"]))

    def write_checkpoint(payload):
        if args.checkpoint is None or payload["elapsed_s"] < args.checkpoint_after_s:
            return
        geometry = payload["geometry"]
        state = payload["surface_state"]
        arrays = {
            "phi": geometry.phi,
            "material_id": geometry.material_id,
        }
        layer_ids = []
        if geometry.material_levelsets is not None:
            for index, (material_id, levelset) in enumerate(
                    sorted(geometry.material_levelsets.items())):
                arrays[f"material_levelset_{index}"] = levelset
                layer_ids.append(int(material_id))
        state_names = []
        upper_bounds = {}
        remap_modes = {}
        if state is not None:
            upper_bounds = dict(state.upper_bounds)
            remap_modes = dict(state.remap_modes)
            for index, (name, values) in enumerate(sorted(state.fields.items())):
                arrays[f"surface_state_{index}"] = values
                state_names.append(name)
        metadata = {
            "dx": geometry.dx,
            "mesh_length_unit_m": geometry.mesh_length_unit_m,
            "mesh_origin_m": geometry.mesh_origin_m,
            "material_levelset_ids": layer_ids,
            "surface_state_names": state_names,
            "surface_state_upper_bounds": upper_bounds,
            "surface_state_remap_modes": remap_modes,
            "surface_state_mesh_fingerprint": payload[
                "surface_state_mesh_fingerprint"],
            "elapsed_s": payload["elapsed_s"],
            "accepted_steps": payload["accepted_steps"],
            "step_duration_s": payload["step_duration_s"],
        }
        arrays["metadata_json"] = np.asarray(json.dumps(metadata, sort_keys=True))
        args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.checkpoint, **arrays)

    if args.resume_checkpoint is not None:
        with np.load(args.resume_checkpoint, allow_pickle=False) as stored:
            metadata = json.loads(str(stored["metadata_json"]))
            layers = {
                int(material_id): stored[f"material_levelset_{index}"]
                for index, material_id in enumerate(
                    metadata["material_levelset_ids"])
            }
            geometry = FeatureGeometry3D(
                stored["phi"], stored["material_id"], metadata["dx"],
                metadata["mesh_length_unit_m"],
                mesh_origin_m=tuple(metadata["mesh_origin_m"]),
                material_levelsets=layers)
            state = None
            if metadata["surface_state_names"]:
                fields = {
                    name: stored[f"surface_state_{index}"]
                    for index, name in enumerate(
                        metadata["surface_state_names"])
                }
                state = MaterialSurfaceState3D(
                    fields, metadata["surface_state_upper_bounds"],
                    metadata["surface_state_remap_modes"])
        boundary = _boundary(matching[0], preregistration)
        router = _router(
            scenario_name=matching[0]["name"],
            tio2_rate_nm_min=max(rates),
            selectivity=args.selectivity,
            density_kg_m3=density)
        pitch_nm = float(preregistration["inferred_geometry_board"]["pitch_nm"])
        step = advance_feature_step_3d(
            geometry, boundary,
            {matching[0]["name"]: "energetic_bombardment"}, router,
            etchable_material_ids=(FILM_MATERIAL, MASK_MATERIAL),
            duration_s=float(metadata["step_duration_s"]),
            source_bounds=(0.0, pitch_nm * 1.0e-3, 0.0, pitch_nm * 1.0e-3),
            source_z=SOURCE_Z_UM,
            surface_state=state,
            surface_state_mesh_fingerprint=metadata[
                "surface_state_mesh_fingerprint"],
            ballistic_transport="face_gather",
            ballistic_periodic_lateral=True,
            ballistic_face_quadrature_points=int(
                preregistration["deterministic_feature_transport"][
                    "triangle_quadrature_points"]),
            profile_periodic_lateral=True,
            topology_change_policy="continue_gas_cavity_and_material_extinction",
            surface_state_remap_backend="indexed_knn",
            reinitialization_method="cr2",
            transport_device=args.transport_device)
        remap_materials = step.state_remap_diagnostics["materials"]
        print(json.dumps({
            "status": "checkpoint step passed",
            "elapsed_s_before_step": metadata["elapsed_s"],
            "accepted_steps_before_step": metadata["accepted_steps"],
            "step_duration_s": metadata["step_duration_s"],
            "next_surface_state_mesh_fingerprint": (
                step.next_surface_state_mesh_fingerprint),
            "surface_state_transfer_fingerprint": (
                step.state_remap_diagnostics["transfer_fingerprint"]),
            "maximum_state_remap_relative_conservation_residual": max(
                float(item["max_relative_conservation_residual"])
                for item in remap_materials.values()),
            "reassigned_unresolved_material_nodes": step.diagnostics[
                "reassigned_unresolved_material_nodes"],
        }, indent=2, sort_keys=True))
        return

    profiles = _run_trajectory(
        width_nm=args.width_nm,
        scenario=matching[0],
        rates_nm_min=rates,
        selectivity=args.selectivity,
        duration_s=args.duration_s,
        dx_nm=args.dx_nm,
        preregistration=preregistration,
        transport_device=args.transport_device,
        pre_step_callback=write_checkpoint,
    )
    print(json.dumps(profiles, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
