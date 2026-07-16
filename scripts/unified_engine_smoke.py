#!/usr/bin/env python3
"""Execute the assembled public common-engine paths on bounded manufactured cases.

This is an operational smoke, not an experimental validation.  It exercises, in one command:

* independently routed mask/substrate mechanisms and conservative same-material redeposition;
* public quasi-static charged profile evolution;
* safe checkpoint serialization and continuation; and
* public finite-arrival ensemble evolution with distinct seeds.

Run from the repository root with ``python scripts/unified_engine_smoke.py``.
"""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import audit_charging_coevolution_c3 as c3
from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.charging_checkpoint_3d import PhysicalChargingCheckpoint3D
from petch.feature_step_3d import make_rectangular_trench_geometry_3d
from petch.material_mechanism_3d import MaterialMechanismRouter3D, MaterialSurfaceState3D
from petch.physical_api import (
    PhysicalChargingEnsembleProcess, PhysicalChargingProcess, PhysicalProcess,
)
from petch.physical_sputtering import PhysicalSputterMechanism, PhysicalSputterParameters
from petch.surface_kinetics import EnergeticYield, ParameterEvidence
from petch.surface_product_redeposition_3d import (
    SurfaceProductRedepositionContract3D, SurfaceProductRedepositionLaw3D,
)


def _sputter(material, inventory, product, reference_yield):
    evidence = {
        name: ParameterEvidence(
            "manufactured unified-engine smoke", "analytic",
            supports_prediction_within_declared_domain=True)
        for name in (
            "bulk_material_unit_density_m3", "sputter_yield",
            "emitted_product_mass_amu", "emission_angular_model",
            "emission_energy_model")}
    return PhysicalSputterMechanism(PhysicalSputterParameters(
        material_name=material, material_inventory_name=inventory,
        projectile_species=("Ar+",), bulk_material_unit_density_m3=1e28,
        sputter_yield=EnergeticYield(reference_yield, 20.0, 100.0),
        emitted_product_name=product, emitted_product_mass_amu=28.0,
        emitted_material_units_per_particle=1.0,
        emission_angular_model="diffuse_cosine", emission_energy_model="thompson",
        emission_energy_parameters={
            "surface_binding_energy_eV": 4.0, "maximum_energy_eV": 100.0},
        evidence=evidence))


def _redeposition_law(name, material_id):
    return SurfaceProductRedepositionLaw3D(
        name, material_id,
        {1: float(material_id == 1), 2: float(material_id == 2)}, 1e28,
        parameter_sources={
            "sticking_probability_by_material": "manufactured same-material smoke",
            "bulk_material_unit_density_m3": "manufactured same-material smoke"},
        parameter_bounds={
            "sticking_probability_by_material": (0.0, 1.0),
            "bulk_material_unit_density_m3": (0.9e28, 1.1e28)})


def _run_material_path():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.2, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.3,
        substrate_top=1.0, etched_depth=0.2)
    before = {
        material_id: np.asarray(field).copy()
        for material_id, field in geometry.material_levelsets.items()}
    boundary = PlasmaBoundaryState((SpeciesBoundaryState(
        "Ar+", 1, 39.948, 1e21, [[0.0, 0.0, 10.0]], [1.0],
        provenance={"source": "manufactured unified-engine smoke"}),),
        reference_plane_m=1.8e-6,
        provenance={"source": "manufactured unified-engine smoke"})
    router = MaterialMechanismRouter3D(
        {1: _sputter("substrate", "substrate_units", "substrate_product", 0.2),
         2: _sputter("mask", "mask_units", "mask_product", 0.05)},
        provenance={
            1: {"source": "manufactured substrate smoke", "bounds": "test-only"},
            2: {"source": "manufactured mask smoke", "bounds": "test-only"}})
    process = PhysicalProcess(
        geometry, boundary, {"Ar+": "energetic_bombardment"}, router, (1, 2),
        1.0, 1, (-0.1, 1.1, -0.1, 0.3), 1.8,
        solver_options=dict(
            ballistic_transport="face_gather", ballistic_face_quadrature_points=3,
            surface_product_redeposition_options={
                "contract": SurfaceProductRedepositionContract3D((
                    _redeposition_law("material_1:substrate_product", 1),
                    _redeposition_law("material_2:mask_product", 2))),
                "rays_per_face": 8, "seed": 11},
            cfl_number=0.3, reinitialize=False, transport_device="cpu"))
    result = process.run()
    step = result.steps[0]
    layer_changed = {
        str(material_id): bool(not np.array_equal(
            result.geometry.material_levelsets[material_id], before[material_id]))
        for material_id in before}
    if (not isinstance(result.surface_state, MaterialSurfaceState3D)
            or not all(layer_changed.values())
            or step.surface_product_redeposition is None
            or step.diagnostics["product_redeposition_relative_balance_error"] > 1e-10):
        raise RuntimeError("assembled multi-material/redeposition smoke failed")
    return dict(
        layer_changed=layer_changed,
        product_population_count=step.diagnostics["product_population_count"],
        redeposition_balance_error=step.diagnostics[
            "product_redeposition_relative_balance_error"],
        manifest_schema=result.run_manifest["schema"],
        mechanism_provenance=result.run_manifest["surface_mechanism"])


def _charged_process(*, finite_arrivals=False):
    charging_options = c3.options()
    solver_options = dict(
        n_position=32 if finite_arrivals else 128, seed=81,
        trajectory_fixed_dt=0.005, trajectory_max_steps=1000,
        reinitialize=False, transport_device="cpu")
    duration = 1e-9 if finite_arrivals else 0.01
    if finite_arrivals:
        charging_options["physical_arrival_statistics"] = "poisson"
        solver_options["bias_mode"] = "physical_time_resolved"
    return PhysicalChargingProcess(
        geometry=c3.geometry(), boundary=c3.boundary(),
        species_role={"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism=c3.mechanism(), charging_system_builder=c3.poisson_system,
        etchable_material_ids=(1,), duration_s=duration, n_steps=1,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=0.25,
        charging_options=charging_options, solver_options=solver_options)


def _run_charging_restart_path():
    process = _charged_process()
    first = process.run()
    with TemporaryDirectory(prefix="petch-unified-smoke-") as directory:
        path = Path(directory) / "checkpoint.npz"
        checkpoint = PhysicalChargingCheckpoint3D.from_result(first.solve)
        checkpoint.save(path)
        restored = PhysicalChargingCheckpoint3D.load(path)
        resumed = process.continue_from_checkpoint(
            restored, duration_s=0.01, n_steps=1, continuation_seed_stride=101).run()
    if (not first.steps[0].charging.converged
            or not resumed.steps[0].charging.converged
            or resumed.run_manifest["restart_source_manifest_sha256"]
            != checkpoint.source_manifest_sha256):
        raise RuntimeError("assembled charging checkpoint/restart smoke failed")
    return dict(
        initial_converged=first.steps[0].charging.converged,
        resumed_converged=resumed.steps[0].charging.converged,
        manifest_schema=first.run_manifest["schema"],
        restart_bound_to_source_manifest=True,
        exact_transport_reused=(
            first.steps[0].feature.transport is first.steps[0].charging.final_step.transport),
        charge_remap_balance_error=first.steps[0].charge_remap.relative_charge_balance_error,
        recovery_and_error_budget=first.run_manifest["recovery_and_error_budget"])


def _run_finite_arrival_path():
    process = _charged_process(finite_arrivals=True)
    ensemble = PhysicalChargingEnsembleProcess(
        process, realization_count=2, seed_stride=101).run()
    if ensemble.seeds != (81, 182) or ensemble.realization_count != 2:
        raise RuntimeError("assembled finite-arrival ensemble smoke failed")
    return dict(
        seeds=list(ensemble.seeds), realization_count=ensemble.realization_count,
        geometry_shape=list(ensemble.mean_levelset.shape),
        nonnegative_geometry_variance=bool(
            np.all(ensemble.standard_deviation_levelset >= 0.0)),
        statistical_claim_ready=ensemble.statistical_claim_ready)


def main():
    report = dict(
        schema="petch-unified-engine-operational-smoke-v1",
        status="pass",
        material_path=_run_material_path(),
        charging_restart_path=_run_charging_restart_path(),
        finite_arrival_path=_run_finite_arrival_path(),
        claim="operational integration only; no experimental validation claim")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
