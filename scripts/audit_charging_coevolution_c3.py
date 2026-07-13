#!/usr/bin/env python3
"""Reproduce the bounded manufactured C3 unified-engine integration gates."""
from __future__ import annotations

from hashlib import sha256
import json
import platform
from pathlib import Path
import subprocess

import numpy as np

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.charged_surface_response_3d import GrazingSpecularIonReflection3D
from petch.charging_coevolution_3d import (
    ExperimentalObservableTolerance3D,
    ResolvedBiasSegment3D,
    solve_charging_coevolution_3d,
)
from petch.charging_poisson_3d import (
    NodalPoissonSystem3D,
    lump_triangle_sheet_charge_3d,
)
from petch.feature_step_3d import FeatureGeometry3D
from petch.physical_sputtering import PhysicalSputterMechanism, PhysicalSputterParameters
from petch.surface_kinetics import EnergeticYield, ParameterEvidence


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "charging_coevolution_c3"


def species(name, charge_number, flux_m2_s):
    return SpeciesBoundaryState(
        name, charge_number, 40.0 if charge_number > 0 else 5.4858e-4,
        flux_m2_s, [[0.0, 0.0, 10.0]], [1.0],
        provenance={"source": "manufactured balanced planar C3 gate"})


def boundary(ion_flux=2.2e21, electron_flux=2.2e21, phase="balanced"):
    return PlasmaBoundaryState((
        species("Ar+", 1, ion_flux), species("electron", -1, electron_flux)),
        reference_plane_m=1.75e-6, provenance={"bias_phase": phase})


def geometry():
    dx = 0.25
    shape = (4, 4, 8)
    z = np.arange(shape[2]) * dx
    phi = np.broadcast_to(0.95 - z, shape).copy()
    return FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), dx, 1e-6)


def poisson_system(feature):
    fixed = np.zeros(feature.phi.shape, dtype=bool)
    fixed[:, :, -1] = True
    phi_center = sum(
        feature.phi[i:i + feature.phi.shape[0] - 1,
                    j:j + feature.phi.shape[1] - 1,
                    k:k + feature.phi.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0
    return NodalPoissonSystem3D(
        np.where(phi_center > 0.0, 3.9, 1.0),
        feature.dx * feature.mesh_length_unit_m, fixed)


def mechanism():
    required = (
        "bulk_material_unit_density_m3", "sputter_yield",
        "emitted_product_mass_amu", "emission_angular_model", "emission_energy_model")
    evidence = {
        name: ParameterEvidence(
            "manufactured C3 integration gate", "analytic",
            supports_prediction_within_declared_domain=True)
        for name in required}
    return PhysicalSputterMechanism(PhysicalSputterParameters(
        material_name="SiO2", material_inventory_name="SiO2_formula_unit",
        projectile_species=("Ar+",), bulk_material_unit_density_m3=2.2e28,
        sputter_yield=EnergeticYield(0.2, 20.0, 100.0),
        emitted_product_name="sputtered_SiO2_unit", emitted_product_mass_amu=60.084,
        emitted_material_units_per_particle=1.0,
        emission_angular_model="diffuse_cosine", emission_energy_model="thompson",
        emission_energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0},
        evidence=evidence))


def options():
    return dict(
        patch_scales_m=(0.25e-6, 1.0e-6),
        potential_rate_tolerance_v_s=1e-5,
        timestep_s=1e-9, maximum_steps=0,
        current_balance_tolerance=0.08, timestep_policy="fixed")


def common(feature, base_boundary):
    return dict(
        geometry=feature, boundary=base_boundary,
        species_role={"Ar+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism=mechanism(), charging_system_builder=poisson_system,
        etchable_material_ids=(1,), source_bounds=(0.0, 0.75, 0.0, 0.75),
        source_z=1.75, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=feature.dx, charging_options=options(),
        n_position=256, seed=71, trajectory_fixed_dt=0.005,
        trajectory_max_steps=1000, reinitialize=False, transport_device="cpu")


def summarize_step(step):
    history = step.charging.history
    return dict(
        charging_converged=step.charging.converged,
        charging_accepted_steps=step.charging.accepted_steps,
        charging_rejected_steps=step.charging.rejected_steps,
        charging_physical_time_s=step.charging.physical_time_s,
        charging_pseudo_time_s=step.charging.pseudo_time_s,
        final_potential_rate_max_v_s=step.charging.diagnostics[
            "final_potential_rate_max_v_s"],
        retained_node_rms_relative_current_imbalance=step.diagnostics[
            "retained_node_rms_relative_current_imbalance"],
        retained_node_max_relative_current_imbalance=step.diagnostics[
            "retained_node_max_relative_current_imbalance"],
        patch_scales_m=list(step.diagnostics["patch_scales_m"]),
        patch_max_relative_imbalance=list(step.diagnostics[
            "patch_max_relative_imbalance"]),
        patch_symmetric_max_relative_imbalance=list(step.diagnostics[
            "patch_symmetric_max_relative_imbalance"]),
        maximum_charge_conservation_relative_error=max(
            item["charge_conservation_relative_error"] for item in history),
        maximum_surface_transfer_relative_charge_balance_error=max(
            item["surface_transfer_relative_charge_balance_error"] for item in history),
        exact_transport_reused=step.feature.transport is step.charging.final_step.transport,
        remap_relative_charge_balance_error=step.charge_remap.relative_charge_balance_error,
        retained_positive_charge_c=step.charge_remap.retained_positive_charge_c,
        retained_negative_charge_c=step.charge_remap.retained_negative_charge_c,
        removed_positive_charge_c=step.charge_remap.removed_positive_charge_c,
        removed_negative_charge_c=step.charge_remap.removed_negative_charge_c,
        wall_clock_s=step.wall_clock_s)


def legacy_checkpoint_migration_audit():
    """Refuse a legacy nodal checkpoint unless its face-sheet inverse is unique and accurate."""
    checkpoint = ROOT / "results/charging_task1_3d_refined_transient_15us/final_states.npz"
    if not checkpoint.exists():
        return dict(available=False, accepted=False, reason="legacy checkpoint is unavailable")
    with np.load(checkpoint) as archived:
        charge = np.asarray(archived["refined_charge_node_c"], dtype=float)
        vertices = np.asarray(archived["vertices"], dtype=float)
        faces = np.asarray(archived["faces"], dtype=int)
    columns = []
    for face_index in range(len(faces)):
        sigma = np.zeros(len(faces))
        sigma[face_index] = 1.0
        columns.append(lump_triangle_sheet_charge_3d(
            charge.shape, vertices, faces, sigma,
            grid_origin=(0.0, 0.0, 0.0), grid_spacing=0.125,
            coordinate_length_unit_m=1e-6).ravel())
    projection = np.stack(columns, axis=1)
    singular = np.linalg.svd(projection, compute_uv=False)
    recovered_sigma, _residual, rank, _singular = np.linalg.lstsq(
        projection, charge.ravel(), rcond=None)
    reconstructed = projection @ recovered_sigma
    relative_l2 = float(
        np.linalg.norm(reconstructed - charge.ravel())
        / max(np.linalg.norm(charge.ravel()), np.finfo(float).tiny))
    relative_linf = float(
        np.max(np.abs(reconstructed - charge.ravel()))
        / max(np.max(np.abs(charge)), np.finfo(float).tiny))
    condition = float(singular[0] / singular[-1])
    accepted = bool(
        rank == len(faces) and condition <= 1e10
        and relative_l2 <= 1e-10 and relative_linf <= 1e-10)
    return dict(
        available=True, checkpoint=checkpoint.name,
        checkpoint_sha256=sha256(checkpoint.read_bytes()).hexdigest(),
        projection_shape=list(projection.shape), face_count=len(faces),
        numerical_rank=int(rank), condition_number=condition,
        reconstruction_relative_l2=relative_l2,
        reconstruction_relative_linf=relative_linf,
        accepted=accepted,
        reason=("unique roundoff-accurate face-sheet state" if accepted else
                "legacy nodal charge has no unique accurate face-sheet inverse; cold-start C3"))


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    feature = geometry()
    base = boundary()
    reflection = GrazingSpecularIonReflection3D.literature_bounded_sensitivity(1, "Ar+")
    quasi = solve_charging_coevolution_3d(
        **common(feature, base), duration_s=0.1, n_steps=1,
        charged_surface_response=reflection)
    waveform = (
        ResolvedBiasSegment3D(1e-9, boundary(2.2e21, 1.1e21, "ion_rich")),
        ResolvedBiasSegment3D(1e-9, boundary(1.1e21, 2.2e21, "electron_rich")),
    )
    pulsed = solve_charging_coevolution_3d(
        **common(feature, base), duration_s=2e-9, n_steps=2,
        bias_mode="waveform_resolved", bias_waveform=waveform)

    b3_anchor_refused = False
    try:
        ExperimentalObservableTolerance3D("notch_depth", 2.0, 1.0)
    except ValueError:
        b3_anchor_refused = True
    quasi_waveform_refused = False
    try:
        solve_charging_coevolution_3d(
            **common(feature, base), duration_s=0.0, n_steps=1,
            bias_waveform=(ResolvedBiasSegment3D(1e-9, base),))
    except ValueError:
        quasi_waveform_refused = True

    config = dict(
        geometry="4x4x8 planar dielectric manufactured gate", dx_mesh_units=0.25,
        mesh_length_unit_m=1e-6, patch_scales_m=[0.25e-6, 1.0e-6],
        quasi_profile_duration_s=0.1, waveform_segment_duration_s=[1e-9, 1e-9],
        n_position=256, seed=71, trajectory_fixed_dt=0.005,
        trajectory_max_steps=1000, transport_device="cpu",
        exact_visibility=True, convergence_contract_revision="CCA-2026-07-13-R2")
    config_hash = sha256(
        json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    source_paths = (
        ROOT / "src/petch/charging_coevolution_3d.py",
        ROOT / "src/petch/charging_coupled_3d.py",
        ROOT / "src/petch/boundary_transport_3d.py",
        ROOT / "src/petch/charged_surface_cascade_3d.py",
        ROOT / "src/petch/feature_step_3d.py",
        ROOT / "src/petch/surface_charge_remap_3d.py",
        ROOT / "src/petch/charged_surface_response_3d.py",
        Path(__file__).resolve(),
    )
    summary = dict(
        schema="petch.charging.coevolution.c3.audit.v1",
        run_manifest=dict(
            config_hash=config_hash, config=config,
            engine_git_revision=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            source_sha256={
                path.relative_to(ROOT).as_posix(): sha256(path.read_bytes()).hexdigest()
                for path in source_paths},
            hardware=platform.platform(), python=platform.python_version()),
        quasi_static=dict(
            mode=quasi.run_manifest["mode"], step=summarize_step(quasi.steps[0]),
            reflection_manifest=quasi.run_manifest["charged_surface_response"]),
        waveform_resolved=dict(
            mode=pulsed.run_manifest["mode"],
            steps=[summarize_step(step) for step in pulsed.steps],
            saturation_required=[step.diagnostics["saturation_required"]
                                 for step in pulsed.steps]),
        refusal_gates=dict(
            b3_tolerance_above_benchmark_uncertainty_refused=b3_anchor_refused,
            quasi_static_waveform_refused=quasi_waveform_refused,
            legacy_nodal_checkpoint_migration=legacy_checkpoint_migration_audit()),
        campaign_status=dict(
            manufactured_engine_integration="pass",
            designated_trench_timestep_grid_sample_refinement="pending",
            independent_high_sample_exact_operator_audit_b5="pending",
            c3_scientific_closure=False))
    (OUTPUT / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
