#!/usr/bin/env python3
"""Evolve Mahorowala oxide-mask features from the reactor/wafer receipt.

The operator is deterministic throughout: fixed RF-phase quadrature, fixed
Gaussian ion-angle quadrature, fixed cosine neutral quadrature, adjoint
face-gather visibility, Chang Eq. 3.13 reaction/specular reflection, measured
species-resolved surface laws, and a level-set moving boundary.  Neither the
observed depth nor an optimizer enters the simulation.

This is still evidence-gated.  The target thesis gives source estimates, not
per-run IEAD/IAD measurements; Cl2+ angular response is transferred from the
measured Cl+ curve; and the RF waveform is an equipment-class sensitivity.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

import numpy as np

from petch.boundary_state import (
    PlasmaBoundaryState,
    SpeciesBoundaryState,
)
from petch.chlorine_species_resolved_si import (
    SpeciesResolvedChlorineSiMechanism,
)
from petch.feature_step_3d import (
    make_rectangular_trench_geometry_3d,
    solve_feature_3d,
)
from petch.reactor_global import DiagnosticConditionedRFSheathTransfer


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "mahorowala_1998_diagnostic_conditioned_depth_projection.json"
)
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "mahorowala_1998_deterministic_feature_depth.json"
)

PITCH_UM = 0.560
OPENING_UM = 0.310
MASK_THICKNESS_UM = 0.200
SUBSTRATE_TOP_UM = 0.650
DOMAIN_HEIGHT_UM = 1.000
CELL_LENGTH_UM = 0.040
SOURCE_Z_UM = 0.980
PLASMA_POTENTIAL_EV = 20.0
ELECTRODE_AREA_M2 = 0.04
ION_FWHM_DEG = 10.0
ETCH_DURATION_S = 75.0


def _gaussian_ion_species(
        name, mass_amu, flux_m2_s, energy_eV, energy_weight, *,
        fwhm_deg=ION_FWHM_DEG, angle_order=9):
    node, weight = np.polynomial.hermite.hermgauss(int(angle_order))
    sigma = np.deg2rad(float(fwhm_deg)) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    angle = np.sqrt(2.0) * sigma * node
    angle_weight = weight / np.sqrt(np.pi)
    energy = np.asarray(energy_eV, dtype=float)
    energy_weight = np.asarray(energy_weight, dtype=float)
    ie, ia = np.meshgrid(
        np.arange(energy.size), np.arange(angle.size), indexing="ij")
    selected_energy = energy[ie.ravel()]
    selected_angle = angle[ia.ravel()]
    speed = np.sqrt(selected_energy)
    velocity = np.column_stack((
        speed * np.sin(selected_angle),
        np.zeros(speed.size),
        speed * np.cos(selected_angle),
    ))
    joint_weight = (
        energy_weight[ie.ravel()] * angle_weight[ia.ravel()]
    )
    return SpeciesBoundaryState(
        name,
        1,
        float(mass_amu),
        float(flux_m2_s),
        velocity,
        joint_weight,
        provenance={
            "energy": "power-closed deterministic RF sheath quadrature",
            "angle": (
                "Mahorowala source estimate: Gaussian IAD, 10 degree FWHM; "
                "fixed Gauss-Hermite quadrature"
            ),
        },
    )


def _cosine_neutral_species(
        name, mass_amu, flux_m2_s, *, energy_eV=0.05,
        polar_order=6, azimuth_order=12):
    node, weight = np.polynomial.legendre.leggauss(int(polar_order))
    mu = 0.5 * (node + 1.0)
    # A plane-crossing isotropic gas has p(mu)=2mu.  The [0,1]
    # transformation contributes 1/2, hence quadrature weight mu*w.
    mu_weight = mu * weight
    azimuth = 2.0 * np.pi * (
        np.arange(int(azimuth_order), dtype=float) + 0.5
    ) / int(azimuth_order)
    azimuth_weight = np.full(azimuth.size, 1.0 / azimuth.size)
    im, ip = np.meshgrid(
        np.arange(mu.size), np.arange(azimuth.size), indexing="ij")
    selected_mu = mu[im.ravel()]
    selected_phi = azimuth[ip.ravel()]
    transverse = np.sqrt(np.maximum(1.0 - selected_mu ** 2, 0.0))
    speed = np.sqrt(float(energy_eV))
    velocity = speed * np.column_stack((
        transverse * np.cos(selected_phi),
        transverse * np.sin(selected_phi),
        selected_mu,
    ))
    joint_weight = mu_weight[im.ravel()] * azimuth_weight[ip.ravel()]
    return SpeciesBoundaryState(
        name,
        0,
        float(mass_amu),
        float(flux_m2_s),
        velocity,
        joint_weight,
        provenance={
            "angle": (
                "analytic plane-crossing isotropic flux; fixed "
                "Gauss-Legendre/azimuth quadrature"
            )
        },
    )


def _boundary(row, *, product_wall):
    sheath = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=ELECTRODE_AREA_M2,
        plasma_potential_eV=PLASMA_POTENTIAL_EV,
        frequency_hz=float(row["rf_sheath_frequency_hz"]),
        collapse_fraction=1.0,
        phase_count=96,
        steps_per_period=128,
        steps_per_transit=128,
        source=(
            "replay of diagnostic-conditioned RF sheath from reactor-depth "
            "receipt"
        ),
    ).predict(
        positive_ion_flux_m2_s={
            "Cl+": float(row["wafer_clplus_flux_m2_s"]),
            "Cl2+": float(row["wafer_cl2plus_flux_m2_s"]),
        },
        electron_temperature_eV=(
            (2.0 / 3.0) * float(row["reactor_mean_electron_energy_eV"])
        ),
        electron_density_m3=float(row["reactor_electron_density_m3"]),
        delivered_bias_power_W=float(row["rf_bias_power_W"]),
    )
    species = [
        _gaussian_ion_species(
            "Cl+",
            35.45,
            row["wafer_clplus_flux_m2_s"],
            sheath.distributions["Cl+"].energy_eV,
            sheath.distributions["Cl+"].weight,
        ),
        _gaussian_ion_species(
            "Cl2+",
            70.90,
            row["wafer_cl2plus_flux_m2_s"],
            sheath.distributions["Cl2+"].energy_eV,
            sheath.distributions["Cl2+"].weight,
        ),
        _cosine_neutral_species(
            "Cl", 35.45, row["wafer_atomic_chlorine_flux_m2_s"]
        ),
    ]
    if product_wall != "none":
        product_key = (
            f"table4_product_{product_wall.removeprefix('table4-')}_wall_"
            "sicl2_flux_m2_s"
            if product_wall.startswith("table4-")
            else f"sicl2_{product_wall}_wall_flux_m2_s"
        )
        product_flux = float(row[product_key])
        species.append(_cosine_neutral_species(
            "SiCl2", 98.991, product_flux
        ))
    return PlasmaBoundaryState(
        tuple(species),
        reference_plane_m=SOURCE_Z_UM * 1.0e-6,
        provenance={
            "reactor_receipt_run": int(row["run"]),
            "product_wall_limit": product_wall,
            "observed_depth_used": False,
        },
    )


def _linear_zero(x0, y0, x1, y1):
    return float(x0 - y0 * (x1 - x0) / (y1 - y0))


def _center_depth_nm(geometry):
    x, _, z = geometry.coordinate_arrays
    center = int(np.argmin(np.abs(x - 0.5 * PITCH_UM)))
    rows = range(1, geometry.phi.shape[1] - 1)
    if geometry.phi.shape[1] <= 2:
        rows = range(geometry.phi.shape[1])
    floors = []
    for j in rows:
        values = geometry.phi[center, j]
        crossing = np.flatnonzero(
            (values[:-1] >= 0.0) & (values[1:] < 0.0)
        )
        for k in crossing:
            root = _linear_zero(z[k], values[k], z[k + 1], values[k + 1])
            if root <= SUBSTRATE_TOP_UM + 0.5 * geometry.dx:
                floors.append(root)
    if not floors:
        raise RuntimeError("final feature has no resolved center floor")
    return float((SUBSTRATE_TOP_UM - np.median(floors)) * 1.0e3)


def simulate(row, *, dx_um, product_wall, reflection):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=PITCH_UM,
        cell_length=max(CELL_LENGTH_UM, 2.0 * dx_um),
        domain_height=DOMAIN_HEIGHT_UM,
        dx=dx_um,
        opening_width=OPENING_UM,
        mask_thickness=MASK_THICKNESS_UM,
        substrate_top=SUBSTRATE_TOP_UM,
        etched_depth=0.0,
        mesh_length_unit_m=1.0e-6,
        substrate_material_id=1,
        mask_material_id=2,
    )
    boundary = _boundary(row, product_wall=product_wall)
    role = {
        species.name: (
            "energetic_bombardment"
            if species.charge_number else "neutral_reactant"
        )
        for species in boundary.species
    }
    mechanism = SpeciesResolvedChlorineSiMechanism(
        strict_by_default=False
    )
    reflection_options = (
        {"model": "chang_sawin_chlorine", "max_bounces": 16,
         "minimum_weight": 1.0e-6}
        if reflection else None
    )
    domain_size = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    solved = solve_feature_3d(
        geometry,
        boundary,
        role,
        mechanism,
        etchable_material_ids=(1,),
        duration_s=ETCH_DURATION_S,
        n_steps=30,
        source_bounds=(0.0, PITCH_UM, 0.0, domain_size[1]),
        source_z=SOURCE_Z_UM,
        n_position=1,
        seed=0,
        cfl_number=0.25,
        reinitialize=True,
        reinitialization_method="cr2",
        profile_periodic_lateral=True,
        ballistic_periodic_lateral=True,
        transport_device="cpu",
        ballistic_transport="face_gather",
        ballistic_face_quadrature_points=3,
        grazing_ion_reflection=reflection_options,
        surface_state_remap_backend="common_refinement",
        topology_change_policy="refuse",
        adaptive_timestep_options={
            "initial_step_duration_s": 1.0,
            "minimum_step_duration_s": 0.01,
            "maximum_step_duration_s": 4.0,
            "target_displacement_cells": 0.30,
            "maximum_displacement_cells": 0.65,
            "maximum_retries_per_step": 20,
            "maximum_accepted_steps": 1000,
        },
    )
    depth = _center_depth_nm(solved.geometry)
    observed = row["observed_feature_depth_nm"]
    return {
        "run": int(row["run"]),
        "observed_feature_depth_nm": observed,
        "predicted_center_feature_depth_nm": depth,
        "signed_error_percent": (
            None if observed is None else 100.0 * (depth / observed - 1.0)
        ),
        "surface_plane_depth_nm": float(row[
            "joined_species_resolved_surface_plane_depth_nm_75s"
            if product_wall == "none"
            else (
                f"table4_product_{product_wall.removeprefix('table4-')}_wall_"
                "surface_plane_depth_nm_75s"
                if product_wall.startswith("table4-")
                else f"sicl2_{product_wall}_wall_surface_plane_depth_nm_75s"
            )
        ]),
        "accepted_feature_steps": len(solved.steps),
        "maximum_step_displacement_cells": float(max(
            step.diagnostics["max_displacement_mesh_units"] / dx_um
            for step in solved.steps
        )),
        "within_declared_scope": solved.validity.within_declared_scope,
        "known_limitations": list(solved.validity.known_limitations),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--runs", default="1")
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument(
        "--product-wall",
        choices=(
            "none", "reflective", "reactive",
            "table4-reflective", "table4-reactive",
        ),
        default="table4-reflective",
    )
    parser.add_argument("--no-reflection", action="store_true")
    arguments = parser.parse_args()
    receipt = json.loads(arguments.input.read_text(encoding="utf-8"))
    selected = {
        int(value) for value in str(arguments.runs).split(",") if value
    }
    rows = [row for row in receipt["rows"] if int(row["run"]) in selected]
    if not rows:
        raise ValueError("no selected runs exist in the reactor receipt")
    jobs = int(arguments.jobs)
    if jobs <= 0:
        raise ValueError("jobs must be positive")
    call_arguments = [
        (
            row,
            float(arguments.dx_um),
            str(arguments.product_wall),
            not arguments.no_reflection,
        )
        for row in rows
    ]
    if jobs == 1:
        results = [
            simulate(
                row,
                dx_um=dx_um,
                product_wall=wall,
                reflection=reflection,
            )
            for row, dx_um, wall, reflection in call_arguments
        ]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = [
                executor.submit(
                    simulate,
                    row,
                    dx_um=dx_um,
                    product_wall=wall,
                    reflection=reflection,
                )
                for row, dx_um, wall, reflection in call_arguments
            ]
            results = [future.result() for future in futures]
    usable = [
        row for row in results if row["signed_error_percent"] is not None
    ]
    payload = {
        "schema": "petch.mahorowala_1998_deterministic_feature_depth.v1",
        "reactor_receipt": str(arguments.input),
        "dx_um": float(arguments.dx_um),
        "product_wall_limit": str(arguments.product_wall),
        "chlorine_specular_reflection": not arguments.no_reflection,
        "observed_depth_used_for_conditioning": False,
        "geometry": {
            "pitch_um": PITCH_UM,
            "opening_um": OPENING_UM,
            "oxide_mask_thickness_um": MASK_THICKNESS_UM,
            "etch_duration_s": ETCH_DURATION_S,
        },
        "rows": results,
        "mape_percent": (
            None if not usable else float(np.mean(np.abs([
                row["signed_error_percent"] for row in usable
            ])))
        ),
        "formal_feature_depth_pass": False,
        "evidence_blockers": [
            "target-tool per-run species-resolved IEAD/IAD measurements",
            "measured Cl2+ angular response",
            "target-tool absorbed bias-power fraction and waveform",
            "self-consistent SiClx collision-power/base-chlorine feedback",
            "target-tool conditioned chamber-wall Si/Cl coverage",
            "measured SiClx+ surface deposition probabilities",
        ],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
