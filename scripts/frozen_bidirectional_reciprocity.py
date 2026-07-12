"""Reproduce the high-sample forward/adjoint audit on a frozen nodal charging state.

The artifact supplies only geometry and voltage. Source physics is the documented AR4 campaign
boundary: finite-transit Ar ions plus half-Maxwellian electrons. Proposal mixtures are numerical
importance distributions and are scored exactly against that unchanged physical boundary density.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from petch.boundary_state import (
    MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    folded_normal_tangential_proposal,
    maxwellian_electron_boundary_state,
    mixture_boundary_proposal,
    qmc_boundary_proposal,
    qmc_boundary_proposal_with_auxiliary,
)
from petch.boundary_transport import (
    adjoint_boundary_state_face_flux,
    bidirectional_boundary_state_cell_flux,
)
from petch.sheath import CollisionlessRFSheath


def _campaign_boundary_and_proposals(use_grazing):
    sheath = CollisionlessRFSheath(
        40.0, 10.0, 2e6, 4.0, 40.0, thickness_m=5e-4)
    ion = collisionless_sheath_boundary_state(
        sheath, 1e19, n_phase=16, tangential_temperature_eV=0.2,
        n_transverse=3, normal_energy_bins=64).get("ion")
    electron = maxwellian_electron_boundary_state(
        4.0, 1e19, n_transverse=3, n_normal=6).get("electron")
    boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)
    tails = {
        "ion": SpeciesBoundaryState(
            "ion_tail", 1, 40.0, 1.0, [[0.0, 0.0, np.sqrt(40.0)]], [1.0],
            density_model=MaxwellianFluxVelocityDensity(40.0)),
        "electron": SpeciesBoundaryState(
            "electron_tail", -1, electron.mass_amu, 1.0,
            [[0.0, 0.0, np.sqrt(16.0)]], [1.0],
            density_model=MaxwellianFluxVelocityDensity(16.0)),
    }
    physical = {"ion": ion, "electron": electron}
    proposals = {}
    for index, name in enumerate(("ion", "electron")):
        source = physical[name]; tail = tails[name]
        if use_grazing:
            components = (
                source,
                folded_normal_tangential_proposal(source, +1),
                folded_normal_tangential_proposal(source, -1),
                tail,
            )
            proposals[name] = mixture_boundary_proposal(
                components, (0.4, 0.25, 0.25, 0.1), name=f"{name}_grazing_proposal")
        else:
            proposals[name] = mixture_boundary_proposal((
                qmc_boundary_proposal(source, 8, seed=101 + 100 * index),
                qmc_boundary_proposal(tail, 8, seed=102 + 100 * index),
            ), (0.8, 0.2), name=f"{name}_proposal")
    return boundary, proposals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact")
    parser.add_argument("--species", choices=("ion", "electron"), default="ion")
    parser.add_argument("--fixed-dt", type=float, required=True)
    parser.add_argument("--face-offset", type=float, default=1e-3)
    parser.add_argument("--source-offset", type=float, default=1e-3)
    parser.add_argument("--adjoint-base", type=int, default=12)
    parser.add_argument("--adjoint-max", type=int, default=18)
    parser.add_argument("--forward-base", type=int, default=14)
    parser.add_argument("--forward-max", type=int, default=19)
    parser.add_argument("--element-absolute", type=float, default=0.01)
    parser.add_argument("--element-relative", type=float, default=0.15)
    parser.add_argument("--cell-i", type=int)
    parser.add_argument("--cell-j", type=int)
    parser.add_argument("--raw-adjoint", action="store_true")
    parser.add_argument("--grazing", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with np.load(args.artifact) as saved:
        solid = saved["solid"]
        potential = (saved["failed_potential"]
                     if "failed_potential" in saved else saved["potential"])
        cells = saved["cells"]; normals = saved["normals"]
    if args.cell_i is not None or args.cell_j is not None:
        if args.cell_i is None or args.cell_j is None:
            raise ValueError("both cell coordinates are required")
        selected = np.all(cells == np.array([args.cell_i, args.cell_j]), axis=1)
        if not np.any(selected):
            raise ValueError("requested cell has no exposed face")
        cells = cells[selected]; normals = normals[selected]

    boundary, proposals = _campaign_boundary_and_proposals(args.grazing)
    proposal = proposals[args.species]
    if args.raw_adjoint:
        for replicate in range(4):
            sampled, auxiliary = qmc_boundary_proposal_with_auxiliary(
                proposal, args.adjoint_base, 1, seed=replicate,
                name=f"{args.species}-adaptive-proposal")
            raw = adjoint_boundary_state_face_flux(
                boundary, args.species, potential, solid, cells, normals,
                proposal_species=sampled, face_position_samples=auxiliary[:, 0],
                fixed_dt=args.fixed_dt, face_offset=args.face_offset)
            print(
                "replicate", replicate, "flux", raw["per_face"].tolist(),
                "ess", raw["effective_sample_size"].tolist(),
                "max_fraction", raw["max_sample_fraction"].tolist(),
                "dominant_surface_velocity", raw["dominant_surface_velocity"].tolist(),
                "dominant_exit_velocity", raw["dominant_exit_velocity"].tolist())
        return

    common = dict(
        n_replicates=4, absolute_tolerance=0.01, relative_tolerance=0.05,
        element_absolute_tolerance=args.element_absolute,
        element_relative_tolerance=args.element_relative, fixed_dt=args.fixed_dt)
    start = time.perf_counter()
    result = bidirectional_boundary_state_cell_flux(
        boundary, args.species, potential, solid, cells, normals,
        proposal_species=proposal,
        adjoint_options=dict(
            common, base_log2=args.adjoint_base, max_log2=args.adjoint_max,
            face_offset=args.face_offset),
        forward_options=dict(
            common, base_log2=args.forward_base, max_log2=args.forward_max,
            source_offset=args.source_offset),
        element_absolute_tolerance=args.element_absolute,
        element_relative_tolerance=args.element_relative, switch_factor=2.0)
    elapsed = time.perf_counter() - start
    worst = int(np.argmax(result["estimator_discrepancy_sigma"]))
    print("seconds", elapsed, "fixed_dt", args.fixed_dt)
    print("worst", tuple(result["unique_cells"][worst]))
    print("discrepancy_sigma", result["estimator_discrepancy_sigma"][worst])
    print("forward", result["forward_cell_mean"][worst],
          result["forward_cell_stderr"][worst])
    print("adjoint", result["adjoint_cell_mean"][worst],
          result["adjoint_cell_stderr"][worst])
    print("all_consistent", bool(np.all(result["estimator_consistent"])))
    print("all_converged", bool(np.all(result["cell_converged"])))
    np.savez(
        args.output, fixed_dt=args.fixed_dt, elapsed_s=elapsed,
        unique_cells=result["unique_cells"],
        forward_cell_mean=result["forward_cell_mean"],
        forward_cell_stderr=result["forward_cell_stderr"],
        adjoint_cell_mean=result["adjoint_cell_mean"],
        adjoint_cell_stderr=result["adjoint_cell_stderr"],
        discrepancy_sigma=result["estimator_discrepancy_sigma"],
        estimator_consistent=result["estimator_consistent"],
        cell_converged=result["cell_converged"], method=result["method"])


if __name__ == "__main__":
    main()
