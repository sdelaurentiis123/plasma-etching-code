#!/usr/bin/env python3
"""Static aspect-ratio-dependent-etching (ARDE) flux gate through the COMMON feature-3d engine.

This measures the first-principles physics behind radical-transport-limited ARDE (de Boer / Blauw
SF6/O2 Si etching) WITHOUT any evolution, chemistry parameter fit, or aspect-ratio-shaped formula.

Mechanism under test (no shortcut): a neutral radical arrives at the mask plane with the analytic
half-Maxwellian *flux* angular law (kinetic cosine measure, not a fitted cone). It undergoes
collisionless ballistic first flight (`gather_boundary_state_ballistic_3d`) and diffuse molecular-flow
re-emission (`solve_diffuse_neutral_radiosity_3d`, H = D + B(1-s)H with area reciprocity
B[i,j]=A[j]F[j->i]/A[i]). On every wall/floor hit it reacts/recombines with ONE physical sticking
coefficient ``s`` -- the single declared surface property, with provenance, not a fudge. The floor
incident flux, normalized to the flat open-wafer flux, is the ARDE transmission the etch rate follows
for a radical-limited process.

The only honest question this gate answers: does the common engine's transport produce an ARDE floor-
flux COLLAPSE, and does that collapse CONVERGE under mesh / form-factor-ray / source-quadrature
refinement (the error-driven AMR requirement)? If it does not converge, no downstream de Boer number
means anything.

Run:  python scripts/deboer_arde_static.py
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from petch.boundary_state import (
    MaxwellianFluxVelocityDensity, PlasmaBoundaryState, SpeciesBoundaryState,
)
from petch.boundary_transport_3d import (
    estimate_diffuse_form_factors_3d, gather_boundary_state_ballistic_3d,
    trace_boundary_state_first_hit_3d,
)
from petch.feature_step_3d import (
    _face_material_ids, _surface_gas_normals, make_rectangular_trench_geometry_3d,
)
from petch.neutral_radiosity_3d import solve_diffuse_neutral_radiosity_3d
from petch.threed import extract_mesh_3d


def thermal_neutral_boundary_state(name, mass_amu, temperature_eV, flux_m2_s, *,
                                   n_transverse=5, n_normal=8, reference_plane_m=0.0):
    """Neutral radical half-Maxwellian flux quadrature (charge 0).

    Identical construction to ``maxwellian_electron_boundary_state`` -- Gauss-Hermite over the two
    tangential Maxwellians and Gauss-Laguerre over the normal-energy exponential -- so the discrete
    incident measure integrates the analytic kinetic flux density exactly. Nodes/weights are numerical
    quadrature only; the physical law is the cosine flux measure, no fitted angular closure.
    """
    temperature = float(temperature_eV)
    hermite_node, hermite_weight = np.polynomial.hermite.hermgauss(int(n_transverse))
    laguerre_node, laguerre_weight = np.polynomial.laguerre.laggauss(int(n_normal))
    ix, iy, iz = np.meshgrid(
        np.arange(hermite_node.size), np.arange(hermite_node.size),
        np.arange(laguerre_node.size), indexing="ij")
    velocity = np.column_stack((
        np.sqrt(temperature) * hermite_node[ix.ravel()],
        np.sqrt(temperature) * hermite_node[iy.ravel()],
        np.sqrt(temperature * laguerre_node[iz.ravel()]),
    ))
    weight = (hermite_weight[ix.ravel()] * hermite_weight[iy.ravel()]
              * laguerre_weight[iz.ravel()] / np.pi)
    species = SpeciesBoundaryState(
        name=name, charge_number=0, mass_amu=float(mass_amu),
        flux_m2_s=float(flux_m2_s), velocity_sqrt_eV=velocity, weight=weight,
        density_model=MaxwellianFluxVelocityDensity(temperature),
        provenance={"model": "analytic_half_maxwellian_flux_neutral"})
    return PlasmaBoundaryState(species=(species,), reference_plane_m=float(reference_plane_m))


def thermal_neutral_qmc_boundary_state(name, mass_amu, temperature_eV, flux_m2_s, *,
                                       log2_samples=16, seed=0, reference_plane_m=0.0):
    """Neutral radical source as N = 2**log2_samples scrambled-Sobol samples of the SAME analytic
    half-Maxwellian flux density.

    Unlike the tensor Gauss-Hermite quadrature, QMC samples concentrate by the physical cosine flux
    measure, so the near-vertical acceptance cone at high AR is resolved by raising N (error-driven,
    the Monte-Carlo AMR the ARDE literature prescribes) rather than by an A-scaled tensor grid. Equal
    weights; the physical law is unchanged.
    """
    from scipy.stats import qmc
    density = MaxwellianFluxVelocityDensity(float(temperature_eV))
    unit = qmc.Sobol(3, scramble=True, seed=int(seed)).random_base2(int(log2_samples))
    velocity = density.sample_flux_velocity(unit)
    weight = np.full(velocity.shape[0], 1.0 / velocity.shape[0])
    species = SpeciesBoundaryState(
        name=name, charge_number=0, mass_amu=float(mass_amu),
        flux_m2_s=float(flux_m2_s), velocity_sqrt_eV=velocity, weight=weight,
        density_model=density,
        provenance={"model": "qmc_half_maxwellian_flux_neutral", "log2_samples": int(log2_samples)})
    return PlasmaBoundaryState(species=(species,), reference_plane_m=float(reference_plane_m))


def floor_transmission(aspect_ratio, sticking, *, opening_um, dx_um, rays_per_face,
                       n_transverse, n_normal, mask_um=0.05, seed=0,
                       transport_method="adjoint", n_position=64, device=None,
                       source_method="quadrature", log2_samples=16):
    """Return the area-weighted floor incident flux / source flux for one (AR, s) point.

    Pure transport + conservation. ``s`` (sticking) is the only physical surface input; it is applied
    uniformly on every gas-facing face (substrate walls, floor, and mask), the single-sticking-
    coefficient Coburn-Winters molecular-flow model.

    ``transport_method``: ``"adjoint"`` uses the deterministic per-face gather (one Warp launch per
    angular atom -- cost-bound at high angular resolution); ``"forward"`` uses the first-hit tracer
    which batches every angular atom x source position into ONE Warp kernel (GPU-ready, so high
    angular resolution needed at high AR is affordable). ``device`` threads to Warp (e.g. "cuda").
    """
    etched = aspect_ratio * opening_um
    substrate_top = etched + max(4.0 * dx_um, 0.05)
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=2.0 * opening_um, cell_length=max(6.0 * dx_um, 0.06),
        domain_height=substrate_top + mask_um + max(6.0 * dx_um, 0.06),
        dx=dx_um, opening_width=opening_um, mask_thickness=mask_um,
        substrate_top=substrate_top, etched_depth=etched)
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    gas_normals = _surface_gas_normals(verts, faces, centroids, geometry)
    domain_size = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(domain_size[2])
    source_flux = 1.0e20
    ref_m = source_z * geometry.mesh_length_unit_m
    if source_method == "qmc":
        boundary = thermal_neutral_qmc_boundary_state(
            "F", 19.0, 0.05, source_flux, log2_samples=log2_samples, seed=seed,
            reference_plane_m=ref_m)
    else:
        boundary = thermal_neutral_boundary_state(
            "F", 19.0, 0.05, source_flux, n_transverse=n_transverse, n_normal=n_normal,
            reference_plane_m=ref_m)

    source_bounds = (0.0, float(domain_size[0]), 0.0, float(domain_size[1]))
    if transport_method == "forward":
        transport = trace_boundary_state_first_hit_3d(
            boundary, {"F": "neutral_reactant"}, verts, faces, areas,
            source_bounds=source_bounds, source_z=source_z,
            mesh_length_unit_m=geometry.mesh_length_unit_m,
            mesh_origin_m=geometry.mesh_origin_m, n_position=n_position,
            periodic_lateral=True, domain_size=domain_size, device=device)
    else:
        transport = gather_boundary_state_ballistic_3d(
            boundary, {"F": "neutral_reactant"}, verts, faces, areas, centroids, gas_normals,
            source_bounds=source_bounds, source_z=source_z,
            mesh_length_unit_m=geometry.mesh_length_unit_m,
            mesh_origin_m=geometry.mesh_origin_m, face_quadrature_points=3,
            periodic_lateral=True, domain_size=domain_size, ray_offset=1e-3 * geometry.dx,
            device=device)
    direct = np.asarray(transport.surface_fluxes.neutral_flux_m2_s["F"], dtype=float)

    factors = estimate_diffuse_form_factors_3d(
        verts, faces, centroids, gas_normals, rays_per_face=rays_per_face, seed=seed,
        domain_size=domain_size, periodic_lateral=True, ray_offset=1e-3 * geometry.dx)
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    reaction = np.full(len(faces), float(sticking))
    solution = solve_diffuse_neutral_radiosity_3d(
        direct, physical_area, factors.source_face, factors.target_face,
        factors.transfer_fraction, factors.escape_fraction, reaction)
    incident = solution.incident_flux_m2_s

    # Floor faces: the etchable-substrate faces at the bottom of the trench (upward gas normal).
    material = _face_material_ids(centroids, geometry)
    is_substrate = material == 1
    upward = gas_normals[:, 2] > 0.5
    if not np.any(is_substrate & upward):
        raise RuntimeError("no upward substrate floor faces found")
    floor_z = centroids[is_substrate & upward, 2].min()
    floor = is_substrate & upward & (centroids[:, 2] <= floor_z + geometry.dx)
    floor_flux = float(np.dot(incident[floor], physical_area[floor]) / physical_area[floor].sum())
    return {
        "aspect_ratio": float(aspect_ratio),
        "sticking": float(sticking),
        "floor_transmission": floor_flux / source_flux,
        "relative_balance_error": float(solution.relative_balance_error),
        "n_faces": int(len(faces)),
        "n_floor_faces": int(floor.sum()),
    }


def geometric_slot_transmission(aspect_ratio):
    """Analytic s=1 target: cosine mouth -> absorbing-wall floor line-of-sight transmission.

    For an infinite 2D slot of width w and depth L (A = L/w), the Hottel crossed-strings view factor
    between the mouth strip and the directly opposed floor strip is ``sqrt(1+A^2) - A``. With s=1 every
    wall hit is absorbed, so the floor flux is exactly this direct view factor. Asymptote ~1/(2A).
    This is a first-principles geometric identity, not a fit.
    """
    a = float(aspect_ratio)
    return float(np.sqrt(1.0 + a * a) - a)


def converged_floor_transmission(aspect_ratio, sticking, *, opening_um, dx_um, rays_per_face,
                                 nt_schedule=(8, 12, 16, 24, 32, 48), n_normal=10, rel_tol=0.03):
    """Adaptive angular refinement (AMR in phase space).

    The floor-reaching acceptance cone has half-angle ~arctan(1/A), so a fixed angular quadrature
    aliases it once the node spacing exceeds 1/A (Coburn-Winters conductance limit; confirmed by the
    flatline artifact this replaces). Refine the transverse-velocity quadrature order until the floor
    transmission stops moving by more than ``rel_tol`` (error-driven stop, not an AR-shaped rule).
    """
    prev = None
    history = []
    for nt in nt_schedule:
        point = floor_transmission(
            aspect_ratio, sticking, opening_um=opening_um, dx_um=dx_um,
            rays_per_face=rays_per_face, n_transverse=nt, n_normal=n_normal)
        t = point["floor_transmission"]
        history.append((nt, t))
        if prev is not None and abs(t - prev) <= rel_tol * max(t, 1e-9):
            return {"transmission": t, "n_transverse_used": nt, "converged": True,
                    "history": history, "n_faces": point["n_faces"]}
        prev = t
    return {"transmission": prev, "n_transverse_used": nt_schedule[-1], "converged": False,
            "history": history, "n_faces": history[-1] if history else None}


def validate_geometric(aspect_ratios=(1.0, 2.0, 4.0, 8.0), *, opening_um=0.10, mask_um=0.05,
                       dx_um=0.02, nt_schedule=(8, 16, 24, 32, 48, 64, 96)):
    """s=1 pure-shadowing gate: converged engine transmission must approach the slot view factor.

    The aperture is the mask top, so the flux-limiting slot runs mask+trench deep: the effective
    aspect ratio is A_eff = A + mask/opening. The analytic target is sqrt(1+A_eff^2) - A_eff.
    """
    print("s=1 pure geometric shadowing: converged engine vs analytic view factor of the")
    print(f"mask+trench slot (A_eff = A + mask/opening = A + {mask_um/opening_um:.2f})")
    print(f"{'AR':>5} {'A_eff':>6} {'engine':>9} {'analytic':>9} {'ratio':>7} {'nt':>4} {'conv':>5}")
    rows = []
    for ar in aspect_ratios:
        r = converged_floor_transmission(
            ar, 1.0, opening_um=opening_um, dx_um=dx_um, rays_per_face=64,
            nt_schedule=nt_schedule)
        a_eff = ar + mask_um / opening_um
        analytic = geometric_slot_transmission(a_eff)
        ratio = r["transmission"] / analytic
        rows.append((ar, a_eff, r["transmission"], analytic, ratio,
                     r["n_transverse_used"], r["converged"]))
        print(f"{ar:>5.1f} {a_eff:>6.2f} {r['transmission']:>9.4f} {analytic:>9.4f} "
              f"{ratio:>7.3f} {r['n_transverse_used']:>4d} {str(r['converged']):>5}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-geometric", action="store_true",
                        help="run the s=1 analytic view-factor gate and exit")
    parser.add_argument("--aspect-ratios", type=float, nargs="+",
                        default=[0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    parser.add_argument("--sticking", type=float, nargs="+", default=[0.01, 0.05, 0.2])
    parser.add_argument("--opening-um", type=float, default=0.10)
    parser.add_argument("--dx-um", type=float, default=0.0125)
    parser.add_argument("--rays-per-face", type=int, default=256)
    parser.add_argument("--n-transverse", type=int, default=5)
    parser.add_argument("--n-normal", type=int, default=8)
    args = parser.parse_args()

    if args.validate_geometric:
        validate_geometric()
        return

    print(f"{'AR':>6}", end="")
    for s in args.sticking:
        print(f"  s={s:<6.3g}", end="")
    print("   (floor transmission = floor flux / source flux)")
    curves = {s: [] for s in args.sticking}
    for ar in args.aspect_ratios:
        print(f"{ar:>6.2f}", end="", flush=True)
        for s in args.sticking:
            point = floor_transmission(
                ar, s, opening_um=args.opening_um, dx_um=args.dx_um,
                rays_per_face=args.rays_per_face, n_transverse=args.n_transverse,
                n_normal=args.n_normal)
            curves[s].append(point)
            print(f"  {point['floor_transmission']:<8.4f}", end="", flush=True)
        print()
    out = {
        "campaign": "deboer_arde_static_flux_gate",
        "engine": "feature-3d (ballistic gather + diffuse radiosity)",
        "opening_um": args.opening_um, "dx_um": args.dx_um,
        "rays_per_face": args.rays_per_face,
        "n_transverse": args.n_transverse, "n_normal": args.n_normal,
        "curves": {str(s): curves[s] for s in args.sticking},
    }
    print("\n" + json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
