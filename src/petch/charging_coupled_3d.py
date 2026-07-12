"""One conservative physical-time update of 3-D dielectric feature charging."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import PlasmaBoundaryState
from .boundary_transport_3d import BoundaryTransport3DResult, trace_boundary_state_field_3d
from .charging_poisson_3d import (
    NodalPoissonSystem3D,
    PoissonDiagnostics3D,
    lump_triangle_sheet_charge_3d,
)
from .sheath import ECHARGE
from .surface_kinetics import FaceResolvedEnergeticFlux


@dataclass(frozen=True)
class DielectricChargingStep3DResult:
    charge_node_c: np.ndarray
    charge_increment_node_c: np.ndarray
    potential_before_v: np.ndarray
    potential_after_v: np.ndarray
    face_current_density_a_m2: np.ndarray
    transport: BoundaryTransport3DResult
    poisson_before: PoissonDiagnostics3D
    poisson_after: PoissonDiagnostics3D
    diagnostics: Mapping[str, float]
    known_limitations: tuple[str, ...]

    def __post_init__(self):
        for name in (
                "charge_node_c", "charge_increment_node_c", "potential_before_v",
                "potential_after_v", "face_current_density_a_m2"):
            array = np.asarray(getattr(self, name), dtype=float).copy()
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))
        object.__setattr__(self, "known_limitations", tuple(self.known_limitations))


def advance_dielectric_charging_3d(
        poisson_system: NodalPoissonSystem3D, charge_node_c, boundary: PlasmaBoundaryState,
        verts, faces, areas, *, source_bounds, source_z, potential_origin,
        potential_spacing, duration_s, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), n_position=256, seed=0,
        trajectory_fixed_dt=0.01, trajectory_max_steps=10000, transport_device=None):
    """Advance stored dielectric charge by the signed incident-particle current.

    The sequence is charge -> Q1 Poisson voltage -> collisionless charged-particle trajectories ->
    signed face current -> compatible Q1 charge projection -> updated Poisson voltage. Every supplied
    triangle is treated as a charge-storing dielectric surface. Dirichlet nodes are external reservoirs,
    so depositing surface charge onto one is refused instead of silently discarding it.

    This is a physical-time forward-Euler update, not the accelerated steady current-balance solve.
    ``duration_s`` must therefore resolve the charging transient selected by the caller.
    """
    if not isinstance(poisson_system, NodalPoissonSystem3D):
        raise TypeError("poisson_system must be a NodalPoissonSystem3D")
    if not np.isfinite(duration_s) or duration_s <= 0.0:
        raise ValueError("duration_s must be finite and positive")
    charge = np.asarray(charge_node_c, dtype=float)
    if charge.shape != poisson_system.shape or not np.all(np.isfinite(charge)):
        raise ValueError("charge_node_c must be a finite grid matching poisson_system")
    if np.any(np.abs(charge[poisson_system.dirichlet_mask]) > 0.0):
        raise ValueError("stored dielectric charge cannot be assigned to Dirichlet reservoir nodes")
    coordinate_spacing = np.asarray(potential_spacing, dtype=float)
    if coordinate_spacing.ndim == 0:
        coordinate_spacing = np.full(3, float(coordinate_spacing))
    expected_spacing_m = coordinate_spacing * float(mesh_length_unit_m)
    if (coordinate_spacing.shape != (3,) or np.any(~np.isfinite(coordinate_spacing))
            or np.any(coordinate_spacing <= 0.0)
            or not np.allclose(
                poisson_system.spacing_m, expected_spacing_m, rtol=1e-12, atol=0.0)):
        raise ValueError(
            "Poisson physical spacing must equal potential_spacing * mesh_length_unit_m")

    charged_species = tuple(species for species in boundary.species if species.charge_number != 0)
    if not charged_species:
        raise ValueError("dielectric charging requires at least one charged boundary species")
    charged_boundary = PlasmaBoundaryState(
        charged_species, boundary.reference_plane_m, provenance=boundary.provenance)
    species_role = {species.name: "energetic_bombardment" for species in charged_species}
    potential_before, poisson_before = poisson_system.solve(charge)
    transport = trace_boundary_state_field_3d(
        charged_boundary, species_role, verts, faces, areas,
        source_bounds=source_bounds, source_z=source_z,
        nodal_potential_v=potential_before, potential_origin=potential_origin,
        potential_spacing=coordinate_spacing,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        n_position=n_position, seed=seed, fixed_dt=trajectory_fixed_dt,
        max_steps=trajectory_max_steps, device=transport_device)

    population_by_name = {
        population.name: population for population in transport.surface_fluxes.energetic_fluxes}
    if set(population_by_name) != set(species_role):
        raise RuntimeError("charged transport did not preserve every species current measure")
    face_count = np.asarray(faces).shape[0]
    face_current = np.zeros(face_count)
    for species in charged_species:
        population = population_by_name[species.name]
        if not isinstance(population, FaceResolvedEnergeticFlux):
            raise RuntimeError("3-D charging requires face-resolved incident events")
        face_current += ECHARGE * float(species.charge_number) * population.flux_m2_s
    sheet_charge_increment = face_current * float(duration_s)
    charge_increment = lump_triangle_sheet_charge_3d(
        poisson_system.shape, verts, faces, sheet_charge_increment,
        grid_origin=potential_origin, grid_spacing=coordinate_spacing,
        coordinate_length_unit_m=mesh_length_unit_m)
    fixed_increment = float(np.sum(np.abs(charge_increment[poisson_system.dirichlet_mask])))
    total_increment = float(np.sum(np.abs(charge_increment)))
    if fixed_increment > 1e-13 * max(total_increment, 1e-300):
        raise ValueError(
            "incident charge projects onto a Dirichlet reservoir; mixed dielectric/conductor "
            "surface handling must be specified explicitly")
    updated_charge = charge + charge_increment
    potential_after, poisson_after = poisson_system.solve(updated_charge)

    areas = np.asarray(areas, dtype=float)
    incident_charge = float(np.dot(
        face_current, areas * float(mesh_length_unit_m) ** 2) * float(duration_s))
    deposited_charge = float(np.sum(charge_increment))
    conservation_residual = deposited_charge - incident_charge
    return DielectricChargingStep3DResult(
        charge_node_c=updated_charge,
        charge_increment_node_c=charge_increment,
        potential_before_v=potential_before,
        potential_after_v=potential_after,
        face_current_density_a_m2=face_current,
        transport=transport,
        poisson_before=poisson_before,
        poisson_after=poisson_after,
        diagnostics=dict(
            duration_s=float(duration_s),
            incident_charge_c=incident_charge,
            deposited_charge_c=deposited_charge,
            charge_conservation_residual_c=conservation_residual,
            maximum_abs_face_current_density_a_m2=float(np.max(np.abs(face_current)))),
        known_limitations=(
            "all supplied surface triangles are treated as charge-storing dielectric",
            "physical-time forward-Euler charge update requires timestep convergence",
            "no secondary-electron emission, reflection, leakage, or surface conduction",
            "no floating-conductor circuit equations",
        ) + tuple(transport.known_limitations))
