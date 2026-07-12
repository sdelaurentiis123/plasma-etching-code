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


@dataclass(frozen=True)
class SteadyDielectricCharging3DResult:
    charge_node_c: np.ndarray
    potential_v: np.ndarray
    positive_current_node_a: np.ndarray
    negative_current_node_a: np.ndarray
    transport: BoundaryTransport3DResult
    poisson: PoissonDiagnostics3D
    history: tuple[Mapping[str, float], ...]
    converged: bool
    rejected_steps: int
    known_limitations: tuple[str, ...]

    def __post_init__(self):
        for name in (
                "charge_node_c", "potential_v", "positive_current_node_a",
                "negative_current_node_a"):
            array = np.asarray(getattr(self, name), dtype=float).copy()
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        object.__setattr__(
            self, "history", tuple(MappingProxyType(dict(item)) for item in self.history))
        object.__setattr__(self, "known_limitations", tuple(self.known_limitations))


class DielectricChargingConvergenceError(RuntimeError):
    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


def _coupled_transport_limitations(transport):
    # The low-level trajectory call correctly says its voltage was supplied. At this coupling level
    # that voltage came from the current charge state, so retaining that line would misreport scope.
    return tuple(
        limitation for limitation in transport.known_limitations
        if limitation != "nodal potential is supplied rather than self-consistently charged")


def _coordinate_spacing_3d(poisson_system, potential_spacing, mesh_length_unit_m):
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
    return coordinate_spacing


def _evaluate_incident_current_3d(
        poisson_system, charge, boundary, verts, faces, areas, *, source_bounds, source_z,
        potential_origin, coordinate_spacing, mesh_length_unit_m, mesh_origin_m,
        n_position, seed, trajectory_fixed_dt, trajectory_max_steps, transport_device):
    charged_species = tuple(species for species in boundary.species if species.charge_number != 0)
    if not charged_species:
        raise ValueError("dielectric charging requires at least one charged boundary species")
    charged_boundary = PlasmaBoundaryState(
        charged_species, boundary.reference_plane_m, provenance=boundary.provenance)
    species_role = {species.name: "energetic_bombardment" for species in charged_species}
    potential, poisson = poisson_system.solve(charge)
    transport = trace_boundary_state_field_3d(
        charged_boundary, species_role, verts, faces, areas,
        source_bounds=source_bounds, source_z=source_z,
        nodal_potential_v=potential, potential_origin=potential_origin,
        potential_spacing=coordinate_spacing,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        n_position=n_position, seed=seed, fixed_dt=trajectory_fixed_dt,
        max_steps=trajectory_max_steps, device=transport_device)

    population_by_name = {
        population.name: population for population in transport.surface_fluxes.energetic_fluxes}
    if set(population_by_name) != set(species_role):
        raise RuntimeError("charged transport did not preserve every species current measure")
    face_count = np.asarray(faces).shape[0]
    positive_face_current = np.zeros(face_count)
    negative_face_current = np.zeros(face_count)
    for species in charged_species:
        population = population_by_name[species.name]
        if not isinstance(population, FaceResolvedEnergeticFlux):
            raise RuntimeError("3-D charging requires face-resolved incident events")
        current = (ECHARGE * abs(float(species.charge_number))
                   * population.flux_m2_s)
        if species.charge_number > 0:
            positive_face_current += current
        else:
            negative_face_current += current
    projection_arguments = dict(
        shape=poisson_system.shape, vertices=verts, faces=faces,
        grid_origin=potential_origin, grid_spacing=coordinate_spacing,
        coordinate_length_unit_m=mesh_length_unit_m)
    positive_node_current = lump_triangle_sheet_charge_3d(
        sigma_c_per_m2=positive_face_current, **projection_arguments)
    negative_node_current = lump_triangle_sheet_charge_3d(
        sigma_c_per_m2=negative_face_current, **projection_arguments)
    return dict(
        potential=potential, poisson=poisson, transport=transport,
        positive_face_current=positive_face_current,
        negative_face_current=negative_face_current,
        positive_node_current=positive_node_current,
        negative_node_current=negative_node_current)


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
    coordinate_spacing = _coordinate_spacing_3d(
        poisson_system, potential_spacing, mesh_length_unit_m)
    evaluated = _evaluate_incident_current_3d(
        poisson_system, charge, boundary, verts, faces, areas,
        source_bounds=source_bounds, source_z=source_z,
        potential_origin=potential_origin, coordinate_spacing=coordinate_spacing,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        n_position=n_position, seed=seed, trajectory_fixed_dt=trajectory_fixed_dt,
        trajectory_max_steps=trajectory_max_steps, transport_device=transport_device)
    face_current = evaluated["positive_face_current"] - evaluated["negative_face_current"]
    current_node = evaluated["positive_node_current"] - evaluated["negative_node_current"]
    charge_increment = current_node * float(duration_s)
    incident_node_current = (
        evaluated["positive_node_current"] + evaluated["negative_node_current"])
    fixed_increment = float(np.sum(
        incident_node_current[poisson_system.dirichlet_mask])) * float(duration_s)
    total_increment = float(np.sum(incident_node_current)) * float(duration_s)
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
        potential_before_v=evaluated["potential"],
        potential_after_v=potential_after,
        face_current_density_a_m2=face_current,
        transport=evaluated["transport"],
        poisson_before=evaluated["poisson"],
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
        ) + _coupled_transport_limitations(evaluated["transport"]))


def solve_dielectric_charging_steady_3d(
        poisson_system: NodalPoissonSystem3D, initial_charge_node_c,
        boundary: PlasmaBoundaryState, verts, faces, areas, *, source_bounds, source_z,
        potential_origin, potential_spacing, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), n_position=256, seed=0,
        trajectory_fixed_dt=0.01, trajectory_max_steps=10000, transport_device=None,
        max_iter=30, min_iter=2, current_balance_tol=1e-3,
        beta=0.5, response_energy_eV=4.0, maximum_voltage_step=8.0,
        trust_growth_tolerance=0.02, minimum_beta=1e-4, require_converged=True):
    """Solve local steady dielectric current balance on the compatible 3-D charge basis.

    The physical residual is ``abs(I+ - I-)/(I+ + I-)`` at every active surface node. Diagonal
    capacitance, ``beta``, and ``response_energy_eV`` precondition the nonlinear solve but do not alter
    that root. Trial steps that increase the RMS physical residual are rejected and retried at half gain.
    """
    if not isinstance(poisson_system, NodalPoissonSystem3D):
        raise TypeError("poisson_system must be a NodalPoissonSystem3D")
    charge = np.asarray(initial_charge_node_c, dtype=float).copy()
    if charge.shape != poisson_system.shape or not np.all(np.isfinite(charge)):
        raise ValueError("initial_charge_node_c must be a finite grid matching poisson_system")
    if np.any(np.abs(charge[poisson_system.dirichlet_mask]) > 0.0):
        raise ValueError("stored dielectric charge cannot be assigned to Dirichlet reservoir nodes")
    if (int(max_iter) != max_iter or int(min_iter) != min_iter or max_iter <= 0 or min_iter <= 0
            or min_iter > max_iter or not np.isfinite(current_balance_tol)
            or current_balance_tol <= 0.0 or not np.isfinite(beta) or beta <= 0.0
            or not np.isfinite(response_energy_eV) or response_energy_eV <= 0.0
            or not np.isfinite(maximum_voltage_step) or maximum_voltage_step <= 0.0
            or not np.isfinite(trust_growth_tolerance) or trust_growth_tolerance < 0.0
            or not np.isfinite(minimum_beta) or minimum_beta <= 0.0 or minimum_beta > beta):
        raise ValueError("invalid steady charging solver controls")
    if (not any(species.charge_number > 0 for species in boundary.species)
            or not any(species.charge_number < 0 for species in boundary.species)):
        raise ValueError("steady dielectric charging requires positive and negative incident species")
    coordinate_spacing = _coordinate_spacing_3d(
        poisson_system, potential_spacing, mesh_length_unit_m)
    support = lump_triangle_sheet_charge_3d(
        poisson_system.shape, verts, faces, np.ones(np.asarray(faces).shape[0]),
        grid_origin=potential_origin, grid_spacing=coordinate_spacing,
        coordinate_length_unit_m=mesh_length_unit_m)
    support_mask = np.abs(support) > 1e-14 * float(np.max(np.abs(support)))
    if np.any(support_mask & poisson_system.dirichlet_mask):
        raise ValueError(
            "dielectric surface projects onto a Dirichlet reservoir; mixed surface equations required")
    support_nodes = np.column_stack(np.where(support_mask))
    if support_nodes.size == 0:
        raise ValueError("dielectric surface has no supported Poisson nodes")
    capacitance = poisson_system.diagonal_capacitance(support_nodes)
    evaluate_arguments = dict(
        poisson_system=poisson_system, boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=source_bounds, source_z=source_z,
        potential_origin=potential_origin, coordinate_spacing=coordinate_spacing,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        n_position=n_position, seed=seed, trajectory_fixed_dt=trajectory_fixed_dt,
        trajectory_max_steps=trajectory_max_steps, transport_device=transport_device)

    beta_current = float(beta); rejected_steps = 0; history = []

    def assess(state_charge):
        evaluated = _evaluate_incident_current_3d(charge=state_charge, **evaluate_arguments)
        positive = evaluated["positive_node_current"][tuple(support_nodes.T)]
        negative = evaluated["negative_node_current"][tuple(support_nodes.T)]
        total = positive + negative
        scale = float(np.max(total)) if total.size else 0.0
        active = total > max(1e-15 * scale, 1e-300)
        if not np.any(active):
            raise RuntimeError("steady charging has no resolved incident current on its surface")
        relative = np.zeros_like(total)
        relative[active] = np.abs(positive[active] - negative[active]) / total[active]
        current_floor = max(1e-15 * scale, 1e-300)
        log_ratio = np.log(
            np.maximum(positive, current_floor) / np.maximum(negative, current_floor))
        merit = float(np.sqrt(np.mean(relative[active] ** 2)))
        maximum = float(np.max(relative[active]))
        return evaluated, positive, negative, active, log_ratio, merit, maximum

    evaluated, positive, negative, active, log_ratio, merit, maximum = assess(charge)
    history.append(dict(
        iteration=1, rms_relative_current_imbalance=merit,
        max_relative_current_imbalance=maximum, beta=beta_current,
        mean_surface_voltage_v=float(np.mean(
            evaluated["potential"][tuple(support_nodes.T)]))))
    while len(history) < int(max_iter) and not (
            len(history) >= int(min_iter) and maximum <= float(current_balance_tol)):
        voltage_step = np.clip(
            beta_current * float(response_energy_eV) * log_ratio,
            -float(maximum_voltage_step), float(maximum_voltage_step))
        voltage_step[~active] = 0.0
        trial_charge = charge.copy()
        trial_charge[tuple(support_nodes.T)] += capacitance * voltage_step
        trial = assess(trial_charge)
        trial_merit = trial[5]
        if trial_merit > merit * (1.0 + float(trust_growth_tolerance)):
            beta_current *= 0.5; rejected_steps += 1
            if beta_current < float(minimum_beta):
                break
            continue
        accepted_beta = beta_current
        if trial_merit < 0.8 * merit:
            beta_current = min(float(beta), 1.2 * beta_current)
        charge = trial_charge
        evaluated, positive, negative, active, log_ratio, merit, maximum = trial
        history.append(dict(
            iteration=len(history) + 1, rms_relative_current_imbalance=merit,
            max_relative_current_imbalance=maximum, beta=accepted_beta,
            mean_surface_voltage_v=float(np.mean(
                evaluated["potential"][tuple(support_nodes.T)]))))

    converged = bool(
        len(history) >= int(min_iter) and maximum <= float(current_balance_tol))
    positive_grid = np.zeros(poisson_system.shape)
    negative_grid = np.zeros(poisson_system.shape)
    positive_grid[tuple(support_nodes.T)] = positive
    negative_grid[tuple(support_nodes.T)] = negative
    result = SteadyDielectricCharging3DResult(
        charge_node_c=charge, potential_v=evaluated["potential"],
        positive_current_node_a=positive_grid, negative_current_node_a=negative_grid,
        transport=evaluated["transport"], poisson=evaluated["poisson"],
        history=tuple(history), converged=converged, rejected_steps=rejected_steps,
        known_limitations=(
            "all supplied surface triangles are treated as charge-storing dielectric",
            "fixed deterministic launch quadrature requires an external sample-refinement ladder",
            "no secondary-electron emission, reflection, leakage, or surface conduction",
            "no floating-conductor circuit equations",
        ) + _coupled_transport_limitations(evaluated["transport"]))
    if require_converged and not converged:
        raise DielectricChargingConvergenceError(
            f"3-D dielectric current balance did not converge in {len(history)} accepted iterations; "
            f"max relative imbalance={maximum:.6g}", result)
    return result
