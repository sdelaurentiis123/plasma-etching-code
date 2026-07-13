"""One conservative physical-time update of 3-D dielectric feature charging."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import PlasmaBoundaryState
from .boundary_transport_3d import (
    BoundaryTransport3DResult, gather_boundary_state_field_adjoint_3d,
    merge_boundary_transport_results_3d, trace_boundary_state_bidirectional_field_3d,
    trace_boundary_state_field_3d,
)
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
    net_current_stderr_node_a: np.ndarray
    transport: BoundaryTransport3DResult
    poisson: PoissonDiagnostics3D
    history: tuple[Mapping[str, float], ...]
    converged: bool
    rejected_steps: int
    known_limitations: tuple[str, ...]

    def __post_init__(self):
        for name in (
                "charge_node_c", "potential_v", "positive_current_node_a",
                "negative_current_node_a", "net_current_stderr_node_a"):
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


class BidirectionalCurrentCertificationError(RuntimeError):
    """A sampled current map failed its per-face precision or consistency contract."""


def _anderson_step(x, residual, x_history, residual_history, gain, depth):
    """Type-II Anderson acceleration for a preconditioned fixed-point residual."""
    x = np.asarray(x, dtype=float); residual = np.asarray(residual, dtype=float)
    if x.shape != residual.shape or x.ndim != 1:
        raise ValueError("Anderson state and residual must be matching vectors")
    if int(depth) != depth or depth <= 0 or not np.isfinite(gain) or gain <= 0.0:
        raise ValueError("Anderson depth and gain must be positive")
    x_history.append(x.copy()); residual_history.append(residual.copy())
    if len(x_history) > int(depth) + 1:
        x_history.pop(0); residual_history.pop(0)
    step = float(gain) * residual
    if len(residual_history) >= 2:
        delta_residual = np.stack([
            residual_history[index + 1] - residual_history[index]
            for index in range(len(residual_history) - 1)], axis=1)
        delta_x = np.stack([
            x_history[index + 1] - x_history[index]
            for index in range(len(x_history) - 1)], axis=1)
        gamma, *_ = np.linalg.lstsq(delta_residual, residual, rcond=1e-8)
        step = step - (delta_x + float(gain) * delta_residual) @ gamma
    return step


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


def _validate_transport_estimators_3d(
        boundary, faces, transport_estimator, face_centroids, face_gas_normals):
    charged_names = {species.name for species in boundary.species if species.charge_number != 0}
    estimator_by_name = (
        {name: transport_estimator for name in charged_names}
        if isinstance(transport_estimator, str) else dict(transport_estimator))
    if (not charged_names or set(estimator_by_name) != charged_names
            or any(value not in {"forward", "adjoint", "bidirectional"}
                   for value in estimator_by_name.values())):
        raise ValueError(
            "transport_estimator must select forward, adjoint, or bidirectional for every charged species")
    if {"adjoint", "bidirectional"} & set(estimator_by_name.values()):
        centroids = np.asarray(face_centroids, dtype=float)
        normals = np.asarray(face_gas_normals, dtype=float)
        if (centroids.shape != (len(faces), 3) or normals.shape != centroids.shape
                or np.any(~np.isfinite(centroids)) or np.any(~np.isfinite(normals))):
            raise ValueError("adjoint charging requires finite centroid and gas-normal arrays per face")
    return estimator_by_name


def _validate_adjoint_proposal_frames_3d(estimator_by_name, adjoint_proposal_frames):
    adjoint_names = {
        name for name, estimator in estimator_by_name.items()
        if estimator in {"adjoint", "bidirectional"}}
    proposal_frames = (
        {name: adjoint_proposal_frames for name in adjoint_names}
        if isinstance(adjoint_proposal_frames, str) else dict(adjoint_proposal_frames))
    if (set(proposal_frames) != adjoint_names
            or any(value not in {"surface_local", "source_aligned"}
                   for value in proposal_frames.values())):
        raise ValueError(
            "adjoint_proposal_frames must select surface_local or source_aligned for adjoint species")
    return proposal_frames


def _evaluate_incident_current_3d(
        poisson_system, charge, boundary, verts, faces, areas, *, source_bounds, source_z,
        potential_origin, coordinate_spacing, mesh_length_unit_m, mesh_origin_m,
        n_position, seed, trajectory_fixed_dt, trajectory_max_steps,
        phase_space_log2_samples, periodic_lateral, transport_estimator,
        face_centroids, face_gas_normals, adjoint_face_quadrature_points,
        adjoint_ray_offset, adjoint_proposals, adjoint_proposal_frames,
        bidirectional_options, transport_device):
    charged_species = tuple(species for species in boundary.species if species.charge_number != 0)
    if not charged_species:
        raise ValueError("dielectric charging requires at least one charged boundary species")
    charged_boundary = PlasmaBoundaryState(
        charged_species, boundary.reference_plane_m, provenance=boundary.provenance)
    species_role = {species.name: "energetic_bombardment" for species in charged_species}
    potential, poisson = poisson_system.solve(charge)
    estimator_by_name = (
        {species.name: str(transport_estimator) for species in charged_species}
        if isinstance(transport_estimator, str) else dict(transport_estimator))
    transports = []; bidirectional_method_hint = {}
    for estimator in ("forward", "adjoint", "bidirectional"):
        selected = tuple(
            species for species in charged_species
            if estimator_by_name.get(species.name) == estimator)
        if not selected:
            continue
        selected_boundary = PlasmaBoundaryState(
            selected, charged_boundary.reference_plane_m, provenance=boundary.provenance)
        selected_role = {species.name: species_role[species.name] for species in selected}
        if estimator == "bidirectional":
            options = {} if bidirectional_options is None else dict(bidirectional_options)
            proposal_subset = (None if adjoint_proposals is None else {
                species.name: adjoint_proposals[species.name] for species in selected
                if species.name in adjoint_proposals})
            frame_subset = (
                adjoint_proposal_frames if isinstance(adjoint_proposal_frames, str) else {
                    species.name: adjoint_proposal_frames[species.name] for species in selected})
            bidirectional = trace_boundary_state_bidirectional_field_3d(
                selected_boundary, selected_role, verts, faces, areas,
                face_centroids, face_gas_normals,
                source_bounds=source_bounds, source_z=source_z,
                nodal_potential_v=potential, potential_origin=potential_origin,
                potential_spacing=coordinate_spacing,
                mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
                seed=seed, fixed_dt=trajectory_fixed_dt, max_steps=trajectory_max_steps,
                periodic_lateral=periodic_lateral, proposal_by_species=proposal_subset,
                proposal_frame_by_species=frame_subset, device=transport_device, **options)
            if not all(item.converged for item in bidirectional.selection_by_species.values()):
                failed = [name for name, item in bidirectional.selection_by_species.items()
                          if not item.converged]
                raise BidirectionalCurrentCertificationError(
                    f"bidirectional current estimator did not certify species {failed}")
            bidirectional_method_hint.update({
                name: item.method.copy()
                for name, item in bidirectional.selection_by_species.items()})
            transports.append(bidirectional.transport)
        elif estimator == "adjoint":
            proposal_subset = (None if adjoint_proposals is None else {
                species.name: adjoint_proposals[species.name] for species in selected
                if species.name in adjoint_proposals})
            transports.append(gather_boundary_state_field_adjoint_3d(
                selected_boundary, selected_role, verts, faces, areas,
                face_centroids, face_gas_normals,
                source_bounds=source_bounds, source_z=source_z,
                nodal_potential_v=potential, potential_origin=potential_origin,
                potential_spacing=coordinate_spacing,
                mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
                face_quadrature_points=adjoint_face_quadrature_points,
                ray_offset=adjoint_ray_offset, fixed_dt=trajectory_fixed_dt,
                max_steps=trajectory_max_steps, periodic_lateral=periodic_lateral,
                proposal_by_species=proposal_subset,
                proposal_frame_by_species=(
                    adjoint_proposal_frames if isinstance(adjoint_proposal_frames, str) else {
                        species.name: adjoint_proposal_frames[species.name]
                        for species in selected}),
                device=transport_device))
        else:
            transports.append(trace_boundary_state_field_3d(
                selected_boundary, selected_role, verts, faces, areas,
                source_bounds=source_bounds, source_z=source_z,
                nodal_potential_v=potential, potential_origin=potential_origin,
                potential_spacing=coordinate_spacing,
                mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
                n_position=n_position, seed=seed, fixed_dt=trajectory_fixed_dt,
                max_steps=trajectory_max_steps,
                phase_space_log2_samples=phase_space_log2_samples,
                periodic_lateral=periodic_lateral, device=transport_device))
    transport = (transports[0] if len(transports) == 1
                 else merge_boundary_transport_results_3d(*transports))

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
        negative_node_current=negative_node_current,
        bidirectional_method_hint=bidirectional_method_hint)


def advance_dielectric_charging_3d(
        poisson_system: NodalPoissonSystem3D, charge_node_c, boundary: PlasmaBoundaryState,
        verts, faces, areas, *, source_bounds, source_z, potential_origin,
        potential_spacing, duration_s, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), n_position=256, seed=0,
        trajectory_fixed_dt=0.01, trajectory_max_steps=10000,
        phase_space_log2_samples=None, periodic_lateral=False,
        transport_estimator="forward", face_centroids=None, face_gas_normals=None,
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-5,
        adjoint_proposals=None, adjoint_proposal_frames="surface_local",
        bidirectional_options=None,
        transport_device=None):
    """Advance stored dielectric charge by the signed incident-particle current.

    The sequence is charge -> Q1 Poisson voltage -> collisionless charged-particle trajectories ->
    signed face current -> compatible Q1 charge projection -> updated Poisson voltage. Every supplied
    triangle is treated as a charge-storing dielectric surface. Dirichlet nodes are external reservoirs,
    so depositing surface charge onto one is refused instead of silently discarding it.

    This is a physical-time forward-Euler update, not the accelerated steady current-balance solve.
    ``duration_s`` must therefore resolve the charging transient selected by the caller.
    ``transport_estimator`` may select ``"forward"`` or ``"adjoint"`` independently for each charged
    species; adjoint species require per-face centroids and gas normals. Adjoint proposal frames may
    likewise be selected per species for broad local incidence or narrow source-aligned incidence.
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
    estimator_by_name = _validate_transport_estimators_3d(
        boundary, faces, transport_estimator, face_centroids, face_gas_normals)
    _validate_adjoint_proposal_frames_3d(estimator_by_name, adjoint_proposal_frames)
    evaluated = _evaluate_incident_current_3d(
        poisson_system, charge, boundary, verts, faces, areas,
        source_bounds=source_bounds, source_z=source_z,
        potential_origin=potential_origin, coordinate_spacing=coordinate_spacing,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        n_position=n_position, seed=seed, trajectory_fixed_dt=trajectory_fixed_dt,
        trajectory_max_steps=trajectory_max_steps,
        phase_space_log2_samples=phase_space_log2_samples,
        periodic_lateral=periodic_lateral,
        transport_estimator=transport_estimator,
        face_centroids=face_centroids, face_gas_normals=face_gas_normals,
        adjoint_face_quadrature_points=adjoint_face_quadrature_points,
        adjoint_ray_offset=adjoint_ray_offset,
        adjoint_proposals=adjoint_proposals,
        adjoint_proposal_frames=adjoint_proposal_frames,
        bidirectional_options=bidirectional_options,
        transport_device=transport_device)
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
        trajectory_fixed_dt=0.01, trajectory_max_steps=10000,
        phase_space_log2_samples=None, periodic_lateral=False,
        transport_estimator="forward", face_centroids=None, face_gas_normals=None,
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-5,
        adjoint_proposals=None, adjoint_proposal_frames="surface_local",
        bidirectional_options=None,
        transport_device=None,
        max_iter=30, min_iter=2, current_balance_tol=1e-3,
        beta=0.5, response_energy_eV=4.0, maximum_voltage_step=8.0,
        trust_growth_tolerance=0.02, minimum_beta=1e-4,
        phase_space_replicates=1, current_confidence_sigma=2.0,
        phase_space_max_log2_samples=None, current_estimator_relative_tol=None,
        nonlinear_update="picard", anderson_depth=4,
        require_converged=True):
    """Solve local steady dielectric current balance on the compatible 3-D charge basis.

    The physical residual is ``abs(I+ - I-)/(I+ + I-)`` at every active surface node. The exact dense
    support-node Poisson response, ``beta``, and ``response_energy_eV`` precondition the nonlinear solve
    but do not alter that root. Trial steps that increase the RMS physical residual are rejected and
    retried at half gain. Forward or reversible-adjoint transport can be selected per charged species.
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
            or not np.isfinite(minimum_beta) or minimum_beta <= 0.0 or minimum_beta > beta
            or int(phase_space_replicates) != phase_space_replicates
            or phase_space_replicates <= 0 or not np.isfinite(current_confidence_sigma)
            or current_confidence_sigma <= 0.0
            or nonlinear_update not in {"picard", "anderson"}
            or int(anderson_depth) != anderson_depth or anderson_depth <= 0
            or int(adjoint_face_quadrature_points) not in {1, 3, 7}
            or not np.isfinite(adjoint_ray_offset) or adjoint_ray_offset <= 0.0):
        raise ValueError("invalid steady charging solver controls")
    estimator_by_name = _validate_transport_estimators_3d(
        boundary, faces, transport_estimator, face_centroids, face_gas_normals)
    _validate_adjoint_proposal_frames_3d(estimator_by_name, adjoint_proposal_frames)
    if phase_space_replicates > 1 and phase_space_log2_samples is None:
        raise ValueError(
            "multiple current replicates require joint continuous-density phase-space sampling")
    if (phase_space_max_log2_samples is None) != (current_estimator_relative_tol is None):
        raise ValueError(
            "adaptive current estimation requires both a maximum phase-space level and tolerance")
    if phase_space_max_log2_samples is not None:
        if (phase_space_log2_samples is None or phase_space_replicates < 2
                or int(phase_space_max_log2_samples) != phase_space_max_log2_samples
                or phase_space_max_log2_samples < phase_space_log2_samples
                or not np.isfinite(current_estimator_relative_tol)
                or current_estimator_relative_tol <= 0.0):
            raise ValueError("invalid adaptive current-estimator controls")
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
    voltage_response = poisson_system.voltage_response(support_nodes)
    evaluate_arguments = dict(
        poisson_system=poisson_system, boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=source_bounds, source_z=source_z,
        potential_origin=potential_origin, coordinate_spacing=coordinate_spacing,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        n_position=n_position, seed=seed, trajectory_fixed_dt=trajectory_fixed_dt,
        trajectory_max_steps=trajectory_max_steps,
        phase_space_log2_samples=phase_space_log2_samples,
        periodic_lateral=periodic_lateral,
        transport_estimator=transport_estimator,
        face_centroids=face_centroids, face_gas_normals=face_gas_normals,
        adjoint_face_quadrature_points=adjoint_face_quadrature_points,
        adjoint_ray_offset=adjoint_ray_offset,
        adjoint_proposals=adjoint_proposals,
        adjoint_proposal_frames=adjoint_proposal_frames,
        bidirectional_options=bidirectional_options,
        transport_device=transport_device)

    beta_current = float(beta); rejected_steps = 0; history = []

    def assess(state_charge):
        base_level = phase_space_log2_samples
        maximum_level = (base_level if phase_space_max_log2_samples is None
                         else int(phase_space_max_log2_samples))
        level_sequence = (None,) if base_level is None else range(int(base_level), maximum_level + 1)
        estimator_converged = phase_space_max_log2_samples is None
        for level in level_sequence:
            evaluations = []
            positive_replicates = []; negative_replicates = []
            for replicate in range(int(phase_space_replicates)):
                arguments = dict(evaluate_arguments)
                arguments["seed"] = int(seed) + 104729 * replicate
                arguments["phase_space_log2_samples"] = level
                item = _evaluate_incident_current_3d(charge=state_charge, **arguments)
                evaluations.append(item)
                positive_replicates.append(
                    item["positive_node_current"][tuple(support_nodes.T)])
                negative_replicates.append(
                    item["negative_node_current"][tuple(support_nodes.T)])
            positive_replicates = np.stack(positive_replicates)
            negative_replicates = np.stack(negative_replicates)
            positive = positive_replicates.mean(axis=0)
            negative = negative_replicates.mean(axis=0)
            signed_replicates = positive_replicates - negative_replicates
            net_stderr = (
                signed_replicates.std(axis=0, ddof=1) / np.sqrt(int(phase_space_replicates))
                if int(phase_space_replicates) > 1 else np.zeros_like(positive))
            total = positive + negative
            scale = float(np.max(total)) if total.size else 0.0
            active = total > max(1e-15 * scale, 1e-300)
            if not np.any(active):
                raise RuntimeError("steady charging has no resolved incident current on its surface")
            uncertainty_width = np.zeros_like(total)
            uncertainty_width[active] = (
                float(current_confidence_sigma) * net_stderr[active] / total[active])
            maximum_uncertainty = float(np.max(uncertainty_width[active]))
            if (phase_space_max_log2_samples is None
                    or maximum_uncertainty <= float(current_estimator_relative_tol)):
                estimator_converged = True
                break
        relative = np.zeros_like(total)
        relative[active] = np.abs(positive[active] - negative[active]) / total[active]
        confidence_envelope = relative + uncertainty_width
        current_floor = max(1e-15 * scale, 1e-300)
        log_ratio = np.log(
            np.maximum(positive, current_floor) / np.maximum(negative, current_floor))
        merit = float(np.sqrt(np.mean(relative[active] ** 2)))
        maximum = float(np.max(relative[active]))
        maximum_confidence = float(np.max(confidence_envelope[active]))
        return (evaluations[0], positive, negative, net_stderr, active, log_ratio,
                merit, maximum, maximum_confidence, maximum_uncertainty,
                estimator_converged, level)

    (evaluated, positive, negative, net_stderr, active, log_ratio,
     merit, maximum, maximum_confidence, maximum_uncertainty,
     estimator_converged, estimator_level) = assess(charge)
    history.append(dict(
        iteration=1, rms_relative_current_imbalance=merit,
        max_relative_current_imbalance=maximum, beta=beta_current,
        confidence_envelope_max_relative_current_imbalance=maximum_confidence,
        current_estimator_max_relative_uncertainty=maximum_uncertainty,
        current_estimator_converged=bool(estimator_converged),
        phase_space_log2_samples=(-1 if estimator_level is None else int(estimator_level)),
        mean_surface_voltage_v=float(np.mean(
            evaluated["potential"][tuple(support_nodes.T)]))))
    if evaluated["bidirectional_method_hint"]:
        frozen_options = ({} if bidirectional_options is None
                          else dict(bidirectional_options))
        discovered_hint = evaluated["bidirectional_method_hint"]
        if "method_hint" in frozen_options:
            supplied_hint = dict(frozen_options["method_hint"])
            if (set(supplied_hint) != set(discovered_hint)
                    or any(not np.array_equal(supplied_hint[name], discovered_hint[name])
                           for name in discovered_hint)):
                raise RuntimeError(
                    "supplied bidirectional method map differs from the certified initial map")
        else:
            frozen_options["method_hint"] = discovered_hint
        frozen_options["require_certification"] = False
        evaluate_arguments["bidirectional_options"] = frozen_options
    anderson_x_history = []; anderson_residual_history = []; cached_voltage_step = None
    while len(history) < int(max_iter) and not (
            len(history) >= int(min_iter)
            and estimator_converged
            and maximum_confidence <= float(current_balance_tol)):
        if cached_voltage_step is None:
            residual = float(response_energy_eV) * log_ratio
            residual[~active] = 0.0
            if nonlinear_update == "anderson":
                surface_voltage = evaluated["potential"][tuple(support_nodes.T)]
                cached_voltage_step = _anderson_step(
                    surface_voltage, residual, anderson_x_history,
                    anderson_residual_history, beta_current, int(anderson_depth))
            else:
                cached_voltage_step = beta_current * residual
        voltage_step = np.clip(
            cached_voltage_step,
            -float(maximum_voltage_step), float(maximum_voltage_step))
        voltage_step[~active] = 0.0
        trial_charge = charge.copy()
        # Invert the exact support-node Poisson response so the proposed charge increment produces
        # the requested voltage step despite strong electrostatic coupling between trench surfaces.
        # This is a preconditioner only: acceptance and convergence still use physical current balance.
        charge_step = np.linalg.solve(voltage_response, voltage_step)
        trial_charge[tuple(support_nodes.T)] += charge_step
        try:
            trial = assess(trial_charge)
        except BidirectionalCurrentCertificationError:
            beta_current *= 0.5; rejected_steps += 1
            if beta_current < float(minimum_beta):
                break
            cached_voltage_step *= 0.5
            continue
        trial_merit = trial[6]
        if trial_merit > merit * (1.0 + float(trust_growth_tolerance)):
            beta_current *= 0.5; rejected_steps += 1
            if beta_current < float(minimum_beta):
                break
            cached_voltage_step *= 0.5
            continue
        accepted_beta = beta_current
        if trial_merit < 0.8 * merit:
            beta_current = min(float(beta), 1.2 * beta_current)
        charge = trial_charge
        cached_voltage_step = None
        (evaluated, positive, negative, net_stderr, active, log_ratio,
         merit, maximum, maximum_confidence, maximum_uncertainty,
         estimator_converged, estimator_level) = trial
        history.append(dict(
            iteration=len(history) + 1, rms_relative_current_imbalance=merit,
            max_relative_current_imbalance=maximum, beta=accepted_beta,
            confidence_envelope_max_relative_current_imbalance=maximum_confidence,
            current_estimator_max_relative_uncertainty=maximum_uncertainty,
            current_estimator_converged=bool(estimator_converged),
            phase_space_log2_samples=(-1 if estimator_level is None else int(estimator_level)),
            mean_surface_voltage_v=float(np.mean(
                evaluated["potential"][tuple(support_nodes.T)]))))

    converged = bool(
        len(history) >= int(min_iter)
        and estimator_converged
        and maximum_confidence <= float(current_balance_tol))
    positive_grid = np.zeros(poisson_system.shape)
    negative_grid = np.zeros(poisson_system.shape)
    net_stderr_grid = np.zeros(poisson_system.shape)
    positive_grid[tuple(support_nodes.T)] = positive
    negative_grid[tuple(support_nodes.T)] = negative
    net_stderr_grid[tuple(support_nodes.T)] = net_stderr
    result = SteadyDielectricCharging3DResult(
        charge_node_c=charge, potential_v=evaluated["potential"],
        positive_current_node_a=positive_grid, negative_current_node_a=negative_grid,
        net_current_stderr_node_a=net_stderr_grid,
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
            f"confidence-envelope max relative imbalance={maximum_confidence:.6g}", result)
    return result
