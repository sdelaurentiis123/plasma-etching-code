"""Physical-time surface-charge state and saturation gates for 3-D co-evolution.

Face sheet charge is authoritative. Every accepted physical or pseudo-time step updates that sheet
charge with the exact kinetic face current and projects it through the same compatible Q1 operator
used by Poisson. The independently updated nodal charge must agree with that projection to roundoff.
This avoids an ill-posed nodal-to-face inverse when a moving surface is remeshed.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from time import perf_counter
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import PlasmaBoundaryState
from .boundary_transport_3d import (
    merge_boundary_transport_results_3d,
    trace_boundary_state_field_3d,
)
from .charging_coupled_3d import (
    DielectricChargingStep3DResult,
    _freeze_certified_bidirectional_options,
    advance_dielectric_charging_3d,
    current_balance_metrics_3d,
)
from .charging_poisson_3d import NodalPoissonSystem3D, lump_triangle_sheet_charge_3d
from .feature_step_3d import (
    FeatureGeometry3D,
    FeatureStep3DResult,
    FeatureStepValidity,
    _face_material_ids,
    _surface_gas_normals,
    advance_feature_step_3d,
)
from .surface_charge_remap_3d import SurfaceChargeRemap3DResult, remap_surface_charge_3d
from .threed import extract_mesh_3d


@dataclass(frozen=True)
class PhysicalPatchBalance3D:
    """Current-balance result at one fixed physical patch scale.

    The first four imbalance fields retain the symmetric ``(Ji-Je)/(Ji+Je)`` diagnostic used by
    the historical node report. ``b2_*`` implements the signed contract's stricter
    ``abs(Ji-Je)/Ji`` patch statistic and is the only measure used as the B2 acceptance gate.
    """

    patch_scale_m: float
    group: np.ndarray
    rms_relative_imbalance: float
    maximum_relative_imbalance: float
    throughput_weighted_rms_relative_imbalance: float
    global_relative_imbalance: float
    b2_rms_ion_normalized_imbalance: float
    b2_maximum_ion_normalized_imbalance: float
    b2_global_ion_normalized_imbalance: float
    active_patch_count: int

    def __post_init__(self):
        group = np.asarray(self.group, dtype=int).copy()
        values = np.asarray([
            self.patch_scale_m, self.rms_relative_imbalance,
            self.maximum_relative_imbalance,
            self.throughput_weighted_rms_relative_imbalance,
            self.global_relative_imbalance,
            self.b2_rms_ion_normalized_imbalance,
            self.b2_maximum_ion_normalized_imbalance,
            self.b2_global_ion_normalized_imbalance], dtype=float)
        if (group.ndim != 1 or np.any(group < 0) or np.any(np.isnan(values))
                or np.any(values < 0.0) or not np.isfinite(self.patch_scale_m)
                or self.patch_scale_m <= 0.0 or int(self.active_patch_count) < 0):
            raise ValueError("invalid physical patch-balance result")
        group.setflags(write=False)
        object.__setattr__(self, "group", group)
        object.__setattr__(self, "active_patch_count", int(self.active_patch_count))


@dataclass(frozen=True)
class ExperimentalObservableTolerance3D:
    """B3 tolerance anchored to combined benchmark and digitization uncertainty."""

    observable: str
    tolerance: float
    benchmark_uncertainty_including_digitization: float
    feature_extent_m: float | None = None

    def __post_init__(self):
        values = np.asarray([
            self.tolerance, self.benchmark_uncertainty_including_digitization], dtype=float)
        if (not self.observable or np.any(~np.isfinite(values)) or np.any(values < 0.0)
                or self.tolerance > self.benchmark_uncertainty_including_digitization
                or (self.feature_extent_m is not None
                    and (not np.isfinite(self.feature_extent_m)
                         or self.feature_extent_m <= 0.0))):
            raise ValueError(
                "an experimental-claim observable tolerance must not exceed its combined "
                "experimental and digitization uncertainty, and feature extent must be positive")


@dataclass(frozen=True)
class ResolvedBiasSegment3D:
    """One explicitly resolved boundary state and physical duration in a pulsed waveform."""

    duration_s: float
    boundary: PlasmaBoundaryState

    def __post_init__(self):
        if (not np.isfinite(self.duration_s) or self.duration_s <= 0.0
                or not isinstance(self.boundary, PlasmaBoundaryState)):
            raise ValueError("a resolved bias segment needs positive duration and a boundary state")


@dataclass(frozen=True)
class SurfaceChargingSaturation3DResult:
    """Saturated or explicitly exhausted face-charge trajectory on one fixed geometry."""

    sigma_c_per_m2: np.ndarray
    face_charge_c: np.ndarray
    charge_node_c: np.ndarray
    potential_v: np.ndarray
    final_step: DielectricChargingStep3DResult
    patch_balance: tuple[PhysicalPatchBalance3D, ...]
    history: tuple[Mapping[str, object], ...]
    converged: bool
    accepted_steps: int
    rejected_steps: int
    physical_time_s: float
    pseudo_time_s: float
    timestep_policy: str
    diagnostics: Mapping[str, object]

    def __post_init__(self):
        sigma = np.asarray(self.sigma_c_per_m2, dtype=float).copy()
        face_charge = np.asarray(self.face_charge_c, dtype=float).copy()
        node = np.asarray(self.charge_node_c, dtype=float).copy()
        potential = np.asarray(self.potential_v, dtype=float).copy()
        if (sigma.ndim != 1 or face_charge.shape != sigma.shape
                or node.shape != potential.shape
                or np.any(~np.isfinite(sigma)) or np.any(~np.isfinite(face_charge))
                or np.any(~np.isfinite(node)) or np.any(~np.isfinite(potential))
                or not isinstance(self.final_step, DielectricChargingStep3DResult)
                or len(self.patch_balance) < 2
                or self.timestep_policy not in {"fixed", "ser"}
                or int(self.accepted_steps) < 0 or int(self.rejected_steps) < 0
                or not np.isfinite(self.physical_time_s) or self.physical_time_s < 0.0
                or not np.isfinite(self.pseudo_time_s) or self.pseudo_time_s < 0.0):
            raise ValueError("invalid surface-charging saturation result")
        for value in (sigma, face_charge, node, potential):
            value.setflags(write=False)
        object.__setattr__(self, "sigma_c_per_m2", sigma)
        object.__setattr__(self, "face_charge_c", face_charge)
        object.__setattr__(self, "charge_node_c", node)
        object.__setattr__(self, "potential_v", potential)
        object.__setattr__(self, "patch_balance", tuple(self.patch_balance))
        object.__setattr__(
            self, "history", tuple(MappingProxyType(dict(item)) for item in self.history))
        object.__setattr__(self, "converged", bool(self.converged))
        object.__setattr__(self, "accepted_steps", int(self.accepted_steps))
        object.__setattr__(self, "rejected_steps", int(self.rejected_steps))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


class SurfaceChargingSaturationError(RuntimeError):
    """A C3 fixed-geometry charging trajectory failed with replayable state/history."""

    def __init__(
            self, message, sigma_c_per_m2, history, accepted_steps, rejected_steps,
            physical_time_s=0.0, pseudo_time_s=0.0):
        super().__init__(message)
        sigma = np.asarray(sigma_c_per_m2, dtype=float).copy()
        clocks = np.asarray([physical_time_s, pseudo_time_s], dtype=float)
        if np.any(~np.isfinite(clocks)) or np.any(clocks < 0.0):
            raise ValueError("failure checkpoint clocks must be finite and nonnegative")
        sigma.setflags(write=False)
        self.sigma_c_per_m2 = sigma
        self.history = tuple(MappingProxyType(dict(item)) for item in history)
        self.accepted_steps = int(accepted_steps)
        self.rejected_steps = int(rejected_steps)
        self.physical_time_s = float(physical_time_s)
        self.pseudo_time_s = float(pseudo_time_s)


def physical_surface_patch_groups_3d(
        face_centroids, face_gas_normals, face_material_id, patch_scale_m, *,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0), patch_origin_m=None):
    """Assign fixed physical surface patches independent of triangle/grid indices.

    Patches are physical Cartesian boxes anchored at ``patch_origin_m`` and split by material and
    dominant oriented surface-normal class. The latter prevents a wall and floor sharing one spatial
    box from being silently merged into a single balance equation.
    """
    centroid = np.asarray(face_centroids, dtype=float)
    normal = np.asarray(face_gas_normals, dtype=float)
    material = np.asarray(face_material_id, dtype=int)
    origin = np.asarray(mesh_origin_m, dtype=float)
    patch_origin = origin if patch_origin_m is None else np.asarray(patch_origin_m, dtype=float)
    if (centroid.ndim != 2 or centroid.shape[1] != 3
            or normal.shape != centroid.shape or material.shape != (len(centroid),)
            or len(centroid) == 0 or np.any(~np.isfinite(centroid))
            or np.any(~np.isfinite(normal)) or np.any(material <= 0)
            or not np.allclose(np.linalg.norm(normal, axis=1), 1.0, rtol=0.0, atol=2e-6)
            or not np.isfinite(patch_scale_m) or patch_scale_m <= 0.0
            or origin.shape != (3,) or patch_origin.shape != (3,)
            or np.any(~np.isfinite(origin)) or np.any(~np.isfinite(patch_origin))
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0):
        raise ValueError("invalid physical surface-patch inputs")
    physical = origin + centroid * float(mesh_length_unit_m)
    cell = np.floor((physical - patch_origin) / float(patch_scale_m) + 1e-12).astype(int)
    dominant_axis = np.argmax(np.abs(normal), axis=1)
    dominant_sign = (normal[np.arange(len(normal)), dominant_axis] >= 0.0).astype(int)
    key = np.column_stack((material, dominant_axis, dominant_sign, cell))
    _, group = np.unique(key, axis=0, return_inverse=True)
    return group.astype(int)


def _patch_balances(
        positive_face_current_density_a_m2, negative_face_current_density_a_m2,
        physical_face_area_m2, patch_groups, patch_scales_m):
    positive = np.asarray(positive_face_current_density_a_m2) * physical_face_area_m2
    negative = np.asarray(negative_face_current_density_a_m2) * physical_face_area_m2
    results = []
    for scale, group in zip(patch_scales_m, patch_groups):
        metrics = current_balance_metrics_3d(positive, negative, group=group)
        active = metrics.active
        ion_normalized = np.divide(
            np.abs(metrics.positive_current_a - metrics.negative_current_a),
            metrics.positive_current_a,
            out=np.full(metrics.positive_current_a.shape, np.inf),
            where=metrics.positive_current_a > 0.0)
        if np.any(active):
            b2_rms = float(np.sqrt(np.mean(ion_normalized[active] ** 2)))
            b2_maximum = float(np.max(ion_normalized[active]))
        else:
            b2_rms = b2_maximum = float("inf")
        total_positive = float(np.sum(metrics.positive_current_a))
        b2_global = (float(abs(
            np.sum(metrics.positive_current_a) - np.sum(metrics.negative_current_a))
            / total_positive) if total_positive > 0.0 else float("inf"))
        results.append(PhysicalPatchBalance3D(
            float(scale), group, metrics.rms_relative_imbalance,
            metrics.maximum_relative_imbalance,
            metrics.throughput_weighted_rms_relative_imbalance,
            metrics.global_relative_imbalance, b2_rms, b2_maximum, b2_global,
            metrics.active_count))
    return tuple(results)


def _ser_candidate_acceptance(
        activated_ser, residual_norm_a, last_residual_norm_a, allowed_residual_growth):
    """Safeguard PTC with the dimensional ODE residual, never a normalized gate ratio."""
    if not activated_ser or last_residual_norm_a is None:
        return True, None
    if residual_norm_a > (1.0 + allowed_residual_growth) * last_residual_norm_a:
        return False, "absolute_current_residual_growth"
    return True, None


def integrate_surface_charging_to_saturation_3d(
        poisson_system: NodalPoissonSystem3D, initial_sigma_c_per_m2,
        boundary: PlasmaBoundaryState, verts, faces, areas, *,
        face_centroids, face_gas_normals, face_material_id,
        source_bounds, source_z, potential_origin, potential_spacing,
        patch_scales_m, potential_rate_tolerance_v_s,
        timestep_s, maximum_steps, current_balance_tolerance=0.08,
        timestep_policy="fixed", maximum_timestep_s=None, minimum_timestep_s=None,
        ser_activation_rms=0.5, ser_maximum_growth=2.0,
        ser_allowed_residual_growth=0.005,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0),
        n_position=256, seed=0, trajectory_fixed_dt=0.01,
        trajectory_max_steps=10000, phase_space_log2_samples=None,
        periodic_lateral=False, transport_estimator="forward",
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-5,
        adjoint_proposals=None, adjoint_proposal_frames="surface_local",
        bidirectional_options=None, transport_device=None,
        charged_surface_response=None, surface_material_state=None,
        response_launch_offset=1e-5, response_fixed_dt=None,
        response_max_bounces=16, response_relative_tail_tolerance=0.0,
        stop_on_saturation=True, scramble_mode="frozen", sampling_seed_stride=1000003,
        fresh_adjoint_proposal_factory=None):
    """Integrate one fixed geometry to B1/B2 saturation with fixed time or safeguarded SER.

    SER follows the residual-ratio rule ``dt[n+1] = dt[n] * ||F[n]||/||F[n+1]||`` with declared
    minimum/maximum steps and an absolute-current-residual rejection safeguard. The signed B2 patch
    ratio remains an acceptance gate for the final state, but cannot reject a dynamical step because
    its local ion denominator can change non-monotonically. SER changes only the explicit step size
    of the same conservative charge ODE. The returned current and all gates are evaluated on the
    exact caller-supplied kinetic surface operator. ``fresh`` scrambling is restricted to fixed
    physical time: each accepted update receives a reproducible independent seed epoch, while the
    final state receives the next epoch for an honest current diagnostic. If adjoint proposals are
    supplied, callers must also supply a factory that regenerates every proposal from that epoch's
    seed. Fresh-scramble SER is deliberately refused because stochastic residual changes cannot
    safely drive its accept/reject controller.
    """
    if not isinstance(poisson_system, NodalPoissonSystem3D):
        raise TypeError("poisson_system must be NodalPoissonSystem3D")
    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=int)
    areas = np.asarray(areas, dtype=float)
    centroid = np.asarray(face_centroids, dtype=float)
    normal = np.asarray(face_gas_normals, dtype=float)
    material = np.asarray(face_material_id, dtype=int)
    sigma = np.asarray(initial_sigma_c_per_m2, dtype=float).copy()
    scales = tuple(float(value) for value in patch_scales_m)
    if (faces.ndim != 2 or faces.shape[1] != 3 or areas.shape != (len(faces),)
            or centroid.shape != (len(faces), 3) or normal.shape != centroid.shape
            or material.shape != (len(faces),) or sigma.shape != (len(faces),)
            or np.any(~np.isfinite(sigma)) or np.any(areas <= 0.0)
            or len(scales) < 2 or len(set(scales)) != len(scales)
            or any(not np.isfinite(value) or value <= 0.0 for value in scales)
            or not np.isfinite(potential_rate_tolerance_v_s)
            or potential_rate_tolerance_v_s <= 0.0
            or not np.isfinite(current_balance_tolerance)
            or current_balance_tolerance <= 0.0
            or not np.isfinite(timestep_s) or timestep_s <= 0.0
            or int(maximum_steps) != maximum_steps or maximum_steps < 0
            or timestep_policy not in {"fixed", "ser"}
            or scramble_mode not in {"frozen", "fresh"}
            or int(sampling_seed_stride) != sampling_seed_stride
            or sampling_seed_stride <= 0
            or not np.isfinite(response_relative_tail_tolerance)
            or not 0.0 <= response_relative_tail_tolerance < 1.0
            or not isinstance(stop_on_saturation, (bool, np.bool_))):
        raise ValueError("invalid C3 surface-charging integration inputs")
    if ((fresh_adjoint_proposal_factory is not None
         and not callable(fresh_adjoint_proposal_factory))
            or (scramble_mode == "frozen" and fresh_adjoint_proposal_factory is not None)
            or (scramble_mode == "fresh" and timestep_policy != "fixed")
            or (scramble_mode == "fresh" and adjoint_proposals is not None
                and fresh_adjoint_proposal_factory is None)):
        raise ValueError("invalid fresh-scramble proposal or timestep controls")
    if timestep_policy == "ser":
        maximum_timestep_s = (float(timestep_s) if maximum_timestep_s is None
                              else float(maximum_timestep_s))
        minimum_timestep_s = (float(timestep_s) / 1024.0 if minimum_timestep_s is None
                              else float(minimum_timestep_s))
        if (not 0.0 < minimum_timestep_s <= timestep_s <= maximum_timestep_s
                or not np.isfinite(ser_activation_rms) or ser_activation_rms <= 0.0
                or not np.isfinite(ser_maximum_growth) or ser_maximum_growth < 1.0
                or not np.isfinite(ser_allowed_residual_growth)
                or ser_allowed_residual_growth < 0.0):
            raise ValueError("invalid safeguarded SER controls")
    else:
        maximum_timestep_s = minimum_timestep_s = float(timestep_s)

    patch_groups = tuple(
        physical_surface_patch_groups_3d(
            centroid, normal, material, scale,
            mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m)
        for scale in scales)
    physical_area = areas * float(mesh_length_unit_m) ** 2
    projection = dict(
        shape=poisson_system.shape, vertices=verts, faces=faces,
        grid_origin=potential_origin, grid_spacing=potential_spacing,
        coordinate_length_unit_m=mesh_length_unit_m)
    charge = lump_triangle_sheet_charge_3d(sigma_c_per_m2=sigma, **projection)
    if np.any(np.abs(charge[poisson_system.dirichlet_mask]) > 0.0):
        raise ValueError("initial surface charge projects onto a Dirichlet reservoir")
    dt = float(timestep_s)
    history = []
    accepted = 0
    rejected = 0
    physical_time = 0.0
    pseudo_time = 0.0
    last_residual_norm = None
    pending_trial = None

    common = dict(
        poisson_system=poisson_system, boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=source_bounds, source_z=source_z, potential_origin=potential_origin,
        potential_spacing=potential_spacing, mesh_length_unit_m=mesh_length_unit_m,
        mesh_origin_m=mesh_origin_m, n_position=n_position, seed=seed,
        trajectory_fixed_dt=trajectory_fixed_dt, trajectory_max_steps=trajectory_max_steps,
        phase_space_log2_samples=phase_space_log2_samples,
        periodic_lateral=periodic_lateral, transport_estimator=transport_estimator,
        face_centroids=centroid, face_gas_normals=normal,
        adjoint_face_quadrature_points=adjoint_face_quadrature_points,
        adjoint_ray_offset=adjoint_ray_offset, adjoint_proposals=adjoint_proposals,
        adjoint_proposal_frames=adjoint_proposal_frames,
        bidirectional_options=bidirectional_options, transport_device=transport_device,
        charged_surface_response=charged_surface_response, face_material_id=material,
        surface_material_state=surface_material_state,
        response_launch_offset=response_launch_offset, response_fixed_dt=response_fixed_dt,
        response_max_bounces=response_max_bounces,
        response_relative_tail_tolerance=response_relative_tail_tolerance)

    final_step = None
    final_patch = None
    converged = False
    attempt = 0
    while True:
        attempt += 1
        sampling_epoch = int(
            0 if scramble_mode == "frozen"
            else accepted + (1 if pending_trial is not None else 0))
        sampling_seed = int(seed) + int(sampling_seed_stride) * sampling_epoch
        common["seed"] = sampling_seed
        if fresh_adjoint_proposal_factory is not None:
            common["adjoint_proposals"] = fresh_adjoint_proposal_factory(sampling_seed)
        try:
            step = advance_dielectric_charging_3d(
                charge_node_c=charge, duration_s=dt, **common)
        except Exception as error:
            raise SurfaceChargingSaturationError(
                f"C3 charging evaluation failed after {accepted} accepted steps: {error}",
                sigma, history, accepted, rejected, physical_time, pseudo_time) from error
        if step.bidirectional_method_hint:
            # Freeze the separately certified method map at its measured sample levels. This is the
            # same replay contract used by the lower physical-time engine and prevents estimator
            # reselection from using the samples that subsequently score an iteration.
            common["bidirectional_options"] = _freeze_certified_bidirectional_options(
                common["bidirectional_options"], step.bidirectional_method_hint,
                step.bidirectional_sampling_provenance)

        node_metrics = current_balance_metrics_3d(
            step.positive_current_node_a, step.negative_current_node_a)
        patch = _patch_balances(
            step.positive_face_current_density_a_m2,
            step.negative_face_current_density_a_m2,
            physical_area, patch_groups, scales)
        potential_rate = float(np.max(np.abs(
            step.potential_after_v - step.potential_before_v)) / dt)
        residual_norm = float(np.linalg.norm(
            step.positive_current_node_a - step.negative_current_node_a))
        charge_conservation_scale = max(
            float(step.diagnostics["absolute_incident_charge_c"]),
            abs(float(step.diagnostics["deposited_charge_c"])), np.finfo(float).tiny)
        charge_conservation_relative_error = abs(float(
            step.diagnostics["charge_conservation_residual_c"])) / charge_conservation_scale
        patch_maximum = max(
            item.b2_maximum_ion_normalized_imbalance for item in patch)
        merit = node_metrics.rms_relative_imbalance
        activated_ser = bool(timestep_policy == "ser" and merit <= ser_activation_rms)
        accept, rejection_reason = _ser_candidate_acceptance(
            activated_ser and pending_trial is not None,
            residual_norm, last_residual_norm,
            float(ser_allowed_residual_growth))
        prospective_accepted_steps = int(
            accepted + (1 if accept and pending_trial is not None else 0))
        item = dict(
            evaluation=int(attempt), accepted=bool(accept),
            scramble_mode=scramble_mode, sampling_epoch=sampling_epoch,
            sampling_seed=sampling_seed,
            accepted_steps=prospective_accepted_steps,
            rejected_steps=int(rejected), timestep_s=float(dt),
            physical_time_s=float(physical_time), pseudo_time_s=float(pseudo_time),
            ser_activated=activated_ser, rejection_reason=rejection_reason,
            potential_rate_max_v_s=potential_rate,
            rms_relative_current_imbalance_node=merit,
            max_relative_current_imbalance_node=node_metrics.maximum_relative_imbalance,
            rms_relative_current_imbalance_face=current_balance_metrics_3d(
                step.positive_face_current_density_a_m2 * physical_area,
                step.negative_face_current_density_a_m2 * physical_area
            ).rms_relative_imbalance,
            max_relative_current_imbalance_face=current_balance_metrics_3d(
                step.positive_face_current_density_a_m2 * physical_area,
                step.negative_face_current_density_a_m2 * physical_area
            ).maximum_relative_imbalance,
            residual_current_norm_a=residual_norm,
            maximum_patch_relative_imbalance=patch_maximum,
            patch_scales_m=scales,
            patch_rms_relative_imbalance=tuple(
                value.b2_rms_ion_normalized_imbalance for value in patch),
            patch_max_relative_imbalance=tuple(
                value.b2_maximum_ion_normalized_imbalance for value in patch),
            patch_symmetric_rms_relative_imbalance=tuple(
                value.rms_relative_imbalance for value in patch),
            patch_symmetric_max_relative_imbalance=tuple(
                value.maximum_relative_imbalance for value in patch),
            incident_charge_c=float(step.diagnostics["incident_charge_c"]),
            positive_incident_charge_c=float(
                step.diagnostics["positive_incident_charge_c"]),
            negative_incident_charge_c=float(
                step.diagnostics["negative_incident_charge_c"]),
            absolute_incident_charge_c=float(
                step.diagnostics["absolute_incident_charge_c"]),
            deposited_charge_c=float(step.diagnostics["deposited_charge_c"]),
            charge_conservation_residual_c=float(
                step.diagnostics["charge_conservation_residual_c"]),
            charge_conservation_relative_error=charge_conservation_relative_error,
            response_tail_closure_relative_absolute_charge_rate=float(getattr(
                step.surface_transfer,
                "tail_closure_relative_absolute_charge_rate", 0.0)),
            response_tail_closure_l1_current_error_bound_relative=float(getattr(
                step.surface_transfer,
                "tail_closure_l1_current_error_bound_relative", 0.0)),
            transport_lineage_replay_count=int(
                step.diagnostics["transport_lineage_replay_count"]),
            transport_lineage_replay_eligible_count=int(
                step.diagnostics["transport_lineage_replay_eligible_count"]),
            transport_lineage_replay_fraction=float(
                step.diagnostics["transport_lineage_replay_fraction"]),
            surface_transfer_relative_charge_balance_error=(
                step.surface_transfer.relative_charge_balance_error))
        history.append(item)
        if (charge_conservation_relative_error > 5e-13
                or step.surface_transfer.relative_charge_balance_error > 5e-13):
            raise SurfaceChargingSaturationError(
                "C3 charge-deposition ledger failed roundoff conservation",
                sigma, history, accepted, rejected, physical_time, pseudo_time)

        if not accept:
            rejected += 1
            sigma = pending_trial["sigma_c_per_m2"]
            charge = pending_trial["charge_node_c"]
            physical_time = pending_trial["physical_time_s"]
            pseudo_time = pending_trial["pseudo_time_s"]
            dt = 0.5 * pending_trial["timestep_s"]
            pending_trial = None
            if dt < float(minimum_timestep_s):
                raise SurfaceChargingSaturationError(
                    "safeguarded SER exhausted its minimum timestep",
                    sigma, history, accepted, rejected, physical_time, pseudo_time)
            continue

        if pending_trial is not None:
            accepted += 1
            pending_trial = None
        final_step = step
        final_patch = patch
        b1 = potential_rate <= float(potential_rate_tolerance_v_s)
        b2 = all(value.b2_maximum_ion_normalized_imbalance <= current_balance_tolerance
                 for value in patch)
        if b1 and b2:
            converged = True
            if stop_on_saturation:
                break
        if accepted >= int(maximum_steps):
            break

        if activated_ser:
            if last_residual_norm is not None and residual_norm > 0.0:
                ratio = last_residual_norm / residual_norm
                growth = min(float(ser_maximum_growth), max(0.5, ratio))
                dt = float(np.clip(dt * growth, minimum_timestep_s, maximum_timestep_s))
            last_residual_norm = residual_norm
        else:
            last_residual_norm = None
            dt = float(timestep_s)

        candidate_sigma = sigma + step.face_current_density_a_m2 * dt
        projected_candidate = lump_triangle_sheet_charge_3d(
            sigma_c_per_m2=candidate_sigma, **projection)
        nodal_candidate = charge + (
            step.positive_current_node_a - step.negative_current_node_a) * dt
        difference = projected_candidate - nodal_candidate
        projection_scale = max(
            float(np.sum(np.abs(projected_candidate))),
            float(np.sum(np.abs(nodal_candidate))), np.finfo(float).tiny)
        consistency = float(np.sum(np.abs(difference)) / projection_scale)
        history[-1]["face_to_node_update_relative_error"] = consistency
        if consistency > 5e-13:
            raise SurfaceChargingSaturationError(
                "face-charge and compatible-Q1 nodal updates diverged",
                sigma, history, accepted, rejected, physical_time, pseudo_time)
        pending_trial = dict(
            sigma_c_per_m2=sigma.copy(), charge_node_c=charge.copy(),
            physical_time_s=float(physical_time), pseudo_time_s=float(pseudo_time),
            timestep_s=float(dt))
        sigma = candidate_sigma
        charge = projected_candidate
        if activated_ser:
            pseudo_time += dt
        else:
            physical_time += dt

    if final_step is None or final_patch is None:  # pragma: no cover - guarded by first evaluation
        raise RuntimeError("surface charging produced no current evaluation")
    face_charge = sigma * physical_area
    return SurfaceChargingSaturation3DResult(
        sigma, face_charge, charge, final_step.potential_before_v,
        final_step, final_patch, tuple(history), converged, accepted, rejected,
        physical_time, pseudo_time, timestep_policy,
        diagnostics=dict(
            potential_rate_tolerance_v_s=float(potential_rate_tolerance_v_s),
            current_balance_tolerance=float(current_balance_tolerance),
            scramble_mode=scramble_mode,
            sampling_seed_stride=int(sampling_seed_stride),
            patch_scales_m=scales,
            exact_operator_statement=(
                "caller-supplied hard-visibility kinetic response; no smoothed residual; "
                f"response-tail tolerance={float(response_relative_tail_tolerance):.17g}"),
            response_relative_tail_tolerance=float(response_relative_tail_tolerance),
            final_response_tail_closure_l1_current_error_bound_relative=float(
                history[-1]["response_tail_closure_l1_current_error_bound_relative"]),
            retained_node_rms_relative_current_imbalance=(
                history[-1]["rms_relative_current_imbalance_node"]),
            retained_node_max_relative_current_imbalance=(
                history[-1]["max_relative_current_imbalance_node"]),
            final_potential_rate_max_v_s=(
                history[-1]["potential_rate_max_v_s"]),
            final_maximum_patch_relative_imbalance=(
                history[-1]["maximum_patch_relative_imbalance"]),
            maximum_transport_lineage_replay_count=max(
                item["transport_lineage_replay_count"] for item in history),
            maximum_transport_lineage_replay_fraction=max(
                item["transport_lineage_replay_fraction"] for item in history)))


@dataclass(frozen=True)
class ChargingCoevolutionStep3DResult:
    """One saturated-geometry transport/profile/remap transaction."""

    charging: SurfaceChargingSaturation3DResult
    feature: FeatureStep3DResult
    charge_remap: SurfaceChargeRemap3DResult
    wall_clock_s: float
    diagnostics: Mapping[str, object]

    def __post_init__(self):
        if (not isinstance(self.charging, SurfaceChargingSaturation3DResult)
                or not isinstance(self.feature, FeatureStep3DResult)
                or not isinstance(self.charge_remap, SurfaceChargeRemap3DResult)
                or not np.isfinite(self.wall_clock_s) or self.wall_clock_s < 0.0):
            raise ValueError("invalid charging co-evolution step")
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


@dataclass(frozen=True)
class ChargingCoevolution3DResult:
    """Quasi-static C3 profile trajectory with remapped surface-charge memory."""

    geometry: FeatureGeometry3D
    surface_state: object
    surface_state_mesh_fingerprint: str
    sigma_c_per_m2: np.ndarray
    steps: tuple[ChargingCoevolutionStep3DResult, ...]
    duration_s: float
    validity: FeatureStepValidity
    run_manifest: Mapping[str, object]

    def __post_init__(self):
        sigma = np.asarray(self.sigma_c_per_m2, dtype=float).copy()
        steps = tuple(self.steps)
        if (not isinstance(self.geometry, FeatureGeometry3D)
                or sigma.ndim != 1 or np.any(~np.isfinite(sigma))
                or not steps
                or any(not isinstance(item, ChargingCoevolutionStep3DResult) for item in steps)
                or not np.isfinite(self.duration_s) or self.duration_s < 0.0
                or not isinstance(self.validity, FeatureStepValidity)):
            raise ValueError("invalid charging co-evolution result")
        sigma.setflags(write=False)
        object.__setattr__(self, "sigma_c_per_m2", sigma)
        object.__setattr__(self, "steps", steps)
        object.__setattr__(self, "run_manifest", MappingProxyType(dict(self.run_manifest)))


def _response_manifest(response):
    if response is None:
        return None
    provenance = getattr(response, "provenance", None)
    if provenance is None:
        raise ValueError(
            "charged surface response requires a provenance manifest with parameters and bounds")
    manifest = dict(provenance)
    parameters = dict(manifest.get("parameters", {}))
    bounds = dict(manifest.get("bounds", {}))
    evidence = dict(manifest.get("evidence", {}))
    if (not parameters or set(parameters) != set(bounds) or set(parameters) != set(evidence)
            or any(len(tuple(bounds[name])) != 2 for name in parameters)):
        raise ValueError(
            "every charged-response parameter must have a source and declared bounds")
    return MappingProxyType(manifest)


def _array_manifest(value):
    array = np.ascontiguousarray(np.asarray(value))
    return dict(
        shape=tuple(int(item) for item in array.shape), dtype=str(array.dtype),
        sha256=sha256(array.view(np.uint8)).hexdigest())


def _manifest_value(value, *, path="manifest"):
    """Return a JSON-compatible provenance value or refuse an opaque runtime object."""
    if value is None or isinstance(value, (str, bool, int, float)):
        if isinstance(value, float) and not np.isfinite(value):
            raise ValueError(f"{path} contains a non-finite value")
        return value
    if isinstance(value, np.generic):
        return _manifest_value(value.item(), path=path)
    if isinstance(value, np.ndarray):
        return _array_manifest(value)
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            if not isinstance(key, (str, int)):
                raise ValueError(f"{path} contains a non-serializable mapping key")
            output[str(key)] = _manifest_value(item, path=f"{path}.{key}")
        return output
    if isinstance(value, (tuple, list)):
        return [_manifest_value(item, path=f"{path}[{index}]")
                for index, item in enumerate(value)]
    provenance = getattr(value, "provenance", None)
    if provenance is not None:
        return dict(
            type=type(value).__name__,
            provenance=_manifest_value(provenance, path=f"{path}.provenance"))
    raise ValueError(
        f"{path} contains opaque {type(value).__name__}; provide a provenance-bearing input")


def _boundary_manifest(boundary):
    return dict(
        reference_plane_m=float(boundary.reference_plane_m),
        provenance=_manifest_value(boundary.provenance, path="boundary.provenance"),
        species=[dict(
            name=species.name, charge_number=int(species.charge_number),
            mass_amu=float(species.mass_amu), flux_m2_s=float(species.flux_m2_s),
            velocity_sqrt_eV=_array_manifest(species.velocity_sqrt_eV),
            weight=_array_manifest(species.weight),
            phase_rad=(None if species.phase_rad is None
                       else _array_manifest(species.phase_rad)),
            position_m=(None if species.position_m is None
                        else _array_manifest(species.position_m)),
            density_model=(None if species.density_model is None
                           else type(species.density_model).__name__),
            provenance=_manifest_value(
                species.provenance, path=f"boundary.species.{species.name}.provenance"))
            for species in boundary.species])


def _merge_final_neutral_transport(
        charged_transport, geometry, boundary, species_role, verts, faces, areas,
        face_gas_normals, potential_v, *,
        source_bounds, source_z, potential_origin, potential_spacing,
        n_position, seed, trajectory_fixed_dt, trajectory_max_steps,
        periodic_lateral, transport_device):
    neutral_species = tuple(
        species for species in boundary.species if species.charge_number == 0)
    if not neutral_species:
        return charged_transport
    neutral_boundary = PlasmaBoundaryState(
        neutral_species, boundary.reference_plane_m, provenance=boundary.provenance)
    neutral_role = {species.name: species_role[species.name] for species in neutral_species}
    neutral = trace_boundary_state_field_3d(
        neutral_boundary, neutral_role, verts, faces, areas,
        source_bounds=source_bounds, source_z=source_z,
        nodal_potential_v=potential_v, potential_origin=potential_origin,
        potential_spacing=potential_spacing,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        mesh_origin_m=geometry.mesh_origin_m,
        n_position=n_position, seed=seed, fixed_dt=trajectory_fixed_dt,
        max_steps=trajectory_max_steps, periodic_lateral=periodic_lateral,
        face_gas_normals=face_gas_normals, device=transport_device)
    return merge_boundary_transport_results_3d(charged_transport, neutral)


def solve_charging_coevolution_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        charging_system_builder, etchable_material_ids,
        duration_s, n_steps, source_bounds, source_z,
        potential_origin, potential_spacing, charging_options,
        charged_surface_response=None, initial_sigma_c_per_m2=None,
        n_position=256, seed=0, trajectory_fixed_dt=0.01,
        trajectory_max_steps=10000, periodic_lateral=False,
        neutral_radiosity_options=None, cfl_number=0.3, reinitialize=True,
        reinitialization_method="skfmm", transport_device=None,
        bias_mode="quasi_static", bias_waveform=None,
        experimental_claim=False, observable_tolerances=()):
    """Run the signed C3 charge/profile co-evolution path.

    At every geometry, surface charge reaches the signed B1/B2 stationary state through the physical
    charge ODE (optionally explicit safeguarded SER), and the final exact charged/re-impact transport
    is reused by surface chemistry. Profile motion is followed by conservative C1 charge remap and a
    rebuilt Poisson operator. In ``waveform_resolved`` mode each declared segment instead advances
    that same physical charge ODE once for the segment duration and co-evolves the profile from the
    exact endpoint transport; no saturation or cycle-averaged bias is assumed. ``quasi_static``
    refuses a supplied waveform. Experimental claims must carry B3 tolerances already anchored to
    the combined experimental and digitization uncertainty.
    """
    if not isinstance(geometry, FeatureGeometry3D):
        raise TypeError("geometry must be FeatureGeometry3D")
    if not isinstance(boundary, PlasmaBoundaryState):
        raise TypeError("boundary must be PlasmaBoundaryState")
    if not callable(charging_system_builder):
        raise TypeError("charging_system_builder must be callable")
    if int(n_steps) != n_steps or n_steps <= 0:
        raise ValueError("n_steps must be a positive integer")
    if not np.isfinite(duration_s) or duration_s < 0.0:
        raise ValueError("duration_s must be finite and nonnegative")
    if bias_mode not in {"quasi_static", "waveform_resolved"}:
        raise ValueError("bias_mode must be 'quasi_static' or 'waveform_resolved'")
    if bias_mode == "quasi_static" and bias_waveform is not None:
        raise ValueError(
            "quasi-static charging refuses pulsed bias; select waveform_resolved co-simulation")
    tolerance_contracts = tuple(observable_tolerances)
    if any(not isinstance(item, ExperimentalObservableTolerance3D)
           for item in tolerance_contracts):
        raise TypeError(
            "observable_tolerances must contain ExperimentalObservableTolerance3D values")
    if bool(experimental_claim) and not tolerance_contracts:
        raise ValueError(
            "an experimental claim requires at least one uncertainty-anchored B3 tolerance")
    if not bool(experimental_claim) and tolerance_contracts:
        raise ValueError("observable tolerances require experimental_claim=True")
    response_manifest = _response_manifest(charged_surface_response)
    role = dict(species_role)
    if set(role) != {species.name for species in boundary.species}:
        raise ValueError("species_role must classify every boundary species")
    if bias_mode == "waveform_resolved":
        if bias_waveform is None:
            raise ValueError("waveform_resolved mode requires explicit bias segments")
        waveform = tuple(bias_waveform)
        if (len(waveform) != int(n_steps)
                or any(not isinstance(item, ResolvedBiasSegment3D) for item in waveform)):
            raise ValueError("n_steps must equal the number of resolved bias segments")
        if not np.isclose(
                sum(item.duration_s for item in waveform), float(duration_s),
                rtol=2e-14, atol=np.finfo(float).tiny):
            raise ValueError("resolved bias-segment durations must sum to duration_s")
        names = tuple(species.name for species in boundary.species)
        charge_and_mass = {
            species.name: (species.charge_number, species.mass_amu)
            for species in boundary.species}
        for item in waveform:
            if (tuple(species.name for species in item.boundary.species) != names
                    or item.boundary.reference_plane_m != boundary.reference_plane_m
                    or any((species.charge_number, species.mass_amu)
                           != charge_and_mass[species.name]
                           for species in item.boundary.species)):
                raise ValueError(
                    "resolved waveform segments must preserve species order, charge, mass, "
                    "and reference plane")
        step_boundaries = tuple(item.boundary for item in waveform)
        step_durations = tuple(item.duration_s for item in waveform)
    else:
        waveform = ()
        step_boundaries = (boundary,) * int(n_steps)
        step_durations = (float(duration_s) / int(n_steps),) * int(n_steps)
    options = dict(charging_options)
    if len({item.observable for item in tolerance_contracts}) != len(tolerance_contracts):
        raise ValueError("experimental observable contracts must have unique names")
    if any(item.feature_extent_m is not None for item in tolerance_contracts):
        scales = tuple(float(value) for value in options.get("patch_scales_m", ()))
        if not scales:
            raise ValueError("feature claims require declared physical patch scales")
        failed = [item.observable for item in tolerance_contracts
                  if item.feature_extent_m is not None
                  and min(scales) > item.feature_extent_m]
        if failed:
            raise ValueError(
                "at least one current-balance patch scale must not exceed the claimed feature "
                "extent for: " + ", ".join(failed))
    forbidden = {
        "poisson_system", "initial_sigma_c_per_m2", "boundary", "verts", "faces", "areas",
        "face_centroids", "face_gas_normals", "face_material_id", "source_bounds", "source_z",
        "potential_origin", "potential_spacing", "mesh_length_unit_m", "mesh_origin_m",
        "n_position", "seed", "trajectory_fixed_dt", "trajectory_max_steps",
        "periodic_lateral", "transport_device", "charged_surface_response",
        "stop_on_saturation"}
    overlap = set(options) & forbidden
    if overlap:
        raise ValueError(
            "charging_options duplicates driver-owned inputs: " + ", ".join(sorted(overlap)))

    current_geometry = geometry
    current_state = None
    current_fingerprint = None
    current_sigma = None if initial_sigma_c_per_m2 is None else np.asarray(
        initial_sigma_c_per_m2, dtype=float).copy()
    results = []
    for step_index, (step_boundary, step_duration) in enumerate(zip(
            step_boundaries, step_durations)):
        wall_start = perf_counter()
        verts, faces, centroids, areas = extract_mesh_3d(
            current_geometry.phi, current_geometry.dx)
        material = _face_material_ids(centroids, current_geometry)
        normals = _surface_gas_normals(verts, faces, centroids, current_geometry)
        if current_sigma is None:
            current_sigma = np.zeros(len(faces))
        if current_sigma.shape != (len(faces),):
            raise ValueError(
                f"co-evolution step {step_index + 1} charge state does not match surface mesh")
        poisson = charging_system_builder(current_geometry)
        if not isinstance(poisson, NodalPoissonSystem3D):
            raise TypeError("charging_system_builder must return NodalPoissonSystem3D")
        if poisson.shape != current_geometry.phi.shape:
            raise ValueError("rebuilt Poisson grid must match the current feature geometry")
        step_charging_options = dict(options)
        if bias_mode == "waveform_resolved":
            # One conservative physical update per explicitly resolved segment. Saturation fields
            # remain in the call only as reported diagnostics and are never interpreted as gates.
            step_charging_options.update(
                timestep_s=float(step_duration), maximum_steps=1, timestep_policy="fixed",
                stop_on_saturation=False)
        charging = integrate_surface_charging_to_saturation_3d(
            poisson, current_sigma, step_boundary, verts, faces, areas,
            face_centroids=centroids, face_gas_normals=normals,
            face_material_id=material, source_bounds=source_bounds, source_z=source_z,
            potential_origin=potential_origin, potential_spacing=potential_spacing,
            mesh_length_unit_m=current_geometry.mesh_length_unit_m,
            mesh_origin_m=current_geometry.mesh_origin_m,
            n_position=n_position, seed=int(seed) + step_index,
            trajectory_fixed_dt=trajectory_fixed_dt,
            trajectory_max_steps=trajectory_max_steps,
            periodic_lateral=periodic_lateral,
            transport_device=transport_device,
            charged_surface_response=charged_surface_response,
            surface_material_state=current_state,
            **step_charging_options)
        if bias_mode == "quasi_static" and not charging.converged:
            raise SurfaceChargingSaturationError(
                f"co-evolution step {step_index + 1} failed signed B1/B2 saturation gates",
                charging.sigma_c_per_m2, charging.history,
                charging.accepted_steps, charging.rejected_steps,
                charging.physical_time_s, charging.pseudo_time_s)
        transport = _merge_final_neutral_transport(
            charging.final_step.transport, current_geometry, step_boundary, role,
            verts, faces, areas, normals, charging.potential_v,
            source_bounds=source_bounds, source_z=source_z,
            potential_origin=potential_origin, potential_spacing=potential_spacing,
            n_position=n_position, seed=int(seed) + step_index,
            trajectory_fixed_dt=trajectory_fixed_dt,
            trajectory_max_steps=trajectory_max_steps,
            periodic_lateral=periodic_lateral, transport_device=transport_device)
        feature = advance_feature_step_3d(
            current_geometry, step_boundary, role, mechanism,
            etchable_material_ids=etchable_material_ids, duration_s=step_duration,
            source_bounds=source_bounds, source_z=source_z,
            surface_state=current_state,
            surface_state_mesh_fingerprint=current_fingerprint,
            n_position=n_position, seed=int(seed) + step_index,
            precomputed_transport=transport,
            neutral_radiosity_options=neutral_radiosity_options,
            cfl_number=cfl_number, reinitialize=reinitialize,
            reinitialization_method=reinitialization_method,
            transport_device=transport_device)
        next_verts, next_faces, next_centroids, _next_areas = extract_mesh_3d(
            feature.geometry.phi, feature.geometry.dx)
        next_material = _face_material_ids(next_centroids, feature.geometry)
        normal_displacement = -feature.face_velocity_mesh_units_s * step_duration
        remap = remap_surface_charge_3d(
            verts, faces, charging.sigma_c_per_m2, material,
            normal_displacement,
            next_verts, next_faces, next_material,
            mesh_length_unit_m=current_geometry.mesh_length_unit_m,
            maximum_distance=max(
                float(np.max(np.abs(normal_displacement))) + 1.5 * current_geometry.dx,
                np.finfo(float).tiny))
        wall_clock = perf_counter() - wall_start
        result = ChargingCoevolutionStep3DResult(
            charging, feature, remap, wall_clock,
            diagnostics=dict(
                step=step_index + 1, profile_duration_s=step_duration,
                bias_mode=bias_mode,
                saturation_required=bool(bias_mode == "quasi_static"),
                charging_saturated=charging.converged,
                charging_accepted_steps=charging.accepted_steps,
                charging_rejected_steps=charging.rejected_steps,
                charging_physical_time_s=charging.physical_time_s,
                charging_pseudo_time_s=charging.pseudo_time_s,
                retained_node_rms_relative_current_imbalance=(
                    charging.diagnostics["retained_node_rms_relative_current_imbalance"]),
                retained_node_max_relative_current_imbalance=(
                    charging.diagnostics["retained_node_max_relative_current_imbalance"]),
                potential_rate_max_v_s=(
                    charging.diagnostics["final_potential_rate_max_v_s"]),
                patch_scales_m=charging.diagnostics["patch_scales_m"],
                patch_max_relative_imbalance=tuple(
                    item.b2_maximum_ion_normalized_imbalance
                    for item in charging.patch_balance),
                patch_symmetric_max_relative_imbalance=tuple(
                    item.maximum_relative_imbalance for item in charging.patch_balance),
                maximum_transport_lineage_replay_count=charging.diagnostics[
                    "maximum_transport_lineage_replay_count"],
                maximum_transport_lineage_replay_fraction=charging.diagnostics[
                    "maximum_transport_lineage_replay_fraction"],
                remap_relative_charge_balance_error=remap.relative_charge_balance_error,
                retained_positive_charge_c=remap.retained_positive_charge_c,
                retained_negative_charge_c=remap.retained_negative_charge_c,
                removed_positive_charge_c=remap.removed_positive_charge_c,
                removed_negative_charge_c=remap.removed_negative_charge_c,
                removed_charge_c=remap.removed_net_charge_c,
                charged_surface_response_manifest=response_manifest))
        results.append(result)
        current_geometry = feature.geometry
        current_state = feature.next_surface_state
        current_fingerprint = feature.next_surface_state_mesh_fingerprint
        current_sigma = remap.sigma_c_per_m2.copy()

    reasons = tuple(
        reason for result in results for reason in result.feature.validity.reasons)
    limitations = tuple(dict.fromkeys(
        limitation for result in results
        for limitation in result.feature.validity.known_limitations)) + ((
        "quasi-static charge/profile separation requires charge saturation faster than profile motion",
    ) if bias_mode == "quasi_static" else (
        "waveform-resolved C3 is smoke-tested co-simulation, not experimentally validated",
    )) + (
        "C3 refinement and independent final-audit evidence are external campaign gates",
    )
    nonpredictive = tuple(dict.fromkeys(
        name for result in results
        for name in result.feature.validity.nonpredictive_parameters))
    manifest = dict(
        mode=bias_mode, n_steps=int(n_steps), duration_s=float(duration_s),
        seed=int(seed), n_position=int(n_position),
        trajectory_fixed_dt=float(trajectory_fixed_dt),
        trajectory_max_steps=int(trajectory_max_steps),
        periodic_lateral=bool(periodic_lateral),
        charging_options=_manifest_value(options, path="charging_options"),
        boundary=_boundary_manifest(boundary),
        waveform=(None if bias_mode == "quasi_static" else [dict(
            duration_s=item.duration_s, boundary=_boundary_manifest(item.boundary))
            for item in waveform]),
        charged_surface_response=(None if response_manifest is None else
                                  _manifest_value(response_manifest)),
        experimental_claim=bool(experimental_claim),
        observable_tolerances=[dict(
            observable=item.observable, tolerance=item.tolerance,
            benchmark_uncertainty_including_digitization=(
                item.benchmark_uncertainty_including_digitization),
            feature_extent_m=item.feature_extent_m)
            for item in tolerance_contracts],
        exact_operator="hard visibility; physical-time charge ODE; compatible Q1 projection",
        convergence_contract_revision="CCA-2026-07-13-R2")
    return ChargingCoevolution3DResult(
        current_geometry, current_state, current_fingerprint, current_sigma,
        tuple(results), float(duration_s),
        FeatureStepValidity(
            not reasons, reasons, limitations, not nonpredictive, nonpredictive),
        manifest)
