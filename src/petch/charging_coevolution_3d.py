"""Physical-time surface-charge state and saturation gates for 3-D co-evolution.

The exact kinetic current is measured on surface triangles and conservatively coupled to the Q1
Poisson load.  ``compatible_q1_charge_state`` makes the field-resolved nodal load authoritative and
stores its unique minimum-density-norm face representative.  This prevents exact P0 face modes that
are invisible to Q1 Poisson from accumulating without electrostatic feedback.  The legacy raw-face
state remains available for historical replay while the compatibility campaign is audited.
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
    average_boundary_transport_results_3d,
    merge_boundary_transport_results_3d,
    trace_boundary_state_field_3d,
)
from .charging_coupled_3d import (
    DielectricChargingStep3DResult,
    _freeze_certified_bidirectional_options,
    _validate_periodic_topology_3d,
    advance_dielectric_charging_3d,
    current_balance_metrics_3d,
)
from .charging_poisson_3d import (
    CompatibleQ1SurfaceChargeProjector3D,
    NodalPoissonSystem3D,
    lump_mixed_surface_density_3d,
    lump_triangle_sheet_charge_3d,
)
from .charging_stationarity_3d import (
    ProfileChargingStationarity3DResult,
    ProfileChargingStationarityBlock3D,
    ProfileChargingStationarityContract3D,
    assess_profile_charging_stationarity_3d,
)
from .feature_step_3d import (
    FeatureGeometry3D,
    FeatureStep3DResult,
    FeatureStepValidity,
    _face_material_ids,
    _surface_gas_normals,
    advance_feature_step_3d,
)
from .surface_charge_remap_3d import SurfaceChargeRemap3DResult, remap_surface_charge_3d
from .surface_kinetics import FaceResolvedEnergeticFlux
from .threed import extract_mesh_3d


CHARGING_RUN_MANIFEST_SCHEMA = "petch-charging-run-manifest-3d-v1"


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
    terminal_window_positive_face_current_density_a_m2: np.ndarray | None = None
    terminal_window_negative_face_current_density_a_m2: np.ndarray | None = None

    def __post_init__(self):
        sigma = np.asarray(self.sigma_c_per_m2, dtype=float).copy()
        face_charge = np.asarray(self.face_charge_c, dtype=float).copy()
        node = np.asarray(self.charge_node_c, dtype=float).copy()
        potential = np.asarray(self.potential_v, dtype=float).copy()
        window_positive = (
            None if self.terminal_window_positive_face_current_density_a_m2 is None
            else np.asarray(
                self.terminal_window_positive_face_current_density_a_m2, dtype=float).copy())
        window_negative = (
            None if self.terminal_window_negative_face_current_density_a_m2 is None
            else np.asarray(
                self.terminal_window_negative_face_current_density_a_m2, dtype=float).copy())
        if (sigma.ndim != 1 or face_charge.shape != sigma.shape
                or node.shape != potential.shape
                or np.any(~np.isfinite(sigma)) or np.any(~np.isfinite(face_charge))
                or np.any(~np.isfinite(node)) or np.any(~np.isfinite(potential))
                or not isinstance(self.final_step, DielectricChargingStep3DResult)
                or len(self.patch_balance) < 2
                or self.timestep_policy not in {"fixed", "ser", "decreasing_gain"}
                or int(self.accepted_steps) < 0 or int(self.rejected_steps) < 0
                or not np.isfinite(self.physical_time_s) or self.physical_time_s < 0.0
                or not np.isfinite(self.pseudo_time_s) or self.pseudo_time_s < 0.0):
            raise ValueError("invalid surface-charging saturation result")
        if ((window_positive is None) != (window_negative is None)
                or (window_positive is not None
                    and (window_positive.shape != sigma.shape
                         or window_negative.shape != sigma.shape
                         or np.any(~np.isfinite(window_positive))
                         or np.any(~np.isfinite(window_negative))
                         or np.any(window_positive < 0.0)
                         or np.any(window_negative < 0.0)))):
            raise ValueError("invalid terminal-window face currents")
        for value in (sigma, face_charge, node, potential):
            value.setflags(write=False)
        if window_positive is not None:
            window_positive.setflags(write=False)
            window_negative.setflags(write=False)
        object.__setattr__(self, "sigma_c_per_m2", sigma)
        object.__setattr__(self, "face_charge_c", face_charge)
        object.__setattr__(self, "charge_node_c", node)
        object.__setattr__(self, "potential_v", potential)
        object.__setattr__(
            self, "terminal_window_positive_face_current_density_a_m2", window_positive)
        object.__setattr__(
            self, "terminal_window_negative_face_current_density_a_m2", window_negative)
        object.__setattr__(self, "patch_balance", tuple(self.patch_balance))
        object.__setattr__(
            self, "history", tuple(MappingProxyType(dict(item)) for item in self.history))
        object.__setattr__(self, "converged", bool(self.converged))
        object.__setattr__(self, "accepted_steps", int(self.accepted_steps))
        object.__setattr__(self, "rejected_steps", int(self.rejected_steps))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


@dataclass(frozen=True)
class CompatibleQ1PseudoTimeProposal3D:
    """One globally compatible projective/PTC proposal awaiting an exact kinetic audit."""

    sigma_c_per_m2: np.ndarray
    face_charge_c: np.ndarray
    charge_node_c: np.ndarray
    potential_v: np.ndarray
    projected_face_current_a: np.ndarray
    pseudo_timestep_s: float
    diagnostics: Mapping[str, object]

    def __post_init__(self):
        sigma = np.asarray(self.sigma_c_per_m2, dtype=float).copy()
        face = np.asarray(self.face_charge_c, dtype=float).copy()
        node = np.asarray(self.charge_node_c, dtype=float).copy()
        potential = np.asarray(self.potential_v, dtype=float).copy()
        current = np.asarray(self.projected_face_current_a, dtype=float).copy()
        if (sigma.ndim != 1 or face.shape != sigma.shape or current.shape != sigma.shape
                or node.shape != potential.shape or np.any(~np.isfinite(sigma))
                or np.any(~np.isfinite(face)) or np.any(~np.isfinite(current))
                or np.any(~np.isfinite(node)) or np.any(~np.isfinite(potential))
                or not np.isfinite(self.pseudo_timestep_s)
                or self.pseudo_timestep_s < 0.0):
            raise ValueError("invalid compatible-Q1 pseudo-time proposal")
        for value in (sigma, face, node, potential, current):
            value.setflags(write=False)
        object.__setattr__(self, "sigma_c_per_m2", sigma)
        object.__setattr__(self, "face_charge_c", face)
        object.__setattr__(self, "charge_node_c", node)
        object.__setattr__(self, "potential_v", potential)
        object.__setattr__(self, "projected_face_current_a", current)
        object.__setattr__(self, "pseudo_timestep_s", float(self.pseudo_timestep_s))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


class SurfaceChargingSaturationError(RuntimeError):
    """A C3 fixed-geometry charging trajectory failed with replayable state/history."""

    def __init__(
            self, message, sigma_c_per_m2, history, accepted_steps, rejected_steps,
            physical_time_s=0.0, pseudo_time_s=0.0, *, state_updates=None,
            resume_sampling_epoch=0, resume_stochastic_gain_age_steps=0):
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
        self.state_updates = int(accepted_steps if state_updates is None else state_updates)
        self.resume_sampling_epoch = int(resume_sampling_epoch)
        self.resume_stochastic_gain_age_steps = int(
            resume_stochastic_gain_age_steps)
        if (self.state_updates < self.accepted_steps
                or self.resume_sampling_epoch < 0
                or self.resume_stochastic_gain_age_steps < 0):
            raise ValueError("failure checkpoint restart metadata must be nonnegative and ordered")


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


def propose_compatible_q1_pseudo_time_step_3d(
        poisson_system, vertices, faces, areas, sigma_c_per_m2,
        mean_face_current_density_a_m2, pseudo_timestep_s, *,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=1.0,
        mesh_length_unit_m=1.0, maximum_potential_jump_v=5.0):
    """Project a terminal-window current forward in pseudo-time without changing its zeros.

    This is the coarse projective/PTC *proposal* layer around the exact kinetic micro-integrator.
    It removes field-invisible P0 face modes, advances the compatible Q1 load with the measured
    terminal-window mean current, and refuses a jump larger than ``maximum_potential_jump_v``.
    The result is never a convergence claim: callers must score candidates with fresh exact
    hard-visibility samples and confirm accepted endpoints using fixed physical time.
    """
    if not isinstance(poisson_system, NodalPoissonSystem3D):
        raise TypeError("poisson_system must be NodalPoissonSystem3D")
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    areas = np.asarray(areas, dtype=float)
    sigma = np.asarray(sigma_c_per_m2, dtype=float)
    current_density = np.asarray(mean_face_current_density_a_m2, dtype=float)
    if (faces.ndim != 2 or faces.shape[1] != 3 or areas.shape != (len(faces),)
            or sigma.shape != (len(faces),) or current_density.shape != sigma.shape
            or np.any(~np.isfinite(sigma)) or np.any(~np.isfinite(current_density))
            or np.any(~np.isfinite(areas)) or np.any(areas <= 0.0)
            or not np.isfinite(pseudo_timestep_s) or pseudo_timestep_s < 0.0
            or not np.isfinite(maximum_potential_jump_v)
            or maximum_potential_jump_v <= 0.0):
        raise ValueError("invalid compatible-Q1 pseudo-time inputs")
    projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        poisson_system, vertices, faces,
        grid_origin=potential_origin, grid_spacing=potential_spacing,
        coordinate_length_unit_m=mesh_length_unit_m)
    physical_area = areas * float(mesh_length_unit_m) ** 2
    area_scale = max(float(np.max(physical_area)), np.finfo(float).tiny)
    if np.max(np.abs(projector.physical_face_area_m2 - physical_area)) / area_scale > 2e-12:
        raise ValueError("declared face areas disagree with the Q1 charge coupling")

    initial_face_charge = sigma * physical_area
    compatible_initial = projector.project_face_charge(initial_face_charge)
    raw_face_current = current_density * physical_area
    projected_face_current = projector.project_face_charge(raw_face_current)
    candidate_face_charge = projector.project_face_charge(
        compatible_initial + projected_face_current * float(pseudo_timestep_s))
    initial_reduced_node = projector.node_charge_from_face_charge(compatible_initial)
    candidate_reduced_node = projector.node_charge_from_face_charge(candidate_face_charge)
    initial_node = poisson_system.canonicalize_reduced_charge(initial_reduced_node)
    candidate_node = poisson_system.canonicalize_reduced_charge(candidate_reduced_node)
    initial_potential, _ = poisson_system.solve(initial_node)
    candidate_potential, _ = poisson_system.solve(candidate_node)
    maximum_jump = float(np.max(np.abs(candidate_potential - initial_potential)))
    if maximum_jump > float(maximum_potential_jump_v):
        raise ValueError(
            f"pseudo-time proposal potential jump {maximum_jump:g} V exceeds "
            f"the declared {float(maximum_potential_jump_v):g} V safeguard")

    expected_reduced_node = initial_reduced_node + projector.node_charge_from_face_charge(
        raw_face_current) * float(pseudo_timestep_s)
    node_scale = max(
        float(np.sum(np.abs(candidate_reduced_node))), np.finfo(float).tiny)
    node_consistency = float(np.sum(np.abs(
        candidate_reduced_node - expected_reduced_node)) / node_scale)
    if node_consistency > 5e-13:
        raise RuntimeError("compatible pseudo-time proposal changed the Q1 current update")
    expected_global = float(
        compatible_initial.sum() + raw_face_current.sum() * float(pseudo_timestep_s))
    charge_scale = max(
        float(np.sum(np.abs(compatible_initial))),
        float(np.sum(np.abs(raw_face_current))) * float(pseudo_timestep_s),
        np.finfo(float).tiny)
    global_error = float(candidate_face_charge.sum() - expected_global)
    if abs(global_error) / charge_scale > 5e-13:
        raise RuntimeError("compatible pseudo-time proposal lost global charge")
    return CompatibleQ1PseudoTimeProposal3D(
        sigma_c_per_m2=candidate_face_charge / physical_area,
        face_charge_c=candidate_face_charge,
        charge_node_c=candidate_node,
        potential_v=candidate_potential,
        projected_face_current_a=projected_face_current,
        pseudo_timestep_s=float(pseudo_timestep_s),
        diagnostics=dict(
            algorithm="compatible-Q1 coarse projective/PTC proposal; exact audit required",
            q1_face_coupling_rank=projector.rank,
            q1_face_coupling_nullity=projector.nullity,
            q1_face_coupling_condition_number=projector.condition_number,
            initial_q1_invisible_fraction=projector.unresolved_fraction(initial_face_charge),
            current_q1_invisible_fraction=projector.unresolved_fraction(raw_face_current),
            maximum_potential_jump_v=maximum_jump,
            maximum_allowed_potential_jump_v=float(maximum_potential_jump_v),
            q1_node_update_relative_l1_error=node_consistency,
            global_charge_error_c=global_error,
            physical_time_advanced_s=0.0,
            exact_audit_required=True))


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


def _terminal_window_current_mean(samples, field):
    """Average nonnegative interval currents without a subtractive sliding sum.

    A terminal window stores ``N + 1`` state evaluations.  The first ``N`` current maps
    advance the state across the window and therefore define its integrated balance.  Reusing a
    running sum and subtracting expired sparse maps can leave tiny negative entries after many
    fresh-scramble updates; those entries are numerical cancellation, but they violate the
    nonnegative-current contract and can abort an otherwise valid unattended transient.  A direct
    reduction is cheap relative to kinetic transport and preserves nonnegativity by construction.
    """
    interval_samples = tuple(samples[:-1])
    if not interval_samples:
        raise ValueError("a terminal current window requires at least one interval")
    arrays = tuple(np.asarray(item[field], dtype=float) for item in interval_samples)
    shape = arrays[0].shape
    if (any(value.shape != shape for value in arrays)
            or any(np.any(~np.isfinite(value)) or np.any(value < 0.0) for value in arrays)):
        raise ValueError("terminal-window currents must be finite nonnegative matching arrays")
    return np.add.reduce(arrays) / float(len(arrays))


def _q1_patch_balance_diagnostics(
        positive_face_current_density_a_m2, negative_face_current_density_a_m2,
        physical_face_area_m2, patch_groups, patch_scales_m, resolved_face_net_current_a,
        functional_null_sensitivity_max):
    """Report field-compatible balance and the unresolved raw numerator at each patch scale.

    These are change-control diagnostics for CCA-R3.  They never replace the signed raw B2 gate
    inside this routine.  The ion denominator remains the exact raw kinetic ion current; only the
    net-current numerator is decomposed into Q1-visible and Q1-null components.
    """
    positive = np.asarray(positive_face_current_density_a_m2) * physical_face_area_m2
    negative = np.asarray(negative_face_current_density_a_m2) * physical_face_area_m2
    resolved = np.asarray(resolved_face_net_current_a, dtype=float)
    if (positive.shape != negative.shape or resolved.shape != positive.shape
            or len(patch_groups) != len(patch_scales_m)
            or len(functional_null_sensitivity_max) != len(patch_scales_m)):
        raise ValueError("invalid Q1 patch-balance diagnostic inputs")
    raw_net = positive - negative
    output = []
    for scale, group, sensitivity in zip(
            patch_scales_m, patch_groups, functional_null_sensitivity_max):
        metrics = current_balance_metrics_3d(positive, negative, group=group)
        group_count = len(metrics.positive_current_a)
        resolved_group = np.bincount(
            group, weights=resolved, minlength=group_count)[:group_count]
        raw_group = metrics.positive_current_a - metrics.negative_current_a
        unresolved_group = raw_group - resolved_group
        resolved_ratio = np.divide(
            np.abs(resolved_group), metrics.positive_current_a,
            out=np.full(group_count, np.inf), where=metrics.positive_current_a > 0.0)
        unresolved_ratio = np.divide(
            np.abs(unresolved_group), metrics.positive_current_a,
            out=np.full(group_count, np.inf), where=metrics.positive_current_a > 0.0)
        active = metrics.active
        if np.any(active):
            resolved_rms = float(np.sqrt(np.mean(resolved_ratio[active] ** 2)))
            resolved_maximum = float(np.max(resolved_ratio[active]))
            unresolved_rms = float(np.sqrt(np.mean(unresolved_ratio[active] ** 2)))
            unresolved_maximum = float(np.max(unresolved_ratio[active]))
        else:
            resolved_rms = resolved_maximum = float("inf")
            unresolved_rms = unresolved_maximum = float("inf")
        total_positive = float(np.sum(metrics.positive_current_a))
        output.append(dict(
            patch_scale_m=float(scale),
            q1_resolved_rms_ion_normalized_imbalance=resolved_rms,
            q1_resolved_maximum_ion_normalized_imbalance=resolved_maximum,
            q1_resolved_global_ion_normalized_imbalance=(
                float(abs(np.sum(resolved_group)) / total_positive)
                if total_positive > 0.0 else float("inf")),
            q1_unresolved_rms_ion_normalized_imbalance=unresolved_rms,
            q1_unresolved_maximum_ion_normalized_imbalance=unresolved_maximum,
            q1_unresolved_global_ion_normalized_imbalance=(
                float(abs(np.sum(unresolved_group)) / total_positive)
                if total_positive > 0.0 else float("inf")),
            q1_patch_functional_null_sensitivity_max=float(sensitivity)))
    unresolved_total = float(np.sum(raw_net - resolved))
    regrouped_unresolved_total = float(np.sum(unresolved_group))
    unresolved_scale = max(
        float(np.sum(np.abs(raw_net - resolved))), np.finfo(float).tiny)
    current_roundoff = 64.0 * np.finfo(float).eps * max(
        float(np.sum(np.abs(raw_net))), float(np.sum(np.abs(resolved))),
        np.finfo(float).tiny)
    if abs(unresolved_total - regrouped_unresolved_total) > max(
            5e-13 * unresolved_scale, current_roundoff):
        # The last patch scale still partitions every face, so its group sum is a global check.
        raise RuntimeError("Q1 patch decomposition lost net current")
    return tuple(output)


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
        stochastic_gain_exponent=0.75, stochastic_gain_offset_steps=16,
        initial_stochastic_gain_age_steps=0,
        ser_activation_rms=0.5, ser_maximum_growth=2.0,
        ser_allowed_residual_growth=0.005,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0),
        n_position=256, seed=0, trajectory_fixed_dt=0.01,
        trajectory_max_steps=10000, trajectory_adaptive_horizon=False,
        trajectory_emergency_max_steps=None, phase_space_log2_samples=None,
        periodic_lateral=False, transport_estimator="forward",
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-5,
        adjoint_proposals=None, adjoint_proposal_frames="surface_local",
        bidirectional_options=None, transport_device=None,
        charged_surface_response=None, surface_material_state=None,
        response_launch_offset=1e-5, response_fixed_dt=None,
        response_max_bounces=16, response_relative_tail_tolerance=0.0,
        response_adaptive_bounce_extension=False,
        response_emergency_max_bounces=None,
        conductor_terminal=None,
        terminal_window_s=None,
        stop_on_saturation=True, scramble_mode="frozen", sampling_seed_stride=1000003,
        initial_sampling_epoch=0, fresh_adjoint_proposal_factory=None,
        progress_callback=None, compatible_q1_charge_state=False,
        physical_arrival_statistics="mean_flux"):
    """Integrate one fixed geometry with physical time, safeguarded SER, or decreasing gain.

    SER follows the residual-ratio rule ``dt[n+1] = dt[n] * ||F[n]||/||F[n+1]||`` with declared
    minimum/maximum steps and an absolute-current-residual rejection safeguard. The signed B2 patch
    ratio remains an acceptance gate for the final state, but cannot reject a dynamical step because
    its local ion denominator can change non-monotonically. SER changes only the explicit step size
    of the same conservative charge ODE. The returned current and all gates are evaluated on the
    exact caller-supplied kinetic surface operator. ``fresh`` scrambling is restricted to fixed
    physical time: each accepted update receives a reproducible independent seed epoch, while the
    final state receives the next epoch for an honest current diagnostic. A resumed run must pass
    that final diagnostic epoch as ``initial_sampling_epoch``: it is reused for the first resumed
    update because it scored, but did not create, the checkpoint state. If adjoint proposals are
    supplied, callers must also supply a factory that regenerates every proposal from that epoch's
    seed. Fresh-scramble SER is deliberately refused because stochastic residual changes cannot
    safely drive its accept/reject controller.

    The decreasing-gain policy is the late-stage stochastic-approximation path for an already
    established fresh-scramble stationary cloud. It uses
    dt[k] = timestep_s * (offset / (offset + age[k]))**exponent on the same conservative charge
    ODE. Requiring 0.5 < exponent <= 1 gives the ideal infinite schedule divergent total gain and
    summable squared gain. This path records pseudo-time only, requires stop_on_saturation=False,
    and must be followed by an independent exact-operator audit and fixed-physical-time
    confirmation. Its checkpointed age prevents a restart from silently resetting the schedule.

    ``terminal_window_s`` activates the signed CCA-R2 terminal-window gate for fixed physical
    time. B1 is the exact endpoint potential drift divided by the declared window duration. B2
    first time-integrates positive and negative patch currents over that same window and only then
    forms the ion-normalized ratio. Instantaneous RMS, worst-node, B1, and B2 diagnostics remain
    separate; a noisy single scramble can neither earn nor veto window saturation.

    When ``trajectory_adaptive_horizon`` is enabled, an incomplete primary, adjoint, or
    surface-reimpact flight is replayed from its identical launch state and sample epoch at the
    same fixed timestep with a doubled work horizon. This is inline numerical recovery, not a
    changed physical operator. An explicit emergency horizon remains a hard integrity stop for
    genuinely trapped or otherwise non-closing trajectories.

    ``progress_callback``, when supplied, is invoked after every fully certified current
    evaluation with read-only views of the current accepted state and a copy of its diagnostic
    record. It is the durability hook for atomic heartbeats and periodic checkpoints; callback
    failure stops with a replayable :class:`SurfaceChargingSaturationError` rather than allowing
    an unattended run to continue without its declared persistence contract.

    ``compatible_q1_charge_state`` projects the initial and every accepted face-charge state onto
    the area-weighted minimum-norm representative of the identical Q1 nodal load.  It changes no
    potential, field, trajectory, resolved current, or global charge.  Exact face modes in the
    null space of the declared Poisson discretization are removed and reported explicitly because
    no self-consistent field at that resolution can respond to them.
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
            or timestep_policy not in {"fixed", "ser", "decreasing_gain"}
            or scramble_mode not in {"frozen", "fresh"}
            or physical_arrival_statistics not in {"mean_flux", "poisson"}
            or int(sampling_seed_stride) != sampling_seed_stride
            or sampling_seed_stride <= 0
            or int(initial_sampling_epoch) != initial_sampling_epoch
            or initial_sampling_epoch < 0
            or int(stochastic_gain_offset_steps) != stochastic_gain_offset_steps
            or stochastic_gain_offset_steps <= 0
            or int(initial_stochastic_gain_age_steps)
                != initial_stochastic_gain_age_steps
            or initial_stochastic_gain_age_steps < 0
            or not np.isfinite(stochastic_gain_exponent)
            or not 0.5 < stochastic_gain_exponent <= 1.0
            or (scramble_mode == "frozen" and initial_sampling_epoch != 0)
            or not np.isfinite(response_relative_tail_tolerance)
            or not 0.0 <= response_relative_tail_tolerance < 1.0
            or not isinstance(response_adaptive_bounce_extension, (bool, np.bool_))
            or (response_emergency_max_bounces is not None
                and (int(response_emergency_max_bounces) != response_emergency_max_bounces
                     or response_emergency_max_bounces <= 0))
            or (response_adaptive_bounce_extension
                and response_emergency_max_bounces is None)
            or not isinstance(trajectory_adaptive_horizon, (bool, np.bool_))
            or (trajectory_emergency_max_steps is not None
                and (int(trajectory_emergency_max_steps) != trajectory_emergency_max_steps
                     or trajectory_emergency_max_steps < trajectory_max_steps))
            or (trajectory_adaptive_horizon
                and trajectory_emergency_max_steps is None)
            or (terminal_window_s is not None
                and (not np.isfinite(terminal_window_s) or terminal_window_s <= 0.0))
            or not isinstance(stop_on_saturation, (bool, np.bool_))
            or not isinstance(compatible_q1_charge_state, (bool, np.bool_))
            or (progress_callback is not None and not callable(progress_callback))):
        raise ValueError("invalid C3 surface-charging integration inputs")
    if ((fresh_adjoint_proposal_factory is not None
         and not callable(fresh_adjoint_proposal_factory))
            or (scramble_mode == "frozen" and fresh_adjoint_proposal_factory is not None)
            or (scramble_mode == "fresh"
                and timestep_policy not in {"fixed", "decreasing_gain"})
            or (physical_arrival_statistics == "poisson"
                and (scramble_mode != "fresh" or stop_on_saturation))
            or (scramble_mode == "frozen" and timestep_policy == "decreasing_gain")
            or (terminal_window_s is not None and timestep_policy != "fixed")
            or (timestep_policy == "decreasing_gain" and stop_on_saturation)
            or (timestep_policy != "decreasing_gain"
                and initial_stochastic_gain_age_steps != 0)
            or (scramble_mode == "fresh" and adjoint_proposals is not None
                and fresh_adjoint_proposal_factory is None)):
        raise ValueError("invalid fresh-scramble proposal or timestep controls")
    _validate_periodic_topology_3d(poisson_system, periodic_lateral)
    terminal_window_steps = None
    if terminal_window_s is not None:
        terminal_window_ratio = float(terminal_window_s) / float(timestep_s)
        terminal_window_steps = int(round(terminal_window_ratio))
        if (terminal_window_steps <= 0
                or not np.isclose(
                    terminal_window_ratio, terminal_window_steps,
                    rtol=2e-13, atol=2e-13)):
            raise ValueError(
                "terminal_window_s must be an integer multiple of the fixed physical timestep")
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
    face_conductor_id = poisson_system.classify_surface_floating_conductors(
        centroid, normal, grid_origin=potential_origin,
        grid_spacing=potential_spacing)
    projection = dict(
        shape=poisson_system.shape, vertices=verts, faces=faces,
        grid_origin=potential_origin, grid_spacing=potential_spacing,
        coordinate_length_unit_m=mesh_length_unit_m)
    charge_projector = None
    initial_unresolved_face_charge_fraction = 0.0
    initial_unresolved_face_charge_l1_c = 0.0
    original_charge = (
        lump_mixed_surface_density_3d(
            poisson_system, verts, faces, sigma, face_conductor_id,
            grid_origin=potential_origin, grid_spacing=potential_spacing,
            coordinate_length_unit_m=mesh_length_unit_m)
        if poisson_system.has_floating_conductors else
        poisson_system.canonicalize_charge(
            lump_triangle_sheet_charge_3d(sigma_c_per_m2=sigma, **projection)))
    if compatible_q1_charge_state:
        charge_projector = (
            CompatibleQ1SurfaceChargeProjector3D.from_mixed_poisson_system(
                poisson_system, verts, faces, face_conductor_id,
                grid_origin=potential_origin, grid_spacing=potential_spacing,
                coordinate_length_unit_m=mesh_length_unit_m)
            if poisson_system.has_floating_conductors else
            CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
                poisson_system, verts, faces,
                grid_origin=potential_origin, grid_spacing=potential_spacing,
                coordinate_length_unit_m=mesh_length_unit_m))
        area_error = np.max(np.abs(
            charge_projector.physical_face_area_m2 - physical_area))
        area_scale = max(float(np.max(physical_area)), np.finfo(float).tiny)
        if area_error / area_scale > 2e-12:
            raise ValueError("declared face areas disagree with the Q1 charge coupling")
        original_face_charge = sigma * physical_area
        compatible_face_charge = charge_projector.project_face_charge(original_face_charge)
        initial_unresolved_face_charge_fraction = charge_projector.unresolved_fraction(
            original_face_charge)
        initial_unresolved_face_charge_l1_c = float(np.sum(np.abs(
            original_face_charge - compatible_face_charge)))
        sigma = compatible_face_charge / physical_area
        charge = (
            lump_mixed_surface_density_3d(
                poisson_system, verts, faces, sigma, face_conductor_id,
                grid_origin=potential_origin, grid_spacing=potential_spacing,
                coordinate_length_unit_m=mesh_length_unit_m)
            if poisson_system.has_floating_conductors else
            poisson_system.canonicalize_charge(
                lump_triangle_sheet_charge_3d(
                    sigma_c_per_m2=sigma, **projection)))
        compatibility_scale = max(
            float(np.sum(np.abs(original_charge))), np.finfo(float).tiny)
        if np.sum(np.abs(charge - original_charge)) / compatibility_scale > 5e-13:
            raise RuntimeError("compatible-Q1 projection changed the electrostatic nodal load")
    else:
        charge = original_charge
    patch_functional_null_sensitivity_max = tuple(
        (0.0 if charge_projector is None else max(
            charge_projector.unresolved_linear_functional_fraction(
                (groups == group).astype(float))
            for group in np.unique(groups)))
        for groups in patch_groups)
    if np.any(np.abs(charge[poisson_system.dirichlet_mask]) > 0.0):
        raise ValueError("initial surface charge projects onto a Dirichlet reservoir")
    def decreasing_gain(age):
        return float(timestep_s) * (
            float(stochastic_gain_offset_steps)
            / (float(stochastic_gain_offset_steps) + float(age))) ** float(
                stochastic_gain_exponent)

    def resume_gain_age(state_updates):
        return int(
            initial_stochastic_gain_age_steps + state_updates
            if timestep_policy == "decreasing_gain" else 0)

    dt = (
        decreasing_gain(initial_stochastic_gain_age_steps)
        if timestep_policy == "decreasing_gain" else float(timestep_s))
    history = []
    terminal_samples = []
    accepted = 0
    rejected = 0
    physical_time = 0.0
    pseudo_time = 0.0
    last_residual_norm = None
    pending_trial = None
    cumulative_unresolved_face_current_l1_c = 0.0

    common = dict(
        poisson_system=poisson_system, boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=source_bounds, source_z=source_z, potential_origin=potential_origin,
        potential_spacing=potential_spacing, mesh_length_unit_m=mesh_length_unit_m,
        mesh_origin_m=mesh_origin_m, n_position=n_position, seed=seed,
        trajectory_fixed_dt=trajectory_fixed_dt, trajectory_max_steps=trajectory_max_steps,
        trajectory_adaptive_horizon=trajectory_adaptive_horizon,
        trajectory_emergency_max_steps=trajectory_emergency_max_steps,
        phase_space_log2_samples=phase_space_log2_samples,
        periodic_lateral=periodic_lateral, transport_estimator=transport_estimator,
        face_centroids=centroid, face_gas_normals=normal,
        adjoint_face_quadrature_points=adjoint_face_quadrature_points,
        adjoint_ray_offset=adjoint_ray_offset, adjoint_proposals=adjoint_proposals,
        adjoint_proposal_frames=adjoint_proposal_frames,
        bidirectional_options=bidirectional_options, transport_device=transport_device,
        charged_surface_response=charged_surface_response, face_material_id=material,
        face_conductor_id=face_conductor_id,
        surface_material_state=surface_material_state,
        response_launch_offset=response_launch_offset, response_fixed_dt=response_fixed_dt,
        response_max_bounces=response_max_bounces,
        response_relative_tail_tolerance=response_relative_tail_tolerance,
        response_adaptive_bounce_extension=response_adaptive_bounce_extension,
        response_emergency_max_bounces=response_emergency_max_bounces,
        conductor_terminal=conductor_terminal,
        physical_arrival_statistics=physical_arrival_statistics)

    final_step = None
    final_patch = None
    final_window_positive_face_current = None
    final_window_negative_face_current = None
    converged = False
    attempt = 0
    while True:
        attempt += 1
        sampling_epoch = int(
            0 if scramble_mode == "frozen"
            else initial_sampling_epoch + accepted + (1 if pending_trial is not None else 0))
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
                sigma, history, accepted, rejected, physical_time, pseudo_time,
                state_updates=accepted + int(pending_trial is not None),
                resume_sampling_epoch=sampling_epoch,
                resume_stochastic_gain_age_steps=resume_gain_age(
                    accepted + int(pending_trial is not None))) from error
        if step.bidirectional_method_hint:
            # Freeze the separately certified method map at its measured sample levels. This is the
            # same replay contract used by the lower physical-time engine and prevents estimator
            # reselection from using the samples that subsequently score an iteration.
            common["bidirectional_options"] = _freeze_certified_bidirectional_options(
                common["bidirectional_options"], step.bidirectional_method_hint,
                step.bidirectional_sampling_provenance)

        positive_independent_node_current_a = poisson_system.reduce_charge(
            step.positive_current_node_a)
        negative_independent_node_current_a = poisson_system.reduce_charge(
            step.negative_current_node_a)
        node_metrics = current_balance_metrics_3d(
            positive_independent_node_current_a,
            negative_independent_node_current_a)
        patch = _patch_balances(
            step.positive_face_current_density_a_m2,
            step.negative_face_current_density_a_m2,
            physical_area, patch_groups, scales)
        potential_rate = float(np.max(np.abs(
            step.potential_after_v - step.potential_before_v)) / dt)
        residual_norm = float(np.linalg.norm(
            positive_independent_node_current_a - negative_independent_node_current_a))
        charge_conservation_scale = max(
            float(step.diagnostics["absolute_incident_charge_c"]),
            abs(float(step.diagnostics["deposited_charge_c"])), np.finfo(float).tiny)
        charge_conservation_relative_error = abs(float(
            step.diagnostics["charge_conservation_residual_c"])) / charge_conservation_scale
        patch_maximum = max(
            item.b2_maximum_ion_normalized_imbalance for item in patch)
        face_net_current_a = step.face_current_density_a_m2 * physical_area
        if charge_projector is None:
            resolved_face_net_current_a = face_net_current_a
            unresolved_face_current_fraction = 0.0
            unresolved_face_current_l1_a = 0.0
            unresolved_face_current_net_a = 0.0
        else:
            resolved_face_net_current_a = charge_projector.project_face_charge(
                face_net_current_a)
            unresolved_face_current_a = face_net_current_a - resolved_face_net_current_a
            unresolved_face_current_fraction = charge_projector.unresolved_fraction(
                face_net_current_a)
            unresolved_face_current_l1_a = float(np.sum(np.abs(unresolved_face_current_a)))
            unresolved_face_current_net_a = float(np.sum(unresolved_face_current_a))
        q1_patch = _q1_patch_balance_diagnostics(
            step.positive_face_current_density_a_m2,
            step.negative_face_current_density_a_m2,
            physical_area, patch_groups, scales, resolved_face_net_current_a,
            patch_functional_null_sensitivity_max)
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
            stochastic_gain_age_steps=int(
                initial_stochastic_gain_age_steps + prospective_accepted_steps),
            stochastic_gain_exponent=(
                float(stochastic_gain_exponent)
                if timestep_policy == "decreasing_gain" else None),
            stochastic_gain_offset_steps=(
                int(stochastic_gain_offset_steps)
                if timestep_policy == "decreasing_gain" else None),
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
            compatible_q1_charge_state=bool(compatible_q1_charge_state),
            unresolved_face_current_fraction=unresolved_face_current_fraction,
            unresolved_face_current_l1_a=unresolved_face_current_l1_a,
            unresolved_face_current_net_a=unresolved_face_current_net_a,
            unresolved_face_current_projection_l1_c=0.0,
            cumulative_unresolved_face_current_projection_l1_c=(
                cumulative_unresolved_face_current_l1_c),
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
            patch_q1_resolved_rms_ion_normalized_imbalance=tuple(
                value["q1_resolved_rms_ion_normalized_imbalance"]
                for value in q1_patch),
            patch_q1_resolved_max_ion_normalized_imbalance=tuple(
                value["q1_resolved_maximum_ion_normalized_imbalance"]
                for value in q1_patch),
            patch_q1_unresolved_rms_ion_normalized_imbalance=tuple(
                value["q1_unresolved_rms_ion_normalized_imbalance"]
                for value in q1_patch),
            patch_q1_unresolved_max_ion_normalized_imbalance=tuple(
                value["q1_unresolved_maximum_ion_normalized_imbalance"]
                for value in q1_patch),
            patch_q1_functional_null_sensitivity_max=tuple(
                patch_functional_null_sensitivity_max),
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
            response_initial_bounce_budget=int(
                step.diagnostics["response_initial_bounce_budget"]),
            response_final_bounce_budget=int(
                step.diagnostics["response_final_bounce_budget"]),
            response_emergency_bounce_limit=int(
                step.diagnostics["response_emergency_bounce_limit"]),
            response_bounce_budget_extension_count=int(
                step.diagnostics["response_bounce_budget_extension_count"]),
            response_derived_bounce_budget=int(
                step.diagnostics["response_derived_bounce_budget"]),
            transport_lineage_replay_count=int(
                step.diagnostics["transport_lineage_replay_count"]),
            transport_lineage_replay_eligible_count=int(
                step.diagnostics["transport_lineage_replay_eligible_count"]),
            transport_lineage_replay_fraction=float(
                step.diagnostics["transport_lineage_replay_fraction"]),
            transport_edge_launch_inset_count=int(
                step.diagnostics["transport_edge_launch_inset_count"]),
            transport_trajectory_horizon_extension_count=int(
                step.diagnostics["transport_trajectory_horizon_extension_count"]),
            transport_trajectory_initial_max_steps=int(
                step.diagnostics["transport_trajectory_initial_max_steps"]),
            transport_trajectory_final_max_steps=int(
                step.diagnostics["transport_trajectory_final_max_steps"]),
            transport_trajectory_emergency_max_steps=int(
                step.diagnostics["transport_trajectory_emergency_max_steps"]),
            physical_arrival_statistics=step.diagnostics.get(
                "physical_arrival_statistics", "mean_flux"),
            physical_arrival_expected_primary_count=step.diagnostics.get(
                "physical_arrival_expected_primary_count", 0.0),
            physical_arrival_realized_primary_count=step.diagnostics.get(
                "physical_arrival_realized_primary_count", 0),
            conductor_terminal_active=bool(
                step.diagnostics.get("conductor_terminal_active", False)),
            conductor_terminal_signed_current_a=float(
                step.diagnostics.get("conductor_terminal_signed_current_a", 0.0)),
            conductor_terminal_absolute_current_a=float(
                step.diagnostics.get("conductor_terminal_absolute_current_a", 0.0)),
            conductor_terminal_face_patch_b2_exclusion=bool(
                step.diagnostics.get(
                    "conductor_terminal_face_patch_b2_exclusion", False)),
            surface_transfer_relative_charge_balance_error=(
                step.surface_transfer.relative_charge_balance_error))
        history.append(item)
        if (charge_conservation_relative_error > 5e-13
                or step.surface_transfer.relative_charge_balance_error > 5e-13):
            raise SurfaceChargingSaturationError(
                "C3 charge-deposition ledger failed roundoff conservation",
                sigma, history, accepted, rejected, physical_time, pseudo_time,
                state_updates=accepted + int(pending_trial is not None),
                resume_sampling_epoch=sampling_epoch,
                resume_stochastic_gain_age_steps=resume_gain_age(
                    accepted + int(pending_trial is not None)))

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
                    sigma, history, accepted, rejected, physical_time, pseudo_time,
                    state_updates=accepted, resume_sampling_epoch=sampling_epoch,
                    resume_stochastic_gain_age_steps=resume_gain_age(accepted))
            continue

        if pending_trial is not None:
            accepted += 1
            pending_trial = None
        gate_patch = patch
        gate_q1_patch = q1_patch
        gate_potential_rate = potential_rate
        terminal_window_ready = False
        if terminal_window_steps is not None:
            sample = dict(
                potential_v=np.asarray(step.potential_before_v, dtype=float).copy(),
                positive_face_current_density_a_m2=np.asarray(
                    step.positive_face_current_density_a_m2, dtype=float).copy(),
                negative_face_current_density_a_m2=np.asarray(
                    step.negative_face_current_density_a_m2, dtype=float).copy())
            terminal_samples.append(sample)
            if len(terminal_samples) > terminal_window_steps + 1:
                terminal_samples.pop(0)
            terminal_window_ready = len(terminal_samples) == terminal_window_steps + 1
            if terminal_window_ready:
                gate_potential_rate = float(np.max(np.abs(
                    terminal_samples[-1]["potential_v"]
                    - terminal_samples[0]["potential_v"])) / float(terminal_window_s))
                positive_window_mean = _terminal_window_current_mean(
                    terminal_samples, "positive_face_current_density_a_m2")
                negative_window_mean = _terminal_window_current_mean(
                    terminal_samples, "negative_face_current_density_a_m2")
                final_window_positive_face_current = positive_window_mean
                final_window_negative_face_current = negative_window_mean
                gate_patch = _patch_balances(
                    positive_window_mean, negative_window_mean,
                    physical_area, patch_groups, scales)
                window_net_current_a = (
                    positive_window_mean - negative_window_mean) * physical_area
                resolved_window_net_current_a = (
                    window_net_current_a if charge_projector is None
                    else charge_projector.project_face_charge(window_net_current_a))
                gate_q1_patch = _q1_patch_balance_diagnostics(
                    positive_window_mean, negative_window_mean,
                    physical_area, patch_groups, scales, resolved_window_net_current_a,
                    patch_functional_null_sensitivity_max)
            else:
                gate_potential_rate = None
        gate_patch_maximum = max(
            value.b2_maximum_ion_normalized_imbalance for value in gate_patch)
        item.update(
            gate_evaluation_mode=(
                "terminal_window" if terminal_window_steps is not None else "instantaneous"),
            gate_potential_rate_max_v_s=gate_potential_rate,
            gate_maximum_patch_relative_imbalance=gate_patch_maximum,
            gate_patch_rms_relative_imbalance=tuple(
                value.b2_rms_ion_normalized_imbalance for value in gate_patch),
            gate_patch_max_relative_imbalance=tuple(
                value.b2_maximum_ion_normalized_imbalance for value in gate_patch),
            gate_patch_q1_resolved_rms_ion_normalized_imbalance=tuple(
                value["q1_resolved_rms_ion_normalized_imbalance"]
                for value in gate_q1_patch),
            gate_patch_q1_resolved_max_ion_normalized_imbalance=tuple(
                value["q1_resolved_maximum_ion_normalized_imbalance"]
                for value in gate_q1_patch),
            gate_patch_q1_unresolved_rms_ion_normalized_imbalance=tuple(
                value["q1_unresolved_rms_ion_normalized_imbalance"]
                for value in gate_q1_patch),
            gate_patch_q1_unresolved_max_ion_normalized_imbalance=tuple(
                value["q1_unresolved_maximum_ion_normalized_imbalance"]
                for value in gate_q1_patch),
            terminal_window_ready=bool(terminal_window_ready),
            terminal_window_duration_s=(
                None if terminal_window_steps is None else float(terminal_window_s)),
            terminal_window_state_count=(
                0 if terminal_window_steps is None else len(terminal_samples)),
            terminal_window_required_state_count=(
                0 if terminal_window_steps is None else terminal_window_steps + 1))
        final_step = step
        final_patch = gate_patch
        gate_ready = bool(
            terminal_window_steps is None or terminal_window_ready)
        b1 = bool(
            gate_ready and gate_potential_rate <= float(potential_rate_tolerance_v_s))
        b2 = bool(
            gate_ready and all(
                value.b2_maximum_ion_normalized_imbalance <= current_balance_tolerance
                for value in gate_patch))
        converged = bool(
            b1 and b2 and timestep_policy != "decreasing_gain")
        item["b1_potential_saturation_satisfied"] = bool(b1)
        item["b2_patch_balance_satisfied"] = bool(b2)
        item["saturation_gates_satisfied"] = converged
        if progress_callback is not None:
            sigma_view = sigma.view()
            charge_view = charge.view()
            potential_view = step.potential_before_v.view()
            sigma_view.setflags(write=False)
            charge_view.setflags(write=False)
            potential_view.setflags(write=False)
            try:
                progress_callback(
                    sigma_c_per_m2=sigma_view,
                    charge_node_c=charge_view,
                    potential_v=potential_view,
                    history_item=MappingProxyType(dict(item)),
                    accepted_steps=int(accepted), rejected_steps=int(rejected),
                    physical_time_s=float(physical_time),
                    pseudo_time_s=float(pseudo_time),
                    resume_sampling_epoch=int(sampling_epoch))
            except Exception as error:
                raise SurfaceChargingSaturationError(
                    f"C3 progress persistence failed after {accepted} accepted steps: {error}",
                    sigma, history, accepted, rejected, physical_time, pseudo_time,
                    state_updates=accepted,
                    resume_sampling_epoch=sampling_epoch,
                    resume_stochastic_gain_age_steps=resume_gain_age(accepted)) from error
        if converged:
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
        elif timestep_policy == "decreasing_gain":
            last_residual_norm = None
            dt = decreasing_gain(
                initial_stochastic_gain_age_steps + accepted)
        else:
            last_residual_norm = None
            dt = float(timestep_s)

        candidate_face_charge = sigma * physical_area + face_net_current_a * dt
        if charge_projector is not None:
            compatible_candidate_face_charge = charge_projector.project_face_charge(
                candidate_face_charge)
            removed_l1_c = float(np.sum(np.abs(
                candidate_face_charge - compatible_candidate_face_charge)))
            cumulative_unresolved_face_current_l1_c += removed_l1_c
            candidate_face_charge = compatible_candidate_face_charge
            history[-1]["unresolved_face_current_projection_l1_c"] = removed_l1_c
            history[-1]["cumulative_unresolved_face_current_projection_l1_c"] = (
                cumulative_unresolved_face_current_l1_c)
        else:
            history[-1]["unresolved_face_current_projection_l1_c"] = 0.0
            history[-1]["cumulative_unresolved_face_current_projection_l1_c"] = 0.0
        candidate_sigma = candidate_face_charge / physical_area
        projected_candidate = (
            lump_mixed_surface_density_3d(
                poisson_system, verts, faces, candidate_sigma, face_conductor_id,
                grid_origin=potential_origin, grid_spacing=potential_spacing,
                coordinate_length_unit_m=mesh_length_unit_m)
            if poisson_system.has_floating_conductors else
            poisson_system.canonicalize_charge(
                lump_triangle_sheet_charge_3d(
                    sigma_c_per_m2=candidate_sigma, **projection)))
        nodal_candidate = poisson_system.canonicalize_charge(
            charge + (step.positive_current_node_a - step.negative_current_node_a) * dt)
        difference = projected_candidate - nodal_candidate
        projection_scale = max(
            float(np.sum(np.abs(projected_candidate))),
            float(np.sum(np.abs(nodal_candidate))), np.finfo(float).tiny)
        consistency = float(np.sum(np.abs(difference)) / projection_scale)
        history[-1]["face_to_node_update_relative_error"] = consistency
        if consistency > 5e-13:
            raise SurfaceChargingSaturationError(
                "face-charge and compatible-Q1 nodal updates diverged",
                sigma, history, accepted, rejected, physical_time, pseudo_time,
                state_updates=accepted, resume_sampling_epoch=sampling_epoch,
                resume_stochastic_gain_age_steps=resume_gain_age(accepted))
        pending_trial = dict(
            sigma_c_per_m2=sigma.copy(), charge_node_c=charge.copy(),
            physical_time_s=float(physical_time), pseudo_time_s=float(pseudo_time),
            timestep_s=float(dt))
        sigma = candidate_sigma
        charge = projected_candidate
        if activated_ser or timestep_policy == "decreasing_gain":
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
            physical_arrival_statistics=physical_arrival_statistics,
            sampling_seed_stride=int(sampling_seed_stride),
            initial_sampling_epoch=int(initial_sampling_epoch),
            resume_sampling_epoch=int(history[-1]["sampling_epoch"]),
            stochastic_gain_exponent=(
                float(stochastic_gain_exponent)
                if timestep_policy == "decreasing_gain" else None),
            stochastic_gain_offset_steps=(
                int(stochastic_gain_offset_steps)
                if timestep_policy == "decreasing_gain" else None),
            initial_stochastic_gain_age_steps=int(
                initial_stochastic_gain_age_steps),
            resume_stochastic_gain_age_steps=int(
                initial_stochastic_gain_age_steps + accepted),
            compatible_q1_charge_state=bool(compatible_q1_charge_state),
            poisson_periodic_axes=tuple(poisson_system.periodic_axes),
            poisson_independent_node_shape=tuple(poisson_system.reduced_shape),
            q1_face_coupling_rank=(
                len(faces) if charge_projector is None else charge_projector.rank),
            q1_face_coupling_nullity=(
                0 if charge_projector is None else charge_projector.nullity),
            q1_face_coupling_condition_number=(
                1.0 if charge_projector is None else charge_projector.condition_number),
            initial_unresolved_face_charge_fraction=(
                initial_unresolved_face_charge_fraction),
            initial_unresolved_face_charge_l1_c=initial_unresolved_face_charge_l1_c,
            maximum_unresolved_face_current_fraction=max(
                item["unresolved_face_current_fraction"] for item in history),
            maximum_unresolved_face_current_l1_a=max(
                item["unresolved_face_current_l1_a"] for item in history),
            maximum_absolute_unresolved_face_current_net_a=max(
                abs(item["unresolved_face_current_net_a"]) for item in history),
            cumulative_unresolved_face_current_projection_l1_c=(
                cumulative_unresolved_face_current_l1_c),
            patch_scales_m=scales,
            gate_evaluation_mode=history[-1]["gate_evaluation_mode"],
            terminal_window_s=(
                None if terminal_window_s is None else float(terminal_window_s)),
            terminal_window_steps=(
                None if terminal_window_steps is None else int(terminal_window_steps)),
            terminal_window_ready=bool(history[-1]["terminal_window_ready"]),
            terminal_window_state_count=int(history[-1]["terminal_window_state_count"]),
            exact_operator_statement=(
                "caller-supplied hard-visibility kinetic response; no smoothed residual; "
                f"response-tail tolerance={float(response_relative_tail_tolerance):.17g}; "
                "terminal-window gates, when enabled, integrate exact physical-time currents "
                "without changing the kinetic operator"),
            response_relative_tail_tolerance=float(response_relative_tail_tolerance),
            response_adaptive_bounce_extension=bool(
                response_adaptive_bounce_extension),
            response_emergency_max_bounces=(
                None if response_emergency_max_bounces is None
                else int(response_emergency_max_bounces)),
            trajectory_adaptive_horizon=bool(trajectory_adaptive_horizon),
            trajectory_initial_max_steps=int(trajectory_max_steps),
            trajectory_emergency_max_steps=(
                None if trajectory_emergency_max_steps is None
                else int(trajectory_emergency_max_steps)),
            final_response_tail_closure_l1_current_error_bound_relative=float(
                history[-1]["response_tail_closure_l1_current_error_bound_relative"]),
            maximum_response_bounce_budget=max(
                item["response_final_bounce_budget"] for item in history),
            maximum_response_bounce_budget_extension_count=max(
                item["response_bounce_budget_extension_count"] for item in history),
            retained_node_rms_relative_current_imbalance=(
                history[-1]["rms_relative_current_imbalance_node"]),
            retained_node_max_relative_current_imbalance=(
                history[-1]["max_relative_current_imbalance_node"]),
            final_instantaneous_potential_rate_max_v_s=(
                history[-1]["potential_rate_max_v_s"]),
            final_potential_rate_max_v_s=(
                history[-1]["gate_potential_rate_max_v_s"]),
            final_instantaneous_maximum_patch_relative_imbalance=(
                history[-1]["maximum_patch_relative_imbalance"]),
            final_maximum_patch_relative_imbalance=(
                history[-1]["gate_maximum_patch_relative_imbalance"]),
            final_patch_q1_resolved_rms_ion_normalized_imbalance=(
                history[-1]["gate_patch_q1_resolved_rms_ion_normalized_imbalance"]),
            final_patch_q1_resolved_max_ion_normalized_imbalance=(
                history[-1]["gate_patch_q1_resolved_max_ion_normalized_imbalance"]),
            final_patch_q1_unresolved_rms_ion_normalized_imbalance=(
                history[-1]["gate_patch_q1_unresolved_rms_ion_normalized_imbalance"]),
            final_patch_q1_unresolved_max_ion_normalized_imbalance=(
                history[-1]["gate_patch_q1_unresolved_max_ion_normalized_imbalance"]),
            patch_q1_functional_null_sensitivity_max=(
                patch_functional_null_sensitivity_max),
            maximum_transport_lineage_replay_count=max(
                item["transport_lineage_replay_count"] for item in history),
            maximum_transport_lineage_replay_fraction=max(
                item["transport_lineage_replay_fraction"] for item in history),
            maximum_transport_edge_launch_inset_count=max(
                item["transport_edge_launch_inset_count"] for item in history),
            maximum_transport_trajectory_horizon_extension_count=max(
                item["transport_trajectory_horizon_extension_count"] for item in history),
            maximum_transport_trajectory_final_max_steps=max(
                item["transport_trajectory_final_max_steps"] for item in history)),
        terminal_window_positive_face_current_density_a_m2=(
            final_window_positive_face_current),
        terminal_window_negative_face_current_density_a_m2=(
            final_window_negative_face_current))


@dataclass(frozen=True)
class ChargingCoevolutionStep3DResult:
    """One saturated-geometry transport/profile/remap transaction."""

    charging: SurfaceChargingSaturation3DResult
    feature: FeatureStep3DResult
    charge_remap: SurfaceChargeRemap3DResult
    wall_clock_s: float
    diagnostics: Mapping[str, object]
    profile_stationarity: ProfileChargingStationarity3DResult | None = None

    def __post_init__(self):
        if (not isinstance(self.charging, SurfaceChargingSaturation3DResult)
                or not isinstance(self.feature, FeatureStep3DResult)
                or not isinstance(self.charge_remap, SurfaceChargeRemap3DResult)
                or (self.profile_stationarity is not None
                    and not isinstance(
                        self.profile_stationarity, ProfileChargingStationarity3DResult))
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


def _geometry_manifest(geometry):
    layers = geometry.material_levelsets
    return dict(
        phi=_array_manifest(geometry.phi),
        material_id=_array_manifest(geometry.material_id),
        dx=float(geometry.dx),
        mesh_length_unit_m=float(geometry.mesh_length_unit_m),
        mesh_origin_m=[float(value) for value in geometry.mesh_origin_m],
        material_levelsets=(None if layers is None else {
            str(material_id): _array_manifest(layers[material_id])
            for material_id in sorted(layers)}))


def _surface_mechanism_manifest(mechanism):
    """Expose the material law in the run record without fabricating provenance.

    Older mechanism classes predate the common-engine manifest contract.  They remain runnable, but
    their absence of a machine-readable provenance block is made explicit.  New composite mechanisms
    such as ``MaterialMechanismRouter3D`` carry their complete declared routing evidence here.
    """
    provenance = getattr(mechanism, "provenance", None)
    return dict(
        type=type(mechanism).__name__,
        provenance=(None if provenance is None else
                    _manifest_value(provenance, path="surface_mechanism.provenance")),
        machine_readable_provenance=bool(provenance is not None))


def _merge_final_neutral_transport(
        charged_transport, geometry, boundary, species_role, verts, faces, areas,
        face_gas_normals, potential_v, *,
        source_bounds, source_z, potential_origin, potential_spacing,
        n_position, seed, trajectory_fixed_dt, trajectory_max_steps,
        trajectory_adaptive_horizon, trajectory_emergency_max_steps,
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
        face_gas_normals=face_gas_normals, device=transport_device,
        adaptive_horizon=trajectory_adaptive_horizon,
        emergency_max_steps=trajectory_emergency_max_steps)
    return merge_boundary_transport_results_3d(charged_transport, neutral)


def _transport_face_flux_map_3d(transport, face_count):
    """Return every delivered species as one face-flux vector for stationarity scoring."""
    output = {
        name: np.asarray(value, dtype=float).copy()
        for name, value in transport.surface_fluxes.neutral_flux_m2_s.items()}
    for population in transport.surface_fluxes.energetic_fluxes:
        if not isinstance(population, FaceResolvedEnergeticFlux):
            raise TypeError("profile stationarity requires face-resolved energetic transport")
        if population.name in output:
            raise ValueError("one species cannot be both neutral and energetic")
        output[population.name] = np.asarray(population.flux_m2_s, dtype=float).copy()
    if (not output or any(value.shape != (face_count,) for value in output.values())
            or any(np.any(~np.isfinite(value)) or np.any(value < 0.0)
                   for value in output.values())):
        raise ValueError("stationarity transport must resolve every species on the fixed mesh")
    return output


def _replicate_mean_and_standard_error(values):
    stack = np.stack([np.asarray(value, dtype=float) for value in values])
    if len(stack) < 2:
        raise ValueError("independent stationarity scoring requires at least two replicates")
    return np.mean(stack, axis=0), np.std(stack, axis=0, ddof=1) / np.sqrt(len(stack))


def solve_charging_coevolution_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        charging_system_builder, etchable_material_ids,
        duration_s, n_steps, source_bounds, source_z,
        potential_origin, potential_spacing, charging_options,
        charged_surface_response=None, initial_sigma_c_per_m2=None,
        initial_surface_state=None, initial_surface_state_mesh_fingerprint=None,
        restart_source_manifest_sha256=None,
        n_position=256, seed=0, trajectory_fixed_dt=0.01,
        trajectory_max_steps=10000, trajectory_adaptive_horizon=False,
        trajectory_emergency_max_steps=None, periodic_lateral=False,
        neutral_radiosity_options=None, neutral_forward_scatter=None,
        neutral_forward_scatter_options=None,
        surface_product_redeposition_options=None,
        cfl_number=0.3, reinitialize=True,
        reinitialization_method="skfmm", transport_device=None,
        profile_motion_enabled=True,
        bias_mode="quasi_static", bias_waveform=None,
        experimental_claim=False, observable_tolerances=(),
        charging_acceptance="signed_r2", stationarity_contract=None,
        stationarity_block_steps=None, stationarity_scoring_replicates=2):
    """Run the signed C3 charge/profile co-evolution path.

    At every geometry, surface charge reaches the signed B1/B2 stationary state through the physical
    charge ODE (optionally explicit safeguarded SER), and the final exact charged/re-impact transport
    is reused by surface chemistry. Profile motion is followed by conservative C1 charge remap and a
    rebuilt Poisson operator. In ``waveform_resolved`` mode each declared segment instead advances
    that same physical charge ODE once for the segment duration and co-evolves the profile from the
    exact endpoint transport; no saturation or cycle-averaged bias is assumed. ``quasi_static``
    refuses a supplied waveform. Experimental claims must carry B3 tolerances already anchored to
    the combined experimental and digitization uncertainty.

    ``charging_acceptance='profile_stationary'`` is an explicit draft alternative for numerical
    product development. It runs two consecutive fresh-scramble physical-time blocks and scores
    both endpoints with disjoint independent replicate ensembles. Profile motion is accepted only
    when the second block no longer materially changes the potential, kinetic currents, delivered
    species fluxes, or predicted profile increment. This path preserves hard visibility and the
    exact kinetic operator, reports signed-R2 diagnostics forever, and cannot support an
    experimental claim until its separately versioned contract is signed.
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
    if ((initial_surface_state is None)
            != (initial_surface_state_mesh_fingerprint is None)
            or (initial_surface_state_mesh_fingerprint is not None
                and (not isinstance(initial_surface_state_mesh_fingerprint, str)
                     or not initial_surface_state_mesh_fingerprint))):
        raise ValueError(
            "initial surface state and its nonempty mesh fingerprint must be supplied together")
    if (restart_source_manifest_sha256 is not None
            and (not isinstance(restart_source_manifest_sha256, str)
                 or len(restart_source_manifest_sha256) != 64
                 or any(character not in "0123456789abcdef"
                        for character in restart_source_manifest_sha256))):
        raise ValueError("restart source manifest must be a lowercase SHA-256 digest")
    if (int(trajectory_max_steps) != trajectory_max_steps or trajectory_max_steps <= 0
            or not isinstance(trajectory_adaptive_horizon, (bool, np.bool_))
            or (trajectory_emergency_max_steps is not None
                and (int(trajectory_emergency_max_steps) != trajectory_emergency_max_steps
                     or trajectory_emergency_max_steps < trajectory_max_steps))
            or (trajectory_adaptive_horizon and trajectory_emergency_max_steps is None)):
        raise ValueError("invalid co-evolution trajectory-horizon controls")
    if not isinstance(profile_motion_enabled, (bool, np.bool_)):
        raise TypeError("profile_motion_enabled must be boolean")
    if bias_mode not in {"quasi_static", "waveform_resolved", "physical_time_resolved"}:
        raise ValueError(
            "bias_mode must be 'quasi_static', 'waveform_resolved', or "
            "'physical_time_resolved'")
    if bias_mode == "quasi_static" and bias_waveform is not None:
        raise ValueError(
            "quasi-static charging refuses pulsed bias; select waveform_resolved co-simulation")
    if bias_mode == "physical_time_resolved" and bias_waveform is not None:
        raise ValueError(
            "constant-boundary physical-time co-evolution does not accept a bias waveform")
    if bias_mode == "physical_time_resolved" and experimental_claim:
        raise ValueError(
            "physical-time stochastic co-evolution is an ensemble mode, not a single-run "
            "experimental claim")
    if not profile_motion_enabled and experimental_claim:
        raise ValueError(
            "a fixed-geometry diagnostic cannot support an experimental profile claim")
    if charging_acceptance not in {"signed_r2", "profile_stationary"}:
        raise ValueError("charging_acceptance must be 'signed_r2' or 'profile_stationary'")
    if charging_acceptance == "signed_r2":
        if stationarity_contract is not None or stationarity_block_steps is not None:
            raise ValueError(
                "profile-stationarity controls require charging_acceptance='profile_stationary'")
    else:
        if bias_mode != "quasi_static":
            raise ValueError("profile stationarity applies only to quasi-static charging")
        if not profile_motion_enabled:
            raise ValueError("profile stationarity requires profile motion to be enabled")
        if not isinstance(stationarity_contract, ProfileChargingStationarityContract3D):
            raise TypeError("profile stationarity requires an explicit stationarity_contract")
        if (stationarity_block_steps is None
                or int(stationarity_block_steps) != stationarity_block_steps
                or stationarity_block_steps <= 0
                or int(stationarity_scoring_replicates) != stationarity_scoring_replicates
                or stationarity_scoring_replicates
                < stationarity_contract.minimum_independent_replicates):
            raise ValueError("invalid stationarity block length or scoring replicate count")
        if experimental_claim:
            raise ValueError(
                "the draft profile-stationarity contract cannot authorize experimental claims")
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
    forward_scatter_manifest = _response_manifest(neutral_forward_scatter)
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
    elif bias_mode == "physical_time_resolved":
        waveform = ()
        step_boundaries = (boundary,) * int(n_steps)
        step_durations = (float(duration_s) / int(n_steps),) * int(n_steps)
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
        "trajectory_adaptive_horizon", "trajectory_emergency_max_steps",
        "periodic_lateral", "transport_device", "charged_surface_response",
        "neutral_forward_scatter", "neutral_forward_scatter_options",
        "surface_product_redeposition_options",
        "stop_on_saturation", "initial_surface_state",
        "initial_surface_state_mesh_fingerprint", "restart_source_manifest_sha256"}
    overlap = set(options) & forbidden
    if overlap:
        raise ValueError(
            "charging_options duplicates driver-owned inputs: " + ", ".join(sorted(overlap)))

    current_geometry = geometry
    current_state = initial_surface_state
    current_fingerprint = initial_surface_state_mesh_fingerprint
    current_sigma = None if initial_sigma_c_per_m2 is None else np.asarray(
        initial_sigma_c_per_m2, dtype=float).copy()
    results = []
    for step_index, (step_boundary, step_duration) in enumerate(zip(
            step_boundaries, step_durations)):
        profile_step_duration = (
            float(step_duration) if profile_motion_enabled else 0.0)
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
        face_conductor_id = poisson.classify_surface_floating_conductors(
            centroids, normals, grid_origin=potential_origin,
            grid_spacing=potential_spacing)
        def integrate_charging(sigma, call_options):
            return integrate_surface_charging_to_saturation_3d(
                poisson, sigma, step_boundary, verts, faces, areas,
                face_centroids=centroids, face_gas_normals=normals,
                face_material_id=material, source_bounds=source_bounds, source_z=source_z,
                potential_origin=potential_origin, potential_spacing=potential_spacing,
                mesh_length_unit_m=current_geometry.mesh_length_unit_m,
                mesh_origin_m=current_geometry.mesh_origin_m,
                n_position=n_position, seed=int(seed) + step_index,
                trajectory_fixed_dt=trajectory_fixed_dt,
                trajectory_max_steps=trajectory_max_steps,
                trajectory_adaptive_horizon=trajectory_adaptive_horizon,
                trajectory_emergency_max_steps=trajectory_emergency_max_steps,
                periodic_lateral=periodic_lateral,
                transport_device=transport_device,
                charged_surface_response=charged_surface_response,
                surface_material_state=current_state,
                **call_options)

        def complete_transport(charged_transport, potential, sampling_seed):
            return _merge_final_neutral_transport(
                charged_transport, current_geometry, step_boundary, role,
                verts, faces, areas, normals, potential,
                source_bounds=source_bounds, source_z=source_z,
                potential_origin=potential_origin, potential_spacing=potential_spacing,
                n_position=n_position, seed=int(sampling_seed),
                trajectory_fixed_dt=trajectory_fixed_dt,
                trajectory_max_steps=trajectory_max_steps,
                trajectory_adaptive_horizon=trajectory_adaptive_horizon,
                trajectory_emergency_max_steps=trajectory_emergency_max_steps,
                periodic_lateral=periodic_lateral, transport_device=transport_device)

        def advance_profile(precomputed_transport, sampling_seed):
            return advance_feature_step_3d(
                current_geometry, step_boundary, role, mechanism,
                etchable_material_ids=etchable_material_ids,
                duration_s=profile_step_duration,
                source_bounds=source_bounds, source_z=source_z,
                surface_state=current_state,
                surface_state_mesh_fingerprint=current_fingerprint,
                n_position=n_position, seed=int(sampling_seed),
                precomputed_transport=precomputed_transport,
                trajectory_max_steps=trajectory_max_steps,
                trajectory_adaptive_horizon=trajectory_adaptive_horizon,
                trajectory_emergency_max_steps=trajectory_emergency_max_steps,
                profile_periodic_lateral=periodic_lateral,
                neutral_radiosity_options=neutral_radiosity_options,
                neutral_forward_scatter=neutral_forward_scatter,
                neutral_forward_scatter_options=neutral_forward_scatter_options,
                surface_product_redeposition_options=(
                    surface_product_redeposition_options),
                cfl_number=cfl_number, reinitialize=reinitialize,
                reinitialization_method=reinitialization_method,
                transport_device=transport_device)

        profile_stationarity = None
        step_charging_options = dict(options)
        if step_charging_options.get("timestep_policy") == "decreasing_gain":
            # Decreasing gain is a bounded stochastic warm-start path. It is deliberately
            # incapable of self-certifying B1/B2, so the public driver owns and supplies the
            # non-stopping integration control rather than accepting it through charging_options.
            step_charging_options["stop_on_saturation"] = False
        if bias_mode in {"waveform_resolved", "physical_time_resolved"}:
            # One conservative physical update per explicitly resolved segment. Saturation fields
            # remain in the call only as reported diagnostics and are never interpreted as gates.
            step_charging_options.update(
                timestep_s=float(step_duration), maximum_steps=1, timestep_policy="fixed",
                stop_on_saturation=False)
            if step_charging_options.get("physical_arrival_statistics") == "poisson":
                step_charging_options["scramble_mode"] = "fresh"
        if charging_acceptance == "signed_r2":
            charging = integrate_charging(current_sigma, step_charging_options)
            if bias_mode == "quasi_static" and not charging.converged:
                raise SurfaceChargingSaturationError(
                    f"co-evolution step {step_index + 1} failed signed B1/B2 saturation gates",
                    charging.sigma_c_per_m2, charging.history,
                    charging.accepted_steps, charging.rejected_steps,
                    charging.physical_time_s, charging.pseudo_time_s,
                    state_updates=charging.accepted_steps,
                    resume_sampling_epoch=charging.diagnostics[
                        "resume_sampling_epoch"],
                    resume_stochastic_gain_age_steps=charging.diagnostics[
                        "resume_stochastic_gain_age_steps"])
            sampling_seed = charging.history[-1]["sampling_seed"]
            transport = complete_transport(
                charging.final_step.transport, charging.potential_v, sampling_seed)
            feature = advance_profile(transport, sampling_seed)
        else:
            # The profile-relevant path deliberately uses fixed physical time and fresh independent
            # scrambles.  It may not inherit stochastic SER decisions or a terminal-window gate.
            block_options = dict(step_charging_options)
            initial_epoch = int(block_options.pop("initial_sampling_epoch", 0))
            block_options.pop("terminal_window_s", None)
            block_options.update(
                maximum_steps=int(stationarity_block_steps), timestep_policy="fixed",
                stop_on_saturation=False, scramble_mode="fresh",
                compatible_q1_charge_state=True,
                initial_sampling_epoch=initial_epoch)
            initial_charge = (
                lump_mixed_surface_density_3d(
                    poisson, verts, faces, current_sigma, face_conductor_id,
                    grid_origin=potential_origin, grid_spacing=potential_spacing,
                    coordinate_length_unit_m=current_geometry.mesh_length_unit_m)
                if poisson.has_floating_conductors else
                poisson.canonicalize_charge(lump_triangle_sheet_charge_3d(
                    shape=poisson.shape, vertices=verts, faces=faces,
                    sigma_c_per_m2=current_sigma, grid_origin=potential_origin,
                    grid_spacing=potential_spacing,
                    coordinate_length_unit_m=current_geometry.mesh_length_unit_m)))
            initial_potential, _ = poisson.solve(initial_charge)
            charging_first = integrate_charging(current_sigma, block_options)
            block_options["initial_sampling_epoch"] = int(
                charging_first.diagnostics["resume_sampling_epoch"])
            charging = integrate_charging(charging_first.sigma_c_per_m2, block_options)

            first_scoring = []
            second_scoring = []
            next_epoch = int(charging.diagnostics["resume_sampling_epoch"]) + 1
            audit_options = dict(block_options)
            audit_options.update(maximum_steps=0, stop_on_saturation=False)
            for replicate in range(int(stationarity_scoring_replicates)):
                audit_options["initial_sampling_epoch"] = next_epoch + replicate
                first_scoring.append(integrate_charging(
                    charging_first.sigma_c_per_m2, audit_options))
            next_epoch += int(stationarity_scoring_replicates)
            for replicate in range(int(stationarity_scoring_replicates)):
                audit_options["initial_sampling_epoch"] = next_epoch + replicate
                second_scoring.append(integrate_charging(
                    charging.sigma_c_per_m2, audit_options))

            def complete_scoring(scoring):
                return tuple(complete_transport(
                    item.final_step.transport, item.potential_v,
                    item.history[-1]["sampling_seed"]) for item in scoring)

            first_transports = complete_scoring(first_scoring)
            second_transports = complete_scoring(second_scoring)
            first_transport = average_boundary_transport_results_3d(*first_transports)
            transport = average_boundary_transport_results_3d(*second_transports)
            first_feature = advance_profile(
                first_transport, first_scoring[0].history[-1]["sampling_seed"])
            feature = advance_profile(
                transport, second_scoring[0].history[-1]["sampling_seed"])
            first_replica_features = tuple(advance_profile(
                value, item.history[-1]["sampling_seed"])
                for value, item in zip(first_transports, first_scoring))
            second_replica_features = tuple(advance_profile(
                value, item.history[-1]["sampling_seed"])
                for value, item in zip(second_transports, second_scoring))

            first_positive, _ = _replicate_mean_and_standard_error([
                item.final_step.positive_face_current_density_a_m2
                for item in first_scoring])
            first_negative, _ = _replicate_mean_and_standard_error([
                item.final_step.negative_face_current_density_a_m2
                for item in first_scoring])
            second_positive, _ = _replicate_mean_and_standard_error([
                item.final_step.positive_face_current_density_a_m2
                for item in second_scoring])
            second_negative, _ = _replicate_mean_and_standard_error([
                item.final_step.negative_face_current_density_a_m2
                for item in second_scoring])
            _, first_net_error = _replicate_mean_and_standard_error([
                item.final_step.face_current_density_a_m2 for item in first_scoring])
            _, second_net_error = _replicate_mean_and_standard_error([
                item.final_step.face_current_density_a_m2 for item in second_scoring])

            first_flux_replicates = [
                _transport_face_flux_map_3d(value, len(faces))
                for value in first_transports]
            second_flux_replicates = [
                _transport_face_flux_map_3d(value, len(faces))
                for value in second_transports]
            species = set(first_flux_replicates[0])
            if (any(set(value) != species for value in first_flux_replicates)
                    or any(set(value) != species for value in second_flux_replicates)):
                raise RuntimeError("stationarity scoring changed delivered species inventory")
            first_flux = _transport_face_flux_map_3d(first_transport, len(faces))
            second_flux = _transport_face_flux_map_3d(transport, len(faces))
            first_flux_error = {}
            second_flux_error = {}
            for name in species:
                _, first_flux_error[name] = _replicate_mean_and_standard_error([
                    value[name] for value in first_flux_replicates])
                _, second_flux_error[name] = _replicate_mean_and_standard_error([
                    value[name] for value in second_flux_replicates])

            length_unit = current_geometry.mesh_length_unit_m
            first_actual_velocity = (
                np.asarray(first_feature.face_velocity_mesh_units_s) * length_unit)
            second_actual_velocity = (
                np.asarray(feature.face_velocity_mesh_units_s) * length_unit)
            first_velocity_mean, first_velocity_error = _replicate_mean_and_standard_error([
                np.asarray(item.face_velocity_mesh_units_s) * length_unit
                for item in first_replica_features])
            second_velocity_mean, second_velocity_error = _replicate_mean_and_standard_error([
                np.asarray(item.face_velocity_mesh_units_s) * length_unit
                for item in second_replica_features])
            # Price any nonlinearity between chemistry(mean transport) and mean chemistry as an
            # additional deterministic uncertainty bound rather than silently assuming linearity.
            first_velocity_error += np.abs(first_actual_velocity - first_velocity_mean)
            second_velocity_error += np.abs(second_actual_velocity - second_velocity_mean)
            physical_area = np.asarray(areas, dtype=float) * length_unit ** 2
            first_block = ProfileChargingStationarityBlock3D(
                potential_start_v=initial_potential,
                potential_end_v=charging_first.potential_v,
                positive_face_current_density_a_m2=first_positive,
                negative_face_current_density_a_m2=first_negative,
                net_face_current_standard_error_a_m2=first_net_error,
                face_area_m2=physical_area,
                species_face_flux_m2_s=first_flux,
                species_face_flux_standard_error_m2_s=first_flux_error,
                profile_velocity_m_s=first_actual_velocity,
                profile_velocity_standard_error_m_s=first_velocity_error,
                profile_increment_m=-first_actual_velocity * profile_step_duration,
                independent_replicates=len(first_scoring),
                scoring_sampling_epochs=tuple(
                    item.history[-1]["sampling_epoch"] for item in first_scoring),
                duration_s=profile_step_duration, exact_hard_visibility=True)
            second_block = ProfileChargingStationarityBlock3D(
                potential_start_v=charging_first.potential_v,
                potential_end_v=charging.potential_v,
                positive_face_current_density_a_m2=second_positive,
                negative_face_current_density_a_m2=second_negative,
                net_face_current_standard_error_a_m2=second_net_error,
                face_area_m2=physical_area,
                species_face_flux_m2_s=second_flux,
                species_face_flux_standard_error_m2_s=second_flux_error,
                profile_velocity_m_s=second_actual_velocity,
                profile_velocity_standard_error_m_s=second_velocity_error,
                profile_increment_m=-second_actual_velocity * profile_step_duration,
                independent_replicates=len(second_scoring),
                scoring_sampling_epochs=tuple(
                    item.history[-1]["sampling_epoch"] for item in second_scoring),
                duration_s=profile_step_duration, exact_hard_visibility=True)
            profile_stationarity = assess_profile_charging_stationarity_3d(
                first_block, second_block, stationarity_contract)
            if not profile_stationarity.passed:
                raise SurfaceChargingSaturationError(
                    f"co-evolution step {step_index + 1} failed profile-relevant "
                    f"stationarity: {'; '.join(profile_stationarity.reasons)}",
                    charging.sigma_c_per_m2, charging.history,
                    charging.accepted_steps, charging.rejected_steps,
                    charging.physical_time_s, charging.pseudo_time_s,
                    resume_sampling_epoch=charging.diagnostics["resume_sampling_epoch"],
                    resume_stochastic_gain_age_steps=charging.diagnostics[
                        "resume_stochastic_gain_age_steps"])
        next_verts, next_faces, next_centroids, _next_areas = extract_mesh_3d(
            feature.geometry.phi, feature.geometry.dx)
        next_material = _face_material_ids(next_centroids, feature.geometry)
        normal_displacement = (
            -feature.face_velocity_mesh_units_s * profile_step_duration)
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
                step=step_index + 1, profile_duration_s=profile_step_duration,
                profile_motion_enabled=bool(profile_motion_enabled),
                bias_mode=bias_mode,
                saturation_required=bool(bias_mode == "quasi_static"),
                charging_acceptance=charging_acceptance,
                charging_saturated=bool(
                    charging.converged if profile_stationarity is None
                    else profile_stationarity.passed),
                signed_r2_satisfied=charging.converged,
                profile_stationarity_satisfied=(
                    None if profile_stationarity is None else profile_stationarity.passed),
                profile_stationarity_potential_drift_upper_v=(
                    None if profile_stationarity is None
                    else profile_stationarity.potential_drift_upper_v),
                profile_stationarity_current_relative_l1_upper=(
                    None if profile_stationarity is None
                    else profile_stationarity.current_relative_l1_upper),
                profile_stationarity_transported_flux_relative_l1_upper=(
                    None if profile_stationarity is None
                    else profile_stationarity.transported_flux_relative_l1_upper),
                profile_stationarity_profile_velocity_relative_l1_upper=(
                    None if profile_stationarity is None
                    else profile_stationarity.profile_velocity_relative_l1_upper),
                profile_stationarity_profile_increment_difference_upper_m=(
                    None if profile_stationarity is None
                    else profile_stationarity.profile_increment_difference_upper_m),
                charging_accepted_steps=charging.accepted_steps,
                charging_rejected_steps=charging.rejected_steps,
                charging_physical_time_s=charging.physical_time_s,
                charging_pseudo_time_s=charging.pseudo_time_s,
                retained_node_rms_relative_current_imbalance=(
                    charging.diagnostics["retained_node_rms_relative_current_imbalance"]),
                retained_node_max_relative_current_imbalance=(
                    charging.diagnostics["retained_node_max_relative_current_imbalance"]),
                instantaneous_potential_rate_max_v_s=(
                    charging.diagnostics["final_instantaneous_potential_rate_max_v_s"]),
                potential_rate_max_v_s=(
                    charging.diagnostics["final_potential_rate_max_v_s"]),
                gate_evaluation_mode=charging.diagnostics["gate_evaluation_mode"],
                terminal_window_s=charging.diagnostics["terminal_window_s"],
                terminal_window_ready=charging.diagnostics["terminal_window_ready"],
                patch_scales_m=charging.diagnostics["patch_scales_m"],
                instantaneous_patch_max_relative_imbalance=tuple(
                    charging.history[-1]["patch_max_relative_imbalance"]),
                patch_max_relative_imbalance=tuple(
                    item.b2_maximum_ion_normalized_imbalance
                    for item in charging.patch_balance),
                patch_symmetric_max_relative_imbalance=tuple(
                    item.maximum_relative_imbalance for item in charging.patch_balance),
                maximum_transport_lineage_replay_count=charging.diagnostics[
                    "maximum_transport_lineage_replay_count"],
                maximum_transport_lineage_replay_fraction=charging.diagnostics[
                    "maximum_transport_lineage_replay_fraction"],
                maximum_transport_edge_launch_inset_count=charging.diagnostics[
                    "maximum_transport_edge_launch_inset_count"],
                maximum_transport_trajectory_horizon_extension_count=charging.diagnostics[
                    "maximum_transport_trajectory_horizon_extension_count"],
                maximum_transport_trajectory_final_max_steps=charging.diagnostics[
                    "maximum_transport_trajectory_final_max_steps"],
                response_tail_closure_l1_current_error_bound_relative=(
                    charging.diagnostics[
                        "final_response_tail_closure_l1_current_error_bound_relative"]),
                maximum_response_bounce_budget=charging.diagnostics[
                    "maximum_response_bounce_budget"],
                maximum_response_bounce_budget_extension_count=charging.diagnostics[
                    "maximum_response_bounce_budget_extension_count"],
                surface_product_redeposition_relative_balance_error=(
                    feature.diagnostics[
                        "product_redeposition_relative_balance_error"]),
                neutral_forward_scatter_rate_s=feature.diagnostics[
                    "neutral_forward_scatter_rate_s"],
                neutral_forward_scatter_landed_rate_s=feature.diagnostics[
                    "neutral_forward_scatter_landed_rate_s"],
                neutral_forward_scatter_particle_balance_error=feature.diagnostics[
                    "neutral_forward_scatter_particle_balance_error"],
                neutral_forward_scatter_energy_balance_error=feature.diagnostics[
                    "neutral_forward_scatter_energy_balance_error"],
                remap_relative_charge_balance_error=remap.relative_charge_balance_error,
                retained_positive_charge_c=remap.retained_positive_charge_c,
                retained_negative_charge_c=remap.retained_negative_charge_c,
                removed_positive_charge_c=remap.removed_positive_charge_c,
                removed_negative_charge_c=remap.removed_negative_charge_c,
                removed_charge_c=remap.removed_net_charge_c,
                charged_surface_response_manifest=response_manifest,
                neutral_forward_scatter_manifest=forward_scatter_manifest),
            profile_stationarity=profile_stationarity)
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
        ("waveform-resolved C3 is smoke-tested co-simulation, not experimentally validated"
         if bias_mode == "waveform_resolved" else
         "physical-time finite-count co-evolution is an ensemble prediction; one realization "
         "is not a deterministic twist forecast"),
    )) + (
        "C3 refinement and independent final-audit evidence are external campaign gates",
    ) + ((
        "profile-relevant stationarity contract is a draft numerical acceptance path and does not "
        "authorize experimental claims",
    ) if charging_acceptance == "profile_stationary" else ())
    nonpredictive = tuple(dict.fromkeys(
        name for result in results
        for name in result.feature.validity.nonpredictive_parameters))
    manifest_charging_options = dict(options)
    # Runtime persistence hooks affect neither the physical operator nor its numerical controls
    # and are intentionally not serializable.  Their artifacts carry their own schema/checksums.
    manifest_charging_options.pop("progress_callback", None)
    manifest = dict(
        schema=CHARGING_RUN_MANIFEST_SCHEMA,
        mode=bias_mode, n_steps=int(n_steps), duration_s=float(duration_s),
        seed=int(seed), n_position=int(n_position),
        initial_geometry=_geometry_manifest(geometry),
        species_role={str(name): str(value) for name, value in sorted(role.items())},
        etchable_material_ids=[
            int(value) for value in sorted(set(etchable_material_ids))],
        source_bounds=[float(value) for value in source_bounds],
        source_z=float(source_z),
        potential_origin=[float(value) for value in potential_origin],
        potential_spacing=_manifest_value(
            np.asarray(potential_spacing, dtype=float), path="potential_spacing"),
        trajectory_fixed_dt=float(trajectory_fixed_dt),
        trajectory_max_steps=int(trajectory_max_steps),
        trajectory_adaptive_horizon=bool(trajectory_adaptive_horizon),
        trajectory_emergency_max_steps=(
            None if trajectory_emergency_max_steps is None
            else int(trajectory_emergency_max_steps)),
        periodic_lateral=bool(periodic_lateral),
        transport_device=(None if transport_device is None else str(transport_device)),
        cfl_number=float(cfl_number),
        reinitialize=bool(reinitialize),
        reinitialization_method=str(reinitialization_method),
        profile_motion_enabled=bool(profile_motion_enabled),
        charging_acceptance=charging_acceptance,
        stationarity_contract=(None if stationarity_contract is None else dict(
            revision=stationarity_contract.revision,
            potential_drift_tolerance_v=(
                stationarity_contract.potential_drift_tolerance_v),
            current_relative_l1_tolerance=(
                stationarity_contract.current_relative_l1_tolerance),
            transported_flux_relative_l1_tolerance=(
                stationarity_contract.transported_flux_relative_l1_tolerance),
            profile_velocity_relative_l1_tolerance=(
                stationarity_contract.profile_velocity_relative_l1_tolerance),
            profile_increment_tolerance_m=(
                stationarity_contract.profile_increment_tolerance_m),
            minimum_independent_replicates=(
                stationarity_contract.minimum_independent_replicates),
            confidence_multiplier=stationarity_contract.confidence_multiplier,
            experimental_claim_authorized=(
                stationarity_contract.authorizes_experimental_claims))),
        stationarity_block_steps=(
            None if stationarity_block_steps is None else int(stationarity_block_steps)),
        stationarity_scoring_replicates=int(stationarity_scoring_replicates),
        charging_options=_manifest_value(
            manifest_charging_options, path="charging_options"),
        initial_surface_state_supplied=bool(initial_surface_state is not None),
        initial_surface_state_mesh_fingerprint=initial_surface_state_mesh_fingerprint,
        restart_source_manifest_sha256=restart_source_manifest_sha256,
        boundary=_boundary_manifest(boundary),
        surface_mechanism=_surface_mechanism_manifest(mechanism),
        neutral_radiosity=(None if neutral_radiosity_options is None else
                           _manifest_value(
                               neutral_radiosity_options, path="neutral_radiosity")),
        neutral_forward_scatter=(
            None if forward_scatter_manifest is None else
            _manifest_value(forward_scatter_manifest)),
        neutral_forward_scatter_options=(
            None if neutral_forward_scatter_options is None else
            _manifest_value(
                neutral_forward_scatter_options,
                path="neutral_forward_scatter_options")),
        waveform=(None if not waveform else [dict(
            duration_s=item.duration_s, boundary=_boundary_manifest(item.boundary))
            for item in waveform]),
        charged_surface_response=(None if response_manifest is None else
                                  _manifest_value(response_manifest)),
        surface_product_redeposition=(
            None if surface_product_redeposition_options is None else
            _manifest_value(
                surface_product_redeposition_options,
                path="surface_product_redeposition")),
        recovery_and_error_budget=dict(
            inline_recovery=dict(
                float64_lineage_replay_count=sum(
                    int(item.diagnostics["maximum_transport_lineage_replay_count"])
                    for item in results),
                trajectory_horizon_extension_count=sum(
                    int(item.diagnostics[
                        "maximum_transport_trajectory_horizon_extension_count"])
                    for item in results),
                charged_cascade_bounce_budget_extension_count=sum(
                    int(item.diagnostics[
                        "maximum_response_bounce_budget_extension_count"])
                    for item in results)),
            priced_absorption=dict(
                maximum_charged_cascade_tail_l1_current_error_bound_relative=max(
                    float(item.diagnostics[
                        "response_tail_closure_l1_current_error_bound_relative"])
                    for item in results)),
            conservation=dict(
                maximum_charge_remap_relative_balance_error=max(
                    float(item.diagnostics["remap_relative_charge_balance_error"])
                    for item in results),
                maximum_surface_product_redeposition_relative_balance_error=max(
                    0.0 if item.diagnostics[
                        "surface_product_redeposition_relative_balance_error"] is None
                    else float(item.diagnostics[
                        "surface_product_redeposition_relative_balance_error"])
                    for item in results)),
            hard_refusal_policy=(
                "incomplete trajectory/cascade, conservation breach, cross-material deposition, "
                "topology event, or corrupt restart state refuses rather than being discarded")),
        experimental_claim=bool(experimental_claim),
        observable_tolerances=[dict(
            observable=item.observable, tolerance=item.tolerance,
            benchmark_uncertainty_including_digitization=(
                item.benchmark_uncertainty_including_digitization),
            feature_extent_m=item.feature_extent_m)
            for item in tolerance_contracts],
        exact_operator=(
            "hard visibility; physical-time charge ODE; compatible Q1 projection"
            + ("; neutralized hard-visibility SiO2 forward scatter"
               if neutral_forward_scatter is not None else "")
            + ("; fixed-geometry transport diagnostic"
               if not profile_motion_enabled else "")),
        convergence_contract_revision=(
            "CCA-2026-07-13-R2" if charging_acceptance == "signed_r2"
            else stationarity_contract.revision))
    return ChargingCoevolution3DResult(
        current_geometry, current_state, current_fingerprint, current_sigma,
        tuple(results), float(duration_s),
        FeatureStepValidity(
            not reasons, reasons, limitations, not nonpredictive, nonpredictive),
        manifest)
