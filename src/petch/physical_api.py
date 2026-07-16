"""Explicit public entry point for the dimensional common 3-D feature engine.

This API intentionally accepts physical contract objects rather than recipe-name shortcuts.  It lives
beside the legacy ViennaPS-shaped compatibility API while mechanisms are re-earned; the two paths never
silently substitute for one another.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from time import perf_counter
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import PlasmaBoundaryState
from .charging_coevolution_3d import (
    ChargingCoevolution3DResult,
    _boundary_manifest,
    _geometry_manifest,
    _manifest_value,
    _surface_mechanism_manifest,
    solve_charging_coevolution_3d,
)
from .charging_checkpoint_3d import PhysicalChargingCheckpoint3D
from .feature_step_3d import FeatureGeometry3D, FeatureSolve3DResult, solve_feature_3d


COMMON_FEATURE_ENGINE = "feature-3d-v1"
COMMON_CHARGING_ENGINE = "feature-charging-coevolution-3d-v1"
COMMON_CHARGING_ENSEMBLE_ENGINE = "feature-charging-ensemble-3d-v1"
COMMON_FEATURE_MANIFEST_SCHEMA = "petch-feature-run-manifest-3d-v1"

_REQUIRED_SOLVER_ARGUMENTS = frozenset({
    "geometry", "boundary", "species_role", "mechanism", "etchable_material_ids",
    "duration_s", "n_steps", "source_bounds", "source_z",
})
_REQUIRED_CHARGING_SOLVER_ARGUMENTS = _REQUIRED_SOLVER_ARGUMENTS | frozenset({
    "charging_system_builder", "potential_origin", "potential_spacing", "charging_options",
})


def _contract_values_equal(first, second):
    """Equality for immutable public-process fields that may contain NumPy arrays."""
    if first is second:
        return True
    if isinstance(first, np.ndarray) or isinstance(second, np.ndarray):
        try:
            return np.array_equal(np.asarray(first), np.asarray(second), equal_nan=True)
        except (TypeError, ValueError):
            return False
    if isinstance(first, Mapping) and isinstance(second, Mapping):
        return (set(first) == set(second)
                and all(_contract_values_equal(first[key], second[key]) for key in first))
    if (isinstance(first, (tuple, list)) and isinstance(second, (tuple, list))
            and len(first) == len(second)):
        return all(_contract_values_equal(a, b) for a, b in zip(first, second))
    try:
        result = first == second
    except (TypeError, ValueError):
        return False
    return bool(result) if isinstance(result, (bool, np.bool_)) else False


@dataclass(frozen=True)
class PhysicalResult:
    """Result from the common dimensional engine, with unambiguous provenance."""

    solve: FeatureSolve3DResult
    wall_time_s: float
    run_manifest: Mapping[str, object]
    engine: str = COMMON_FEATURE_ENGINE

    def __post_init__(self):
        if not isinstance(self.solve, FeatureSolve3DResult):
            raise TypeError("solve must be a FeatureSolve3DResult")
        if not np.isfinite(self.wall_time_s) or self.wall_time_s < 0.0:
            raise ValueError("wall_time_s must be finite and nonnegative")
        if self.engine != COMMON_FEATURE_ENGINE:
            raise ValueError("invalid common feature engine identifier")
        object.__setattr__(self, "run_manifest", MappingProxyType(dict(self.run_manifest)))

    @property
    def geometry(self):
        return self.solve.geometry

    @property
    def surface_state(self):
        return self.solve.surface_state

    @property
    def steps(self):
        return self.solve.steps

    @property
    def duration_s(self):
        return self.solve.duration_s

    @property
    def validity(self):
        return self.solve.validity

    def measure_trench_profile(self, **measurement_contract):
        """Measure declared notch/bow bands on the final common-engine geometry."""
        from .profile_observables_3d import measure_trench_profile_observables_3d
        return measure_trench_profile_observables_3d(
            self.geometry, **measurement_contract)


@dataclass(frozen=True)
class PhysicalProcess:
    """One fully explicit common-engine feature process.

    ``solver_options`` may select numerical accuracy and supported physical operators, but may not
    replace the required physical case fields.  Unknown chemistry parameters remain properties of the
    supplied mechanism, where their provenance and validity are reported.
    """

    geometry: FeatureGeometry3D
    boundary: PlasmaBoundaryState
    species_role: Mapping[str, str]
    mechanism: object
    etchable_material_ids: tuple[int, ...]
    duration_s: float
    n_steps: int
    source_bounds: tuple[float, float, float, float]
    source_z: float
    solver_options: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.geometry, FeatureGeometry3D):
            raise TypeError("geometry must be FeatureGeometry3D")
        if not isinstance(self.boundary, PlasmaBoundaryState):
            raise TypeError("boundary must be PlasmaBoundaryState")
        role = dict(self.species_role)
        if not role or any(not isinstance(name, str) or not isinstance(value, str)
                           or not name or not value for name, value in role.items()):
            raise ValueError("species_role must map nonempty species names to nonempty roles")
        species = {item.name for item in self.boundary.species}
        if set(role) != species:
            raise ValueError("species_role must cover every and only boundary species")
        etchable = tuple(int(value) for value in self.etchable_material_ids)
        if not etchable or any(value <= 0 for value in etchable) or len(set(etchable)) != len(etchable):
            raise ValueError("etchable_material_ids must be unique positive ids")
        bounds = tuple(float(value) for value in self.source_bounds)
        if (len(bounds) != 4 or np.any(~np.isfinite(bounds))
                or bounds[1] <= bounds[0] or bounds[3] <= bounds[2]):
            raise ValueError("source_bounds must be finite (xmin,xmax,ymin,ymax) bounds")
        if (not np.isfinite(self.duration_s) or self.duration_s < 0.0
                or int(self.n_steps) != self.n_steps or int(self.n_steps) <= 0
                or not np.isfinite(self.source_z)):
            raise ValueError("invalid process duration, step count, or source height")
        options = dict(self.solver_options)
        overlap = _REQUIRED_SOLVER_ARGUMENTS & set(options)
        if overlap:
            raise ValueError(
                "solver_options cannot override required case fields: " + ", ".join(sorted(overlap)))
        object.__setattr__(self, "species_role", MappingProxyType(role))
        object.__setattr__(self, "etchable_material_ids", etchable)
        object.__setattr__(self, "duration_s", float(self.duration_s))
        object.__setattr__(self, "n_steps", int(self.n_steps))
        object.__setattr__(self, "source_bounds", bounds)
        object.__setattr__(self, "source_z", float(self.source_z))
        object.__setattr__(self, "solver_options", MappingProxyType(options))

    @property
    def engine(self):
        return COMMON_FEATURE_ENGINE

    def run(self):
        manifest = dict(
            schema=COMMON_FEATURE_MANIFEST_SCHEMA,
            engine=COMMON_FEATURE_ENGINE,
            mode="uncharged_physical_profile",
            initial_geometry=_geometry_manifest(self.geometry),
            boundary=_boundary_manifest(self.boundary),
            species_role={
                str(name): str(value) for name, value in sorted(self.species_role.items())},
            surface_mechanism=_surface_mechanism_manifest(self.mechanism),
            etchable_material_ids=[int(value) for value in self.etchable_material_ids],
            duration_s=self.duration_s, n_steps=self.n_steps,
            source_bounds=[float(value) for value in self.source_bounds],
            source_z=self.source_z,
            solver_options=_manifest_value(
                self.solver_options, path="physical_process.solver_options"),
            exact_operator="hard visibility; material-routed surface law; level-set motion")
        start = perf_counter()
        result = solve_feature_3d(
            self.geometry, self.boundary, self.species_role, self.mechanism,
            etchable_material_ids=self.etchable_material_ids,
            duration_s=self.duration_s, n_steps=self.n_steps,
            source_bounds=self.source_bounds, source_z=self.source_z,
            **self.solver_options)
        return PhysicalResult(result, perf_counter() - start, manifest)


@dataclass(frozen=True)
class PhysicalChargingResult:
    """Result from the unified charge/transport/chemistry/profile engine."""

    solve: ChargingCoevolution3DResult
    wall_time_s: float
    engine: str = COMMON_CHARGING_ENGINE

    def __post_init__(self):
        if not isinstance(self.solve, ChargingCoevolution3DResult):
            raise TypeError("solve must be a ChargingCoevolution3DResult")
        if not np.isfinite(self.wall_time_s) or self.wall_time_s < 0.0:
            raise ValueError("wall_time_s must be finite and nonnegative")
        if self.engine != COMMON_CHARGING_ENGINE:
            raise ValueError("invalid common charging engine identifier")

    @property
    def geometry(self):
        return self.solve.geometry

    @property
    def surface_state(self):
        return self.solve.surface_state

    @property
    def surface_charge_c_per_m2(self):
        return self.solve.sigma_c_per_m2

    @property
    def steps(self):
        return self.solve.steps

    @property
    def duration_s(self):
        return self.solve.duration_s

    @property
    def validity(self):
        return self.solve.validity

    @property
    def run_manifest(self):
        return self.solve.run_manifest

    def measure_trench_profile(self, **measurement_contract):
        """Measure declared notch/bow bands on the final charged geometry."""
        from .profile_observables_3d import measure_trench_profile_observables_3d
        return measure_trench_profile_observables_3d(
            self.geometry, **measurement_contract)


@dataclass(frozen=True)
class PhysicalChargingProcess:
    """One explicit run of the unified 3-D charging co-evolution engine.

    This public process does not silently enable charging on :class:`PhysicalProcess`.  Callers must
    supply the Poisson builder, charge-integration controls, and electrostatic coordinate contract.
    Optional reflection/SEE, waveform, acceptance, and transport controls remain explicit
    ``solver_options`` and are provenance-checked by the underlying C3 engine.
    """

    geometry: FeatureGeometry3D
    boundary: PlasmaBoundaryState
    species_role: Mapping[str, str]
    mechanism: object
    charging_system_builder: object
    etchable_material_ids: tuple[int, ...]
    duration_s: float
    n_steps: int
    source_bounds: tuple[float, float, float, float]
    source_z: float
    potential_origin: tuple[float, float, float]
    potential_spacing: object
    charging_options: Mapping[str, object]
    solver_options: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        # Reuse the uncharged public contract for the common geometry/boundary/profile fields.
        validated = PhysicalProcess(
            self.geometry, self.boundary, self.species_role, self.mechanism,
            self.etchable_material_ids, self.duration_s, self.n_steps,
            self.source_bounds, self.source_z)
        if not callable(self.charging_system_builder):
            raise TypeError("charging_system_builder must be callable")
        origin = np.asarray(self.potential_origin, dtype=float)
        spacing = np.asarray(self.potential_spacing, dtype=float)
        if (origin.shape != (3,) or np.any(~np.isfinite(origin))
                or spacing.size not in {1, 3} or np.any(~np.isfinite(spacing))
                or np.any(spacing <= 0.0)):
            raise ValueError("invalid electrostatic coordinate origin or spacing")
        charging_options = dict(self.charging_options)
        if not charging_options:
            raise ValueError("charging_options must explicitly declare the charging integration")
        solver_options = dict(self.solver_options)
        overlap = _REQUIRED_CHARGING_SOLVER_ARGUMENTS & set(solver_options)
        if overlap:
            raise ValueError(
                "solver_options cannot override required charged case fields: "
                + ", ".join(sorted(overlap)))
        object.__setattr__(self, "species_role", validated.species_role)
        object.__setattr__(self, "etchable_material_ids", validated.etchable_material_ids)
        object.__setattr__(self, "duration_s", validated.duration_s)
        object.__setattr__(self, "n_steps", validated.n_steps)
        object.__setattr__(self, "source_bounds", validated.source_bounds)
        object.__setattr__(self, "source_z", validated.source_z)
        object.__setattr__(self, "potential_origin", tuple(float(value) for value in origin))
        object.__setattr__(
            self, "potential_spacing",
            float(spacing) if spacing.ndim == 0 else tuple(float(value) for value in spacing))
        object.__setattr__(self, "charging_options", MappingProxyType(charging_options))
        object.__setattr__(self, "solver_options", MappingProxyType(solver_options))

    @property
    def engine(self):
        return COMMON_CHARGING_ENGINE

    def run(self):
        start = perf_counter()
        result = solve_charging_coevolution_3d(
            self.geometry, self.boundary, self.species_role, self.mechanism,
            charging_system_builder=self.charging_system_builder,
            etchable_material_ids=self.etchable_material_ids,
            duration_s=self.duration_s, n_steps=self.n_steps,
            source_bounds=self.source_bounds, source_z=self.source_z,
            potential_origin=self.potential_origin,
            potential_spacing=self.potential_spacing,
            charging_options=self.charging_options,
            **self.solver_options)
        return PhysicalChargingResult(result, perf_counter() - start)

    def continue_from(
            self, previous: PhysicalChargingResult, *, duration_s, n_steps,
            seed=None, continuation_seed_stride=1000003):
        """Build the next profile batch without dropping remapped charge or surface state.

        A distinct seed is derived by default so a batch boundary cannot silently replay the prior
        batch's quadrature/arrival stream.  The returned process is immutable and has not run yet.
        """
        if not isinstance(previous, PhysicalChargingResult):
            raise TypeError("previous must be a PhysicalChargingResult")
        if (int(continuation_seed_stride) != continuation_seed_stride
                or continuation_seed_stride <= 0):
            raise ValueError("continuation_seed_stride must be a positive integer")
        if (not np.isfinite(duration_s) or duration_s < 0.0
                or int(n_steps) != n_steps or n_steps <= 0):
            raise ValueError("continuation duration and step count are invalid")
        options = dict(self.solver_options)
        overlap = {
            "initial_sigma_c_per_m2", "initial_surface_state",
            "initial_surface_state_mesh_fingerprint"} & set(options)
        if overlap:
            raise ValueError(
                "continuation refuses pre-existing restart state in solver_options: "
                + ", ".join(sorted(overlap)))
        if seed is None:
            seed = (int(options.get("seed", 0))
                    + int(continuation_seed_stride) * len(previous.steps))
        if int(seed) != seed or seed < 0:
            raise ValueError("continuation seed must be a nonnegative integer")
        options.update(
            seed=int(seed),
            initial_sigma_c_per_m2=previous.surface_charge_c_per_m2,
            initial_surface_state=previous.surface_state,
            initial_surface_state_mesh_fingerprint=(
                previous.solve.surface_state_mesh_fingerprint))
        return replace(
            self, geometry=previous.geometry,
            duration_s=float(duration_s), n_steps=int(n_steps),
            solver_options=options)

    def continue_from_checkpoint(
            self, checkpoint: PhysicalChargingCheckpoint3D, *, duration_s, n_steps,
            seed=None, continuation_seed_stride=1000003):
        """Build a continuation process from a safe versioned disk/in-memory checkpoint."""
        if not isinstance(checkpoint, PhysicalChargingCheckpoint3D):
            raise TypeError("checkpoint must be PhysicalChargingCheckpoint3D")
        if (not np.isfinite(duration_s) or duration_s < 0.0
                or int(n_steps) != n_steps or n_steps <= 0
                or int(continuation_seed_stride) != continuation_seed_stride
                or continuation_seed_stride <= 0):
            raise ValueError("invalid checkpoint-continuation controls")
        options = dict(self.solver_options)
        overlap = {
            "initial_sigma_c_per_m2", "initial_surface_state",
            "initial_surface_state_mesh_fingerprint", "restart_source_manifest_sha256",
        } & set(options)
        if overlap:
            raise ValueError(
                "continuation refuses pre-existing restart state in solver_options: "
                + ", ".join(sorted(overlap)))
        if seed is None:
            seed = (int(options.get("seed", 0))
                    + int(continuation_seed_stride) * checkpoint.completed_steps)
        if int(seed) != seed or seed < 0:
            raise ValueError("continuation seed must be a nonnegative integer")
        options.update(
            seed=int(seed), initial_sigma_c_per_m2=checkpoint.sigma_c_per_m2,
            initial_surface_state=checkpoint.restore_surface_state(),
            initial_surface_state_mesh_fingerprint=(
                checkpoint.surface_state_mesh_fingerprint),
            restart_source_manifest_sha256=checkpoint.source_manifest_sha256)
        return replace(
            self, geometry=checkpoint.geometry,
            duration_s=float(duration_s), n_steps=int(n_steps),
            solver_options=options)


@dataclass(frozen=True)
class PhysicalChargingEnsembleResult:
    """Finite-count charging/profile realizations and their grid-level geometry moments."""

    realizations: tuple[PhysicalChargingResult, ...]
    seeds: tuple[int, ...]
    mean_levelset: np.ndarray
    standard_deviation_levelset: np.ndarray
    wall_time_s: float
    engine: str = COMMON_CHARGING_ENSEMBLE_ENGINE

    def __post_init__(self):
        realizations = tuple(self.realizations)
        seeds = tuple(int(value) for value in self.seeds)
        mean = np.asarray(self.mean_levelset, dtype=float).copy()
        deviation = np.asarray(self.standard_deviation_levelset, dtype=float).copy()
        if (len(realizations) < 2 or len(seeds) != len(realizations)
                or len(set(seeds)) != len(seeds)
                or any(not isinstance(item, PhysicalChargingResult) for item in realizations)
                or mean.shape != deviation.shape or mean.size == 0
                or np.any(~np.isfinite(mean)) or np.any(~np.isfinite(deviation))
                or np.any(deviation < 0.0)
                or not np.isfinite(self.wall_time_s) or self.wall_time_s < 0.0
                or self.engine != COMMON_CHARGING_ENSEMBLE_ENGINE):
            raise ValueError("invalid physical charging ensemble result")
        mean.setflags(write=False)
        deviation.setflags(write=False)
        object.__setattr__(self, "realizations", realizations)
        object.__setattr__(self, "seeds", seeds)
        object.__setattr__(self, "mean_levelset", mean)
        object.__setattr__(self, "standard_deviation_levelset", deviation)

    @property
    def realization_count(self):
        return len(self.realizations)

    @property
    def statistical_claim_ready(self):
        # C5 additionally requires N doubling, sample-level doubling, and the isotropy control.
        return False

    def measure_centerline(self, **measurement_contract):
        """Measure twist/tilt statistics from the realized geometries and a declared opening ROI."""
        from .profile_observables_3d import measure_feature_centerline_ensemble_3d
        return measure_feature_centerline_ensemble_3d(
            (item.geometry for item in self.realizations), **measurement_contract)

    def measure_trench_profile(self, **measurement_contract):
        """Measure notch/bow uncertainty from the independent realized geometries."""
        from .profile_observables_3d import measure_trench_profile_ensemble_3d
        return measure_trench_profile_ensemble_3d(
            (item.geometry for item in self.realizations), **measurement_contract)


@dataclass(frozen=True)
class PhysicalChargingEnsembleProcess:
    """Reproducible finite-count ensemble wrapper around one physical-time C3 process.

    This wrapper intentionally accepts only ``physical_time_resolved`` runs with dimensional
    Poisson arrivals.  It does not relabel fresh quadrature scrambles as shot noise, and it never
    promotes a single realization to a deterministic twisting prediction.
    """

    process: PhysicalChargingProcess
    realization_count: int
    seed_stride: int = 1000003

    def __post_init__(self):
        if not isinstance(self.process, PhysicalChargingProcess):
            raise TypeError("process must be PhysicalChargingProcess")
        if (int(self.realization_count) != self.realization_count
                or self.realization_count < 2
                or int(self.seed_stride) != self.seed_stride or self.seed_stride <= 0):
            raise ValueError("an ensemble requires at least two realizations and a seed stride")
        if self.process.solver_options.get("bias_mode") != "physical_time_resolved":
            raise ValueError("charging ensembles require physical_time_resolved co-evolution")
        if self.process.charging_options.get("physical_arrival_statistics") != "poisson":
            raise ValueError("charging ensembles require dimensional Poisson physical arrivals")
        object.__setattr__(self, "realization_count", int(self.realization_count))
        object.__setattr__(self, "seed_stride", int(self.seed_stride))

    @property
    def engine(self):
        return COMMON_CHARGING_ENSEMBLE_ENGINE

    def run_twist_refinement(
            self, sample_refined_process, *, aspect_ratio, measurement_contract,
            refinement_contract):
        """Run one paired C5 condition at base N, doubled N, and doubled ray samples.

        ``self.realization_count`` is the doubled-N population.  Its first half is the nested base-N
        population, so those particles are never recomputed or confounded with a different seed set.
        The sample-refined process must be the identical physical case with exactly twice
        ``n_position`` and the same realization seeds.
        """
        from .profile_observables_3d import measure_feature_centerline_ensemble_3d
        from .twist_campaign_3d import (
            TwistEnsembleRefinementContract3D,
            score_twist_condition_campaign_3d,
        )
        if not isinstance(sample_refined_process, PhysicalChargingEnsembleProcess):
            raise TypeError("sample_refined_process must be a PhysicalChargingEnsembleProcess")
        if not isinstance(refinement_contract, TwistEnsembleRefinementContract3D):
            raise TypeError("refinement_contract must be TwistEnsembleRefinementContract3D")
        if (self.realization_count % 2
                or self.realization_count // 2 < refinement_contract.minimum_realizations
                or sample_refined_process.realization_count != self.realization_count
                or sample_refined_process.seed_stride != self.seed_stride):
            raise ValueError(
                "twist refinement requires nested base N, doubled N, and paired realization seeds")
        base = self.process
        refined = sample_refined_process.process
        identity_fields = (
            "geometry", "boundary", "mechanism", "charging_system_builder")
        if any(getattr(base, field) is not getattr(refined, field) for field in identity_fields):
            raise ValueError("sample refinement must reuse the identical physical case objects")
        equal_fields = (
            "species_role", "etchable_material_ids", "duration_s", "n_steps",
            "source_bounds", "source_z", "potential_origin", "potential_spacing",
            "charging_options")
        if any(not _contract_values_equal(
                getattr(base, field), getattr(refined, field)) for field in equal_fields):
            raise ValueError("sample refinement changed a physical case field")
        base_options = dict(base.solver_options)
        refined_options = dict(refined.solver_options)
        base_n_position = int(base_options.pop("n_position", 256))
        refined_n_position = int(refined_options.pop("n_position", 256))
        if (refined_n_position != 2 * base_n_position
                or not _contract_values_equal(base_options, refined_options)):
            raise ValueError(
                "sample refinement may only double n_position on the identical operator")
        doubled_result = self.run()
        refined_result = sample_refined_process.run()
        if doubled_result.seeds != refined_result.seeds:
            raise RuntimeError("sample refinement did not preserve paired realization seeds")
        base_count = self.realization_count // 2
        base_centerline = measure_feature_centerline_ensemble_3d(
            (item.geometry for item in doubled_result.realizations[:base_count]),
            **dict(measurement_contract))
        doubled_centerline = doubled_result.measure_centerline(**dict(measurement_contract))
        refined_centerline = refined_result.measure_centerline(**dict(measurement_contract))
        return score_twist_condition_campaign_3d(
            base_centerline, doubled_centerline, refined_centerline,
            aspect_ratio=aspect_ratio,
            base_transport_sample_count=base_n_position,
            refined_transport_sample_count=refined_n_position,
            base_seeds=doubled_result.seeds[:base_count],
            doubled_seeds=doubled_result.seeds,
            doubled_sample_seeds=refined_result.seeds,
            contract=refinement_contract)

    def run(self):
        start = perf_counter()
        base_seed = int(self.process.solver_options.get("seed", 0))
        seeds = tuple(
            base_seed + self.seed_stride * index for index in range(self.realization_count))
        realizations = []
        for seed in seeds:
            options = dict(self.process.solver_options)
            options["seed"] = seed
            realizations.append(replace(
                self.process, solver_options=options).run())
        levelsets = np.stack([
            item.geometry.phi for item in realizations])
        return PhysicalChargingEnsembleResult(
            tuple(realizations), seeds,
            np.mean(levelsets, axis=0), np.std(levelsets, axis=0, ddof=1),
            perf_counter() - start)
