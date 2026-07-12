"""Explicit public entry point for the dimensional common 3-D feature engine.

This API intentionally accepts physical contract objects rather than recipe-name shortcuts.  It lives
beside the legacy ViennaPS-shaped compatibility API while mechanisms are re-earned; the two paths never
silently substitute for one another.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from time import perf_counter
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import PlasmaBoundaryState
from .feature_step_3d import FeatureGeometry3D, FeatureSolve3DResult, solve_feature_3d


COMMON_FEATURE_ENGINE = "feature-3d-v1"

_REQUIRED_SOLVER_ARGUMENTS = frozenset({
    "geometry", "boundary", "species_role", "mechanism", "etchable_material_ids",
    "duration_s", "n_steps", "source_bounds", "source_z",
})


@dataclass(frozen=True)
class PhysicalResult:
    """Result from the common dimensional engine, with unambiguous provenance."""

    solve: FeatureSolve3DResult
    wall_time_s: float
    engine: str = COMMON_FEATURE_ENGINE

    def __post_init__(self):
        if not isinstance(self.solve, FeatureSolve3DResult):
            raise TypeError("solve must be a FeatureSolve3DResult")
        if not np.isfinite(self.wall_time_s) or self.wall_time_s < 0.0:
            raise ValueError("wall_time_s must be finite and nonnegative")
        if self.engine != COMMON_FEATURE_ENGINE:
            raise ValueError("invalid common feature engine identifier")

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
        start = perf_counter()
        result = solve_feature_3d(
            self.geometry, self.boundary, self.species_role, self.mechanism,
            etchable_material_ids=self.etchable_material_ids,
            duration_s=self.duration_s, n_steps=self.n_steps,
            source_bounds=self.source_bounds, source_z=self.source_z,
            **self.solver_options)
        return PhysicalResult(result, perf_counter() - start)
