"""Finite-count physical arrival sampling for stochastic feature charging.

Randomized QMC scrambles quantify/integrate a continuum boundary distribution; changing their seed
does not by itself define physical shot noise.  This module supplies the missing dimensional bridge.
For one sparse face-resolved transport measure and a declared physical time window,

``lambda_j = event_flux_j * face_area[event_face_j] * duration``

is the expected number of physical particles represented by event ``j``.  Independent Poisson
counts are drawn and converted back to a sparse flux measure.  Therefore the realized landed charge
is exactly ``q e sum(N_j)`` for a perfect absorber, while ensemble mean and variance both equal the
declared arrival count.  No arbitrary "particles per simulation step" scale enters the physics.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .sheath import ECHARGE
from .surface_kinetics import FaceResolvedEnergeticFlux


@dataclass(frozen=True)
class PhysicalArrivalSample3D:
    """One reproducible finite-count realization of a sparse incident population."""

    population: FaceResolvedEnergeticFlux
    source_event_count: np.ndarray
    source_expected_event_count: np.ndarray
    realized_face_count: np.ndarray
    expected_face_count: np.ndarray
    duration_s: float
    seed: int
    diagnostics: Mapping[str, object]

    def __post_init__(self):
        if not isinstance(self.population, FaceResolvedEnergeticFlux):
            raise TypeError("population must be FaceResolvedEnergeticFlux")
        count = np.asarray(self.source_event_count, dtype=np.int64).copy()
        expected = np.asarray(self.source_expected_event_count, dtype=float).copy()
        face_count = np.asarray(self.realized_face_count, dtype=np.int64).copy()
        face_expected = np.asarray(self.expected_face_count, dtype=float).copy()
        if (count.ndim != 1 or expected.shape != count.shape or np.any(count < 0)
                or np.any(~np.isfinite(expected)) or np.any(expected < 0.0)
                or face_count.shape != (self.population.face_count,)
                or face_expected.shape != face_count.shape or np.any(face_count < 0)
                or np.any(~np.isfinite(face_expected)) or np.any(face_expected < 0.0)
                or int(self.seed) != self.seed or self.seed < 0
                or not np.isfinite(self.duration_s) or self.duration_s <= 0.0
                or int(np.sum(count)) != int(np.sum(face_count))
                or not np.isclose(
                    float(np.sum(expected)), float(np.sum(face_expected)),
                    rtol=3e-15, atol=np.finfo(float).tiny)):
            raise ValueError("invalid physical-arrival sample")
        for value in (count, expected, face_count, face_expected):
            value.setflags(write=False)
        object.__setattr__(self, "source_event_count", count)
        object.__setattr__(self, "source_expected_event_count", expected)
        object.__setattr__(self, "realized_face_count", face_count)
        object.__setattr__(self, "expected_face_count", face_expected)
        object.__setattr__(self, "duration_s", float(self.duration_s))
        object.__setattr__(self, "seed", int(self.seed))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))

    @property
    def realized_particle_count(self):
        return int(np.sum(self.realized_face_count))

    @property
    def expected_particle_count(self):
        return float(np.sum(self.expected_face_count))

    def landed_charge_c(self, charge_number):
        if int(charge_number) != charge_number:
            raise ValueError("charge_number must be an integer")
        return float(int(charge_number) * ECHARGE * self.realized_particle_count)


def sample_physical_poisson_arrivals_3d(
        population: FaceResolvedEnergeticFlux, face_area_m2, duration_s, *, seed):
    """Draw physical particle counts and return their exact equivalent sparse flux measure.

    Transport events with zero realized count are omitted from the returned population.  Events with
    count greater than one stay compact: their flux contribution is ``N/(area*duration)`` and the
    energy/angle/lineage data remain attached once.  This preserves downstream yield integration
    without allocating one record per physical particle.
    """
    if not isinstance(population, FaceResolvedEnergeticFlux):
        raise TypeError("population must be FaceResolvedEnergeticFlux")
    area = np.asarray(face_area_m2, dtype=float)
    if (area.shape != (population.face_count,) or np.any(~np.isfinite(area))
            or np.any(area <= 0.0) or not np.isfinite(duration_s) or duration_s <= 0.0
            or int(seed) != seed or seed < 0):
        raise ValueError("invalid physical-arrival area, duration, or seed")
    expected = (np.asarray(population.event_flux_m2_s, dtype=float)
                * area[population.event_face] * float(duration_s))
    # NumPy's Poisson sampler protects its signed-64 implementation only near the upper tail.  A
    # refusal here asks callers to shorten the physical window rather than overflow particle counts.
    maximum_safe_mean = np.iinfo(np.int64).max - 16.0 * np.sqrt(np.iinfo(np.int64).max)
    if np.any(expected > maximum_safe_mean):
        raise ValueError("physical arrival window is too long for finite-count sampling")
    count = np.random.default_rng(int(seed)).poisson(expected).astype(np.int64, copy=False)
    selected = count > 0
    selected_face = population.event_face[selected]
    selected_flux = (
        count[selected].astype(float)
        / (area[selected_face] * float(duration_s)))

    def selected_optional(name):
        value = getattr(population, name)
        return None if value is None else np.asarray(value)[selected]

    realized_population = FaceResolvedEnergeticFlux(
        population.name, population.face_count,
        selected_face, selected_flux,
        population.event_energy_eV[selected],
        population.event_cosine_incidence[selected],
        event_position=selected_optional("event_position"),
        event_incident_direction=selected_optional("event_incident_direction"))
    face_count = np.bincount(
        population.event_face, weights=count, minlength=population.face_count).astype(np.int64)
    expected_face_count = np.bincount(
        population.event_face, weights=expected, minlength=population.face_count)
    realized_identity = np.bincount(
        realized_population.event_face,
        weights=(realized_population.event_flux_m2_s
                 * area[realized_population.event_face] * float(duration_s)),
        minlength=population.face_count)
    identity_tolerance = 32.0 * np.finfo(float).eps * max(
        1.0, float(np.max(face_count, initial=0)))
    if not np.allclose(
            realized_identity, face_count, rtol=0.0, atol=identity_tolerance):
        raise RuntimeError("finite-count flux conversion lost a physical arrival")
    return PhysicalArrivalSample3D(
        population=realized_population,
        source_event_count=count,
        source_expected_event_count=expected,
        realized_face_count=face_count,
        expected_face_count=expected_face_count,
        duration_s=float(duration_s), seed=int(seed),
        diagnostics=dict(
            statistics="independent Poisson physical arrivals",
            dimensional_rate="event flux density * physical face area * physical duration",
            compact_multiple_arrivals=True,
            primary_arrivals_only=True,
            response_branching_statistics=(
                "conditional-mean unless the charged surface response supplies a stochastic law")))
