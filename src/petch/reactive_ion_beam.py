"""Measured, species-resolved SiO2 reactive-ion beam closures.

This module intentionally distinguishes a data-anchored surface operator from
an independently validated predictive mechanism.  Karahashi's mass-selected
beam measurements resolve ion identity and energy directly, but only at normal
incidence, without simultaneous neutral radicals, and only from 250--2000 eV.
Interpolation inside that support is useful and auditable.  Extrapolation or
application to an unpublished reactor ion mixture is refused.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from types import MappingProxyType

import numpy as np


@dataclass(frozen=True)
class ReactiveIonYieldSeries:
    species: str
    energy_eV: np.ndarray
    yield_sio2_per_ion: np.ndarray
    digitization_yield_uncertainty: np.ndarray

    def __post_init__(self):
        energy = np.asarray(self.energy_eV, dtype=float).copy()
        yield_value = np.asarray(self.yield_sio2_per_ion, dtype=float).copy()
        uncertainty = np.asarray(
            self.digitization_yield_uncertainty, dtype=float).copy()
        if (not self.species
                or energy.ndim != 1 or energy.size < 2
                or yield_value.shape != energy.shape
                or uncertainty.shape != energy.shape
                or np.any(~np.isfinite(energy))
                or np.any(~np.isfinite(yield_value))
                or np.any(~np.isfinite(uncertainty))
                or np.any(energy <= 0.0)
                or np.any(np.diff(energy) <= 0.0)
                or np.any(yield_value < 0.0)
                or np.any(uncertainty <= 0.0)):
            raise ValueError("invalid reactive-ion beam yield series")
        for value in (energy, yield_value, uncertainty):
            value.setflags(write=False)
        object.__setattr__(self, "energy_eV", energy)
        object.__setattr__(self, "yield_sio2_per_ion", yield_value)
        object.__setattr__(self, "digitization_yield_uncertainty", uncertainty)


class Karahashi2007ReactiveIonYieldTable:
    """Piecewise-linear normal-incidence closure over measured Figure 4 points.

    This object is fitted evidence: predictions at the tabulated points are
    reproductions, not validation.  It is suitable for direct-beam regression
    and bounded interpolation studies.  It is not a substitute for a reactor's
    species-resolved ion spectrum.
    """

    def __init__(self, series, *, source_table_sha256: str,
                 cosine_tolerance: float = 1e-5):
        by_name = {}
        for item in series:
            if not isinstance(item, ReactiveIonYieldSeries):
                raise TypeError("series must contain ReactiveIonYieldSeries")
            if item.species in by_name:
                raise ValueError(f"duplicate reactive-ion series {item.species}")
            by_name[item.species] = item
        if set(by_name) != {"F+", "CF+", "CF2+", "CF3+"}:
            raise ValueError("Karahashi table requires F+, CF+, CF2+, and CF3+")
        if (len(source_table_sha256) != 64
                or any(char not in "0123456789abcdef"
                       for char in source_table_sha256)
                or not np.isfinite(cosine_tolerance)
                or cosine_tolerance < 0.0):
            raise ValueError("invalid reactive-ion table provenance or tolerance")
        self.series = MappingProxyType(by_name)
        self.source_table_sha256 = source_table_sha256
        self.cosine_tolerance = float(cosine_tolerance)
        digest = hashlib.sha256()
        digest.update(b"karahashi-2007-normal-incidence-linear-v1")
        digest.update(source_table_sha256.encode("ascii"))
        digest.update(np.float64(self.cosine_tolerance).tobytes())
        for name in sorted(by_name):
            item = by_name[name]
            digest.update(name.encode("ascii"))
            digest.update(item.energy_eV.tobytes())
            digest.update(item.yield_sio2_per_ion.tobytes())
            digest.update(item.digitization_yield_uncertainty.tobytes())
        self.fingerprint = digest.hexdigest()

    @classmethod
    def from_observations(cls, observations, *, source_table_sha256: str,
                          cosine_tolerance: float = 1e-5):
        grouped = {}
        for row in observations:
            grouped.setdefault(row.species, []).append(row)
        series = []
        for species, rows in sorted(grouped.items()):
            rows = sorted(rows, key=lambda row: row.energy_eV)
            series.append(ReactiveIonYieldSeries(
                species=species,
                energy_eV=np.array([row.energy_eV for row in rows]),
                yield_sio2_per_ion=np.array(
                    [row.yield_sio2_per_ion for row in rows]),
                digitization_yield_uncertainty=np.array(
                    [row.digitization_yield_uncertainty for row in rows]),
            ))
        return cls(
            series, source_table_sha256=source_table_sha256,
            cosine_tolerance=cosine_tolerance)

    @property
    def supported_species(self):
        return tuple(self.series)

    def energy_domain_eV(self, species):
        item = self.series.get(species)
        if item is None:
            return None
        return float(item.energy_eV[0]), float(item.energy_eV[-1])

    def supports(self, species, energy_eV, cosine_incidence):
        item = self.series.get(species)
        energy = np.asarray(energy_eV, dtype=float)
        cosine = np.asarray(cosine_incidence, dtype=float)
        energy, cosine = np.broadcast_arrays(energy, cosine)
        if item is None:
            return np.zeros(energy.shape, dtype=bool)
        return (
            (energy >= item.energy_eV[0])
            & (energy <= item.energy_eV[-1])
            & (np.abs(cosine - 1.0) <= self.cosine_tolerance)
        )

    def evaluate(self, species, energy_eV, cosine_incidence):
        item = self.series.get(species)
        if item is None:
            raise ValueError(f"no Karahashi reactive-ion series for {species}")
        energy = np.asarray(energy_eV, dtype=float)
        cosine = np.asarray(cosine_incidence, dtype=float)
        energy, cosine = np.broadcast_arrays(energy, cosine)
        if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
            raise ValueError("reactive-ion energies must be finite and nonnegative")
        if (np.any(~np.isfinite(cosine))
                or np.any((cosine < 0.0) | (cosine > 1.0))):
            raise ValueError("reactive-ion incidence cosines must lie in [0, 1]")
        supported = self.supports(species, energy, cosine)
        if not np.all(supported):
            lo, hi = self.energy_domain_eV(species)
            raise ValueError(
                f"{species} leaves Karahashi's normal-incidence "
                f"{lo:g}--{hi:g} eV measured support")
        return np.interp(
            energy, item.energy_eV, item.yield_sio2_per_ion)

    def evaluate_uncertainty(self, species, energy_eV, cosine_incidence):
        item = self.series.get(species)
        # Reuse the same support checks as the central value.
        energy, cosine = np.broadcast_arrays(
            np.asarray(energy_eV, dtype=float),
            np.asarray(cosine_incidence, dtype=float))
        self.evaluate(species, energy, cosine)
        return np.interp(
            energy,
            item.energy_eV,
            item.digitization_yield_uncertainty)
