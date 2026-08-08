"""Deterministic electron-energy finite volumes and linear EEPF moments.

The isotropic electron energy probability function (EEPF) is normalized as

``integral sqrt(epsilon) f0(epsilon) d epsilon = 1``.

This is the convention used by Hagelaar and Pitchford (2005).  Values stored
here are piecewise constant over fixed energy cells, so normalization and mean
energy weights are exact for that finite-volume representation.  Collision
cross sections are piecewise linear on their source knots and their rate and
incident-energy moments are integrated exactly inside each energy cell.

The Scharfetter--Gummel class is the conservative reference operator on which
the collision Boltzmann solve will be assembled.  It is not by itself a swarm,
reactor, wafer-boundary, or feature-depth model; those evidence boundaries are
exposed on every returned solution.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy import sparse
from scipy.integrate import simpson
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from .electron_collision_deck import ElectronCollisionProcess
from .network import E_CHARGE_C, ELECTRON_MASS_KG


ELECTRON_SPEED_PER_SQRT_EV_M_S = math.sqrt(
    2.0 * E_CHARGE_C / ELECTRON_MASS_KG)
_MOMENTUM_KINDS = frozenset({"EFFECTIVE", "ELASTIC", "MOMENTUM"})


def _readonly_float_array(values, *, name: str) -> np.ndarray:
    array = np.array(values, dtype=float, copy=True)
    if array.size == 0 or np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must be finite and non-empty")
    array.setflags(write=False)
    return array


def _scalar_or_array(value: np.ndarray | float):
    array = np.asarray(value, dtype=float)
    return float(array) if array.ndim == 0 else array


def _collision_support_nodes(
    process: ElectronCollisionProcess,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return interpolation nodes and the lower resolved-population bound."""
    source_energy = np.asarray(process.electron_energy_eV, dtype=float)
    source_sigma = np.asarray(process.cross_section_m2, dtype=float)
    if process.kind in _MOMENTUM_KINDS:
        return source_energy, source_sigma, float(source_energy[0])
    threshold = max(0.0, float(process.energy_loss_eV or 0.0))
    first = int(np.searchsorted(source_energy, threshold, side="left"))
    if first < source_energy.size and source_energy[first] == threshold:
        energy = source_energy[first:]
        sigma = source_sigma[first:]
    else:
        energy = np.concatenate((np.array([threshold]), source_energy[first:]))
        sigma = np.concatenate((np.array([0.0]), source_sigma[first:]))
    if energy.size < 2:
        raise ValueError("collision process has insufficient resolved support")
    return energy, sigma, 0.0


def _integrate_energy_cross_section(
    energy: np.ndarray,
    cross_section: np.ndarray,
    lower_eV: float,
    upper_eV: float,
    *,
    energy_power: int,
) -> float:
    """Integrate ``E**energy_power * sigma(E)`` on source-linear knots."""
    return _PiecewiseLinearCrossSectionIntegral(
        energy, cross_section, powers=(energy_power,)
    ).integrate(lower_eV, upper_eV, energy_power=energy_power)


class _PiecewiseLinearCrossSectionIntegral:
    """Exact O(log n) antiderivatives for a source-linear cross section."""

    def __init__(self, energy, cross_section, *, powers=(1, 2)):
        self.energy = np.asarray(energy, dtype=float)
        self.cross_section = np.asarray(cross_section, dtype=float)
        self.powers = tuple(sorted({int(value) for value in powers}))
        if (
            self.energy.ndim != 1
            or self.cross_section.shape != self.energy.shape
            or self.energy.size < 2
            or np.any(~np.isfinite(self.energy))
            or np.any(~np.isfinite(self.cross_section))
            or np.any(np.diff(self.energy) <= 0.0)
            or np.any(self.cross_section < 0.0)
            or not self.powers
            or any(value < 0 for value in self.powers)
        ):
            raise ValueError("invalid piecewise-linear cross-section integral")
        energy_extended = self.energy.astype(np.longdouble)
        cross_section_extended = self.cross_section.astype(np.longdouble)
        self.energy_extended = energy_extended
        self.cross_section_extended = cross_section_extended
        self.slope = (
            np.diff(cross_section_extended) / np.diff(energy_extended))
        self.cumulative = {}
        for power in self.powers:
            segment = self._segment_integral(
                self.energy[:-1], self.energy[1:], power)
            self.cumulative[power] = np.concatenate((
                np.array([0.0], dtype=np.longdouble),
                np.cumsum(segment, dtype=np.longdouble),
            ))

    def _segment_integral(self, lower, upper, power):
        """Integrate each native segment in a local energy coordinate.

        A global ``m E + b`` primitive catastrophically cancels when a tiny
        cross section changes over a narrow interval at high energy. Expanding
        ``E**p`` around each source knot keeps the same exact polynomial while
        evaluating only non-cancelling local powers.
        """
        lower = np.asarray(lower, dtype=np.longdouble)
        upper = np.asarray(upper, dtype=np.longdouble)
        base = self.energy_extended[:-1]
        t_lower = lower - base
        t_upper = upper - base
        result = np.zeros_like(base)
        for exponent in range(power + 1):
            coefficient = (
                math.comb(power, exponent)
                * base ** (power - exponent)
            )
            result += coefficient * (
                self.cross_section_extended[:-1]
                * (t_upper ** (exponent + 1)
                   - t_lower ** (exponent + 1))
                / (exponent + 1)
                + self.slope
                * (t_upper ** (exponent + 2)
                   - t_lower ** (exponent + 2))
                / (exponent + 2)
            )
        return result

    def _partial_segment_integral(
        self, segment: int, upper: np.longdouble, power: int
    ) -> np.longdouble:
        base = self.energy_extended[segment]
        span = upper - base
        result = np.longdouble(0.0)
        for exponent in range(power + 1):
            coefficient = (
                math.comb(power, exponent)
                * base ** (power - exponent)
            )
            result += coefficient * (
                self.cross_section_extended[segment]
                * span ** (exponent + 1) / (exponent + 1)
                + self.slope[segment]
                * span ** (exponent + 2) / (exponent + 2)
            )
        return result

    def _primitive(self, value: float, power: int) -> np.longdouble:
        x = min(max(float(value), self.energy[0]), self.energy[-1])
        if x <= self.energy[0]:
            return np.longdouble(0.0)
        if x >= self.energy[-1]:
            return self.cumulative[power][-1]
        segment = int(np.searchsorted(self.energy, x, side="right") - 1)
        x = np.longdouble(x)
        partial = self._partial_segment_integral(segment, x, power)
        return self.cumulative[power][segment] + partial

    def integrate(
        self,
        lower_eV: float,
        upper_eV: float,
        *,
        energy_power: int,
    ) -> float:
        power = int(energy_power)
        if power not in self.cumulative or power != energy_power:
            raise ValueError("energy power was not compiled")
        lower = max(float(lower_eV), float(self.energy[0]))
        upper = min(float(upper_eV), float(self.energy[-1]))
        if upper <= lower:
            return 0.0
        upper_primitive = self._primitive(upper, power)
        lower_primitive = self._primitive(lower, power)
        value = upper_primitive - lower_primitive
        # Cumulative subtraction can leave a signed ulp for vanishing tail
        # intervals. A physically negative resolved cross-section integral is
        # still rejected at a scale larger than roundoff.
        scale = max(
            abs(upper_primitive),
            abs(lower_primitive),
            np.finfo(np.longdouble).tiny,
        )
        if value < -32.0 * np.finfo(np.longdouble).eps * scale:
            raise FloatingPointError("negative cross-section antiderivative")
        return float(max(value, 0.0))


@dataclass(frozen=True)
class ElectronEnergyGrid:
    """A fixed-topology electron-energy finite-volume grid."""

    boundaries_eV: tuple[float, ...]

    def __post_init__(self):
        boundaries = tuple(float(value) for value in self.boundaries_eV)
        if (
            len(boundaries) < 3
            or any(not math.isfinite(value) for value in boundaries)
            or boundaries[0] < 0.0
            or any(
                right <= left
                for left, right in zip(boundaries, boundaries[1:])
            )
        ):
            raise ValueError("invalid electron-energy grid")
        object.__setattr__(self, "boundaries_eV", boundaries)

    @classmethod
    def linear(
        cls,
        maximum_energy_eV: float,
        cell_count: int,
        *,
        minimum_energy_eV: float = 0.0,
    ) -> "ElectronEnergyGrid":
        maximum = float(maximum_energy_eV)
        minimum = float(minimum_energy_eV)
        if (
            not math.isfinite(maximum)
            or not math.isfinite(minimum)
            or minimum < 0.0
            or maximum <= minimum
            or int(cell_count) != cell_count
            or cell_count < 2
        ):
            raise ValueError("invalid linear electron-energy grid")
        return cls(tuple(np.linspace(minimum, maximum, int(cell_count) + 1)))

    @classmethod
    def piecewise_linear(
        cls,
        breakpoints_eV: tuple[float, ...],
        cells_per_segment: tuple[int, ...],
        *,
        inserted_boundaries_eV: tuple[float, ...] = (),
    ) -> "ElectronEnergyGrid":
        """Build a fixed nonuniform grid with exact physical boundaries.

        Each consecutive breakpoint pair is divided linearly by its declared
        cell count. Additional boundaries, such as collision thresholds, are
        inserted exactly without perturbing the surrounding fixed topology.
        """
        breakpoints = tuple(float(value) for value in breakpoints_eV)
        raw_counts = tuple(cells_per_segment)
        counts = tuple(int(value) for value in raw_counts)
        inserted = tuple(float(value) for value in inserted_boundaries_eV)
        if (
            len(breakpoints) < 2
            or len(counts) != len(breakpoints) - 1
            or any(not math.isfinite(value) for value in breakpoints)
            or breakpoints[0] < 0.0
            or any(
                right <= left
                for left, right in zip(breakpoints, breakpoints[1:])
            )
            or any(value < 1 for value in counts)
            or any(raw != converted for raw, converted in zip(
                raw_counts, counts))
            or any(not math.isfinite(value) for value in inserted)
            or any(
                value <= breakpoints[0] or value >= breakpoints[-1]
                for value in inserted
            )
        ):
            raise ValueError("invalid piecewise-linear electron-energy grid")
        pieces = [
            np.linspace(left, right, count + 1)
            for left, right, count in zip(
                breakpoints[:-1], breakpoints[1:], counts)
        ]
        boundaries = np.unique(np.concatenate((
            *pieces,
            np.asarray(inserted, dtype=float),
        )))
        return cls(tuple(boundaries))

    @property
    def cell_count(self) -> int:
        return len(self.boundaries_eV) - 1

    @property
    def boundaries(self) -> np.ndarray:
        return np.asarray(self.boundaries_eV, dtype=float)

    @property
    def cell_centers_eV(self) -> np.ndarray:
        boundaries = self.boundaries
        return 0.5 * (boundaries[:-1] + boundaries[1:])

    @property
    def cell_widths_eV(self) -> np.ndarray:
        return np.diff(self.boundaries)

    @property
    def normalization_weights(self) -> np.ndarray:
        """Exact cell weights for ``integral sqrt(E) f0 dE``."""
        boundaries = self.boundaries
        return (2.0 / 3.0) * (
            boundaries[1:] ** 1.5 - boundaries[:-1] ** 1.5
        )

    @property
    def mean_energy_weights_eV(self) -> np.ndarray:
        """Exact cell weights for ``integral E^(3/2) f0 dE``."""
        boundaries = self.boundaries
        return (2.0 / 5.0) * (
            boundaries[1:] ** 2.5 - boundaries[:-1] ** 2.5
        )


@dataclass(frozen=True)
class ElectronEnergyDistribution:
    """One or a batch of normalized piecewise-constant EEPFs."""

    grid: ElectronEnergyGrid
    eepf_eV_minus_3_over_2: np.ndarray
    normalization_tolerance: float = 5.0e-12

    def __post_init__(self):
        if not isinstance(self.grid, ElectronEnergyGrid):
            raise TypeError("an ElectronEnergyGrid is required")
        values = _readonly_float_array(
            self.eepf_eV_minus_3_over_2, name="EEPF")
        if values.ndim < 1 or values.shape[-1] != self.grid.cell_count:
            raise ValueError("EEPF shape does not match energy grid")
        if np.any(values < 0.0):
            raise ValueError("EEPF must be nonnegative")
        tolerance = float(self.normalization_tolerance)
        if not math.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("invalid EEPF normalization tolerance")
        normalization = np.sum(
            values * self.grid.normalization_weights, axis=-1)
        if np.any(np.abs(normalization - 1.0) > tolerance):
            raise ValueError("EEPF is not normalized on the declared grid")
        object.__setattr__(self, "eepf_eV_minus_3_over_2", values)
        object.__setattr__(self, "normalization_tolerance", tolerance)

    @classmethod
    def from_unnormalized(
        cls,
        grid: ElectronEnergyGrid,
        values,
        *,
        normalization_tolerance: float = 5.0e-12,
    ) -> "ElectronEnergyDistribution":
        normalized = normalize_eepf(grid, values)
        return cls(grid, normalized, normalization_tolerance)

    @classmethod
    def maxwellian(
        cls,
        grid: ElectronEnergyGrid,
        electron_temperature_eV,
    ) -> "ElectronEnergyDistribution":
        temperatures = np.asarray(electron_temperature_eV, dtype=float)
        if (
            np.any(~np.isfinite(temperatures))
            or np.any(temperatures <= 0.0)
        ):
            raise ValueError("electron temperature must be positive")
        centers = grid.cell_centers_eV
        raw = np.exp(-centers / temperatures[..., np.newaxis])
        return cls.from_unnormalized(grid, raw)

    @property
    def batch_shape(self) -> tuple[int, ...]:
        return self.eepf_eV_minus_3_over_2.shape[:-1]

    @property
    def normalization(self):
        return _scalar_or_array(np.sum(
            self.eepf_eV_minus_3_over_2
            * self.grid.normalization_weights,
            axis=-1,
        ))

    @property
    def mean_energy_eV(self):
        return _scalar_or_array(np.sum(
            self.eepf_eV_minus_3_over_2
            * self.grid.mean_energy_weights_eV,
            axis=-1,
        ))

    def population_fraction_in_last_cells(self, cell_count: int = 1):
        if int(cell_count) != cell_count or not 1 <= cell_count <= self.grid.cell_count:
            raise ValueError("invalid tail cell count")
        weights = self.grid.normalization_weights[-int(cell_count):]
        value = np.sum(
            self.eepf_eV_minus_3_over_2[..., -int(cell_count):] * weights,
            axis=-1,
        )
        return _scalar_or_array(value)

    def convergence_receipt(
        self,
        *,
        tail_cell_count: int = 4,
        maximum_tail_population_fraction: float = 1.0e-8,
    ) -> dict[str, object]:
        maximum_tail = float(maximum_tail_population_fraction)
        if not math.isfinite(maximum_tail) or not 0.0 < maximum_tail < 1.0:
            raise ValueError("invalid tail population tolerance")
        tail = np.asarray(
            self.population_fraction_in_last_cells(tail_cell_count))
        return {
            "normalization": self.normalization,
            "mean_energy_eV": self.mean_energy_eV,
            "tail_cell_count": int(tail_cell_count),
            "tail_population_fraction": _scalar_or_array(tail),
            "maximum_tail_population_fraction": maximum_tail,
            "normalization_passed": bool(np.all(
                np.abs(np.asarray(self.normalization) - 1.0)
                <= self.normalization_tolerance
            )),
            "positivity_passed": bool(np.all(
                self.eepf_eV_minus_3_over_2 >= 0.0)),
            "tail_passed": bool(np.all(tail <= maximum_tail)),
            "supports_swarm_validation": False,
            "supports_reactor_state_prediction": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
        }


def normalize_eepf(
    grid: ElectronEnergyGrid,
    unnormalized_values,
) -> np.ndarray:
    """Normalize a batch-first EEPF without changing its fixed topology."""
    if not isinstance(grid, ElectronEnergyGrid):
        raise TypeError("an ElectronEnergyGrid is required")
    values = np.asarray(unnormalized_values, dtype=float)
    if (
        values.ndim < 1
        or values.shape[-1] != grid.cell_count
        or np.any(~np.isfinite(values))
        or np.any(values < 0.0)
    ):
        raise ValueError("invalid unnormalized EEPF")
    normalization = np.sum(
        values * grid.normalization_weights, axis=-1, keepdims=True)
    if np.any(normalization <= 0.0):
        raise ValueError("EEPF normalization must be positive")
    return values / normalization


def normalize_eepf_jvp(
    grid: ElectronEnergyGrid,
    unnormalized_values,
    unnormalized_tangent,
) -> np.ndarray:
    """Analytic JVP of :func:`normalize_eepf`."""
    values = np.asarray(unnormalized_values, dtype=float)
    tangent = np.asarray(unnormalized_tangent, dtype=float)
    if values.shape != tangent.shape:
        raise ValueError("EEPF value and tangent shapes must match")
    normalized = normalize_eepf(grid, values)
    normalization = np.sum(
        values * grid.normalization_weights, axis=-1, keepdims=True)
    normalization_tangent = np.sum(
        tangent * grid.normalization_weights, axis=-1, keepdims=True)
    return (
        tangent / normalization
        - normalized * normalization_tangent / normalization
    )


def normalize_eepf_vjp(
    grid: ElectronEnergyGrid,
    unnormalized_values,
    normalized_cotangent,
) -> np.ndarray:
    """Analytic VJP of :func:`normalize_eepf`."""
    values = np.asarray(unnormalized_values, dtype=float)
    cotangent = np.asarray(normalized_cotangent, dtype=float)
    if values.shape != cotangent.shape:
        raise ValueError("EEPF value and cotangent shapes must match")
    normalized = normalize_eepf(grid, values)
    normalization = np.sum(
        values * grid.normalization_weights, axis=-1, keepdims=True)
    projection = np.sum(
        cotangent * normalized, axis=-1, keepdims=True)
    return (
        cotangent
        - grid.normalization_weights * projection
    ) / normalization


@dataclass(frozen=True)
class ElectronCollisionMoments:
    """Rate and incident-energy moments for one collision process."""

    rate_coefficient_m3_s: float | np.ndarray
    incident_energy_moment_eV_m3_s: float | np.ndarray
    collision_weighted_mean_incident_energy_eV: float | np.ndarray
    unresolved_population_fraction: float | np.ndarray
    process_kind: str
    target: str
    product: str | None
    supports_swarm_validation: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


@dataclass(frozen=True)
class ElectronCollisionMomentKernel:
    """Exact linear moments for a piecewise-linear collision cross section."""

    grid: ElectronEnergyGrid
    process: ElectronCollisionProcess
    rate_weights_m3_s: tuple[float, ...]
    incident_energy_weights_eV_m3_s: tuple[float, ...]
    unresolved_population_weights: tuple[float, ...]
    resolved_energy_interval_eV: tuple[float, float]

    def __post_init__(self):
        if not isinstance(self.grid, ElectronEnergyGrid):
            raise TypeError("an ElectronEnergyGrid is required")
        if not isinstance(self.process, ElectronCollisionProcess):
            raise TypeError("an ElectronCollisionProcess is required")
        arrays = (
            self.rate_weights_m3_s,
            self.incident_energy_weights_eV_m3_s,
            self.unresolved_population_weights,
        )
        if any(len(item) != self.grid.cell_count for item in arrays):
            raise ValueError("collision-kernel shape does not match grid")
        if any(
            not math.isfinite(value) or value < 0.0
            for item in arrays
            for value in item
        ):
            raise ValueError("invalid collision-kernel weight")
        lower, upper = self.resolved_energy_interval_eV
        if not 0.0 <= lower < upper:
            raise ValueError("invalid resolved collision support")

    @classmethod
    def from_process(
        cls,
        grid: ElectronEnergyGrid,
        process: ElectronCollisionProcess,
    ) -> "ElectronCollisionMomentKernel":
        if not isinstance(grid, ElectronEnergyGrid):
            raise TypeError("an ElectronEnergyGrid is required")
        if not isinstance(process, ElectronCollisionProcess):
            raise TypeError("an ElectronCollisionProcess is required")

        energy, sigma, resolved_lower = _collision_support_nodes(process)
        resolved_upper = float(energy[-1])
        integral = _PiecewiseLinearCrossSectionIntegral(
            energy, sigma, powers=(1, 2))

        boundaries = grid.boundaries
        rate_weights = np.zeros(grid.cell_count)
        energy_weights = np.zeros(grid.cell_count)
        resolved_population_weights = np.zeros(grid.cell_count)

        for cell in range(grid.cell_count):
            cell_lower = boundaries[cell]
            cell_upper = boundaries[cell + 1]
            overlap_lower = max(cell_lower, float(energy[0]))
            overlap_upper = min(cell_upper, resolved_upper)
            if overlap_upper > overlap_lower:
                rate_weights[cell] += ELECTRON_SPEED_PER_SQRT_EV_M_S * (
                    integral.integrate(
                        overlap_lower,
                        overlap_upper,
                        energy_power=1,
                    )
                )
                energy_weights[cell] += ELECTRON_SPEED_PER_SQRT_EV_M_S * (
                    integral.integrate(
                        overlap_lower,
                        overlap_upper,
                        energy_power=2,
                    )
                )

            population_lower = max(cell_lower, resolved_lower)
            population_upper = min(cell_upper, resolved_upper)
            if population_upper > population_lower:
                resolved_population_weights[cell] = (2.0 / 3.0) * (
                    population_upper ** 1.5 - population_lower ** 1.5
                )

        unresolved = np.maximum(
            grid.normalization_weights - resolved_population_weights, 0.0)
        if np.any(rate_weights < -1.0e-30) or np.any(energy_weights < -1.0e-30):
            raise FloatingPointError("negative integrated collision weight")
        rate_weights = np.maximum(rate_weights, 0.0)
        energy_weights = np.maximum(energy_weights, 0.0)
        return cls(
            grid=grid,
            process=process,
            rate_weights_m3_s=tuple(rate_weights),
            incident_energy_weights_eV_m3_s=tuple(energy_weights),
            unresolved_population_weights=tuple(unresolved),
            resolved_energy_interval_eV=(resolved_lower, resolved_upper),
        )

    def _validated_values(
        self,
        distribution: ElectronEnergyDistribution,
    ) -> np.ndarray:
        if not isinstance(distribution, ElectronEnergyDistribution):
            raise TypeError("an ElectronEnergyDistribution is required")
        if distribution.grid != self.grid:
            raise ValueError("collision kernel and EEPF grids differ")
        return distribution.eepf_eV_minus_3_over_2

    def evaluate(
        self,
        distribution: ElectronEnergyDistribution,
        *,
        maximum_unresolved_population_fraction: float = 1.0e-8,
    ) -> ElectronCollisionMoments:
        values = self._validated_values(distribution)
        maximum_unresolved = float(maximum_unresolved_population_fraction)
        if (
            not math.isfinite(maximum_unresolved)
            or not 0.0 <= maximum_unresolved < 1.0
        ):
            raise ValueError("invalid unresolved-population tolerance")
        unresolved = np.sum(
            values * np.asarray(self.unresolved_population_weights), axis=-1)
        if np.any(unresolved > maximum_unresolved):
            raise ValueError(
                "EEPF population outside collision cross-section support "
                "exceeds the declared tolerance"
            )
        rate = np.sum(
            values * np.asarray(self.rate_weights_m3_s), axis=-1)
        incident = np.sum(
            values
            * np.asarray(self.incident_energy_weights_eV_m3_s),
            axis=-1,
        )
        mean = np.divide(
            incident,
            rate,
            out=np.full_like(np.asarray(rate), np.nan, dtype=float),
            where=np.asarray(rate) > 0.0,
        )
        return ElectronCollisionMoments(
            rate_coefficient_m3_s=_scalar_or_array(rate),
            incident_energy_moment_eV_m3_s=_scalar_or_array(incident),
            collision_weighted_mean_incident_energy_eV=(
                _scalar_or_array(mean)),
            unresolved_population_fraction=_scalar_or_array(unresolved),
            process_kind=self.process.kind,
            target=self.process.target,
            product=self.process.product,
        )

    def rate_jvp(self, eepf_tangent) -> float | np.ndarray:
        tangent = np.asarray(eepf_tangent, dtype=float)
        if tangent.shape[-1:] != (self.grid.cell_count,):
            raise ValueError("EEPF tangent shape does not match grid")
        return _scalar_or_array(np.sum(
            tangent * np.asarray(self.rate_weights_m3_s), axis=-1))

    def rate_vjp(self, rate_cotangent) -> np.ndarray:
        cotangent = np.asarray(rate_cotangent, dtype=float)
        return (
            cotangent[..., np.newaxis]
            * np.asarray(self.rate_weights_m3_s)
        )

    def incident_energy_jvp(self, eepf_tangent) -> float | np.ndarray:
        tangent = np.asarray(eepf_tangent, dtype=float)
        if tangent.shape[-1:] != (self.grid.cell_count,):
            raise ValueError("EEPF tangent shape does not match grid")
        return _scalar_or_array(np.sum(
            tangent
            * np.asarray(self.incident_energy_weights_eV_m3_s),
            axis=-1,
        ))

    def incident_energy_vjp(self, moment_cotangent) -> np.ndarray:
        cotangent = np.asarray(moment_cotangent, dtype=float)
        return (
            cotangent[..., np.newaxis]
            * np.asarray(self.incident_energy_weights_eV_m3_s)
        )


def _bernoulli_function(value: np.ndarray) -> np.ndarray:
    """Stable Bernoulli function ``x / (exp(x) - 1)``."""
    x = np.asarray(value, dtype=float)
    result = np.empty_like(x)
    small = np.abs(x) < 1.0e-5
    positive = x > 50.0
    negative = x < -50.0
    regular = ~(small | positive | negative)
    x_small = x[small]
    result[small] = (
        1.0 - 0.5 * x_small + x_small ** 2 / 12.0
        - x_small ** 4 / 720.0
    )
    result[positive] = x[positive] * np.exp(-x[positive]) / (
        1.0 - np.exp(-x[positive]))
    result[negative] = -x[negative] / (1.0 - np.exp(x[negative]))
    result[regular] = x[regular] / np.expm1(x[regular])
    return result


@dataclass(frozen=True)
class ScharfetterGummelSolution:
    """A normalized conservative energy-space solution and its ledgers."""

    distribution: ElectronEnergyDistribution
    compatibility_multiplier: float | np.ndarray
    energy_flux_faces: np.ndarray
    maximum_augmented_residual: float
    maximum_physical_residual: float
    source_integral_sum: float | np.ndarray
    supports_collision_boltzmann_solve: bool = False
    supports_swarm_validation: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


@dataclass(frozen=True)
class ScharfetterGummelEnergyOperator:
    """Batch-first conservative reference operator with zero-flux boundaries.

    For a face with locally constant drift ``W`` and diffusion ``D``, the
    exact face flux between adjacent cell centers is used.  A bordered
    normalization equation removes the zero-flux nullspace.  The returned
    compatibility multiplier is zero only when the supplied cell-integrated
    source is compatible with the closed energy domain.
    """

    grid: ElectronEnergyGrid
    drift_eV_s: np.ndarray
    diffusion_eV2_s: np.ndarray

    def __post_init__(self):
        if not isinstance(self.grid, ElectronEnergyGrid):
            raise TypeError("an ElectronEnergyGrid is required")
        drift = _readonly_float_array(self.drift_eV_s, name="energy drift")
        diffusion = _readonly_float_array(
            self.diffusion_eV2_s, name="energy diffusion")
        expected = self.grid.cell_count - 1
        if (
            drift.ndim < 1
            or diffusion.shape != drift.shape
            or drift.shape[-1] != expected
            or np.any(diffusion <= 0.0)
        ):
            raise ValueError("invalid Scharfetter-Gummel face coefficients")
        object.__setattr__(self, "drift_eV_s", drift)
        object.__setattr__(self, "diffusion_eV2_s", diffusion)

    @property
    def batch_shape(self) -> tuple[int, ...]:
        return self.drift_eV_s.shape[:-1]

    def _face_coefficients(
        self,
        flat_batch_index: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        centers = self.grid.cell_centers_eV
        spacing = np.diff(centers)
        drift = self.drift_eV_s.reshape(
            (-1, self.grid.cell_count - 1))[flat_batch_index]
        diffusion = self.diffusion_eV2_s.reshape(
            (-1, self.grid.cell_count - 1))[flat_batch_index]
        peclet = drift * spacing / diffusion
        left = diffusion / spacing * _bernoulli_function(-peclet)
        right = -diffusion / spacing * _bernoulli_function(peclet)
        return left, right

    def _matrix(self, flat_batch_index: int) -> sparse.csr_matrix:
        cell_count = self.grid.cell_count
        left, right = self._face_coefficients(flat_batch_index)
        diagonal = np.empty(cell_count)
        diagonal[0] = left[0]
        diagonal[1:-1] = left[1:] - right[:-1]
        diagonal[-1] = -right[-1]
        return sparse.diags(
            diagonals=(-left, diagonal, right),
            offsets=(-1, 0, 1),
            shape=(cell_count, cell_count),
            format="csr",
        )

    def _augmented_matrix(self, flat_batch_index: int) -> sparse.csr_matrix:
        matrix = self._matrix(flat_batch_index)
        compatibility_direction = sparse.csr_matrix(
            self.grid.cell_widths_eV[:, np.newaxis])
        normalization = sparse.csr_matrix(
            self.grid.normalization_weights[np.newaxis, :])
        return sparse.bmat(
            (
                (matrix, compatibility_direction),
                (normalization, sparse.csr_matrix((1, 1))),
            ),
            format="csr",
        )

    def _broadcast_cells(self, values, *, name: str) -> np.ndarray:
        array = np.asarray(values, dtype=float)
        target_shape = self.batch_shape + (self.grid.cell_count,)
        try:
            output = np.broadcast_to(array, target_shape)
        except ValueError as exc:
            raise ValueError(f"{name} cannot broadcast to operator batch") from exc
        if np.any(~np.isfinite(output)):
            raise ValueError(f"{name} must be finite")
        return np.asarray(output)

    def solve(
        self,
        source_integrals,
        *,
        maximum_compatibility_multiplier: float = 1.0e-10,
    ) -> ScharfetterGummelSolution:
        source = self._broadcast_cells(
            source_integrals, name="energy-space source")
        maximum_compatibility = float(maximum_compatibility_multiplier)
        if not math.isfinite(maximum_compatibility) or maximum_compatibility < 0.0:
            raise ValueError("invalid compatibility tolerance")
        flat_source = source.reshape((-1, self.grid.cell_count))
        flat_distribution = np.empty_like(flat_source)
        multipliers = np.empty(flat_source.shape[0])
        augmented_residuals = []
        physical_residuals = []
        fluxes = np.zeros((flat_source.shape[0], self.grid.cell_count + 1))

        for batch_index, source_row in enumerate(flat_source):
            augmented = self._augmented_matrix(batch_index)
            rhs = np.concatenate((source_row, np.array([1.0])))
            state = np.asarray(spsolve(augmented, rhs), dtype=float)
            if np.any(~np.isfinite(state)):
                raise FloatingPointError("nonfinite energy-space solution")
            distribution = state[:-1]
            multiplier = float(state[-1])
            minimum = float(np.min(distribution))
            if minimum < 0.0:
                raise FloatingPointError("negative Scharfetter-Gummel solution")
            flat_distribution[batch_index] = distribution
            multipliers[batch_index] = multiplier

            left, right = self._face_coefficients(batch_index)
            fluxes[batch_index, 1:-1] = (
                left * distribution[:-1] + right * distribution[1:])
            matrix = self._matrix(batch_index)
            compatibility_direction = self.grid.cell_widths_eV
            augmented_residuals.append(float(np.max(np.abs(
                matrix @ distribution
                + compatibility_direction * multiplier
                - source_row
            ))))
            physical_residuals.append(float(np.max(np.abs(
                matrix @ distribution - source_row))))

        if np.any(np.abs(multipliers) > maximum_compatibility):
            raise ValueError(
                "energy-space source is incompatible with zero-flux boundaries"
            )
        distribution_shape = self.batch_shape + (self.grid.cell_count,)
        flux_shape = self.batch_shape + (self.grid.cell_count + 1,)
        distribution = ElectronEnergyDistribution(
            self.grid, flat_distribution.reshape(distribution_shape))
        output_multipliers = multipliers.reshape(self.batch_shape)
        source_sum = np.sum(source, axis=-1)
        return ScharfetterGummelSolution(
            distribution=distribution,
            compatibility_multiplier=_scalar_or_array(output_multipliers),
            energy_flux_faces=fluxes.reshape(flux_shape),
            maximum_augmented_residual=max(augmented_residuals),
            maximum_physical_residual=max(physical_residuals),
            source_integral_sum=_scalar_or_array(source_sum),
        )

    def implicit_source_jvp(self, source_tangent):
        """Exact implicit JVP of normalized state with respect to source."""
        tangent = self._broadcast_cells(
            source_tangent, name="energy-space source tangent")
        flat = tangent.reshape((-1, self.grid.cell_count))
        state_tangent = np.empty_like(flat)
        multiplier_tangent = np.empty(flat.shape[0])
        for batch_index, tangent_row in enumerate(flat):
            rhs = np.concatenate((tangent_row, np.array([0.0])))
            state = np.asarray(spsolve(
                self._augmented_matrix(batch_index), rhs), dtype=float)
            state_tangent[batch_index] = state[:-1]
            multiplier_tangent[batch_index] = state[-1]
        return (
            state_tangent.reshape(
                self.batch_shape + (self.grid.cell_count,)),
            _scalar_or_array(multiplier_tangent.reshape(self.batch_shape)),
        )

    def implicit_source_vjp(
        self,
        distribution_cotangent,
        compatibility_cotangent=0.0,
    ) -> np.ndarray:
        """Exact implicit VJP of normalized state with respect to source."""
        cotangent = self._broadcast_cells(
            distribution_cotangent, name="distribution cotangent")
        multiplier_cotangent = np.asarray(
            compatibility_cotangent, dtype=float)
        try:
            multiplier_cotangent = np.broadcast_to(
                multiplier_cotangent, self.batch_shape)
        except ValueError as exc:
            raise ValueError(
                "compatibility cotangent cannot broadcast to operator batch"
            ) from exc
        flat = cotangent.reshape((-1, self.grid.cell_count))
        flat_multiplier = multiplier_cotangent.reshape(-1)
        source_cotangent = np.empty_like(flat)
        for batch_index, cotangent_row in enumerate(flat):
            rhs = np.concatenate((
                cotangent_row,
                np.array([flat_multiplier[batch_index]]),
            ))
            adjoint = np.asarray(spsolve(
                self._augmented_matrix(batch_index).T, rhs), dtype=float)
            source_cotangent[batch_index] = adjoint[:-1]
        return source_cotangent.reshape(
            self.batch_shape + (self.grid.cell_count,))


K_BOLTZMANN_J_K = 1.380649e-23


@dataclass(frozen=True)
class TwoTermBoltzmannCondition:
    """One local-field gas state for the deterministic two-term solve.

    ``angular_field_frequency_over_density_m3_s`` is omega/N. Zero selects
    the DC equation. A positive value selects the high-frequency two-term
    heating operator with ``E/N`` interpreted as the RMS field.
    """

    reduced_electric_field_Td: float
    gas_temperature_K: float
    target_mole_fractions: Mapping[str, float]
    growth_model: str = "no_growth"
    inelastic_momentum_closure: str = "isotropic_source_reproduction"
    ionization_energy_sharing: str = "equal_sharing"
    initial_electron_temperature_eV: float = 2.0
    angular_field_frequency_over_density_m3_s: float = 0.0
    electron_electron_coulomb_model: str = "none"
    electron_to_neutral_density_ratio: float = 0.0
    gas_number_density_m3: float | None = None

    def __post_init__(self):
        field = float(self.reduced_electric_field_Td)
        temperature = float(self.gas_temperature_K)
        reduced_frequency = float(
            self.angular_field_frequency_over_density_m3_s)
        ionization_degree = float(self.electron_to_neutral_density_ratio)
        gas_density = (
            None if self.gas_number_density_m3 is None
            else float(self.gas_number_density_m3))
        fractions = {
            str(name).strip(): float(value)
            for name, value in dict(self.target_mole_fractions).items()
        }
        if (
            not math.isfinite(field)
            or field < 0.0
            or not math.isfinite(temperature)
            or temperature <= 0.0
            or not math.isfinite(reduced_frequency)
            or reduced_frequency < 0.0
            or not fractions
            or any(not name for name in fractions)
            or any(
                not math.isfinite(value) or value < 0.0
                for value in fractions.values()
            )
            or not math.isclose(
                sum(fractions.values()), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
            or self.growth_model not in {"no_growth", "temporal_growth"}
            or self.inelastic_momentum_closure
            != "isotropic_source_reproduction"
            or self.ionization_energy_sharing != "equal_sharing"
            or self.electron_electron_coulomb_model not in {
                "none", "isotropic_classical_debye"
            }
            or not math.isfinite(ionization_degree)
            or ionization_degree < 0.0
            or (
                self.electron_electron_coulomb_model == "none"
                and (ionization_degree != 0.0 or gas_density is not None)
            )
            or (
                self.electron_electron_coulomb_model
                == "isotropic_classical_debye"
                and (
                    ionization_degree <= 0.0
                    or gas_density is None
                    or not math.isfinite(gas_density)
                    or gas_density <= 0.0
                )
            )
            or not math.isfinite(self.initial_electron_temperature_eV)
            or self.initial_electron_temperature_eV <= 0.0
        ):
            raise ValueError("invalid two-term Boltzmann condition")
        object.__setattr__(self, "reduced_electric_field_Td", field)
        object.__setattr__(self, "gas_temperature_K", temperature)
        object.__setattr__(
            self, "target_mole_fractions", MappingProxyType(fractions))
        object.__setattr__(
            self,
            "initial_electron_temperature_eV",
            float(self.initial_electron_temperature_eV),
        )
        object.__setattr__(
            self,
            "angular_field_frequency_over_density_m3_s",
            reduced_frequency,
        )
        object.__setattr__(
            self, "electron_to_neutral_density_ratio", ionization_degree)
        object.__setattr__(self, "gas_number_density_m3", gas_density)


@dataclass(frozen=True)
class TwoTermElectronTransportMoments:
    """Local two-term flux moments with deliberately bounded evidence."""

    flux_reduced_mobility_m_inv_V_inv_s_inv: float
    scalar_reduced_diffusion_m_inv_s_inv: float
    dissipative_reduced_mobility_m_inv_V_inv_s_inv: float
    reduced_field_power_gain_eV_m3_s: float
    mean_electron_speed_m_s: float
    isotropic_wall_flux_coefficient_m_s: float
    mean_wall_loss_electron_energy_eV: float
    supports_flux_transport_moments: bool = True
    supports_direct_swarm_grade: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        values = np.asarray((
            self.flux_reduced_mobility_m_inv_V_inv_s_inv,
            self.scalar_reduced_diffusion_m_inv_s_inv,
            self.dissipative_reduced_mobility_m_inv_V_inv_s_inv,
            self.reduced_field_power_gain_eV_m3_s,
            self.mean_electron_speed_m_s,
            self.isotropic_wall_flux_coefficient_m_s,
            self.mean_wall_loss_electron_energy_eV,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values[:3] <= 0.0)
            or values[3] < 0.0
            or np.any(values[4:] <= 0.0)
        ):
            raise ValueError("invalid two-term electron transport moments")
        for name, value in zip((
            "flux_reduced_mobility_m_inv_V_inv_s_inv",
            "scalar_reduced_diffusion_m_inv_s_inv",
            "dissipative_reduced_mobility_m_inv_V_inv_s_inv",
            "reduced_field_power_gain_eV_m3_s",
            "mean_electron_speed_m_s",
            "isotropic_wall_flux_coefficient_m_s",
            "mean_wall_loss_electron_energy_eV",
        ), values):
            object.__setattr__(self, name, float(value))


@dataclass(frozen=True)
class TwoTermBoltzmannSolution:
    """Converged local-field two-term EEPF with explicit evidence limits."""

    distribution: ElectronEnergyDistribution
    collision_moments: tuple[ElectronCollisionMoments, ...]
    transport_moments: TwoTermElectronTransportMoments
    reduced_electric_field_Td: float
    gas_temperature_K: float
    angular_field_frequency_over_density_m3_s: float
    growth_model: str
    net_growth_rate_coefficient_m3_s: float
    iteration_count: int
    weighted_iteration_residual: float
    maximum_equation_residual_m3_s: float
    particle_source_from_collision_operator_m3_s: float
    particle_growth_closure_error_m3_s: float
    minimum_raw_eepf_before_roundoff_projection: float
    roundoff_negative_population_fraction: float
    energy_flux_faces_reduced: np.ndarray
    collision_deck_sha256: str
    electron_electron_coulomb_model: str = "none"
    coulomb_logarithm: float | None = None
    nonlinear_iteration_count: int = 0
    nonlinear_weighted_residual: float = 0.0
    growth_root_evaluations: int = 0
    solver_id: str = "petch-two-term-sg-v1"
    supports_collision_boltzmann_solve: bool = True
    supports_direct_swarm_grade: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


class DeterministicTwoTermBoltzmannSolver:
    """Original conservative two-term EEPF solver for local reactor rates.

    The physical coefficients follow Hagelaar--Pitchford equations 7, 27,
    and 39--43 without electron--electron collisions.  Excitation,
    attachment, and equal-sharing ionization are assembled as conservative
    finite-volume scattering-out/in maps.  Temporal growth uses the paper's
    normalization source and effective momentum-transfer correction.

    This solver intentionally stops short of the density-gradient and SST
    hierarchies required to grade mean-arrival-time drift, longitudinal
    diffusion, and effective ionization in attaching chlorine.
    """

    def __init__(self, grid: ElectronEnergyGrid, collision_deck):
        from .electron_collision_deck import ElectronCollisionDeck

        if not isinstance(grid, ElectronEnergyGrid):
            raise TypeError("an ElectronEnergyGrid is required")
        if not isinstance(collision_deck, ElectronCollisionDeck):
            raise TypeError("an ElectronCollisionDeck is required")
        self.grid = grid
        self.collision_deck = collision_deck
        self._collision_kernels = tuple(
            ElectronCollisionMomentKernel.from_process(self.grid, process)
            for process in self.collision_deck.processes
        )
        self._collision_source_components = tuple(
            self._build_collision_source_component(process, kernel)
            for process, kernel in zip(
                self.collision_deck.processes, self._collision_kernels)
        )
        self._momentum_bases = {
            "internal_faces": self._compile_momentum_basis(
                self.grid.boundaries[1:-1]),
            "upper_boundaries": self._compile_momentum_basis(
                self.grid.boundaries[1:]),
            "centers": self._compile_momentum_basis(
                self.grid.cell_centers_eV),
        }

    def _validate_condition(self, condition: TwoTermBoltzmannCondition):
        if not isinstance(condition, TwoTermBoltzmannCondition):
            raise TypeError("a TwoTermBoltzmannCondition is required")
        deck_targets = set(self.collision_deck.targets)
        fraction_targets = {
            name
            for name, fraction in condition.target_mole_fractions.items()
            if fraction > 0.0
        }
        if fraction_targets != deck_targets:
            raise ValueError(
                "positive target mole fractions must exactly match deck targets"
            )
        if any(
            process.kind == "EFFECTIVE"
            for process in self.collision_deck.processes
        ):
            raise ValueError(
                "effective momentum-transfer rows require an explicit "
                "elastic/inelastic deconvolution before this solver"
            )
        for target in deck_targets:
            momentum_count = sum(
                process.target == target
                and process.kind in {"ELASTIC", "MOMENTUM"}
                for process in self.collision_deck.processes
            )
            if momentum_count != 1:
                raise ValueError(
                    "each target requires exactly one elastic or momentum row"
                )

    @staticmethod
    def _cross_section_at(
        process: ElectronCollisionProcess,
        energies_eV: np.ndarray,
    ) -> np.ndarray:
        nodes, cross_sections, resolved_lower = _collision_support_nodes(
            process)
        if energies_eV[-1] > nodes[-1]:
            raise ValueError(
                "energy grid exceeds collision cross-section support"
            )
        if process.kind in _MOMENTUM_KINDS and energies_eV[0] < resolved_lower:
            raise ValueError(
                "energy grid extends below momentum cross-section support"
            )
        values = np.interp(
            np.clip(energies_eV, nodes[0], nodes[-1]), nodes, cross_sections)
        if process.kind not in _MOMENTUM_KINDS:
            values = np.where(energies_eV < nodes[0], 0.0, values)
        return values

    def _transport_coefficients(
        self,
        condition: TwoTermBoltzmannCondition,
        growth_rate_coefficient_m3_s: float,
        electron_electron_coefficients=None,
    ) -> tuple[np.ndarray, np.ndarray]:
        face_energy = self.grid.boundaries[1:-1]
        sigma_m, sigma_epsilon = self._weighted_momentum_basis(
            condition, "internal_faces")
        growth_correction = np.divide(
            growth_rate_coefficient_m3_s,
            ELECTRON_SPEED_PER_SQRT_EV_M_S * np.sqrt(face_energy),
        )
        sigma_tilde = sigma_m + growth_correction
        if np.any(sigma_m <= 0.0) or np.any(sigma_tilde <= 0.0):
            raise ValueError(
                "nonpositive effective momentum-transfer cross section"
            )
        reduced_field = condition.reduced_electric_field_Td * 1.0e-21
        gas_thermal_energy_eV = (
            K_BOLTZMANN_J_K * condition.gas_temperature_K / E_CHARGE_C)
        reduced_frequency = (
            condition.angular_field_frequency_over_density_m3_s)
        if reduced_frequency == 0.0:
            # Preserve the established DC arithmetic path bit-for-bit.
            field_heating = (
                (reduced_field ** 2 / 3.0) * face_energy / sigma_tilde)
        else:
            equivalent_ac_cross_section = np.divide(
                reduced_frequency,
                ELECTRON_SPEED_PER_SQRT_EV_M_S * np.sqrt(face_energy),
            )
            field_heating_inverse_cross_section = np.divide(
                sigma_tilde,
                sigma_tilde ** 2 + equivalent_ac_cross_section ** 2,
            )
            field_heating = (
                (reduced_field ** 2 / 3.0)
                * face_energy * field_heating_inverse_cross_section)
        drift = (
            -ELECTRON_SPEED_PER_SQRT_EV_M_S
            * face_energy ** 2 * sigma_epsilon
        )
        diffusion = ELECTRON_SPEED_PER_SQRT_EV_M_S * (
            field_heating
            + gas_thermal_energy_eV
            * face_energy ** 2 * sigma_epsilon
        )
        if electron_electron_coefficients is not None:
            drift = (
                drift
                + np.asarray(
                    electron_electron_coefficients.drift_eV_m3_s,
                    dtype=float,
                )
            )
            diffusion = (
                diffusion
                + np.asarray(
                    electron_electron_coefficients.diffusion_eV2_m3_s,
                    dtype=float,
                )
            )
            if drift.shape != face_energy.shape or diffusion.shape != (
                face_energy.shape
            ):
                raise ValueError(
                    "electron-electron coefficients differ from solver grid")
        if np.any(diffusion <= 0.0) or np.any(~np.isfinite(diffusion)):
            raise ValueError("nonpositive electron energy diffusion")
        return drift, diffusion

    def _momentum_cross_sections(
        self,
        condition: TwoTermBoltzmannCondition,
        energies_eV: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        energies = np.asarray(energies_eV, dtype=float)
        sigma_m = np.zeros_like(energies)
        sigma_epsilon = np.zeros_like(energies)
        for process in self.collision_deck.processes:
            fraction = condition.target_mole_fractions[process.target]
            cross_section = self._cross_section_at(process, energies)
            sigma_m += fraction * cross_section
            if process.kind in {"ELASTIC", "MOMENTUM"}:
                sigma_epsilon += (
                    2.0 * float(process.mass_ratio)
                    * fraction * cross_section
                )
        return sigma_m, sigma_epsilon

    def _compile_momentum_basis(
        self,
        energies_eV: np.ndarray,
    ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        """Compile target-resolved momentum arrays on one fixed solver grid."""
        energies = np.asarray(energies_eV, dtype=float)
        basis = {
            target: (np.zeros_like(energies), np.zeros_like(energies))
            for target in self.collision_deck.targets
        }
        for process in self.collision_deck.processes:
            momentum, energy_transfer = basis[process.target]
            cross_section = self._cross_section_at(process, energies)
            momentum += cross_section
            if process.kind in {"ELASTIC", "MOMENTUM"}:
                energy_transfer += (
                    2.0 * float(process.mass_ratio) * cross_section)
        return basis

    def _weighted_momentum_basis(
        self,
        condition: TwoTermBoltzmannCondition,
        location: str,
    ) -> tuple[np.ndarray, np.ndarray]:
        basis = self._momentum_bases[location]
        first = next(iter(basis.values()))
        sigma_m = np.zeros_like(first[0])
        sigma_epsilon = np.zeros_like(first[1])
        for target, (momentum, energy_transfer) in basis.items():
            fraction = condition.target_mole_fractions[target]
            sigma_m += fraction * momentum
            sigma_epsilon += fraction * energy_transfer
        return sigma_m, sigma_epsilon

    def _build_collision_source_component(
        self,
        process: ElectronCollisionProcess,
        kernel: ElectronCollisionMomentKernel,
    ) -> sparse.csr_matrix:
        """Compile one unit-target-fraction scattering operator."""
        cell_count = self.grid.cell_count
        boundaries = self.grid.boundaries
        source = np.zeros((cell_count, cell_count))
        if process.kind in _MOMENTUM_KINDS:
            return sparse.csr_matrix(source)
        out_weights = np.asarray(kernel.rate_weights_m3_s)
        source[np.arange(cell_count), np.arange(cell_count)] -= out_weights
        if process.kind == "ATTACHMENT":
            return sparse.csr_matrix(source)
        if process.kind in {"EXCITATION", "ROTATION"}:
            energy_loss = float(process.energy_loss_eV or 0.0)
            parent_lower_from_child = lambda value: value + energy_loss
            parent_upper_from_child = lambda value: value + energy_loss
            in_factor = 1.0
        elif process.kind == "IONIZATION":
            energy_loss = float(process.energy_loss_eV or 0.0)
            outgoing_electron_count = 1 + process.electron_number_change
            parent_lower_from_child = (
                lambda value: outgoing_electron_count * value + energy_loss)
            parent_upper_from_child = (
                lambda value: outgoing_electron_count * value + energy_loss)
            in_factor = float(outgoing_electron_count)
        else:
            raise ValueError("unsupported inelastic collision kind")

        energy, cross_section, _ = _collision_support_nodes(process)
        integral = _PiecewiseLinearCrossSectionIntegral(
            energy, cross_section, powers=(1,))
        for destination in range(cell_count):
            parent_lower = parent_lower_from_child(boundaries[destination])
            parent_upper = parent_upper_from_child(
                boundaries[destination + 1])
            parent_lower = max(parent_lower, boundaries[0], energy[0])
            parent_upper = min(parent_upper, boundaries[-1], energy[-1])
            if parent_upper <= parent_lower:
                continue
            first_parent = max(
                0,
                int(np.searchsorted(
                    boundaries, parent_lower, side="right") - 1),
            )
            last_parent = min(
                cell_count - 1,
                int(np.searchsorted(
                    boundaries, parent_upper, side="left")),
            )
            for parent in range(first_parent, last_parent + 1):
                lower = max(parent_lower, boundaries[parent])
                upper = min(parent_upper, boundaries[parent + 1])
                if upper <= lower:
                    continue
                event_integral = (
                    ELECTRON_SPEED_PER_SQRT_EV_M_S
                    * integral.integrate(
                        lower,
                        upper,
                        energy_power=1,
                    )
                )
                source[destination, parent] += in_factor * event_integral
        return sparse.csr_matrix(source)

    def _collision_source_matrix(
        self,
        condition: TwoTermBoltzmannCondition,
    ) -> tuple[sparse.csr_matrix, tuple[ElectronCollisionMomentKernel, ...]]:
        source = sparse.csr_matrix(
            (self.grid.cell_count, self.grid.cell_count))
        for process, component in zip(
            self.collision_deck.processes,
            self._collision_source_components,
        ):
            fraction = condition.target_mole_fractions[process.target]
            if fraction > 0.0 and component.nnz:
                source = source + fraction * component
        return source, self._collision_kernels

    def _transport_moments(
        self,
        condition: TwoTermBoltzmannCondition,
        values: np.ndarray,
        growth_rate_coefficient_m3_s: float,
    ) -> TwoTermElectronTransportMoments:
        """Evaluate Hagelaar--Pitchford flux moments on the solved EEPF."""
        boundaries = self.grid.boundaries
        centers = self.grid.cell_centers_eV
        momentum_boundary, _ = self._weighted_momentum_basis(
            condition, "upper_boundaries")
        effective_boundary = momentum_boundary + np.divide(
            growth_rate_coefficient_m3_s,
            ELECTRON_SPEED_PER_SQRT_EV_M_S * np.sqrt(boundaries[1:]),
        )
        momentum_center, _ = self._weighted_momentum_basis(
            condition, "centers")
        effective_center = momentum_center + np.divide(
            growth_rate_coefficient_m3_s,
            ELECTRON_SPEED_PER_SQRT_EV_M_S * np.sqrt(centers),
        )
        if (
            np.any(effective_boundary <= 0.0)
            or np.any(effective_center <= 0.0)
        ):
            raise ValueError(
                "nonpositive effective momentum cross section in moments"
            )

        derivative = np.concatenate((
            np.array([0.0]),
            np.diff(values) / np.diff(centers),
            np.array([0.0]),
        ))
        mobility_integrand = np.zeros_like(boundaries)
        mobility_integrand[1:] = (
            derivative[1:]
            * boundaries[1:]
            / effective_boundary
        )
        reduced_mobility = -(
            ELECTRON_SPEED_PER_SQRT_EV_M_S / 3.0
        ) * simpson(mobility_integrand, x=boundaries)
        reduced_frequency = (
            condition.angular_field_frequency_over_density_m3_s)
        if reduced_frequency == 0.0:
            dissipative_reduced_mobility = reduced_mobility
        else:
            equivalent_ac_cross_section = np.zeros_like(boundaries)
            equivalent_ac_cross_section[1:] = np.divide(
                reduced_frequency,
                ELECTRON_SPEED_PER_SQRT_EV_M_S * np.sqrt(boundaries[1:]),
            )
            dissipative_integrand = np.zeros_like(boundaries)
            dissipative_integrand[1:] = (
                derivative[1:]
                * boundaries[1:]
                * effective_boundary
                / (
                    effective_boundary ** 2
                    + equivalent_ac_cross_section[1:] ** 2
                )
            )
            dissipative_reduced_mobility = -(
                ELECTRON_SPEED_PER_SQRT_EV_M_S / 3.0
            ) * simpson(dissipative_integrand, x=boundaries)

        # The EEPF is piecewise constant in this solver.  Midpoint integration
        # is therefore the representation-consistent scalar-diffusion moment.
        reduced_diffusion = (
            ELECTRON_SPEED_PER_SQRT_EV_M_S / 3.0
        ) * float(np.sum(
            values
            * centers
            / effective_center
            * self.grid.cell_widths_eV
        ))
        reduced_field_V_m2 = condition.reduced_electric_field_Td * 1.0e-21
        speed_weights = 0.5 * (
            boundaries[1:] ** 2 - boundaries[:-1] ** 2)
        wall_energy_weights = (1.0 / 3.0) * (
            boundaries[1:] ** 3 - boundaries[:-1] ** 3)
        reduced_speed_moment = float(np.dot(values, speed_weights))
        mean_speed = (
            ELECTRON_SPEED_PER_SQRT_EV_M_S * reduced_speed_moment)
        mean_wall_energy = float(
            np.dot(values, wall_energy_weights) / reduced_speed_moment)
        return TwoTermElectronTransportMoments(
            flux_reduced_mobility_m_inv_V_inv_s_inv=reduced_mobility,
            scalar_reduced_diffusion_m_inv_s_inv=reduced_diffusion,
            dissipative_reduced_mobility_m_inv_V_inv_s_inv=(
                dissipative_reduced_mobility),
            reduced_field_power_gain_eV_m3_s=(
                dissipative_reduced_mobility * reduced_field_V_m2 ** 2
            ),
            mean_electron_speed_m_s=mean_speed,
            isotropic_wall_flux_coefficient_m_s=0.25 * mean_speed,
            mean_wall_loss_electron_energy_eV=mean_wall_energy,
        )

    @staticmethod
    def _solve_normalized_nullspace(
        matrix: sparse.csr_matrix,
        normalization_weights: np.ndarray,
        compatibility_direction: np.ndarray,
        *,
        allow_physical_negative: bool = False,
    ) -> tuple[np.ndarray, float, float, float]:
        augmented = sparse.bmat(
            (
                (matrix, sparse.csr_matrix(
                    compatibility_direction[:, np.newaxis])),
                (sparse.csr_matrix(normalization_weights[np.newaxis, :]),
                 sparse.csr_matrix((1, 1))),
            ),
            format="csr",
        )
        state = np.asarray(spsolve(
            augmented,
            np.concatenate((np.zeros(matrix.shape[0]), np.array([1.0]))),
        ), dtype=float)
        if np.any(~np.isfinite(state)):
            raise FloatingPointError("nonfinite two-term Boltzmann state")
        values = state[:-1]
        minimum = float(np.min(values))
        scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
        if minimum < -1.0e-12 * scale and not allow_physical_negative:
            raise FloatingPointError("negative two-term Boltzmann EEPF")
        negative_population = float(np.dot(
            np.maximum(-values, 0.0), normalization_weights))
        if minimum < 0.0 and not allow_physical_negative:
            # Sparse direct solves can leave signed roundoff in an
            # exponentially small tail.  The projected mass and raw minimum
            # are returned as receipts; a physical-scale negative solution
            # failed above and can never be repaired here.
            values = np.maximum(values, 0.0)
            values /= np.dot(values, normalization_weights)
        return values, float(state[-1]), minimum, negative_population

    def _temporal_growth_eigenstate(
        self,
        condition: TwoTermBoltzmannCondition,
        collision_source: sparse.csr_matrix,
        *,
        relative_tolerance: float,
        maximum_iterations: int,
        electron_electron_coefficients=None,
        initial_growth_rate_coefficient_m3_s: float | None = None,
    ) -> tuple[
        np.ndarray,
        float,
        int,
        float,
        float,
        float,
        float,
        ScharfetterGummelEnergyOperator,
        int,
    ]:
        """Solve the PT growth coefficient as a bordered scalar eigen-root."""
        normalization_weights = self.grid.normalization_weights
        compatibility_direction = self.grid.cell_widths_eV
        column_growth = np.asarray(collision_source.sum(axis=0)).ravel()
        population_specific_growth = np.divide(
            column_growth,
            normalization_weights,
        )
        lower = float(np.min(population_specific_growth))
        upper = float(np.max(population_specific_growth))
        face_energy = self.grid.boundaries[1:-1]
        sigma_m, _ = self._weighted_momentum_basis(
            condition, "internal_faces")
        physical_lower = float(np.max(
            -sigma_m
            * ELECTRON_SPEED_PER_SQRT_EV_M_S
            * np.sqrt(face_energy)
        ))
        rate_scale = max(abs(lower), abs(upper), 1.0e-30)
        lower = max(lower, physical_lower + 1.0e-12 * rate_scale)
        if upper <= lower:
            raise ValueError("temporal-growth eigenvalue domain is empty")

        def state_at(growth: float):
            drift, diffusion = self._transport_coefficients(
                condition, growth, electron_electron_coefficients)
            transport = ScharfetterGummelEnergyOperator(
                self.grid, drift, diffusion)
            matrix = (
                transport._matrix(0)
                - collision_source
                + sparse.diags(
                    growth * normalization_weights, format="csr")
            )
            state = self._solve_normalized_nullspace(
                matrix,
                normalization_weights,
                compatibility_direction,
                allow_physical_negative=True,
            )
            return state, transport, matrix

        state_by_growth = {}

        def evaluate_growth(growth: float):
            growth = float(growth)
            if growth in state_by_growth:
                cached = state_by_growth[growth]
                return None if cached is None else float(cached[0][1])
            try:
                evaluated_state = state_at(growth)
            except (FloatingPointError, ValueError):
                state_by_growth[growth] = None
                return None
            compatibility = evaluated_state[0][1]
            if math.isfinite(compatibility):
                state_by_growth[growth] = evaluated_state
                return float(compatibility)
            state_by_growth[growth] = None
            return None

        brackets = []
        if initial_growth_rate_coefficient_m3_s is not None:
            initial = float(initial_growth_rate_coefficient_m3_s)
            if math.isfinite(initial):
                center = min(max(initial, lower), upper)
                center_value = evaluate_growth(center)
                span = max(
                    1.0e-8 * (upper - lower),
                    1.0e-6 * max(abs(center), rate_scale),
                    np.finfo(float).tiny,
                )
                for _ in range(14):
                    left = max(lower, center - span)
                    right = min(upper, center + span)
                    left_value = evaluate_growth(left)
                    right_value = evaluate_growth(right)
                    candidates = (
                        (left, left_value, center, center_value),
                        (center, center_value, right, right_value),
                        (left, left_value, right, right_value),
                    )
                    for x0, y0, x1, y1 in candidates:
                        if y0 == 0.0:
                            brackets = [(x0, x0)]
                            break
                        if y1 == 0.0:
                            brackets = [(x1, x1)]
                            break
                        if (
                            y0 is not None
                            and y1 is not None
                            and y0 * y1 < 0.0
                        ):
                            brackets = [(x0, x1)]
                            break
                    if brackets or (left == lower and right == upper):
                        break
                    span *= 4.0

        if not brackets:
            sample_count = min(
                max(33, int(math.sqrt(maximum_iterations)) * 4), 129)
            samples = np.linspace(lower, upper, sample_count)
            if lower < 0.0 < upper:
                # Attaching gases can place the physical eigen-root many
                # orders of magnitude closer to zero than either pure column
                # bound. Preserve the exhaustive signed-log fallback.
                logarithmic_fraction = np.geomspace(1.0e-12, 1.0, 49)
                samples = np.unique(np.concatenate((
                    samples,
                    -abs(lower) * logarithmic_fraction,
                    np.array([0.0]),
                    upper * logarithmic_fraction,
                )))
                samples = samples[(samples >= lower) & (samples <= upper)]
            for growth in samples:
                evaluate_growth(float(growth))
        evaluated = sorted(
            (growth, float(state[0][1]))
            for growth, state in state_by_growth.items()
            if state is not None
        )
        if not brackets and len(evaluated) < 2:
            raise RuntimeError("could not evaluate temporal-growth eigen-root")

        if not brackets:
            for left_item, right_item in zip(evaluated, evaluated[1:]):
                left_growth, left_value = left_item
                right_growth, right_value = right_item
                if left_value == 0.0:
                    brackets.append((left_growth, left_growth))
                elif left_value * right_value < 0.0:
                    brackets.append((left_growth, right_growth))
            if evaluated[-1][1] == 0.0:
                brackets.append((evaluated[-1][0], evaluated[-1][0]))
        if not brackets:
            raise RuntimeError("temporal-growth compatibility root is unbracketed")

        candidates = []
        total_iterations = 0
        # Solve the scalar eigen-root to floating-point precision.  The
        # user tolerance is a convergence gate on the resulting conservation
        # closure, not permission to stop the root solve early.
        absolute_tolerance = np.finfo(float).tiny
        for bracket_lower, bracket_upper in brackets:
            if bracket_lower == bracket_upper:
                root = bracket_lower
                iterations = 0
            else:
                root, result = brentq(
                    lambda growth: evaluate_growth(growth),
                    bracket_lower,
                    bracket_upper,
                    xtol=absolute_tolerance,
                    rtol=4.0 * np.finfo(float).eps,
                    maxiter=maximum_iterations,
                    full_output=True,
                    disp=True,
                )
                iterations = int(result.iterations)
            total_iterations += iterations
            (
                raw_values,
                compatibility,
                minimum,
                negative_population,
            ), transport, matrix = state_by_growth.get(root) or state_at(root)
            scale = max(
                float(np.max(np.abs(raw_values))),
                np.finfo(float).tiny,
            )
            if minimum < -1.0e-10 * scale:
                continue
            values = np.maximum(raw_values, 0.0)
            values /= np.dot(values, normalization_weights)
            collision_growth = float(np.sum(collision_source @ values))
            closure_error = collision_growth - root
            candidates.append((
                abs(closure_error) / rate_scale,
                values,
                root,
                total_iterations,
                compatibility,
                minimum,
                negative_population,
                transport,
                matrix,
            ))
        if not candidates:
            raise RuntimeError(
                "temporal-growth roots have no positive normalized eigenstate"
            )
        candidates.sort(key=lambda item: item[0])
        (
            relative_closure,
            values,
            growth,
            iteration_count,
            compatibility,
            minimum,
            negative_population,
            transport,
            _matrix,
        ) = candidates[0]
        return (
            values,
            growth,
            iteration_count,
            relative_closure,
            compatibility,
            minimum,
            negative_population,
            transport,
            len(state_by_growth),
        )

    def solve(
        self,
        condition: TwoTermBoltzmannCondition,
        *,
        initial_solution: TwoTermBoltzmannSolution | None = None,
        relative_tolerance: float = 1.0e-9,
        maximum_iterations: int = 200,
        damping: float = 0.7,
        maximum_tail_population_fraction: float = 1.0e-8,
    ) -> TwoTermBoltzmannSolution:
        self._validate_condition(condition)
        if initial_solution is not None:
            if not isinstance(initial_solution, TwoTermBoltzmannSolution):
                raise TypeError("initial_solution must be a two-term solution")
            if (
                initial_solution.collision_deck_sha256
                != self.collision_deck.payload_sha256
                or not np.array_equal(
                    initial_solution.distribution.grid.boundaries,
                    self.grid.boundaries,
                )
            ):
                raise ValueError(
                    "initial_solution must use the same grid and collision deck"
                )
        tolerance = float(relative_tolerance)
        if (
            not math.isfinite(tolerance)
            or tolerance <= 0.0
            or int(maximum_iterations) != maximum_iterations
            or maximum_iterations < 1
            or not math.isfinite(damping)
            or not 0.0 < damping <= 1.0
        ):
            raise ValueError("invalid two-term solver controls")
        collision_source, kernels = self._collision_source_matrix(condition)
        number_changing = any(
            process.kind in {"ATTACHMENT", "IONIZATION"}
            for process in self.collision_deck.processes
        )
        if condition.growth_model == "no_growth" and number_changing:
            raise ValueError(
                "number-changing collisions require temporal_growth"
            )

        normalization_weights = self.grid.normalization_weights
        compatibility_direction = self.grid.cell_widths_eV
        coulomb_coefficients = None
        coulomb_logarithm = None
        nonlinear_iteration_count = 0
        nonlinear_weighted_residual = 0.0
        growth_root_evaluations = 0
        if condition.growth_model == "temporal_growth":
            growth_seed = (
                None
                if initial_solution is None
                else initial_solution.net_growth_rate_coefficient_m3_s
            )

            def solve_temporal(frozen_coulomb=None):
                nonlocal growth_seed, growth_root_evaluations
                result = self._temporal_growth_eigenstate(
                    condition,
                    collision_source,
                    relative_tolerance=tolerance,
                    maximum_iterations=int(maximum_iterations),
                    electron_electron_coefficients=frozen_coulomb,
                    initial_growth_rate_coefficient_m3_s=growth_seed,
                )
                growth_seed = result[1]
                growth_root_evaluations += result[-1]
                return result

            if condition.electron_electron_coulomb_model != "none":
                from dataclasses import replace

                from .electron_coulomb import (
                    IsotropicElectronElectronCoulombKernel,
                )

                coulomb_kernel = IsotropicElectronElectronCoulombKernel(
                    self.grid)
                previous_coefficients = None
                if initial_solution is None:
                    state = solve_temporal()
                    previous_values = state[0]
                else:
                    state = None
                    previous_values = np.asarray(
                        initial_solution.distribution
                        .eepf_eV_minus_3_over_2,
                        dtype=float,
                    )
                for nonlinear_iteration_count in range(
                    1, int(maximum_iterations) + 1
                ):
                    evaluated = coulomb_kernel.evaluate(
                        ElectronEnergyDistribution(
                            self.grid, previous_values),
                        electron_to_neutral_density_ratio=(
                            condition.electron_to_neutral_density_ratio),
                        gas_number_density_m3=(
                            condition.gas_number_density_m3),
                    )
                    if previous_coefficients is None:
                        frozen = evaluated
                    else:
                        frozen = replace(
                            evaluated,
                            drift_eV_m3_s=(
                                damping * evaluated.drift_eV_m3_s
                                + (1.0 - damping)
                                * previous_coefficients.drift_eV_m3_s),
                            diffusion_eV2_m3_s=(
                                damping * evaluated.diffusion_eV2_m3_s
                                + (1.0 - damping)
                                * previous_coefficients.diffusion_eV2_m3_s),
                        )
                    new_state = solve_temporal(frozen)
                    new_values = new_state[0]
                    nonlinear_weighted_residual = float(np.sum(
                        np.abs(new_values - previous_values)
                        * normalization_weights
                    ))
                    state = new_state
                    previous_values = new_values
                    previous_coefficients = frozen
                    if nonlinear_weighted_residual <= tolerance:
                        break
                else:
                    raise RuntimeError(
                        "electron-electron Coulomb fixed point did not converge"
                    )
                coulomb_coefficients = coulomb_kernel.evaluate(
                    ElectronEnergyDistribution(self.grid, state[0]),
                    electron_to_neutral_density_ratio=(
                        condition.electron_to_neutral_density_ratio),
                    gas_number_density_m3=condition.gas_number_density_m3,
                )
                coulomb_logarithm = float(
                    coulomb_coefficients.coulomb_logarithm)
            else:
                state = solve_temporal()
            (
                final_values,
                final_growth,
                iteration_count,
                weighted_residual,
                compatibility,
                minimum_raw_eepf,
                negative_population,
                transport,
                _final_growth_root_evaluations,
            ) = state
            # The bordered solve can carry signed roundoff in a vanishing
            # high-energy tail.  Projection makes the returned EEPF physical;
            # close the reported growth and its final operator on that exact
            # returned state rather than on the pre-projection vector.
            eigen_growth = final_growth
            final_growth = float(np.sum(collision_source @ final_values))
            growth_scale = max(abs(eigen_growth), abs(final_growth), 1.0e-30)
            weighted_residual = max(
                weighted_residual,
                abs(final_growth - eigen_growth) / growth_scale,
            )
            if weighted_residual > tolerance:
                raise RuntimeError(
                    "temporal-growth conservation closure did not converge"
                )
            drift, diffusion = self._transport_coefficients(
                condition, final_growth, coulomb_coefficients)
            transport = ScharfetterGummelEnergyOperator(
                self.grid, drift, diffusion)
            final_matrix = (
                transport._matrix(0)
                - collision_source
                + sparse.diags(
                    final_growth * normalization_weights,
                    format="csr",
                )
            )
        else:
            if condition.electron_electron_coulomb_model != "none":
                raise ValueError(
                    "electron-electron Coulomb coupling currently requires "
                    "temporal_growth")
            drift, diffusion = self._transport_coefficients(condition, 0.0)
            transport = ScharfetterGummelEnergyOperator(
                self.grid, drift, diffusion)
            final_matrix = transport._matrix(0) - collision_source
            (
                final_values,
                compatibility,
                minimum_raw_eepf,
                negative_population,
            ) = self._solve_normalized_nullspace(
                final_matrix,
                normalization_weights,
                compatibility_direction,
            )
            final_growth = 0.0
            iteration_count = 1
            weighted_residual = 0.0
        equation_residual = np.asarray(final_matrix @ final_values)
        distribution = ElectronEnergyDistribution(self.grid, final_values)
        tail_receipt = distribution.convergence_receipt(
            tail_cell_count=min(4, self.grid.cell_count),
            maximum_tail_population_fraction=(
                maximum_tail_population_fraction),
        )
        if not tail_receipt["tail_passed"]:
            raise ValueError(
                "electron-energy upper boundary failed tail convergence"
            )

        moments = tuple(
            kernel.evaluate(
                distribution,
                maximum_unresolved_population_fraction=(
                    maximum_tail_population_fraction),
            )
            for kernel in kernels
        )
        transport_moments = self._transport_moments(
            condition, final_values, final_growth)
        left, right = transport._face_coefficients(0)
        fluxes = np.zeros(self.grid.cell_count + 1)
        fluxes[1:-1] = (
            left * final_values[:-1] + right * final_values[1:])
        particle_source = float(np.sum(collision_source @ final_values))
        expected_particle_source = (
            final_growth
            if condition.growth_model == "temporal_growth"
            else 0.0
        )
        return TwoTermBoltzmannSolution(
            distribution=distribution,
            collision_moments=moments,
            transport_moments=transport_moments,
            reduced_electric_field_Td=condition.reduced_electric_field_Td,
            gas_temperature_K=condition.gas_temperature_K,
            angular_field_frequency_over_density_m3_s=(
                condition.angular_field_frequency_over_density_m3_s),
            growth_model=condition.growth_model,
            net_growth_rate_coefficient_m3_s=final_growth,
            iteration_count=iteration_count,
            weighted_iteration_residual=weighted_residual,
            maximum_equation_residual_m3_s=float(
                np.max(np.abs(equation_residual))),
            particle_source_from_collision_operator_m3_s=particle_source,
            particle_growth_closure_error_m3_s=(
                particle_source - expected_particle_source),
            minimum_raw_eepf_before_roundoff_projection=minimum_raw_eepf,
            roundoff_negative_population_fraction=negative_population,
            energy_flux_faces_reduced=fluxes,
            collision_deck_sha256=self.collision_deck.payload_sha256,
            electron_electron_coulomb_model=(
                condition.electron_electron_coulomb_model),
            coulomb_logarithm=coulomb_logarithm,
            nonlinear_iteration_count=nonlinear_iteration_count,
            nonlinear_weighted_residual=nonlinear_weighted_residual,
            growth_root_evaluations=growth_root_evaluations,
        )
