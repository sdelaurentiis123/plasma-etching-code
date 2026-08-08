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
    lower = max(float(lower_eV), float(energy[0]))
    upper = min(float(upper_eV), float(energy[-1]))
    if upper <= lower:
        return 0.0
    knots = np.concatenate((
        np.array([lower]),
        energy[(energy > lower) & (energy < upper)],
        np.array([upper]),
    ))
    sigma = np.interp(knots, energy, cross_section)
    integral = 0.0
    for left, right, sigma_left, sigma_right in zip(
        knots[:-1], knots[1:], sigma[:-1], sigma[1:]
    ):
        slope = (sigma_right - sigma_left) / (right - left)
        intercept = sigma_left - slope * left
        integral += (
            intercept / (energy_power + 1.0)
            * (right ** (energy_power + 1) - left ** (energy_power + 1))
            + slope / (energy_power + 2.0)
            * (right ** (energy_power + 2) - left ** (energy_power + 2))
        )
    return float(integral)


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
                    _integrate_energy_cross_section(
                        energy,
                        sigma,
                        overlap_lower,
                        overlap_upper,
                        energy_power=1,
                    )
                )
                energy_weights[cell] += ELECTRON_SPEED_PER_SQRT_EV_M_S * (
                    _integrate_energy_cross_section(
                        energy,
                        sigma,
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
    """One local-field gas state for the deterministic two-term solve."""

    reduced_electric_field_Td: float
    gas_temperature_K: float
    target_mole_fractions: Mapping[str, float]
    growth_model: str = "no_growth"
    inelastic_momentum_closure: str = "isotropic_source_reproduction"
    ionization_energy_sharing: str = "equal_sharing"
    initial_electron_temperature_eV: float = 2.0

    def __post_init__(self):
        field = float(self.reduced_electric_field_Td)
        temperature = float(self.gas_temperature_K)
        fractions = {
            str(name).strip(): float(value)
            for name, value in dict(self.target_mole_fractions).items()
        }
        if (
            not math.isfinite(field)
            or field < 0.0
            or not math.isfinite(temperature)
            or temperature <= 0.0
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


@dataclass(frozen=True)
class TwoTermBoltzmannSolution:
    """Converged local-field two-term EEPF with explicit evidence limits."""

    distribution: ElectronEnergyDistribution
    collision_moments: tuple[ElectronCollisionMoments, ...]
    reduced_electric_field_Td: float
    gas_temperature_K: float
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
    ) -> tuple[np.ndarray, np.ndarray]:
        face_energy = self.grid.boundaries[1:-1]
        sigma_m = np.zeros_like(face_energy)
        sigma_epsilon = np.zeros_like(face_energy)
        for process in self.collision_deck.processes:
            fraction = condition.target_mole_fractions[process.target]
            cross_section = self._cross_section_at(process, face_energy)
            sigma_m += fraction * cross_section
            if process.kind in {"ELASTIC", "MOMENTUM"}:
                sigma_epsilon += (
                    2.0 * float(process.mass_ratio)
                    * fraction * cross_section
                )
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
        drift = (
            -ELECTRON_SPEED_PER_SQRT_EV_M_S
            * face_energy ** 2 * sigma_epsilon
        )
        diffusion = ELECTRON_SPEED_PER_SQRT_EV_M_S * (
            (reduced_field ** 2 / 3.0) * face_energy / sigma_tilde
            + gas_thermal_energy_eV
            * face_energy ** 2 * sigma_epsilon
        )
        if np.any(diffusion <= 0.0) or np.any(~np.isfinite(diffusion)):
            raise ValueError("nonpositive electron energy diffusion")
        return drift, diffusion

    def _collision_source_matrix(
        self,
        condition: TwoTermBoltzmannCondition,
    ) -> tuple[sparse.csr_matrix, tuple[ElectronCollisionMomentKernel, ...]]:
        cell_count = self.grid.cell_count
        boundaries = self.grid.boundaries
        source = np.zeros((cell_count, cell_count))
        kernels = []
        for process in self.collision_deck.processes:
            kernel = ElectronCollisionMomentKernel.from_process(
                self.grid, process)
            kernels.append(kernel)
            if process.kind in _MOMENTUM_KINDS:
                continue
            fraction = condition.target_mole_fractions[process.target]
            out_weights = fraction * np.asarray(kernel.rate_weights_m3_s)
            source[np.arange(cell_count), np.arange(cell_count)] -= out_weights
            if process.kind == "ATTACHMENT":
                continue
            if process.kind in {"EXCITATION", "ROTATION"}:
                energy_loss = float(process.energy_loss_eV or 0.0)
                parent_lower_from_child = lambda value: value + energy_loss
                parent_upper_from_child = lambda value: value + energy_loss
                in_factor = 1.0
            elif process.kind == "IONIZATION":
                energy_loss = float(process.energy_loss_eV or 0.0)
                parent_lower_from_child = lambda value: 2.0 * value + energy_loss
                parent_upper_from_child = lambda value: 2.0 * value + energy_loss
                in_factor = 2.0
            else:
                raise ValueError("unsupported inelastic collision kind")

            energy, cross_section, _ = _collision_support_nodes(process)
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
                        * fraction
                        * _integrate_energy_cross_section(
                            energy,
                            cross_section,
                            lower,
                            upper,
                            energy_power=1,
                        )
                    )
                    source[destination, parent] += in_factor * event_integral
        return sparse.csr_matrix(source), tuple(kernels)

    @staticmethod
    def _solve_normalized_nullspace(
        matrix: sparse.csr_matrix,
        normalization_weights: np.ndarray,
        compatibility_direction: np.ndarray,
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
        scale = max(float(np.max(values)), 1.0)
        if minimum < -1.0e-12 * scale:
            raise FloatingPointError("negative two-term Boltzmann EEPF")
        negative_population = float(np.dot(
            np.maximum(-values, 0.0), normalization_weights))
        if minimum < 0.0:
            # Sparse direct solves can leave signed roundoff in an
            # exponentially small tail.  The projected mass and raw minimum
            # are returned as receipts; a physical-scale negative solution
            # failed above and can never be repaired here.
            values = np.maximum(values, 0.0)
            values /= np.dot(values, normalization_weights)
        return values, float(state[-1]), minimum, negative_population

    def solve(
        self,
        condition: TwoTermBoltzmannCondition,
        *,
        relative_tolerance: float = 1.0e-9,
        maximum_iterations: int = 200,
        damping: float = 0.7,
        maximum_tail_population_fraction: float = 1.0e-8,
    ) -> TwoTermBoltzmannSolution:
        self._validate_condition(condition)
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

        values = ElectronEnergyDistribution.maxwellian(
            self.grid,
            condition.initial_electron_temperature_eV,
        ).eepf_eV_minus_3_over_2
        normalization_weights = self.grid.normalization_weights
        compatibility_direction = self.grid.cell_widths_eV
        weighted_residual = math.inf
        compatibility = math.inf
        minimum_raw_eepf = 0.0
        negative_population = 0.0
        iteration_count = 0

        for iteration_count in range(1, int(maximum_iterations) + 1):
            collision_source_values = collision_source @ values
            growth = float(np.sum(collision_source_values))
            if condition.growth_model == "no_growth":
                growth = 0.0
            drift, diffusion = self._transport_coefficients(condition, growth)
            transport = ScharfetterGummelEnergyOperator(
                self.grid, drift, diffusion)
            matrix = transport._matrix(0) - collision_source
            if condition.growth_model == "temporal_growth":
                matrix = matrix + sparse.diags(
                    growth * normalization_weights, format="csr")
            (
                candidate,
                compatibility,
                minimum_raw_eepf,
                negative_population,
            ) = self._solve_normalized_nullspace(
                matrix,
                normalization_weights,
                compatibility_direction,
            )
            weighted_residual = float(np.sum(
                np.abs(candidate - values) * normalization_weights))
            values = damping * candidate + (1.0 - damping) * values
            if weighted_residual <= tolerance:
                break
        else:
            raise RuntimeError("two-term Boltzmann iteration did not converge")

        collision_source_values = collision_source @ values
        growth = float(np.sum(collision_source_values))
        if condition.growth_model == "no_growth":
            growth = 0.0
        drift, diffusion = self._transport_coefficients(condition, growth)
        transport = ScharfetterGummelEnergyOperator(self.grid, drift, diffusion)
        matrix = transport._matrix(0) - collision_source
        if condition.growth_model == "temporal_growth":
            matrix = matrix + sparse.diags(
                growth * normalization_weights, format="csr")
        (
            final_values,
            compatibility,
            minimum_raw_eepf,
            negative_population,
        ) = self._solve_normalized_nullspace(
            matrix,
            normalization_weights,
            compatibility_direction,
        )
        final_growth = float(np.sum(collision_source @ final_values))
        if condition.growth_model == "no_growth":
            final_growth = 0.0
        final_matrix = transport._matrix(0) - collision_source
        if condition.growth_model == "temporal_growth":
            final_matrix = final_matrix + sparse.diags(
                final_growth * normalization_weights, format="csr")
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
            reduced_electric_field_Td=condition.reduced_electric_field_Td,
            gas_temperature_K=condition.gas_temperature_K,
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
        )
