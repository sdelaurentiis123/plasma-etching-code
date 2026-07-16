"""Unified plasma-to-feature boundary state.

Every analytic source, sheath model, reactor solver, diagnostic reconstruction, or learned surrogate must
produce this representation. Transport engines consume it without knowing how it was generated.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping, Protocol

import numpy as np
from scipy.special import betaincinv
from scipy.stats import norm, qmc

from .sheath import CollisionlessRFSheath, ECHARGE


class BoundaryDensityModel(Protocol):
    """Normalized incident flux-density evaluator used by adjoint transport."""
    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None): ...
    @property
    def sampling_dimension(self): ...
    def sample_flux_velocity(self, unit_interval): ...


class BoundaryDensityModel2D(Protocol):
    """Normalized incident density projected into a declared 2-D model plane."""
    @property
    def sampling_dimension(self): ...
    def sample_flux_velocity_2d(self, unit_interval): ...


def _unit_interval_samples(value, dimension):
    samples = np.asarray(value, dtype=float)
    if (samples.ndim != 2 or samples.shape[1] != int(dimension)
            or np.any(~np.isfinite(samples)) or np.any(samples < 0.0)
            or np.any(samples > 1.0)):
        raise ValueError(f"unit_interval must have shape (n,{int(dimension)}) in [0,1]")
    # Inverse normal/exponential maps are infinite at the endpoints. Deterministic digital nets may
    # contain zero exactly, so move endpoints by one representable value without perturbing interiors.
    return np.clip(samples, np.nextafter(0.0, 1.0), np.nextafter(1.0, 0.0))


@dataclass(frozen=True)
class RectilinearVelocityHistogramDensity:
    """Normalized piecewise-constant joint velocity flux density.

    This representation can carry reactor/PIC output or diagnostic histograms without assuming that
    energy and angle factor. The third velocity coordinate is incident-normal and must be nonnegative.
    """
    edges: tuple[np.ndarray, np.ndarray, np.ndarray]
    probability_mass: np.ndarray

    def __post_init__(self):
        edges = tuple(np.asarray(edge, dtype=float).copy() for edge in self.edges)
        if len(edges) != 3 or any(edge.ndim != 1 or edge.size < 2 or np.any(np.diff(edge) <= 0) for edge in edges):
            raise ValueError("three strictly increasing velocity edges are required")
        if edges[2][0] < 0.0:
            raise ValueError("incident-normal velocity support must be nonnegative")
        mass = np.asarray(self.probability_mass, dtype=float).copy()
        expected = tuple(edge.size - 1 for edge in edges)
        if mass.shape != expected or np.any(mass < 0.0) or not np.all(np.isfinite(mass)):
            raise ValueError("probability_mass has invalid shape or values")
        total = float(mass.sum())
        if total <= 0.0:
            raise ValueError("histogram must have positive probability mass")
        mass /= total
        for edge in edges: edge.setflags(write=False)
        mass.setflags(write=False)
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "probability_mass", mass)

    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
        velocity = np.asarray(velocity_sqrt_eV, dtype=float)
        if velocity.shape[-1] != 3:
            raise ValueError("velocity must end in three components")
        flat = velocity.reshape(-1, 3)
        index = [np.searchsorted(self.edges[d], flat[:, d], side="right") - 1 for d in range(3)]
        valid = np.ones(flat.shape[0], dtype=bool)
        for d in range(3): valid &= (index[d] >= 0) & (index[d] < self.probability_mass.shape[d])
        result = np.full(flat.shape[0], -np.inf)
        if valid.any():
            iv = tuple(idx[valid] for idx in index)
            volume = np.ones(valid.sum())
            for d in range(3):
                volume *= np.diff(self.edges[d])[index[d][valid]]
            density = self.probability_mass[iv] / volume
            positive = density > 0.0
            selected = np.where(valid)[0]
            result[selected[positive]] = np.log(density[positive])
        return result.reshape(velocity.shape[:-1])

    @property
    def sampling_dimension(self):
        return 4

    def sample_flux_velocity(self, unit_interval):
        samples = _unit_interval_samples(unit_interval, self.sampling_dimension)
        cumulative = np.cumsum(self.probability_mass.ravel())
        flat_index = np.searchsorted(cumulative, samples[:, 0], side="right")
        flat_index = np.minimum(flat_index, cumulative.size - 1)
        index = np.column_stack(np.unravel_index(flat_index, self.probability_mass.shape))
        velocity = np.empty((samples.shape[0], 3))
        for axis in range(3):
            lower = self.edges[axis][index[:, axis]]
            upper = self.edges[axis][index[:, axis] + 1]
            velocity[:, axis] = lower + samples[:, axis + 1] * (upper - lower)
        return velocity


@dataclass(frozen=True)
class IonEnergyTransverseMaxwellianDensity:
    """Joint ion flux density from normal-energy mass and transverse ion temperature."""
    normal_energy_edges_eV: np.ndarray
    probability_mass: np.ndarray
    tangential_temperature_eV: float

    def __post_init__(self):
        edge = np.asarray(self.normal_energy_edges_eV, dtype=float).copy()
        mass = np.asarray(self.probability_mass, dtype=float).copy()
        if edge.ndim != 1 or edge.size < 2 or np.any(np.diff(edge) <= 0.0) or edge[0] < 0.0:
            raise ValueError("normal energy edges must be nonnegative and strictly increasing")
        if mass.shape != (edge.size - 1,) or np.any(mass < 0.0) or mass.sum() <= 0.0:
            raise ValueError("normal-energy probability mass has invalid shape or values")
        if self.tangential_temperature_eV <= 0.0:
            raise ValueError("tangential ion temperature must be positive")
        mass /= mass.sum(); edge.setflags(write=False); mass.setflags(write=False)
        object.__setattr__(self, "normal_energy_edges_eV", edge)
        object.__setattr__(self, "probability_mass", mass)

    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
        velocity = np.asarray(velocity_sqrt_eV, dtype=float)
        flat = velocity.reshape(-1, 3); vz = flat[:, 2]; energy = vz * vz
        index = np.searchsorted(self.normal_energy_edges_eV, energy, side="right") - 1
        valid = (vz > 0.0) & (index >= 0) & (index < self.probability_mass.size)
        result = np.full(flat.shape[0], -np.inf)
        if valid.any():
            width = np.diff(self.normal_energy_edges_eV)[index[valid]]
            p_energy = self.probability_mass[index[valid]] / width
            temperature = float(self.tangential_temperature_eV)
            transverse_log = (-np.log(np.pi * temperature)
                              - (flat[valid, 0] ** 2 + flat[valid, 1] ** 2) / temperature)
            # Transform normal energy E=vz^2 to normal speed: p(vz)=2*vz*p(E).
            positive = p_energy > 0.0
            selected = np.where(valid)[0]
            result[selected[positive]] = (transverse_log[positive]
                                          + np.log(2.0 * vz[selected[positive]] * p_energy[positive]))
        return result.reshape(velocity.shape[:-1])

    @property
    def sampling_dimension(self):
        return 3

    def sample_flux_velocity(self, unit_interval):
        samples = _unit_interval_samples(unit_interval, self.sampling_dimension)
        temperature = float(self.tangential_temperature_eV)
        velocity = np.empty((samples.shape[0], 3))
        velocity[:, :2] = np.sqrt(temperature / 2.0) * norm.ppf(samples[:, :2])
        cumulative = np.cumsum(self.probability_mass)
        index = np.searchsorted(cumulative, samples[:, 2], side="right")
        index = np.minimum(index, self.probability_mass.size - 1)
        previous = np.where(index > 0, cumulative[np.maximum(index - 1, 0)], 0.0)
        local = (samples[:, 2] - previous) / self.probability_mass[index]
        energy = (self.normal_energy_edges_eV[index]
                  + local * (self.normal_energy_edges_eV[index + 1]
                             - self.normal_energy_edges_eV[index]))
        velocity[:, 2] = np.sqrt(np.maximum(energy, 0.0))
        return velocity


def _sample_histogram_energy(unit, edges, probability_mass):
    cumulative = np.cumsum(probability_mass)
    index = np.searchsorted(cumulative, unit, side="right")
    index = np.minimum(index, probability_mass.size - 1)
    previous = np.where(index > 0, cumulative[np.maximum(index - 1, 0)], 0.0)
    local = (unit - previous) / probability_mass[index]
    return (edges[index]
            + local * (edges[index + 1] - edges[index]))


@dataclass(frozen=True)
class IonEnergyTransverseDensity2D:
    """A 2-D ion projection with one transverse thermal degree of freedom.

    This is the direct representation of a two-dimensional source model:
    normal energy is sampled from the declared histogram and the single
    in-plane transverse velocity is Maxwellian.  It avoids folding an
    unmodeled out-of-plane energy component back into the 2-D trajectory.
    """
    normal_energy_edges_eV: np.ndarray
    probability_mass: np.ndarray
    tangential_temperature_eV: float

    def __post_init__(self):
        edge = np.asarray(self.normal_energy_edges_eV, dtype=float).copy()
        mass = np.asarray(self.probability_mass, dtype=float).copy()
        if (edge.ndim != 1 or edge.size < 2 or np.any(np.diff(edge) <= 0.0)
                or edge[0] < 0.0
                or mass.shape != (edge.size - 1,) or np.any(mass < 0.0)
                or mass.sum() <= 0.0
                or not np.isfinite(self.tangential_temperature_eV)
                or self.tangential_temperature_eV <= 0.0):
            raise ValueError("invalid 2-D ion energy/transverse density")
        mass /= mass.sum()
        edge.setflags(write=False); mass.setflags(write=False)
        object.__setattr__(self, "normal_energy_edges_eV", edge)
        object.__setattr__(self, "probability_mass", mass)

    @property
    def sampling_dimension(self):
        return 2

    def sample_flux_velocity_2d(self, unit_interval):
        samples = _unit_interval_samples(unit_interval, self.sampling_dimension)
        energy = _sample_histogram_energy(
            samples[:, 1], self.normal_energy_edges_eV, self.probability_mass)
        transverse = (
            np.sqrt(float(self.tangential_temperature_eV) / 2.0)
            * norm.ppf(samples[:, 0]))
        return np.column_stack((
            transverse, np.sqrt(np.maximum(energy, 0.0))))


@dataclass(frozen=True)
class EnergyCosineAngleDensity2D:
    """Factorized 2-D energy and signed-angle incident distribution.

    The angle probability density is proportional to ``cos(theta)**p`` on
    ``[-pi/2, pi/2]``.  The beta inverse below samples that law exactly after
    the change of variables ``x=(sin(theta)+1)/2``.
    """
    energy_edges_eV: np.ndarray
    probability_mass: np.ndarray
    cosine_power: float

    def __post_init__(self):
        edge = np.asarray(self.energy_edges_eV, dtype=float).copy()
        mass = np.asarray(self.probability_mass, dtype=float).copy()
        if (edge.ndim != 1 or edge.size < 2 or np.any(np.diff(edge) <= 0.0)
                or edge[0] < 0.0
                or mass.shape != (edge.size - 1,) or np.any(mass < 0.0)
                or mass.sum() <= 0.0 or not np.isfinite(self.cosine_power)
                or self.cosine_power <= -1.0):
            raise ValueError("invalid 2-D energy/cosine-angle density")
        mass /= mass.sum()
        edge.setflags(write=False); mass.setflags(write=False)
        object.__setattr__(self, "energy_edges_eV", edge)
        object.__setattr__(self, "probability_mass", mass)

    @property
    def sampling_dimension(self):
        return 2

    def sample_flux_velocity_2d(self, unit_interval):
        samples = _unit_interval_samples(unit_interval, self.sampling_dimension)
        energy = _sample_histogram_energy(
            samples[:, 0], self.energy_edges_eV, self.probability_mass)
        beta_shape = 0.5 * (float(self.cosine_power) + 1.0)
        sine = 2.0 * betaincinv(
            beta_shape, beta_shape, samples[:, 1]) - 1.0
        sine = np.clip(sine, -1.0, 1.0)
        cosine = np.sqrt(np.maximum(1.0 - sine * sine, 0.0))
        speed = np.sqrt(np.maximum(energy, 0.0))
        return np.column_stack((speed * sine, speed * cosine))


@dataclass(frozen=True)
class MaxwellianFluxVelocityDensity:
    """Normalized half-space Maxwellian flux density in energy-scaled velocity coordinates.

    With ``|v|^2`` measured in eV, the two tangential components have density
    ``exp(-v_t^2/T)/sqrt(pi*T)`` and the positive incident-normal component has density
    ``2*vz*exp(-vz^2/T)/T``. This is the kinetic flux measure, not a fitted angular law.
    """
    temperature_eV: float

    def __post_init__(self):
        if not np.isfinite(self.temperature_eV) or self.temperature_eV <= 0.0:
            raise ValueError("Maxwellian temperature must be positive and finite")

    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
        velocity = np.asarray(velocity_sqrt_eV, dtype=float)
        if velocity.shape[-1] != 3:
            raise ValueError("velocity must end in three components")
        flat = velocity.reshape(-1, 3); vz = flat[:, 2]
        result = np.full(flat.shape[0], -np.inf)
        valid = vz > 0.0
        if valid.any():
            temperature = float(self.temperature_eV)
            energy = np.sum(flat[valid] ** 2, axis=1)
            result[valid] = (np.log(2.0 * vz[valid]) - np.log(np.pi)
                             - 2.0 * np.log(temperature) - energy / temperature)
        return result.reshape(velocity.shape[:-1])

    @property
    def sampling_dimension(self):
        return 3

    def sample_flux_velocity(self, unit_interval):
        samples = _unit_interval_samples(unit_interval, self.sampling_dimension)
        temperature = float(self.temperature_eV)
        velocity = np.empty((samples.shape[0], 3))
        velocity[:, :2] = np.sqrt(temperature / 2.0) * norm.ppf(samples[:, :2])
        velocity[:, 2] = np.sqrt(-temperature * np.log1p(-samples[:, 2]))
        return velocity


@dataclass(frozen=True)
class FoldedNormalTangentialDensity:
    """Source density rotated into a grazing-incidence surface proposal.

    The source coordinates are ``(transverse_x, transverse_y, incident_normal)``.  A vertical
    material face instead needs local coordinates ``(tangent, out_of_plane, inward_normal)``.
    This normalized pushforward maps source normal speed to signed surface tangent and folds the
    signed source transverse speed into positive surface-normal speed.  It changes only the numerical
    importance distribution; the physical density is still evaluated at the plasma boundary.
    """
    source: BoundaryDensityModel
    tangent_sign: int = 1

    def __post_init__(self):
        if int(self.tangent_sign) not in (-1, 1):
            raise ValueError("tangent_sign must be -1 or +1")
        object.__setattr__(self, "tangent_sign", int(self.tangent_sign))

    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
        velocity = np.asarray(velocity_sqrt_eV, dtype=float)
        if velocity.shape[-1] != 3:
            raise ValueError("velocity must end in three components")
        flat = velocity.reshape(-1, 3)
        tangent = self.tangent_sign * flat[:, 0]
        normal = flat[:, 2]
        valid = (tangent > 0.0) & (normal > 0.0)
        result = np.full(flat.shape[0], -np.inf)
        if valid.any():
            positive = np.column_stack((normal[valid], flat[valid, 1], tangent[valid]))
            negative = positive.copy(); negative[:, 0] *= -1.0
            phase = (None if phase_rad is None
                     else np.asarray(phase_rad).reshape(-1)[valid])
            position = (None if position_m is None
                        else np.asarray(position_m).reshape(-1, 2)[valid])
            result[valid] = np.logaddexp(
                self.source.log_flux_density(positive, phase, position),
                self.source.log_flux_density(negative, phase, position))
        return result.reshape(velocity.shape[:-1])

    @property
    def sampling_dimension(self):
        if not hasattr(self.source, "sampling_dimension"):
            raise TypeError("source density must implement deterministic flux sampling")
        return int(self.source.sampling_dimension)

    def sample_flux_velocity(self, unit_interval):
        if not hasattr(self.source, "sample_flux_velocity"):
            raise TypeError("source density must implement deterministic flux sampling")
        source = self.source.sample_flux_velocity(unit_interval)
        return np.column_stack((
            self.tangent_sign * source[:, 2], source[:, 1], np.abs(source[:, 0])))


@dataclass(frozen=True)
class MixtureBoundaryDensity:
    """Normalized mixture of boundary densities used for support-complete numerical proposals."""
    components: tuple[BoundaryDensityModel, ...]
    weight: np.ndarray

    def __post_init__(self):
        components = tuple(self.components)
        weight = np.asarray(self.weight, dtype=float).copy()
        if (not components or weight.shape != (len(components),) or np.any(weight < 0.0)
                or not np.all(np.isfinite(weight)) or weight.sum() <= 0.0):
            raise ValueError("mixture requires matching nonnegative component weights")
        weight /= weight.sum(); weight.setflags(write=False)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "weight", weight)

    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
        terms = np.stack([
            np.log(mixture_weight) + component.log_flux_density(
                velocity_sqrt_eV, phase_rad, position_m)
            for component, mixture_weight in zip(self.components, self.weight)
            if mixture_weight > 0.0
        ])
        return np.logaddexp.reduce(terms, axis=0)

    @property
    def sampling_dimension(self):
        dimensions = []
        for component in self.components:
            if not hasattr(component, "sampling_dimension"):
                raise TypeError("every mixture component must implement deterministic flux sampling")
            dimensions.append(int(component.sampling_dimension))
        return 1 + max(dimensions)

    def sample_flux_velocity(self, unit_interval):
        samples = _unit_interval_samples(unit_interval, self.sampling_dimension)
        cumulative = np.cumsum(self.weight)
        selection = np.searchsorted(cumulative, samples[:, 0], side="right")
        selection = np.minimum(selection, len(self.components) - 1)
        velocity = np.empty((samples.shape[0], 3))
        for index, component in enumerate(self.components):
            selected = selection == index
            if not np.any(selected):
                continue
            if not hasattr(component, "sample_flux_velocity"):
                raise TypeError("every mixture component must implement deterministic flux sampling")
            dimension = int(component.sampling_dimension)
            velocity[selected] = component.sample_flux_velocity(
                samples[selected, 1:1 + dimension])
        return velocity


def _readonly_array(value, shape_tail=()):
    array = np.asarray(value, dtype=float).copy()
    if array.ndim != 1 + len(shape_tail) or (shape_tail and array.shape[1:] != shape_tail):
        raise ValueError(f"expected array shape (n,{','.join(map(str, shape_tail))})")
    if not np.all(np.isfinite(array)):
        raise ValueError("boundary arrays must be finite")
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class SpeciesBoundaryState:
    """Weighted joint phase-space measure for one incident species.

    ``velocity_sqrt_eV`` has shape `(n,3)` and follows the feature-engine convention: squaring and
    summing components gives kinetic energy in eV. Component 2 is positive toward the feature.
    """
    name: str
    charge_number: int
    mass_amu: float
    flux_m2_s: float
    velocity_sqrt_eV: np.ndarray
    weight: np.ndarray
    phase_rad: np.ndarray | None = None
    position_m: np.ndarray | None = None
    density_model: BoundaryDensityModel | None = None
    provenance: Mapping[str, object] = field(default_factory=dict)
    density_model_2d: BoundaryDensityModel2D | None = None

    def __post_init__(self):
        if not self.name:
            raise ValueError("species name is required")
        if self.mass_amu <= 0.0 or self.flux_m2_s < 0.0:
            raise ValueError("mass must be positive and flux nonnegative")
        velocity = _readonly_array(self.velocity_sqrt_eV, (3,))
        weight = np.asarray(self.weight, dtype=float).copy()
        if weight.shape != (velocity.shape[0],) or np.any(weight < 0.0) or not np.all(np.isfinite(weight)):
            raise ValueError("weights must be finite, nonnegative, and match sample count")
        total = float(weight.sum())
        if total <= 0.0:
            raise ValueError("weights must have positive mass")
        weight /= total; weight.setflags(write=False)
        if np.any(velocity[:, 2] < 0.0):
            raise ValueError("incident normal velocity coordinate must be nonnegative")
        phase = None if self.phase_rad is None else _readonly_array(self.phase_rad)
        position = None if self.position_m is None else _readonly_array(self.position_m, (2,))
        if phase is not None and phase.shape[0] != velocity.shape[0]:
            raise ValueError("phase must match sample count")
        if position is not None and position.shape[0] != velocity.shape[0]:
            raise ValueError("position must match sample count")
        if (self.density_model_2d is not None
                and (not hasattr(self.density_model_2d, "sampling_dimension")
                     or not hasattr(self.density_model_2d, "sample_flux_velocity_2d"))):
            raise ValueError("2-D density model must provide sampling dimension and sampler")
        object.__setattr__(self, "velocity_sqrt_eV", velocity)
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "phase_rad", phase)
        object.__setattr__(self, "position_m", position)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def kinetic_energy_eV(self):
        return np.sum(self.velocity_sqrt_eV ** 2, axis=1)

    @property
    def mean_energy_eV(self):
        return float(np.dot(self.weight, self.kinetic_energy_eV))

    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
        if self.density_model is None:
            raise ValueError(f"species {self.name!r} has no continuous boundary density model")
        return self.density_model.log_flux_density(velocity_sqrt_eV, phase_rad, position_m)

    @property
    def flux_sampling_dimension(self):
        if self.density_model is None or not hasattr(self.density_model, "sampling_dimension"):
            raise ValueError(f"species {self.name!r} has no deterministic continuous-density sampler")
        return int(self.density_model.sampling_dimension)

    def sample_flux_velocity(self, unit_interval):
        if self.density_model is None or not hasattr(self.density_model, "sample_flux_velocity"):
            raise ValueError(f"species {self.name!r} has no deterministic continuous-density sampler")
        velocity = np.asarray(self.density_model.sample_flux_velocity(unit_interval), dtype=float)
        if (velocity.ndim != 2 or velocity.shape[1] != 3 or np.any(~np.isfinite(velocity))
                or np.any(velocity[:, 2] < 0.0)):
            raise RuntimeError("boundary density sampler returned invalid incident velocities")
        return velocity

    @property
    def flux_sampling_dimension_2d(self):
        model = (
            self.density_model
            if self.density_model_2d is None else self.density_model_2d)
        if model is None or not hasattr(model, "sampling_dimension"):
            raise ValueError(f"species {self.name!r} has no deterministic 2-D density sampler")
        return int(model.sampling_dimension)

    def sample_flux_velocity_2d(self, unit_interval):
        if self.density_model_2d is not None:
            velocity = np.asarray(
                self.density_model_2d.sample_flux_velocity_2d(unit_interval),
                dtype=float)
        else:
            velocity_3d = self.sample_flux_velocity(unit_interval)
            # Default projection for genuinely 3-D sources: retain the sampled
            # x/z direction while preserving total kinetic energy.
            energy = np.einsum("rc,rc->r", velocity_3d, velocity_3d)
            theta = np.arctan2(velocity_3d[:, 0], velocity_3d[:, 2])
            speed = np.sqrt(np.maximum(energy, 0.0))
            velocity = np.column_stack((
                speed * np.sin(theta), speed * np.cos(theta)))
        if (velocity.ndim != 2 or velocity.shape[1] != 2
                or np.any(~np.isfinite(velocity))
                or np.any(velocity[:, 1] < 0.0)):
            raise RuntimeError("2-D boundary density sampler returned invalid velocities")
        return velocity


@dataclass(frozen=True)
class PlasmaBoundaryState:
    species: tuple[SpeciesBoundaryState, ...]
    reference_plane_m: float
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        species = tuple(self.species)
        if not species or len({item.name for item in species}) != len(species):
            raise ValueError("boundary state requires uniquely named species")
        if not np.isfinite(self.reference_plane_m):
            raise ValueError("reference_plane_m must be finite")
        object.__setattr__(self, "species", species)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))

    def get(self, name):
        for item in self.species:
            if item.name == name:
                return item
        raise KeyError(name)

    @property
    def current_density_A_m2(self):
        return float(ECHARGE * sum(item.charge_number * item.flux_m2_s for item in self.species))


def maxwellian_electron_boundary_state(temperature_eV, flux_m2_s, *, n_transverse=5, n_normal=8,
                                        electron_name="electron", reference_plane_m=0.0):
    """Construct a deterministic electron half-Maxwellian flux quadrature.

    Gauss-Hermite integrates each tangential Maxwellian and Gauss-Laguerre integrates the exponential
    normal-energy distribution. The nodes and weights are numerical quadrature only; the physical law
    is the analytic kinetic flux density above.
    """
    temperature = float(temperature_eV)
    if temperature <= 0.0 or int(n_transverse) <= 0 or int(n_normal) <= 0:
        raise ValueError("positive temperature and quadrature orders are required")
    hermite_node, hermite_weight = np.polynomial.hermite.hermgauss(int(n_transverse))
    laguerre_node, laguerre_weight = np.polynomial.laguerre.laggauss(int(n_normal))
    ix, iy, iz = np.meshgrid(
        np.arange(hermite_node.size), np.arange(hermite_node.size),
        np.arange(laguerre_node.size), indexing="ij")
    velocity = np.column_stack((
        np.sqrt(temperature) * hermite_node[ix.ravel()],
        np.sqrt(temperature) * hermite_node[iy.ravel()],
        np.sqrt(temperature * laguerre_node[iz.ravel()]),
    ))
    weight = (hermite_weight[ix.ravel()] * hermite_weight[iy.ravel()]
              * laguerre_weight[iz.ravel()] / np.pi)
    electron = SpeciesBoundaryState(
        name=electron_name, charge_number=-1, mass_amu=5.485799e-4,
        flux_m2_s=float(flux_m2_s), velocity_sqrt_eV=velocity, weight=weight,
        density_model=MaxwellianFluxVelocityDensity(temperature),
        provenance={"model": "analytic_half_maxwellian_flux"},
    )
    return PlasmaBoundaryState(
        species=(electron,), reference_plane_m=float(reference_plane_m),
        provenance={"source": "kinetic_maxwellian"},
    )


def mixture_boundary_proposal(components, mixture_weight=None, *, name="proposal"):
    """Combine species quadratures into an exactly scored multiple-importance proposal.

    This object is numerical, not a plasma source: component weights may change estimator variance but
    cannot change the physical density scored by adjoint transport. Components must provide continuous
    densities. Phase/position labels are retained when every component supplies them.
    """
    components = tuple(components)
    if not components or any(item.density_model is None for item in components):
        raise ValueError("proposal components require continuous density models")
    if mixture_weight is None:
        mixture_weight = np.ones(len(components))
    mixture_weight = np.asarray(mixture_weight, dtype=float)
    if (mixture_weight.shape != (len(components),) or np.any(mixture_weight < 0.0)
            or mixture_weight.sum() <= 0.0):
        raise ValueError("invalid mixture weights")
    mixture_weight = mixture_weight / mixture_weight.sum()
    velocity = np.concatenate([item.velocity_sqrt_eV for item in components])
    weight = np.concatenate([
        fraction * item.weight for fraction, item in zip(mixture_weight, components)])
    phase = (np.concatenate([item.phase_rad for item in components])
             if all(item.phase_rad is not None for item in components) else None)
    position = (np.concatenate([item.position_m for item in components])
                if all(item.position_m is not None for item in components) else None)
    first = components[0]
    return SpeciesBoundaryState(
        name=name, charge_number=first.charge_number, mass_amu=first.mass_amu, flux_m2_s=1.0,
        velocity_sqrt_eV=velocity, weight=weight, phase_rad=phase, position_m=position,
        density_model=MixtureBoundaryDensity(
            tuple(item.density_model for item in components), mixture_weight),
        provenance={"role": "numerical_multiple_importance_proposal",
                    "components": tuple(item.name for item in components)},
    )


def folded_normal_tangential_proposal(template, tangent_sign, *, name=None):
    """Return a template for the normalized grazing-incidence pushforward of a source density."""
    if template.density_model is None:
        raise ValueError("grazing-incidence proposal requires a continuous source density")
    sign = int(tangent_sign)
    if sign not in (-1, 1):
        raise ValueError("tangent_sign must be -1 or +1")
    velocity = template.velocity_sqrt_eV
    transformed = np.column_stack((
        sign * velocity[:, 2], velocity[:, 1], np.abs(velocity[:, 0])))
    return SpeciesBoundaryState(
        name=(f"{template.name}-grazing-{sign:+d}" if name is None else name),
        charge_number=template.charge_number, mass_amu=template.mass_amu, flux_m2_s=1.0,
        velocity_sqrt_eV=transformed, weight=template.weight,
        density_model=FoldedNormalTangentialDensity(template.density_model, sign),
        provenance={"role": "normal_tangential_grazing_importance_proposal",
                    "source": template.name, "tangent_sign": sign},
    )


def qmc_boundary_proposal(template: SpeciesBoundaryState, log2_samples, seed=0, *, name=None):
    """Scrambled-Sobol proposal sampled from a supported analytic/tabulated density.

    The returned equal weights are numerical Monte Carlo weights. Repeated seeds give independent
    randomized replicates for error estimation. Mixtures stratify every component separately, avoiding
    random loss of a low-probability support-completion component.
    """
    level = int(log2_samples)
    if level < 0 or template.density_model is None:
        raise ValueError("nonnegative sample level and a continuous density are required")
    density = template.density_model
    proposal_name = template.name if name is None else name
    if isinstance(density, MixtureBoundaryDensity):
        components = []
        for index, component_density in enumerate(density.components):
            component = SpeciesBoundaryState(
                name=f"{proposal_name}:{index}", charge_number=template.charge_number,
                mass_amu=template.mass_amu, flux_m2_s=1.0,
                velocity_sqrt_eV=[[0.0, 0.0, 1.0]], weight=[1.0],
                density_model=component_density)
            components.append(qmc_boundary_proposal(
                component, level, seed=int(seed) + 104729 * index,
                name=component.name))
        return mixture_boundary_proposal(components, density.weight, name=proposal_name)

    if isinstance(density, FoldedNormalTangentialDensity):
        source_template = SpeciesBoundaryState(
            name=f"{proposal_name}:source", charge_number=template.charge_number,
            mass_amu=template.mass_amu, flux_m2_s=1.0,
            velocity_sqrt_eV=[[0.0, 0.0, 1.0]], weight=[1.0],
            density_model=density.source)
        source = qmc_boundary_proposal(
            source_template, level, seed=int(seed), name=source_template.name)
        velocity = source.velocity_sqrt_eV
        transformed = np.column_stack((
            density.tangent_sign * velocity[:, 2], velocity[:, 1], np.abs(velocity[:, 0])))
        return SpeciesBoundaryState(
            name=proposal_name, charge_number=template.charge_number,
            mass_amu=template.mass_amu, flux_m2_s=1.0,
            velocity_sqrt_eV=transformed, weight=source.weight,
            density_model=density,
            provenance={"role": "scrambled_sobol_grazing_proposal",
                        "log2_samples": level, "seed": int(seed),
                        "tangent_sign": density.tangent_sign})

    if not hasattr(density, "sampling_dimension") or not hasattr(density, "sample_flux_velocity"):
        raise TypeError(f"QMC sampling is not implemented for {type(density).__name__}")
    n = 2 ** level
    u = qmc.Sobol(
        int(density.sampling_dimension), scramble=True, seed=int(seed)).random_base2(level)
    velocity = density.sample_flux_velocity(u)
    return SpeciesBoundaryState(
        name=proposal_name, charge_number=template.charge_number, mass_amu=template.mass_amu,
        flux_m2_s=1.0, velocity_sqrt_eV=velocity, weight=np.ones(n), density_model=density,
        provenance={"role": "scrambled_sobol_proposal", "log2_samples": level, "seed": int(seed)},
    )


def qmc_boundary_proposal_with_auxiliary(
        template: SpeciesBoundaryState, log2_samples, auxiliary_dimension=1, seed=0, *, name=None):
    """Sample velocity and auxiliary coordinates from one joint scrambled-Sobol rule.

    Surface position is coupled to the phase-space acceptance map, so independently permuting a
    one-dimensional position stratification against a velocity QMC rule is not a joint low-discrepancy
    integration. Mixtures remain component-stratified and every component receives its own joint rule.
    """
    level = int(log2_samples); auxiliary_dimension = int(auxiliary_dimension)
    if level < 0 or auxiliary_dimension <= 0 or template.density_model is None:
        raise ValueError("require a nonnegative level, positive auxiliary dimension, and density")
    density = template.density_model
    proposal_name = template.name if name is None else name
    if isinstance(density, MixtureBoundaryDensity):
        components = []; auxiliary = []
        for index, component_density in enumerate(density.components):
            component = SpeciesBoundaryState(
                name=f"{proposal_name}:{index}", charge_number=template.charge_number,
                mass_amu=template.mass_amu, flux_m2_s=1.0,
                velocity_sqrt_eV=[[0.0, 0.0, 1.0]], weight=[1.0],
                density_model=component_density)
            sampled, component_auxiliary = qmc_boundary_proposal_with_auxiliary(
                component, level, auxiliary_dimension,
                seed=int(seed) + 104729 * index, name=component.name)
            components.append(sampled); auxiliary.append(component_auxiliary)
        return (mixture_boundary_proposal(components, density.weight, name=proposal_name),
                np.concatenate(auxiliary, axis=0))
    if not hasattr(density, "sampling_dimension") or not hasattr(density, "sample_flux_velocity"):
        raise TypeError(f"QMC sampling is not implemented for {type(density).__name__}")
    physical_dimension = int(density.sampling_dimension)
    unit = qmc.Sobol(
        physical_dimension + auxiliary_dimension, scramble=True,
        seed=int(seed)).random_base2(level)
    velocity = density.sample_flux_velocity(unit[:, :physical_dimension])
    proposal = SpeciesBoundaryState(
        name=proposal_name, charge_number=template.charge_number, mass_amu=template.mass_amu,
        flux_m2_s=1.0, velocity_sqrt_eV=velocity, weight=np.ones(2 ** level),
        density_model=density,
        provenance={"role": "joint_scrambled_sobol_proposal", "log2_samples": level,
                    "seed": int(seed), "auxiliary_dimension": auxiliary_dimension})
    return proposal, unit[:, physical_dimension:].copy()


def collisionless_sheath_boundary_state(sheath: CollisionlessRFSheath, flux_m2_s, *, n_phase=256,
                                         ion_name="ion", reference_plane_m=0.0,
                                         tangential_temperature_eV=None, n_transverse=3,
                                         normal_energy_bins=64, density_phase_count=None):
    """Construct the common boundary state from the finite-transit collisionless sheath."""
    phase = 2.0 * np.pi * (np.arange(int(n_phase)) + 0.5) / int(n_phase)
    energy = sheath.ion_impact_energies(phase)
    density_model = None
    if tangential_temperature_eV is None:
        velocity = np.zeros((phase.size, 3)); velocity[:, 2] = np.sqrt(energy)
        weight = np.ones(phase.size)
        sample_phase = phase
    else:
        nodes, gh_weight = np.polynomial.hermite.hermgauss(int(n_transverse))
        transverse = np.sqrt(float(tangential_temperature_eV)) * nodes
        transverse_weight = gh_weight / np.sqrt(np.pi)
        ex, ey, ep = np.meshgrid(np.arange(nodes.size), np.arange(nodes.size),
                                 np.arange(phase.size), indexing="ij")
        velocity = np.column_stack((transverse[ex.ravel()], transverse[ey.ravel()],
                                    np.sqrt(energy[ep.ravel()])))
        weight = (transverse_weight[ex.ravel()] * transverse_weight[ey.ravel()])
        sample_phase = phase[ep.ravel()]
        # The transport quadrature and continuous-density representation have independent accuracy
        # requirements. Binning only the output phase nodes creates artificial zero-density holes in
        # the continuous pushforward of uniform RF phase, which causes catastrophic adjoint weights
        # after even small electrostatic energy shifts. Densely integrate the same sheath map instead.
        if density_phase_count is None:
            density_phase_count = max(4096, 64 * int(normal_energy_bins))
        density_phase = (2.0 * np.pi
                         * (np.arange(int(density_phase_count)) + 0.5)
                         / int(density_phase_count))
        density_energy = sheath.ion_impact_energies(density_phase)
        span = max(float(np.ptp(density_energy)), 1e-6)
        lo = max(0.0, float(density_energy.min()) - 0.01 * span)
        hi = float(density_energy.max()) + 0.01 * span
        edges = np.linspace(lo, hi, int(normal_energy_bins) + 1)
        mass, _ = np.histogram(density_energy, bins=edges)
        density_model = IonEnergyTransverseMaxwellianDensity(
            edges, mass.astype(float), float(tangential_temperature_eV))
    ion = SpeciesBoundaryState(
        name=ion_name, charge_number=1, mass_amu=sheath.ion_mass_amu,
        flux_m2_s=float(flux_m2_s), velocity_sqrt_eV=velocity,
        weight=weight, phase_rad=sample_phase, density_model=density_model,
        provenance={"model": "collisionless_finite_transit_child_sheath"},
    )
    return PlasmaBoundaryState(
        species=(ion,), reference_plane_m=float(reference_plane_m),
        provenance={"source": "CollisionlessRFSheath"},
    )


def instantaneous_sinusoidal_ion_boundary_state(V_dc, V_rf, Te_eV, ion_mass_amu, flux_m2_s, *,
                                                 n_phase=256, ion_name="ion", reference_plane_m=0.0):
    """Named instantaneous/zero-transit limiting constructor; not universal production physics."""
    phase = 2.0 * np.pi * (np.arange(int(n_phase)) + 0.5) / int(n_phase)
    energy = 0.5 * float(Te_eV) + float(V_dc) + float(V_rf) * np.sin(phase)
    if np.any(energy < 0.0):
        raise ValueError("instantaneous sheath energy became negative")
    velocity = np.zeros((phase.size, 3)); velocity[:, 2] = np.sqrt(energy)
    ion = SpeciesBoundaryState(
        name=ion_name, charge_number=1, mass_amu=float(ion_mass_amu), flux_m2_s=float(flux_m2_s),
        velocity_sqrt_eV=velocity, weight=np.ones(phase.size), phase_rad=phase,
        provenance={"model": "instantaneous_sinusoidal_limit"},
    )
    return PlasmaBoundaryState(species=(ion,), reference_plane_m=float(reference_plane_m),
                               provenance={"source": "analytic_limit"})
