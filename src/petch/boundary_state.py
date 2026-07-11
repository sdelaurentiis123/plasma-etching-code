"""Unified plasma-to-feature boundary state.

Every analytic source, sheath model, reactor solver, diagnostic reconstruction, or learned surrogate must
produce this representation. Transport engines consume it without knowing how it was generated.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping, Protocol

import numpy as np

from .sheath import CollisionlessRFSheath, ECHARGE


class BoundaryDensityModel(Protocol):
    """Normalized incident flux-density evaluator used by adjoint transport."""
    def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None): ...


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


def collisionless_sheath_boundary_state(sheath: CollisionlessRFSheath, flux_m2_s, *, n_phase=256,
                                         ion_name="ion", reference_plane_m=0.0,
                                         tangential_temperature_eV=None, n_transverse=3,
                                         normal_energy_bins=64):
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
        span = max(float(np.ptp(energy)), 1e-6)
        lo = max(0.0, float(energy.min()) - 0.01 * span)
        hi = float(energy.max()) + 0.01 * span
        edges = np.linspace(lo, hi, int(normal_energy_bins) + 1)
        mass, _ = np.histogram(energy, bins=edges)
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
