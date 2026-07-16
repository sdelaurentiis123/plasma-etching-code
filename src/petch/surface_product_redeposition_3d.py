"""Conservative same-material redeposition for explicitly resolved surface products.

The transport operator already knows how an emitted diffuse population moves between triangles.  This
module adds the deliberately narrow material closure needed to feed a *declared* sticking event back to
the moving interface.  It never guesses a product identity, sticking probability, density, or film
material.  Version 1 permits only same-material growth: a product declared as material ``m`` may react
only on faces already owned by ``m``.  Cross-material coatings require a new level-set material layer
and are refused instead of being hidden inside an effective etch rate.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .neutral_radiosity_3d import (
    DiffuseFormFactors3D,
    DiffuseSurfaceEmissionSolve3D,
    transport_surface_product_population_3d,
)
from .surface_exchange import SurfaceProductPopulation


_PARAMETERS = frozenset({
    "sticking_probability_by_material", "bulk_material_unit_density_m3",
})


@dataclass(frozen=True)
class SurfaceProductRedepositionLaw3D:
    """One population's bounded sticking and same-material volume-conversion law."""

    population_name: str
    deposited_material_id: int
    sticking_probability_by_material: Mapping[int, float]
    bulk_material_unit_density_m3: float
    parameter_sources: Mapping[str, str]
    parameter_bounds: Mapping[str, tuple[float, float]]

    def __post_init__(self):
        probability = {
            int(material): float(value)
            for material, value in dict(self.sticking_probability_by_material).items()}
        sources = dict(self.parameter_sources)
        bounds = {
            str(name): tuple(float(value) for value in supplied)
            for name, supplied in dict(self.parameter_bounds).items()}
        material_id = int(self.deposited_material_id)
        density = float(self.bulk_material_unit_density_m3)
        if (not isinstance(self.population_name, str) or not self.population_name
                or material_id <= 0 or not probability
                or any(key <= 0 for key in probability)
                or any(not np.isfinite(value) or value < 0.0 or value > 1.0
                       for value in probability.values())
                or not np.isfinite(density) or density <= 0.0
                or set(sources) != _PARAMETERS or set(bounds) != _PARAMETERS
                or any(not isinstance(value, str) or not value for value in sources.values())
                or any(len(value) != 2 or not np.all(np.isfinite(value))
                       or value[0] > value[1] for value in bounds.values())):
            raise ValueError("invalid surface-product redeposition law")
        probability_bounds = bounds["sticking_probability_by_material"]
        density_bounds = bounds["bulk_material_unit_density_m3"]
        if (probability_bounds[0] < 0.0 or probability_bounds[1] > 1.0
                or any(not probability_bounds[0] <= value <= probability_bounds[1]
                       for value in probability.values())
                or not density_bounds[0] <= density <= density_bounds[1]
                or density_bounds[0] <= 0.0):
            raise ValueError("redeposition parameters must lie inside their declared bounds")
        object.__setattr__(self, "deposited_material_id", material_id)
        object.__setattr__(self, "sticking_probability_by_material", MappingProxyType(probability))
        object.__setattr__(self, "bulk_material_unit_density_m3", density)
        object.__setattr__(self, "parameter_sources", MappingProxyType(sources))
        object.__setattr__(self, "parameter_bounds", MappingProxyType(bounds))

    @property
    def provenance(self):
        return MappingProxyType(dict(
            model="same-material-diffuse-redeposition-v1",
            parameters=dict(
                sticking_probability_by_material={
                    str(key): value
                    for key, value in self.sticking_probability_by_material.items()},
                bulk_material_unit_density_m3=self.bulk_material_unit_density_m3,
                deposited_material_id=self.deposited_material_id),
            sources=dict(self.parameter_sources),
            bounds={name: list(value) for name, value in self.parameter_bounds.items()}))


@dataclass(frozen=True)
class SurfaceProductRedepositionContract3D:
    """Exact set of product laws enabled for one feature run."""

    laws: tuple[SurfaceProductRedepositionLaw3D, ...]

    def __post_init__(self):
        laws = tuple(self.laws)
        if (not laws or any(not isinstance(item, SurfaceProductRedepositionLaw3D)
                            for item in laws)
                or len({item.population_name for item in laws}) != len(laws)):
            raise ValueError("redeposition contract requires unique product laws")
        object.__setattr__(self, "laws", laws)

    @property
    def by_population(self):
        return MappingProxyType({item.population_name: item for item in self.laws})

    @property
    def provenance(self):
        return MappingProxyType(dict(
            model="surface-product-redeposition-3d-v1",
            closure="same-material growth only; cross-material coatings refused",
            laws={item.population_name: dict(item.provenance) for item in self.laws}))


@dataclass(frozen=True)
class SurfaceProductRedeposition3DResult:
    """Conserved emitted/deposited/escaped material and signed growth contribution."""

    transport_by_population: Mapping[str, DiffuseSurfaceEmissionSolve3D]
    deposited_units_m2: Mapping[str, np.ndarray]
    normal_growth_velocity_m_s: np.ndarray
    emitted_material_units_s: Mapping[str, float]
    deposited_material_units_s: Mapping[str, float]
    escaped_material_units_s: Mapping[str, float]
    maximum_relative_balance_error: float
    contract_provenance: Mapping[str, object]

    def __post_init__(self):
        transport = dict(self.transport_by_population)
        deposited = {}
        for name, supplied in dict(self.deposited_units_m2).items():
            value = np.asarray(supplied, dtype=float).copy()
            if not name or value.ndim != 1 or np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid redeposited material inventory")
            value.setflags(write=False)
            deposited[name] = value
        growth = np.asarray(self.normal_growth_velocity_m_s, dtype=float).copy()
        rates = [dict(self.emitted_material_units_s), dict(self.deposited_material_units_s),
                 dict(self.escaped_material_units_s)]
        if (not transport or any(not isinstance(value, DiffuseSurfaceEmissionSolve3D)
                                 for value in transport.values())
                or not deposited or growth.ndim != 1 or np.any(~np.isfinite(growth))
                or np.any(growth < 0.0)
                or any(set(value) != set(transport) for value in rates)
                or any(not np.isfinite(item) or item < 0.0
                       for value in rates for item in value.values())
                or not np.isfinite(self.maximum_relative_balance_error)
                or self.maximum_relative_balance_error < 0.0):
            raise ValueError("invalid surface-product redeposition result")
        growth.setflags(write=False)
        object.__setattr__(self, "transport_by_population", MappingProxyType(transport))
        object.__setattr__(self, "deposited_units_m2", MappingProxyType(deposited))
        object.__setattr__(self, "normal_growth_velocity_m_s", growth)
        object.__setattr__(self, "emitted_material_units_s", MappingProxyType(rates[0]))
        object.__setattr__(self, "deposited_material_units_s", MappingProxyType(rates[1]))
        object.__setattr__(self, "escaped_material_units_s", MappingProxyType(rates[2]))
        object.__setattr__(self, "contract_provenance", MappingProxyType(
            dict(self.contract_provenance)))


def transport_surface_product_redeposition_3d(
        populations, duration_s, face_area_m2, form_factors: DiffuseFormFactors3D,
        face_material_id, evolving_face_mask,
        contract: SurfaceProductRedepositionContract3D, *,
        relative_tolerance=1e-10, maximum_iterations=500):
    """Transport every declared product and convert same-material capture to interface growth."""
    populations = tuple(populations)
    area = np.asarray(face_area_m2, dtype=float)
    material = np.asarray(face_material_id, dtype=int)
    evolving = np.asarray(evolving_face_mask, dtype=bool)
    if (not isinstance(contract, SurfaceProductRedepositionContract3D)
            or not isinstance(form_factors, DiffuseFormFactors3D)
            or form_factors.face_count != area.size
            or area.ndim != 1 or material.shape != area.shape or evolving.shape != area.shape
            or np.any(~np.isfinite(area)) or np.any(area <= 0.0) or np.any(material <= 0)
            or not np.isfinite(duration_s) or duration_s <= 0.0
            or any(not isinstance(item, SurfaceProductPopulation) for item in populations)
            or len({item.name for item in populations}) != len(populations)):
        raise ValueError("invalid surface-product redeposition inputs")
    population_by_name = {item.name: item for item in populations}
    laws = contract.by_population
    if set(population_by_name) != set(laws):
        raise ValueError(
            "redeposition contract must cover every and only emitted surface-product population")
    unique_material = set(int(value) for value in np.unique(material))
    deposited = {}
    growth = np.zeros(area.shape)
    solutions = {}; emitted_rates = {}; deposited_rates = {}; escaped_rates = {}
    maximum_balance = 0.0
    for name in sorted(population_by_name):
        population = population_by_name[name]
        law = laws[name]
        if not unique_material.issubset(law.sticking_probability_by_material):
            missing = unique_material - set(law.sticking_probability_by_material)
            raise ValueError(
                f"redeposition law {name!r} lacks target materials {sorted(missing)}")
        reaction = np.asarray([
            law.sticking_probability_by_material[int(value)] for value in material], dtype=float)
        solution = transport_surface_product_population_3d(
            population, float(duration_s), area, form_factors, reaction,
            relative_tolerance=relative_tolerance, maximum_iterations=maximum_iterations)
        reacted_units = (solution.reacted_flux_m2_s * float(duration_s)
                         * population.material_units_per_particle)
        incompatible = (material != law.deposited_material_id) & (reacted_units > 0.0)
        if np.any(incompatible):
            raise ValueError(
                f"product {name!r} would form a cross-material coating; v1 requires an explicit "
                "new material layer")
        if np.any((~evolving) & (reacted_units > 0.0)):
            raise ValueError(
                f"product {name!r} deposits on a pinned surface; include that material in the "
                "evolving mechanism router")
        inventory = population.source_inventory
        deposited[inventory] = deposited.get(inventory, np.zeros(area.shape)) + reacted_units
        growth += reacted_units / law.bulk_material_unit_density_m3 / float(duration_s)
        emitted = solution.emitted_rate_s * population.material_units_per_particle
        captured = solution.reacted_rate_s * population.material_units_per_particle
        escaped = ((solution.escaped_without_impact_rate_s
                    + solution.escaped_after_reflection_rate_s)
                   * population.material_units_per_particle)
        scale = max(emitted, np.finfo(float).tiny)
        balance = abs(emitted - captured - escaped) / scale
        solutions[name] = solution
        emitted_rates[name] = emitted
        deposited_rates[name] = captured
        escaped_rates[name] = escaped
        maximum_balance = max(maximum_balance, balance)
    if maximum_balance > max(20.0 * float(relative_tolerance), 2e-12):
        raise RuntimeError(
            f"surface-product material balance failed: {maximum_balance:.3e}")
    return SurfaceProductRedeposition3DResult(
        solutions, deposited, growth, emitted_rates, deposited_rates, escaped_rates,
        maximum_balance, contract.provenance)
