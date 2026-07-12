"""Material-conserving physical sputtering mechanism for the common feature engine."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import (
    SurfaceMaterialExchange, SurfaceProductPopulation, validate_surface_product_routing,
)
from .surface_kinetics import EnergeticYield, MechanismValidity, ParameterEvidence, SurfaceFluxes


@dataclass(frozen=True)
class PhysicalSputterState:
    removed_material_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        value = np.asarray(self.removed_material_units_m2, dtype=float).copy()
        if np.any(~np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError("removed material inventory must be finite and nonnegative")
        value.setflags(write=False)
        object.__setattr__(self, "removed_material_units_m2", value)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.zeros(shape))

    def conservative_surface_fields(self):
        return {"removed_material_units_m2": self.removed_material_units_m2}

    def conservative_surface_upper_bounds(self):
        return {"removed_material_units_m2": None}

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != {"removed_material_units_m2"}:
            raise ValueError("physical-sputter remap fields do not match its state contract")
        return type(self)(fields["removed_material_units_m2"])


@dataclass(frozen=True)
class PhysicalSputterParameters:
    material_name: str
    material_inventory_name: str
    projectile_species: tuple[str, ...]
    bulk_material_unit_density_m3: float
    sputter_yield: EnergeticYield
    emitted_product_name: str
    emitted_product_mass_amu: float
    emitted_material_units_per_particle: float
    emission_angular_model: str
    emission_energy_model: str
    emission_energy_parameters: Mapping[str, float]
    evidence: Mapping[str, ParameterEvidence]
    known_omissions: tuple[str, ...] = (
        "collision-cascade damage and implantation state are not represented",
        "emitted energy affects downstream interactions only when their material law consumes it",
        "redeposited material state and resputtering require a material mechanism",
    )

    def __post_init__(self):
        projectiles = tuple(self.projectile_species)
        if (not self.material_name or not self.material_inventory_name or not projectiles
                or any(not item for item in projectiles) or len(set(projectiles)) != len(projectiles)
                or not np.isfinite(self.bulk_material_unit_density_m3)
                or self.bulk_material_unit_density_m3 <= 0.0
                or not isinstance(self.sputter_yield, EnergeticYield)
                or not self.emitted_product_name
                or not np.isfinite(self.emitted_product_mass_amu)
                or self.emitted_product_mass_amu <= 0.0
                or not np.isfinite(self.emitted_material_units_per_particle)
                or self.emitted_material_units_per_particle <= 0.0
                or self.emission_angular_model not in {"diffuse_cosine"}):
            raise ValueError("invalid physical-sputter parameters")
        if not self.emission_energy_model:
            raise ValueError("physical sputtering requires an emitted energy model")
        energy_parameters = dict(self.emission_energy_parameters)
        # Reuse the product contract as the authoritative energy-law validator.
        SurfaceProductPopulation(
            self.emitted_product_name, self.material_inventory_name, np.asarray(0.0),
            self.emitted_material_units_per_particle, self.emitted_product_mass_amu,
            angular_model=self.emission_angular_model,
            energy_model=self.emission_energy_model,
            energy_parameters=energy_parameters)
        evidence = dict(self.evidence)
        required = {
            "bulk_material_unit_density_m3", "sputter_yield",
            "emitted_product_mass_amu", "emission_angular_model", "emission_energy_model",
        }
        if set(evidence) != required or any(
                not isinstance(item, ParameterEvidence) for item in evidence.values()):
            raise ValueError("physical-sputter evidence must cover every physical input")
        object.__setattr__(self, "projectile_species", projectiles)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "emission_energy_parameters", MappingProxyType(energy_parameters))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))


@dataclass(frozen=True)
class PhysicalSputterStepResult:
    state: PhysicalSputterState
    etch_velocity_m_s: np.ndarray
    removed_material_units_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity

    def __post_init__(self):
        for name in ("etch_velocity_m_s", "removed_material_units_m2"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            value.setflags(write=False); object.__setattr__(self, name, value)
        products = validate_surface_product_routing(
            self.material_exchange, self.product_populations)
        object.__setattr__(self, "product_populations", products)


class PhysicalSputterMechanism:
    """Physical sputtering with explicit target removal and outgoing material populations."""

    def __init__(self, parameters: PhysicalSputterParameters):
        if not isinstance(parameters, PhysicalSputterParameters):
            raise TypeError("parameters must be PhysicalSputterParameters")
        self.parameters = parameters

    @staticmethod
    def initial_state(shape=()):
        return PhysicalSputterState.bare(shape)

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in par.projectile_species
            and np.any(np.asarray(population.flux_m2_s) > 0.0)}))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no physical-sputter channel")
        nonpredictive = tuple(sorted(
            name for name, evidence in par.evidence.items()
            if not evidence.supports_prediction_within_declared_domain))
        return MechanismValidity(
            within_declared_scope=not reasons, reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=par.known_omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive)

    def advance(self, state, fluxes: SurfaceFluxes, duration_s: float, *, strict=True):
        if not isinstance(state, PhysicalSputterState):
            raise TypeError("physical sputtering requires PhysicalSputterState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))
        shape = state.removed_material_units_m2.shape
        removal_rate = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name in self.parameters.projectile_species:
                removal_rate = removal_rate + np.broadcast_to(
                    population.yield_rate_m2_s(self.parameters.sputter_yield), shape)
        removed = removal_rate * float(duration_s)
        updated = PhysicalSputterState(state.removed_material_units_m2 + removed)
        inventory = self.parameters.material_inventory_name
        exchange = SurfaceMaterialExchange(
            removed_units_m2={inventory: removed}, outgoing_units_m2={inventory: removed},
            unresolved_units_m2={}, deposited_units_m2={},
            known_limitations=(
                "outgoing physical-sputter material is not redeposited unless product transport is enabled",
            ))
        product = SurfaceProductPopulation(
            name=self.parameters.emitted_product_name, source_inventory=inventory,
            integrated_particle_count_m2=(
                removed / self.parameters.emitted_material_units_per_particle),
            material_units_per_particle=self.parameters.emitted_material_units_per_particle,
            mass_amu=self.parameters.emitted_product_mass_amu,
            angular_model=self.parameters.emission_angular_model,
            energy_model=self.parameters.emission_energy_model,
            energy_parameters=self.parameters.emission_energy_parameters,
            provenance={
                "angular_source": self.parameters.evidence["emission_angular_model"].source,
                "energy_source": self.parameters.evidence["emission_energy_model"].source,
            },
            relative_standard_uncertainty=max(
                value for value in (
                    self.parameters.evidence["emission_angular_model"].relative_standard_uncertainty,
                    self.parameters.evidence["emission_energy_model"].relative_standard_uncertainty,
                ) if value is not None) if any(
                    value is not None for value in (
                        self.parameters.evidence["emission_angular_model"].relative_standard_uncertainty,
                        self.parameters.evidence["emission_energy_model"].relative_standard_uncertainty,
                    )) else None)
        return PhysicalSputterStepResult(
            state=updated,
            etch_velocity_m_s=removal_rate / self.parameters.bulk_material_unit_density_m3,
            removed_material_units_m2=removed, material_exchange=exchange,
            product_populations=(product,), validity=validity)
