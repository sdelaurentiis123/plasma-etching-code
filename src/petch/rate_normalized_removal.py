"""Explicitly conditional surface removal normalized by a blanket rate.

This mechanism is a bridge between verified feature transport and an observed
or cross-machine blanket etch rate.  It is intentionally *not* a microscopic
surface law: the supplied blanket velocity sets the normal-incidence scale,
while the deterministic feature transport supplies the local ion-dose ratio.
Keeping this bridge as a separate mechanism prevents a process rate from being
disguised as a universal sputter yield.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .physical_sputtering import PhysicalSputterState
from .surface_exchange import unresolved_surface_exchange
from .surface_kinetics import (
    MechanismValidity,
    ParameterEvidence,
    SurfaceFluxes,
)


@dataclass(frozen=True)
class RateNormalizedRemovalParameters:
    """Inputs for a conditional ion-dose-to-removal conversion.

    ``blanket_removal_velocity_m_s`` is the velocity that the declared
    projectile flux would produce on an unobstructed horizontal witness.
    Feature faces inherit only the local dose ratio.  Energy and incidence
    response beyond what is already represented in transport are omitted and
    reported in the validity record.
    """

    material_name: str
    material_inventory_name: str
    projectile_species: tuple[str, ...]
    reference_projectile_flux_m2_s: float
    blanket_removal_velocity_m_s: float
    bulk_material_unit_density_m3: float
    evidence: Mapping[str, ParameterEvidence]
    declared_inert_neutral_species: tuple[str, ...] = ()
    known_omissions: tuple[str, ...] = (
        "blanket rate is transferred without a species-resolved surface reaction law",
        "ion-energy and incidence-angle yield response beyond geometric projected flux is omitted",
        "neutral-radical, polymer, passivation, charging, and redeposition responses are omitted",
        "removed-product identity and emission distribution are unresolved",
    )

    def __post_init__(self):
        projectiles = tuple(self.projectile_species)
        if (
            not self.material_name
            or not self.material_inventory_name
            or not projectiles
            or any(not item for item in projectiles)
            or len(set(projectiles)) != len(projectiles)
            or not np.isfinite(self.reference_projectile_flux_m2_s)
            or self.reference_projectile_flux_m2_s <= 0.0
            or not np.isfinite(self.blanket_removal_velocity_m_s)
            or self.blanket_removal_velocity_m_s <= 0.0
            or not np.isfinite(self.bulk_material_unit_density_m3)
            or self.bulk_material_unit_density_m3 <= 0.0
        ):
            raise ValueError("invalid rate-normalized removal parameters")
        evidence = dict(self.evidence)
        inert = tuple(str(name) for name in self.declared_inert_neutral_species)
        required = {
            "reference_projectile_flux_m2_s",
            "blanket_removal_velocity_m_s",
            "bulk_material_unit_density_m3",
        }
        if set(evidence) != required or any(
            not isinstance(item, ParameterEvidence) for item in evidence.values()
        ):
            raise ValueError(
                "rate-normalized removal evidence must cover every physical input"
            )
        if (
            any(not name for name in inert)
            or len(set(inert)) != len(inert)
            or set(inert) & set(projectiles)
        ):
            raise ValueError("invalid declared inert neutral species")
        object.__setattr__(self, "projectile_species", projectiles)
        object.__setattr__(self, "declared_inert_neutral_species", inert)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))


@dataclass(frozen=True)
class RateNormalizedRemovalStepResult:
    state: PhysicalSputterState
    etch_velocity_m_s: np.ndarray
    removed_material_units_m2: np.ndarray
    material_exchange: object
    product_populations: tuple[()] = ()
    validity: MechanismValidity | None = None

    def __post_init__(self):
        velocity = np.asarray(self.etch_velocity_m_s, dtype=float).copy()
        removed = np.asarray(self.removed_material_units_m2, dtype=float).copy()
        if (
            not isinstance(self.state, PhysicalSputterState)
            or velocity.shape != removed.shape
            or np.any(~np.isfinite(velocity))
            or np.any(velocity < 0.0)
            or np.any(~np.isfinite(removed))
            or np.any(removed < 0.0)
            or not isinstance(self.validity, MechanismValidity)
        ):
            raise ValueError("invalid rate-normalized removal result")
        velocity.setflags(write=False)
        removed.setflags(write=False)
        object.__setattr__(self, "etch_velocity_m_s", velocity)
        object.__setattr__(self, "removed_material_units_m2", removed)


class RateNormalizedRemovalMechanism:
    """Map local deterministic projectile dose to a conditional etch speed."""

    def __init__(self, parameters: RateNormalizedRemovalParameters):
        if not isinstance(parameters, RateNormalizedRemovalParameters):
            raise TypeError("parameters must be RateNormalizedRemovalParameters")
        self.parameters = parameters

    @staticmethod
    def initial_state(shape=()):
        return PhysicalSputterState.bare(shape)

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        unsupported_neutral = tuple(sorted(
            name
            for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in par.declared_inert_neutral_species
            and np.any(np.asarray(value) > 0.0)
        ))
        unsupported_energetic = tuple(sorted({
            population.name
            for population in fluxes.energetic_fluxes
            if population.name not in par.projectile_species
            and np.any(np.asarray(population.flux_m2_s) > 0.0)
        }))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no rate-normalized removal channel")
        nonpredictive = tuple(sorted(
            name
            for name, evidence in par.evidence.items()
            if not evidence.supports_prediction_within_declared_domain
        ))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=par.known_omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive,
        )

    def advance(self, state, fluxes: SurfaceFluxes, duration_s: float, *, strict=True):
        if not isinstance(state, PhysicalSputterState):
            raise TypeError("rate-normalized removal requires PhysicalSputterState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons)
            )
        shape = state.removed_material_units_m2.shape
        local_flux = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name in self.parameters.projectile_species:
                local_flux = local_flux + np.broadcast_to(
                    np.asarray(population.flux_m2_s, dtype=float), shape
                )
        velocity = (
            self.parameters.blanket_removal_velocity_m_s
            * local_flux
            / self.parameters.reference_projectile_flux_m2_s
        )
        removed = (
            velocity
            * self.parameters.bulk_material_unit_density_m3
            * float(duration_s)
        )
        updated = PhysicalSputterState(
            state.removed_material_units_m2 + removed
        )
        exchange = unresolved_surface_exchange(
            removed_units_m2={
                self.parameters.material_inventory_name: removed,
            },
            limitations=(
                "rate-normalized removal does not identify emitted products",
            ),
        )
        return RateNormalizedRemovalStepResult(
            state=updated,
            etch_velocity_m_s=velocity,
            removed_material_units_m2=removed,
            material_exchange=exchange,
            product_populations=(),
            validity=validity,
        )
