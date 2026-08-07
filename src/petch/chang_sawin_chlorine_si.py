"""Species- and atom-balanced Chang--Sawin Ar+/chlorine etching of poly-Si.

This module implements the deliberately small steady surface model in J. P.
Chang's MIT thesis (1998), Chapter 3, Tables 3.3--3.4 and Eqs. (3.9)--(3.11).
At a fixed ion energy and normal incidence,

``theta = s R / (s R + 4 beta)``

``Y = Y0 (1 - theta) + beta theta``

where ``R`` is neutral-particle flux divided by Ar+ flux.  The coefficient
``s`` is retained under the source's careful name, *surface chlorination
coefficient*.  It is not relabelled as a universal sticking probability.

The implementation makes the source's implied elemental bookkeeping explicit:
the chemical branch removes one Si and consumes four surface Cl atoms per
assumed SiCl4 product; the physical branch emits elemental Si.  Both product
energy-angle distributions are unknown, so neither is transport-ready.

This is a mechanistic beam-regressed closure, not a first-principles potential.
In particular, the SiCl4 branch is the source model's low-energy assumption,
and the molecular-Cl2 angular response was not measured.  Validity reporting
keeps those limitations visible.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import (
    SurfaceMaterialExchange,
    SurfaceProductPopulation,
    validate_surface_product_routing,
)
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    MechanismValidity,
    ParameterEvidence,
    SurfaceFluxes,
)


def _chang_sawin_chemical_angular_factor(cosine_incidence):
    """Minimal source-constrained interpolation of the measured angular curve.

    The measurements are flat through 40 degrees, are approximately 0.7 and
    0.5 of normal at 60 and 70 degrees, and assume zero at grazing incidence.
    The MCFPM-lineage class-2 interpolation introduces no extra fitted
    coefficient: unity through 45 degrees, then projected-flux roll-off.
    """
    cosine = np.asarray(cosine_incidence, dtype=float)
    if (np.any(~np.isfinite(cosine))
            or np.any((cosine < 0.0) | (cosine > 1.0))):
        raise ValueError("incidence cosines must lie in [0, 1]")
    return np.minimum(cosine / np.cos(np.pi / 4.0), 1.0)


@dataclass(frozen=True)
class ChangSawinArClSiState:
    """Cumulative elemental inventories carried across geometry remaps."""

    removed_si_atoms_m2: np.ndarray | float = 0.0
    consumed_chlorine_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        values = np.broadcast_arrays(
            np.asarray(self.removed_si_atoms_m2, dtype=float),
            np.asarray(self.consumed_chlorine_atoms_m2, dtype=float),
        )
        for name, supplied in zip(
                ("removed_si_atoms_m2", "consumed_chlorine_atoms_m2"), values):
            value = np.array(supplied, copy=True)
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("chlorine/Si inventories must be finite and nonnegative")
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.zeros(shape), np.zeros(shape))

    def conservative_surface_fields(self):
        return {
            "removed_si_atoms_m2": self.removed_si_atoms_m2,
            "consumed_chlorine_atoms_m2": self.consumed_chlorine_atoms_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {
            "removed_si_atoms_m2": None,
            "consumed_chlorine_atoms_m2": None,
        }

    def surface_field_remap_modes(self):
        return {
            "removed_si_atoms_m2": "conservative",
            "consumed_chlorine_atoms_m2": "conservative",
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("Chang--Sawin Si remap fields do not match its state")
        return type(self)(
            fields["removed_si_atoms_m2"],
            fields["consumed_chlorine_atoms_m2"],
        )


@dataclass(frozen=True)
class ChangSawinArClSiParameters:
    """One source-fixed Ar+/neutral-chlorine beam condition."""

    ion_species: str
    neutral_species: str
    ion_energy_eV: float
    neutral_chlorine_atoms_per_particle: int
    surface_chlorination_coefficient: float
    physical_sputter_yield_si_per_ion: float
    ion_enhanced_yield_si_per_ion_at_full_chlorination: float
    product_chlorine_atoms_per_si: int
    bulk_si_atom_density_m3: float
    maximum_measured_incidence_angle_deg: float
    angular_response_measured_for_neutral: bool
    evidence: Mapping[str, ParameterEvidence]
    energy_tolerance_eV: float = 1e-6

    def __post_init__(self):
        evidence = dict(self.evidence)
        required = {
            "surface_chlorination_coefficient",
            "physical_sputter_yield_si_per_ion",
            "ion_enhanced_yield_si_per_ion_at_full_chlorination",
            "product_chlorine_atoms_per_si",
            "bulk_si_atom_density_m3",
            "angular_response",
        }
        if (not self.ion_species or not self.neutral_species
                or not np.isfinite(self.ion_energy_eV)
                or self.ion_energy_eV <= 0.0
                or int(self.neutral_chlorine_atoms_per_particle)
                != self.neutral_chlorine_atoms_per_particle
                or self.neutral_chlorine_atoms_per_particle not in {1, 2}
                or not np.isfinite(self.surface_chlorination_coefficient)
                or self.surface_chlorination_coefficient < 0.0
                or self.surface_chlorination_coefficient
                > self.neutral_chlorine_atoms_per_particle
                or not np.isfinite(self.physical_sputter_yield_si_per_ion)
                or self.physical_sputter_yield_si_per_ion < 0.0
                or not np.isfinite(
                    self.ion_enhanced_yield_si_per_ion_at_full_chlorination)
                or self.ion_enhanced_yield_si_per_ion_at_full_chlorination < 0.0
                or int(self.product_chlorine_atoms_per_si)
                != self.product_chlorine_atoms_per_si
                or self.product_chlorine_atoms_per_si <= 0
                or not np.isfinite(self.bulk_si_atom_density_m3)
                or self.bulk_si_atom_density_m3 <= 0.0
                or not np.isfinite(self.maximum_measured_incidence_angle_deg)
                or not 0.0 <= self.maximum_measured_incidence_angle_deg < 90.0
                or not isinstance(self.angular_response_measured_for_neutral, bool)
                or not np.isfinite(self.energy_tolerance_eV)
                or self.energy_tolerance_eV < 0.0
                or set(evidence) != required
                or any(not isinstance(item, ParameterEvidence)
                       for item in evidence.values())):
            raise ValueError("invalid Chang--Sawin Ar+/chlorine Si parameters")
        object.__setattr__(
            self, "neutral_chlorine_atoms_per_particle",
            int(self.neutral_chlorine_atoms_per_particle))
        object.__setattr__(
            self, "product_chlorine_atoms_per_si",
            int(self.product_chlorine_atoms_per_si))
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @classmethod
    def molecular_chlorine_100eV(cls):
        """Table-3.4 Ar+ / Cl2 card, regressed only at 100 eV."""
        thesis = (
            "chang-thesis, "
            "Chapter 3, Tables 3.3--3.4 and Eqs. 3.9--3.11"
        )
        angular = (
            "chang-thesis, Chapter 3, Figures 3.7--3.8: 100 eV Ar+/atomic-Cl "
            "angular response; class-2 interpolation checked against Chapter 5 "
            "Figure 5.9"
        )
        return cls(
            ion_species="Ar+",
            neutral_species="Cl2",
            ion_energy_eV=100.0,
            neutral_chlorine_atoms_per_particle=2,
            surface_chlorination_coefficient=0.07,
            physical_sputter_yield_si_per_ion=0.07,
            ion_enhanced_yield_si_per_ion_at_full_chlorination=0.83,
            product_chlorine_atoms_per_si=4,
            bulk_si_atom_density_m3=5.0e28,
            maximum_measured_incidence_angle_deg=60.0,
            angular_response_measured_for_neutral=False,
            evidence={
                "surface_chlorination_coefficient": ParameterEvidence(
                    thesis, "beam-yield regression", note=(
                        "Table 3.4 s=0.07 for molecular chlorine; retained as "
                        "an effective Cl-site-fill coefficient, not promoted "
                        "to a universal sticking probability"),
                    supports_prediction_within_declared_domain=True),
                "physical_sputter_yield_si_per_ion": ParameterEvidence(
                    thesis, "controlled-beam measurement/regression",
                    note="Table 3.4 Y0=0.07 at 100 eV",
                    supports_prediction_within_declared_domain=True),
                "ion_enhanced_yield_si_per_ion_at_full_chlorination":
                    ParameterEvidence(
                        thesis, "beam-yield regression",
                        note="Table 3.4 beta=0.83 at 100 eV",
                        supports_prediction_within_declared_domain=True),
                "product_chlorine_atoms_per_si": ParameterEvidence(
                    thesis, "published model assumption",
                    note="Table 3.3 assumes SiCl4 dominates in this low-energy model"),
                "bulk_si_atom_density_m3": ParameterEvidence(
                    "crystalline-Si density and molar mass; rounded 5.0e28 atoms/m3",
                    "derived physical constant",
                    supports_prediction_within_declared_domain=True),
                "angular_response": ParameterEvidence(
                    angular, "cross-neutral model transfer", note=(
                        "Ar+/Cl2 angular yields were not separately reported")),
            },
        )

    @classmethod
    def atomic_chlorine_100eV(cls):
        """Table-3.4 Ar+ / atomic-Cl card at the measured 100 eV condition."""
        base = cls.molecular_chlorine_100eV()
        thesis = (
            "chang-thesis, "
            "Chapter 3, Tables 3.3--3.4 and Eqs. 3.9--3.11"
        )
        evidence = dict(base.evidence)
        evidence.update({
            "surface_chlorination_coefficient": ParameterEvidence(
                thesis, "beam-yield regression", note="Table 3.4 s=0.30",
                supports_prediction_within_declared_domain=True),
            "ion_enhanced_yield_si_per_ion_at_full_chlorination":
                ParameterEvidence(
                    thesis, "beam-yield regression",
                    note="Table 3.4 beta=3.59 at 100 eV",
                    supports_prediction_within_declared_domain=True),
            "angular_response": ParameterEvidence(
                "chang-thesis, Chapter 3, Figures 3.7--3.8",
                "controlled-beam angular measurement",
                note="100 eV Ar+ with atomic Cl, measured through 60 degrees",
                supports_prediction_within_declared_domain=True),
        })
        return cls(
            ion_species="Ar+",
            neutral_species="Cl",
            ion_energy_eV=100.0,
            neutral_chlorine_atoms_per_particle=1,
            surface_chlorination_coefficient=0.30,
            physical_sputter_yield_si_per_ion=0.07,
            ion_enhanced_yield_si_per_ion_at_full_chlorination=3.59,
            product_chlorine_atoms_per_si=4,
            bulk_si_atom_density_m3=base.bulk_si_atom_density_m3,
            maximum_measured_incidence_angle_deg=60.0,
            angular_response_measured_for_neutral=True,
            evidence=evidence,
        )


@dataclass(frozen=True)
class ChangSawinArClSiStepResult:
    state: ChangSawinArClSiState
    etch_velocity_m_s: np.ndarray
    chlorination_fraction: np.ndarray
    physical_removed_si_atoms_m2: np.ndarray
    chemical_removed_si_atoms_m2: np.ndarray
    consumed_chlorine_atoms_m2: np.ndarray
    consumed_neutral_particles_m2: np.ndarray
    steady_site_balance_residual_cl_atoms_m2_s: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity

    def __post_init__(self):
        names = (
            "etch_velocity_m_s", "chlorination_fraction",
            "physical_removed_si_atoms_m2", "chemical_removed_si_atoms_m2",
            "consumed_chlorine_atoms_m2", "consumed_neutral_particles_m2",
            "steady_site_balance_residual_cl_atoms_m2_s",
        )
        arrays = {}
        for name in names:
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Chang--Sawin step result")
            value.setflags(write=False)
            arrays[name] = value
            object.__setattr__(self, name, value)
        if (np.any(arrays["chlorination_fraction"] > 1.0)
                or not isinstance(self.state, ChangSawinArClSiState)
                or not isinstance(self.material_exchange, SurfaceMaterialExchange)
                or not isinstance(self.validity, MechanismValidity)):
            raise ValueError("invalid Chang--Sawin step-result contract")
        populations = validate_surface_product_routing(
            self.material_exchange, tuple(self.product_populations))
        object.__setattr__(self, "product_populations", populations)
        removed = (
            arrays["physical_removed_si_atoms_m2"]
            + arrays["chemical_removed_si_atoms_m2"])
        scale = np.maximum(removed, 1.0)
        if np.any(
            np.abs(
                self.material_exchange.removed_units_m2["Si_atom"] - removed)
            > 64.0 * np.finfo(float).eps * scale
        ):
            raise ValueError("Chang--Sawin Si branches do not close")


class ChangSawinArClSiMechanism:
    """Fixed-condition Ar+/Cl or Ar+/Cl2 poly-Si surface closure."""

    def __init__(self, parameters: ChangSawinArClSiParameters | None = None):
        self.parameters = (
            ChangSawinArClSiParameters.molecular_chlorine_100eV()
            if parameters is None else parameters)
        if not isinstance(self.parameters, ChangSawinArClSiParameters):
            raise TypeError("parameters must be ChangSawinArClSiParameters")

    @staticmethod
    def initial_state(shape=()):
        return ChangSawinArClSiState.bare(shape)

    @property
    def provenance(self):
        par = self.parameters
        return MappingProxyType({
            "model": "Chang--Sawin steady Ar+/chlorine/poly-Si site balance",
            "equations": ["Chang thesis Eq. 3.9", "Eq. 3.10", "Eq. 3.11"],
            "ion_species": par.ion_species,
            "neutral_species": par.neutral_species,
            "fixed_condition": {
                "ion_energy_eV": par.ion_energy_eV,
                "source_substrate_temperature_K": 313.15,
            },
            "parameters": {
                "neutral_chlorine_atoms_per_particle":
                    par.neutral_chlorine_atoms_per_particle,
                "surface_chlorination_coefficient":
                    par.surface_chlorination_coefficient,
                "physical_sputter_yield_si_per_ion":
                    par.physical_sputter_yield_si_per_ion,
                "ion_enhanced_yield_si_per_ion_at_full_chlorination":
                    par.ion_enhanced_yield_si_per_ion_at_full_chlorination,
                "product_chlorine_atoms_per_si":
                    par.product_chlorine_atoms_per_si,
                "bulk_si_atom_density_m3": par.bulk_si_atom_density_m3,
            },
            "evidence": {
                name: {
                    "source": item.source,
                    "evidence_type": item.evidence_type,
                    "supports_prediction_within_declared_domain":
                        item.supports_prediction_within_declared_domain,
                    "note": item.note,
                }
                for name, item in par.evidence.items()
            },
            "claim": (
                "beam-regressed mechanistic closure at one ion energy; "
                "not a first-principles interatomic model"
            ),
        })

    @staticmethod
    def _population_measures(population, angular_factor, shape):
        if isinstance(population, FaceResolvedEnergeticFlux):
            total = np.broadcast_to(population.flux_m2_s, shape)
            weighted = np.bincount(
                population.event_face,
                weights=(
                    population.event_flux_m2_s
                    * angular_factor(population.event_cosine_incidence)
                ),
                minlength=population.face_count,
            )
            return total, np.broadcast_to(weighted, shape)
        if isinstance(population, EnergeticFlux):
            total = np.broadcast_to(
                np.asarray(population.flux_m2_s, dtype=float), shape)
            mean_factor = float(np.dot(
                population.weight,
                angular_factor(population.cosine_incidence),
            ))
            return total, total * mean_factor
        raise TypeError(type(population).__name__)  # pragma: no cover

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name != par.neutral_species
            and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != par.ion_species
            and np.any(np.asarray(item.flux_m2_s) > 0.0)
        }))
        wrong_energy = False
        exceeds_measured_angle = False
        uses_unmeasured_molecular_angle = False
        for population in fluxes.energetic_fluxes:
            if population.name != par.ion_species:
                continue
            if isinstance(population, FaceResolvedEnergeticFlux):
                selected = population.event_flux_m2_s > 0.0
                energy = population.event_energy_eV[selected]
                cosine = population.event_cosine_incidence[selected]
            else:
                selected = population.weight > 0.0
                energy = population.energy_eV[selected]
                cosine = population.cosine_incidence[selected]
            wrong_energy |= bool(np.any(
                np.abs(energy - par.ion_energy_eV) > par.energy_tolerance_eV))
            angle = np.rad2deg(np.arccos(cosine))
            exceeds_measured_angle |= bool(np.any(
                angle > par.maximum_measured_incidence_angle_deg + 1e-10))
            uses_unmeasured_molecular_angle |= bool(
                not par.angular_response_measured_for_neutral
                and np.any(angle > 1e-10))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no declared Ar+/chlorine/Si channel")
        if wrong_energy:
            reasons.append(
                f"ion energy leaves the fixed {par.ion_energy_eV:g} eV beam card")
        if exceeds_measured_angle:
            reasons.append("ion incidence exceeds the measured angular range")
        if uses_unmeasured_molecular_angle:
            reasons.append("molecular-Cl2 angular response was not measured")
        nonpredictive = [
            name for name, item in par.evidence.items()
            if not item.supports_prediction_within_declared_domain
            and name != "angular_response"
        ]
        if uses_unmeasured_molecular_angle:
            nonpredictive.append("angular_response")
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "quasi-steady chlorination has no measured transient site density",
                "SiCl4 is an assumed dominant low-energy product branch",
                "product energy-angle distributions are absent",
                "surface roughness, crystallographic orientation, doping, and damage memory are unresolved",
                "off-normal source algebra multiplies yield after its normal-incidence site balance",
            ),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=tuple(sorted(set(nonpredictive))),
        )

    def advance(
            self, state: ChangSawinArClSiState, fluxes: SurfaceFluxes,
            duration_s: float, *, strict=True):
        if not isinstance(state, ChangSawinArClSiState):
            raise TypeError("Chang--Sawin mechanism requires ChangSawinArClSiState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))
        par = self.parameters
        shape = state.removed_si_atoms_m2.shape
        neutral_flux = np.broadcast_to(
            np.asarray(
                fluxes.neutral_flux_m2_s.get(par.neutral_species, 0.0),
                dtype=float),
            shape,
        )
        ion_flux = np.zeros(shape)
        angular_ion_flux = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name != par.ion_species:
                continue
            total, weighted = self._population_measures(
                population, _chang_sawin_chemical_angular_factor, shape)
            ion_flux = ion_flux + total
            angular_ion_flux = angular_ion_flux + weighted
        adsorption_drive = par.surface_chlorination_coefficient * neutral_flux
        removal_drive = (
            par.product_chlorine_atoms_per_si
            * par.ion_enhanced_yield_si_per_ion_at_full_chlorination
            * ion_flux
        )
        denominator = adsorption_drive + removal_drive
        chlorination = np.zeros(shape)
        active = denominator > 0.0
        chlorination[active] = adsorption_drive[active] / denominator[active]
        physical_rate = (
            angular_ion_flux
            * par.physical_sputter_yield_si_per_ion
            * (1.0 - chlorination)
        )
        chemical_rate = (
            angular_ion_flux
            * par.ion_enhanced_yield_si_per_ion_at_full_chlorination
            * chlorination
        )
        duration = float(duration_s)
        physical_removed = physical_rate * duration
        chemical_removed = chemical_rate * duration
        removed = physical_removed + chemical_removed
        consumed_cl = (
            par.product_chlorine_atoms_per_si * chemical_removed)
        consumed_neutral = (
            consumed_cl / par.neutral_chlorine_atoms_per_particle)
        site_balance_residual = np.maximum(
            adsorption_drive * (1.0 - chlorination)
            - par.product_chlorine_atoms_per_si * chemical_rate,
            0.0,
        )
        updated = ChangSawinArClSiState(
            state.removed_si_atoms_m2 + removed,
            state.consumed_chlorine_atoms_m2 + consumed_cl,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed},
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=(
                "Si versus SiCl4 branching follows the published simplified model",
                "neither product has a measured launch energy-angle distribution",
            ),
        )
        products = (
            SurfaceProductPopulation(
                name="Si_physical",
                source_inventory="Si_atom",
                integrated_particle_count_m2=physical_removed,
                material_units_per_particle=1.0,
                mass_amu=28.085,
                provenance={
                    "source": par.evidence[
                        "physical_sputter_yield_si_per_ion"].source,
                    "branch": "published physical-sputter branch",
                    "missing": "differential emitted energy-angle distribution",
                },
            ),
            SurfaceProductPopulation(
                name="SiCl4_chemical",
                source_inventory="Si_atom",
                integrated_particle_count_m2=chemical_removed,
                material_units_per_particle=1.0,
                mass_amu=169.885,
                provenance={
                    "source": par.evidence[
                        "product_chlorine_atoms_per_si"].source,
                    "branch": "published low-energy product assumption",
                    "missing": "measured branching and differential emission",
                },
            ),
        )
        return ChangSawinArClSiStepResult(
            state=updated,
            etch_velocity_m_s=(
                physical_rate + chemical_rate
            ) / par.bulk_si_atom_density_m3,
            chlorination_fraction=chlorination,
            physical_removed_si_atoms_m2=physical_removed,
            chemical_removed_si_atoms_m2=chemical_removed,
            consumed_chlorine_atoms_m2=consumed_cl,
            consumed_neutral_particles_m2=consumed_neutral,
            steady_site_balance_residual_cl_atoms_m2_s=site_balance_residual,
            material_exchange=exchange,
            product_populations=products,
            validity=validity,
        )
