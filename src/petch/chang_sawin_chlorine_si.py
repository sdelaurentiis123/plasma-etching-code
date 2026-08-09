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


@dataclass(frozen=True)
class ChangSawinClIonSiParameters:
    """Measured Cl+/Cl/poly-Si cards from Chang thesis Table 5.2.

    The tabulated coefficients are interpolated linearly in ``sqrt(E)`` only
    between measured cards.  That is the source's stated energy coordinate;
    no threshold, prefactor, or feature-rate target is fitted here.
    """

    ion_energy_eV: np.ndarray
    surface_chlorination_coefficient: np.ndarray
    ion_enhanced_yield_si_per_ion_at_full_chlorination: np.ndarray
    product_chlorine_atoms_per_si: int
    bulk_si_atom_density_m3: float
    minimum_measured_neutral_to_ion_ratio: float
    maximum_measured_incidence_angle_deg: float
    evidence: Mapping[str, ParameterEvidence]

    def __post_init__(self):
        energy = np.asarray(self.ion_energy_eV, dtype=float).copy()
        chlorination = np.asarray(
            self.surface_chlorination_coefficient, dtype=float).copy()
        enhanced = np.asarray(
            self.ion_enhanced_yield_si_per_ion_at_full_chlorination,
            dtype=float,
        ).copy()
        evidence = dict(self.evidence)
        required = {
            "surface_chlorination_coefficient",
            "ion_enhanced_yield_si_per_ion_at_full_chlorination",
            "product_chlorine_atoms_per_si",
            "bulk_si_atom_density_m3",
            "angular_response",
            "energy_interpolation",
        }
        if (
            energy.ndim != 1
            or energy.size < 2
            or chlorination.shape != energy.shape
            or enhanced.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(np.diff(energy) <= 0.0)
            or np.any(energy <= 0.0)
            or np.any(~np.isfinite(chlorination))
            or np.any((chlorination < 0.0) | (chlorination > 1.0))
            or np.any(~np.isfinite(enhanced))
            or np.any(enhanced < 0.0)
            or int(self.product_chlorine_atoms_per_si)
            != self.product_chlorine_atoms_per_si
            or self.product_chlorine_atoms_per_si <= 0
            or not np.isfinite(self.bulk_si_atom_density_m3)
            or self.bulk_si_atom_density_m3 <= 0.0
            or not np.isfinite(self.minimum_measured_neutral_to_ion_ratio)
            or self.minimum_measured_neutral_to_ion_ratio < 0.0
            or not np.isfinite(self.maximum_measured_incidence_angle_deg)
            or not 0.0 < self.maximum_measured_incidence_angle_deg < 90.0
            or set(evidence) != required
            or any(not isinstance(item, ParameterEvidence)
                   for item in evidence.values())
        ):
            raise ValueError("invalid Chang--Sawin Cl+/Cl/Si parameters")
        for value in (energy, chlorination, enhanced):
            value.setflags(write=False)
        object.__setattr__(self, "ion_energy_eV", energy)
        object.__setattr__(
            self, "surface_chlorination_coefficient", chlorination)
        object.__setattr__(
            self,
            "ion_enhanced_yield_si_per_ion_at_full_chlorination",
            enhanced,
        )
        object.__setattr__(
            self, "product_chlorine_atoms_per_si",
            int(self.product_chlorine_atoms_per_si))
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @classmethod
    def chang_thesis_table_5_2(cls):
        thesis = (
            "chang-thesis, Chapter 5, Table 5.2 and Eqs. 5.1--5.5"
        )
        return cls(
            ion_energy_eV=np.asarray((35.0, 60.0, 100.0)),
            surface_chlorination_coefficient=np.asarray((0.18, 0.32, 0.45)),
            ion_enhanced_yield_si_per_ion_at_full_chlorination=np.asarray(
                (1.14, 2.42, 3.61)),
            product_chlorine_atoms_per_si=4,
            bulk_si_atom_density_m3=5.0e28,
            minimum_measured_neutral_to_ion_ratio=5.0,
            maximum_measured_incidence_angle_deg=70.0,
            evidence={
                "surface_chlorination_coefficient": ParameterEvidence(
                    thesis,
                    "controlled-beam regression",
                    note="Table 5.2 s at 35, 60, and 100 eV",
                    supports_prediction_within_declared_domain=True,
                ),
                "ion_enhanced_yield_si_per_ion_at_full_chlorination":
                    ParameterEvidence(
                        thesis,
                        "controlled-beam regression",
                        note="Table 5.2 beta at 35, 60, and 100 eV",
                        supports_prediction_within_declared_domain=True,
                    ),
                "product_chlorine_atoms_per_si": ParameterEvidence(
                    "chang-thesis, Chapter 5, Table 5.1",
                    "published low-energy product assumption",
                    note=(
                        "SiCl4 is assumed with simultaneous Cl+/Cl; the "
                        "source warns that unsaturated SiClx matters at "
                        "neutral-to-ion ratios below five"
                    ),
                ),
                "bulk_si_atom_density_m3": ParameterEvidence(
                    "crystalline-Si density and molar mass; rounded 5.0e28 atoms/m3",
                    "derived physical constant",
                    supports_prediction_within_declared_domain=True,
                ),
                "angular_response": ParameterEvidence(
                    "chang-thesis, Chapter 5, Figures 5.8--5.9",
                    "controlled-beam angular measurement",
                    note=(
                        "35 eV Cl+/Cl: flat through about 40 deg, 30% and "
                        "50% reductions at 60 and 70 deg"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "energy_interpolation": ParameterEvidence(
                    "chang-thesis, Chapter 5, Figure 5.13 and Table 5.2",
                    "source-stated interpolation coordinate",
                    note=(
                        "piecewise linear in sqrt(E) through the three "
                        "printed cards; no extrapolation is permitted"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
            },
        )

    @property
    def energy_domain_eV(self):
        return float(self.ion_energy_eV[0]), float(self.ion_energy_eV[-1])

    def coefficients(self, energy_eV, *, allow_extrapolation=False):
        energy = np.asarray(energy_eV, dtype=float)
        lower, upper = self.energy_domain_eV
        if (
            np.any(~np.isfinite(energy))
            or np.any(energy <= 0.0)
            or (
                not allow_extrapolation
                and np.any((energy < lower) | (energy > upper))
            )
        ):
            raise ValueError(
                "Cl+/Cl beam energy leaves the measured 35--100 eV domain")
        coordinate = np.sqrt(energy)
        source_coordinate = np.sqrt(self.ion_energy_eV)

        def interpolate(values):
            result = np.interp(coordinate, source_coordinate, values)
            if allow_extrapolation:
                result = np.asarray(result)
                below = coordinate < source_coordinate[0]
                above = coordinate > source_coordinate[-1]
                low_slope = (
                    (values[1] - values[0])
                    / (source_coordinate[1] - source_coordinate[0]))
                high_slope = (
                    (values[-1] - values[-2])
                    / (source_coordinate[-1] - source_coordinate[-2]))
                result = np.where(
                    below,
                    values[0] + low_slope * (
                        coordinate - source_coordinate[0]),
                    result,
                )
                result = np.where(
                    above,
                    values[-1] + high_slope * (
                        coordinate - source_coordinate[-1]),
                    result,
                )
                result = np.maximum(result, 0.0)
            return result

        return (
            interpolate(self.surface_chlorination_coefficient),
            interpolate(
                self.ion_enhanced_yield_si_per_ion_at_full_chlorination),
        )


@dataclass(frozen=True)
class ChangSawinClIonSiStepResult:
    state: ChangSawinArClSiState
    etch_velocity_m_s: np.ndarray
    chlorination_fraction: np.ndarray
    removed_si_atoms_m2: np.ndarray
    chlorine_atoms_supplied_by_neutrals_m2: np.ndarray
    chlorine_atoms_supplied_by_ions_m2: np.ndarray
    steady_site_balance_residual_cl_atoms_m2_s: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity

    def __post_init__(self):
        names = (
            "etch_velocity_m_s",
            "chlorination_fraction",
            "removed_si_atoms_m2",
            "chlorine_atoms_supplied_by_neutrals_m2",
            "chlorine_atoms_supplied_by_ions_m2",
            "steady_site_balance_residual_cl_atoms_m2_s",
        )
        arrays = {}
        for name in names:
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Chang--Sawin Cl+ step result")
            value.setflags(write=False)
            arrays[name] = value
            object.__setattr__(self, name, value)
        if (
            np.any(arrays["chlorination_fraction"] > 1.0)
            or not isinstance(self.state, ChangSawinArClSiState)
            or not isinstance(self.material_exchange, SurfaceMaterialExchange)
            or not isinstance(self.validity, MechanismValidity)
        ):
            raise ValueError("invalid Chang--Sawin Cl+ result contract")
        populations = validate_surface_product_routing(
            self.material_exchange, tuple(self.product_populations))
        object.__setattr__(self, "product_populations", populations)
        supplied = (
            arrays["chlorine_atoms_supplied_by_neutrals_m2"]
            + arrays["chlorine_atoms_supplied_by_ions_m2"]
        )
        required = 4.0 * arrays["removed_si_atoms_m2"]
        if np.any(
            np.abs(supplied - required)
            > 128.0 * np.finfo(float).eps * np.maximum(required, 1.0)
        ):
            raise ValueError("Chang--Sawin Cl+ chlorine inventory does not close")


class ChangSawinClIonSiMechanism:
    """Absolute Cl+/Cl/poly-Si beam closure across the measured cards.

    For a monoenergetic normal beam this evaluates Chang Eqs. 5.1--5.3
    directly.  For a resolved narrow distribution, the source terms are
    quadrature-integrated before the one steady site balance is solved.  The
    latter is an explicit deterministic extension and remains listed as a
    model-form limitation because the source did not measure broad IEADs.
    """

    def __init__(self, parameters: ChangSawinClIonSiParameters | None = None):
        self.parameters = (
            ChangSawinClIonSiParameters.chang_thesis_table_5_2()
            if parameters is None else parameters)
        if not isinstance(self.parameters, ChangSawinClIonSiParameters):
            raise TypeError("parameters must be ChangSawinClIonSiParameters")

    @staticmethod
    def initial_state(shape=()):
        return ChangSawinArClSiState.bare(shape)

    @property
    def provenance(self):
        par = self.parameters
        return MappingProxyType({
            "model": "Chang--Sawin steady Cl+/Cl/poly-Si site balance",
            "equations": ["Chang thesis Eq. 5.1", "Eq. 5.2", "Eq. 5.3"],
            "ion_species": "Cl+",
            "neutral_species": "Cl",
            "measured_energy_cards_eV": par.ion_energy_eV.tolist(),
            "surface_chlorination_coefficients": (
                par.surface_chlorination_coefficient.tolist()),
            "full_chlorination_yields": (
                par.ion_enhanced_yield_si_per_ion_at_full_chlorination.tolist()),
            "energy_interpolation": "piecewise linear in sqrt(E); no extrapolation",
            "minimum_measured_neutral_to_ion_ratio": (
                par.minimum_measured_neutral_to_ion_ratio),
            "maximum_measured_incidence_angle_deg": (
                par.maximum_measured_incidence_angle_deg),
            "evidence": {
                name: {
                    "source": item.source,
                    "evidence_type": item.evidence_type,
                    "supports_prediction_within_declared_domain": (
                        item.supports_prediction_within_declared_domain),
                    "note": item.note,
                }
                for name, item in par.evidence.items()
            },
            "claim": (
                "beam-regressed absolute-rate closure inside the measured "
                "energy, angle, and flux-ratio domain"
            ),
        })

    @staticmethod
    def _events(population):
        if isinstance(population, FaceResolvedEnergeticFlux):
            return (
                population.event_face,
                population.event_flux_m2_s,
                population.event_energy_eV,
                population.event_cosine_incidence,
                population.face_count,
            )
        if isinstance(population, EnergeticFlux):
            flux = float(np.asarray(population.flux_m2_s))
            return (
                np.zeros(population.weight.size, dtype=int),
                flux * population.weight,
                population.energy_eV,
                population.cosine_incidence,
                1,
            )
        raise TypeError(type(population).__name__)  # pragma: no cover

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name != "Cl" and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != "Cl+" and np.any(np.asarray(item.flux_m2_s) > 0.0)
        }))
        energy_outside = False
        angle_outside = False
        total_ion = None
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl+":
                continue
            if isinstance(population, FaceResolvedEnergeticFlux):
                selected = population.event_flux_m2_s > 0.0
                energy = population.event_energy_eV[selected]
                cosine = population.event_cosine_incidence[selected]
                ion = population.flux_m2_s
            else:
                selected = population.weight > 0.0
                energy = population.energy_eV[selected]
                cosine = population.cosine_incidence[selected]
                ion = np.asarray(population.flux_m2_s)
            lower, upper = par.energy_domain_eV
            energy_outside |= bool(np.any((energy < lower) | (energy > upper)))
            angle_outside |= bool(np.any(
                np.rad2deg(np.arccos(cosine))
                > par.maximum_measured_incidence_angle_deg + 1.0e-10))
            total_ion = ion if total_ion is None else total_ion + ion
        ratio_outside = False
        if total_ion is not None:
            ion = np.asarray(total_ion, dtype=float)
            neutral = np.asarray(
                fluxes.neutral_flux_m2_s.get("Cl", 0.0), dtype=float)
            ion, neutral = np.broadcast_arrays(ion, neutral)
            active = ion > 0.0
            neutral_to_ion = np.full(ion.shape, np.inf, dtype=float)
            np.divide(neutral, ion, out=neutral_to_ion, where=active)
            ratio_outside = bool(np.any(
                active
                & (neutral_to_ion < par.minimum_measured_neutral_to_ion_ratio)))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no declared Cl+/Cl/Si channel")
        if energy_outside:
            reasons.append("Cl+ energy leaves the measured 35--100 eV cards")
        if angle_outside:
            reasons.append("Cl+ incidence exceeds the measured 70 degree range")
        if ratio_outside:
            reasons.append("Cl/Cl+ flux ratio falls below the measured model domain of five")
        nonpredictive = tuple(sorted(
            name for name, item in par.evidence.items()
            if not item.supports_prediction_within_declared_domain))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "Cl2+ surface kinetics were outside Chang's measured mechanism",
                "broad-IEAD coefficients are flux-quadrature extensions of monoenergetic cards",
                "SiCl4 is the low-energy simultaneous-Cl+/Cl product assumption",
                "SiClx branching below neutral-to-ion ratio five is unresolved",
                "etch-product return, implantation, damage, and transient site density are unresolved",
            ),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive,
        )

    def _integrated_drives(self, fluxes, shape, *, allow_extrapolation=False):
        ion_flux = np.zeros(shape)
        ionic_chlorination_drive = np.zeros(shape)
        chlorination_coefficient_flux = np.zeros(shape)
        enhanced_removal_capacity = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl+":
                continue
            face, event_flux, energy, cosine, face_count = self._events(population)
            s_coefficient, beta = self.parameters.coefficients(
                energy, allow_extrapolation=allow_extrapolation)
            angular = _chang_sawin_chemical_angular_factor(cosine)
            if isinstance(population, FaceResolvedEnergeticFlux):
                integrated = lambda values: np.bincount(  # noqa: E731
                    face, weights=values, minlength=face_count)
                ion_flux = ion_flux + np.broadcast_to(
                    integrated(event_flux), shape)
                ionic_chlorination_drive = (
                    ionic_chlorination_drive + np.broadcast_to(
                        integrated(event_flux * angular), shape))
                chlorination_coefficient_flux = (
                    chlorination_coefficient_flux + np.broadcast_to(
                        integrated(event_flux * s_coefficient), shape))
                enhanced_removal_capacity = (
                    enhanced_removal_capacity + np.broadcast_to(
                        integrated(event_flux * angular * beta), shape))
            else:
                ion_flux = ion_flux + np.broadcast_to(
                    np.sum(event_flux), shape)
                ionic_chlorination_drive = (
                    ionic_chlorination_drive + np.broadcast_to(
                        np.sum(event_flux * angular), shape))
                chlorination_coefficient_flux = (
                    chlorination_coefficient_flux + np.broadcast_to(
                        np.sum(event_flux * s_coefficient), shape))
                enhanced_removal_capacity = (
                    enhanced_removal_capacity + np.broadcast_to(
                        np.sum(event_flux * angular * beta), shape))
        effective_s = np.zeros(shape)
        active = ion_flux > 0.0
        effective_s[active] = (
            chlorination_coefficient_flux[active] / ion_flux[active])
        return (
            ion_flux,
            ionic_chlorination_drive,
            effective_s,
            enhanced_removal_capacity,
        )

    def advance(
        self,
        state: ChangSawinArClSiState,
        fluxes: SurfaceFluxes,
        duration_s: float,
        *,
        strict=True,
    ):
        if not isinstance(state, ChangSawinArClSiState):
            raise TypeError("Chang--Sawin Cl+ mechanism requires ChangSawinArClSiState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))
        shape = state.removed_si_atoms_m2.shape
        _, ionic_drive, effective_s, removal_capacity = (
            self._integrated_drives(
                fluxes, shape, allow_extrapolation=not strict))
        neutral_flux = np.broadcast_to(
            np.asarray(fluxes.neutral_flux_m2_s.get("Cl", 0.0), dtype=float),
            shape,
        )
        neutral_drive = effective_s * neutral_flux
        supply_drive = neutral_drive + ionic_drive
        denominator = (
            supply_drive
            + self.parameters.product_chlorine_atoms_per_si * removal_capacity)
        chlorination = np.zeros(shape)
        active = denominator > 0.0
        chlorination[active] = supply_drive[active] / denominator[active]
        removal_rate = removal_capacity * chlorination
        neutral_supply_rate = neutral_drive * (1.0 - chlorination)
        ionic_supply_rate = ionic_drive * (1.0 - chlorination)
        balance_residual = np.abs(
            neutral_supply_rate + ionic_supply_rate
            - self.parameters.product_chlorine_atoms_per_si * removal_rate)
        duration = float(duration_s)
        removed = removal_rate * duration
        neutral_supplied = neutral_supply_rate * duration
        ionic_supplied = ionic_supply_rate * duration
        updated = ChangSawinArClSiState(
            state.removed_si_atoms_m2 + removed,
            state.consumed_chlorine_atoms_m2
            + neutral_supplied + ionic_supplied,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed},
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=(
                "SiCl4 is the source's simultaneous-Cl+/Cl low-energy product assumption",
                "product energy-angle distribution is not measured",
            ),
        )
        products = (
            SurfaceProductPopulation(
                name="SiCl4_chemical",
                source_inventory="Si_atom",
                integrated_particle_count_m2=removed,
                material_units_per_particle=1.0,
                mass_amu=169.885,
                provenance={
                    "source": self.parameters.evidence[
                        "product_chlorine_atoms_per_si"].source,
                    "branch": "published low-energy simultaneous-Cl+/Cl assumption",
                    "missing": "measured differential emission distribution",
                },
            ),
        )
        return ChangSawinClIonSiStepResult(
            state=updated,
            etch_velocity_m_s=(
                removal_rate / self.parameters.bulk_si_atom_density_m3),
            chlorination_fraction=chlorination,
            removed_si_atoms_m2=removed,
            chlorine_atoms_supplied_by_neutrals_m2=neutral_supplied,
            chlorine_atoms_supplied_by_ions_m2=ionic_supplied,
            steady_site_balance_residual_cl_atoms_m2_s=balance_residual,
            material_exchange=exchange,
            product_populations=products,
            validity=validity,
        )


@dataclass(frozen=True)
class ChangSawinSiCl2SuppressionParameters:
    """Chang Chapter-5 etch-product suppression constants.

    ``delta`` is the independently measured stable-film SiCl2 sticking
    coefficient.  ``eta`` was regressed by Chang to the separate SiCl2/Cl+/Cl
    beam sweep in Figure 5.10 and used in the printed Eq. 5.6/Figure 5.14.
    The latter is therefore constitutive evidence for transfer to another
    reactor, not a first-principles reaction probability.
    """

    sicl2_sticking_coefficient: float
    chlorinated_sicl2_reaction_coefficient: float
    evidence: Mapping[str, ParameterEvidence]

    def __post_init__(self):
        evidence = dict(self.evidence)
        if (
            not np.isfinite(self.sicl2_sticking_coefficient)
            or not 0.0 <= self.sicl2_sticking_coefficient <= 1.0
            or not np.isfinite(self.chlorinated_sicl2_reaction_coefficient)
            or self.chlorinated_sicl2_reaction_coefficient < 0.0
            or set(evidence) != {
                "sicl2_sticking_coefficient",
                "chlorinated_sicl2_reaction_coefficient",
                "single_coverage_closure",
            }
            or any(not isinstance(item, ParameterEvidence)
                   for item in evidence.values())
        ):
            raise ValueError("invalid Chang SiCl2-suppression parameters")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @classmethod
    def chang_thesis_equation_5_6(cls):
        return cls(
            sicl2_sticking_coefficient=0.3,
            chlorinated_sicl2_reaction_coefficient=10.0,
            evidence={
                "sicl2_sticking_coefficient": ParameterEvidence(
                    "chang-thesis, Chapter 5 section 5.4",
                    "controlled-beam deposition measurement",
                    note=(
                        "delta=0.3 from DCS-beam flux and laser-"
                        "interferometry stable-film growth"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "chlorinated_sicl2_reaction_coefficient": ParameterEvidence(
                    "chang-thesis, Chapter 5 Eq. 5.6 and Figure 5.14",
                    "controlled-beam regression",
                    note=(
                        "eta=10 fits the independent 35 eV, Cl/Cl+=200 "
                        "SiCl2 suppression sweep"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "single_coverage_closure": ParameterEvidence(
                    "chang-thesis, Chapter 5 Eq. 5.6",
                    "published semi-empirical model form",
                    note=(
                        "one steady coverage predicts etch suppression but "
                        "does not resolve a conservative SiClx film inventory"
                    ),
                ),
            },
        )


@dataclass(frozen=True)
class ChangSawinSiCl2SuppressionStepResult:
    state: ChangSawinArClSiState
    etch_velocity_m_s: np.ndarray
    chlorination_fraction: np.ndarray
    removed_si_atoms_m2: np.ndarray
    sicl2_to_clplus_flux_ratio: np.ndarray
    site_balance_residual_sites_m2_s: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity

    def __post_init__(self):
        arrays = {}
        for name in (
            "etch_velocity_m_s",
            "chlorination_fraction",
            "removed_si_atoms_m2",
            "sicl2_to_clplus_flux_ratio",
            "site_balance_residual_sites_m2_s",
        ):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Chang SiCl2-suppression result")
            value.setflags(write=False)
            arrays[name] = value
            object.__setattr__(self, name, value)
        if (
            np.any(arrays["chlorination_fraction"] > 1.0)
            or not isinstance(self.state, ChangSawinArClSiState)
            or not isinstance(self.material_exchange, SurfaceMaterialExchange)
            or not isinstance(self.validity, MechanismValidity)
        ):
            raise ValueError("invalid Chang SiCl2-suppression contract")
        populations = validate_surface_product_routing(
            self.material_exchange, tuple(self.product_populations))
        object.__setattr__(self, "product_populations", populations)


class ChangSawinClIonSiCl2SuppressionMechanism(ChangSawinClIonSiMechanism):
    """Printed Chang Eq. 5.6 for Cl+/Cl/Si with incident SiCl2.

    The ion-energy and angular drives are integrated exactly as in the base
    mechanism.  Incident SiCl2 then competes in the *one* common steady site
    balance before the substrate-removal rate is evaluated.  This is the
    strongest rate law the source supports.  It deliberately does not invent
    a transport-ready deposited-film population because Eq. 5.6 has no
    independent SiClx coverage or conservative byproduct-state equation.
    """

    def __init__(
        self,
        parameters: ChangSawinClIonSiParameters | None = None,
        suppression_parameters: ChangSawinSiCl2SuppressionParameters | None = None,
    ):
        super().__init__(parameters)
        self.suppression_parameters = (
            ChangSawinSiCl2SuppressionParameters.chang_thesis_equation_5_6()
            if suppression_parameters is None else suppression_parameters
        )
        if not isinstance(
            self.suppression_parameters,
            ChangSawinSiCl2SuppressionParameters,
        ):
            raise TypeError(
                "suppression_parameters must be "
                "ChangSawinSiCl2SuppressionParameters"
            )

    @property
    def provenance(self):
        base = dict(super().provenance)
        suppression = self.suppression_parameters
        base.update({
            "model": "Chang--Sawin Cl+/Cl/Si with SiCl2 suppression",
            "equations": [
                "Chang thesis Eq. 5.2",
                "Chang thesis Eq. 5.6",
            ],
            "depositing_species": "SiCl2",
            "sicl2_sticking_coefficient": (
                suppression.sicl2_sticking_coefficient),
            "chlorinated_sicl2_reaction_coefficient": (
                suppression.chlorinated_sicl2_reaction_coefficient),
            "claim": (
                "beam-regressed absolute etch suppression inside the "
                "measured 35 eV, normal-incidence SiCl2 sweep"
            ),
        })
        return MappingProxyType(base)

    def validity(self, fluxes: SurfaceFluxes):
        base_fluxes = SurfaceFluxes(
            {
                name: value
                for name, value in fluxes.neutral_flux_m2_s.items()
                if name != "SiCl2"
            },
            fluxes.energetic_fluxes,
        )
        base = super().validity(base_fluxes)
        reasons = list(base.reasons)
        # The suppression experiment is only at 35 eV and normal incidence.
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl+":
                continue
            _, event_flux, energy, cosine, _ = self._events(population)
            selected = event_flux > 0.0
            if np.any(np.abs(energy[selected] - 35.0) > 1.0e-10):
                reasons.append(
                    "SiCl2 suppression was measured only at 35 eV Cl+"
                )
            if np.any(np.abs(cosine[selected] - 1.0) > 1.0e-10):
                reasons.append(
                    "SiCl2 suppression was measured only at normal incidence"
                )
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in {"Cl", "SiCl2"}
            and np.any(np.asarray(value) > 0.0)
        ))
        if unsupported_neutral:
            reasons.append(
                "positive incident flux has no declared Cl+/Cl/SiCl2/Si channel"
            )
        nonpredictive = tuple(sorted(set(base.nonpredictive_parameters) | {
            name for name, item in self.suppression_parameters.evidence.items()
            if not item.supports_prediction_within_declared_domain
        }))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(dict.fromkeys(reasons)),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=tuple(dict.fromkeys(
                base.known_model_form_omissions + (
                    "Eq. 5.6 has no conservative SiClx-film state equation",
                    "SiCl2 suppression energy and angular response beyond the 35 eV normal beam is unresolved",
                )
            )),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive,
        )

    def advance(
        self,
        state: ChangSawinArClSiState,
        fluxes: SurfaceFluxes,
        duration_s: float,
        *,
        strict=True,
    ):
        if not isinstance(state, ChangSawinArClSiState):
            raise TypeError(
                "Chang SiCl2 suppression requires ChangSawinArClSiState"
            )
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons)
            )
        shape = state.removed_si_atoms_m2.shape
        ion_flux, ionic_drive, effective_s, removal_capacity = (
            self._integrated_drives(
                fluxes, shape, allow_extrapolation=not strict
            )
        )
        cl_flux = np.broadcast_to(
            np.asarray(fluxes.neutral_flux_m2_s.get("Cl", 0.0), dtype=float),
            shape,
        )
        sicl2_flux = np.broadcast_to(
            np.asarray(
                fluxes.neutral_flux_m2_s.get("SiCl2", 0.0), dtype=float
            ),
            shape,
        )
        delta = self.suppression_parameters.sicl2_sticking_coefficient
        eta = (
            self.suppression_parameters
            .chlorinated_sicl2_reaction_coefficient
        )
        chlorine_drive = effective_s * cl_flux + ionic_drive
        sicl2_drive = delta * sicl2_flux
        denominator = (
            chlorine_drive
            + sicl2_drive
            + self.parameters.product_chlorine_atoms_per_si
            * removal_capacity
            + 3.0 * eta * sicl2_drive
        )
        chlorination = np.zeros(shape)
        active = denominator > 0.0
        chlorination[active] = (
            (chlorine_drive[active] + sicl2_drive[active])
            / denominator[active]
        )
        removal_rate = removal_capacity * chlorination
        # This is the exact algebraic residual of printed Eq. 5.6 in
        # dimensional site-rate form.  It is not promoted to an atom balance.
        site_supply = (
            (chlorine_drive + sicl2_drive) * (1.0 - chlorination)
        )
        site_loss = chlorination * (
            self.parameters.product_chlorine_atoms_per_si
            * removal_capacity
            + 3.0 * eta * sicl2_drive
        )
        residual = np.abs(site_supply - site_loss)
        duration = float(duration_s)
        removed = removal_rate * duration
        # Preserve the base state's explicitly declared substrate inventory.
        # Returned SiCl2 is not counted as substrate material.
        consumed_substrate_cl = (
            self.parameters.product_chlorine_atoms_per_si * removed
        )
        updated = ChangSawinArClSiState(
            state.removed_si_atoms_m2 + removed,
            state.consumed_chlorine_atoms_m2 + consumed_substrate_cl,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed},
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=(
                "Eq. 5.6 predicts substrate etch suppression but has no conservative SiClx-film inventory",
                "returned SiCl2 fate and SiCl4 differential emission are unresolved",
            ),
        )
        products = (
            SurfaceProductPopulation(
                name="SiCl4_substrate_with_SiCl2_suppression",
                source_inventory="Si_atom",
                integrated_particle_count_m2=removed,
                material_units_per_particle=1.0,
                mass_amu=169.885,
                provenance={
                    "source": "chang-thesis, Chapter 5 Eqs. 5.2 and 5.6",
                    "branch": "substrate-Si etch under empirical SiCl2 suppression",
                    "missing": "returned-SiCl2 inventory and differential emission",
                },
            ),
        )
        ratio = np.zeros(shape)
        active_ion = ion_flux > 0.0
        ratio[active_ion] = sicl2_flux[active_ion] / ion_flux[active_ion]
        return ChangSawinSiCl2SuppressionStepResult(
            state=updated,
            etch_velocity_m_s=(
                removal_rate / self.parameters.bulk_si_atom_density_m3
            ),
            chlorination_fraction=chlorination,
            removed_si_atoms_m2=removed,
            sicl2_to_clplus_flux_ratio=ratio,
            site_balance_residual_sites_m2_s=residual,
            material_exchange=exchange,
            product_populations=products,
            validity=validity,
        )


@dataclass(frozen=True)
class BaloochCl2IonSiParameters:
    """Cl2+ removal of highly chlorinated poly-Si from Chang Figure 5.7.

    Chang prints the fitted slope and plots the intercept.  The intercept is
    recovered by the checksum-bound 600-dpi pixel audit in
    ``data/experimental/chang_1998_figure5_7``.  The underlying Cl2+ points
    are attributed to Balooch et al.; Chang states that their approximately
    1e-4 Torr chlorine background produced a highly chlorinated surface.
    """

    slope_si_per_ion_per_sqrt_eV: float
    threshold_energy_eV: float
    minimum_energy_eV: float
    maximum_energy_eV: float
    minimum_high_chlorination_fraction: float
    bulk_si_atom_density_m3: float
    evidence: Mapping[str, ParameterEvidence]

    def __post_init__(self):
        evidence = dict(self.evidence)
        required = {
            "slope_si_per_ion_per_sqrt_eV",
            "threshold_energy_eV",
            "energy_domain_eV",
            "surface_state",
            "bulk_si_atom_density_m3",
        }
        if (
            not np.isfinite(self.slope_si_per_ion_per_sqrt_eV)
            or self.slope_si_per_ion_per_sqrt_eV <= 0.0
            or not np.isfinite(self.threshold_energy_eV)
            or self.threshold_energy_eV <= 0.0
            or not np.isfinite(self.minimum_energy_eV)
            or not np.isfinite(self.maximum_energy_eV)
            or self.minimum_energy_eV < self.threshold_energy_eV - 1.0e-8
            or self.maximum_energy_eV <= self.minimum_energy_eV
            or not np.isfinite(self.minimum_high_chlorination_fraction)
            or not 0.0 < self.minimum_high_chlorination_fraction <= 1.0
            or not np.isfinite(self.bulk_si_atom_density_m3)
            or self.bulk_si_atom_density_m3 <= 0.0
            or set(evidence) != required
            or any(not isinstance(item, ParameterEvidence)
                   for item in evidence.values())
        ):
            raise ValueError("invalid Balooch Cl2+/poly-Si parameters")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @classmethod
    def chang_figure5_7(cls):
        source = "chang-thesis, Figure 5.7, Balooch Cl2+ series"
        return cls(
            slope_si_per_ion_per_sqrt_eV=0.22,
            threshold_energy_eV=25.998846756576185,
            minimum_energy_eV=25.998846756576185,
            maximum_energy_eV=625.0,
            minimum_high_chlorination_fraction=0.85,
            bulk_si_atom_density_m3=5.0e28,
            evidence={
                "slope_si_per_ion_per_sqrt_eV": ParameterEvidence(
                    source,
                    "source-printed controlled-beam fit",
                    note=(
                        "Chang prints slopes 0.06, 0.22, and 0.57 for Cl+, "
                        "Cl2+, and Cl+/Cl, respectively"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "threshold_energy_eV": ParameterEvidence(
                    (
                        "data/experimental/chang_1998_figure5_7/"
                        "digitization_manifest.json"
                    ),
                    "checksum-bound 600-dpi fit-line digitization",
                    note=(
                        "PIL/NumPy replay recovers 25.999 eV with a declared "
                        "2 eV digitization bound"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "energy_domain_eV": ParameterEvidence(
                    source,
                    "plotted controlled-beam support",
                    note="prediction is bounded to the plotted 26--625 eV support",
                    supports_prediction_within_declared_domain=True,
                ),
                "surface_state": ParameterEvidence(
                    "chang-thesis, Chapter 5 text immediately above Figure 5.7",
                    "source interpretation",
                    note=(
                        "Balooch used about 1e-4 Torr chlorine; Chang attributes "
                        "the enhanced yield to saturated chlorine coverage. The "
                        "0.85 guard labels high chlorination and never rescales yield"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "bulk_si_atom_density_m3": ParameterEvidence(
                    "crystalline-Si density and molar mass; rounded 5.0e28 atoms/m3",
                    "derived physical constant",
                    supports_prediction_within_declared_domain=True,
                ),
            },
        )

    @property
    def energy_domain_eV(self):
        return self.minimum_energy_eV, self.maximum_energy_eV

    def yield_si_per_ion(self, energy_eV, *, allow_extrapolation=False):
        energy = np.asarray(energy_eV, dtype=float)
        if (
            np.any(~np.isfinite(energy))
            or np.any(energy < 0.0)
            or (
                not allow_extrapolation
                and np.any(
                    (energy < self.minimum_energy_eV)
                    | (energy > self.maximum_energy_eV)
                )
            )
        ):
            raise ValueError(
                "Cl2+ energy leaves the measured 26--625 eV domain")
        return self.slope_si_per_ion_per_sqrt_eV * np.maximum(
            np.sqrt(energy) - np.sqrt(self.threshold_energy_eV), 0.0)


@dataclass(frozen=True)
class BaloochCl2IonSiStepResult:
    state: ChangSawinArClSiState
    etch_velocity_m_s: np.ndarray
    removed_si_atoms_m2: np.ndarray
    incident_cl2plus_ions_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity

    def __post_init__(self):
        arrays = {}
        for name in (
            "etch_velocity_m_s", "removed_si_atoms_m2",
            "incident_cl2plus_ions_m2",
        ):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Balooch Cl2+ step result")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
            arrays[name] = value
        if (
            not isinstance(self.state, ChangSawinArClSiState)
            or not isinstance(self.material_exchange, SurfaceMaterialExchange)
            or not isinstance(self.validity, MechanismValidity)
        ):
            raise ValueError("invalid Balooch Cl2+ result contract")
        populations = validate_surface_product_routing(
            self.material_exchange, tuple(self.product_populations))
        object.__setattr__(self, "product_populations", populations)


class BaloochCl2IonSiMechanism:
    """Normal-incidence Cl2+ yield on highly chlorinated poly-Si."""

    def __init__(self, parameters: BaloochCl2IonSiParameters | None = None):
        self.parameters = (
            BaloochCl2IonSiParameters.chang_figure5_7()
            if parameters is None else parameters)
        if not isinstance(self.parameters, BaloochCl2IonSiParameters):
            raise TypeError("parameters must be BaloochCl2IonSiParameters")

    @staticmethod
    def initial_state(shape=()):
        return ChangSawinArClSiState.bare(shape)

    @property
    def provenance(self):
        par = self.parameters
        return MappingProxyType({
            "model": "Chang Figure 5.7 / Balooch Cl2+ poly-Si yield",
            "ion_species": "Cl2+",
            "surface_state": "highly chlorinated poly-Si",
            "yield_law": (
                "Y=0.22*max(sqrt(E_eV)-sqrt(25.998846756576185),0)"
            ),
            "measured_energy_domain_eV": list(par.energy_domain_eV),
            "incidence_angle_domain_deg": [0.0, 0.0],
            "feature_depth_used": False,
            "evidence": {
                name: {
                    "source": item.source,
                    "evidence_type": item.evidence_type,
                    "supports_prediction_within_declared_domain": (
                        item.supports_prediction_within_declared_domain),
                    "note": item.note,
                }
                for name, item in par.evidence.items()
            },
        })

    @staticmethod
    def _events(population):
        if isinstance(population, FaceResolvedEnergeticFlux):
            return (
                population.event_face,
                population.event_flux_m2_s,
                population.event_energy_eV,
                population.event_cosine_incidence,
                population.face_count,
            )
        if isinstance(population, EnergeticFlux):
            flux = float(np.asarray(population.flux_m2_s))
            return (
                np.zeros(population.weight.size, dtype=int),
                flux * population.weight,
                population.energy_eV,
                population.cosine_incidence,
                1,
            )
        raise TypeError(type(population).__name__)  # pragma: no cover

    def validity(self, fluxes: SurfaceFluxes, surface_chlorination_fraction):
        par = self.parameters
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name != "Cl" and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != "Cl2+" and np.any(np.asarray(item.flux_m2_s) > 0.0)
        }))
        outside_energy = False
        outside_angle = False
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl2+":
                continue
            _, event_flux, energy, cosine, _ = self._events(population)
            selected = event_flux > 0.0
            lower, upper = par.energy_domain_eV
            outside_energy |= bool(np.any(
                (energy[selected] < lower) | (energy[selected] > upper)))
            outside_angle |= bool(np.any(cosine[selected] < 1.0 - 1.0e-10))
        chlorination = np.asarray(surface_chlorination_fraction, dtype=float)
        outside_surface = bool(
            np.any(~np.isfinite(chlorination))
            or np.any(chlorination < par.minimum_high_chlorination_fraction)
            or np.any(chlorination > 1.0)
        )
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no declared Cl2+/chlorinated-Si channel")
        if outside_energy:
            reasons.append("Cl2+ energy leaves the measured 26--625 eV domain")
        if outside_angle:
            reasons.append("Cl2+ incidence leaves the measured normal-incidence domain")
        if outside_surface:
            reasons.append("poly-Si surface is not inside the declared high-chlorination scope")
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "Cl2+ incidence-angle response was not measured",
                "emitted SiClx product identity and differential distribution are unresolved",
                "high-chlorination transfer is a scope gate and does not rescale yield",
                "implantation, damage, and transient site density are unresolved",
            ),
            parameter_evidence_supports_prediction=True,
            nonpredictive_parameters=(),
        )

    def advance(
        self,
        state: ChangSawinArClSiState,
        fluxes: SurfaceFluxes,
        duration_s: float,
        *,
        surface_chlorination_fraction,
        strict=True,
    ):
        if not isinstance(state, ChangSawinArClSiState):
            raise TypeError("Balooch Cl2+ mechanism requires ChangSawinArClSiState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes, surface_chlorination_fraction)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))
        shape = state.removed_si_atoms_m2.shape
        ion_rate = np.zeros(shape)
        removal_rate = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl2+":
                continue
            face, event_flux, energy, _, face_count = self._events(population)
            yield_per_ion = self.parameters.yield_si_per_ion(
                energy, allow_extrapolation=not strict)
            if isinstance(population, FaceResolvedEnergeticFlux):
                ion_rate = ion_rate + np.broadcast_to(
                    np.bincount(face, weights=event_flux, minlength=face_count),
                    shape,
                )
                removal_rate = removal_rate + np.broadcast_to(
                    np.bincount(
                        face, weights=event_flux * yield_per_ion,
                        minlength=face_count,
                    ),
                    shape,
                )
            else:
                ion_rate = ion_rate + np.broadcast_to(np.sum(event_flux), shape)
                removal_rate = removal_rate + np.broadcast_to(
                    np.sum(event_flux * yield_per_ion), shape)
        duration = float(duration_s)
        removed = removal_rate * duration
        incident = ion_rate * duration
        updated = ChangSawinArClSiState(
            state.removed_si_atoms_m2 + removed,
            state.consumed_chlorine_atoms_m2,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed},
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=(
                "Figure 5.7 does not resolve the emitted SiClx product identity",
                "product energy-angle distribution is not measured",
            ),
        )
        products = (
            SurfaceProductPopulation(
                name="SiClx_from_Cl2plus_unresolved",
                source_inventory="Si_atom",
                integrated_particle_count_m2=removed,
                material_units_per_particle=1.0,
                mass_amu=28.0855,
                provenance={
                    "source": "chang-thesis Figure 5.7 / Balooch Cl2+ series",
                    "branch": "measured total Si removal",
                    "missing": "chlorine stoichiometry and differential emission",
                },
            ),
        )
        return BaloochCl2IonSiStepResult(
            state=updated,
            etch_velocity_m_s=(
                removal_rate / self.parameters.bulk_si_atom_density_m3),
            removed_si_atoms_m2=removed,
            incident_cl2plus_ions_m2=incident,
            material_exchange=exchange,
            product_populations=products,
            validity=validity,
        )
