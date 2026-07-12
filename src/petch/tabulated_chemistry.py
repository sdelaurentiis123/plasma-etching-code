"""Stateful surface mechanisms driven by sourced interaction tables."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .surface_interaction_table import SurfaceInteractionTable
from .surface_kinetics import (
    EnergeticFlux, FaceResolvedEnergeticFlux, MechanismValidity, ParameterEvidence,
    SurfaceFluxes,
)
from .surface_exchange import (
    SurfaceMaterialExchange, SurfaceProductPopulation, unresolved_surface_exchange,
    validate_surface_product_routing,
)


@dataclass(frozen=True)
class TabulatedSiSurfaceState:
    removed_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        value = np.asarray(self.removed_atoms_m2, dtype=float).copy()
        if np.any(~np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError("removed Si inventory must be finite and nonnegative")
        value.setflags(write=False); object.__setattr__(self, "removed_atoms_m2", value)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.zeros(shape))

    def conservative_surface_fields(self):
        return {"removed_atoms_m2": self.removed_atoms_m2}

    def conservative_surface_upper_bounds(self):
        return {"removed_atoms_m2": None}

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != {"removed_atoms_m2"}:
            raise ValueError("tabulated Si remap fields do not match its state contract")
        return type(self)(fields["removed_atoms_m2"])


@dataclass(frozen=True)
class TabulatedSiSurfaceStepResult:
    state: TabulatedSiSurfaceState
    etch_velocity_m_s: np.ndarray
    etch_velocity_standard_uncertainty_m_s: np.ndarray
    removed_atoms_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    table_fingerprint: str
    validity: MechanismValidity

    def __post_init__(self):
        for name in (
                "etch_velocity_m_s", "etch_velocity_standard_uncertainty_m_s",
                "removed_atoms_m2"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            value.setflags(write=False); object.__setattr__(self, name, value)


@dataclass(frozen=True)
class TabulatedSiPhysicalSputterStepResult:
    state: TabulatedSiSurfaceState
    etch_velocity_m_s: np.ndarray
    etch_velocity_standard_uncertainty_m_s: np.ndarray
    removed_atoms_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    table_fingerprint: str
    validity: MechanismValidity

    def __post_init__(self):
        for name in (
                "etch_velocity_m_s", "etch_velocity_standard_uncertainty_m_s",
                "removed_atoms_m2"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            value.setflags(write=False); object.__setattr__(self, name, value)
        populations = validate_surface_product_routing(
            self.material_exchange, self.product_populations)
        object.__setattr__(self, "product_populations", populations)


class TabulatedSiPhysicalSputterMechanism:
    """Normal-incidence Ar+ physical sputtering of Si from the pinned DeepMD table.

    The table supplies energy-dependent total Si yield, not a differential energy/angular distribution
    for emitted Si. Removed material is therefore routed exactly to neutral Si products but those products
    remain ineligible for transport until that missing distribution is supplied from evidence.
    """

    def __init__(
            self, interaction_table: SurfaceInteractionTable, bulk_atom_density_m3: float,
            bulk_density_evidence: ParameterEvidence, *, ion_species="Ar+",
            cosine_tolerance=1e-5):
        table = interaction_table
        if (table.material != "Si(100)" or table.incident_species != ("Ar+",)
                or len(table.axes) != 1 or table.axes[0].name != "ion_energy"
                or set(table.outputs) != {
                    "physical_sputter_yield", "amorphous_layer_thickness"}
                or table.output_units["physical_sputter_yield"] != "Si/Ar+"):
            raise ValueError("interaction table does not implement the Si Ar+ sputter contract")
        conditions = dict(table.provenance.get("conditions", {}))
        if conditions.get("incidence_angle_deg") != 0.0:
            raise ValueError("Si sputter table must declare normal incidence")
        if (not np.isfinite(bulk_atom_density_m3) or bulk_atom_density_m3 <= 0.0
                or not isinstance(bulk_density_evidence, ParameterEvidence)
                or not ion_species or not np.isfinite(cosine_tolerance)
                or cosine_tolerance < 0.0):
            raise ValueError("invalid tabulated physical-sputter inputs")
        self.table = table
        self.bulk_atom_density_m3 = float(bulk_atom_density_m3)
        self.bulk_density_evidence = bulk_density_evidence
        self.ion_species = str(ion_species)
        self.cosine_tolerance = float(cosine_tolerance)

    @staticmethod
    def initial_state(shape=()):
        return TabulatedSiSurfaceState.bare(shape)

    def validity(self, fluxes: SurfaceFluxes):
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != self.ion_species and np.any(np.asarray(item.flux_m2_s) > 0.0)}))
        wrong_angle = False
        for population in fluxes.energetic_fluxes:
            if population.name != self.ion_species:
                continue
            cosine = (population.event_cosine_incidence
                      if isinstance(population, FaceResolvedEnergeticFlux)
                      else population.cosine_incidence)
            wrong_angle |= bool(np.any(np.abs(cosine - 1.0) > self.cosine_tolerance))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no physical-sputter table channel")
        if wrong_angle:
            reasons.append("ion incidence leaves the fixed normal-incidence sputter condition")
        nonpredictive = []
        if not self.bulk_density_evidence.supports_prediction_within_declared_domain:
            nonpredictive.append("bulk_atom_density_m3")
        if self.table.provenance.get("supports_prediction_within_declared_domain") is not True:
            nonpredictive.append("interaction_table")
        return MechanismValidity(
            within_declared_scope=not reasons, reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "physical-sputter table has no incidence-angle sweep",
                "amorphous-layer thickness is reported but not evolved as a damage state",
                "emitted Si energy and angular distributions are absent from the source table",
            ),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=tuple(nonpredictive))

    def _rate_and_uncertainty(self, population, shape):
        if isinstance(population, FaceResolvedEnergeticFlux):
            evaluated = self.table.evaluate({"ion_energy": population.event_energy_eV})
            rate = np.bincount(
                population.event_face,
                weights=(population.event_flux_m2_s
                         * evaluated.values["physical_sputter_yield"]),
                minlength=population.face_count)
            # Source-table uncertainty is epistemic; use the fully correlated linear propagation.
            uncertainty = np.bincount(
                population.event_face,
                weights=(population.event_flux_m2_s
                         * evaluated.standard_uncertainty["physical_sputter_yield"]),
                minlength=population.face_count)
            return np.broadcast_to(rate, shape), np.broadcast_to(uncertainty, shape)
        evaluated = self.table.evaluate({"ion_energy": population.energy_eV})
        mean_yield = float(np.dot(
            population.weight, evaluated.values["physical_sputter_yield"]))
        mean_uncertainty = float(np.dot(
            population.weight,
            evaluated.standard_uncertainty["physical_sputter_yield"]))
        flux = np.broadcast_to(np.asarray(population.flux_m2_s, dtype=float), shape)
        return flux * mean_yield, flux * mean_uncertainty

    def advance(self, state, fluxes: SurfaceFluxes, duration_s: float, *, strict=True):
        if not isinstance(state, TabulatedSiSurfaceState):
            raise TypeError("Si physical sputtering requires TabulatedSiSurfaceState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))
        shape = state.removed_atoms_m2.shape
        removal_rate = np.zeros(shape); removal_uncertainty = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name == self.ion_species:
                rate, uncertainty = self._rate_and_uncertainty(population, shape)
                removal_rate += rate; removal_uncertainty += uncertainty
        removed = removal_rate * float(duration_s)
        updated = TabulatedSiSurfaceState(state.removed_atoms_m2 + removed)
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed}, unresolved_units_m2={},
            deposited_units_m2={}, known_limitations=(
                "outgoing Si is materially identified but lacks a sourced launch distribution",))
        products = (SurfaceProductPopulation(
            "Si", "Si_atom", removed, 1.0, 28.085,
            provenance={
                "source": self.table.provenance["source"],
                "table_fingerprint": self.table.fingerprint,
                "evidence_type": self.table.provenance["evidence_type"],
                "missing": "differential emitted energy-angle distribution",
            }),)
        return TabulatedSiPhysicalSputterStepResult(
            state=updated,
            etch_velocity_m_s=removal_rate / self.bulk_atom_density_m3,
            etch_velocity_standard_uncertainty_m_s=(
                removal_uncertainty / self.bulk_atom_density_m3),
            removed_atoms_m2=removed, material_exchange=exchange,
            product_populations=products, table_fingerprint=self.table.fingerprint,
            validity=validity)


class TabulatedSiClArMechanism:
    """Si-Cl2-Ar+ RIE at the exact fixed conditions released in OSTI 2589032.

    The archived RIE table varies only Cl2:Ar+ incident flux ratio. The mechanism therefore refuses
    non-100-eV ions, non-normal impacts, and ratios outside 10--200. It does not invent angle, energy,
    coverage, or temperature dependence absent from the source.
    """

    def __init__(
            self, interaction_table: SurfaceInteractionTable, bulk_atom_density_m3: float,
            bulk_density_evidence: ParameterEvidence, *, ion_species="Ar+", neutral_species="Cl2",
            energy_tolerance_eV=1e-6, cosine_tolerance=1e-5):
        table = interaction_table
        if (table.material != "Si(100)" or table.incident_species != ("Ar+", "Cl2")
                or len(table.axes) != 1 or table.axes[0].name != "cl2_to_ar_flux_ratio"
                or set(table.outputs) != {"reactive_etch_yield"}
                or table.output_units["reactive_etch_yield"] != "Si/Ar+"):
            raise ValueError("interaction table does not implement the Si-Cl2-Ar+ RIE contract")
        conditions = dict(table.provenance.get("conditions", {}))
        if (conditions.get("ar_ion_energy_eV") != 100.0
                or conditions.get("incidence_angle_deg") != 0.0):
            raise ValueError("Si-Cl2-Ar+ RIE table must declare its fixed energy and incidence")
        if (not np.isfinite(bulk_atom_density_m3) or bulk_atom_density_m3 <= 0.0
                or not isinstance(bulk_density_evidence, ParameterEvidence)
                or not ion_species or not neutral_species
                or not np.isfinite(energy_tolerance_eV) or energy_tolerance_eV < 0.0
                or not np.isfinite(cosine_tolerance) or cosine_tolerance < 0.0):
            raise ValueError("invalid tabulated Si mechanism inputs")
        self.table = table
        self.bulk_atom_density_m3 = float(bulk_atom_density_m3)
        self.bulk_density_evidence = bulk_density_evidence
        self.ion_species = str(ion_species); self.neutral_species = str(neutral_species)
        self.energy_tolerance_eV = float(energy_tolerance_eV)
        self.cosine_tolerance = float(cosine_tolerance)

    @staticmethod
    def initial_state(shape=()):
        return TabulatedSiSurfaceState.bare(shape)

    @staticmethod
    def _population_measure(population, shape):
        flux = np.asarray(population.flux_m2_s, dtype=float)
        return np.broadcast_to(flux, shape)

    def validity(self, fluxes: SurfaceFluxes):
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name != self.neutral_species and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != self.ion_species and np.any(np.asarray(item.flux_m2_s) > 0.0)}))
        positive_ion = False; wrong_energy = False; wrong_angle = False
        for population in fluxes.energetic_fluxes:
            if population.name != self.ion_species:
                continue
            positive_ion |= bool(np.any(np.asarray(population.flux_m2_s) > 0.0))
            if isinstance(population, FaceResolvedEnergeticFlux):
                selected = population.event_flux_m2_s > 0.0
                energy = population.event_energy_eV[selected]
                cosine = population.event_cosine_incidence[selected]
            elif isinstance(population, EnergeticFlux):
                selected = population.weight > 0.0
                energy = population.energy_eV[selected]
                cosine = population.cosine_incidence[selected]
            else:  # pragma: no cover - SurfaceFluxes validates population types
                raise TypeError(type(population).__name__)
            wrong_energy |= bool(np.any(np.abs(energy - 100.0) > self.energy_tolerance_eV))
            wrong_angle |= bool(np.any(np.abs(cosine - 1.0) > self.cosine_tolerance))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no Si-Cl2-Ar+ table channel")
        positive_neutral = bool(np.any(np.asarray(
            fluxes.neutral_flux_m2_s.get(self.neutral_species, 0.0)) > 0.0))
        if positive_ion != positive_neutral:
            reasons.append(
                "archived RIE table requires simultaneous positive Cl2 and Ar+ incident flux")
        if wrong_energy:
            reasons.append("ion energy leaves the fixed 100 eV RIE table condition")
        if wrong_angle:
            reasons.append("ion incidence leaves the fixed normal-incidence RIE table condition")
        nonpredictive = []
        if not self.bulk_density_evidence.supports_prediction_within_declared_domain:
            nonpredictive.append("bulk_atom_density_m3")
        if self.table.provenance.get("supports_prediction_within_declared_domain") is not True:
            nonpredictive.append("interaction_table")
        return MechanismValidity(
            within_declared_scope=not reasons, reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "RIE table has no incidence-angle or ion-energy sweep",
                f"normal incidence is accepted only within cosine tolerance {self.cosine_tolerance:g}",
                "surface coverage and damage memory are implicit in archived steady RIE yields",
                "etch-product branching is available only for the separate 80 eV ALE table",
            ),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=tuple(nonpredictive))

    def advance(self, state, fluxes: SurfaceFluxes, duration_s: float, *, strict=True):
        if not isinstance(state, TabulatedSiSurfaceState):
            raise TypeError("Si-Cl2-Ar+ mechanism requires TabulatedSiSurfaceState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))
        shape = state.removed_atoms_m2.shape
        ion_flux = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name == self.ion_species:
                ion_flux = ion_flux + self._population_measure(population, shape)
        neutral_flux = np.broadcast_to(
            np.asarray(fluxes.neutral_flux_m2_s.get(self.neutral_species, 0.0), dtype=float), shape)
        active = ion_flux > 0.0
        ratio = np.zeros(shape); ratio[active] = neutral_flux[active] / ion_flux[active]
        yield_value = np.zeros(shape); yield_uncertainty = np.zeros(shape)
        if np.any(active):
            evaluated = self.table.evaluate({
                "cl2_to_ar_flux_ratio": ratio[active]})
            yield_value[active] = evaluated.values["reactive_etch_yield"]
            yield_uncertainty[active] = evaluated.standard_uncertainty[
                "reactive_etch_yield"]
        removal_rate = ion_flux * yield_value
        removal_uncertainty = ion_flux * yield_uncertainty
        removed = removal_rate * float(duration_s)
        updated = TabulatedSiSurfaceState(state.removed_atoms_m2 + removed)
        exchange = unresolved_surface_exchange(
            removed_units_m2={"Si_atom": removed},
            limitations=(
                "reactive Si-Cl2-Ar+ product branching is absent from the released RIE table",
                "unresolved removed Si is not eligible for redeposition transport",
            ))
        return TabulatedSiSurfaceStepResult(
            state=updated,
            etch_velocity_m_s=removal_rate / self.bulk_atom_density_m3,
            etch_velocity_standard_uncertainty_m_s=(
                removal_uncertainty / self.bulk_atom_density_m3),
            removed_atoms_m2=removed, material_exchange=exchange,
            table_fingerprint=self.table.fingerprint,
            validity=validity)
