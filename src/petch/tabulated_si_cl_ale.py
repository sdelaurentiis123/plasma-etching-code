"""Dose-resolved, element-balanced Si--Cl2--Ar+ ALE product replay.

The Kounis-Melas et al. DeepMD dataset releases five product-yield averages
for consecutive Ar+ dose windows during the ion-only half-cycle of a
prechlorinated Si surface.  Figures 13--14 of the accepted manuscript bind
those rows to 215 eV, normal-incidence Ar+; Supplementary Figures S1--S2
provide the SiCl2 and Cl branches.

This module treats the tabulated values as *window averages*.  It integrates
the overlap of an incident ion dose with each released window exactly instead
of inventing a smooth dose law between bin centres.  Si, SiCl, SiCl2, and Cl
products are retained separately, and both Si and Cl elemental ledgers close.
The caller must provide the initial retained-Cl inventory.  The mechanism
refuses a dose that would emit more Cl than that inventory or extend beyond
the released sequence.

This is an atom-balanced replay of a DeepMD result, not an experimental
validation or a first-principles feature model.  Product energy-angle
distributions, the preceding Cl2 loading transient, and continuation beyond
the released dose remain unknown.
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
from .surface_interaction_table import (
    SurfaceInteractionDomainError,
    SurfaceInteractionTable,
)
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    MechanismValidity,
    ParameterEvidence,
    SurfaceFluxes,
)


_PRODUCT_STOICHIOMETRY = MappingProxyType({
    "si_yield": ("Si", 1.0, 0.0, 28.085),
    "sicl_yield": ("SiCl", 1.0, 1.0, 63.538),
    "sicl2_yield": ("SiCl2", 1.0, 2.0, 98.991),
    "cl_yield": ("Cl", 0.0, 1.0, 35.453),
})


def _immutable_nonnegative(value, name):
    array = np.asarray(value, dtype=float).copy()
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    array.setflags(write=False)
    return array


def _balance_tolerance(*values):
    scale = np.ones(())
    for value in values:
        scale = np.maximum(scale, np.abs(np.asarray(value, dtype=float)))
    return 128.0 * np.finfo(float).eps * scale


@dataclass(frozen=True)
class TabulatedSiClAleState:
    """Local dose and elemental inventories per square metre."""

    ar_ion_dosage_m2: np.ndarray | float
    loaded_chlorine_atoms_m2: np.ndarray | float
    retained_chlorine_atoms_m2: np.ndarray | float
    emitted_chlorine_atoms_m2: np.ndarray | float = 0.0
    removed_si_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        supplied = np.broadcast_arrays(
            np.asarray(self.ar_ion_dosage_m2, dtype=float),
            np.asarray(self.loaded_chlorine_atoms_m2, dtype=float),
            np.asarray(self.retained_chlorine_atoms_m2, dtype=float),
            np.asarray(self.emitted_chlorine_atoms_m2, dtype=float),
            np.asarray(self.removed_si_atoms_m2, dtype=float),
        )
        names = (
            "ar_ion_dosage_m2",
            "loaded_chlorine_atoms_m2",
            "retained_chlorine_atoms_m2",
            "emitted_chlorine_atoms_m2",
            "removed_si_atoms_m2",
        )
        for name, value in zip(names, supplied):
            object.__setattr__(
                self, name, _immutable_nonnegative(value, name))
        residual = (
            self.loaded_chlorine_atoms_m2
            - self.retained_chlorine_atoms_m2
            - self.emitted_chlorine_atoms_m2
        )
        if np.any(np.abs(residual) > _balance_tolerance(
                self.loaded_chlorine_atoms_m2,
                self.retained_chlorine_atoms_m2,
                self.emitted_chlorine_atoms_m2)):
            raise ValueError("ALE chlorine inventory does not close")

    @classmethod
    def chlorinated(cls, retained_chlorine_atoms_m2, shape=None):
        retained = np.asarray(retained_chlorine_atoms_m2, dtype=float)
        if shape is not None:
            retained = np.broadcast_to(retained, shape)
        zero = np.zeros(retained.shape, dtype=float)
        return cls(zero, retained, retained, zero, zero)

    @property
    def shape(self):
        return self.ar_ion_dosage_m2.shape

    def conservative_surface_fields(self):
        return {
            "ar_ion_dosage_m2": self.ar_ion_dosage_m2,
            "loaded_chlorine_atoms_m2": self.loaded_chlorine_atoms_m2,
            "retained_chlorine_atoms_m2": self.retained_chlorine_atoms_m2,
            "emitted_chlorine_atoms_m2": self.emitted_chlorine_atoms_m2,
            "removed_si_atoms_m2": self.removed_si_atoms_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {name: None for name in self.conservative_surface_fields()}

    def surface_field_remap_modes(self):
        return {
            "ar_ion_dosage_m2": "intensive",
            "loaded_chlorine_atoms_m2": "conservative",
            "retained_chlorine_atoms_m2": "conservative",
            "emitted_chlorine_atoms_m2": "conservative",
            "removed_si_atoms_m2": "conservative",
        }

    def with_conservative_surface_fields(self, supplied):
        supplied = dict(supplied)
        if set(supplied) != set(self.conservative_surface_fields()):
            raise ValueError("Si--Cl ALE remap fields do not match its state")
        return type(self)(**supplied)


@dataclass(frozen=True)
class TabulatedSiClAleStepResult:
    state: TabulatedSiClAleState
    etch_velocity_m_s: np.ndarray
    removed_si_atoms_m2: np.ndarray
    emitted_chlorine_atoms_m2: np.ndarray
    product_counts_m2: Mapping[str, np.ndarray]
    product_count_standard_uncertainty_m2: Mapping[str, np.ndarray]
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    dose_bin_edges_m2: np.ndarray
    table_fingerprint: str
    validity: MechanismValidity

    def __post_init__(self):
        for name in (
                "etch_velocity_m_s",
                "removed_si_atoms_m2",
                "emitted_chlorine_atoms_m2",
                "dose_bin_edges_m2"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Si--Cl ALE step array")
            value.setflags(write=False)
            object.__setattr__(self, name, value)

        def freeze(mapping):
            output = {}
            for name, value in dict(mapping).items():
                array = _immutable_nonnegative(value, name)
                output[name] = array
            return MappingProxyType(output)

        object.__setattr__(
            self, "product_counts_m2", freeze(self.product_counts_m2))
        object.__setattr__(
            self, "product_count_standard_uncertainty_m2",
            freeze(self.product_count_standard_uncertainty_m2))
        populations = validate_surface_product_routing(
            self.material_exchange, self.product_populations)
        object.__setattr__(self, "product_populations", populations)


class TabulatedSiClAleProductMechanism:
    """Exact bin-integrated replay of the released 215 eV ALE products."""

    def __init__(
            self,
            interaction_table: SurfaceInteractionTable,
            bulk_si_atom_density_m3: float,
            bulk_density_evidence: ParameterEvidence,
            *,
            ion_species="Ar+",
            energy_tolerance_eV=1e-6,
            cosine_tolerance=1e-5):
        table = interaction_table
        expected_outputs = set(_PRODUCT_STOICHIOMETRY)
        conditions = dict(table.provenance.get("conditions", {}))
        if (table.material != "Si(100)"
                or table.incident_species != ("Ar+", "Cl2")
                or len(table.axes) != 1
                or table.axes[0].name != "ar_ion_dosage"
                or table.axes[0].unit != "1e15 cm^-2"
                or set(table.outputs) != expected_outputs
                or conditions.get("ar_ion_energy_eV") != 215.0
                or conditions.get("incidence_angle_deg") != 0.0):
            raise ValueError(
                "interaction table does not implement the 215 eV "
                "Si--Cl2--Ar+ ALE product contract")
        if (not np.isfinite(bulk_si_atom_density_m3)
                or bulk_si_atom_density_m3 <= 0.0
                or not isinstance(bulk_density_evidence, ParameterEvidence)
                or not ion_species
                or not np.isfinite(energy_tolerance_eV)
                or energy_tolerance_eV < 0.0
                or not np.isfinite(cosine_tolerance)
                or cosine_tolerance < 0.0):
            raise ValueError("invalid tabulated Si--Cl ALE inputs")

        centres_m2 = table.axes[0].values * 1.0e19
        spacing = np.diff(centres_m2)
        if (not np.allclose(spacing, spacing[0], rtol=1e-12, atol=0.0)
                or not np.isclose(
                    centres_m2[0], 0.5 * spacing[0],
                    rtol=1e-12, atol=0.0)):
            raise ValueError(
                "ALE product coordinates do not define contiguous "
                "equal-width dose windows starting at zero")
        edges = np.concatenate((
            np.array([0.0]),
            0.5 * (centres_m2[:-1] + centres_m2[1:]),
            np.array([centres_m2[-1] + 0.5 * spacing[-1]]),
        ))
        edges.setflags(write=False)
        self.table = table
        self.dose_bin_edges_m2 = edges
        self.bulk_si_atom_density_m3 = float(bulk_si_atom_density_m3)
        self.bulk_density_evidence = bulk_density_evidence
        self.ion_species = str(ion_species)
        self.energy_tolerance_eV = float(energy_tolerance_eV)
        self.cosine_tolerance = float(cosine_tolerance)

    @staticmethod
    def initial_state(retained_chlorine_atoms_m2, shape=None):
        return TabulatedSiClAleState.chlorinated(
            retained_chlorine_atoms_m2, shape=shape)

    @property
    def provenance(self):
        return MappingProxyType({
            "model": "Kounis-Melas 215 eV dose-window product replay",
            "table_fingerprint": self.table.fingerprint,
            "dataset_doi": self.table.provenance["dataset_doi"],
            "paper_doi": self.table.provenance["paper_doi"],
            "dose_interpretation": (
                "piecewise-constant released window averages; no interpolation"),
            "initial_condition": (
                "caller-supplied retained Cl inventory; no fitted loading dose"),
            "claim": (
                "atom-balanced DeepMD replay, not experiment or "
                "first-principles feature prediction"),
        })

    def validity(self, fluxes: SurfaceFluxes):
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != self.ion_species
            and np.any(np.asarray(item.flux_m2_s) > 0.0)
        }))
        wrong_energy = False
        wrong_angle = False
        for population in fluxes.energetic_fluxes:
            if population.name != self.ion_species:
                continue
            if isinstance(population, FaceResolvedEnergeticFlux):
                selected = population.event_flux_m2_s > 0.0
                energy = population.event_energy_eV[selected]
                cosine = population.event_cosine_incidence[selected]
            elif isinstance(population, EnergeticFlux):
                selected = population.weight > 0.0
                energy = population.energy_eV[selected]
                cosine = population.cosine_incidence[selected]
            else:  # pragma: no cover
                raise TypeError(type(population).__name__)
            wrong_energy |= bool(np.any(
                np.abs(energy - 215.0) > self.energy_tolerance_eV))
            wrong_angle |= bool(np.any(
                np.abs(cosine - 1.0) > self.cosine_tolerance))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append(
                "positive incident flux has no released 215 eV ALE "
                "bombardment-stage channel")
        if wrong_energy:
            reasons.append("ion energy leaves the fixed 215 eV product table")
        if wrong_angle:
            reasons.append(
                "ion incidence leaves the fixed normal-incidence product table")
        nonpredictive = []
        if not self.bulk_density_evidence.supports_prediction_within_declared_domain:
            nonpredictive.append("bulk_si_atom_density_m3")
        if self.table.provenance.get(
                "supports_prediction_within_declared_domain") is not True:
            nonpredictive.append("interaction_table")
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "the preceding Cl2 loading transient is not in Products.csv",
                "initial retained-Cl inventory must be supplied independently",
                "dose-window averages do not resolve within-window transients",
                "product energy-angle distributions are absent",
                "the sequence cannot continue beyond the released final dose",
                "DeepMD evidence is not a direct experimental measurement",
            ),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=tuple(nonpredictive))

    def _integrated_counts(self, start_m2, stop_m2):
        start = np.asarray(start_m2, dtype=float)
        stop = np.asarray(stop_m2, dtype=float)
        counts = {
            name: np.zeros(start.shape, dtype=float)
            for name in _PRODUCT_STOICHIOMETRY
        }
        variance = {
            name: np.zeros(start.shape, dtype=float)
            for name in _PRODUCT_STOICHIOMETRY
        }
        for index, (left, right) in enumerate(zip(
                self.dose_bin_edges_m2[:-1],
                self.dose_bin_edges_m2[1:])):
            overlap = np.maximum(
                np.minimum(stop, right) - np.maximum(start, left), 0.0)
            for name in counts:
                counts[name] += overlap * self.table.outputs[name][index]
                uncertainty = self.table.standard_uncertainty[name][index]
                variance[name] += (overlap * uncertainty) ** 2
        return counts, {
            name: np.sqrt(value) for name, value in variance.items()
        }

    @staticmethod
    def _ion_flux(population, shape):
        return np.broadcast_to(
            np.asarray(population.flux_m2_s, dtype=float), shape)

    def advance(
            self,
            state: TabulatedSiClAleState,
            fluxes: SurfaceFluxes,
            duration_s: float,
            *,
            strict=True):
        if not isinstance(state, TabulatedSiClAleState):
            raise TypeError(
                "Si--Cl ALE mechanism requires TabulatedSiClAleState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))
        ion_flux = np.zeros(state.shape, dtype=float)
        for population in fluxes.energetic_fluxes:
            if population.name == self.ion_species:
                ion_flux += self._ion_flux(population, state.shape)
        added_dose = ion_flux * float(duration_s)
        final_dose = state.ar_ion_dosage_m2 + added_dose
        maximum = self.dose_bin_edges_m2[-1]
        tolerance = 128.0 * np.finfo(float).eps * max(1.0, maximum)
        if np.any(final_dose > maximum + tolerance):
            raise SurfaceInteractionDomainError(
                "Ar+ dose exceeds the released ALE product sequence")
        final_dose = np.minimum(final_dose, maximum)
        counts, count_uncertainty = self._integrated_counts(
            state.ar_ion_dosage_m2, final_dose)
        removed_si = sum(
            counts[name] * stoichiometry[1]
            for name, stoichiometry in _PRODUCT_STOICHIOMETRY.items())
        emitted_cl = sum(
            counts[name] * stoichiometry[2]
            for name, stoichiometry in _PRODUCT_STOICHIOMETRY.items())
        if np.any(emitted_cl > (
                state.retained_chlorine_atoms_m2
                + _balance_tolerance(
                    emitted_cl, state.retained_chlorine_atoms_m2))):
            raise ValueError(
                "released products require more Cl than the supplied "
                "retained-layer inventory")
        retained_cl = np.maximum(
            state.retained_chlorine_atoms_m2 - emitted_cl, 0.0)
        updated = TabulatedSiClAleState(
            final_dose,
            state.loaded_chlorine_atoms_m2,
            retained_cl,
            state.emitted_chlorine_atoms_m2 + emitted_cl,
            state.removed_si_atoms_m2 + removed_si,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={
                "Si_atom": removed_si,
                "Cl_atom": emitted_cl,
            },
            outgoing_units_m2={
                "Si_atom": removed_si,
                "Cl_atom": emitted_cl,
            },
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=(
                "product identities are known but launch distributions are not",
                "Cl is drawn only from the caller-supplied preloaded inventory",
            ),
        )
        products = []
        for output_name, (
                species_name, si_per_particle, cl_per_particle, mass_amu
        ) in _PRODUCT_STOICHIOMETRY.items():
            if si_per_particle:
                source = "Si_atom"
                primary_units = si_per_particle
                additional = (
                    {"Cl_atom": cl_per_particle}
                    if cl_per_particle else {})
            else:
                source = "Cl_atom"
                primary_units = cl_per_particle
                additional = {}
            products.append(SurfaceProductPopulation(
                name=species_name,
                source_inventory=source,
                integrated_particle_count_m2=counts[output_name],
                material_units_per_particle=primary_units,
                mass_amu=mass_amu,
                provenance={
                    "source": self.table.provenance["source"],
                    "dataset_doi": self.table.provenance["dataset_doi"],
                    "table_fingerprint": self.table.fingerprint,
                    "ion_energy_eV": 215.0,
                    "dose_semantics": "released window-average yield",
                    "missing": (
                        "differential emitted energy-angle distribution"),
                },
                additional_source_inventories_per_particle=additional,
            ))
        duration = float(duration_s)
        velocity = (
            removed_si / self.bulk_si_atom_density_m3 / duration
            if duration > 0.0 else np.zeros(state.shape, dtype=float))
        return TabulatedSiClAleStepResult(
            state=updated,
            etch_velocity_m_s=velocity,
            removed_si_atoms_m2=removed_si,
            emitted_chlorine_atoms_m2=emitted_cl,
            product_counts_m2=counts,
            product_count_standard_uncertainty_m2=count_uncertainty,
            material_exchange=exchange,
            product_populations=tuple(products),
            dose_bin_edges_m2=self.dose_bin_edges_m2,
            table_fingerprint=self.table.fingerprint,
            validity=validity,
        )
