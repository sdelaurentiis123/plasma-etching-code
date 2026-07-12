"""Dimensional, stateful surface kinetics for feature-scale etching.

This module separates the reusable numerical/physical kernel from chemistry parameter sets.  A reaction
topology may be sourced from literature while its probabilities remain measured, calculated, calibrated,
or unknown inputs.  Nothing in this module silently promotes a fitted coefficient to a universal constant.

The energetic law and reduced reaction topology follow F. Krueger, *Modeling and Optimization of High
Aspect Ratio Plasma Etching* (2024), Eq. 2.40, Sec. 6.4, and Appendix B, DOI 10.7302/23106.  The named
angular forms refer to J. P. Chang and H. H. Sawin, DOI 10.1116/1.580692, and J. D. Kress et al., DOI
10.1116/1.581948.  Their numerical coefficients remain explicit inputs.

The first mechanism is a deliberately reduced SiO2/fluorocarbon network.  It resolves oxide-fluorocarbon
complex coverage and a finite polymer inventory, and it conservatively accounts for complex formation,
bare/complex energetic removal, polymer deposition, and O/energetic polymer removal.  Crosslinking,
species-resolved complex stoichiometry, redeposition, and mask chemistry remain outside this reduced
mechanism and are reported as model-form omissions rather than hidden in an effective rate scale.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np


@dataclass(frozen=True)
class ParameterEvidence:
    """Provenance and uncertainty for one physical chemistry input."""

    source: str
    evidence_type: str
    relative_standard_uncertainty: float | None = None
    note: str = ""

    def __post_init__(self):
        if not self.source or not self.evidence_type:
            raise ValueError("parameter evidence requires source and evidence_type")
        uncertainty = self.relative_standard_uncertainty
        if uncertainty is not None and (not np.isfinite(uncertainty) or uncertainty < 0.0):
            raise ValueError("relative standard uncertainty must be finite and nonnegative")


@dataclass(frozen=True)
class EnergeticYield:
    """Threshold energy-angle yield law used for ions and energetic neutrals.

    The energy factor follows Krueger's MCFPM Eq. (2.40):
    ``Y = Y_ref * max((E-E_th)/(E_ref-E_th), 0)**n * f(theta)``.
    It is a yield, not necessarily a Bernoulli probability, so values are not clipped to one.
    """

    reference_yield: float
    threshold_energy_eV: float
    reference_energy_eV: float
    energy_exponent: float = 1.0
    angular_model: str = "none"
    angular_parameter: float | None = None

    def __post_init__(self):
        if (not np.isfinite(self.reference_yield) or self.reference_yield < 0.0
                or not np.isfinite(self.threshold_energy_eV) or self.threshold_energy_eV < 0.0
                or not np.isfinite(self.reference_energy_eV)
                or self.reference_energy_eV <= self.threshold_energy_eV
                or not np.isfinite(self.energy_exponent) or self.energy_exponent <= 0.0):
            raise ValueError("invalid energetic-yield parameters")
        if self.angular_model not in {"none", "chang_sawin_1997", "kress_1999"}:
            raise ValueError(f"unknown angular yield model: {self.angular_model}")
        if self.angular_model == "kress_1999":
            if (self.angular_parameter is None or not np.isfinite(self.angular_parameter)
                    or self.angular_parameter < 0.0):
                raise ValueError("kress_1999 requires a nonnegative angular_parameter")

    def evaluate(self, energy_eV, cosine_incidence):
        energy = np.asarray(energy_eV, dtype=float)
        cosine = np.asarray(cosine_incidence, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
            raise ValueError("incident energies must be finite and nonnegative")
        if np.any(~np.isfinite(cosine)) or np.any((cosine < 0.0) | (cosine > 1.0)):
            raise ValueError("incidence cosines must lie in [0, 1]")
        scaled_energy = np.maximum(
            (energy - self.threshold_energy_eV)
            / (self.reference_energy_eV - self.threshold_energy_eV), 0.0)
        if self.angular_model == "none":
            angular = np.ones(np.broadcast_shapes(energy.shape, cosine.shape))
        elif self.angular_model == "chang_sawin_1997":
            theta = np.arccos(cosine)
            angular = np.where(
                cosine >= 0.5, 1.0, np.maximum(3.0 - 6.0 * theta / np.pi, 0.0))
        else:
            angular = ((1.0 + float(self.angular_parameter) * (1.0 - cosine * cosine))
                       * cosine)
        return self.reference_yield * scaled_energy ** self.energy_exponent * angular


@dataclass(frozen=True)
class EnergeticFlux:
    """One incident ion/hot-neutral population with a normalized joint energy-angle rule."""

    name: str
    flux_m2_s: np.ndarray | float
    energy_eV: np.ndarray
    cosine_incidence: np.ndarray
    weight: np.ndarray

    def __post_init__(self):
        if not self.name:
            raise ValueError("energetic flux requires a species name")
        flux = np.asarray(self.flux_m2_s, dtype=float).copy()
        energy = np.asarray(self.energy_eV, dtype=float).copy()
        cosine = np.asarray(self.cosine_incidence, dtype=float).copy()
        weight = np.asarray(self.weight, dtype=float).copy()
        if (np.any(~np.isfinite(flux)) or np.any(flux < 0.0)
                or energy.ndim != 1 or cosine.shape != energy.shape or weight.shape != energy.shape
                or np.any(~np.isfinite(energy)) or np.any(energy < 0.0)
                or np.any(~np.isfinite(cosine)) or np.any((cosine < 0.0) | (cosine > 1.0))
                or np.any(~np.isfinite(weight)) or np.any(weight < 0.0) or weight.sum() <= 0.0):
            raise ValueError("invalid energetic flux")
        weight /= weight.sum()
        for array in (flux, energy, cosine, weight):
            array.setflags(write=False)
        object.__setattr__(self, "flux_m2_s", flux)
        object.__setattr__(self, "energy_eV", energy)
        object.__setattr__(self, "cosine_incidence", cosine)
        object.__setattr__(self, "weight", weight)

    def mean_yield(self, yield_law: EnergeticYield):
        values = yield_law.evaluate(self.energy_eV, self.cosine_incidence)
        return float(np.dot(self.weight, values))


@dataclass(frozen=True)
class SurfaceFluxes:
    """Dimensional incident fluxes at one or many surface elements."""

    neutral_flux_m2_s: Mapping[str, np.ndarray | float]
    energetic_fluxes: tuple[EnergeticFlux, ...] = ()

    def __post_init__(self):
        neutral = {}
        for name, value in self.neutral_flux_m2_s.items():
            if not name:
                raise ValueError("neutral species names must be nonempty")
            array = np.asarray(value, dtype=float).copy()
            if np.any(~np.isfinite(array)) or np.any(array < 0.0):
                raise ValueError(f"invalid neutral flux for {name}")
            array.setflags(write=False); neutral[name] = array
        object.__setattr__(self, "neutral_flux_m2_s", MappingProxyType(neutral))
        energetic = tuple(self.energetic_fluxes)
        if any(not isinstance(item, EnergeticFlux) for item in energetic):
            raise TypeError("energetic_fluxes must contain EnergeticFlux objects")
        object.__setattr__(self, "energetic_fluxes", energetic)


@dataclass(frozen=True)
class SiO2SurfaceState:
    """Reduced per-area surface state.

    ``complex_fraction`` is the fraction of accessible oxide sites converted to an
    oxide-fluorocarbon complex. ``polymer_units_m2`` is a physical areal inventory, not a coverage
    fit parameter. ``removed_formula_units_m2`` integrates removed SiO2 formula units.
    """

    complex_fraction: np.ndarray | float
    polymer_units_m2: np.ndarray | float
    removed_formula_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        arrays = np.broadcast_arrays(
            np.asarray(self.complex_fraction, dtype=float),
            np.asarray(self.polymer_units_m2, dtype=float),
            np.asarray(self.removed_formula_units_m2, dtype=float))
        complex_fraction, polymer, removed = [np.array(item, copy=True) for item in arrays]
        if (np.any(~np.isfinite(complex_fraction))
                or np.any((complex_fraction < 0.0) | (complex_fraction > 1.0))
                or np.any(~np.isfinite(polymer)) or np.any(polymer < 0.0)
                or np.any(~np.isfinite(removed)) or np.any(removed < 0.0)):
            raise ValueError("invalid SiO2 surface state")
        for array in (complex_fraction, polymer, removed): array.setflags(write=False)
        object.__setattr__(self, "complex_fraction", complex_fraction)
        object.__setattr__(self, "polymer_units_m2", polymer)
        object.__setattr__(self, "removed_formula_units_m2", removed)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.zeros(shape), np.zeros(shape), np.zeros(shape))


@dataclass(frozen=True)
class ReducedSiO2FluorocarbonParameters:
    """Physical inputs to the reduced SiO2/fluorocarbon mechanism.

    Species maps make chemistry data inputs rather than branches in the transport/interface engine.
    Every nonzero incident neutral must be assigned a physical channel or the validity check refuses it.
    """

    site_density_m2: float
    bulk_formula_density_m3: float
    polymer_monolayer_density_m2: float
    complex_formation_probability: Mapping[str, float]
    polymer_deposition_probability_on_substrate: Mapping[str, float]
    polymer_deposition_probability_on_polymer: Mapping[str, float]
    oxygen_species: str
    oxygen_polymer_etch_probability: float
    bare_sio2_yield: EnergeticYield
    complex_sio2_yield: EnergeticYield
    polymer_sputter_yield: EnergeticYield
    evidence: Mapping[str, ParameterEvidence] = field(default_factory=dict)
    known_omissions: tuple[str, ...] = (
        "polymer_crosslinking",
        "species_resolved_oxide_fluorocarbon_complexes",
        "surface_fluorination_sequence",
        "etch_product_redeposition",
        "amorphous_carbon_mask_chemistry",
    )

    def __post_init__(self):
        for name, value in (
                ("site_density_m2", self.site_density_m2),
                ("bulk_formula_density_m3", self.bulk_formula_density_m3),
                ("polymer_monolayer_density_m2", self.polymer_monolayer_density_m2)):
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if (not self.oxygen_species or not np.isfinite(self.oxygen_polymer_etch_probability)
                or not 0.0 <= self.oxygen_polymer_etch_probability <= 1.0):
            raise ValueError("invalid oxygen polymer-etch input")
        maps = {}
        for name in (
                "complex_formation_probability",
                "polymer_deposition_probability_on_substrate",
                "polymer_deposition_probability_on_polymer"):
            values = dict(getattr(self, name))
            if any((not species or not np.isfinite(value) or value < 0.0 or value > 1.0)
                   for species, value in values.items()):
                raise ValueError(f"invalid probability map: {name}")
            maps[name] = MappingProxyType(values)
        for name, values in maps.items(): object.__setattr__(self, name, values)
        evidence = dict(self.evidence)
        if any(not isinstance(item, ParameterEvidence) for item in evidence.values()):
            raise TypeError("parameter evidence values must be ParameterEvidence objects")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))


@dataclass(frozen=True)
class MechanismValidity:
    within_declared_scope: bool
    reasons: tuple[str, ...]
    unsupported_neutral_species: tuple[str, ...]
    known_model_form_omissions: tuple[str, ...]


@dataclass(frozen=True)
class SurfaceStepResult:
    state: SiO2SurfaceState
    etch_velocity_m_s: np.ndarray
    formed_complex_units_m2: np.ndarray
    removed_complex_units_m2: np.ndarray
    removed_bare_formula_units_m2: np.ndarray
    deposited_polymer_units_m2: np.ndarray
    removed_polymer_units_m2: np.ndarray
    validity: MechanismValidity


class ReducedSiO2FluorocarbonMechanism:
    """Vectorized conservative kernel for the declared reduced reaction network."""

    def __init__(self, parameters: ReducedSiO2FluorocarbonParameters):
        self.parameters = parameters

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        supported = (set(par.complex_formation_probability)
                     | set(par.polymer_deposition_probability_on_substrate)
                     | set(par.polymer_deposition_probability_on_polymer)
                     | {par.oxygen_species})
        unsupported = tuple(sorted(
            name for name, flux in fluxes.neutral_flux_m2_s.items()
            if np.any(np.asarray(flux) > 0.0) and name not in supported))
        required_evidence = {
            "site_density_m2", "bulk_formula_density_m3", "polymer_monolayer_density_m2",
            "complex_formation_probability", "polymer_deposition_probability_on_substrate",
            "polymer_deposition_probability_on_polymer", "oxygen_polymer_etch_probability",
            "bare_sio2_yield", "complex_sio2_yield", "polymer_sputter_yield",
        }
        missing_evidence = tuple(sorted(required_evidence - set(par.evidence)))
        reasons = []
        if unsupported:
            reasons.append("positive incident neutral flux has no declared reaction channel")
        if missing_evidence:
            reasons.append("missing parameter evidence: " + ", ".join(missing_evidence))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported,
            known_model_form_omissions=par.known_omissions)

    @staticmethod
    def _broadcast(value, shape):
        return np.broadcast_to(np.asarray(value, dtype=float), shape)

    def _neutral_weighted_rate(self, fluxes, probability, shape):
        total = np.zeros(shape)
        for species, reaction_probability in probability.items():
            if species in fluxes.neutral_flux_m2_s:
                total = total + self._broadcast(
                    fluxes.neutral_flux_m2_s[species], shape) * reaction_probability
        return total

    def _energetic_rate(self, fluxes, yield_law, shape):
        total = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            total = total + self._broadcast(population.flux_m2_s, shape) * population.mean_yield(
                yield_law)
        return total

    def _polymer_step(self, inventory, fluxes, duration_s, shape):
        """Exact constant-flux solution and exact integrated deposition/removal bookkeeping."""
        par = self.parameters; monolayer = par.polymer_monolayer_density_m2
        deposit_substrate = self._neutral_weighted_rate(
            fluxes, par.polymer_deposition_probability_on_substrate, shape)
        deposit_polymer = self._neutral_weighted_rate(
            fluxes, par.polymer_deposition_probability_on_polymer, shape)
        oxygen_flux = self._broadcast(
            fluxes.neutral_flux_m2_s.get(par.oxygen_species, 0.0), shape)
        removal_capacity = (oxygen_flux * par.oxygen_polymer_etch_probability
                            + self._energetic_rate(fluxes, par.polymer_sputter_yield, shape))
        if not np.any(deposit_substrate + deposit_polymer + removal_capacity > 0.0):
            return np.array(inventory, copy=True), np.zeros(shape), np.zeros(shape)

        # dN/dt = D_sub + (D_poly-D_sub-R) * (1-exp(-N/N_mono)).  With y=exp(N/N_mono)
        # this is a linear ODE, so positivity and the N=0 boundary are preserved without clipping.
        coefficient = deposit_polymer - removal_capacity
        transition = deposit_polymer - deposit_substrate - removal_capacity
        exponent = coefficient * duration_s / monolayer
        y0 = np.exp(inventory / monolayer)
        small = np.abs(coefficient) <= 1e-14 * np.maximum(
            deposit_substrate + deposit_polymer + removal_capacity, 1.0)
        y1 = np.empty(shape)
        y1[small] = y0[small] - transition[small] * duration_s / monolayer
        if np.any(~small):
            em1 = np.expm1(exponent[~small])
            y1[~small] = (y0[~small] * (1.0 + em1)
                          - transition[~small] / coefficient[~small] * em1)
        if np.any(y1 < 1.0 - 1e-12):
            raise RuntimeError("polymer inventory integrator violated positivity")
        y1 = np.maximum(y1, 1.0)
        updated = monolayer * np.log(y1)

        delta = updated - inventory
        nonzero_transition = np.abs(transition) > 1e-14 * np.maximum(
            deposit_substrate + deposit_polymer + removal_capacity, 1.0)
        coverage_integral = np.empty(shape)
        coverage_integral[nonzero_transition] = (
            delta[nonzero_transition]
            - deposit_substrate[nonzero_transition] * duration_s
        ) / transition[nonzero_transition]
        remaining = ~nonzero_transition
        if np.any(remaining):
            ds = deposit_substrate[remaining]; n0 = inventory[remaining]
            growing = ds > 0.0
            local = np.empty(ds.shape)
            local[~growing] = duration_s * (1.0 - np.exp(-n0[~growing] / monolayer))
            if np.any(growing):
                local[growing] = (
                    duration_s
                    - np.exp(-n0[growing] / monolayer) * monolayer / ds[growing]
                    * (-np.expm1(-ds[growing] * duration_s / monolayer)))
            coverage_integral[remaining] = local
        coverage_integral = np.clip(coverage_integral, 0.0, duration_s)
        deposited = (deposit_substrate * duration_s
                     + (deposit_polymer - deposit_substrate) * coverage_integral)
        removed = removal_capacity * coverage_integral
        polymer_balance = (updated - inventory) - (deposited - removed)
        polymer_scale = np.maximum.reduce((
            np.abs(updated - inventory), np.abs(deposited), np.abs(removed), np.ones(shape)))
        if np.any(np.abs(polymer_balance) > 5e-13 * polymer_scale):
            raise RuntimeError("polymer reaction bookkeeping failed conservation")
        return updated, deposited, removed

    def _substrate_step(self, complex_fraction, polymer_inventory, fluxes, duration_s, shape):
        par = self.parameters
        access = np.exp(-polymer_inventory / par.polymer_monolayer_density_m2)
        formation_event_rate = self._neutral_weighted_rate(
            fluxes, par.complex_formation_probability, shape) * access
        complex_removal_event_rate = self._energetic_rate(
            fluxes, par.complex_sio2_yield, shape) * access
        bare_removal_event_rate = self._energetic_rate(
            fluxes, par.bare_sio2_yield, shape) * access
        formation_hazard = formation_event_rate / par.site_density_m2
        removal_hazard = complex_removal_event_rate / par.site_density_m2
        total_hazard = formation_hazard + removal_hazard
        equilibrium = np.divide(
            formation_hazard, total_hazard, out=np.zeros(shape), where=total_hazard > 0.0)
        decay = np.exp(-total_hazard * duration_s)
        updated = equilibrium + (complex_fraction - equilibrium) * decay
        integral_complex = np.empty(shape)
        active = total_hazard > 0.0
        integral_complex[active] = (
            equilibrium[active] * duration_s
            + (complex_fraction[active] - equilibrium[active])
            * (-np.expm1(-total_hazard[active] * duration_s)) / total_hazard[active])
        integral_complex[~active] = complex_fraction[~active] * duration_s
        removed_complex = complex_removal_event_rate * integral_complex
        formed_complex = formation_event_rate * (duration_s - integral_complex)
        removed_bare = bare_removal_event_rate * (duration_s - integral_complex)
        site_change = (updated - complex_fraction) * par.site_density_m2
        site_balance = site_change - (formed_complex - removed_complex)
        site_scale = np.maximum.reduce((
            np.abs(site_change), np.abs(formed_complex), np.abs(removed_complex), np.ones(shape)))
        if np.any(np.abs(site_balance) > 5e-13 * site_scale):
            raise RuntimeError("oxide-complex reaction bookkeeping failed conservation")
        return updated, formed_complex, removed_complex, removed_bare

    def advance(self, state: SiO2SurfaceState, fluxes: SurfaceFluxes, duration_s: float, *,
                max_step_s: float | None = None, strict=True):
        """Advance constant incident fluxes using conservative Strang-split exact sub-operators."""
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        if max_step_s is not None and (not np.isfinite(max_step_s) or max_step_s <= 0.0):
            raise ValueError("max_step_s must be positive and finite")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))
        shape = state.complex_fraction.shape
        n_steps = (1 if duration_s == 0.0 or max_step_s is None else
                   max(1, int(np.ceil(duration_s / max_step_s))))
        step = duration_s / n_steps
        complex_fraction = np.array(state.complex_fraction, copy=True)
        polymer = np.array(state.polymer_units_m2, copy=True)
        removed_total = np.array(state.removed_formula_units_m2, copy=True)
        formed_complex = np.zeros(shape); removed_complex = np.zeros(shape)
        removed_bare = np.zeros(shape); deposited_polymer = np.zeros(shape)
        removed_polymer = np.zeros(shape)
        for _ in range(n_steps):
            polymer, deposited, removed = self._polymer_step(
                polymer, fluxes, 0.5 * step, shape)
            deposited_polymer += deposited; removed_polymer += removed
            complex_fraction, formed, removed_c, removed_b = self._substrate_step(
                complex_fraction, polymer, fluxes, step, shape)
            formed_complex += formed; removed_complex += removed_c; removed_bare += removed_b
            removed_total += removed_c + removed_b
            polymer, deposited, removed = self._polymer_step(
                polymer, fluxes, 0.5 * step, shape)
            deposited_polymer += deposited; removed_polymer += removed
        new_state = SiO2SurfaceState(complex_fraction, polymer, removed_total)
        velocity = ((removed_complex + removed_bare)
                    / self.parameters.bulk_formula_density_m3
                    / duration_s if duration_s > 0.0 else np.zeros(shape))
        return SurfaceStepResult(
            state=new_state,
            etch_velocity_m_s=np.asarray(velocity),
            formed_complex_units_m2=formed_complex,
            removed_complex_units_m2=removed_complex,
            removed_bare_formula_units_m2=removed_bare,
            deposited_polymer_units_m2=deposited_polymer,
            removed_polymer_units_m2=removed_polymer,
            validity=validity)
