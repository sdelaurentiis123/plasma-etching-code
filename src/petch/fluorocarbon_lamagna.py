"""La Magna--Garozzo fluorocarbon etching on the common petch engine.

This module reproduces the three-coverage SiO2 fluorocarbon model implemented by
ViennaPS 4.6.1 while retaining petch's dimensional fluxes, material ledgers, and
neutral/surface fixed-point contract.  The algebraic surface state is

    theta_pe = J_e S_pe / (J_e S_pe + J_i Y_p)
    theta_p  = J_p S_p / (J_i Y_p theta_pe + delta_p)
    theta_e  = J_e S_e (1-theta_p)
               / (k_ie J_i Y_ie + k_ev J_ev + J_e S_e)

and exposed SiO2 is removed by thermal chemical, ion-enhanced, and physical
sputtering channels.  A saturated polymer coverage grows a finite surface-film
inventory; an existing film must be removed before substrate recession.  Thus
the ViennaPS chemistry can use the same moving-surface path as charging,
reflection, and material routing without introducing a second feature solver.

Primary scientific source: A. La Magna and G. Garozzo, Journal of The
Electrochemical Society 150, F178--F185 (2003), DOI 10.1149/1.1602084.
Behavioral parity source: ViennaTools/ViennaPS, ``psFluorocarbonEtching.hpp``,
master commit 2956ed587984c6dc38be24c6e2390e10c9b2f0a7 (GPL-3.0).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange
from .surface_kinetics import (
    MechanismValidity, ParameterEvidence, SteinbruchelYield, SurfaceFluxes,
)


_SOURCE_COMMIT = "2956ed587984c6dc38be24c6e2390e10c9b2f0a7"
_SOURCE_URL = (
    "https://github.com/ViennaTools/ViennaPS/blob/"
    f"{_SOURCE_COMMIT}/include/viennaps/models/psFluorocarbonEtching.hpp"
)
_PAPER_DOI = "https://doi.org/10.1149/1.1602084"

_REQUIRED_INPUTS = frozenset({
    "bulk_formula_density_m3",
    "polymer_unit_density_m3",
    "substrate_etchant_sticking_probability",
    "substrate_polymer_sticking_probability",
    "polymer_etchant_sticking_probability",
    "polymer_polymer_sticking_probability",
    "reference_etchant_flux_m2_s",
    "polymer_loss_rate_m2_s",
    "temperature_K",
    "chemical_rate_coefficient",
    "chemical_activation_energy_eV",
    "ion_enhanced_coverage_loss_factor",
    "chemical_coverage_loss_factor",
    "physical_sputter_yield",
    "ion_enhanced_yield",
    "polymer_removal_yield",
})


def _yield_manifest(law: SteinbruchelYield):
    if not isinstance(law, SteinbruchelYield):
        raise TypeError("La Magna fluorocarbon channels require SteinbruchelYield laws")
    return {
        "type": "SteinbruchelYield",
        "prefactor_per_sqrt_eV": float(law.prefactor_per_sqrt_eV),
        "threshold_energy_eV": float(law.threshold_energy_eV),
        "angular_model": law.angular_model,
        "angular_parameter": (
            None if law.angular_parameter is None else float(law.angular_parameter)),
    }


def _evidence_manifest(evidence):
    return {
        name: {
            "source": item.source,
            "evidence_type": item.evidence_type,
            "relative_standard_uncertainty": item.relative_standard_uncertainty,
            "note": item.note,
            "supports_prediction_within_declared_domain": (
                item.supports_prediction_within_declared_domain),
        }
        for name, item in evidence.items()
    }


@dataclass(frozen=True)
class LaMagnaFluorocarbonState:
    """Warm-start coverages plus finite film and cumulative oxide inventories."""

    etchant_coverage: np.ndarray | float = 0.0
    polymer_coverage: np.ndarray | float = 0.0
    etchant_on_polymer_coverage: np.ndarray | float = 0.0
    polymer_film_units_m2: np.ndarray | float = 0.0
    removed_formula_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        values = np.broadcast_arrays(
            np.asarray(self.etchant_coverage, dtype=float),
            np.asarray(self.polymer_coverage, dtype=float),
            np.asarray(self.etchant_on_polymer_coverage, dtype=float),
            np.asarray(self.polymer_film_units_m2, dtype=float),
            np.asarray(self.removed_formula_units_m2, dtype=float),
        )
        e, p, pe, film, removed = [np.array(value, copy=True) for value in values]
        if (any(np.any(~np.isfinite(value)) for value in values)
                or np.any((e < 0.0) | (e > 1.0))
                or np.any((p < 0.0) | (p > 1.0))
                or np.any((pe < 0.0) | (pe > 1.0))
                or np.any(film < 0.0) or np.any(removed < 0.0)):
            raise ValueError("invalid La Magna fluorocarbon surface state")
        for value in (e, p, pe, film, removed):
            value.setflags(write=False)
        object.__setattr__(self, "etchant_coverage", e)
        object.__setattr__(self, "polymer_coverage", p)
        object.__setattr__(self, "etchant_on_polymer_coverage", pe)
        object.__setattr__(self, "polymer_film_units_m2", film)
        object.__setattr__(self, "removed_formula_units_m2", removed)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape)
        return cls(zero, zero, zero, zero, zero)

    def conservative_surface_fields(self):
        return {
            "etchant_coverage": self.etchant_coverage,
            "polymer_coverage": self.polymer_coverage,
            "etchant_on_polymer_coverage": self.etchant_on_polymer_coverage,
            "polymer_film_units_m2": self.polymer_film_units_m2,
            "removed_formula_units_m2": self.removed_formula_units_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {
            "etchant_coverage": 1.0,
            "polymer_coverage": 1.0,
            "etchant_on_polymer_coverage": 1.0,
            "polymer_film_units_m2": None,
            "removed_formula_units_m2": None,
        }

    def surface_field_remap_modes(self):
        """Separate algebraic intensive coverages from physical areal inventories."""
        return {
            "etchant_coverage": "intensive",
            "polymer_coverage": "intensive",
            "etchant_on_polymer_coverage": "intensive",
            "polymer_film_units_m2": "conservative",
            "removed_formula_units_m2": "conservative",
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("La Magna fluorocarbon remap fields do not match its state contract")
        return type(self)(
            fields["etchant_coverage"], fields["polymer_coverage"],
            fields["etchant_on_polymer_coverage"], fields["polymer_film_units_m2"],
            fields["removed_formula_units_m2"],
        )


@dataclass(frozen=True)
class LaMagnaFluorocarbonParameters:
    material_name: str
    material_inventory_name: str
    etchant_species: tuple[str, ...]
    polymer_species: tuple[str, ...]
    projectile_species: tuple[str, ...]
    bulk_formula_density_m3: float
    polymer_unit_density_m3: float
    substrate_etchant_sticking_probability: float
    substrate_polymer_sticking_probability: float
    polymer_etchant_sticking_probability: float
    polymer_polymer_sticking_probability: float
    reference_etchant_flux_m2_s: float
    polymer_loss_rate_m2_s: float
    temperature_K: float
    chemical_rate_coefficient: float
    chemical_activation_energy_eV: float
    ion_enhanced_coverage_loss_factor: float
    chemical_coverage_loss_factor: float
    physical_sputter_yield: SteinbruchelYield
    ion_enhanced_yield: SteinbruchelYield
    polymer_removal_yield: SteinbruchelYield
    evidence: Mapping[str, ParameterEvidence]
    neutral_transport_mode: str = "viennaps_4_6_1"
    coverage_zero_tolerance_m2_s: float = 0.0
    known_omissions: tuple[str, ...] = (
        "reactor-to-feature fluxes must be supplied externally",
        "volatile SiO2 fluorocarbon product identities and return transport are unresolved",
        "the polymer film is an areal surface state rather than a separately labeled level set",
        "surface charging is supplied by the common charging engine rather than this chemistry law",
    )

    def __post_init__(self):
        etchant = tuple(self.etchant_species)
        polymer = tuple(self.polymer_species)
        projectiles = tuple(self.projectile_species)
        probabilities = (
            self.substrate_etchant_sticking_probability,
            self.substrate_polymer_sticking_probability,
            self.polymer_etchant_sticking_probability,
            self.polymer_polymer_sticking_probability,
        )
        positive = (
            self.bulk_formula_density_m3, self.polymer_unit_density_m3,
            self.reference_etchant_flux_m2_s, self.temperature_K,
            self.chemical_rate_coefficient, self.ion_enhanced_coverage_loss_factor,
            self.chemical_coverage_loss_factor,
        )
        if (not self.material_name or not self.material_inventory_name
                or not etchant or not polymer or not projectiles
                or any(not name for name in etchant + polymer + projectiles)
                or len(set(etchant + polymer)) != len(etchant + polymer)
                or set(etchant + polymer) & set(projectiles)
                or any(not np.isfinite(value) or value <= 0.0 for value in positive)
                or any(not np.isfinite(value) or not 0.0 <= value <= 1.0
                       for value in probabilities)
                or not np.isfinite(self.polymer_loss_rate_m2_s)
                or self.polymer_loss_rate_m2_s < 0.0
                or not np.isfinite(self.chemical_activation_energy_eV)
                or self.chemical_activation_energy_eV < 0.0
                or not np.isfinite(self.coverage_zero_tolerance_m2_s)
                or self.coverage_zero_tolerance_m2_s < 0.0
                or self.neutral_transport_mode not in {
                    "viennaps_4_6_1", "species_specific"}):
            raise ValueError("invalid La Magna fluorocarbon parameters")
        for law in (
                self.physical_sputter_yield, self.ion_enhanced_yield,
                self.polymer_removal_yield):
            _yield_manifest(law)
        evidence = dict(self.evidence)
        if (set(evidence) != _REQUIRED_INPUTS
                or any(not isinstance(item, ParameterEvidence) for item in evidence.values())):
            raise ValueError("La Magna evidence must cover every physical input")
        object.__setattr__(self, "etchant_species", etchant)
        object.__setattr__(self, "polymer_species", polymer)
        object.__setattr__(self, "projectile_species", projectiles)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))

    @classmethod
    def viennaps_4_6_1_reference(
            cls, *, reference_etchant_flux_m2_s,
            etchant_species=("FC_etchant",), polymer_species=("FC_polymer",),
            projectile_species=("Ar+",), neutral_transport_mode="viennaps_4_6_1"):
        """Documented ViennaPS reference parameters in SI units.

        The polymer material overrides (density 2.0, beta_e 0.6, doubled A_ie)
        are the ones used by the official ``stackEtching`` example.  They remain
        transferred defaults, not Jeong-calibrated values.
        """
        evidence = {
            name: ParameterEvidence(
                _PAPER_DOI if name not in {
                    "polymer_unit_density_m3", "polymer_etchant_sticking_probability",
                    "polymer_polymer_sticking_probability", "polymer_removal_yield",
                } else _SOURCE_URL,
                "published_model_parameter" if name not in {
                    "polymer_unit_density_m3", "polymer_etchant_sticking_probability",
                    "polymer_polymer_sticking_probability", "polymer_removal_yield",
                } else "ViennaPS_4.6.1_example_parameter",
                note="Transferred reference value; not calibrated to the target experiment.",
                supports_prediction_within_declared_domain=False,
            )
            for name in _REQUIRED_INPUTS
        }
        return cls(
            material_name="SiO2", material_inventory_name="SiO2_formula_unit",
            etchant_species=tuple(etchant_species), polymer_species=tuple(polymer_species),
            projectile_species=tuple(projectile_species),
            bulk_formula_density_m3=2.2e28, polymer_unit_density_m3=2.0e28,
            substrate_etchant_sticking_probability=0.9,
            substrate_polymer_sticking_probability=0.26,
            polymer_etchant_sticking_probability=0.6,
            polymer_polymer_sticking_probability=0.26,
            reference_etchant_flux_m2_s=float(reference_etchant_flux_m2_s),
            polymer_loss_rate_m2_s=1.0e19, temperature_K=300.0,
            chemical_rate_coefficient=0.002789491704544977,
            chemical_activation_energy_eV=0.168,
            ion_enhanced_coverage_loss_factor=2.0,
            chemical_coverage_loss_factor=2.0,
            physical_sputter_yield=SteinbruchelYield(
                0.0139, 18.0, angular_model="kress_1999", angular_parameter=9.3),
            ion_enhanced_yield=SteinbruchelYield(
                0.0361, 4.0, angular_model="kress_1999", angular_parameter=0.0),
            polymer_removal_yield=SteinbruchelYield(
                0.0722, 4.0, angular_model="kress_1999", angular_parameter=0.0),
            evidence=evidence, neutral_transport_mode=neutral_transport_mode,
        )


@dataclass(frozen=True)
class LaMagnaFluorocarbonStepResult:
    state: LaMagnaFluorocarbonState
    etch_velocity_m_s: np.ndarray
    normal_growth_velocity_m_s: np.ndarray
    etchant_coverage: np.ndarray
    polymer_coverage: np.ndarray
    etchant_on_polymer_coverage: np.ndarray
    chemical_removal_rate_m2_s: np.ndarray
    ion_enhanced_removal_rate_m2_s: np.ndarray
    physical_sputter_rate_m2_s: np.ndarray
    polymer_deposition_rate_m2_s: np.ndarray
    polymer_removal_rate_m2_s: np.ndarray
    removed_formula_units_m2: np.ndarray
    deposited_polymer_units_m2: np.ndarray
    removed_polymer_units_m2: np.ndarray
    etchant_site_balance_residual_m2_s: np.ndarray
    polymer_site_balance_residual_m2_s: np.ndarray
    etchant_on_polymer_balance_residual_m2_s: np.ndarray
    transport_fixed_point_change: np.ndarray
    material_exchange: SurfaceMaterialExchange
    validity: MechanismValidity
    product_populations: tuple = ()

    def __post_init__(self):
        if not isinstance(self.state, LaMagnaFluorocarbonState):
            raise TypeError("invalid La Magna fluorocarbon result state")
        nonnegative = {
            "etch_velocity_m_s", "normal_growth_velocity_m_s", "etchant_coverage",
            "polymer_coverage", "etchant_on_polymer_coverage",
            "chemical_removal_rate_m2_s", "ion_enhanced_removal_rate_m2_s",
            "physical_sputter_rate_m2_s", "polymer_deposition_rate_m2_s",
            "polymer_removal_rate_m2_s", "removed_formula_units_m2",
            "deposited_polymer_units_m2", "removed_polymer_units_m2",
        }
        for name in (
                *nonnegative, "etchant_site_balance_residual_m2_s",
                "polymer_site_balance_residual_m2_s",
                "etchant_on_polymer_balance_residual_m2_s",
                "transport_fixed_point_change"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or (name in nonnegative and np.any(value < 0.0)):
                raise ValueError(f"invalid La Magna result field: {name}")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if not isinstance(self.material_exchange, SurfaceMaterialExchange):
            raise TypeError("La Magna result requires a material-exchange ledger")
        if not isinstance(self.validity, MechanismValidity):
            raise TypeError("La Magna result requires mechanism validity")
        object.__setattr__(self, "product_populations", tuple(self.product_populations))
        if self.product_populations:
            raise ValueError("La Magna v1 leaves volatile products explicitly unresolved")


class LaMagnaGarozzoFluorocarbonMechanism:
    """ViennaPS-compatible quasi-steady fluorocarbon law on the common engine."""

    quasi_steady_surface_state = True

    def __init__(self, parameters: LaMagnaFluorocarbonParameters):
        if not isinstance(parameters, LaMagnaFluorocarbonParameters):
            raise TypeError("parameters must be LaMagnaFluorocarbonParameters")
        self.parameters = parameters
        par = parameters
        self.provenance = MappingProxyType({
            "model": "lamagna-garozzo-fluorocarbon-common-engine-v1",
            "scientific_source": _PAPER_DOI,
            "behavioral_parity_source": _SOURCE_URL,
            "viennaps_source_commit": _SOURCE_COMMIT,
            "neutral_transport_mode": par.neutral_transport_mode,
            "parameters": {
                "material_name": par.material_name,
                "material_inventory_name": par.material_inventory_name,
                "etchant_species": list(par.etchant_species),
                "polymer_species": list(par.polymer_species),
                "projectile_species": list(par.projectile_species),
                "bulk_formula_density_m3": par.bulk_formula_density_m3,
                "polymer_unit_density_m3": par.polymer_unit_density_m3,
                "substrate_etchant_sticking_probability": (
                    par.substrate_etchant_sticking_probability),
                "substrate_polymer_sticking_probability": (
                    par.substrate_polymer_sticking_probability),
                "polymer_etchant_sticking_probability": (
                    par.polymer_etchant_sticking_probability),
                "polymer_polymer_sticking_probability": (
                    par.polymer_polymer_sticking_probability),
                "reference_etchant_flux_m2_s": par.reference_etchant_flux_m2_s,
                "polymer_loss_rate_m2_s": par.polymer_loss_rate_m2_s,
                "temperature_K": par.temperature_K,
                "chemical_rate_coefficient": par.chemical_rate_coefficient,
                "chemical_activation_energy_eV": par.chemical_activation_energy_eV,
                "ion_enhanced_coverage_loss_factor": (
                    par.ion_enhanced_coverage_loss_factor),
                "chemical_coverage_loss_factor": par.chemical_coverage_loss_factor,
                "physical_sputter_yield": _yield_manifest(par.physical_sputter_yield),
                "ion_enhanced_yield": _yield_manifest(par.ion_enhanced_yield),
                "polymer_removal_yield": _yield_manifest(par.polymer_removal_yield),
            },
            "sources": _evidence_manifest(par.evidence),
            "known_omissions": list(par.known_omissions),
        })

    @staticmethod
    def initial_state(shape=()):
        return LaMagnaFluorocarbonState.bare(shape)

    @staticmethod
    def _broadcast(value, shape, label):
        try:
            return np.broadcast_to(np.asarray(value, dtype=float), shape)
        except ValueError as error:
            raise ValueError(f"{label} does not match the fluorocarbon state shape") from error

    def _sum_neutral(self, fluxes, names, shape):
        total = np.zeros(shape)
        for name in names:
            total = total + self._broadcast(
                fluxes.neutral_flux_m2_s.get(name, 0.0), shape, name)
        return total

    def _energetic_rate(self, fluxes, law, shape):
        total = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name in self.parameters.projectile_species:
                total = total + self._broadcast(
                    population.yield_rate_m2_s(law), shape, population.name)
        return total

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        supported_neutral = set(par.etchant_species) | set(par.polymer_species)
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in supported_neutral and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in par.projectile_species
            and np.any(np.asarray(population.flux_m2_s) > 0.0)}))
        reasons = []
        if unsupported_neutral:
            reasons.append("positive neutral flux has no declared fluorocarbon reaction channel")
        if unsupported_energetic:
            reasons.append(
                "positive energetic flux has no declared fluorocarbon channel: "
                + ", ".join(unsupported_energetic))
        nonpredictive = tuple(sorted(
            name for name, evidence in par.evidence.items()
            if not evidence.supports_prediction_within_declared_domain))
        return MechanismValidity(
            within_declared_scope=not reasons, reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=par.known_omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive,
        )

    def neutral_reaction_probability(self, state: LaMagnaFluorocarbonState):
        if not isinstance(state, LaMagnaFluorocarbonState):
            raise TypeError("La Magna neutral probabilities require LaMagnaFluorocarbonState")
        par = self.parameters
        on_polymer = state.polymer_film_units_m2 > 0.0
        beta_e = np.where(
            on_polymer, par.polymer_etchant_sticking_probability,
            par.substrate_etchant_sticking_probability)
        if par.neutral_transport_mode == "viennaps_4_6_1":
            beta_p = beta_e
        else:
            beta_p = np.where(
                on_polymer, par.polymer_polymer_sticking_probability,
                par.substrate_polymer_sticking_probability)
        available = np.maximum(
            1.0 - state.etchant_coverage - state.polymer_coverage, 0.0)
        result = {
            name: np.minimum(beta_e * available, 1.0)
            for name in par.etchant_species}
        result.update({
            name: np.minimum(beta_p * available, 1.0)
            for name in par.polymer_species})
        return MappingProxyType(result)

    def _coverage_solution(self, fluxes, shape, *, on_polymer):
        par = self.parameters
        etchant = self._sum_neutral(fluxes, par.etchant_species, shape)
        polymer = self._sum_neutral(fluxes, par.polymer_species, shape)
        ion_sp = self._energetic_rate(fluxes, par.physical_sputter_yield, shape)
        ion_ie = self._energetic_rate(fluxes, par.ion_enhanced_yield, shape)
        ion_polymer = self._energetic_rate(fluxes, par.polymer_removal_yield, shape)
        beta_e = (par.polymer_etchant_sticking_probability if on_polymer
                  else par.substrate_etchant_sticking_probability)
        beta_p = (par.polymer_polymer_sticking_probability if on_polymer
                  else par.substrate_polymer_sticking_probability)
        zero = par.coverage_zero_tolerance_m2_s

        pe_denominator = etchant * par.polymer_etchant_sticking_probability + ion_polymer
        pe = np.divide(
            etchant * par.polymer_etchant_sticking_probability, pe_denominator,
            out=np.zeros(shape), where=pe_denominator > zero)

        p_denominator = ion_polymer * pe + par.polymer_loss_rate_m2_s
        raw_p = np.divide(
            polymer * beta_p, p_denominator,
            out=np.zeros(shape), where=p_denominator > zero)
        forced_saturation = (
            (polymer > zero) & ((pe <= zero) | (ion_polymer <= zero)))
        saturated = forced_saturation | (raw_p >= 1.0)
        p = np.where(saturated, 1.0, raw_p)

        boltzmann_eV_K = 8.617333262145e-5
        chemical_reference = (
            par.chemical_rate_coefficient * par.reference_etchant_flux_m2_s
            * np.exp(-par.chemical_activation_energy_eV
                     / (boltzmann_eV_K * par.temperature_K)))
        e_supply = etchant * beta_e
        e_denominator = (
            par.ion_enhanced_coverage_loss_factor * ion_ie
            + par.chemical_coverage_loss_factor * chemical_reference + e_supply)
        e = np.divide(
            e_supply * (1.0 - p), e_denominator,
            out=np.zeros(shape), where=e_denominator > zero)

        chemical = chemical_reference * e
        enhanced = ion_ie * e
        sputter = ion_sp * (1.0 - e)
        polymer_net = (
            polymer * par.polymer_polymer_sticking_probability - ion_polymer * pe)
        return {
            "etchant": etchant, "polymer": polymer, "ion_sp": ion_sp,
            "ion_ie": ion_ie, "ion_polymer": ion_polymer,
            "e": e, "p": p, "pe": pe, "saturated": saturated,
            "chemical_reference": np.broadcast_to(chemical_reference, shape),
            "chemical": chemical, "enhanced": enhanced, "sputter": sputter,
            "substrate_removal": chemical + enhanced + sputter,
            "polymer_net": polymer_net, "beta_e": beta_e, "beta_p": beta_p,
        }

    @staticmethod
    def _select(mask, left, right):
        return np.where(mask, left, right)

    def advance(
            self, state: LaMagnaFluorocarbonState, fluxes: SurfaceFluxes,
            duration_s: float, *, strict=True):
        if not isinstance(state, LaMagnaFluorocarbonState):
            raise TypeError("La Magna advance requires LaMagnaFluorocarbonState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))

        par = self.parameters
        shape = state.etchant_coverage.shape
        substrate = self._coverage_solution(fluxes, shape, on_polymer=False)
        film = self._coverage_solution(fluxes, shape, on_polymer=True)
        film_present = state.polymer_film_units_m2 > 0.0
        selected = {
            name: self._select(film_present, film[name], substrate[name])
            for name in substrate}

        dt = float(duration_s)
        growth_rate = np.where(
            selected["saturated"], np.maximum(selected["polymer_net"], 0.0), 0.0)
        film_removal_rate = np.where(
            film_present, np.maximum(-selected["polymer_net"], 0.0), 0.0)
        deposited_polymer = growth_rate * dt
        removed_polymer = np.minimum(
            state.polymer_film_units_m2, film_removal_rate * dt)
        final_film = (
            state.polymer_film_units_m2 + deposited_polymer - removed_polymer)

        depletion_time = np.divide(
            state.polymer_film_units_m2, film_removal_rate,
            out=np.full(shape, np.inf), where=film_removal_rate > 0.0)
        substrate_time = np.where(
            ~film_present & ~substrate["saturated"], dt, 0.0)
        substrate_time = np.where(
            film_present & (depletion_time < dt) & ~substrate["saturated"],
            dt - depletion_time, substrate_time)
        removed_formula = substrate["substrate_removal"] * substrate_time

        if dt > 0.0:
            etch_velocity = (
                removed_formula / par.bulk_formula_density_m3
                + removed_polymer / par.polymer_unit_density_m3) / dt
            growth_velocity = deposited_polymer / par.polymer_unit_density_m3 / dt
        else:
            etch_velocity = np.where(
                film_present,
                np.maximum(-selected["polymer_net"], 0.0) / par.polymer_unit_density_m3,
                np.where(substrate["saturated"], 0.0,
                         substrate["substrate_removal"] / par.bulk_formula_density_m3))
            growth_velocity = np.where(
                selected["saturated"], np.maximum(selected["polymer_net"], 0.0)
                / par.polymer_unit_density_m3, 0.0)

        end_on_polymer = final_film > 0.0
        end = {
            name: self._select(end_on_polymer, film[name], substrate[name])
            for name in substrate}
        updated = LaMagnaFluorocarbonState(
            end["e"], end["p"], end["pe"], final_film,
            state.removed_formula_units_m2 + removed_formula)
        fixed_point_change = np.maximum.reduce((
            np.abs(updated.etchant_coverage - state.etchant_coverage),
            np.abs(updated.polymer_coverage - state.polymer_coverage),
            np.abs(updated.etchant_on_polymer_coverage
                   - state.etchant_on_polymer_coverage),
        ))

        etchant_balance = (
            end["etchant"] * end["beta_e"]
            * (1.0 - end["p"] - end["e"])
            - (par.ion_enhanced_coverage_loss_factor * end["ion_ie"]
               + par.chemical_coverage_loss_factor * end["chemical_reference"])
            * end["e"])
        polymer_balance = (
            end["polymer"] * end["beta_p"]
            - (end["ion_polymer"] * end["pe"] + par.polymer_loss_rate_m2_s)
            * end["p"])
        pe_balance = (
            end["etchant"] * par.polymer_etchant_sticking_probability
            * (1.0 - end["pe"]) - end["ion_polymer"] * end["pe"])

        exchange = unresolved_surface_exchange(
            removed_units_m2={
                par.material_inventory_name: removed_formula,
                "fluorocarbon_film_unit": removed_polymer,
            },
            deposited_units_m2={"fluorocarbon_film_unit": deposited_polymer},
            limitations=(
                "volatile SiO2/fluorocarbon product branching and return transport are unresolved",
                "incident etchant/polymer site balance is reported separately from material inventory",
            ))
        return LaMagnaFluorocarbonStepResult(
            state=updated, etch_velocity_m_s=etch_velocity,
            normal_growth_velocity_m_s=growth_velocity,
            etchant_coverage=end["e"], polymer_coverage=end["p"],
            etchant_on_polymer_coverage=end["pe"],
            chemical_removal_rate_m2_s=end["chemical"],
            ion_enhanced_removal_rate_m2_s=end["enhanced"],
            physical_sputter_rate_m2_s=end["sputter"],
            polymer_deposition_rate_m2_s=growth_rate,
            polymer_removal_rate_m2_s=film_removal_rate,
            removed_formula_units_m2=removed_formula,
            deposited_polymer_units_m2=deposited_polymer,
            removed_polymer_units_m2=removed_polymer,
            etchant_site_balance_residual_m2_s=etchant_balance,
            polymer_site_balance_residual_m2_s=polymer_balance,
            etchant_on_polymer_balance_residual_m2_s=pe_balance,
            transport_fixed_point_change=fixed_point_change,
            material_exchange=exchange, validity=validity,
        )
