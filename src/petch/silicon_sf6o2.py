"""Provenance-bearing Belen SF6/O2 silicon kinetics for the common feature engine.

This is the common-engine form of the surface law used by Belen et al. and ViennaPS.  It
contains the channel absent from the earlier de-Boer adapter: fluorinated silicon can leave by
the direct chemical term even without an ion impact.  The quasi-steady site balances are

    0 = gamma_F Gamma_F (1-theta_F-theta_O)
        - (k_sigma + nu_ie Y_ie Gamma_i) theta_F
    0 = gamma_O Gamma_O (1-theta_F-theta_O)
        - (beta_sigma + Y_O Gamma_i) theta_O

and the silicon removal rate is

    R_Si = k_sigma theta_F / nu_F + Y_sp Gamma_i + theta_F Y_ie Gamma_i.

The neutral radiosity sink uses the same ``gamma * vacant-site fraction`` as the site balance.
Consequently transport and chemistry form a small physical fixed point.  ``advance(..., 0)``
updates that algebraic surface state without moving material, so an orchestrator can converge the
same operator before applying a profile step.  No photon channel is present: the declared SF6/O2
model and the de Boer evidence do not provide a photon flux or photon-assisted yield.

Primary model sources: R. J. Belen et al., JVST A 23, 99 (2005),
DOI 10.1116/1.1830495; R. J. Belen et al., JVST A 23, 1430 (2005),
DOI 10.1116/1.2013317.  The square-root energetic law follows Steinbruechel,
Appl. Phys. Lett. 55, 1960 (1989).
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange
from .surface_kinetics import (
    EnergeticYield, MechanismValidity, ParameterEvidence, SteinbruchelYield,
    SurfaceFluxes,
)


_REQUIRED_INPUTS = frozenset({
    "site_density_m2",
    "bulk_si_atom_density_m3",
    "fluorine_sticking_probability",
    "oxygen_sticking_probability",
    "spontaneous_fluorine_removal_rate_m2_s",
    "oxygen_desorption_rate_m2_s",
    "physical_sputter_yield",
    "ion_enhanced_yield",
    "oxygen_sputter_yield",
    "fluorine_atoms_per_removed_si",
    "ion_enhanced_fluorine_release_per_si",
})


def _freeze_bound(value, path):
    if isinstance(value, Mapping):
        supplied = dict(value)
        if not supplied:
            raise ValueError(f"{path} must not be empty")
        return MappingProxyType({
            str(name): _freeze_bound(item, f"{path}.{name}")
            for name, item in supplied.items()})
    try:
        lower, upper = value
    except (TypeError, ValueError) as error:
        raise ValueError(f"{path} must be a numeric [lower, upper] pair") from error
    lower = float(lower); upper = float(upper)
    if not np.isfinite(lower) or not np.isfinite(upper) or upper < lower:
        raise ValueError(f"{path} has invalid bounds")
    return (lower, upper)


def _yield_manifest(law):
    if isinstance(law, SteinbruchelYield):
        return {
            "type": "SteinbruchelYield",
            "prefactor_per_sqrt_eV": float(law.prefactor_per_sqrt_eV),
            "threshold_energy_eV": float(law.threshold_energy_eV),
            "angular_model": law.angular_model,
            "angular_parameter": (
                None if law.angular_parameter is None else float(law.angular_parameter)),
        }
    if isinstance(law, EnergeticYield):
        return {
            "type": "EnergeticYield",
            "reference_yield": float(law.reference_yield),
            "threshold_energy_eV": float(law.threshold_energy_eV),
            "reference_energy_eV": float(law.reference_energy_eV),
            "energy_exponent": float(law.energy_exponent),
            "angular_model": law.angular_model,
            "angular_parameter": (
                None if law.angular_parameter is None else float(law.angular_parameter)),
        }
    raise TypeError("Belen silicon channels require a declared energetic-yield law")


def _yield_numeric_parameters(law):
    manifest = _yield_manifest(law)
    return {
        name: value for name, value in manifest.items()
        if name not in {"type", "angular_model"} and value is not None}


def _require_inside_bounds(actual, bounds, path):
    if isinstance(bounds, Mapping):
        actual = dict(actual)
        if set(actual) != set(bounds):
            raise ValueError(f"{path} bounds do not cover every numeric parameter")
        for name, value in actual.items():
            _require_inside_bounds(value, bounds[name], f"{path}.{name}")
        return
    lower, upper = bounds
    value = float(actual)
    if value < lower or value > upper:
        raise ValueError(f"{path}={value:g} lies outside declared bounds [{lower:g}, {upper:g}]")


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
class BelenSiliconState:
    """Warm-start state plus authoritative cumulative removed-Si inventory."""

    available_site_fraction: np.ndarray | float = 1.0
    removed_si_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        available, removed = np.broadcast_arrays(
            np.asarray(self.available_site_fraction, dtype=float),
            np.asarray(self.removed_si_atoms_m2, dtype=float))
        available = np.array(available, copy=True); removed = np.array(removed, copy=True)
        if (np.any(~np.isfinite(available))
                or np.any((available < 0.0) | (available > 1.0))
                or np.any(~np.isfinite(removed)) or np.any(removed < 0.0)):
            raise ValueError("invalid Belen silicon surface state")
        available.setflags(write=False); removed.setflags(write=False)
        object.__setattr__(self, "available_site_fraction", available)
        object.__setattr__(self, "removed_si_atoms_m2", removed)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.ones(shape), np.zeros(shape))

    def conservative_surface_fields(self):
        return {
            "available_site_fraction": self.available_site_fraction,
            "removed_si_atoms_m2": self.removed_si_atoms_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {"available_site_fraction": 1.0, "removed_si_atoms_m2": None}

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != {"available_site_fraction", "removed_si_atoms_m2"}:
            raise ValueError("Belen silicon remap fields do not match its state contract")
        return type(self)(fields["available_site_fraction"], fields["removed_si_atoms_m2"])


@dataclass(frozen=True)
class BelenSiliconParameters:
    material_name: str
    material_inventory_name: str
    fluorine_species: str
    oxygen_species: str
    projectile_species: tuple[str, ...]
    site_density_m2: float
    bulk_si_atom_density_m3: float
    fluorine_sticking_probability: float
    oxygen_sticking_probability: float
    spontaneous_fluorine_removal_rate_m2_s: float
    oxygen_desorption_rate_m2_s: float
    physical_sputter_yield: EnergeticYield | SteinbruchelYield
    ion_enhanced_yield: EnergeticYield | SteinbruchelYield
    oxygen_sputter_yield: EnergeticYield | SteinbruchelYield
    fluorine_atoms_per_removed_si: float
    ion_enhanced_fluorine_release_per_si: float
    evidence: Mapping[str, ParameterEvidence]
    parameter_bounds: Mapping[str, object]
    known_omissions: tuple[str, ...] = (
        "crystal-orientation-dependent silicon etching is not represented",
        "SiOxFy composition and temperature-dependent passivation are reduced to oxygen coverage",
        "volatile SiFx product identity and return transport are unresolved",
        "surface roughness and micro-masking are not represented",
        "no photon-assisted reaction channel is declared",
    )

    def __post_init__(self):
        projectiles = tuple(self.projectile_species)
        yield_types = (EnergeticYield, SteinbruchelYield)
        if (not self.material_name or not self.material_inventory_name
                or not self.fluorine_species or not self.oxygen_species
                or self.fluorine_species == self.oxygen_species
                or not projectiles or any(not name for name in projectiles)
                or len(set(projectiles)) != len(projectiles)
                or not np.isfinite(self.site_density_m2) or self.site_density_m2 <= 0.0
                or not np.isfinite(self.bulk_si_atom_density_m3)
                or self.bulk_si_atom_density_m3 <= 0.0
                or not np.isfinite(self.fluorine_sticking_probability)
                or not 0.0 <= self.fluorine_sticking_probability <= 1.0
                or not np.isfinite(self.oxygen_sticking_probability)
                or not 0.0 <= self.oxygen_sticking_probability <= 1.0
                or not np.isfinite(self.spontaneous_fluorine_removal_rate_m2_s)
                or self.spontaneous_fluorine_removal_rate_m2_s <= 0.0
                or not np.isfinite(self.oxygen_desorption_rate_m2_s)
                or self.oxygen_desorption_rate_m2_s <= 0.0
                or not isinstance(self.physical_sputter_yield, yield_types)
                or not isinstance(self.ion_enhanced_yield, yield_types)
                or not isinstance(self.oxygen_sputter_yield, yield_types)
                or not np.isfinite(self.fluorine_atoms_per_removed_si)
                or self.fluorine_atoms_per_removed_si <= 0.0
                or not np.isfinite(self.ion_enhanced_fluorine_release_per_si)
                or self.ion_enhanced_fluorine_release_per_si < 0.0):
            raise ValueError("invalid Belen silicon parameters")
        evidence = dict(self.evidence)
        bounds = dict(self.parameter_bounds)
        if set(evidence) != _REQUIRED_INPUTS or any(
                not isinstance(item, ParameterEvidence) for item in evidence.values()):
            raise ValueError("Belen silicon evidence must cover every physical input")
        if set(bounds) != _REQUIRED_INPUTS:
            raise ValueError("Belen silicon bounds must cover every physical input")
        bounds = {
            name: _freeze_bound(value, f"parameter_bounds.{name}")
            for name, value in bounds.items()}
        actual = {
            "site_density_m2": self.site_density_m2,
            "bulk_si_atom_density_m3": self.bulk_si_atom_density_m3,
            "fluorine_sticking_probability": self.fluorine_sticking_probability,
            "oxygen_sticking_probability": self.oxygen_sticking_probability,
            "spontaneous_fluorine_removal_rate_m2_s": (
                self.spontaneous_fluorine_removal_rate_m2_s),
            "oxygen_desorption_rate_m2_s": self.oxygen_desorption_rate_m2_s,
            "physical_sputter_yield": _yield_numeric_parameters(self.physical_sputter_yield),
            "ion_enhanced_yield": _yield_numeric_parameters(self.ion_enhanced_yield),
            "oxygen_sputter_yield": _yield_numeric_parameters(self.oxygen_sputter_yield),
            "fluorine_atoms_per_removed_si": self.fluorine_atoms_per_removed_si,
            "ion_enhanced_fluorine_release_per_si": (
                self.ion_enhanced_fluorine_release_per_si),
        }
        for name in _REQUIRED_INPUTS:
            _require_inside_bounds(actual[name], bounds[name], f"parameter_bounds.{name}")
        object.__setattr__(self, "projectile_species", projectiles)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "parameter_bounds", MappingProxyType(bounds))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))


@dataclass(frozen=True)
class BelenSiliconStepResult:
    state: BelenSiliconState
    etch_velocity_m_s: np.ndarray
    fluorine_coverage: np.ndarray
    oxygen_coverage: np.ndarray
    available_site_fraction: np.ndarray
    chemical_removal_rate_m2_s: np.ndarray
    physical_sputter_rate_m2_s: np.ndarray
    ion_enhanced_removal_rate_m2_s: np.ndarray
    removed_si_atoms_m2: np.ndarray
    fluorine_site_balance_residual_m2_s: np.ndarray
    oxygen_site_balance_residual_m2_s: np.ndarray
    transport_fixed_point_change: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple = ()
    validity: MechanismValidity | None = None

    def __post_init__(self):
        if not isinstance(self.state, BelenSiliconState):
            raise TypeError("invalid Belen silicon step state")
        for name in (
                "etch_velocity_m_s", "fluorine_coverage", "oxygen_coverage",
                "available_site_fraction", "chemical_removal_rate_m2_s",
                "physical_sputter_rate_m2_s", "ion_enhanced_removal_rate_m2_s",
                "removed_si_atoms_m2", "fluorine_site_balance_residual_m2_s",
                "oxygen_site_balance_residual_m2_s", "transport_fixed_point_change"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)):
                raise ValueError(f"non-finite Belen silicon result field: {name}")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        object.__setattr__(self, "product_populations", tuple(self.product_populations))
        if self.product_populations:
            raise ValueError("Belen silicon v1 leaves volatile products explicitly unresolved")
        if not isinstance(self.material_exchange, SurfaceMaterialExchange):
            raise TypeError("Belen silicon result requires a material-exchange ledger")
        if not isinstance(self.validity, MechanismValidity):
            raise TypeError("Belen silicon result requires mechanism validity")


class BelenSiliconSF6O2Mechanism:
    """Quasi-steady Belen SF6/O2 law on the common dimensional engine."""

    quasi_steady_surface_state = True

    def __init__(self, parameters: BelenSiliconParameters):
        if not isinstance(parameters, BelenSiliconParameters):
            raise TypeError("parameters must be BelenSiliconParameters")
        self.parameters = parameters
        par = parameters
        self.provenance = MappingProxyType({
            "model": "belen-silicon-sf6o2-common-engine-v1",
            "equations": "Belen 2005 pseudo-steady coupled F/O coverage and Si removal",
            "primary_sources": [
                "https://doi.org/10.1116/1.1830495",
                "https://doi.org/10.1116/1.2013317",
                "https://doi.org/10.1063/1.102336",
            ],
            "parameters": {
                "material_name": par.material_name,
                "material_inventory_name": par.material_inventory_name,
                "fluorine_species": par.fluorine_species,
                "oxygen_species": par.oxygen_species,
                "projectile_species": list(par.projectile_species),
                "site_density_m2": par.site_density_m2,
                "bulk_si_atom_density_m3": par.bulk_si_atom_density_m3,
                "fluorine_sticking_probability": par.fluorine_sticking_probability,
                "oxygen_sticking_probability": par.oxygen_sticking_probability,
                "spontaneous_fluorine_removal_rate_m2_s": (
                    par.spontaneous_fluorine_removal_rate_m2_s),
                "oxygen_desorption_rate_m2_s": par.oxygen_desorption_rate_m2_s,
                "physical_sputter_yield": _yield_manifest(par.physical_sputter_yield),
                "ion_enhanced_yield": _yield_manifest(par.ion_enhanced_yield),
                "oxygen_sputter_yield": _yield_manifest(par.oxygen_sputter_yield),
                "fluorine_atoms_per_removed_si": par.fluorine_atoms_per_removed_si,
                "ion_enhanced_fluorine_release_per_si": (
                    par.ion_enhanced_fluorine_release_per_si),
            },
            "sources": _evidence_manifest(par.evidence),
            "bounds": par.parameter_bounds,
            "known_omissions": list(par.known_omissions),
        })

    @staticmethod
    def initial_state(shape=()):
        return BelenSiliconState.bare(shape)

    def neutral_reaction_probability(self, state: BelenSiliconState):
        if not isinstance(state, BelenSiliconState):
            raise TypeError("Belen silicon neutral probabilities require BelenSiliconState")
        available = state.available_site_fraction
        return MappingProxyType({
            self.parameters.fluorine_species:
                self.parameters.fluorine_sticking_probability * available,
            self.parameters.oxygen_species:
                self.parameters.oxygen_sticking_probability * available,
        })

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        supported_neutral = {par.fluorine_species, par.oxygen_species}
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in supported_neutral and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in par.projectile_species
            and np.any(np.asarray(population.flux_m2_s) > 0.0)}))
        reasons = []
        if unsupported_neutral:
            reasons.append("positive incident neutral flux has no declared silicon reaction channel")
        if unsupported_energetic:
            reasons.append(
                "positive energetic flux has no declared silicon channel: "
                + ", ".join(unsupported_energetic))
        nonpredictive = tuple(sorted(
            name for name, evidence in par.evidence.items()
            if not evidence.supports_prediction_within_declared_domain))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=par.known_omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive)

    @staticmethod
    def _broadcast(value, shape, label):
        try:
            return np.broadcast_to(np.asarray(value, dtype=float), shape)
        except ValueError as error:
            raise ValueError(f"{label} does not match the Belen silicon state shape") from error

    def _energetic_rate(self, fluxes, law, shape):
        total = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name in self.parameters.projectile_species:
                total = total + self._broadcast(
                    population.yield_rate_m2_s(law), shape, population.name)
        return total

    def advance(
            self, state: BelenSiliconState, fluxes: SurfaceFluxes,
            duration_s: float, *, strict=True):
        if not isinstance(state, BelenSiliconState):
            raise TypeError("Belen silicon advance requires BelenSiliconState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError("surface mechanism outside declared scope: " + "; ".join(validity.reasons))

        par = self.parameters
        shape = state.available_site_fraction.shape
        fluorine_incident = self._broadcast(
            fluxes.neutral_flux_m2_s.get(par.fluorine_species, 0.0),
            shape, par.fluorine_species)
        oxygen_incident = self._broadcast(
            fluxes.neutral_flux_m2_s.get(par.oxygen_species, 0.0),
            shape, par.oxygen_species)
        fluorine_supply = par.fluorine_sticking_probability * fluorine_incident
        oxygen_supply = par.oxygen_sticking_probability * oxygen_incident
        physical = self._energetic_rate(fluxes, par.physical_sputter_yield, shape)
        ion_enhanced = self._energetic_rate(fluxes, par.ion_enhanced_yield, shape)
        oxygen_sputter = self._energetic_rate(fluxes, par.oxygen_sputter_yield, shape)

        fluorine_loss = (
            par.spontaneous_fluorine_removal_rate_m2_s
            + par.ion_enhanced_fluorine_release_per_si * ion_enhanced)
        oxygen_loss = par.oxygen_desorption_rate_m2_s + oxygen_sputter
        fluorine_ratio = fluorine_supply / fluorine_loss
        oxygen_ratio = oxygen_supply / oxygen_loss
        denominator = 1.0 + fluorine_ratio + oxygen_ratio
        available = 1.0 / denominator
        fluorine_coverage = fluorine_ratio * available
        oxygen_coverage = oxygen_ratio * available

        chemical = (
            par.spontaneous_fluorine_removal_rate_m2_s * fluorine_coverage
            / par.fluorine_atoms_per_removed_si)
        enhanced_removal = ion_enhanced * fluorine_coverage
        removal_rate = chemical + physical + enhanced_removal
        removed = removal_rate * float(duration_s)
        updated = BelenSiliconState(
            available, state.removed_si_atoms_m2 + removed)

        fluorine_balance = (
            fluorine_supply * available - fluorine_loss * fluorine_coverage)
        oxygen_balance = oxygen_supply * available - oxygen_loss * oxygen_coverage
        fixed_point_change = available - state.available_site_fraction
        exchange = unresolved_surface_exchange(
            removed_units_m2={par.material_inventory_name: removed},
            limitations=(
                "volatile SiFx product branching and return transport are unresolved",
                "incident F/O site balance is reported separately from the target-material ledger",
            ))
        return BelenSiliconStepResult(
            state=updated,
            etch_velocity_m_s=removal_rate / par.bulk_si_atom_density_m3,
            fluorine_coverage=fluorine_coverage,
            oxygen_coverage=oxygen_coverage,
            available_site_fraction=available,
            chemical_removal_rate_m2_s=chemical,
            physical_sputter_rate_m2_s=physical,
            ion_enhanced_removal_rate_m2_s=enhanced_removal,
            removed_si_atoms_m2=removed,
            fluorine_site_balance_residual_m2_s=fluorine_balance,
            oxygen_site_balance_residual_m2_s=oxygen_balance,
            transport_fixed_point_change=fixed_point_change,
            material_exchange=exchange,
            validity=validity)
