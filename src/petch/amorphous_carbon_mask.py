"""Finite-film amorphous-carbon mask chemistry for the common feature engine.

The v1 state is deliberately smaller than Krüger's MCFPM mask lattice. It retains the channels
needed to test necking/clogging and mask erosion without creating a second feature solver:

* fluorocarbon deposition on exposed amorphous carbon and on an existing film,
* oxygen removal and energetic sputtering of that film,
* oxygen and energetic removal of exposed amorphous carbon.

The finite-film update reuses the conservative analytic operator in
``ReducedSiO2FluorocarbonMechanism``.  Krüger's resolved lattice distinguishes fresh and crosslinked
film: the latter has fewer open radical-attachment sites.  This reduced model exposes one bounded,
fixed effective crosslinked-growth fraction that blends those two *published* attachment laws.  It
does not pretend to resolve bond creation/scission or polymer identity.  Consequently the Krüger
parameter factory remains a development replay until calibration and held-out transfer establish
adequacy.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .material_mechanism_3d import MaterialMechanismRouter3D
from .surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange
from .surface_kinetics import (
    EnergeticYield,
    MechanismValidity,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
    SurfaceFluxes,
)


@dataclass(frozen=True)
class AmorphousCarbonMaskState:
    """Conservative per-area mask-film and cumulative carbon-removal inventories."""

    polymer_units_m2: np.ndarray | float = 0.0
    removed_carbon_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        polymer, removed = [
            np.array(value, copy=True)
            for value in np.broadcast_arrays(
                np.asarray(self.polymer_units_m2, dtype=float),
                np.asarray(self.removed_carbon_atoms_m2, dtype=float))
        ]
        if (np.any(~np.isfinite(polymer)) or np.any(polymer < 0.0)
                or np.any(~np.isfinite(removed)) or np.any(removed < 0.0)):
            raise ValueError("invalid amorphous-carbon mask state")
        polymer.setflags(write=False)
        removed.setflags(write=False)
        object.__setattr__(self, "polymer_units_m2", polymer)
        object.__setattr__(self, "removed_carbon_atoms_m2", removed)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape)
        return cls(zero, zero)

    def conservative_surface_fields(self):
        return {
            "polymer_units_m2": self.polymer_units_m2,
            "removed_carbon_atoms_m2": self.removed_carbon_atoms_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {
            "polymer_units_m2": None,
            "removed_carbon_atoms_m2": None,
        }

    def surface_field_remap_modes(self):
        return {
            "polymer_units_m2": "conservative",
            "removed_carbon_atoms_m2": "conservative",
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("amorphous-carbon mask remap fields do not match its state contract")
        return type(self)(
            fields["polymer_units_m2"], fields["removed_carbon_atoms_m2"])


@dataclass(frozen=True)
class AmorphousCarbonMaskParameters:
    """Physical inputs for a finite fluorocarbon film on amorphous carbon."""

    bulk_carbon_atom_density_m3: float
    polymer_unit_density_m3: float
    polymer_monolayer_density_m2: float
    polymer_deposition_probability_on_carbon: Mapping[str, float]
    polymer_deposition_probability_on_polymer: Mapping[str, float]
    polymer_deposition_probability_on_crosslinked_polymer: Mapping[str, float]
    effective_crosslinked_growth_fraction: float
    oxygen_species: str
    oxygen_polymer_etch_probability: float
    oxygen_carbon_etch_probability: float
    projectile_species: tuple[str, ...]
    polymer_sputter_yield: EnergeticYield
    carbon_sputter_yield: EnergeticYield
    declared_inert_neutral_species: tuple[str, ...]
    evidence: Mapping[str, ParameterEvidence]
    known_omissions: tuple[str, ...] = (
        "polymer identities and crosslink formation/scission are not resolved",
        "polymer carbonization and atomic-F chemistry are not resolved",
        "volatile mask and polymer product identities and return transport are unresolved",
    )

    def __post_init__(self):
        positive = (
            self.bulk_carbon_atom_density_m3,
            self.polymer_unit_density_m3,
            self.polymer_monolayer_density_m2,
        )
        probabilities = (
            self.oxygen_polymer_etch_probability,
            self.oxygen_carbon_etch_probability,
            self.effective_crosslinked_growth_fraction,
        )
        if (any(not np.isfinite(value) or value <= 0.0 for value in positive)
                or any(not np.isfinite(value) or not 0.0 <= value <= 1.0
                       for value in probabilities)
                or not self.oxygen_species
                or not isinstance(self.polymer_sputter_yield, EnergeticYield)
                or not isinstance(self.carbon_sputter_yield, EnergeticYield)):
            raise ValueError("invalid amorphous-carbon mask parameters")
        maps = {}
        for name in (
                "polymer_deposition_probability_on_carbon",
                "polymer_deposition_probability_on_polymer",
                "polymer_deposition_probability_on_crosslinked_polymer"):
            values = dict(getattr(self, name))
            if any((not species or not np.isfinite(value)
                    or value < 0.0 or value > 1.0)
                   for species, value in values.items()):
                raise ValueError(f"invalid mask probability map: {name}")
            maps[name] = MappingProxyType(values)
            object.__setattr__(self, name, maps[name])
        projectiles = tuple(str(name) for name in self.projectile_species)
        inert = tuple(str(name) for name in self.declared_inert_neutral_species)
        reactive = (
            set(maps["polymer_deposition_probability_on_carbon"])
            | set(maps["polymer_deposition_probability_on_polymer"])
            | {self.oxygen_species})
        if (not projectiles or any(not name for name in projectiles)
                or len(set(projectiles)) != len(projectiles)
                or any(not name for name in inert)
                or len(set(inert)) != len(inert)
                or set(inert) & reactive):
            raise ValueError("invalid mask projectile or inert-species declaration")
        evidence = dict(self.evidence)
        required = {
            "bulk_carbon_atom_density_m3",
            "polymer_unit_density_m3",
            "polymer_monolayer_density_m2",
            "polymer_deposition_probability_on_carbon",
            "polymer_deposition_probability_on_polymer",
            "polymer_deposition_probability_on_crosslinked_polymer",
            "effective_crosslinked_growth_fraction",
            "oxygen_polymer_etch_probability",
            "oxygen_carbon_etch_probability",
            "projectile_species",
            "polymer_sputter_yield",
            "carbon_sputter_yield",
            "declared_inert_neutral_species",
        }
        if (set(evidence) != required
                or any(not isinstance(item, ParameterEvidence)
                       for item in evidence.values())):
            raise ValueError("mask evidence must cover every physical input")
        object.__setattr__(self, "projectile_species", projectiles)
        object.__setattr__(self, "declared_inert_neutral_species", inert)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))

    @classmethod
    def krueger_2024_reduced_projection(
            cls, *, projectile_species=("ions",),
            effective_crosslinked_growth_fraction=0.0):
        """Return the four-feature calibrated Krüger mask projection.

        ``projectile_species`` names the explicit kinetic population used to represent the
        paper's aggregate ion/hot-neutral mixture. The bundled Krüger evidence has no
        species-resolved IEAD, so this remains nonpredictive until that boundary is supplied.

        ``effective_crosslinked_growth_fraction`` is a base-SEM calibration closure in [0, 1].
        Zero reproduces the former all-fresh-film projection.  One uses Appendix B's 0.02
        attachment probability for a crosslinked film.  Intermediate values blend collision
        probabilities, not outputs or profile velocities.  The fraction is fixed for every
        geometry/process in a validation campaign and does not replace a future dynamic bond model.
        """
        article = "https://doi.org/10.1116/6.0003554"
        thesis = "https://doi.org/10.7302/23106"
        site_source = "https://doi.org/10.1088/1361-6463/aa6f40"

        def evidence(source, evidence_type, note, *, supports=False):
            return ParameterEvidence(
                source, evidence_type, note=note,
                supports_prediction_within_declared_domain=supports)

        parameter_evidence = {
            "bulk_carbon_atom_density_m3": evidence(
                "assumed 1.8 g/cm3 amorphous-carbon density",
                "declared_material_assumption",
                "The experimental mask density is unreported; 9.03e28 C atoms/m3 "
                "must be sensitivity-bounded before a predictive mask-erosion claim."),
            "polymer_unit_density_m3": evidence(
                "https://github.com/ViennaTools/ViennaPS/tree/v4.6.1",
                "transferred_feature_model_parameter",
                "2e28 film units/m3 is a geometry conversion, not a Krüger fit."),
            "polymer_monolayer_density_m2": evidence(
                site_source, "transferred_surface_site_density",
                "1e15 cm^-2 converts the voxel film into an areal inventory."),
            "polymer_deposition_probability_on_carbon": evidence(
                article, "published_calibrated_parameter",
                "Table V four-feature pd,poly-AC=0.0842, applied to the four "
                "Appendix-B depositing radicals."),
            "polymer_deposition_probability_on_polymer": evidence(
                thesis, "published_model_reduced_projection",
                "Appendix-B un-crosslinked growth is 0.1 for CF/CF2/CF3 and "
                "0.03 for C2F3; film identity is collapsed."),
            "polymer_deposition_probability_on_crosslinked_polymer": evidence(
                thesis, "published_model_reduced_projection",
                "Appendix-B growth on crosslinked CF/CF2/CF3 sites is 0.02 for "
                "each of CF, CF2, CF3, and C2F3."),
            "effective_crosslinked_growth_fraction": evidence(
                article, "base_sem_calibration_parameter",
                "A single fixed blend between the published fresh/crosslinked growth "
                "laws replaces unresolved bond topology; it must be frozen before transfer."),
            "oxygen_polymer_etch_probability": evidence(
                article, "published_calibrated_parameter",
                "Table V four-feature pe,poly=0.0628."),
            "oxygen_carbon_etch_probability": evidence(
                thesis, "published_model_parameter",
                "Appendix-B AC + O reaction probability 1e-5."),
            "projectile_species": evidence(
                article, "unresolved_published_boundary_mixture",
                "Table I publishes only aggregate ion flux and Fig. 4 a combined IEAD."),
            "polymer_sputter_yield": evidence(
                thesis, "published_model_parameter",
                "Appendix-B un-crosslinked polymer p0=0.9, Eth=20 eV, n=0.5, "
                "Er=500 eV, Kress angular response."),
            "carbon_sputter_yield": evidence(
                thesis, "published_model_parameter",
                "Appendix-B AC p0=0.001, Eth=200 eV, n=0.4, Er=250 eV, "
                "Kress angular response."),
            "declared_inert_neutral_species": evidence(
                thesis, "published_absent_channel",
                "C3F4 is present in the HPEM flux table but absent from Appendix B."),
        }
        return cls(
            bulk_carbon_atom_density_m3=9.03e28,
            polymer_unit_density_m3=2.0e28,
            polymer_monolayer_density_m2=1.0e19,
            polymer_deposition_probability_on_carbon={
                "CF": 0.0842, "CF2": 0.0842, "CF3": 0.0842, "C2F3": 0.0842,
            },
            polymer_deposition_probability_on_polymer={
                "CF": 0.1, "CF2": 0.1, "CF3": 0.1, "C2F3": 0.03,
            },
            polymer_deposition_probability_on_crosslinked_polymer={
                "CF": 0.02, "CF2": 0.02, "CF3": 0.02, "C2F3": 0.02,
            },
            effective_crosslinked_growth_fraction=float(
                effective_crosslinked_growth_fraction),
            oxygen_species="O",
            oxygen_polymer_etch_probability=0.0628,
            oxygen_carbon_etch_probability=1.0e-5,
            projectile_species=tuple(projectile_species),
            polymer_sputter_yield=EnergeticYield(
                0.9, 20.0, 500.0, energy_exponent=0.5,
                angular_model="kress_1999", angular_parameter=9.3),
            carbon_sputter_yield=EnergeticYield(
                0.001, 200.0, 250.0, energy_exponent=0.4,
                angular_model="kress_1999", angular_parameter=9.3),
            declared_inert_neutral_species=("C3F4",),
            evidence=parameter_evidence,
            known_omissions=(
                "polymer identities and crosslink formation/scission are collapsed",
                "crosslinked growth is a fixed effective fraction rather than a dynamic bond state",
                "crosslinked-film sputter resistance is not represented",
                "polymer carbonization and atomic-F chemistry are not represented",
                "volatile mask and film product identities and return transport are unresolved",
                "the aggregate ion/hot-neutral composition and full IEAD require a boundary closure",
                "the amorphous-carbon density is assumed because the experiment does not report it",
            ),
        )


@dataclass(frozen=True)
class AmorphousCarbonMaskStepResult:
    state: AmorphousCarbonMaskState
    etch_velocity_m_s: np.ndarray
    normal_growth_velocity_m_s: np.ndarray
    deposited_polymer_units_m2: np.ndarray
    removed_polymer_units_m2: np.ndarray
    removed_carbon_atoms_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    validity: MechanismValidity
    product_populations: tuple = ()

    def __post_init__(self):
        if not isinstance(self.state, AmorphousCarbonMaskState):
            raise TypeError("invalid amorphous-carbon mask result state")
        for name in (
                "etch_velocity_m_s", "normal_growth_velocity_m_s",
                "deposited_polymer_units_m2", "removed_polymer_units_m2",
                "removed_carbon_atoms_m2"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError(f"invalid amorphous-carbon mask result field: {name}")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if not isinstance(self.material_exchange, SurfaceMaterialExchange):
            raise TypeError("mask result requires a material-exchange ledger")
        if not isinstance(self.validity, MechanismValidity):
            raise TypeError("mask result requires mechanism validity")
        object.__setattr__(self, "product_populations", tuple(self.product_populations))
        if self.product_populations:
            raise ValueError("mask v1 leaves volatile products explicitly unresolved")


class AmorphousCarbonMaskMechanism:
    """Finite-film mask law sharing the common engine and exact film integrator."""

    def __init__(self, parameters: AmorphousCarbonMaskParameters):
        if not isinstance(parameters, AmorphousCarbonMaskParameters):
            raise TypeError("parameters must be AmorphousCarbonMaskParameters")
        self.parameters = parameters
        par = parameters

        def yield_manifest(law):
            return {
                "reference_yield": float(law.reference_yield),
                "threshold_energy_eV": float(law.threshold_energy_eV),
                "reference_energy_eV": float(law.reference_energy_eV),
                "energy_exponent": float(law.energy_exponent),
                "angular_model": law.angular_model,
                "angular_parameter": (
                    None if law.angular_parameter is None
                    else float(law.angular_parameter)),
            }

        self.provenance = MappingProxyType({
            "model": "finite-film-amorphous-carbon-mask-common-engine-v1",
            "parameters": {
                "bulk_carbon_atom_density_m3": par.bulk_carbon_atom_density_m3,
                "polymer_unit_density_m3": par.polymer_unit_density_m3,
                "polymer_monolayer_density_m2": par.polymer_monolayer_density_m2,
                "polymer_deposition_probability_on_carbon": dict(
                    par.polymer_deposition_probability_on_carbon),
                "polymer_deposition_probability_on_polymer": dict(
                    par.polymer_deposition_probability_on_polymer),
                "polymer_deposition_probability_on_crosslinked_polymer": dict(
                    par.polymer_deposition_probability_on_crosslinked_polymer),
                "effective_crosslinked_growth_fraction": float(
                    par.effective_crosslinked_growth_fraction),
                "oxygen_species": par.oxygen_species,
                "oxygen_polymer_etch_probability": par.oxygen_polymer_etch_probability,
                "oxygen_carbon_etch_probability": par.oxygen_carbon_etch_probability,
                "projectile_species": list(par.projectile_species),
                "polymer_sputter_yield": yield_manifest(par.polymer_sputter_yield),
                "carbon_sputter_yield": yield_manifest(par.carbon_sputter_yield),
                "declared_inert_neutral_species": list(
                    par.declared_inert_neutral_species),
            },
            "sources": {
                name: {
                    "source": item.source,
                    "evidence_type": item.evidence_type,
                    "relative_standard_uncertainty": item.relative_standard_uncertainty,
                    "note": item.note,
                    "supports_prediction_within_declared_domain": (
                        item.supports_prediction_within_declared_domain),
                }
                for name, item in par.evidence.items()
            },
            "known_omissions": list(par.known_omissions),
        })
        # Reuse the exact finite-film inventory operator. Dummy oxide-removal channels have
        # identically zero yield; mask-carbon removal remains authoritative in this class.
        analytic = ParameterEvidence(
            "analytic zero channel in amorphous-carbon adapter", "analytic",
            supports_prediction_within_declared_domain=True)
        crosslinked_fraction = float(par.effective_crosslinked_growth_fraction)
        effective_polymer_growth = {
            species: (
                (1.0 - crosslinked_fraction)
                * par.polymer_deposition_probability_on_polymer.get(species, 0.0)
                + crosslinked_fraction
                * par.polymer_deposition_probability_on_crosslinked_polymer.get(
                    species,
                    par.polymer_deposition_probability_on_polymer.get(species, 0.0)))
            for species in (
                set(par.polymer_deposition_probability_on_polymer)
                | set(par.polymer_deposition_probability_on_crosslinked_polymer))
        }
        self._effective_polymer_growth_probability = MappingProxyType(
            effective_polymer_growth)
        internal_evidence = {
            "site_density_m2": analytic,
            "bulk_formula_density_m3": analytic,
            "polymer_monolayer_density_m2": par.evidence[
                "polymer_monolayer_density_m2"],
            "complex_formation_probability": analytic,
            "polymer_deposition_probability_on_substrate": par.evidence[
                "polymer_deposition_probability_on_carbon"],
            "polymer_deposition_probability_on_polymer": par.evidence[
                "effective_crosslinked_growth_fraction"],
            "oxygen_polymer_etch_probability": par.evidence[
                "oxygen_polymer_etch_probability"],
            "bare_sio2_yield": analytic,
            "complex_sio2_yield": analytic,
            "polymer_sputter_yield": par.evidence["polymer_sputter_yield"],
            "declared_inert_neutral_species": par.evidence[
                "declared_inert_neutral_species"],
        }
        self._film = ReducedSiO2FluorocarbonMechanism(
            ReducedSiO2FluorocarbonParameters(
                site_density_m2=par.polymer_monolayer_density_m2,
                bulk_formula_density_m3=1.0,
                polymer_monolayer_density_m2=par.polymer_monolayer_density_m2,
                complex_formation_probability={},
                polymer_deposition_probability_on_substrate=(
                    par.polymer_deposition_probability_on_carbon),
                polymer_deposition_probability_on_polymer=(
                    effective_polymer_growth),
                oxygen_species=par.oxygen_species,
                oxygen_polymer_etch_probability=par.oxygen_polymer_etch_probability,
                bare_sio2_yield=EnergeticYield(0.0, 0.0, 1.0),
                complex_sio2_yield=EnergeticYield(0.0, 0.0, 1.0),
                polymer_sputter_yield=par.polymer_sputter_yield,
                declared_inert_neutral_species=par.declared_inert_neutral_species,
                evidence=internal_evidence,
                known_omissions=par.known_omissions))

    @staticmethod
    def initial_state(shape=()):
        return AmorphousCarbonMaskState.bare(shape)

    @staticmethod
    def _broadcast(value, shape, label):
        try:
            return np.broadcast_to(np.asarray(value, dtype=float), shape)
        except ValueError as error:
            raise ValueError(f"{label} does not match the mask state shape") from error

    def _selected_fluxes(self, fluxes):
        selected = tuple(
            population for population in fluxes.energetic_fluxes
            if population.name in self.parameters.projectile_species)
        return SurfaceFluxes(fluxes.neutral_flux_m2_s, selected)

    def _energetic_rate(self, fluxes, law, shape):
        total = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name in self.parameters.projectile_species:
                total = total + self._broadcast(
                    population.yield_rate_m2_s(law), shape, population.name)
        return total

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        supported_neutral = (
            set(par.polymer_deposition_probability_on_carbon)
            | set(par.polymer_deposition_probability_on_polymer)
            | {par.oxygen_species}
            | set(par.declared_inert_neutral_species))
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in supported_neutral and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in par.projectile_species
            and np.any(np.asarray(population.flux_m2_s) > 0.0)}))
        reasons = []
        if unsupported_neutral:
            reasons.append("positive neutral flux has no declared mask reaction channel")
        if unsupported_energetic:
            reasons.append(
                "positive energetic flux has no declared mask channel: "
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

    def neutral_reaction_probability(self, state: AmorphousCarbonMaskState):
        if not isinstance(state, AmorphousCarbonMaskState):
            raise TypeError("mask neutral probabilities require AmorphousCarbonMaskState")
        internal = SiO2SurfaceState(
            np.zeros(state.polymer_units_m2.shape), state.polymer_units_m2)
        return self._film.neutral_reaction_probability(internal)

    def _instantaneous_film_rates(self, state, fluxes):
        par = self.parameters
        shape = state.polymer_units_m2.shape
        access = np.exp(
            -state.polymer_units_m2 / par.polymer_monolayer_density_m2)
        coverage = 1.0 - access
        deposition = np.zeros(shape)
        for species in (
                set(par.polymer_deposition_probability_on_carbon)
                | set(self._effective_polymer_growth_probability)):
            incident = self._broadcast(
                fluxes.neutral_flux_m2_s.get(species, 0.0), shape, species)
            deposition = deposition + incident * (
                par.polymer_deposition_probability_on_carbon.get(species, 0.0) * access
                + self._effective_polymer_growth_probability.get(species, 0.0)
                * coverage)
        oxygen = self._broadcast(
            fluxes.neutral_flux_m2_s.get(par.oxygen_species, 0.0),
            shape, par.oxygen_species)
        removal = (
            oxygen * par.oxygen_polymer_etch_probability
            + self._energetic_rate(fluxes, par.polymer_sputter_yield, shape)
        ) * coverage
        return deposition, removal

    def advance(
            self, state: AmorphousCarbonMaskState, fluxes: SurfaceFluxes,
            duration_s: float, *, strict=True):
        if not isinstance(state, AmorphousCarbonMaskState):
            raise TypeError("mask advance requires AmorphousCarbonMaskState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: " + "; ".join(validity.reasons))
        par = self.parameters
        shape = state.polymer_units_m2.shape
        selected = self._selected_fluxes(fluxes)
        initial_film = SiO2SurfaceState(
            np.zeros(shape), state.polymer_units_m2)
        half = self._film.advance(initial_film, selected, 0.5 * float(duration_s))
        midpoint_access = np.exp(
            -half.state.polymer_units_m2 / par.polymer_monolayer_density_m2)
        oxygen = self._broadcast(
            fluxes.neutral_flux_m2_s.get(par.oxygen_species, 0.0),
            shape, par.oxygen_species)
        carbon_removal_rate = midpoint_access * (
            oxygen * par.oxygen_carbon_etch_probability
            + self._energetic_rate(fluxes, par.carbon_sputter_yield, shape))
        removed_carbon = carbon_removal_rate * float(duration_s)
        second = self._film.advance(half.state, selected, 0.5 * float(duration_s))
        deposited_polymer = (
            half.deposited_polymer_units_m2 + second.deposited_polymer_units_m2)
        removed_polymer = (
            half.removed_polymer_units_m2 + second.removed_polymer_units_m2)
        updated = AmorphousCarbonMaskState(
            second.state.polymer_units_m2,
            state.removed_carbon_atoms_m2 + removed_carbon)

        if duration_s > 0.0:
            deposition_rate = deposited_polymer / float(duration_s)
            film_removal_rate = removed_polymer / float(duration_s)
        else:
            deposition_rate, film_removal_rate = self._instantaneous_film_rates(
                state, selected)
        etch_velocity = (
            carbon_removal_rate / par.bulk_carbon_atom_density_m3
            + film_removal_rate / par.polymer_unit_density_m3)
        growth_velocity = deposition_rate / par.polymer_unit_density_m3
        exchange = unresolved_surface_exchange(
            removed_units_m2={
                "amorphous_carbon_atom": removed_carbon,
                "fluorocarbon_film_unit": removed_polymer,
            },
            deposited_units_m2={
                "fluorocarbon_film_unit": deposited_polymer,
            },
            limitations=(
                "volatile amorphous-carbon and fluorocarbon products are unresolved",
                "unresolved removed material is not eligible for return transport",
                "gross film recession and growth are reported separately so their "
                "difference drives the common level-set velocity",
            ))
        return AmorphousCarbonMaskStepResult(
            state=updated,
            etch_velocity_m_s=etch_velocity,
            normal_growth_velocity_m_s=growth_velocity,
            deposited_polymer_units_m2=deposited_polymer,
            removed_polymer_units_m2=removed_polymer,
            removed_carbon_atoms_m2=removed_carbon,
            material_exchange=exchange,
            validity=validity)


def build_krueger_2024_material_router_3d(
        *, oxide_material_id=1, mask_material_id=2, projectile_species=("ions",),
        effective_mask_crosslinked_growth_fraction=0.0,
        oxide_etch_yield_scale=1.0):
    """Build one material router for the reduced Krüger oxide/mask development replay."""
    oxide_id = int(oxide_material_id)
    mask_id = int(mask_material_id)
    if oxide_id <= 0 or mask_id <= 0 or oxide_id == mask_id:
        raise ValueError("Krüger oxide and mask material ids must be distinct and positive")
    oxide = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.krueger_2024_reduced_projection(
            oxide_etch_yield_scale=oxide_etch_yield_scale))
    mask = AmorphousCarbonMaskMechanism(
        AmorphousCarbonMaskParameters.krueger_2024_reduced_projection(
            projectile_species=projectile_species,
            effective_crosslinked_growth_fraction=(
                effective_mask_crosslinked_growth_fraction)))
    return MaterialMechanismRouter3D(
        {oxide_id: oxide, mask_id: mask},
        provenance={
            oxide_id: {
                "role": "SiO2 substrate",
                "model": dict(oxide.provenance),
                "claim_status": "development_replay",
            },
            mask_id: {
                "role": "amorphous-carbon mask",
                "model": dict(mask.provenance),
                "claim_status": "development_replay",
            },
        })
