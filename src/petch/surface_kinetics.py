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

from .surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange


@dataclass(frozen=True)
class ParameterEvidence:
    """Provenance and uncertainty for one physical chemistry input."""

    source: str
    evidence_type: str
    relative_standard_uncertainty: float | None = None
    note: str = ""
    supports_prediction_within_declared_domain: bool = False

    def __post_init__(self):
        if not self.source or not self.evidence_type:
            raise ValueError("parameter evidence requires source and evidence_type")
        uncertainty = self.relative_standard_uncertainty
        if uncertainty is not None and (not np.isfinite(uncertainty) or uncertainty < 0.0):
            raise ValueError("relative standard uncertainty must be finite and nonnegative")
        if not isinstance(self.supports_prediction_within_declared_domain, bool):
            raise TypeError("prediction-support flag must be boolean")


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

    def yield_rate_m2_s(self, yield_law: EnergeticYield):
        return np.asarray(self.flux_m2_s) * self.mean_yield(yield_law)


@dataclass(frozen=True)
class FaceResolvedEnergeticFlux:
    """Exact sparse energetic hit events on a surface mesh.

    ``event_flux_m2_s`` is each event's contribution to the dimensional flux density of its hit
    face.  Retaining events avoids energy/angle histogram binning and mean-energy yield bias.  A
    production transport backend may compact identical events, but it must preserve this measure.
    """

    name: str
    face_count: int
    event_face: np.ndarray
    event_flux_m2_s: np.ndarray
    event_energy_eV: np.ndarray
    event_cosine_incidence: np.ndarray
    event_position: np.ndarray | None = None
    event_incident_direction: np.ndarray | None = None

    def __post_init__(self):
        if not self.name or int(self.face_count) <= 0:
            raise ValueError("face-resolved energetic flux requires a name and positive face_count")
        face = np.asarray(self.event_face, dtype=int).copy()
        flux = np.asarray(self.event_flux_m2_s, dtype=float).copy()
        energy = np.asarray(self.event_energy_eV, dtype=float).copy()
        cosine = np.asarray(self.event_cosine_incidence, dtype=float).copy()
        if (face.ndim != 1 or flux.shape != face.shape or energy.shape != face.shape
                or cosine.shape != face.shape or np.any(face < 0) or np.any(face >= self.face_count)
                or np.any(~np.isfinite(flux)) or np.any(flux < 0.0)
                or np.any(~np.isfinite(energy)) or np.any(energy < 0.0)
                or np.any(~np.isfinite(cosine)) or np.any((cosine < 0.0) | (cosine > 1.0))):
            raise ValueError("invalid face-resolved energetic events")
        for array in (face, flux, energy, cosine): array.setflags(write=False)
        object.__setattr__(self, "face_count", int(self.face_count))
        object.__setattr__(self, "event_face", face)
        object.__setattr__(self, "event_flux_m2_s", flux)
        object.__setattr__(self, "event_energy_eV", energy)
        object.__setattr__(self, "event_cosine_incidence", cosine)
        for name in ("event_position", "event_incident_direction"):
            supplied = getattr(self, name)
            if supplied is None:
                continue
            value = np.asarray(supplied, dtype=float).copy()
            if value.shape != (face.size, 3) or np.any(~np.isfinite(value)):
                raise ValueError(f"{name} must be a finite (event_count, 3) array")
            if (name == "event_incident_direction" and value.size
                    and not np.allclose(
                        np.linalg.norm(value, axis=1), 1.0, rtol=0.0, atol=2e-6)):
                raise ValueError("event incident directions must be unit vectors")
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    @property
    def flux_m2_s(self):
        return np.bincount(
            self.event_face, weights=self.event_flux_m2_s, minlength=self.face_count)

    def yield_rate_m2_s(self, yield_law: EnergeticYield):
        event_yield = yield_law.evaluate(self.event_energy_eV, self.event_cosine_incidence)
        return np.bincount(
            self.event_face, weights=self.event_flux_m2_s * event_yield,
            minlength=self.face_count)


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
        if any(not isinstance(item, (EnergeticFlux, FaceResolvedEnergeticFlux))
               for item in energetic):
            raise TypeError(
                "energetic_fluxes must contain EnergeticFlux or FaceResolvedEnergeticFlux objects")
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

    def conservative_surface_fields(self):
        return {
            "complex_fraction": self.complex_fraction,
            "polymer_units_m2": self.polymer_units_m2,
            "removed_formula_units_m2": self.removed_formula_units_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {
            "complex_fraction": 1.0,
            "polymer_units_m2": None,
            "removed_formula_units_m2": None,
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("SiO2 remap fields do not match its state contract")
        return type(self)(
            fields["complex_fraction"], fields["polymer_units_m2"],
            fields["removed_formula_units_m2"])


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
    parameter_evidence_supports_prediction: bool
    nonpredictive_parameters: tuple[str, ...]


@dataclass(frozen=True)
class SurfaceStepResult:
    state: SiO2SurfaceState
    etch_velocity_m_s: np.ndarray
    formed_complex_units_m2: np.ndarray
    removed_complex_units_m2: np.ndarray
    removed_bare_formula_units_m2: np.ndarray
    deposited_polymer_units_m2: np.ndarray
    removed_polymer_units_m2: np.ndarray
    material_exchange: SurfaceMaterialExchange
    validity: MechanismValidity


class ReducedSiO2FluorocarbonMechanism:
    """Vectorized conservative kernel for the declared reduced reaction network."""

    def __init__(self, parameters: ReducedSiO2FluorocarbonParameters):
        self.parameters = parameters

    @staticmethod
    def initial_state(shape=()):
        return SiO2SurfaceState.bare(shape)

    def neutral_reaction_probability(self, state: SiO2SurfaceState):
        """Per-collision probability that each neutral leaves the ballistic population.

        The loss is the sum of the same competing channels advanced by this mechanism: complex
        formation on accessible uncomplexed oxide, deposition on bare/polymer-covered surface, and
        oxygen removal of polymer. Transport may diffusely re-emit the remaining probability.
        """
        if not isinstance(state, SiO2SurfaceState):
            raise TypeError("neutral reaction probabilities require SiO2SurfaceState")
        par = self.parameters
        access = np.exp(-state.polymer_units_m2 / par.polymer_monolayer_density_m2)
        polymer_coverage = 1.0 - access
        probability = {}

        def add(species, value):
            probability[species] = probability.get(species, np.zeros_like(access)) + value

        for species, value in par.complex_formation_probability.items():
            add(species, value * access * (1.0 - state.complex_fraction))
        all_deposition_species = (
            set(par.polymer_deposition_probability_on_substrate)
            | set(par.polymer_deposition_probability_on_polymer))
        for species in all_deposition_species:
            add(species,
                par.polymer_deposition_probability_on_substrate.get(species, 0.0) * access
                + par.polymer_deposition_probability_on_polymer.get(species, 0.0)
                * polymer_coverage)
        add(par.oxygen_species, par.oxygen_polymer_etch_probability * polymer_coverage)
        if any(np.any(value > 1.0 + 5e-14) for value in probability.values()):
            raise ValueError(
                "competing neutral reaction probabilities exceed one; chemistry inputs are invalid")
        return MappingProxyType({
            species: np.minimum(np.maximum(value, 0.0), 1.0)
            for species, value in probability.items()})

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
        nonpredictive = tuple(sorted(
            name for name in required_evidence
            if name not in par.evidence
            or not par.evidence[name].supports_prediction_within_declared_domain))
        reasons = []
        if unsupported:
            reasons.append("positive incident neutral flux has no declared reaction channel")
        if missing_evidence:
            reasons.append("missing parameter evidence: " + ", ".join(missing_evidence))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported,
            known_model_form_omissions=par.known_omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive)

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
            total = total + self._broadcast(population.yield_rate_m2_s(yield_law), shape)
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
        active_reaction = deposit_substrate + deposit_polymer + removal_capacity > 0.0
        if not np.any(active_reaction):
            return np.array(inventory, copy=True), np.zeros(shape), np.zeros(shape)

        # dN/dt = D_sub + (D_poly-D_sub-R) * (1-exp(-N/N_mono)).  With y=exp(N/N_mono)
        # this is a linear ODE, so positivity and the N=0 boundary are preserved without clipping.
        coefficient = deposit_polymer - removal_capacity
        transition = deposit_polymer - deposit_substrate - removal_capacity
        exponent = coefficient * duration_s / monolayer
        log_y0 = inventory / monolayer
        small = np.abs(coefficient) <= 1e-14 * np.maximum(
            deposit_substrate + deposit_polymer + removal_capacity, 1.0)
        delta_log_y = np.empty(shape)
        if np.any(small):
            # coefficient -> 0: y1=y0-transition*t/M. Evaluate as log(y0)+log1p(delta/y0)
            # so an already thick polymer never forms exp(N/M).
            delta = -transition[small] * duration_s / monolayer
            ratio = delta * np.exp(-log_y0[small])
            if np.any(1.0 + ratio <= 0.0):
                raise RuntimeError("polymer inventory integrator violated positivity")
            delta_log_y[small] = np.log1p(ratio)
        if np.any(~small):
            # y1/y0 = exp(x)*d + (1-d), where d=1-a/y0 and
            # a=transition/coefficient.  Evaluate log(y1/y0) directly so a thick film's small
            # increment never comes from subtracting two O(N/M) logarithms.
            selected_log_y0 = log_y0[~small]
            a = transition[~small] / coefficient[~small]
            difference = 1.0 - a * np.exp(-selected_log_y0)
            selected_exponent = exponent[~small]
            local_delta = np.empty(selected_exponent.shape)
            growing = selected_exponent >= 0.0
            if np.any(growing):
                argument = (difference[growing]
                            + (1.0 - difference[growing])
                            * np.exp(-selected_exponent[growing]))
                if np.any(argument <= 0.0):
                    raise RuntimeError("polymer inventory integrator violated positivity")
                local_delta[growing] = selected_exponent[growing] + np.log(argument)
            if np.any(~growing):
                # Here coefficient<0 implies a>1.  When y0 is much larger than a and x is very
                # negative, log1p(expm1(x)*d) rounds to log(0) even though the exact limiting ratio
                # is a/y0.  Resolve the two positive terms in log space; use a factored difference
                # only when y0<a and the exponential term is subtractive.
                local_a = a[~growing]
                local_log_y0 = selected_log_y0[~growing]
                local_difference = difference[~growing]
                nonnegative = local_difference >= 0.0
                negative_delta = np.empty(local_difference.shape)
                if np.any(nonnegative):
                    log_first = np.where(
                        local_difference[nonnegative] > 0.0,
                        selected_exponent[~growing][nonnegative]
                        + np.log(np.maximum(
                            local_difference[nonnegative], np.finfo(float).tiny)),
                        -np.inf)
                    log_second = np.log(local_a[nonnegative]) - local_log_y0[nonnegative]
                    negative_delta[nonnegative] = np.logaddexp(log_first, log_second)
                if np.any(~nonnegative):
                    log_b = np.log(local_a[~nonnegative]) - local_log_y0[~nonnegative]
                    subtract_fraction = (
                        np.exp(selected_exponent[~growing][~nonnegative])
                        * (-np.expm1(-log_b)))
                    if np.any(subtract_fraction >= 1.0):
                        raise RuntimeError("polymer inventory integrator violated positivity")
                    negative_delta[~nonnegative] = log_b + np.log1p(-subtract_fraction)
                local_delta[~growing] = negative_delta
            if np.any(~np.isfinite(local_delta)):
                raise RuntimeError("polymer inventory integrator violated positivity")
            delta_log_y[~small] = local_delta
        log_y1 = log_y0 + delta_log_y
        if np.any(log_y1 < -1e-12):
            raise RuntimeError("polymer inventory integrator violated positivity")
        # Carry the analytic increment directly.  Reconstructing it later as N1-N0 loses several
        # digits once a thick film (O(1e21) m^-2) changes by only O(1e17) m^-2 in one substep.
        # That cancellation is numerical bookkeeping error, not a physical non-conservation.
        at_inventory_floor = log_y1 <= 0.0
        delta = np.where(
            at_inventory_floor, -inventory, monolayer * delta_log_y)
        updated = np.where(at_inventory_floor, 0.0, inventory + delta)
        # Preserve the exact no-flux identity face-by-face. Otherwise exp(log(N)) changes a large
        # inactive inventory by an ULP when a different face is active, creating a false local source.
        updated = np.where(active_reaction, updated, inventory)
        delta = np.where(active_reaction, delta, 0.0)
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
        polymer_balance = delta - (deposited - removed)
        polymer_scale = np.maximum.reduce((
            np.abs(delta), np.abs(deposited), np.abs(removed), np.ones(shape)))
        relative_balance = np.abs(polymer_balance) / polymer_scale
        if np.any(relative_balance > 5e-13):
            failed = np.unravel_index(np.argmax(relative_balance), relative_balance.shape)
            raise RuntimeError(
                "polymer reaction bookkeeping failed conservation: "
                f"max relative residual={float(relative_balance[failed]):.3e}, "
                f"Dsub={float(deposit_substrate[failed]):.6e}, "
                f"Dpoly={float(deposit_polymer[failed]):.6e}, "
                f"R={float(removal_capacity[failed]):.6e}, "
                f"N0={float(inventory[failed]):.6e}, N1={float(updated[failed]):.6e}, "
                f"deposited={float(deposited[failed]):.6e}, "
                f"removed={float(removed[failed]):.6e}")
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
        active = total_hazard > 0.0
        updated = np.where(active, updated, complex_fraction)
        integral_complex = np.empty(shape)
        integral_complex[active] = (
            equilibrium[active] * duration_s
            + (complex_fraction[active] - equilibrium[active])
            * (-np.expm1(-total_hazard[active] * duration_s)) / total_hazard[active])
        integral_complex[~active] = complex_fraction[~active] * duration_s
        removed_complex = complex_removal_event_rate * integral_complex
        formed_from_rate = formation_event_rate * (duration_s - integral_complex)
        removed_bare = bare_removal_event_rate * (duration_s - integral_complex)
        site_change = (updated - complex_fraction) * par.site_density_m2
        # This identity is the conservative state equation and is better conditioned near saturation
        # than subtracting two large rate integrals. Keep the independent rate integral as an audit.
        formed_complex = site_change + removed_complex
        rate_consistency = formed_from_rate - formed_complex
        site_scale = np.maximum.reduce((
            np.abs(site_change), np.abs(formed_from_rate), np.abs(removed_complex), np.ones(shape)))
        relative_consistency = np.abs(rate_consistency) / site_scale
        # A fractional coverage cannot represent fewer than O(eps*site_density) surface units.
        # Use that absolute floating-point floor in addition to the relative analytic audit; this is
        # relevant only on essentially unilluminated faces, where both physical rates are sub-ULP.
        roundoff_floor = (8.0 * np.finfo(float).eps * par.site_density_m2
                          * np.maximum(np.abs(complex_fraction), 1.0))
        failed_mask = np.abs(rate_consistency) > np.maximum(5e-12 * site_scale, roundoff_floor)
        if np.any(failed_mask):
            failed_score = np.where(failed_mask, relative_consistency, -np.inf)
            failed = np.unravel_index(np.argmax(failed_score), relative_consistency.shape)
            raise RuntimeError(
                "oxide-complex analytic rate integral disagrees with conserved state update: "
                f"max relative residual={float(relative_consistency[failed]):.3e}, "
                f"formation={float(formation_event_rate[failed]):.6e}, "
                f"removal={float(complex_removal_event_rate[failed]):.6e}, "
                f"theta0={float(complex_fraction[failed]):.6e}, "
                f"theta1={float(updated[failed]):.6e}, "
                f"formed_rate={float(formed_from_rate[failed]):.6e}, "
                f"formed_conserved={float(formed_complex[failed]):.6e}, "
                f"removed={float(removed_complex[failed]):.6e}")
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
        removed_sio2 = removed_complex + removed_bare
        exchange = unresolved_surface_exchange(
            removed_units_m2={
                "SiO2_formula_unit": removed_sio2,
                "fluorocarbon_film_unit": removed_polymer,
            },
            deposited_units_m2={"fluorocarbon_film_unit": deposited_polymer},
            limitations=(
                "reactive SiO2 and fluorocarbon product identities and branching are unresolved",
                "unresolved removed material is not eligible for redeposition transport",
            ))
        return SurfaceStepResult(
            state=new_state,
            etch_velocity_m_s=np.asarray(velocity),
            formed_complex_units_m2=formed_complex,
            removed_complex_units_m2=removed_complex,
            removed_bare_formula_units_m2=removed_bare,
            deposited_polymer_units_m2=deposited_polymer,
            removed_polymer_units_m2=removed_polymer,
            material_exchange=exchange,
            validity=validity)
