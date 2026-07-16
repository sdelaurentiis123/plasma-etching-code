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

The optional second-order complex-removal activation follows the nearest-neighbour mixing-layer
construction in W. Guo, MIT PhD thesis (2009), Sec. 4.3: the COF2-forming event probability contains a
C--O factor and two C--F bond factors.  It reuses the same bounded complex state; it does not introduce an
unobserved film state or alter transport.  Order one remains the default for backward compatibility.
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


def _angular_yield_factor(cosine, angular_model, angular_parameter):
    """Shared, dimensionless incidence-angle factor for energetic surface yields."""
    cosine = np.asarray(cosine, dtype=float)
    if np.any(~np.isfinite(cosine)) or np.any((cosine < 0.0) | (cosine > 1.0)):
        raise ValueError("incidence cosines must lie in [0, 1]")
    if angular_model == "none":
        return np.ones(cosine.shape)
    if angular_model == "chang_sawin_1997":
        theta = np.arccos(cosine)
        return np.where(
            cosine >= 0.5, 1.0, np.maximum(3.0 - 6.0 * theta / np.pi, 0.0))
    if angular_model == "kress_1999":
        return ((1.0 + float(angular_parameter) * (1.0 - cosine * cosine))
                * cosine)
    raise ValueError(f"unknown angular yield model: {angular_model}")


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
        scaled_energy = np.maximum(
            (energy - self.threshold_energy_eV)
            / (self.reference_energy_eV - self.threshold_energy_eV), 0.0)
        angular = _angular_yield_factor(
            cosine, self.angular_model, self.angular_parameter)
        return self.reference_yield * scaled_energy ** self.energy_exponent * angular


@dataclass(frozen=True)
class LowEnergyActivationYield:
    """Low-energy activation law used by Huang--Kushner surface states.

    Huang et al. (JVST A 37, 031304, 2019), Eq. (2), use

    ``p(E, theta) = p0 * max(0, 1 - E / E_max) * f(theta)``.

    Their reaction table also lists a 5 eV lower cutoff: energetic partners below that
    energy are reclassified as thermal neutrals and do not activate the surface.  Keeping
    this law separate from :class:`EnergeticYield` prevents a high-energy sputter law from
    being accidentally reused as a low-energy activation law.
    """

    zero_energy_yield: float
    minimum_energy_eV: float
    maximum_energy_eV: float
    angular_model: str = "none"
    angular_parameter: float | None = None

    def __post_init__(self):
        if (not np.isfinite(self.zero_energy_yield) or self.zero_energy_yield < 0.0
                or not np.isfinite(self.minimum_energy_eV)
                or self.minimum_energy_eV < 0.0
                or not np.isfinite(self.maximum_energy_eV)
                or self.maximum_energy_eV <= self.minimum_energy_eV):
            raise ValueError("invalid low-energy activation-yield parameters")
        if self.angular_model not in {"none", "chang_sawin_1997", "kress_1999"}:
            raise ValueError(f"unknown angular yield model: {self.angular_model}")
        if self.angular_model == "kress_1999" and (
                self.angular_parameter is None
                or not np.isfinite(self.angular_parameter)
                or self.angular_parameter < 0.0):
            raise ValueError("kress_1999 requires a nonnegative angular_parameter")

    def evaluate(self, energy_eV, cosine_incidence):
        energy = np.asarray(energy_eV, dtype=float)
        cosine = np.asarray(cosine_incidence, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
            raise ValueError("incident energies must be finite and nonnegative")
        energy_factor = np.maximum(1.0 - energy / self.maximum_energy_eV, 0.0)
        energy_factor = np.where(energy >= self.minimum_energy_eV, energy_factor, 0.0)
        angular = _angular_yield_factor(
            cosine, self.angular_model, self.angular_parameter)
        return self.zero_energy_yield * energy_factor * angular


@dataclass(frozen=True)
class SteinbruchelYield:
    """Square-root threshold yield used by the Belen/ViennaPS silicon model.

    ``Y(E, theta) = A * max(sqrt(E) - sqrt(E_th), 0) * f(theta)``.
    This is kept separate from :class:`EnergeticYield`: the two energy laws are not
    algebraically interchangeable away from their reference energy.
    """

    prefactor_per_sqrt_eV: float
    threshold_energy_eV: float
    angular_model: str = "none"
    angular_parameter: float | None = None

    def __post_init__(self):
        if (not np.isfinite(self.prefactor_per_sqrt_eV)
                or self.prefactor_per_sqrt_eV < 0.0
                or not np.isfinite(self.threshold_energy_eV)
                or self.threshold_energy_eV < 0.0):
            raise ValueError("invalid Steinbruchel-yield parameters")
        if self.angular_model not in {"none", "chang_sawin_1997", "kress_1999"}:
            raise ValueError(f"unknown angular yield model: {self.angular_model}")
        if self.angular_model == "kress_1999" and (
                self.angular_parameter is None
                or not np.isfinite(self.angular_parameter)
                or self.angular_parameter < 0.0):
            raise ValueError("kress_1999 requires a nonnegative angular_parameter")

    def evaluate(self, energy_eV, cosine_incidence):
        energy = np.asarray(energy_eV, dtype=float)
        cosine = np.asarray(cosine_incidence, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
            raise ValueError("incident energies must be finite and nonnegative")
        energy_factor = np.maximum(
            np.sqrt(energy) - np.sqrt(self.threshold_energy_eV), 0.0)
        angular = _angular_yield_factor(
            cosine, self.angular_model, self.angular_parameter)
        return self.prefactor_per_sqrt_eV * energy_factor * angular


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
    ``activated_complex_fraction`` is the absolute fraction of oxide sites that are both
    complexed and low-energy activated; it may not exceed ``complex_fraction``.
    ``activated_polymer_fraction`` is the conditional activated fraction of the exposed
    fluorocarbon-polymer surface.  The two states remain separate because Huang et al. assign
    different activation energy windows to complex (5--70 eV) and polymer (5--30 eV) sites.
    """

    complex_fraction: np.ndarray | float
    polymer_units_m2: np.ndarray | float
    removed_formula_units_m2: np.ndarray | float = 0.0
    activated_complex_fraction: np.ndarray | float = 0.0
    activated_polymer_fraction: np.ndarray | float = 0.0

    def __post_init__(self):
        arrays = np.broadcast_arrays(
            np.asarray(self.complex_fraction, dtype=float),
            np.asarray(self.polymer_units_m2, dtype=float),
            np.asarray(self.removed_formula_units_m2, dtype=float),
            np.asarray(self.activated_complex_fraction, dtype=float),
            np.asarray(self.activated_polymer_fraction, dtype=float))
        complex_fraction, polymer, removed, activated_complex, activated_polymer = [
            np.array(item, copy=True) for item in arrays]
        if (np.any(~np.isfinite(complex_fraction))
                or np.any((complex_fraction < 0.0) | (complex_fraction > 1.0))
                or np.any(~np.isfinite(polymer)) or np.any(polymer < 0.0)
                or np.any(~np.isfinite(removed)) or np.any(removed < 0.0)
                or np.any(~np.isfinite(activated_complex))
                or np.any(activated_complex < 0.0)
                or np.any(activated_complex > complex_fraction + 16.0 * np.finfo(float).eps)
                or np.any(~np.isfinite(activated_polymer))
                or np.any((activated_polymer < 0.0) | (activated_polymer > 1.0))):
            raise ValueError("invalid SiO2 surface state")
        # Erase only representational excursions at the dependent upper bound.
        activated_complex = np.minimum(activated_complex, complex_fraction)
        activated_polymer = np.where(polymer > 0.0, activated_polymer, 0.0)
        for array in (
                complex_fraction, polymer, removed,
                activated_complex, activated_polymer):
            array.setflags(write=False)
        object.__setattr__(self, "complex_fraction", complex_fraction)
        object.__setattr__(self, "polymer_units_m2", polymer)
        object.__setattr__(self, "removed_formula_units_m2", removed)
        object.__setattr__(self, "activated_complex_fraction", activated_complex)
        object.__setattr__(self, "activated_polymer_fraction", activated_polymer)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape)
        return cls(zero, zero, zero, zero, zero)

    def conservative_surface_fields(self):
        return {
            "complex_fraction": self.complex_fraction,
            "polymer_units_m2": self.polymer_units_m2,
            "removed_formula_units_m2": self.removed_formula_units_m2,
            "activated_complex_fraction": self.activated_complex_fraction,
            "activated_polymer_fraction": self.activated_polymer_fraction,
        }

    def conservative_surface_upper_bounds(self):
        return {
            "complex_fraction": 1.0,
            "polymer_units_m2": None,
            "removed_formula_units_m2": None,
            "activated_complex_fraction": 1.0,
            "activated_polymer_fraction": 1.0,
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("SiO2 remap fields do not match its state contract")
        return type(self)(
            fields["complex_fraction"], fields["polymer_units_m2"],
            fields["removed_formula_units_m2"],
            np.minimum(
                fields["activated_complex_fraction"], fields["complex_fraction"]),
            np.where(
                np.asarray(fields["polymer_units_m2"]) > 0.0,
                fields["activated_polymer_fraction"], 0.0))


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
    complex_removal_reaction_order: int = 1
    activated_polymer_deposition_probability_on_substrate: Mapping[str, float] = field(
        default_factory=dict)
    activated_polymer_deposition_probability_on_polymer: Mapping[str, float] = field(
        default_factory=dict)
    complex_activation_yield: LowEnergyActivationYield | None = None
    polymer_activation_yield: LowEnergyActivationYield | None = None
    activation_energetic_species: tuple[str, ...] = ()
    energetic_polymer_deposition_yield: LowEnergyActivationYield | None = None
    energetic_polymer_deposition_species: tuple[str, ...] = ()
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
        if (int(self.complex_removal_reaction_order) != self.complex_removal_reaction_order
                or int(self.complex_removal_reaction_order) not in {1, 2}):
            raise ValueError("complex removal reaction order must be one or two")
        object.__setattr__(
            self, "complex_removal_reaction_order", int(self.complex_removal_reaction_order))
        if (self.complex_activation_yield is not None
                and not isinstance(self.complex_activation_yield, LowEnergyActivationYield)):
            raise TypeError("complex_activation_yield must be LowEnergyActivationYield or None")
        if (self.polymer_activation_yield is not None
                and not isinstance(self.polymer_activation_yield, LowEnergyActivationYield)):
            raise TypeError("polymer_activation_yield must be LowEnergyActivationYield or None")
        if (self.energetic_polymer_deposition_yield is not None
                and not isinstance(
                    self.energetic_polymer_deposition_yield,
                    LowEnergyActivationYield)):
            raise TypeError(
                "energetic_polymer_deposition_yield must be "
                "LowEnergyActivationYield or None")
        if (self.complex_activation_yield is not None
                and self.complex_removal_reaction_order != 1):
            raise ValueError(
                "activated complex sites currently require first-order complex removal")
        maps = {}
        for name in (
                "complex_formation_probability",
                "polymer_deposition_probability_on_substrate",
                "polymer_deposition_probability_on_polymer",
                "activated_polymer_deposition_probability_on_substrate",
                "activated_polymer_deposition_probability_on_polymer"):
            values = dict(getattr(self, name))
            if any((not species or not np.isfinite(value) or value < 0.0 or value > 1.0)
                   for species, value in values.items()):
                raise ValueError(f"invalid probability map: {name}")
            maps[name] = MappingProxyType(values)
        for name, values in maps.items(): object.__setattr__(self, name, values)
        if (maps["activated_polymer_deposition_probability_on_substrate"]
                and self.complex_activation_yield is None):
            raise ValueError(
                "activated substrate deposition requires a complex activation yield")
        if (maps["activated_polymer_deposition_probability_on_polymer"]
                and self.polymer_activation_yield is None):
            raise ValueError(
                "activated polymer deposition requires a polymer activation yield")
        activation_species = tuple(str(name) for name in self.activation_energetic_species)
        if any(not name for name in activation_species) or len(set(activation_species)) != len(
                activation_species):
            raise ValueError("activation energetic species must be unique nonempty names")
        deposition_species = tuple(
            str(name) for name in self.energetic_polymer_deposition_species)
        if (any(not name for name in deposition_species)
                or len(set(deposition_species)) != len(deposition_species)
                or bool(self.energetic_polymer_deposition_yield)
                != bool(deposition_species)):
            raise ValueError(
                "energetic polymer deposition requires a nonempty unique species list "
                "and a yield law")
        object.__setattr__(self, "activation_energetic_species", activation_species)
        object.__setattr__(
            self, "energetic_polymer_deposition_species", deposition_species)
        evidence = dict(self.evidence)
        if any(not isinstance(item, ParameterEvidence) for item in evidence.values()):
            raise TypeError("parameter evidence values must be ParameterEvidence objects")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))
        object.__setattr__(self, "known_omissions", tuple(self.known_omissions))

    @classmethod
    def huang_kushner_2019_reduced_projection(
            cls, *, energetic_response_scale=1.0):
        """Project the published MCFPM Ar/C4F8/O2 table onto this reduced state.

        Huang et al. resolve activated sites and several oxide--fluorocarbon complexes.  This
        common-engine mechanism has one bounded complex fraction and one polymer inventory, so the
        constructor deliberately preserves only the reactions representable by those states:
        CF/CF2/CF3 and reactive C3F6 passivate oxide; every plotted fluorocarbon radical can grow
        polymer; ions remove bare oxide, complex, and polymer with the published threshold laws.

        ``energetic_response_scale`` is not a new yield shape.  It is the single experiment-adapter
        calibration allowed when the source reports electron density and self-bias but not the
        species-resolved ion/hot-neutral flux or IEAD.  A held-out validation must report the scale
        and may not retune it.  A value of one reproduces the published MCFPM p0 values exactly.
        """
        scale = float(energetic_response_scale)
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError("energetic_response_scale must be positive and finite")
        huang = "https://doi.org/10.1116/1.5090606"
        kaler = "https://doi.org/10.1088/1361-6463/aa6f40"

        def evidence(source, evidence_type, note="", *, supports=False):
            return ParameterEvidence(
                source, evidence_type, note=note,
                supports_prediction_within_declared_domain=supports)

        parameter_evidence = {
            "site_density_m2": evidence(
                kaler, "source_model_assumption",
                "CF2 uptake site density 1e15 cm^-2; transferred to the reduced state."),
            "bulk_formula_density_m3": evidence(
                "fused SiO2 density 2200 kg/m3 and molar mass 60.0843 g/mol",
                "material_constant_derived", supports=True),
            "polymer_monolayer_density_m2": evidence(
                kaler, "source_model_assumption",
                "CF2 uptake site density 1e15 cm^-2; transferred to the reduced state."),
            "complex_formation_probability": evidence(
                huang, "published_model_reduced_projection",
                "Table I p0 values; C3F6 uses the paper's reactive CxFy class."),
            "polymer_deposition_probability_on_substrate": evidence(
                huang, "published_model_reduced_projection",
                "Table I activated-complex deposition projected onto accessible substrate."),
            "polymer_deposition_probability_on_polymer": evidence(
                huang, "published_model_parameter", "Table I polymer-growth p0 values."),
            "activated_polymer_deposition_probability_on_substrate": evidence(
                huang, "published_model_parameter",
                "Table I polymer deposition on activated SiO2CmFn sites."),
            "activated_polymer_deposition_probability_on_polymer": evidence(
                huang, "published_model_parameter",
                "Table I P* sticking probabilities, ten times the unactivated P values."),
            "complex_activation_yield": evidence(
                huang, "published_model_parameter",
                "Table I and Eq. (2): p0=0.1, 5--70 eV low-energy window."),
            "polymer_activation_yield": evidence(
                huang, "published_model_parameter",
                "Table I and Eq. (2): p0=0.3, 5--30 eV low-energy window."),
            "activation_energetic_species": evidence(
                huang, "published_model_scope",
                "Empty tuple means every chemistry-facing ion/hot-neutral population; "
                "feature transport filters out non-bombarding species before this kernel."),
            "energetic_polymer_deposition_yield": evidence(
                huang, "published_model_parameter",
                "Table I p0=0.1, 5--70 eV for CFx+ and CxFy+ deposition on "
                "oxide-fluorocarbon complexes."),
            "energetic_polymer_deposition_species": evidence(
                huang, "published_model_scope",
                "Table I fluorocarbon positive-ion family; Ar+ and other non-FC ions "
                "do not use the direct polymer-deposition channel."),
            "oxygen_polymer_etch_probability": evidence(
                huang, "declared_absent_channel",
                "Jeong used Ar/C4F8 without an incident O channel; probability is zero."),
            "bare_sio2_yield": evidence(
                huang, "published_model_times_single_anchor_scale",
                f"Table I p0=0.9, Eth=70 eV, Er=140 eV; common scale={scale:.17g}."),
            "complex_sio2_yield": evidence(
                huang, "published_model_times_single_anchor_scale",
                f"Table I p0=0.75, Eth=35 eV, Er=140 eV; common scale={scale:.17g}."),
            "polymer_sputter_yield": evidence(
                huang, "published_model_times_single_anchor_scale",
                f"Table I p0=0.3, Eth=30 eV, Er=140 eV; common scale={scale:.17g}."),
        }
        polymer_probability = {
            "CF": 0.002, "CF2": 0.0015,
            "FC_complex_02": 0.001, "FC_polymer_heavy": 0.001,
        }
        activated_polymer_probability = {
            species: 10.0 * value for species, value in polymer_probability.items()}
        return cls(
            site_density_m2=1.0e19,
            bulk_formula_density_m3=2.205e28,
            polymer_monolayer_density_m2=1.0e19,
            complex_formation_probability={
                "CF": 0.4, "CF2": 0.3, "FC_complex_02": 0.2,
            },
            # Table I only deposits polymer on activated oxide-complex sites.  The previous
            # projection applied that rate to every accessible substrate site.
            polymer_deposition_probability_on_substrate={},
            polymer_deposition_probability_on_polymer=polymer_probability,
            oxygen_species="O", oxygen_polymer_etch_probability=0.0,
            bare_sio2_yield=EnergeticYield(
                0.9 * scale, 70.0, 140.0, energy_exponent=0.5,
                angular_model="kress_1999", angular_parameter=9.3),
            complex_sio2_yield=EnergeticYield(
                0.75 * scale, 35.0, 140.0, energy_exponent=0.5,
                angular_model="chang_sawin_1997"),
            polymer_sputter_yield=EnergeticYield(
                0.3 * scale, 30.0, 140.0, energy_exponent=0.5,
                angular_model="kress_1999", angular_parameter=9.3),
            activated_polymer_deposition_probability_on_substrate=polymer_probability,
            activated_polymer_deposition_probability_on_polymer=(
                activated_polymer_probability),
            complex_activation_yield=LowEnergyActivationYield(
                0.1, 5.0, 70.0,
                angular_model="kress_1999", angular_parameter=9.3),
            polymer_activation_yield=LowEnergyActivationYield(
                0.3, 5.0, 30.0,
                angular_model="kress_1999", angular_parameter=9.3),
            activation_energetic_species=(),
            energetic_polymer_deposition_yield=LowEnergyActivationYield(
                0.1, 5.0, 70.0,
                angular_model="kress_1999", angular_parameter=9.3),
            energetic_polymer_deposition_species=(
                "CF+", "CF2+", "CF3+", "C2F3+", "C2F4+", "C2F5+",
                "C3F5+", "C3F6+", "C3F7+", "C4F7+"),
            evidence=parameter_evidence,
            known_omissions=(
                "the published bare-SiO2 activated state is not yet resolved",
                "species-resolved oxide-fluorocarbon complexes are collapsed",
                "Jeong Figure 6 omits the atomic-F boundary flux",
                "species-resolved positive-ion and hot-neutral boundary fluxes are unreported",
                "etch-product redeposition is unresolved",
                "amorphous-carbon-mask chemistry is unresolved",
            ),
        )


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
        if not isinstance(parameters, ReducedSiO2FluorocarbonParameters):
            raise TypeError("parameters must be ReducedSiO2FluorocarbonParameters")
        self.parameters = parameters
        par = parameters

        def yield_manifest(law):
            if law is None:
                return None
            if isinstance(law, LowEnergyActivationYield):
                return {
                    "type": type(law).__name__,
                    "zero_energy_yield": float(law.zero_energy_yield),
                    "minimum_energy_eV": float(law.minimum_energy_eV),
                    "maximum_energy_eV": float(law.maximum_energy_eV),
                    "angular_model": law.angular_model,
                    "angular_parameter": (
                        None if law.angular_parameter is None
                        else float(law.angular_parameter)),
                }
            return {
                "type": type(law).__name__,
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
            "model": "reduced-sio2-fluorocarbon-common-engine-v1",
            "parameters": {
                "site_density_m2": float(par.site_density_m2),
                "bulk_formula_density_m3": float(par.bulk_formula_density_m3),
                "polymer_monolayer_density_m2": float(par.polymer_monolayer_density_m2),
                "complex_formation_probability": dict(
                    par.complex_formation_probability),
                "polymer_deposition_probability_on_substrate": dict(
                    par.polymer_deposition_probability_on_substrate),
                "polymer_deposition_probability_on_polymer": dict(
                    par.polymer_deposition_probability_on_polymer),
                "activated_polymer_deposition_probability_on_substrate": dict(
                    par.activated_polymer_deposition_probability_on_substrate),
                "activated_polymer_deposition_probability_on_polymer": dict(
                    par.activated_polymer_deposition_probability_on_polymer),
                "oxygen_species": par.oxygen_species,
                "oxygen_polymer_etch_probability": float(
                    par.oxygen_polymer_etch_probability),
                "complex_removal_reaction_order": int(
                    par.complex_removal_reaction_order),
                "bare_sio2_yield": yield_manifest(par.bare_sio2_yield),
                "complex_sio2_yield": yield_manifest(par.complex_sio2_yield),
                "polymer_sputter_yield": yield_manifest(par.polymer_sputter_yield),
                "complex_activation_yield": yield_manifest(
                    par.complex_activation_yield),
                "polymer_activation_yield": yield_manifest(
                    par.polymer_activation_yield),
                "activation_energetic_species": list(
                    par.activation_energetic_species),
                "energetic_polymer_deposition_yield": yield_manifest(
                    par.energetic_polymer_deposition_yield),
                "energetic_polymer_deposition_species": list(
                    par.energetic_polymer_deposition_species),
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

    @staticmethod
    def initial_state(shape=()):
        return SiO2SurfaceState.bare(shape)

    @staticmethod
    def _state_weighted_probability(base_probability, activated_probability, fraction):
        """Interpolate an unactivated/activated probability at a bounded site fraction."""
        base = dict(base_probability)
        activated = dict(activated_probability)
        return {
            species: (
                base.get(species, 0.0)
                + (activated.get(species, base.get(species, 0.0))
                   - base.get(species, 0.0)) * fraction)
            for species in set(base) | set(activated)
        }

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
        substrate_probability = self._state_weighted_probability(
            par.polymer_deposition_probability_on_substrate,
            par.activated_polymer_deposition_probability_on_substrate,
            state.activated_complex_fraction)
        polymer_probability = self._state_weighted_probability(
            par.polymer_deposition_probability_on_polymer,
            par.activated_polymer_deposition_probability_on_polymer,
            state.activated_polymer_fraction)
        all_deposition_species = set(substrate_probability) | set(polymer_probability)
        for species in all_deposition_species:
            add(species,
                substrate_probability.get(species, 0.0) * access
                + polymer_probability.get(species, 0.0) * polymer_coverage)
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
                     | set(par.activated_polymer_deposition_probability_on_substrate)
                     | set(par.activated_polymer_deposition_probability_on_polymer)
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
        if par.complex_activation_yield is not None:
            required_evidence.update({
                "activated_polymer_deposition_probability_on_substrate",
                "complex_activation_yield", "activation_energetic_species"})
        if par.polymer_activation_yield is not None:
            required_evidence.update({
                "activated_polymer_deposition_probability_on_polymer",
                "polymer_activation_yield", "activation_energetic_species"})
        if par.energetic_polymer_deposition_yield is not None:
            required_evidence.update({
                "energetic_polymer_deposition_yield",
                "energetic_polymer_deposition_species"})
        if par.complex_removal_reaction_order != 1:
            required_evidence.add("complex_removal_reaction_order")
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

    def _activation_rate(self, fluxes, yield_law, shape):
        if yield_law is None:
            return np.zeros(shape)
        selected_species = set(self.parameters.activation_energetic_species)
        total = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if selected_species and population.name not in selected_species:
                continue
            total = total + self._broadcast(
                population.yield_rate_m2_s(yield_law), shape)
        return total

    def _energetic_polymer_deposition_rate(self, fluxes, shape):
        law = self.parameters.energetic_polymer_deposition_yield
        if law is None:
            return np.zeros(shape)
        selected_species = set(
            self.parameters.energetic_polymer_deposition_species)
        total = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name not in selected_species:
                continue
            total = total + self._broadcast(
                population.yield_rate_m2_s(law), shape)
        return total

    def _activation_step(
            self, activated_complex, activated_polymer, complex_fraction,
            polymer_inventory, fluxes, duration_s, shape):
        """Exact bounded activation update at fixed complex/polymer inventories."""
        if duration_s == 0.0:
            return (
                np.array(activated_complex, copy=True),
                np.array(activated_polymer, copy=True))
        par = self.parameters
        complex_hazard = (
            self._activation_rate(fluxes, par.complex_activation_yield, shape)
            / par.site_density_m2)
        polymer_hazard = (
            self._activation_rate(fluxes, par.polymer_activation_yield, shape)
            / par.polymer_monolayer_density_m2)
        complex_updated = (
            complex_fraction
            - (complex_fraction - activated_complex)
            * np.exp(-complex_hazard * duration_s))
        polymer_present = polymer_inventory > 0.0
        polymer_updated = np.where(
            polymer_present,
            1.0 - (1.0 - activated_polymer)
            * np.exp(-polymer_hazard * duration_s),
            0.0)
        # Preserve exact identity where the corresponding activation channel is absent.
        complex_updated = np.where(
            complex_hazard > 0.0, complex_updated, activated_complex)
        polymer_updated = np.where(
            polymer_present & (polymer_hazard > 0.0),
            polymer_updated, np.where(polymer_present, activated_polymer, 0.0))
        tolerance = 32.0 * np.finfo(float).eps
        if (np.any(complex_updated < -tolerance)
                or np.any(complex_updated > complex_fraction + tolerance)
                or np.any(polymer_updated < -tolerance)
                or np.any(polymer_updated > 1.0 + tolerance)):
            raise RuntimeError("low-energy activation update violated boundedness")
        return (
            np.minimum(np.maximum(complex_updated, 0.0), complex_fraction),
            np.minimum(np.maximum(polymer_updated, 0.0), 1.0))

    def _polymer_step(
            self, inventory, complex_fraction, activated_complex, activated_polymer,
            fluxes, duration_s, shape):
        """Exact constant-flux solution and exact integrated deposition/removal bookkeeping."""
        par = self.parameters; monolayer = par.polymer_monolayer_density_m2
        substrate_probability = self._state_weighted_probability(
            par.polymer_deposition_probability_on_substrate,
            par.activated_polymer_deposition_probability_on_substrate,
            activated_complex)
        polymer_probability = self._state_weighted_probability(
            par.polymer_deposition_probability_on_polymer,
            par.activated_polymer_deposition_probability_on_polymer,
            activated_polymer)
        deposit_substrate = self._neutral_weighted_rate(
            fluxes, substrate_probability, shape)
        # Huang--Kushner Table I also permits a low-energy fluorocarbon ion
        # to deposit one polymer unit directly on an oxide-fluorocarbon
        # complex.  This is distinct from activation: it consumes no thermal
        # neutral and therefore must enter the conservative polymer source
        # explicitly.  ``complex_fraction`` is the accessible substrate-site
        # fraction; the analytic polymer ODE applies the overlayer coverage.
        deposit_substrate = (
            deposit_substrate
            + self._energetic_polymer_deposition_rate(fluxes, shape)
            * complex_fraction)
        deposit_polymer = self._neutral_weighted_rate(
            fluxes, polymer_probability, shape)
        oxygen_flux = self._broadcast(
            fluxes.neutral_flux_m2_s.get(par.oxygen_species, 0.0), shape)
        removal_capacity = (oxygen_flux * par.oxygen_polymer_etch_probability
                            + self._energetic_rate(fluxes, par.polymer_sputter_yield, shape))
        active_reaction = deposit_substrate + deposit_polymer + removal_capacity > 0.0
        # N=0 with no substrate-nucleation channel is an exact invariant.  Polymer-growth and
        # removal rates are multiplied by the polymer coverage, which is identically zero there;
        # evaluating the transformed analytic solution can otherwise manufacture an O(ULP)
        # inventory when its large growth/removal coefficients nearly cancel.
        active_reaction &= (deposit_substrate > 0.0) | (inventory > 0.0)
        # A float64 inventory cannot represent a change below O(eps*N_mono).  On nearly dark
        # faces, evaluating the exact log solution below that floor can turn a sub-unit physical
        # change into an O(10--100) cancellation artifact and then fail its own ledger.  Price that
        # representational limit explicitly: if even the sum of all incident reaction capacities
        # over this substep is below the same 32-eps surface-unit floor used by the complex-state
        # integrator, the authoritative state and every associated ledger remain exactly unchanged.
        roundoff_floor = (
            32.0 * np.finfo(float).eps * monolayer
            * np.maximum(np.abs(np.asarray(inventory, dtype=float)) / monolayer, 1.0))
        throughput_bound = duration_s * (
            deposit_substrate + deposit_polymer + removal_capacity)
        active_reaction &= throughput_bound > roundoff_floor
        if not np.any(active_reaction):
            return np.array(inventory, copy=True), np.zeros(shape), np.zeros(shape)
        deposit_substrate = np.where(active_reaction, deposit_substrate, 0.0)
        deposit_polymer = np.where(active_reaction, deposit_polymer, 0.0)
        removal_capacity = np.where(active_reaction, removal_capacity, 0.0)

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
                    # y1/y0 = exp(x)*(1-a/y0) + a/y0.  Factoring exp(x) gives
                    # x + log1p((a/y0)*expm1(-x)).  This preserves the small positive
                    # correction to x when a thick film is removed slowly.  A direct logaddexp
                    # can round to the wrong side of x at N/M~30, which then makes the implied
                    # coverage exceed the step duration and creates a false ledger failure.
                    local_x = selected_exponent[~growing][nonnegative]
                    a_over_y0 = (
                        local_a[nonnegative] * np.exp(-local_log_y0[nonnegative]))
                    mild_decay = local_x > -1.0
                    local_negative_delta = np.empty(local_x.shape)
                    if np.any(mild_decay):
                        local_negative_delta[mild_decay] = (
                            local_x[mild_decay]
                            + np.log1p(
                                a_over_y0[mild_decay]
                                * np.expm1(-local_x[mild_decay])))
                    if np.any(~mild_decay):
                        # For a large negative exponent the factored expm1 can overflow and its
                        # logarithm nearly cancels x.  The original log-sum representation is
                        # well-conditioned in that floor-reaching regime.
                        strong_difference = local_difference[nonnegative][~mild_decay]
                        log_first = np.where(
                            strong_difference > 0.0,
                            local_x[~mild_decay]
                            + np.log(np.maximum(
                                strong_difference, np.finfo(float).tiny)),
                            -np.inf)
                        log_second = (
                            np.log(local_a[nonnegative][~mild_decay])
                            - local_log_y0[nonnegative][~mild_decay])
                        local_negative_delta[~mild_decay] = np.logaddexp(
                            log_first, log_second)
                    negative_delta[nonnegative] = local_negative_delta
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
        if par.complex_removal_reaction_order == 2:
            return self._quadratic_substrate_step(
                complex_fraction, formation_event_rate, complex_removal_event_rate,
                bare_removal_event_rate, formation_hazard, removal_hazard,
                duration_s, shape)
        total_hazard = formation_hazard + removal_hazard
        equilibrium = np.divide(
            formation_hazard, total_hazard, out=np.zeros(shape), where=total_hazard > 0.0)
        bare_equilibrium = np.divide(
            removal_hazard, total_hazard, out=np.zeros(shape), where=total_hazard > 0.0)
        decay = np.exp(-total_hazard * duration_s)
        active = total_hazard > 0.0
        # Evolve both complementary coverages analytically.  Near theta=1, forming the
        # bare-site exposure as ``duration - integral(theta)`` catastrophically cancels;
        # near theta=0 the converse is true.  Direct theta/q solutions and integrals retain
        # the small physical channel on both limits while describing the identical ODE.
        bare_fraction = 1.0 - complex_fraction
        updated_complex = equilibrium + (complex_fraction - equilibrium) * decay
        updated_bare = bare_equilibrium + (bare_fraction - bare_equilibrium) * decay
        updated = np.where(updated_complex >= 0.5, 1.0 - updated_bare, updated_complex)
        updated = np.where(active, updated, complex_fraction)
        integral_kernel = np.empty(shape)
        integral_kernel[active] = (
            -np.expm1(-total_hazard[active] * duration_s) / total_hazard[active])
        integral_kernel[~active] = duration_s
        integral_complex = np.empty(shape)
        integral_bare = np.empty(shape)
        integral_complex[active] = (
            equilibrium[active] * duration_s
            + (complex_fraction[active] - equilibrium[active]) * integral_kernel[active])
        integral_bare[active] = (
            bare_equilibrium[active] * duration_s
            + (bare_fraction[active] - bare_equilibrium[active]) * integral_kernel[active])
        integral_complex[~active] = complex_fraction[~active] * duration_s
        integral_bare[~active] = bare_fraction[~active] * duration_s
        removed_complex = complex_removal_event_rate * integral_complex
        formed_from_rate = formation_event_rate * integral_bare
        removed_bare = bare_removal_event_rate * integral_bare
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

    def _quadratic_substrate_step(
            self, complex_fraction, formation_event_rate, complex_removal_event_rate,
            bare_removal_event_rate, formation_hazard, removal_hazard, duration_s, shape):
        """Exact constant-flux update for ``dtheta/dt=a(1-theta)-b theta**2``.

        The positive and negative roots of the Riccati equation give a closed-form state update.
        Its analytic time integral supplies both bare-site exposure and the formation ledger; the
        complex-removal ledger follows from the conserved site equation.  Thus the nonlinear
        activation changes the physical rate law without weakening any bookkeeping identity.
        """
        theta0 = np.asarray(complex_fraction, dtype=float)
        a = np.asarray(formation_hazard, dtype=float)
        b = np.asarray(removal_hazard, dtype=float)
        updated = np.array(theta0, copy=True)
        integral_theta = np.array(theta0 * duration_s, copy=True)

        formation_only = (a > 0.0) & (b == 0.0)
        if np.any(formation_only):
            local_a = a[formation_only]
            change = ((1.0 - theta0[formation_only])
                      * (-np.expm1(-local_a * duration_s)))
            updated[formation_only] = theta0[formation_only] + change
            integral_theta[formation_only] = duration_s - change / local_a

        removal_only = (a == 0.0) & (b > 0.0)
        if np.any(removal_only):
            local_b = b[removal_only]
            scaled = local_b * theta0[removal_only] * duration_s
            updated[removal_only] = theta0[removal_only] / (1.0 + scaled)
            integral_theta[removal_only] = np.log1p(scaled) / local_b

        coupled = (a > 0.0) & (b > 0.0)
        if np.any(coupled):
            local_a = a[coupled]
            local_b = b[coupled]
            local_theta0 = theta0[coupled]
            discriminant = np.sqrt(local_a * local_a + 4.0 * local_a * local_b)
            positive_root = 2.0 * local_a / (local_a + discriminant)
            negative_root = -local_a / (local_b * positive_root)
            root_ratio = ((local_theta0 - positive_root)
                          / (local_theta0 - negative_root))
            evolved_ratio = root_ratio * np.exp(-discriminant * duration_s)
            local_updated = ((positive_root - evolved_ratio * negative_root)
                             / (1.0 - evolved_ratio))
            updated[coupled] = local_updated
            log_ratio = np.log1p(
                (local_theta0 - local_updated) / (local_updated - negative_root))
            integral_theta[coupled] = (
                positive_root * duration_s + log_ratio / local_b)

        state_tolerance = 32.0 * np.finfo(float).eps
        if (np.any(updated < -state_tolerance) or np.any(updated > 1.0 + state_tolerance)
                or np.any(integral_theta < -state_tolerance * max(duration_s, 1.0))
                or np.any(integral_theta > duration_s * (1.0 + state_tolerance))):
            raise RuntimeError("quadratic oxide-complex integrator violated boundedness")
        # Only erase representational excursions at an exact physical boundary.
        updated = np.where(updated < 0.0, 0.0, np.where(updated > 1.0, 1.0, updated))
        integral_theta = np.where(
            integral_theta < 0.0, 0.0,
            np.where(integral_theta > duration_s, duration_s, integral_theta))

        formed_complex = formation_event_rate * (duration_s - integral_theta)
        site_change = (updated - theta0) * self.parameters.site_density_m2
        removed_complex = formed_complex - site_change
        removed_bare = bare_removal_event_rate * (duration_s - integral_theta)
        scale = np.maximum.reduce((
            np.abs(formed_complex), np.abs(removed_complex), np.abs(site_change),
            np.ones(shape)))
        tolerance = np.maximum(
            5e-12 * scale,
            32.0 * np.finfo(float).eps * self.parameters.site_density_m2)
        if (np.any(formed_complex < -tolerance)
                or np.any(removed_complex < -tolerance)
                or np.any(np.abs(site_change - (formed_complex - removed_complex)) > tolerance)):
            raise RuntimeError("quadratic oxide-complex bookkeeping failed conservation")
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
        activated_complex = np.array(state.activated_complex_fraction, copy=True)
        activated_polymer = np.array(state.activated_polymer_fraction, copy=True)
        formed_complex = np.zeros(shape); removed_complex = np.zeros(shape)
        removed_bare = np.zeros(shape); deposited_polymer = np.zeros(shape)
        removed_polymer = np.zeros(shape)
        for _ in range(n_steps):
            activated_complex, activated_polymer = self._activation_step(
                activated_complex, activated_polymer, complex_fraction,
                polymer, fluxes, 0.5 * step, shape)
            polymer, deposited, removed = self._polymer_step(
                polymer, complex_fraction, activated_complex, activated_polymer,
                fluxes, 0.5 * step, shape)
            deposited_polymer += deposited; removed_polymer += removed
            activated_polymer *= np.exp(
                -removed / self.parameters.polymer_monolayer_density_m2)
            activated_polymer = np.where(polymer > 0.0, activated_polymer, 0.0)
            access = np.exp(
                -polymer / self.parameters.polymer_monolayer_density_m2)
            complex_removal_hazard = (
                self._energetic_rate(
                    fluxes, self.parameters.complex_sio2_yield, shape)
                * access / self.parameters.site_density_m2)
            complex_fraction, formed, removed_c, removed_b = self._substrate_step(
                complex_fraction, polymer, fluxes, step, shape)
            activated_complex *= np.exp(-complex_removal_hazard * step)
            activated_complex = np.minimum(activated_complex, complex_fraction)
            formed_complex += formed; removed_complex += removed_c; removed_bare += removed_b
            removed_total += removed_c + removed_b
            polymer, deposited, removed = self._polymer_step(
                polymer, complex_fraction, activated_complex, activated_polymer,
                fluxes, 0.5 * step, shape)
            deposited_polymer += deposited; removed_polymer += removed
            activated_polymer *= np.exp(
                -removed / self.parameters.polymer_monolayer_density_m2)
            activated_polymer = np.where(polymer > 0.0, activated_polymer, 0.0)
            activated_complex, activated_polymer = self._activation_step(
                activated_complex, activated_polymer, complex_fraction,
                polymer, fluxes, 0.5 * step, shape)
        new_state = SiO2SurfaceState(
            complex_fraction, polymer, removed_total,
            activated_complex, activated_polymer)
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
