"""Atom-balanced Gray Ar+/F beam closure for silicon.

Gray's 1993 MIT thesis, Chapter 5, reduces the normally incident Ar+/atomic-F
beam measurements to a two-site steady balance.  For one ion energy,

``theta_F = s0 R / (s0 R + 2 beta2 (1 + b))``

``Y = p0 (1 - theta_F) + beta2 (1 + b) theta_F``

where ``R`` is the F/Ar+ particle-flux ratio.  The source supplies independent
square-root energy laws for bare-Si sputtering and the fluorinated-surface
coefficient, plus a mass-spectrometric fragment branch:

``p0 = 0.0337 max(sqrt(E) - sqrt(20 eV), 0)``

``beta2 = 0.687 max(sqrt(E) - sqrt(4 eV), 0)``

``b = 0.009 sqrt(E/eV)``.

This implementation integrates the joint ion-energy measure before solving the
site balance, routes removed Si into Si, SiF2, and SiF4 populations, and checks
the required F-atom inventory.  It deliberately does not invent an off-normal
angular law or transfer Ar+ coefficients to SFx+ ions.

The closure is beam-regressed and phenomenological, not an elementary
first-principles potential.  It is nevertheless substantially more fundamental
than a profile-depth scale: species, energy, surface occupation, product
identity, and elemental bookkeeping are explicit, and no feature datum enters.
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
    SteinbruchelYield,
    SurfaceFluxes,
)


@dataclass(frozen=True)
class GrayArFSiState:
    """Quasi-steady fluorination plus cumulative elemental inventories."""

    fluorinated_fraction: np.ndarray | float = 0.0
    removed_si_atoms_m2: np.ndarray | float = 0.0
    consumed_f_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        supplied = np.broadcast_arrays(
            np.asarray(self.fluorinated_fraction, dtype=float),
            np.asarray(self.removed_si_atoms_m2, dtype=float),
            np.asarray(self.consumed_f_atoms_m2, dtype=float),
        )
        fluorinated, removed, consumed = [
            np.array(item, copy=True) for item in supplied
        ]
        if (np.any(~np.isfinite(fluorinated))
                or np.any((fluorinated < 0.0) | (fluorinated > 1.0))
                or np.any(~np.isfinite(removed)) or np.any(removed < 0.0)
                or np.any(~np.isfinite(consumed)) or np.any(consumed < 0.0)):
            raise ValueError("invalid Gray Ar+/F silicon state")
        for value in (fluorinated, removed, consumed):
            value.setflags(write=False)
        object.__setattr__(self, "fluorinated_fraction", fluorinated)
        object.__setattr__(self, "removed_si_atoms_m2", removed)
        object.__setattr__(self, "consumed_f_atoms_m2", consumed)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape)
        return cls(zero, zero, zero)

    def conservative_surface_fields(self):
        return {
            "fluorinated_fraction": self.fluorinated_fraction,
            "removed_si_atoms_m2": self.removed_si_atoms_m2,
            "consumed_f_atoms_m2": self.consumed_f_atoms_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {
            "fluorinated_fraction": 1.0,
            "removed_si_atoms_m2": None,
            "consumed_f_atoms_m2": None,
        }

    def surface_field_remap_modes(self):
        return {
            "fluorinated_fraction": "intensive",
            "removed_si_atoms_m2": "conservative",
            "consumed_f_atoms_m2": "conservative",
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("Gray Ar+/F silicon remap fields do not match")
        return type(self)(
            fields["fluorinated_fraction"],
            fields["removed_si_atoms_m2"],
            fields["consumed_f_atoms_m2"],
        )


@dataclass(frozen=True)
class GrayArFSiParameters:
    """Source-fixed parameters and the declared beam domain."""

    ion_species: str
    neutral_species: str
    initial_f_sticking_coefficient: float
    physical_sputter_yield: SteinbruchelYield
    fluorinated_surface_yield: SteinbruchelYield
    fragment_branch_per_sqrt_eV: float
    bulk_si_atom_density_m3: float
    minimum_measured_energy_eV: float
    maximum_measured_energy_eV: float
    minimum_regression_flux_ratio: float
    normal_incidence_tolerance_deg: float
    evidence: Mapping[str, ParameterEvidence]

    def __post_init__(self):
        evidence = dict(self.evidence)
        required = {
            "initial_f_sticking_coefficient",
            "physical_sputter_yield",
            "fluorinated_surface_yield",
            "fragment_branch",
            "bulk_si_atom_density_m3",
        }
        if (not self.ion_species or not self.neutral_species
                or self.ion_species == self.neutral_species
                or not np.isfinite(self.initial_f_sticking_coefficient)
                or not 0.0 <= self.initial_f_sticking_coefficient <= 1.0
                or not isinstance(
                    self.physical_sputter_yield, SteinbruchelYield)
                or not isinstance(
                    self.fluorinated_surface_yield, SteinbruchelYield)
                or not np.isfinite(self.fragment_branch_per_sqrt_eV)
                or self.fragment_branch_per_sqrt_eV < 0.0
                or not np.isfinite(self.bulk_si_atom_density_m3)
                or self.bulk_si_atom_density_m3 <= 0.0
                or not np.isfinite(self.minimum_measured_energy_eV)
                or self.minimum_measured_energy_eV <= 0.0
                or not np.isfinite(self.maximum_measured_energy_eV)
                or self.maximum_measured_energy_eV
                <= self.minimum_measured_energy_eV
                or not np.isfinite(self.minimum_regression_flux_ratio)
                or self.minimum_regression_flux_ratio <= 0.0
                or not np.isfinite(self.normal_incidence_tolerance_deg)
                or not 0.0 <= self.normal_incidence_tolerance_deg < 90.0
                or set(evidence) != required
                or any(not isinstance(item, ParameterEvidence)
                       for item in evidence.values())):
            raise ValueError("invalid Gray Ar+/F silicon parameters")
        if (self.physical_sputter_yield.angular_model != "none"
                or self.fluorinated_surface_yield.angular_model != "none"):
            raise ValueError(
                "Gray normal-incidence cards may not hide an angular model")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @classmethod
    def beam_energy_law(cls):
        thesis = (
            "gray-1993-thesis, Chapter 5: Tables 5-1 and 5-9; "
            "Eqs. 5-27, 5-30, 5-31, and 5-34"
        )
        return cls(
            ion_species="Ar+",
            neutral_species="F",
            initial_f_sticking_coefficient=0.2,
            physical_sputter_yield=SteinbruchelYield(0.0337, 20.0),
            fluorinated_surface_yield=SteinbruchelYield(0.687, 4.0),
            fragment_branch_per_sqrt_eV=0.009,
            bulk_si_atom_density_m3=5.0e28,
            minimum_measured_energy_eV=20.0,
            maximum_measured_energy_eV=1000.0,
            minimum_regression_flux_ratio=15.0,
            normal_incidence_tolerance_deg=1.0e-8,
            evidence={
                "initial_f_sticking_coefficient": ParameterEvidence(
                    thesis,
                    "controlled-beam yield regression",
                    note=(
                        "Table 5-9 constant-s0 fit; apparent initial F "
                        "sticking on ion-renewed clean Si sites"),
                    supports_prediction_within_declared_domain=True,
                ),
                "physical_sputter_yield": ParameterEvidence(
                    thesis,
                    "controlled-beam and cross-laboratory sputter regression",
                    note="Table 5-1 and Eq. 5-7: A=0.0337, Eth about 20 eV",
                    supports_prediction_within_declared_domain=True,
                ),
                "fluorinated_surface_yield": ParameterEvidence(
                    thesis,
                    "controlled-beam yield regression",
                    note=(
                        "Eq. 5-34: beta2=0.687(sqrt(E)-sqrt(4)); "
                        "Table 5-9 spans 20--1000 eV"),
                    supports_prediction_within_declared_domain=True,
                ),
                "fragment_branch": ParameterEvidence(
                    thesis,
                    "modulated-beam mass-spectrometric regression",
                    note=(
                        "Eq. 5-27: summed unsaturated SiFx/SiF4 branch "
                        "represented by SiF2 in the source's reduced model"),
                    supports_prediction_within_declared_domain=True,
                ),
                "bulk_si_atom_density_m3": ParameterEvidence(
                    "crystalline-Si density and molar mass; rounded",
                    "derived physical constant",
                    supports_prediction_within_declared_domain=True,
                ),
            },
        )


@dataclass(frozen=True)
class GrayArFSiStepResult:
    state: GrayArFSiState
    etch_velocity_m_s: np.ndarray
    fluorinated_fraction: np.ndarray
    physical_removed_si_atoms_m2: np.ndarray
    sif2_removed_si_atoms_m2: np.ndarray
    sif4_removed_si_atoms_m2: np.ndarray
    consumed_f_atoms_m2: np.ndarray
    incident_f_atoms_m2: np.ndarray
    steady_site_balance_abs_residual_f_atoms_m2_s: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity

    def __post_init__(self):
        names = (
            "etch_velocity_m_s",
            "fluorinated_fraction",
            "physical_removed_si_atoms_m2",
            "sif2_removed_si_atoms_m2",
            "sif4_removed_si_atoms_m2",
            "consumed_f_atoms_m2",
            "incident_f_atoms_m2",
            "steady_site_balance_abs_residual_f_atoms_m2_s",
        )
        values = {}
        for name in names:
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Gray Ar+/F silicon step result")
            value.setflags(write=False)
            values[name] = value
            object.__setattr__(self, name, value)
        if (np.any(values["fluorinated_fraction"] > 1.0)
                or np.any(
                    values["consumed_f_atoms_m2"]
                    > values["incident_f_atoms_m2"]
                    + 64.0 * np.finfo(float).eps
                    * np.maximum(values["incident_f_atoms_m2"], 1.0))
                or not isinstance(self.state, GrayArFSiState)
                or not isinstance(self.material_exchange, SurfaceMaterialExchange)
                or not isinstance(self.validity, MechanismValidity)):
            raise ValueError("invalid Gray Ar+/F silicon result contract")
        products = validate_surface_product_routing(
            self.material_exchange, tuple(self.product_populations))
        object.__setattr__(self, "product_populations", products)
        removed = (
            values["physical_removed_si_atoms_m2"]
            + values["sif2_removed_si_atoms_m2"]
            + values["sif4_removed_si_atoms_m2"]
        )
        expected = self.material_exchange.removed_units_m2["Si_atom"]
        if np.any(
            np.abs(expected - removed)
            > 64.0 * np.finfo(float).eps * np.maximum(removed, 1.0)
        ):
            raise ValueError("Gray Ar+/F silicon branches do not close")


class GrayArFSiMechanism:
    """Normally incident, species-resolved Ar+/F/Si beam closure."""

    quasi_steady_surface_state = True

    def __init__(self, parameters: GrayArFSiParameters | None = None):
        self.parameters = (
            GrayArFSiParameters.beam_energy_law()
            if parameters is None else parameters
        )
        if not isinstance(self.parameters, GrayArFSiParameters):
            raise TypeError("parameters must be GrayArFSiParameters")

    @staticmethod
    def initial_state(shape=()):
        return GrayArFSiState.bare(shape)

    @property
    def provenance(self):
        par = self.parameters
        return MappingProxyType({
            "model": "Gray-1993 normally-incident Ar+/F/Si site model",
            "equations": [
                "5-7", "5-27", "5-30", "5-31", "5-34", "5-36",
            ],
            "parameters": {
                "initial_f_sticking_coefficient":
                    par.initial_f_sticking_coefficient,
                "physical_sputter_yield": {
                    "prefactor_per_sqrt_eV":
                        par.physical_sputter_yield.prefactor_per_sqrt_eV,
                    "threshold_energy_eV":
                        par.physical_sputter_yield.threshold_energy_eV,
                },
                "fluorinated_surface_yield": {
                    "prefactor_per_sqrt_eV":
                        par.fluorinated_surface_yield.prefactor_per_sqrt_eV,
                    "threshold_energy_eV":
                        par.fluorinated_surface_yield.threshold_energy_eV,
                },
                "fragment_branch_per_sqrt_eV":
                    par.fragment_branch_per_sqrt_eV,
                "bulk_si_atom_density_m3": par.bulk_si_atom_density_m3,
            },
            "declared_domain": {
                "ion_species": par.ion_species,
                "neutral_species": par.neutral_species,
                "energy_eV": [
                    par.minimum_measured_energy_eV,
                    par.maximum_measured_energy_eV,
                ],
                "incidence": "normal only",
                "fit_used_f_to_ar_ratio_greater_than":
                    par.minimum_regression_flux_ratio,
            },
            "evidence": {
                name: {
                    "source": item.source,
                    "evidence_type": item.evidence_type,
                    "note": item.note,
                    "supports_prediction_within_declared_domain":
                        item.supports_prediction_within_declared_domain,
                }
                for name, item in par.evidence.items()
            },
            "claim": (
                "beam-regressed species/energy/site/product closure; "
                "not an elementary first-principles potential and not a "
                "feature-depth calibration"
            ),
        })

    def neutral_reaction_probability(self, state: GrayArFSiState):
        if not isinstance(state, GrayArFSiState):
            raise TypeError("Gray neutral probabilities require GrayArFSiState")
        return MappingProxyType({
            self.parameters.neutral_species: (
                self.parameters.initial_f_sticking_coefficient
                * (1.0 - state.fluorinated_fraction)
            ),
        })

    @staticmethod
    def _positive_energy_angle(population):
        if isinstance(population, FaceResolvedEnergeticFlux):
            selected = population.event_flux_m2_s > 0.0
            return (
                population.event_energy_eV[selected],
                population.event_cosine_incidence[selected],
            )
        if isinstance(population, EnergeticFlux):
            if not np.any(np.asarray(population.flux_m2_s) > 0.0):
                return np.empty(0), np.empty(0)
            selected = population.weight > 0.0
            return (
                population.energy_eV[selected],
                population.cosine_incidence[selected],
            )
        raise TypeError(type(population).__name__)  # pragma: no cover

    @staticmethod
    def _reduce_population(population, value_per_event, shape):
        if isinstance(population, FaceResolvedEnergeticFlux):
            reduced = np.bincount(
                population.event_face,
                weights=population.event_flux_m2_s * value_per_event(
                    population.event_energy_eV),
                minlength=population.face_count,
            )
            return np.broadcast_to(reduced, shape)
        if isinstance(population, EnergeticFlux):
            mean = float(np.dot(
                population.weight,
                value_per_event(population.energy_eV),
            ))
            return np.broadcast_to(
                np.asarray(population.flux_m2_s, dtype=float), shape
            ) * mean
        raise TypeError(type(population).__name__)  # pragma: no cover

    @staticmethod
    def _reduce_flux(population, shape):
        return np.broadcast_to(
            np.asarray(population.flux_m2_s, dtype=float), shape)

    def _incident_measures(self, fluxes, shape):
        par = self.parameters
        ion_flux = np.zeros(shape)
        physical = np.zeros(shape)
        sif4_capacity = np.zeros(shape)
        sif2_capacity = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name != par.ion_species:
                continue
            ion_flux = ion_flux + self._reduce_flux(population, shape)
            physical = physical + self._reduce_population(
                population,
                lambda energy: par.physical_sputter_yield.evaluate(
                    energy, np.ones(np.asarray(energy).shape)),
                shape,
            )
            sif4_capacity = sif4_capacity + self._reduce_population(
                population,
                lambda energy: par.fluorinated_surface_yield.evaluate(
                    energy, np.ones(np.asarray(energy).shape)),
                shape,
            )
            sif2_capacity = sif2_capacity + self._reduce_population(
                population,
                lambda energy: (
                    par.fluorinated_surface_yield.evaluate(
                        energy, np.ones(np.asarray(energy).shape))
                    * par.fragment_branch_per_sqrt_eV * np.sqrt(energy)
                ),
                shape,
            )
        return ion_flux, physical, sif2_capacity, sif4_capacity

    def validity(self, state, fluxes):
        if not isinstance(state, GrayArFSiState):
            raise TypeError("Gray validity requires GrayArFSiState")
        par = self.parameters
        shape = state.fluorinated_fraction.shape
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name != par.neutral_species
            and np.any(np.asarray(value) > 0.0)
        ))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != par.ion_species
            and np.any(np.asarray(item.flux_m2_s) > 0.0)
        }))
        leaves_energy = False
        off_normal = False
        for population in fluxes.energetic_fluxes:
            if population.name != par.ion_species:
                continue
            energy, cosine = self._positive_energy_angle(population)
            leaves_energy |= bool(np.any(
                (energy < par.minimum_measured_energy_eV)
                | (energy > par.maximum_measured_energy_eV)
            ))
            if cosine.size:
                angle = np.rad2deg(np.arccos(cosine))
                off_normal |= bool(np.any(
                    angle > par.normal_incidence_tolerance_deg))
        f_flux = np.broadcast_to(np.asarray(
            fluxes.neutral_flux_m2_s.get(par.neutral_species, 0.0),
            dtype=float), shape)
        ion_flux, _, _, _ = self._incident_measures(fluxes, shape)
        ratio = np.zeros(shape)
        active = ion_flux > 0.0
        ratio[active] = f_flux[active] / ion_flux[active]
        below_fit_ratio = bool(np.any(
            active & (ratio > 0.0)
            & (ratio <= par.minimum_regression_flux_ratio)))
        thermal_only = bool(np.any((~active) & (f_flux > 0.0)))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append(
                "positive incident flux has no declared Ar+/F/Si channel")
        if leaves_energy:
            reasons.append("ion energy leaves the measured 20--1000 eV board")
        if off_normal:
            reasons.append(
                "Gray measured normal incidence; no angular law is declared")
        if below_fit_ratio:
            reasons.append(
                "positive F/Ar+ ratio is at or below the R>15 regression board")
        if thermal_only:
            reasons.append(
                "F-only thermal silicon etching is not included in this beam card")
        nonpredictive = tuple(sorted(
            name for name, item in par.evidence.items()
            if not item.supports_prediction_within_declared_domain
        ))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "Gray's reduced CSTR/site model is phenomenological rather than elementary",
                "thermal F-only etching is omitted from this ion-enhanced card",
                "the summed unsaturated SiFx signal is represented as SiF2",
                "no product launch energy-angle distribution was measured",
                "no off-normal ion response was measured",
                "the Ar+ card is not transferred to reactive SFx+ ions",
                "transient site density and collision-cascade depth are unresolved",
            ),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive,
        )

    def advance(
            self, state: GrayArFSiState, fluxes: SurfaceFluxes,
            duration_s: float, *, strict=True):
        if not isinstance(state, GrayArFSiState):
            raise TypeError("Gray mechanism requires GrayArFSiState")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(state, fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))
        par = self.parameters
        shape = state.fluorinated_fraction.shape
        f_flux = np.broadcast_to(np.asarray(
            fluxes.neutral_flux_m2_s.get(par.neutral_species, 0.0),
            dtype=float), shape)
        _, physical_capacity, sif2_capacity, sif4_capacity = (
            self._incident_measures(fluxes, shape))
        chemical_capacity = sif2_capacity + sif4_capacity
        adsorption_drive = par.initial_f_sticking_coefficient * f_flux
        denominator = adsorption_drive + 2.0 * chemical_capacity
        theta = np.zeros(shape)
        active = denominator > 0.0
        theta[active] = adsorption_drive[active] / denominator[active]

        physical_rate = physical_capacity * (1.0 - theta)
        sif2_rate = sif2_capacity * theta
        sif4_rate = sif4_capacity * theta
        consumed_f_rate = 2.0 * sif2_rate + 4.0 * sif4_rate
        if np.any(
            consumed_f_rate
            > f_flux + 64.0 * np.finfo(float).eps * np.maximum(f_flux, 1.0)
        ):
            raise RuntimeError(
                "Gray product routing demands more F atoms than incident")
        site_residual = np.abs(
            adsorption_drive * (1.0 - theta)
            - 2.0 * chemical_capacity * theta
        )
        duration = float(duration_s)
        physical_removed = physical_rate * duration
        sif2_removed = sif2_rate * duration
        sif4_removed = sif4_rate * duration
        removed = physical_removed + sif2_removed + sif4_removed
        consumed_f = consumed_f_rate * duration
        incident_f = f_flux * duration
        updated = GrayArFSiState(
            theta,
            state.removed_si_atoms_m2 + removed,
            state.consumed_f_atoms_m2 + consumed_f,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed},
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=(
                "summed unsaturated SiFx mass-spectral branch is reduced to SiF2",
                "emission energy-angle distributions are not measured",
            ),
        )
        products = (
            SurfaceProductPopulation(
                "Si_physical", "Si_atom", physical_removed, 1.0, 28.085,
                provenance={
                    "branch": "Gray Table 5-7 physical bare-Si sputter",
                    "missing": "differential emission energy-angle law",
                },
            ),
            SurfaceProductPopulation(
                "SiF2_fragment", "Si_atom", sif2_removed, 1.0, 66.082,
                provenance={
                    "branch": (
                        "Gray Eq. 5-27 summed unsaturated branch, represented "
                        "as SiF2 by the reduced model"),
                    "missing": "resolved SiFx distribution and differential emission",
                },
            ),
            SurfaceProductPopulation(
                "SiF4_chemical", "Si_atom", sif4_removed, 1.0, 104.079,
                provenance={
                    "branch": "Gray Table 5-7 ballistic-mixing SiF4 branch",
                    "missing": "differential emission energy-angle law",
                },
            ),
        )
        return GrayArFSiStepResult(
            state=updated,
            etch_velocity_m_s=(
                physical_rate + sif2_rate + sif4_rate
            ) / par.bulk_si_atom_density_m3,
            fluorinated_fraction=theta,
            physical_removed_si_atoms_m2=physical_removed,
            sif2_removed_si_atoms_m2=sif2_removed,
            sif4_removed_si_atoms_m2=sif4_removed,
            consumed_f_atoms_m2=consumed_f,
            incident_f_atoms_m2=incident_f,
            steady_site_balance_abs_residual_f_atoms_m2_s=site_residual,
            material_exchange=exchange,
            product_populations=products,
            validity=validity,
        )
