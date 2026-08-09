"""Deterministic species-resolved Cl+/Cl2+/Cl poly-Si surface closure.

This module joins two independently measured beam laws without averaging ion
identity: Chang's Cl+/Cl steady site balance and the Balooch Cl2+ yield line
reproduced in Chang Figure 5.7.  A returned-SiCl2 flux optionally activates
Chang Eq. 5.6 in the same site balance.  No feature depth enters the join.

The data do not provide the Cl2+ angular response or a conservative SiClx film
state.  Strict mode therefore admits Cl2+ only at normal incidence and the
measured high-chlorination state.  Non-strict mode transfers the measured Cl+
angular roll-off to Cl2+ as an explicitly reported sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from .chang_sawin_chlorine_si import (
    BaloochCl2IonSiParameters,
    ChangSawinArClSiState,
    ChangSawinClIonSiMechanism,
    ChangSawinClIonSiParameters,
    ChangSawinSiCl2SuppressionParameters,
    _chang_sawin_chemical_angular_factor,
)
from .surface_exchange import (
    SurfaceMaterialExchange,
    SurfaceProductPopulation,
    validate_surface_product_routing,
)
from .surface_kinetics import (
    FaceResolvedEnergeticFlux,
    MechanismValidity,
    SurfaceFluxes,
)


@dataclass(frozen=True)
class SpeciesResolvedChlorineSiStepResult:
    state: ChangSawinArClSiState
    etch_velocity_m_s: np.ndarray
    chlorination_fraction: np.ndarray
    clplus_removal_rate_si_m2_s: np.ndarray
    cl2plus_removal_rate_si_m2_s: np.ndarray
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
            "clplus_removal_rate_si_m2_s",
            "cl2plus_removal_rate_si_m2_s",
            "removed_si_atoms_m2",
            "sicl2_to_clplus_flux_ratio",
            "site_balance_residual_sites_m2_s",
        ):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid species-resolved chlorine/Si result")
            value.setflags(write=False)
            arrays[name] = value
            object.__setattr__(self, name, value)
        if (
            np.any(arrays["chlorination_fraction"] > 1.0)
            or not isinstance(self.state, ChangSawinArClSiState)
            or not isinstance(self.material_exchange, SurfaceMaterialExchange)
            or not isinstance(self.validity, MechanismValidity)
        ):
            raise ValueError("invalid species-resolved chlorine/Si contract")
        populations = validate_surface_product_routing(
            self.material_exchange, tuple(self.product_populations)
        )
        object.__setattr__(self, "product_populations", populations)


class SpeciesResolvedChlorineSiMechanism:
    """One common surface step for Cl+, Cl2+, Cl, and optional SiCl2."""

    quasi_steady_surface_state = True

    def __init__(
        self,
        clplus_parameters: ChangSawinClIonSiParameters | None = None,
        cl2plus_parameters: BaloochCl2IonSiParameters | None = None,
        suppression_parameters: ChangSawinSiCl2SuppressionParameters | None = None,
        cl2plus_coverage_mode: str = "saturated_card",
        strict_by_default: bool = True,
    ):
        self.clplus = ChangSawinClIonSiMechanism(clplus_parameters)
        self.cl2plus_parameters = (
            BaloochCl2IonSiParameters.chang_figure5_7()
            if cl2plus_parameters is None else cl2plus_parameters
        )
        self.suppression_parameters = (
            ChangSawinSiCl2SuppressionParameters.chang_thesis_equation_5_6()
            if suppression_parameters is None else suppression_parameters
        )
        if not isinstance(self.cl2plus_parameters, BaloochCl2IonSiParameters):
            raise TypeError("cl2plus_parameters must be BaloochCl2IonSiParameters")
        if not isinstance(
            self.suppression_parameters,
            ChangSawinSiCl2SuppressionParameters,
        ):
            raise TypeError(
                "suppression_parameters must be "
                "ChangSawinSiCl2SuppressionParameters"
            )
        if not isinstance(strict_by_default, (bool, np.bool_)):
            raise TypeError("strict_by_default must be boolean")
        if cl2plus_coverage_mode not in {
            "saturated_card", "source_bounded_linear"
        }:
            raise ValueError("unknown Cl2+ coverage mode")
        self.cl2plus_coverage_mode = str(cl2plus_coverage_mode)
        self.strict_by_default = bool(strict_by_default)

    @staticmethod
    def initial_state(shape=()):
        return ChangSawinArClSiState.bare(shape)

    @property
    def provenance(self):
        return MappingProxyType({
            "model": "species-resolved Chang/Balooch chlorine/poly-Si join",
            "clplus": dict(self.clplus.provenance),
            "cl2plus": {
                "source": "chang-thesis Figure 5.7; Balooch series",
                "slope_si_per_ion_per_sqrt_eV": (
                    self.cl2plus_parameters
                    .slope_si_per_ion_per_sqrt_eV
                ),
                "threshold_energy_eV": (
                    self.cl2plus_parameters.threshold_energy_eV
                ),
                "coverage_mode": self.cl2plus_coverage_mode,
                "source_bounded_zero_coverage_fraction": (
                    np.sqrt(2.0) * 0.06
                    / self.cl2plus_parameters.slope_si_per_ion_per_sqrt_eV
                ),
            },
            "sicl2": {
                "equation": "chang-thesis Eq. 5.6",
                "delta": (
                    self.suppression_parameters.sicl2_sticking_coefficient
                ),
                "eta": (
                    self.suppression_parameters
                    .chlorinated_sicl2_reaction_coefficient
                ),
            },
            "feature_depth_used_for_conditioning": False,
            "strict_by_default": self.strict_by_default,
        })

    @staticmethod
    def _positive_events(population):
        if isinstance(population, FaceResolvedEnergeticFlux):
            return (
                population.event_face,
                population.event_flux_m2_s,
                population.event_energy_eV,
                population.event_cosine_incidence,
                population.face_count,
            )
        return ChangSawinClIonSiMechanism._events(population)

    def validity(self, fluxes: SurfaceFluxes):
        neutral_unknown = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in {"Cl", "SiCl2"}
            and np.any(np.asarray(value) > 0.0)
        ))
        energetic_unknown = tuple(sorted({
            population.name
            for population in fluxes.energetic_fluxes
            if population.name not in {"Cl+", "Cl2+"}
            and np.any(np.asarray(population.flux_m2_s) > 0.0)
        }))
        reasons = []
        if neutral_unknown or energetic_unknown:
            reasons.append("positive incident flux has no declared chlorine/Si channel")
        clplus_only = SurfaceFluxes(
            {"Cl": fluxes.neutral_flux_m2_s.get("Cl", 0.0)},
            tuple(
                population for population in fluxes.energetic_fluxes
                if population.name == "Cl+"
            ),
        )
        base = self.clplus.validity(clplus_only)
        reasons.extend(base.reasons)
        cl2_energy_outside = False
        cl2_off_normal = False
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl2+":
                continue
            _, event_flux, energy, cosine, _ = self._positive_events(population)
            selected = event_flux > 0.0
            lower, upper = self.cl2plus_parameters.energy_domain_eV
            cl2_energy_outside |= bool(np.any(
                (energy[selected] < lower) | (energy[selected] > upper)
            ))
            cl2_off_normal |= bool(np.any(
                np.abs(cosine[selected] - 1.0) > 1.0e-10
            ))
        if cl2_energy_outside:
            reasons.append("Cl2+ energy leaves the measured 26--625 eV domain")
        if cl2_off_normal:
            reasons.append("Cl2+ angular yield was measured only at normal incidence")
        if self.cl2plus_coverage_mode == "source_bounded_linear":
            reasons.append(
                "Cl2+ coverage interpolation is bounded by printed slopes but unmeasured"
            )
        nonpredictive = tuple(sorted(set(base.nonpredictive_parameters) | {
            "cl2plus_angular_response",
            "single_coverage_closure",
            *(
                {"cl2plus_coverage_interpolation"}
                if self.cl2plus_coverage_mode == "source_bounded_linear"
                else set()
            ),
        }))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(dict.fromkeys(reasons)),
            unsupported_neutral_species=neutral_unknown,
            known_model_form_omissions=tuple(dict.fromkeys(
                base.known_model_form_omissions + (
                    "Cl2+ angular response is unresolved; non-strict mode transfers the measured Cl+ roll-off",
                    "Cl2+ product identity and chlorine consumption are unresolved",
                    "Eq. 5.6 has no conservative returned-SiClx film state",
                )
            )),
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive,
        )

    def _cl2plus_capacity(self, fluxes, shape, *, allow_extrapolation):
        capacity = np.zeros(shape)
        for population in fluxes.energetic_fluxes:
            if population.name != "Cl2+":
                continue
            face, event_flux, energy, cosine, face_count = (
                self._positive_events(population)
            )
            yield_si = self.cl2plus_parameters.yield_si_per_ion(
                energy, allow_extrapolation=allow_extrapolation
            )
            angular = _chang_sawin_chemical_angular_factor(cosine)
            if isinstance(population, FaceResolvedEnergeticFlux):
                local = np.bincount(
                    face,
                    weights=event_flux * yield_si * angular,
                    minlength=face_count,
                )
                capacity += np.broadcast_to(local, shape)
            else:
                capacity += np.broadcast_to(
                    np.sum(event_flux * yield_si * angular), shape
                )
        return capacity

    def advance(
        self,
        state: ChangSawinArClSiState,
        fluxes: SurfaceFluxes,
        duration_s: float,
        *,
        strict=None,
    ):
        if not isinstance(state, ChangSawinArClSiState):
            raise TypeError("species-resolved chlorine/Si state mismatch")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        strict = self.strict_by_default if strict is None else bool(strict)
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons)
            )
        shape = state.removed_si_atoms_m2.shape
        clplus_flux, ionic_drive, effective_s, clplus_capacity = (
            self.clplus._integrated_drives(
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
            + self.clplus.parameters.product_chlorine_atoms_per_si
            * clplus_capacity
            + 3.0 * eta * sicl2_drive
        )
        theta = np.zeros(shape)
        active = denominator > 0.0
        theta[active] = (
            (chlorine_drive[active] + sicl2_drive[active])
            / denominator[active]
        )
        clplus_rate = clplus_capacity * theta
        cl2plus_capacity = self._cl2plus_capacity(
            fluxes, shape, allow_extrapolation=not strict
        )
        if self.cl2plus_coverage_mode == "saturated_card":
            cl2plus_rate = cl2plus_capacity
        else:
            # Chang prints slopes 0.06 for Cl+ physical sputtering and 0.22
            # for Balooch Cl2+ in a saturated chlorine background, and states
            # that momentum scaling alone predicts sqrt(2) times the Cl+
            # value.  Those two endpoints bound this explicit sensitivity;
            # no feature observation selects the interpolation.
            bare_fraction = (
                np.sqrt(2.0) * 0.06
                / self.cl2plus_parameters.slope_si_per_ion_per_sqrt_eV
            )
            coverage_scale = bare_fraction + (1.0 - bare_fraction) * theta
            cl2plus_rate = cl2plus_capacity * coverage_scale
        if strict and np.any(
            (cl2plus_rate > 0.0)
            & (theta < self.cl2plus_parameters.minimum_high_chlorination_fraction)
        ):
            raise ValueError(
                "Cl2+ yield requires the measured high-chlorination surface state"
            )
        site_supply = (
            (chlorine_drive + sicl2_drive) * (1.0 - theta)
        )
        site_loss = theta * (
            self.clplus.parameters.product_chlorine_atoms_per_si
            * clplus_capacity
            + 3.0 * eta * sicl2_drive
        )
        site_residual = np.abs(site_supply - site_loss)
        duration = float(duration_s)
        removed_clplus = clplus_rate * duration
        removed_cl2plus = cl2plus_rate * duration
        removed = removed_clplus + removed_cl2plus
        updated = ChangSawinArClSiState(
            state.removed_si_atoms_m2 + removed,
            state.consumed_chlorine_atoms_m2
            + self.clplus.parameters.product_chlorine_atoms_per_si
            * removed_clplus,
        )
        exchange = SurfaceMaterialExchange(
            removed_units_m2={"Si_atom": removed},
            outgoing_units_m2={"Si_atom": removed},
            unresolved_units_m2={},
            deposited_units_m2={},
            known_limitations=validity.known_model_form_omissions,
        )
        products = (
            SurfaceProductPopulation(
                "SiCl4_from_Clplus",
                "Si_atom",
                removed_clplus,
                1.0,
                169.885,
                provenance={
                    "source": "chang-thesis Chapter 5",
                    "missing": "differential emission",
                },
            ),
            SurfaceProductPopulation(
                "SiClx_from_Cl2plus_unresolved",
                "Si_atom",
                removed_cl2plus,
                1.0,
                98.991,
                provenance={
                    "source": "chang-thesis Figure 5.7 Balooch series",
                    "missing": "product identity and differential emission",
                },
            ),
        )
        ratio = np.zeros(shape)
        active_clplus = clplus_flux > 0.0
        ratio[active_clplus] = (
            sicl2_flux[active_clplus] / clplus_flux[active_clplus]
        )
        return SpeciesResolvedChlorineSiStepResult(
            state=updated,
            etch_velocity_m_s=(
                (clplus_rate + cl2plus_rate)
                / self.clplus.parameters.bulk_si_atom_density_m3
            ),
            chlorination_fraction=theta,
            clplus_removal_rate_si_m2_s=clplus_rate,
            cl2plus_removal_rate_si_m2_s=cl2plus_rate,
            removed_si_atoms_m2=removed,
            sicl2_to_clplus_flux_ratio=ratio,
            site_balance_residual_sites_m2_s=site_residual,
            material_exchange=exchange,
            product_populations=products,
            validity=validity,
        )
