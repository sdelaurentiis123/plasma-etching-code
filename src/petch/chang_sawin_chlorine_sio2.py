"""Beam-measured Chang--Sawin Ar+/Cl etching of silicon dioxide.

Chang's MIT thesis, Chapter 4, reduces the oxide surface to one chlorinated
fraction and reports two absolute normal-incidence cards (70 and 100 eV), a
measured incidence-angle response, and the steady laws

``theta_Cl = s R / (s R + beta)``

``Y = c(phi) [Y0 + (beta - Y0) theta_Cl]``.

Here ``R`` is the atomic-Cl/Ar+ flux ratio and the yield unit is one removed
SiO2 formula unit per incident Ar+.  The implementation integrates a resolved
energy-angle quadrature deterministically.  It does not silently relabel Cl+
or Cl2+ as Ar+: projectile transfer is a separate, evidence-gated operation in
the reactor-to-mask diagnostic.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import unresolved_surface_exchange
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    MechanismValidity,
    ParameterEvidence,
    SurfaceFluxes,
)


# Original-pixel PIL digitization of Chang Figure 4.14.  Values are absolute
# SiO2-formula/Ar+ yields at 100 eV and Cl/Ar+=90; only their ratio is used as
# c(phi), preserving the independently printed Table-4.2 normalization.
FIGURE_4_14_ANGLE_DEG = np.asarray(
    (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0))
FIGURE_4_14_OXIDE_YIELD = np.asarray(
    (0.0735511714, 0.0737361282, 0.0758323058, 0.0907521578,
     0.1413070284, 0.1214549938, 0.0))


def chang_sawin_sio2_angular_factor(cosine_incidence):
    """Piecewise-linear replay of the measured Figure-4.14 oxide markers."""
    cosine = np.asarray(cosine_incidence, dtype=float)
    if (
        np.any(~np.isfinite(cosine))
        or np.any((cosine < 0.0) | (cosine > 1.0))
    ):
        raise ValueError("incidence cosines must lie in [0, 1]")
    angle = np.rad2deg(np.arccos(cosine))
    factor = np.interp(
        angle,
        FIGURE_4_14_ANGLE_DEG,
        FIGURE_4_14_OXIDE_YIELD / FIGURE_4_14_OXIDE_YIELD[0],
    )
    return np.maximum(factor, 0.0)


@dataclass(frozen=True)
class ChangSawinArClSiO2State:
    removed_sio2_formula_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        value = np.asarray(self.removed_sio2_formula_units_m2, dtype=float).copy()
        if np.any(~np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError("oxide removal inventory must be finite and nonnegative")
        value.setflags(write=False)
        object.__setattr__(self, "removed_sio2_formula_units_m2", value)

    @classmethod
    def bare(cls, shape=()):
        return cls(np.zeros(shape))

    def conservative_surface_fields(self):
        return {
            "removed_sio2_formula_units_m2": (
                self.removed_sio2_formula_units_m2),
        }

    def conservative_surface_upper_bounds(self):
        return {"removed_sio2_formula_units_m2": None}

    def surface_field_remap_modes(self):
        return {"removed_sio2_formula_units_m2": "conservative"}

    def with_conservative_surface_fields(self, fields):
        if set(fields) != {"removed_sio2_formula_units_m2"}:
            raise ValueError("Chang--Sawin SiO2 remap fields do not match")
        return type(self)(fields["removed_sio2_formula_units_m2"])


@dataclass(frozen=True)
class ChangSawinArClSiO2Parameters:
    ion_energy_eV: np.ndarray
    physical_yield_sio2_per_ion: np.ndarray
    full_chlorination_yield_sio2_per_ion: np.ndarray
    surface_chlorination_coefficient: np.ndarray
    bulk_sio2_formula_density_m3: float
    maximum_neutral_to_ion_ratio: float
    evidence: Mapping[str, ParameterEvidence]

    def __post_init__(self):
        energy = np.asarray(self.ion_energy_eV, dtype=float).copy()
        y0 = np.asarray(self.physical_yield_sio2_per_ion, dtype=float).copy()
        beta = np.asarray(
            self.full_chlorination_yield_sio2_per_ion, dtype=float).copy()
        sticking = np.asarray(
            self.surface_chlorination_coefficient, dtype=float).copy()
        evidence = dict(self.evidence)
        required = {
            "beam_cards", "angular_response", "bulk_sio2_formula_density_m3",
        }
        if (
            energy.ndim != 1
            or energy.size < 2
            or y0.shape != energy.shape
            or beta.shape != energy.shape
            or sticking.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(np.diff(energy) <= 0.0)
            or np.any(~np.isfinite(y0))
            or np.any(y0 < 0.0)
            or np.any(~np.isfinite(beta))
            or np.any(beta < y0)
            or np.any(~np.isfinite(sticking))
            or np.any(sticking < 0.0)
            or not np.isfinite(self.bulk_sio2_formula_density_m3)
            or self.bulk_sio2_formula_density_m3 <= 0.0
            or not np.isfinite(self.maximum_neutral_to_ion_ratio)
            or self.maximum_neutral_to_ion_ratio <= 0.0
            or set(evidence) != required
            or any(not isinstance(item, ParameterEvidence)
                   for item in evidence.values())
        ):
            raise ValueError("invalid Chang--Sawin Ar+/Cl/SiO2 parameters")
        for value in (energy, y0, beta, sticking):
            value.setflags(write=False)
        object.__setattr__(self, "ion_energy_eV", energy)
        object.__setattr__(self, "physical_yield_sio2_per_ion", y0)
        object.__setattr__(
            self, "full_chlorination_yield_sio2_per_ion", beta)
        object.__setattr__(self, "surface_chlorination_coefficient", sticking)
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @classmethod
    def chang_thesis_table_4_2(cls):
        thesis = (
            "chang-thesis, Chapter 4, Eqs. 4.3--4.6 and Table 4.2"
        )
        return cls(
            ion_energy_eV=np.asarray((70.0, 100.0)),
            physical_yield_sio2_per_ion=np.asarray((0.01, 0.02)),
            full_chlorination_yield_sio2_per_ion=np.asarray((0.04, 0.08)),
            surface_chlorination_coefficient=np.asarray((0.001, 0.005)),
            bulk_sio2_formula_density_m3=2.2042e28,
            maximum_neutral_to_ion_ratio=120.0,
            evidence={
                "beam_cards": ParameterEvidence(
                    thesis,
                    "controlled-beam measurement/regression",
                    note=(
                        "Table 4.2 reports (Y0,beta,s)=(0.01,0.04,0.001) "
                        "at 70 eV and (0.02,0.08,0.005) at 100 eV"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "angular_response": ParameterEvidence(
                    "chang-thesis, Chapter 4, Figures 4.7 and 4.14",
                    "controlled-beam angular measurement",
                    note=(
                        "100 eV Ar+/Cl at Cl/Ar+=90; original-pixel PIL "
                        "digitization with maximum near 60 degrees"
                    ),
                    supports_prediction_within_declared_domain=True,
                ),
                "bulk_sio2_formula_density_m3": ParameterEvidence(
                    "2.20 g/cm3 fused SiO2 and 60.0843 g/mol molar mass",
                    "derived physical constant",
                    supports_prediction_within_declared_domain=True,
                ),
            },
        )

    @property
    def energy_domain_eV(self):
        return float(self.ion_energy_eV[0]), float(self.ion_energy_eV[-1])

    def coefficients(self, energy_eV, *, allow_extrapolation=False):
        energy = np.asarray(energy_eV, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
            raise ValueError("ion energies must be finite and nonnegative")
        lower, upper = self.energy_domain_eV
        if not allow_extrapolation and np.any((energy < lower) | (energy > upper)):
            raise ValueError("oxide energy leaves the measured 70--100 eV cards")
        coordinate = np.sqrt(energy)
        source = np.sqrt(self.ion_energy_eV)

        def interpolate(values):
            result = np.interp(coordinate, source, values)
            if allow_extrapolation:
                low_slope = (values[1] - values[0]) / (source[1] - source[0])
                high_slope = (values[-1] - values[-2]) / (
                    source[-1] - source[-2])
                result = np.where(
                    coordinate < source[0],
                    values[0] + low_slope * (coordinate - source[0]),
                    result,
                )
                result = np.where(
                    coordinate > source[-1],
                    values[-1] + high_slope * (coordinate - source[-1]),
                    result,
                )
            return np.maximum(result, 0.0)

        return (
            interpolate(self.physical_yield_sio2_per_ion),
            interpolate(self.full_chlorination_yield_sio2_per_ion),
            interpolate(self.surface_chlorination_coefficient),
        )


@dataclass(frozen=True)
class ChangSawinArClSiO2StepResult:
    state: ChangSawinArClSiO2State
    etch_velocity_m_s: np.ndarray
    chlorination_fraction: np.ndarray
    removal_rate_sio2_formula_m2_s: np.ndarray
    mean_yield_sio2_formula_per_ion: np.ndarray
    removed_sio2_formula_units_m2: np.ndarray
    material_exchange: object
    validity: MechanismValidity

    def __post_init__(self):
        for name in (
            "etch_velocity_m_s", "chlorination_fraction",
            "removal_rate_sio2_formula_m2_s",
            "mean_yield_sio2_formula_per_ion",
            "removed_sio2_formula_units_m2",
        ):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("invalid Chang--Sawin oxide result")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if np.any(self.chlorination_fraction > 1.0):
            raise ValueError("oxide chlorination fraction exceeds unity")


class ChangSawinArClSiO2Mechanism:
    """Absolute Ar+/Cl/SiO2 rate law inside the Chapter-4 beam domain."""

    quasi_steady_surface_state = True

    def __init__(self, parameters: ChangSawinArClSiO2Parameters | None = None):
        self.parameters = (
            ChangSawinArClSiO2Parameters.chang_thesis_table_4_2()
            if parameters is None else parameters
        )
        if not isinstance(self.parameters, ChangSawinArClSiO2Parameters):
            raise TypeError("parameters must be ChangSawinArClSiO2Parameters")

    @staticmethod
    def initial_state(shape=()):
        return ChangSawinArClSiO2State.bare(shape)

    @property
    def provenance(self):
        par = self.parameters
        return MappingProxyType({
            "model": "Chang--Sawin steady Ar+/Cl/SiO2 site balance",
            "equations": ["Chang thesis Eq. 4.4", "Eq. 4.5"],
            "measured_energy_cards_eV": par.ion_energy_eV.tolist(),
            "angular_data": "checksum-backed Figure 4.14 PIL digitization",
            "projectile": "Ar+ only",
            "feature_depth_used_for_conditioning": False,
        })

    @staticmethod
    def _events(population):
        if isinstance(population, FaceResolvedEnergeticFlux):
            return (
                population.event_face,
                population.event_flux_m2_s,
                population.event_energy_eV,
                population.event_cosine_incidence,
                population.face_count,
            )
        if isinstance(population, EnergeticFlux):
            flux = float(np.asarray(population.flux_m2_s))
            return (
                np.zeros(population.weight.size, dtype=int),
                flux * population.weight,
                population.energy_eV,
                population.cosine_incidence,
                1,
            )
        raise TypeError(type(population).__name__)

    def validity(self, fluxes: SurfaceFluxes):
        par = self.parameters
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name != "Cl" and np.any(np.asarray(value) > 0.0)
        ))
        unsupported_energetic = tuple(sorted({
            item.name for item in fluxes.energetic_fluxes
            if item.name != "Ar+" and np.any(np.asarray(item.flux_m2_s) > 0.0)
        }))
        total_ion = None
        energy_outside = False
        for population in fluxes.energetic_fluxes:
            if population.name != "Ar+":
                continue
            if isinstance(population, FaceResolvedEnergeticFlux):
                selected = population.event_flux_m2_s > 0.0
                energy = population.event_energy_eV[selected]
                ion = population.flux_m2_s
            else:
                selected = population.weight > 0.0
                energy = population.energy_eV[selected]
                ion = np.asarray(population.flux_m2_s)
            lower, upper = par.energy_domain_eV
            energy_outside |= bool(np.any((energy < lower) | (energy > upper)))
            total_ion = ion if total_ion is None else total_ion + ion
        ratio_outside = False
        if total_ion is not None:
            neutral, ion = np.broadcast_arrays(
                np.asarray(fluxes.neutral_flux_m2_s.get("Cl", 0.0)),
                np.asarray(total_ion),
            )
            active = ion > 0.0
            ratio = np.zeros(ion.shape)
            np.divide(neutral, ion, out=ratio, where=active)
            ratio_outside = bool(np.any(
                active & (ratio > par.maximum_neutral_to_ion_ratio)
            ))
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive flux has no declared Ar+/Cl/SiO2 channel")
        if energy_outside:
            reasons.append("Ar+ energy leaves the measured 70--100 eV cards")
        if ratio_outside:
            reasons.append("Cl/Ar+ ratio leaves the measured 0--120 domain")
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "Cl+ and Cl2+ projectile transfer is not supplied by this Ar+ beam law",
                "broad IEADs are deterministic quadrature extensions of monoenergetic cards",
                "product identity and differential emission are unresolved",
                "surface charging, implantation depth, and transient oxide stoichiometry are unresolved",
            ),
            parameter_evidence_supports_prediction=True,
            nonpredictive_parameters=(),
        )

    def advance(self, state, fluxes, duration_s, *, strict=True):
        if not isinstance(state, ChangSawinArClSiO2State):
            raise TypeError("Chang--Sawin oxide state mismatch")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons)
            )
        shape = state.removed_sio2_formula_units_m2.shape
        ion_flux = np.zeros(shape)
        beta_flux = np.zeros(shape)
        s_flux = np.zeros(shape)
        events = []
        for population in fluxes.energetic_fluxes:
            if population.name != "Ar+":
                continue
            face, event_flux, energy, cosine, face_count = self._events(population)
            y0, beta, sticking = self.parameters.coefficients(
                energy, allow_extrapolation=not strict
            )
            if isinstance(population, FaceResolvedEnergeticFlux):
                integrate = lambda values: np.bincount(  # noqa: E731
                    face, weights=values, minlength=face_count)
                ion_flux += np.broadcast_to(integrate(event_flux), shape)
                beta_flux += np.broadcast_to(integrate(event_flux * beta), shape)
                s_flux += np.broadcast_to(integrate(event_flux * sticking), shape)
            else:
                ion_flux += np.broadcast_to(np.sum(event_flux), shape)
                beta_flux += np.broadcast_to(np.sum(event_flux * beta), shape)
                s_flux += np.broadcast_to(np.sum(event_flux * sticking), shape)
            events.append((population, face, event_flux, y0, beta, cosine, face_count))

        effective_s = np.zeros(shape)
        active_ion = ion_flux > 0.0
        effective_s[active_ion] = s_flux[active_ion] / ion_flux[active_ion]
        neutral = np.broadcast_to(
            np.asarray(fluxes.neutral_flux_m2_s.get("Cl", 0.0), dtype=float),
            shape,
        )
        supply = effective_s * neutral
        denominator = supply + beta_flux
        chlorination = np.zeros(shape)
        active = denominator > 0.0
        chlorination[active] = supply[active] / denominator[active]

        removal_rate = np.zeros(shape)
        for population, face, event_flux, y0, beta, cosine, face_count in events:
            angular = chang_sawin_sio2_angular_factor(cosine)
            if isinstance(population, FaceResolvedEnergeticFlux):
                local_theta = chlorination[face]
                event_yield = angular * (y0 + (beta - y0) * local_theta)
                removal_rate += np.broadcast_to(np.bincount(
                    face,
                    weights=event_flux * event_yield,
                    minlength=face_count,
                ), shape)
            else:
                event_yield = angular * (
                    y0 + (beta - y0) * float(chlorination)
                )
                removal_rate += np.broadcast_to(
                    np.sum(event_flux * event_yield), shape)

        mean_yield = np.zeros(shape)
        mean_yield[active_ion] = removal_rate[active_ion] / ion_flux[active_ion]
        removed = removal_rate * float(duration_s)
        updated = ChangSawinArClSiO2State(
            state.removed_sio2_formula_units_m2 + removed
        )
        exchange = unresolved_surface_exchange(
            removed_units_m2={"SiO2_formula": removed},
            limitations=(
                "Chang assumes Si, O2, and SiCl2 major products but does not measure branching",
                "product energy-angle distributions are unresolved",
            ),
        )
        return ChangSawinArClSiO2StepResult(
            state=updated,
            etch_velocity_m_s=(
                removal_rate / self.parameters.bulk_sio2_formula_density_m3
            ),
            chlorination_fraction=chlorination,
            removal_rate_sio2_formula_m2_s=removal_rate,
            mean_yield_sio2_formula_per_ion=mean_yield,
            removed_sio2_formula_units_m2=removed,
            material_exchange=exchange,
            validity=validity,
        )
