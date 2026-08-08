"""Published-model chlorine ion mobility and global wall transport.

Lymberopoulos and Economou (1995) publish reduced mobilities for ``Cl2+``,
``Cl+``, and ``Cl-``.  The paper states that the underlying constant collision
cross sections came from a private communication and does not print them.
Ramamurthi and Economou (2002) reuse the same values at a different ion
temperature.  The exact values are therefore useful for source-model
reproduction but are quarantined from predictive evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from petch.sheath import bohm_speed

from .chlorine import CHLORINE_ATOM_MASS_AMU, CHLORINE_MOLECULE_MASS_AMU
from .chlorine_particle_model import (
    REACTOR_SCALAR_EVIDENCE_KINDS,
    ChlorineChargedTransportState,
    ChlorineFixedPressureCondition,
    PositiveIonWallTransport,
    ReactorScalarInput,
)
from .network import E_CHARGE_C
from .transport import ATOMIC_MASS_UNIT_KG

LYMBEROPOULOS_1995_REDUCED_CL2PLUS_MOBILITY_M_INV_V_INV_S_INV = 5.62e21
LYMBEROPOULOS_1995_REDUCED_CLPLUS_MOBILITY_M_INV_V_INV_S_INV = 6.48e21
LYMBEROPOULOS_1995_REDUCED_CLMINUS_MOBILITY_M_INV_V_INV_S_INV = 6.48e21

_PREDICTIVE_EVIDENCE_KINDS = frozenset({"measured", "validated_model"})
_POSITIVE_ION_MASS_AMU = MappingProxyType({
    "Cl2+": CHLORINE_MOLECULE_MASS_AMU,
    "Cl+": CHLORINE_ATOM_MASS_AMU,
})


@dataclass(frozen=True)
class ReducedIonMobility:
    """Density-reduced ion mobility with a strict ion-temperature domain."""

    reduced_mobility_m_inv_V_inv_s_inv: float
    reference_ion_temperature_eV: float
    valid_ion_temperature_eV: tuple[float, float]
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        reduced_mobility = float(
            self.reduced_mobility_m_inv_V_inv_s_inv)
        reference_temperature = float(self.reference_ion_temperature_eV)
        try:
            lower, upper = (
                float(value) for value in self.valid_ion_temperature_eV)
        except (TypeError, ValueError):
            raise ValueError(
                "ion-temperature domain must contain two numbers") from None
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            not np.isfinite(reduced_mobility)
            or reduced_mobility <= 0.0
            or not np.isfinite(reference_temperature)
            or reference_temperature <= 0.0
            or not np.isfinite(lower)
            or not np.isfinite(upper)
            or lower <= 0.0
            or upper < lower
            or not lower <= reference_temperature <= upper
            or not str(self.source).strip()
            or self.evidence_kind not in REACTOR_SCALAR_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid reduced ion mobility")
        object.__setattr__(
            self,
            "reduced_mobility_m_inv_V_inv_s_inv",
            reduced_mobility,
        )
        object.__setattr__(
            self, "reference_ion_temperature_eV", reference_temperature)
        object.__setattr__(
            self, "valid_ion_temperature_eV", (lower, upper))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS
            and self.relative_uncertainty is not None
        )

    def evaluate(
        self,
        *,
        total_neutral_density_m3: float,
        ion_temperature_eV: float,
    ) -> "IonMobilityState":
        density = float(total_neutral_density_m3)
        temperature = float(ion_temperature_eV)
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError("total neutral density must be positive")
        if (
            not np.isfinite(temperature)
            or not self.valid_ion_temperature_eV[0]
            <= temperature
            <= self.valid_ion_temperature_eV[1]
        ):
            raise ValueError(
                "ion temperature is outside the mobility evidence domain")
        return IonMobilityState(
            mobility_m2_V_s=(
                self.reduced_mobility_m_inv_V_inv_s_inv / density),
            total_neutral_density_m3=density,
            ion_temperature_eV=temperature,
            source=self.source,
            evidence_kind=self.evidence_kind,
            relative_uncertainty=self.relative_uncertainty,
            provenance={
                **self.provenance,
                "reduced_mobility_m_inv_V_inv_s_inv": (
                    self.reduced_mobility_m_inv_V_inv_s_inv),
                "reference_ion_temperature_eV": (
                    self.reference_ion_temperature_eV),
                "valid_ion_temperature_eV": self.valid_ion_temperature_eV,
            },
        )


@dataclass(frozen=True)
class IonMobilityState:
    """One evaluated ion mobility and its source conditions."""

    mobility_m2_V_s: float
    total_neutral_density_m3: float
    ion_temperature_eV: float
    source: str
    evidence_kind: str
    relative_uncertainty: float | None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        values = np.asarray([
            self.mobility_m2_V_s,
            self.total_neutral_density_m3,
            self.ion_temperature_eV,
        ], dtype=float)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not str(self.source).strip()
            or self.evidence_kind not in REACTOR_SCALAR_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid ion mobility state")
        object.__setattr__(self, "mobility_m2_V_s", float(values[0]))
        object.__setattr__(
            self, "total_neutral_density_m3", float(values[1]))
        object.__setattr__(self, "ion_temperature_eV", float(values[2]))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS
            and self.relative_uncertainty is not None
        )


def lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities(
) -> Mapping[str, ReducedIonMobility]:
    """Return the exact 1995 Table-text reduced mobilities in SI."""
    values = {
        "Cl2+": LYMBEROPOULOS_1995_REDUCED_CL2PLUS_MOBILITY_M_INV_V_INV_S_INV,
        "Cl+": LYMBEROPOULOS_1995_REDUCED_CLPLUS_MOBILITY_M_INV_V_INV_S_INV,
        "Cl-": LYMBEROPOULOS_1995_REDUCED_CLMINUS_MOBILITY_M_INV_V_INV_S_INV,
    }
    return MappingProxyType({
        species: ReducedIonMobility(
            reduced_mobility_m_inv_V_inv_s_inv=value,
            reference_ion_temperature_eV=0.12,
            valid_ion_temperature_eV=(0.12, 0.12),
            source=(
                "lymberopoulos-economou-1995-chlorine Table-text "
                f"{species} reduced mobility"
            ),
            evidence_kind="published_model",
            relative_uncertainty=None,
            provenance={
                "doi": "10.1109/27.467977",
                "source_value_cm-1_V-1_s-1": value / 100.0,
                "source_collision_data": (
                    "constant cross sections from private communication "
                    "with M. J. Kushner; values not printed"
                ),
                "coefficient_selection_target": None,
            },
        )
        for species, value in values.items()
    })


def ramamurthi_economou_2002_chlorine_reduced_ion_mobilities(
) -> Mapping[str, ReducedIonMobility]:
    """Return the 2002 reuse at 300 K, quarantined from prediction."""
    temperature_eV = 300.0 * 1.380649e-23 / E_CHARGE_C
    values = {
        "Cl2+": LYMBEROPOULOS_1995_REDUCED_CL2PLUS_MOBILITY_M_INV_V_INV_S_INV,
        "Cl+": LYMBEROPOULOS_1995_REDUCED_CLPLUS_MOBILITY_M_INV_V_INV_S_INV,
        "Cl-": LYMBEROPOULOS_1995_REDUCED_CLMINUS_MOBILITY_M_INV_V_INV_S_INV,
    }
    return MappingProxyType({
        species: ReducedIonMobility(
            reduced_mobility_m_inv_V_inv_s_inv=value,
            reference_ion_temperature_eV=temperature_eV,
            valid_ion_temperature_eV=(temperature_eV, temperature_eV),
            source=(
                "ramamurthi-economou-2002-chlorine Table III "
                f"{species} reduced mobility"
            ),
            evidence_kind="published_model",
            relative_uncertainty=None,
            provenance={
                "doi": "10.1116/1.1450581",
                "source_value_cm-1_V-1_s-1": value / 100.0,
                "source_reference": (
                    "Lymberopoulos-Economou 1995, DOI 10.1109/27.467977"),
                "temperature_conflict": (
                    "same reduced mobility appears at 0.12 eV in 1995"),
                "coefficient_selection_target": None,
            },
        )
        for species, value in values.items()
    })


@dataclass(frozen=True)
class LeeEconomouChlorineChargedTransportProvider:
    """Compose published ion mobilities with Lee global edge factors.

    Mobility supplies ``nu_m = e/(m mu)`` and ``D_a = mu Te``.  The momentum
    mean free path uses the Maxwellian mean ion speed at the declared ion
    temperature.  Species-specific edge factors are retained; this relaxes
    Lee's common-edge-factor approximation while keeping every input visible.
    """

    reduced_mobilities: Mapping[str, ReducedIonMobility]
    ion_temperature: ReactorScalarInput
    name: str = "lee_economou_chlorine_charged_transport"
    version: str = "1"

    def __post_init__(self):
        mobilities = dict(self.reduced_mobilities)
        if (
            set(mobilities) != set(_POSITIVE_ION_MASS_AMU)
            or any(
                not isinstance(item, ReducedIonMobility)
                for item in mobilities.values()
            )
            or not isinstance(self.ion_temperature, ReactorScalarInput)
            or self.ion_temperature.unit != "eV"
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid chlorine mobility transport provider")
        object.__setattr__(
            self, "reduced_mobilities", MappingProxyType(mobilities))

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
    ) -> ChlorineChargedTransportState:
        if not isinstance(condition, ChlorineFixedPressureCondition):
            raise TypeError("chlorine fixed-pressure condition is required")
        densities = {
            str(name): float(value) for name, value in densities_m3.items()
        }
        required = {"e", "Cl2", "Cl", "Cl2+", "Cl+", "Cl-"}
        if (
            set(densities) != required
            or any(
                not np.isfinite(value) or value <= 0.0
                for value in densities.values()
            )
        ):
            raise ValueError("invalid chlorine density state")
        electron_temperature = condition.electron_temperature.value
        ion_temperature = self.ion_temperature.value
        electronegativity = densities["Cl-"] / densities["e"]
        total_neutral_density = densities["Cl2"] + densities["Cl"]
        transport = {}
        for species, mass_amu in _POSITIVE_ION_MASS_AMU.items():
            mobility = self.reduced_mobilities[species].evaluate(
                total_neutral_density_m3=total_neutral_density,
                ion_temperature_eV=ion_temperature,
            )
            mass_kg = mass_amu * ATOMIC_MASS_UNIT_KG
            momentum_collision_frequency = (
                E_CHARGE_C / (mass_kg * mobility.mobility_m2_V_s))
            mean_ion_speed = np.sqrt(
                8.0 * E_CHARGE_C * ion_temperature
                / (np.pi * mass_kg)
            )
            momentum_mean_free_path = (
                mean_ion_speed / momentum_collision_frequency)
            ambipolar_diffusivity = (
                mobility.mobility_m2_V_s * electron_temperature)
            speed = bohm_speed(electron_temperature, mass_amu)
            edge = condition.geometry.electronegative_edge_factors(
                electronegativity=electronegativity,
                electron_to_ion_temperature_ratio=(
                    electron_temperature / ion_temperature),
                ion_mean_free_path_m=momentum_mean_free_path,
                bohm_speed_m_s=speed,
                ambipolar_diffusion_m2_s=ambipolar_diffusivity,
            )
            transport[species] = PositiveIonWallTransport(
                axial_flux_velocity_m_s=edge.axial * speed,
                radial_flux_velocity_m_s=edge.radial * speed,
                source=(
                    f"{mobility.source}; lee-lieberman-1994-global "
                    "Eqs. 13-14"
                ),
                evidence_kind="published_model",
                relative_uncertainty=None,
                provenance={
                    **mobility.provenance,
                    "mobility_m2_V_s": mobility.mobility_m2_V_s,
                    "total_neutral_density_m3": total_neutral_density,
                    "momentum_collision_frequency_s_inv": (
                        momentum_collision_frequency),
                    "mean_ion_speed_m_s": mean_ion_speed,
                    "momentum_mean_free_path_m": momentum_mean_free_path,
                    "ambipolar_diffusivity_m2_s": ambipolar_diffusivity,
                    "bohm_speed_m_s": speed,
                    "electronegativity": electronegativity,
                    "edge_factor_axial": edge.axial,
                    "edge_factor_radial": edge.radial,
                    "mean_free_path_velocity_convention": (
                        "Maxwellian mean ion speed sqrt(8 e Ti/(pi m))"),
                    "edge_factor_species_assumption": (
                        "species-specific; relaxes Lee common-edge factor"),
                },
            )
        return ChlorineChargedTransportState(
            geometry=condition.geometry,
            positive_ion_transport=transport,
            negative_ion_confinement_source=(
                "lee-lieberman-1994-global parabolic negative-ion profile "
                "with zero sheath-edge density"
            ),
            negative_ion_confinement_evidence="published_model",
            negative_ion_confinement_relative_uncertainty=None,
        )
