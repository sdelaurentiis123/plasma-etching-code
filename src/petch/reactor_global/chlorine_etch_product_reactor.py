"""Deterministic 0-D transfer from wafer Si removal to SiClx fluxes.

At fixed chlorine/electron state the Lee--Graves--Lieberman Table-IV network
is linear in every Si-containing species.  This module exploits that structure
instead of time marching or Monte Carlo: one dense linear solve closes product
fragmentation, ionization, mutual neutralization, pumping, and the paper's two
wall limits.  It is intended to sit inside a scalar wafer/surface fixed point.

The calculation is one-way with respect to the base chlorine plasma.  It
therefore reports, but does not hide, the missing feedback of product collision
power and released chlorine into the EEPF/chlorine balances.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from petch.sheath import bohm_speed

from .chlorine_etch_products import (
    build_lee_graves_lieberman_etch_product_network,
)
from .geometry import CylindricalReactor, ElectropositiveEdgeFactors
from .network import RateContext
from .network import E_CHARGE_C, ElectronArrheniusRateCoefficient
from .wafer_power_transfer import isotropic_thermal_particle_flux_m2_s


_NEUTRALS = ("Si", "SiCl", "SiCl2", "SiCl3", "SiCl4")
_IONS = tuple(f"{name}+" for name in _NEUTRALS)
_PRODUCTS = _NEUTRALS + _IONS
_MASS_AMU = MappingProxyType({
    species.name: species.mass_amu
    for species in build_lee_graves_lieberman_etch_product_network().species
    if species.name in _PRODUCTS
})


@dataclass(frozen=True)
class EtchProductWallBoundary:
    """One explicit Lee wall limit for the silicon-containing species."""

    name: str
    neutral_sticking_probability: Mapping[str, float]
    positive_ion_sticking_probability: float
    source: str

    def __post_init__(self):
        neutral = {
            str(name): float(value)
            for name, value in self.neutral_sticking_probability.items()
        }
        ion = float(self.positive_ion_sticking_probability)
        if (
            not str(self.name).strip()
            or set(neutral) != set(_NEUTRALS)
            or any(
                not math.isfinite(value) or not 0.0 <= value <= 1.0
                for value in neutral.values()
            )
            or not math.isfinite(ion)
            or not 0.0 <= ion <= 1.0
            or not str(self.source).strip()
        ):
            raise ValueError("invalid etch-product wall boundary")
        object.__setattr__(
            self,
            "neutral_sticking_probability",
            MappingProxyType(neutral),
        )
        object.__setattr__(self, "positive_ion_sticking_probability", ion)


def lee_1995_reflective_product_wall() -> EtchProductWallBoundary:
    """Perfect product reflection/ion neutralization, Section 2.3."""
    return EtchProductWallBoundary(
        name="lee_1995_nonreactive_reflective",
        neutral_sticking_probability={name: 0.0 for name in _NEUTRALS},
        positive_ion_sticking_probability=0.0,
        source=(
            "lee-graves-lieberman-1995-etch-products Section 2.3 "
            "nonreactive-wall limit"
        ),
    )


def lee_1995_reactive_product_wall() -> EtchProductWallBoundary:
    """Silicon-covered reactive-wall endpoint from Section 2.3.1.

    Only Si and SiCl neutral redeposition is retained, as required by item 1.
    Their probabilities are the paper's silicon-passivated endpoint.  All
    positive ions use the paper's assumed probability 0.5; the remainder is
    neutralized and returned to the gas.
    """
    return EtchProductWallBoundary(
        name="lee_1995_reactive_silicon_passivated_endpoint",
        neutral_sticking_probability={
            "Si": 1.0,
            "SiCl": 0.4,
            "SiCl2": 0.0,
            "SiCl3": 0.0,
            "SiCl4": 0.0,
        },
        positive_ion_sticking_probability=0.5,
        source=(
            "lee-graves-lieberman-1995-etch-products Section 2.3.1 "
            "items 1, 4, and 5"
        ),
    )


@dataclass(frozen=True)
class EtchProductPlasmaCondition:
    geometry: CylindricalReactor
    neutral_control_volume_m3: float
    electron_density_m3: float
    chlorine_atom_density_m3: float
    chlorine_negative_ion_density_m3: float
    electron_temperature_eV: float
    gas_temperature_K: float
    exhaust_loss_frequency_s_inv: float
    common_edge_factors: ElectropositiveEdgeFactors
    wall_boundary: EtchProductWallBoundary
    source: str

    def __post_init__(self):
        values = (
            self.neutral_control_volume_m3,
            self.electron_density_m3,
            self.chlorine_atom_density_m3,
            self.chlorine_negative_ion_density_m3,
            self.electron_temperature_eV,
            self.gas_temperature_K,
            self.exhaust_loss_frequency_s_inv,
        )
        if (
            not isinstance(self.geometry, CylindricalReactor)
            or any(not math.isfinite(float(value)) or value <= 0.0 for value in values)
            or self.neutral_control_volume_m3 < self.geometry.volume_m3
            or not isinstance(self.common_edge_factors, ElectropositiveEdgeFactors)
            or not isinstance(self.wall_boundary, EtchProductWallBoundary)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid etch-product plasma condition")

    @property
    def active_volume_fraction(self) -> float:
        return float(self.geometry.volume_m3 / self.neutral_control_volume_m3)


@dataclass(frozen=True)
class EtchProductPlasmaSolution:
    densities_m3: Mapping[str, float]
    wafer_neutral_flux_m2_s: Mapping[str, float]
    wafer_positive_ion_flux_m2_s: Mapping[str, float]
    pumping_loss_molecule_s: Mapping[str, float]
    wall_deposition_loss_molecule_s: Mapping[str, float]
    positive_ion_wall_loss_molecule_s: Mapping[str, float]
    chlorine_atom_source_m3_s: float
    table4_threshold_power_lower_bound_W_m3: float
    injected_si_atom_s: float
    silicon_inventory_relative_residual: float
    linear_balance_maximum_relative_residual: float
    matrix_condition_number: float
    missing_feedback_closures: tuple[str, ...] = (
        "SiClx collision energy loss in the self-consistent EEPF",
        "released chlorine in the base chlorine particle balance",
        "product-ion mobility or species-specific edge factors",
        "conditioned chamber-wall Si/Cl coverages between limiting endpoints",
    )

    def __post_init__(self):
        mappings = (
            (self.densities_m3, _PRODUCTS),
            (self.wafer_neutral_flux_m2_s, _NEUTRALS),
            (self.wafer_positive_ion_flux_m2_s, _IONS),
            (self.pumping_loss_molecule_s, _NEUTRALS),
            (self.wall_deposition_loss_molecule_s, _PRODUCTS),
            (self.positive_ion_wall_loss_molecule_s, _IONS),
        )
        for raw, names in mappings:
            values = {str(name): float(value) for name, value in raw.items()}
            if set(values) != set(names) or any(
                not math.isfinite(value) or value < 0.0
                for value in values.values()
            ):
                raise ValueError("invalid etch-product solution mapping")
            object.__setattr__(
                self,
                next(
                    field for field in (
                        "densities_m3", "wafer_neutral_flux_m2_s",
                        "wafer_positive_ion_flux_m2_s",
                        "pumping_loss_molecule_s",
                        "wall_deposition_loss_molecule_s",
                        "positive_ion_wall_loss_molecule_s",
                    ) if getattr(self, field) is raw
                ),
                MappingProxyType(values),
            )
        scalars = (
            self.injected_si_atom_s,
            self.silicon_inventory_relative_residual,
            self.linear_balance_maximum_relative_residual,
            self.matrix_condition_number,
            self.chlorine_atom_source_m3_s,
            self.table4_threshold_power_lower_bound_W_m3,
        )
        if (
            any(not math.isfinite(float(value)) for value in scalars)
            or self.injected_si_atom_s < 0.0
            or abs(self.silicon_inventory_relative_residual) > 2.0e-10
            or self.linear_balance_maximum_relative_residual > 2.0e-10
            or self.matrix_condition_number <= 0.0
            or self.chlorine_atom_source_m3_s < 0.0
            or self.table4_threshold_power_lower_bound_W_m3 < 0.0
        ):
            raise ValueError("etch-product solution failed conservation")

    @property
    def total_neutral_density_m3(self) -> float:
        return float(sum(self.densities_m3[name] for name in _NEUTRALS))

    @property
    def total_positive_ion_density_m3(self) -> float:
        return float(sum(self.densities_m3[name] for name in _IONS))

    def chlorine_feedback_lower_bound(self):
        """Return conserved pressure/charge feedback and minimum power loss.

        The Table-IV Arrhenius activation energies are lower bounds on the
        corresponding inelastic loss.  The inverse-polynomial channels and
        the Table-V excitation/momentum ledger are deliberately assigned zero
        here, so this object cannot claim a predictive electron-power closure.
        """
        from .chlorine_eedf_model import EEDFChlorineFixedFeedback

        return EEDFChlorineFixedFeedback(
            chlorine_species_source_m3_s={
                "Cl2": 0.0,
                "Cl": self.chlorine_atom_source_m3_s,
            },
            extra_neutral_density_m3=self.total_neutral_density_m3,
            extra_positive_charge_density_m3=(
                self.total_positive_ion_density_m3),
            extra_collisional_power_density_W_m3=(
                self.table4_threshold_power_lower_bound_W_m3),
            extra_charged_wall_power_density_W_m3=0.0,
            source=(
                "Lee-Graves-Lieberman Table-IV particle feedback with "
                "Arrhenius-threshold-only product collision power lower bound"
            ),
            supports_prediction=False,
        )


class LeeEtchProductLinearReactor:
    """One exact linear Table-IV product solve at a fixed base plasma state."""

    def __init__(self):
        self.network = build_lee_graves_lieberman_etch_product_network()
        self._network_index = {
            name: index for index, name in enumerate(self.network.species_names)
        }

    @staticmethod
    def _ion_wall_frequency_s_inv(
        condition: EtchProductPlasmaCondition,
        species: str,
    ) -> float:
        speed = bohm_speed(condition.electron_temperature_eV, _MASS_AMU[species])
        axial_area = 2.0 * math.pi * condition.geometry.radius_m ** 2
        radial_area = (
            2.0 * math.pi
            * condition.geometry.radius_m
            * condition.geometry.length_m
        )
        return float(
            speed
            * (
                condition.common_edge_factors.axial * axial_area
                + condition.common_edge_factors.radial * radial_area
            )
            / condition.geometry.volume_m3
        )

    def solve(
        self,
        condition: EtchProductPlasmaCondition,
        *,
        gross_si_removal_flux_m2_s: float,
        exposed_silicon_area_m2: float,
        sicl2_product_fraction: float = 0.6,
        sicl4_product_fraction: float = 0.4,
    ) -> EtchProductPlasmaSolution:
        if not isinstance(condition, EtchProductPlasmaCondition):
            raise TypeError("etch-product plasma condition is required")
        removal = float(gross_si_removal_flux_m2_s)
        exposed_area = float(exposed_silicon_area_m2)
        fractions = (float(sicl2_product_fraction), float(sicl4_product_fraction))
        if (
            not math.isfinite(removal) or removal < 0.0
            or not math.isfinite(exposed_area) or exposed_area <= 0.0
            or any(not math.isfinite(value) or value < 0.0 for value in fractions)
            or not math.isclose(sum(fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-14)
        ):
            raise ValueError("invalid wafer product source")

        context = RateContext(
            electron_temperature_eV=condition.electron_temperature_eV,
            gas_temperature_K=condition.gas_temperature_K,
        )
        fixed = {
            name: 0.0 for name in self.network.species_names
        }
        fixed.update({
            "e": condition.electron_density_m3,
            "Cl": condition.chlorine_atom_density_m3,
            "Cl-": condition.chlorine_negative_ion_density_m3,
        })
        matrix = np.zeros((len(_PRODUCTS), len(_PRODUCTS)))
        for column, species in enumerate(_PRODUCTS):
            densities = dict(fixed)
            densities[species] = 1.0
            source = self.network.source_vector_m3_s(densities, context)
            for row, product in enumerate(_PRODUCTS):
                matrix[row, column] = source[self._network_index[product]]

        active_fraction = condition.active_volume_fraction
        matrix *= active_fraction
        neutral_wall_frequency = {}
        ion_wall_frequency = {}
        for neutral in _NEUTRALS:
            neutral_wall_frequency[neutral] = (
                isotropic_thermal_particle_flux_m2_s(
                    1.0,
                    condition.gas_temperature_K,
                    mass_amu=_MASS_AMU[neutral],
                )
                * condition.geometry.physical_area_m2
                / condition.geometry.volume_m3
                * condition.wall_boundary.neutral_sticking_probability[neutral]
            )
            index = _PRODUCTS.index(neutral)
            matrix[index, index] -= (
                condition.exhaust_loss_frequency_s_inv
                + active_fraction * neutral_wall_frequency[neutral]
            )
        ion_sticking = condition.wall_boundary.positive_ion_sticking_probability
        for ion in _IONS:
            neutral = ion[:-1]
            frequency = self._ion_wall_frequency_s_inv(condition, ion)
            ion_wall_frequency[ion] = frequency
            ion_index = _PRODUCTS.index(ion)
            neutral_index = _PRODUCTS.index(neutral)
            matrix[ion_index, ion_index] -= active_fraction * frequency
            matrix[neutral_index, ion_index] += (
                active_fraction * (1.0 - ion_sticking) * frequency
            )

        injected_si_atom_s = removal * exposed_area
        external = np.zeros(len(_PRODUCTS))
        external[_PRODUCTS.index("SiCl2")] = (
            sicl2_product_fraction
            * injected_si_atom_s
            / condition.neutral_control_volume_m3
        )
        external[_PRODUCTS.index("SiCl4")] = (
            sicl4_product_fraction
            * injected_si_atom_s
            / condition.neutral_control_volume_m3
        )
        densities_vector = np.linalg.solve(matrix, -external)
        scale = max(float(np.max(densities_vector)), 1.0)
        if np.min(densities_vector) < -1.0e-11 * scale:
            raise RuntimeError("etch-product linear solve produced negative density")
        densities_vector = np.maximum(densities_vector, 0.0)
        residual = matrix @ densities_vector + external
        residual_scale = np.maximum(
            np.abs(matrix) @ densities_vector + np.abs(external),
            1.0,
        )
        maximum_relative_residual = float(np.max(np.abs(residual) / residual_scale))
        densities = dict(zip(_PRODUCTS, densities_vector))

        all_densities = dict(fixed)
        all_densities.update(densities)
        network_source = self.network.source_vector_m3_s(
            all_densities, context)
        chlorine_atom_source = float(
            active_fraction
            * network_source[self._network_index["Cl"]]
        )
        event_rates = self.network.event_rates_m3_s(all_densities, context)
        threshold_loss_eV_m3_s = 0.0
        for reaction, event_rate in zip(self.network.reactions, event_rates):
            coefficient = reaction.rate_coefficient
            if isinstance(coefficient, ElectronArrheniusRateCoefficient):
                threshold_loss_eV_m3_s += (
                    coefficient.activation_eV * event_rate)
        threshold_power_lower_bound = float(
            active_fraction * E_CHARGE_C * threshold_loss_eV_m3_s)

        pumping_loss = {
            name: (
                condition.exhaust_loss_frequency_s_inv
                * densities[name]
                * condition.neutral_control_volume_m3
            )
            for name in _NEUTRALS
        }
        wall_loss = {
            name: (
                active_fraction
                * neutral_wall_frequency[name]
                * densities[name]
                * condition.neutral_control_volume_m3
            )
            for name in _NEUTRALS
        }
        wall_loss.update({
            ion: (
                active_fraction
                * ion_sticking
                * ion_wall_frequency[ion]
                * densities[ion]
                * condition.neutral_control_volume_m3
            )
            for ion in _IONS
        })
        positive_ion_wall_loss = {
            ion: (
                active_fraction
                * ion_wall_frequency[ion]
                * densities[ion]
                * condition.neutral_control_volume_m3
            )
            for ion in _IONS
        }
        silicon_loss = sum(pumping_loss.values()) + sum(wall_loss.values())
        silicon_residual = (
            injected_si_atom_s - silicon_loss
        ) / max(injected_si_atom_s, silicon_loss, 1.0)
        neutral_flux = {
            name: isotropic_thermal_particle_flux_m2_s(
                densities[name],
                condition.gas_temperature_K,
                mass_amu=_MASS_AMU[name],
            )
            for name in _NEUTRALS
        }
        ion_flux = {
            ion: (
                densities[ion]
                * condition.common_edge_factors.axial
                * bohm_speed(
                    condition.electron_temperature_eV,
                    _MASS_AMU[ion],
                )
            )
            for ion in _IONS
        }
        return EtchProductPlasmaSolution(
            densities_m3=densities,
            wafer_neutral_flux_m2_s=neutral_flux,
            wafer_positive_ion_flux_m2_s=ion_flux,
            pumping_loss_molecule_s=pumping_loss,
            wall_deposition_loss_molecule_s=wall_loss,
            positive_ion_wall_loss_molecule_s=positive_ion_wall_loss,
            chlorine_atom_source_m3_s=chlorine_atom_source,
            table4_threshold_power_lower_bound_W_m3=(
                threshold_power_lower_bound),
            injected_si_atom_s=injected_si_atom_s,
            silicon_inventory_relative_residual=silicon_residual,
            linear_balance_maximum_relative_residual=maximum_relative_residual,
            matrix_condition_number=float(np.linalg.cond(matrix)),
        )
