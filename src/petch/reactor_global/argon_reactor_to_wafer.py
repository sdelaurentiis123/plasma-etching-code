"""Deterministic pure-Ar reactor-knob to wafer-boundary stack.

This module composes existing, separately audited operators without hiding an
equipment calibration:

``absorbed bulk power, p, Tg, geometry`` -> global Ar state and Bohm flux
``bias power, frequency`` -> power-closed collisionless RF sheath
``p, Tg`` -> implicit deterministic elastic/CX discrete-ordinates lift

The input powers are absorbed/delivered powers, not generator forward powers.
The ion collision order is closed by a bounded absorbing-system solve, but the
Maxwellian floating-potential relation is a published-model closure, not a
measurement of a driven CCP plasma potential.  Fast neutrals after their next
neutral--neutral collision also remain unresolved.  Consequently this stack
is a source-grounded pure-Ar sensitivity path for the resolved ion boundary,
not a full reactor/sheath model, and is not promoted to equipment or
feature-depth prediction.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .argon import ARGON_MASS_AMU
from .collisional_sheath import (
    ArgonBornMayerPhelpsCollisionModel,
    DeterministicArgonCollisionalSheathTransfer,
    DeterministicCollisionalSheathSolution,
    PowerClosedArgonCollisionalSheathProjection,
)
from .current_driven_rf_sheath import (
    PeriodicCurrentDensity,
    TurnerChabertCurrentDrivenSheath,
)
from .model import (
    ArgonGlobalCondition,
    ArgonGlobalSolution,
    LeeLiebermanArgonGlobalModel,
)
from .network import ELECTRON_MASS_KG
from .moving_collisional_sheath_discrete_ordinates import (
    DeterministicMovingCollisionalRFSheath,
)
from .transport import ATOMIC_MASS_UNIT_KG
from .wafer_sheath_transfer import DiagnosticConditionedRFSheathTransfer
from ..sheath import bohm_speed


def maxwellian_floating_sheath_potential_eV(
    electron_temperature_eV: float,
    ion_mass_amu: float = ARGON_MASS_AMU,
) -> float:
    """Return the collisionless Maxwellian electron/ion current-balance drop."""
    temperature = float(electron_temperature_eV)
    mass = float(ion_mass_amu)
    if (
        not math.isfinite(temperature)
        or temperature <= 0.0
        or not math.isfinite(mass)
        or mass <= 0.0
    ):
        raise ValueError("positive electron temperature and ion mass required")
    ratio = mass * ATOMIC_MASS_UNIT_KG / (
        2.0 * np.pi * ELECTRON_MASS_KG)
    return float(0.5 * temperature * np.log(ratio))


@dataclass(frozen=True)
class ArgonReactorToWaferCondition:
    global_condition: ArgonGlobalCondition
    delivered_bias_power_W: float
    bias_frequency_hz: float
    sheath_collapse_fraction: float = 1.0
    plasma_potential_source: str = (
        "collisionless Maxwellian electron/Bohm-ion floating current balance")
    plasma_potential_evidence: str = "published_model"

    def __post_init__(self):
        values = np.asarray([
            self.delivered_bias_power_W,
            self.bias_frequency_hz,
            self.sheath_collapse_fraction,
        ], dtype=float)
        if (
            not isinstance(self.global_condition, ArgonGlobalCondition)
            or np.any(~np.isfinite(values))
            or self.delivered_bias_power_W < 0.0
            or self.bias_frequency_hz <= 0.0
            or not 0.0 <= self.sheath_collapse_fraction <= 1.0
            or not str(self.plasma_potential_source).strip()
            or self.plasma_potential_evidence not in {
                "measured", "validated_model", "published_model",
                "assumed", "sensitivity",
            }
        ):
            raise ValueError("invalid Ar reactor-to-wafer condition")


@dataclass(frozen=True)
class ArgonReactorToWaferSolution:
    condition_id: str
    global_plasma: ArgonGlobalSolution
    plasma_potential_eV: float
    wafer: PowerClosedArgonCollisionalSheathProjection
    maximum_conservation_residual: float
    provenance: Mapping[str, object]
    supports_equipment_prediction: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.global_plasma, ArgonGlobalSolution)
            or not isinstance(
                self.wafer, PowerClosedArgonCollisionalSheathProjection)
            or not math.isfinite(self.plasma_potential_eV)
            or self.plasma_potential_eV <= 0.0
            or not math.isfinite(self.maximum_conservation_residual)
            or self.maximum_conservation_residual < 0.0
            or self.maximum_conservation_residual > 2.0e-8
            or self.supports_equipment_prediction
            or self.supports_feature_depth
        ):
            raise ValueError("invalid Ar reactor-to-wafer solution")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def boundary(self):
        return self.wafer.to_boundary_state()


@dataclass(frozen=True)
class DeterministicArgonReactorToWaferModel:
    global_model: LeeLiebermanArgonGlobalModel
    collisional_transfer: DeterministicArgonCollisionalSheathTransfer
    sheath_phase_count: int = 48
    sheath_steps_per_period: int = 256
    sheath_steps_per_transit: int = 256

    def __post_init__(self):
        if (
            not isinstance(self.global_model, LeeLiebermanArgonGlobalModel)
            or not isinstance(
                self.collisional_transfer,
                DeterministicArgonCollisionalSheathTransfer,
            )
            or int(self.sheath_phase_count) < 16
            or int(self.sheath_steps_per_period) < 32
            or int(self.sheath_steps_per_transit) < 32
        ):
            raise ValueError("invalid deterministic Ar reactor-to-wafer model")

    def solve(
        self,
        condition: ArgonReactorToWaferCondition,
    ) -> ArgonReactorToWaferSolution:
        if not isinstance(condition, ArgonReactorToWaferCondition):
            raise TypeError("an Ar reactor-to-wafer condition is required")
        plasma = self.global_model.solve(condition.global_condition)
        potential = maxwellian_floating_sheath_potential_eV(
            plasma.electron_temperature_eV, ARGON_MASS_AMU)
        electrode_area = np.pi * condition.global_condition.geometry.radius_m ** 2
        collisionless = DiagnosticConditionedRFSheathTransfer(
            ion_mass_amu={"Ar+": ARGON_MASS_AMU},
            electrode_area_m2=electrode_area,
            plasma_potential_eV=potential,
            frequency_hz=condition.bias_frequency_hz,
            collapse_fraction=condition.sheath_collapse_fraction,
            phase_count=int(self.sheath_phase_count),
            steps_per_period=int(self.sheath_steps_per_period),
            steps_per_transit=int(self.sheath_steps_per_transit),
            source=(
                f"{condition.plasma_potential_source}; bias power closure"),
        ).predict(
            positive_ion_flux_m2_s={"Ar+": plasma.axial_ion_flux_m2_s},
            electron_temperature_eV=plasma.electron_temperature_eV,
            electron_density_m3=plasma.electron_density_m3,
            delivered_bias_power_W=condition.delivered_bias_power_W,
        )
        wafer = self.collisional_transfer.project(
            collisionless,
            pressure_Pa=condition.global_condition.pressure_Pa,
            gas_temperature_K=condition.global_condition.gas_temperature_K,
        )
        maximum_residual = max(
            plasma.maximum_normalized_residual,
            abs(collisionless.power_closure_relative_residual),
            wafer.collisional.probability_ledger_relative_residual,
            wafer.collisional.maximum_resolved_energy_ledger_relative_residual,
        )
        return ArgonReactorToWaferSolution(
            condition_id=condition.global_condition.condition_id,
            global_plasma=plasma,
            plasma_potential_eV=potential,
            wafer=wafer,
            maximum_conservation_residual=maximum_residual,
            provenance={
                "stack": "absorbed_knobs_to_collisional_Ar_wafer_v2",
                "absorbed_bulk_power_boundary_kind": (
                    condition.global_condition.absorbed_power_boundary_kind),
                "plasma_potential_source": (
                    condition.plasma_potential_source),
                "plasma_potential_evidence": (
                    condition.plasma_potential_evidence),
                "generator_forward_power_inversion_closed": False,
                "ion_collision_order_closed": wafer.collisional.provenance[
                    "ion_collision_order_closed"],
                "fast_neutral_wafer_transport_closed": False,
                "feature_depth_used": False,
            },
            supports_equipment_prediction=False,
            supports_feature_depth=False,
        )


@dataclass(frozen=True)
class ArgonCurrentDrivenReactorToWaferCondition:
    """Pure-Ar bulk condition plus a de-embedded sheath-current waveform.

    This is the physically identifiable seam between a plasma/electrical
    reactor model and the moving sheath.  A generator-power setpoint is not a
    substitute: its mapping to electrode current requires the matching
    network, stray impedances, and plasma load.
    """

    global_condition: ArgonGlobalCondition
    sheath_current_density: PeriodicCurrentDensity

    def __post_init__(self):
        if (
            not isinstance(self.global_condition, ArgonGlobalCondition)
            or not isinstance(
                self.sheath_current_density, PeriodicCurrentDensity)
        ):
            raise ValueError("invalid current-driven Ar reactor condition")


@dataclass(frozen=True)
class ArgonCurrentDrivenReactorToWaferSolution:
    """Audited global-plasma -> moving-sheath -> collisional-ion boundary."""

    condition_id: str
    global_plasma: ArgonGlobalSolution
    sheath_edge_density_m3: float
    sheath: TurnerChabertCurrentDrivenSheath
    wafer: DeterministicCollisionalSheathSolution
    bohm_flux_seam_relative_residual: float
    maximum_conservation_residual: float
    provenance: Mapping[str, object]
    supports_resolved_ion_boundary_prediction: bool
    supports_equipment_prediction: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        scalars = np.asarray([
            self.sheath_edge_density_m3,
            self.bohm_flux_seam_relative_residual,
            self.maximum_conservation_residual,
        ])
        expected_support = (
            self.global_plasma.supports_prediction
            and self.sheath.current.supports_predictive_boundary
            and self.maximum_conservation_residual <= 2.0e-8
        )
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.global_plasma, ArgonGlobalSolution)
            or not isinstance(self.sheath, TurnerChabertCurrentDrivenSheath)
            or not isinstance(
                self.wafer, DeterministicCollisionalSheathSolution)
            or np.any(~np.isfinite(scalars))
            or np.any(scalars < 0.0)
            or self.sheath_edge_density_m3 <= 0.0
            or self.maximum_conservation_residual > 2.0e-8
            or self.supports_resolved_ion_boundary_prediction
            != expected_support
            or self.supports_equipment_prediction
            or self.supports_feature_depth
        ):
            raise ValueError("invalid current-driven Ar reactor solution")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def boundary(self):
        return self.wafer.to_boundary_state(
            ion_name="Ar+", ion_mass_amu=ARGON_MASS_AMU)


@dataclass(frozen=True)
class DeterministicCurrentDrivenArgonReactorToWaferModel:
    """Compose the global Ar model with a moving collisional RF sheath.

    RF phase, ion direction, position, energy, and transverse-energy fraction
    are all kinetic coordinates in the wafer operator.  The sparse absorbing
    solve closes every ion elastic/CX collision order without particle Monte
    Carlo.  The remaining fast-neutral lineage and generator/network inverse
    stay explicit open gates.
    """

    global_model: LeeLiebermanArgonGlobalModel
    collision_model: ArgonBornMayerPhelpsCollisionModel
    sheath_phase_quadrature_count: int = 4096
    phase_node_count: int = 12
    position_node_count: int = 7
    total_energy_node_count: int = 7
    transverse_fraction_node_count: int = 9
    initial_thermal_radial_order: int = 2
    output_azimuth_order: int = 4
    impact_quadrature_order: int = 2
    collision_azimuth_order: int = 4
    collision_event_quadrature_order: int = 3
    steps_per_period: int = 192
    steps_per_transit: int = 192
    maximum_transit_periods: float = 30.0

    def __post_init__(self):
        orders = (
            self.sheath_phase_quadrature_count,
            self.phase_node_count,
            self.position_node_count,
            self.total_energy_node_count,
            self.transverse_fraction_node_count,
            self.initial_thermal_radial_order,
            self.output_azimuth_order,
            self.impact_quadrature_order,
            self.collision_azimuth_order,
            self.collision_event_quadrature_order,
            self.steps_per_period,
            self.steps_per_transit,
        )
        if (
            not isinstance(self.global_model, LeeLiebermanArgonGlobalModel)
            or not isinstance(
                self.collision_model, ArgonBornMayerPhelpsCollisionModel)
            or any(int(value) < 1 for value in orders)
            or int(self.sheath_phase_quadrature_count) < 256
            or int(self.phase_node_count) < 4
            or int(self.position_node_count) < 3
            or int(self.total_energy_node_count) < 3
            or int(self.transverse_fraction_node_count) < 3
            or int(self.steps_per_period) < 16
            or int(self.steps_per_transit) < 16
            or not math.isfinite(float(self.maximum_transit_periods))
            or self.maximum_transit_periods <= 0.0
        ):
            raise ValueError("invalid current-driven Ar reactor model")

    def solve(
        self,
        condition: ArgonCurrentDrivenReactorToWaferCondition,
    ) -> ArgonCurrentDrivenReactorToWaferSolution:
        if not isinstance(
            condition, ArgonCurrentDrivenReactorToWaferCondition
        ):
            raise TypeError("a current-driven Ar condition is required")
        plasma = self.global_model.solve(condition.global_condition)
        speed = bohm_speed(
            plasma.electron_temperature_eV, ARGON_MASS_AMU)
        edge_density = plasma.axial_ion_flux_m2_s / speed
        seam_residual = abs(
            edge_density * speed - plasma.axial_ion_flux_m2_s
        ) / plasma.axial_ion_flux_m2_s
        sheath = TurnerChabertCurrentDrivenSheath(
            current=condition.sheath_current_density,
            electron_temperature_eV=plasma.electron_temperature_eV,
            ion_mass_amu=ARGON_MASS_AMU,
            sheath_edge_density_m3=edge_density,
            phase_quadrature_count=int(
                self.sheath_phase_quadrature_count),
            provenance={
                "sheath_edge_density_source": (
                    "global axial Bohm flux divided by Bohm speed"),
                "global_condition_id": plasma.condition_id,
            },
        )
        wafer = DeterministicMovingCollisionalRFSheath(
            sheath=sheath,
            collision_model=self.collision_model,
            gas_number_density_m3=(
                condition.global_condition.neutral_ground_density_m3),
            neutral_gas_temperature_K=(
                condition.global_condition.gas_temperature_K),
            source_ion_flux_m2_s=plasma.axial_ion_flux_m2_s,
            phase_node_count=int(self.phase_node_count),
            position_node_count=int(self.position_node_count),
            total_energy_node_count=int(self.total_energy_node_count),
            transverse_fraction_node_count=int(
                self.transverse_fraction_node_count),
            initial_thermal_radial_order=int(
                self.initial_thermal_radial_order),
            output_azimuth_order=int(self.output_azimuth_order),
            impact_quadrature_order=int(self.impact_quadrature_order),
            collision_azimuth_order=int(self.collision_azimuth_order),
            collision_event_quadrature_order=int(
                self.collision_event_quadrature_order),
            steps_per_period=int(self.steps_per_period),
            steps_per_transit=int(self.steps_per_transit),
            maximum_transit_periods=float(self.maximum_transit_periods),
            provenance={
                "stack": "global_Ar_to_moving_collisional_wafer_v1",
                "feature_depth_used": False,
            },
        ).solve()
        maximum_residual = max(
            plasma.maximum_normalized_residual,
            seam_residual,
            sheath.child_current_relative_residual,
            sheath.charge_voltage_relative_residual,
            wafer.probability_ledger_relative_residual,
            wafer.maximum_resolved_energy_ledger_relative_residual,
        )
        resolved_support = (
            plasma.supports_prediction
            and condition.sheath_current_density.supports_predictive_boundary
            and maximum_residual <= 2.0e-8
        )
        return ArgonCurrentDrivenReactorToWaferSolution(
            condition_id=plasma.condition_id,
            global_plasma=plasma,
            sheath_edge_density_m3=edge_density,
            sheath=sheath,
            wafer=wafer,
            bohm_flux_seam_relative_residual=seam_residual,
            maximum_conservation_residual=maximum_residual,
            provenance={
                "stack": "global_Ar_to_moving_collisional_wafer_v1",
                "bulk_absorbed_power_source": plasma.absorbed_power_source,
                "bulk_absorbed_power_evidence": (
                    plasma.absorbed_power_evidence),
                "sheath_current_source": (
                    condition.sheath_current_density.source),
                "sheath_current_evidence": (
                    condition.sheath_current_density.evidence_kind),
                "sheath_edge_density_source": (
                    "global axial Bohm flux divided by Bohm speed"),
                "moving_electron_front_resolved": True,
                "time_dependent_Poisson_field_resolved": True,
                "RF_phase_is_kinetic_state": True,
                "ion_collision_order_closed": True,
                "fast_neutral_transport_closed": False,
                "low_energy_ion_angular_scattering_closed": True,
                "low_energy_ion_angular_model": (
                    "Phelps LXCat isotropic/backscatter decomposition"),
                "generator_matching_network_inversion_closed": False,
                "molecular_chemistry_closed": False,
                "feature_depth_used": False,
            },
            supports_resolved_ion_boundary_prediction=resolved_support,
            supports_equipment_prediction=False,
            supports_feature_depth=False,
        )
