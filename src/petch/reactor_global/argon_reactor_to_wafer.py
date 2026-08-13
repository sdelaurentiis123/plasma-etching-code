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
is a physics-complete pure-Ar sensitivity path to the resolved ion boundary,
but it is not promoted to equipment or feature-depth prediction.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .argon import ARGON_MASS_AMU
from .collisional_sheath import (
    DeterministicArgonCollisionalSheathTransfer,
    PowerClosedArgonCollisionalSheathProjection,
)
from .model import (
    ArgonGlobalCondition,
    ArgonGlobalSolution,
    LeeLiebermanArgonGlobalModel,
)
from .network import ELECTRON_MASS_KG
from .transport import ATOMIC_MASS_UNIT_KG
from .wafer_sheath_transfer import DiagnosticConditionedRFSheathTransfer


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
