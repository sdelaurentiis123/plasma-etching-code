"""Quasi-steady sheath propagation of evidence-gated self-bias histories.

Chamber conditioning evolves on seconds-to-minutes while the 13.56 MHz sheath
period is tens of nanoseconds.  This scale separation permits deterministic
RF-sheath solutions at fixed slow-time nodes.  The present operator is an
explicit collisionless baseline: it refuses collisional or absolute-depth
status until molecular ion-neutral cross sections and target fluxes are
supplied.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .oxford80_self_bias import Oxford80SelfBiasTransfer
from .wafer_sheath_transfer import (
    DiagnosticConditionedRFSheathTransfer,
    PowerClosedRFSheathProjection,
)


BOLTZMANN_J_K = 1.380649e-23


@dataclass(frozen=True)
class QuasiSteadySelfBiasSheathSnapshot:
    history_name: str
    slow_time_s: float
    normalized_time: float
    bias_magnitude_V: float
    projection: PowerClosedRFSheathProjection

    def __post_init__(self):
        if (
            not str(self.history_name).strip()
            or not math.isfinite(self.slow_time_s)
            or self.slow_time_s < 0.0
            or not math.isfinite(self.normalized_time)
            or not 0.0 <= self.normalized_time <= 1.0
            or not math.isfinite(self.bias_magnitude_V)
            or self.bias_magnitude_V <= 0.0
            or not isinstance(self.projection, PowerClosedRFSheathProjection)
            or not math.isclose(
                self.projection.bias_dc_component_v,
                self.bias_magnitude_V,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ):
            raise ValueError("invalid quasi-steady self-bias sheath snapshot")


@dataclass(frozen=True)
class CollisionlessSelfBiasSheathEnsemble:
    snapshots_by_history: Mapping[
        str, tuple[QuasiSteadySelfBiasSheathSnapshot, ...]
    ]
    target_pressure_Pa: float
    neutral_gas_temperature_K: float
    neutral_number_density_m3: float
    time_scale_separation_ratio: float
    molecular_collision_cross_sections_supplied: bool = False
    species_resolved_target_fluxes_measured_or_validated: bool = False
    supports_collisional_target_iead: bool = False
    supports_absolute_depth_prediction: bool = False

    def __post_init__(self):
        snapshots = {
            str(name): tuple(values)
            for name, values in self.snapshots_by_history.items()
        }
        if (
            not snapshots
            or any(not name or not values for name, values in snapshots.items())
            or any(
                item.history_name != name
                for name, values in snapshots.items()
                for item in values
            )
            or any(
                np.any(np.diff([item.slow_time_s for item in values]) <= 0.0)
                for values in snapshots.values()
            )
            or not math.isfinite(self.target_pressure_Pa)
            or self.target_pressure_Pa <= 0.0
            or not math.isfinite(self.neutral_gas_temperature_K)
            or self.neutral_gas_temperature_K <= 0.0
            or not math.isfinite(self.neutral_number_density_m3)
            or self.neutral_number_density_m3 <= 0.0
            or not math.isfinite(self.time_scale_separation_ratio)
            or self.time_scale_separation_ratio <= 1.0e6
            or self.molecular_collision_cross_sections_supplied
            or self.species_resolved_target_fluxes_measured_or_validated
            or self.supports_collisional_target_iead
            or self.supports_absolute_depth_prediction
        ):
            raise ValueError("invalid collisionless self-bias sheath ensemble")
        object.__setattr__(
            self, "snapshots_by_history", MappingProxyType(snapshots))


def build_collisionless_self_bias_sheath_ensemble(
    *,
    bias_transfer: Oxford80SelfBiasTransfer,
    sheath_transfer: DiagnosticConditionedRFSheathTransfer,
    positive_ion_flux_m2_s: Mapping[str, float],
    electron_temperature_eV: float,
    electron_density_m3: float,
    neutral_gas_temperature_K: float,
    normalized_time_nodes: np.ndarray | tuple[float, ...] = (
        0.0, 0.25, 0.5, 0.75, 1.0),
) -> CollisionlessSelfBiasSheathEnsemble:
    """Project every voltage history at deterministic slow-time nodes.

    The species flux argument is an upstream reactor boundary, not inferred by
    this function.  Passing a development-model flux therefore keeps the
    returned validation flags false even though energy and power ledgers close.
    """

    nodes = np.asarray(normalized_time_nodes, dtype=float)
    duration = bias_transfer.target.duration_s
    if (
        duration is None
        or nodes.ndim != 1
        or nodes.size < 2
        or np.any(~np.isfinite(nodes))
        or nodes[0] != 0.0
        or nodes[-1] != 1.0
        or np.any(np.diff(nodes) <= 0.0)
    ):
        raise ValueError("normalized time nodes must increase from zero to one")
    temperature = float(neutral_gas_temperature_K)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("neutral gas temperature must be positive")

    slow_time = nodes * duration
    snapshots = {}
    for history in bias_transfer.histories:
        values = []
        for normalized, time_s in zip(nodes, slow_time):
            bias = history.at(float(time_s))
            projection = sheath_transfer.project_from_bias_dc_component(
                positive_ion_flux_m2_s=positive_ion_flux_m2_s,
                electron_temperature_eV=electron_temperature_eV,
                electron_density_m3=electron_density_m3,
                bias_dc_component_v=bias,
            )
            values.append(QuasiSteadySelfBiasSheathSnapshot(
                history_name=history.name,
                slow_time_s=float(time_s),
                normalized_time=float(normalized),
                bias_magnitude_V=bias,
                projection=projection,
            ))
        snapshots[history.name] = tuple(values)

    pressure_pa = 0.133322368 * bias_transfer.target.pressure_mTorr
    number_density = pressure_pa / (BOLTZMANN_J_K * temperature)
    rf_period_s = 1.0 / sheath_transfer.frequency_hz
    return CollisionlessSelfBiasSheathEnsemble(
        snapshots_by_history=snapshots,
        target_pressure_Pa=pressure_pa,
        neutral_gas_temperature_K=temperature,
        neutral_number_density_m3=number_density,
        time_scale_separation_ratio=duration / rf_period_s,
    )
