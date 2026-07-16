"""External current terminals for floating conductors outside the feature mesh.

Feature electrostatics resolves micron-scale surfaces. Some experiments connect those surfaces to
millimeter-scale collectors that remain outside the modeled Poisson domain. This module represents
that missing circuit element as a conservative current contribution to an existing floating
conductor inventory. It does not add volume plasma charge, alter local particle trajectories, or
pretend that an unresolved pad is part of the feature geometry.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .charging_poisson_3d import NodalPoissonSystem3D


def _machine_value(value, path):
    if value is None or isinstance(value, (str, bool, int, float)):
        if isinstance(value, float) and not np.isfinite(value):
            raise ValueError(f"{path} contains a non-finite value")
        return value
    if isinstance(value, np.generic):
        return _machine_value(value.item(), path)
    if isinstance(value, Mapping):
        return {
            str(key): _machine_value(item, f"{path}.{key}")
            for key, item in value.items()
        }
    if isinstance(value, (tuple, list)):
        return [
            _machine_value(item, f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    raise ValueError(f"{path} contains opaque {type(value).__name__}")


@dataclass(frozen=True)
class ConductorTerminalCurrent3D:
    """Positive/negative nodal current inventories from one external terminal."""

    positive_node_current_a: np.ndarray
    negative_node_current_a: np.ndarray
    signed_current_a_by_conductor: Mapping[int, float]
    provenance: Mapping[str, object]

    def __post_init__(self):
        positive = np.asarray(self.positive_node_current_a, dtype=float).copy()
        negative = np.asarray(self.negative_node_current_a, dtype=float).copy()
        signed = {int(key): float(value)
                  for key, value in dict(self.signed_current_a_by_conductor).items()}
        if (
            positive.shape != negative.shape
            or positive.ndim != 3
            or np.any(~np.isfinite(positive))
            or np.any(~np.isfinite(negative))
            or np.any(positive < 0.0)
            or np.any(negative < 0.0)
            or not signed
            or any(key <= 0 or not np.isfinite(value) for key, value in signed.items())
        ):
            raise ValueError("invalid conductor-terminal current contribution")
        represented = float(np.sum(positive) - np.sum(negative))
        declared = float(sum(signed.values()))
        scale = max(float(np.sum(positive) + np.sum(negative)), np.finfo(float).tiny)
        if abs(represented - declared) > 2e-14 * scale:
            raise ValueError("conductor-terminal nodal and component ledgers disagree")
        positive.setflags(write=False)
        negative.setflags(write=False)
        object.__setattr__(self, "positive_node_current_a", positive)
        object.__setattr__(self, "negative_node_current_a", negative)
        object.__setattr__(
            self, "signed_current_a_by_conductor", MappingProxyType(signed))
        object.__setattr__(
            self, "provenance",
            MappingProxyType(_machine_value(dict(self.provenance), "provenance")))

    @property
    def signed_total_current_a(self):
        return float(sum(self.signed_current_a_by_conductor.values()))

    @property
    def absolute_total_current_a(self):
        return float(np.sum(self.positive_node_current_a)
                     + np.sum(self.negative_node_current_a))


@dataclass(frozen=True)
class RemotePadElectronCollector3D:
    """R1-dominant constant-electron-current limit of the Nozawa pad circuit.

    Nozawa et al. report that pad-to-plasma resistance dominates polysilicon line resistance and
    that notch depth is insensitive to a line distance of about 27 mm. In that experimentally
    stated limit, collected electron current scales with pad perimeter and is delivered to each
    connected equipotential conductor without a resolved voltage drop along the line.

    The collection coefficient remains an input with explicit source and bounds. This class is an
    engine mechanism, not a value for the Nozawa held-out cases: those cases may use it only after
    the pad perimeters and one legal calibration of the coefficient are committed.
    """

    collector_perimeter_m_by_conductor: Mapping[int, float]
    electron_current_per_perimeter_a_m: float
    coefficient_bounds_a_m: tuple[float, float]
    source: str
    coefficient_evidence: str
    topology_evidence: str

    def __post_init__(self):
        perimeter = {
            int(key): float(value)
            for key, value in dict(self.collector_perimeter_m_by_conductor).items()
        }
        bounds = tuple(float(value) for value in self.coefficient_bounds_a_m)
        coefficient = float(self.electron_current_per_perimeter_a_m)
        if (
            not perimeter
            or any(key <= 0 or not np.isfinite(value) or value <= 0.0
                   for key, value in perimeter.items())
            or not np.isfinite(coefficient)
            or coefficient <= 0.0
            or len(bounds) != 2
            or any(not np.isfinite(value) or value <= 0.0 for value in bounds)
            or not bounds[0] <= coefficient <= bounds[1]
            or bounds[1] <= bounds[0]
            or not str(self.source).strip()
            or not str(self.coefficient_evidence).strip()
            or not str(self.topology_evidence).strip()
        ):
            raise ValueError("invalid remote-pad electron collector")
        object.__setattr__(
            self, "collector_perimeter_m_by_conductor",
            MappingProxyType(perimeter),
        )
        object.__setattr__(self, "coefficient_bounds_a_m", bounds)

    @property
    def provenance(self):
        return MappingProxyType({
            "model": "remote-pad-electron-collector-r1-dominant-v1",
            "source": self.source,
            "physical_limit": (
                "constant electron source; pad-plasma resistance R1 dominates "
                "polysilicon-line resistance R2"),
            "parameters": {
                "electron_current_per_perimeter_a_m": (
                    float(self.electron_current_per_perimeter_a_m)),
                "collector_perimeter_m_by_conductor": {
                    str(key): value
                    for key, value in self.collector_perimeter_m_by_conductor.items()
                },
            },
            "bounds": {
                "electron_current_per_perimeter_a_m": list(
                    self.coefficient_bounds_a_m),
            },
            "evidence": {
                "electron_current_per_perimeter_a_m": self.coefficient_evidence,
                "collector_topology": self.topology_evidence,
            },
            "resolved_line_resistance": False,
            "external_pad_capacitance": False,
            "volume_plasma_charge_added": False,
        })

    def current_contribution(self, poisson_system: NodalPoissonSystem3D):
        """Route remote electron current into exact floating-conductor inventories."""
        if not isinstance(poisson_system, NodalPoissonSystem3D):
            raise TypeError("poisson_system must be NodalPoissonSystem3D")
        declared = set(self.collector_perimeter_m_by_conductor)
        available = set(poisson_system.floating_conductor_ids)
        if not declared.issubset(available):
            raise ValueError(
                "remote-pad conductor ids are absent from the Poisson system: "
                + ", ".join(map(str, sorted(declared - available))))
        positive = np.zeros(poisson_system.shape)
        negative = np.zeros(poisson_system.shape)
        signed = {}
        for conductor_id, perimeter_m in self.collector_perimeter_m_by_conductor.items():
            electron_current = (
                float(self.electron_current_per_perimeter_a_m) * perimeter_m)
            node = poisson_system.floating_conductor_representative_node(conductor_id)
            negative[node] += electron_current
            signed[conductor_id] = -electron_current
        return ConductorTerminalCurrent3D(
            positive, negative, signed, self.provenance)
