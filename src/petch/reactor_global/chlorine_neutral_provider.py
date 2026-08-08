"""State-dependent chlorine neutral-wall transport composition."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Protocol

import numpy as np

from .chlorine_particle_model import ChlorineFixedPressureCondition
from .chlorine_transport import (
    ChlorineNeutralWallTransport,
    solve_chlorine_neutral_wall_transport,
)
from .chlorine_wall import (
    ChlorineIncidentVelocityState,
    ChlorineWallRecombinationBoundary,
)
from .neutral_transport import (
    ChapmanEnskogBinaryDiffusivity,
    ReducedNeutralDiffusivity,
)

_CHLORINE_SPECIES = frozenset({"e", "Cl2", "Cl", "Cl2+", "Cl+", "Cl-"})


class ChlorineWallBoundaryProvider(Protocol):
    """Condition-scoped wall-response law used by neutral transport."""

    name: str
    version: str

    def predict(
        self,
        *,
        cl_to_cl2_ratio: float,
        pressure_Pa: float,
        icp_power_W: float,
        gas_temperature_K: float,
    ) -> ChlorineWallRecombinationBoundary:
        ...


@dataclass(frozen=True)
class StateDependentChlorineNeutralTransportProvider:
    """Recompute the exact neutral Robin loss from the current Cl/Cl2 state."""

    wall_recombination_provider: ChlorineWallBoundaryProvider
    incident_velocity_state: ChlorineIncidentVelocityState
    diffusivity_model: (
        ReducedNeutralDiffusivity | ChapmanEnskogBinaryDiffusivity)
    name: str = "state_dependent_chlorine_neutral_transport"
    version: str = "1"
    neutral_density_basis: str = "state_total_neutral_particles"

    def __post_init__(self):
        wall_provider = self.wall_recombination_provider
        if (
            not hasattr(wall_provider, "predict")
            or not str(getattr(wall_provider, "name", "")).strip()
            or not str(getattr(wall_provider, "version", "")).strip()
            or not isinstance(
                self.incident_velocity_state,
                ChlorineIncidentVelocityState,
            )
            or not isinstance(
                self.diffusivity_model,
                (ReducedNeutralDiffusivity, ChapmanEnskogBinaryDiffusivity),
            )
            or not str(self.name).strip()
            or not str(self.version).strip()
            or self.neutral_density_basis != "state_total_neutral_particles"
        ):
            raise ValueError("invalid state-dependent neutral provider")

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
    ) -> ChlorineNeutralWallTransport:
        if not isinstance(condition, ChlorineFixedPressureCondition):
            raise TypeError("chlorine fixed-pressure condition is required")
        densities = {
            str(name): float(value) for name, value in densities_m3.items()
        }
        if (
            set(densities) != _CHLORINE_SPECIES
            or any(
                not np.isfinite(value) or value <= 0.0
                for value in densities.values()
            )
        ):
            raise ValueError("invalid chlorine density state")
        ratio = densities["Cl"] / densities["Cl2"]
        boundary = self.wall_recombination_provider.predict(
            cl_to_cl2_ratio=ratio,
            pressure_Pa=condition.pressure.value,
            icp_power_W=condition.source_power.value,
            gas_temperature_K=condition.gas_temperature.value,
        )
        if not isinstance(boundary, ChlorineWallRecombinationBoundary):
            raise TypeError("wall provider returned an invalid boundary")
        return solve_chlorine_neutral_wall_transport(
            geometry=condition.geometry,
            wall_boundary=boundary,
            incident_velocity_state=self.incident_velocity_state,
            diffusivity_model=self.diffusivity_model,
            total_neutral_density_m3=(densities["Cl2"] + densities["Cl"]),
            gas_temperature_K=condition.gas_temperature.value,
            cl_to_cl2_ratio=ratio,
            pressure_Pa=condition.pressure.value,
            icp_power_W=condition.source_power.value,
        )
