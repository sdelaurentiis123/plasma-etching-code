"""Reusable chlorine global-state to axisymmetric wafer-ion provider.

The provider lifts a declared volume-average ``Cl+``/``Cl2+``/``Cl-`` state
through the deterministic quasineutral drift-diffusion tier.  It uses public
Lee--Economou reduced mobilities, free-ion Einstein diffusion, Bohm loss for
positive ions, an absorbing negative-ion boundary, and an explicitly supplied
source moment.  The inferred volumetric source amplitudes and all particle
ledgers remain visible.

This is deliberately a *conditional spatial provider*: it predicts how a
global inventory partitions to a finite wafer, but it does not independently
predict that inventory or the ICP power-deposition field.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from petch.sheath import bohm_speed

from .axisymmetric_drift_diffusion import (
    DeterministicQuasineutralInventoryLift,
    QuasineutralInventoryLiftResult,
    QuasineutralInventoryLiftTangent,
)
from .axisymmetric_reaction_diffusion import AxisymmetricFiniteVolumeGrid
from .chlorine import CHLORINE_ATOM_MASS_AMU, CHLORINE_MOLECULE_MASS_AMU
from .chlorine_charged_transport import (
    ReducedIonMobility,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
)
from .geometry import CylindricalReactor


CHLORINE_AXISYMMETRIC_CHARGED_SPECIES = ("Cl+", "Cl2+", "Cl-")
CHLORINE_AXISYMMETRIC_POSITIVE_IONS = ("Cl+", "Cl2+")
_CHARGE_NUMBER = np.asarray([1.0, 1.0, -1.0])
_POSITIVE_ION_MASS_AMU = MappingProxyType({
    "Cl+": CHLORINE_ATOM_MASS_AMU,
    "Cl2+": CHLORINE_MOLECULE_MASS_AMU,
})


def _exact_mapping(
    values: Mapping[str, float], expected: tuple[str, ...], *, positive: bool
) -> MappingProxyType:
    converted = {str(name): float(value) for name, value in values.items()}
    if (
        set(converted) != set(expected)
        or any(
            not math.isfinite(value)
            or (value <= 0.0 if positive else value < 0.0)
            for value in converted.values()
        )
    ):
        raise ValueError("invalid chlorine axisymmetric mapping")
    return MappingProxyType(converted)


@dataclass(frozen=True)
class ChlorineAxisymmetricTransportInput:
    """Minimum global state needed by the charged spatial tier."""

    condition_id: str
    geometry: CylindricalReactor
    charged_density_m3: Mapping[str, float]
    total_neutral_density_m3: float
    mean_electron_energy_eV: float
    source: str
    global_state_supports_prediction: bool = False

    def __post_init__(self):
        density = _exact_mapping(
            self.charged_density_m3,
            CHLORINE_AXISYMMETRIC_CHARGED_SPECIES,
            positive=True,
        )
        neutral = float(self.total_neutral_density_m3)
        energy = float(self.mean_electron_energy_eV)
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.geometry, CylindricalReactor)
            or not math.isfinite(neutral)
            or neutral <= 0.0
            or not math.isfinite(energy)
            or energy <= 0.0
            or not str(self.source).strip()
        ):
            raise ValueError("invalid chlorine axisymmetric transport input")
        object.__setattr__(self, "charged_density_m3", density)
        object.__setattr__(self, "total_neutral_density_m3", neutral)
        object.__setattr__(self, "mean_electron_energy_eV", energy)

    @property
    def electron_temperature_eV(self) -> float:
        return (2.0 / 3.0) * self.mean_electron_energy_eV

    @classmethod
    def from_eedf_solution(cls, condition, solution):
        """Adapt an EEDF chlorine solve without weakening its evidence flags."""
        if condition.condition_id != solution.condition_id:
            raise ValueError("chlorine EEDF condition/solution mismatch")
        return cls(
            condition_id=solution.condition_id,
            geometry=condition.geometry,
            charged_density_m3={
                name: solution.densities_m3[name]
                for name in CHLORINE_AXISYMMETRIC_CHARGED_SPECIES
            },
            total_neutral_density_m3=(
                solution.densities_m3["Cl"] + solution.densities_m3["Cl2"]
            ),
            mean_electron_energy_eV=solution.mean_electron_energy_eV,
            source=(
                "EEDF chlorine global solution " + solution.condition_id
            ),
            global_state_supports_prediction=bool(
                solution.supports_reactor_state_prediction
            ),
        )


@dataclass(frozen=True)
class ChlorineAxisymmetricWaferIonResult:
    input_state: ChlorineAxisymmetricTransportInput
    lift_result: QuasineutralInventoryLiftResult
    source_moment_name: str
    source_moment_provenance: str
    source_moment_measured_or_validated: bool
    wafer_radius_m: float
    wafer_positive_ion_flux_m2_s: Mapping[str, float]
    full_lower_endcap_positive_ion_flux_m2_s: Mapping[str, float]
    mobility_m2_V_s: Mapping[str, float]
    bohm_wall_velocity_m_s: Mapping[str, float]
    supports_reactor_state_prediction: bool = False
    supports_absolute_wafer_flux_prediction: bool = False
    supports_implicit_differentiation: bool = True

    def __post_init__(self):
        if (
            not isinstance(self.input_state, ChlorineAxisymmetricTransportInput)
            or not isinstance(self.lift_result, QuasineutralInventoryLiftResult)
            or not str(self.source_moment_name).strip()
            or not str(self.source_moment_provenance).strip()
            or not math.isfinite(float(self.wafer_radius_m))
            or not 0.0 < float(self.wafer_radius_m) <= self.input_state.geometry.radius_m
            or bool(self.supports_reactor_state_prediction)
            or bool(self.supports_absolute_wafer_flux_prediction)
            or not bool(self.supports_implicit_differentiation)
        ):
            raise ValueError("invalid chlorine axisymmetric wafer result")
        wafer = _exact_mapping(
            self.wafer_positive_ion_flux_m2_s,
            CHLORINE_AXISYMMETRIC_POSITIVE_IONS,
            positive=False,
        )
        endcap = _exact_mapping(
            self.full_lower_endcap_positive_ion_flux_m2_s,
            CHLORINE_AXISYMMETRIC_POSITIVE_IONS,
            positive=False,
        )
        mobility = _exact_mapping(
            self.mobility_m2_V_s,
            CHLORINE_AXISYMMETRIC_CHARGED_SPECIES,
            positive=True,
        )
        bohm = _exact_mapping(
            self.bohm_wall_velocity_m_s,
            CHLORINE_AXISYMMETRIC_POSITIVE_IONS,
            positive=True,
        )
        object.__setattr__(self, "wafer_positive_ion_flux_m2_s", wafer)
        object.__setattr__(
            self, "full_lower_endcap_positive_ion_flux_m2_s", endcap)
        object.__setattr__(self, "mobility_m2_V_s", mobility)
        object.__setattr__(self, "bohm_wall_velocity_m_s", bohm)
        object.__setattr__(self, "wafer_radius_m", float(self.wafer_radius_m))

    @property
    def total_wafer_positive_ion_flux_m2_s(self) -> float:
        return float(sum(self.wafer_positive_ion_flux_m2_s.values()))


class DeterministicChlorineAxisymmetricTransport:
    """Conditional global-inventory to finite-wafer charged transport."""

    def __init__(
        self,
        *,
        grid: AxisymmetricFiniteVolumeGrid,
        ion_temperature_eV: float,
        positive_negative_recombination_m3_s: float,
        reduced_mobilities: Mapping[str, ReducedIonMobility] | None = None,
    ):
        temperature = float(ion_temperature_eV)
        recombination = float(positive_negative_recombination_m3_s)
        mobilities = dict(
            lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities()
            if reduced_mobilities is None else reduced_mobilities
        )
        if (
            not isinstance(grid, AxisymmetricFiniteVolumeGrid)
            or not math.isfinite(temperature)
            or temperature <= 0.0
            or not math.isfinite(recombination)
            or recombination < 0.0
            or set(mobilities) != set(CHLORINE_AXISYMMETRIC_CHARGED_SPECIES)
            or any(not isinstance(value, ReducedIonMobility)
                   for value in mobilities.values())
        ):
            raise ValueError("invalid deterministic chlorine spatial provider")
        self.grid = grid
        self.ion_temperature_eV = temperature
        self.positive_negative_recombination_m3_s = recombination
        self.reduced_mobilities = MappingProxyType(mobilities)

    def _transport_arrays(self, state: ChlorineAxisymmetricTransportInput):
        mobility = np.asarray([
            self.reduced_mobilities[name].evaluate(
                total_neutral_density_m3=state.total_neutral_density_m3,
                ion_temperature_eV=self.ion_temperature_eV,
            ).mobility_m2_V_s
            for name in CHLORINE_AXISYMMETRIC_CHARGED_SPECIES
        ])
        bohm = np.asarray([
            bohm_speed(
                state.electron_temperature_eV,
                _POSITIVE_ION_MASS_AMU[name],
            )
            for name in CHLORINE_AXISYMMETRIC_POSITIVE_IONS
        ])
        return mobility, bohm

    def _build_lift(
        self,
        state: ChlorineAxisymmetricTransportInput,
        source_shape: np.ndarray,
        *,
        source_moment_name: str,
        source_moment_provenance: str,
    ):
        if (
            not isinstance(state, ChlorineAxisymmetricTransportInput)
            or state.geometry != self.grid.geometry
            or not str(source_moment_name).strip()
            or not str(source_moment_provenance).strip()
        ):
            raise ValueError("chlorine spatial input/grid mismatch")
        shape = np.asarray(source_shape, dtype=float)
        expected = (
            self.grid.radial_cell_count, self.grid.axial_cell_count)
        if shape.shape != expected:
            raise ValueError("chlorine source moment has wrong grid shape")
        mobility, bohm = self._transport_arrays(state)
        lift = DeterministicQuasineutralInventoryLift(
            grid=self.grid,
            species_names=CHLORINE_AXISYMMETRIC_CHARGED_SPECIES,
            charge_number=_CHARGE_NUMBER,
            mobility_m2_V_s=mobility,
            ion_temperature_eV=np.full(3, self.ion_temperature_eV),
            electron_temperature_eV=state.electron_temperature_eV,
            wall_velocity_m_s=np.vstack((
                np.repeat(bohm[:, None], 3, axis=1),
                np.full((1, 3), np.inf),
            )),
            source_shape=np.stack((shape, shape, shape)),
            source=(
                state.source + "; source moment " + source_moment_name
                + ": " + source_moment_provenance
            ),
            positive_negative_recombination_m3_s=(
                self.positive_negative_recombination_m3_s),
        )
        return lift, mobility, bohm

    def predict(
        self,
        state: ChlorineAxisymmetricTransportInput,
        source_shape: np.ndarray,
        *,
        source_moment_name: str,
        source_moment_provenance: str,
        source_moment_measured_or_validated: bool,
        wafer_radius_m: float,
        initial_electrostatic_potential_V: np.ndarray | None = None,
        relative_tolerance: float = 5.0e-8,
        maximum_iterations: int = 500,
    ) -> ChlorineAxisymmetricWaferIonResult:
        lift, mobility, bohm = self._build_lift(
            state,
            source_shape,
            source_moment_name=source_moment_name,
            source_moment_provenance=source_moment_provenance,
        )
        target = np.asarray([
            state.charged_density_m3[name]
            for name in CHLORINE_AXISYMMETRIC_CHARGED_SPECIES
        ])
        result = lift.solve(
            target,
            relative_tolerance=relative_tolerance,
            maximum_iterations=maximum_iterations,
            initial_electrostatic_potential_V=(
                initial_electrostatic_potential_V),
        )
        radius = float(wafer_radius_m)
        wafer = {
            name: result.solution.lower_endcap_area_average_flux_m2_s(
                name, wafer_radius_m=radius)
            for name in CHLORINE_AXISYMMETRIC_POSITIVE_IONS
        }
        endcap = {
            name: result.solution.lower_endcap_area_average_flux_m2_s(
                name, wafer_radius_m=self.grid.geometry.radius_m)
            for name in CHLORINE_AXISYMMETRIC_POSITIVE_IONS
        }
        return ChlorineAxisymmetricWaferIonResult(
            input_state=state,
            lift_result=result,
            source_moment_name=source_moment_name,
            source_moment_provenance=source_moment_provenance,
            source_moment_measured_or_validated=bool(
                source_moment_measured_or_validated),
            wafer_radius_m=radius,
            wafer_positive_ion_flux_m2_s=wafer,
            full_lower_endcap_positive_ion_flux_m2_s=endcap,
            mobility_m2_V_s=dict(zip(
                CHLORINE_AXISYMMETRIC_CHARGED_SPECIES, mobility)),
            bohm_wall_velocity_m_s=dict(zip(
                CHLORINE_AXISYMMETRIC_POSITIVE_IONS, bohm)),
        )

    def target_inventory_jvp(
        self,
        result: ChlorineAxisymmetricWaferIonResult,
        target_volume_average_density_tangent_m3: np.ndarray,
    ) -> tuple[QuasineutralInventoryLiftTangent, Mapping[str, float]]:
        """Return exact implicit field and finite-wafer flux derivatives."""
        if not isinstance(result, ChlorineAxisymmetricWaferIonResult):
            raise TypeError("a chlorine axisymmetric wafer result is required")
        source_shape = result.lift_result.solution.condition.source_rate_m3_s
        source_shape = source_shape / (
            result.lift_result.inferred_source_amplitude_m3_s[:, None, None]
        )
        # All three species share the declared source moment by construction.
        lift, _, _ = self._build_lift(
            result.input_state,
            source_shape[0],
            source_moment_name=result.source_moment_name,
            source_moment_provenance=result.source_moment_provenance,
        )
        tangent = lift.implicit_target_inventory_jvp(
            result.lift_result,
            target_volume_average_density_tangent_m3,
        )
        solution = result.lift_result.solution
        radius = result.wafer_radius_m
        outer = np.minimum(self.grid.radial_edges_m[1:], radius)
        inner = np.minimum(self.grid.radial_edges_m[:-1], radius)
        area = np.pi * np.maximum(outer ** 2 - inner ** 2, 0.0)
        wafer_area = np.pi * radius ** 2
        flux_tangent = {}
        for index, name in enumerate(CHLORINE_AXISYMMETRIC_POSITIVE_IONS):
            lower_density = solution.density_m3[index, :, 0]
            lower_flux = solution.lower_endcap_flux_m2_s[index]
            effective_velocity = np.divide(
                lower_flux,
                lower_density,
                out=np.zeros_like(lower_flux),
                where=lower_density > 0.0,
            )
            flux_tangent[name] = float(np.dot(
                effective_velocity * tangent.density_tangent_m3[index, :, 0],
                area,
            ) / wafer_area)
        return tangent, MappingProxyType(flux_tangent)
