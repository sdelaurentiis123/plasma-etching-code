"""Deterministic Oxford-NPG80 radial ion-flux closure for the Zhu condition.

The target reactor is a 13.56 MHz parallel-plate CCP with a 240 mm powered
electrode.  Zhao et al., Phys. Rev. Lett. 122, 185002 (2019), measured nearly
uniform fundamental current and nearly identical center/edge sheath motion at
13.56 MHz in a comparably sized 3 Pa CCP.  Oxford specifies a full-area gas
showerhead.  The first spatial rung therefore uses a volume-uniform source
moment and lets ambipolar diffusion, electronegative sheath-edge depletion,
and finite sidewall loss predict the radial partition.

The operator does not refit the 0-D axial flux.  It consumes the conserved
positive-ion inventory, electron state, and Lee--Lieberman inputs, then grades
its finite-volume lower-endcap flux against the independently computed 0-D
flux.  Ion mobility is anchored to Basurto and de Urquijo's mass-resolved
CHF2+ in CHF3 measurement.  This is a deterministic, differentiable spatial
closure; it is not a serial-number-specific chamber calibration.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from petch.sheath import bohm_speed

from .axisymmetric_reaction_diffusion import (
    AxisymmetricFiniteVolumeGrid,
    AxisymmetricInventoryLiftResult,
    DeterministicAxisymmetricInventoryLift,
)
from .chf3_ion_mobility import (
    CHF2MobilityState,
    load_basurto_2002_chf2_chf3_mobility_model,
)
from .geometry import CylindricalReactor


@dataclass(frozen=True)
class PositiveIonIdentity:
    mass_amu: float
    charge_number: int = 1

    def __post_init__(self):
        if (
            not math.isfinite(float(self.mass_amu))
            or float(self.mass_amu) <= 0.0
            or int(self.charge_number) != self.charge_number
            or int(self.charge_number) <= 0
        ):
            raise ValueError("invalid positive-ion identity")
        object.__setattr__(self, "mass_amu", float(self.mass_amu))
        object.__setattr__(self, "charge_number", int(self.charge_number))


ZHU_POSITIVE_ION_IDENTITIES = MappingProxyType({
    "H+": PositiveIonIdentity(1.0),
    "H2+": PositiveIonIdentity(2.0),
    "CH+": PositiveIonIdentity(13.0),
    "O+": PositiveIonIdentity(16.0),
    "F+": PositiveIonIdentity(19.0),
    "CF+": PositiveIonIdentity(31.0),
    "CHF+": PositiveIonIdentity(32.0),
    "O2+": PositiveIonIdentity(32.0),
    "S+": PositiveIonIdentity(32.0),
    "F2+": PositiveIonIdentity(38.0),
    "CF2+": PositiveIonIdentity(50.0),
    "CHF2+": PositiveIonIdentity(51.0),
    "SF+": PositiveIonIdentity(51.0),
    "CF3+": PositiveIonIdentity(69.0),
    "SF2+": PositiveIonIdentity(70.0),
    "SF2++": PositiveIonIdentity(70.0, 2),
    "SF3+": PositiveIonIdentity(89.0),
    "SF4+": PositiveIonIdentity(108.0),
    "SF4++": PositiveIonIdentity(108.0, 2),
    "SF5+": PositiveIonIdentity(127.0),
})


def _positive_mapping(values: Mapping[str, float], *, name: str):
    converted = {str(species): float(value) for species, value in values.items()}
    if (
        not converted
        or not set(converted).issubset(ZHU_POSITIVE_ION_IDENTITIES)
        or any(
            not species or not math.isfinite(value) or value < 0.0
            for species, value in converted.items()
        )
        or sum(converted.values()) <= 0.0
    ):
        raise ValueError(f"invalid {name}")
    return MappingProxyType(converted)


@dataclass(frozen=True)
class ZhuAxisymmetricTransportInput:
    condition_id: str
    geometry: CylindricalReactor
    positive_ion_density_m3: Mapping[str, float]
    global_axial_positive_ion_flux_m2_s: Mapping[str, float]
    electron_density_m3: float
    electronegativity: float
    mean_electron_energy_eV: float
    total_neutral_density_m3: float
    ion_temperature_eV: float
    ion_momentum_mean_free_path_m: float
    source: str

    def __post_init__(self):
        density = _positive_mapping(
            self.positive_ion_density_m3, name="positive-ion density")
        flux = _positive_mapping(
            self.global_axial_positive_ion_flux_m2_s,
            name="global axial positive-ion flux",
        )
        values = np.asarray([
            self.electron_density_m3,
            self.electronegativity,
            self.mean_electron_energy_eV,
            self.total_neutral_density_m3,
            self.ion_temperature_eV,
            self.ion_momentum_mean_free_path_m,
        ], dtype=float)
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.geometry, CylindricalReactor)
            or set(density) != set(flux)
            or np.any(~np.isfinite(values))
            or self.electron_density_m3 <= 0.0
            or self.electronegativity < 0.0
            or np.any(values[2:] <= 0.0)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid Zhu axisymmetric transport input")
        object.__setattr__(self, "positive_ion_density_m3", density)
        object.__setattr__(self, "global_axial_positive_ion_flux_m2_s", flux)

    @property
    def electron_temperature_eV(self) -> float:
        return (2.0 / 3.0) * float(self.mean_electron_energy_eV)

    @property
    def total_positive_ion_density_m3(self) -> float:
        return float(sum(self.positive_ion_density_m3.values()))

    @property
    def effective_bohm_mass_amu(self) -> float:
        """Mass reproducing the density-weighted multi-ion Bohm current."""
        total = self.total_positive_ion_density_m3
        inverse_sqrt_mass = sum(
            density
            * math.sqrt(ZHU_POSITIVE_ION_IDENTITIES[name].charge_number)
            / math.sqrt(ZHU_POSITIVE_ION_IDENTITIES[name].mass_amu)
            for name, density in self.positive_ion_density_m3.items()
        ) / total
        return 1.0 / inverse_sqrt_mass ** 2


@dataclass(frozen=True)
class ZhuAxisymmetricWaferIonResult:
    input_state: ZhuAxisymmetricTransportInput
    lift_result: AxisymmetricInventoryLiftResult
    reference_mobility: CHF2MobilityState
    radial_center_m: np.ndarray
    total_lower_endcap_flux_m2_s: np.ndarray
    species_lower_endcap_flux_m2_s: Mapping[str, np.ndarray]
    full_electrode_average_flux_m2_s: float
    global_model_average_flux_m2_s: float
    global_to_spatial_relative_residual: float
    optic_radius_m: float
    optic_average_flux_m2_s: float
    optic_to_full_electrode_flux_ratio: float
    source_moment: str
    source_moment_provenance: str
    supports_conditioned_radial_partition: bool = True
    supports_absolute_target_wafer_flux: bool = False
    supports_implicit_differentiation: bool = True

    def __post_init__(self):
        if (
            not isinstance(self.input_state, ZhuAxisymmetricTransportInput)
            or not isinstance(self.lift_result, AxisymmetricInventoryLiftResult)
            or not isinstance(self.reference_mobility, CHF2MobilityState)
            or not str(self.source_moment).strip()
            or not str(self.source_moment_provenance).strip()
            or not self.supports_conditioned_radial_partition
            or self.supports_absolute_target_wafer_flux
            or not self.supports_implicit_differentiation
        ):
            raise ValueError("invalid Zhu axisymmetric wafer-ion result")
        radial = np.asarray(self.radial_center_m, dtype=float).copy()
        total = np.asarray(
            self.total_lower_endcap_flux_m2_s, dtype=float).copy()
        species = {
            str(name): np.asarray(value, dtype=float).copy()
            for name, value in self.species_lower_endcap_flux_m2_s.items()
        }
        scalars = np.asarray([
            self.full_electrode_average_flux_m2_s,
            self.global_model_average_flux_m2_s,
            abs(self.global_to_spatial_relative_residual),
            self.optic_radius_m,
            self.optic_average_flux_m2_s,
            self.optic_to_full_electrode_flux_ratio,
        ])
        if (
            radial.ndim != 1
            or total.shape != radial.shape
            or set(species) != set(self.input_state.positive_ion_density_m3)
            or any(value.shape != radial.shape for value in species.values())
            or np.any(~np.isfinite(radial))
            or np.any(~np.isfinite(total))
            or np.any(total < 0.0)
            or any(np.any(~np.isfinite(value)) or np.any(value < 0.0)
                   for value in species.values())
            or np.any(~np.isfinite(scalars))
            or np.any(scalars[:2] <= 0.0)
            or not 0.0 <= abs(self.global_to_spatial_relative_residual) < 0.2
            or not 0.0 < self.optic_radius_m <= self.input_state.geometry.radius_m
            or self.optic_average_flux_m2_s <= 0.0
            or self.optic_to_full_electrode_flux_ratio <= 0.0
            or not np.allclose(
                sum(species.values()), total, rtol=2.0e-14, atol=0.0)
        ):
            raise ValueError("Zhu axisymmetric wafer-ion conservation gate failed")
        radial.setflags(write=False)
        total.setflags(write=False)
        for value in species.values():
            value.setflags(write=False)
        object.__setattr__(self, "radial_center_m", radial)
        object.__setattr__(self, "total_lower_endcap_flux_m2_s", total)
        object.__setattr__(
            self, "species_lower_endcap_flux_m2_s",
            MappingProxyType(species),
        )


class DeterministicZhuAxisymmetricCCPTransport:
    """Lift the conserved 0-D Zhu state into a radial lower-electrode flux."""

    def __init__(
        self,
        *,
        grid: AxisymmetricFiniteVolumeGrid,
        mobility_reduced_field_Td: float = 50.0,
    ):
        if not isinstance(grid, AxisymmetricFiniteVolumeGrid):
            raise TypeError("an axisymmetric grid is required")
        field = float(mobility_reduced_field_Td)
        if not math.isfinite(field) or field <= 0.0:
            raise ValueError("invalid mobility reference field")
        self.grid = grid
        self.mobility_reduced_field_Td = field
        self.mobility_model = load_basurto_2002_chf2_chf3_mobility_model()

    def predict(
        self,
        state: ZhuAxisymmetricTransportInput,
        *,
        optic_radius_m: float,
    ) -> ZhuAxisymmetricWaferIonResult:
        if state.geometry != self.grid.geometry:
            raise ValueError("Zhu state and axisymmetric grid geometry differ")
        mobility = self.mobility_model.evaluate(
            reduced_field_Td=self.mobility_reduced_field_Td,
            total_neutral_density_m3=state.total_neutral_density_m3,
        )
        temperature = state.electron_temperature_eV
        diffusion = mobility.actual_mobility_m2_V_s * temperature
        edge = state.geometry.electronegative_edge_factors(
            electronegativity=state.electronegativity,
            electron_to_ion_temperature_ratio=(
                temperature / state.ion_temperature_eV),
            ion_mean_free_path_m=state.ion_momentum_mean_free_path_m,
            include_high_pressure_diffusion=False,
        )
        bohm = bohm_speed(temperature, state.effective_bohm_mass_amu)
        wall_velocity = np.asarray([[
            edge.axial * bohm,
            edge.axial * bohm,
            edge.radial * bohm,
        ]])
        shape = np.ones((
            1,
            self.grid.radial_cell_count,
            self.grid.axial_cell_count,
        ))
        lift = DeterministicAxisymmetricInventoryLift(
            grid=self.grid,
            species_names=("positive-ion mixture",),
            diffusion_coefficient_m2_s=np.asarray([diffusion]),
            volume_reaction_matrix_s_inv=np.zeros((1, 1)),
            wall_velocity_m_s=wall_velocity,
            source_shape=shape,
            source=(
                state.source
                + "; uniform 13.56 MHz CCP source moment; Basurto CHF2+ "
                "mobility; Lee-Lieberman electronegative edge factors"
            ),
        )
        result = lift.solve(np.asarray([
            state.total_positive_ion_density_m3]))
        solution = result.solution
        total_profile = solution.lower_endcap_flux_m2_s[0]
        full_average = solution.lower_endcap_area_average_flux_m2_s(
            "positive-ion mixture",
            wafer_radius_m=state.geometry.radius_m,
        )
        optic_average = solution.lower_endcap_area_average_flux_m2_s(
            "positive-ion mixture",
            wafer_radius_m=float(optic_radius_m),
        )
        global_average = float(sum(
            state.global_axial_positive_ion_flux_m2_s.values()))
        fractions = {
            name: flux / global_average
            for name, flux in state.global_axial_positive_ion_flux_m2_s.items()
        }
        return ZhuAxisymmetricWaferIonResult(
            input_state=state,
            lift_result=result,
            reference_mobility=mobility,
            radial_center_m=self.grid.radial_centers_m,
            total_lower_endcap_flux_m2_s=total_profile,
            species_lower_endcap_flux_m2_s={
                name: fraction * total_profile
                for name, fraction in fractions.items()
            },
            full_electrode_average_flux_m2_s=full_average,
            global_model_average_flux_m2_s=global_average,
            global_to_spatial_relative_residual=(
                (full_average - global_average) / global_average),
            optic_radius_m=float(optic_radius_m),
            optic_average_flux_m2_s=optic_average,
            optic_to_full_electrode_flux_ratio=optic_average / full_average,
            source_moment="volume-uniform fundamental CCP source",
            source_moment_provenance=(
                "Zhao et al. PRL 122, 185002 measured radially uniform "
                "13.56 MHz fundamental current and nearly identical radial "
                "sheath motion in a comparable 3 Pa parallel-plate CCP; "
                "Oxford specifies a full-area showerhead"
            ),
        )
