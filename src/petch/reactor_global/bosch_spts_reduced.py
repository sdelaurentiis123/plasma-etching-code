"""Deterministic measured-waveform reactor-to-wafer tier for SPTS Bosch etching.

The Sayyed et al. dataset records one SPTS Omega i2L DSi Rapier run at 5 Hz.
The tool has independent ICP source and platen-bias channels; its second ICP
source is present in the trace schema but is zero throughout this experiment.
This module turns the measured gas, pressure, net source power, and platen Vpp
into a compact phase-resolved plasma state and then uses the common
axisymmetric finite-volume inventory lift to predict the radial wafer flux.

The reduced state contains effective atomic-F, C4F8-derived film precursor,
and positive-ion populations.  Each obeys an exactly integrated linear balance
over a DAQ interval,

    dn/dt = S(machine waveform) - n/tau(pressure),

with production capped by both absorbed RF power and inlet particle supply.
This is a deterministic differentiable equipment-transfer closure, not a
replacement for the repository's product-resolved SF6 collision deck.  Its
effective lifetimes and absorbed-power coupling are explicitly bounded
tool-transfer parameters and must be fitted only on the preregistered
calibration lots.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
import math

import numpy as np

from ..bosch_process_data import (
    BoschProcessTrace, C4F8_FLOW_CHANNEL, SF6_FLOW_CHANNEL,
)
from .axisymmetric_reaction_diffusion import (
    AxisymmetricFiniteVolumeGrid, DeterministicAxisymmetricInventoryLift,
    normalized_annular_skin_source, normalized_exponential_skin_source,
)
from .geometry import CylindricalReactor


_ELECTRON_CHARGE_C = 1.602176634e-19
_BOLTZMANN_J_K = 1.380649e-23
_ATOMIC_MASS_KG = 1.66053906660e-27
_SCCM_PARTICLES_S = 4.477962e17
_PRESSURE_CHANNEL = "Stat3_Etch_MV_Pressure"
_SOURCE_LOAD_CHANNEL = "Stat3_Etch_MV_SourceRFLoadPower"
_SOURCE_REFLECTED_CHANNEL = "Stat3_Etch_MV_SourceRFReflectedPower"
_SOURCE2_LOAD_CHANNEL = "Stat3_Etch_MV_SourceRF2LoadPower"
_PLATEN_VPP_CHANNEL = "Stat3_Etch_MV_PlatenRFPeakToPeak"

_SPECIES = ("F", "C4F8_film_precursor", "positive_ion")


@dataclass(frozen=True)
class BoschSPTSReducedParameters:
    """Physical and equipment-transfer inputs for one Rapier configuration."""

    reactor_radius_m: float = 0.16
    reactor_length_m: float = 0.18
    wafer_radius_m: float = 0.10
    gas_temperature_K: float = 300.0
    electron_temperature_eV: float = 4.0
    effective_positive_ion_mass_amu: float = 127.0
    absorbed_source_power_fraction: float = 0.75
    sf6_dissociation_energy_cost_eV: float = 80.0
    c4f8_fragment_energy_cost_eV: float = 120.0
    ion_pair_energy_cost_eV: float = 250.0
    f_atoms_per_power_limited_dissociation: float = 2.0
    film_units_per_power_limited_fragmentation: float = 1.0
    f_reference_lifetime_s: float = 0.004
    film_precursor_reference_lifetime_s: float = 0.003
    positive_ion_reference_lifetime_s: float = 2.0e-5
    lifetime_reference_pressure_torr: float = 0.04
    neutral_lifetime_pressure_exponent: float = 1.0
    pressure_channel_to_torr: float = 1.0
    sheath_bias_fraction_of_vpp: float = 0.25
    collisional_ion_energy_transmission: float = 0.70
    plasma_potential_per_electron_temperature: float = 4.0
    radial_cell_count: int = 24
    axial_cell_count: int = 24
    source_axial_skin_depth_m: float = 0.045
    source_ring_radius_m: float = 0.085
    source_radial_width_m: float = 0.055
    source_central_fraction: tuple[float, float, float] = (0.0, 0.0, 0.0)
    central_source_radial_scale_m: float = 0.060
    central_source_radial_power: float = 2.0
    diffusion_coefficient_m2_s: tuple[float, float, float] = (0.12, 0.06, 0.20)
    lower_wall_velocity_m_s: tuple[float, float, float] = (120.0, 45.0, 1700.0)
    upper_wall_velocity_m_s: tuple[float, float, float] = (45.0, 20.0, 500.0)
    side_wall_velocity_m_s: tuple[float, float, float] = (55.0, 25.0, 650.0)

    def __post_init__(self):
        positive = (
            "reactor_radius_m", "reactor_length_m", "wafer_radius_m",
            "gas_temperature_K", "electron_temperature_eV",
            "effective_positive_ion_mass_amu", "sf6_dissociation_energy_cost_eV",
            "c4f8_fragment_energy_cost_eV", "ion_pair_energy_cost_eV",
            "f_atoms_per_power_limited_dissociation",
            "film_units_per_power_limited_fragmentation",
            "f_reference_lifetime_s", "film_precursor_reference_lifetime_s",
            "positive_ion_reference_lifetime_s",
            "lifetime_reference_pressure_torr", "pressure_channel_to_torr",
            "plasma_potential_per_electron_temperature",
            "source_axial_skin_depth_m", "source_radial_width_m",
            "central_source_radial_scale_m", "central_source_radial_power",
        )
        if (any(not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0.0
                for name in positive)
                or self.wafer_radius_m > self.reactor_radius_m
                or not 0.0 <= self.source_ring_radius_m <= self.reactor_radius_m
                or not 0.0 < self.absorbed_source_power_fraction <= 1.0
                or not 0.0 < self.sheath_bias_fraction_of_vpp <= 1.0
                or not 0.0 < self.collisional_ion_energy_transmission <= 1.0
                or not math.isfinite(self.neutral_lifetime_pressure_exponent)
                or not 0.0 <= self.neutral_lifetime_pressure_exponent <= 2.0
                or int(self.radial_cell_count) != self.radial_cell_count
                or int(self.axial_cell_count) != self.axial_cell_count
                or min(self.radial_cell_count, self.axial_cell_count) < 4):
            raise ValueError("invalid SPTS reduced-reactor parameters")
        for name in (
                "diffusion_coefficient_m2_s", "lower_wall_velocity_m_s",
                "upper_wall_velocity_m_s", "side_wall_velocity_m_s"):
            value = tuple(float(item) for item in getattr(self, name))
            if len(value) != len(_SPECIES) or any(
                    not math.isfinite(item) or item <= 0.0 for item in value):
                raise ValueError(f"invalid SPTS axisymmetric parameter: {name}")
            object.__setattr__(self, name, value)
        central = tuple(float(value) for value in self.source_central_fraction)
        if (len(central) != len(_SPECIES)
                or any(not math.isfinite(value) or not 0.0 <= value <= 1.0
                       for value in central)):
            raise ValueError("invalid species-resolved central-source fractions")
        object.__setattr__(self, "source_central_fraction", central)

    @property
    def reactor_volume_m3(self):
        return math.pi * self.reactor_radius_m ** 2 * self.reactor_length_m

    def manifest(self):
        return {
            "schema": "petch-spts-bosch-reduced-parameters-v1",
            **{name: getattr(self, name) for name in self.__dataclass_fields__},
            "reactor_volume_m3": self.reactor_volume_m3,
            "equipment_evidence": {
                "process_record": "https://zenodo.org/records/17122442",
                "tool_architecture": (
                    "SPTS Rapier dual ICP sources and dual gas inlets; second source "
                    "is measured zero in this dataset"),
                "tool_architecture_source": (
                    "https://www.epfl.ch/research/facilities/cmi/equipment/etching/spts-rapier/"),
                "reactor_model_precedent": (
                    "AVS 2017 Quantemol/SPTS Rapier HPEM equipment and feature simulation"),
            },
            "calibration_status": (
                "defaults are physically bounded initialization values; tool-transfer "
                "parameters are not certified until calibration receipt is sealed"),
        }


@dataclass(frozen=True)
class BoschSPTSReducedReactorSolution:
    trace_key: str
    interval_midpoint_s: np.ndarray
    interval_duration_s: np.ndarray
    volume_average_density_m3: np.ndarray
    source_rate_m3_s: np.ndarray
    ion_energy_eV: np.ndarray
    final_density_m3: np.ndarray
    integrated_production_m3: np.ndarray
    integrated_loss_m3: np.ndarray
    maximum_inventory_ledger_relative_residual: float
    provenance: dict

    def __post_init__(self):
        time = np.asarray(self.interval_midpoint_s, dtype=float).copy()
        duration = np.asarray(self.interval_duration_s, dtype=float).copy()
        density = np.asarray(self.volume_average_density_m3, dtype=float).copy()
        source = np.asarray(self.source_rate_m3_s, dtype=float).copy()
        energy = np.asarray(self.ion_energy_eV, dtype=float).copy()
        final = np.asarray(self.final_density_m3, dtype=float).copy()
        produced = np.asarray(self.integrated_production_m3, dtype=float).copy()
        lost = np.asarray(self.integrated_loss_m3, dtype=float).copy()
        n = time.size
        if (not self.trace_key or duration.shape != (n,)
                or density.shape != (n, len(_SPECIES))
                or source.shape != density.shape or energy.shape != (n,)
                or final.shape != (len(_SPECIES),)
                or produced.shape != final.shape or lost.shape != final.shape
                or any(np.any(~np.isfinite(value)) for value in (
                    time, duration, density, source, energy, final, produced, lost))
                or np.any(duration <= 0.0) or np.any(density < 0.0)
                or np.any(source < 0.0) or np.any(energy < 0.0)
                or np.any(final < 0.0) or np.any(produced < 0.0)
                or np.any(lost < 0.0)
                or not math.isfinite(self.maximum_inventory_ledger_relative_residual)
                or not 0.0 <= self.maximum_inventory_ledger_relative_residual < 1.0e-10):
            raise ValueError("invalid SPTS reduced-reactor solution")
        for value in (time, duration, density, source, energy, final, produced, lost):
            value.setflags(write=False)
        object.__setattr__(self, "interval_midpoint_s", time)
        object.__setattr__(self, "interval_duration_s", duration)
        object.__setattr__(self, "volume_average_density_m3", density)
        object.__setattr__(self, "source_rate_m3_s", source)
        object.__setattr__(self, "ion_energy_eV", energy)
        object.__setattr__(self, "final_density_m3", final)
        object.__setattr__(self, "integrated_production_m3", produced)
        object.__setattr__(self, "integrated_loss_m3", lost)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


@dataclass(frozen=True)
class BoschSPTSWaferBoundaryTrace:
    reactor: BoschSPTSReducedReactorSolution
    species_names: tuple[str, ...]
    radial_centers_m: np.ndarray
    radial_flux_m2_s: np.ndarray
    wafer_area_average_flux_m2_s: np.ndarray
    maximum_axisymmetric_species_ledger_relative_residual: float
    inventory_lift_condition_number: float
    source_jvp_supported: bool

    def __post_init__(self):
        names = tuple(self.species_names)
        radial = np.asarray(self.radial_centers_m, dtype=float).copy()
        flux = np.asarray(self.radial_flux_m2_s, dtype=float).copy()
        average = np.asarray(self.wafer_area_average_flux_m2_s, dtype=float).copy()
        n = self.reactor.interval_midpoint_s.size
        if (names != _SPECIES or radial.ndim != 1 or radial.size < 4
                or flux.shape != (n, len(names), radial.size)
                or average.shape != (n, len(names))
                or np.any(~np.isfinite(radial)) or np.any(radial <= 0.0)
                or np.any(~np.isfinite(flux)) or np.any(flux < 0.0)
                or np.any(~np.isfinite(average)) or np.any(average < 0.0)
                or not math.isfinite(
                    self.maximum_axisymmetric_species_ledger_relative_residual)
                or self.maximum_axisymmetric_species_ledger_relative_residual >= 1.0e-8
                or not math.isfinite(self.inventory_lift_condition_number)
                or self.inventory_lift_condition_number < 1.0
                or not bool(self.source_jvp_supported)):
            raise ValueError("invalid SPTS wafer-boundary trace")
        for value in (radial, flux, average):
            value.setflags(write=False)
        object.__setattr__(self, "species_names", names)
        object.__setattr__(self, "radial_centers_m", radial)
        object.__setattr__(self, "radial_flux_m2_s", flux)
        object.__setattr__(self, "wafer_area_average_flux_m2_s", average)


def _relaxation_interval(initial, source_rate, lifetime_s, duration_s):
    """Return exact interval-mean/end density and integrated first-order loss."""
    x = duration_s / lifetime_s
    decay = math.exp(-x)
    steady = source_rate * lifetime_s
    end = steady + (initial - steady) * decay
    phi = -math.expm1(-x) / x
    mean = steady + (initial - steady) * phi
    produced = source_rate * duration_s
    lost = produced - (end - initial)
    scale = max(abs(produced), abs(end - initial), 1.0)
    if end < -1.0e-12 * max(abs(steady), abs(initial), 1.0) or lost < -1.0e-12 * scale:
        raise RuntimeError("exact reduced-reactor relaxation lost positivity")
    return max(mean, 0.0), max(end, 0.0), max(lost, 0.0)


def solve_bosch_spts_reduced_reactor(
        trace: BoschProcessTrace,
        parameters: BoschSPTSReducedParameters,
        *, initial_density_m3=(0.0, 0.0, 0.0)):
    if not isinstance(trace, BoschProcessTrace):
        raise TypeError("trace must be BoschProcessTrace")
    if not isinstance(parameters, BoschSPTSReducedParameters):
        raise TypeError("parameters must be BoschSPTSReducedParameters")
    initial = np.asarray(initial_density_m3, dtype=float)
    if (initial.shape != (len(_SPECIES),) or np.any(~np.isfinite(initial))
            or np.any(initial < 0.0)):
        raise ValueError("invalid reduced-reactor initial density")

    elapsed = trace.elapsed_s
    dt = np.diff(elapsed)
    interval_count = dt.size
    density = np.zeros((interval_count, len(_SPECIES)))
    source = np.zeros_like(density)
    energy = np.zeros(interval_count)
    produced = np.zeros(len(_SPECIES))
    lost = np.zeros(len(_SPECIES))
    state = initial.copy()
    volume = parameters.reactor_volume_m3
    channels = trace.channels

    for index, duration in enumerate(dt):
        sf6_flow = max(float(channels[SF6_FLOW_CHANNEL][index]), 0.0)
        c4f8_flow = max(float(channels[C4F8_FLOW_CHANNEL][index]), 0.0)
        total_flow = sf6_flow + c4f8_flow
        source2 = max(float(channels[_SOURCE2_LOAD_CHANNEL][index]), 0.0)
        if source2 > 1.0e-9:
            raise ValueError(
                "SPTS reduced v1 is frozen for the measured source-2-off experiment")
        net_source_power = max(
            float(channels[_SOURCE_LOAD_CHANNEL][index])
            - float(channels[_SOURCE_REFLECTED_CHANNEL][index]), 0.0)
        absorbed_power = parameters.absorbed_source_power_fraction * net_source_power

        sf6_supply_s = sf6_flow * _SCCM_PARTICLES_S
        c4f8_supply_s = c4f8_flow * _SCCM_PARTICLES_S
        total_supply_s = total_flow * _SCCM_PARTICLES_S
        f_power_s = (
            absorbed_power
            / (parameters.sf6_dissociation_energy_cost_eV * _ELECTRON_CHARGE_C)
            * parameters.f_atoms_per_power_limited_dissociation
            if sf6_flow > 0.0 else 0.0)
        film_power_s = (
            absorbed_power
            / (parameters.c4f8_fragment_energy_cost_eV * _ELECTRON_CHARGE_C)
            * parameters.film_units_per_power_limited_fragmentation
            if c4f8_flow > 0.0 else 0.0)
        ion_power_s = (
            absorbed_power
            / (parameters.ion_pair_energy_cost_eV * _ELECTRON_CHARGE_C)
            if total_flow > 0.0 else 0.0)
        source[index] = (
            min(f_power_s, 6.0 * sf6_supply_s) / volume,
            min(film_power_s, c4f8_supply_s) / volume,
            min(ion_power_s, total_supply_s) / volume,
        )

        pressure_torr = max(
            float(channels[_PRESSURE_CHANNEL][index])
            * parameters.pressure_channel_to_torr, 1.0e-6)
        neutral_scale = np.clip(
            (pressure_torr / parameters.lifetime_reference_pressure_torr)
            ** parameters.neutral_lifetime_pressure_exponent,
            0.25, 4.0)
        lifetimes = (
            parameters.f_reference_lifetime_s * neutral_scale,
            parameters.film_precursor_reference_lifetime_s * neutral_scale,
            parameters.positive_ion_reference_lifetime_s,
        )
        for species, lifetime in enumerate(lifetimes):
            mean, end, interval_loss = _relaxation_interval(
                state[species], source[index, species], lifetime, duration)
            density[index, species] = mean
            state[species] = end
            produced[species] += source[index, species] * duration
            lost[species] += interval_loss

        vpp = abs(float(channels[_PLATEN_VPP_CHANNEL][index]))
        energy[index] = (
            parameters.plasma_potential_per_electron_temperature
            * parameters.electron_temperature_eV
            + parameters.collisional_ion_energy_transmission
            * parameters.sheath_bias_fraction_of_vpp * vpp)

    ledger = np.abs(initial + produced - lost - state)
    ledger_scale = np.maximum.reduce((
        np.abs(initial) + np.abs(produced), np.abs(lost) + np.abs(state),
        np.ones(len(_SPECIES))))
    maximum_ledger = float(np.max(ledger / ledger_scale))
    return BoschSPTSReducedReactorSolution(
        trace_key=trace.experiment_key,
        interval_midpoint_s=0.5 * (elapsed[:-1] + elapsed[1:]),
        interval_duration_s=dt,
        volume_average_density_m3=density,
        source_rate_m3_s=source,
        ion_energy_eV=energy,
        final_density_m3=state,
        integrated_production_m3=produced,
        integrated_loss_m3=lost,
        maximum_inventory_ledger_relative_residual=maximum_ledger,
        provenance={
            "model": "spts-bosch-measured-waveform-reduced-reactor-v1",
            "species_order": list(_SPECIES),
            "equation": "exact piecewise-constant dn/dt=S-n/tau",
            "production_limit": "minimum of absorbed-power and inlet-particle supply",
            "source2_measured_off_required": True,
            "parameters": parameters.manifest(),
        })


class DeterministicBoschSPTSReactorToWafer:
    """Reduced 0-D state plus deterministic axisymmetric radial transfer."""

    def __init__(self, parameters: BoschSPTSReducedParameters):
        if not isinstance(parameters, BoschSPTSReducedParameters):
            raise TypeError("parameters must be BoschSPTSReducedParameters")
        self.parameters = parameters
        geometry = CylindricalReactor(
            radius_m=parameters.reactor_radius_m,
            length_m=parameters.reactor_length_m)
        grid = AxisymmetricFiniteVolumeGrid.uniform(
            geometry,
            radial_cell_count=parameters.radial_cell_count,
            axial_cell_count=parameters.axial_cell_count)
        annular_source = normalized_annular_skin_source(
            grid,
            axial_skin_depth_m=parameters.source_axial_skin_depth_m,
            ring_radius_m=parameters.source_ring_radius_m,
            radial_width_m=parameters.source_radial_width_m)
        central_source = normalized_exponential_skin_source(
            grid,
            axial_skin_depth_m=parameters.source_axial_skin_depth_m,
            radial_scale_m=parameters.central_source_radial_scale_m,
            radial_power=parameters.central_source_radial_power)
        source_shapes = np.stack(tuple(
            central_fraction * central_source
            + (1.0 - central_fraction) * annular_source
            for central_fraction in parameters.source_central_fraction))
        wall = np.column_stack((
            parameters.lower_wall_velocity_m_s,
            parameters.upper_wall_velocity_m_s,
            parameters.side_wall_velocity_m_s))
        self._lift = DeterministicAxisymmetricInventoryLift(
            grid=grid, species_names=_SPECIES,
            diffusion_coefficient_m2_s=np.asarray(
                parameters.diffusion_coefficient_m2_s),
            volume_reaction_matrix_s_inv=np.zeros((len(_SPECIES), len(_SPECIES))),
            wall_velocity_m_s=wall,
            source_shape=source_shapes,
            source=(
                "SPTS Rapier species-resolved central/annular source-1 skin "
                "moments; source-2 is measured off"))
        unit = self._lift.solve(np.ones(len(_SPECIES)))
        self._unit_lower_flux_per_density_m_s = (
            unit.solution.lower_endcap_flux_m2_s.copy())
        self._maximum_axisymmetric_ledger = (
            unit.solution.maximum_species_ledger_relative_residual)
        self._grid = grid

    def solve(
            self, trace: BoschProcessTrace, *, initial_density_m3=(0.0, 0.0, 0.0)):
        reactor = solve_bosch_spts_reduced_reactor(
            trace, self.parameters, initial_density_m3=initial_density_m3)
        radial_flux = (
            reactor.volume_average_density_m3[:, :, None]
            * self._unit_lower_flux_per_density_m_s[None, :, :])
        grid = self._grid
        clipped_outer = np.minimum(
            grid.radial_edges_m[1:], self.parameters.wafer_radius_m)
        clipped_inner = np.minimum(
            grid.radial_edges_m[:-1], self.parameters.wafer_radius_m)
        area = np.pi * np.maximum(clipped_outer ** 2 - clipped_inner ** 2, 0.0)
        average = np.einsum("tsr,r->ts", radial_flux, area) / (
            np.pi * self.parameters.wafer_radius_m ** 2)
        return BoschSPTSWaferBoundaryTrace(
            reactor=reactor,
            species_names=_SPECIES,
            radial_centers_m=grid.radial_centers_m,
            radial_flux_m2_s=radial_flux,
            wafer_area_average_flux_m2_s=average,
            maximum_axisymmetric_species_ledger_relative_residual=(
                self._maximum_axisymmetric_ledger),
            inventory_lift_condition_number=(
                self._lift.source_response_condition_number),
            source_jvp_supported=True)

    def density_to_radial_flux_jvp(self, density_tangent_m3):
        tangent = np.asarray(density_tangent_m3, dtype=float)
        if (tangent.ndim != 2 or tangent.shape[1] != len(_SPECIES)
                or np.any(~np.isfinite(tangent))):
            raise ValueError("density tangent must have shape (time, 3)")
        return (
            tangent[:, :, None]
            * self._unit_lower_flux_per_density_m_s[None, :, :])
