"""Deterministic ion Boltzmann transport in a moving RF sheath.

The earlier implicit sheath solver closes all ion collision orders, but does
so independently at each entry phase in an effective static Child field.  In
this module RF phase is part of the kinetic state.  Between collisions ions
follow the time-dependent Turner--Chabert Poisson field; after an elastic or
charge-exchange event, the new ion is deposited at its collision position,
energy, direction, *and RF phase*.  The complete collision series is again
summed by the sparse absorbing solve

``(I - Q.T) visits = source``.

No particle Monte Carlo is used.  First-event probabilities are integrated on
a fixed time mesh; collision impact parameter and azimuth use fixed Gaussian
and periodic quadrature.  The density derivative is the exact derivative of
the assembled discrete operator.

This closes the moving-field ion problem.  It does not close repeated
fast-neutral collisions, generator/matching-network inversion, molecular-ion
cross sections, or an etch surface law, so feature-depth support remains
false.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu

from .collisional_sheath import (
    BOLTZMANN_J_K,
    E_CHARGE_C,
    ArgonBornMayerPhelpsCollisionModel,
    CollisionalIonEnergyAngleDistribution,
    CollisionalSheathDensityTangent,
    DeterministicCollisionalSheathSolution,
    _canonical_axisymmetric_velocity,
    _equal_mass_collision_velocities,
)
from .current_driven_rf_sheath import TurnerChabertCurrentDrivenSheath


@dataclass(frozen=True)
class _MovingCollisionBranch:
    ion_velocity: np.ndarray
    neutral_velocity: np.ndarray
    position_fraction: float
    phase_rad: float
    coefficient: float
    density_derivative: float
    charge_exchange: bool
    below_angular_support: bool


@dataclass(frozen=True)
class _MovingStateKernel:
    endpoint_kind: str
    endpoint_velocity: np.ndarray
    endpoint_phase_rad: float
    survival: float
    survival_density_derivative: float
    collision_probability: float
    collision_probability_density_derivative: float
    branches: tuple[_MovingCollisionBranch, ...]
    unit_optical_depth: float
    maximum_collision_energy_residual: float
    maximum_verlet_energy_defect_eV: float


@dataclass(frozen=True)
class MovingSheathConvergenceReceipt:
    coarse: DeterministicCollisionalSheathSolution
    fine: DeterministicCollisionalSheathSolution
    mean_energy_relative_change: float
    rms_angle_relative_change: float
    collision_count_relative_change: float
    arrival_probability_relative_change: float
    probability_ledger_residual: float
    passed: bool
    limits: Mapping[str, float]

    def __post_init__(self):
        values = np.asarray([
            self.mean_energy_relative_change,
            self.rms_angle_relative_change,
            self.collision_count_relative_change,
            self.arrival_probability_relative_change,
            self.probability_ledger_residual,
        ])
        if (
            np.any(~np.isfinite(values))
            or np.any(values < 0.0)
            or self.passed != all((
                self.mean_energy_relative_change
                <= self.limits["mean_energy_relative_change"],
                self.rms_angle_relative_change
                <= self.limits["rms_angle_relative_change"],
                self.collision_count_relative_change
                <= self.limits["collision_count_relative_change"],
                self.arrival_probability_relative_change
                <= self.limits["arrival_probability_relative_change"],
                self.probability_ledger_residual
                <= self.limits["probability_ledger_residual"],
            ))
        ):
            raise ValueError("invalid moving-sheath convergence receipt")
        object.__setattr__(self, "limits", MappingProxyType(dict(self.limits)))


@dataclass(frozen=True)
class DeterministicMovingCollisionalRFSheath:
    """Phase-space discrete ordinates for a time-dependent collisional sheath."""

    sheath: TurnerChabertCurrentDrivenSheath
    collision_model: ArgonBornMayerPhelpsCollisionModel
    gas_number_density_m3: float
    neutral_gas_temperature_K: float
    source_ion_flux_m2_s: float
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
    provenance: Mapping[str, object] = field(default_factory=lambda: {
        "feature_depth_used": False,
        "state_coordinates": (
            "RF phase x position x signed direction x total energy x "
            "transverse-energy fraction"
        ),
    })

    def __post_init__(self):
        values = np.asarray([
            self.gas_number_density_m3,
            self.neutral_gas_temperature_K,
            self.source_ion_flux_m2_s,
            self.maximum_transit_periods,
        ])
        orders = (
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
            not isinstance(self.sheath, TurnerChabertCurrentDrivenSheath)
            or not isinstance(
                self.collision_model, ArgonBornMayerPhelpsCollisionModel)
            or np.any(~np.isfinite(values))
            or self.gas_number_density_m3 < 0.0
            or self.neutral_gas_temperature_K <= 0.0
            or self.source_ion_flux_m2_s <= 0.0
            or self.maximum_transit_periods <= 0.0
            or any(int(value) < 1 for value in orders)
            or int(self.phase_node_count) < 4
            or int(self.position_node_count) < 3
            or int(self.total_energy_node_count) < 3
            or int(self.transverse_fraction_node_count) < 3
            or int(self.output_azimuth_order) < 2
            or int(self.collision_azimuth_order) % 2 != 0
            or int(self.steps_per_period) < 16
            or int(self.steps_per_transit) < 16
        ):
            raise ValueError("invalid moving collisional RF sheath")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def tangential_temperature_eV(self) -> float:
        return (
            BOLTZMANN_J_K * self.neutral_gas_temperature_K / E_CHARGE_C)

    @property
    def phase_axis_rad(self) -> np.ndarray:
        return (
            2.0
            * np.pi
            * (np.arange(int(self.phase_node_count), dtype=float) + 0.5)
            / int(self.phase_node_count)
        )

    def _initial_velocity_quadrature(self):
        node, weight = np.polynomial.laguerre.laggauss(
            int(self.initial_thermal_radial_order))
        bohm_energy = 0.5 * float(self.sheath.Te_eV)
        for radial, coefficient in zip(node, weight):
            yield (
                np.array([
                    math.sqrt(self.tangential_temperature_eV * float(radial)),
                    0.0,
                    math.sqrt(bohm_energy),
                ]),
                float(coefficient),
            )

    @staticmethod
    def _linear_deposit(axis: np.ndarray, value: float):
        if value <= axis[0]:
            return ((0, 1.0),)
        if value >= axis[-1]:
            return ((len(axis) - 1, 1.0),)
        upper = int(np.searchsorted(axis, value, side="right"))
        lower = upper - 1
        fraction = (value - axis[lower]) / (axis[upper] - axis[lower])
        return ((lower, 1.0 - fraction), (upper, fraction))

    def _periodic_phase_deposit(self, phase_rad: float):
        count = int(self.phase_node_count)
        coordinate = (
            np.mod(float(phase_rad), 2.0 * np.pi)
            * count / (2.0 * np.pi)
            - 0.5
        )
        lower_unwrapped = math.floor(coordinate)
        fraction = coordinate - lower_unwrapped
        lower = int(lower_unwrapped % count)
        upper = int((lower_unwrapped + 1) % count)
        if fraction <= 2.0e-15:
            return ((lower, 1.0),)
        return ((lower, 1.0 - fraction), (upper, fraction))

    def _axes(self):
        initial = tuple(self._initial_velocity_quadrature())
        initial_max = max(float(np.dot(v, v)) for v, _ in initial)
        energy_cap = initial_max + self.sheath.maximum_voltage_v
        position = np.linspace(0.0, 1.0, int(self.position_node_count))
        # Quadratic energy spacing resolves slow CX ions while linear
        # deposition still preserves the energy moment exactly.
        unit = np.linspace(0.0, 1.0, int(self.total_energy_node_count))
        energy = energy_cap * unit ** 2
        mu = np.linspace(0.0, 1.0, int(self.transverse_fraction_node_count))
        transverse = 1.0 - mu[::-1] ** 2
        return initial, position, energy, transverse

    @staticmethod
    def _velocity_from_state(
        direction_index: int,
        energy_eV: float,
        transverse_fraction: float,
    ) -> np.ndarray:
        transverse = math.sqrt(max(energy_eV * transverse_fraction, 0.0))
        normal = math.sqrt(max(energy_eV * (1.0 - transverse_fraction), 0.0))
        if direction_index == 0:
            normal = -normal
        return np.array([transverse, 0.0, normal])

    def _state_deposition(
        self,
        velocity: np.ndarray,
        position_fraction: float,
        phase_rad: float,
        position_axis: np.ndarray,
        energy_axis: np.ndarray,
        transverse_axis: np.ndarray,
    ):
        total = float(np.dot(velocity, velocity))
        fraction = (
            0.0 if total <= np.finfo(float).tiny
            else float(np.dot(velocity[:2], velocity[:2])) / total
        )
        direction = 0 if float(velocity[2]) < 0.0 else 1
        shape = (
            len(self.phase_axis_rad), 2, len(position_axis),
            len(energy_axis), len(transverse_axis),
        )
        for ip, wp in self._periodic_phase_deposit(phase_rad):
            for ix, wx in self._linear_deposit(
                position_axis, float(position_fraction)
            ):
                for ie, we in self._linear_deposit(energy_axis, total):
                    for it, wt in self._linear_deposit(
                        transverse_axis, fraction
                    ):
                        flat = int(np.ravel_multi_index(
                            (ip, direction, ix, ie, it), shape))
                        yield flat, float(wp * wx * we * wt)

    def _output_deposition(
        self,
        velocity: np.ndarray,
        phase_rad: float,
        energy_axis: np.ndarray,
        transverse_axis: np.ndarray,
    ):
        total = float(np.dot(velocity, velocity))
        fraction = (
            0.0 if total <= np.finfo(float).tiny
            else float(np.dot(velocity[:2], velocity[:2])) / total
        )
        shape = (len(self.phase_axis_rad), len(energy_axis), len(transverse_axis))
        for ip, wp in self._periodic_phase_deposit(phase_rad):
            for ie, we in self._linear_deposit(energy_axis, total):
                for it, wt in self._linear_deposit(transverse_axis, fraction):
                    flat = int(np.ravel_multi_index((ip, ie, it), shape))
                    yield flat, float(wp * we * wt)

    def _time_step_s(self) -> float:
        bohm = math.sqrt(0.5 * self.sheath.Te_eV)
        maximum_speed_sqrt_eV = math.sqrt(
            0.5 * self.sheath.Te_eV + self.sheath.maximum_voltage_v)
        # sqrt(eV) velocity coordinates convert to physical velocity with the
        # same factor in numerator/denominator; use a physical transit from
        # the common sheath method for clarity.
        mass_kg = self.sheath.ion_mass_amu * 1.66053906660e-27
        physical_bohm = math.sqrt(2.0 * E_CHARGE_C / mass_kg) * bohm
        physical_maximum = (
            math.sqrt(2.0 * E_CHARGE_C / mass_kg)
            * maximum_speed_sqrt_eV
        )
        transit = (
            2.0 * self.sheath.maximum_width_m
            / max(physical_bohm + physical_maximum, 1.0e-30)
        )
        return min(
            self.sheath.period_s / int(self.steps_per_period),
            transit / int(self.steps_per_transit),
        )

    def _trajectory_kernel(
        self,
        velocity_sqrt_eV: np.ndarray,
        position_fraction: float,
        phase_rad: float,
        impact_cache: dict[float, tuple[np.ndarray, np.ndarray]],
    ) -> _MovingStateKernel:
        velocity = np.asarray(velocity_sqrt_eV, dtype=float).copy()
        position = float(position_fraction) * self.sheath.maximum_width_m
        time = (
            np.mod(float(phase_rad), 2.0 * np.pi)
            * self.sheath.period_s / (2.0 * np.pi)
        )
        # A zero-normal-speed ordinate on the plasma-edge boundary is not an
        # incoming Bohm characteristic.  Tangential energy cannot carry it
        # into the sheath because E(x=0)=0, so it belongs to the plasma-edge
        # outflow sink just like a negative-normal-speed ordinate.
        if position <= 2.0e-15 and velocity[2] <= 1.0e-14:
            return _MovingStateKernel(
                endpoint_kind="escaped",
                endpoint_velocity=velocity,
                endpoint_phase_rad=float(phase_rad),
                survival=1.0,
                survival_density_derivative=0.0,
                collision_probability=0.0,
                collision_probability_density_derivative=0.0,
                branches=(),
                unit_optical_depth=0.0,
                maximum_collision_energy_residual=0.0,
                maximum_verlet_energy_defect_eV=0.0,
            )
        mass_kg = self.sheath.ion_mass_amu * 1.66053906660e-27
        velocity_scale = math.sqrt(2.0 * E_CHARGE_C / mass_kg)
        dt = self._time_step_s()
        maximum_steps = int(math.ceil(
            self.maximum_transit_periods * self.sheath.period_s / dt))
        density = float(self.gas_number_density_m3)
        accumulated_unit_hazard = 0.0
        trajectory_segments = []
        maximum_verlet_defect = 0.0
        endpoint_kind = ""

        for _ in range(maximum_steps):
            physical_velocity = velocity * velocity_scale
            a0 = (
                E_CHARGE_C
                * self.sheath.electric_field_V_m(position, time)
                / mass_kg
            )
            normal_half = physical_velocity[2] + 0.5 * a0 * dt
            next_position = position + normal_half * dt
            next_time = time + dt
            a1 = (
                E_CHARGE_C
                * self.sheath.electric_field_V_m(next_position, next_time)
                / mass_kg
            )
            next_physical = physical_velocity.copy()
            next_physical[2] = normal_half + 0.5 * a1 * dt
            next_velocity = next_physical / velocity_scale

            endpoint = None
            if next_position >= self.sheath.maximum_width_m:
                endpoint = self.sheath.maximum_width_m
                endpoint_kind = "arrived"
            elif next_position <= 0.0 and next_physical[2] < 0.0:
                endpoint = 0.0
                endpoint_kind = "escaped"
            fraction = 1.0
            if endpoint is not None:
                denominator = next_position - position
                fraction = float(np.clip(
                    (endpoint - position) / denominator
                    if abs(denominator) > 1.0e-30 else 0.0,
                    0.0,
                    1.0,
                ))
            segment_dt = dt * fraction
            midpoint_time = time + 0.5 * segment_dt
            midpoint_position = position + 0.5 * normal_half * segment_dt
            midpoint_physical = 0.5 * (
                physical_velocity + (
                    physical_velocity
                    + fraction * (next_physical - physical_velocity)
                )
            )
            midpoint_velocity = midpoint_physical / velocity_scale
            event_energy = float(np.dot(midpoint_velocity, midpoint_velocity))
            speed = float(np.linalg.norm(midpoint_physical))
            unit_segment_hazard = (
                0.0 if speed <= 0.0 or event_energy <= 0.0
                else self.collision_model.total_cross_section_m2(event_energy)
                * speed * segment_dt
            )
            if unit_segment_hazard > 0.0:
                trajectory_segments.append((
                    accumulated_unit_hazard,
                    accumulated_unit_hazard + unit_segment_hazard,
                    position,
                    position + fraction * (next_position - position),
                    time,
                    time + segment_dt,
                    velocity.copy(),
                    (
                        physical_velocity
                        + fraction * (next_physical - physical_velocity)
                    ) / velocity_scale,
                ))
            accumulated_unit_hazard += unit_segment_hazard

            # In the zero-density case this is an independent check on the
            # time integrator: kinetic gain plus instantaneous potential is
            # not conserved in an RF field, but a static snapshot is.  Store
            # only the local Verlet work defect over one step.
            work_eV = 0.5 * (
                self.sheath.electric_field_V_m(position, time)
                + self.sheath.electric_field_V_m(next_position, next_time)
            ) * (next_position - position)
            kinetic_change = (
                float(np.dot(next_velocity, next_velocity))
                - float(np.dot(velocity, velocity))
            )
            maximum_verlet_defect = max(
                maximum_verlet_defect, abs(kinetic_change - work_eV))
            if endpoint is not None:
                velocity = (
                    physical_velocity
                    + fraction * (next_physical - physical_velocity)
                ) / velocity_scale
                position = endpoint
                time += segment_dt
                break
            position = next_position
            time = next_time
            velocity = next_velocity
        else:
            # Exactly stationary zero-energy ordinate at the mean sheath edge
            # is not a kinetic ion state; the presheath would restore Bohm
            # entry.  Treat it as an escape sink rather than fabricate a long
            # residence collision source.
            if position <= 2.0e-15 and np.linalg.norm(velocity) <= 1.0e-14:
                endpoint_kind = "escaped"
            else:
                raise RuntimeError(
                    "moving collisional trajectory exceeded transit horizon")

        survival = math.exp(-density * accumulated_unit_hazard)
        survival_derivative = -accumulated_unit_hazard * survival
        collision_probability = 1.0 - survival
        collision_derivative = -survival_derivative
        event_records = []
        if accumulated_unit_hazard > 0.0:
            # The first-event density in unit optical depth tau is
            # n_g exp(-n_g tau).  A fixed Gaussian rule over [0, tau_total],
            # normalized to the exact collision probability, avoids one
            # scattering expansion per Verlet step.  Its n_g -> 0 limit is
            # uniform in optical depth and retains the exact boundary JVP.
            node, weight = np.polynomial.legendre.leggauss(
                int(self.collision_event_quadrature_order))
            optical_node = 0.5 * accumulated_unit_hazard * (node + 1.0)
            optical_weight = 0.5 * accumulated_unit_hazard * weight
            raw = optical_weight * np.exp(-density * optical_node)
            raw_derivative = -optical_node * raw
            raw_total = float(np.sum(raw))
            raw_total_derivative = float(np.sum(raw_derivative))
            if raw_total <= 0.0 or not math.isfinite(raw_total):
                raise RuntimeError("first-event optical quadrature lost mass")
            segment_ends = np.asarray([
                segment[1] for segment in trajectory_segments])
            for tau, raw_value, raw_dn in zip(
                optical_node, raw, raw_derivative
            ):
                segment_index = min(
                    int(np.searchsorted(segment_ends, tau, side="right")),
                    len(trajectory_segments) - 1,
                )
                segment = trajectory_segments[segment_index]
                tau0, tau1, x0, x1, t0, t1, v0, v1 = segment
                local = float(np.clip(
                    (float(tau) - tau0) / max(tau1 - tau0, 1.0e-300),
                    0.0,
                    1.0,
                ))
                conditional = float(raw_value / raw_total)
                conditional_derivative = float(
                    (raw_dn * raw_total
                     - raw_value * raw_total_derivative)
                    / raw_total ** 2
                )
                event_probability = collision_probability * conditional
                event_derivative = (
                    collision_derivative * conditional
                    + collision_probability * conditional_derivative
                )
                event_position = x0 + local * (x1 - x0)
                event_time = t0 + local * (t1 - t0)
                event_velocity = v0 + local * (v1 - v0)
                event_records.append((
                    event_position / self.sheath.maximum_width_m,
                    np.mod(
                        2.0 * np.pi * event_time / self.sheath.period_s,
                        2.0 * np.pi,
                    ),
                    event_velocity,
                    event_probability,
                    event_derivative,
                ))
        azimuth = (
            2.0
            * np.pi
            * (np.arange(int(self.collision_azimuth_order)) + 0.5)
            / int(self.collision_azimuth_order)
        )
        branches = []
        maximum_collision_residual = 0.0
        for x, phase, event_velocity, coefficient, derivative in event_records:
            event_energy = float(np.dot(event_velocity, event_velocity))
            cache_key = round(event_energy, 8)
            if cache_key not in impact_cache:
                impact_cache[cache_key] = self.collision_model.impact_quadrature(
                    event_energy, int(self.impact_quadrature_order))
            angles, weights = impact_cache[cache_key]
            below = (
                event_energy
                < self.collision_model.born_mayer_minimum_lab_energy_eV)
            elastic_probability, charge_exchange_probability = (
                self.collision_model.channel_probabilities(event_energy))
            for angle, impact_weight in zip(angles, weights):
                for phi in azimuth:
                    projectile, target = _equal_mass_collision_velocities(
                        event_velocity, float(angle), float(phi))
                    cx_projectile, cx_target = (
                        (projectile, target)
                        if not below
                        else _equal_mass_collision_velocities(
                            event_velocity, 0.0, float(phi))
                    )
                    residual = abs(
                        float(np.dot(projectile, projectile))
                        + float(np.dot(target, target))
                        - event_energy
                    ) / max(event_energy, 1.0)
                    maximum_collision_residual = max(
                        maximum_collision_residual, residual)
                    common = float(impact_weight) / len(azimuth)
                    for ion, neutral, channel, charge_exchange in (
                        (
                            projectile, target,
                            elastic_probability, False,
                        ),
                        (
                            cx_target, cx_projectile,
                            charge_exchange_probability,
                            True,
                        ),
                    ):
                        fraction = channel * common
                        branches.append(_MovingCollisionBranch(
                            ion_velocity=_canonical_axisymmetric_velocity(ion),
                            neutral_velocity=_canonical_axisymmetric_velocity(
                                neutral),
                            position_fraction=float(x),
                            phase_rad=float(phase),
                            coefficient=coefficient * fraction,
                            density_derivative=derivative * fraction,
                            charge_exchange=charge_exchange,
                            below_angular_support=below,
                        ))
        if (
            abs(sum(branch.coefficient for branch in branches)
                - collision_probability) > 3.0e-11
            or abs(sum(branch.density_derivative for branch in branches)
                   - collision_derivative) > 3.0e-25
        ):
            raise RuntimeError("moving collision quadrature lost probability")
        return _MovingStateKernel(
            endpoint_kind=endpoint_kind,
            endpoint_velocity=_canonical_axisymmetric_velocity(velocity),
            endpoint_phase_rad=float(np.mod(
                2.0 * np.pi * time / self.sheath.period_s,
                2.0 * np.pi,
            )),
            survival=survival,
            survival_density_derivative=survival_derivative,
            collision_probability=collision_probability,
            collision_probability_density_derivative=collision_derivative,
            branches=tuple(branches),
            unit_optical_depth=accumulated_unit_hazard,
            maximum_collision_energy_residual=maximum_collision_residual,
            maximum_verlet_energy_defect_eV=maximum_verlet_defect,
        )

    def _neutral_lower_bound(self, branch: _MovingCollisionBranch):
        velocity = branch.neutral_velocity
        speed = float(np.linalg.norm(velocity))
        normal = float(velocity[2])
        if speed <= 0.0 or normal <= 0.0:
            return "escaped", branch.coefficient, branch.density_derivative
        unit_hazard = (
            self.sheath.maximum_width_m
            * (1.0 - branch.position_fraction)
            * speed / normal
            * self.collision_model.total_cross_section_m2(speed * speed)
        )
        survival = math.exp(-self.gas_number_density_m3 * unit_hazard)
        survival_derivative = -unit_hazard * survival
        arrival = branch.coefficient * survival
        arrival_derivative = (
            branch.density_derivative * survival
            + branch.coefficient * survival_derivative)
        unresolved = branch.coefficient * (1.0 - survival)
        unresolved_derivative = (
            branch.density_derivative * (1.0 - survival)
            - branch.coefficient * survival_derivative)
        flight_time = (
            self.sheath.maximum_width_m * (1.0 - branch.position_fraction)
            / (normal * math.sqrt(
                2.0 * E_CHARGE_C
                / (self.sheath.ion_mass_amu * 1.66053906660e-27)))
        )
        arrival_phase = np.mod(
            branch.phase_rad
            + 2.0 * np.pi * flight_time / self.sheath.period_s,
            2.0 * np.pi,
        )
        return (
            "arrived", arrival, arrival_derivative,
            unresolved, unresolved_derivative, float(arrival_phase),
        )

    def _assemble_and_solve(self, density_direction: float):
        initial, position_axis, energy_axis, transverse_axis = self._axes()
        shape = (
            len(self.phase_axis_rad), 2, len(position_axis),
            len(energy_axis), len(transverse_axis),
        )
        output_shape = (
            len(self.phase_axis_rad), len(energy_axis), len(transverse_axis))
        state_count = int(np.prod(shape))
        output_count = int(np.prod(output_shape))
        impact_cache: dict[float, tuple[np.ndarray, np.ndarray]] = {}

        q_row, q_col, q_data, dq_data = [], [], [], []
        b_row, b_col, b_data, db_data = [], [], [], []
        n_row, n_col, n_data, dn_data = [], [], [], []
        arrival = np.zeros(state_count)
        arrival_dn = np.zeros(state_count)
        escaped = np.zeros(state_count)
        escaped_dn = np.zeros(state_count)
        collision = np.zeros(state_count)
        collision_dn = np.zeros(state_count)
        charge_exchange = np.zeros(state_count)
        charge_exchange_dn = np.zeros(state_count)
        birth_energy = np.zeros(state_count)
        birth_energy_dn = np.zeros(state_count)
        neutral_arrival = np.zeros(state_count)
        neutral_arrival_dn = np.zeros(state_count)
        neutral_unresolved = np.zeros(state_count)
        neutral_unresolved_dn = np.zeros(state_count)
        neutral_escape = np.zeros(state_count)
        neutral_escape_dn = np.zeros(state_count)
        below_support = np.zeros(state_count)
        below_support_dn = np.zeros(state_count)
        maximum_row_residual = 0.0
        maximum_collision_energy_residual = 0.0
        maximum_verlet_defect = 0.0

        for flat in range(state_count):
            ip, direction, ix, ie, it = np.unravel_index(flat, shape)
            velocity = self._velocity_from_state(
                direction, energy_axis[ie], transverse_axis[it])
            kernel = self._trajectory_kernel(
                velocity, position_axis[ix], self.phase_axis_rad[ip],
                impact_cache)
            maximum_collision_energy_residual = max(
                maximum_collision_energy_residual,
                kernel.maximum_collision_energy_residual)
            maximum_verlet_defect = max(
                maximum_verlet_defect,
                kernel.maximum_verlet_energy_defect_eV)
            collision[flat] = kernel.collision_probability
            collision_dn[flat] = (
                kernel.collision_probability_density_derivative)
            if kernel.endpoint_kind == "arrived":
                arrival[flat] = kernel.survival
                arrival_dn[flat] = kernel.survival_density_derivative
                for output, weight in self._output_deposition(
                    kernel.endpoint_velocity, kernel.endpoint_phase_rad,
                    energy_axis, transverse_axis,
                ):
                    b_row.append(flat)
                    b_col.append(output)
                    b_data.append(kernel.survival * weight)
                    db_data.append(
                        kernel.survival_density_derivative * weight)
            else:
                escaped[flat] = kernel.survival
                escaped_dn[flat] = kernel.survival_density_derivative
            for branch in kernel.branches:
                if branch.charge_exchange:
                    charge_exchange[flat] += branch.coefficient
                    charge_exchange_dn[flat] += branch.density_derivative
                if branch.below_angular_support:
                    below_support[flat] += branch.coefficient
                    below_support_dn[flat] += branch.density_derivative
                energy = float(np.dot(
                    branch.neutral_velocity, branch.neutral_velocity))
                birth_energy[flat] += branch.coefficient * energy
                birth_energy_dn[flat] += branch.density_derivative * energy
                for destination, weight in self._state_deposition(
                    branch.ion_velocity, branch.position_fraction,
                    branch.phase_rad, position_axis, energy_axis,
                    transverse_axis,
                ):
                    q_row.append(flat)
                    q_col.append(destination)
                    q_data.append(branch.coefficient * weight)
                    dq_data.append(branch.density_derivative * weight)
                neutral = self._neutral_lower_bound(branch)
                if neutral[0] == "escaped":
                    neutral_escape[flat] += neutral[1]
                    neutral_escape_dn[flat] += neutral[2]
                else:
                    _, value, derivative, unresolved, unresolved_dn, phase = neutral
                    neutral_arrival[flat] += value
                    neutral_arrival_dn[flat] += derivative
                    neutral_unresolved[flat] += unresolved
                    neutral_unresolved_dn[flat] += unresolved_dn
                    for output, weight in self._output_deposition(
                        branch.neutral_velocity, phase,
                        energy_axis, transverse_axis,
                    ):
                        n_row.append(flat)
                        n_col.append(output)
                        n_data.append(value * weight)
                        dn_data.append(derivative * weight)
            maximum_row_residual = max(
                maximum_row_residual,
                abs(arrival[flat] + escaped[flat]
                    + sum(branch.coefficient for branch in kernel.branches)
                    - 1.0),
            )

        q = coo_matrix(
            (q_data, (q_row, q_col)), shape=(state_count, state_count)).tocsc()
        dq = coo_matrix(
            (dq_data, (q_row, q_col)), shape=(state_count, state_count)).tocsc()
        output = coo_matrix(
            (b_data, (b_row, b_col)), shape=(state_count, output_count)).tocsc()
        output_dn = coo_matrix(
            (db_data, (b_row, b_col)), shape=(state_count, output_count)).tocsc()
        neutral_output = coo_matrix(
            (n_data, (n_row, n_col)), shape=(state_count, output_count)).tocsc()
        neutral_output_dn = coo_matrix(
            (dn_data, (n_row, n_col)), shape=(state_count, output_count)).tocsc()

        source = np.zeros(state_count)
        source_dn = np.zeros(state_count)
        direct_output = np.zeros(output_count)
        direct_output_dn = np.zeros(output_count)
        direct_neutral = np.zeros(output_count)
        direct_neutral_dn = np.zeros(output_count)
        names = (
            "arrival", "escaped", "collision", "charge_exchange",
            "birth_energy", "neutral_arrival", "neutral_unresolved",
            "neutral_escape", "below_support", "optical_depth",
        )
        direct = {name: np.zeros(2) for name in names}
        phase_weight = 1.0 / len(self.phase_axis_rad)
        for phase in self.phase_axis_rad:
            for velocity, thermal_weight in initial:
                weight = phase_weight * thermal_weight
                kernel = self._trajectory_kernel(
                    velocity, 0.0, phase, impact_cache)
                maximum_collision_energy_residual = max(
                    maximum_collision_energy_residual,
                    kernel.maximum_collision_energy_residual)
                maximum_verlet_defect = max(
                    maximum_verlet_defect,
                    kernel.maximum_verlet_energy_defect_eV)
                direct["collision"] += weight * np.array([
                    kernel.collision_probability,
                    kernel.collision_probability_density_derivative,
                ])
                direct["optical_depth"][0] += (
                    weight * self.gas_number_density_m3
                    * kernel.unit_optical_depth)
                if kernel.endpoint_kind == "arrived":
                    direct["arrival"] += weight * np.array([
                        kernel.survival,
                        kernel.survival_density_derivative,
                    ])
                    for destination, coefficient in self._output_deposition(
                        kernel.endpoint_velocity, kernel.endpoint_phase_rad,
                        energy_axis, transverse_axis,
                    ):
                        direct_output[destination] += (
                            weight * kernel.survival * coefficient)
                        direct_output_dn[destination] += (
                            weight * kernel.survival_density_derivative
                            * coefficient)
                else:
                    direct["escaped"] += weight * np.array([
                        kernel.survival,
                        kernel.survival_density_derivative,
                    ])
                for branch in kernel.branches:
                    if branch.charge_exchange:
                        direct["charge_exchange"] += weight * np.array([
                            branch.coefficient, branch.density_derivative])
                    if branch.below_angular_support:
                        direct["below_support"] += weight * np.array([
                            branch.coefficient, branch.density_derivative])
                    neutral_energy = float(np.dot(
                        branch.neutral_velocity, branch.neutral_velocity))
                    direct["birth_energy"] += weight * neutral_energy * np.array([
                        branch.coefficient, branch.density_derivative])
                    for destination, coefficient in self._state_deposition(
                        branch.ion_velocity, branch.position_fraction,
                        branch.phase_rad, position_axis, energy_axis,
                        transverse_axis,
                    ):
                        source[destination] += (
                            weight * branch.coefficient * coefficient)
                        source_dn[destination] += (
                            weight * branch.density_derivative * coefficient)
                    neutral = self._neutral_lower_bound(branch)
                    if neutral[0] == "escaped":
                        direct["neutral_escape"] += weight * np.array([
                            neutral[1], neutral[2]])
                    else:
                        _, value, derivative, unresolved, dn_unresolved, nphase = neutral
                        direct["neutral_arrival"] += weight * np.array([
                            value, derivative])
                        direct["neutral_unresolved"] += weight * np.array([
                            unresolved, dn_unresolved])
                        for destination, coefficient in self._output_deposition(
                            branch.neutral_velocity, nphase,
                            energy_axis, transverse_axis,
                        ):
                            direct_neutral[destination] += (
                                weight * value * coefficient)
                            direct_neutral_dn[destination] += (
                                weight * derivative * coefficient)

        operator = eye(state_count, format="csc") - q.T
        factor = splu(operator)
        visits = factor.solve(source)
        visits_dn = factor.solve(source_dn + dq.T @ visits)
        solve_residual = float(np.linalg.norm(
            operator @ visits - source, ord=np.inf)
            / max(np.linalg.norm(source, ord=np.inf), 1.0))
        tangent_residual = float(np.linalg.norm(
            operator @ visits_dn - (source_dn + dq.T @ visits), ord=np.inf)
            / max(np.linalg.norm(source_dn + dq.T @ visits, ord=np.inf), 1.0))

        bins = direct_output + np.asarray(output.T @ visits).ravel()
        bins_dn = (
            direct_output_dn
            + np.asarray(output_dn.T @ visits).ravel()
            + np.asarray(output.T @ visits_dn).ravel())
        neutral_bins = (
            direct_neutral + np.asarray(neutral_output.T @ visits).ravel())
        neutral_bins_dn = (
            direct_neutral_dn
            + np.asarray(neutral_output_dn.T @ visits).ravel()
            + np.asarray(neutral_output.T @ visits_dn).ravel())

        def total(name, array, derivative):
            return (
                float(direct[name][0] + np.dot(visits, array)),
                float(direct[name][1] + np.dot(visits_dn, array)
                      + np.dot(visits, derivative)),
            )

        totals = {
            "arrival": total("arrival", arrival, arrival_dn),
            "escaped": total("escaped", escaped, escaped_dn),
            "collision": total("collision", collision, collision_dn),
            "charge_exchange": total(
                "charge_exchange", charge_exchange, charge_exchange_dn),
            "birth_energy": total(
                "birth_energy", birth_energy, birth_energy_dn),
            "neutral_arrival": total(
                "neutral_arrival", neutral_arrival, neutral_arrival_dn),
            "neutral_unresolved": total(
                "neutral_unresolved", neutral_unresolved,
                neutral_unresolved_dn),
            "neutral_escape": total(
                "neutral_escape", neutral_escape, neutral_escape_dn),
            "below_support": total(
                "below_support", below_support, below_support_dn),
        }
        return {
            "energy_axis": energy_axis,
            "transverse_axis": transverse_axis,
            "bins": bins,
            "bins_dn": bins_dn,
            "neutral_bins": neutral_bins,
            "neutral_bins_dn": neutral_bins_dn,
            "totals": totals,
            "direct_uncollided": float(direct["arrival"][0]),
            "direct_uncollided_dn": float(direct["arrival"][1]),
            "mean_initial_optical_depth": float(direct["optical_depth"][0]),
            "maximum_row_residual": maximum_row_residual,
            "maximum_collision_energy_residual": (
                maximum_collision_energy_residual),
            "maximum_verlet_defect_eV": maximum_verlet_defect,
            "solve_residual": solve_residual,
            "tangent_residual": tangent_residual,
            "density_direction": density_direction,
        }

    def _solve(self, density_direction: float | None):
        direction = 0.0 if density_direction is None else float(density_direction)
        if not math.isfinite(direction):
            raise ValueError("density tangent must be finite")
        assembled = self._assemble_and_solve(direction)
        totals = assembled["totals"]
        arrived, arrived_dn_per_density = totals["arrival"]
        escaped, escaped_dn_per_density = totals["escaped"]
        unresolved = max(0.0, 1.0 - arrived - escaped)
        unresolved_dn_per_density = -arrived_dn_per_density - escaped_dn_per_density
        probability_residual = abs(arrived + escaped + unresolved - 1.0)

        output_shape = (
            len(self.phase_axis_rad), len(assembled["energy_axis"]),
            len(assembled["transverse_axis"]),
        )
        velocity = []
        raw_weight = []
        raw_tangent = []
        phase_value = []
        neutral_velocity = []
        neutral_weight = []
        neutral_phase = []
        azimuth = (
            2.0
            * np.pi
            * (np.arange(int(self.output_azimuth_order)) + 0.5)
            / int(self.output_azimuth_order)
        )
        for flat, probability in enumerate(assembled["bins"]):
            if probability <= 0.0:
                continue
            ip, ie, it = np.unravel_index(flat, output_shape)
            energy = float(assembled["energy_axis"][ie])
            fraction = float(assembled["transverse_axis"][it])
            vt = math.sqrt(max(energy * fraction, 0.0))
            vn = math.sqrt(max(energy * (1.0 - fraction), 0.0))
            for phi in azimuth:
                velocity.append(np.array([
                    vt * math.cos(phi), vt * math.sin(phi), vn]))
                raw_weight.append(probability / len(azimuth))
                raw_tangent.append(
                    assembled["bins_dn"][flat] * direction / len(azimuth))
                phase_value.append(float(self.phase_axis_rad[ip]))
        raw_weight_array = np.asarray(raw_weight)
        raw_tangent_array = np.asarray(raw_tangent)
        normalized = raw_weight_array / arrived
        arrived_tangent = arrived_dn_per_density * direction
        normalized_tangent = (
            raw_tangent_array * arrived
            - raw_weight_array * arrived_tangent
        ) / arrived ** 2
        distribution = CollisionalIonEnergyAngleDistribution(
            velocity_sqrt_eV=np.asarray(velocity),
            weight=normalized,
            entry_phase_rad=np.asarray(phase_value),
        )

        neutral_arrival = totals["neutral_arrival"][0]
        neutral_distribution = None
        if neutral_arrival > 0.0:
            for flat, probability in enumerate(assembled["neutral_bins"]):
                if probability <= 0.0:
                    continue
                ip, ie, it = np.unravel_index(flat, output_shape)
                energy = float(assembled["energy_axis"][ie])
                fraction = float(assembled["transverse_axis"][it])
                vt = math.sqrt(max(energy * fraction, 0.0))
                vn = math.sqrt(max(energy * (1.0 - fraction), 0.0))
                for phi in azimuth:
                    neutral_velocity.append(np.array([
                        vt * math.cos(phi), vt * math.sin(phi), vn]))
                    neutral_weight.append(probability / len(azimuth))
                    neutral_phase.append(float(self.phase_axis_rad[ip]))
            neutral_distribution = CollisionalIonEnergyAngleDistribution(
                velocity_sqrt_eV=np.asarray(neutral_velocity),
                weight=np.asarray(neutral_weight) / neutral_arrival,
                entry_phase_rad=np.asarray(neutral_phase),
            )
        neutral_birth = totals["collision"][0]
        neutral_unresolved = totals["neutral_unresolved"][0]
        neutral_escape = totals["neutral_escape"][0]
        neutral_lineage_residual = abs(
            neutral_arrival + neutral_unresolved + neutral_escape
            - neutral_birth)
        reference = self.sheath.ion_impact_energies(
            self.phase_axis_rad,
            steps_per_period=int(self.steps_per_period),
            steps_per_transit=int(self.steps_per_transit),
            max_periods=float(self.maximum_transit_periods),
        )
        solution = DeterministicCollisionalSheathSolution(
            distribution=distribution,
            resolved_fast_neutral_distribution=neutral_distribution,
            source_ion_flux_m2_s=float(self.source_ion_flux_m2_s),
            arriving_ion_flux_m2_s=self.source_ion_flux_m2_s * arrived,
            resolved_fast_neutral_flux_m2_s=(
                self.source_ion_flux_m2_s * neutral_arrival),
            ion_arrival_probability=arrived,
            unresolved_probability=unresolved,
            escaped_probability=escaped,
            uncollided_arrival_probability=assembled["direct_uncollided"],
            expected_collision_count_lower_bound=neutral_birth,
            expected_charge_exchange_count_lower_bound=(
                totals["charge_exchange"][0]),
            expected_fast_neutral_birth_count_lower_bound=neutral_birth,
            expected_fast_neutral_birth_energy_lower_bound_eV_per_source_ion=(
                totals["birth_energy"][0]),
            resolved_fast_neutral_arrivals_per_source_ion=neutral_arrival,
            unresolved_fast_neutral_collisions_per_source_ion=(
                neutral_unresolved),
            escaped_fast_neutrals_per_source_ion=neutral_escape,
            fast_neutral_lineage_ledger_relative_residual=(
                neutral_lineage_residual),
            maximum_resolved_energy_ledger_relative_residual=(
                assembled["maximum_collision_energy_residual"]),
            probability_ledger_relative_residual=probability_residual,
            collisionless_reference_mean_normal_energy_eV=float(
                np.mean(reference)),
            mean_total_optical_depth=assembled["mean_initial_optical_depth"],
            maximum_total_optical_depth=assembled["mean_initial_optical_depth"],
            below_born_mayer_support_collision_probability_lower_bound=(
                totals["below_support"][0]),
            model_source=(
                self.collision_model.source + "; " + self.sheath.source),
            provenance={
                **dict(self.collision_model.provenance),
                **dict(self.sheath.certification_receipt()),
                **dict(self.provenance),
                "solver": "moving-field implicit discrete ordinates",
                "RF_phase_is_kinetic_state": True,
                "charge_exchange_birth_phase_resolved": True,
                "moving_sheath_self_consistency_closed": True,
                "ion_collision_order_closed": True,
                "fast_neutral_transport_closed": False,
                "phase_node_count": int(self.phase_node_count),
                "position_node_count": int(self.position_node_count),
                "total_energy_node_count": int(self.total_energy_node_count),
                "transverse_fraction_node_count": int(
                    self.transverse_fraction_node_count),
                "steps_per_period": int(self.steps_per_period),
                "steps_per_transit": int(self.steps_per_transit),
                "collision_event_quadrature_order": int(
                    self.collision_event_quadrature_order),
                "linear_solve_relative_residual": assembled["solve_residual"],
                "tangent_linear_solve_relative_residual": (
                    assembled["tangent_residual"]),
                "maximum_row_probability_residual": (
                    assembled["maximum_row_residual"]),
                "maximum_verlet_work_defect_eV": (
                    assembled["maximum_verlet_defect_eV"]),
                "supports_generator_power_inversion": False,
                "supports_feature_depth": False,
                "feature_depth_used": False,
            },
        )
        tangent = None
        if density_direction is not None:
            neutral_lineage_tangent_residual = direction * (
                totals["neutral_arrival"][1]
                + totals["neutral_unresolved"][1]
                + totals["neutral_escape"][1]
                - totals["collision"][1]
            )
            tangent = CollisionalSheathDensityTangent(
                gas_number_density_tangent_m3=direction,
                distribution_weight_tangent=normalized_tangent,
                ion_arrival_probability_tangent=arrived_tangent,
                unresolved_probability_tangent=(
                    unresolved_dn_per_density * direction),
                escaped_probability_tangent=escaped_dn_per_density * direction,
                uncollided_arrival_probability_tangent=(
                    assembled["direct_uncollided_dn"] * direction),
                expected_collision_count_tangent=(
                    totals["collision"][1] * direction),
                expected_charge_exchange_count_tangent=(
                    totals["charge_exchange"][1] * direction),
                expected_fast_neutral_birth_count_tangent=(
                    totals["collision"][1] * direction),
                expected_fast_neutral_birth_energy_tangent_eV=(
                    totals["birth_energy"][1] * direction),
                resolved_fast_neutral_arrivals_tangent=(
                    totals["neutral_arrival"][1] * direction),
                unresolved_fast_neutral_collisions_tangent=(
                    totals["neutral_unresolved"][1] * direction),
                escaped_fast_neutrals_tangent=(
                    totals["neutral_escape"][1] * direction),
                fast_neutral_lineage_ledger_tangent_residual=(
                    neutral_lineage_tangent_residual),
                mean_impact_energy_tangent_eV=float(np.dot(
                    normalized_tangent, distribution.energy_eV)),
                probability_ledger_tangent_residual=(
                    direction * (
                        arrived_dn_per_density + escaped_dn_per_density
                        + unresolved_dn_per_density)),
            )
        return solution, tangent

    def solve(self) -> DeterministicCollisionalSheathSolution:
        solution, _ = self._solve(None)
        return solution

    def density_jvp(self, gas_number_density_tangent_m3: float):
        solution, tangent = self._solve(float(gas_number_density_tangent_m3))
        assert tangent is not None
        return solution, tangent


def certify_moving_sheath_convergence(
    coarse: DeterministicCollisionalSheathSolution,
    fine: DeterministicCollisionalSheathSolution,
    *,
    mean_energy_relative_limit: float = 2.0e-2,
    rms_angle_relative_limit: float = 3.0e-2,
    collision_count_relative_limit: float = 3.0e-2,
    arrival_probability_relative_limit: float = 1.0e-2,
    probability_ledger_limit: float = 2.0e-11,
) -> MovingSheathConvergenceReceipt:
    for value in (coarse, fine):
        if (
            value.provenance.get("RF_phase_is_kinetic_state") is not True
            or value.provenance.get("ion_collision_order_closed") is not True
        ):
            raise ValueError("both inputs must be moving-sheath solutions")

    def relative(a: float, b: float) -> float:
        return float(abs(b - a) / max(abs(a), abs(b), np.finfo(float).tiny))

    values = {
        "mean_energy_relative_change": relative(
            coarse.distribution.mean_energy_eV,
            fine.distribution.mean_energy_eV),
        "rms_angle_relative_change": relative(
            math.sqrt(coarse.distribution.mean_squared_polar_angle_rad2),
            math.sqrt(fine.distribution.mean_squared_polar_angle_rad2)),
        "collision_count_relative_change": relative(
            coarse.expected_collision_count_lower_bound,
            fine.expected_collision_count_lower_bound),
        "arrival_probability_relative_change": relative(
            coarse.ion_arrival_probability, fine.ion_arrival_probability),
        "probability_ledger_residual": max(
            coarse.probability_ledger_relative_residual,
            fine.probability_ledger_relative_residual),
    }
    limits = {
        "mean_energy_relative_change": float(mean_energy_relative_limit),
        "rms_angle_relative_change": float(rms_angle_relative_limit),
        "collision_count_relative_change": float(collision_count_relative_limit),
        "arrival_probability_relative_change": float(
            arrival_probability_relative_limit),
        "probability_ledger_residual": float(probability_ledger_limit),
    }
    return MovingSheathConvergenceReceipt(
        coarse=coarse, fine=fine, **values,
        passed=all(values[name] <= limit for name, limit in limits.items()),
        limits=limits,
    )
