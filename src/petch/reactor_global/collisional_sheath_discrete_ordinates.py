"""Bounded deterministic discrete-ordinates transport through an RF sheath.

This module closes the ion collision-order expansion in
``collisional_sheath`` without a particle Monte Carlo and without an
exponentially growing path tree.  For each RF entry phase it constructs the
finite-state linear Boltzmann operator on fixed potential-coordinate,
tangential-energy, normal-energy, and direction ordinates.  The complete
Neumann series is evaluated by the absorbing-system solve

``(I - Q.T) x = s``.

Linear deposition preserves probability, kinetic-energy components, and the
electrostatic invariant in expectation.  Backscattered ions are followed to
the plasma edge or through their electrostatic turning point.  The exact JVP
of the *discrete* implicit solve with respect to gas density is evaluated from

``(I - Q.T) dx = ds + dQ.T @ x``.

The fast-neutral boundary remains deliberately conservative: every neutral
born in an ion collision is transported through its no-further-collision
branch, while subsequent neutral--neutral collisions remain an explicit
unresolved ledger.  Consequently this provider closes ion transport but does
not promote itself to a feature-depth boundary.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu

from ..sheath import CollisionlessWaveformSheath
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


_NEUTRAL_TRANSPORT_RESULT = tuple[
    str, float, float, float, float
] | tuple[str, float, float]


@dataclass(frozen=True)
class _CollisionBranch:
    ion_velocity: np.ndarray
    neutral_velocity: np.ndarray
    position_fraction: float
    coefficient: float
    density_derivative: float
    charge_exchange: bool
    below_angular_support: bool


@dataclass(frozen=True)
class _StateKernel:
    endpoint_kind: str
    endpoint_velocity: np.ndarray
    survival: float
    survival_density_derivative: float
    collision_probability: float
    collision_probability_density_derivative: float
    branches: tuple[_CollisionBranch, ...]
    unit_optical_depth: float
    maximum_energy_residual: float


@dataclass(frozen=True)
class _PhaseSolution:
    total_energy_axis_eV: np.ndarray
    transverse_fraction_axis: np.ndarray
    arrival_bin_probability: np.ndarray
    arrival_bin_density_derivative: np.ndarray
    neutral_arrival_bin_probability: np.ndarray
    neutral_arrival_bin_density_derivative: np.ndarray
    arrival_probability: float
    arrival_probability_density_derivative: float
    escape_probability: float
    escape_probability_density_derivative: float
    uncollided_arrival_probability: float
    uncollided_arrival_probability_density_derivative: float
    collision_count: float
    collision_count_density_derivative: float
    charge_exchange_count: float
    charge_exchange_count_density_derivative: float
    neutral_birth_count: float
    neutral_birth_count_density_derivative: float
    neutral_birth_energy_eV: float
    neutral_birth_energy_density_derivative_eV: float
    neutral_arrival: float
    neutral_arrival_density_derivative: float
    neutral_unresolved: float
    neutral_unresolved_density_derivative: float
    neutral_escape: float
    neutral_escape_density_derivative: float
    below_support_probability: float
    below_support_probability_density_derivative: float
    mean_initial_optical_depth: float
    maximum_energy_residual: float
    maximum_row_probability_residual: float
    linear_solve_relative_residual: float
    tangent_linear_solve_relative_residual: float


@dataclass(frozen=True)
class DeterministicDiscreteOrdinatesRFSheath:
    """Implicit, converged ion transport for a phase-conditioned RF sheath."""

    sheath: CollisionlessWaveformSheath
    collision_model: ArgonBornMayerPhelpsCollisionModel
    gas_number_density_m3: float
    neutral_gas_temperature_K: float
    source_ion_flux_m2_s: float
    phase_count: int = 24
    initial_thermal_radial_order: int = 2
    initial_thermal_azimuth_order: int = 4
    potential_node_count: int = 9
    total_energy_node_count: int = 7
    transverse_fraction_node_count: int = 13
    position_quadrature_order: int = 3
    hazard_quadrature_order: int = 5
    impact_quadrature_order: int = 2
    collision_azimuth_order: int = 4
    steps_per_period: int = 256
    steps_per_transit: int = 256
    provenance: Mapping[str, object] = field(default_factory=lambda: {
        "state_coordinate": (
            "signed direction x Child-potential coordinate x tangential "
            "energy x normal energy"
        ),
        "linear_system": "absorbing discrete-ordinates Boltzmann operator",
        "feature_depth_used": False,
    })

    def __post_init__(self):
        if not isinstance(self.sheath, CollisionlessWaveformSheath):
            raise TypeError("a CollisionlessWaveformSheath is required")
        if not isinstance(
            self.collision_model, ArgonBornMayerPhelpsCollisionModel
        ):
            raise TypeError("an Ar Born-Mayer collision model is required")
        values = np.asarray([
            self.gas_number_density_m3,
            self.neutral_gas_temperature_K,
            self.source_ion_flux_m2_s,
        ], dtype=float)
        orders = (
            self.phase_count,
            self.initial_thermal_radial_order,
            self.initial_thermal_azimuth_order,
            self.potential_node_count,
            self.total_energy_node_count,
            self.transverse_fraction_node_count,
            self.position_quadrature_order,
            self.hazard_quadrature_order,
            self.impact_quadrature_order,
            self.collision_azimuth_order,
            self.steps_per_period,
            self.steps_per_transit,
        )
        if (
            np.any(~np.isfinite(values))
            or self.gas_number_density_m3 < 0.0
            or self.neutral_gas_temperature_K <= 0.0
            or self.source_ion_flux_m2_s <= 0.0
            or any(int(value) < 1 for value in orders)
            or int(self.phase_count) < 4
            or int(self.initial_thermal_azimuth_order) < 2
            or int(self.potential_node_count) < 3
            or int(self.total_energy_node_count) < 3
            or int(self.transverse_fraction_node_count) < 3
            or int(self.collision_azimuth_order) % 2 != 0
        ):
            raise ValueError("invalid discrete-ordinates collisional sheath")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def tangential_temperature_eV(self) -> float:
        return (
            BOLTZMANN_J_K * self.neutral_gas_temperature_K / E_CHARGE_C)

    @staticmethod
    def _potential_shape(position_fraction: float) -> float:
        return float(np.clip(position_fraction, 0.0, 1.0) ** (4.0 / 3.0))

    @staticmethod
    def _position_from_potential_shape(potential_shape: float) -> float:
        return float(np.clip(potential_shape, 0.0, 1.0) ** 0.75)

    def _initial_velocity_quadrature(self):
        radial_node, radial_weight = np.polynomial.laguerre.laggauss(
            int(self.initial_thermal_radial_order))
        bohm_energy = 0.5 * float(self.sheath.Te_eV)
        for node, weight in zip(radial_node, radial_weight):
            yield (
                np.array([
                    math.sqrt(self.tangential_temperature_eV * float(node)),
                    0.0,
                    math.sqrt(bohm_energy),
                ]),
                float(weight),
            )

    def _velocity_on_leg(
        self,
        velocity: np.ndarray,
        start_fraction: float,
        position_fraction: float,
        phase_gain_eV: float,
        direction_sign: float,
    ) -> np.ndarray:
        value = np.asarray(velocity, dtype=float)
        normal_energy = (
            value[2] ** 2
            + phase_gain_eV * (
                self._potential_shape(position_fraction)
                - self._potential_shape(start_fraction)
            )
        )
        if normal_energy < -2.0e-10:
            raise RuntimeError("trajectory crossed an inaccessible potential")
        return np.array([
            float(value[0]),
            float(value[1]),
            math.copysign(math.sqrt(max(normal_energy, 0.0)), direction_sign),
        ])

    def _trajectory(
        self,
        velocity: np.ndarray,
        start_fraction: float,
        phase_gain_eV: float,
    ) -> tuple[tuple[tuple[float, float, float], ...], str, np.ndarray]:
        """Return ordered monotone legs and the collisionless endpoint."""
        value = np.asarray(velocity, dtype=float)
        start = float(start_fraction)
        if value[2] >= 0.0:
            endpoint = self._velocity_on_leg(
                value, start, 1.0, phase_gain_eV, 1.0)
            return ((start, 1.0, 1.0),), "arrived", endpoint

        start_u = self._potential_shape(start)
        normal_energy = float(value[2] ** 2)
        edge_energy = normal_energy - phase_gain_eV * start_u
        if phase_gain_eV <= 0.0 or edge_energy >= -2.0e-12:
            endpoint = self._velocity_on_leg(
                value, start, 0.0, phase_gain_eV, -1.0)
            return ((start, 0.0, -1.0),), "escaped", endpoint

        turning_u = np.clip(
            start_u - normal_energy / phase_gain_eV, 0.0, start_u)
        turning = self._position_from_potential_shape(float(turning_u))
        endpoint = self._velocity_on_leg(
            value, start, 1.0, phase_gain_eV, 1.0)
        legs = []
        if start - turning > 2.0e-15:
            legs.append((start, turning, -1.0))
        legs.append((turning, 1.0, 1.0))
        return tuple(legs), "arrived", endpoint

    def _unit_hazard_rate(self, velocity: np.ndarray) -> float:
        normal = abs(float(velocity[2]))
        speed = float(np.linalg.norm(velocity))
        if speed <= 0.0:
            return math.inf
        if normal <= np.finfo(float).tiny:
            return math.inf
        return float(
            self.sheath.thickness
            * self.collision_model.total_cross_section_m2(speed * speed)
            * speed / normal
        )

    def _integrated_leg_unit_hazard(
        self,
        velocity: np.ndarray,
        start_fraction: float,
        leg_start: float,
        leg_end: float,
        phase_gain_eV: float,
        direction_sign: float,
    ) -> float:
        if abs(leg_end - leg_start) <= 2.0e-15:
            return 0.0
        node, weight = np.polynomial.legendre.leggauss(
            int(self.hazard_quadrature_order))
        half = 0.5 * (leg_end - leg_start)
        center = 0.5 * (leg_end + leg_start)
        total = 0.0
        for local, quadrature_weight in zip(node, weight):
            position = center + half * float(local)
            propagated = self._velocity_on_leg(
                velocity,
                start_fraction,
                position,
                phase_gain_eV,
                direction_sign,
            )
            total += (
                float(quadrature_weight)
                * self._unit_hazard_rate(propagated)
            )
        return float(abs(half) * total)

    def _state_kernel(
        self,
        velocity: np.ndarray,
        start_fraction: float,
        phase_gain_eV: float,
        impact_cache: dict[float, tuple[np.ndarray, np.ndarray]],
    ) -> _StateKernel:
        density = float(self.gas_number_density_m3)
        legs, endpoint_kind, endpoint_velocity = self._trajectory(
            velocity, start_fraction, phase_gain_eV)
        leg_hazard = [
            self._integrated_leg_unit_hazard(
                velocity,
                start_fraction,
                leg_start,
                leg_end,
                phase_gain_eV,
                direction,
            )
            for leg_start, leg_end, direction in legs
        ]
        unit_optical_depth = float(np.sum(leg_hazard))
        if not math.isfinite(unit_optical_depth):
            raise RuntimeError("non-finite ion optical depth")
        survival = math.exp(-density * unit_optical_depth)
        survival_derivative = -unit_optical_depth * survival
        collision_probability = -math.expm1(-density * unit_optical_depth)
        collision_derivative = -survival_derivative
        if unit_optical_depth <= np.finfo(float).tiny:
            return _StateKernel(
                endpoint_kind=endpoint_kind,
                endpoint_velocity=endpoint_velocity,
                survival=1.0,
                survival_density_derivative=0.0,
                collision_probability=0.0,
                collision_probability_density_derivative=0.0,
                branches=(),
                unit_optical_depth=0.0,
                maximum_energy_residual=0.0,
            )

        position_node, position_weight = np.polynomial.legendre.leggauss(
            int(self.position_quadrature_order))
        event_records = []
        raw_total = 0.0
        raw_total_derivative = 0.0
        completed_hazard = 0.0
        for (
            leg_start,
            leg_end,
            direction,
        ), full_leg_hazard in zip(legs, leg_hazard):
            half = 0.5 * (leg_end - leg_start)
            center = 0.5 * (leg_end + leg_start)
            for local, quadrature_weight in zip(
                position_node, position_weight
            ):
                position = center + half * float(local)
                event_velocity = self._velocity_on_leg(
                    velocity,
                    start_fraction,
                    position,
                    phase_gain_eV,
                    direction,
                )
                partial_hazard = self._integrated_leg_unit_hazard(
                    velocity,
                    start_fraction,
                    leg_start,
                    position,
                    phase_gain_eV,
                    direction,
                )
                cumulative_hazard = completed_hazard + partial_hazard
                raw = (
                    abs(half)
                    * float(quadrature_weight)
                    * self._unit_hazard_rate(event_velocity)
                    * math.exp(-density * cumulative_hazard)
                )
                raw_derivative = -cumulative_hazard * raw
                event_records.append((
                    position,
                    event_velocity,
                    raw,
                    raw_derivative,
                ))
                raw_total += raw
                raw_total_derivative += raw_derivative
            completed_hazard += full_leg_hazard
        if raw_total <= 0.0 or not math.isfinite(raw_total):
            raise RuntimeError("first-collision quadrature lost its mass")

        azimuths = (
            2.0 * np.pi
            * (np.arange(int(self.collision_azimuth_order)) + 0.5)
            / int(self.collision_azimuth_order)
        )
        branches: list[_CollisionBranch] = []
        maximum_energy_residual = 0.0
        for position, event_velocity, raw, raw_derivative in event_records:
            conditional = raw / raw_total
            conditional_derivative = (
                raw_derivative * raw_total
                - raw * raw_total_derivative
            ) / raw_total ** 2
            event_coefficient = collision_probability * conditional
            event_derivative = (
                collision_derivative * conditional
                + collision_probability * conditional_derivative
            )
            event_energy = float(np.dot(event_velocity, event_velocity))
            cache_key = round(event_energy, 8)
            if cache_key not in impact_cache:
                impact_cache[cache_key] = self.collision_model.impact_quadrature(
                    event_energy, int(self.impact_quadrature_order))
            angles, impact_weight = impact_cache[cache_key]
            below_support = (
                event_energy
                < self.collision_model.born_mayer_minimum_lab_energy_eV
            )
            elastic_probability, charge_exchange_probability = (
                self.collision_model.channel_probabilities(event_energy))
            for angle, angular_weight in zip(angles, impact_weight):
                for azimuth in azimuths:
                    projectile, target = _equal_mass_collision_velocities(
                        event_velocity, float(angle), float(azimuth))
                    cx_projectile, cx_target = (
                        (projectile, target)
                        if not below_support
                        else _equal_mass_collision_velocities(
                            event_velocity, 0.0, float(azimuth))
                    )
                    residual = abs(
                        float(np.dot(projectile, projectile))
                        + float(np.dot(target, target))
                        - event_energy
                    ) / max(event_energy, 1.0)
                    maximum_energy_residual = max(
                        maximum_energy_residual, residual)
                    common = float(angular_weight) / len(azimuths)
                    for (
                        ion_velocity,
                        neutral_velocity,
                        channel_probability,
                        charge_exchange,
                    ) in (
                        (
                            projectile,
                            target,
                            elastic_probability,
                            False,
                        ),
                        (
                            cx_target,
                            cx_projectile,
                            charge_exchange_probability,
                            True,
                        ),
                    ):
                        fraction = channel_probability * common
                        branches.append(_CollisionBranch(
                            ion_velocity=_canonical_axisymmetric_velocity(
                                ion_velocity),
                            neutral_velocity=_canonical_axisymmetric_velocity(
                                neutral_velocity),
                            position_fraction=float(position),
                            coefficient=event_coefficient * fraction,
                            density_derivative=event_derivative * fraction,
                            charge_exchange=charge_exchange,
                            below_angular_support=below_support,
                        ))
        branch_probability = float(sum(
            branch.coefficient for branch in branches))
        branch_derivative = float(sum(
            branch.density_derivative for branch in branches))
        if (
            abs(branch_probability - collision_probability) > 2.0e-12
            or abs(branch_derivative - collision_derivative) > 2.0e-28
        ):
            raise RuntimeError("collision channel quadrature lost probability")
        return _StateKernel(
            endpoint_kind=endpoint_kind,
            endpoint_velocity=endpoint_velocity,
            survival=survival,
            survival_density_derivative=survival_derivative,
            collision_probability=collision_probability,
            collision_probability_density_derivative=collision_derivative,
            branches=tuple(branches),
            unit_optical_depth=unit_optical_depth,
            maximum_energy_residual=maximum_energy_residual,
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

    def _state_deposition(
        self,
        velocity: np.ndarray,
        position_fraction: float,
        potential_axis: np.ndarray,
        normalized_energy_axis: np.ndarray,
        transverse_fraction_axis: np.ndarray,
        *,
        phase_gain_eV: float,
        initial_energy_cap_eV: float,
    ):
        sign_index = 0 if float(velocity[2]) < 0.0 else 1
        total_energy = float(np.dot(velocity, velocity))
        transverse_fraction = (
            0.0 if total_energy <= np.finfo(float).tiny
            else float(np.dot(velocity[:2], velocity[:2])) / total_energy
        )
        for iu, wu in self._linear_deposit(
            potential_axis, self._potential_shape(position_fraction)
        ):
            local_cap = (
                initial_energy_cap_eV
                + phase_gain_eV * float(potential_axis[iu])
            )
            normalized_energy = (
                0.0 if local_cap <= np.finfo(float).tiny
                else total_energy / local_cap
            )
            for it, wt in self._linear_deposit(
                normalized_energy_axis, normalized_energy
            ):
                for iz, wz in self._linear_deposit(
                    transverse_fraction_axis, transverse_fraction
                ):
                    flat = int(np.ravel_multi_index(
                        (sign_index, iu, it, iz),
                        (
                            2,
                            len(potential_axis),
                            len(normalized_energy_axis),
                            len(transverse_fraction_axis),
                        ),
                    ))
                    yield flat, float(wu * wt * wz)

    def _output_deposition(
        self,
        velocity: np.ndarray,
        total_energy_axis: np.ndarray,
        transverse_fraction_axis: np.ndarray,
    ):
        total_energy = float(np.dot(velocity, velocity))
        transverse_fraction = (
            0.0 if total_energy <= np.finfo(float).tiny
            else float(np.dot(velocity[:2], velocity[:2])) / total_energy
        )
        for it, wt in self._linear_deposit(
            total_energy_axis, total_energy
        ):
            for iz, wz in self._linear_deposit(
                transverse_fraction_axis, transverse_fraction
            ):
                flat = int(np.ravel_multi_index(
                    (it, iz), (
                        len(total_energy_axis),
                        len(transverse_fraction_axis),
                    )))
                yield flat, float(wt * wz)

    def _neutral_no_further_collision(
        self,
        branch: _CollisionBranch,
    ) -> _NEUTRAL_TRANSPORT_RESULT:
        velocity = branch.neutral_velocity
        speed = float(np.linalg.norm(velocity))
        normal = float(velocity[2])
        if speed <= 0.0 or normal <= 0.0:
            return "escaped", branch.coefficient, branch.density_derivative
        unit_hazard = (
            self.sheath.thickness
            * (1.0 - branch.position_fraction)
            * speed / normal
            * self.collision_model.total_cross_section_m2(speed * speed)
        )
        survival = math.exp(-self.gas_number_density_m3 * unit_hazard)
        survival_derivative = -unit_hazard * survival
        arrival = branch.coefficient * survival
        arrival_derivative = (
            branch.density_derivative * survival
            + branch.coefficient * survival_derivative
        )
        unresolved = branch.coefficient * (1.0 - survival)
        unresolved_derivative = (
            branch.density_derivative * (1.0 - survival)
            - branch.coefficient * survival_derivative
        )
        return (
            "arrived",
            arrival,
            arrival_derivative,
            unresolved,
            unresolved_derivative,
        )

    def _solve_phase(self, phase_gain_eV: float) -> _PhaseSolution:
        gain = float(phase_gain_eV)
        initial_nodes = tuple(self._initial_velocity_quadrature())
        initial_max_energy = max(
            float(np.dot(velocity, velocity)) for velocity, _ in initial_nodes)
        energy_cap = gain + initial_max_energy
        potential_axis = np.linspace(0.0, 1.0, int(self.potential_node_count))
        # At fixed potential coordinate the electrostatic invariant bounds the
        # local kinetic energy by ``E_initial + gain*u``.  Normalizing total
        # energy by that local cap both enforces the invariant and avoids the
        # severe low-u resolution loss of one global energy axis.
        normalized_energy_axis = np.linspace(
            0.0, 1.0, int(self.total_energy_node_count))
        wafer_energy_axis = energy_cap * normalized_energy_axis
        # ``1-mu^2`` is the transverse-energy fraction.  A cosine ordinate
        # resolves the narrow normal-incidence core far better than a uniform
        # fraction grid while retaining both grazing and normal endpoints.
        mu_axis = np.linspace(
            0.0, 1.0, int(self.transverse_fraction_node_count))
        transverse_fraction_axis = 1.0 - mu_axis[::-1] ** 2
        shape = (
            2,
            len(potential_axis),
            len(normalized_energy_axis),
            len(transverse_fraction_axis),
        )
        state_count = int(np.prod(shape))
        output_count = (
            len(wafer_energy_axis) * len(transverse_fraction_axis))
        impact_cache: dict[float, tuple[np.ndarray, np.ndarray]] = {}

        q_row: list[int] = []
        q_col: list[int] = []
        q_data: list[float] = []
        dq_data: list[float] = []
        b_row: list[int] = []
        b_col: list[int] = []
        b_data: list[float] = []
        db_data: list[float] = []
        n_row: list[int] = []
        n_col: list[int] = []
        n_data: list[float] = []
        dn_data: list[float] = []
        arrival = np.zeros(state_count)
        arrival_dn = np.zeros(state_count)
        escaped = np.zeros(state_count)
        escaped_dn = np.zeros(state_count)
        collision = np.zeros(state_count)
        collision_dn = np.zeros(state_count)
        charge_exchange = np.zeros(state_count)
        charge_exchange_dn = np.zeros(state_count)
        neutral_birth_energy = np.zeros(state_count)
        neutral_birth_energy_dn = np.zeros(state_count)
        neutral_arrival = np.zeros(state_count)
        neutral_arrival_dn = np.zeros(state_count)
        neutral_unresolved = np.zeros(state_count)
        neutral_unresolved_dn = np.zeros(state_count)
        neutral_escape = np.zeros(state_count)
        neutral_escape_dn = np.zeros(state_count)
        below_support = np.zeros(state_count)
        below_support_dn = np.zeros(state_count)
        maximum_energy_residual = 0.0
        maximum_row_residual = 0.0

        for flat in range(state_count):
            sign_index, iu, it, iz = np.unravel_index(flat, shape)
            sign = -1.0 if sign_index == 0 else 1.0
            position = self._position_from_potential_shape(
                float(potential_axis[iu]))
            local_cap = initial_max_energy + gain * float(potential_axis[iu])
            total_energy = local_cap * float(normalized_energy_axis[it])
            transverse_fraction = float(transverse_fraction_axis[iz])
            velocity = np.array([
                math.sqrt(total_energy * transverse_fraction),
                0.0,
                sign * math.sqrt(total_energy * (1.0 - transverse_fraction)),
            ])
            kernel = self._state_kernel(
                velocity, position, gain, impact_cache)
            maximum_energy_residual = max(
                maximum_energy_residual, kernel.maximum_energy_residual)
            if kernel.endpoint_kind == "arrived":
                arrival[flat] = kernel.survival
                arrival_dn[flat] = kernel.survival_density_derivative
                for output, fraction in self._output_deposition(
                    kernel.endpoint_velocity,
                    wafer_energy_axis,
                    transverse_fraction_axis,
                ):
                    b_row.append(flat)
                    b_col.append(output)
                    b_data.append(kernel.survival * fraction)
                    db_data.append(
                        kernel.survival_density_derivative * fraction)
            else:
                escaped[flat] = kernel.survival
                escaped_dn[flat] = kernel.survival_density_derivative

            collision[flat] = kernel.collision_probability
            collision_dn[flat] = (
                kernel.collision_probability_density_derivative)
            row_transition_probability = 0.0
            for branch in kernel.branches:
                if branch.charge_exchange:
                    charge_exchange[flat] += branch.coefficient
                    charge_exchange_dn[flat] += branch.density_derivative
                born_energy = float(np.dot(
                    branch.neutral_velocity, branch.neutral_velocity))
                neutral_birth_energy[flat] += branch.coefficient * born_energy
                neutral_birth_energy_dn[flat] += (
                    branch.density_derivative * born_energy)
                if branch.below_angular_support:
                    below_support[flat] += branch.coefficient
                    below_support_dn[flat] += branch.density_derivative
                for destination, fraction in self._state_deposition(
                    branch.ion_velocity,
                    branch.position_fraction,
                    potential_axis,
                    normalized_energy_axis,
                    transverse_fraction_axis,
                    phase_gain_eV=gain,
                    initial_energy_cap_eV=initial_max_energy,
                ):
                    q_row.append(flat)
                    q_col.append(destination)
                    q_data.append(branch.coefficient * fraction)
                    dq_data.append(branch.density_derivative * fraction)
                    row_transition_probability += branch.coefficient * fraction
                neutral_result = self._neutral_no_further_collision(branch)
                kind = neutral_result[0]
                if kind == "escaped":
                    neutral_escape[flat] += neutral_result[1]
                    neutral_escape_dn[flat] += neutral_result[2]
                else:
                    _, value, derivative, unresolved, unresolved_derivative = (
                        neutral_result)
                    neutral_arrival[flat] += value
                    neutral_arrival_dn[flat] += derivative
                    neutral_unresolved[flat] += unresolved
                    neutral_unresolved_dn[flat] += unresolved_derivative
                    for output, fraction in self._output_deposition(
                        branch.neutral_velocity,
                        wafer_energy_axis,
                        transverse_fraction_axis,
                    ):
                        n_row.append(flat)
                        n_col.append(output)
                        n_data.append(value * fraction)
                        dn_data.append(derivative * fraction)
            maximum_row_residual = max(
                maximum_row_residual,
                abs(
                    arrival[flat] + escaped[flat]
                    + row_transition_probability
                    - 1.0
                ),
            )

        matrix_shape = (state_count, state_count)
        q = coo_matrix(
            (q_data, (q_row, q_col)), shape=matrix_shape).tocsr()
        dq = coo_matrix(
            (dq_data, (q_row, q_col)), shape=matrix_shape).tocsr()
        b = coo_matrix(
            (b_data, (b_row, b_col)),
            shape=(state_count, output_count),
        ).tocsr()
        db = coo_matrix(
            (db_data, (b_row, b_col)),
            shape=(state_count, output_count),
        ).tocsr()
        neutral_output = coo_matrix(
            (n_data, (n_row, n_col)),
            shape=(state_count, output_count),
        ).tocsr()
        neutral_output_dn = coo_matrix(
            (dn_data, (n_row, n_col)),
            shape=(state_count, output_count),
        ).tocsr()

        source = np.zeros(state_count)
        source_dn = np.zeros(state_count)
        direct_arrival = np.zeros(output_count)
        direct_arrival_dn = np.zeros(output_count)
        direct_neutral = np.zeros(output_count)
        direct_neutral_dn = np.zeros(output_count)
        direct_scalars = {
            name: [0.0, 0.0]
            for name in (
                "arrival", "escape", "uncollided", "collision", "cx",
                "birth_energy", "neutral_arrival", "neutral_unresolved",
                "neutral_escape", "below_support",
            )
        }
        initial_optical_depth = []
        initial_optical_weight = []
        for initial_velocity, thermal_weight in initial_nodes:
            kernel = self._state_kernel(
                initial_velocity, 0.0, gain, impact_cache)
            maximum_energy_residual = max(
                maximum_energy_residual, kernel.maximum_energy_residual)
            initial_optical_depth.append(
                self.gas_number_density_m3 * kernel.unit_optical_depth)
            initial_optical_weight.append(thermal_weight)
            weight = float(thermal_weight)
            if kernel.endpoint_kind == "arrived":
                direct_scalars["arrival"][0] += weight * kernel.survival
                direct_scalars["arrival"][1] += (
                    weight * kernel.survival_density_derivative)
                direct_scalars["uncollided"][0] += weight * kernel.survival
                direct_scalars["uncollided"][1] += (
                    weight * kernel.survival_density_derivative)
                for output, fraction in self._output_deposition(
                    kernel.endpoint_velocity,
                    wafer_energy_axis,
                    transverse_fraction_axis,
                ):
                    direct_arrival[output] += (
                        weight * kernel.survival * fraction)
                    direct_arrival_dn[output] += (
                        weight
                        * kernel.survival_density_derivative
                        * fraction
                    )
            else:
                direct_scalars["escape"][0] += weight * kernel.survival
                direct_scalars["escape"][1] += (
                    weight * kernel.survival_density_derivative)
            direct_scalars["collision"][0] += (
                weight * kernel.collision_probability)
            direct_scalars["collision"][1] += (
                weight * kernel.collision_probability_density_derivative)
            for branch in kernel.branches:
                coefficient = weight * branch.coefficient
                derivative = weight * branch.density_derivative
                if branch.charge_exchange:
                    direct_scalars["cx"][0] += coefficient
                    direct_scalars["cx"][1] += derivative
                born_energy = float(np.dot(
                    branch.neutral_velocity, branch.neutral_velocity))
                direct_scalars["birth_energy"][0] += (
                    coefficient * born_energy)
                direct_scalars["birth_energy"][1] += derivative * born_energy
                if branch.below_angular_support:
                    direct_scalars["below_support"][0] += coefficient
                    direct_scalars["below_support"][1] += derivative
                for destination, fraction in self._state_deposition(
                    branch.ion_velocity,
                    branch.position_fraction,
                    potential_axis,
                    normalized_energy_axis,
                    transverse_fraction_axis,
                    phase_gain_eV=gain,
                    initial_energy_cap_eV=initial_max_energy,
                ):
                    source[destination] += coefficient * fraction
                    source_dn[destination] += derivative * fraction
                neutral_result = self._neutral_no_further_collision(branch)
                if neutral_result[0] == "escaped":
                    direct_scalars["neutral_escape"][0] += (
                        weight * neutral_result[1])
                    direct_scalars["neutral_escape"][1] += (
                        weight * neutral_result[2])
                else:
                    _, value, dn_value, unresolved, dn_unresolved = (
                        neutral_result)
                    direct_scalars["neutral_arrival"][0] += weight * value
                    direct_scalars["neutral_arrival"][1] += weight * dn_value
                    direct_scalars["neutral_unresolved"][0] += (
                        weight * unresolved)
                    direct_scalars["neutral_unresolved"][1] += (
                        weight * dn_unresolved)
                    for output, fraction in self._output_deposition(
                        branch.neutral_velocity,
                        wafer_energy_axis,
                        transverse_fraction_axis,
                    ):
                        direct_neutral[output] += weight * value * fraction
                        direct_neutral_dn[output] += (
                            weight * dn_value * fraction)

        operator = (eye(state_count, format="csc") - q.T.tocsc())
        factor = splu(operator)
        visits = factor.solve(source)
        tangent_rhs = source_dn + dq.T @ visits
        visits_dn = factor.solve(np.asarray(tangent_rhs))
        solve_residual = float(np.linalg.norm(
            operator @ visits - source, ord=np.inf
        ) / max(np.linalg.norm(source, ord=np.inf), 1.0))
        tangent_solve_residual = float(np.linalg.norm(
            operator @ visits_dn - tangent_rhs, ord=np.inf
        ) / max(np.linalg.norm(tangent_rhs, ord=np.inf), 1.0))
        if (
            np.any(~np.isfinite(visits))
            or np.any(visits < -2.0e-12)
            or solve_residual > 2.0e-10
            or tangent_solve_residual > 2.0e-28
        ):
            raise RuntimeError("discrete-ordinates absorbing solve failed")
        visits = np.maximum(visits, 0.0)

        arrival_bins = direct_arrival + np.asarray(b.T @ visits).ravel()
        arrival_bins_dn = (
            direct_arrival_dn
            + np.asarray(db.T @ visits).ravel()
            + np.asarray(b.T @ visits_dn).ravel()
        )
        neutral_bins = (
            direct_neutral
            + np.asarray(neutral_output.T @ visits).ravel())
        neutral_bins_dn = (
            direct_neutral_dn
            + np.asarray(neutral_output_dn.T @ visits).ravel()
            + np.asarray(neutral_output.T @ visits_dn).ravel()
        )

        def total(name: str, value: np.ndarray, derivative: np.ndarray):
            return (
                direct_scalars[name][0] + float(np.dot(visits, value)),
                direct_scalars[name][1]
                + float(np.dot(visits_dn, value))
                + float(np.dot(visits, derivative)),
            )

        arrival_total = total("arrival", arrival, arrival_dn)
        escape_total = total("escape", escaped, escaped_dn)
        collision_total = total("collision", collision, collision_dn)
        cx_total = total("cx", charge_exchange, charge_exchange_dn)
        birth_energy_total = total(
            "birth_energy", neutral_birth_energy, neutral_birth_energy_dn)
        neutral_arrival_total = total(
            "neutral_arrival", neutral_arrival, neutral_arrival_dn)
        neutral_unresolved_total = total(
            "neutral_unresolved", neutral_unresolved, neutral_unresolved_dn)
        neutral_escape_total = total(
            "neutral_escape", neutral_escape, neutral_escape_dn)
        below_total = total(
            "below_support", below_support, below_support_dn)
        return _PhaseSolution(
            total_energy_axis_eV=wafer_energy_axis,
            transverse_fraction_axis=transverse_fraction_axis,
            arrival_bin_probability=arrival_bins,
            arrival_bin_density_derivative=arrival_bins_dn,
            neutral_arrival_bin_probability=neutral_bins,
            neutral_arrival_bin_density_derivative=neutral_bins_dn,
            arrival_probability=arrival_total[0],
            arrival_probability_density_derivative=arrival_total[1],
            escape_probability=escape_total[0],
            escape_probability_density_derivative=escape_total[1],
            uncollided_arrival_probability=(
                direct_scalars["uncollided"][0]),
            uncollided_arrival_probability_density_derivative=(
                direct_scalars["uncollided"][1]),
            collision_count=collision_total[0],
            collision_count_density_derivative=collision_total[1],
            charge_exchange_count=cx_total[0],
            charge_exchange_count_density_derivative=cx_total[1],
            neutral_birth_count=collision_total[0],
            neutral_birth_count_density_derivative=collision_total[1],
            neutral_birth_energy_eV=birth_energy_total[0],
            neutral_birth_energy_density_derivative_eV=birth_energy_total[1],
            neutral_arrival=neutral_arrival_total[0],
            neutral_arrival_density_derivative=neutral_arrival_total[1],
            neutral_unresolved=neutral_unresolved_total[0],
            neutral_unresolved_density_derivative=neutral_unresolved_total[1],
            neutral_escape=neutral_escape_total[0],
            neutral_escape_density_derivative=neutral_escape_total[1],
            below_support_probability=below_total[0],
            below_support_probability_density_derivative=below_total[1],
            mean_initial_optical_depth=float(np.average(
                initial_optical_depth, weights=initial_optical_weight)),
            maximum_energy_residual=maximum_energy_residual,
            maximum_row_probability_residual=maximum_row_residual,
            linear_solve_relative_residual=solve_residual,
            tangent_linear_solve_relative_residual=tangent_solve_residual,
        )

    def _solve(
        self,
        *,
        gas_number_density_tangent_m3: float | None,
    ) -> tuple[
        DeterministicCollisionalSheathSolution,
        CollisionalSheathDensityTangent | None,
    ]:
        direction = (
            0.0 if gas_number_density_tangent_m3 is None
            else float(gas_number_density_tangent_m3)
        )
        if not math.isfinite(direction):
            raise ValueError("neutral-density tangent must be finite")
        phases = (
            2.0 * np.pi
            * (np.arange(int(self.phase_count), dtype=float) + 0.5)
            / int(self.phase_count)
        )
        collisionless_energy = self.sheath.ion_impact_energies(
            phases,
            steps_per_period=int(self.steps_per_period),
            steps_per_transit=int(self.steps_per_transit),
        )
        bohm_energy = 0.5 * float(self.sheath.Te_eV)
        gains = np.maximum(collisionless_energy - bohm_energy, 0.0)
        phase_cache: dict[float, _PhaseSolution] = {}
        phase_solutions = []
        for gain in gains:
            key = round(float(gain), 10)
            if key not in phase_cache:
                phase_cache[key] = self._solve_phase(float(gain))
            phase_solutions.append(phase_cache[key])

        arriving_velocity = []
        arriving_raw_weight = []
        arriving_raw_tangent = []
        arriving_phase = []
        neutral_velocity = []
        neutral_raw_weight = []
        neutral_phase = []
        phase_weight = 1.0 / int(self.phase_count)
        output_azimuth = (
            2.0 * np.pi
            * (np.arange(int(self.initial_thermal_azimuth_order)) + 0.5)
            / int(self.initial_thermal_azimuth_order)
        )
        for phase, phase_solution in zip(phases, phase_solutions):
            output_shape = (
                len(phase_solution.total_energy_axis_eV),
                len(phase_solution.transverse_fraction_axis),
            )
            for flat, probability in enumerate(
                phase_solution.arrival_bin_probability
            ):
                if probability <= 0.0:
                    continue
                it, iz = np.unravel_index(flat, output_shape)
                total_energy = float(phase_solution.total_energy_axis_eV[it])
                fraction = float(
                    phase_solution.transverse_fraction_axis[iz])
                transverse = math.sqrt(total_energy * fraction)
                normal = math.sqrt(total_energy * (1.0 - fraction))
                raw = phase_weight * probability
                raw_tangent = (
                    phase_weight
                    * phase_solution.arrival_bin_density_derivative[flat]
                    * direction
                )
                for azimuth in output_azimuth:
                    arriving_velocity.append(np.array([
                        transverse * math.cos(float(azimuth)),
                        transverse * math.sin(float(azimuth)),
                        normal,
                    ]))
                    arriving_raw_weight.append(raw / len(output_azimuth))
                    arriving_raw_tangent.append(
                        raw_tangent / len(output_azimuth))
                    arriving_phase.append(float(phase))
            for flat, probability in enumerate(
                phase_solution.neutral_arrival_bin_probability
            ):
                if probability <= 0.0:
                    continue
                it, iz = np.unravel_index(flat, output_shape)
                total_energy = float(phase_solution.total_energy_axis_eV[it])
                fraction = float(
                    phase_solution.transverse_fraction_axis[iz])
                transverse = math.sqrt(total_energy * fraction)
                normal = math.sqrt(total_energy * (1.0 - fraction))
                raw = phase_weight * probability
                for azimuth in output_azimuth:
                    neutral_velocity.append(np.array([
                        transverse * math.cos(float(azimuth)),
                        transverse * math.sin(float(azimuth)),
                        normal,
                    ]))
                    neutral_raw_weight.append(raw / len(output_azimuth))
                    neutral_phase.append(float(phase))

        def average(name: str) -> float:
            return float(np.mean([
                getattr(solution, name) for solution in phase_solutions]))

        def tangent_average(name: str) -> float:
            return direction * average(name)

        arrived = average("arrival_probability")
        arrived_tangent = tangent_average(
            "arrival_probability_density_derivative")
        escaped = average("escape_probability")
        escaped_tangent = tangent_average(
            "escape_probability_density_derivative")
        unresolved = max(0.0, 1.0 - arrived - escaped)
        unresolved_tangent = -arrived_tangent - escaped_tangent
        probability_residual = abs(arrived + escaped + unresolved - 1.0)
        raw_weight = np.asarray(arriving_raw_weight)
        raw_tangent = np.asarray(arriving_raw_tangent)
        normalized = raw_weight / arrived
        normalized_tangent = (
            raw_tangent * arrived - raw_weight * arrived_tangent
        ) / arrived ** 2
        distribution = CollisionalIonEnergyAngleDistribution(
            velocity_sqrt_eV=np.asarray(arriving_velocity),
            weight=normalized,
            entry_phase_rad=np.asarray(arriving_phase),
        )
        neutral_arrival = average("neutral_arrival")
        neutral_distribution = None
        if neutral_arrival > 0.0:
            neutral_distribution = CollisionalIonEnergyAngleDistribution(
                velocity_sqrt_eV=np.asarray(neutral_velocity),
                weight=np.asarray(neutral_raw_weight) / neutral_arrival,
                entry_phase_rad=np.asarray(neutral_phase),
            )
        neutral_birth = average("neutral_birth_count")
        neutral_unresolved = average("neutral_unresolved")
        neutral_escape = average("neutral_escape")
        neutral_lineage_residual = abs(
            neutral_arrival + neutral_unresolved + neutral_escape
            - neutral_birth
        )
        neutral_lineage_tangent_residual = direction * (
            average("neutral_arrival_density_derivative")
            + average("neutral_unresolved_density_derivative")
            + average("neutral_escape_density_derivative")
            - average("neutral_birth_count_density_derivative")
        )
        max_linear_residual = max(
            phase.linear_solve_relative_residual
            for phase in phase_solutions)
        max_tangent_residual = max(
            phase.tangent_linear_solve_relative_residual
            for phase in phase_solutions)
        max_row_residual = max(
            phase.maximum_row_probability_residual
            for phase in phase_solutions)
        max_energy_residual = max(
            phase.maximum_energy_residual for phase in phase_solutions)
        solution = DeterministicCollisionalSheathSolution(
            distribution=distribution,
            resolved_fast_neutral_distribution=neutral_distribution,
            source_ion_flux_m2_s=float(self.source_ion_flux_m2_s),
            arriving_ion_flux_m2_s=float(self.source_ion_flux_m2_s) * arrived,
            resolved_fast_neutral_flux_m2_s=(
                float(self.source_ion_flux_m2_s) * neutral_arrival),
            ion_arrival_probability=arrived,
            unresolved_probability=unresolved,
            escaped_probability=escaped,
            uncollided_arrival_probability=average(
                "uncollided_arrival_probability"),
            expected_collision_count_lower_bound=average("collision_count"),
            expected_charge_exchange_count_lower_bound=average(
                "charge_exchange_count"),
            expected_fast_neutral_birth_count_lower_bound=neutral_birth,
            expected_fast_neutral_birth_energy_lower_bound_eV_per_source_ion=(
                average("neutral_birth_energy_eV")),
            resolved_fast_neutral_arrivals_per_source_ion=neutral_arrival,
            unresolved_fast_neutral_collisions_per_source_ion=(
                neutral_unresolved),
            escaped_fast_neutrals_per_source_ion=neutral_escape,
            fast_neutral_lineage_ledger_relative_residual=(
                neutral_lineage_residual),
            maximum_resolved_energy_ledger_relative_residual=(
                max_energy_residual),
            probability_ledger_relative_residual=probability_residual,
            collisionless_reference_mean_normal_energy_eV=float(
                np.mean(collisionless_energy)),
            mean_total_optical_depth=average("mean_initial_optical_depth"),
            maximum_total_optical_depth=max(
                phase.mean_initial_optical_depth
                for phase in phase_solutions),
            below_born_mayer_support_collision_probability_lower_bound=(
                average("below_support_probability")),
            model_source=self.collision_model.source,
            provenance={
                **dict(self.collision_model.provenance),
                **dict(self.provenance),
                "solver": "implicit_discrete_ordinates_absorbing_system",
                "phase_count": int(self.phase_count),
                "potential_node_count": int(self.potential_node_count),
                "total_energy_node_count": int(
                    self.total_energy_node_count),
                "transverse_fraction_node_count": int(
                    self.transverse_fraction_node_count),
                "position_quadrature_order": int(
                    self.position_quadrature_order),
                "hazard_quadrature_order": int(
                    self.hazard_quadrature_order),
                "impact_quadrature_order": int(
                    self.impact_quadrature_order),
                "collision_azimuth_order": int(
                    self.collision_azimuth_order),
                "ion_collision_order_closed": True,
                "ion_backscatter_turning_resolved": True,
                "linear_solve_relative_residual": max_linear_residual,
                "tangent_linear_solve_relative_residual": (
                    max_tangent_residual),
                "maximum_row_probability_residual": max_row_residual,
                "fast_neutral_transport_closed": False,
                "fast_neutral_no_further_collision_branch_resolved": True,
                "moving_sheath_self_consistency_closed": False,
                "feature_depth_used": False,
            },
        )
        tangent = None
        if gas_number_density_tangent_m3 is not None:
            tangent = CollisionalSheathDensityTangent(
                gas_number_density_tangent_m3=direction,
                distribution_weight_tangent=normalized_tangent,
                ion_arrival_probability_tangent=arrived_tangent,
                unresolved_probability_tangent=unresolved_tangent,
                escaped_probability_tangent=escaped_tangent,
                uncollided_arrival_probability_tangent=tangent_average(
                    "uncollided_arrival_probability_density_derivative"),
                expected_collision_count_tangent=tangent_average(
                    "collision_count_density_derivative"),
                expected_charge_exchange_count_tangent=tangent_average(
                    "charge_exchange_count_density_derivative"),
                expected_fast_neutral_birth_count_tangent=tangent_average(
                    "neutral_birth_count_density_derivative"),
                expected_fast_neutral_birth_energy_tangent_eV=tangent_average(
                    "neutral_birth_energy_density_derivative_eV"),
                resolved_fast_neutral_arrivals_tangent=tangent_average(
                    "neutral_arrival_density_derivative"),
                unresolved_fast_neutral_collisions_tangent=tangent_average(
                    "neutral_unresolved_density_derivative"),
                escaped_fast_neutrals_tangent=tangent_average(
                    "neutral_escape_density_derivative"),
                fast_neutral_lineage_ledger_tangent_residual=(
                    neutral_lineage_tangent_residual),
                mean_impact_energy_tangent_eV=float(np.dot(
                    normalized_tangent, distribution.energy_eV)),
                probability_ledger_tangent_residual=(
                    arrived_tangent + unresolved_tangent + escaped_tangent),
            )
        return solution, tangent

    def solve(self) -> DeterministicCollisionalSheathSolution:
        solution, _ = self._solve(gas_number_density_tangent_m3=None)
        return solution

    def density_jvp(
        self,
        gas_number_density_tangent_m3: float,
    ) -> tuple[
        DeterministicCollisionalSheathSolution,
        CollisionalSheathDensityTangent,
    ]:
        """Exact density JVP of the implicit discrete transport operator."""
        solution, tangent = self._solve(
            gas_number_density_tangent_m3=gas_number_density_tangent_m3)
        assert tangent is not None
        return solution, tangent


@dataclass(frozen=True)
class DiscreteOrdinatesConvergenceReceipt:
    coarse: DeterministicCollisionalSheathSolution
    fine: DeterministicCollisionalSheathSolution
    mean_energy_relative_change: float
    rms_angle_relative_change: float
    collision_count_relative_change: float
    neutral_arrival_relative_change: float
    ion_probability_ledger_residual: float
    ion_transport_closed: bool
    passed: bool
    limits: Mapping[str, float]

    def __post_init__(self):
        values = np.asarray([
            self.mean_energy_relative_change,
            self.rms_angle_relative_change,
            self.collision_count_relative_change,
            self.neutral_arrival_relative_change,
            self.ion_probability_ledger_residual,
        ], dtype=float)
        if (
            not isinstance(self.coarse, DeterministicCollisionalSheathSolution)
            or not isinstance(
                self.fine, DeterministicCollisionalSheathSolution)
            or np.any(~np.isfinite(values))
            or np.any(values < 0.0)
            or self.ion_transport_closed is not True
            or self.passed != all((
                self.mean_energy_relative_change
                <= self.limits["mean_energy_relative_change"],
                self.rms_angle_relative_change
                <= self.limits["rms_angle_relative_change"],
                self.collision_count_relative_change
                <= self.limits["collision_count_relative_change"],
                self.neutral_arrival_relative_change
                <= self.limits["neutral_arrival_relative_change"],
                self.ion_probability_ledger_residual
                <= self.limits["ion_probability_ledger_residual"],
            ))
        ):
            raise ValueError("invalid discrete-ordinates convergence receipt")
        object.__setattr__(
            self, "limits", MappingProxyType(dict(self.limits)))


def certify_discrete_ordinates_convergence(
    coarse: DeterministicCollisionalSheathSolution,
    fine: DeterministicCollisionalSheathSolution,
    *,
    mean_energy_relative_limit: float = 5.0e-3,
    rms_angle_relative_limit: float = 2.0e-2,
    collision_count_relative_limit: float = 1.0e-2,
    neutral_arrival_relative_limit: float = 1.0e-2,
    probability_ledger_limit: float = 2.0e-11,
) -> DiscreteOrdinatesConvergenceReceipt:
    """Grade two successively refined deterministic transport solutions."""
    if (
        coarse.provenance.get("ion_collision_order_closed") is not True
        or fine.provenance.get("ion_collision_order_closed") is not True
    ):
        raise ValueError("both solutions must close the ion collision order")

    def relative(coarse_value: float, fine_value: float) -> float:
        return float(abs(fine_value - coarse_value) / max(
            abs(fine_value), abs(coarse_value), np.finfo(float).tiny))

    coarse_angle = math.sqrt(
        coarse.distribution.mean_squared_polar_angle_rad2)
    fine_angle = math.sqrt(fine.distribution.mean_squared_polar_angle_rad2)
    limits = {
        "mean_energy_relative_change": float(mean_energy_relative_limit),
        "rms_angle_relative_change": float(rms_angle_relative_limit),
        "collision_count_relative_change": float(
            collision_count_relative_limit),
        "neutral_arrival_relative_change": float(
            neutral_arrival_relative_limit),
        "ion_probability_ledger_residual": float(probability_ledger_limit),
    }
    values = {
        "mean_energy_relative_change": relative(
            coarse.distribution.mean_energy_eV,
            fine.distribution.mean_energy_eV,
        ),
        "rms_angle_relative_change": relative(coarse_angle, fine_angle),
        "collision_count_relative_change": relative(
            coarse.expected_collision_count_lower_bound,
            fine.expected_collision_count_lower_bound,
        ),
        "neutral_arrival_relative_change": relative(
            coarse.resolved_fast_neutral_arrivals_per_source_ion,
            fine.resolved_fast_neutral_arrivals_per_source_ion,
        ),
        "ion_probability_ledger_residual": max(
            coarse.probability_ledger_relative_residual,
            fine.probability_ledger_relative_residual,
        ),
    }
    return DiscreteOrdinatesConvergenceReceipt(
        coarse=coarse,
        fine=fine,
        **values,
        ion_transport_closed=True,
        passed=all(values[name] <= limit for name, limit in limits.items()),
        limits=limits,
    )
