"""Deterministic axisymmetric inductive-field and power-deposition tier.

For a harmonic azimuthal field ``E_theta exp(i omega t)`` in an axisymmetric
conducting plasma, the quasistatic bulk equation is

``-[(1/r)d_r(r d_r E_theta) + d_z^2 E_theta - E_theta/r^2]
  + i mu0 omega sigma E_theta = 0``.

The finite-volume operator below is conservative and complex symmetric.  A
declared coil-side tangential field is imposed on the upper plasma boundary;
the axis, lower endcap, and sidewall use homogeneous tangential-field
conditions.  Its discrete complex-power identity closes the time-average
Ohmic absorption to roundoff.  Exact conductivity and boundary-field JVPs
reuse the same sparse factorization.

This module predicts a deposited-power *shape* conditional on conductivity
and the coil-side field.  It does not infer generator-to-plasma coupling or a
coil field from source-power setpoint alone.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .axisymmetric_reaction_diffusion import AxisymmetricFiniteVolumeGrid
from .network import E_CHARGE_C


VACUUM_PERMEABILITY_H_M = 4.0e-7 * np.pi


@dataclass(frozen=True)
class AxisymmetricInductiveFieldCondition:
    grid: AxisymmetricFiniteVolumeGrid
    frequency_hz: float
    conductivity_S_m: np.ndarray
    upper_boundary_electric_field_V_m: np.ndarray
    source: str
    conductivity_source: str
    upper_boundary_source: str

    def __post_init__(self):
        if not isinstance(self.grid, AxisymmetricFiniteVolumeGrid):
            raise TypeError("an axisymmetric finite-volume grid is required")
        nr = self.grid.radial_cell_count
        nz = self.grid.axial_cell_count
        frequency = float(self.frequency_hz)
        conductivity = np.asarray(
            self.conductivity_S_m, dtype=complex).copy()
        upper = np.asarray(
            self.upper_boundary_electric_field_V_m, dtype=complex).copy()
        if (
            not math.isfinite(frequency)
            or frequency <= 0.0
            or conductivity.shape != (nr, nz)
            or upper.shape != (nr,)
            or np.any(~np.isfinite(conductivity.real))
            or np.any(~np.isfinite(conductivity.imag))
            or np.any(conductivity.real <= 0.0)
            or np.any(~np.isfinite(upper.real))
            or np.any(~np.isfinite(upper.imag))
            or float(np.max(np.abs(upper))) <= 0.0
            or not str(self.source).strip()
            or not str(self.conductivity_source).strip()
            or not str(self.upper_boundary_source).strip()
        ):
            raise ValueError("invalid axisymmetric inductive-field condition")
        conductivity.setflags(write=False)
        upper.setflags(write=False)
        object.__setattr__(self, "frequency_hz", frequency)
        object.__setattr__(self, "conductivity_S_m", conductivity)
        object.__setattr__(
            self, "upper_boundary_electric_field_V_m", upper)

    @property
    def angular_frequency_rad_s(self) -> float:
        return 2.0 * np.pi * self.frequency_hz


@dataclass(frozen=True)
class AxisymmetricInductiveFieldSolution:
    condition: AxisymmetricInductiveFieldCondition
    electric_field_V_m: np.ndarray
    absorbed_power_density_W_m3: np.ndarray
    normalized_absorbed_power_source_moment: np.ndarray
    total_absorbed_power_W: float
    boundary_complex_power_V2_m: complex
    complex_power_ledger_relative_residual: float
    linear_system_relative_residual: float
    supports_implicit_differentiation: bool = True
    supports_generator_to_plasma_power_prediction: bool = False

    def __post_init__(self):
        if not isinstance(self.condition, AxisymmetricInductiveFieldCondition):
            raise TypeError("inductive-field solution condition mismatch")
        shape = self.condition.conductivity_S_m.shape
        field = np.asarray(self.electric_field_V_m, dtype=complex).copy()
        power = np.asarray(
            self.absorbed_power_density_W_m3, dtype=float).copy()
        moment = np.asarray(
            self.normalized_absorbed_power_source_moment, dtype=float).copy()
        if (
            field.shape != shape
            or power.shape != shape
            or moment.shape != shape
            or np.any(~np.isfinite(field.real))
            or np.any(~np.isfinite(field.imag))
            or np.any(~np.isfinite(power))
            or np.any(~np.isfinite(moment))
            or np.any(power < 0.0)
            or np.any(moment < 0.0)
            or not math.isfinite(float(self.total_absorbed_power_W))
            or self.total_absorbed_power_W <= 0.0
            or not 0.0 <= self.complex_power_ledger_relative_residual < 1.0e-10
            or not 0.0 <= self.linear_system_relative_residual < 1.0e-10
            or not bool(self.supports_implicit_differentiation)
            or bool(self.supports_generator_to_plasma_power_prediction)
        ):
            raise ValueError("inductive-field conservation gate failed")
        volume = self.condition.grid.cell_volume_m3
        average = float(
            np.sum(moment * volume) / self.condition.grid.geometry.volume_m3)
        if not math.isclose(average, 1.0, rel_tol=2.0e-12, abs_tol=2.0e-12):
            raise ValueError("inductive power source moment is not normalized")
        for array in (field, power, moment):
            array.setflags(write=False)
        object.__setattr__(self, "electric_field_V_m", field)
        object.__setattr__(self, "absorbed_power_density_W_m3", power)
        object.__setattr__(
            self, "normalized_absorbed_power_source_moment", moment)
        object.__setattr__(
            self, "total_absorbed_power_W", float(self.total_absorbed_power_W))


@dataclass(frozen=True)
class AxisymmetricInductiveFieldTangent:
    conductivity_tangent_S_m: np.ndarray
    upper_boundary_electric_field_tangent_V_m: np.ndarray
    electric_field_tangent_V_m: np.ndarray
    absorbed_power_density_tangent_W_m3: np.ndarray
    normalized_source_moment_tangent: np.ndarray
    total_absorbed_power_tangent_W: float
    linearized_system_relative_residual: float

    def __post_init__(self):
        conductivity = np.asarray(
            self.conductivity_tangent_S_m, dtype=complex).copy()
        upper = np.asarray(
            self.upper_boundary_electric_field_tangent_V_m,
            dtype=complex,
        ).copy()
        field = np.asarray(
            self.electric_field_tangent_V_m, dtype=complex).copy()
        power = np.asarray(
            self.absorbed_power_density_tangent_W_m3, dtype=float).copy()
        moment = np.asarray(
            self.normalized_source_moment_tangent, dtype=float).copy()
        if (
            conductivity.shape != field.shape
            or power.shape != field.shape
            or moment.shape != field.shape
            or upper.shape != (field.shape[0],)
            or any(
                np.any(~np.isfinite(array.real))
                or np.any(~np.isfinite(array.imag))
                for array in (conductivity, upper, field)
            )
            or np.any(~np.isfinite(power))
            or np.any(~np.isfinite(moment))
            or not math.isfinite(float(self.total_absorbed_power_tangent_W))
            or not 0.0 <= self.linearized_system_relative_residual < 1.0e-9
        ):
            raise ValueError("invalid inductive-field tangent")
        for array in (conductivity, upper, field, power, moment):
            array.setflags(write=False)
        object.__setattr__(self, "conductivity_tangent_S_m", conductivity)
        object.__setattr__(
            self, "upper_boundary_electric_field_tangent_V_m", upper)
        object.__setattr__(self, "electric_field_tangent_V_m", field)
        object.__setattr__(
            self, "absorbed_power_density_tangent_W_m3", power)
        object.__setattr__(
            self, "normalized_source_moment_tangent", moment)


class DeterministicAxisymmetricInductiveField:
    """Sparse harmonic azimuthal-field solve with exact parameter JVPs."""

    def __init__(self, condition: AxisymmetricInductiveFieldCondition):
        if not isinstance(condition, AxisymmetricInductiveFieldCondition):
            raise TypeError("an inductive-field condition is required")
        self.condition = condition
        self._operator, self._boundary_map = self._assemble()
        try:
            self._factorization = splu(self._operator)
        except RuntimeError as error:
            raise ValueError("inductive-field operator is singular") from error

    def _index(self, radial: int, axial: int) -> int:
        return radial * self.condition.grid.axial_cell_count + axial

    def _assemble(self) -> tuple[csc_matrix, np.ndarray]:
        condition = self.condition
        grid = condition.grid
        nr, nz = grid.radial_cell_count, grid.axial_cell_count
        operator = lil_matrix((nr * nz, nr * nz), dtype=complex)
        upper_map = np.zeros(nr)
        radial_center = grid.radial_centers_m
        axial_center = grid.axial_centers_m
        volume = grid.cell_volume_m3
        for radial in range(nr):
            for axial in range(nz):
                row = self._index(radial, axial)
                operator[row, row] += (
                    volume[radial, axial] / radial_center[radial] ** 2
                    + 1j
                    * VACUUM_PERMEABILITY_H_M
                    * condition.angular_frequency_rad_s
                    * condition.conductivity_S_m[radial, axial]
                    * volume[radial, axial]
                )
                if radial > 0:
                    area = (
                        2.0 * np.pi * grid.radial_edges_m[radial]
                        * (grid.axial_edges_m[axial + 1]
                           - grid.axial_edges_m[axial])
                    )
                    conductance = area / (
                        radial_center[radial] - radial_center[radial - 1])
                    neighbor = self._index(radial - 1, axial)
                    operator[row, row] += conductance
                    operator[row, neighbor] -= conductance
                if radial < nr - 1:
                    area = (
                        2.0 * np.pi * grid.radial_edges_m[radial + 1]
                        * (grid.axial_edges_m[axial + 1]
                           - grid.axial_edges_m[axial])
                    )
                    conductance = area / (
                        radial_center[radial + 1] - radial_center[radial])
                    neighbor = self._index(radial + 1, axial)
                    operator[row, row] += conductance
                    operator[row, neighbor] -= conductance
                else:
                    conductance = (
                        grid.outer_radial_face_area_m2[axial]
                        / (grid.radial_edges_m[-1] - radial_center[-1])
                    )
                    operator[row, row] += conductance
                if axial > 0:
                    conductance = (
                        grid.axial_face_area_m2[radial]
                        / (axial_center[axial] - axial_center[axial - 1])
                    )
                    neighbor = self._index(radial, axial - 1)
                    operator[row, row] += conductance
                    operator[row, neighbor] -= conductance
                else:
                    conductance = (
                        grid.axial_face_area_m2[radial]
                        / (axial_center[0] - grid.axial_edges_m[0])
                    )
                    operator[row, row] += conductance
                if axial < nz - 1:
                    conductance = (
                        grid.axial_face_area_m2[radial]
                        / (axial_center[axial + 1] - axial_center[axial])
                    )
                    neighbor = self._index(radial, axial + 1)
                    operator[row, row] += conductance
                    operator[row, neighbor] -= conductance
                else:
                    conductance = (
                        grid.axial_face_area_m2[radial]
                        / (grid.axial_edges_m[-1] - axial_center[-1])
                    )
                    operator[row, row] += conductance
                    upper_map[radial] = conductance
        return operator.tocsc(), upper_map

    def _rhs(self, upper_boundary_field: np.ndarray) -> np.ndarray:
        rhs = np.zeros(self._operator.shape[0], dtype=complex)
        nz = self.condition.grid.axial_cell_count
        for radial, (conductance, value) in enumerate(zip(
            self._boundary_map, upper_boundary_field
        )):
            rhs[radial * nz + nz - 1] = conductance * value
        return rhs

    def solve(self) -> AxisymmetricInductiveFieldSolution:
        condition = self.condition
        rhs = self._rhs(condition.upper_boundary_electric_field_V_m)
        flat_field = self._factorization.solve(rhs)
        field = flat_field.reshape(condition.conductivity_S_m.shape)
        residual = float(
            np.linalg.norm(self._operator @ flat_field - rhs)
            / max(np.linalg.norm(rhs), 1.0)
        )
        power_density = (
            0.5 * condition.conductivity_S_m.real * np.abs(field) ** 2)
        volume = condition.grid.cell_volume_m3
        total_power = float(np.sum(power_density * volume))
        average_power = total_power / condition.grid.geometry.volume_m3
        moment = power_density / average_power
        boundary_complex_power = np.vdot(flat_field, rhs)
        ledger_power = float(
            boundary_complex_power.imag
            / (
                2.0
                * VACUUM_PERMEABILITY_H_M
                * condition.angular_frequency_rad_s
            )
        )
        ledger_residual = abs(ledger_power - total_power) / max(
            abs(ledger_power), total_power, 1.0e-300)
        return AxisymmetricInductiveFieldSolution(
            condition=condition,
            electric_field_V_m=field,
            absorbed_power_density_W_m3=power_density,
            normalized_absorbed_power_source_moment=moment,
            total_absorbed_power_W=total_power,
            boundary_complex_power_V2_m=complex(boundary_complex_power),
            complex_power_ledger_relative_residual=float(ledger_residual),
            linear_system_relative_residual=residual,
        )

    def parameter_jvp(
        self,
        solution: AxisymmetricInductiveFieldSolution,
        *,
        conductivity_tangent_S_m: np.ndarray,
        upper_boundary_electric_field_tangent_V_m: np.ndarray,
    ) -> AxisymmetricInductiveFieldTangent:
        if (
            not isinstance(solution, AxisymmetricInductiveFieldSolution)
            or solution.condition is not self.condition
        ):
            raise ValueError("inductive-field solution/model mismatch")
        condition = self.condition
        conductivity_tangent = np.asarray(
            conductivity_tangent_S_m, dtype=complex)
        upper_tangent = np.asarray(
            upper_boundary_electric_field_tangent_V_m, dtype=complex)
        if (
            conductivity_tangent.shape != condition.conductivity_S_m.shape
            or upper_tangent.shape
            != condition.upper_boundary_electric_field_V_m.shape
            or np.any(~np.isfinite(conductivity_tangent.real))
            or np.any(~np.isfinite(conductivity_tangent.imag))
            or np.any(~np.isfinite(upper_tangent.real))
            or np.any(~np.isfinite(upper_tangent.imag))
        ):
            raise ValueError("invalid inductive-field parameter tangent")
        volume = condition.grid.cell_volume_m3
        operator_tangent_field = (
            1j
            * VACUUM_PERMEABILITY_H_M
            * condition.angular_frequency_rad_s
            * conductivity_tangent
            * volume
            * solution.electric_field_V_m
        ).ravel()
        rhs_tangent = self._rhs(upper_tangent)
        field_tangent_flat = self._factorization.solve(
            rhs_tangent - operator_tangent_field)
        field_tangent = field_tangent_flat.reshape(
            condition.conductivity_S_m.shape)
        power_tangent = 0.5 * (
            conductivity_tangent.real
            * np.abs(solution.electric_field_V_m) ** 2
            + 2.0
            * condition.conductivity_S_m.real
            * np.real(
                np.conj(solution.electric_field_V_m) * field_tangent)
        )
        total_tangent = float(np.sum(power_tangent * volume))
        average = (
            solution.total_absorbed_power_W
            / condition.grid.geometry.volume_m3)
        average_tangent = total_tangent / condition.grid.geometry.volume_m3
        moment_tangent = (
            power_tangent / average
            - solution.normalized_absorbed_power_source_moment
            * average_tangent / average
        )
        linearized_residual = float(
            np.linalg.norm(
                self._operator @ field_tangent_flat
                + operator_tangent_field
                - rhs_tangent
            )
            / max(
                np.linalg.norm(rhs_tangent),
                np.linalg.norm(operator_tangent_field),
                1.0,
            )
        )
        return AxisymmetricInductiveFieldTangent(
            conductivity_tangent_S_m=conductivity_tangent,
            upper_boundary_electric_field_tangent_V_m=upper_tangent,
            electric_field_tangent_V_m=field_tangent,
            absorbed_power_density_tangent_W_m3=power_tangent,
            normalized_source_moment_tangent=moment_tangent,
            total_absorbed_power_tangent_W=total_tangent,
            linearized_system_relative_residual=linearized_residual,
        )


def drude_conductivity_from_two_term_mobilities(
    *,
    electron_density_m3: float,
    neutral_density_m3: float,
    frequency_hz: float,
    flux_reduced_mobility_m_inv_V_inv_s_inv: float,
    dissipative_reduced_mobility_m_inv_V_inv_s_inv: float,
) -> complex:
    """Recover a one-pole complex conductivity from EEPF mobility moments.

    The two-term solver supplies the DC/flux mobility and the RF-dissipative
    mobility.  Their ratio identifies the Drude collision-to-frequency ratio;
    the returned real part exactly reproduces the EEPF power-gain moment.
    """
    values = np.asarray((
        electron_density_m3,
        neutral_density_m3,
        frequency_hz,
        flux_reduced_mobility_m_inv_V_inv_s_inv,
        dissipative_reduced_mobility_m_inv_V_inv_s_inv,
    ), dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("invalid Drude conductivity inputs")
    electron_density, neutral_density, frequency, flux_reduced, dissipative = values
    ratio = dissipative / flux_reduced
    if not 0.0 < ratio <= 1.0 + 1.0e-10:
        raise ValueError("dissipative mobility exceeds DC flux mobility")
    ratio = min(float(ratio), 1.0)
    dc_mobility = flux_reduced / neutral_density
    real_conductivity = (
        E_CHARGE_C * electron_density * dissipative / neutral_density)
    if 1.0 - ratio < 1.0e-14:
        return complex(real_conductivity, 0.0)
    omega_over_collision_frequency = math.sqrt((1.0 - ratio) / ratio)
    return complex(
        real_conductivity,
        -real_conductivity * omega_over_collision_frequency,
    )
