"""Quasineutral deterministic charged-species reactor transport.

This is the charged-particle companion to ``axisymmetric_reaction_diffusion``.
Its bulk/sheath split and quasineutral Boltzmann-electron closure follow
``wise-1996-rapid-2d-cl``.
For a fixed bulk potential it discretizes Nernst--Planck drift diffusion with
Scharfetter--Gummel face fluxes on the same axisymmetric finite-volume grid.
The exponential fitting is positivity preserving and resolves drift-dominated
faces without particle Monte Carlo.

``DeterministicQuasineutralInventoryLift`` closes the potential by the
Boltzmann-electron relation and quasineutrality.  Positive and negative ion
inventories remain declared 0-D constraints; their distributed source
amplitudes are inferred and exposed.  Consequently this tier predicts the
spatial wafer partition conditional on a global chemistry state, not the
global state itself.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import LinearOperator, gmres, splu
from scipy.optimize import root

from .axisymmetric_reaction_diffusion import (
    AXISYMMETRIC_WALL_LABELS,
    AxisymmetricFiniteVolumeGrid,
    _effective_robin_velocity,
)


def _bernoulli(value: float) -> float:
    """Stable Bernoulli function x/(exp(x)-1)."""
    x = float(value)
    if abs(x) < 1.0e-5:
        return 1.0 - 0.5 * x + x * x / 12.0 - x ** 4 / 720.0
    if x > 50.0:
        return x * math.exp(-x)
    if x < -50.0:
        return -x
    return x / math.expm1(x)


def _bernoulli_derivative(value: float) -> float:
    """Stable derivative of the Bernoulli function."""
    x = float(value)
    if abs(x) < 1.0e-4:
        return -0.5 + x / 6.0 - x ** 3 / 180.0 + x ** 5 / 5040.0
    if x > 50.0:
        return (1.0 - x) * math.exp(-x)
    if x < -50.0:
        return -1.0
    denominator = math.expm1(x)
    return (denominator - x * math.exp(x)) / denominator ** 2


@dataclass(frozen=True)
class AxisymmetricDriftDiffusionCondition:
    grid: AxisymmetricFiniteVolumeGrid
    species_names: tuple[str, ...]
    charge_number: np.ndarray
    mobility_m2_V_s: np.ndarray
    temperature_eV: np.ndarray
    electrostatic_potential_V: np.ndarray
    volume_reaction_matrix_s_inv: np.ndarray
    source_rate_m3_s: np.ndarray
    wall_velocity_m_s: np.ndarray
    source: str

    def __post_init__(self):
        if not isinstance(self.grid, AxisymmetricFiniteVolumeGrid):
            raise TypeError("an axisymmetric finite-volume grid is required")
        names = tuple(str(name) for name in self.species_names)
        count = len(names)
        nr = self.grid.radial_cell_count
        nz = self.grid.axial_cell_count
        charge = np.asarray(self.charge_number, dtype=float).copy()
        mobility = np.asarray(self.mobility_m2_V_s, dtype=float).copy()
        temperature = np.asarray(self.temperature_eV, dtype=float).copy()
        potential = np.asarray(
            self.electrostatic_potential_V, dtype=float).copy()
        reaction = np.asarray(
            self.volume_reaction_matrix_s_inv, dtype=float).copy()
        source_rate = np.asarray(self.source_rate_m3_s, dtype=float).copy()
        wall = np.asarray(self.wall_velocity_m_s, dtype=float).copy()
        off_diagonal = reaction.copy()
        if count:
            np.fill_diagonal(off_diagonal, 0.0)
        scale = max(1.0, float(np.max(np.abs(reaction))) if reaction.size else 1.0)
        if (
            count < 1
            or any(not name for name in names)
            or len(set(names)) != count
            or charge.shape != (count,)
            or mobility.shape != (count,)
            or temperature.shape != (count,)
            or potential.shape != (nr, nz)
            or reaction.shape != (count, count)
            or source_rate.shape != (count, nr, nz)
            or wall.shape != (count, len(AXISYMMETRIC_WALL_LABELS))
            or np.any(~np.isfinite(charge))
            or np.any(charge == 0.0)
            or np.any(~np.isfinite(mobility))
            or np.any(mobility <= 0.0)
            or np.any(~np.isfinite(temperature))
            or np.any(temperature <= 0.0)
            or np.any(~np.isfinite(potential))
            or np.any(~np.isfinite(reaction))
            or np.any(np.diag(reaction) < 0.0)
            or np.any(off_diagonal > 1.0e-14 * scale)
            or np.any(np.sum(reaction, axis=0) < -1.0e-12 * scale)
            or np.any(~np.isfinite(source_rate))
            or np.any(source_rate < 0.0)
            or np.any(np.isnan(wall))
            or np.any(wall < 0.0)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid axisymmetric drift-diffusion condition")
        for array in (
            charge, mobility, temperature, potential, reaction, source_rate, wall
        ):
            array.setflags(write=False)
        object.__setattr__(self, "species_names", names)
        object.__setattr__(self, "charge_number", charge)
        object.__setattr__(self, "mobility_m2_V_s", mobility)
        object.__setattr__(self, "temperature_eV", temperature)
        object.__setattr__(self, "electrostatic_potential_V", potential)
        object.__setattr__(self, "volume_reaction_matrix_s_inv", reaction)
        object.__setattr__(self, "source_rate_m3_s", source_rate)
        object.__setattr__(self, "wall_velocity_m_s", wall)

    @property
    def diffusion_coefficient_m2_s(self) -> np.ndarray:
        return self.mobility_m2_V_s * self.temperature_eV


@dataclass(frozen=True)
class AxisymmetricDriftDiffusionSolution:
    condition: AxisymmetricDriftDiffusionCondition
    density_m3: np.ndarray
    lower_endcap_flux_m2_s: np.ndarray
    upper_endcap_flux_m2_s: np.ndarray
    sidewall_flux_m2_s: np.ndarray
    integrated_source_rate_s: np.ndarray
    integrated_volume_reaction_rate_s: np.ndarray
    integrated_wall_loss_rate_s: np.ndarray
    maximum_species_ledger_relative_residual: float
    linear_system_relative_residual: float

    def __post_init__(self):
        if not isinstance(self.condition, AxisymmetricDriftDiffusionCondition):
            raise TypeError("drift-diffusion solution condition mismatch")
        count = len(self.condition.species_names)
        nr = self.condition.grid.radial_cell_count
        nz = self.condition.grid.axial_cell_count
        shapes = {
            "density_m3": (count, nr, nz),
            "lower_endcap_flux_m2_s": (count, nr),
            "upper_endcap_flux_m2_s": (count, nr),
            "sidewall_flux_m2_s": (count, nz),
            "integrated_source_rate_s": (count,),
            "integrated_volume_reaction_rate_s": (count,),
            "integrated_wall_loss_rate_s": (count,),
        }
        for name, shape in shapes.items():
            array = np.asarray(getattr(self, name), dtype=float).copy()
            if array.shape != shape or np.any(~np.isfinite(array)):
                raise ValueError("invalid drift-diffusion solution")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if (
            np.any(self.density_m3 < 0.0)
            or np.any(self.lower_endcap_flux_m2_s < 0.0)
            or np.any(self.upper_endcap_flux_m2_s < 0.0)
            or np.any(self.sidewall_flux_m2_s < 0.0)
            or not 0.0 <= self.maximum_species_ledger_relative_residual < 1.0e-8
            or not 0.0 <= self.linear_system_relative_residual < 1.0e-8
        ):
            raise ValueError("drift-diffusion conservation gate failed")

    def lower_endcap_area_average_flux_m2_s(
        self, species_name: str, *, wafer_radius_m: float
    ) -> float:
        try:
            index = self.condition.species_names.index(species_name)
        except ValueError as error:
            raise ValueError("unknown drift-diffusion species") from error
        radius = float(wafer_radius_m)
        grid = self.condition.grid
        if not math.isfinite(radius) or not 0.0 < radius <= grid.geometry.radius_m:
            raise ValueError("wafer radius lies outside the reactor grid")
        outer = np.minimum(grid.radial_edges_m[1:], radius)
        inner = np.minimum(grid.radial_edges_m[:-1], radius)
        area = np.pi * np.maximum(outer ** 2 - inner ** 2, 0.0)
        return float(
            np.dot(self.lower_endcap_flux_m2_s[index], area)
            / (np.pi * radius ** 2)
        )


class DeterministicAxisymmetricDriftDiffusion:
    """Fixed-potential Scharfetter--Gummel finite-volume solve."""

    def __init__(self, condition: AxisymmetricDriftDiffusionCondition):
        if not isinstance(condition, AxisymmetricDriftDiffusionCondition):
            raise TypeError("a drift-diffusion condition is required")
        self.condition = condition
        self._operator, self._boundary_velocity = self._assemble_operator()
        try:
            self._factorization = splu(self._operator)
        except RuntimeError as error:
            raise ValueError("drift-diffusion operator is singular") from error
        pivot = np.abs(self._factorization.U.diagonal())
        if pivot.size == 0 or float(np.min(pivot)) <= 1.0e-13 * max(
            float(np.max(pivot)), 1.0
        ):
            raise ValueError("drift-diffusion operator is singular")

    def _index(self, species: int, radial: int, axial: int) -> int:
        grid = self.condition.grid
        return (
            species * grid.radial_cell_count * grid.axial_cell_count
            + radial * grid.axial_cell_count + axial
        )

    def _add_face(
        self,
        operator: lil_matrix,
        *,
        species: int,
        left: tuple[int, int],
        right: tuple[int, int],
        area_m2: float,
        distance_m: float,
    ):
        condition = self.condition
        charge = condition.charge_number[species]
        temperature = condition.temperature_eV[species]
        potential = condition.electrostatic_potential_V
        psi = charge * (
            potential[right] - potential[left]
        ) / temperature
        conductance = (
            condition.diffusion_coefficient_m2_s[species]
            * area_m2 / distance_m
        )
        b_forward = _bernoulli(psi)
        b_reverse = _bernoulli(-psi)
        left_index = self._index(species, *left)
        right_index = self._index(species, *right)
        operator[left_index, left_index] += conductance * b_forward
        operator[left_index, right_index] -= conductance * b_reverse
        operator[right_index, right_index] += conductance * b_reverse
        operator[right_index, left_index] -= conductance * b_forward

    def _assemble_operator(self) -> tuple[csc_matrix, np.ndarray]:
        condition = self.condition
        grid = condition.grid
        count = len(condition.species_names)
        nr, nz = grid.radial_cell_count, grid.axial_cell_count
        size = count * nr * nz
        operator = lil_matrix((size, size), dtype=float)
        radial_center = grid.radial_centers_m
        axial_center = grid.axial_centers_m
        wall_effective = np.zeros((count, 3))
        for species in range(count):
            diffusion = condition.diffusion_coefficient_m2_s[species]
            wall_effective[species] = (
                _effective_robin_velocity(
                    diffusion, condition.wall_velocity_m_s[species, 0],
                    axial_center[0] - grid.axial_edges_m[0]),
                _effective_robin_velocity(
                    diffusion, condition.wall_velocity_m_s[species, 1],
                    grid.axial_edges_m[-1] - axial_center[-1]),
                _effective_robin_velocity(
                    diffusion, condition.wall_velocity_m_s[species, 2],
                    grid.radial_edges_m[-1] - radial_center[-1]),
            )
            for radial in range(nr):
                for axial in range(nz):
                    row = self._index(species, radial, axial)
                    volume = grid.cell_volume_m3[radial, axial]
                    for coupled in range(count):
                        coefficient = (
                            condition.volume_reaction_matrix_s_inv[
                                species, coupled] * volume
                        )
                        if coefficient != 0.0:
                            operator[row, self._index(
                                coupled, radial, axial)] += coefficient
                    if axial == 0:
                        operator[row, row] += (
                            wall_effective[species, 0]
                            * grid.axial_face_area_m2[radial]
                        )
                    if axial == nz - 1:
                        operator[row, row] += (
                            wall_effective[species, 1]
                            * grid.axial_face_area_m2[radial]
                        )
                    if radial == nr - 1:
                        operator[row, row] += (
                            wall_effective[species, 2]
                            * grid.outer_radial_face_area_m2[axial]
                        )
            for radial in range(nr - 1):
                face_radius = grid.radial_edges_m[radial + 1]
                for axial in range(nz):
                    area = 2.0 * np.pi * face_radius * (
                        grid.axial_edges_m[axial + 1]
                        - grid.axial_edges_m[axial]
                    )
                    self._add_face(
                        operator,
                        species=species,
                        left=(radial, axial),
                        right=(radial + 1, axial),
                        area_m2=area,
                        distance_m=radial_center[radial + 1] - radial_center[radial],
                    )
            for radial in range(nr):
                for axial in range(nz - 1):
                    self._add_face(
                        operator,
                        species=species,
                        left=(radial, axial),
                        right=(radial, axial + 1),
                        area_m2=grid.axial_face_area_m2[radial],
                        distance_m=axial_center[axial + 1] - axial_center[axial],
                    )
        return operator.tocsc(), wall_effective

    def _source_vector(self, source_rate_m3_s: np.ndarray) -> np.ndarray:
        return np.asarray(
            source_rate_m3_s
            * self.condition.grid.cell_volume_m3[None, :, :]
        ).ravel()

    def source_jvp(self, source_rate_tangent_m3_s: np.ndarray) -> np.ndarray:
        tangent = np.asarray(source_rate_tangent_m3_s, dtype=float)
        if (
            tangent.shape != self.condition.source_rate_m3_s.shape
            or np.any(~np.isfinite(tangent))
        ):
            raise ValueError("source tangent shape or values are invalid")
        return self._factorization.solve(
            self._source_vector(tangent)).reshape(tangent.shape)

    def potential_operator_jvp(
        self,
        density_m3: np.ndarray,
        electrostatic_potential_tangent_V: np.ndarray,
    ) -> np.ndarray:
        """Return ``dA/dphi[dphi] @ density`` for the SG operator.

        This is the exact fixed-topology directional derivative of every
        interior Scharfetter--Gummel face.  Wall and volume-reaction terms do
        not depend on the bulk potential in this closure.
        """
        condition = self.condition
        grid = condition.grid
        density = np.asarray(density_m3, dtype=float)
        tangent = np.asarray(electrostatic_potential_tangent_V, dtype=float)
        if (
            density.shape != condition.source_rate_m3_s.shape
            or tangent.shape != condition.electrostatic_potential_V.shape
            or np.any(~np.isfinite(density))
            or np.any(~np.isfinite(tangent))
        ):
            raise ValueError("invalid drift-diffusion potential JVP request")
        result = np.zeros_like(density)

        def add_face(
            species: int,
            left: tuple[int, int],
            right: tuple[int, int],
            area_m2: float,
            distance_m: float,
        ):
            charge = condition.charge_number[species]
            temperature = condition.temperature_eV[species]
            psi = charge * (
                condition.electrostatic_potential_V[right]
                - condition.electrostatic_potential_V[left]
            ) / temperature
            psi_tangent = charge * (
                tangent[right] - tangent[left]
            ) / temperature
            conductance = (
                condition.diffusion_coefficient_m2_s[species]
                * area_m2 / distance_m
            )
            flux_tangent = conductance * psi_tangent * (
                _bernoulli_derivative(psi) * density[(species, *left)]
                + _bernoulli_derivative(-psi) * density[(species, *right)]
            )
            result[(species, *left)] += flux_tangent
            result[(species, *right)] -= flux_tangent

        nr, nz = grid.radial_cell_count, grid.axial_cell_count
        for species in range(len(condition.species_names)):
            for radial in range(nr - 1):
                for axial in range(nz):
                    area = (
                        2.0 * np.pi * grid.radial_edges_m[radial + 1]
                        * (
                            grid.axial_edges_m[axial + 1]
                            - grid.axial_edges_m[axial]
                        )
                    )
                    add_face(
                        species,
                        (radial, axial),
                        (radial + 1, axial),
                        area,
                        grid.radial_centers_m[radial + 1]
                        - grid.radial_centers_m[radial],
                    )
            for radial in range(nr):
                for axial in range(nz - 1):
                    add_face(
                        species,
                        (radial, axial),
                        (radial, axial + 1),
                        grid.axial_face_area_m2[radial],
                        grid.axial_centers_m[axial + 1]
                        - grid.axial_centers_m[axial],
                    )
        return result

    def fixed_source_potential_jvp(
        self,
        density_m3: np.ndarray,
        electrostatic_potential_tangent_V: np.ndarray,
    ) -> np.ndarray:
        """Differentiate the fixed-source density through the sparse solve."""
        operator_tangent_density = self.potential_operator_jvp(
            density_m3, electrostatic_potential_tangent_V)
        return -self._factorization.solve(
            operator_tangent_density.ravel()
        ).reshape(operator_tangent_density.shape)

    def fixed_source_parameter_jvp(
        self,
        density_m3: np.ndarray,
        *,
        electrostatic_potential_tangent_V: np.ndarray,
        volume_reaction_matrix_tangent_s_inv: np.ndarray,
    ) -> np.ndarray:
        """Differentiate density for fixed source through potential/reaction.

        The caller supplies the directional derivative of the local linear
        reaction matrix.  This keeps the sparse solve generic while allowing
        nonlinear chemistry closures to use exact implicit derivatives.
        """
        density = np.asarray(density_m3, dtype=float)
        reaction_tangent = np.asarray(
            volume_reaction_matrix_tangent_s_inv, dtype=float)
        count = len(self.condition.species_names)
        if (
            density.shape != self.condition.source_rate_m3_s.shape
            or reaction_tangent.shape != (count, count)
            or np.any(~np.isfinite(density))
            or np.any(~np.isfinite(reaction_tangent))
        ):
            raise ValueError("invalid drift-diffusion parameter JVP request")
        operator_tangent_density = self.potential_operator_jvp(
            density, electrostatic_potential_tangent_V)
        operator_tangent_density += np.einsum(
            "st,tij,ij->sij",
            reaction_tangent,
            density,
            self.condition.grid.cell_volume_m3,
        )
        return -self._factorization.solve(
            operator_tangent_density.ravel()
        ).reshape(operator_tangent_density.shape)

    def solve(self) -> AxisymmetricDriftDiffusionSolution:
        rhs = self._source_vector(self.condition.source_rate_m3_s)
        raw = self._factorization.solve(rhs)
        # Two deterministic refinement steps materially improve the ledger for
        # drift-dominated matrices whose Scharfetter--Gummel coefficients span
        # many decades.
        for _ in range(2):
            raw += self._factorization.solve(rhs - self._operator @ raw)
        scale = max(1.0, float(np.max(np.abs(raw))))
        if float(np.min(raw)) < -1.0e-10 * scale:
            raise RuntimeError("drift-diffusion solve lost positivity")
        raw[np.abs(raw) < 1.0e-13 * scale] = 0.0
        density = np.maximum(raw, 0.0).reshape(
            self.condition.source_rate_m3_s.shape)
        residual = float(
            np.linalg.norm(self._operator @ density.ravel() - rhs)
            / max(np.linalg.norm(rhs), 1.0)
        )
        grid = self.condition.grid
        lower = self._boundary_velocity[:, 0, None] * density[:, :, 0]
        upper = self._boundary_velocity[:, 1, None] * density[:, :, -1]
        side = self._boundary_velocity[:, 2, None] * density[:, -1, :]
        source = np.sum(
            self.condition.source_rate_m3_s
            * grid.cell_volume_m3[None, :, :], axis=(1, 2))
        reaction = np.einsum(
            "st,tij,ij->s",
            self.condition.volume_reaction_matrix_s_inv,
            density,
            grid.cell_volume_m3,
        )
        wall = (
            np.einsum("si,i->s", lower, grid.axial_face_area_m2)
            + np.einsum("si,i->s", upper, grid.axial_face_area_m2)
            + np.einsum("sj,j->s", side, grid.outer_radial_face_area_m2)
        )
        ledger = np.abs(source - reaction - wall)
        ledger_scale = np.maximum.reduce((
            np.abs(source), np.abs(reaction) + np.abs(wall),
            np.ones_like(source)))
        return AxisymmetricDriftDiffusionSolution(
            condition=self.condition,
            density_m3=density,
            lower_endcap_flux_m2_s=lower,
            upper_endcap_flux_m2_s=upper,
            sidewall_flux_m2_s=side,
            integrated_source_rate_s=source,
            integrated_volume_reaction_rate_s=reaction,
            integrated_wall_loss_rate_s=wall,
            maximum_species_ledger_relative_residual=float(
                np.max(ledger / ledger_scale)),
            linear_system_relative_residual=residual,
        )


@dataclass(frozen=True)
class QuasineutralInventoryLiftResult:
    solution: AxisymmetricDriftDiffusionSolution
    electron_density_m3: np.ndarray
    electrostatic_potential_V: np.ndarray
    target_volume_average_density_m3: np.ndarray
    recovered_volume_average_density_m3: np.ndarray
    inferred_source_amplitude_m3_s: np.ndarray
    nonlinear_solver_evaluations: int
    maximum_potential_fixed_point_relative_residual: float
    maximum_inventory_relative_residual: float
    minimum_electron_density_m3: float
    supports_reactor_state_prediction: bool = False
    supports_implicit_differentiation: bool = True

    def __post_init__(self):
        if not isinstance(self.solution, AxisymmetricDriftDiffusionSolution):
            raise TypeError("quasineutral lift requires a spatial solution")
        count = len(self.solution.condition.species_names)
        grid_shape = (
            self.solution.condition.grid.radial_cell_count,
            self.solution.condition.grid.axial_cell_count,
        )
        specifications = {
            "electron_density_m3": grid_shape,
            "electrostatic_potential_V": grid_shape,
            "target_volume_average_density_m3": (count,),
            "recovered_volume_average_density_m3": (count,),
            "inferred_source_amplitude_m3_s": (count,),
        }
        for name, shape in specifications.items():
            array = np.asarray(getattr(self, name), dtype=float).copy()
            if array.shape != shape or np.any(~np.isfinite(array)):
                raise ValueError("invalid quasineutral inventory lift")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if (
            np.any(self.electron_density_m3 <= 0.0)
            or np.any(self.inferred_source_amplitude_m3_s < 0.0)
            or int(self.nonlinear_solver_evaluations) < 1
            or not 0.0 <= self.maximum_potential_fixed_point_relative_residual < 1.0e-7
            or not 0.0 <= self.maximum_inventory_relative_residual < 1.0e-10
            or self.minimum_electron_density_m3 <= 0.0
            or bool(self.supports_reactor_state_prediction)
            or not bool(self.supports_implicit_differentiation)
        ):
            raise ValueError("quasineutral inventory-lift certification failed")


@dataclass(frozen=True)
class QuasineutralInventoryLiftTangent:
    """Implicit derivative of a converged quasineutral inventory lift."""

    target_volume_average_density_tangent_m3: np.ndarray
    density_tangent_m3: np.ndarray
    electron_density_tangent_m3: np.ndarray
    electrostatic_potential_tangent_V: np.ndarray
    inferred_source_amplitude_tangent_m3_s: np.ndarray
    maximum_linearized_fixed_point_relative_residual: float
    maximum_inventory_tangent_relative_residual: float

    def __post_init__(self):
        target = np.asarray(
            self.target_volume_average_density_tangent_m3, dtype=float).copy()
        density = np.asarray(self.density_tangent_m3, dtype=float).copy()
        electron = np.asarray(self.electron_density_tangent_m3, dtype=float).copy()
        potential = np.asarray(
            self.electrostatic_potential_tangent_V, dtype=float).copy()
        source = np.asarray(
            self.inferred_source_amplitude_tangent_m3_s, dtype=float).copy()
        if (
            density.ndim != 3
            or target.shape != (density.shape[0],)
            or source.shape != target.shape
            or electron.shape != density.shape[1:]
            or potential.shape != density.shape[1:]
            or any(np.any(~np.isfinite(array)) for array in (
                target, density, electron, potential, source
            ))
            or not 0.0 <= self.maximum_linearized_fixed_point_relative_residual < 1.0e-7
            or not 0.0 <= self.maximum_inventory_tangent_relative_residual < 1.0e-8
        ):
            raise ValueError("invalid quasineutral inventory-lift tangent")
        for array in (target, density, electron, potential, source):
            array.setflags(write=False)
        object.__setattr__(
            self, "target_volume_average_density_tangent_m3", target)
        object.__setattr__(self, "density_tangent_m3", density)
        object.__setattr__(self, "electron_density_tangent_m3", electron)
        object.__setattr__(
            self, "electrostatic_potential_tangent_V", potential)
        object.__setattr__(
            self, "inferred_source_amplitude_tangent_m3_s", source)


class DeterministicQuasineutralInventoryLift:
    """Gummel-iterated Boltzmann-electron/quasineutral spatial lift."""

    def __init__(
        self,
        *,
        grid: AxisymmetricFiniteVolumeGrid,
        species_names: tuple[str, ...],
        charge_number: np.ndarray,
        mobility_m2_V_s: np.ndarray,
        ion_temperature_eV: np.ndarray,
        electron_temperature_eV: float,
        wall_velocity_m_s: np.ndarray,
        source_shape: np.ndarray,
        source: str,
        positive_negative_recombination_m3_s: float = 0.0,
    ):
        names = tuple(species_names)
        charge = np.asarray(charge_number, dtype=float).copy()
        if (
            charge.shape != (len(names),)
            or np.count_nonzero(charge > 0.0) < 1
            or np.count_nonzero(charge < 0.0) < 1
            or np.any(np.abs(charge) != 1.0)
        ):
            raise ValueError(
                "quasineutral lift requires singly charged positive and negative ions")
        shape = np.asarray(source_shape, dtype=float).copy()
        expected = (len(names), grid.radial_cell_count, grid.axial_cell_count)
        if shape.shape != expected or np.any(~np.isfinite(shape)) or np.any(shape < 0.0):
            raise ValueError("invalid quasineutral source moments")
        volume = grid.cell_volume_m3
        means = np.sum(shape * volume[None, :, :], axis=(1, 2)) / grid.geometry.volume_m3
        if np.any(np.abs(means - 1.0) > 1.0e-10):
            raise ValueError("quasineutral source moments must average to one")
        electron_temperature = float(electron_temperature_eV)
        if not math.isfinite(electron_temperature) or electron_temperature <= 0.0:
            raise ValueError("invalid electron temperature")
        self.grid = grid
        self.species_names = names
        self.charge_number = charge
        self.mobility_m2_V_s = np.asarray(mobility_m2_V_s, dtype=float).copy()
        self.ion_temperature_eV = np.asarray(ion_temperature_eV, dtype=float).copy()
        self.electron_temperature_eV = electron_temperature
        self.wall_velocity_m_s = np.asarray(wall_velocity_m_s, dtype=float).copy()
        self.source_shape = shape
        self.source = str(source).strip()
        recombination = float(positive_negative_recombination_m3_s)
        if not math.isfinite(recombination) or recombination < 0.0:
            raise ValueError("invalid positive-negative recombination coefficient")
        self.positive_negative_recombination_m3_s = recombination
        # Run the shared condition validator once at zero potential/source.
        AxisymmetricDriftDiffusionCondition(
            grid=grid,
            species_names=names,
            charge_number=charge,
            mobility_m2_V_s=self.mobility_m2_V_s,
            temperature_eV=self.ion_temperature_eV,
            electrostatic_potential_V=np.zeros((grid.radial_cell_count, grid.axial_cell_count)),
            volume_reaction_matrix_s_inv=np.zeros((len(names), len(names))),
            source_rate_m3_s=np.zeros(expected),
            wall_velocity_m_s=self.wall_velocity_m_s,
            source=self.source,
        )

    def _fixed_potential_lift(
        self,
        potential_V: np.ndarray,
        target_density_m3: np.ndarray,
    ) -> tuple[AxisymmetricDriftDiffusionSolution, np.ndarray]:
        count = len(self.species_names)
        reaction = self._reaction_matrix(target_density_m3)
        zero = np.zeros_like(self.source_shape)
        condition = AxisymmetricDriftDiffusionCondition(
            grid=self.grid,
            species_names=self.species_names,
            charge_number=self.charge_number,
            mobility_m2_V_s=self.mobility_m2_V_s,
            temperature_eV=self.ion_temperature_eV,
            electrostatic_potential_V=potential_V,
            volume_reaction_matrix_s_inv=reaction,
            source_rate_m3_s=zero,
            wall_velocity_m_s=self.wall_velocity_m_s,
            source=self.source,
        )
        transport = DeterministicAxisymmetricDriftDiffusion(condition)
        volume = self.grid.cell_volume_m3
        unit_average = np.empty(count)
        for species in range(count):
            tangent = np.zeros_like(self.source_shape)
            tangent[species] = self.source_shape[species]
            density = transport.source_jvp(tangent)
            unit_average[species] = (
                np.sum(density[species] * volume) / self.grid.geometry.volume_m3
            )
        amplitude = target_density_m3 / unit_average
        source_rate = amplitude[:, None, None] * self.source_shape
        solved_condition = AxisymmetricDriftDiffusionCondition(
            grid=self.grid,
            species_names=self.species_names,
            charge_number=self.charge_number,
            mobility_m2_V_s=self.mobility_m2_V_s,
            temperature_eV=self.ion_temperature_eV,
            electrostatic_potential_V=potential_V,
            volume_reaction_matrix_s_inv=reaction,
            source_rate_m3_s=source_rate,
            wall_velocity_m_s=self.wall_velocity_m_s,
            source=self.source + "; source amplitudes inferred from 0-D inventory",
        )
        return DeterministicAxisymmetricDriftDiffusion(solved_condition).solve(), amplitude

    def _reaction_matrix(self, target_density_m3: np.ndarray) -> np.ndarray:
        count = len(self.species_names)
        positive_density = float(np.sum(
            target_density_m3[self.charge_number > 0.0]))
        negative_density = float(np.sum(
            target_density_m3[self.charge_number < 0.0]))
        reaction = np.zeros((count, count))
        reaction[np.diag_indices(count)] = np.where(
            self.charge_number > 0.0,
            self.positive_negative_recombination_m3_s * negative_density,
            self.positive_negative_recombination_m3_s * positive_density,
        )
        return reaction

    def _reaction_matrix_jvp(
        self, target_density_tangent_m3: np.ndarray
    ) -> np.ndarray:
        count = len(self.species_names)
        positive_tangent = float(np.sum(
            target_density_tangent_m3[self.charge_number > 0.0]))
        negative_tangent = float(np.sum(
            target_density_tangent_m3[self.charge_number < 0.0]))
        tangent = np.zeros((count, count))
        tangent[np.diag_indices(count)] = np.where(
            self.charge_number > 0.0,
            self.positive_negative_recombination_m3_s * negative_tangent,
            self.positive_negative_recombination_m3_s * positive_tangent,
        )
        return tangent

    def solve(
        self,
        target_volume_average_density_m3: np.ndarray,
        *,
        relative_tolerance: float = 2.0e-8,
        maximum_iterations: int = 400,
        relaxation_fraction: float = 0.08,
        initial_electrostatic_potential_V: np.ndarray | None = None,
    ) -> QuasineutralInventoryLiftResult:
        target = np.asarray(target_volume_average_density_m3, dtype=float).copy()
        count = len(self.species_names)
        if (
            target.shape != (count,)
            or np.any(~np.isfinite(target))
            or np.any(target <= 0.0)
            or np.dot(self.charge_number, target) <= 0.0
            or not 0.0 < relative_tolerance < 1.0e-4
            or int(maximum_iterations) < 1
            or not 0.0 < relaxation_fraction <= 1.0
        ):
            raise ValueError("invalid quasineutral inventory-lift request")
        # ``relaxation_fraction`` is retained as an API-compatibility guard for
        # the superseded scalar Gummel iteration.  The deterministic Newton-
        # Krylov solver selects its steps from fixed residual/JVP probes.
        del relaxation_fraction
        grid_shape = (
            self.grid.radial_cell_count, self.grid.axial_cell_count)
        if initial_electrostatic_potential_V is None:
            initial_potential = np.zeros(grid_shape)
        else:
            initial_potential = np.asarray(
                initial_electrostatic_potential_V, dtype=float).copy()
            if (
                initial_potential.shape != grid_shape
                or np.any(~np.isfinite(initial_potential))
            ):
                raise ValueError("invalid initial electrostatic potential")
        volume = self.grid.cell_volume_m3
        volume_total = self.grid.geometry.volume_m3
        evaluation_count = 0

        def fixed_point_residual(flat_potential: np.ndarray) -> np.ndarray:
            nonlocal evaluation_count
            evaluation_count += 1
            potential = np.asarray(flat_potential, dtype=float).reshape(
                grid_shape)
            solution, _ = self._fixed_potential_lift(potential, target)
            electron_density = np.tensordot(
                self.charge_number, solution.density_m3, axes=(0, 0))
            electron_scale = max(
                float(np.sum(electron_density * volume) / volume_total), 1.0)
            # A zero-potential first Gummel iterate can locally put the more
            # strongly confined negative ion above total positive charge.
            # Continue from a strictly positive barrier field; acceptance
            # still requires the unmodified quasineutral density to be
            # positive everywhere.
            barrier_floor = 1.0e-8 * electron_scale
            positive_electron_density = np.maximum(
                electron_density, barrier_floor)
            candidate = self.electron_temperature_eV * np.log(
                positive_electron_density / electron_scale)
            candidate -= float(np.sum(candidate * volume) / volume_total)
            return ((candidate - potential)
                    / max(self.electron_temperature_eV, 1.0)).ravel()

        nonlinear = root(
            fixed_point_residual,
            initial_potential.ravel(),
            method="krylov",
            options={
                "maxiter": int(maximum_iterations),
                "fatol": float(relative_tolerance),
            },
        )
        fixed_point_error = float(np.max(np.abs(nonlinear.fun)))
        if (
            not nonlinear.success
            or not math.isfinite(fixed_point_error)
            or fixed_point_error >= relative_tolerance
        ):
            raise RuntimeError(
                "quasineutral Newton-Krylov solve failed to converge: "
                f"relative potential residual={fixed_point_error}; "
                f"message={nonlinear.message}"
            )
        potential = np.asarray(nonlinear.x, dtype=float).reshape(grid_shape)
        solution, amplitude = self._fixed_potential_lift(potential, target)
        electron_density = np.tensordot(
            self.charge_number, solution.density_m3, axes=(0, 0))
        if float(np.min(electron_density)) <= 0.0:
            raise RuntimeError(
                "converged quasineutral lift has nonpositive electron density")
        recovered = np.sum(
            solution.density_m3 * volume[None, :, :], axis=(1, 2)
        ) / volume_total
        inventory_residual = float(np.max(
            np.abs(recovered - target) / np.maximum(target, 1.0)))
        return QuasineutralInventoryLiftResult(
            solution=solution,
            electron_density_m3=electron_density,
            electrostatic_potential_V=potential,
            target_volume_average_density_m3=target,
            recovered_volume_average_density_m3=recovered,
            inferred_source_amplitude_m3_s=amplitude,
            nonlinear_solver_evaluations=evaluation_count,
            maximum_potential_fixed_point_relative_residual=(
                fixed_point_error),
            maximum_inventory_relative_residual=inventory_residual,
            minimum_electron_density_m3=float(np.min(electron_density)),
        )

    def implicit_target_inventory_jvp(
        self,
        result: QuasineutralInventoryLiftResult,
        target_volume_average_density_tangent_m3: np.ndarray,
        *,
        relative_tolerance: float = 2.0e-9,
        maximum_iterations: int = 400,
    ) -> QuasineutralInventoryLiftTangent:
        """Differentiate a converged lift by the implicit-function theorem.

        The Scharfetter--Gummel face derivative, target-dependent mutual
        neutralization, source-amplitude inventory constraint, Boltzmann
        electron closure, and potential gauge all enter the same matrix-free
        bordered solve.  No nonlinear iterations are unrolled and no finite
        differences enter the released derivative.
        """
        if not isinstance(result, QuasineutralInventoryLiftResult):
            raise TypeError("a converged quasineutral lift result is required")
        solution = result.solution
        if (
            solution.condition.grid is not self.grid
            or solution.condition.species_names != self.species_names
            or not np.array_equal(
                solution.condition.charge_number, self.charge_number)
        ):
            raise ValueError("quasineutral tangent/result model mismatch")
        target_tangent = np.asarray(
            target_volume_average_density_tangent_m3, dtype=float).copy()
        count = len(self.species_names)
        if (
            target_tangent.shape != (count,)
            or np.any(~np.isfinite(target_tangent))
            or not 0.0 < float(relative_tolerance) < 1.0e-5
            or int(maximum_iterations) < 1
        ):
            raise ValueError("invalid quasineutral inventory tangent request")
        transport = DeterministicAxisymmetricDriftDiffusion(
            solution.condition)
        density = solution.density_m3
        target = result.target_volume_average_density_m3
        volume = self.grid.cell_volume_m3
        volume_total = self.grid.geometry.volume_m3
        electron_density = result.electron_density_m3
        potential_scale = max(self.electron_temperature_eV, 1.0)
        grid_shape = electron_density.shape

        def density_response(
            potential_tangent_V: np.ndarray,
            inventory_tangent_m3: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray]:
            fixed_source = transport.fixed_source_parameter_jvp(
                density,
                electrostatic_potential_tangent_V=potential_tangent_V,
                volume_reaction_matrix_tangent_s_inv=(
                    self._reaction_matrix_jvp(inventory_tangent_m3)
                ),
            )
            fixed_average = np.sum(
                fixed_source * volume[None, :, :], axis=(1, 2)
            ) / volume_total
            amplitude_fraction_tangent = (
                inventory_tangent_m3 - fixed_average
            ) / target
            adjusted = (
                fixed_source
                + density * amplitude_fraction_tangent[:, None, None]
            )
            source_amplitude_tangent = (
                result.inferred_source_amplitude_m3_s
                * amplitude_fraction_tangent
            )
            return adjusted, source_amplitude_tangent

        def closure_tangent(
            potential_tangent_V: np.ndarray,
            inventory_tangent_m3: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            density_tangent, source_tangent = density_response(
                potential_tangent_V, inventory_tangent_m3)
            electron_tangent = np.tensordot(
                self.charge_number, density_tangent, axes=(0, 0))
            electron_log_tangent = electron_tangent / electron_density
            candidate_tangent = (
                self.electron_temperature_eV * electron_log_tangent)
            candidate_tangent -= float(
                np.sum(candidate_tangent * volume) / volume_total)
            fixed_point_tangent = (
                candidate_tangent - potential_tangent_V
            ) / potential_scale
            return (
                fixed_point_tangent,
                density_tangent,
                electron_tangent,
                source_tangent,
            )

        zeros_target = np.zeros(count)

        def matvec(flat_potential_tangent: np.ndarray) -> np.ndarray:
            potential_tangent = np.asarray(
                flat_potential_tangent, dtype=float).reshape(grid_shape)
            return closure_tangent(
                potential_tangent, zeros_target)[0].ravel()

        target_closure_tangent = closure_tangent(
            np.zeros(grid_shape), target_tangent)[0]
        operator = LinearOperator(
            (target_closure_tangent.size, target_closure_tangent.size),
            matvec=matvec,
            dtype=float,
        )
        flat_potential_tangent, info = gmres(
            operator,
            -target_closure_tangent.ravel(),
            rtol=float(relative_tolerance),
            atol=0.0,
            maxiter=int(maximum_iterations),
            restart=min(40, target_closure_tangent.size),
        )
        if info != 0 or np.any(~np.isfinite(flat_potential_tangent)):
            raise RuntimeError(
                "quasineutral implicit JVP failed to converge: "
                f"gmres_info={info}"
            )
        potential_tangent = flat_potential_tangent.reshape(grid_shape)
        (
            fixed_point_tangent,
            density_tangent,
            electron_tangent,
            source_tangent,
        ) = closure_tangent(potential_tangent, target_tangent)
        recovered_tangent = np.sum(
            density_tangent * volume[None, :, :], axis=(1, 2)
        ) / volume_total
        inventory_residual = float(np.max(
            np.abs(recovered_tangent - target_tangent)
            / max(float(np.max(np.abs(target_tangent))), 1.0)
        ))
        return QuasineutralInventoryLiftTangent(
            target_volume_average_density_tangent_m3=target_tangent,
            density_tangent_m3=density_tangent,
            electron_density_tangent_m3=electron_tangent,
            electrostatic_potential_tangent_V=potential_tangent,
            inferred_source_amplitude_tangent_m3_s=source_tangent,
            maximum_linearized_fixed_point_relative_residual=float(
                np.max(np.abs(fixed_point_tangent))
            ),
            maximum_inventory_tangent_relative_residual=inventory_residual,
        )
