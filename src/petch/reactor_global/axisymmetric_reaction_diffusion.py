"""Conservative deterministic axisymmetric reactor-transport tier.

This module lifts a well-mixed chemistry state into a spatial wafer boundary
without particle Monte Carlo.  Multiple species diffuse on a fixed
axisymmetric finite-volume grid, exchange through a local linear reaction
matrix, and leave through independent lower-endcap, upper-endcap, and
sidewall Robin velocities.  The sparse operator is an M-matrix, so the solve
is positivity preserving; every returned species carries an integrated
source/reaction/wall ledger.

The topology is fixed and the residual is linear.  Source JVPs therefore use
the same sparse factorization as the primal solve, providing an exact
deterministic differentiability contract for batch recipe calculations.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .geometry import CylindricalReactor


AXISYMMETRIC_WALL_LABELS = (
    "lower_endcap_wafer_plane",
    "upper_endcap",
    "cylindrical_sidewall",
)


@dataclass(frozen=True)
class AxisymmetricFiniteVolumeGrid:
    radial_edges_m: np.ndarray
    axial_edges_m: np.ndarray

    def __post_init__(self):
        radial = np.asarray(self.radial_edges_m, dtype=float).copy()
        axial = np.asarray(self.axial_edges_m, dtype=float).copy()
        if (
            radial.ndim != 1
            or axial.ndim != 1
            or min(radial.size, axial.size) < 3
            or radial[0] != 0.0
            or axial[0] != 0.0
            or np.any(~np.isfinite(radial))
            or np.any(~np.isfinite(axial))
            or np.any(np.diff(radial) <= 0.0)
            or np.any(np.diff(axial) <= 0.0)
        ):
            raise ValueError("invalid axisymmetric finite-volume grid")
        radial.setflags(write=False)
        axial.setflags(write=False)
        object.__setattr__(self, "radial_edges_m", radial)
        object.__setattr__(self, "axial_edges_m", axial)

    @classmethod
    def uniform(
        cls,
        geometry: CylindricalReactor,
        *,
        radial_cell_count: int,
        axial_cell_count: int,
    ) -> "AxisymmetricFiniteVolumeGrid":
        if not isinstance(geometry, CylindricalReactor):
            raise TypeError("a cylindrical reactor geometry is required")
        nr = int(radial_cell_count)
        nz = int(axial_cell_count)
        if (
            nr != radial_cell_count
            or nz != axial_cell_count
            or min(nr, nz) < 2
        ):
            raise ValueError("axisymmetric grid requires at least 2x2 cells")
        return cls(
            radial_edges_m=np.linspace(0.0, geometry.radius_m, nr + 1),
            axial_edges_m=np.linspace(0.0, geometry.length_m, nz + 1),
        )

    @property
    def geometry(self) -> CylindricalReactor:
        return CylindricalReactor(
            radius_m=float(self.radial_edges_m[-1]),
            length_m=float(self.axial_edges_m[-1]),
        )

    @property
    def radial_cell_count(self) -> int:
        return self.radial_edges_m.size - 1

    @property
    def axial_cell_count(self) -> int:
        return self.axial_edges_m.size - 1

    @property
    def radial_centers_m(self) -> np.ndarray:
        return 0.5 * (self.radial_edges_m[:-1] + self.radial_edges_m[1:])

    @property
    def axial_centers_m(self) -> np.ndarray:
        return 0.5 * (self.axial_edges_m[:-1] + self.axial_edges_m[1:])

    @property
    def cell_volume_m3(self) -> np.ndarray:
        annular_area = np.pi * (
            self.radial_edges_m[1:] ** 2
            - self.radial_edges_m[:-1] ** 2)
        return annular_area[:, None] * np.diff(self.axial_edges_m)[None, :]

    @property
    def axial_face_area_m2(self) -> np.ndarray:
        return np.pi * (
            self.radial_edges_m[1:] ** 2
            - self.radial_edges_m[:-1] ** 2)

    @property
    def outer_radial_face_area_m2(self) -> np.ndarray:
        return (
            2.0 * np.pi * self.geometry.radius_m
            * np.diff(self.axial_edges_m))


@dataclass(frozen=True)
class AxisymmetricReactionDiffusionCondition:
    grid: AxisymmetricFiniteVolumeGrid
    species_names: tuple[str, ...]
    diffusion_coefficient_m2_s: np.ndarray
    volume_reaction_matrix_s_inv: np.ndarray
    source_rate_m3_s: np.ndarray
    wall_velocity_m_s: np.ndarray
    source: str

    def __post_init__(self):
        if not isinstance(self.grid, AxisymmetricFiniteVolumeGrid):
            raise TypeError("an axisymmetric finite-volume grid is required")
        names = tuple(str(name) for name in self.species_names)
        diffusion = np.asarray(
            self.diffusion_coefficient_m2_s, dtype=float).copy()
        reaction = np.asarray(
            self.volume_reaction_matrix_s_inv, dtype=float).copy()
        source_rate = np.asarray(self.source_rate_m3_s, dtype=float).copy()
        wall = np.asarray(self.wall_velocity_m_s, dtype=float).copy()
        species_count = len(names)
        expected_field = (
            species_count,
            self.grid.radial_cell_count,
            self.grid.axial_cell_count,
        )
        off_diagonal = reaction.copy()
        if species_count:
            np.fill_diagonal(off_diagonal, 0.0)
        column_sum = np.sum(reaction, axis=0) if species_count else np.array([])
        scale = max(1.0, float(np.max(np.abs(reaction))) if reaction.size else 1.0)
        if (
            species_count < 1
            or any(not name for name in names)
            or len(set(names)) != species_count
            or diffusion.shape != (species_count,)
            or reaction.shape != (species_count, species_count)
            or source_rate.shape != expected_field
            or wall.shape != (species_count, len(AXISYMMETRIC_WALL_LABELS))
            or np.any(~np.isfinite(diffusion))
            or np.any(diffusion <= 0.0)
            or np.any(~np.isfinite(reaction))
            or np.any(np.diag(reaction) < 0.0)
            or np.any(off_diagonal > 1.0e-14 * scale)
            or np.any(column_sum < -1.0e-12 * scale)
            or np.any(~np.isfinite(source_rate))
            or np.any(source_rate < 0.0)
            or np.any(np.isnan(wall))
            or np.any(wall < 0.0)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid axisymmetric reaction-diffusion condition")
        for array in (diffusion, reaction, source_rate, wall):
            array.setflags(write=False)
        object.__setattr__(self, "species_names", names)
        object.__setattr__(self, "diffusion_coefficient_m2_s", diffusion)
        object.__setattr__(self, "volume_reaction_matrix_s_inv", reaction)
        object.__setattr__(self, "source_rate_m3_s", source_rate)
        object.__setattr__(self, "wall_velocity_m_s", wall)


@dataclass(frozen=True)
class AxisymmetricReactionDiffusionSolution:
    condition: AxisymmetricReactionDiffusionCondition
    density_m3: np.ndarray
    lower_endcap_flux_m2_s: np.ndarray
    upper_endcap_flux_m2_s: np.ndarray
    sidewall_flux_m2_s: np.ndarray
    integrated_source_rate_s: np.ndarray
    integrated_volume_reaction_rate_s: np.ndarray
    integrated_wall_loss_rate_s: np.ndarray
    maximum_species_ledger_relative_residual: float
    linear_system_relative_residual: float
    minimum_density_m3: float
    source_jvp_supported: bool = True

    def __post_init__(self):
        condition = self.condition
        if not isinstance(condition, AxisymmetricReactionDiffusionCondition):
            raise TypeError("reaction-diffusion solution condition mismatch")
        species_count = len(condition.species_names)
        nr = condition.grid.radial_cell_count
        nz = condition.grid.axial_cell_count
        expected = {
            "density_m3": (species_count, nr, nz),
            "lower_endcap_flux_m2_s": (species_count, nr),
            "upper_endcap_flux_m2_s": (species_count, nr),
            "sidewall_flux_m2_s": (species_count, nz),
            "integrated_source_rate_s": (species_count,),
            "integrated_volume_reaction_rate_s": (species_count,),
            "integrated_wall_loss_rate_s": (species_count,),
        }
        for name, shape in expected.items():
            array = np.asarray(getattr(self, name), dtype=float).copy()
            if array.shape != shape or np.any(~np.isfinite(array)):
                raise ValueError("invalid axisymmetric reaction-diffusion solution")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if (
            np.any(self.density_m3 < 0.0)
            or np.any(self.lower_endcap_flux_m2_s < 0.0)
            or np.any(self.upper_endcap_flux_m2_s < 0.0)
            or np.any(self.sidewall_flux_m2_s < 0.0)
            or not math.isfinite(self.maximum_species_ledger_relative_residual)
            or not 0.0 <= self.maximum_species_ledger_relative_residual < 1.0e-8
            or not math.isfinite(self.linear_system_relative_residual)
            or not 0.0 <= self.linear_system_relative_residual < 1.0e-10
            or not math.isfinite(self.minimum_density_m3)
            or self.minimum_density_m3 < 0.0
            or not bool(self.source_jvp_supported)
        ):
            raise ValueError("reaction-diffusion conservation gate failed")

    def lower_endcap_area_average_flux_m2_s(
        self,
        species_name: str,
        *,
        wafer_radius_m: float,
    ) -> float:
        """Average a piecewise-annular lower-endcap flux over a wafer disk."""
        try:
            species_index = self.condition.species_names.index(species_name)
        except ValueError as error:
            raise ValueError("unknown reaction-diffusion species") from error
        radius = float(wafer_radius_m)
        grid = self.condition.grid
        if (
            not math.isfinite(radius)
            or not 0.0 < radius <= grid.geometry.radius_m
        ):
            raise ValueError("wafer radius lies outside the reactor grid")
        clipped_outer = np.minimum(grid.radial_edges_m[1:], radius)
        clipped_inner = np.minimum(grid.radial_edges_m[:-1], radius)
        area = np.pi * np.maximum(
            clipped_outer ** 2 - clipped_inner ** 2, 0.0)
        return float(
            np.dot(self.lower_endcap_flux_m2_s[species_index], area)
            / (np.pi * radius ** 2))


@dataclass(frozen=True)
class AxisymmetricInventoryLiftResult:
    """Spatial state conditioned only on a declared global inventory.

    This is a dimensional lift, not an independent chemistry prediction.  The
    inferred source amplitudes are returned explicitly so a downstream audit
    cannot confuse a volume-average constraint with a measured local source.
    """

    solution: AxisymmetricReactionDiffusionSolution
    target_volume_average_density_m3: np.ndarray
    recovered_volume_average_density_m3: np.ndarray
    inferred_source_amplitude_m3_s: np.ndarray
    source_response_condition_number: float
    maximum_inventory_relative_residual: float
    supports_reactor_state_prediction: bool = False

    def __post_init__(self):
        if not isinstance(self.solution, AxisymmetricReactionDiffusionSolution):
            raise TypeError("inventory lift requires a spatial solution")
        species_count = len(self.solution.condition.species_names)
        for name in (
            "target_volume_average_density_m3",
            "recovered_volume_average_density_m3",
            "inferred_source_amplitude_m3_s",
        ):
            array = np.asarray(getattr(self, name), dtype=float).copy()
            if (
                array.shape != (species_count,)
                or np.any(~np.isfinite(array))
                or np.any(array < 0.0)
            ):
                raise ValueError("invalid inventory-lift state")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if (
            not math.isfinite(self.source_response_condition_number)
            or self.source_response_condition_number < 1.0
            or not math.isfinite(self.maximum_inventory_relative_residual)
            or not 0.0 <= self.maximum_inventory_relative_residual < 1.0e-10
            or bool(self.supports_reactor_state_prediction)
        ):
            raise ValueError("inventory-lift certification failed")


def _effective_robin_velocity(
    diffusion_m2_s: float,
    wall_velocity_m_s: float,
    half_cell_distance_m: float,
) -> float:
    if math.isinf(wall_velocity_m_s):
        return diffusion_m2_s / half_cell_distance_m
    if wall_velocity_m_s == 0.0:
        return 0.0
    return 1.0 / (
        1.0 / wall_velocity_m_s
        + half_cell_distance_m / diffusion_m2_s)


class DeterministicAxisymmetricReactionDiffusion:
    """Sparse factorized solver with exact source-direction derivatives."""

    def __init__(self, condition: AxisymmetricReactionDiffusionCondition):
        if not isinstance(condition, AxisymmetricReactionDiffusionCondition):
            raise TypeError("a reaction-diffusion condition is required")
        self.condition = condition
        self._operator, self._boundary_velocity = self._assemble_operator()
        try:
            self._factorization = splu(self._operator)
        except RuntimeError as error:
            raise ValueError(
                "reaction-diffusion operator is singular; add a volume or wall loss"
            ) from error
        pivot_magnitude = np.abs(self._factorization.U.diagonal())
        pivot_scale = max(float(np.max(pivot_magnitude)), 1.0)
        if (
            pivot_magnitude.size == 0
            or float(np.min(pivot_magnitude)) <= 1.0e-13 * pivot_scale
        ):
            raise ValueError(
                "reaction-diffusion operator is singular; add a volume or wall loss"
            )

    def _index(self, species: int, radial: int, axial: int) -> int:
        grid = self.condition.grid
        return (
            species * grid.radial_cell_count * grid.axial_cell_count
            + radial * grid.axial_cell_count
            + axial)

    def _assemble_operator(self) -> tuple[csc_matrix, np.ndarray]:
        condition = self.condition
        grid = condition.grid
        species_count = len(condition.species_names)
        nr = grid.radial_cell_count
        nz = grid.axial_cell_count
        cell_count = nr * nz
        operator = lil_matrix((species_count * cell_count,) * 2, dtype=float)
        volume = grid.cell_volume_m3
        radial_center = grid.radial_centers_m
        axial_center = grid.axial_centers_m
        axial_area = grid.axial_face_area_m2
        wall_effective = np.zeros((species_count, 3))
        for species in range(species_count):
            diffusion = condition.diffusion_coefficient_m2_s[species]
            lower_velocity = _effective_robin_velocity(
                diffusion,
                condition.wall_velocity_m_s[species, 0],
                axial_center[0] - grid.axial_edges_m[0],
            )
            upper_velocity = _effective_robin_velocity(
                diffusion,
                condition.wall_velocity_m_s[species, 1],
                grid.axial_edges_m[-1] - axial_center[-1],
            )
            side_velocity = _effective_robin_velocity(
                diffusion,
                condition.wall_velocity_m_s[species, 2],
                grid.radial_edges_m[-1] - radial_center[-1],
            )
            wall_effective[species] = (
                lower_velocity, upper_velocity, side_velocity)
            for radial in range(nr):
                for axial in range(nz):
                    row = self._index(species, radial, axial)
                    diagonal = 0.0
                    if radial > 0:
                        face_radius = grid.radial_edges_m[radial]
                        face_area = (
                            2.0 * np.pi * face_radius
                            * (grid.axial_edges_m[axial + 1]
                               - grid.axial_edges_m[axial]))
                        conductance = (
                            diffusion * face_area
                            / (radial_center[radial]
                               - radial_center[radial - 1]))
                        diagonal += conductance
                        operator[row, self._index(
                            species, radial - 1, axial)] -= conductance
                    if radial < nr - 1:
                        face_radius = grid.radial_edges_m[radial + 1]
                        face_area = (
                            2.0 * np.pi * face_radius
                            * (grid.axial_edges_m[axial + 1]
                               - grid.axial_edges_m[axial]))
                        conductance = (
                            diffusion * face_area
                            / (radial_center[radial + 1]
                               - radial_center[radial]))
                        diagonal += conductance
                        operator[row, self._index(
                            species, radial + 1, axial)] -= conductance
                    else:
                        conductance = (
                            side_velocity
                            * grid.outer_radial_face_area_m2[axial])
                        diagonal += conductance
                    if axial > 0:
                        conductance = (
                            diffusion * axial_area[radial]
                            / (axial_center[axial]
                               - axial_center[axial - 1]))
                        diagonal += conductance
                        operator[row, self._index(
                            species, radial, axial - 1)] -= conductance
                    else:
                        diagonal += lower_velocity * axial_area[radial]
                    if axial < nz - 1:
                        conductance = (
                            diffusion * axial_area[radial]
                            / (axial_center[axial + 1]
                               - axial_center[axial]))
                        diagonal += conductance
                        operator[row, self._index(
                            species, radial, axial + 1)] -= conductance
                    else:
                        diagonal += upper_velocity * axial_area[radial]
                    for coupled_species in range(species_count):
                        coefficient = condition.volume_reaction_matrix_s_inv[
                            species, coupled_species] * volume[radial, axial]
                        if coupled_species == species:
                            diagonal += coefficient
                        elif coefficient != 0.0:
                            operator[row, self._index(
                                coupled_species, radial, axial)] += coefficient
                    operator[row, row] += diagonal
        return operator.tocsc(), wall_effective

    def _source_vector(self, source_rate_m3_s: np.ndarray) -> np.ndarray:
        return np.asarray(
            source_rate_m3_s * self.condition.grid.cell_volume_m3[None, :, :],
            dtype=float,
        ).ravel()

    def solve(self) -> AxisymmetricReactionDiffusionSolution:
        rhs = self._source_vector(self.condition.source_rate_m3_s)
        raw_density = self._factorization.solve(rhs)
        maximum = max(1.0, float(np.max(np.abs(raw_density))))
        if np.min(raw_density) < -1.0e-10 * maximum:
            raise RuntimeError("reaction-diffusion solve lost positivity")
        raw_density[np.abs(raw_density) < 1.0e-13 * maximum] = 0.0
        density = np.maximum(raw_density, 0.0).reshape(
            self.condition.source_rate_m3_s.shape)
        linear_residual = float(
            np.linalg.norm(self._operator @ density.ravel() - rhs)
            / max(np.linalg.norm(rhs), 1.0))
        grid = self.condition.grid
        lower_flux = self._boundary_velocity[:, 0, None] * density[:, :, 0]
        upper_flux = self._boundary_velocity[:, 1, None] * density[:, :, -1]
        side_flux = self._boundary_velocity[:, 2, None] * density[:, -1, :]
        integrated_source = np.sum(
            self.condition.source_rate_m3_s
            * grid.cell_volume_m3[None, :, :],
            axis=(1, 2),
        )
        integrated_reaction = np.einsum(
            "st,tij,ij->s",
            self.condition.volume_reaction_matrix_s_inv,
            density,
            grid.cell_volume_m3,
        )
        integrated_wall = (
            np.einsum("si,i->s", lower_flux, grid.axial_face_area_m2)
            + np.einsum("si,i->s", upper_flux, grid.axial_face_area_m2)
            + np.einsum(
                "sj,j->s", side_flux, grid.outer_radial_face_area_m2)
        )
        ledger = np.abs(
            integrated_source - integrated_reaction - integrated_wall)
        ledger_scale = np.maximum.reduce((
            np.abs(integrated_source),
            np.abs(integrated_reaction) + np.abs(integrated_wall),
            np.ones_like(integrated_source),
        ))
        maximum_ledger = float(np.max(ledger / ledger_scale))
        return AxisymmetricReactionDiffusionSolution(
            condition=self.condition,
            density_m3=density,
            lower_endcap_flux_m2_s=lower_flux,
            upper_endcap_flux_m2_s=upper_flux,
            sidewall_flux_m2_s=side_flux,
            integrated_source_rate_s=integrated_source,
            integrated_volume_reaction_rate_s=integrated_reaction,
            integrated_wall_loss_rate_s=integrated_wall,
            maximum_species_ledger_relative_residual=maximum_ledger,
            linear_system_relative_residual=linear_residual,
            minimum_density_m3=float(np.min(density)),
        )

    def source_jvp(self, source_rate_tangent_m3_s: np.ndarray) -> np.ndarray:
        """Exact JVP of density with respect to the distributed source."""
        tangent = np.asarray(source_rate_tangent_m3_s, dtype=float)
        if (
            tangent.shape != self.condition.source_rate_m3_s.shape
            or np.any(~np.isfinite(tangent))
        ):
            raise ValueError("source tangent shape or values are invalid")
        return self._factorization.solve(
            self._source_vector(tangent)).reshape(tangent.shape)


class DeterministicAxisymmetricInventoryLift:
    """Infer nonnegative spatial source amplitudes from global inventories.

    One nonnegative, volume-average-one source moment is declared per species.
    Because the transport operator is linear, a small dense response matrix
    maps those source amplitudes to global inventories.  Its inversion enforces
    the supplied 0-D state exactly while the finite-volume solve predicts the
    corresponding local wafer/sidewall partition.  Negative inferred sources
    fail closed: they prove the chosen moments/reaction network are
    incompatible with the supplied inventory.
    """

    def __init__(
        self,
        *,
        grid: AxisymmetricFiniteVolumeGrid,
        species_names: tuple[str, ...],
        diffusion_coefficient_m2_s: np.ndarray,
        volume_reaction_matrix_s_inv: np.ndarray,
        wall_velocity_m_s: np.ndarray,
        source_shape: np.ndarray,
        source: str,
    ):
        names = tuple(species_names)
        shape = np.asarray(source_shape, dtype=float).copy()
        expected_shape = (
            len(names), grid.radial_cell_count, grid.axial_cell_count)
        if (
            shape.shape != expected_shape
            or np.any(~np.isfinite(shape))
            or np.any(shape < 0.0)
        ):
            raise ValueError("invalid inventory-lift source moments")
        volume = grid.cell_volume_m3
        means = np.sum(shape * volume[None, :, :], axis=(1, 2)) / (
            grid.geometry.volume_m3)
        if np.any(np.abs(means - 1.0) > 1.0e-10):
            raise ValueError(
                "inventory-lift source moments must have volume average one")
        shape.setflags(write=False)
        self.source_shape = shape
        self.source = str(source).strip()
        if not self.source:
            raise ValueError("inventory-lift source provenance is required")
        zero_condition = AxisymmetricReactionDiffusionCondition(
            grid=grid,
            species_names=names,
            diffusion_coefficient_m2_s=diffusion_coefficient_m2_s,
            volume_reaction_matrix_s_inv=volume_reaction_matrix_s_inv,
            source_rate_m3_s=np.zeros(expected_shape),
            wall_velocity_m_s=wall_velocity_m_s,
            source=self.source,
        )
        self._transport = DeterministicAxisymmetricReactionDiffusion(
            zero_condition)
        species_count = len(names)
        response = np.zeros((species_count, species_count))
        for source_species in range(species_count):
            tangent = np.zeros(expected_shape)
            tangent[source_species] = shape[source_species]
            density = self._transport.source_jvp(tangent)
            response[:, source_species] = np.sum(
                density * volume[None, :, :], axis=(1, 2)
            ) / grid.geometry.volume_m3
        condition_number = float(np.linalg.cond(response))
        if not math.isfinite(condition_number) or condition_number > 1.0e12:
            raise ValueError(
                "inventory-lift source response is singular or ill-conditioned")
        self._source_response_s_m3 = response
        self.source_response_condition_number = condition_number

    @property
    def condition(self) -> AxisymmetricReactionDiffusionCondition:
        return self._transport.condition

    def _source_amplitudes(self, target_density_m3: np.ndarray) -> np.ndarray:
        target = np.asarray(target_density_m3, dtype=float)
        species_count = len(self.condition.species_names)
        if (
            target.shape != (species_count,)
            or np.any(~np.isfinite(target))
            or np.any(target < 0.0)
        ):
            raise ValueError("invalid global inventory target")
        amplitude = np.linalg.solve(self._source_response_s_m3, target)
        scale = max(1.0, float(np.max(np.abs(amplitude))))
        if float(np.min(amplitude)) < -1.0e-11 * scale:
            raise ValueError(
                "global inventory requires a negative source amplitude; "
                "the spatial moments or reaction network are incompatible")
        return np.maximum(amplitude, 0.0)

    def solve(
        self,
        target_volume_average_density_m3: np.ndarray,
    ) -> AxisymmetricInventoryLiftResult:
        target = np.asarray(
            target_volume_average_density_m3, dtype=float).copy()
        amplitude = self._source_amplitudes(target)
        source_rate = amplitude[:, None, None] * self.source_shape
        base = self.condition
        condition = AxisymmetricReactionDiffusionCondition(
            grid=base.grid,
            species_names=base.species_names,
            diffusion_coefficient_m2_s=base.diffusion_coefficient_m2_s,
            volume_reaction_matrix_s_inv=base.volume_reaction_matrix_s_inv,
            source_rate_m3_s=source_rate,
            wall_velocity_m_s=base.wall_velocity_m_s,
            source=(self.source + "; amplitudes inferred from declared 0-D inventory"),
        )
        solution = DeterministicAxisymmetricReactionDiffusion(condition).solve()
        volume = base.grid.cell_volume_m3
        recovered = np.sum(
            solution.density_m3 * volume[None, :, :], axis=(1, 2)
        ) / base.grid.geometry.volume_m3
        residual = float(np.max(
            np.abs(recovered - target) / np.maximum(np.abs(target), 1.0)))
        return AxisymmetricInventoryLiftResult(
            solution=solution,
            target_volume_average_density_m3=target,
            recovered_volume_average_density_m3=recovered,
            inferred_source_amplitude_m3_s=amplitude,
            source_response_condition_number=(
                self.source_response_condition_number),
            maximum_inventory_relative_residual=residual,
        )

    def target_inventory_jvp(
        self,
        target_density_tangent_m3: np.ndarray,
    ) -> np.ndarray:
        """Exact density JVP with respect to a global-inventory direction."""
        tangent = np.asarray(target_density_tangent_m3, dtype=float)
        species_count = len(self.condition.species_names)
        if tangent.shape != (species_count,) or np.any(~np.isfinite(tangent)):
            raise ValueError("invalid inventory tangent")
        amplitude_tangent = np.linalg.solve(
            self._source_response_s_m3, tangent)
        return self._transport.source_jvp(
            amplitude_tangent[:, None, None] * self.source_shape)


def normalized_exponential_skin_source(
    grid: AxisymmetricFiniteVolumeGrid,
    *,
    axial_skin_depth_m: float,
    radial_scale_m: float,
    radial_power: float = 2.0,
) -> np.ndarray:
    """Return a volume-average-one ICP-like top/center source moment."""
    if not isinstance(grid, AxisymmetricFiniteVolumeGrid):
        raise TypeError("an axisymmetric finite-volume grid is required")
    skin = float(axial_skin_depth_m)
    radial_scale = float(radial_scale_m)
    power = float(radial_power)
    if (
        not math.isfinite(skin)
        or skin <= 0.0
        or not math.isfinite(radial_scale)
        or radial_scale <= 0.0
        or not math.isfinite(power)
        or power <= 0.0
    ):
        raise ValueError("invalid axisymmetric skin-source scale")
    radial = grid.radial_centers_m[:, None]
    axial = grid.axial_centers_m[None, :]
    raw = np.exp(
        -(grid.geometry.length_m - axial) / skin
        - (radial / radial_scale) ** power
    )
    weighted_mean = float(np.sum(raw * grid.cell_volume_m3)
                          / grid.geometry.volume_m3)
    if weighted_mean <= 0.0:
        raise RuntimeError("axisymmetric source normalization underflowed")
    return raw / weighted_mean


def normalized_annular_skin_source(
    grid: AxisymmetricFiniteVolumeGrid,
    *,
    axial_skin_depth_m: float,
    ring_radius_m: float,
    radial_width_m: float,
) -> np.ndarray:
    """Return a volume-average-one top/annular ICP source moment.

    Planar induction coils deposit power in a toroidal region rather than at
    the symmetry axis.  This fixed Gaussian-ring/exponential-skin moment is a
    deterministic sensitivity interface for an EM solution or a measured
    deposition map; it is not itself an equipment calibration.
    """
    if not isinstance(grid, AxisymmetricFiniteVolumeGrid):
        raise TypeError("an axisymmetric finite-volume grid is required")
    skin = float(axial_skin_depth_m)
    ring = float(ring_radius_m)
    width = float(radial_width_m)
    if (
        not math.isfinite(skin)
        or skin <= 0.0
        or not math.isfinite(ring)
        or not 0.0 <= ring <= grid.geometry.radius_m
        or not math.isfinite(width)
        or width <= 0.0
    ):
        raise ValueError("invalid annular skin-source scale")
    radial = grid.radial_centers_m[:, None]
    axial = grid.axial_centers_m[None, :]
    raw = np.exp(
        -(grid.geometry.length_m - axial) / skin
        - 0.5 * ((radial - ring) / width) ** 2
    )
    weighted_mean = float(
        np.sum(raw * grid.cell_volume_m3) / grid.geometry.volume_m3)
    if weighted_mean <= 0.0:
        raise RuntimeError("annular source normalization underflowed")
    return raw / weighted_mean
