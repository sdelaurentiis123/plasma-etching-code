"""Conservative deterministic cylindrical 3-D inventory-to-wafer transport.

This tier resolves ``(r, phi, z)`` with periodic azimuth.  Each effective
species diffuses independently after the product-resolved 0-D chemistry has
set its volume-average inventory.  Robin loss is applied at the lower/upper
endcaps and cylindrical sidewall.  The finite-volume operator is an M-matrix,
so a nonnegative source produces a nonnegative density.

The source-to-inventory map is linear.  One sparse factorization per species
therefore supplies the primal unit-inventory response and an exact JVP for any
time-dependent inventory perturbation.  This is the deterministic 3-D upgrade
of the axisymmetric inventory lift; no particle Monte Carlo is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .axisymmetric_reaction_diffusion import _effective_robin_velocity
from .geometry import CylindricalReactor


@dataclass(frozen=True)
class CylindricalFiniteVolumeGrid:
    radial_edges_m: np.ndarray
    azimuthal_edges_rad: np.ndarray
    axial_edges_m: np.ndarray

    def __post_init__(self):
        radial = np.asarray(self.radial_edges_m, dtype=float).copy()
        azimuthal = np.asarray(self.azimuthal_edges_rad, dtype=float).copy()
        axial = np.asarray(self.axial_edges_m, dtype=float).copy()
        if (
            radial.ndim != 1 or azimuthal.ndim != 1 or axial.ndim != 1
            or min(radial.size, axial.size) < 3 or azimuthal.size < 5
            or radial[0] != 0.0 or axial[0] != 0.0
            or azimuthal[0] != 0.0
            or not np.isclose(azimuthal[-1], 2.0 * np.pi)
            or any(np.any(~np.isfinite(value)) for value in (
                radial, azimuthal, axial))
            or any(np.any(np.diff(value) <= 0.0) for value in (
                radial, azimuthal, axial))
            or not np.allclose(
                np.diff(azimuthal), np.diff(azimuthal)[0],
                rtol=0.0, atol=2.0e-14)
        ):
            raise ValueError("invalid cylindrical finite-volume grid")
        for value in (radial, azimuthal, axial):
            value.setflags(write=False)
        object.__setattr__(self, "radial_edges_m", radial)
        object.__setattr__(self, "azimuthal_edges_rad", azimuthal)
        object.__setattr__(self, "axial_edges_m", axial)

    @classmethod
    def uniform(cls, geometry, *, radial_cell_count, azimuthal_cell_count,
                axial_cell_count):
        if not isinstance(geometry, CylindricalReactor):
            raise TypeError("a cylindrical reactor is required")
        nr = int(radial_cell_count)
        nphi = int(azimuthal_cell_count)
        nz = int(axial_cell_count)
        if (nr != radial_cell_count or nphi != azimuthal_cell_count
                or nz != axial_cell_count or min(nr, nz) < 2 or nphi < 4):
            raise ValueError("cylindrical grid requires at least 2x4x2 cells")
        return cls(
            np.linspace(0.0, geometry.radius_m, nr + 1),
            np.linspace(0.0, 2.0 * np.pi, nphi + 1),
            np.linspace(0.0, geometry.length_m, nz + 1))

    @property
    def geometry(self):
        return CylindricalReactor(
            radius_m=float(self.radial_edges_m[-1]),
            length_m=float(self.axial_edges_m[-1]))

    @property
    def radial_cell_count(self):
        return self.radial_edges_m.size - 1

    @property
    def azimuthal_cell_count(self):
        return self.azimuthal_edges_rad.size - 1

    @property
    def axial_cell_count(self):
        return self.axial_edges_m.size - 1

    @property
    def radial_centers_m(self):
        return 0.5 * (self.radial_edges_m[:-1] + self.radial_edges_m[1:])

    @property
    def azimuthal_centers_rad(self):
        return 0.5 * (
            self.azimuthal_edges_rad[:-1] + self.azimuthal_edges_rad[1:])

    @property
    def axial_centers_m(self):
        return 0.5 * (self.axial_edges_m[:-1] + self.axial_edges_m[1:])

    @property
    def cell_volume_m3(self):
        sector_area = 0.5 * (
            self.radial_edges_m[1:] ** 2
            - self.radial_edges_m[:-1] ** 2)[:, None]
        dphi = np.diff(self.azimuthal_edges_rad)[None, :]
        dz = np.diff(self.axial_edges_m)[None, None, :]
        return sector_area[:, :, None] * dphi[:, :, None] * dz

    @property
    def endcap_cell_area_m2(self):
        sector_area = 0.5 * (
            self.radial_edges_m[1:] ** 2
            - self.radial_edges_m[:-1] ** 2)[:, None]
        return sector_area * np.diff(self.azimuthal_edges_rad)[None, :]

    @property
    def sidewall_cell_area_m2(self):
        return (
            self.geometry.radius_m
            * np.diff(self.azimuthal_edges_rad)[:, None]
            * np.diff(self.axial_edges_m)[None, :])


@dataclass(frozen=True)
class CylindricalInventoryLiftResult:
    species_names: tuple[str, ...]
    grid: CylindricalFiniteVolumeGrid
    target_volume_average_density_m3: np.ndarray
    density_m3: np.ndarray
    lower_endcap_flux_m2_s: np.ndarray
    upper_endcap_flux_m2_s: np.ndarray
    sidewall_flux_m2_s: np.ndarray
    inferred_source_amplitude_m3_s: np.ndarray
    maximum_species_ledger_relative_residual: float
    maximum_inventory_relative_residual: float
    maximum_linear_system_relative_residual: float
    minimum_density_m3: float
    source_jvp_supported: bool = True

    def __post_init__(self):
        names = tuple(self.species_names)
        grid = self.grid
        ns = len(names)
        shapes = {
            "target_volume_average_density_m3": (ns,),
            "density_m3": (
                ns, grid.radial_cell_count, grid.azimuthal_cell_count,
                grid.axial_cell_count),
            "lower_endcap_flux_m2_s": (
                ns, grid.radial_cell_count, grid.azimuthal_cell_count),
            "upper_endcap_flux_m2_s": (
                ns, grid.radial_cell_count, grid.azimuthal_cell_count),
            "sidewall_flux_m2_s": (
                ns, grid.azimuthal_cell_count, grid.axial_cell_count),
            "inferred_source_amplitude_m3_s": (ns,),
        }
        if not names or len(set(names)) != ns or any(not name for name in names):
            raise ValueError("invalid cylindrical inventory species")
        for name, shape in shapes.items():
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if (value.shape != shape or np.any(~np.isfinite(value))
                    or np.any(value < 0.0)):
                raise ValueError("invalid cylindrical inventory result")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if (not 0.0 <= self.maximum_species_ledger_relative_residual < 1.0e-8
                or not 0.0 <= self.maximum_inventory_relative_residual < 1.0e-10
                or not 0.0 <= self.maximum_linear_system_relative_residual < 1.0e-10
                or not math.isfinite(self.minimum_density_m3)
                or self.minimum_density_m3 < 0.0
                or not bool(self.source_jvp_supported)):
            raise ValueError("cylindrical inventory certification failed")


class DeterministicCylindricalIndependentInventoryLift:
    """Factorized independent-species 3-D diffusion and wall-loss lift."""

    def __init__(self, *, grid, species_names, diffusion_coefficient_m2_s,
                 wall_velocity_m_s, source_shape, source):
        if not isinstance(grid, CylindricalFiniteVolumeGrid):
            raise TypeError("a cylindrical finite-volume grid is required")
        names = tuple(str(name) for name in species_names)
        diffusion = np.asarray(diffusion_coefficient_m2_s, dtype=float).copy()
        wall = np.asarray(wall_velocity_m_s, dtype=float).copy()
        shape = np.asarray(source_shape, dtype=float).copy()
        expected = (
            len(names), grid.radial_cell_count, grid.azimuthal_cell_count,
            grid.axial_cell_count)
        if (not names or len(set(names)) != len(names)
                or diffusion.shape != (len(names),)
                or wall.shape != (len(names), 3)
                or shape.shape != expected
                or np.any(~np.isfinite(diffusion)) or np.any(diffusion <= 0.0)
                or np.any(~np.isfinite(wall)) or np.any(wall <= 0.0)
                or np.any(~np.isfinite(shape)) or np.any(shape < 0.0)
                or not str(source).strip()):
            raise ValueError("invalid cylindrical inventory-lift inputs")
        volume = grid.cell_volume_m3
        means = np.sum(shape * volume[None], axis=(1, 2, 3)) / grid.geometry.volume_m3
        if np.any(np.abs(means - 1.0) > 1.0e-10):
            raise ValueError("cylindrical source moments must average to one")
        self.grid = grid
        self.species_names = names
        self.diffusion_coefficient_m2_s = diffusion
        self.wall_velocity_m_s = wall
        self.source_shape = shape
        self.source = str(source)
        self._factorizations = []
        self._operators = []
        self._boundary_velocity = []
        unit_density = []
        unit_lower = []
        unit_upper = []
        unit_side = []
        unit_amplitude = []
        unit_ledger = []
        unit_linear = []
        for species in range(len(names)):
            operator, effective = self._assemble_species_operator(species)
            factor = splu(operator)
            rhs = (shape[species] * volume).ravel()
            raw = factor.solve(rhs)
            scale = max(1.0, float(np.max(np.abs(raw))))
            if np.min(raw) < -1.0e-10 * scale:
                raise RuntimeError("cylindrical inventory solve lost positivity")
            raw[np.abs(raw) < 1.0e-13 * scale] = 0.0
            density = np.maximum(raw, 0.0).reshape(shape[species].shape)
            response = float(np.sum(density * volume) / grid.geometry.volume_m3)
            if not math.isfinite(response) or response <= 0.0:
                raise ValueError("cylindrical inventory response is singular")
            amplitude = 1.0 / response
            density *= amplitude
            lower = effective[0] * density[:, :, 0]
            upper = effective[1] * density[:, :, -1]
            side = effective[2] * density[-1]
            source_rate = amplitude * float(np.sum(shape[species] * volume))
            wall_rate = (
                float(np.sum(lower * grid.endcap_cell_area_m2))
                + float(np.sum(upper * grid.endcap_cell_area_m2))
                + float(np.sum(side * grid.sidewall_cell_area_m2)))
            ledger = abs(source_rate - wall_rate) / max(source_rate, wall_rate, 1.0)
            linear = float(
                np.linalg.norm(operator @ density.ravel() - amplitude * rhs)
                / max(np.linalg.norm(amplitude * rhs), 1.0))
            self._operators.append(operator)
            self._factorizations.append(factor)
            self._boundary_velocity.append(effective)
            unit_density.append(density)
            unit_lower.append(lower)
            unit_upper.append(upper)
            unit_side.append(side)
            unit_amplitude.append(amplitude)
            unit_ledger.append(ledger)
            unit_linear.append(linear)
        self._unit_density_per_m3 = np.stack(unit_density)
        self._unit_lower_flux_per_density_m_s = np.stack(unit_lower)
        self._unit_upper_flux_per_density_m_s = np.stack(unit_upper)
        self._unit_side_flux_per_density_m_s = np.stack(unit_side)
        self._unit_source_amplitude_s_inv = np.asarray(unit_amplitude)
        self.maximum_unit_ledger_relative_residual = float(max(unit_ledger))
        self.maximum_unit_linear_system_relative_residual = float(max(unit_linear))

    def _index(self, radial, azimuthal, axial):
        nphi = self.grid.azimuthal_cell_count
        nz = self.grid.axial_cell_count
        return (radial * nphi + azimuthal) * nz + axial

    @property
    def unit_lower_flux_per_density_m_s(self):
        """Read-only lower-wall flux response to unit average inventory."""
        view = self._unit_lower_flux_per_density_m_s.view()
        view.setflags(write=False)
        return view

    def _assemble_species_operator(self, species):
        grid = self.grid
        nr = grid.radial_cell_count
        nphi = grid.azimuthal_cell_count
        nz = grid.axial_cell_count
        operator = lil_matrix((nr * nphi * nz,) * 2, dtype=float)
        D = float(self.diffusion_coefficient_m2_s[species])
        radial_center = grid.radial_centers_m
        axial_center = grid.axial_centers_m
        dr = np.diff(grid.radial_edges_m)
        dphi = float(np.diff(grid.azimuthal_edges_rad)[0])
        dz = np.diff(grid.axial_edges_m)
        lower_velocity = _effective_robin_velocity(
            D, self.wall_velocity_m_s[species, 0],
            axial_center[0] - grid.axial_edges_m[0])
        upper_velocity = _effective_robin_velocity(
            D, self.wall_velocity_m_s[species, 1],
            grid.axial_edges_m[-1] - axial_center[-1])
        side_velocity = _effective_robin_velocity(
            D, self.wall_velocity_m_s[species, 2],
            grid.radial_edges_m[-1] - radial_center[-1])
        for radial in range(nr):
            for azimuthal in range(nphi):
                left_phi = (azimuthal - 1) % nphi
                right_phi = (azimuthal + 1) % nphi
                for axial in range(nz):
                    row = self._index(radial, azimuthal, axial)
                    diagonal = 0.0
                    if radial > 0:
                        area = grid.radial_edges_m[radial] * dphi * dz[axial]
                        conductance = D * area / (
                            radial_center[radial] - radial_center[radial - 1])
                        diagonal += conductance
                        operator[row, self._index(
                            radial - 1, azimuthal, axial)] -= conductance
                    if radial < nr - 1:
                        area = grid.radial_edges_m[radial + 1] * dphi * dz[axial]
                        conductance = D * area / (
                            radial_center[radial + 1] - radial_center[radial])
                        diagonal += conductance
                        operator[row, self._index(
                            radial + 1, azimuthal, axial)] -= conductance
                    else:
                        diagonal += (
                            side_velocity
                            * grid.sidewall_cell_area_m2[azimuthal, axial])
                    azimuth_area = dr[radial] * dz[axial]
                    azimuth_distance = radial_center[radial] * dphi
                    azimuth_conductance = D * azimuth_area / azimuth_distance
                    diagonal += 2.0 * azimuth_conductance
                    operator[row, self._index(
                        radial, left_phi, axial)] -= azimuth_conductance
                    operator[row, self._index(
                        radial, right_phi, axial)] -= azimuth_conductance
                    endcap_area = grid.endcap_cell_area_m2[radial, azimuthal]
                    if axial > 0:
                        conductance = D * endcap_area / (
                            axial_center[axial] - axial_center[axial - 1])
                        diagonal += conductance
                        operator[row, self._index(
                            radial, azimuthal, axial - 1)] -= conductance
                    else:
                        diagonal += lower_velocity * endcap_area
                    if axial < nz - 1:
                        conductance = D * endcap_area / (
                            axial_center[axial + 1] - axial_center[axial])
                        diagonal += conductance
                        operator[row, self._index(
                            radial, azimuthal, axial + 1)] -= conductance
                    else:
                        diagonal += upper_velocity * endcap_area
                    operator[row, row] += diagonal
        return operator.tocsc(), np.asarray(
            (lower_velocity, upper_velocity, side_velocity))

    def solve(self, target_volume_average_density_m3):
        target = np.asarray(target_volume_average_density_m3, dtype=float)
        if (target.shape != (len(self.species_names),)
                or np.any(~np.isfinite(target)) or np.any(target < 0.0)):
            raise ValueError("invalid cylindrical target inventory")
        density = target[:, None, None, None] * self._unit_density_per_m3
        volume = self.grid.cell_volume_m3
        recovered = np.sum(density * volume[None], axis=(1, 2, 3)) / (
            self.grid.geometry.volume_m3)
        inventory_residual = float(np.max(
            np.abs(recovered - target) / np.maximum(target, 1.0)))
        return CylindricalInventoryLiftResult(
            species_names=self.species_names,
            grid=self.grid,
            target_volume_average_density_m3=target,
            density_m3=density,
            lower_endcap_flux_m2_s=(
                target[:, None, None] * self._unit_lower_flux_per_density_m_s),
            upper_endcap_flux_m2_s=(
                target[:, None, None] * self._unit_upper_flux_per_density_m_s),
            sidewall_flux_m2_s=(
                target[:, None, None] * self._unit_side_flux_per_density_m_s),
            inferred_source_amplitude_m3_s=(
                target * self._unit_source_amplitude_s_inv),
            maximum_species_ledger_relative_residual=(
                self.maximum_unit_ledger_relative_residual),
            maximum_inventory_relative_residual=inventory_residual,
            maximum_linear_system_relative_residual=(
                self.maximum_unit_linear_system_relative_residual),
            minimum_density_m3=float(np.min(density)))

    def target_inventory_jvp(self, target_density_tangent_m3):
        tangent = np.asarray(target_density_tangent_m3, dtype=float)
        if (tangent.shape != (len(self.species_names),)
                or np.any(~np.isfinite(tangent))):
            raise ValueError("invalid cylindrical inventory tangent")
        return tangent[:, None, None, None] * self._unit_density_per_m3

    def density_to_lower_flux_jvp(self, density_tangent_m3):
        tangent = np.asarray(density_tangent_m3, dtype=float)
        if (tangent.ndim != 2 or tangent.shape[1] != len(self.species_names)
                or np.any(~np.isfinite(tangent))):
            raise ValueError("density tangent must have shape (time, species)")
        return tangent[:, :, None, None] * self._unit_lower_flux_per_density_m_s

    def source_shape_to_unit_lower_flux(self, source_shape):
        """Solve a new positive source moment without refactorizing transport.

        Diffusion and Robin wall losses define the sparse operator; source
        moments enter only through its right-hand side.  Calibration can
        therefore vary a physical source map while reusing the exact same
        factorization.  The returned response is normalized to one unit of
        volume-average density for each species.
        """
        shape = np.asarray(source_shape, dtype=float)
        expected = (
            len(self.species_names), self.grid.radial_cell_count,
            self.grid.azimuthal_cell_count, self.grid.axial_cell_count)
        volume = self.grid.cell_volume_m3
        if (shape.shape != expected or np.any(~np.isfinite(shape))
                or np.any(shape < 0.0)):
            raise ValueError("invalid cylindrical source moments")
        means = np.sum(shape * volume[None], axis=(1, 2, 3)) / (
            self.grid.geometry.volume_m3)
        if np.any(np.abs(means - 1.0) > 1.0e-10):
            raise ValueError("cylindrical source moments must average to one")

        unit_lower = []
        ledgers = []
        linears = []
        for species, (factor, operator, effective) in enumerate(zip(
                self._factorizations, self._operators, self._boundary_velocity)):
            rhs = (shape[species] * volume).ravel()
            raw = factor.solve(rhs)
            scale = max(1.0, float(np.max(np.abs(raw))))
            if np.min(raw) < -1.0e-10 * scale:
                raise RuntimeError("cylindrical inventory solve lost positivity")
            raw[np.abs(raw) < 1.0e-13 * scale] = 0.0
            density = np.maximum(raw, 0.0).reshape(shape[species].shape)
            response = float(
                np.sum(density * volume) / self.grid.geometry.volume_m3)
            if not math.isfinite(response) or response <= 0.0:
                raise ValueError("cylindrical inventory response is singular")
            amplitude = 1.0 / response
            density *= amplitude
            lower = effective[0] * density[:, :, 0]
            upper = effective[1] * density[:, :, -1]
            side = effective[2] * density[-1]
            source_rate = amplitude * float(np.sum(shape[species] * volume))
            wall_rate = (
                float(np.sum(lower * self.grid.endcap_cell_area_m2))
                + float(np.sum(upper * self.grid.endcap_cell_area_m2))
                + float(np.sum(side * self.grid.sidewall_cell_area_m2)))
            ledgers.append(
                abs(source_rate - wall_rate) / max(source_rate, wall_rate, 1.0))
            linears.append(float(
                np.linalg.norm(operator @ density.ravel() - amplitude * rhs)
                / max(np.linalg.norm(amplitude * rhs), 1.0)))
            unit_lower.append(lower)
        output = np.stack(unit_lower)
        output.setflags(write=False)
        return output, float(max(ledgers)), float(max(linears))


def normalized_cylindrical_annular_skin_source(
        grid, *, axial_skin_depth_m, ring_radius_m, radial_width_m,
        cosine_coefficients=(), sine_coefficients=()):
    """Positive volume-average-one annular source with azimuthal harmonics."""
    if not isinstance(grid, CylindricalFiniteVolumeGrid):
        raise TypeError("a cylindrical finite-volume grid is required")
    skin = float(axial_skin_depth_m)
    ring = float(ring_radius_m)
    width = float(radial_width_m)
    cosine = tuple(float(value) for value in cosine_coefficients)
    sine = tuple(float(value) for value in sine_coefficients)
    if (not math.isfinite(skin) or skin <= 0.0
            or not math.isfinite(ring) or not 0.0 <= ring <= grid.geometry.radius_m
            or not math.isfinite(width) or width <= 0.0
            or len(cosine) != len(sine)
            or any(not math.isfinite(value) or abs(value) > 5.0
                   for value in cosine + sine)):
        raise ValueError("invalid cylindrical source parameters")
    radial = grid.radial_centers_m[:, None, None]
    azimuthal = grid.azimuthal_centers_rad[None, :, None]
    axial = grid.axial_centers_m[None, None, :]
    angular_log = np.zeros((1, grid.azimuthal_cell_count, 1))
    for order, (cosine_value, sine_value) in enumerate(
            zip(cosine, sine), start=1):
        angular_log += (
            cosine_value * np.cos(order * azimuthal)
            + sine_value * np.sin(order * azimuthal))
    raw = np.exp(
        -(grid.geometry.length_m - axial) / skin
        - 0.5 * ((radial - ring) / width) ** 2
        + angular_log)
    mean = float(
        np.sum(raw * grid.cell_volume_m3) / grid.geometry.volume_m3)
    if not math.isfinite(mean) or mean <= 0.0:
        raise RuntimeError("cylindrical source normalization failed")
    return raw / mean
