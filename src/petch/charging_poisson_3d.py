"""Compatible 3-D electrostatics for feature-scale dielectric charging.

The unknown is nodal voltage and the source is physical free charge in coulombs.  Trilinear Q1
elements assemble ``integral epsilon_r grad(Na).grad(Nb) dV`` on a uniform rectilinear grid, so
material-interface displacement continuity follows from the weak form rather than a fitted field law.
The resulting voltage array is directly consumable by :mod:`petch.boundary_transport_3d`.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix, vstack
from scipy.sparse.linalg import splu

from .charging_poisson import EPS0


_LOCAL_OFFSETS = np.asarray(tuple(product((0, 1), repeat=3)), dtype=int)
# Positive-weight Dunavant degree-four triangle rule.  Q1 basis functions restricted to an
# arbitrarily oriented planar triangle are degree at most three, so this is exact while retaining
# the physical nonnegativity of every nodal load.  The older four-point degree-three rule used a
# negative centroid weight; clipping endpoint roundoff to the owning Q1 cell could then turn a
# mathematically nonnegative integral negative at grid-aligned marching-cubes faces.
_TRIANGLE_BARYCENTRIC = np.asarray([
    [0.4459484909159648863, 0.4459484909159648863, 0.1081030181680702274],
    [0.4459484909159648863, 0.1081030181680702274, 0.4459484909159648863],
    [0.1081030181680702274, 0.4459484909159648863, 0.4459484909159648863],
    [0.0915762135097707435, 0.0915762135097707435, 0.8168475729804585131],
    [0.0915762135097707435, 0.8168475729804585131, 0.0915762135097707435],
    [0.8168475729804585131, 0.0915762135097707435, 0.0915762135097707435],
])
_TRIANGLE_QUADRATURE_WEIGHT = np.asarray([
    0.2233815896780114657, 0.2233815896780114657, 0.2233815896780114657,
    0.1099517436553218676, 0.1099517436553218676, 0.1099517436553218676,
])
_TRIANGLE_BARYCENTRIC /= np.sum(_TRIANGLE_BARYCENTRIC, axis=1)[:, None]
_TRIANGLE_QUADRATURE_WEIGHT /= np.sum(_TRIANGLE_QUADRATURE_WEIGHT)


def _spacing_vector(spacing_m):
    spacing = np.asarray(spacing_m, dtype=float)
    if spacing.ndim == 0:
        spacing = np.full(3, float(spacing))
    if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("spacing_m must be one positive finite value or a length-three vector")
    return spacing


def _local_q1_stiffness_3d(spacing_m):
    """Return the physical-cell trilinear Laplacian matrix for unit relative permittivity."""
    spacing = _spacing_vector(spacing_m)
    gauss = 0.5 + np.asarray((-1.0, 1.0)) / (2.0 * np.sqrt(3.0))
    local = np.zeros((8, 8))
    for coordinate in product(gauss, repeat=3):
        coordinate = np.asarray(coordinate)
        gradient = np.empty((8, 3))
        for node_index, offset in enumerate(_LOCAL_OFFSETS):
            for axis in range(3):
                derivative = 1.0 if offset[axis] else -1.0
                for other_axis in range(3):
                    if other_axis == axis:
                        continue
                    derivative *= (coordinate[other_axis] if offset[other_axis]
                                   else 1.0 - coordinate[other_axis])
                gradient[node_index, axis] = derivative / spacing[axis]
        # Two-point Gauss weights on [0,1] are 1/2 in each direction.
        local += gradient @ gradient.T * np.prod(spacing) / 8.0
    return 0.5 * (local + local.T)


def assemble_q1_stiffness_3d(epsilon_r, spacing_m):
    """Assemble ``integral epsilon_r grad(Na).grad(Nb) dV`` on a 3-D cell grid."""
    epsilon_r = np.asarray(epsilon_r, dtype=float)
    if (epsilon_r.ndim != 3 or epsilon_r.size == 0 or not np.all(np.isfinite(epsilon_r))
            or np.any(epsilon_r <= 0.0)):
        raise ValueError("epsilon_r must be a nonempty finite positive 3-D cell grid")
    local_unit = _local_q1_stiffness_3d(spacing_m)
    node_shape = tuple(np.asarray(epsilon_r.shape) + 1)
    rows = []; columns = []; values = []
    for cell in np.ndindex(epsilon_r.shape):
        node_coordinates = _LOCAL_OFFSETS + np.asarray(cell)
        nodes = np.ravel_multi_index(node_coordinates.T, node_shape)
        local = epsilon_r[cell] * local_unit
        rows.extend(np.repeat(nodes, 8))
        columns.extend(np.tile(nodes, 8))
        values.extend(local.ravel())
    node_count = int(np.prod(node_shape))
    return coo_matrix(
        (values, (rows, columns)), shape=(node_count, node_count)).tocsc()


@dataclass(frozen=True)
class PoissonDiagnostics3D:
    max_abs_free_charge_residual_c: float
    rms_free_charge_residual_c: float
    free_nodes: int
    electrostatic_energy_j: float
    specified_charge_c: float
    dirichlet_reaction_charge_c: float
    charge_balance_c: float
    maximum_floating_conductor_voltage_spread_v: float = 0.0
    floating_conductor_ids: tuple[int, ...] = ()
    floating_conductor_voltage_v: tuple[float, ...] = ()
    floating_conductor_charge_c: tuple[float, ...] = ()


class NodalPoissonSystem3D:
    """Reusable SI-unit Poisson factorization for one fixed 3-D material grid.

    Cell-centered ``epsilon_r`` may vary arbitrarily. ``charge_node_c`` is already projected onto the
    Q1 nodal basis. Boundaries not present in ``dirichlet_mask`` carry the natural zero-normal-
    displacement condition. At least one Dirichlet node is required to fix the voltage gauge.

    ``periodic_axes`` identifies the two endpoint planes of each selected axis *before* the linear
    system is factorized.  The public voltage grid retains both endpoint planes because trajectory
    interpolation needs the closed cell, but those entries are prolongations of one independent
    unknown and are therefore bitwise identical.  Full-grid charge on identified endpoints is
    summed onto that unknown; :meth:`canonicalize_charge` returns an equal-share full-grid
    representative without changing the physical periodic load.
    """

    def __init__(
            self, epsilon_r, spacing_m, dirichlet_mask, dirichlet_voltage=None, *,
            periodic_axes=(), floating_conductor_node_ids=None):
        self.epsilon_r = np.asarray(epsilon_r, dtype=float).copy()
        if self.epsilon_r.ndim != 3:
            raise ValueError("epsilon_r must be a 3-D cell grid")
        self.spacing_m = _spacing_vector(spacing_m).copy()
        self.shape = tuple(np.asarray(self.epsilon_r.shape) + 1)
        raw_periodic_axes = tuple(periodic_axes)
        if any(isinstance(axis, (bool, np.bool_)) or int(axis) != axis
               for axis in raw_periodic_axes):
            raise ValueError("periodic_axes must contain unique integer axes from 0, 1, 2")
        periodic_axes = tuple(sorted(int(axis) for axis in raw_periodic_axes))
        if (len(set(periodic_axes)) != len(periodic_axes)
                or any(axis < 0 or axis >= 3 for axis in periodic_axes)):
            raise ValueError("periodic_axes must contain unique axes from 0, 1, 2")
        self.periodic_axes = periodic_axes
        mask = np.asarray(dirichlet_mask, dtype=bool)
        if mask.shape != self.shape or not np.any(mask):
            raise ValueError("dirichlet_mask must match the nodal grid and fix at least one node")
        if dirichlet_voltage is None:
            voltage = np.zeros(self.shape)
        else:
            voltage = np.asarray(dirichlet_voltage, dtype=float)
            if voltage.shape != self.shape or not np.all(np.isfinite(voltage)):
                raise ValueError("dirichlet_voltage must be a finite nodal grid")
        self.dirichlet_mask = mask.copy()
        self.dirichlet_voltage = voltage.copy()
        if floating_conductor_node_ids is None:
            conductor_node_ids = np.zeros(self.shape, dtype=int)
        else:
            conductor_node_ids = np.asarray(
                floating_conductor_node_ids, dtype=int).copy()
            if (conductor_node_ids.shape != self.shape
                    or np.any(conductor_node_ids < 0)):
                raise ValueError(
                    "floating_conductor_node_ids must be a nonnegative nodal grid")
        if np.any((conductor_node_ids > 0) & self.dirichlet_mask):
            raise ValueError(
                "a floating-conductor node cannot also be a Dirichlet reservoir")
        self.floating_conductor_node_ids = conductor_node_ids
        self.stiffness = assemble_q1_stiffness_3d(self.epsilon_r, self.spacing_m)
        flat_mask = self.dirichlet_mask.ravel()
        if not periodic_axes:
            self.reduced_shape = self.shape
            self._full_to_reduced = np.arange(flat_mask.size, dtype=int)
            self._periodic_multiplicity = np.ones(flat_mask.size, dtype=int)
            self._nodal_charge_reduction = None
            self._nodal_voltage_prolongation = None
            self._operator = self.stiffness
            reduced_mask = flat_mask
            reduced_voltage = self.dirichlet_voltage.ravel()
            reduced_conductor_node_ids = conductor_node_ids.ravel().copy()
        else:
            reduced_shape = np.asarray(self.shape, dtype=int)
            reduced_shape[list(periodic_axes)] -= 1
            if np.any(reduced_shape < 1):
                raise ValueError("a periodic axis needs at least one physical cell")
            self.reduced_shape = tuple(int(value) for value in reduced_shape)
            coordinates = np.indices(self.shape).reshape(3, -1).T
            for axis in periodic_axes:
                coordinates[:, axis] %= reduced_shape[axis]
            self._full_to_reduced = np.ravel_multi_index(
                coordinates.T, self.reduced_shape)
            reduced_count = int(np.prod(self.reduced_shape))
            full_count = int(np.prod(self.shape))
            prolongation = coo_matrix(
                (np.ones(full_count),
                 (np.arange(full_count), self._full_to_reduced)),
                shape=(full_count, reduced_count)).tocsc()
            reduction = prolongation.T.tocsc()
            self._nodal_charge_reduction = reduction
            self._nodal_voltage_prolongation = prolongation
            self._periodic_multiplicity = np.bincount(
                self._full_to_reduced, minlength=reduced_count)
            self._operator = csc_matrix(reduction @ self.stiffness @ prolongation)
            fixed_count = np.bincount(
                self._full_to_reduced, weights=flat_mask.astype(int),
                minlength=reduced_count).astype(int)
            if np.any((fixed_count != 0) & (fixed_count != self._periodic_multiplicity)):
                raise ValueError(
                    "periodically identified Dirichlet masks must agree")
            reduced_mask = fixed_count > 0
            reduced_voltage = np.zeros(reduced_count)
            flat_voltage = self.dirichlet_voltage.ravel()
            for reduced_index in np.flatnonzero(reduced_mask):
                selected = flat_mask & (self._full_to_reduced == reduced_index)
                values = flat_voltage[selected]
                if values.size == 0 or not np.allclose(
                        values, values[0], rtol=0.0, atol=2e-14):
                    raise ValueError(
                        "periodically identified Dirichlet voltages must agree")
                reduced_voltage[reduced_index] = values[0]
            reduced_conductor_node_ids = np.zeros(reduced_count, dtype=int)
            flat_conductor_node_ids = conductor_node_ids.ravel()
            for reduced_index in range(reduced_count):
                values = np.unique(flat_conductor_node_ids[
                    self._full_to_reduced == reduced_index])
                if len(values) != 1:
                    raise ValueError(
                        "periodically identified floating-conductor ids must agree")
                reduced_conductor_node_ids[reduced_index] = int(values[0])
        self._reduced_dirichlet_mask = reduced_mask.copy()
        self._reduced_dirichlet_voltage = reduced_voltage.copy()
        self.fixed = np.flatnonzero(reduced_mask)
        self.free = np.flatnonzero(~reduced_mask)
        if self.free.size == 0:
            raise ValueError("dirichlet_mask must leave at least one free node")
        self._free_lookup = np.full(reduced_mask.size, -1, dtype=int)
        self._free_lookup[self.free] = np.arange(self.free.size)
        self._kff = csc_matrix(self._operator[self.free][:, self.free])
        self._kfc = csc_matrix(self._operator[self.free][:, self.fixed])
        self._factor = splu(self._kff)
        self._reduced_conductor_node_ids = reduced_conductor_node_ids
        self.floating_conductor_ids = tuple(
            int(value) for value in sorted(set(reduced_conductor_node_ids) - {0}))
        self._floating_conductor_representative_full_node = {}
        self._conductor_condensation = None
        self._conductor_factor = None
        if self.floating_conductor_ids:
            free_component = reduced_conductor_node_ids[self.free]
            if any(not np.any(free_component == component)
                   for component in self.floating_conductor_ids):
                raise ValueError("every floating conductor must contain a free Q1 node")
            nonconductor = np.flatnonzero(free_component == 0)
            group_column = {
                component: len(nonconductor) + index
                for index, component in enumerate(self.floating_conductor_ids)}
            column = np.empty(len(self.free), dtype=int)
            column[nonconductor] = np.arange(len(nonconductor))
            for component in self.floating_conductor_ids:
                column[free_component == component] = group_column[component]
            condensation = coo_matrix(
                (np.ones(len(self.free)), (np.arange(len(self.free)), column)),
                shape=(len(self.free),
                       len(nonconductor) + len(self.floating_conductor_ids))).tocsc()
            condensed_operator = csc_matrix(
                condensation.T @ self._kff @ condensation)
            self._conductor_condensation = condensation
            self._conductor_factor = splu(condensed_operator)
            flat_full_component = conductor_node_ids.ravel()
            for component in self.floating_conductor_ids:
                full_index = int(np.flatnonzero(
                    flat_full_component == component)[0])
                self._floating_conductor_representative_full_node[component] = (
                    tuple(int(value) for value in np.unravel_index(
                        full_index, self.shape)))

    @property
    def nodal_charge_reduction(self):
        """Sparse map from duplicated full nodes to independent periodic charge nodes."""
        return self._nodal_charge_reduction

    def reduce_charge(self, charge_node_c):
        """Sum full-grid nodal charge onto independent periodic degrees of freedom."""
        charge = np.asarray(charge_node_c, dtype=float)
        if charge.shape != self.shape or not np.all(np.isfinite(charge)):
            raise ValueError("charge_node_c must be a finite nodal grid")
        if self._nodal_charge_reduction is None:
            return charge.copy()
        reduced = np.asarray(self._nodal_charge_reduction @ charge.ravel()).ravel()
        return reduced.reshape(self.reduced_shape)

    def canonicalize_reduced_charge(self, reduced_charge_c):
        """Return the periodic/equipotential canonical full-grid charge representative.

        For each floating conductor, only the charge distribution is changed: its exact total
        free charge is preserved while the Q1 electrostatic solution is constrained to one
        equipotential.  Non-conductor nodal loads are preserved.  This is a linear electrostatic
        redistribution, not a nonlinear current-balance or root solve.
        """
        reduced = np.asarray(reduced_charge_c, dtype=float)
        if reduced.shape != self.reduced_shape or not np.all(np.isfinite(reduced)):
            raise ValueError("reduced_charge_c must match the independent periodic grid")
        reduced = self._redistribute_floating_conductor_charge(reduced.ravel())
        if self._nodal_charge_reduction is None:
            return reduced.reshape(self.shape).copy()
        flat = reduced[self._full_to_reduced]
        flat = flat / self._periodic_multiplicity[self._full_to_reduced]
        return flat.reshape(self.shape)

    def canonicalize_charge(self, charge_node_c):
        """Remove only periodic duplicates and internal floating-conductor redistribution."""
        return self.canonicalize_reduced_charge(self.reduce_charge(charge_node_c))

    @property
    def has_floating_conductors(self):
        return bool(self.floating_conductor_ids)

    def floating_conductor_representative_node(self, component_id):
        """Return one full-grid node used only to inject a component's conserved total current."""
        component = int(component_id)
        try:
            return self._floating_conductor_representative_full_node[component]
        except KeyError as error:
            raise KeyError(f"unknown floating conductor component {component}") from error

    def classify_surface_floating_conductors(
            self, face_centroids, face_gas_normals, *, grid_origin, grid_spacing):
        """Classify interface triangles by probing from the interface into the solid.

        Coordinates use the same declared mesh unit as ``grid_origin`` and ``grid_spacing``.
        Zero denotes an ordinary dielectric face.  The probe never infers electrical connectivity;
        it only transfers the caller-declared nodal component ids to the extracted surface.
        """
        centroid = np.asarray(face_centroids, dtype=float)
        normal = np.asarray(face_gas_normals, dtype=float)
        origin = np.asarray(grid_origin, dtype=float)
        spacing = np.asarray(grid_spacing, dtype=float)
        if spacing.ndim == 0:
            spacing = np.full(3, float(spacing))
        if (centroid.ndim != 2 or centroid.shape[1] != 3
                or normal.shape != centroid.shape or origin.shape != (3,)
                or spacing.shape != (3,) or np.any(~np.isfinite(centroid))
                or np.any(~np.isfinite(normal)) or np.any(~np.isfinite(origin))
                or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0)
                or not np.allclose(
                    np.linalg.norm(normal, axis=1), 1.0, rtol=0.0, atol=2e-6)):
            raise ValueError("invalid surface conductor-classification inputs")
        if not self.has_floating_conductors:
            return np.zeros(len(centroid), dtype=int)
        result = np.zeros(len(centroid), dtype=int)
        grid_maximum = origin + (np.asarray(self.shape) - 1) * spacing
        for fraction in (0.15, 0.35, 0.65, 0.95):
            selected = np.flatnonzero(result == 0)
            if not len(selected):
                break
            probe = centroid[selected] - fraction * float(np.min(spacing)) * normal[selected]
            probe = np.minimum(np.maximum(probe, origin), grid_maximum)
            index = np.rint((probe - origin) / spacing).astype(int)
            component = self.floating_conductor_node_ids[tuple(index.T)]
            found = component > 0
            result[selected[found]] = component[found]
        return result

    def _redistribute_floating_conductor_charge(self, reduced_charge_flat):
        charge = np.asarray(reduced_charge_flat, dtype=float).copy()
        if charge.shape != (int(np.prod(self.reduced_shape)),):
            raise ValueError("reduced charge has the wrong flattened shape")
        if not self.has_floating_conductors:
            return charge
        fixed_voltage = self._reduced_dirichlet_voltage[self.fixed]
        free_charge = charge[self.free]
        right_hand_side = self._conductor_condensation.T @ (
            free_charge / EPS0 - self._kfc @ fixed_voltage)
        condensed_voltage = self._conductor_factor.solve(
            np.asarray(right_hand_side).ravel())
        free_voltage = np.asarray(
            self._conductor_condensation @ condensed_voltage).ravel()
        redistributed = EPS0 * np.asarray(
            self._kff @ free_voltage + self._kfc @ fixed_voltage).ravel()
        free_component = self._reduced_conductor_node_ids[self.free]
        redistributed[free_component == 0] = free_charge[free_component == 0]
        for component in self.floating_conductor_ids:
            selected = free_component == component
            correction = float(
                np.sum(free_charge[selected]) - np.sum(redistributed[selected]))
            redistributed[np.flatnonzero(selected)[0]] += correction
        charge[self.free] = redistributed
        return charge

    def solve(self, charge_node_c=None):
        if charge_node_c is None:
            charge = np.zeros(self.shape)
        else:
            charge = np.asarray(charge_node_c, dtype=float)
            if charge.shape != self.shape or not np.all(np.isfinite(charge)):
                raise ValueError("charge_node_c must be a finite nodal grid")
        reduced_charge = self._redistribute_floating_conductor_charge(
            self.reduce_charge(charge).ravel())
        fixed_voltage = self._reduced_dirichlet_voltage[self.fixed]
        load = reduced_charge / EPS0
        right_hand_side = load[self.free] - self._kfc @ fixed_voltage
        reduced_voltage = self._reduced_dirichlet_voltage.copy()
        reduced_voltage[self.free] = self._factor.solve(right_hand_side)
        if self._nodal_voltage_prolongation is None:
            flat_voltage = reduced_voltage
        else:
            flat_voltage = np.asarray(
                self._nodal_voltage_prolongation @ reduced_voltage).ravel()
        voltage = flat_voltage.reshape(self.shape)
        residual_charge = EPS0 * (self._operator @ reduced_voltage) - reduced_charge
        free_residual = residual_charge[self.free]
        specified_charge = float(np.sum(reduced_charge))
        dirichlet_reaction = float(np.sum(residual_charge[self.fixed]))
        conductor_voltage = []
        conductor_charge = []
        maximum_conductor_spread = 0.0
        for component in self.floating_conductor_ids:
            selected = self._reduced_conductor_node_ids == component
            values = reduced_voltage[selected]
            conductor_voltage.append(float(np.mean(values)))
            conductor_charge.append(float(np.sum(reduced_charge[selected])))
            maximum_conductor_spread = max(
                maximum_conductor_spread,
                float(np.max(values) - np.min(values)))
        diagnostics = PoissonDiagnostics3D(
            max_abs_free_charge_residual_c=(float(np.max(np.abs(free_residual)))
                                            if free_residual.size else 0.0),
            rms_free_charge_residual_c=(float(np.sqrt(np.mean(free_residual ** 2)))
                                        if free_residual.size else 0.0),
            free_nodes=int(self.free.size),
            electrostatic_energy_j=(
                0.5 * EPS0 * float(reduced_voltage @ (
                    self._operator @ reduced_voltage))),
            specified_charge_c=specified_charge,
            dirichlet_reaction_charge_c=dirichlet_reaction,
            charge_balance_c=specified_charge + dirichlet_reaction,
            maximum_floating_conductor_voltage_spread_v=maximum_conductor_spread,
            floating_conductor_ids=self.floating_conductor_ids,
            floating_conductor_voltage_v=tuple(conductor_voltage),
            floating_conductor_charge_c=tuple(conductor_charge))
        return voltage, diagnostics

    def diagonal_capacitance(self, nodes):
        """Return exact diagonal nodal charge-to-voltage capacitances in farads."""
        response = self.voltage_response(nodes)
        diagonal = np.diag(response)
        if np.any(diagonal <= 0.0) or not np.all(np.isfinite(diagonal)):
            raise RuntimeError("Poisson response diagonal must be finite and positive")
        return 1.0 / diagonal

    def voltage_response(self, nodes):
        """Return the exact support-node voltage response in volts per coulomb.

        Column ``j`` is the voltage on every requested node after depositing one coulomb on node
        ``j`` with all Dirichlet voltages held fixed.  The dense matrix is a nonlinear-solver
        preconditioner for modest feature-surface supports; the volume Poisson operator remains sparse.
        """
        nodes = np.asarray(nodes, dtype=int)
        if nodes.ndim != 2 or nodes.shape[1] != 3 or nodes.size == 0:
            raise ValueError("nodes must have shape (n,3)")
        if np.any(nodes < 0) or np.any(nodes >= np.asarray(self.shape)):
            raise ValueError("response nodes lie outside the nodal grid")
        flat = np.ravel_multi_index(nodes.T, self.shape)
        reduced_flat = self._full_to_reduced[flat]
        if np.unique(reduced_flat).size != reduced_flat.size:
            raise ValueError("response nodes must be unique modulo periodic identification")
        free_rows = self._free_lookup[reduced_flat]
        if np.any(free_rows < 0):
            raise ValueError("response is defined only for free nodes")
        right_hand_side = np.zeros((self.free.size, len(flat)))
        right_hand_side[free_rows, np.arange(len(flat))] = 1.0 / EPS0
        response = self._factor.solve(right_hand_side)
        support_response = response[free_rows]
        if (not np.all(np.isfinite(support_response))
                or not np.allclose(support_response, support_response.T, rtol=2e-10, atol=0.0)):
            raise RuntimeError("Poisson support response must be finite and symmetric")
        return support_response


def lump_triangle_sheet_charge_3d(
        shape, vertices, faces, sigma_c_per_m2, *, grid_origin=(0.0, 0.0, 0.0),
        grid_spacing=1.0, coordinate_length_unit_m=1.0):
    """Project piecewise-constant triangle sheet charge onto a matching Q1 volume grid.

    Each triangle must lie inside one grid cell, as marching-cubes triangles do. A degree-three exact
    triangle rule integrates the trilinear basis; total charge and its first spatial moments are thus
    conserved without nearest-node assignment. Coordinates may use any declared length unit.
    """
    shape = tuple(int(value) for value in shape)
    input_vertices = np.asarray(vertices)
    vertices = np.asarray(vertices, dtype=float); faces = np.asarray(faces, dtype=int)
    sigma = np.asarray(sigma_c_per_m2, dtype=float)
    origin = np.asarray(grid_origin, dtype=float)
    spacing = np.asarray(grid_spacing, dtype=float)
    if spacing.ndim == 0:
        spacing = np.full(3, float(spacing))
    if (len(shape) != 3 or min(shape) < 2 or vertices.ndim != 2 or vertices.shape[1] != 3
            or faces.ndim != 2 or faces.shape[1] != 3 or sigma.shape != (len(faces),)
            or origin.shape != (3,) or spacing.shape != (3,)
            or np.any(~np.isfinite(vertices)) or np.any(~np.isfinite(sigma))
            or np.any(~np.isfinite(origin)) or np.any(~np.isfinite(spacing))
            or np.any(spacing <= 0.0) or np.any(faces < 0) or np.any(faces >= len(vertices))
            or not np.isfinite(coordinate_length_unit_m) or coordinate_length_unit_m <= 0.0):
        raise ValueError("invalid triangle sheet-charge projection inputs")
    cell_shape = np.asarray(shape) - 1
    normalized_vertices = (vertices - origin) / spacing
    # Marching-cubes vertices are float32.  A vertex on a physical endpoint can therefore land a few
    # ulps beyond the matching nodal-grid endpoint after division by ``spacing`` (for example,
    # float32(0.3) / 0.01 = 30.00000119).  Scale the admission tolerance by the source precision and
    # grid extent, then clamp only the admitted roundoff so downstream cell selection remains exact.
    source_epsilon = (np.finfo(input_vertices.dtype).eps
                      if np.issubdtype(input_vertices.dtype, np.floating)
                      else np.finfo(float).eps)
    # The common driver validates/casts mesh coordinates before they reach this routine. Preserve
    # the marching-cubes precision contract even after that harmless float32 -> float64 promotion:
    # an exact float32 round trip proves that the larger source ulp is still the relevant bound.
    if (np.issubdtype(input_vertices.dtype, np.floating)
            and input_vertices.dtype.itemsize > np.dtype(np.float32).itemsize
            and np.array_equal(
                input_vertices, input_vertices.astype(np.float32).astype(input_vertices.dtype))):
        source_epsilon = max(source_epsilon, np.finfo(np.float32).eps)
    tolerance = max(1e-10, 8.0 * source_epsilon * max(1.0, float(np.max(cell_shape))))
    if (np.any(normalized_vertices < -tolerance)
            or np.any(normalized_vertices > cell_shape + tolerance)):
        raise ValueError("triangle vertices lie outside the nodal grid")
    normalized_vertices = np.clip(normalized_vertices, 0.0, cell_shape)

    charge = np.zeros(shape)
    for face_index, face in enumerate(faces):
        triangle = vertices[face]
        normalized_triangle = normalized_vertices[face]
        lower = np.floor(np.min(normalized_triangle, axis=0) + tolerance).astype(int)
        lower = np.minimum(lower, cell_shape - 1)
        if np.any(np.max(normalized_triangle, axis=0) > lower + 1.0 + tolerance):
            raise ValueError("each triangle must lie within one potential-grid cell")
        physical_edge_a = (triangle[1] - triangle[0]) * float(coordinate_length_unit_m)
        physical_edge_b = (triangle[2] - triangle[0]) * float(coordinate_length_unit_m)
        area_m2 = 0.5 * float(np.linalg.norm(np.cross(physical_edge_a, physical_edge_b)))
        if area_m2 <= 0.0:
            raise ValueError("sheet-charge triangles must have positive area")
        quadrature_point = _TRIANGLE_BARYCENTRIC @ normalized_triangle
        for point, weight in zip(
                quadrature_point, _TRIANGLE_QUADRATURE_WEIGHT):
            fraction = np.clip(point - lower, 0.0, 1.0)
            nodal_weight = np.prod(np.where(
                _LOCAL_OFFSETS == 1, fraction[None, :], 1.0 - fraction[None, :]), axis=1)
            contribution = float(sigma[face_index]) * area_m2 * float(weight) * nodal_weight
            nodes = _LOCAL_OFFSETS + lower
            np.add.at(charge, tuple(nodes.T), contribution)
    return charge


def lump_mixed_surface_density_3d(
        poisson_system, vertices, faces, surface_density, face_conductor_id, *,
        grid_origin=(0.0, 0.0, 0.0), grid_spacing=1.0,
        coordinate_length_unit_m=1.0, canonicalize=True):
    """Conservatively couple dielectric faces and pooled floating conductors to Q1.

    ``surface_density`` may carry any surface-integrated rate or inventory unit (for example
    C/m2 or A/m2). Ordinary dielectric triangles use the exact Q1 sheet coupling. Every triangle
    tagged with a positive conductor id instead contributes its *complete integrated amount* to
    that electrical component; :meth:`NodalPoissonSystem3D.canonicalize_charge` then computes the
    unique equipotential surface redistribution while preserving the component total.  With
    ``canonicalize=False`` the component total remains at its representative node.  That
    nonnegative representation is the correct one for separate ion/electron arrival diagnostics;
    signed induced redistribution belongs only to the electrostatic state update.
    """
    if not isinstance(poisson_system, NodalPoissonSystem3D):
        raise TypeError("poisson_system must be NodalPoissonSystem3D")
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    density = np.asarray(surface_density, dtype=float)
    component = np.asarray(face_conductor_id, dtype=int)
    if (faces.ndim != 2 or faces.shape[1] != 3
            or density.shape != (len(faces),) or component.shape != density.shape
            or np.any(~np.isfinite(density)) or np.any(component < 0)
            or not isinstance(canonicalize, (bool, np.bool_))
            or not set(np.unique(component)).issubset(
                {0, *poisson_system.floating_conductor_ids})):
        raise ValueError("invalid mixed surface-density or conductor classification")
    dielectric_density = np.where(component == 0, density, 0.0)
    coupled = lump_triangle_sheet_charge_3d(
        poisson_system.shape, vertices, faces, dielectric_density,
        grid_origin=grid_origin, grid_spacing=grid_spacing,
        coordinate_length_unit_m=coordinate_length_unit_m)
    triangle = vertices[faces]
    physical_area = 0.5 * np.linalg.norm(
        np.cross(triangle[:, 1] - triangle[:, 0],
                 triangle[:, 2] - triangle[:, 0]), axis=1)
    physical_area *= float(coordinate_length_unit_m) ** 2
    if np.any(~np.isfinite(physical_area)) or np.any(physical_area <= 0.0):
        raise ValueError("mixed surface-density triangles must have positive area")
    for conductor_id in poisson_system.floating_conductor_ids:
        selected = component == conductor_id
        if np.any(selected):
            node = poisson_system.floating_conductor_representative_node(conductor_id)
            coupled[node] += float(np.dot(density[selected], physical_area[selected]))
    return poisson_system.canonicalize_charge(coupled) if canonicalize else coupled


def triangle_sheet_face_charge_coupling_3d(
        shape, vertices, faces, *, grid_origin=(0.0, 0.0, 0.0), grid_spacing=1.0,
        coordinate_length_unit_m=1.0):
    """Return the conservative map from triangle charge to Q1 nodal charge.

    Column ``f`` contains the Q1 nodal load produced by one coulomb distributed uniformly over
    triangle ``f``.  Consequently every column sums to one.  The map exposes an important
    compatibility property that a direct face-to-node projection can otherwise hide: a P0
    triangle charge space can contain modes in the exact null space of the Q1 field solve.
    Carrying those modes as authoritative charge lets them accumulate without electrostatic
    feedback.  :class:`CompatibleQ1SurfaceChargeProjector3D` removes only that unresolvable
    component while preserving the complete nodal load and global charge.
    """
    shape = tuple(int(value) for value in shape)
    input_vertices = np.asarray(vertices)
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    origin = np.asarray(grid_origin, dtype=float)
    spacing = np.asarray(grid_spacing, dtype=float)
    if spacing.ndim == 0:
        spacing = np.full(3, float(spacing))
    if (len(shape) != 3 or min(shape) < 2 or vertices.ndim != 2 or vertices.shape[1] != 3
            or faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0
            or origin.shape != (3,) or spacing.shape != (3,)
            or np.any(~np.isfinite(vertices)) or np.any(~np.isfinite(origin))
            or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0)
            or np.any(faces < 0) or np.any(faces >= len(vertices))
            or not np.isfinite(coordinate_length_unit_m)
            or coordinate_length_unit_m <= 0.0):
        raise ValueError("invalid triangle sheet-charge coupling inputs")

    cell_shape = np.asarray(shape) - 1
    normalized_vertices = (vertices - origin) / spacing
    source_epsilon = (np.finfo(input_vertices.dtype).eps
                      if np.issubdtype(input_vertices.dtype, np.floating)
                      else np.finfo(float).eps)
    if (np.issubdtype(input_vertices.dtype, np.floating)
            and input_vertices.dtype.itemsize > np.dtype(np.float32).itemsize
            and np.array_equal(
                input_vertices, input_vertices.astype(np.float32).astype(input_vertices.dtype))):
        source_epsilon = max(source_epsilon, np.finfo(np.float32).eps)
    tolerance = max(1e-10, 8.0 * source_epsilon * max(1.0, float(np.max(cell_shape))))
    if (np.any(normalized_vertices < -tolerance)
            or np.any(normalized_vertices > cell_shape + tolerance)):
        raise ValueError("triangle vertices lie outside the nodal grid")
    normalized_vertices = np.clip(normalized_vertices, 0.0, cell_shape)

    rows = []
    columns = []
    values = []
    physical_area_m2 = np.empty(len(faces))
    for face_index, face in enumerate(faces):
        triangle = vertices[face]
        normalized_triangle = normalized_vertices[face]
        lower = np.floor(np.min(normalized_triangle, axis=0) + tolerance).astype(int)
        lower = np.minimum(lower, cell_shape - 1)
        if np.any(np.max(normalized_triangle, axis=0) > lower + 1.0 + tolerance):
            raise ValueError("each triangle must lie within one potential-grid cell")
        physical_edge_a = (triangle[1] - triangle[0]) * float(coordinate_length_unit_m)
        physical_edge_b = (triangle[2] - triangle[0]) * float(coordinate_length_unit_m)
        physical_area_m2[face_index] = 0.5 * float(
            np.linalg.norm(np.cross(physical_edge_a, physical_edge_b)))
        if physical_area_m2[face_index] <= 0.0:
            raise ValueError("sheet-charge triangles must have positive area")
        quadrature_point = _TRIANGLE_BARYCENTRIC @ normalized_triangle
        for point, weight in zip(
                quadrature_point, _TRIANGLE_QUADRATURE_WEIGHT):
            fraction = np.clip(point - lower, 0.0, 1.0)
            nodal_weight = np.prod(np.where(
                _LOCAL_OFFSETS == 1, fraction[None, :], 1.0 - fraction[None, :]), axis=1)
            nodes = _LOCAL_OFFSETS + lower
            rows.extend(np.ravel_multi_index(nodes.T, shape))
            columns.extend(np.full(8, face_index, dtype=int))
            values.extend(float(weight) * nodal_weight)
    coupling = coo_matrix(
        (values, (rows, columns)), shape=(int(np.prod(shape)), len(faces))).tocsc()
    column_sum = np.asarray(coupling.sum(axis=0)).ravel()
    if not np.allclose(column_sum, 1.0, rtol=0.0, atol=2e-14):
        raise RuntimeError("Q1 triangle charge coupling lost partition of unity")
    return coupling, physical_area_m2


class CompatibleQ1SurfaceChargeProjector3D:
    """Represent surface charge in the part of P0 face space visible to Q1 Poisson.

    The electrostatic source actually consumed by the engine is the Q1 nodal load.  If the P0
    triangle space has more independent modes than that nodal trace, retaining an arbitrary face
    representative creates exact field-invisible modes.  This projector uses the unique face
    charge with minimum area-weighted charge-density norm for a given nodal load.  Projection:

    * preserves the Q1 nodal charge load to roundoff;
    * preserves exact global charge because the coupling is a partition of unity;
    * is idempotent; and
    * removes no field-resolved information.

    The dense SVD is constructed once per modest feature surface.  It is deliberately diagnostic
    and explicit: rank, nullity, and the smallest resolved singular value are available to run
    manifests instead of allowing a hidden discretization null space.
    """

    def __init__(self, coupling, physical_face_area_m2, *, relative_rank_tolerance=None):
        coupling = csc_matrix(coupling, dtype=float)
        area = np.asarray(physical_face_area_m2, dtype=float)
        if (coupling.ndim != 2 or coupling.shape[1] == 0
                or area.shape != (coupling.shape[1],)
                or np.any(~np.isfinite(area)) or np.any(area <= 0.0)):
            raise ValueError("invalid compatible-Q1 surface-charge projector inputs")
        dense = coupling.toarray()
        if np.any(~np.isfinite(dense)):
            raise ValueError("surface-charge coupling must be finite")
        square_root_area = np.sqrt(area)
        scaled = dense * square_root_area[None, :]
        left, singular, right_transpose = np.linalg.svd(scaled, full_matrices=False)
        if singular.size == 0 or singular[0] <= 0.0:
            raise ValueError("surface-charge coupling has no resolved Q1 mode")
        if relative_rank_tolerance is None:
            relative_rank_tolerance = max(scaled.shape) * np.finfo(float).eps
        if (not np.isfinite(relative_rank_tolerance)
                or not 0.0 < relative_rank_tolerance < 1.0):
            raise ValueError("relative_rank_tolerance must lie strictly between zero and one")
        threshold = float(relative_rank_tolerance) * float(singular[0])
        rank = int(np.sum(singular > threshold))
        if rank == 0:
            raise ValueError("surface-charge coupling has no resolved Q1 mode")

        self.coupling = coupling
        self.physical_face_area_m2 = area.copy()
        self.rank = rank
        self.nullity = int(coupling.shape[1] - rank)
        self.relative_rank_tolerance = float(relative_rank_tolerance)
        self.singular_values = singular.copy()
        self.condition_number = float(singular[0] / singular[rank - 1])
        self._node_shape = None
        self._sqrt_area = square_root_area
        self._resolved_right = right_transpose[:rank].T.copy()
        self._resolved_left = left[:, :rank].copy()
        self._resolved_singular = singular[:rank].copy()
        self._field_coupling = coupling
        self._inventory_constraint_count = 0

    @classmethod
    def from_triangles(
            cls, shape, vertices, faces, *, grid_origin=(0.0, 0.0, 0.0),
            grid_spacing=1.0, coordinate_length_unit_m=1.0,
            relative_rank_tolerance=None):
        coupling, area = triangle_sheet_face_charge_coupling_3d(
            shape, vertices, faces, grid_origin=grid_origin, grid_spacing=grid_spacing,
            coordinate_length_unit_m=coordinate_length_unit_m)
        result = cls(
            coupling, area, relative_rank_tolerance=relative_rank_tolerance)
        result._node_shape = tuple(int(value) for value in shape)
        return result

    @classmethod
    def from_poisson_system(
            cls, poisson_system, vertices, faces, *, grid_origin=(0.0, 0.0, 0.0),
            grid_spacing=1.0, coordinate_length_unit_m=1.0,
            relative_rank_tolerance=None):
        """Construct the face projector for the Poisson operator's actual charge space.

        A periodic Q1 system first sums duplicated endpoint loads onto independent periodic
        unknowns.  Building compatibility against the unreduced full grid would therefore retain
        face modes that are invisible to the *actual* field operator.  This constructor composes
        the conservative triangle coupling with that exact reduction before computing its rank and
        minimum-norm representative.  Nonperiodic systems reduce to :meth:`from_triangles`.
        """
        if not isinstance(poisson_system, NodalPoissonSystem3D):
            raise TypeError("poisson_system must be NodalPoissonSystem3D")
        if poisson_system.has_floating_conductors:
            raise ValueError(
                "the minimum-norm face projector does not yet define a unique conductor-face "
                "representative; use the mixed surface-density coupling")
        coupling, area = triangle_sheet_face_charge_coupling_3d(
            poisson_system.shape, vertices, faces, grid_origin=grid_origin,
            grid_spacing=grid_spacing,
            coordinate_length_unit_m=coordinate_length_unit_m)
        if poisson_system.nodal_charge_reduction is not None:
            coupling = csc_matrix(poisson_system.nodal_charge_reduction @ coupling)
        result = cls(
            coupling, area, relative_rank_tolerance=relative_rank_tolerance)
        result._node_shape = tuple(int(value) for value in poisson_system.reduced_shape)
        return result

    @classmethod
    def from_mixed_poisson_system(
            cls, poisson_system, vertices, faces, face_conductor_id, *,
            grid_origin=(0.0, 0.0, 0.0), grid_spacing=1.0,
            coordinate_length_unit_m=1.0, relative_rank_tolerance=None):
        """Construct the compatible state for dielectric faces and floating conductors.

        Dielectric triangles retain the same Q1-visible minimum-density representative used by
        :meth:`from_poisson_system`.  Every face on one floating conductor couples only through
        that conductor's integrated free charge, matching :func:`lump_mixed_surface_density_3d`.
        Extra linear constraints preserve each conductor inventory separately, so projection may
        redistribute charge *within* an equipotential component but can never move charge between
        a dielectric and a conductor or between distinct conductors.

        In area-weighted density coordinates, the minimum-norm representative of one fixed
        conductor total is uniform surface density.  This is not a conductivity model added to a
        dielectric: it is the unique stored representative of the equipotential conductor already
        declared in the Poisson operator.
        """
        if not isinstance(poisson_system, NodalPoissonSystem3D):
            raise TypeError("poisson_system must be NodalPoissonSystem3D")
        if not poisson_system.has_floating_conductors:
            raise ValueError("mixed compatible coupling requires floating conductors")
        vertices = np.asarray(vertices, dtype=float)
        faces = np.asarray(faces, dtype=int)
        component = np.asarray(face_conductor_id, dtype=int)
        if (component.shape != (len(faces),)
                or np.any(component < 0)
                or not set(np.unique(component)).issubset(
                    {0, *poisson_system.floating_conductor_ids})):
            raise ValueError("face_conductor_id does not match the Poisson components")

        triangle_coupling, area = triangle_sheet_face_charge_coupling_3d(
            poisson_system.shape, vertices, faces, grid_origin=grid_origin,
            grid_spacing=grid_spacing,
            coordinate_length_unit_m=coordinate_length_unit_m)
        triangle_coo = triangle_coupling.tocoo()
        keep = component[triangle_coo.col] == 0
        rows = triangle_coo.row[keep].tolist()
        columns = triangle_coo.col[keep].tolist()
        values = triangle_coo.data[keep].tolist()
        for face_index in np.flatnonzero(component > 0):
            representative = poisson_system.floating_conductor_representative_node(
                component[face_index])
            rows.append(int(np.ravel_multi_index(
                representative, poisson_system.shape)))
            columns.append(int(face_index))
            values.append(1.0)
        field_coupling = coo_matrix(
            (values, (rows, columns)),
            shape=(int(np.prod(poisson_system.shape)), len(faces))).tocsc()
        if poisson_system.nodal_charge_reduction is not None:
            field_coupling = csc_matrix(
                poisson_system.nodal_charge_reduction @ field_coupling)
        # Rows untouched by this particular surface cannot change its row space. Removing them
        # keeps the one-time dense SVD bounded for realistic feature grids.
        active_field_rows = np.flatnonzero(field_coupling.getnnz(axis=1))
        field_coupling = field_coupling[active_field_rows]

        inventory_rows = []
        for conductor_id in poisson_system.floating_conductor_ids:
            inventory_rows.append(
                csc_matrix((component == conductor_id).astype(float)[None, :]))
        constrained_coupling = vstack(
            (field_coupling, *inventory_rows), format="csc")
        result = cls(
            constrained_coupling, area,
            relative_rank_tolerance=relative_rank_tolerance)
        result._field_coupling = field_coupling
        result._inventory_constraint_count = len(inventory_rows)
        result._mixed_conductor_ids = poisson_system.floating_conductor_ids
        result._node_shape = None
        return result

    def node_charge_from_face_charge(self, face_charge_c):
        face_charge = np.asarray(face_charge_c, dtype=float)
        if (face_charge.shape != (self.coupling.shape[1],)
                or np.any(~np.isfinite(face_charge))):
            raise ValueError("face_charge_c must be one finite value per triangle")
        node = np.asarray(self._field_coupling @ face_charge).ravel()
        return node if self._node_shape is None else node.reshape(self._node_shape)

    def project_face_charge(self, face_charge_c):
        """Return the area-weighted minimum-norm face representative of the same nodal load."""
        face_charge = np.asarray(face_charge_c, dtype=float)
        if (face_charge.shape != (self.coupling.shape[1],)
                or np.any(~np.isfinite(face_charge))):
            raise ValueError("face_charge_c must be one finite value per triangle")
        if self.nullity == 0:
            return face_charge.copy()
        density_coordinate = face_charge / self._sqrt_area
        resolved = self._resolved_right @ (self._resolved_right.T @ density_coordinate)
        return self._sqrt_area * resolved

    def face_charge_from_node_charge(self, charge_node_c, *, relative_residual_tolerance=1e-11):
        """Reconstruct the compatible face representative for an in-range Q1 nodal load."""
        if self._inventory_constraint_count:
            raise ValueError(
                "mixed surface reconstruction also requires floating-conductor inventories")
        node = np.asarray(charge_node_c, dtype=float)
        if self._node_shape is not None and node.shape == self._node_shape:
            node = node.ravel()
        if (node.shape != (self.coupling.shape[0],) or np.any(~np.isfinite(node))
                or not np.isfinite(relative_residual_tolerance)
                or relative_residual_tolerance <= 0.0):
            raise ValueError("invalid nodal charge reconstruction inputs")
        coefficient = (self._resolved_left.T @ node) / self._resolved_singular
        face_charge = self._sqrt_area * (self._resolved_right @ coefficient)
        reconstructed = np.asarray(self.coupling @ face_charge).ravel()
        scale = max(float(np.linalg.norm(node)), np.finfo(float).tiny)
        if np.linalg.norm(reconstructed - node) / scale > relative_residual_tolerance:
            raise ValueError("nodal charge contains a component outside the surface coupling range")
        return face_charge

    def unresolved_fraction(self, face_charge_c):
        """Return the fraction of area-weighted charge-density norm invisible to Q1 Poisson."""
        face_charge = np.asarray(face_charge_c, dtype=float)
        projected = self.project_face_charge(face_charge)
        coordinate = face_charge / self._sqrt_area
        scale = np.linalg.norm(coordinate)
        if scale == 0.0:
            return 0.0
        return float(np.linalg.norm((face_charge - projected) / self._sqrt_area) / scale)

    def unresolved_linear_functional_fraction(self, face_charge_functional):
        """Return how strongly a face-charge statistic responds to Q1-null modes.

        ``face_charge_functional`` defines a scalar statistic ``l @ q`` on integrated face
        charges.  A value of zero means that statistic is completely determined by the Q1 nodal
        load.  A nonzero value means two face states with identical electrostatic fields can report
        different values of the statistic.  Patch-summed current balance is such a functional and
        must be checked before it is used as a convergence gate.
        """
        functional = np.asarray(face_charge_functional, dtype=float)
        if (functional.shape != (self.coupling.shape[1],)
                or np.any(~np.isfinite(functional))):
            raise ValueError("face_charge_functional must be one finite value per triangle")
        density_dual = functional * self._sqrt_area
        scale = np.linalg.norm(density_dual)
        if scale == 0.0 or self.nullity == 0:
            return 0.0
        resolved = self._resolved_right @ (self._resolved_right.T @ density_dual)
        return float(np.linalg.norm(density_dual - resolved) / scale)
