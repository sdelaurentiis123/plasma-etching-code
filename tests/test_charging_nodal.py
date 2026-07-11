import numpy as np
from scipy.stats import gamma as gamma_dist, norm, qmc

from petch.charging_backward import (
    backward_electron_floor_liouville,
    backward_electron_gather,
    backward_ion_gather,
)
from petch.charging_nodal import nodal_domain, nodal_field_at, solve_nodal_laplace, trace_nodal


def _flat_gap(nx=12, gas_rows=10, voltage=7.0):
    solid = np.zeros((nx, gas_rows + 1), dtype=bool)
    solid[:, gas_rows] = True
    surface = np.zeros_like(solid, dtype=float)
    surface[:, gas_rows] = voltage
    return solid, surface


def test_nodal_domain_places_voltage_on_exact_material_face():
    solid, surface = _flat_gap(voltage=7.0)
    active, fixed, value = nodal_domain(solid, surface)
    assert np.all(active[:, 10])
    assert np.all(fixed[:, 10])
    assert np.all(value[:, 10] == 7.0)
    assert np.all(fixed[:, 0])
    assert np.all(value[:, 0] == 0.0)


def test_nodal_laplace_reproduces_linear_parallel_plate_solution():
    solid, surface = _flat_gap(nx=12, gas_rows=10, voltage=7.0)
    V, diag = solve_nodal_laplace(solid, surface, sweeps=2000, omega=1.7, tolerance=1e-11)
    expected = 7.0 * np.arange(11) / 10.0
    assert np.allclose(V[:, :11], expected[None, :], atol=2e-9)
    assert diag["max_abs"] < 1e-11
    ex, ez = nodal_field_at(V, 5.25, 4.75)
    assert abs(ex) < 1e-10
    assert np.isclose(ez, -0.7, atol=2e-9)


def test_nodal_laplace_constant_boundary_is_constant():
    solid, surface = _flat_gap(nx=8, gas_rows=6, voltage=0.0)
    V, diag = solve_nodal_laplace(solid, surface, sweeps=20)
    assert np.all(V == 0.0)
    assert diag["max_abs"] == 0.0


def test_nodal_tracer_cannot_tunnel_through_one_cell_wall():
    nx, nz = 12, 8
    solid = np.zeros((nx, nz), dtype=bool)
    solid[6, :] = True
    V = np.zeros((nx + 1, nz + 1))
    x = np.array([2.25]); z = np.array([3.4])
    vx = np.array([30.0]); vz = np.array([0.0])
    hit_x, hit_z, *_ = trace_nodal(V, solid, x, z, vx, vz, 1.0, nx, nz, 500, 0.4, 0.1)
    assert hit_x[0] == 6
    assert hit_z[0] == 3


def test_nodal_tracer_checks_second_crossed_face_on_diagonal_step():
    nx, nz = 10, 10
    solid = np.zeros((nx, nz), dtype=bool)
    solid[:, 5] = True
    V = np.zeros((nx + 1, nz + 1))
    # Cross x=4 first into gas, then z=5 into the solid floor during the same step.
    x = np.array([3.95]); z = np.array([4.8])
    vx = np.array([1.0]); vz = np.array([3.0])
    hit_x, hit_z, *_ = trace_nodal(V, solid, x, z, vx, vz, 1.0, nx, nz, 20, 0.4, 0.1)
    assert hit_x[0] == 4
    assert hit_z[0] == 5


def test_nodal_tracer_lateral_reflection_keeps_remaining_step():
    nx, nz = 8, 6
    solid = np.zeros((nx, nz), dtype=bool)
    # A vertical detector catches the particle only if reflection preserves the remaining distance.
    solid[1, :] = True
    V = np.zeros((nx + 1, nz + 1))
    x = np.array([0.1]); z = np.array([2.5])
    vx = np.array([-4.0]); vz = np.array([0.0])
    hit_x, hit_z, *_ = trace_nodal(V, solid, x, z, vx, vz, 1.0, nx, nz, 20, 0.4, 0.1,
                                   fixed_dt=0.2)
    assert hit_x[0] == 1
    assert hit_z[0] == 2


def test_nodal_tracer_uniform_field_energy_error_is_small():
    nx, nz = 8, 16
    solid = np.zeros((nx, nz), dtype=bool)
    solid[:, 12:] = True
    V = np.broadcast_to(-np.arange(nz + 1, dtype=float), (nx + 1, nz + 1)).copy()
    x = np.array([4.0]); z = np.array([1.0])
    vx = np.array([0.0]); vz = np.array([1.0])
    hit_x, hit_z, impact, *_ = trace_nodal(V, solid, x, z, vx, vz, 1.0, nx, nz, 1000, 0.15, 0.1)
    assert hit_x[0] >= 0 and hit_z[0] == 12
    # H=v^2+qV: initial H=0, so K at the exact z=12 face is 12.
    assert np.isclose(impact[0], 12.0, atol=2e-3)


def _nodal_trench(floor_potential):
    nx, nz = 96, 72
    left, right, floor_z = 12, 82, 64
    solid = np.zeros((nx, nz), dtype=bool)
    solid[left, 1:floor_z + 1] = True
    solid[right, 1:floor_z + 1] = True
    solid[left:right + 1, floor_z] = True
    nodal_V = np.broadcast_to(
        floor_potential * np.arange(nz + 1) / floor_z, (nx + 1, nz + 1)).copy()
    surface = np.zeros((nx, nz))
    surface[:, floor_z] = floor_potential
    cells = [(x, floor_z) for x in range(left + 1, right)]
    normals = [(0.0, -1.0)] * len(cells)
    return solid, nodal_V, surface, cells, normals, left, right, floor_z


def test_nodal_electron_forward_backward_reciprocity_in_linear_field():
    solid, V, surface, cells, normals, left, right, floor_z = _nodal_trench(7.0)
    nx, nz = solid.shape; zero = np.zeros_like(surface)
    backward = backward_electron_gather(
        solid, zero, zero, surface, cells, normals, n_log2=10, n_scramble=3,
        seed=71, nodal_potential=V).mean()
    u = qmc.Sobol(d=4, scramble=True, seed=73).random_base2(15)
    energy = gamma_dist.ppf(u[:, 0], a=2.0, scale=4.0)
    ct = np.sqrt(u[:, 1])
    vx = np.sqrt(energy) * np.sqrt(1.0 - ct * ct) * np.cos(2.0 * np.pi * u[:, 2])
    vz = np.sqrt(energy) * ct
    x = left + 1.0 + u[:, 3] * (right - left - 1.0)
    z = np.full_like(x, 1e-3)
    hx, hz, *_ = trace_nodal(V, solid, x, z, vx, vz, -1.0, nx, nz, 200 * nz, 0.15, 0.1)
    forward = np.mean((hz == floor_z) & (hx > left) & (hx < right))
    assert np.isclose(backward, forward, rtol=0.04, atol=0.006), (backward, forward)
    for shifted_fraction in (0.2, 0.8):
        liouville = backward_electron_floor_liouville(
            solid, V, surface, cells, n_log2=10, n_scramble=3, seed=89,
            shifted_fraction=shifted_fraction).mean()
        assert np.isclose(liouville, forward, rtol=0.04, atol=0.006), (liouville, forward)

def test_nodal_ion_forward_backward_reciprocity_in_linear_field():
    solid, V, surface, cells, normals, left, right, floor_z = _nodal_trench(5.0)
    nx, nz = solid.shape; zero = np.zeros_like(surface)
    backward = backward_ion_gather(
        solid, zero, zero, surface, cells, normals, n_log2=10, n_scramble=3,
        seed=79, nodal_potential=V).mean()
    u = qmc.Sobol(d=3, scramble=True, seed=83).random_base2(15)
    phase = 2.0 * np.pi * u[:, 0]
    vx = np.sqrt(0.25) * norm.ppf(np.clip(u[:, 1], 1e-6, 1.0 - 1e-6))
    vz = np.sqrt(2.0 + 37.0 + 30.0 * np.sin(phase))
    x = left + 1.0 + u[:, 2] * (right - left - 1.0)
    z = np.full_like(x, 1e-3)
    hx, hz, *_ = trace_nodal(V, solid, x, z, vx, vz, 1.0, nx, nz, 200 * nz, 0.15, 0.1)
    forward = np.mean((hz == floor_z) & (hx > left) & (hx < right))
    assert np.isclose(backward, forward, rtol=0.04, atol=0.006), (backward, forward)
