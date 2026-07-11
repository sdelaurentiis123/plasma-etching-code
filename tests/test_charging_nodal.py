import numpy as np

from petch.charging_nodal import nodal_domain, nodal_field_at, solve_nodal_laplace


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
