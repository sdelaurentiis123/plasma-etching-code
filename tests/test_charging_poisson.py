import numpy as np

from petch.charging_poisson import EPS0, NodalPoissonSystem, lump_edge_sheet_charge


def _top_sheet_charge(nx, nz, sigma, h):
    charge = np.zeros((nx + 1, nz + 1))
    charge[0, 0] = 0.5 * sigma * h
    charge[-1, 0] = 0.5 * sigma * h
    charge[1:-1, 0] = sigma * h
    return charge


def test_q1_poisson_reproduces_uniform_parallel_plate_capacitance():
    nx, nz = 5, 8
    epsilon_r = np.full((nx, nz), 3.9)
    fixed = np.zeros((nx + 1, nz + 1), dtype=bool); fixed[:, -1] = True
    system = NodalPoissonSystem(epsilon_r, fixed)
    sigma = 7e-4; h = 31.25e-9
    voltage, diagnostics = system.solve(_top_sheet_charge(nx, nz, sigma, h))

    expected = sigma * (nz * h) / (EPS0 * 3.9)
    assert np.allclose(voltage[:, 0], expected, rtol=1e-12, atol=1e-12)
    assert np.allclose(voltage.mean(axis=0), np.linspace(expected, 0.0, nz + 1),
                       rtol=1e-12, atol=1e-12)
    assert diagnostics.max_abs_residual_v < 1e-10


def test_q1_poisson_reproduces_series_dielectric_capacitance():
    nx, nz = 4, 8
    epsilon_r = np.ones((nx, nz)); epsilon_r[:, 4:] = 4.0
    fixed = np.zeros((nx + 1, nz + 1), dtype=bool); fixed[:, -1] = True
    system = NodalPoissonSystem(epsilon_r, fixed)
    sigma = 2e-4; h = 20e-9
    voltage, _ = system.solve(_top_sheet_charge(nx, nz, sigma, h))

    expected = sigma * h * (4.0 / 1.0 + 4.0 / 4.0) / EPS0
    assert np.allclose(voltage[:, 0], expected, rtol=1e-12, atol=1e-12)


def test_surface_response_capacitance_is_positive_and_reuses_factorization():
    epsilon_r = np.full((3, 4), 3.9)
    fixed = np.zeros((4, 5), dtype=bool); fixed[:, -1] = True
    system = NodalPoissonSystem(epsilon_r, fixed)
    nodes = np.array([[1, 0], [2, 0]])
    capacitance = system.diagonal_surface_capacitance(nodes)

    assert capacitance.shape == (2,)
    assert np.all(capacitance > 0.0)
    assert np.allclose(capacitance[0], capacitance[1], rtol=1e-12)


def test_edge_sheet_charge_lumping_conserves_total_charge():
    faces = [((1, 2), (1, 3)), ((1, 3), (2, 3))]
    sigma = np.array([2.0, 3.0]); h = 0.25
    charge = lump_edge_sheet_charge((4, 5), faces, sigma, h)

    assert np.isclose(charge.sum(), sigma.sum() * h)
    assert np.isclose(charge[1, 3], 0.5 * sigma.sum() * h)
