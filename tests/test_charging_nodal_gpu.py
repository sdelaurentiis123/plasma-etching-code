import numpy as np
import pytest

from petch.charging_nodal import trace_nodal_cpu

wp = pytest.importorskip("warp")
from petch.charging_nodal_gpu import trace_nodal_warp


def _compare(cpu, warp, atol=2e-12):
    for index in (0, 1, 4, 7, 8):
        assert np.array_equal(cpu[index], warp[index]), index
    for index in (2, 3, 5, 6, 9, 10):
        assert np.allclose(cpu[index], warp[index], rtol=2e-12, atol=atol), index


def test_warp_nodal_matches_exact_face_hits_exits_and_reflections():
    nx, nz = 9, 8
    solid = np.zeros((nx, nz), dtype=bool)
    solid[6, :] = True
    solid[2:6, 6:] = True
    potential = np.zeros((nx + 1, nz + 1))
    x = np.array([2.25, 3.95, 1.2, 0.05])
    z = np.array([3.4, 5.8, 0.2, 2.0])
    vx = np.array([30.0, 1.0, 0.1, -2.0])
    vz = np.array([0.0, 3.0, -2.0, -0.2])
    arguments = (potential, solid, x, z, vx, vz, 1.0, nx, nz, 500, 0.4, 0.1)

    cpu = trace_nodal_cpu(*arguments)
    warp = trace_nodal_warp(*arguments, device="cpu")

    _compare(cpu, warp)


@pytest.mark.parametrize("charge", [-1.0, 1.0])
def test_warp_nodal_matches_nonuniform_q1_orbits(charge):
    nx, nz = 18, 14
    solid = np.zeros((nx, nz), dtype=bool)
    solid[:4, 1:] = True; solid[14:, 1:] = True; solid[4:14, 11:] = True
    ii, jj = np.meshgrid(np.arange(nx + 1), np.arange(nz + 1), indexing="ij")
    potential = 0.025 * jj + 0.0012 * (ii - nx / 2.0) * jj
    rng = np.random.default_rng(314159)
    count = 48
    x = rng.uniform(4.1, 13.9, count)
    z = rng.uniform(0.05, 9.8, count)
    vx = rng.normal(0.0, 1.2, count)
    vz = rng.uniform(0.4, 4.0, count)
    arguments = (
        potential, solid, x, z, vx, vz, charge, nx, nz, 1000, 0.15, 0.1, 0.02)

    cpu = trace_nodal_cpu(*arguments)
    warp = trace_nodal_warp(*arguments, device="cpu")

    _compare(cpu, warp, atol=2e-10)
