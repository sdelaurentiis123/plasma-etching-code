import numpy as np

from petch.charging_backward import (
    backward_electron_gather,
    backward_ion_gather,
)


def _open_flat(V_surface=0.0):
    """Collisionless horizontal wafer with an unobstructed source plane."""
    nx, nz = 24, 16
    solid = np.zeros((nx, nz), dtype=bool)
    solid[:, -1] = True
    field = np.zeros((nx, nz), dtype=float)
    cell = (nx // 2, nz - 1)
    potential = np.zeros((nx, nz), dtype=float)
    potential[cell] = V_surface
    return solid, field, potential, [cell], [(0.0, -1.0)]


def _electron_flux(V_surface):
    solid, field, potential, cells, normals = _open_flat(V_surface)
    return backward_electron_gather(
        solid,
        field,
        field,
        potential,
        cells,
        normals,
        Te=4.0,
        n_log2=11,
        n_scramble=2,
        seed=19,
    )[0]


def test_backward_electron_gather_obeys_langmuir_gate():
    """Gate B: exp(V/Te) retardation below zero and saturation above zero."""
    flux_zero = _electron_flux(0.0)
    flux_retarded = _electron_flux(-4.0)
    flux_positive = _electron_flux(4.0)

    assert np.isclose(flux_zero, 1.0, atol=0.015)
    assert np.isclose(flux_retarded, np.exp(-1.0), atol=0.015)
    assert np.isclose(flux_positive, 1.0, atol=0.015)


def test_backward_ion_gather_is_unity_on_open_uncharged_wafer():
    solid, field, potential, cells, normals = _open_flat(0.0)
    flux = backward_ion_gather(
        solid,
        field,
        field,
        potential,
        cells,
        normals,
        n_log2=11,
        n_scramble=2,
        seed=23,
    )[0]

    assert np.isclose(flux, 1.0, atol=0.01)


def test_backward_gather_is_reproducible_for_fixed_seed():
    first = _electron_flux(-4.0)
    second = _electron_flux(-4.0)

    assert first == second
