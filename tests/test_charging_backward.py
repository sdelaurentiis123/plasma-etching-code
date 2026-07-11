import numpy as np
import pytest
from scipy.stats import gamma as gamma_dist, norm, qmc

from petch.charging_backward import (
    _current_balance_diagnostics,
    _interval_current_balance_diagnostics,
    _laplace_residual,
    adaptive_backward_ion_gather,
    backward_electron_gather,
    backward_ion_gather,
)
from petch.charging_general import _trace_general


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


def test_current_balance_diagnostics_match_solver_component_physics():
    Gi = np.array([1.0, 2.0, 2.0, 1e-8])
    Ge = np.array([1.0, 1.0, 3.0, 1e-8])
    comp = np.array([0, 1, 1, 0])

    result = _current_balance_diagnostics(Gi, Ge, comp)

    assert result['active_count'] == 3
    assert result['inactive_count'] == 1
    assert np.isclose(result['log_ratio'][0], 0.0)
    assert np.allclose(result['log_ratio'][1:3], np.log(4.0 / 4.0))
    assert np.isclose(result['pooled'][1]['Gi'], 4.0)
    assert np.isclose(result['pooled'][1]['Ge'], 4.0)
    assert np.isclose(result['max_abs_log_ratio'], 0.0)


def test_current_balance_pools_multiple_faces_of_one_insulator_cell():
    Gi = np.array([1.0, 0.25, 0.5])
    Ge = np.array([0.25, 1.0, 0.5])
    comp = np.zeros(3, dtype=int)
    cells = [(4, 7), (4, 7), (8, 9)]

    result = _current_balance_diagnostics(Gi, Ge, comp, cells)

    assert np.allclose(result['log_ratio'], 0.0)
    assert result['active_count'] == 3
    assert np.isclose(result['max_abs_log_ratio'], 0.0)


def test_interval_current_balance_updates_only_certified_imbalance():
    result = _interval_current_balance_diagnostics(
        Gi=np.array([1.0, 2.0, 0.01]), Ge=np.array([0.9, 1.0, 0.02]),
        Gi_stderr=np.array([0.1, 0.1, 0.01]), Ge_stderr=np.array([0.1, 0.1, 0.01]),
        comp=np.array([0, 1, 0]), cells=[(0, 0), (1, 0), (2, 0)],
        active_flux=0.05, confidence_sigma=2.0)
    assert result['log_ratio'][0] == 0.0
    assert result['log_ratio'][1] > 0.0
    assert not result['active'][2]


def test_laplace_residual_is_zero_for_constant_harmonic_field():
    potential = np.full((12, 9), 7.25)
    gas = np.ones_like(potential, dtype=bool)

    result = _laplace_residual(potential, gas)

    assert result['max_abs'] == 0.0
    assert result['rms'] == 0.0


@pytest.mark.parametrize("floor_potential", [0.0, 9.2])
def test_backward_forward_electron_reciprocity_in_frozen_field_trench(floor_potential):
    """Independent forward launch and backward gather must score the same frozen geometry."""
    nx, nz = 160, 104
    left, right, floor_z = 20, 132, 92
    solid = np.zeros((nx, nz), dtype=bool)
    solid[left, 2:floor_z + 1] = True
    solid[right, 2:floor_z + 1] = True
    solid[left:right + 1, floor_z] = True
    potential_slope = floor_potential / floor_z
    Ex = np.zeros((nx, nz), dtype=float)
    Ez = np.full((nx, nz), -potential_slope, dtype=float)
    potential = np.broadcast_to(potential_slope * np.arange(nz), (nx, nz)).copy()
    floor_cells = [(x, floor_z) for x in range(left + 1, right)]
    floor_normals = [(0.0, -1.0)] * len(floor_cells)

    backward = backward_electron_gather(
        solid, Ex, Ez, potential, floor_cells, floor_normals,
        n_log2=11, n_scramble=3, seed=31,
    ).mean()

    sampler = qmc.Sobol(d=4, scramble=True, seed=47)
    u = sampler.random_base2(16)
    energy = gamma_dist.ppf(u[:, 0], a=2.0, scale=4.0)
    cos_theta = np.sqrt(u[:, 1])
    sin_theta = np.sqrt(1.0 - cos_theta ** 2)
    vx = np.sqrt(energy) * sin_theta * np.cos(2.0 * np.pi * u[:, 2])
    vz = np.sqrt(energy) * cos_theta
    x = left + 1.0 + u[:, 3] * (right - left - 1.0)
    z = np.full_like(x, 0.51)
    hit_x, hit_z, *_ = _trace_general(
        Ex, Ez, solid, x, z, vx, vz, -1.0, nx, nz,
        200 * nz, 0.15, 0.10,
    )
    forward = np.mean((hit_z == floor_z) & (hit_x > left) & (hit_x < right))

    assert np.isclose(backward, forward, rtol=0.04, atol=0.005), (backward, forward)


def test_trace_general_bounds_energy_error_in_uniform_field_at_production_step():
    """Manufactured 1-D orbit: v^2 + qV is constant with V=-z and q=+1."""
    nx, nz = 8, 16
    solid = np.zeros((nx, nz), dtype=bool)
    solid[:, 12:] = True
    Ex = np.zeros((nx, nz), dtype=float)
    Ez = np.ones((nx, nz), dtype=float)
    x = np.array([4.0]); z = np.array([1.0])
    vx = np.array([0.0]); vz = np.array([1.0])

    hit_x, hit_z, impact_energy, *_ = _trace_general(
        Ex, Ez, solid, x, z, vx, vz, 1.0, nx, nz,
        1000, 0.15, 0.10,
    )

    assert hit_x[0] >= 0 and hit_z[0] == 12
    # Cell-crossing impact detection limits this gate; the production step stays below 0.7% error.
    assert np.isclose(impact_energy[0], 12.0, atol=0.08)


@pytest.mark.parametrize("exit_state_weight,exit_energy_mixture", [(False, 0.0), (True, 0.0), (True, 0.2)])
@pytest.mark.parametrize("floor_potential", [0.0, 5.0])
def test_backward_forward_ion_reciprocity_in_frozen_field_trench(
        floor_potential, exit_state_weight, exit_energy_mixture):
    nx, nz = 160, 104
    left, right, floor_z = 20, 132, 92
    solid = np.zeros((nx, nz), dtype=bool)
    solid[left, 2:floor_z + 1] = True
    solid[right, 2:floor_z + 1] = True
    solid[left:right + 1, floor_z] = True
    slope = floor_potential / floor_z
    Ex = np.zeros((nx, nz), dtype=float)
    Ez = np.full((nx, nz), -slope, dtype=float)
    potential = np.broadcast_to(slope * np.arange(nz), (nx, nz)).copy()
    floor_cells = [(x, floor_z) for x in range(left + 1, right)]
    floor_normals = [(0.0, -1.0)] * len(floor_cells)

    backward = backward_ion_gather(
        solid, Ex, Ez, potential, floor_cells, floor_normals,
        n_log2=11, n_scramble=3, seed=53, exit_state_weight=exit_state_weight,
        exit_energy_mixture=exit_energy_mixture,
    ).mean()

    sampler = qmc.Sobol(d=3, scramble=True, seed=59)
    u = sampler.random_base2(16)
    phase = 2.0 * np.pi * u[:, 0]
    sheath_energy = 37.0 + 30.0 * np.sin(phase)
    weights = np.ones_like(sheath_energy)
    vx = np.sqrt(0.25) * norm.ppf(np.clip(u[:, 1], 1e-6, 1.0 - 1e-6))
    vz = np.sqrt(2.0 + sheath_energy)
    x = left + 1.0 + u[:, 2] * (right - left - 1.0)
    z = np.full_like(x, 0.51)
    hit_x, hit_z, *_ = _trace_general(
        Ex, Ez, solid, x, z, vx, vz, 1.0, nx, nz,
        200 * nz, 0.15, 0.10,
    )
    floor_hit = (hit_z == floor_z) & (hit_x > left) & (hit_x < right)
    forward = weights[floor_hit].sum() / weights.sum()

    assert np.isclose(backward, forward, rtol=0.04, atol=0.005), (backward, forward)


def test_ion_exit_state_weight_rejects_benchmark_shaped_phase_law():
    solid, field, potential, cells, normals = _open_flat(0.0)
    with pytest.raises(ValueError, match="uniform RF phase"):
        backward_ion_gather(
            solid, field, field, potential, cells, normals,
            n_log2=5, n_scramble=1, ied_phase_exponent=0.35, exit_state_weight=True,
        )


def test_adaptive_ion_gather_is_geometry_agnostic_on_open_surface():
    solid, field, potential, cells, normals = _open_flat(0.0)
    result = adaptive_backward_ion_gather(
        solid, field, field, potential, cells, normals,
        base_log2=5, max_log2=8, n_replicates=3,
        absolute_tolerance=5e-3, relative_tolerance=0.0,
        ied_phase_exponent=0.0,
    )
    assert result.converged
    assert np.isclose(result.total_mean, 1.0, atol=5e-3)
