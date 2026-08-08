import numpy as np
import pytest

from petch.reactor_global import (
    ElectronEnergyDistribution,
    ElectronEnergyGrid,
    IsotropicElectronElectronCoulombKernel,
    ScharfetterGummelEnergyOperator,
)


def test_coulomb_kernel_integrals_are_exact_for_constant_eepf():
    maximum = 20.0
    grid = ElectronEnergyGrid.linear(maximum, 80)
    distribution = ElectronEnergyDistribution.from_unnormalized(
        grid, np.ones(grid.cell_count))
    state = IsotropicElectronElectronCoulombKernel(grid).evaluate(
        distribution,
        electron_to_neutral_density_ratio=7.0e-4,
        gas_number_density_m3=6.0e19,
    )
    constant = 3.0 / (2.0 * maximum ** 1.5)
    faces = grid.boundaries[1:-1]
    np.testing.assert_allclose(
        state.lower_population_integral,
        constant * (2.0 / 3.0) * faces ** 1.5,
        rtol=2.0e-14,
        atol=2.0e-15,
    )
    np.testing.assert_allclose(
        state.lower_energy_integral_eV,
        constant * (2.0 / 5.0) * faces ** 2.5,
        rtol=2.0e-14,
        atol=2.0e-15,
    )
    np.testing.assert_allclose(
        state.upper_eepf_integral_eV_minus_1_over_2,
        constant * (maximum - faces),
        rtol=2.0e-14,
        atol=2.0e-15,
    )
    assert state.kinetic_temperature_eV == pytest.approx(
        2.0 * distribution.mean_energy_eV / 3.0)
    assert state.coulomb_logarithm > 1.0
    assert np.all(state.drift_eV_m3_s < 0.0)
    assert np.all(state.diffusion_eV2_m3_s > 0.0)
    assert not state.supports_reactor_state_prediction
    assert not state.supports_feature_depth


def test_coulomb_maxwellian_drift_diffusion_identity_converges():
    temperature = 3.0
    errors = []
    for cells in (400, 1600):
        grid = ElectronEnergyGrid.linear(40.0, cells)
        distribution = ElectronEnergyDistribution.maxwellian(
            grid, temperature)
        state = IsotropicElectronElectronCoulombKernel(grid).evaluate(
            distribution,
            electron_to_neutral_density_ratio=7.0e-4,
            gas_number_density_m3=6.0e19,
        )
        drift_to_diffusion = (
            state.drift_eV_m3_s / state.diffusion_eV2_m3_s)
        errors.append(float(np.max(np.abs(
            temperature * drift_to_diffusion + 1.0))))
    # Equation 38 has exactly zero flux for any Maxwellian. The remaining
    # error is the declared piecewise-constant EEPF representation and must
    # fall under refinement at close to first order.
    assert errors[1] < 1.0e-3
    assert errors[1] < 0.27 * errors[0]


def test_coulomb_operator_conserves_particles_and_energy_under_refinement():
    relative_energy_defects = []
    for cells in (200, 800):
        grid = ElectronEnergyGrid.linear(200.0, cells)
        cold = ElectronEnergyDistribution.maxwellian(
            grid, 1.5).eepf_eV_minus_3_over_2
        hot = ElectronEnergyDistribution.maxwellian(
            grid, 6.0).eepf_eV_minus_3_over_2
        distribution = ElectronEnergyDistribution(
            grid, 0.4 * cold + 0.6 * hot)
        state = IsotropicElectronElectronCoulombKernel(grid).evaluate(
            distribution,
            electron_to_neutral_density_ratio=7.0e-4,
            gas_number_density_m3=6.0e19,
        )
        operator = ScharfetterGummelEnergyOperator(
            grid, state.drift_eV_m3_s, state.diffusion_eV2_m3_s)
        source = operator._matrix(0) @ (
            distribution.eepf_eV_minus_3_over_2)
        assert np.sum(source) == pytest.approx(0.0, abs=2.0e-28)
        signed_energy = float(np.dot(grid.cell_centers_eV, source))
        absolute_energy = float(np.sum(np.abs(
            grid.cell_centers_eV * source)))
        relative_energy_defects.append(abs(signed_energy) / absolute_energy)

    # The continuous Landau term conserves electron energy. The finite-volume
    # defect from representing the nonlinear EEPF as cell constants must
    # collapse under refinement rather than acting as a hidden power source.
    assert relative_energy_defects[1] < 3.0e-3
    assert relative_energy_defects[1] < 0.12 * relative_energy_defects[0]


def test_coulomb_jvp_vjp_are_exact_batch_adjoint_pair():
    rng = np.random.default_rng(9417)
    grid = ElectronEnergyGrid.linear(30.0, 120)
    distribution = ElectronEnergyDistribution.maxwellian(
        grid, np.array([2.2, 3.4]))
    tangent = rng.normal(size=distribution.eepf_eV_minus_3_over_2.shape)
    # The derivative is defined on raw EEPF values; use a normalization-
    # preserving direction so it is also tangent to the public distribution.
    tangent -= (
        np.sum(tangent * grid.normalization_weights, axis=-1)
        / np.sum(grid.normalization_weights ** 2)
    )[..., np.newaxis] * grid.normalization_weights
    ratio = np.array([3.0e-4, 8.0e-4])
    ratio_tangent = np.array([2.0e-5, -3.0e-5])
    gas_density = np.array([4.0e19, 7.0e19])
    gas_tangent = np.array([2.0e18, -1.0e18])
    kernel = IsotropicElectronElectronCoulombKernel(grid)
    drift_tangent, diffusion_tangent = kernel.jvp(
        distribution,
        tangent,
        electron_to_neutral_density_ratio=ratio,
        electron_to_neutral_density_ratio_tangent=ratio_tangent,
        gas_number_density_m3=gas_density,
        gas_number_density_m3_tangent=gas_tangent,
    )
    drift_bar = rng.normal(size=drift_tangent.shape)
    diffusion_bar = rng.normal(size=diffusion_tangent.shape)
    eepf_bar, ratio_bar, gas_bar = kernel.vjp(
        distribution,
        drift_bar,
        diffusion_bar,
        electron_to_neutral_density_ratio=ratio,
        gas_number_density_m3=gas_density,
    )
    forward_dot = float(np.sum(drift_tangent * drift_bar) + np.sum(
        diffusion_tangent * diffusion_bar))
    reverse_dot = float(
        np.sum(tangent * eepf_bar)
        + np.sum(ratio_tangent * ratio_bar)
        + np.sum(gas_tangent * gas_bar)
    )
    assert forward_dot == pytest.approx(reverse_dot, rel=2.0e-13, abs=1e-27)


def test_coulomb_kernel_fails_closed_outside_classical_log_domain():
    grid = ElectronEnergyGrid.linear(1.0e-5, 20)
    distribution = ElectronEnergyDistribution.maxwellian(grid, 1.0e-6)
    with pytest.raises(ValueError, match="outside its domain"):
        IsotropicElectronElectronCoulombKernel(grid).evaluate(
            distribution,
            electron_to_neutral_density_ratio=1.0,
            gas_number_density_m3=1.0e30,
        )
