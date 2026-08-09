import numpy as np
import pytest

from petch.reactor_global import (
    AxisymmetricDriftDiffusionCondition,
    AxisymmetricFiniteVolumeGrid,
    AxisymmetricReactionDiffusionCondition,
    CylindricalReactor,
    DeterministicAxisymmetricDriftDiffusion,
    DeterministicAxisymmetricReactionDiffusion,
    DeterministicQuasineutralInventoryLift,
    normalized_exponential_skin_source,
)


def _grid(nr=10, nz=8):
    return AxisymmetricFiniteVolumeGrid.uniform(
        CylindricalReactor(radius_m=0.15, length_m=0.08),
        radial_cell_count=nr,
        axial_cell_count=nz,
    )


def test_zero_potential_reduces_exactly_to_reaction_diffusion():
    grid = _grid()
    source = np.full((1, 10, 8), 3.0e18)
    wall = np.array([[200.0, 300.0, 40.0]])
    diffusion = np.array([0.24])
    reaction = np.array([[4.0]])
    drift_condition = AxisymmetricDriftDiffusionCondition(
        grid=grid,
        species_names=("ion",),
        charge_number=np.array([1.0]),
        mobility_m2_V_s=np.array([1.2]),
        temperature_eV=np.array([0.2]),
        electrostatic_potential_V=np.zeros((10, 8)),
        volume_reaction_matrix_s_inv=reaction,
        source_rate_m3_s=source,
        wall_velocity_m_s=wall,
        source="zero-potential equivalence test",
    )
    reaction_condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=("ion",),
        diffusion_coefficient_m2_s=diffusion,
        volume_reaction_matrix_s_inv=reaction,
        source_rate_m3_s=source,
        wall_velocity_m_s=wall,
        source="zero-potential equivalence test",
    )
    drift = DeterministicAxisymmetricDriftDiffusion(drift_condition).solve()
    scalar = DeterministicAxisymmetricReactionDiffusion(
        reaction_condition).solve()
    assert drift.density_m3 == pytest.approx(
        scalar.density_m3, rel=2.0e-13)
    assert drift.lower_endcap_flux_m2_s == pytest.approx(
        scalar.lower_endcap_flux_m2_s, rel=2.0e-13)
    assert drift.maximum_species_ledger_relative_residual < 1.0e-12


def test_scharfetter_gummel_preserves_discrete_boltzmann_equilibrium():
    grid = _grid(12, 10)
    radial = grid.radial_centers_m[:, None]
    axial = grid.axial_centers_m[None, :]
    potential = 0.7 * radial / grid.geometry.radius_m + 0.3 * axial / grid.geometry.length_m
    temperature = 0.25
    equilibrium = 2.0e16 * np.exp(-potential / temperature)
    condition = AxisymmetricDriftDiffusionCondition(
        grid=grid,
        species_names=("positive",),
        charge_number=np.array([1.0]),
        mobility_m2_V_s=np.array([0.8]),
        temperature_eV=np.array([temperature]),
        electrostatic_potential_V=potential,
        volume_reaction_matrix_s_inv=np.array([[3.0]]),
        source_rate_m3_s=(3.0 * equilibrium)[None, :, :],
        wall_velocity_m_s=np.zeros((1, 3)),
        source="discrete Boltzmann equilibrium",
    )
    solution = DeterministicAxisymmetricDriftDiffusion(condition).solve()
    assert solution.density_m3[0] == pytest.approx(
        equilibrium, rel=2.0e-12)
    assert solution.integrated_wall_loss_rate_s[0] == pytest.approx(0.0)


def test_quasineutral_lift_conserves_inventory_and_charge_everywhere():
    grid = _grid(14, 12)
    source = normalized_exponential_skin_source(
        grid, axial_skin_depth_m=0.025, radial_scale_m=0.12)
    lift = DeterministicQuasineutralInventoryLift(
        grid=grid,
        species_names=("Cl+", "Cl2+", "Cl-"),
        charge_number=np.array([1.0, 1.0, -1.0]),
        mobility_m2_V_s=np.array([18.0, 14.0, 20.0]),
        ion_temperature_eV=np.array([0.15, 0.15, 0.15]),
        electron_temperature_eV=3.0,
        wall_velocity_m_s=np.array([
            [2500.0, 2500.0, 2500.0],
            [1800.0, 1800.0, 1800.0],
            [np.inf, np.inf, np.inf],
        ]),
        source_shape=np.stack((source, source, source)),
        source="manufactured electronegative chlorine lift",
        positive_negative_recombination_m3_s=5.0e-14,
    )
    target = np.array([3.0e16, 8.0e16, 7.0e16])
    result = lift.solve(
        target,
        relative_tolerance=5.0e-8,
        maximum_iterations=800,
        relaxation_fraction=0.05,
    )
    assert result.recovered_volume_average_density_m3 == pytest.approx(
        target, rel=2.0e-12)
    assert np.min(result.electron_density_m3) > 0.0
    assert result.solution.maximum_species_ledger_relative_residual < 1.0e-11
    assert result.maximum_inventory_relative_residual < 1.0e-12
    assert result.nonlinear_solver_evaluations < 800
    assert result.supports_implicit_differentiation is True


def test_quasineutral_lift_implicit_inventory_jvp_matches_centered_difference():
    grid = _grid(6, 5)
    source = normalized_exponential_skin_source(
        grid, axial_skin_depth_m=0.030, radial_scale_m=0.13)
    lift = DeterministicQuasineutralInventoryLift(
        grid=grid,
        species_names=("Cl+", "Cl2+", "Cl-"),
        charge_number=np.array([1.0, 1.0, -1.0]),
        mobility_m2_V_s=np.array([18.0, 14.0, 20.0]),
        ion_temperature_eV=np.array([0.15, 0.15, 0.15]),
        electron_temperature_eV=3.0,
        wall_velocity_m_s=np.array([
            [2500.0, 2500.0, 2500.0],
            [1800.0, 1800.0, 1800.0],
            [np.inf, np.inf, np.inf],
        ]),
        source_shape=np.stack((source, source, source)),
        source="manufactured implicit-JVP chlorine lift",
        positive_negative_recombination_m3_s=5.0e-14,
    )
    target = np.array([3.0e16, 8.0e16, 7.0e16])
    tangent = target * np.array([0.11, -0.07, 0.05])
    result = lift.solve(
        target, relative_tolerance=2.0e-9, maximum_iterations=800)
    derivative = lift.implicit_target_inventory_jvp(
        result, tangent, relative_tolerance=2.0e-9)
    step = 2.0e-4
    plus = lift.solve(
        target + step * tangent,
        relative_tolerance=2.0e-9,
        maximum_iterations=800,
        initial_electrostatic_potential_V=result.electrostatic_potential_V,
    )
    minus = lift.solve(
        target - step * tangent,
        relative_tolerance=2.0e-9,
        maximum_iterations=800,
        initial_electrostatic_potential_V=result.electrostatic_potential_V,
    )
    finite_density = (
        plus.solution.density_m3 - minus.solution.density_m3
    ) / (2.0 * step)
    finite_potential = (
        plus.electrostatic_potential_V - minus.electrostatic_potential_V
    ) / (2.0 * step)
    assert derivative.density_tangent_m3 == pytest.approx(
        finite_density, rel=3.0e-5, abs=2.0e8)
    assert derivative.electrostatic_potential_tangent_V == pytest.approx(
        finite_potential, rel=2.0e-5, abs=2.0e-7)
    assert derivative.maximum_linearized_fixed_point_relative_residual < 1.0e-8
    assert derivative.maximum_inventory_tangent_relative_residual < 1.0e-12


def test_drift_diffusion_condition_rejects_neutral_charge_and_bad_reaction():
    grid = _grid(4, 4)
    base = dict(
        grid=grid,
        species_names=("ion",),
        mobility_m2_V_s=np.ones(1),
        temperature_eV=np.ones(1),
        electrostatic_potential_V=np.zeros((4, 4)),
        volume_reaction_matrix_s_inv=np.ones((1, 1)),
        source_rate_m3_s=np.ones((1, 4, 4)),
        wall_velocity_m_s=np.ones((1, 3)),
        source="invalid charge gate",
    )
    with pytest.raises(ValueError):
        AxisymmetricDriftDiffusionCondition(
            charge_number=np.array([0.0]), **base)
