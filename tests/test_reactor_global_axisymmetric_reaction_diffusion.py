import numpy as np
import pytest
from scipy.special import j0, jn_zeros

from petch.reactor_global import (
    AxisymmetricFiniteVolumeGrid,
    AxisymmetricReactionDiffusionCondition,
    CylindricalReactor,
    DeterministicAxisymmetricInventoryLift,
    DeterministicAxisymmetricReactionDiffusion,
    normalized_exponential_skin_source,
    normalized_annular_skin_source,
)


def _uniform_grid(n=12):
    return AxisymmetricFiniteVolumeGrid.uniform(
        CylindricalReactor(radius_m=0.2, length_m=0.1),
        radial_cell_count=n,
        axial_cell_count=n,
    )


def test_reflecting_uniform_volume_loss_recovers_exact_well_mixed_limit():
    grid = _uniform_grid()
    condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=("Cl",),
        diffusion_coefficient_m2_s=np.array([0.12]),
        volume_reaction_matrix_s_inv=np.array([[5.0]]),
        source_rate_m3_s=np.full((1, 12, 12), 10.0),
        wall_velocity_m_s=np.zeros((1, 3)),
        source="manufactured reflecting well-mixed limit",
    )
    solution = DeterministicAxisymmetricReactionDiffusion(condition).solve()
    assert solution.density_m3 == pytest.approx(2.0, rel=2.0e-13)
    assert solution.integrated_wall_loss_rate_s[0] == pytest.approx(0.0)
    assert solution.integrated_source_rate_s[0] == pytest.approx(
        solution.integrated_volume_reaction_rate_s[0], rel=2.0e-13)
    assert solution.maximum_species_ledger_relative_residual < 1.0e-12


def _bessel_sine_error(cell_count):
    geometry = CylindricalReactor(radius_m=0.2, length_m=0.1)
    grid = AxisymmetricFiniteVolumeGrid.uniform(
        geometry,
        radial_cell_count=cell_count,
        axial_cell_count=cell_count,
    )
    diffusion = 0.08
    loss = 3.0
    alpha = float(jn_zeros(0, 1)[0])
    radial = grid.radial_centers_m[:, None]
    axial = grid.axial_centers_m[None, :]
    exact = (
        j0(alpha * radial / geometry.radius_m)
        * np.sin(np.pi * axial / geometry.length_m)
    )
    eigenvalue = (
        (alpha / geometry.radius_m) ** 2
        + (np.pi / geometry.length_m) ** 2)
    source = (diffusion * eigenvalue + loss) * exact
    condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=("ion",),
        diffusion_coefficient_m2_s=np.array([diffusion]),
        volume_reaction_matrix_s_inv=np.array([[loss]]),
        source_rate_m3_s=source[None, :, :],
        wall_velocity_m_s=np.full((1, 3), np.inf),
        source="Bessel-sine continuum manufactured solution",
    )
    numerical = DeterministicAxisymmetricReactionDiffusion(
        condition).solve().density_m3[0]
    weight = grid.cell_volume_m3
    return float(np.sqrt(
        np.sum(weight * (numerical - exact) ** 2)
        / np.sum(weight * exact ** 2)))


def test_absorbing_bessel_sine_mode_converges_second_order():
    coarse = _bessel_sine_error(12)
    fine = _bessel_sine_error(24)
    assert fine < 0.32 * coarse
    assert fine < 3.0e-3


def test_coupled_conversion_matrix_conserves_and_remains_positive():
    grid = _uniform_grid(10)
    source = np.zeros((2, 10, 10))
    source[0] = 12.0
    condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=("A", "B"),
        diffusion_coefficient_m2_s=np.array([0.1, 0.05]),
        volume_reaction_matrix_s_inv=np.array([
            [4.0, 0.0],
            [-4.0, 3.0],
        ]),
        source_rate_m3_s=source,
        wall_velocity_m_s=np.zeros((2, 3)),
        source="manufactured conservative A to B conversion",
    )
    solution = DeterministicAxisymmetricReactionDiffusion(condition).solve()
    assert solution.density_m3[0] == pytest.approx(3.0, rel=3.0e-13)
    assert solution.density_m3[1] == pytest.approx(4.0, rel=3.0e-13)
    assert solution.maximum_species_ledger_relative_residual < 1.0e-12


def test_wall_flux_ledger_partial_wafer_and_exact_source_jvp():
    grid = _uniform_grid(9)
    shape = normalized_exponential_skin_source(
        grid,
        axial_skin_depth_m=0.025,
        radial_scale_m=0.15,
    )
    assert np.sum(shape * grid.cell_volume_m3) == pytest.approx(
        grid.geometry.volume_m3)
    condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=("Cl+",),
        diffusion_coefficient_m2_s=np.array([0.2]),
        volume_reaction_matrix_s_inv=np.array([[0.0]]),
        source_rate_m3_s=(2.0e20 * shape)[None, :, :],
        wall_velocity_m_s=np.array([[600.0, 300.0, 50.0]]),
        source="ICP-like ion source and asymmetric wall collection",
    )
    model = DeterministicAxisymmetricReactionDiffusion(condition)
    solution = model.solve()
    wafer_flux = solution.lower_endcap_area_average_flux_m2_s(
        "Cl+", wafer_radius_m=0.1)
    full_flux = solution.lower_endcap_area_average_flux_m2_s(
        "Cl+", wafer_radius_m=0.2)
    assert wafer_flux > 0.0
    assert full_flux > 0.0
    assert solution.integrated_source_rate_s[0] == pytest.approx(
        solution.integrated_wall_loss_rate_s[0], rel=2.0e-12)

    tangent = (0.3e20 * shape)[None, :, :]
    jvp = model.source_jvp(tangent)
    epsilon = 1.0e-5
    perturbed_condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=condition.species_names,
        diffusion_coefficient_m2_s=condition.diffusion_coefficient_m2_s,
        volume_reaction_matrix_s_inv=condition.volume_reaction_matrix_s_inv,
        source_rate_m3_s=condition.source_rate_m3_s + epsilon * tangent,
        wall_velocity_m_s=condition.wall_velocity_m_s,
        source="source-JVP finite perturbation check",
    )
    perturbed = DeterministicAxisymmetricReactionDiffusion(
        perturbed_condition).solve()
    finite_difference = (
        perturbed.density_m3 - solution.density_m3) / epsilon
    assert jvp == pytest.approx(finite_difference, rel=2.0e-9, abs=1.0e-3)


def test_non_m_matrix_and_unanchored_system_fail_closed():
    grid = _uniform_grid(4)
    with pytest.raises(ValueError):
        AxisymmetricReactionDiffusionCondition(
            grid=grid,
            species_names=("A", "B"),
            diffusion_coefficient_m2_s=np.ones(2),
            volume_reaction_matrix_s_inv=np.array([[1.0, 0.2], [0.0, 1.0]]),
            source_rate_m3_s=np.ones((2, 4, 4)),
            wall_velocity_m_s=np.ones((2, 3)),
            source="invalid positive offdiagonal",
        )
    condition = AxisymmetricReactionDiffusionCondition(
        grid=grid,
        species_names=("A",),
        diffusion_coefficient_m2_s=np.ones(1),
        volume_reaction_matrix_s_inv=np.zeros((1, 1)),
        source_rate_m3_s=np.ones((1, 4, 4)),
        wall_velocity_m_s=np.zeros((1, 3)),
        source="unanchored reflecting diffusion",
    )
    with pytest.raises(ValueError, match="singular"):
        DeterministicAxisymmetricReactionDiffusion(condition)


def test_inventory_lift_recovers_global_state_and_local_wafer_partition():
    grid = _uniform_grid(11)
    top_center = normalized_exponential_skin_source(
        grid,
        axial_skin_depth_m=0.018,
        radial_scale_m=0.09,
    )
    uniform = np.ones_like(top_center)
    lift = DeterministicAxisymmetricInventoryLift(
        grid=grid,
        species_names=("Cl+", "Cl2+"),
        diffusion_coefficient_m2_s=np.array([0.18, 0.11]),
        volume_reaction_matrix_s_inv=np.diag([2.0, 1.0]),
        wall_velocity_m_s=np.array([
            [500.0, 300.0, 60.0],
            [350.0, 250.0, 40.0],
        ]),
        source_shape=np.stack((top_center, uniform)),
        source="manufactured global-inventory spatial lift",
    )
    target = np.array([2.5e16, 1.2e16])
    result = lift.solve(target)
    assert result.recovered_volume_average_density_m3 == pytest.approx(
        target, rel=3.0e-13)
    assert result.maximum_inventory_relative_residual < 1.0e-12
    assert np.all(result.inferred_source_amplitude_m3_s > 0.0)
    center_wafer_flux = (
        result.solution.lower_endcap_area_average_flux_m2_s(
            "Cl+", wafer_radius_m=0.1))
    full_endcap_flux = (
        result.solution.lower_endcap_area_average_flux_m2_s(
            "Cl+", wafer_radius_m=0.2))
    assert center_wafer_flux > full_endcap_flux
    assert result.supports_reactor_state_prediction is False

    tangent = np.array([0.3e16, -0.1e16])
    jvp = lift.target_inventory_jvp(tangent)
    epsilon = 1.0e-5
    perturbed = lift.solve(target + epsilon * tangent)
    finite_difference = (
        perturbed.solution.density_m3 - result.solution.density_m3
    ) / epsilon
    assert jvp == pytest.approx(finite_difference, rel=2.0e-9, abs=5.0e5)


def test_inventory_lift_rejects_incompatible_conversion_inventory():
    grid = _uniform_grid(6)
    lift = DeterministicAxisymmetricInventoryLift(
        grid=grid,
        species_names=("A", "B"),
        diffusion_coefficient_m2_s=np.array([0.1, 0.1]),
        volume_reaction_matrix_s_inv=np.array([
            [5.0, 0.0],
            [-5.0, 1.0],
        ]),
        wall_velocity_m_s=np.full((2, 3), 0.1),
        source_shape=np.ones((2, 6, 6)),
        source="manufactured incompatible A-to-B inventory",
    )
    with pytest.raises(ValueError, match="negative source"):
        lift.solve(np.array([2.0e16, 1.0e10]))
