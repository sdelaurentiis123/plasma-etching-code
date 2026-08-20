import numpy as np
import pytest

from petch.reactor_global.cylindrical_inventory_lift import (
    CylindricalFiniteVolumeGrid,
    DeterministicCylindricalIndependentInventoryLift,
    normalized_cylindrical_annular_skin_source,
)
from petch.reactor_global.geometry import CylindricalReactor


@pytest.fixture(scope="module")
def grid():
    return CylindricalFiniteVolumeGrid.uniform(
        CylindricalReactor(radius_m=0.16, length_m=0.18),
        radial_cell_count=10, azimuthal_cell_count=12, axial_cell_count=9)


def _lift(grid, source):
    return DeterministicCylindricalIndependentInventoryLift(
        grid=grid, species_names=("F", "ion"),
        diffusion_coefficient_m2_s=(0.4, 5.0),
        wall_velocity_m_s=((120.0, 45.0, 55.0), (1700.0, 500.0, 650.0)),
        source_shape=np.stack((source, source)), source="manufactured source")


def test_axisymmetric_source_recovers_azimuthally_invariant_positive_solution(grid):
    source = normalized_cylindrical_annular_skin_source(
        grid, axial_skin_depth_m=0.045, ring_radius_m=0.10,
        radial_width_m=0.04)
    lift = _lift(grid, source)
    result = lift.solve(np.array([2.0e18, 4.0e16]))

    assert result.maximum_species_ledger_relative_residual < 1.0e-12
    assert result.maximum_inventory_relative_residual < 1.0e-12
    assert result.maximum_linear_system_relative_residual < 1.0e-12
    assert np.min(result.density_m3) >= 0.0
    assert np.ptp(result.lower_endcap_flux_m2_s, axis=2) == pytest.approx(
        0.0, abs=2.0e-12 * np.max(result.lower_endcap_flux_m2_s))


def test_positive_harmonic_source_produces_conserved_nonaxisymmetric_wafer_flux(grid):
    source = normalized_cylindrical_annular_skin_source(
        grid, axial_skin_depth_m=0.045, ring_radius_m=0.12,
        radial_width_m=0.025, cosine_coefficients=(0.3, -0.2),
        sine_coefficients=(-0.1, 0.15))
    lift = _lift(grid, source)
    result = lift.solve(np.array([1.0e18, 1.0e17]))

    assert result.maximum_species_ledger_relative_residual < 1.0e-12
    assert np.any(np.ptp(result.lower_endcap_flux_m2_s, axis=2) > 0.0)
    volume = grid.cell_volume_m3
    recovered = np.sum(result.density_m3 * volume[None], axis=(1, 2, 3)) / (
        grid.geometry.volume_m3)
    assert recovered == pytest.approx([1.0e18, 1.0e17], rel=2.0e-13)


def test_inventory_jvp_is_exact(grid):
    source = normalized_cylindrical_annular_skin_source(
        grid, axial_skin_depth_m=0.05, ring_radius_m=0.09,
        radial_width_m=0.05, cosine_coefficients=(0.1,),
        sine_coefficients=(0.05,))
    lift = _lift(grid, source)
    base = np.array([2.0e18, 3.0e17])
    tangent = np.array([-3.0e16, 8.0e15])
    epsilon = 1.0e-4
    plus = lift.solve(base + epsilon * tangent).density_m3
    minus = lift.solve(base - epsilon * tangent).density_m3
    finite_difference = (plus - minus) / (2.0 * epsilon)

    assert lift.target_inventory_jvp(tangent) == pytest.approx(
        finite_difference, rel=2.0e-10, abs=1.0e3)
    time_tangent = np.stack((tangent, -0.5 * tangent))
    assert lift.density_to_lower_flux_jvp(time_tangent).shape == (
        2, 2, grid.radial_cell_count, grid.azimuthal_cell_count)


def test_new_source_moment_reuses_operator_with_exact_conservation(grid):
    base_source = normalized_cylindrical_annular_skin_source(
        grid, axial_skin_depth_m=0.05, ring_radius_m=0.09,
        radial_width_m=0.05)
    lift = _lift(grid, base_source)
    replay, replay_ledger, replay_linear = (
        lift.source_shape_to_unit_lower_flux(
            np.stack((base_source, base_source))))
    changed = normalized_cylindrical_annular_skin_source(
        grid, axial_skin_depth_m=0.05, ring_radius_m=0.11,
        radial_width_m=0.025, cosine_coefficients=(0.2, -0.1),
        sine_coefficients=(-0.15, 0.05))
    changed_flux, changed_ledger, changed_linear = (
        lift.source_shape_to_unit_lower_flux(np.stack((changed, changed))))

    assert replay == pytest.approx(
        lift._unit_lower_flux_per_density_m_s, rel=2.0e-14, abs=1.0e-16)
    assert max(replay_ledger, changed_ledger) < 1.0e-12
    assert max(replay_linear, changed_linear) < 1.0e-12
    assert np.all(changed_flux >= 0.0)
    assert np.any(np.ptp(changed_flux, axis=2) > 0.0)
