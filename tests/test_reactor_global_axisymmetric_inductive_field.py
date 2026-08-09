import numpy as np
import pytest
from scipy.special import j1, jn_zeros

from petch.reactor_global import (
    AxisymmetricFiniteVolumeGrid,
    AxisymmetricInductiveFieldCondition,
    CylindricalReactor,
    DeterministicAxisymmetricInductiveField,
    VACUUM_PERMEABILITY_H_M,
    drude_conductivity_from_two_term_mobilities,
)


def _condition(nr, nz, *, conductivity=18.0 + 0.0j, upper_scale=1.0):
    geometry = CylindricalReactor(radius_m=0.18, length_m=0.10)
    grid = AxisymmetricFiniteVolumeGrid.uniform(
        geometry, radial_cell_count=nr, axial_cell_count=nz)
    alpha = float(jn_zeros(1, 1)[0] / geometry.radius_m)
    upper = upper_scale * j1(alpha * grid.radial_centers_m)
    return AxisymmetricInductiveFieldCondition(
        grid=grid,
        frequency_hz=13.56e6,
        conductivity_S_m=np.full((nr, nz), conductivity, dtype=complex),
        upper_boundary_electric_field_V_m=upper,
        source="manufactured separated cylindrical ICP mode",
        conductivity_source="uniform manufactured plasma conductivity",
        upper_boundary_source="first J1 coil-side boundary mode",
    )


def _exact(condition):
    alpha = float(jn_zeros(1, 1)[0] / condition.grid.geometry.radius_m)
    q = np.sqrt(
        alpha ** 2
        + 1j
        * VACUUM_PERMEABILITY_H_M
        * condition.angular_frequency_rad_s
        * condition.conductivity_S_m[0, 0]
    )
    radial = j1(alpha * condition.grid.radial_centers_m)[:, None]
    axial = np.sinh(q * condition.grid.axial_centers_m[None, :]) / np.sinh(
        q * condition.grid.geometry.length_m)
    return radial * axial


def _weighted_error(solution):
    exact = _exact(solution.condition)
    volume = solution.condition.grid.cell_volume_m3
    return float(np.sqrt(
        np.sum(np.abs(solution.electric_field_V_m - exact) ** 2 * volume)
        / np.sum(np.abs(exact) ** 2 * volume)
    ))


def test_inductive_field_recovers_bessel_skin_mode_and_power_ledger():
    coarse = DeterministicAxisymmetricInductiveField(
        _condition(18, 14)).solve()
    fine = DeterministicAxisymmetricInductiveField(
        _condition(36, 28)).solve()
    assert _weighted_error(fine) < 8.0e-3
    assert _weighted_error(coarse) / _weighted_error(fine) > 3.2
    assert fine.complex_power_ledger_relative_residual < 2.0e-13
    assert fine.linear_system_relative_residual < 2.0e-13
    assert fine.total_absorbed_power_W > 0.0
    volume = fine.condition.grid.cell_volume_m3
    assert (
        np.sum(fine.normalized_absorbed_power_source_moment * volume)
        / fine.condition.grid.geometry.volume_m3
    ) == pytest.approx(1.0, abs=2.0e-13)


def test_inductive_field_exact_parameter_jvp_matches_centered_difference():
    condition = _condition(12, 10, conductivity=14.0 - 3.0j)
    model = DeterministicAxisymmetricInductiveField(condition)
    solution = model.solve()
    radial = condition.grid.radial_centers_m[:, None]
    axial = condition.grid.axial_centers_m[None, :]
    conductivity_tangent = condition.conductivity_S_m * (
        0.08 * radial / condition.grid.geometry.radius_m
        - 0.03j * axial / condition.grid.geometry.length_m
    )
    upper_tangent = (
        0.07 - 0.02j
    ) * condition.upper_boundary_electric_field_V_m
    tangent = model.parameter_jvp(
        solution,
        conductivity_tangent_S_m=conductivity_tangent,
        upper_boundary_electric_field_tangent_V_m=upper_tangent,
    )
    step = 2.0e-5

    def shifted(sign):
        changed = AxisymmetricInductiveFieldCondition(
            grid=condition.grid,
            frequency_hz=condition.frequency_hz,
            conductivity_S_m=(
                condition.conductivity_S_m
                + sign * step * conductivity_tangent),
            upper_boundary_electric_field_V_m=(
                condition.upper_boundary_electric_field_V_m
                + sign * step * upper_tangent),
            source=condition.source,
            conductivity_source=condition.conductivity_source,
            upper_boundary_source=condition.upper_boundary_source,
        )
        return DeterministicAxisymmetricInductiveField(changed).solve()

    plus = shifted(1.0)
    minus = shifted(-1.0)
    finite_field = (
        plus.electric_field_V_m - minus.electric_field_V_m) / (2.0 * step)
    finite_moment = (
        plus.normalized_absorbed_power_source_moment
        - minus.normalized_absorbed_power_source_moment
    ) / (2.0 * step)
    assert tangent.electric_field_tangent_V_m == pytest.approx(
        finite_field, rel=2.0e-8, abs=2.0e-9)
    assert tangent.normalized_source_moment_tangent == pytest.approx(
        finite_moment, rel=3.0e-7, abs=2.0e-8)
    assert tangent.linearized_system_relative_residual < 2.0e-13


def test_drude_conductivity_reproduces_dissipative_mobility():
    electron_density = 2.4e17
    neutral_density = 3.1e20
    flux_reduced = 7.0e24
    dissipative = 5.6e24
    conductivity = drude_conductivity_from_two_term_mobilities(
        electron_density_m3=electron_density,
        neutral_density_m3=neutral_density,
        frequency_hz=13.56e6,
        flux_reduced_mobility_m_inv_V_inv_s_inv=flux_reduced,
        dissipative_reduced_mobility_m_inv_V_inv_s_inv=dissipative,
    )
    expected_real = (
        1.602176634e-19 * electron_density * dissipative / neutral_density)
    assert conductivity.real == pytest.approx(expected_real, rel=2.0e-15)
    assert conductivity.imag < 0.0
    assert abs(conductivity.imag / conductivity.real) == pytest.approx(0.5)
