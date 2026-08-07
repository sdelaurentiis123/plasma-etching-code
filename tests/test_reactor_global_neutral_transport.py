import numpy as np
import pytest
from scipy.special import jn_zeros

from petch.reactor_global import CylindricalReactor
from petch.reactor_global.neutral_transport import (
    ReducedNeutralDiffusivity,
    solve_cylindrical_neutral_wall_loss,
)


def test_exact_cylindrical_robin_mode_satisfies_both_wall_boundaries():
    geometry = CylindricalReactor(radius_m=0.14, length_m=0.10)
    state = solve_cylindrical_neutral_wall_loss(
        geometry=geometry,
        diffusivity_m2_s=0.8,
        mean_thermal_speed_m_s=420.0,
        wall_reaction_probability=0.02,
    )
    assert state.numerical_closure_passes
    assert state.maximum_eigen_residual < 1.0e-12
    assert 0.0 < state.radial_dimensionless_root < jn_zeros(0, 1)[0]
    assert 0.0 < state.axial_dimensionless_root < 0.5 * np.pi
    assert (
        0.0
        < state.exact_loss_frequency_s_inv
        < state.surface_limit_frequency_s_inv
    )
    assert (
        state.exact_loss_frequency_s_inv
        < state.absorbing_wall_frequency_s_inv
    )


def test_exact_mode_reaches_absorbing_and_surface_reaction_limits():
    geometry = CylindricalReactor(radius_m=0.14, length_m=0.10)
    absorbing = solve_cylindrical_neutral_wall_loss(
        geometry=geometry,
        diffusivity_m2_s=1.0e-10,
        mean_thermal_speed_m_s=1000.0,
        wall_reaction_probability=1.0,
    )
    assert absorbing.exact_loss_frequency_s_inv == pytest.approx(
        absorbing.absorbing_wall_frequency_s_inv,
        rel=2.0e-7,
    )

    reaction_limited = solve_cylindrical_neutral_wall_loss(
        geometry=geometry,
        diffusivity_m2_s=1.0e5,
        mean_thermal_speed_m_s=400.0,
        wall_reaction_probability=0.01,
    )
    assert reaction_limited.exact_loss_frequency_s_inv == pytest.approx(
        reaction_limited.surface_limit_frequency_s_inv,
        rel=2.0e-6,
    )


def test_nonreactive_wall_has_zero_loss_and_infinite_residence_time():
    state = solve_cylindrical_neutral_wall_loss(
        geometry=CylindricalReactor(radius_m=0.14, length_m=0.10),
        diffusivity_m2_s=1.0,
        mean_thermal_speed_m_s=400.0,
        wall_reaction_probability=0.0,
    )
    assert state.exact_loss_frequency_s_inv == 0.0
    assert state.chantry_loss_frequency_s_inv == 0.0
    assert np.isinf(state.extrapolation_length_m)
    assert np.isinf(state.residence_time_s)
    assert state.numerical_closure_passes


@pytest.mark.parametrize("probability", (0.001, 0.01, 0.1, 1.0))
@pytest.mark.parametrize("diffusivity", (0.01, 0.1, 1.0, 10.0))
def test_chantry_resistance_sum_stays_within_its_published_cylinder_error(
    probability,
    diffusivity,
):
    state = solve_cylindrical_neutral_wall_loss(
        geometry=CylindricalReactor(radius_m=0.14, length_m=0.10),
        diffusivity_m2_s=diffusivity,
        mean_thermal_speed_m_s=420.0,
        wall_reaction_probability=probability,
    )
    relative_error = abs(
        state.chantry_loss_frequency_s_inv
        / state.exact_loss_frequency_s_inv
        - 1.0
    )
    assert relative_error <= 0.11


def test_reduced_diffusivity_refuses_temperature_extrapolation():
    model = ReducedNeutralDiffusivity(
        reduced_diffusivity_m_inv_s=6.21e20,
        reference_temperature_K=500.0,
        valid_temperature_K=(500.0, 500.0),
        source="primary published model",
        evidence_kind="published_model",
        relative_uncertainty=None,
        provenance={"fit_target": None},
    )
    state = model.evaluate(
        total_neutral_density_m3=1.0e21,
        gas_temperature_K=500.0,
    )
    assert state.diffusivity_m2_s == pytest.approx(0.621)
    assert state.provenance["fit_target"] is None
    assert not state.supports_prediction
    with pytest.raises(ValueError, match="temperature"):
        model.evaluate(
            total_neutral_density_m3=1.0e21,
            gas_temperature_K=300.0,
        )
