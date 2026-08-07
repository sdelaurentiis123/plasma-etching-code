import numpy as np
import pytest

from petch.reactor_global import (
    CylindricalReactor,
    ElectronegativeEdgeFactors,
    ElectropositiveEdgeFactors,
)


def test_cylinder_geometry_matches_analytic_definitions():
    reactor = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    expected_volume = np.pi * 0.1525 ** 2 * 0.075
    expected_area = 2.0 * np.pi * 0.1525 ** 2 + 2.0 * np.pi * 0.1525 * 0.075
    expected_diffusion_length = 1.0 / np.sqrt(
        (np.pi / 0.075) ** 2 + (2.405 / 0.1525) ** 2)
    assert reactor.volume_m3 == expected_volume
    assert reactor.physical_area_m2 == expected_area
    assert reactor.diffusion_length_m == expected_diffusion_length


def test_effective_area_retains_axial_and_radial_edge_factors():
    reactor = CylindricalReactor(radius_m=0.2, length_m=0.1)
    factors = ElectropositiveEdgeFactors(axial=0.25, radial=0.5)
    expected = (
        0.25 * 2.0 * np.pi * 0.2 ** 2
        + 0.5 * 2.0 * np.pi * 0.2 * 0.1
    )
    assert reactor.effective_loss_area_m2(factors) == expected


def test_low_pressure_edge_factors_match_electropositive_equations():
    reactor = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    mean_free_path = 0.03
    factors = reactor.electropositive_edge_factors(
        ion_mean_free_path_m=mean_free_path,
        include_high_pressure_diffusion=False,
    )
    expected_axial = 0.86 / np.sqrt(3.0 + 0.075 / (2.0 * mean_free_path))
    expected_radial = 0.8 / np.sqrt(4.0 + 0.1525 / mean_free_path)
    assert factors.axial == expected_axial
    assert factors.radial == expected_radial


def test_high_pressure_terms_reduce_edge_density_and_require_complete_inputs():
    reactor = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    low = reactor.electropositive_edge_factors(
        ion_mean_free_path_m=0.01,
        include_high_pressure_diffusion=False,
    )
    full = reactor.electropositive_edge_factors(
        ion_mean_free_path_m=0.01,
        bohm_speed_m_s=2500.0,
        ambipolar_diffusion_m2_s=4.0,
    )
    assert full.axial < low.axial
    assert full.radial < low.radial
    with pytest.raises(ValueError, match="require Bohm speed and diffusivity"):
        reactor.electropositive_edge_factors(
            ion_mean_free_path_m=0.01,
            bohm_speed_m_s=2500.0,
        )


def test_zero_electronegativity_recovers_electropositive_edge_factors_exactly():
    reactor = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    electropositive = reactor.electropositive_edge_factors(
        ion_mean_free_path_m=0.01,
        bohm_speed_m_s=2500.0,
        ambipolar_diffusion_m2_s=4.0,
    )
    electronegative = reactor.electronegative_edge_factors(
        electronegativity=0.0,
        electron_to_ion_temperature_ratio=100.0,
        ion_mean_free_path_m=0.01,
        bohm_speed_m_s=2500.0,
        ambipolar_diffusion_m2_s=4.0,
    )
    assert isinstance(electronegative, ElectronegativeEdgeFactors)
    assert electronegative.electronegative_correction == 1.0
    assert electronegative.axial == electropositive.axial
    assert electronegative.radial == electropositive.radial


def test_electronegative_correction_matches_lee_lieberman_equations_13_14():
    reactor = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    alpha = 10.0
    temperature_ratio = 100.0
    baseline = reactor.electropositive_edge_factors(
        ion_mean_free_path_m=0.03,
        include_high_pressure_diffusion=False,
    )
    factors = reactor.electronegative_edge_factors(
        electronegativity=alpha,
        electron_to_ion_temperature_ratio=temperature_ratio,
        ion_mean_free_path_m=0.03,
        include_high_pressure_diffusion=False,
    )
    expected_correction = (
        1.0 + 3.0 * alpha / temperature_ratio
    ) / (1.0 + alpha)
    assert factors.electronegative_correction == expected_correction
    assert factors.axial == expected_correction * baseline.axial
    assert factors.radial == expected_correction * baseline.radial
    assert factors.axial < baseline.axial
    assert factors.radial < baseline.radial
    assert reactor.effective_loss_area_m2(factors) < (
        reactor.effective_loss_area_m2(baseline))


@pytest.mark.parametrize(
    ("electronegativity", "temperature_ratio", "message"),
    [
        (-1.0, 100.0, "electronegativity"),
        (0.0, 0.0, "temperature ratio"),
        (0.0, np.inf, "temperature ratio"),
    ],
)
def test_electronegative_edge_factors_reject_invalid_plasma_ratios(
        electronegativity, temperature_ratio, message):
    reactor = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    with pytest.raises(ValueError, match=message):
        reactor.electronegative_edge_factors(
            electronegativity=electronegativity,
            electron_to_ion_temperature_ratio=temperature_ratio,
            ion_mean_free_path_m=0.03,
            include_high_pressure_diffusion=False,
        )
