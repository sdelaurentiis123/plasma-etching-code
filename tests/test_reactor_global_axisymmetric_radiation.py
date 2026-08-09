import numpy as np
import pytest
from scipy.constants import atomic_mass

from petch.reactor_global import (
    AxisymmetricRadiationMomentField,
    CylindricalReactor,
    ResonanceLineData,
    certify_deterministic_axisymmetric_resonance_escape,
    deterministic_axisymmetric_resonance_escape,
    deterministic_cylinder_resonance_escape,
)


LINE = ResonanceLineData(
    wavelength_nm=106.665980,
    transition_probability_s_inv=1.32e8,
    lower_statistical_weight=1.0,
    upper_statistical_weight=3.0,
    absorber_mass_kg=39.948 * atomic_mass,
    source="nist-argon-persistent-lines",
)


def _field(temperature, absorber, emitter):
    return AxisymmetricRadiationMomentField(
        radial_nodes_m=np.array([0.0, 0.05625, 0.1125]),
        axial_nodes_m=np.array([0.0, 0.06, 0.12]),
        gas_temperature_K=np.asarray(temperature, dtype=float),
        absorber_density_m3=np.asarray(absorber, dtype=float),
        emitter_density_m3=np.asarray(emitter, dtype=float),
        source="synthetic deterministic moment test",
    )


def test_axisymmetric_uniform_field_matches_homogeneous_cylinder_limit():
    shape = (3, 3)
    field = _field(
        np.full(shape, 500.0),
        np.full(shape, 2.0e20),
        np.ones(shape),
    )
    moment = deterministic_axisymmetric_resonance_escape(
        field, LINE, source_quadrature_order=6,
        direction_quadrature_order=6, path_quadrature_order=12,
        frequency_quadrature_order=24)
    homogeneous = deterministic_cylinder_resonance_escape(
        CylindricalReactor(radius_m=0.1125, length_m=0.12),
        LINE, absorber_density_m3=2.0e20, gas_temperature_K=500.0,
        geometry_quadrature_order=8, frequency_quadrature_order=24,
        maximum_doppler_standard_deviations=1024.0)
    assert moment.trapping_factor == pytest.approx(
        homogeneous.trapping_factor, rel=2.0e-14)


def test_cold_dense_shell_and_line_specific_emitter_change_escape():
    uniform = _field(
        np.full((3, 3), 700.0),
        np.full((3, 3), 1.5e20),
        np.ones((3, 3)),
    )
    cold_shell_temperature = np.array([
        [700.0, 700.0, 700.0],
        [600.0, 650.0, 600.0],
        [350.0, 350.0, 350.0],
    ])
    cold_shell_density = 1.5e20 * 700.0 / cold_shell_temperature
    core_emitter = np.array([
        [1.0, 3.0, 1.0],
        [0.3, 1.0, 0.3],
        [0.01, 0.01, 0.01],
    ])
    structured = _field(
        cold_shell_temperature, cold_shell_density, core_emitter)
    first = deterministic_axisymmetric_resonance_escape(
        uniform, LINE, source_quadrature_order=6,
        direction_quadrature_order=6, path_quadrature_order=12,
        frequency_quadrature_order=24)
    second = deterministic_axisymmetric_resonance_escape(
        structured, LINE, source_quadrature_order=6,
        direction_quadrature_order=6, path_quadrature_order=12,
        frequency_quadrature_order=24)
    repeated = deterministic_axisymmetric_resonance_escape(
        structured, LINE, source_quadrature_order=6,
        direction_quadrature_order=6, path_quadrature_order=12,
        frequency_quadrature_order=24)
    assert second.trapping_factor != pytest.approx(first.trapping_factor)
    assert repeated == second
    assert second.frequency_profile_normalization_minimum == pytest.approx(
        1.0, abs=2.0e-8)


def test_axisymmetric_field_is_immutable_and_fails_on_negative_density():
    shape = (3, 3)
    field = _field(
        np.full(shape, 500.0), np.full(shape, 2.0e20), np.ones(shape))
    with pytest.raises(ValueError):
        field.absorber_density_m3[0, 0] = 0.0
    bad = np.full(shape, 2.0e20)
    bad[0, 0] = -1.0
    with pytest.raises(ValueError):
        _field(np.full(shape, 500.0), bad, np.ones(shape))


def test_axisymmetric_transport_exposes_a_coarse_refined_gate():
    shape = (3, 3)
    field = _field(
        np.full(shape, 500.0), np.full(shape, 2.0e18), np.ones(shape))
    receipt = certify_deterministic_axisymmetric_resonance_escape(
        field,
        LINE,
        source_quadrature_order=4,
        direction_quadrature_order=4,
        path_quadrature_order=4,
        frequency_quadrature_order=8,
        relative_tolerance=0.2,
    )
    assert receipt.converged
    assert receipt.refined.source_quadrature_order == 8
    assert receipt.refined.direction_quadrature_order == 8
    assert receipt.refined.path_quadrature_order == 8
    assert receipt.refined.frequency_quadrature_order == 16
    assert receipt.relative_escape_change >= 0.0
    with pytest.raises(ValueError):
        certify_deterministic_axisymmetric_resonance_escape(
            field, LINE, refinement_factor=1)
