import numpy as np
import pytest

from petch.reactor_global import (
    PeriodicCurrentDensity,
    TurnerChabertCurrentDrivenSheath,
)


def _current(scale=1.0, evidence_kind="assumed"):
    return PeriodicCurrentDensity(
        fundamental_frequency_hz=13.56e6,
        harmonic_number=np.array([1]),
        sine_A_m2=np.array([-35.0 * scale]),
        cosine_A_m2=np.array([0.0]),
        source="manufactured Turner-Chabert single-frequency gate",
        evidence_kind=evidence_kind,
    )


def _sheath(scale=1.0, evidence_kind="assumed"):
    return TurnerChabertCurrentDrivenSheath(
        current=_current(scale, evidence_kind=evidence_kind),
        electron_temperature_eV=4.0,
        ion_mass_amu=39.948,
        sheath_edge_density_m3=2.0e16,
        phase_quadrature_count=2048,
    )


def test_single_frequency_recovers_turner_chabert_shape_and_xi():
    sheath = _sheath()
    phase = sheath.phase_grid_rad
    time = phase / (2.0 * np.pi * sheath.current.fundamental_frequency_hz)
    normalized_charge = np.asarray(sheath.normalized_charge(time))
    # A sinusoidal RF current produces a normalized integrated charge spanning
    # exactly [0, 1], and source equation (18) gives xi = 163/384.
    assert normalized_charge.min() == pytest.approx(0.0, abs=2.0e-15)
    assert normalized_charge.max() == pytest.approx(1.0, abs=2.0e-15)
    assert sheath.xi == pytest.approx(163.0 / 384.0, abs=2.0e-12)
    voltage = np.asarray(sheath.voltage(time))
    assert voltage.min() == pytest.approx(0.0, abs=2.0e-12)
    assert voltage.max() == pytest.approx(
        sheath.maximum_voltage_v, rel=0.0, abs=2.0e-12)


def test_charge_current_poisson_and_child_ledgers_close():
    sheath = _sheath()
    phase = sheath.phase_grid_rad
    time = phase / (2.0 * np.pi * sheath.current.fundamental_frequency_hz)
    dt = sheath.current.period_s / len(time)
    charge = np.asarray(sheath.current.charge_primitive_C_m2(time))
    current = np.asarray(sheath.current.current_density_A_m2(time))
    # Periodic centered difference independently checks dQ/dt = J.
    reconstructed = (np.roll(charge, -1) - np.roll(charge, 1)) / (2.0 * dt)
    assert np.max(np.abs(reconstructed - current)) / np.max(
        np.abs(current)) < 2.0e-6
    assert sheath.charge_voltage_relative_residual < 2.0e-14
    assert sheath.child_current_relative_residual < 2.0e-14


def test_moving_front_field_obeys_instantaneous_poisson_geometry():
    sheath = _sheath()
    frequency = sheath.current.fundamental_frequency_hz
    phase = sheath.phase_grid_rad
    time = phase / (2.0 * np.pi * frequency)
    normalized_charge = np.asarray(sheath.normalized_charge(time))
    selected = int(np.argmin(np.abs(normalized_charge - 0.4)))
    t = float(time[selected])
    front = float(sheath.electron_front_fraction(t))
    width = sheath.maximum_width_m
    assert sheath.electric_field_V_m(0.99 * front * width, t) == 0.0
    assert sheath.potential_drop_v(0.99 * front * width, t) == 0.0
    assert sheath.electric_field_V_m(width, t) > 0.0
    assert sheath.potential_drop_v(width, t) == pytest.approx(
        sheath.voltage(t), rel=2.0e-14, abs=2.0e-12)


def test_current_scale_power_laws_and_exact_jvp():
    base = _sheath()
    scale = 1.0001
    perturbed = _sheath(scale)
    assert perturbed.maximum_width_m / base.maximum_width_m == pytest.approx(
        scale ** 3, rel=2.0e-12)
    assert perturbed.maximum_voltage_v / base.maximum_voltage_v == pytest.approx(
        scale ** 4, rel=2.0e-12)
    direction = 0.37
    tangent = base.current_scale_jvp(direction)
    assert tangent.maximum_width_tangent_m == pytest.approx(
        3.0 * base.maximum_width_m * direction)
    assert tangent.maximum_voltage_tangent_v == pytest.approx(
        4.0 * base.maximum_voltage_v * direction)


def test_moving_sheath_ion_trajectory_is_deterministic_and_physical():
    sheath = _sheath()
    phases = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
    first = sheath.ion_impact_energies(
        phases, steps_per_period=128, steps_per_transit=192)
    second = sheath.ion_impact_energies(
        phases, steps_per_period=128, steps_per_transit=192)
    np.testing.assert_array_equal(first, second)
    assert np.all(np.isfinite(first))
    assert np.all(first >= 0.5 * sheath.electron_temperature_eV)
    assert np.ptp(first) > 0.0
    # A resolved moving front does not imply a generator-power or depth claim.
    receipt = sheath.certification_receipt()
    assert receipt["moving_electron_front_resolved"] is True
    assert receipt["time_dependent_poisson_field_resolved"] is True
    assert receipt["supports_generator_power_inversion"] is False
    assert receipt["supports_feature_depth"] is False


def test_fourier_export_preserves_voltage_and_evidence_gate():
    sheath = _sheath(evidence_kind="measured_sheath_current")
    waveform = sheath.voltage_fourier_projection(harmonic_count=24)
    phase = sheath.phase_grid_rad
    time = phase / (2.0 * np.pi * sheath.current.fundamental_frequency_hz)
    relative = np.linalg.norm(
        waveform.voltage(time) - sheath.voltage(time)
    ) / np.linalg.norm(sheath.voltage(time))
    assert relative < 2.0e-5
    assert waveform.evidence_kind == "validated_reactor_model"
    assert sheath.supports_predictive_boundary


def test_dual_frequency_waveform_is_deterministic_and_bounded():
    current = PeriodicCurrentDensity(
        fundamental_frequency_hz=1.0e6,
        harmonic_number=np.array([1, 40]),
        sine_A_m2=np.array([-25.0, -12.0]),
        cosine_A_m2=np.array([0.0, 4.0]),
        source="manufactured dual-frequency current",
    )
    sheath = TurnerChabertCurrentDrivenSheath(
        current=current,
        electron_temperature_eV=3.5,
        ion_mass_amu=39.948,
        sheath_edge_density_m3=1.0e16,
        phase_quadrature_count=4096,
    )
    time = sheath.phase_grid_rad / (2.0 * np.pi * 1.0e6)
    normalized = np.asarray(sheath.normalized_charge(time))
    voltage = np.asarray(sheath.voltage(time))
    assert np.all((normalized >= 0.0) & (normalized <= 1.0))
    assert np.all(voltage >= -2.0e-12)
    assert sheath.maximum_frequency_hz == pytest.approx(40.0e6)
    assert 0.0 < sheath.xi < 1.0


def test_invalid_current_evidence_and_empty_waveforms_fail_closed():
    with pytest.raises(ValueError, match="invalid periodic"):
        PeriodicCurrentDensity(
            fundamental_frequency_hz=1.0e6,
            harmonic_number=np.array([1]),
            sine_A_m2=np.array([0.0]),
            cosine_A_m2=np.array([0.0]),
            source="zero",
        )
    with pytest.raises(ValueError, match="invalid periodic"):
        PeriodicCurrentDensity(
            fundamental_frequency_hz=1.0e6,
            harmonic_number=np.array([1]),
            sine_A_m2=np.array([1.0]),
            cosine_A_m2=np.array([0.0]),
            source="bad evidence",
            evidence_kind="feature_depth_fit",
        )
