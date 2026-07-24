"""Rung-0 gates for the LER metrology twin: synthesis <-> estimation closure."""

import numpy as np
import pytest

from petch.ler_metrology import (
    EdgeStatistics,
    averaged_psd_nm3,
    fit_edge_statistics,
    height_height_correlation_nm2,
    mack_noise_floor_subtraction,
    palasantzas_psd_nm3,
    periodogram_psd_nm3,
    sigma_from_psd_nm,
    synthesize_edge_nm,
)

_DECLARED = EdgeStatistics(sigma_nm=2.0, correlation_length_nm=25.0,
                           roughness_exponent=0.75)


def test_palasantzas_normalization_integrates_to_sigma_squared():
    f = np.linspace(0.0, 5.0, 400000)
    psd = palasantzas_psd_nm3(f, _DECLARED)
    assert sigma_from_psd_nm(f, psd) == pytest.approx(_DECLARED.sigma_nm, rel=2e-3)


def test_synthesis_estimation_round_trip_recovers_declared_statistics():
    edges = np.stack([
        synthesize_edge_nm(_DECLARED, n_points=4096, spacing_nm=1.0, seed=s)
        for s in range(64)])
    frequencies, psd = averaged_psd_nm3(edges, spacing_nm=1.0)
    sigma = sigma_from_psd_nm(frequencies, psd)
    assert sigma == pytest.approx(_DECLARED.sigma_nm, rel=0.05)
    recovered = fit_edge_statistics(frequencies, psd)
    assert recovered.sigma_nm == pytest.approx(_DECLARED.sigma_nm, rel=0.05)
    assert recovered.roughness_exponent == pytest.approx(
        _DECLARED.roughness_exponent, abs=0.15)
    assert recovered.correlation_length_nm == pytest.approx(
        _DECLARED.correlation_length_nm, rel=0.5)


def test_synthesis_is_deterministic_per_seed():
    a = synthesize_edge_nm(_DECLARED, n_points=1024, spacing_nm=1.0, seed=7)
    b = synthesize_edge_nm(_DECLARED, n_points=1024, spacing_nm=1.0, seed=7)
    c = synthesize_edge_nm(_DECLARED, n_points=1024, spacing_nm=1.0, seed=8)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_hhcf_saturates_at_twice_sigma_squared():
    edges = [synthesize_edge_nm(_DECLARED, n_points=8192, spacing_nm=1.0, seed=s)
             for s in range(8)]
    saturations = []
    for edge in edges:
        lags, hhcf = height_height_correlation_nm2(edge, spacing_nm=1.0)
        tail = hhcf[lags > 8.0 * _DECLARED.correlation_length_nm]
        saturations.append(float(np.mean(tail)))
    assert float(np.mean(saturations)) == pytest.approx(
        2.0 * _DECLARED.sigma_nm ** 2, rel=0.25)


def test_noise_floor_subtraction_recovers_clean_sigma():
    edges = np.stack([
        synthesize_edge_nm(_DECLARED, n_points=4096, spacing_nm=1.0, seed=s)
        for s in range(64)])
    frequencies, clean = averaged_psd_nm3(edges, spacing_nm=1.0)
    noisy = clean + 0.05  # flat CD-SEM-style noise floor (nm^3)
    corrected, floor = mack_noise_floor_subtraction(frequencies, noisy)
    assert floor == pytest.approx(0.05, rel=0.3)
    assert sigma_from_psd_nm(frequencies, corrected) == pytest.approx(
        sigma_from_psd_nm(frequencies, clean), rel=0.1)


def test_invalid_statistics_refused():
    with pytest.raises(ValueError):
        EdgeStatistics(-1.0, 25.0, 0.75)
    with pytest.raises(ValueError):
        EdgeStatistics(2.0, 0.0, 0.75)
    with pytest.raises(ValueError):
        EdgeStatistics(2.0, 25.0, 1.5)
