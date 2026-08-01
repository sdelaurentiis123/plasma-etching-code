"""Rung-1 gates for the LER spectral transfer layer.

Every gate is a round trip through a KNOWN analytic transfer function: push a
synthetic Palasantzas ensemble through T(f) (plus known additive intrinsic
noise), then recover T(f), the coherence and the intrinsic PSD from the pairs
alone. Nothing is fitted anywhere in the module under test.
"""

import numpy as np
import pytest

from petch.ler_metrology import (
    EdgeStatistics,
    periodogram_psd_nm3,
    sigma_from_psd_nm,
    synthesize_edge_nm,
)
from petch.ler_transfer import (
    apply_transfer_to_edges,
    cross_psd_nm3,
    estimate_transfer,
    predict_output_psd_nm3,
    predict_output_statistics,
    spectral_noise_edge_nm,
)

_DECLARED = EdgeStatistics(sigma_nm=2.0, correlation_length_nm=25.0,
                           roughness_exponent=0.75)
_N_POINTS = 1024
_SPACING = 1.0
_REALIZATIONS = 256
_TRANSFER_XI = 12.0          # nm: the known low-pass length
_INTRINSIC_PSD = 0.02        # nm^3: known flat etch-added floor


def _known_transfer(frequency_per_nm):
    """Analytic low-pass: T(f) = 1 / (1 + (2 pi f xi_t)^2), T(0) = 1."""
    return 1.0 / (1.0 + (2.0 * np.pi * np.asarray(frequency_per_nm)
                         * _TRANSFER_XI) ** 2)


def _frequencies():
    return np.fft.rfftfreq(_N_POINTS, d=_SPACING)


def _input_ensemble(n=_REALIZATIONS, offset=0):
    return np.stack([
        synthesize_edge_nm(_DECLARED, n_points=_N_POINTS,
                           spacing_nm=_SPACING, seed=offset + s)
        for s in range(n)])


@pytest.fixture(scope="module")
def pairs():
    """Input ensemble, its noisy output, and the estimate recovered from them."""
    frequencies = _frequencies()
    intrinsic = np.full(frequencies.size, _INTRINSIC_PSD)
    intrinsic[0] = 0.0
    edges_in = _input_ensemble()
    edges_out = apply_transfer_to_edges(
        edges_in, spacing_nm=_SPACING, transfer=_known_transfer,
        intrinsic_psd_nm3=intrinsic, seed=99)
    estimate = estimate_transfer(edges_in, edges_out, spacing_nm=_SPACING)
    return edges_in, edges_out, estimate


def test_known_transfer_and_intrinsic_are_recovered(pairs):
    """Self-consistency: recover both |T|^2 and the added PSD from pairs only.

    |T|^2 is gated where the input still explains the output (coherence > 0.9,
    33 bins): max 6%, mean 2%. The intrinsic floor is gated in the complement
    (bins where the added power exceeds the transferred power tenfold): band
    mean within 2% of truth, per-bin mean within 10% (per-bin scatter is the
    irreducible chi-square spread of a 256-realization periodogram).
    """
    _, _, estimate = pairs
    frequencies = estimate.frequencies_per_nm
    truth = _known_transfer(frequencies) ** 2

    coherent = estimate.coherence > 0.9
    assert coherent.sum() >= 20
    relative = np.abs(estimate.magnitude_squared[coherent]
                      - truth[coherent]) / truth[coherent]
    assert relative.max() < 0.06
    assert relative.mean() < 0.02

    transferred = truth * estimate.input_psd_nm3
    tail = (frequencies > 0.0) & (_INTRINSIC_PSD > 10.0 * transferred)
    assert tail.sum() >= 100
    recovered = estimate.intrinsic_psd_nm3[tail]
    assert recovered.mean() == pytest.approx(_INTRINSIC_PSD, rel=0.02)
    per_bin = np.abs(recovered - _INTRINSIC_PSD) / _INTRINSIC_PSD
    assert per_bin.mean() < 0.10


def test_unity_transfer_limit_is_exact():
    """Output identical to input: T == 1, gamma^2 == 1, intrinsic == 0."""
    edges = _input_ensemble(n=16)
    estimate = estimate_transfer(edges, edges, spacing_nm=_SPACING)
    powered = estimate.input_psd_nm3 > 0.0
    assert np.abs(estimate.transfer[powered] - 1.0).max() < 1e-12
    assert np.array_equal(estimate.coherence[powered],
                          np.ones(int(powered.sum())))
    scale = estimate.output_psd_nm3.max()
    assert np.abs(estimate.intrinsic_psd_nm3).max() < 1e-12 * scale


def test_zero_transfer_limit_attributes_everything_to_intrinsic():
    """Output independent of input: |T|^2 PSD_in vanishes, intrinsic == PSD_out."""
    frequencies = _frequencies()
    edges_in = _input_ensemble()
    noise = np.stack([
        spectral_noise_edge_nm(np.full(frequencies.size, 0.05),
                               n_points=_N_POINTS, spacing_nm=_SPACING,
                               rng=np.random.default_rng(1000 + i))
        for i in range(_REALIZATIONS)])
    estimate = estimate_transfer(edges_in, noise, spacing_nm=_SPACING)
    powered = (estimate.input_psd_nm3 > 0.0) & (estimate.output_psd_nm3 > 0.0)
    explained = (estimate.magnitude_squared[powered]
                 * estimate.input_psd_nm3[powered]
                 / estimate.output_psd_nm3[powered])
    assert explained.mean() < 0.02
    residual = (estimate.intrinsic_psd_nm3[powered]
                / estimate.output_psd_nm3[powered])
    assert residual.mean() == pytest.approx(1.0, rel=0.02)


def test_coherence_separates_transferred_from_intrinsic_power():
    """gamma^2 -> 1 with no added noise, gamma^2 -> 0 for pure noise."""
    frequencies = _frequencies()
    edges_in = _input_ensemble(n=64)
    clean_out = apply_transfer_to_edges(edges_in, spacing_nm=_SPACING,
                                        transfer=_known_transfer)
    clean = estimate_transfer(edges_in, clean_out, spacing_nm=_SPACING)
    powered = clean.input_psd_nm3 > 0.0
    assert clean.coherence[powered].min() > 0.99
    # ... and the noiseless transfer is recovered essentially exactly.
    truth = _known_transfer(frequencies) ** 2
    assert np.abs(clean.magnitude_squared[powered]
                  - truth[powered]).max() < 1e-9

    noise = np.stack([
        spectral_noise_edge_nm(np.full(frequencies.size, 0.05),
                               n_points=_N_POINTS, spacing_nm=_SPACING,
                               rng=np.random.default_rng(2000 + i))
        for i in range(64)])
    incoherent = estimate_transfer(edges_in, noise, spacing_nm=_SPACING)
    assert incoherent.coherence[powered].mean() < 0.02
    assert incoherent.coherence[powered].max() < 0.20


def test_parseval_consistency_of_predicted_sigma(pairs):
    """Predicted PSD must carry the realized ensemble's power.

    Against the measured output spectrum the prediction closes to 1e-3 (it is
    the same power, reassembled); against the realized edge ensemble it closes
    to 3% (trapezoid quadrature on the discrete grid plus the M/(M-1) residual
    correction).
    """
    _, edges_out, estimate = pairs
    frequencies = estimate.frequencies_per_nm
    predicted = predict_output_psd_nm3(frequencies, estimate.input_psd_nm3,
                                       estimate)
    sigma_predicted = sigma_from_psd_nm(frequencies, predicted)
    sigma_measured = sigma_from_psd_nm(frequencies, estimate.output_psd_nm3)
    assert sigma_predicted == pytest.approx(sigma_measured, rel=1e-3)
    sigma_realized = float(np.std(edges_out, axis=1).mean())
    assert sigma_predicted == pytest.approx(sigma_realized, rel=0.03)


def test_preregistered_transfer_signs(pairs):
    """Design-doc Rung-1 signs: T -> 1 at low f, sigma shrinks, xi grows."""
    _, _, estimate = pairs
    assert estimate.magnitude_squared[1] == pytest.approx(1.0, abs=0.02)
    _, _, predicted = predict_output_statistics(_DECLARED, estimate,
                                                include_intrinsic=False)
    assert predicted.sigma_nm < _DECLARED.sigma_nm
    assert (predicted.correlation_length_nm
            > _DECLARED.correlation_length_nm)


def test_naive_ratio_is_biased_where_the_cross_estimator_is_not(pairs):
    """The bias the H1 estimator exists to avoid, measured explicitly.

    In the noise-dominated band the true |T|^2 is ~4e-6: the naive periodogram
    ratio overstates it by ~2e5 (it attributes the whole etch-added floor to
    transfer, by construction), the cross-spectral estimate by ~9e2 — two
    orders of magnitude less biased, and it credits under 1% of the output
    power to the input instead of all of it.
    """
    _, _, estimate = pairs
    frequencies = estimate.frequencies_per_nm
    truth = _known_transfer(frequencies) ** 2
    powered = estimate.input_psd_nm3 > 0.0
    naive = np.zeros_like(truth)
    naive[powered] = (estimate.output_psd_nm3[powered]
                      / estimate.input_psd_nm3[powered])
    noise_band = powered & (frequencies > 0.2)
    assert noise_band.sum() >= 100
    naive_bias = naive[noise_band].mean() / truth[noise_band].mean()
    cross_bias = (estimate.magnitude_squared[noise_band].mean()
                  / truth[noise_band].mean())
    assert naive_bias > 1.0e4
    assert cross_bias < naive_bias / 100.0
    explained = (estimate.magnitude_squared[noise_band]
                 * estimate.input_psd_nm3[noise_band]
                 / estimate.output_psd_nm3[noise_band])
    assert explained.mean() < 0.01


def test_cross_psd_matches_the_metrology_normalization():
    """S_xx must reproduce the Rung-0 periodogram, so ratios are dimensionless."""
    edge = synthesize_edge_nm(_DECLARED, n_points=_N_POINTS,
                              spacing_nm=_SPACING, seed=3)
    frequencies, auto = cross_psd_nm3(edge, edge, spacing_nm=_SPACING)
    reference_f, reference = periodogram_psd_nm3(edge, spacing_nm=_SPACING)
    assert np.array_equal(frequencies, reference_f)
    assert auto.real == pytest.approx(reference, rel=1e-12, abs=1e-18)
    assert np.abs(auto.imag).max() < 1e-18


def test_invalid_inputs_are_refused():
    edges = _input_ensemble(n=4)
    with pytest.raises(ValueError):
        cross_psd_nm3(edges, edges[:2], spacing_nm=_SPACING)
    with pytest.raises(ValueError):
        cross_psd_nm3(edges, edges, spacing_nm=0.0)
    with pytest.raises(ValueError):
        cross_psd_nm3(np.zeros((2, 4)), np.zeros((2, 4)), spacing_nm=1.0)
    estimate = estimate_transfer(edges, edges, spacing_nm=_SPACING)
    with pytest.raises(ValueError):
        predict_output_psd_nm3(np.linspace(0.0, 0.5, 16),
                               np.ones(8), estimate)
    with pytest.raises(ValueError):
        apply_transfer_to_edges(edges, spacing_nm=_SPACING,
                                transfer=_known_transfer,
                                intrinsic_psd_nm3=np.ones(7))
