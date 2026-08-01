"""LER spectral transfer layer (LER modality Rung 1).

Implements the Naulleau/Gallatin transfer algebra that the whole modality is
built on (RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md sect. 3a):

    PSD_out(f) = |T(f)|^2 * PSD_in(f) + PSD_intrinsic(f)

with two directions of use:

* **estimation** — given paired input/output edge ensembles (from the engine's
  perturbation probe, or from an experiment's ADI/AEI pairs), recover the
  complex transfer function T(f), the magnitude-squared coherence gamma^2(f),
  and the intrinsic (etch-added) PSD;
* **prediction** — given a measured T(f) + intrinsic PSD and a measured or
  Palasantzas input PSD, forward-apply them and report the predicted output
  triplet (sigma, xi, alpha) through the Rung-0 metrology estimators.

Nothing here is fitted. T(f) is *measured* — from simulation pairs or from
experimental pairs — never assumed or parameterized; this module is only the
estimation and prediction machinery around it.

Edge convention (mask edge -> etched edge)
-----------------------------------------
An edge is a 1-D displacement signal u(y) in nm from the nominal straight
line, sampled uniformly along the line direction y with spacing dy in nm;
positive u points outward (feature locally wider). The *input* edge is the
mask/resist edge before etch (ADI); the *output* edge is the etched-feature
edge after etch (AEI), taken at a declared measurement plane (top-down, or a
sidewall height), sampled on the same y grid over the same box length. T(f) is
therefore the mask-to-etched-edge transfer at that plane: |T| < 1 means the
etch smooths that frequency, |T| > 1 means it amplifies it, and arg T carries
any lateral shift of the pattern (zero for a translationally symmetric etch).

Spectral conventions follow ``ler_metrology``: frequency f in 1/nm, one-sided
spectra normalized so sigma^2 = integral PSD(f) df, and the DC bin carries no
roughness by construction.

Estimator notes (the bias handling)
-----------------------------------
The naive periodogram ratio PSD_out/PSD_in is biased *upward* by exactly the
intrinsic term (it attributes etch-added roughness to transfer). The
cross-spectral H1 estimator used here,

    T(f) = S_xy(f) / S_xx(f),

is asymptotically unbiased whenever the intrinsic noise is uncorrelated with
the input edge, because the uncorrelated part averages out of the cross
spectrum. Two finite-ensemble corrections are applied when ``bias_correct``
is set (the default): the magnitude-squared coherence of M realizations is
biased high by ~1/M, corrected as (M*g - 1)/(M - 1); and the residual
non-coherent power is biased low by the one complex degree of freedom spent
fitting T per bin, corrected by M/(M - 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ler_metrology import (
    EdgeStatistics,
    fit_edge_statistics,
    palasantzas_psd_nm3,
    sigma_from_psd_nm,
)


def _as_ensemble(edges_nm):
    """Rows = realizations; a single edge is a one-row ensemble."""
    edges = np.atleast_2d(np.asarray(edges_nm, dtype=float))
    if edges.ndim != 2 or edges.shape[1] < 8:
        raise ValueError("edge ensemble must be rows of >= 8 samples")
    return edges


def cross_psd_nm3(edges_a_nm, edges_b_nm, *, spacing_nm: float):
    """Ensemble-averaged one-sided cross-PSD S_ab(f) = <conj(A) B>.

    Normalized identically to ``ler_metrology.periodogram_psd_nm3`` (same
    2 dy^2 / L scaling, same zeroed DC bin, same halved Nyquist bin), so that
    ``cross_psd_nm3(x, x)`` is the auto-PSD and ratios of these spectra are
    dimensionless.
    """
    a = _as_ensemble(edges_a_nm)
    b = _as_ensemble(edges_b_nm)
    if a.shape != b.shape:
        raise ValueError("paired ensembles must have matching shapes")
    if spacing_nm <= 0.0:
        raise ValueError("spacing must be positive")
    n_points = a.shape[1]
    length_nm = n_points * spacing_nm
    spectrum_a = np.fft.rfft(a - a.mean(axis=1, keepdims=True), axis=1)
    spectrum_b = np.fft.rfft(b - b.mean(axis=1, keepdims=True), axis=1)
    cross = (2.0 * spacing_nm ** 2 / length_nm) * np.conj(spectrum_a) * spectrum_b
    cross = cross.mean(axis=0)
    cross[0] = 0.0
    if n_points % 2 == 0:
        cross[-1] *= 0.5
    frequencies = np.fft.rfftfreq(n_points, d=spacing_nm)
    return frequencies, cross


@dataclass(frozen=True)
class TransferEstimate:
    """A measured mask-to-etched-edge transfer function and its residual.

    ``transfer`` is complex T(f); ``coherence`` is gamma^2(f) in [0, 1] (1 =
    every bit of output power at that frequency is explained by the input, 0 =
    none of it is); ``intrinsic_psd_nm3`` is the etch-added power the input
    cannot explain. The DC bin is carried for grid alignment but holds no
    information (zeroed by the metrology convention).
    """

    frequencies_per_nm: np.ndarray
    transfer: np.ndarray
    coherence: np.ndarray
    intrinsic_psd_nm3: np.ndarray
    input_psd_nm3: np.ndarray
    output_psd_nm3: np.ndarray
    n_realizations: int

    @property
    def magnitude_squared(self):
        """|T(f)|^2 — the quantity that multiplies the input PSD."""
        return np.abs(self.transfer) ** 2

    def resample(self, frequencies_per_nm):
        """Interpolate (|T|^2, intrinsic PSD) onto another frequency grid.

        Linear in f, clamped at the measured band edges; the informationless
        DC bin is dropped first so the low-frequency clamp uses the lowest
        *measured* frequency (the plateau), not a placeholder.
        """
        target = np.asarray(frequencies_per_nm, dtype=float)
        measured = self.frequencies_per_nm[1:]
        magnitude = self.magnitude_squared[1:]
        intrinsic = self.intrinsic_psd_nm3[1:]
        return (np.interp(target, measured, magnitude),
                np.interp(target, measured, intrinsic))


def estimate_transfer(input_edges_nm, output_edges_nm, *, spacing_nm: float,
                      bias_correct: bool = True) -> TransferEstimate:
    """Recover T(f), gamma^2(f) and the intrinsic PSD from paired ensembles.

    ``input_edges_nm`` and ``output_edges_nm`` are matched row-for-row: row m
    of the output is what the etch produced from row m of the input, on the
    same y grid. See the module docstring for the estimator and its bias
    corrections.
    """
    frequencies, s_xx = cross_psd_nm3(input_edges_nm, input_edges_nm,
                                      spacing_nm=spacing_nm)
    _, s_yy = cross_psd_nm3(output_edges_nm, output_edges_nm,
                            spacing_nm=spacing_nm)
    _, s_xy = cross_psd_nm3(input_edges_nm, output_edges_nm,
                            spacing_nm=spacing_nm)
    auto_in = s_xx.real
    auto_out = s_yy.real
    realizations = int(_as_ensemble(input_edges_nm).shape[0])

    powered = auto_in > 0.0
    transfer = np.zeros(frequencies.size, dtype=complex)
    transfer[powered] = s_xy[powered] / auto_in[powered]

    coherence = np.zeros(frequencies.size)
    both = powered & (auto_out > 0.0)
    coherence[both] = np.abs(s_xy[both]) ** 2 / (auto_in[both] * auto_out[both])
    coherence = np.clip(coherence, 0.0, 1.0)

    # Residual (non-coherent) output power: what the input cannot explain.
    intrinsic = auto_out - np.abs(transfer) ** 2 * auto_in

    if bias_correct and realizations > 1:
        coherence = np.clip(
            (realizations * coherence - 1.0) / (realizations - 1.0), 0.0, 1.0)
        intrinsic = intrinsic * realizations / (realizations - 1.0)

    intrinsic = np.maximum(intrinsic, 0.0)
    intrinsic[0] = 0.0
    return TransferEstimate(frequencies, transfer, coherence, intrinsic,
                            auto_in, auto_out, realizations)


def predict_output_psd_nm3(frequencies_per_nm, input_psd_nm3,
                           estimate: TransferEstimate, *,
                           include_intrinsic: bool = True):
    """Forward-apply a measured transfer: |T|^2 PSD_in + PSD_intrinsic."""
    frequencies = np.asarray(frequencies_per_nm, dtype=float)
    psd_in = np.asarray(input_psd_nm3, dtype=float)
    if psd_in.shape != frequencies.shape:
        raise ValueError("input PSD must match the frequency grid")
    magnitude, intrinsic = estimate.resample(frequencies)
    predicted = magnitude * psd_in
    if include_intrinsic:
        predicted = predicted + intrinsic
    return predicted


def predict_output_statistics(input_statistics: EdgeStatistics,
                              estimate: TransferEstimate, *,
                              frequencies_per_nm=None,
                              include_intrinsic: bool = True):
    """Predict the post-etch (sigma, xi, alpha) for a Palasantzas input edge.

    Returns ``(frequencies, predicted_psd, EdgeStatistics)`` — the triplet is
    read back with the Rung-0 estimators, so a prediction and a measurement
    are always reduced by the same code.
    """
    if frequencies_per_nm is None:
        frequencies_per_nm = estimate.frequencies_per_nm
    frequencies = np.asarray(frequencies_per_nm, dtype=float)
    psd_in = palasantzas_psd_nm3(frequencies, input_statistics)
    psd_out = predict_output_psd_nm3(frequencies, psd_in, estimate,
                                     include_intrinsic=include_intrinsic)
    return frequencies, psd_out, fit_edge_statistics(frequencies, psd_out)


def predicted_sigma_nm(frequencies_per_nm, psd_nm3):
    """sigma of a predicted PSD (Parseval), via the Rung-0 integrator."""
    return sigma_from_psd_nm(frequencies_per_nm, psd_nm3)


def _transfer_on_grid(transfer, frequencies_per_nm):
    """Evaluate a callable T(f), or interpolate a tabulated (f, T) pair."""
    frequencies = np.asarray(frequencies_per_nm, dtype=float)
    if callable(transfer):
        return np.asarray(transfer(frequencies), dtype=complex)
    grid, values = transfer
    grid = np.asarray(grid, dtype=float)
    values = np.asarray(values, dtype=complex)
    return (np.interp(frequencies, grid, values.real)
            + 1j * np.interp(frequencies, grid, values.imag))


def spectral_noise_edge_nm(psd_nm3, *, n_points: int, spacing_nm: float,
                           rng: np.random.Generator):
    """Draw an edge realization with a prescribed one-sided PSD.

    Uses exactly the ``ler_metrology.synthesize_edge_nm`` normalization, so a
    noise realization drawn here and a synthetic edge drawn there are on the
    same Parseval footing.
    """
    psd = np.asarray(psd_nm3, dtype=float)
    length_nm = n_points * spacing_nm
    variance = psd * length_nm / (4.0 * spacing_nm ** 2)
    spectrum = np.sqrt(np.maximum(variance, 0.0)) * (
        rng.standard_normal(psd.size) + 1j * rng.standard_normal(psd.size))
    spectrum[0] = 0.0
    if n_points % 2 == 0:
        spectrum[-1] = np.sqrt(2.0) * spectrum[-1].real
    return np.fft.irfft(spectrum, n=n_points)


def apply_transfer_to_edges(input_edges_nm, *, spacing_nm: float, transfer,
                            intrinsic_psd_nm3=None, seed: int | None = None):
    """Push an input ensemble through a known T(f) plus optional added noise.

    ``transfer`` is a callable f -> complex or a tabulated ``(frequencies,
    values)`` pair. This is the forward operator whose inverse
    ``estimate_transfer`` performs; it exists so a known transfer can be
    round-tripped (the self-consistency gate) and so a measured transfer can
    generate output-edge realizations for downstream metrology.
    """
    edges = _as_ensemble(input_edges_nm)
    n_points = edges.shape[1]
    frequencies = np.fft.rfftfreq(n_points, d=spacing_nm)
    kernel = _transfer_on_grid(transfer, frequencies)
    spectrum = np.fft.rfft(edges - edges.mean(axis=1, keepdims=True), axis=1)
    filtered = np.fft.irfft(spectrum * kernel[None, :], n=n_points, axis=1)
    if intrinsic_psd_nm3 is not None:
        rng = np.random.default_rng(0 if seed is None else int(seed))
        psd = np.asarray(intrinsic_psd_nm3, dtype=float)
        if psd.shape != frequencies.shape:
            raise ValueError("intrinsic PSD must be on the edge rfft grid")
        noise = np.stack([
            spectral_noise_edge_nm(psd, n_points=n_points,
                                   spacing_nm=spacing_nm, rng=rng)
            for _ in range(edges.shape[0])])
        filtered = filtered + noise
    return filtered
