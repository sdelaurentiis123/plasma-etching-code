"""Line-edge-roughness metrology twin (LER modality Rung 0).

Implements the measurement side of the PSD-transfer architecture
(RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md): synthesis of self-affine test
edges with declared (sigma, correlation length xi, roughness exponent alpha),
periodogram PSD and height-height correlation estimators, and the Mack
noise-floor subtraction used on CD-SEM data. Everything is deterministic
given an explicit seed; nothing here touches the feature engine.

Conventions: edge displacement u(y) in nm sampled uniformly along the line
(spacing dy in nm); spatial frequency f in 1/nm; one-sided PSD normalized so
sigma^2 = integral PSD(f) df.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class EdgeStatistics:
    """Declared self-affine edge statistics (the standard LER triplet)."""

    sigma_nm: float
    correlation_length_nm: float
    roughness_exponent: float

    def __post_init__(self):
        if (not np.isfinite(self.sigma_nm) or self.sigma_nm < 0.0
                or not np.isfinite(self.correlation_length_nm)
                or self.correlation_length_nm <= 0.0
                or not np.isfinite(self.roughness_exponent)
                or not 0.0 < self.roughness_exponent <= 1.0):
            raise ValueError("invalid edge statistics")


def palasantzas_psd_nm3(frequency_per_nm, statistics: EdgeStatistics):
    """One-sided Palasantzas self-affine PSD (the standard LER model form).

    PSD(f) = C / (1 + (2 pi f xi)^2)^(alpha + 1/2), with C fixed analytically
    so that the integral over f in [0, inf) equals sigma^2:
    integral df / (1+(2 pi f xi)^2)^(a+1/2) = B(1/2, a) / (4 pi xi) * 2
    where B is the Beta function — giving C = sigma^2 * 4 pi xi / (2 B(1/2, a))
    ... implemented via the Beta function directly, no numeric fit.
    """
    from scipy.special import beta as beta_function

    f = np.asarray(frequency_per_nm, dtype=float)
    sigma = statistics.sigma_nm
    xi = statistics.correlation_length_nm
    alpha = statistics.roughness_exponent
    # integral_0^inf dt (1+t^2)^-(alpha+1/2) = B(1/2, alpha)/2, so with
    # t = 2 pi f xi the PSD integral is C * B(1/2, alpha) / (4 pi xi) = sigma^2.
    normalization = (4.0 * np.pi * xi
                     / beta_function(0.5, alpha))
    return (sigma ** 2 * normalization
            / (1.0 + (2.0 * np.pi * f * xi) ** 2) ** (alpha + 0.5))


def synthesize_edge_nm(statistics: EdgeStatistics, *, n_points: int,
                       spacing_nm: float, seed: int):
    """Deterministic spectral synthesis of a periodic self-affine edge.

    Draws independent Gaussian Fourier amplitudes with variance proportional
    to the declared PSD, so ensemble statistics converge to the declared
    triplet exactly; a single realization carries the usual finite-length
    scatter (quantified in tests, not hidden).
    """
    if n_points < 8 or spacing_nm <= 0.0:
        raise ValueError("edge synthesis requires n_points >= 8 and positive spacing")
    rng = np.random.default_rng(int(seed))
    length_nm = n_points * spacing_nm
    frequencies = np.fft.rfftfreq(n_points, d=spacing_nm)
    psd = palasantzas_psd_nm3(frequencies, statistics)
    psd[0] = 0.0  # zero-mean edge; DC carries no roughness
    # Parseval-consistent with periodogram_psd_nm3: the estimator computes
    # 2 dy^2 |rfft(u)|^2 / L, so E|S_k|^2 must equal PSD_k L / (2 dy^2),
    # i.e. per-component variance PSD_k L / (4 dy^2).
    variance = psd * length_nm / (4.0 * spacing_nm ** 2)
    real = rng.standard_normal(frequencies.size)
    imaginary = rng.standard_normal(frequencies.size)
    spectrum = np.sqrt(np.maximum(variance, 0.0)) * (real + 1j * imaginary)
    spectrum[0] = 0.0
    if n_points % 2 == 0:
        spectrum[-1] = np.sqrt(2.0) * spectrum[-1].real
    edge = np.fft.irfft(spectrum, n=n_points)
    return edge


def periodogram_psd_nm3(edge_nm, *, spacing_nm: float):
    """One-sided periodogram PSD estimate of a (near-)periodic sampled edge."""
    edge = np.asarray(edge_nm, dtype=float)
    if edge.ndim != 1 or edge.size < 8 or spacing_nm <= 0.0:
        raise ValueError("PSD estimation requires a 1-D edge of >= 8 samples")
    centered = edge - edge.mean()
    spectrum = np.fft.rfft(centered)
    length_nm = edge.size * spacing_nm
    psd = 2.0 * spacing_nm ** 2 * np.abs(spectrum) ** 2 / length_nm
    psd[0] = 0.0
    if edge.size % 2 == 0:
        psd[-1] *= 0.5
    frequencies = np.fft.rfftfreq(edge.size, d=spacing_nm)
    return frequencies, psd


def averaged_psd_nm3(edges_nm, *, spacing_nm: float):
    """Ensemble-averaged periodogram over many sampled edges (rows)."""
    edges = np.atleast_2d(np.asarray(edges_nm, dtype=float))
    accumulated = None
    for row in edges:
        frequencies, psd = periodogram_psd_nm3(row, spacing_nm=spacing_nm)
        accumulated = psd if accumulated is None else accumulated + psd
    return frequencies, accumulated / edges.shape[0]


def sigma_from_psd_nm(frequencies_per_nm, psd_nm3):
    """sigma = sqrt(integral PSD df) by trapezoid; inverse of the synthesis norm."""
    return float(np.sqrt(np.trapz(np.asarray(psd_nm3, dtype=float),
                                  np.asarray(frequencies_per_nm, dtype=float))))


def height_height_correlation_nm2(edge_nm, *, spacing_nm: float,
                                  max_lag_points: int | None = None):
    """HHCF(r) = <(u(y+r) - u(y))^2> over the periodic sample."""
    edge = np.asarray(edge_nm, dtype=float)
    n = edge.size
    lags = np.arange(1, (n // 2 if max_lag_points is None
                         else min(max_lag_points, n // 2)) + 1)
    values = np.empty(lags.size)
    for index, lag in enumerate(lags):
        difference = np.roll(edge, -int(lag)) - edge
        values[index] = float(np.mean(difference ** 2))
    return lags * spacing_nm, values


def mack_noise_floor_subtraction(frequencies_per_nm, psd_nm3, *,
                                 noise_fraction=0.2):
    """Estimate and subtract the flat CD-SEM noise floor (Mack's method).

    The floor is the median PSD over the highest ``noise_fraction`` of the
    frequency range, where genuine self-affine PSDs have rolled off far below
    it; subtraction is clamped at zero. Returns (corrected_psd, floor).
    """
    frequencies = np.asarray(frequencies_per_nm, dtype=float)
    psd = np.asarray(psd_nm3, dtype=float)
    if not 0.0 < noise_fraction < 1.0:
        raise ValueError("noise_fraction must be in (0, 1)")
    cut = frequencies >= frequencies.max() * (1.0 - noise_fraction)
    if not np.any(cut):
        raise ValueError("no high-frequency band available for the noise floor")
    floor = float(np.median(psd[cut]))
    return np.maximum(psd - floor, 0.0), floor


def fit_edge_statistics(frequencies_per_nm, psd_nm3) -> EdgeStatistics:
    """Recover (sigma, xi, alpha) from a measured PSD.

    sigma from the integral; alpha from the log-log high-frequency slope
    (slope = -(2 alpha + 1)); xi from the half-power point of the declared
    Palasantzas form. Pure closed-form estimators — no iterative fitting.
    """
    frequencies = np.asarray(frequencies_per_nm, dtype=float)
    psd = np.asarray(psd_nm3, dtype=float)
    positive = (frequencies > 0.0) & (psd > 0.0)
    frequencies = frequencies[positive]
    psd = psd[positive]
    if frequencies.size < 8:
        raise ValueError("edge-statistics fit requires >= 8 positive PSD points")
    sigma = float(np.sqrt(np.trapz(psd, frequencies)))
    # High-frequency slope over the top decade.
    high = frequencies >= frequencies.max() / 10.0
    slope = np.polyfit(np.log(frequencies[high]), np.log(psd[high]), 1)[0]
    alpha = float(np.clip((-slope - 1.0) / 2.0, 0.05, 1.0))
    # Half-power point: PSD(f_c) = PSD(0)/2^(alpha+1/2) at 2 pi f_c xi = 1.
    plateau = float(np.median(psd[frequencies <= frequencies.max() * 1e-2])
                    if np.any(frequencies <= frequencies.max() * 1e-2)
                    else psd[0])
    target = plateau / 2.0 ** (alpha + 0.5)
    crossing = np.flatnonzero(psd <= target)
    f_c = float(frequencies[crossing[0]]) if crossing.size else float(
        frequencies[psd.size // 2])
    xi = float(1.0 / (2.0 * np.pi * f_c))
    return EdgeStatistics(sigma, xi, alpha)
