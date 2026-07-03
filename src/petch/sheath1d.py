"""RF-sheath arrival samplers for the 2-D Hwang-Giapis charging solver.

This is a compact source-boundary model, not a reactor model. It exists to remove the most
damaging shortcut in ``charging2d``: independent analytic draws for ion energy, ion angle, electron
energy, electron angle, and RF phase. The public function returns particles already at the feature
mouth, so the in-feature trace still owns all surface charging and Laplace-field physics.
"""
from __future__ import annotations

import numpy as np


def _burst_phase(rng, n, V_dc, V_rf, Te):
    phi = np.linspace(0.0, 2.0 * np.pi, 1440)
    w = np.exp(-(V_dc + V_rf * np.sin(phi)) / Te)
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return np.interp(rng.uniform(0.0, 1.0, n), cdf, phi)


def sample_sheath_source(kind, n, rng, Te=4.0, V_dc=37.0, V_rf=30.0,
                         iadf_hwhm_deg=4.3, cos_power=0.6,
                         boundary_um=3.7, sheath_um=89.0,
                         n_phase=7, ion_spread_floor=0.12):
    """Sample feature-mouth arrivals from a reduced RF-sheath source.

    ``kind`` is ``"ion"`` or ``"electron"``. Returns ``E0, theta, sB`` where ``E0`` is kinetic
    energy in the existing charging2d convention, ``theta`` is the in-plane angle from vertical, and
    ``sB`` is the residual sheath-barrier multiplier used by electrons inside the feature.

    Ions: finite-transit RF averaging over nearby phases gives a joint energy-angle distribution.
    It preserves the published HWHM at the mean sheath voltage but avoids independent marginals.

    Electrons: sheath-crossing vertical residual energy is exponential; transverse energy is kept
    but the feature-mouth angle is represented as a projected 3-D arrival to avoid the earlier
    over-wide purely in-plane source. The residual lower-sheath barrier remains explicit in ``sB``.
    """
    if kind not in ("ion", "electron"):
        raise ValueError(f"unknown sheath source kind: {kind}")

    if kind == "ion":
        phi0 = rng.uniform(0.0, 2.0 * np.pi, n)
        # Ions sample the RF field over a finite transit. The phase window is intentionally set by
        # the sheath voltage ratio only, not fit to HG data: high RF modulation -> broader window.
        span = 0.5 * np.pi * V_rf / max(V_dc + V_rf, 1.0)
        offs = np.linspace(-span, span, int(n_phase))
        weight = np.hanning(int(n_phase))
        if not np.any(weight):
            weight = np.ones(int(n_phase))
        weight = weight / weight.sum()
        vs = V_dc + V_rf * np.sin(phi0[:, None] + offs[None, :])
        E0 = np.maximum((vs * weight[None, :]).sum(axis=1), 0.5)
        # Transverse presheath temperature gives the main anticorrelation, but RF transit smears it.
        sig0 = np.deg2rad(iadf_hwhm_deg) / 1.1774
        smear = ion_spread_floor + (1.0 - ion_spread_floor) * np.sqrt(V_dc / np.maximum(E0, 0.5))
        th = rng.normal(0.0, sig0 * smear, n)
        sB = np.zeros(n)
        return E0, th, sB

    phi_e = _burst_phase(rng, n, V_dc, V_rf, Te)
    Ez_res = rng.exponential(Te, n)
    Et = rng.exponential(Te, n)
    az = rng.uniform(0.0, 2.0 * np.pi, n)
    vx_comp = np.sqrt(Et) * np.cos(az)
    th = np.arctan2(vx_comp, np.sqrt(Ez_res))
    E0 = Ez_res + vx_comp * vx_comp
    frac = (boundary_um / sheath_um) ** (4.0 / 3.0)
    sB = frac * (V_dc + V_rf * np.sin(phi_e))
    return E0, th, sB
