"""Feature charging: the Hwang-Giapis current-balance closure (insulating floor).

Replaces the retired flux-ratio throttle f = 1 - charge_alpha*(1 - Ge/Gi), which was structurally
wrong (no potential state / energy resolution; broke charge conservation; deleted throttled ions;
collimated-Gaussian electrons). Published physics (Hwang & Giapis, JVST B 15, 70 (1997); JAP 82,
566 (1997); Arnold & Sawin, JAP 70, 5314 (1991)):

  - Every insulating surface element charges until ion and electron currents balance.
  - The floor potential V_f DECELERATES ions: only the IEDF slice with E > e*V_f lands (with energy
    E - e*V_f); the sub-threshold slice is deflected to the sidewall foot (NOT destroyed) -- that
    deflected, accelerated flux is the notching driver.
  - Electrons are Maxwellian (T_e), isotropic-flux at the sheath edge (EADF ~ cos^0.6 at the wafer
    for T_e = 4 eV) -- NOT collimated. Their floor flux is geometric shading SOFTENED by
    electrostatic attraction to the positive floor (HG: the reduction is "significantly smaller
    than ... the decrease in the solid angle").
  - V_f saturates near the ion-energy scale (HG: 8 -> 33 V over AR 1 -> 4 with <E_ion> ~ 37 eV;
    Matsui APL 78, 883: 300 eV ions only fully cut off above AR ~ 7).

0-D closure per floor: solve  Q(e*V_f) = F_e(AR, V_f)  for V_f, where Q is the IEDF survival
function and F_e the electron floor-arrival fraction (Monte-Carlo over Maxwellian speeds,
cos-power angles, ballistic flight in the linear in-trench potential). At the sheath boundary the
time-averaged ion and electron fluxes are equal, so the normalized floor flux is Q(e*V_f) itself.

Gate: the 8-point Gamma_i,floor(AR) curve of HG JAP 82, 566 Fig. 4 (0.59 -> 0.22 over AR 1 -> 4),
V_f(AR) in 8-33 V, and the Matsui high-energy asymptote. See scripts/charging_gate.py.
"""
from __future__ import annotations

import numpy as np


def ied_sample(n, V_dc=37.0, V_rf=30.0, rng=None):
    """Bimodal RF-sheath IEDF: E = V_dc + V_rf*cos(phi), phi ~ U(0, pi) (arcsine density,
    horns at V_dc +- V_rf). HG conditions: V_s = 37 + 30 sin(wt) -> E in ~7-67 eV."""
    rng = rng or np.random.default_rng(0)
    phi = rng.uniform(0.0, np.pi, n)
    return V_dc + V_rf * np.cos(phi)


def ied_survival(V_f, V_dc=37.0, V_rf=30.0):
    """Q(e*V_f) = P(E_ion > e*V_f) for the bimodal arcsine IEDF (analytic)."""
    x = (np.asarray(V_f, float) - V_dc) / V_rf
    x = np.clip(x, -1.0, 1.0)
    return np.arccos(x) / np.pi


def electron_floor_fraction(AR, V_f, Te=4.0, cos_power=0.6, n=200000, seed=0):
    """Fraction of mouth-entering electrons that reach the floor of a trench of aspect ratio AR,
    with an attracting linear potential 0 -> +V_f from mouth to floor (ballistic, 2D trench).

    Electrons: speeds from the Maxwellian flux distribution (energy ~ Gamma(2, Te) for flux-weighted
    sampling in 3D), angles from the wafer-level EADF ~ cos^p(theta) (HG fit p = 0.6 at Te = 4 eV).
    Trajectory in the trench: lateral velocity constant, vertical velocity grows by the potential
    drop; the electron reaches the floor iff its lateral excursion stays inside the width. Entry x
    uniform across the mouth. Normalized so V_f = 0, AR -> 0 gives 1."""
    rng = np.random.default_rng(seed)
    # flux-weighted Maxwellian energies: pdf ~ E * exp(-E/Te) (Gamma k=2)
    E = rng.gamma(2.0, Te, n)                          # eV
    # EADF ~ cos^p(theta) FLUX distribution: sample theta with pdf ~ cos^p(th)*sin(th)
    u = rng.uniform(0.0, 1.0, n)
    ct = (1.0 - u) ** (1.0 / (cos_power + 1.0))        # cos(theta) via inverse CDF of cos^p*sin
    st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
    # azimuth folds into the 2D trench plane: lateral speed component = v*st*|cos(az)|
    az = rng.uniform(0.0, 2.0 * np.pi, n)
    v = np.sqrt(E)                                     # units where m/2 = 1: v = sqrt(E[eV])
    vz0 = v * ct
    vx = v * st * np.abs(np.cos(az))
    x0 = rng.uniform(0.0, 1.0, n)                      # entry position across the mouth (W = 1)
    D = AR                                             # depth in units of W
    # linear accelerating potential: vz(z)^2 = vz0^2 + 2*(V_f)*(z/D) in eV units (vz^2 == E_z)
    # time-of-flight to depth z: integrate dz/vz -> analytic for linear field
    a = np.maximum(V_f, 1e-12) / D                     # dE_z/dz
    # vz(z) = sqrt(vz0^2 + 2*a*z); t(D) = (vz(D)-vz0)/a  (from vz dvz = a dz, dz = vz dt)
    vzD = np.sqrt(vz0 ** 2 + 2.0 * a * D)
    t = np.where(V_f > 1e-9, (vzD - vz0) / a, D / np.maximum(vz0, 1e-12))
    xlat = vx * t                                      # lateral excursion during descent
    # reflect off... no: sidewalls ABSORB electrons (they charge negative; HG sidewall e-flux grows
    # with AR). Electron reaches floor iff it stays inside the trench: x0 + xlat within [0,1]
    # lateral direction random sign:
    sgn = np.where(rng.uniform(0, 1, n) < 0.5, 1.0, -1.0)
    xf = x0 + sgn * xlat
    reach = (xf >= 0.0) & (xf <= 1.0)
    return float(np.mean(reach))


def floor_balance(AR, Te=4.0, V_dc=37.0, V_rf=30.0, cos_power=0.6, n=200000, seed=0,
                  v_lo=0.0, v_hi=None, tol=0.25):
    """Solve Q(e*V_f) = F_e(AR, V_f) for the floor potential by bisection.
    Returns (V_f, floor_flux = Q(e*V_f))."""
    if v_hi is None:
        v_hi = V_dc + V_rf
    lo, hi = v_lo, v_hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        q = ied_survival(mid, V_dc, V_rf)
        fe = electron_floor_fraction(AR, mid, Te, cos_power, n, seed)
        if q > fe:            # more ions than electrons -> floor charges more positive
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    vf = 0.5 * (lo + hi)
    return vf, float(ied_survival(vf, V_dc, V_rf))
