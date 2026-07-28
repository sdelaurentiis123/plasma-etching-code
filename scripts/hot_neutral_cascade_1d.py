#!/usr/bin/env python3
"""1-D hot-neutral wall-cascade discriminator (no fitting).

Trench of aspect ratio AR: ions launched from the mouth with a narrow
angular spread; every wall collision applies the MCFPM leftover-probability
selection under a candidate interpretation, and Eq. 2.34 energy retention
(E_ts=100 eV, E_c=10 eV, theta_c=70 deg; Huang thesis). Specular bounces
preserve the polar angle in a straight trench, so the cascade is a survival
product over bounce count. Compares predicted etch-front hot-neutral flux
versus Huang's published values (3.1e15 at AR~0, 8.0e15 peak at AR 4,
1.1e15 at AR 40; ions 2.0e15 -> 0.3e15).
"""
import numpy as np

E_TS, E_C, THETA_C = 100.0, 10.0, 70.0
E_ION = 1500.0
ION_FLUX_OPEN = 2.0e15
HOT_AT_0, HOT_PEAK_AR4, HOT_AT_40 = 3.1e15, 8.0e15, 1.1e15


def kress(cos_t, B=9.3):
    return np.maximum((1.0 + B * (1.0 - cos_t ** 2)) * cos_t, 0.0)


def f_energy(E, p0=0.9, eth=20.0, e0=500.0, q=0.5):
    return p0 * np.maximum(E ** q - eth ** q, 0.0) / (e0 ** q - eth ** q)


def selection_probability(E, cos_t, mode):
    if mode == "B":          # f(E) and angular both in selection
        return np.clip(f_energy(E) * kress(cos_t), 0.0, 1.0)
    if mode == "Bp":         # angular only; f(E) scales removal count
        return np.clip(0.9 * kress(cos_t), 0.0, 1.0)
    raise ValueError(mode)


def retained_energy(E, theta_deg):
    if E > E_TS and theta_deg > THETA_C:
        return E                          # pure specular
    if E < E_C or theta_deg < THETA_C:
        return 0.0                        # diffusive: treat as lost to cascade
    return E * (theta_deg - THETA_C) / (90.0 - THETA_C) \
        * (E - E_C) / (E_TS - E_C)        # Eq. 2.34 interpolation


def front_flux(AR, mode, spread_deg=2.0, n_rays=4000):
    """Hot-neutral flux reaching the floor of a trench of aspect ratio AR."""
    rng = np.random.default_rng(7)
    polar = np.abs(rng.normal(0.0, spread_deg, n_rays))    # deg off vertical
    surviving = 0.0
    direct_ions = 0.0
    for th in polar:
        # Ray geometry: hits floor directly if tan(th) < 1/AR (width=1).
        if np.tan(np.radians(th)) < 1.0 / max(AR, 1e-9):
            direct_ions += 1.0
            continue
        # Wall bounce count down the trench: each traverse covers
        # width/tan(th) of depth; need total depth AR.
        depth_per_bounce = 1.0 / np.tan(np.radians(th))
        bounces = int(np.ceil(AR / depth_per_bounce))
        E = E_ION
        w = 1.0
        theta_wall = 90.0 - th                 # incidence on vertical wall
        cos_wall = np.cos(np.radians(theta_wall))
        alive = True
        for _ in range(bounces):
            p_react = selection_probability(E, cos_wall, mode)
            w *= (1.0 - p_react)
            E = retained_energy(E, theta_wall)
            if w < 1e-6 or E <= E_C:
                alive = False
                break
        if alive:
            surviving += w
    return (ION_FLUX_OPEN * surviving / n_rays,
            ION_FLUX_OPEN * direct_ions / n_rays)


if __name__ == "__main__":
    print(f"{'AR':>4} {'mode':>4} {'hot@front':>12} {'ions@front':>12}")
    for mode in ("B", "Bp"):
        for AR in (1, 4, 8, 20, 40):
            hot, ions = front_flux(AR, mode)
            print(f"{AR:>4} {mode:>4} {hot:12.3e} {ions:12.3e}")
    print("\nHuang published: hot 8.0e15@AR4 (peak), 1.1e15@AR40; "
          "ions 2.0e15->0.3e15")
