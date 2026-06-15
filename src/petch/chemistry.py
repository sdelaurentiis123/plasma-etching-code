"""Surface chemistry -> normal velocity.

`surface_rate` dispatches on flags.chemistry:
  - "langmuir": the original PoC competitive-Langmuir steady state (byte-identical default).
  - "belen":    the exact ViennaPS coupled-coverage model (contributor #1) — added in Step 3.
"""
import numpy as np


def _ied_yield(A, Eth, par):
    """Yield Y = A*max(sqrt(E)-sqrt(Eth),0) integrated over the ion energy distribution (IED).

    par['ied_mode']:
      'mean'    : evaluate at Emean (PoC default; cheapest, ignores distribution).
      'gauss'   : integrate over a Gaussian N(Emean, Esig) -- matches ViennaPS (initNormalDistEnergy).
      'bimodal' : integrate over the BIMODAL arcsine sheath IED (horns at +-dE/2) -- the REAL
                  low-frequency-bias RF-sheath distribution (Kawamura 1999). BEYOND ViennaPS.
    sqrt(E) is concave so <Y(E)> < Y(<E>) (Jensen): using the mean OVER-estimates the yield; the
    effect grows with distribution width and is sharply nonlinear when the low horn nears Eth.
    Integration is a fixed quadrature (~free: done once per step, not per ray)."""
    mode = par.get('ied_mode', 'mean')
    Em = float(par['Emean']); sEth = np.sqrt(max(Eth, 0.0))
    if mode == 'mean':
        return float(A * max(np.sqrt(max(Em, 0.0)) - sEth, 0.0))
    if mode == 'gauss':
        Es = float(par.get('Esig', 10.0))
        z = np.linspace(-3.5, 3.5, 25); w = np.exp(-0.5 * z * z); w /= w.sum()
        E = np.clip(Em + Es * z, 0.0, None)
        return float((w * A * np.maximum(np.sqrt(E) - sEth, 0.0)).sum())
    # 'bimodal' arcsine: E = Em + (dE/2)cos(phi), phi ~ U(0,pi) -> arcsine density (horns at +-dE/2)
    dE = float(par.get('ied_dE', 40.0))
    phi = np.linspace(0.0, np.pi, 49)
    E = np.clip(Em + 0.5 * dE * np.cos(phi), 0.0, None)
    return float(np.mean(A * np.maximum(np.sqrt(E) - sEth, 0.0)))


def _yields(par):
    """Ion-enhanced, physical-sputter, passivation yields, integrated over the ion energy
    distribution (par['ied_mode']; default 'mean' = evaluate at Emean)."""
    return (_ied_yield(par['A_ie'], par['Eth_ie'], par),
            _ied_yield(par['A_sp'], par['Eth_sp'], par),
            _ied_yield(par['A_p'], par['Eth_p'], par))


def angular_factors(cos_i, par, mode):
    """Per-channel ion angular yield factors (f_ie for ion-enhanced/passivation, f_sp for sputter).

    'cosine' (PoC): both = cos(theta) -> reproduces the original single forward-peaked factor.
    'viennaps' (contributor #3): the exact ViennaPS forms — sputter peaks at oblique angles
    (B_sp), ion-enhanced is flat to 60 deg then falls to 0 at grazing.
    """
    c = np.clip(cos_i, 0.0, 1.0)
    if mode == "viennaps":
        B = par.get('B_sp', 9.3)
        f_sp = np.maximum((1.0 + B * (1.0 - c * c)) * c, 0.0)
        theta = np.arccos(c)
        f_ie = np.where(c >= 0.5, 1.0, np.maximum(3.0 - 6.0 * theta / np.pi, 0.0))
        return f_ie, f_sp
    return c, c                                         # cosine (PoC)


def surface_rate_langmuir(m_i, m_F, m_O, cos_i, is_mask, par, flags=None):
    """Competitive-Langmuir steady state. Per-channel angular yields; cosine mode == PoC."""
    Yie, Ysp, Yp = _yields(par)
    mode = "cosine" if flags is None else getattr(flags, "yield_angular", "cosine")
    f_ie, f_sp = angular_factors(cos_i, par, mode)
    Fi = par['ionFlux'] * m_i                           # geometric ion flux (angular is in yields)
    Fev = par['Fflux'] * m_F * par.get('cal_F', 1.0)    # flux-normalization calibration
    Fp = par['Oflux'] * m_O
    eps = 1e-9
    Yie_a = Yie * f_ie; Ysp_a = Ysp * f_sp; Yp_a = Yp * f_ie
    # competitive Langmuir steady-state coverages
    rF = par['s_F'] * Fev / (Yie_a * Fi + eps)          # theta_F / bare
    rO = par['s_O'] * Fp / (Yp_a * Fi + eps)            # theta_O / bare
    bare = 1.0 / (1.0 + rF + rO)
    thF = rF * bare
    V = (1.0 / par['rho']) * (Yie_a * Fi * thF + Ysp_a * Fi * bare)   # Si removal -> velocity
    V = V * par['rate_scale']
    V[is_mask] = 0.0                                    # mask not etched
    return V


def surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=None):
    """Dispatch surface chemistry model. Default = langmuir (PoC behavior)."""
    model = "langmuir" if flags is None else getattr(flags, "chemistry", "langmuir")
    if model == "langmuir":
        return surface_rate_langmuir(m_i, m_F, m_O, cos_i, is_mask, par, flags)
    elif model == "belen":
        from .belen import surface_rate_belen
        return surface_rate_belen(m_i, m_F, m_O, cos_i, is_mask, par, flags)
    else:
        raise ValueError(f"unknown chemistry model: {model}")
