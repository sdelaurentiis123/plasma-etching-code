"""Surface chemistry -> normal velocity.

`surface_rate` dispatches on flags.chemistry:
  - "langmuir": the original PoC competitive-Langmuir steady state (byte-identical default).
  - "belen":    the exact ViennaPS coupled-coverage model (contributor #1) — added in Step 3.
"""
import numpy as np


def _yields(par):
    """Energy-dependent yields evaluated at the mean ion energy (PoC; contributor #2 = 'mean')."""
    E = par['Emean']
    sqrtE = np.sqrt(max(E, 0.0))
    Yie = par['A_ie'] * max(sqrtE - np.sqrt(par['Eth_ie']), 0.0)
    Ysp = par['A_sp'] * max(sqrtE - np.sqrt(par['Eth_sp']), 0.0)
    Yp = par['A_p'] * max(sqrtE - np.sqrt(par['Eth_p']), 0.0)
    return Yie, Ysp, Yp


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
