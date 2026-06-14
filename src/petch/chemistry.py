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


def surface_rate_langmuir(m_i, m_F, m_O, cos_i, is_mask, par):
    """Original PoC: competitive-Langmuir steady state. Verbatim from feature_etch.py."""
    Yie, Ysp, Yp = _yields(par)
    # angular factor for ion yields (forward-peaked; ~cos incidence)
    fang = np.clip(cos_i, 0.0, 1.0)
    Fi = par['ionFlux'] * m_i * fang
    Fev = par['Fflux'] * m_F * par.get('cal_F', 1.0)   # flux-normalization calibration
    Fp = par['Oflux'] * m_O
    eps = 1e-9
    # competitive Langmuir steady-state coverages
    rF = par['s_F'] * Fev / (Yie * Fi + eps)        # theta_F / bare
    rO = par['s_O'] * Fp / (Yp * Fi + eps)          # theta_O / bare
    bare = 1.0 / (1.0 + rF + rO)
    thF = rF * bare
    V = (1.0 / par['rho']) * (Yie * Fi * thF + Ysp * Fi * bare)   # Si removal -> normal velocity
    V = V * par['rate_scale']
    V[is_mask] = 0.0                                    # mask not etched
    return V


def surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=None):
    """Dispatch surface chemistry model. Default = langmuir (PoC behavior)."""
    model = "langmuir" if flags is None else getattr(flags, "chemistry", "langmuir")
    if model == "langmuir":
        return surface_rate_langmuir(m_i, m_F, m_O, cos_i, is_mask, par)
    elif model == "belen":
        from .belen import surface_rate_belen
        return surface_rate_belen(m_i, m_F, m_O, cos_i, is_mask, par, flags)
    else:
        raise ValueError(f"unknown chemistry model: {model}")
