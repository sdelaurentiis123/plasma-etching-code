"""Contributor #1: the exact ViennaPS Belen/Ertl coupled-coverage SF6/O2 model.

Ported from ViennaPS `psPlasmaEtching.hpp`. The differences from the PoC competitive-Langmuir
steady state (`chemistry.surface_rate_langmuir`) are:

  1. TWO coupled coverages with cross terms (not a single bare ratio):
         a = (k_sigma + 2*GY_ie) / Gb_E
         b = (beta_sigma + GY_p) / Gb_P
         theta_F = 1 / (1 + a*(1 + 1/b))
  2. A STANDALONE chemical-etch term `k_sigma*theta_F/4` in the rate that survives even where
     ions are blocked (absent from the PoC — prime suspect for the absolute-rate gap).
  3. ViennaPS sticking betaE=0.7 / betaO=1.0 (the driver runs transport with these for belen).

The yields here are still mean-energy + cosine-angular (the PoC defaults) so that this toggle
isolates the COVERAGE/RATE form. Contributors #2 (IED) and #3 (ViennaPS angular yield) are
separate flags layered on top later.

`par['rate_scale']` plays the role of the ViennaPS `unitConversion` (flux/density -> um/min).
Whether it can be set to a single physical constant (vs the PoC's empirical 0.29) is the
measured question — see scripts/run_phase0.py.
"""
import numpy as np
from .chemistry import _yields, angular_factors


def surface_rate_belen(m_i, m_F, m_O, cos_i, is_mask, par, flags=None):
    Yie, Ysp, Yp = _yields(par)
    mode = "cosine" if flags is None else getattr(flags, "yield_angular", "cosine")
    f_ie, f_sp = angular_factors(cos_i, par, mode)
    Fi = par['ionFlux'] * m_i                          # geometric ion flux (angular in yields)
    Yie_a = Yie * f_ie; Ysp_a = Ysp * f_sp; Yp_a = Yp * f_ie
    eps = 1e-9

    GY_ie = Yie_a * Fi                                # ion-enhanced etchant removal rate
    GY_p = Yp_a * Fi                                  # ion-enhanced passivation removal rate
    # Flux normalization to the ViennaPS convention: our m_F records STUCK flux (open-field ~=
    # betaE), ViennaPS normalizes ARRIVING flux to 1 on open field. Divide by fnorm (=betaE/betaO)
    # to match. fnorm=1.0 recovers the uncorrected PoC-style normalization (for A/B testing).
    fnE = par.get('fnorm_E', 1.0)
    fnO = par.get('fnorm_O', 1.0)
    calF = par.get('cal_F', 1.0)                       # flux-normalization calibration (~12)
    Gb_E = par['Fflux'] * m_F * calF / fnE + eps       # arriving F flux (ViennaPS convention)
    Gb_P = par['Oflux'] * m_O / fnO + eps              # arriving O flux

    a = (par['k_sigma'] + 2.0 * GY_ie) / Gb_E
    b = (par['beta_sigma'] + GY_p) / Gb_P
    thF = 1.0 / (1.0 + a * (1.0 + 1.0 / (b + eps)))   # fluorine coverage (coupled to O via b)

    # ViennaPS rate: chemical etch + physical sputter + ion-enhanced etch.
    # Ysp_scale (default 1.0 = exact ViennaPS) lifts the AR-INDEPENDENT physical-sputter floor: in
    # deep features the chemical/ion-enhanced terms vanish with the neutral coverage, so this directional
    # ion-sputter term sets the high-AR etch-rate floor. The de Boer cryo-DRIE floor (~0.20 of the mouth
    # rate at AR 40) is ~3x the ViennaPS-default sputter floor -> calibrate it here, not via sticking.
    sps = par.get('Ysp_scale', 1.0)
    rate = par['k_sigma'] * thF / 4.0 + sps * Ysp_a * Fi + thF * Yie_a * Fi
    V = (1.0 / par['rho']) * rate * par['rate_scale']
    V[is_mask] = 0.0
    return V
