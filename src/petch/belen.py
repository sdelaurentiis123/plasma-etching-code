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
from .chemistry import _yields


def surface_rate_belen(m_i, m_F, m_O, cos_i, is_mask, par, flags=None):
    Yie, Ysp, Yp = _yields(par)
    fang = np.clip(cos_i, 0.0, 1.0)
    Fi = par['ionFlux'] * m_i * fang                  # ion flux x angular (PoC convention)
    eps = 1e-9

    GY_ie = Yie * Fi                                  # ion-enhanced etchant removal rate
    GY_p = Yp * Fi                                    # ion-enhanced passivation removal rate
    Gb_E = par['Fflux'] * m_F + eps                   # adsorbed F flux (m_F carries betaE sticking)
    Gb_P = par['Oflux'] * m_O + eps                   # adsorbed O flux (m_O carries betaO sticking)

    a = (par['k_sigma'] + 2.0 * GY_ie) / Gb_E
    b = (par['beta_sigma'] + GY_p) / Gb_P
    thF = 1.0 / (1.0 + a * (1.0 + 1.0 / (b + eps)))   # fluorine coverage (coupled to O via b)

    # ViennaPS rate: chemical etch + physical sputter + ion-enhanced etch
    rate = par['k_sigma'] * thF / 4.0 + Ysp * Fi + thF * Yie * Fi
    V = (1.0 / par['rho']) * rate * par['rate_scale']
    V[is_mask] = 0.0
    return V
