"""Cryogenic etch chemistry — temperature-dependent physisorption enhancement.

Below ~0 degC a condensed/physisorbed etchant layer (HF/H2O in CF4/H2 "pseudo-wet" plasmas; C4F8 or HF
in cryo-ALE) builds up on the surface and multiplies the etch rate, because the physisorbed species'
surface residence time grows as exp(E_ads/kT) as T drops. This is the hot 2023-2026 chemistry (3D-NAND
/ DRAM) and no open feature-scale tool has it. Reduced model: a Langmuir physisorption coverage theta(T)
times an ion-activated etch channel, layered on the base (ALE/sputter) rate -- smooth in T and every
parameter, so it drops straight into the differentiable pipeline (see ale_diff.py pattern).

Physics + numbers (honest, benchmark-gated -- see refs):
  - Isotherm:   theta(T) = K(T)*p / (1 + K(T)*p),  K(T)*p = A * exp(E_ads / (kB*T))   [rises as T falls]
  - Rate:       ER(T) = R_base * (1 + gain * theta(T))
  - E_ads ~ 0.4-0.5 eV is a FIRM independently-measured physisorption energy (Antoun et al. Sci.Rep.
    11,357 2021: E_d=0.406 eV, t_d0=1e-11 s; HF cryo-ALE ~0.5 eV) -- fixed as the physical prior, not
    free-floated.
Benchmark anchor (CF4/H2 pseudo-wet, Small Methods 2024, doi 10.1002/smtd.202400090, Fig 1):
  SiO2 ER = 2.3 nm/s for T>0 degC (plateau); 3.76 nm/s at -60 degC = 1.6x the +20 degC value.
  (NOTE: the "rate doubles" folk-anchor is really 1.6x, and from CF4/H2 -- not CHF3. Gated on 1.6x.)
Independent cross-check: HF cryo-ALE SiO2 EPC 0.25 -> 0.79 nm/cycle over +20 -> -60 degC (~3.2x, a
different system/onset). C4F8 cryo-ALE desorption cliff: etches at -120 degC, not -110 (residence-time).
"""
import numpy as np

KB = 8.617333e-5   # eV/K
E_ADS = 0.40       # physisorption adsorption enthalpy, eV (firm: Antoun 0.406 eV / HF ~0.5 eV)
A_KP = 1.34e-9     # K0*p prefactor, calibrated so theta(+20C)~0.01 (warm plateau)
R_BASE = 2.3       # base SiO2 etch rate warm-side plateau, nm/s (Small Methods 2024 Fig 1)
GAIN = 0.80        # physisorbed-channel gain, calibrated to the -60C anchor (3.76 nm/s = 1.6x)


def physisorption_coverage(T_C, E_ads=E_ADS, A=A_KP):
    """Langmuir physisorption coverage theta(T), T in degC. theta->0 warm, ->1 cold; smooth+differentiable.
    K(T)*p = A*exp(E_ads/(kB*T)) grows as T drops (longer surface residence)."""
    T_K = np.asarray(T_C, dtype=float) + 273.15
    Kp = A * np.exp(E_ads / (KB * T_K))
    return Kp / (1.0 + Kp)


def cryo_etch_rate(T_C, R_base=R_BASE, gain=GAIN, E_ads=E_ADS, A=A_KP):
    """Cryo etch rate (nm/s): base rate times the physisorbed condensed-etchant enhancement.
    ER(T) = R_base * (1 + gain*theta_phys(T)). Rises as T drops, saturates cold."""
    return R_base * (1.0 + gain * physisorption_coverage(T_C, E_ads=E_ads, A=A))


def residence_time(T_C, E_d=0.406, t_d0=1e-11):
    """Physisorbed-species surface residence time t_d = t_d0*exp(E_d/kB T) (Antoun et al. 2021). The
    etch-on/off cliff is where t_d crosses the process step time: for C4F8 (E_d=0.406 eV) that is the
    -120 degC (etches) vs -110 degC (does not) threshold."""
    T_K = np.asarray(T_C, dtype=float) + 273.15
    return t_d0 * np.exp(E_d / (KB * T_K))


def enhancement_factor(T_C, **kw):
    """Etch-rate enhancement relative to the warm plateau: ER(T)/R_base = 1 + gain*theta(T)."""
    return cryo_etch_rate(T_C, **kw) / kw.get("R_base", R_BASE)
