"""Plasma / surface parameters and model-selection flags.

The default `PAR` reproduces the original feature_etch.py proof-of-concept exactly.
The extra ViennaPS/Belen constants (k_sigma, beta_sigma, B_sp, betaE, betaO) are used
only by the `belen` chemistry model (contributor #1); they are inert under the default
`langmuir` path.
"""
from dataclasses import dataclass, asdict
import numpy as np

# ----------------------------- plasma / surface parameters -----------------------------
# (ViennaPS SF6O2 defaults; matches the original PoC)
PAR = dict(
    ionFlux=12.0, Fflux=1800.0, Oflux=100.0,   # relative fluxes (ViennaPS SF6O2 defaults)
    Emean=100.0, Esig=10.0,                     # ion energy distribution (eV)
    ion_ang_sigma=np.deg2rad(2.5),              # ion angular spread (near-vertical)
    Eth_ie=15.0, Eth_sp=20.0, Eth_p=10.0,       # yield thresholds (eV)
    A_ie=7.0, A_sp=0.0337, A_p=3.0,             # yield prefactors
    s_F=0.20, s_O=0.30,                         # neutral sticking coefficients (PoC values)
    rho=5.02,                                   # substrate density factor
    rate_scale=1.0,                             # global rate calibration (knob; goal = 1.0)
    # --- ViennaPS Belen/Ertl extras (used by chemistry='belen' only) ---
    k_sigma=300.0,                              # chemical-etch coupling (1e15 cm^-2 s^-1)
    beta_sigma=0.04,                            # passivation coupling (1e15 cm^-2 s^-1)
    B_sp=9.3,                                   # angular sputter coefficient
    betaE=0.7,                                  # ViennaPS fluorine sticking (Si)
    betaO=1.0,                                  # ViennaPS oxygen sticking (Si)
    # Flux-normalization calibration: our open-field-normalized etchant flux underrepresents
    # the absolute ViennaPS flux ratio by ~12x. Fitting this to the ViennaPS ground-truth ARDE
    # gives the minimum at cal_F=12 (ARDE rmse 0.110 -> 0.0165). This is the dominant fidelity
    # fix and it un-sticks 3D HARC. cal_F=1.0 recovers the uncalibrated PoC. See FINDINGS.md.
    cal_F=12.0,
    # Ion energy distribution (IED) for yield integration. 'mean' = evaluate yields at Emean (PoC);
    # 'gauss' = integrate over N(Emean,Esig) (matches ViennaPS); 'bimodal' = arcsine sheath IED of
    # full width ied_dE (the REAL low-freq-bias distribution, beyond ViennaPS). See chemistry._ied_yield.
    ied_mode='mean',
    ied_dE=40.0,            # bimodal IED full width (eV); Delta-E ~ Vs/(f*sqrt(M)) (Kawamura 1999)
    # Redeposition (flags.redeposition): product emitted from each face ~ etch rate, sticks with
    # s_redep, and k_redep couples the redeposited flux back into a velocity reduction (passivation).
    s_redep=0.5,
    k_redep=1.0,
    # Surface charging (flags.surface_charging; BEYOND ViennaPS). Electrons arrive diffusely (cosine,
    # unity-sticking) so they are MORE geometrically shadowed in HARC than the directional ions -> the
    # floor floats positive -> the effective floor ion flux is throttled toward the electron arrival
    # rate (insulator current balance): f_charge = 1 - charge_alpha*(1 - Gamma_e/Gamma_i). charge_alpha
    # in [0,1] = insulator-ness / charging strength (0 = off). Calibrated to Hwang-Giapis 1997 (floor
    # ion current drops ~60% by aspect-ratio ~4). eFlux = electron source flux (open-field normalized).
    charge_alpha=0.0,
    eFlux=1.0,
    # Electron angular spread at the surface (rad). Electrons reach the wafer through the sheath with
    # PARTIAL collimation -- narrower than a full cosine (Lambertian) hemisphere but wider than the
    # near-vertical ions. This sets how much MORE shadowed electrons are than ions in HARC (i.e. the
    # steepness of the charging rolloff). Calibrated to Hwang-Giapis 1997. ~0.5 rad (29 deg).
    e_ang_sigma=0.5,
)


@dataclass
class Flags:
    """Model-selection toggles. Defaults reproduce the original PoC byte-for-byte.

    Each non-default value flips exactly one bias contributor (or speedup) so the
    harness can attribute a measured number to each.
    """
    # Defaults = the ViennaPS-calibrated ACCURATE config (belen + viennaps angular + cal_F=12).
    # Set chemistry='langmuir', yield_angular='cosine', cal_F=1.0 to recover the original PoC.
    chemistry: str = "belen"         # "belen" (accurate, contributor #1) | "langmuir" (PoC)
    yield_energy: str = "mean"       # "mean" (PoC) | "ied"   (contributor #2)
    yield_angular: str = "viennaps"  # "viennaps" (accurate, contributor #3) | "cosine" (PoC)
    ion_reflection: bool = False     # contributor #4
    advection: str = "upwind1"       # "upwind1" (PoC) | "weno_rk2" (contributor #5)
    flux_clips: bool = True          # contributor #6a (remove heuristic flux caps)
    reemit_law: str = "3d_sqrtU"     # "3d_sqrtU" (PoC) | "2d_arcsin" (contributor #6b)
    sampling: str = "pseudo"         # "pseudo" (PoC) | "sobol" (speedup: QMC)
    transport_split: bool = False    # speedup: ion few-ray / neutral radiosity
    ion_reflection: bool = False     # contributor #4: ion specular reflection (3D)
    coverage_sticking: bool = False  # Langmuir coverage-dependent neutral sticking (3D ARDE fix)
    redeposition: bool = False       # BEYOND ViennaPS: etch-product redeposition (sidewall passivation/taper)
    neutral_transport: str = "mc"    # "mc" (Russian-roulette MC) | "radiosity" (deterministic, exact deep-floor)
    warm_start_coverage: bool = False  # speed (accuracy-neutral): seed the coverage fixed point from the
                                       # previous step's coverage -> converges in 1-2 iters not 4 (same
                                       # fixed point). The front moves <1 cell/step so coverage barely moves.
    surface_charging: bool = False   # BEYOND ViennaPS: differential electron/ion charging throttles the
                                     # deep-HARC floor ion flux (Hwang-Giapis). See PAR charge_alpha/eFlux.

    def to_dict(self):
        return asdict(self)


DEFAULT_FLAGS = Flags()
