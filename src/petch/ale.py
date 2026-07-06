"""Directional plasma Atomic Layer Etching (Si / Cl2 / Ar+) — cyclic three-layer site-balance ROM.

Reduced-order model of Vella & Graves (2025), "Si-Cl2-Ar+ Atomic Layer Etching Window", OSTI 2586627
(https://www.osti.gov/biblio/2586627). A cyclic wrapper around the same coupled-coverage kinetics
petch already uses for continuous etch: each cycle is (1) a Cl2 modification dose that chlorinates the
top layer (analytic exponential) followed by (2) an Ar+ bombardment step that removes the chlorinated
Si (two mildly-stiff coupled ODEs). Self-limitation emerges from the yields: below ~20 eV the chemical
removal completes and the physical-sputter leak Y_Si is negligible (flat EPC = ALE window); above
~22.5 eV the sputter leak dominates and EPC rises (loss of self-limitation).

Everything here is smooth in the ion energy E and the parameters, so the whole cycle composes with
reverse-mode autodiff (the differentiable-calibration path). Pure numpy; no box needed.

Benchmark gates (see CHEMISTRY_FRONTIER_ALE.md): EPC ~0.67/0.75/1.0 Å/cycle at 15/17.5/20 eV
(the 15-20 eV window), rise to ~1.7/6 at 22.5/30 eV; Cl uptake plateau 0.8e15 cm^-2; synergy S->~100%.
"""
import numpy as np

# fixed constants (Vella-Graves ROM)
SIGMA1 = 1.0e15      # top-layer site density, sites/cm^2
GAMMA_CL = 0.25      # Cl sticking during modification
K_CL_MIX = 0.45      # top<->mixed Cl exchange coefficient
J_CL2 = 9.8e17       # Cl2 flux during modification, cm^-2 s^-1
J_AR = 3.7e16        # Ar+ flux during bombardment, cm^-2 s^-1
N_SI = 5.0e22        # Si atomic density, cm^-3  (EPC[A] = N_removed[cm^-2]/N_SI*1e8)
E_TH_MIX = 2.81      # mixed-layer-depth threshold energy, eV


def yields(E):
    """Energy-dependent sputter yields and mixed-layer depth ratio sigma2/sigma1. All clamped >=0,
    all smooth in E (linear / sqrt) so gradients are exact. Returns dict."""
    E = np.asarray(E, dtype=float)
    return dict(
        Y_Cl=np.maximum(0.0, 3.2e-3 * E - 4.25e-2),
        Y_SiCl=np.maximum(0.0, 1.0e-4 * E - 1.4e-3),
        Y_SiCl2=np.maximum(0.0, 2.0e-3 * E - 2.9e-2),
        Y_Si=np.maximum(0.0, 4.5e-5 * E - 9.0e-4),
        s2_ratio=np.maximum(1e-6, 0.77 * (np.sqrt(np.maximum(E, 0.0)) - np.sqrt(E_TH_MIX))),
    )


def modification_step(theta1_init, t_dose):
    """Cl2 dose: analytic chlorination of the top layer. Returns (theta1, cl_uptake[cm^-2])."""
    alpha = 2.0 * J_CL2 * GAMMA_CL / SIGMA1
    theta1 = 1.0 - (1.0 - theta1_init) * np.exp(-alpha * t_dose)
    cl_uptake = (theta1 - theta1_init) * SIGMA1   # Cl atoms taken up this step
    return theta1, cl_uptake


def bombardment_step(theta1, theta2, E, t_barr, dt=0.02):
    """Ar+ bombardment: integrate the two coupled coverage ODEs (Vella-Graves Eqs. 6-7), accumulating
    Si removed via ALL THREE product channels (Eqs. 9-11):
      dtheta1/dt = -(J_Ar/s1)(Y_Cl t1 + Y_SiCl t1 + 2 Y_SiCl2 t1^2 + K t1)      [Eq 6]
      dtheta2/dt = -(J_Ar/s2)(Y_Cl t2 + Y_SiCl t2 + 2 Y_SiCl2 t2^2 - K t1)      [Eq 7]
      J_SiCl  = J_Ar Y_SiCl  (t1+t2)   [Eq 9,  1 Si/molecule]
      J_SiCl2 = J_Ar Y_SiCl2 (t1+t2)   [Eq 10, 1 Si/molecule]
      J_Si    = J_Ar Y_Si    (1 - t2)  [Eq 11, bare-Si sputter, from the mixed layer]
    The chemical SiCl/SiCl2 channels drive the 15-20 eV window (Y_Si=0 there); the bare-Si sputter
    leak turns on ~20 eV and breaks self-limitation. Forward Euler, dt=0.02 s resolves the ~0.06 s
    K-mixing transient. Returns (theta1, theta2, N_Si_removed[cm^-2])."""
    y = yields(E)
    s1 = SIGMA1
    s2 = y["s2_ratio"] * SIGMA1
    n = max(1, int(round(t_barr / dt)))
    dt = t_barr / n
    N_si = 0.0
    for _ in range(n):
        r1 = (J_AR / s1) * (y["Y_Cl"] * theta1 + y["Y_SiCl"] * theta1
                            + 2.0 * y["Y_SiCl2"] * theta1 ** 2 + K_CL_MIX * theta1)
        r2 = (J_AR / s2) * (y["Y_Cl"] * theta2 + y["Y_SiCl"] * theta2
                            + 2.0 * y["Y_SiCl2"] * theta2 ** 2 - K_CL_MIX * theta1)
        theta_sum = theta1 + theta2
        j_si = J_AR * (y["Y_SiCl"] * theta_sum + y["Y_SiCl2"] * theta_sum + y["Y_Si"] * (1.0 - theta2))
        theta1 = max(0.0, theta1 - dt * r1)
        theta2 = min(1.0, max(0.0, theta2 - dt * r2))
        N_si += dt * j_si
    return theta1, theta2, N_si


def run_ale(E, n_cycles=8, t_dose=0.112, t_barr=113.5, dt=0.02, return_history=False):
    """Run n_cycles of (modification, bombardment) to cyclic steady state at ion energy E (eV).
    t_dose, t_barr are the step durations (s). Returns dict with EPC (A/cycle, steady tail-mean over
    the last 3 cycles), Cl uptake per mod step, and the modification/removal-only decomposition for
    the ALE synergy metric S = (EPC - alpha - beta)/EPC (Kanarik)."""
    theta1, theta2 = 0.0, 0.0
    epc_hist, cl_hist = [], []
    for c in range(n_cycles):
        theta1, cl_uptake = modification_step(theta1, t_dose)
        theta1, theta2, N_si = bombardment_step(theta1, theta2, E, t_barr, dt=dt)
        epc_hist.append(N_si / N_SI * 1e8)   # Angstrom/cycle
        cl_hist.append(cl_uptake)
    tail = slice(max(0, n_cycles - 3), n_cycles)
    epc = float(np.mean(epc_hist[tail]))
    out = dict(E=float(E), epc=epc, cl_uptake=float(np.mean(cl_hist[tail])),
               epc_history=epc_hist, cl_history=cl_hist)
    if return_history:
        out["theta_final"] = (theta1, theta2)
    return out


def synergy(E, t_dose=0.112, t_barr=113.5, dt=0.02):
    """ALE synergy S = (EPC - alpha - beta)/EPC (Kanarik JVST A 33,020802). alpha = etch from
    modification alone (no Ar+), beta = etch from Ar+ alone (no Cl2 dose). S->~100% = ideal ALE."""
    full = run_ale(E, t_dose=t_dose, t_barr=t_barr, dt=dt)["epc"]
    # beta: bombardment with NO chlorination (theta1=theta2=0 each cycle) -> pure sputter leak
    _, _, N_beta = bombardment_step(0.0, 0.0, E, t_barr, dt=dt)
    beta = N_beta / N_SI * 1e8
    # alpha: modification alone removes nothing (no ion-driven removal) -> ~0
    alpha = 0.0
    S = (full - alpha - beta) / full if full > 0 else 0.0
    return dict(E=float(E), epc=full, alpha=alpha, beta=beta, synergy=float(S))
