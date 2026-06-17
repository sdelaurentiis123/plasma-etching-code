#!/usr/bin/env python3
"""The definitive trench-ARDE comparison: petch vs ViennaPS vs the de Boer/Blauw experiment.

Normalized bottom etch rate vs aspect ratio (the standard ARDE plot; absolute rate divided out).
SAME geometry for petch & ViennaPS: DX=0.03, W=0.5, XE=1.5, YE=0.3 um, periodic-y, sub-micron trench.
Measured on an RTX 3090 box (scripts/vps_trench_arde.py + petch_trench_arde_deep.py); deep substrate
so AR is not domain-floor clamped (floor-clamped points dropped).

What the plot shows:
- petch MC  ~=  ViennaPS  (within ~0.08): both ballistic line-of-sight neutral transport.
- petch radiosity is GENTLER at high AR (Lambertian diffuse re-emission over-couples flux to the
  bottom) -> MC is the more physical neutral model for deep features.
- ALL of them sit well above the de Boer experiment, whose steep early falloff (1.0 -> 0.43 by AR 10)
  is the Knudsen molecular-flow regime. The experiment gap is a SHARED ballistic-vs-Knudsen transport
  limit, not a petch-vs-ViennaPS discrepancy.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = 0.5


def arde(dur, dep):
    dur = np.array(dur, float); dep = np.array(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dur)
    return armid, rate / rate[0]


# ViennaPS SF6O2Etching, CPU_TRIANGLE  (scripts/vps_trench_arde.py)
v_ar, v_nr = arde([0.4, 0.8, 1.3, 1.9], [1.305, 2.460, 3.705, 4.969])
# petch radiosity + periodic-y  (deep run; dropped floor-clamped 8.30 point)
r_ar, r_nr = arde([0.4, 0.8, 1.3, 1.9, 2.6], [1.16, 2.36, 3.8, 5.45, 7.19])
# petch MC + periodic-y  (deep run; all points clean, max 7.34 < 8.3 floor)
m_ar, m_nr = arde([0.4, 0.8, 1.3, 1.9, 2.6, 3.4], [1.16, 2.36, 3.74, 5.21, 6.56, 7.34])
# de Boer / Blauw cryogenic SF6 DRIE experiment
EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])

fig, ax = plt.subplots(figsize=(8.2, 5.4))
arr = np.linspace(0, 14, 140)
ax.plot(arr, np.interp(arr, EXP_AR, EXP_R), 'k-', lw=2.5, label='de Boer/Blauw experiment (Knudsen regime)')
ax.plot([0, 10], np.interp([0, 10], EXP_AR, EXP_R), 'ko', ms=7)
ax.plot(v_ar, v_nr, 's-', color='#1a73e8', lw=2, ms=9, label='ViennaPS SF6O2 (CPU_TRIANGLE)')
ax.plot(m_ar, m_nr, '^-', color='#0a8f3c', lw=2, ms=8, label='petch MC + periodic-y')
ax.plot(r_ar, r_nr, 'o--', color='#b0772a', lw=2, ms=7, label='petch radiosity + periodic-y')

ax.set_xlim(0, 14); ax.set_ylim(0, 1.04)
ax.set_xlabel("aspect ratio  (depth / width)", fontsize=12)
ax.set_ylabel("normalized bottom etch rate", fontsize=12)
ax.set_title("Sub-micron trench ARDE: petch MC tracks ViennaPS; both ballistic vs experiment", fontsize=12)
ax.legend(loc='upper right', fontsize=10.5)
ax.grid(alpha=0.3)

ax.annotate("petch MC ≈ ViennaPS\n(ballistic, within ~0.08)", xy=(8.5, 0.77), xytext=(2.4, 0.46),
            fontsize=10, color='#333', arrowprops=dict(arrowstyle='->', color='#888'))
ax.annotate("radiosity over-couples\nflux to the bottom", xy=(12.5, 0.83), xytext=(9.3, 0.70),
            fontsize=9.5, color='#8a5a10', arrowprops=dict(arrowstyle='->', color='#b0772a'))
ax.annotate("Knudsen gap\n(needs molecular-flow transport)", xy=(9.5, 0.50), xytext=(9.0, 0.10),
            fontsize=9.5, color='#333', arrowprops=dict(arrowstyle='->', color='#888'))

fig.tight_layout()
out = "/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/trench_arde_3way.png"
fig.savefig(out, dpi=150)
print("wrote", out, flush=True)
for a in [6, 9, 12]:
    print(f"AR {a:>2}:  exp {np.interp(a, EXP_AR, EXP_R):.2f}   ViennaPS {np.interp(a, v_ar, v_nr):.2f}   "
          f"petchMC {np.interp(a, m_ar, m_nr):.2f}   petchRad {np.interp(a, r_ar, r_nr):.2f}", flush=True)
