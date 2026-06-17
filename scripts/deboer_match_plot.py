#!/usr/bin/env python3
"""Plot the de Boer match from deboer_final.npz: petch DEFAULT (ViennaPS regime) vs petch with DE-BOER
PROCESS params vs the experiment. Shows that the experiment gap is a PROCESS-PARAMETER calibration
(etchant-starved flux ratio + sub-degree IADF + sputter floor), not a missing transport model."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

npz = sys.argv[1] if len(sys.argv) > 1 else "deboer_final.npz"
d = np.load(npz)
ar_d, nr_d, ar_b, nr_b = d['ar_d'], d['nr_d'], d['ar_b'], d['nr_b']
EXP_AR, EXP_R = d['exp_ar'], d['exp_r']
rmse_b, rmse_d = float(d['rmse_b']), float(d['rmse_d'])


def clean(ar, nr):
    ok = ~np.isnan(nr)
    ar, nr = ar[ok], nr[ok]
    if len(nr) >= 5:                          # light 5-pt smoothing of the MC-noisy gradient curve
        k = np.ones(5) / 5.0
        nr = np.convolve(np.clip(nr, 0, 1.0), k, mode='same')
        nr[:2] = np.clip(nr[:2], 0, 1.0); nr[-2:] = nr[-2:]
    return ar, np.clip(nr, 0, 1.02)


fig, ax = plt.subplots(figsize=(8.6, 5.6))
xr = np.linspace(0, 40, 200)
ax.plot(xr, np.interp(xr, EXP_AR, EXP_R), 'k-', lw=2.6, label='de Boer/Blauw experiment (cryo SF6 DRIE)')
ax.plot(EXP_AR, EXP_R, 'ko', ms=8)
ax.plot(*clean(ar_d, nr_d), '--', color='#1a73e8', lw=2,
        label=f'petch, ViennaPS-default params  (RMSE {rmse_d:.2f})')
ax.plot(*clean(ar_b, nr_b), '-', color='#c0392b', lw=2.2,
        label=f'petch, de-Boer process params  (RMSE {rmse_b:.2f})')

# shade the two regimes: knee region (params help) vs floor region (frontier)
ax.axvspan(8, 20, color='#2e7d32', alpha=0.07)
ax.axvspan(20, 40, color='#9e9e9e', alpha=0.10)

ax.set_xlim(0, 40); ax.set_ylim(0, 1.04)
ax.set_xlabel("aspect ratio  (depth / width)", fontsize=12)
ax.set_ylabel("normalized bottom etch rate", fontsize=12)
ax.set_title("Closing the de Boer gap: de-Boer process params fix the KNEE;\n"
             "the high-AR FLOOR (ion delivery past AR~20) is the remaining frontier", fontsize=11.5)
ax.legend(loc='upper right', fontsize=10)
ax.grid(alpha=0.3)
ax.annotate("ViennaPS regime: etchant-rich,\nbroad IADF -> too gentle",
            xy=(17, np.interp(17, *clean(ar_d, nr_d))), xytext=(20.0, 0.80),
            fontsize=9.2, color='#1a4e8a', arrowprops=dict(arrowstyle='->', color='#1a73e8'))
ax.annotate("knee moved by etchant-starved\nflux ratio + narrower IADF",
            xy=(13, np.interp(13, *clean(ar_b, nr_b))), xytext=(2.4, 0.18),
            fontsize=9.2, color='#7a241a', arrowprops=dict(arrowstyle='->', color='#c0392b'))
ax.annotate("floor frontier: petch (& ViennaPS)\ncollapse past AR~20; experiment holds ~0.25",
            xy=(30, 0.245), xytext=(22.5, 0.42),
            fontsize=9.2, color='#555', arrowprops=dict(arrowstyle='->', color='#888'))

fig.tight_layout()
out = "/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/deboer_match.png"
fig.savefig(out, dpi=150)
print("wrote", out)
