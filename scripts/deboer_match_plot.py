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
        label='petch, ViennaPS-default params (gentle, never knees)')
ax.plot(*clean(ar_b, nr_b), '-', color='#c0392b', lw=2.2,
        label='petch, knee-tightened (matches to AR~10, then stalls)')

arb_c, nrb_c = clean(ar_b, nr_b)
ar_stall = arb_c.max()
# shade: knee region (params nail it) vs the collapse/floor region (one frontier)
ax.axvspan(0, 10.5, color='#2e7d32', alpha=0.07)
ax.axvspan(10.5, 40, color='#9e9e9e', alpha=0.10)

ax.set_xlim(0, 40); ax.set_ylim(0, 1.04)
ax.set_xlabel("aspect ratio  (depth / width)", fontsize=12)
ax.set_ylabel("normalized bottom etch rate", fontsize=12)
ax.set_title("de Boer: petch now matches the KNEE to AR~10 (cal_F=1.5) — but starving for the\n"
             "knee makes the deep rate COLLAPSE and the trench STALL by AR~%d. Knee & floor are one problem."
             % int(round(ar_stall)), fontsize=10.5)
ax.legend(loc='upper right', fontsize=10)
ax.grid(alpha=0.3)
ax.axvline(ar_stall, color='#c0392b', ls=':', lw=1.3, alpha=0.7)
ax.annotate("ViennaPS regime: etchant-rich\n-> too gentle, never knees",
            xy=(17, np.interp(17, *clean(ar_d, nr_d))), xytext=(20.0, 0.82),
            fontsize=9.2, color='#1a4e8a', arrowprops=dict(arrowstyle='->', color='#1a73e8'))
ax.annotate("knee MATCHES experiment\nthrough AR~10",
            xy=(9, np.interp(9, arb_c, nrb_c)), xytext=(1.6, 0.30),
            fontsize=9.2, color='#7a241a', arrowprops=dict(arrowstyle='->', color='#c0392b'))
ax.annotate("then COLLAPSES / trench stalls\n(deep rate -> 0, not a sustained floor)",
            xy=(ar_stall, 0.16), xytext=(20.5, 0.45),
            fontsize=9.2, color='#7a241a', arrowprops=dict(arrowstyle='->', color='#c0392b'))
ax.annotate("experiment sustains a gradual\nconductance decline to ~0.20",
            xy=(32, 0.225), xytext=(24.5, 0.62),
            fontsize=9.2, color='#333', arrowprops=dict(arrowstyle='->', color='#888'))

fig.tight_layout()
out = "/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/deboer_match.png"
fig.savefig(out, dpi=150)
print("wrote", out)
