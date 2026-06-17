"""Floor demo figure: knee config +/- the bimodal ion core, vs the de Boer experiment. Shows the core
lifts the mid-AR floor (the mechanism works) but cannot sustain de Boer's flat 0.25 tail -- the residual
is a REACTION-MODEL limit (the Belen rate's dominant ion-enhanced term dies with coverage at depth;
the only coverage-independent channel, physical sputter, is ~200x smaller -> structurally tiny floor)."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load(sys.argv[1] if len(sys.argv) > 1 else "deboer_floor.npz")
ar_n, nr_n, ar_c, nr_c = d['ar_n'], d['nr_n'], d['ar_c'], d['nr_c']
EXP_AR, EXP_R = d['exp_ar'], d['exp_r']


def cl(ar, nr):
    ok = ~np.isnan(nr); ar, nr = ar[ok], nr[ok]
    if len(nr) >= 5:
        nr = np.convolve(np.clip(nr, 0, 1.05), np.ones(5) / 5, mode='same')
    return ar, np.clip(nr, 0, 1.02)


fig, ax = plt.subplots(figsize=(8.4, 5.4))
xr = np.linspace(0, 40, 200)
ax.plot(xr, np.interp(xr, EXP_AR, EXP_R), 'k-', lw=2.6, label='de Boer/Blauw experiment')
ax.plot(EXP_AR, EXP_R, 'ko', ms=8)
ax.plot(*cl(ar_n, nr_n), '--', color='#9534c0', lw=2, label='petch knee config, single-lobe IADF')
ax.plot(*cl(ar_c, nr_c), '-', color='#0a8f3c', lw=2.2, label='petch knee config + bimodal ion core')
ax.axvspan(13, 23, color='#2e7d32', alpha=0.08)

ax.set_xlim(0, 40); ax.set_ylim(0, 1.04)
ax.set_xlabel("aspect ratio  (depth / width)", fontsize=12)
ax.set_ylabel("normalized bottom etch rate", fontsize=12)
ax.set_title("Bimodal ion core lifts the mid-AR floor (mechanism works) but cannot sustain de Boer's\n"
             "flat tail: the residual is a REACTION-MODEL limit (no large coverage-independent channel)",
             fontsize=10.3)
ax.legend(loc='upper right', fontsize=10)
ax.grid(alpha=0.3)
ax.annotate("core delivers ions to the\nbottom -> floor lifted here",
            xy=(20, 0.21), xytext=(21.5, 0.45),
            fontsize=9.2, color='#1b5e20', arrowprops=dict(arrowstyle='->', color='#0a8f3c'))
ax.annotate("experiment holds a sustained ~0.25 floor:\nimplies a coverage-INDEPENDENT etch channel\n"
            "the Belen rate (petch & ViennaPS) lacks",
            xy=(35, 0.21), xytext=(13.5, 0.66),
            fontsize=9.0, color='#333', arrowprops=dict(arrowstyle='->', color='#888'))

fig.tight_layout()
out = "/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/deboer_floor.png"
fig.savefig(out, dpi=150)
print("wrote", out)
