"""Figure: petch trench ARDE vs ViennaPS-GPU (ground truth), before/after faithful ion reflection."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AR = np.array([3.7, 6.1, 8.6])
vps = np.array([1.0, 0.861, 0.731])
before = np.array([1.0, 0.775, 0.468])          # petch, ions stick on first hit (no reflection)
after = np.array([1.0, 0.913, 0.602])           # petch, faithful ViennaPS ion reflection
after_e = np.array([0.0, 0.010, 0.017])         # seed-averaged std

fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.plot(AR, vps, 'o-', color='black', lw=2.4, ms=9, label='ViennaPS-GPU (ground truth, deterministic)')
ax.plot(AR, before, 's--', color='#c0392b', lw=1.8, ms=7, alpha=0.8,
        label='petch — ions stick on first hit  (RMSE 0.152)')
ax.errorbar(AR, after, yerr=after_e, fmt='^-', color='#2471a3', lw=2.2, ms=9, capsize=4,
            label='petch — faithful ion reflection  (RMSE 0.080)')
ax.set_xlabel('aspect ratio (depth / width)', fontsize=12)
ax.set_ylabel('normalized bottom etch rate', fontsize=12)
ax.set_title('Trench ARDE: petch vs ViennaPS  (SF6/O2, dx=0.04, same RTX 3090)', fontsize=12)
ax.legend(fontsize=10, loc='lower left')
ax.grid(alpha=0.3)
ax.set_ylim(0.3, 1.05)
ax.annotate('ViennaPS ion: sticking=0, coned-cosine reflection,\n'
            'energy loss per bounce until E<Eth  -> funnels ions\n'
            'down the trench to feed the deep floor.\n'
            'petch was missing this; now ported verbatim.',
            xy=(8.6, 0.731), xytext=(4.4, 0.45), fontsize=8.5,
            bbox=dict(boxstyle='round', fc='#fdf2e9', ec='#e67e22', alpha=0.9),
            arrowprops=dict(arrowstyle='->', color='#e67e22'))
fig.tight_layout()
fig.savefig("ionrefl_match.png", dpi=140)
print("wrote ionrefl_match.png")
