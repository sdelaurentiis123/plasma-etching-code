#!/usr/bin/env python3
"""Decisive convergence figure: normalized bottom etch rate at AR~8.6 vs ray count, for petch and
ViennaPS on the SAME sub-micron trench. ViennaPS is FLAT (converged at its default 1000 rays/pt); petch
converges UP to ~0.80. They settle ~0.06 apart -> a real model difference, NOT smoothing (formula is
identical to ViennaRay) and NOT under-sampling (both converged). The early 'petch no-smooth @40k = 0.735
~= ViennaPS' was an under-sampling artifact (open marker)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# nr at AR~8.6 (measured this session)
pe_rays = np.array([40000, 160000, 400000])
pe_full = np.array([0.787, 0.807, 0.795])      # full-smooth (default)
pe_none = np.array([0.735, np.nan, 0.794])     # no-smooth (40k = under-sampling artifact)
vp_rays = np.array([1000, 4000])
vp_nr = np.array([0.732, 0.732])

fig, ax = plt.subplots(figsize=(8, 5))
ax.axhline(0.732, color='#1a73e8', ls=':', lw=1, alpha=0.6)
ax.axhline(0.80, color='#0a8f3c', ls=':', lw=1, alpha=0.6)
ax.plot(vp_rays, vp_nr, 's-', color='#1a73e8', ms=11, lw=2, label='ViennaPS (converged at 1000 rays/pt)')
ax.plot(pe_rays, pe_full, '^-', color='#0a8f3c', ms=10, lw=2, label='petch, full smooth (converges → ~0.80)')
ax.plot([40000, 400000], [0.735, 0.794], 'o--', color='#9534c0', ms=9, lw=1.6, alpha=0.85,
        label='petch, no smooth')
ax.annotate('under-sampling artifact\n(few rays → fake-steep,\ncoincidentally ≈ ViennaPS)',
            xy=(40000, 0.735), xytext=(60000, 0.55), fontsize=9, color='#6a2c8f',
            arrowprops=dict(arrowstyle='->', color='#9534c0'))
ax.annotate('~0.06 real model gap\n(petch gentler at convergence)', xy=(400000, 0.765),
            xytext=(120000, 0.84), fontsize=9.5, color='#333',
            arrowprops=dict(arrowstyle='-[, widthB=2.0', color='#888'))

ax.set_xscale('log')
ax.set_xlabel('ray count  (petch: neutral rays  |  ViennaPS: rays/point)', fontsize=11)
ax.set_ylabel('normalized bottom etch rate at AR ≈ 8.6', fontsize=11)
ax.set_title('Trench ARDE convergence: ViennaPS flat (converged), petch converges ~0.06 gentler\n'
             '→ a real model difference, not smoothing (formula identical) or sampling', fontsize=10.5)
ax.set_ylim(0.5, 0.95); ax.legend(loc='upper left', fontsize=9.5); ax.grid(alpha=0.3, which='both')
fig.tight_layout()
out = "/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/trench_convergence.png"
fig.savefig(out, dpi=150); print("wrote", out)
