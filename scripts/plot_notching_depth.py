#!/usr/bin/env python3
"""Figure: charging-driven notch at the poly/oxide junction. Left: notch profile cross-section
(the localized foot undercut). Right: notch depth vs AR -- petch (charging on/off) and HG's measured
shape. Reads notching_depth_result.npz."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load(os.path.join(os.path.dirname(__file__), "..", "notching_depth_result.npz"))
ar, non, noff = d['ar'], d['notch_on'], d['notch_off']
hg_ar, hg_notch = d['hg_ar'], d['hg_notch']
W = float(d['W']); dx = 0.25; corr = float(d['corr'])
xs = d['xs']; cx = 0.5 * (xs[0] + xs[-1])

fig, (axp, axd) = plt.subplots(1, 2, figsize=(12, 5))

# --- left: notch profile cross-section at the deepest resolved AR (AR4) ---
phi = d['phi_AR4']                          # (nx, nz)
zs = np.arange(phi.shape[1]) * dx           # reconstruct z (grid starts at 0, dx=0.25)
poly_th = 4.0 * W; OX = 1.5; sub_top = poly_th + OX; z_stop = OX
# shade gas (phi<0) vs solid; show only the trench/foot region
XX, ZZ = np.meshgrid(xs - cx, zs, indexing='ij')
axp.contourf(XX, ZZ, (phi < 0).astype(float), levels=[0.5, 1.5], colors=['#7fb3ff'], alpha=0.85)
axp.contour(XX, ZZ, phi, levels=[0.0], colors='k', linewidths=1.2)
axp.axhline(z_stop, color='#b5651d', lw=2.5, label='buried oxide (etch-stop)')
axp.axvline(W / 2, color='0.4', ls='--', lw=1, label='as-drawn wall')
axp.axvline(-W / 2, color='0.4', ls='--', lw=1)
axp.annotate('notch', xy=(non[-1] + W / 2, z_stop + 0.4), xytext=(W / 2 + 1.0, z_stop + 1.6),
             arrowprops=dict(arrowstyle='->', color='crimson'), color='crimson', fontweight='bold')
axp.set_xlabel('x from center (um)'); axp.set_ylabel('z (um)')
axp.set_title(f'Charging notch at poly/oxide foot (AR 4, W={W:g}um, 100% overetch)')
axp.set_ylim(0, sub_top + 0.5); axp.legend(loc='upper right', fontsize=8)

# --- right: notch depth vs AR, petch (twin axis for HG's different feature size) ---
axd.plot(ar, non, 'o-', color='crimson', lw=2, ms=8, label=f'petch, charging ON (W={W:g}um)')
axd.plot(ar, noff, 's--', color='0.5', lw=1.5, ms=6, label='petch, charging OFF')
axd.set_xlabel('aspect ratio at oxide'); axd.set_ylabel('petch notch depth (um)', color='crimson')
axd.tick_params(axis='y', labelcolor='crimson'); axd.set_ylim(0, max(non.max(), 0.2) * 1.4)
ax2 = axd.twinx()
ax2.plot(hg_ar, hg_notch, '^:', color='navy', lw=2, ms=8, label='Hwang-Giapis measured (W~0.5um)')
ax2.set_ylabel('HG measured notch depth (um)', color='navy')
ax2.tick_params(axis='y', labelcolor='navy'); ax2.set_ylim(0, hg_notch.max() * 1.4)
axd.set_title(f'Notch vs AR: monotone rise (Fujiwara), shape r={corr:.2f} vs HG')
l1, la1 = axd.get_legend_handles_labels(); l2, la2 = ax2.get_legend_handles_labels()
axd.legend(l1 + l2, la1 + la2, loc='upper left', fontsize=8)
axd.text(0.5, -0.16, 'Charging OFF -> no notch (perfectly anisotropic). Absolute um uncalibrated; '
         'AR1 notch below dx=0.25 grid resolution.', transform=axd.transAxes, ha='center', fontsize=7.5, color='0.3')

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "..", "viz", "notching_depth.png")
plt.savefig(out, dpi=130, bbox_inches='tight')
print("wrote", out)
