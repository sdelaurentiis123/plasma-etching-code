"""Static instantaneous ARDE nr(AR)=V_floor/V_field on a clean W=0.5 um trench: petch's new
deterministic DDA transport vs petch MC, with the ViennaPS and de Boer-wafer references. Shows
(1) DDA fixes the MC deep-floor under-sampling (single-static-eval collapse), (2) DDA reproduces
ViennaPS's ballistic rolloff (~0.73 @ AR8.6), (3) both ballistic engines sit ABOVE the real wafer.

Measured 2026-06-29 (reconciliation run). Sources:
  petch-MC : arde_static.py on RTX 3090 (CUDA), 200k rays, ion_reflection, dx=0.04, static.
  petch-DDA: arde_dda (local), neutral_transport='dda', 64 dirs, 40k ion rays, ion_reflection, dx=0.05.
  ViennaPS : SF6O2Etching GPU_TRIANGLE deep-AR ~0.73 @ AR~8.6 (documented arde_static reference).
  de Boer  : Jansen/de Boer/Blauw cryo SF6/O2 wafer (W=2 um trench).
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AR = np.array([2, 4, 6, 8, 10], float)
PETCH_MC  = np.array([0.536, 0.186, 0.077, 0.043, 0.029])     # under-samples deep floor (static)
PETCH_DDA = np.array([0.989, 0.941, 0.857, 0.759, 0.652])     # deterministic, tracks ViennaPS
VPS_AR, VPS_NR = 8.6, 0.73                                     # ViennaPS GPU_TRIANGLE reference
DEBOER_AR = np.array([10.0, 20.0, 40.0]); DEBOER_NR = np.array([0.43, 0.29, 0.20])   # W=2 um wafer

fig, ax = plt.subplots(figsize=(8.6, 5.6))
ax.plot(AR, PETCH_DDA, "o-", color="#2471c7", lw=2.6, ms=8, label="petch DDA (deterministic, new)")
ax.plot(AR, PETCH_MC,  "s--", color="#c0392b", lw=2.2, ms=7, label="petch MC (200k rays, static)")
ax.plot([VPS_AR], [VPS_NR], "*", color="#16a085", ms=22, label="ViennaPS-GPU (ballistic ref)", zorder=6)
ax.plot(DEBOER_AR, DEBOER_NR, "kP", ms=12, label="de Boer wafer (measured, W=2 µm)", zorder=6)
ax.plot(DEBOER_AR, DEBOER_NR, "k:", lw=1.2, alpha=0.5)
ax.annotate("MC under-samples the deep floor\n(single static eval) → collapses",
            xy=(8, 0.043), xytext=(4.4, 0.30), color="#c0392b", fontsize=9.5,
            arrowprops=dict(arrowstyle="->", color="#c0392b"))
ax.annotate("DDA reproduces ViennaPS's\nballistic rolloff (≈0.71 vs 0.73)",
            xy=(8.6, 0.73), xytext=(2.3, 0.45), color="#16a085", fontsize=9.5,
            arrowprops=dict(arrowstyle="->", color="#16a085"))
ax.annotate("both ballistic engines sit\nABOVE the real wafer",
            xy=(10, 0.43), xytext=(6.0, 0.16), color="0.25", fontsize=9.5,
            arrowprops=dict(arrowstyle="->", color="0.4"))
ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylabel("normalized floor etch rate  $n_r = V_{floor}/V_{field}$")
ax.set_xlim(1.5, 11); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=10)
ax.set_title("Static ARDE — deterministic DDA fixes petch-MC's deep-floor under-sampling\n"
             "and matches ViennaPS (W = 0.5 µm trench, SF$_6$/O$_2$)", fontsize=12, fontweight="bold")
plt.tight_layout()
for p in [os.path.join(HERE, "viz", "dda_vs_mc_arde.png")]:
    os.makedirs(os.path.dirname(p), exist_ok=True); plt.savefig(p, dpi=150); print("saved", p)
np.savez(os.path.join(HERE, "dda_vs_mc_arde.npz"), ar=AR, petch_mc=PETCH_MC, petch_dda=PETCH_DDA,
         vps_ar=VPS_AR, vps_nr=VPS_NR, deboer_ar=DEBOER_AR, deboer_nr=DEBOER_NR)
print("petch-DDA @ AR8.6 (interp):", round(float(np.interp(8.6, AR, PETCH_DDA)), 3), "vs ViennaPS", VPS_NR)
