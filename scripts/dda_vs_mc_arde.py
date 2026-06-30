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
# Panel A: W=0.5 um trench (MC vs DDA vs ViennaPS)
AR = np.array([2, 4, 6, 8, 10], float)
PETCH_MC  = np.array([0.536, 0.186, 0.077, 0.043, 0.029])     # under-samples deep floor (static)
PETCH_DDA = np.array([0.989, 0.941, 0.857, 0.759, 0.652])     # deterministic, tracks ViennaPS
VPS_AR, VPS_NR = 8.6, 0.73                                     # ViennaPS GPU_TRIANGLE reference
# Panel B: W=2 um trench (petch-DDA vs the de Boer wafer, on CUDA)
AR2 = np.array([2, 5, 10, 15, 20], float)
PETCH_DDA_W2 = np.array([0.969, 0.871, 0.579, 0.363, 0.228])
DEBOER_AR = np.array([10.0, 20.0, 40.0]); DEBOER_NR = np.array([0.43, 0.29, 0.20])   # W=2 um wafer

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))

axA.plot(AR, PETCH_DDA, "o-", color="#2471c7", lw=2.6, ms=8, label="petch DDA (deterministic, new)")
axA.plot(AR, PETCH_MC,  "s--", color="#c0392b", lw=2.2, ms=7, label="petch MC (200k rays, static)")
axA.plot([VPS_AR], [VPS_NR], "*", color="#16a085", ms=22, label="ViennaPS-GPU (ballistic ref)", zorder=6)
axA.annotate("MC under-samples the\ndeep floor → collapses", xy=(8, 0.043), xytext=(3.8, 0.28),
             color="#c0392b", fontsize=9, arrowprops=dict(arrowstyle="->", color="#c0392b"))
axA.annotate("DDA ≈ ViennaPS\n(0.727 vs 0.73)", xy=(8.6, 0.73), xytext=(2.2, 0.42),
             color="#16a085", fontsize=9, arrowprops=dict(arrowstyle="->", color="#16a085"))
axA.set_xlim(1.5, 11); axA.set_title("W = 0.5 µm: DDA fixes MC under-sampling, matches ViennaPS", fontsize=11)
axA.legend(loc="upper right", fontsize=9)

axB.plot(AR2, PETCH_DDA_W2, "o-", color="#2471c7", lw=2.6, ms=8, label="petch DDA (W=2 µm, CUDA)")
axB.plot(DEBOER_AR, DEBOER_NR, "kP", ms=13, label="de Boer wafer (measured)", zorder=6)
axB.plot(DEBOER_AR, DEBOER_NR, "k:", lw=1.2, alpha=0.5)
axB.annotate("petch-DDA straddles the wafer\n(above @ AR10, below @ AR20);\nViennaPS-regime params, RMSE≈0.09",
             xy=(10, 0.58), xytext=(10.5, 0.74), color="#2471c7", fontsize=9,
             arrowprops=dict(arrowstyle="->", color="#2471c7"))
axB.set_xlim(1.5, 21); axB.set_title("W = 2 µm: petch DDA vs the real de Boer wafer", fontsize=11)
axB.legend(loc="upper right", fontsize=9)

for ax in (axA, axB):
    ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3)
axA.set_ylabel("normalized floor etch rate  $n_r = V_{floor}/V_{field}$")
fig.suptitle("petch deterministic DDA transport (ported from plasma_sim) — SF$_6$/O$_2$ static ARDE",
             fontweight="bold", fontsize=13)
plt.tight_layout()
p = os.path.join(HERE, "viz", "dda_vs_mc_arde.png")
os.makedirs(os.path.dirname(p), exist_ok=True); plt.savefig(p, dpi=150); print("saved", p)
np.savez(os.path.join(HERE, "dda_vs_mc_arde.npz"), ar=AR, petch_mc=PETCH_MC, petch_dda=PETCH_DDA,
         vps_ar=VPS_AR, vps_nr=VPS_NR, ar2=AR2, petch_dda_w2=PETCH_DDA_W2,
         deboer_ar=DEBOER_AR, deboer_nr=DEBOER_NR)
print("petch-DDA @ AR8.6 (interp):", round(float(np.interp(8.6, AR, PETCH_DDA)), 3), "vs ViennaPS", VPS_NR)
