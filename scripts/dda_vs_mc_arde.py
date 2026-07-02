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
# All petch curves use the FAITHFUL reflected ion (2026-07-02 fix: the deterministic paths previously
# ignored flags.ion_reflection; earlier published values below in comments were legacy-ion).
# Panel A: W=0.5 um trench (MC vs DDA vs ViennaPS)
AR = np.array([2, 4, 6, 8, 10], float)
PETCH_MC  = np.array([0.536, 0.186, 0.077, 0.043, 0.029])     # under-samples deep floor (static)
PETCH_DDA = np.array([0.993, 0.967, 0.914, 0.836, 0.742])     # faithful ion (legacy-ion was 0.99..0.65)
# ViennaPS MEASURED static reference (2026-07-02 box run, driver 570 GPU_TRIANGLE): pre-carved
# MakeTrench (makeMask=False, carve verified d0==D), 0.15-min window, reflecting ions. This replaces
# the earlier single unverified "0.73 @ AR8.6" reference point.
VPS_STATIC = np.array([0.911, 0.820, 0.728, 0.626, 0.534])
# Panel B: W=2 um trench vs the de Boer wafer
AR2 = np.array([2, 5, 10, 15, 20], float)
PETCH_DDA_W2 = np.array([1.000, 1.005, 0.757, 0.486, 0.233])  # ballistic DDA, faithful ion
KN_AR = np.array([0, 5, 10, 15, 20, 30, 40], float)           # knudsen + faithful ion, wls=1.4
KN_NR = np.array([1.0, 0.722, 0.482, 0.348, 0.271, 0.188, 0.142])   # RMSE 0.040 vs wafer -> PASSES
DEBOER_AR = np.array([0.0, 10.0, 20.0, 40.0]); DEBOER_NR = np.array([1.0, 0.43, 0.29, 0.20])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))

axA.plot(AR, PETCH_DDA, "o-", color="#2471c7", lw=2.6, ms=8, label="petch DDA (deterministic)")
axA.plot(AR, VPS_STATIC, "*-", color="#16a085", lw=2.4, ms=13, label="ViennaPS-GPU (measured static)", zorder=6)
axA.plot(AR, PETCH_MC,  "s--", color="#c0392b", lw=2.2, ms=7, label="petch MC (200k rays, static)")
axA.annotate("MC under-samples the\ndeep floor → collapses", xy=(8, 0.043), xytext=(3.8, 0.25),
             color="#c0392b", fontsize=9, arrowprops=dict(arrowstyle="->", color="#c0392b"))
axA.annotate("DDA is noise-free but sits 0.08–0.21\nGENTLER than ViennaPS (grows with AR):\nits re-emission needs calibration",
             xy=(9, 0.58), xytext=(2.0, 0.44), color="#16a085", fontsize=9,
             arrowprops=dict(arrowstyle="->", color="#16a085"))
axA.set_xlim(1.5, 11); axA.set_title("W = 0.5 µm: DDA fixes MC under-sampling; gap to ViennaPS is its open calibration", fontsize=11)
axA.legend(loc="upper right", fontsize=9)

axB.plot(KN_AR, KN_NR, "o-", color="#1e8449", lw=2.8, ms=8,
         label="petch Knudsen + faithful ion (wls=1.4) — RMSE 0.040 ✓")
axB.plot(AR2, PETCH_DDA_W2, "o--", color="#2471c7", lw=2.0, ms=7, alpha=0.8,
         label="petch DDA (ballistic reference)")
axB.plot(DEBOER_AR, DEBOER_NR, "kP", ms=13, label="de Boer wafer (measured)", zorder=6)
axB.plot(DEBOER_AR, DEBOER_NR, "k:", lw=1.2, alpha=0.5)
axB.annotate("Knudsen conductance + AR-independent\nion floor passes the 0.05 wafer gate",
             xy=(20, 0.271), xytext=(16, 0.55), color="#1e8449", fontsize=9,
             arrowprops=dict(arrowstyle="->", color="#1e8449"))
axB.set_xlim(0, 42); axB.set_title("W = 2 µm: petch vs the real de Boer wafer", fontsize=11)
axB.legend(loc="upper right", fontsize=9)

for ax in (axA, axB):
    ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3)
axA.set_ylabel("normalized floor etch rate  $n_r = V_{floor}/V_{field}$")
fig.suptitle("petch deterministic transport (reconciled from plasma_sim) — SF$_6$/O$_2$ static ARDE, faithful reflected ion",
             fontweight="bold", fontsize=13)
plt.tight_layout()
p = os.path.join(HERE, "viz", "dda_vs_mc_arde.png")
os.makedirs(os.path.dirname(p), exist_ok=True); plt.savefig(p, dpi=150); print("saved", p)
np.savez(os.path.join(HERE, "dda_vs_mc_arde.npz"), ar=AR, petch_mc=PETCH_MC, petch_dda=PETCH_DDA,
         vps_static=VPS_STATIC, ar2=AR2, petch_dda_w2=PETCH_DDA_W2,
         kn_ar=KN_AR, kn_nr=KN_NR, deboer_ar=DEBOER_AR, deboer_nr=DEBOER_NR)
print("petch-DDA vs ViennaPS-static @ AR10:", PETCH_DDA[-1], "vs", VPS_STATIC[-1],
      f"(gap {PETCH_DDA[-1]-VPS_STATIC[-1]:.2f})")
