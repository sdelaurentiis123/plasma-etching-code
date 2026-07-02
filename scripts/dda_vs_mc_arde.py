"""Static instantaneous ARDE nr(AR)=V_floor/V_field: petch's deterministic neutral transports vs
the MEASURED ViennaPS static reference and the de Boer wafer. All petch curves use the faithful
reflected ion.

DATA PROVENANCE (2026-07-02/03, measured — do not edit without re-running):
- VPS_STATIC: ViennaPS GPU_TRIANGLE, driver-570 box, pre-carved MakeTrench (makeMask=False, carve
  verified d0==D), 0.15-min window, reflecting ions. W=0.5um trench.
- RADIO: petch neutral_transport='radiosity' + radiosity_solver='gmres' (mesh form-factor,
  n_ff=128), scripts/dda_static_gate.py -> RMSE 0.043 vs VPS_STATIC = PASSES the 0.05 gate.
- PETCH_DDA: the grid-march DDA AFTER its five transport fixes (sky double-count, area-mean
  radiosity, sub-cell walls, solid-seeking remit, full-sphere quadrature). Too STEEP in the
  passivated-wall (albedo~0.99) regime — a documented structural limit of the remit-field
  representation (per-hop grazing losses compound over ~50-bounce duct cascades; the form-factor
  radiosity handles that regime exactly). EARLIER published DDA values (0.99..0.74, "0.727 vs
  0.73", "0.08-0.21 gentler") were computed with a sky double-counting bug and are RETRACTED.
- PETCH_MC: 200k rays, single static eval (deep-floor under-sampling collapse, unchanged).
- KN: knudsen + faithful ion, wls=1.4 static harness (the wafer-gate lineage; unchanged).
- RADIO_W2: radiosity+gmres ballistic reference at the de Boer width.
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Panel A: W = 0.5 um trench
AR = np.array([2, 4, 6, 8, 10], float)
VPS_STATIC = np.array([0.911, 0.820, 0.728, 0.626, 0.534])    # measured ViennaPS reference
RADIO = np.array([0.908, 0.826, 0.755, 0.677, 0.609])         # petch radiosity+GMRES: RMSE 0.043 PASS
PETCH_DDA = np.array([0.833, 0.635, 0.462, 0.332, 0.234])     # DDA post-fix: passivated-regime limit
PETCH_MC = np.array([0.536, 0.186, 0.077, 0.043, 0.029])      # MC single static eval: starves
# Panel B: W = 2 um (de Boer width)
KN_AR = np.array([0, 5, 10, 15, 20, 30, 40], float)
KN_NR = np.array([1.0, 0.722, 0.482, 0.348, 0.271, 0.188, 0.142])   # knudsen static (wafer lineage)
RADIO_W2_AR = np.array([2, 5, 10, 15, 20], float)
RADIO_W2 = np.array([0.884, 0.787, 0.656, 0.546, 0.445])      # ballistic radiosity reference
DEBOER_AR = np.array([0.0, 10.0, 20.0, 40.0]); DEBOER_NR = np.array([1.0, 0.43, 0.29, 0.20])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))

axA.plot(AR, VPS_STATIC, "*-", color="#16a085", lw=2.4, ms=13, label="ViennaPS-GPU (measured static)", zorder=6)
axA.plot(AR, RADIO, "o-", color="#2471c7", lw=2.6, ms=8, label="petch radiosity+GMRES — RMSE 0.043 ✓")
axA.plot(AR, PETCH_DDA, "^--", color="#8e44ad", lw=1.8, ms=7, alpha=0.85,
         label="petch DDA (post-fix; passivated-regime limit)")
axA.plot(AR, PETCH_MC, "s--", color="#c0392b", lw=1.8, ms=7, alpha=0.85, label="petch MC (200k rays, static)")
axA.annotate("radiosity+GMRES passes the\nmeasured-ViennaPS gate (0.043)",
             xy=(6, 0.755), xytext=(2.2, 0.36), color="#2471c7", fontsize=9,
             arrowprops=dict(arrowstyle="->", color="#2471c7"))
axA.annotate("MC starves; DDA too steep in the\nmirror-wall regime (documented limit)",
             xy=(8, 0.33), xytext=(5.4, 0.12), color="0.35", fontsize=8.5,
             arrowprops=dict(arrowstyle="->", color="0.5"))
axA.set_xlim(1.5, 11); axA.set_title("W = 0.5 µm: deterministic transports vs measured ViennaPS", fontsize=11)
axA.legend(loc="upper right", fontsize=8.5)

axB.plot(KN_AR, KN_NR, "o-", color="#1e8449", lw=2.4, ms=7,
         label="petch Knudsen + reflected ion (static)")
axB.plot(RADIO_W2_AR, RADIO_W2, "o--", color="#2471c7", lw=2.0, ms=7, alpha=0.85,
         label="petch ballistic radiosity (reference)")
axB.plot(DEBOER_AR, DEBOER_NR, "kP", ms=13, label="de Boer wafer (measured)", zorder=6)
axB.plot(DEBOER_AR, DEBOER_NR, "k:", lw=1.2, alpha=0.5)
axB.annotate("ballistic sits above the wafer;\nKnudsen conductance + ion floor\ncarry the wafer physics",
             xy=(20, 0.271), xytext=(16, 0.6), color="#1e8449", fontsize=9,
             arrowprops=dict(arrowstyle="->", color="#1e8449"))
axB.set_xlim(0, 42); axB.set_title("W = 2 µm: petch vs the real de Boer wafer", fontsize=11)
axB.legend(loc="upper right", fontsize=9)

for ax in (axA, axB):
    ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3)
axA.set_ylabel("normalized floor etch rate  $n_r = V_{floor}/V_{field}$")
fig.suptitle("petch deterministic transport — SF$_6$/O$_2$ static ARDE, faithful reflected ion",
             fontweight="bold", fontsize=13)
plt.tight_layout()
p = os.path.join(HERE, "viz", "dda_vs_mc_arde.png")
os.makedirs(os.path.dirname(p), exist_ok=True); plt.savefig(p, dpi=150); print("saved", p)
np.savez(os.path.join(HERE, "dda_vs_mc_arde.npz"), ar=AR, vps_static=VPS_STATIC, radio=RADIO,
         petch_dda=PETCH_DDA, petch_mc=PETCH_MC, kn_ar=KN_AR, kn_nr=KN_NR,
         radio_w2_ar=RADIO_W2_AR, radio_w2=RADIO_W2, deboer_ar=DEBOER_AR, deboer_nr=DEBOER_NR)
print("radiosity RMSE vs measured ViennaPS:", round(float(np.sqrt(np.mean((RADIO - VPS_STATIC) ** 2))), 3))
