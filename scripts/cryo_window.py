"""C5 gate + viz: cryogenic etch enhancement from petch's physisorption model vs published data.
Left: SiO2 etch rate vs temperature (petch curve + CF4/H2 pseudo-wet anchors, Small Methods 2024),
with the physisorption coverage theta(T) that drives it. Right: the etch-rate enhancement factor with
the independent HF cryo-ALE cross-check (~3.2x at -60C, different system). Saves viz/cryo_window.png."""
import sys, os; sys.path.insert(0, 'src')
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from petch.cryo import cryo_etch_rate, physisorption_coverage, enhancement_factor

T = np.linspace(-120, 30, 151)
er = cryo_etch_rate(T)
th = physisorption_coverage(T)
enh = enhancement_factor(T)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.8, 5.0))

# Panel A: etch rate vs T + physisorption coverage
axA.axvspan(-60, -40, color="#cfe3ff", alpha=0.5, label="condensation onset band")
axA.plot(T, er, "-", color="#2471c7", lw=2.6, label="petch cryo model")
axA.plot([20, -60], [2.3, 3.76], "k*", ms=16, zorder=6,
         label="CF4/H2 pseudo-wet (Small Methods 2024)")
axA.annotate("3.76 nm/s\n(1.6× warm)", (-60, 3.76), textcoords="offset points", xytext=(8, -6), fontsize=8.5)
axA.annotate("2.3 nm/s plateau", (20, 2.3), textcoords="offset points", xytext=(-70, 6), fontsize=8.5)
axA.set_xlabel("substrate temperature (°C)"); axA.set_ylabel("SiO₂ etch rate (nm/s)", color="#2471c7")
axA.set_title("Cryo etch enhancement — petch vs CF4/H2 pseudo-wet", fontsize=11)
axA.legend(fontsize=8.5, loc="upper right"); axA.grid(alpha=0.3)
ax2 = axA.twinx()
ax2.plot(T, th, "--", color="#e07b39", lw=1.8)
ax2.set_ylabel("physisorption coverage θ(T)", color="#e07b39"); ax2.set_ylim(0, 1.05)

# Panel B: enhancement factor + independent HF cryo-ALE cross-check
axB.axhline(1.0, color="0.6", ls=":", lw=1.3)
axB.plot(T, enh, "-", color="#7b3fe0", lw=2.6, label="petch ER(T)/R_base")
axB.plot([20, -60], [1.0, 1.6], "k*", ms=14, label="CF4/H2 anchor (1.6× @ −60 °C)")
axB.plot([20, -60], [1.0, 3.2], "^", color="0.45", ms=9, label="HF cryo-ALE cross-check (3.2×, diff. system)")
axB.set_xlabel("substrate temperature (°C)"); axB.set_ylabel("etch-rate enhancement (×)")
axB.set_title("Enhancement vs T (E_ads=0.4 eV fixed from measurement)", fontsize=11)
axB.legend(fontsize=8.5, loc="upper right"); axB.grid(alpha=0.3)

fig.suptitle("C5 — cryogenic etch chemistry (physisorption, differentiable, unoccupied in open source)", fontsize=11.5, y=1.0)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/cryo_window.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/cryo_window.png", dpi=130, bbox_inches="tight")
print("ER(+20)=%.2f  ER(-60)=%.2f  ratio=%.2f" % (cryo_etch_rate(20.0), cryo_etch_rate(-60.0), cryo_etch_rate(-60.0)/cryo_etch_rate(20.0)), flush=True)
print("saved viz/cryo_window.png", flush=True)
print("DONE", flush=True)
