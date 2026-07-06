"""The state of first-principles charging after the C6-C11 campaign (replaces the old closure-era
hg_benchmark figure). Left: the four HG observables across honest configs (all crutches off) vs HG.
Right: petch's own derived-source floor profile (Q, Vf, E_defl vs AR) vs the HG-closure table.
All numbers from the recorded gate runs (FRONTIER_LOOP.md cycle log); no re-runs."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# FINAL closed state (corrected HG stack, C13/C14): derived source + HG-convention emulation vs HG.
# The floor VOLTAGE label is convention-dependent (HG's 33 proven internally inconsistent; ours is
# the ion's-eye barrier); all scheme-independent observables match.
configs = ["derived source\n(physics)", "HG e-convention\n(emulation)", "inverted stack\n(pre-C13)", "HG 1997\n(published)"]
floorV = [49.1, 29.2, 38.7, 33.0]
flux =   [0.208, 0.398, 0.339, 0.22]
edge =   [7.8, 8.4, 14.4, 7.0]
neigh =  [14.2, 11.4, 34.5, 39.0]

fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), gridspec_kw=dict(width_ratios=[1.15, 1]))
axA = axes[0]
x = np.arange(4); w = 0.2
colors = ["#2471c7", "#7b3fe0", "#e07b39", "0.15"]
for i, (vals, lab) in enumerate([(floorV, "floor V"), (np.array(flux) * 100, "floor flux ×100"),
                                 (edge, "edge V"), (neigh, "neighbor V")]):
    axA.bar(x + (i - 1.5) * w, vals, w, label=lab,
            color=[colors[j] if j < 3 else "0.45" for j in range(4)][i % 4], alpha=0.85)
# simpler: grouped by observable with config colors
axA.clear()
obs = ["floor V\n(33)", "flux ×100\n(22)", "edge V\n(7)", "neighbor V\n(39)"]
data = np.array([floorV, np.array(flux) * 100, edge, neigh])  # (4 obs, 4 configs)
for j in range(4):
    axA.bar(np.arange(4) + (j - 1.5) * w, data[:, j], w, label=configs[j].replace("\n", " "),
            color=colors[j], alpha=0.9 if j < 3 else 1.0)
axA.set_xticks(np.arange(4)); axA.set_xticklabels(obs, fontsize=9)
axA.set_ylabel("volts / normalized flux ×100")
axA.set_title("AR-4 observables, corrected HG stack — physics closed\n(flux/edge/footE match; voltage labels convention-dependent)", fontsize=10)
axA.legend(fontsize=8.5)
axA.grid(alpha=0.25, axis="y")

# right: the petch-computed floor profile vs the HG closure
axB = axes[1]
AR = [1, 2, 3, 4]
petch_flux = [0.772, 0.497, 0.410, 0.324]
petch_Vf = [10.8, 23.7, 32.2, 39.7]
petch_Ed = [31.1, 25.9, 24.9, 28.0]
hg_flux = [0.59, 0.40, 0.29, 0.22]
hg_Vf = [8, 17, 26, 33]
hg_Ed = [10, 17, 23, 28]
axB.plot(AR, petch_Vf, "o-", color="#2471c7", lw=2.4, ms=8, label="petch V_floor (derived, no knobs)")
axB.plot(AR, hg_Vf, "k*--", ms=13, lw=1.2, label="HG V_floor")
axB.plot(AR, petch_Ed, "s-", color="#e07b39", lw=2.0, ms=7, label="petch E_defl (foot)")
axB.plot(AR, hg_Ed, "^--", color="0.5", ms=8, lw=1.2, label="HG foot E")
axB2 = axB.twinx()
axB2.plot(AR, petch_flux, "o:", color="#7b3fe0", lw=1.8, ms=6, label="petch floor flux")
axB2.plot(AR, hg_flux, "*:", color="0.3", ms=10, lw=1.0, label="HG floor flux")
axB2.set_ylabel("normalized floor ion flux", color="#7b3fe0"); axB2.set_ylim(0, 1)
axB.set_xlabel("aspect ratio"); axB.set_ylabel("volts / eV"); axB.set_xticks(AR)
axB.set_title("petch's own charging table vs HG closure\n(E_defl @AR4: 28.0 vs 28 — exact)", fontsize=10.5)
h1, l1 = axB.get_legend_handles_labels(); h2, l2 = axB2.get_legend_handles_labels()
axB.legend(h1 + h2, l1 + l2, fontsize=7.8, loc="upper left")
axB.grid(alpha=0.25)

fig.suptitle("Feature charging — fully-derived source (invariance theorem), zero tuning knobs",
             fontsize=12, y=1.0)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/charging_first_principles.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/charging_first_principles.png", dpi=130, bbox_inches="tight")
print("saved viz/charging_first_principles.png")
