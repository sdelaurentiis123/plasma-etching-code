"""C3 gate + viz: the Si/Cl2/Ar+ ALE window from petch's site-balance ROM vs Vella-Graves 2025.
Left: EPC vs Ar+ energy with the 15-20 eV ALE window shaded, petch vs ROM (Fig 11) + MD (Fig 8).
Right: ALE synergy S = (EPC-alpha-beta)/EPC -> ~100% inside the window, collapsing outside.
Saves viz/ale_window.png (+ docs/ copy). Everything from first principles; nothing tuned."""
import sys, os; sys.path.insert(0, 'src')
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from petch.ale import run_ale, synergy

E = np.array([15, 17.5, 20, 22.5, 25, 27.5, 30])
rom = np.array([0.6, 0.7, 0.9, 1.4, 1.9, 4.1, 4.8])   # Vella-Graves ROM, Fig 11 (visual read)
md = np.array([0.6, 0.75, 0.95, 1.8, 2.45, 3.7, 6.2])  # their MD, Fig 8 (visual read)
petch = np.array([run_ale(float(e))["epc"] for e in E])
syn = np.array([synergy(float(e))["synergy"] for e in E])
for e, p, r in zip(E, petch, rom):
    print(f"E={e:>5} eV: petch EPC={p:.3f}  ROM={r:.2f} A/cyc", flush=True)

fig, (axE, axS) = plt.subplots(1, 2, figsize=(12.8, 5.0))
axE.axvspan(15, 20, color="#bfe3c0", alpha=0.5, label="ALE window (15–20 eV)")
axE.axhline(1.3585, color="0.5", ls=":", lw=1.3, label="1 monolayer (1.36 Å)")
axE.plot(E, md, "^--", color="0.45", lw=1.6, ms=7, label="Vella–Graves MD (Fig 8)")
axE.plot(E, rom, "s--", color="#e07b39", lw=1.8, ms=7, label="Vella–Graves ROM (Fig 11)")
axE.plot(E, petch, "o-", color="#2471c7", lw=2.6, ms=8, label="petch ale.py (site-balance ROM)")
axE.set_xlabel("Ar⁺ ion energy (eV)"); axE.set_ylabel("etch per cycle (Å/cycle)")
axE.set_title("Si/Cl₂/Ar⁺ ALE window — petch vs Vella–Graves 2025", fontsize=11)
axE.legend(fontsize=8.5, loc="upper left"); axE.grid(alpha=0.3)

axS.axvspan(15, 20, color="#bfe3c0", alpha=0.5)
axS.plot(E, syn * 100, "o-", color="#7b3fe0", lw=2.6, ms=8)
axS.axhline(100, color="0.5", ls=":", lw=1.3, label="ideal ALE (100%)")
axS.set_xlabel("Ar⁺ ion energy (eV)"); axS.set_ylabel("ALE synergy S (%)")
axS.set_title("Synergy → 100% in the window, collapses outside (Kanarik)", fontsize=11)
axS.set_ylim(0, 110); axS.legend(fontsize=9); axS.grid(alpha=0.3)

fig.suptitle("C3 — open, differentiable, feature-scale ALE (unoccupied in open source)", fontsize=11.5, y=1.0)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/ale_window.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/ale_window.png", dpi=130, bbox_inches="tight")
print("saved viz/ale_window.png", flush=True)
print("DONE", flush=True)
