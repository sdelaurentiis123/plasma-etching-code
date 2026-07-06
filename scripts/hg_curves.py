"""The full HG curve comparison: every channel of Figs 2-6 (JAP 82,566) vs petch's derived-source
solver on the TRUE stack, correct per-surface definitions, zero knobs. HG values are direct figure
reads (+-0.02 flux / +-1 V read tolerance). Requires the c15c_AR*.npz sweep outputs."""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/private/tmp/claude-501/-Users-stanislavdelaurentiis-chip-etch/ab7b4b64-cc35-4337-a2e2-8d822f56aadd/scratchpad"
ARs = [1.0, 2.0, 3.0, 4.0]
p = {k: [] for k in ["e_bottom", "e_pr", "e_poly_in", "e_poly_out", "i_bottom", "i_poly_in",
                     "E_defl_face", "floorV", "corner", "edgeV", "neighV"]}
h = {k: [] for k in p}
for a in ARs:
    z = np.load(f"{D}/c15c_AR{a:g}.npz")
    for k in p: p[k].append(float(z[k]))
    try:
        zh = np.load(f"{D}/c15c_hg_AR{a:g}.npz")
        for k in h: h[k].append(float(zh[k]))
    except FileNotFoundError:
        h = None

# HG digitized (figure reads)
hg = dict(e_bottom=[.59,.40,.29,.22], e_pr=[.28,.41,.49,.57], e_poly_in=[.05,.04,.03,.02],
          e_poly_out=[.19,.20,.19,.19], i_bottom=[.59,.40,.29,.22], i_poly_in=[.15,.14,.14,.13],
          E_defl_face=[15,20,24,28], floorV=[8,17,26,33], edgeV=[1,4,6,7.5], neighV=[6,19,30,39])

fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.6))

axA = axes[0,0]
for k, lab, c in [("e_bottom","trench bottom","#2471c7"), ("e_pr","PR sidewalls","#e07b39"),
                  ("e_poly_out","poly outer","#3aa655"), ("e_poly_in","poly inner","#7b3fe0")]:
    axA.plot(ARs, p[k], "o-", color=c, lw=2.2, ms=7, label=f"petch {lab}")
    axA.plot(ARs, hg[k], "*--", color=c, ms=11, lw=1.1, alpha=0.75)
axA.set_title("Electron fluxes vs AR — solid: petch, stars: HG Fig 3", fontsize=10.5)
axA.set_xlabel("aspect ratio"); axA.set_ylabel("normalized electron flux")
axA.legend(fontsize=8); axA.grid(alpha=0.3)

axB = axes[0,1]
axB.plot(ARs, p["i_bottom"], "o-", color="#2471c7", lw=2.2, ms=7, label="petch bottom")
axB.plot(ARs, hg["i_bottom"], "*--", color="#2471c7", ms=11, lw=1.1, alpha=0.75, label="HG bottom")
axB.plot(ARs, p["i_poly_in"], "s-", color="#c0392b", lw=2.2, ms=7, label="petch poly-inner (notch foot)")
axB.plot(ARs, hg["i_poly_in"], "*--", color="#c0392b", ms=11, lw=1.1, alpha=0.75, label="HG poly-inner")
if h: axB.plot(ARs, h["i_bottom"], "o:", color="#2471c7", lw=1.4, ms=5, alpha=0.6, label="petch, HG e-convention")
axB.set_title("Ion fluxes vs AR — the foot current matches (Fig 4)", fontsize=10.5)
axB.set_xlabel("aspect ratio"); axB.set_ylabel("normalized ion flux")
axB.legend(fontsize=8.5); axB.grid(alpha=0.3)

axC = axes[1,0]
axC.plot(ARs, p["E_defl_face"], "o-", color="#c0392b", lw=2.4, ms=8, label="petch (poly-inner FACE impacts)")
axC.plot(ARs, hg["E_defl_face"], "*--", color="0.2", ms=12, lw=1.2, label="HG Fig 4 (right axis)")
axC.set_title("Deflected-ion energy at the notch foot — RISING trend\n(restored by the face-surface definition)", fontsize=10)
axC.set_xlabel("aspect ratio"); axC.set_ylabel("mean impact energy (eV)")
axC.legend(fontsize=8.5); axC.grid(alpha=0.3)

axD = axes[1,1]
for k, lab, c in [("floorV","floor center","#2471c7"), ("neighV","neighbor line","#7b3fe0"),
                  ("edgeV","edge line","#3aa655")]:
    axD.plot(ARs, p[k], "o-", color=c, lw=2.2, ms=7, label=f"petch {lab}")
    axD.plot(ARs, hg[k], "*--", color=c, ms=11, lw=1.1, alpha=0.75)
if h:
    for k, c in [("floorV","#2471c7"), ("neighV","#7b3fe0"), ("edgeV","#3aa655")]:
        axD.plot(ARs, h[k], "o:", color=c, lw=1.4, ms=5, alpha=0.6)
    axD.plot([], [], "o:", color="0.4", lw=1.4, ms=5, alpha=0.8, label="dotted: HG e-convention mode (brackets HG)")
axD.plot(ARs, p["corner"], "^-", color="#e07b39", lw=1.6, ms=6, label="petch PR corner (HG Fig 2: ≈−4.5)")
axD.axhline(-4.5, color="#e07b39", ls=":", lw=1.2)
axD.set_title("Potentials vs AR — solid: petch, stars: HG Figs 5/6\n(labels convention-sensitive; shapes are the physics)", fontsize=10)
axD.set_xlabel("aspect ratio"); axD.set_ylabel("potential (V)")
axD.legend(fontsize=8); axD.grid(alpha=0.3)

fig.suptitle("Hwang–Giapis 1997, all published curves vs petch's derived-source solver "
             "(true stack, zero knobs)", fontsize=12.5, y=0.995)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/hg_curves.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/hg_curves.png", dpi=130, bbox_inches="tight")
print("saved viz/hg_curves.png")
