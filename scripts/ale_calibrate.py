"""C4 gate + viz: DIFFERENTIABLE ALE. Reverse-mode gradients through the whole cyclic site-balance
chemistry (src/petch/ale_diff.py, torch), and gradient-based inverse design of the ion energy for a
target etch-per-cycle. The moat: no open feature-scale etch tool exposes chemistry gradients.
Left: EPC(E) with autograd tangents (dEPC/dE). Right: Newton inversion path to a target EPC.
Saves viz/ale_calibrate.png (+ docs/ copy)."""
import sys, os; sys.path.insert(0, 'src')
import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from petch.ale_diff import epc_torch, dEPC_dE, invert_energy_for_epc

Es = np.linspace(15, 30, 31)
epc = np.array([float(epc_torch(torch.tensor(e, dtype=torch.float64)).detach()) for e in Es])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.8, 5.0))

# Panel A: EPC curve + autograd tangent lines (process sensitivity dEPC/dE)
axA.axvspan(15, 20, color="#bfe3c0", alpha=0.5, label="ALE window")
axA.plot(Es, epc, "-", color="#2471c7", lw=2.6, label="EPC(E)  (differentiable model)")
for E0 in (17.5, 22.5, 27.5):
    e0, g = dEPC_dE(E0)
    xs = np.array([E0 - 1.6, E0 + 1.6])
    axA.plot(xs, e0 + g * (xs - E0), "--", color="#e07b39", lw=2.0)
    axA.plot([E0], [e0], "o", color="#e07b39", ms=8)
    axA.annotate(f"dEPC/dE={g:.3f}", (E0, e0), textcoords="offset points", xytext=(6, -14), fontsize=8.5, color="#a8500f")
axA.set_xlabel("Ar⁺ ion energy (eV)"); axA.set_ylabel("etch per cycle (Å/cycle)")
axA.set_title("Reverse-mode gradients dEPC/dE (autograd = tangents)", fontsize=11)
axA.legend(fontsize=9, loc="upper left"); axA.grid(alpha=0.3)

# Panel B: inverse design -- Newton descent on the loss to hit a target EPC
target = 1.0
E_sol, epc_sol, hist = invert_energy_for_epc(target, E0=24.0)
loss = (epc - target) ** 2
axB.plot(Es, loss, "-", color="0.5", lw=2.0, label="loss = (EPC−target)²")
hx = [h[0] for h in hist]; hy = [(h[1] - target) ** 2 for h in hist]
axB.plot(hx, hy, "o-", color="#7b3fe0", lw=1.6, ms=6, label=f"Newton path (E₀=24 → {E_sol:.2f} eV)")
axB.axvline(E_sol, color="#7b3fe0", ls=":", lw=1.4)
axB.set_xlabel("Ar⁺ ion energy (eV)"); axB.set_ylabel("calibration loss")
axB.set_title(f"Inverse design: solve E for target EPC={target} Å/cyc → {E_sol:.2f} eV", fontsize=11)
axB.legend(fontsize=9); axB.grid(alpha=0.3)

fig.suptitle("C4 — differentiable feature-scale ALE chemistry (gradients + inverse design, the moat)", fontsize=11.5, y=1.0)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/ale_calibrate.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/ale_calibrate.png", dpi=130, bbox_inches="tight")
print(f"solved E={E_sol:.3f} eV for target EPC={target} (epc there = {epc_sol:.3f})", flush=True)
print("saved viz/ale_calibrate.png", flush=True)
print("DONE", flush=True)
