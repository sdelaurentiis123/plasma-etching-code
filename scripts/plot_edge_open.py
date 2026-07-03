#!/usr/bin/env python3
"""Figure: edge/open-area current balance diagnostics for the HG edge-line cell.

Reads charging_gate_result.npz and notching_gate_result.npz, written by the gate scripts.
The plot makes the distinction explicit: HG Fig. 3 reports gross outer-side electron flux, while
the conductor update uses net open-side electron surplus after the modeled ion counterflux.
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cg = np.load(os.path.join(HERE, "charging_gate_result.npz"))
ng = np.load(os.path.join(HERE, "notching_gate_result.npz"))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.2))

ar = cg["ar"]
axA.plot(ar, cg["edge_hg_electron_gross"], "k*--", ms=12, lw=1.2,
         label="HG Fig. 3 outer electron gross")
axA.plot(ar, cg["edge_electron_gross"], "o-", lw=2.3, label="modeled outer electron gross")
axA.plot(ar, cg["edge_ion_gross"], "s-", lw=2.0, label="modeled outer ion gross")
axA.plot(ar, cg["edge_net_electron"], "^-", lw=2.0, label="net electron surplus")
axA.set_xlabel("aspect ratio")
axA.set_ylabel("flux normalized to trench opening")
axA.grid(alpha=0.3)
axA.legend(fontsize=9)
axA.set_title("Open-side current budget")

v_gap = ng["vpoly"] - ng["vedge"]
axB.plot(ng["ar"], ng["vedge"], "o-", lw=2.3, label="edge-line poly")
axB.plot(ng["ar"], ng["vpoly"], "s-", lw=2.3, label="neighbor poly")
axB.plot(ng["ar"], v_gap, "^-", lw=2.0, label="neighbor - edge")
axB.set_xlabel("aspect ratio")
axB.set_ylabel("potential (V)")
axB.grid(alpha=0.3)
axB.legend(fontsize=9)
axB.set_title("Line-to-line potential split")

model = str(cg["edge_open_model"]) if "edge_open_model" in cg else "unknown"
fig.suptitle(f"HG edge/open diagnostics ({model})", fontweight="bold", fontsize=13)
plt.tight_layout()
p = os.path.join(HERE, "viz", "edge_open_current.png")
plt.savefig(p, dpi=150)
print("saved", p)
