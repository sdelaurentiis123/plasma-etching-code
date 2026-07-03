#!/usr/bin/env python3
"""Figure: petch's 2-D Hwang-Giapis charging solver vs the published gate.
Left: floor ion flux vs AR (model curve vs the 8 digitized HG points + Matsui note).
Right: the steady-state 2-D potential map at AR=4 (the mechanism made visible:
positive floor/foot, negative mask tops, the in-trench field that deflects ions).
Reads charging_gate_result.npz (written by scripts/charging_gate.py)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from petch.charging2d import solve_edge_array_charging, solve_trench_charging

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = np.load(os.path.join(HERE, "charging_gate_result.npz"))

def scalar(key, default):
    return g[key].item() if key in g and getattr(g[key], "shape", ()) == () else default

solve_kwargs = dict(
    see_model=str(scalar("see_model", "none")),
    see_generations=int(scalar("see_generations", 1)),
    source_model=str(scalar("source_model", "analytic")),
    poly_mode=str(scalar("poly_mode", "tied")),
    poly_bias_V=float(scalar("poly_bias_V", 0.0)),
)
charging_geometry = str(scalar("charging_geometry", "trench"))
charging_W = int(scalar("charging_W", 32))
charging_mouth = int(scalar("charging_mouth", 237))
if "edge_open_model" in g:
    solve_kwargs["edge_open_model"] = str(scalar("edge_open_model", "none"))
if "edge_open_electron_flux" in g and float(g["edge_open_electron_flux"]) >= 0.0:
    solve_kwargs["edge_open_electron_flux"] = float(g["edge_open_electron_flux"])

print("solving AR=4 for the potential map (smooth=True, cosmetic only)...", flush=True)
if charging_geometry == "edge_array":
    edge_kwargs = dict(see_model=solve_kwargs["see_model"],
                       see_generations=solve_kwargs["see_generations"],
                       source_model=solve_kwargs["source_model"],
                       edge_open_model=solve_kwargs.get("edge_open_model", "none"))
    if "edge_open_electron_flux" in solve_kwargs:
        edge_kwargs["edge_open_electron_flux"] = solve_kwargs["edge_open_electron_flux"]
    r4 = solve_edge_array_charging(4.0, W=charging_W, mouth=charging_mouth,
                                   n_per_iter=5000, n_iter=120, seed=7, smooth=True, **edge_kwargs)
else:
    r4 = solve_trench_charging(4.0, W=charging_W, mouth=charging_mouth,
                               n_per_iter=8000, n_iter=140, seed=7, smooth=True, **solve_kwargs)
V = r4["V"]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw=dict(width_ratios=[1.15, 1]))

axA.plot(g["ar"], g["model"], "o-", color="#2471c7", lw=2.6, ms=8,
         label=f"petch 2-D charging solver (RMSE {float(g['rmse']):.3f})")
axA.plot(g["ar"], g["hg"], "k*", ms=16, label="Hwang–Giapis 1997 (digitized)", zorder=6)
axA.plot(g["ar"], g["hg"], "k--", lw=1.2, alpha=0.5)
axA.set_xlabel("aspect ratio"); axA.set_ylabel("normalized floor ion flux")
axA.set_ylim(0, 0.8); axA.grid(alpha=0.3); axA.legend(loc="upper right", fontsize=10)
axA.set_title("Floor ion flux vs AR — model vs published data", fontsize=11)
axA.text(1.05, 0.34, f"mechanism config (RMSE {float(g['rmse']):.3f})\n"
         "closure config passes the gate at 0.039\nnothing tuned",
         fontsize=8.5, color="0.35", va="top", ha="left")
axA.annotate("with 300 eV ions the floor stays open\n(0.56 @ AR 4 — the Matsui asymptote)",
             xy=(4.0, 0.22), xytext=(1.6, 0.10), fontsize=9, color="0.3",
             arrowprops=dict(arrowstyle="->", color="0.5"))

# exact solver geometry (aligns the outline to the field)
gm = r4["geom"]; W_, mouth_, nz_, nx_ = gm["W"], gm["mouth"], gm["nz"], gm["nx"]
z0 = max(0, mouth_ - 30)
wall = np.zeros_like(V, bool)
cy = "#4dd0e1"
if charging_geometry == "edge_array":
    edge0, edge1 = gm["edge0"], gm["edge1"]
    trench0, trench1 = gm["trench0"], gm["trench1"]
    neigh0, neigh1 = gm["neigh0"], gm["neigh1"]
    next0 = gm.get("next0", neigh1)
    wall[edge0:edge1, mouth_:] = True
    wall[neigh0:neigh1, mouth_:] = True
    wall[next0:, mouth_:] = True
    outline = [
        ([edge0, edge1], [mouth_, mouth_]), ([neigh0, neigh1], [mouth_, mouth_]),
        ([next0, nx_], [mouth_, mouth_]),
        ([edge0, edge0], [mouth_, nz_ - 1]), ([edge1, edge1], [mouth_, nz_ - 1]),
        ([neigh0, neigh0], [mouth_, nz_ - 1]), ([neigh1, neigh1], [mouth_, nz_ - 1]),
        ([trench0, trench1], [nz_ - 1, nz_ - 1]),
    ]
else:
    pad_ = gm["pad"]
    wall[:pad_, mouth_:] = True
    wall[pad_ + W_:, mouth_:] = True
    outline = [
        ([0, pad_], [mouth_, mouth_]), ([pad_ + W_, nx_], [mouth_, mouth_]),
        ([pad_, pad_], [mouth_, nz_ - 1]), ([pad_ + W_, pad_ + W_], [mouth_, nz_ - 1]),
        ([pad_, pad_ + W_], [nz_ - 1, nz_ - 1]),
    ]
Vshow = np.where(wall, np.nan, V)
Vc = Vshow[:, z0:]
cmap = plt.cm.inferno.copy(); cmap.set_bad("#08040f")
im = axB.imshow(Vc.T, origin="upper", cmap=cmap, aspect="auto", extent=[0, nx_, nz_, z0])
axB.contour(Vc.T, levels=10, colors="w", linewidths=0.4, alpha=0.5,
            extent=[0, nx_, nz_, z0], origin="upper")
for xs, ys in outline:
    axB.plot(xs, ys, color=cy, lw=2)
axB.set_xlim(0, nx_); axB.set_ylim(nz_, z0)
axB.set_title(f"Steady-state potential at AR 4  (floor {r4['V_floor_center']:.0f} V)", fontsize=11)
axB.set_xlabel("x (cells)"); axB.set_ylabel("z (cells, plasma at top)")
plt.colorbar(im, ax=axB, label="V (volts, sheath edge = 0)")

fig.suptitle("Feature charging — petch's 2-D Hwang–Giapis solver vs the published data",
             fontweight="bold", fontsize=13)
plt.tight_layout()
p = os.path.join(HERE, "viz", "charging_hg.png")
plt.savefig(p, dpi=150)
print("saved", p)
