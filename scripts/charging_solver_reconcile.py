"""C1 of the frontier loop: reconcile the two charging solvers and show the floor over-charge is
BRACKETED by the electron launch/source model. Old edge_array solver (sheath launch plane
boundary_um=3.7 + rf_bursts) OVER-delivers floor electrons; the general engine's simple z=1 launch
UNDER-delivers; HG's floorV 33 / floorFlux 0.22 sits between. => the faithful phase-resolved source
(the C2 kinetic engine) is the lever. Saves viz/charging_solver_reconcile.png (+ docs/ copy)."""
import sys, time, os; sys.path.insert(0, 'src')
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from petch.charging2d import _build_edge_array_geometry, solve_edge_array_charging
from petch.charging_general import solve_charging, GAS, INSULATOR, CONDUCTOR

W, MOUTH, AR = 16, 80, 4.0
g = _build_edge_array_geometry(AR, W=W, mouth=MOUTH)
t0, t1 = g['trench0'], g['trench1']; nx = g['nx']; fz = 143
mat = np.where(g['solid'], np.where(g['cond'] > 0, CONDUCTOR, INSULATOR), GAS).astype(np.int64)

NITS = [120, 400, 700]
old_V, old_F, new_V, new_F = [], [], [], []
for nit in NITS:
    ro = solve_edge_array_charging(AR, W=W, mouth=MOUTH, n_per_iter=6000, n_iter=nit, seed=7)
    old_V.append(ro['V_floor_center']); old_F.append(ro['floor_flux'])
    rn = solve_charging(mat, mouth=MOUTH, field_model='laplace', electron_model='trace',
                        electron_open_vf=True, ied_bias=0.25, open_wall_boost=1.4,
                        n_per_iter=6000, n_iter=nit, seed=7)
    d = rn['ntot'] / nx
    new_V.append(rn['Vs'][t0:t1, fz].mean()); new_F.append(rn['ion_counts'][t0:t1, fz].mean() / d)
    print(f"nit={nit}: OLD V={old_V[-1]:.1f} F={old_F[-1]:.3f} | NEW V={new_V[-1]:.1f} F={new_F[-1]:.3f}", flush=True)

HG_V, HG_F = 33.0, 0.22
fig, (axV, axF) = plt.subplots(1, 2, figsize=(12.5, 5.0))
x = np.array(NITS)
axV.plot(x, old_V, "o-", color="#e07b39", lw=2.4, ms=8, label="old edge_array (sheath launch plane)")
axV.plot(x, new_V, "s-", color="#2471c7", lw=2.4, ms=8, label="general engine (z=1 launch)")
axV.axhline(HG_V, color="k", ls="--", lw=1.6, label="Hwang–Giapis target (33 V)")
axV.fill_between(x, old_V, new_V, color="0.85", alpha=0.5, zorder=0)
axV.set_xlabel("charging iterations (convergence)"); axV.set_ylabel("floor-center potential (V)")
axV.set_title("Floor over-charge is BRACKETED by the electron source", fontsize=11)
axV.legend(fontsize=9); axV.grid(alpha=0.3)

axF.plot(x, old_F, "o-", color="#e07b39", lw=2.4, ms=8, label="old edge_array (over-delivers e⁻)")
axF.plot(x, new_F, "s-", color="#2471c7", lw=2.4, ms=8, label="general engine (under-delivers e⁻)")
axF.axhline(HG_F, color="k", ls="--", lw=1.6, label="Hwang–Giapis target (0.22)")
axF.fill_between(x, old_F, new_F, color="0.85", alpha=0.5, zorder=0)
axF.set_xlabel("charging iterations (convergence)"); axF.set_ylabel("normalized floor ion flux")
axF.set_title("HG's floor sits between the two launch models", fontsize=11)
axF.legend(fontsize=9); axF.grid(alpha=0.3)

fig.suptitle("C1 reconcile — the electron launch/source is the lever for the floor over-charge "
             "(→ faithful phase-resolved source = C2 kinetic engine)", fontsize=11.5, y=1.00)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/charging_solver_reconcile.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/charging_solver_reconcile.png", dpi=130, bbox_inches="tight")
print("saved viz/charging_solver_reconcile.png", flush=True)
print("DONE", flush=True)
