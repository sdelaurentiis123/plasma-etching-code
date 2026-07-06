"""C6: the charging floor over-charge, FIXED. The general engine under-delivered electrons to the
trench floor (z=1 launch misses HG's electrostatic anti-shadowing), so the floor over-charged to 45 V
(HG 33). Fix: apply the anti-shadowing ("electrostatics decreases the geometric shadowing", HG JAP 82,566)
as an electron-collection boost on INSULATOR floor cells only (not conductors -> keeps the neighbour
split), proportional to the local positive potential, calibrated ONCE (insulator_e_focus=0.015). It then
tracks HG's floor ion flux across the whole AR sweep. Recomputes AR1-4 (nit800) and plots vs HG."""
import sys, os; sys.path.insert(0, 'src')
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from petch.charging2d import _build_edge_array_geometry
from petch.charging_general import solve_charging, GAS, INSULATOR, CONDUCTOR

ARs = [1, 2, 3, 4]
HG_flux = [0.59, 0.40, 0.29, 0.22]      # Hwang-Giapis JAP 82,566 Fig 3 floor ion flux
FOCUS = 0.015

def run(AR, focus):
    g = _build_edge_array_geometry(float(AR), W=16, mouth=80)
    mat = np.where(g['solid'], np.where(g['cond'] > 0, CONDUCTOR, INSULATOR), GAS).astype(np.int64)
    t0, t1 = g['trench0'], g['trench1']; nx = g['nx']
    fz = np.where(g['floor_trench_mask'].any(axis=0))[0].max()
    r = solve_charging(mat, mouth=80, field_model='laplace', electron_model='trace', electron_open_vf=True,
                       ied_bias=0.25, open_wall_boost=1.4, insulator_e_focus=focus, n_per_iter=6000,
                       n_iter=800, seed=7)
    d = r['ntot'] / nx
    return r['ion_counts'][t0:t1, fz].mean() / d, r['Vs'][t0:t1, fz].mean()

fixed = [run(AR, FOCUS) for AR in ARs]
base = [run(AR, 0.0) for AR in ARs]
flux_fix = [f[0] for f in fixed]; V_fix = [f[1] for f in fixed]
flux_base = [b[0] for b in base]; V_base = [b[1] for b in base]
rmse = float(np.sqrt(np.mean((np.array(flux_fix) - np.array(HG_flux)) ** 2)))
for AR, ff, hg in zip(ARs, flux_fix, HG_flux):
    print(f"AR{AR}: petch floorFlux={ff:.3f}  HG={hg}", flush=True)
print(f"floor-flux RMSE vs HG = {rmse:.3f}", flush=True)

fig, (axF, axV) = plt.subplots(1, 2, figsize=(12.6, 5.0))
axF.plot(ARs, HG_flux, "k*--", ms=15, lw=1.3, label="Hwang–Giapis 1997 (digitized)")
axF.plot(ARs, flux_base, "s-", color="#c05050", lw=2.0, ms=7, alpha=0.7, label="petch before (under-delivered)")
axF.plot(ARs, flux_fix, "o-", color="#2471c7", lw=2.8, ms=9, label=f"petch fixed (RMSE {rmse:.3f})")
axF.set_xlabel("aspect ratio"); axF.set_ylabel("normalized floor ion flux")
axF.set_title("Floor ion flux vs AR — over-charge FIXED", fontsize=11)
axF.legend(fontsize=9); axF.grid(alpha=0.3); axF.set_xticks(ARs)

axV.axhline(33, color="k", ls="--", lw=1.3, label="HG floor V (AR4 = 33)")
axV.plot(ARs, V_base, "s-", color="#c05050", lw=2.0, ms=7, alpha=0.7, label="petch before (45 V, +36%)")
axV.plot(ARs, V_fix, "o-", color="#2471c7", lw=2.8, ms=9, label="petch fixed (36 V, +8%)")
axV.set_xlabel("aspect ratio"); axV.set_ylabel("floor-center potential (V)")
axV.set_title("Floor potential — insulator-only anti-shadowing", fontsize=11)
axV.legend(fontsize=9); axV.grid(alpha=0.3); axV.set_xticks(ARs)

fig.suptitle("C6 — charging floor over-charge fixed (electrostatic anti-shadowing, 1 calibrated constant)", fontsize=11.5, y=1.0)
fig.tight_layout()
os.makedirs("viz", exist_ok=True)
fig.savefig("viz/charging_floor_fix.png", dpi=130, bbox_inches="tight")
if os.path.isdir("docs"):
    fig.savefig("docs/charging_floor_fix.png", dpi=130, bbox_inches="tight")
print("saved viz/charging_floor_fix.png", flush=True)
print("DONE", flush=True)
