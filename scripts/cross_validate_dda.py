"""Numerics divergence stress-test: petch (Monte-Carlo, Warp) vs Craig Xu Chen's plasma_sim
(deterministic DDA discrete-ordinates, Metal) on MATCHED Belen SF6/O2 chemistry, so ONLY the
neutral-transport method differs. We compare the normalized ARDE rate n_r vs aspect ratio for
hole + trench, across petch's {mc, knudsen, radiosity} and Craig's {dda(gmres), knudsen}, plus
a dx-convergence check, and report WHERE the engines diverge.

petch runs on CPU (Warp-CPU); plasma_sim runs on the Apple Metal GPU (mlx). ARDE is normalized
(n_r=1 at AR_REF) so the per-engine absolute rate_scale is irrelevant — only the rolloff shape is
compared. Writes cross_validate_dda.npz + viz/cross_validate_dda.png and prints a divergence table.
"""
import os, sys, time, json
import numpy as np
os.environ.setdefault("PETCH_DEVICE", "cpu")
sys.path.insert(0, "/Users/stanislavdelaurentiis/chip-etch")          # plasma_sim
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import petch
from petch import threed as t3
import plasma_sim as ps
from plasma_sim.solver3d import Solver3D

W = 0.5                       # feature width / diameter (um)
MASK = 0.3
AR_REF = 2.0                  # normalize n_r=1 here
EXTENT = 1.2

# ---- matched Belen SF6/O2 chemistry (identical in both engines; only transport differs) ----
def petch_par(rate_scale):
    p = dict(petch.PAR)
    p.update(ied_mode='gauss', rate_scale=rate_scale, radiosity_solver='gmres',
             ionFlux=12.0, Fflux=1800.0, Oflux=100.0, k_sigma=300.0, beta_sigma=0.04, B_sp=9.3,
             betaE=0.7, betaO=1.0, A_ie=7.0, A_sp=0.0337, A_p=3.0, Eth_ie=15.0, Eth_sp=20.0,
             Eth_p=10.0, rho=5.02, Emean=100.0, Esig=10.0)
    return p

def craig_model(rate_scale):
    return ps.SF6O2.default(ion_flux=12, f_flux=1800, o_flux=100, cal_f=1.0, ied_mode='gauss',
                            mean_energy=100, sigma_energy=10, k_sigma=300, beta_sigma=0.04, b_sp=9.3,
                            beta_f=0.7, beta_o=1.0, a_ie=7, a_sp=0.0337, a_p=3,
                            eth_ie=15, eth_sp=20, eth_p=10, rho=5.02, rate_scale=rate_scale)

# ---- engine runners: return (times, center_depths) ----
def run_petch(shape, transport, dx, ion_refl, rate_scale=0.5, t_end=9.0, n_steps=60):
    par = petch_par(rate_scale)
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     ion_reflection=ion_refl, sampling="sobol", neutral_transport=transport)
    sub = 8.0
    t0 = time.time()
    g = t3.run_etch_3d(Lx=EXTENT, Ly=EXTENT, Lz=sub + MASK + 0.5, dx=dx, trench_width=W, mask_th=MASK,
                       sub_top=sub, t_end=t_end, n_steps=n_steps, hole=(shape == "hole"),
                       par=par, flags=fl, n_ion=8000, n_neu=8000, reinit_method="fsm",
                       verbose=False, record_depth_every=1)
    dt = time.time() - t0
    dh = g['depth_history']
    step_dt = t_end / n_steps
    times = np.array([s * step_dt for s, _ in dh])
    depth = np.array([d for _, d in dh])
    return times, depth, dt

def run_craig(shape, solver_name, dx, rate_scale=0.022, duration=5.0, nframes=16):
    mdl = craig_model(rate_scale)
    if shape == "hole":
        dom = ps.Domain.via(extent=EXTENT, dx=dx, diameter=W, mask=MASK, depth=7.5, headroom=0.5)
    else:
        dom = ps.Domain.trench(extent=EXTENT, dx=dx, width=W, mask=MASK, depth=7.5, headroom=0.5)
    cfg = ps.Process(dom, mdl, duration=duration, use_gpu=True).config()
    cfg.neutral_solver = solver_name
    solver = Solver3D(cfg)
    rr = np.sqrt(solver.X**2 + solver.Y**2) if shape == "hole" else np.abs(solver.X)
    cen = rr < max(cfg.via_radius * 0.5, cfg.dx)
    rec = []
    def on_frame(t, phi):
        m = (phi > 0) & cen
        rec.append((t, max(0.0, -float(solver.Z[m].min())) if m.any() else 0.0))
    t0 = time.time()
    solver.run(duration, list(np.linspace(0.0, duration, nframes)), on_frame)
    dt = time.time() - t0
    rec = np.array(rec)
    return rec[:, 0], rec[:, 1], dt

# ---- normalized ARDE rate n_r vs AR ----
def nr_curve(times, depth):
    """floor rate = d(depth)/dt; AR at interval midpoint; normalize n_r=1 at AR_REF."""
    depth = np.maximum.accumulate(depth)                 # enforce monotone (guard MC jitter)
    rate = np.diff(depth) / np.maximum(np.diff(times), 1e-9)
    ar = 0.5 * (depth[1:] + depth[:-1]) / W
    keep = ar > 0.3
    ar, rate = ar[keep], rate[keep]
    if ar.size < 3 or ar.max() < AR_REF:
        return None
    order = np.argsort(ar); ar, rate = ar[order], rate[order]
    r_ref = np.interp(AR_REF, ar, rate)
    if r_ref <= 0:
        return None
    return ar, rate / r_ref

CASES = [
    # (engine, label, shape, transport/solver, dx, ion_refl)
    ("petch", "petch mc",        "hole", "mc",        0.05, False),
    ("petch", "petch knudsen",   "hole", "knudsen",   0.05, False),
    ("petch", "petch radiosity", "hole", "radiosity", 0.05, False),
    ("petch", "petch mc+ionrefl","hole", "mc",        0.05, True),
    ("craig", "Craig DDA",       "hole", "gmres",     0.05, None),
    ("craig", "Craig knudsen",   "hole", "knudsen",   0.05, None),
    ("petch", "petch mc",        "trench", "mc",        0.05, False),
    ("petch", "petch knudsen",   "trench", "knudsen",   0.05, False),
    ("petch", "petch radiosity", "trench", "radiosity", 0.05, False),
    ("petch", "petch mc+ionrefl","trench", "mc",        0.05, True),
    ("craig", "Craig DDA",       "trench", "gmres",     0.05, None),
    ("craig", "Craig knudsen",   "trench", "knudsen",   0.05, None),
]
# dx-convergence: hole petch-mc + Craig-DDA at coarse/fine dx (0.05 already above)
DX_CONV = [
    ("petch", "petch mc",  "hole", "mc",    0.08, False),
    ("craig", "Craig DDA", "hole", "gmres", 0.08, None),
    ("petch", "petch mc",  "hole", "mc",    0.04, False),
    ("craig", "Craig DDA", "hole", "gmres", 0.04, None),
]

def run_case(engine, label, shape, mode, dx, ion_refl):
    if engine == "petch":
        return run_petch(shape, mode, dx, ion_refl)
    return run_craig(shape, mode, dx)

results = {}      # key -> dict(times, depth, wall, ar, nr)
def do(cases, tag):
    for (engine, label, shape, mode, dx, ion_refl) in cases:
        key = f"{label}|{shape}|dx{dx}"
        try:
            times, depth, wall = run_case(engine, label, shape, mode, dx, ion_refl)
            cur = nr_curve(times, depth)
            ar, nr = (cur if cur is not None else (np.array([]), np.array([])))
            results[key] = dict(engine=engine, label=label, shape=shape, dx=dx,
                                times=times, depth=depth, wall=wall, ar=ar, nr=nr)
            print(f"[{tag}] {key:34s} wall={wall:6.1f}s  AR_end={depth[-1]/W:4.1f}  pts={len(ar)}", flush=True)
        except Exception as e:
            print(f"[{tag}] {key:34s} FAILED: {type(e).__name__}: {e}", flush=True)

print("=== main battery (dx=0.05) ===", flush=True)
do(CASES, "main")
print("=== dx-convergence (hole) ===", flush=True)
do(DX_CONV, "dxconv")

# ---- common AR grid + divergence table ----
def grid_for(shape):
    petch_mc = results.get(f"petch mc|{shape}|dx0.05")
    craig    = results.get(f"Craig DDA|{shape}|dx0.05")
    hi = 6.5
    for r in (petch_mc, craig):
        if r is not None and len(r['ar']):
            hi = min(hi, float(r['ar'].max()))
    return np.arange(AR_REF, hi + 1e-6, 0.5)

def resample(r, grid):
    if r is None or len(r['ar']) < 2:
        return np.full(len(grid), np.nan)
    return np.interp(grid, r['ar'], r['nr'], left=np.nan, right=np.nan)

print("\n================= DIVERGENCE TABLE (normalized n_r vs AR) =================", flush=True)
table = {}
for shape in ("hole", "trench"):
    grid = grid_for(shape)
    table[shape] = dict(grid=grid.tolist())
    labels = ["petch mc", "petch knudsen", "petch radiosity", "petch mc+ionrefl", "Craig DDA", "Craig knudsen"]
    cols = {lab: resample(results.get(f"{lab}|{shape}|dx0.05"), grid) for lab in labels}
    for lab in labels:
        table[shape][lab] = np.where(np.isfinite(cols[lab]), cols[lab], None).tolist()
    delta = np.abs(cols["petch mc"] - cols["Craig DDA"])
    table[shape]["delta_petchmc_vs_dda"] = np.where(np.isfinite(delta), delta, None).tolist()
    print(f"\n--- {shape} ---  (n_r=1 at AR={AR_REF})", flush=True)
    print("  AR  " + "".join(f"{lab:>17s}" for lab in labels) + "   |dMC-DDA|", flush=True)
    for i, ar in enumerate(grid):
        row = f" {ar:4.1f} " + "".join(f"{cols[lab][i]:17.3f}" if np.isfinite(cols[lab][i]) else f"{'-':>17s}" for lab in labels)
        row += f"   {delta[i]:8.3f}" if np.isfinite(delta[i]) else f"   {'-':>8s}"
        print(row, flush=True)
    if np.isfinite(delta).any():
        j = int(np.nanargmax(delta))
        print(f"  >> MAX divergence petch-MC vs Craig-DDA: {delta[j]:.3f} at AR={grid[j]:.1f}", flush=True)
        table[shape]["max_div"] = dict(ar=float(grid[j]), delta=float(delta[j]))

# dx-convergence summary: divergence at a fixed AR vs dx
print("\n================= dx-CONVERGENCE (hole, petch-mc vs Craig-DDA) =================", flush=True)
for dx in (0.08, 0.05, 0.04):
    rp = results.get(f"petch mc|hole|dx{dx}"); rc = results.get(f"Craig DDA|hole|dx{dx}")
    if rp is None or rc is None or not len(rp['ar']) or not len(rc['ar']):
        print(f"  dx={dx}: missing", flush=True); continue
    hi = min(rp['ar'].max(), rc['ar'].max())
    g = np.arange(AR_REF, hi + 1e-6, 0.5)
    d = np.abs(np.interp(g, rp['ar'], rp['nr']) - np.interp(g, rc['ar'], rc['nr']))
    print(f"  dx={dx}: mean|dMC-DDA|={np.nanmean(d):.3f}  max={np.nanmax(d):.3f} at AR={g[np.nanargmax(d)]:.1f}  (AR_max={hi:.1f})", flush=True)

# ---- save + plot ----
np.savez(os.path.join(HERE, "cross_validate_dda.npz"),
         results={k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv)
                      for kk, vv in v.items()} for k, v in results.items()},
         table=table, W=W, AR_REF=AR_REF, allow_pickle=True)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
STYLE = {"petch mc": ("#2471c7", "-"), "petch knudsen": ("#27ae60", "--"),
         "petch radiosity": ("#8e44ad", ":"), "petch mc+ionrefl": ("#16a085", "-."),
         "Craig DDA": ("#c0392b", "-"), "Craig knudsen": ("#e67e22", "--")}
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
for ax, shape in zip(axes, ("hole", "trench")):
    for lab, (col, ls) in STYLE.items():
        r = results.get(f"{lab}|{shape}|dx0.05")
        if r is not None and len(r['ar']):
            ax.plot(r['ar'], r['nr'], ls, color=col, lw=2.2, label=lab)
    ax.set_title(f"{shape}  (W = {W} µm, dx = 0.05)"); ax.set_xlabel("aspect ratio  depth / W")
    ax.grid(alpha=0.3); ax.set_xlim(AR_REF, None); ax.set_ylim(0, 1.15)
    md = table[shape].get("max_div")
    if md:
        ax.axvline(md["ar"], color="0.5", ls=":", lw=1)
        ax.annotate(f"max |petch-MC − DDA|\n= {md['delta']:.2f} @ AR {md['ar']:.1f}",
                    xy=(md["ar"], 0.5), xytext=(AR_REF + 0.2, 0.12), fontsize=9,
                    arrowprops=dict(arrowstyle="->", color="0.4"))
axes[0].set_ylabel("normalized etch rate  $n_r$"); axes[0].legend(fontsize=9, loc="upper right")
fig.suptitle("petch (MC) vs Craig's plasma_sim (DDA) — transport-only numerics cross-validation",
             fontweight="bold")
plt.tight_layout()
os.makedirs(os.path.join(HERE, "viz"), exist_ok=True)
out = os.path.join(HERE, "viz", "cross_validate_dda.png")
plt.savefig(out, dpi=150); print("\nsaved", out, flush=True)
print("saved", os.path.join(HERE, "cross_validate_dda.npz"), flush=True)
print("DONE", flush=True)
