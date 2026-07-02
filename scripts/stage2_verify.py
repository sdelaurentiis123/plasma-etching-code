"""STAGE 2 — verifiable gates, all with the FROZEN calibration (wls=1.4, faithful ion).

TEST 1 (width-sweep PREDICTION): frozen knob, W = 1/2/4 um trenches. Gates:
  (a) nr(AR) collapses across widths (de Boer's empirical finding) -- max cross-width spread;
  (b) each width's RMSE vs the wafer curve <= 0.05 WITHOUT re-calibration.
TEST 2 (Gomez absolute rate): the single rate_scale that gives the open-field 1.3 um/min,
  and whether the SAME constant holds across W (field rate should be W-independent).
TEST 3 (evolving-etch consistency): evolve the W=2 trench (knudsen+ion), depth(t) -> nr(AR),
  gate: RMSE vs the static curve <= 0.05 (the pass is not a static-harness artifact).
"""
import os, time; os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np, warp as wp, skfmm
import petch
from petch import threed as t3

WAFER_AR = np.array([0.0, 10.0, 20.0, 40.0]); WAFER_NR = np.array([1.0, 0.43, 0.29, 0.20])
WLS = 1.4                          # FROZEN (calibrated at W=2; never touched below)

def static_rate(W, D, dx, extra_par=None):
    XE, YE = max(5*W, 6.0), max(2.5*W, 2.0)
    sub_top = D + 2*dx + 0.5; Lz = sub_top + 0.6
    nx, ny, nz = round(XE/dx), round(YE/dx), round(Lz/dx)
    xs, ys, zs = np.arange(nx)*dx, np.arange(ny)*dx, np.arange(nz)*dx
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij'); r = np.abs(X - XE/2)
    phi = skfmm.distance(-np.maximum(Z - sub_top, -(np.maximum(r - W/2, (sub_top - D) - Z))), dx=dx)
    geo = dict(xs=xs, ys=ys, zs=zs, dx=dx, phi=phi, Lx=XE, Ly=YE, Lz=Lz, sub_top=sub_top,
               trench_width=W, hole=False, mask=np.zeros_like(phi, bool))
    verts, faces, cen, areas = t3.extract_mesh_3d(phi, dx)
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=t3.DEVICE),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=t3.DEVICE))
    par = dict(petch.PAR); par['knudsen_wall_loss_scale'] = WLS; par['periodic_y'] = 1
    if extra_par: par.update(extra_par)
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     sampling="sobol", ion_reflection=True, neutral_transport='knudsen')
    m_i, m_F, m_O, cos_i = t3.mc_flux_3d_knudsen(mesh, verts, faces, cen, areas, geo, par,
                                                 n_ion=100000, seed=0, flags=fl, n_fp=3)
    V = t3.surface_rate(m_i, m_F, m_O, cos_i, np.zeros(len(faces), bool), par, flags=fl)
    cx = cen[:, 0] - XE/2; cz = cen[:, 2]
    floor = (np.abs(cx) < 0.15*W) & (cz < (sub_top - D) + 3*dx)
    field = cz > sub_top - 0.5*dx
    return float(np.nanmean(V[floor])), float(np.nanmean(V[field]))

print("=== TEST 1: width-sweep PREDICTION (frozen wls=1.4) ===", flush=True)
curves = {}
for W, dx in [(1.0, 0.125), (2.0, 0.25), (4.0, 0.5)]:
    nr = [1.0]
    for AR in (10, 20, 40):
        f, f0 = static_rate(W, AR*W, dx)
        nr.append(f / f0)
    nr = np.array(nr); curves[W] = nr
    rmse = float(np.sqrt(np.mean((nr - WAFER_NR)**2)))
    print(f"  W={W:3.1f}um: nr={np.round(nr,3)}  RMSE_vs_wafer={rmse:.3f}  [{'PASS' if rmse<=0.05 else 'fail'}]", flush=True)
allc = np.array(list(curves.values()))
spread = float(np.max(allc.max(0) - allc.min(0)))
print(f"  cross-width collapse: max spread = {spread:.3f}  [{'PASS (collapses)' if spread <= 0.08 else 'fail'}]", flush=True)

print("=== TEST 2: Gomez open-field absolute rate (1.3 um/min +-10%) ===", flush=True)
for W, dx in [(2.0, 0.25), (0.5, 0.05)]:
    _, f0 = static_rate(W, 0.5*W, dx)        # shallow feature; field faces give open rate
    # f0 is in petch velocity units * rate_scale(default). rate_scale needed for 1.3 um/min:
    rs_needed = 1.3 / max(60.0 * f0 / petch.PAR['rate_scale'], 1e-30)
    print(f"  W={W}: field V={f0:.4g} um/s (rate_scale={petch.PAR['rate_scale']})  ->  rate_scale for 1.3um/min = {rs_needed:.4g}", flush=True)
print("  gate: the two rate_scale values should agree (one global constant)", flush=True)

print("=== TEST 3: evolving-etch consistency (W=2, knudsen+ion) ===", flush=True)
par = dict(petch.PAR); par['knudsen_wall_loss_scale'] = WLS; par['periodic_y'] = 1
par['rate_scale'] = 0.35
fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                 sampling="sobol", ion_reflection=True, neutral_transport='knudsen',
                 warm_start_coverage=True)
t0 = time.time()
g = t3.run_etch_3d(Lx=10.0, Ly=5.0, Lz=46.0, dx=0.25, trench_width=2.0, mask_th=2.0, sub_top=42.0,
                   t_end=10.0, n_steps=120, hole=False, par=par, flags=fl, n_ion=40000, n_neu=40000,
                   reinit_method="fsm", verbose=False, record_depth_every=2)
h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
tm = st / 120 * 10.0
dd = np.maximum.accumulate(dd)
rate = np.gradient(dd, tm); ar = dd / 2.0
r0 = rate[ar < 2].max() if (ar < 2).any() else rate.max()
nr_ev = np.clip(rate / max(r0, 1e-9), 0, 1.5)
static_nr = np.array([1.0] + [static_rate(2.0, AR*2.0, 0.25)[0] / static_rate(2.0, AR*2.0, 0.25)[1]
                              for AR in ()])  # static already measured in TEST 1 (W=2 row)
st_nr = curves[2.0]
ev_at = np.interp([10.0, 20.0], ar, nr_ev, right=np.nan)
comp = np.array([ev_at[0], ev_at[1]]); ref = np.array([st_nr[1], st_nr[2]])
ok = np.isfinite(comp)
rmse_es = float(np.sqrt(np.nanmean((comp[ok] - ref[ok])**2))) if ok.any() else float('nan')
print(f"  evolving reached AR {ar.max():.1f} in {time.time()-t0:.0f}s; nr@AR10/20 evolving={np.round(comp,3)} vs static={np.round(ref,3)}", flush=True)
print(f"  evolving-vs-static RMSE = {rmse_es:.3f}  [{'PASS' if rmse_es <= 0.05 else 'fail'}]", flush=True)
print("STAGE2 DONE", flush=True)
