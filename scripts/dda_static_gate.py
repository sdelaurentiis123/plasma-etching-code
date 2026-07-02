#!/usr/bin/env python3
"""Deterministic-transport static-ARDE gate vs the MEASURED ViennaPS reference.
PRIMARY (passing): neutral_transport='radiosity' + radiosity_solver='gmres' -> RMSE 0.043.
SECONDARY (reported): the grid-march DDA after its transport fixes -- too steep in the
passivated-wall regime (documented structural limit; see RECONCILIATION). Reference: ViennaPS
static reference (2026-07-02 box run, pre-carved trenches, reflecting ions, GPU_TRIANGLE):
nr = 0.911/0.820/0.728/0.626/0.534 at AR 2/4/6/8/10, W=0.5um trench. Gate: RMSE <= 0.05."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.environ.setdefault("PETCH_DEVICE", "cpu")
import numpy as np, skfmm, warp as wp
import petch
from petch import threed as t3

W, XE, YE, SUB, DX = 0.5, 2.0, 1.0, 8.0, 0.05
VPS_AR = np.array([2, 4, 6, 8, 10], float)
VPS_NR = np.array([0.911, 0.820, 0.728, 0.626, 0.534])

def clean_trench_phi(D):
    Lz = SUB + 0.4
    nx, ny, nz = round(XE/DX), round(YE/DX), round(Lz/DX)
    xs, ys, zs = np.arange(nx)*DX, np.arange(ny)*DX, np.arange(nz)*DX
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij'); r = np.abs(X - XE/2)
    f_sub = Z - SUB; f_box = np.maximum(r - W/2, (SUB - D) - Z)
    phi = skfmm.distance(-np.maximum(f_sub, -f_box), dx=DX)
    return dict(xs=xs, ys=ys, zs=zs, dx=DX, phi=phi, Lx=XE, Ly=YE, Lz=Lz, sub_top=SUB,
                trench_width=W, hole=False, mask=np.zeros_like(phi, bool))

def nr(AR, transport='radiosity', n_dir=64, n_re=12):
    geo = clean_trench_phi(AR * W)
    verts, faces, cen, areas = t3.extract_mesh_3d(geo['phi'], DX)
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=t3.DEVICE),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=t3.DEVICE))
    par = dict(petch.PAR); par['dda_n_dir'] = n_dir; par['dda_n_reemit'] = n_re; par['periodic_y'] = 1
    par['radiosity_solver'] = 'gmres'
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     sampling="sobol", ion_reflection=True, neutral_transport=transport)
    if transport == 'radiosity':
        m_i, m_F, m_O, cos_i = t3.mc_flux_3d_radiosity(mesh, verts, faces, cen, areas, geo, par,
                                                       n_ion=40000, n_ff=128, seed=0, flags=fl, n_fp=4)
    else:
        m_i, m_F, m_O, cos_i = t3.mc_flux_3d_dda(mesh, verts, faces, cen, areas, geo, par,
                                                 n_ion=40000, seed=0, flags=fl, n_fp=3)
    V = t3.surface_rate(m_i, m_F, m_O, cos_i, np.zeros(len(faces), bool), par, flags=fl)
    cx = cen[:, 0] - XE/2; cz = cen[:, 2]
    floor = (np.abs(cx) < 0.15) & (cz < (SUB - AR*W) + 3*DX)
    field = cz > SUB - 0.5*DX
    return float(np.nanmean(V[floor]) / np.nanmean(V[field]))

if __name__ == "__main__":
    for transport in ("radiosity", "dda"):
        vals = []
        for i, ar in enumerate(VPS_AR):
            t0 = time.time(); v = nr(ar, transport=transport); vals.append(v)
            print(f"  [{transport}] AR{ar:4.0f}: petch={v:.3f}  ViennaPS={VPS_NR[i]:.3f}  delta={v-VPS_NR[i]:+.3f}  ({time.time()-t0:.0f}s)", flush=True)
        vals = np.array(vals)
        rmse = float(np.sqrt(np.mean((vals - VPS_NR) ** 2)))
        gate = 'PASS' if rmse <= 0.05 else ('fail (documented passivated-regime limit)' if transport == 'dda' else 'fail')
        print(f"  {transport}-vs-ViennaPS RMSE = {rmse:.3f}  [{gate}] (gate 0.05)", flush=True)
