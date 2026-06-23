"""Localize petch's residual deep-AR gap. Across an AR sweep, recompute the flux on each final mesh
and dump the FLOOR-CENTER rate decomposition (chemical k_sigma*thF/4 vs sputter vs ion-enhanced
thF*GY_ie) plus the floor fluxes/coverage. Shows which term collapses with aspect ratio. cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3
from petch.threed import _belen_coverages

DX, W, XE, YE, SUB = 0.04, 0.5, 1.5, 0.3, 6.0
FLG = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                  warm_start_coverage=True, sampling="sobol", neutral_transport="mc", ion_reflection=True)


def run(dur):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W,
               mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 22)), par=p, flags=FLG,
                       n_ion=80000, n_neu=80000, reinit_method="fsm", verbose=False, **GEO)
    verts, faces, centroids, areas = t3.extract_mesh_3d(g['phi'], DX)
    mesh = t3.wp.Mesh(points=t3.wp.array(verts, dtype=t3.wp.vec3, device=t3.DEVICE),
                      indices=t3.wp.array(faces.flatten(), dtype=t3.wp.int32, device=t3.DEVICE))
    m_i, m_F, m_O, cos_i, bare = t3.mc_flux_3d_coupled(mesh, verts, faces, areas, g, p, n_ion=80000,
                                                       n_neu=80000, seed=7, sampling="sobol", flags=FLG)
    thF, thO = _belen_coverages(m_i, m_F, m_O, cos_i, p, FLG)
    sp_n, ie_n, ip_n = p['_ion_yield']; ionF = p['ionFlux']
    chem = p['k_sigma'] * thF / 4.0; sput = ionF * sp_n; ienh = thF * ionF * ie_n
    tot = chem + sput + ienh
    cz = centroids[:, 2]; cx = centroids[:, 0]
    dep = float(g['sub_top'] - cz.min())
    # floor-center = deepest band, near trench center x
    flo = (cz < cz.min() + 2.5 * DX) & (np.abs(cx - XE / 2) < 0.12)
    if flo.sum() < 1:
        flo = (cz < cz.min() + 2.5 * DX)
    s = flo
    return dict(AR=dep / W, dep=dep, m_i=m_i[s].mean(), m_F=m_F[s].mean(), m_O=m_O[s].mean(),
                thF=thF[s].mean(), chem=chem[s].mean(), sput=sput[s].mean(), ienh=ienh[s].mean(),
                tot=tot[s].mean(), n=int(s.sum()))


print("floor-center decomposition vs aspect ratio (petch, faithful ion)\n", flush=True)
print("  AR   dep   nfc |  m_i   m_F   m_O  thF  |  chem  sput  ienh   tot  | tot/tot(AR~3)", flush=True)
rows = [run(d) for d in [0.7, 1.2, 1.7, 2.2, 2.7]]
tot0 = rows[0]['tot']
for r in rows:
    print(f"  {r['AR']:4.1f} {r['dep']:4.2f} {r['n']:4d} | {r['m_i']:.2f}  {r['m_F']:.2f}  {r['m_O']:.2f}  "
          f"{r['thF']:.3f} | {r['chem']:5.1f} {r['sput']:5.2f} {r['ienh']:5.1f} {r['tot']:6.1f} | "
          f"{r['tot']/max(tot0,1e-9):.3f}", flush=True)
print("\n  m_i=ion arrival (w/ reflection), m_F/m_O=neutral F/O flux, thF=F coverage", flush=True)
print("  tot/tot0 is the normalized floor rate ~ ARDE nr (compare to ViennaPS 1.0,0.861,0.731)", flush=True)
