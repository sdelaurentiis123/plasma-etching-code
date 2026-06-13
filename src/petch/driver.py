"""Main etch driver. Ported from feature_etch.run_etch.

Threads a `Flags` object through chemistry (and, later, transport) so model variants can be
A/B'd without forking the loop. Default flags reproduce the original PoC behavior. Exposes
`n_part_ion`/`n_part_neu` so the convergence harness can sweep ray count.
"""
import time
import numpy as np

from .params import PAR, DEFAULT_FLAGS
from .geometry import (make_trench, extract_surface, orient_normals,
                       seg_in_mask, profile_bottom)
from .transport import mc_flux
from .chemistry import surface_rate
from .levelset import advect, extend_velocity, reinit


def run_etch(W=20.0, H=20.0, dx=0.25, trench_width=8.0, mask_thickness=2.0,
             sub_top=15.0, t_end=3.0, n_steps=60, par=None, flags=None,
             n_part_ion=20000, n_part_neu=20000, seed=0, snapshots=None, verbose=True):
    if par is None:
        par = PAR
    if flags is None:
        flags = DEFAULT_FLAGS
    X, Y, xs, ys, phi, mask, nx, ny = make_trench(W, H, dx, trench_width, mask_thickness, sub_top)
    mask_phi = phi.copy()                  # to re-stamp mask each step
    # belen chemistry uses ViennaPS sticking; keep transport re-emission consistent with it
    mc_par = par
    if getattr(flags, "chemistry", "langmuir") == "belen":
        mc_par = dict(par); mc_par['s_F'] = par['betaE']; mc_par['s_O'] = par['betaO']
    y_src = H - dx
    dt = t_end / n_steps
    band = 5 * dx
    timings = dict(raytrace=0.0, total=0.0, chem=0.0, advect=0.0)
    snaps = {}
    t_total0 = time.time()
    for step in range(n_steps):
        segs, mid, nrm, L = extract_surface(phi, xs, ys, dx)
        if len(segs) == 0:
            break
        nrm = orient_normals(mid, nrm, phi, xs, ys, dx)
        is_mask = seg_in_mask(mid, mask, xs, ys, dx)
        t0 = time.time()
        m_i, m_F, m_O, cos_i = mc_flux(segs, mid, nrm, is_mask, L, y_src, W, mc_par,
                                       n_part_ion=n_part_ion, n_part_neu=n_part_neu, seed=seed)
        timings['raytrace'] += time.time() - t0
        t0 = time.time()
        V = surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=flags)
        timings['chem'] += time.time() - t0
        F = extend_velocity(V, mid, phi, xs, ys, dx, band)
        t0 = time.time()
        Vmax = max(V.max(), 1e-6)
        nsub = int(np.ceil(Vmax * dt / (0.4 * dx)))
        nsub = max(1, min(nsub, 40))
        for _ in range(nsub):
            phi = advect(phi, F, dx, dt / nsub)
            phi[mask] = mask_phi[mask]
        phi = reinit(phi, dx)                   # restore signed distance
        timings['advect'] += time.time() - t0
        if snapshots and step in snapshots:
            snaps[step] = extract_surface(phi, xs, ys, dx)[0]
        if verbose and step % 10 == 0:
            depth = sub_top - profile_bottom(phi, xs, ys, dx, W)
            print(f"  step {step:3d}/{n_steps}  etch depth ~ {depth:5.2f} um  Vmax {Vmax:.3f}")
    timings['total'] = time.time() - t_total0
    final_segs = extract_surface(phi, xs, ys, dx)[0]
    return dict(phi=phi, xs=xs, ys=ys, dx=dx, segs=final_segs, snaps=snaps,
                timings=timings, sub_top=sub_top, mask=mask, X=X, Y=Y)
