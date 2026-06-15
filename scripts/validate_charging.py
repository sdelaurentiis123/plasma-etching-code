#!/usr/bin/env python3
"""Calibrate + validate surface charging vs Hwang-Giapis 1997 (floor ion current vs AR, ~0.5 at AR4).
Electrons reach the wafer with a Gaussian angular spread e_ang_sigma (sheath collimation). Etch holes
to increasing depth (once, charging off), then for each e_ang_sigma compute the floor charging factor
vs AR and compare to HG. Pick the e_ang_sigma whose rolloff matches. Run on a box: PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import warp as wp
import petch
from petch import threed as t3

DEV = t3.DEVICE
DX, W, N = 0.25, 4.0, 50000
GEO = dict(Lx=12, Ly=12, Lz=40, mask_th=2, sub_top=34, hole=True)
HG = {0: 1.0, 2: 0.72, 4: 0.50, 6: 0.42}     # Hwang-Giapis floor ion current vs AR


def mesh_of(phi):
    v, f, c, a = t3.extract_mesh_3d(phi, DX)
    A = np.maximum(a, 0.3 * np.median(a))
    mesh = wp.Mesh(points=wp.array(v, dtype=wp.vec3, device=DEV), indices=wp.array(f.flatten(), dtype=wp.int32, device=DEV))
    return mesh, c, A


def trace(mesh, F, A, kind, sig, sd):
    rng = np.random.default_rng(sd)
    o, d = t3._source3d(kind, N, GEO['Lx'], GEO['Ly'], GEO['Lz']-DX, sig, 'sobol', rng, sd)
    if kind == 'ion':
        fi = wp.zeros(F, dtype=float, device=DEV); ai = wp.zeros(F, dtype=float, device=DEV)
        wp.launch(t3._trace3d, dim=N, device=DEV, inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEV),
                  wp.array(d, dtype=wp.vec3, device=DEV), 1.0, 0, sd, 0, 0.34, 0.8, fi, ai])
        fl = fi.numpy()
    else:
        fl = wp.zeros(F, dtype=float, device=DEV)
        wp.launch(t3._trace3d_cov_rr, dim=N, device=DEV, inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEV),
                  wp.array(d, dtype=wp.vec3, device=DEV), wp.array(np.ones(F, np.float32), dtype=float, device=DEV), 1.0, sd, fl])
        fl = fl.numpy()
    return np.clip((fl / A) / (N / (GEO['Lx']*GEO['Ly'])), 0.0, 8.0)


# etch holes to increasing depth ONCE (charging off -> geometry only)
phis = []
for ns in [6, 14, 22, 32]:
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=ns, par=dict(petch.PAR, n_fp=1),
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, t_end=ns*0.06, **GEO)
    phis.append(g['phi'])

print(f"device={DEV}  calibrating electron angular spread to Hwang-Giapis\n", flush=True)
print(f"  HG floor current vs AR: {HG}\n")
for e_sig in [0.3, 0.5, 0.8, 1.1]:
    row = []
    for phi in phis:
        mesh, c, A = mesh_of(phi); F = len(c)
        m_i = trace(mesh, F, A, 'ion', petch.PAR['ion_ang_sigma'], 11)
        m_e = trace(mesh, F, A, 'ion', e_sig, 13)                  # electrons = ion-source w/ e_sig spread
        ref_i = np.percentile(m_i[m_i > 1e-6], 90); ref_e = np.percentile(m_e[m_e > 1e-6], 90)
        sh_i = np.clip(m_i/ref_i, 0, 1); sh_e = np.clip(m_e/ref_e, 0, 1)
        ratio = np.clip(sh_e/np.maximum(sh_i, 1e-3), 0, 1)
        x, y, z = c[:, 0], c[:, 1], c[:, 2]
        r = np.sqrt((x-GEO['Lx']/2)**2 + (y-GEO['Ly']/2)**2)
        zf = z[(r < W*0.55)].min() if (r < W*0.55).any() else GEO['sub_top']
        floor = (r < W*0.55) & (z < zf + 1.2)
        ar = (GEO['sub_top'] - z[floor].mean()) / W
        row.append((ar, float(ratio[floor].mean())))
    s = "  ".join(f"AR{ar:.0f}:{rr:.2f}" for ar, rr in row)
    err = np.mean([abs(rr - np.interp(ar, list(HG), list(HG.values()))) for ar, rr in row])
    print(f"  e_ang_sigma={e_sig:.1f}: {s}   |HG-err|={err:.3f}", flush=True)
print("\n  pick the e_ang_sigma with smallest |HG-err| -> set par['e_ang_sigma'].")
