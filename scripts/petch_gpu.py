"""petch-GPU (fudge-free, wrap-fixed, full GPU stack) speed vs ViennaPS-GPU in /root/vps_gpu.json, SAME
box, SAME holes. Warmed. PETCH_DEVICE=cuda."""
import os, time, json
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

vps = json.load(open("/root/vps_gpu.json"))
DX, DUR, XE = 0.25, 3.0, 14.0


def run(diam):
    p = dict(petch.PAR); p['rate_scale'] = 1.0
    p.update(flux_smooth_gpu=True, gpu_source=True, gpu_mesh=True, gpu_warmstart=True, device_flux=True, n_fp=1)
    GEO = dict(Lx=XE, Ly=XE, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=DUR)
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=diam, dx=DX, n_steps=40, par=p, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.max_depth_3d(g)


print(f"device={t3.DEVICE}  petch-GPU vs ViennaPS-GPU (OptiX) -- SAME BOX, holes d=4/6/8\n", flush=True)
print(f"  {'d':>3}  {'ViennaPS-GPU':>13}  {'petch-GPU':>11}  {'speedup':>8}", flush=True)
ratios = []
for diam in [4.0, 6.0, 8.0]:
    run(diam)                      # warm (compile + cache)
    w, dep = run(diam)
    vw = vps[str(diam)]['wall']; vdep = vps[str(diam)]['depth']
    r = vw / w; ratios.append(r)
    print(f"  {diam:>3}  {vw:6.2f}s d{vdep:4.1f}  {w:5.2f}s d{dep:4.1f}  {r:6.1f}x", flush=True)
print(f"\n  petch is {np.mean(ratios):.1f}x faster than ViennaPS-GPU (OptiX), same box, fudge-free.", flush=True)
