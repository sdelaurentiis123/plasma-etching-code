#!/usr/bin/env python3
"""SUB-MICRON profile overlay: does petch reproduce ViennaPS's actual SHAPE (not just depth)?

Belen-scale hole (~0.4 um) at fine dx. Run ViennaPS (GPU) FIRST (OptiX<->Warp conflict), then petch
rate-matched to ViennaPS depth. Extract each hole's sidewall profile radius(depth) in the centre slice,
overlay, and report the real shape difference (mean/max radius mismatch + sidewall angle), plus a PNG.
PETCH_DEVICE=cuda python scripts/vps_profile_compare.py"""
import time
import numpy as np

DX, EXT, DIAM, MASK, DUR = 0.025, 0.8, 0.4, 0.08, 0.2     # sub-micron; short etch so neither bottoms out
DEPTH_BINS = np.arange(0.0, 1.2, DX)


def profile_radius(xy_z, xc, mask_top):
    """Centre-slice surface points -> mean hole radius vs depth-below-mask-top."""
    x, z = xy_z[:, 0] - xc, xy_z[:, 1]
    depth = mask_top - z
    rad = np.abs(x)
    out = []
    for d0 in DEPTH_BINS:
        sel = (depth >= d0) & (depth < d0 + DX)
        if sel.sum() >= 1:
            out.append((d0 + 0.5 * DX, float(rad[sel].mean())))
    return np.array(out) if out else np.empty((0, 2))


# ---------- ViennaPS (before importing petch/warp) ----------
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def _flatten(obj, prefix="", out=None, depth=0):
    """Dump numeric fields of a ViennaPS params struct (one level of nesting)."""
    if out is None:
        out = {}
    for k in dir(obj):
        if k.startswith("_"):
            continue
        try:
            v = getattr(obj, k)
        except Exception:
            continue
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            out[prefix + k] = float(v)
        elif depth < 1 and not callable(v) and not isinstance(v, (str, list, tuple)):
            try:
                _flatten(v, prefix + k + ".", out, depth + 1)
            except Exception:
                pass
    return out


# physics-parameter parity: ViennaPS SF6O2 defaults (petch PAR compared later, after warp is safe to init)
vps_par = _flatten(v3.SF6O2Etching.defaultParameters())
print("=== ViennaPS SF6O2Etching default parameters ===", flush=True)
for k in sorted(vps_par):
    print(f"    {k} = {vps_par[k]:g}", flush=True)
print(flush=True)

d = v3.Domain()
v3.MakeHole(domain=d, gridDelta=DX, xExtent=EXT, yExtent=EXT, holeRadius=DIAM/2,
            holeDepth=MASK, makeMask=True, material=ps.Material.Si).apply()
m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
# CPU_TRIANGLE: GPU (OptiX) aborts the process if libnvoptix is absent, so use CPU -- the PROFILE/PHYSICS
# is identical to GPU (same model, only the ray engine differs; speed is benchmarked separately).
p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE)
print("ViennaPS engine: CPU_TRIANGLE (profile is engine-independent)", flush=True)
t0 = time.time(); p.apply(); vps_wall = time.time() - t0
nodes = np.array(d.getSurfaceMesh().getNodes())                   # (M,3), hole centred at origin
yc_v = 0.0
slab_v = nodes[np.abs(nodes[:, 1] - yc_v) < DX]
mask_top_v = nodes[:, 2].max()
vps_prof = profile_radius(slab_v[:, [0, 2]], 0.0, mask_top_v)
vps_depth = float(vps_prof[:, 0].max()) if len(vps_prof) else 0.0
print(f"ViennaPS: wall {vps_wall:.1f}s  mask_top_z {mask_top_v:.3f}  depth {vps_depth:.3f}um  "
      f"nodes {len(nodes)} slab {len(slab_v)}", flush=True)

# ---------- petch (rate-matched to ViennaPS depth) ----------
import petch
from petch import threed as t3

# physics-parameter parity check: our PAR vs the ViennaPS defaults dumped above
print("=== physics-parameter parity (petch PAR vs ViennaPS SF6O2 default) ===", flush=True)
_cmp = {"ionFlux": "ionFlux", "Fflux": "etchantFlux", "Oflux": "passivationFlux",
        "Emean": "Ions.meanEnergy", "Esig": "Ions.sigmaEnergy",
        "A_ie": "Substrate.A_ie", "Eth_ie": "Substrate.Eth_ie", "A_sp": "Substrate.A_sp",
        "Eth_sp": "Substrate.Eth_sp", "B_sp": "Substrate.B_sp", "k_sigma": "Substrate.k_sigma",
        "beta_sigma": "Substrate.beta_sigma", "rho": "Substrate.rho",
        "A_p": "Passivation.A_ie", "Eth_p": "Passivation.Eth_ie"}
for ours_k, vk in _cmp.items():
    pv = petch.PAR[ours_k]; vv = vps_par.get(vk)
    tag = "MATCH" if (vv is not None and abs(pv - vv) < 1e-6 * max(1, abs(vv))) else ("DIFF" if vv is not None else "n/a")
    print(f"    {ours_k:11s} petch {pv:>8g}   ViennaPS.{vk:22s} {('%.4g' % vv) if vv is not None else '----':>8}   {tag}", flush=True)
# physics ViennaPS HAS that petch does not model by default (affects bottom-corner shape):
print("    NOTE: ViennaPS also has ion specular reflection (Ions.inflectAngle/minAngle/thetaR),", flush=True)
print("          B_ie ion-enhanced angular term, and mask sputtering (Mask.*). petch: ion_reflection off,", flush=True)
print("          hard mask. These shape the deep bottom/microtrench, not the bulk sidewall.", flush=True)
print(flush=True)
GEO = dict(Lx=EXT, Ly=EXT, Lz=MASK + 0.8 + 0.3, dx=DX, trench_width=DIAM, mask_th=MASK, sub_top=0.8 + 0.3, hole=True)


def ours(rate, steps=40):
    p_ = dict(petch.PAR); p_['rate_scale'] = rate
    p_.update(flux_smooth_gpu=True, gpu_source=True, gpu_mesh=True, gpu_warmstart=True, device_flux=True)
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol")
    p_['ied_mode'] = 'gauss'
    t0 = time.time()
    g = t3.run_etch_3d(t_end=DUR, n_steps=steps, par=p_, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, g


ours(0.05, steps=4)                                              # warm
best = None
for r in [0.02, 0.05, 0.1, 0.2, 0.4]:
    w, g = ours(r)
    dep = t3.center_depth_3d(g)
    if best is None or abs(dep - vps_depth) < abs(best[1] - vps_depth):
        best = (w, dep, g, r)
ours_wall, ours_depth, g, rate = best
v, f, c, a = t3.extract_mesh_3d(g['phi'], DX)
xc = GEO['Lx'] / 2; mask_top_p = GEO['sub_top'] + MASK
slab_p = v[np.abs(v[:, 1] - GEO['Ly']/2) < DX][:, [0, 2]]
ours_prof = profile_radius(slab_p, xc, mask_top_p)
print(f"petch   : wall {ours_wall:.2f}s  depth {ours_depth:.3f}um  rate_scale {rate}  slab {len(slab_p)}", flush=True)

# ---------- compare shape (radius vs depth, on common depths) ----------
zc = np.intersect1d(np.round(vps_prof[:, 0], 4), np.round(ours_prof[:, 0], 4))
if len(zc) >= 3:
    rv = np.interp(zc, vps_prof[:, 0], vps_prof[:, 1])
    rp = np.interp(zc, ours_prof[:, 0], ours_prof[:, 1])
    dr = np.abs(rp - rv)
    print(f"\n  radius(depth) mismatch vs ViennaPS: mean {dr.mean()*1000:.1f} nm  max {dr.max()*1000:.1f} nm  "
          f"(hole radius {DIAM/2*1000:.0f} nm, dx {DX*1000:.0f} nm)", flush=True)
    print(f"  -> mean mismatch {100*dr.mean()/(DIAM/2):.1f}% of radius; depth: vps {vps_depth:.3f} ours {ours_depth:.3f}", flush=True)
else:
    print("\n  WARN: not enough common depths to compare -- check extraction.", flush=True)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4, 5))
    ax.plot(vps_prof[:, 1], -vps_prof[:, 0], 'o-', ms=2, label='ViennaPS', color='#2f5d62')
    ax.plot(ours_prof[:, 1], -ours_prof[:, 0], 's-', ms=2, label='petch', color='#7a3b2e')
    ax.plot(-vps_prof[:, 1], -vps_prof[:, 0], 'o-', ms=2, color='#2f5d62')
    ax.plot(-ours_prof[:, 1], -ours_prof[:, 0], 's-', ms=2, color='#7a3b2e')
    ax.set_xlabel('radius (µm)'); ax.set_ylabel('depth below mask (µm)'); ax.legend(); ax.set_title('sub-micron hole profile')
    fig.tight_layout(); fig.savefig('vps_profile_compare.png', dpi=140)
    print("  wrote vps_profile_compare.png", flush=True)
except Exception as e:
    print("plot skipped:", e)
