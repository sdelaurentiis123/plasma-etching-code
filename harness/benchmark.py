"""Fidelity benchmark: run a flags/par config and score it on two independent axes.

Accuracy axis  -> center depth (width-8) and ARDE curve vs the cached ViennaPS reference.
Speed axis     -> ray throughput (primary launches / raytrace second) and floor roughness
                  (a Monte-Carlo variance proxy; the thing QMC / the species split attack).

ViennaPS has no arm64 wheel on this Mac, so the accuracy reference is the cached SCALAR
depths in harness/reference/summary.json (no Chamfer/sidewall until a GPU/Linux box).
"""
import json
import os
import numpy as np

import petch

HERE = os.path.dirname(__file__)
REF_PATH = os.path.join(HERE, "reference", "summary.json")

# width-8 benchmark geometry (matches the original run_benchmark.py)
BENCH = dict(W=20.0, H=24.0, dx=0.25, mask_thickness=2.0, sub_top=18.0, t_end=3.0, n_steps=60)
ARDE_WIDTHS = [4.0, 6.0, 8.0, 12.0]


def load_reference():
    with open(REF_PATH) as f:
        return json.load(f)


def floor_roughness(result, frac=0.6):
    """Std of the floor depth across the central trench region — a variance/noise proxy."""
    seg = result['segs']
    if len(seg) == 0:
        return float('nan')
    W = (result['xs'][-1] - result['xs'][0]) + result['dx']
    xc, dep = petch.ours_profile(seg, result['sub_top'], W)
    if len(dep) == 0:
        return float('nan')
    dmax = dep.max()
    # floor = central, near-deepest points
    sel = (np.abs(xc) < 0.6 * W * 0.2) & (dep > frac * dmax)
    return float(np.std(dep[sel])) if sel.sum() >= 3 else float('nan')


def rays_per_sec(result, n_part_ion, n_part_neu, n_steps):
    """Primary particle launches per raytrace-second (throughput proxy)."""
    primaries = (n_part_ion + 2 * n_part_neu) * n_steps
    rt = result['timings']['raytrace']
    return primaries / rt if rt > 0 else float('nan')


def run_case(par, flags, width, n_part_ion=20000, n_part_neu=20000, **over):
    """Run one width and return the result dict plus derived metrics."""
    cfg = dict(BENCH); cfg.update(over)
    r = petch.run_etch(trench_width=width, par=par, flags=flags,
                       n_part_ion=n_part_ion, n_part_neu=n_part_neu, verbose=False, **cfg)
    r['_depth'] = petch.center_depth(r)
    r['_rough'] = floor_roughness(r)
    r['_rays_per_sec'] = rays_per_sec(r, n_part_ion, n_part_neu, cfg['n_steps'])
    return r


def scorecard(par, flags, label="", widths=ARDE_WIDTHS, n_part_ion=20000, n_part_neu=20000):
    """Full accuracy + speed scorecard for one (par, flags) configuration."""
    ref = load_reference()
    depths = {}
    speed = {}
    for w in widths:
        r = run_case(par, flags, w, n_part_ion=n_part_ion, n_part_neu=n_part_neu)
        depths[w] = r['_depth']
        if w == 8.0:
            speed = dict(rays_per_sec=r['_rays_per_sec'], rough=r['_rough'],
                         total_s=r['timings']['total'])
    d8 = depths.get(8.0, float('nan'))
    vps8 = ref['depth_vps']
    # ARDE: normalize to widest trench (matches the figure-(c) convention)
    widest = depths[widths[-1]]
    if widest <= 1e-6:   # saturated/failed run -> ARDE undefined
        arde_rmse = float('nan')
    else:
        our_arde = np.array([depths[w] for w in widths]) / widest
        vps_arde = np.array(ref['vps_depth']) / ref['vps_depth'][-1]
        arde_rmse = float(np.sqrt(np.mean((our_arde - vps_arde) ** 2)))
    return dict(label=label, depths=depths, d8=d8, vps8=vps8,
                d8_abs_err=d8 - vps8, d8_pct_err=100 * (d8 - vps8) / vps8,
                arde_rmse=arde_rmse, speed=speed,
                rate_scale=par.get('rate_scale', 1.0), chemistry=flags.chemistry)


def print_scorecard(sc):
    s = sc['speed']
    print(f"  [{sc['label']}]  chemistry={sc['chemistry']}  rate_scale={sc['rate_scale']}")
    print(f"    width-8 depth : {sc['d8']:.3f} um   vs ViennaPS {sc['vps8']:.3f}   "
          f"(Δ {sc['d8_abs_err']:+.3f} um, {sc['d8_pct_err']:+.1f}%)")
    print(f"    ARDE rmse     : {sc['arde_rmse']:.4f}   "
          f"(depths {', '.join('%.2f'%sc['depths'][w] for w in sorted(sc['depths']))})")
    if s:
        print(f"    speed         : {s['rays_per_sec']/1e6:.2f} M rays/s   "
              f"floor_rough {s['rough']:.3f} um   total {s['total_s']:.1f}s")
