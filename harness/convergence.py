"""Convergence suite: separate discretization error (grid) from sampling noise (rays).

Key finding the harness exposes: at baseline (pseudorandom MC, ~20k rays) the per-run depth
noise is comparable to or larger than the grid-convergence signal. So both studies use
SEED-AVERAGING: at each setting we run several independent MC realizations (via the new
`seed` knob) and report the mean +/- the across-seed standard deviation.

- grid_refinement: refine dx -> 0 on the seed-averaged depth; Richardson order estimate.
- ray_refinement: at each N, std(depth) across seeds IS the Monte-Carlo noise; it should fall
  ~1/sqrt(N) for pseudorandom sampling (the baseline QMC must beat in Step 4).
"""
import numpy as np
import petch

# compact, fast convergence geometry
CONV = dict(W=12.0, H=14.0, trench_width=5.0, mask_thickness=2.0, sub_top=10.0,
            t_end=2.0, n_steps=40)


def depths_over_seeds(par, flags, dx, n_part, seeds):
    """Center depth for each MC seed at a given (dx, n_part)."""
    out = []
    for s in seeds:
        r = petch.run_etch(dx=dx, par=par, flags=flags, n_part_ion=n_part,
                           n_part_neu=n_part, seed=s, verbose=False, **CONV)
        out.append(petch.center_depth(r))
    return np.array(out)


def richardson(d_coarse, d_mid, d_fine, ratio=2.0):
    """Estimate continuum value and observed order from 3 successively halved grids."""
    num = d_coarse - d_mid
    den = d_mid - d_fine
    if abs(den) < 1e-12 or num / den <= 0:
        return d_fine, float('nan')
    p = np.log(num / den) / np.log(ratio)
    d0 = d_fine + (d_fine - d_mid) / (ratio ** p - 1.0)
    return float(d0), float(p)


def grid_refinement(par, flags, dxs=(0.5, 0.25, 0.125), n_part=20000, n_seeds=4):
    seeds = list(range(n_seeds))
    means, sems = [], []
    for dx in dxs:
        d = depths_over_seeds(par, flags, dx, n_part, seeds)
        m, sd = d.mean(), d.std()
        means.append(m); sems.append(sd / np.sqrt(len(seeds)))
        print(f"    dx={dx:<6} depth={m:.3f} +/- {sd/np.sqrt(len(seeds)):.3f} um  (seed std {sd:.3f})")
    out = dict(dxs=list(dxs), means=means, sems=sems)
    if len(means) >= 3:
        d0, p = richardson(means[-3], means[-2], means[-1])
        out['extrapolated'] = d0; out['order'] = p
        signal = abs(means[-1] - means[-3])
        print(f"    Richardson: continuum depth ~ {d0:.3f} um, order p ~ {p:.2f} "
              f"(grid signal {signal:.3f} um vs seed-noise {sems[-1]:.3f} um)")
    return out


def ray_refinement(par, flags, Ns=(5000, 10000, 20000, 40000), dx=0.25, n_seeds=5):
    seeds = list(range(n_seeds))
    means, stds = [], []
    for N in Ns:
        d = depths_over_seeds(par, flags, dx, N, seeds)
        means.append(d.mean()); stds.append(d.std())
        print(f"    N={N:<7} depth={d.mean():.3f} um   MC noise (seed std)={d.std():.3f} um")
    out = dict(Ns=list(Ns), means=means, stds=stds)
    # check noise ~ 1/sqrt(N): std * sqrt(N) should be ~constant
    consts = [s * np.sqrt(N) for s, N in zip(stds, Ns)]
    print(f"    std*sqrt(N) (flat => 1/sqrt(N) scaling): {['%.1f'%c for c in consts]}")
    out['std_sqrtN'] = consts
    return out
