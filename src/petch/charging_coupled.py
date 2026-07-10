"""Dynamic (self-consistent) charging <-> etch coupling.

The static path (threed.py _apply_hg_charging, mode "petch") looks up (Q, Vf, E_defl) from a
precomputed AR-1..4 table (charging_general.petch_floor_profile) and CLAMPS for AR>4 -- it cannot
represent high-aspect-ratio charging. This module re-solves charging on the ACTUAL front geometry at
the CURRENT aspect ratio each etch checkpoint, warm-started from the previous solution, so the floor
throttle and deflected-ion foot energy track the real evolving feature -- the capability no open tool
has, and the one the marquee HAR (3D-NAND) application needs.

`dynamic_floor_profile(AR, ...)` is a drop-in replacement for `petch_floor_profile(AR)`: same
(Q, Vf, E_defl) return, but computed live from the derivation-clean solver (source_model="sheath",
charge_update="log"), with optional warm-start seeding for cheap re-solves across checkpoints.
"""
import numpy as np

from .charging2d import _build_edge_array_geometry
from .charging_general import (solve_charging, _trace_general, sample_sheath_source,
                               GAS, INSULATOR, CONDUCTOR)


def _extract(g, r, arrival):
    """(Q, Vf, E_defl) from a solved charging result on edge-array geometry g.
    Q = floor-center ion flux / field-free arrival; Vf = floor-center potential; E_defl = mean impact
    energy on the poly-inner sidewall FACE (the deflected-ion notch driver, C15 corrected definition)."""
    t0, t1 = g['trench0'], g['trench1']; nx = g['nx']
    fz = np.where(g['floor_trench_mask'].any(axis=0))[0].max()
    e1 = g['edge1']; zp = g['z_poly0']
    d = r['ntot'] / nx
    flux = r['ion_counts'][t0 + 4:t1 - 4, fz].mean() / d
    Q = float(flux / max(arrival, 1e-9))
    Vf = float(r['Vs'][t0 + 4:t1 - 4, fz].mean())
    cEi = r['ion_counts'][e1 - 1, zp:fz]; EEi = r['ion_energy'][e1 - 1, zp:fz]
    E_defl = float(EEi.sum() / max(cEi.sum(), 1e-9))
    return Q, Vf, E_defl


def _field_free_arrival(g):
    """Geometric floor ion arrival with NO charging field (denominator for Q). Derived sheath ions
    through zero field. ~independent of the charging state, so it can be cached per geometry."""
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']
    fz = np.where(g['floor_trench_mask'].any(axis=0))[0].max()
    rng = np.random.default_rng(11)
    x, z, vx, vz = sample_sheath_source(120000, rng, nx, "ion")
    Ex = np.zeros((nx, nz)); Ez = np.zeros((nx, nz))
    hix, hiz, E, _, _, _, _ = _trace_general(Ex, Ez, solid, x, z, vx, vz, 1.0, nx, nz, 160 * nz, 0.15, 0.10)
    cnt = np.zeros((nx, nz)); m = hix >= 0
    np.add.at(cnt, (hix[m], hiz[m]), 1.0)
    return cnt[t0 + 4:t1 - 4, fz].mean() / (120000 / nx)


def dynamic_floor_profile(AR, W=16, mouth=80, n_iter=800, n_per_iter=6000, seed_state=None,
                          return_state=False, arrival=None, leak_rate=0.0):
    """Live charging solve on the front at aspect ratio AR. Drop-in for petch_floor_profile(AR).
    Returns (Q, Vf, E_defl); if return_state, also (state, arrival) for warm-starting the NEXT
    checkpoint (state seeds the solver; arrival caches the geometric denominator).

    Warm-start: pass the previous checkpoint's `state` as seed_state and its `arrival`. The solver
    seeds only when the grid shape matches (same geometry), else it safely cold-solves."""
    g = _build_edge_array_geometry(float(AR), W=W, mouth=mouth)
    mat = np.where(g['solid'], np.where(g['cond'] > 0, CONDUCTOR, INSULATOR), GAS).astype(np.int64)
    if arrival is None:
        arrival = _field_free_arrival(g)
    r = solve_charging(mat, mouth=mouth, field_model='laplace', electron_model='trace',
                       electron_open_vf=False, charge_update="log", source_model="sheath",
                       trace_dt=0.15, trace_dt_field=0.10, trace_steps=120,
                       n_per_iter=n_per_iter, n_iter=n_iter, seed=7, seed_state=seed_state,
                       leak_rate=leak_rate)
    Q, Vf, E_defl = _extract(g, r, arrival)
    if return_state:
        # build the seed for the NEXT checkpoint from the standard return keys (the solver seeds
        # only when the grid shape matches, so this warm-starts a same-geometry re-solve).
        st = {'Vs': r['Vs'], 'V': r['V'], 'sig_sheet': r.get('sigma'), 'Vc': r['Vc']}
        return (Q, Vf, E_defl), (st, arrival)
    return Q, Vf, E_defl
