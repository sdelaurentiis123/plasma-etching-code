"""Deterministic charged-particle transport for feature charging (Route B, the frontier engine).

Replaces the Monte-Carlo ion/electron tracers with a FIXED quadrature over the sheath launch phase
space (position x transverse-velocity x phase/energy), each ordinate integrated through the field and
deposited on impact. Noise-free (no particle count to sweep), scales to all aspect ratios without
sampling starvation (a resolved ordinate fan always covers the narrow near-vertical cone that reaches
a deep floor), and differentiable (smooth quadrature + soft-binning deposition -> autodiff-clean).

This module starts with the ION source and the FIELD-FREE arrival as the correctness gate: the same
distribution as charging_general.sample_sheath_source("ion"), expressed as weighted quadrature nodes
instead of random draws, traced through the SAME _trace_general kernel. If the deterministic arrival
reproduces the MC arrival (and is noise-free + converges under refinement), the quadrature+weighting
is correct and we build the charged (field-coupled) solve on top.

Design notes (per literature + advisor):
- transverse velocity v_perp ~ N(0, 0.5*Ti): GAUSS-HERMITE nodes (exact for the Gaussian, and the
  near-zero-v_perp = most-vertical = deepest-reaching rays are included BY CONSTRUCTION -> no starvation).
- phase -> energy: the instantaneous RF sheath Vs=V_dc+V_rf*sin(phase) with the HG Fig-4a importance
  weight w~Vs^-0.35 (the bimodal-bathtub IEDF asymmetry). Midpoint grid over phase, weighted by w.
- normalization matches the MC arrival EXACTLY: each trench-mouth column emits total weight 1, so the
  floor-band deposition is the per-column transmission probability -- identical estimand to the MC
  cnt[floor].mean()/(N/nx). Direct A/B against charging_general is therefore apples-to-apples.
"""
import numpy as np

from .charging2d import _build_edge_array_geometry
from .charging_general import _trace_general


def ion_source_quadrature(cols, nx, n_vperp=21, n_phase=64, Te=4.0, Ti=0.5, V_dc=37.0, V_rf=30.0):
    """Deterministic node set for the RF-sheath ion source. Returns (x, z, vperp, vz, weight) flat
    arrays; the per-column weight sums to 1 (a probability quadrature over the (phase, vperp) launch
    space), replicated across each launch column in `cols`."""
    two_pi = 2.0 * np.pi
    ph = (np.arange(n_phase) + 0.5) / n_phase * two_pi        # midpoint rule over phase
    Vs = V_dc + V_rf * np.sin(ph)
    w_ph = Vs ** -0.35; w_ph = w_ph / w_ph.sum()              # IEDF asymmetry weight (HG Fig 4a)
    vz_ph = np.sqrt(0.5 * Te + Vs)                            # Bohm entry + instantaneous sheath gain
    sigma = np.sqrt(0.5 * Ti)
    # v_perp quadrature: BOUNDED midpoint grid over [-4sigma, 4sigma] with Gaussian weights -- NOT
    # Gauss-Hermite. GH puts nodes out to +/-7 sigma, spawning unphysical wide-angle ions that
    # mis-integrate the SHARP geometric transmission cutoff (reaches floor iff angle < cutoff) and
    # drift into neighboring features. A bounded fine grid resolves the cutoff with uniform density
    # and stays inside the physical +/-4 sigma the real Gaussian occupies (advisor: grade/resolve the
    # cone, don't let sampling starvation become quadrature starvation).
    span = 4.0 * sigma
    edges = np.linspace(-span, span, n_vperp + 1)
    vperp_n = 0.5 * (edges[:-1] + edges[1:])
    w_vp = np.exp(-0.5 * (vperp_n / sigma) ** 2); w_vp = w_vp / w_vp.sum()
    # tensor (phase x vperp) velocity+weight set for ONE column; weight sums to 1
    VZ = np.repeat(vz_ph, n_vperp)
    VP = np.tile(vperp_n, n_phase)
    W = np.repeat(w_ph, n_vperp) * np.tile(w_vp, n_phase)
    ncol = len(cols)
    X = np.repeat(np.asarray(cols, float), VZ.size)
    Z = np.ones_like(X)
    return X, Z, np.tile(VP, ncol), np.tile(VZ, ncol), np.tile(W, ncol)


def deterministic_ion_transport(g, Ex, Ez, n_vperp=61, n_phase=96, Te=4.0, Ti=0.5,
                                V_dc=37.0, V_rf=30.0, trace_dt=0.15, trace_dt_field=0.10,
                                trace_steps=120):
    """Charged deterministic ion transport: the quadrature fan integrated through the ACTUAL field
    (Ex, Ez) -- trajectories bend. Drop-in for charging_general.trace('ion', ...), but noise-free and
    WEIGHTED. Returns (counts, energy) grids [per-column-normalized: each trench column emits weight 1].
    Same _trace_general kernel and same Boris-style push as the MC path, so identical physics."""
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']
    cols = np.arange(t0, t1)
    X, Z, VP, VZ, W = ion_source_quadrature(cols, nx, n_vperp=n_vperp, n_phase=n_phase,
                                            Te=Te, Ti=Ti, V_dc=V_dc, V_rf=V_rf)
    hix, hiz, E, _, _ = _trace_general(Ex, Ez, solid, X, Z, VP, VZ, 1.0, nx, nz,
                                       int(trace_steps) * nz, trace_dt, trace_dt_field)
    counts = np.zeros((nx, nz)); energy = np.zeros((nx, nz))
    m = hix >= 0
    np.add.at(counts, (hix[m], hiz[m]), W[m])
    np.add.at(energy, (hix[m], hiz[m]), W[m] * E[m])
    return counts, energy


def deterministic_arrival(g, n_vperp=21, n_phase=64):
    """Field-free floor arrival (Q denominator) from the deterministic ion quadrature. Same estimand
    as charging_general MC arrival, but noise-free. Straight-line rays (E=0), same _trace_general."""
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']
    fz = np.where(g['floor_trench_mask'].any(axis=0))[0].max()
    cols = np.arange(t0, t1)
    X, Z, VP, VZ, W = ion_source_quadrature(cols, nx, n_vperp=n_vperp, n_phase=n_phase)
    Ex = np.zeros((nx, nz)); Ez = np.zeros((nx, nz))
    hix, hiz, _, _, _ = _trace_general(Ex, Ez, solid, X, Z, VP, VZ, 1.0, nx, nz, 160 * nz, 0.15, 0.10)
    cnt = np.zeros((nx, nz)); m = hix >= 0
    np.add.at(cnt, (hix[m], hiz[m]), W[m])
    return float(cnt[t0 + 4:t1 - 4, fz].mean())
