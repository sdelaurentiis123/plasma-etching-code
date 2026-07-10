"""Preintegration estimator for the electron FLOOR-band flux fraction -- the scalable, accurate,
noise-controlled replacement for the tensor-grid electron floor (which OOMs at high AR).

Method (Griewank-Kuo-Sloan preintegration / conditional QMC, validated vs MC to within physics tol):
  outer scrambled-Sobol over (E, phase, az, launch-x); inner = fixed multi-band-safe ct-scan over the
  analytically crossed range [sqrt(Vs/E), 1] with the density integrated as int 2ct*fate dct. The fate
  oracle is the SAME _trace_general the MC uses (no reimplemented bookkeeping). Cost = N_outer * n_inner
  traces -- FLAT in aspect ratio (no tensor-product blowup), so it scales to arbitrarily high AR.

Returns the floor-band electron flux in the grid[band].mean()-per-crossed convention (same as
deterministic_electron_transport's floor value and the MC reference 0.0794 at AR6), so the caller can
override the moderate-resolution floor directly: ce[floor band cells] = P_acc * scale.
"""
import numpy as np

from .charging_general import _trace_general

_TWO_PI = 2.0 * np.pi


def preint_floor_fraction(g, Ex, Ez, n_log2=13, n_scramble=4, n_inner=64,
                          Te=4.0, V_dc=37.0, V_rf=30.0, trace_dt=0.15,
                          trace_dt_field=0.10, trace_steps=120, tol=5e-3):
    """Floor-band electron flux fraction on the field (Ex,Ez), via preintegration+QMC. Auto-stops
    when the scramble spread is under `tol`. Returns (value, stderr)."""
    from scipy.stats import qmc, gamma
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']
    fz = np.where(g['floor_trench_mask'].any(axis=0))[0].max()
    b0, b1 = t0 + 4, t1 - 4                    # floor band columns (edge cells excluded, per the ref)
    conv = float(t1 - t0) / float(b1 - b0)     # (ncol/band) convention factor -> grid[band].mean() scale
    msteps = int(trace_steps) * nz

    def fate(E, ct, phase, az, x):
        Vs = V_dc + V_rf * np.sin(phase)
        vz = np.sqrt(np.maximum(E * ct * ct - Vs, 0.0))
        st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
        vp = np.sqrt(E) * st * np.cos(az)
        hx, hz, _, _, _ = _trace_general(Ex, Ez, solid, x.astype(float), np.ones_like(E), vp, vz,
                                         -1.0, nx, nz, msteps, trace_dt, trace_dt_field)
        return ((hz == fz) & (hx >= b0) & (hx < b1)).astype(np.float64)

    def one(seed):
        s = qmc.Sobol(d=4, scramble=True, seed=seed)
        u = s.random_base2(n_log2)
        E = gamma.ppf(u[:, 0], a=2.0, scale=Te)
        phase = u[:, 1] * _TWO_PI; az = u[:, 2] * _TWO_PI; x = t0 + u[:, 3] * (t1 - t0)
        Vs = V_dc + V_rf * np.sin(phase)
        crossed = E > Vs
        ct_lo = np.sqrt(np.clip(Vs / E, 0.0, 1.0))
        span = 1.0 - ct_lo
        num = np.zeros_like(E)
        for j in (np.arange(n_inner) + 0.5) / n_inner:      # multi-band-safe ct-scan
            ctj = ct_lo + j * span
            num += 2.0 * ctj * fate(E, ctj, phase, az, x) * (span / n_inner)
        num = np.where(crossed, num, 0.0)
        den = np.where(crossed, 1.0 - ct_lo ** 2, 0.0)
        return num.sum() / max(den.sum(), 1e-12)

    vals = []
    for k in range(n_scramble):
        vals.append(one(k))
        if len(vals) >= 2 and np.std(vals) / np.sqrt(len(vals)) < tol:
            break
    vals = np.array(vals)
    # P(land|cross) * (ncol/band) = the grid[band].mean()-per-crossed convention (== MC 0.0794 @ AR6)
    return float(vals.mean() * conv), float(vals.std() / np.sqrt(len(vals)) * conv)


def preint_floor_ion_fraction(g, Ex, Ez, n_log2=11, n_scramble=4, n_inner=96,
                              Te=4.0, Ti=0.5, V_dc=37.0, V_rf=30.0, trace_dt=0.15,
                              trace_dt_field=0.10, trace_steps=120, tol=5e-3):
    """Floor-band ION flux fraction on the field (Ex,Ez), via preintegration+QMC -- the ion twin of
    preint_floor_fraction, fixing the deep-AR ci under-resolution (the tensor v_perp step ~1.8 deg
    equals the AR15 acceptance half-cone, so the floor ci was a 1-3 quadrature-node quantity).

    Source = sample_sheath_source(kind='ion') exactly: ions ALL cross (accelerated, no retardation),
    phase ~ uniform importance-weighted w~Vs^-0.35 (HG Fig-4a IEDF asymmetry), vz=sqrt(0.5*Te+Vs),
    v_perp ~ N(0, sqrt(0.5*Ti)) conserved. Outer scrambled-Sobol over (phase, x); inner = full-range
    v_perp scan via the probit map v_perp = sigma*Phi^-1(u) (Gaussian density absorbed EXACTLY, nodes
    auto-concentrated in the near-vertical acceptance cone, multi-band-safe through the warped field --
    at high floor V the acceptance is field-deflected, so we scan the whole crossed range like the
    electron ct-scan does). Fate oracle = _trace_general verbatim (q=+1). Cost flat in AR.

    Same grid[band].mean()-per-launched convention as the deterministic tensor floor, so the caller
    overrides: ci[floor band] = P_acc * scale (shape-preserving). Returns (value, stderr)."""
    from scipy.stats import qmc, norm
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']
    fz = np.where(g['floor_trench_mask'].any(axis=0))[0].max()
    b0, b1 = t0 + 4, t1 - 4
    conv = float(t1 - t0) / float(b1 - b0)
    msteps = int(trace_steps) * nz
    sigma = np.sqrt(0.5 * Ti)
    vp_nodes = sigma * norm.ppf((np.arange(n_inner) + 0.5) / n_inner)   # probit-mapped v_perp scan

    def fate(vp, vz, x):
        hx, hz, _, _, _ = _trace_general(Ex, Ez, solid, x.astype(float), np.ones_like(vz), vp, vz,
                                         1.0, nx, nz, msteps, trace_dt, trace_dt_field)
        return ((hz == fz) & (hx >= b0) & (hx < b1)).astype(np.float64)

    def one(seed):
        s = qmc.Sobol(d=2, scramble=True, seed=seed)
        u = s.random_base2(n_log2)
        phase = u[:, 0] * _TWO_PI
        x = t0 + u[:, 1] * (t1 - t0)
        Vs = V_dc + V_rf * np.sin(phase)
        w = Vs ** -0.35                       # IEDF asymmetry weight (HG Fig 4a), matches MC+tensor
        vz = np.sqrt(0.5 * Te + Vs)           # Bohm entry + instantaneous sheath gain
        num = np.zeros_like(Vs)
        for vpj in vp_nodes:                  # inner: E[fate] over v_perp ~ N(0, sigma), density exact
            num += fate(np.full(Vs.size, vpj), vz, x) / n_inner
        return float((w * num).sum() / w.sum())

    vals = []
    for k in range(n_scramble):
        vals.append(one(k))
        if len(vals) >= 2 and np.std(vals) / np.sqrt(len(vals)) < tol:
            break
    vals = np.array(vals)
    # P(land in band) * (ncol/band) = grid[band].mean() per-launched-per-column convention
    return float(vals.mean() * conv), float(vals.std() / np.sqrt(len(vals)) * conv)
