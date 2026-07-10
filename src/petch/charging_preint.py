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
        hx, hz, _, _, _, _, _ = _trace_general(Ex, Ez, solid, x.astype(float), np.ones_like(E), vp, vz,
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
        hx, hz, _, _, _, _, _ = _trace_general(Ex, Ez, solid, x.astype(float), np.ones_like(vz), vp, vz,
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


def preint_wall_fraction(g, Ex, Ez, n_log2=13, n_scramble=4, n_inner=64,
                         Te=4.0, V_dc=37.0, V_rf=30.0, trace_dt=0.15,
                         trace_dt_field=0.10, trace_steps=120, tol=5e-3):
    """UPPER-SIDEWALL electron flux PROFILE on the field (Ex,Ez), via preintegration+QMC -- the WALL
    analog of preint_floor_fraction, fixing the deep-AR non-monotone floor by forming the physical
    electron-shading DIPOLE. The deterministic down-going quadrature fan structurally under-samples
    side-arriving electrons onto the upper insulator sidewalls, so those walls receive ~0 traced flux
    and charge POSITIVE (+34V) instead of negative -> the dipole never forms -> the floor under-charges.
    This estimator traces electrons through the ACTUAL field (same _trace_general fate oracle, same
    scrambled-Sobol source, same crossing-restricted ct-scan as the floor) and bins the ones that land
    on each wall cell -> the walls charge negative and self-limit.

    Three design points that make it CORRECT (do not "simplify" any):
    0. FULL-WIDTH launch (x over ALL nx columns), NOT the floor's aperture launch x in [t0,t1). The
       wide-angle electrons that physically charge the UPPER WALLS enter obliquely from launch positions
       OVER THE MASK, outside the trench aperture; only near-vertical electrons launched over the aperture
       reach the floor. Copying the floor's [t0,t1) launch truncated ~93% of the wall flux (measured 13.6x
       under-delivery, walls stayed positive). With full-width launch the per-cell normalization factor is
       conv = nx (not t1-t0): the P(cell|cross)*conv*scale count only matches the deterministic deposit
       when conv equals the launch-domain width. (Floor: launch width = aperture = t1-t0, so conv=ncol.)
    1. PER-CELL, not per-band. The dipole IS a depth gradient (less negative at mouth, more negative
       deep), so the accumulator is a per-cell histogram [2 walls, nrows], NOT a single band scalar.
       Convention mirrors the floor: a wall cell is a "band of size 1", so counts[c,r] = P(cell|cross)
       * nx * scale, i.e. flux[c,r]*scale with flux = Pcell*nx. Same Sobol/den/_trace_general as the
       floor => this IS the electron-beam-consistent flux the down-going fan WOULD deposit on those
       cells if it sampled side-arrivals; ion wall flux uses the same `scale`, so ci-ce=0 lands the wall
       on the right floating potential.
    2. NO explicit thr=exp(V/Te) factor. The crossing gate uses Vs = the MOUTH sheath potential
       (positive), NOT the wall surface potential, so it is invariant as the wall charges negative and
       never zeros the flux. Retardation of the negative wall lives ENTIRELY in the traced field inside
       _trace_general. An explicit thr would DOUBLE-COUNT retardation. This sign-discipline is why the
       identical (floor-positive) machinery is sign-correct for a negative-charging wall: sign-agnostic
       at the gate, sign-correct in the trace.

    Returns (mean, stderr), both shape (2, nrows): row 0 = LEFT wall (col trench0-1), row 1 = RIGHT wall
    (col trench1); columns index insulator rows [mouth : z_poly0]. Caller overrides
    counts[cL/cR, mouth:z_poly0] = flux*scale (shape-exact, no borrowed shape)."""
    from scipy.stats import qmc, gamma
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']
    cL, cR = t0 - 1, t1                                 # left/right trench INSULATOR walls
    r0, r1 = int(g['mouth']), int(g['z_poly0'])         # PR (insulator) rows; below r1 the wall is conductor
    nrows = r1 - r0
    conv = float(nx)                                    # FULL-WIDTH launch (see note): per-cell factor = nx
    if nrows <= 0:                                       # no insulator wall (all-conductor line): nothing to do
        return np.zeros((2, 0)), np.zeros((2, 0))
    msteps = int(trace_steps) * nz

    def hits(E, ct, phase, az, x):
        Vs = V_dc + V_rf * np.sin(phase)
        vz = np.sqrt(np.maximum(E * ct * ct - Vs, 0.0))
        st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
        vp = np.sqrt(E) * st * np.cos(az)
        hx, hz, _, _, _, _, _ = _trace_general(Ex, Ez, solid, x.astype(float), np.ones_like(E), vp, vz,
                                         -1.0, nx, nz, msteps, trace_dt, trace_dt_field)
        return hx, hz

    def one(seed):
        s = qmc.Sobol(d=4, scramble=True, seed=seed)
        u = s.random_base2(n_log2)
        E = gamma.ppf(u[:, 0], a=2.0, scale=Te)
        phase = u[:, 1] * _TWO_PI; az = u[:, 2] * _TWO_PI; x = u[:, 3] * (nx - 1)
        Vs = V_dc + V_rf * np.sin(phase)
        crossed = E > Vs
        ct_lo = np.sqrt(np.clip(Vs / E, 0.0, 1.0))
        span = 1.0 - ct_lo
        num = np.zeros((2, nrows))
        for j in (np.arange(n_inner) + 0.5) / n_inner:  # multi-band-safe ct-scan (same as floor)
            ctj = ct_lo + j * span
            hx, hz = hits(E, ctj, phase, az, x)
            w = np.where(crossed, 2.0 * ctj * (span / n_inner), 0.0)
            mL = (hx == cL) & (hz >= r0) & (hz < r1)
            np.add.at(num[0], (hz[mL] - r0).astype(np.intp), w[mL])
            mR = (hx == cR) & (hz >= r0) & (hz < r1)
            np.add.at(num[1], (hz[mR] - r0).astype(np.intp), w[mR])
        den = np.where(crossed, 1.0 - ct_lo ** 2, 0.0).sum()
        return num / max(den, 1e-12)                    # P(land in each wall cell | crossed)

    vals = []
    for k in range(n_scramble):
        vals.append(one(k))
        if len(vals) >= 2:                              # auto-stop on the total wall-flux spread (CRN seeds)
            tot = np.array([v.sum() for v in vals])
            if np.std(tot) / np.sqrt(len(vals)) < tol:
                break
    vals = np.array(vals)                               # (k, 2, nrows)
    return vals.mean(axis=0) * conv, vals.std(axis=0) / np.sqrt(len(vals)) * conv
