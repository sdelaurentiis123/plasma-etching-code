"""Grid-native 1-D Knudsen conductance tail for the deep-AR neutral flux.

Ported from Craig Xu Chen's `plasma_sim/solver3d.py` (`_neutral_flux_knudsen`),
adapted to petch's grid conventions: phi is indexed [ix, iy, iz], z is the
vertical (depth) axis (`geo['zs']`), the substrate top is at `z = sub_top`, the
feature is etched downward into `z < sub_top`. NOTE petch uses phi > 0 = SOLID,
phi < 0 = gas (the opposite sign from plasma_sim's phi < 0 = Si).

Physics: model the open feature as a vertical duct. From the level set, extract
per-depth gas cross-sectional area a(z) and sidewall perimeter p(z); form a
hydraulic radius r = 2a/p and a free-molecular (Knudsen) conductance
C(z) = (2/3) r a / dx. Link adjacent depth slices in series by harmonic-mean
conductance, impose unit neutral concentration at the opening (top BC) and a
coverage-weighted reaction sink at floor surfaces, then Thomas-solve the
tridiagonal system for the neutral concentration profile down the feature.

This is a reduced 1-D transport model intended as a CHEAP high-AR neutral tail
(a diffusion-synthetic-acceleration-style alternative to many-bounce MC), NOT an
accuracy fix: petch's Russian-roulette MC already samples the free-molecular
regime. Validate before claiming any accuracy delta. See memory
[[reconcile-craig-into-petch]] / [[accuracy-yardstick-vs-viennaps]].
"""
from __future__ import annotations

import numpy as np


def _via_slice_geometry(gas, below_iz, dx):
    """Per-depth open area and sidewall perimeter for an axisymmetric (hole) feature.

    `gas` is the boolean gas mask [nx, ny, nz]; `below_iz` are the z-slice indices
    inside the substrate (z < sub_top). Returns area[nz], perimeter[nz]."""
    nz = gas.shape[2]
    area = np.zeros(nz)
    perimeter = np.zeros(nz)
    for iz in below_iz:
        g = gas[:, :, iz]
        if not g.any():
            continue
        area[iz] = float(g.sum()) * dx * dx
        faces = (np.count_nonzero(g[:-1, :] != g[1:, :])
                 + np.count_nonzero(g[:, :-1] != g[:, 1:]))
        perimeter[iz] = max(float(faces) * dx, dx)
    return area, perimeter


def _trench_slice_geometry(gas, below_iz, dx):
    """Per-depth open area and perimeter for a translationally-invariant trench.

    The trench runs along y; the open width is measured in x. Area is per unit
    trench length (mean open x-width); perimeter is the two long sidewalls."""
    nz = gas.shape[2]
    area = np.zeros(nz)
    perimeter = np.zeros(nz)
    for iz in below_iz:
        widths = gas[:, :, iz].sum(axis=0) * dx        # open x-width per y-row
        open_widths = widths[widths > 0.0]
        if open_widths.size == 0:
            continue
        area[iz] = float(open_widths.mean())
        perimeter[iz] = 2.0
    return area, perimeter


def slice_geometry(phi, zs, dx, sub_top, feature_shape):
    """Per-depth gas area a(z) and sidewall perimeter p(z) from the level set."""
    gas = phi < 0.0                                    # petch: phi < 0 = gas, phi > 0 = solid
    below_iz = np.flatnonzero(zs < sub_top)
    if feature_shape == "trench":
        return _trench_slice_geometry(gas, below_iz, dx)
    return _via_slice_geometry(gas, below_iz, dx)


def _solve_tridiagonal(lower, diag, upper, rhs):
    """Thomas algorithm. Verbatim from plasma_sim (numerically guarded)."""
    n = len(diag)
    if n == 0:
        return np.zeros(0)
    a = lower.astype(float).copy()
    b = diag.astype(float).copy()
    c = upper.astype(float).copy()
    d = rhs.astype(float).copy()
    for i in range(1, n):
        w = a[i - 1] / max(b[i - 1], 1.0e-30)
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    x = np.zeros(n)
    x[-1] = d[-1] / max(b[-1], 1.0e-30)
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / max(b[i], 1.0e-30)
    return x


def conductance_profile(phi, zs, dx, sub_top, feature_shape, slice_loss,
                        wall_loss_scale=1.85):
    """Solve the 1-D Knudsen conductance system for the neutral depth profile.

    `slice_loss[nz]` : per-slice kinetic absorption 0.25*Sigma(s_face*A_face)/norm — the
    coverage-weighted consumption of ALL surface faces in that depth slice (geometry-invariant:
    a flat floor reduces to the classic 0.25*s*a floor term; a rounded evolving front
    distributes the same total consumption over its slices instead of double-counting).
    Returns profile[nz] in [0, 1]: neutral concentration vs depth, 1 at/above the
    opening, decaying with depth. Slices outside the feature default to 1.0."""
    nz = phi.shape[2]
    area, perimeter = slice_geometry(phi, zs, dx, sub_top, feature_shape)
    valid = (zs < sub_top) & (area > 0.0) & (perimeter > 0.0)
    profile = np.ones(nz)
    if not valid.any():
        return profile

    idx = np.flatnonzero(valid)[::-1]                  # opening slice first, deepest last
    m = len(idx)
    a = area[idx]
    p = perimeter[idx]
    radius = np.maximum(2.0 * a / np.maximum(p, 1.0e-30), dx)
    conductance = (2.0 / 3.0) * radius * a / dx        # free-molecular (Knudsen) conductance
    # wall_loss_scale is a scalar (default: uniform effective floor loss) OR a per-slice array of
    # length nz (opt-in passivation-linked front loss -- see knudsen_face_flux). Per-slice lets the
    # effective F-consumption enhancement relax with depth as the sidewall column passivates.
    ls = np.asarray(wall_loss_scale, float)
    if ls.ndim == 0:
        loss_scale = np.maximum(float(ls), 0.0)
    else:
        loss_scale = np.maximum(ls[idx], 0.0)
    sink = loss_scale * np.maximum(slice_loss[idx], 0.0)

    link = np.zeros(max(m - 1, 0))
    for i in range(m - 1):
        if abs(idx[i + 1] - idx[i]) != 1:              # non-adjacent slices: no link
            continue
        c0, c1 = conductance[i], conductance[i + 1]
        link[i] = 2.0 * c0 * c1 / max(c0 + c1, 1.0e-30)   # series (harmonic-mean) conductance

    diag = sink.copy()
    lower = np.zeros(max(m - 1, 0))
    upper = np.zeros(max(m - 1, 0))
    rhs = np.zeros(m)
    top = conductance[0]
    diag[0] += top
    rhs[0] += top                                      # unit-concentration opening BC
    for i, g in enumerate(link):
        diag[i] += g
        diag[i + 1] += g
        upper[i] = -g
        lower[i] = -g
    conc = _solve_tridiagonal(lower, np.maximum(diag, 1.0e-30), upper, rhs)
    profile[idx] = np.clip(conc, 0.0, 1.0)
    return profile


def _nearest_valid_slice(valid_iz, iz):
    """For each query slice in `iz`, nearest slice index that is in `valid_iz`."""
    if valid_iz.size == 0:
        return np.full(len(iz), -1, dtype=np.int64)
    pos = np.searchsorted(valid_iz, iz)
    right = np.clip(pos, 0, valid_iz.size - 1)
    left = np.clip(pos - 1, 0, valid_iz.size - 1)
    choose_right = np.abs(valid_iz[right] - iz) < np.abs(valid_iz[left] - iz)
    return np.where(choose_right, valid_iz[right], valid_iz[left])


def _front_loss_scale(zs, dx, sub_top, nz, W, wall_loss_scale, band_W, pass_frac, ar_pass):
    """OPT-IN passivation-linked wall loss (flags.knudsen_front_loss): a per-slice effective
    loss-scale array replacing the uniform scalar. PHYSICS: the >1 wall_loss_scale is a proxy for
    F recombination on the reactive (contested) sidewall band adjacent to the etch front. On a
    SHALLOW feature the whole sidewall column is freshly exposed -> full recombination -> full
    scale. As the feature deepens, only a FIXED fresh band (height band_W feature-widths) near the
    front stays contested; the tall column above it is SiOxFy-passivated and inert to F (Blauw JVST
    B 18,3453; Coburn-Winters). So the depth-integrated recombination the neutral flux sees relaxes
    from the full scale (shallow) toward pass_frac*scale (deep, mostly-passivated walls) on an AR
    decay scale ar_pass. Per-slice: each slice's local AR = (sub_top - z)/W sets its fresh fraction.
        scale(AR) = wls * (pass_frac + (1-pass_frac) * exp(-max(0, AR - band_W)/ar_pass))
    band_W: fresh-band height in feature widths (full strength through AR<=band_W -> preserves the
    knee). pass_frac: passivated-wall loss as a fraction of full (5-20% for a truly inert wall; the
    de Boer floor calibrates nearer 0.5 = the bare kinetic floor term surviving under the fudge)."""
    zc = zs[:nz]
    ar_slice = np.maximum((sub_top - zc) / max(W, 1.0e-9), 0.0)
    g = np.exp(-np.maximum(ar_slice - band_W, 0.0) / max(ar_pass, 1.0e-9))
    return max(float(wall_loss_scale), 0.0) * (pass_frac + (1.0 - pass_frac) * g)


def knudsen_face_flux(phi, zs, dx, sub_top, feature_shape, centroids, gas_nz, face_areas, Ly,
                      s_face, wall_loss_scale=1.85, center=None, half_width=None,
                      front_loss=False, band_W=8.0, pass_frac=0.5, ar_pass=15.0):
    """Per-face neutral arrival multiplier from the 1-D Knudsen conductance solve.

    centroids  : (F,3) face centroids; gas_nz : (F,) into-gas normal z (floor classifier);
    face_areas : (F,) triangle areas; Ly : trench length (per-unit-length norm; ignored for vias);
    s_face     : (F,) per-face sticking (= bare*beta, the coverage physics).
    LITERATURE MODEL (Coburn-Winters APL 55, 2730; Blauw JVST B 18, 3453 "negligible sidewall F
    loss"): the reaction sink lives at the BOTTOM only; sidewalls are elastic diffuse scatterers
    (cryo-passivated). Per-slice sink = 0.25*Sigma(s_face*A_face) over the FLOOR-classified faces
    (gas_nz > 0.7) in that slice — a flat floor reduces exactly to the classic 0.25*s*a term, and
    a rounded evolving front contributes its ACTUAL bottom area once (the old form applied the
    full duct-area sink at every slice the front touched — the evolving over-consumption bug).
    An all-faces sink was tried and rejected: the coupled coverage fixed point then has a
    strongly-attracting collapsed state (F and O both starve -> nothing passivates -> bare=1 ->
    max sticking -> stays starved). Returns m[F]."""
    nz = phi.shape[2]
    z0 = zs[0]
    area, perimeter = slice_geometry(phi, zs, dx, sub_top, feature_shape)
    valid_iz = np.flatnonzero((zs < sub_top) & (area > 0.0) & (perimeter > 0.0))
    iz_face = np.clip(np.round((centroids[:, 2] - z0) / dx).astype(int), 0, nz - 1)
    # Mesh-face centroids sit ON the zero contour, so a face's raw slice index can round into
    # the gas-free slice just below the last valid one -- the sink would land on an invalid slice
    # and never enter the solve (the flat-profile bug). Snap every face to its nearest VALID
    # (gas-carrying) slice before accumulating; band-cell-based codes (plasma_sim) get this for free.
    snap = _nearest_valid_slice(valid_iz, iz_face)
    below = centroids[:, 2] < sub_top
    # Footprint filter (plasma_sim's _feature_mask, restored): field-plane faces sit at sub_top-eps
    # after marching cubes, would count as "below", snap to the duct mouth, and dump the entire
    # field's absorption into the chain (collapses the profile). The field is the reservoir --
    # already represented by the unit-concentration top BC -- so only faces inside the feature
    # footprint belong in the sink.
    inside = np.ones(len(centroids), bool)
    if center is not None and half_width is not None:
        band = half_width + 2.0 * dx
        if feature_shape == "trench":
            inside = np.abs(centroids[:, 0] - center[0]) <= band
        else:
            inside = np.hypot(centroids[:, 0] - center[0], centroids[:, 1] - center[1]) <= band

    slice_loss = np.zeros(nz)
    sel = below & inside & (snap >= 0) & (np.asarray(gas_nz) > 0.7)   # bottom faces only (Blauw)
    if sel.any():
        sA = np.clip(np.asarray(s_face, float)[sel], 0.0, 1.0) * np.asarray(face_areas, float)[sel]
        np.add.at(slice_loss, snap[sel], 0.25 * sA)
        if feature_shape == "trench":
            slice_loss /= max(float(Ly), 1.0e-30)      # per unit trench length (area a is too)

    wls_arg = wall_loss_scale
    if front_loss:
        W = 2.0 * float(half_width) if half_width is not None else dx
        wls_arg = _front_loss_scale(zs, dx, sub_top, nz, W, wall_loss_scale,
                                    band_W, pass_frac, ar_pass)
    profile = conductance_profile(phi, zs, dx, sub_top, feature_shape,
                                  slice_loss, wls_arg)
    m = np.ones(len(centroids))
    has = snap >= 0
    m[has & below] = profile[snap[has & below]]
    return np.clip(m, 0.0, 8.0)
