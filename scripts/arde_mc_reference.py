"""Independent particle Monte Carlo reference for reactive ARDE floor transmission.

Completely separate from the common engine's diffuse-radiosity transport: this launches particles
from the source plane (half-Maxwellian cosine flux), traces straight-line hits against the analytic
periodic box (open region above the mask; the opening tube through mask+trench; the floor), and on
each surface hit either REACTS (absorbed, probability s) or DIFFUSELY RE-EMITS (cosine about the
surface normal, probability 1-s), until absorbed or escaping the top. It tallies floor-incident flux.
Agreement with the engine's radiosity (scripts/deboer_arde_static.floor_transmission) validates the
reactive transport by an independent method.

2D in (x,z): the geometry is translationally invariant along the cell-length axis and periodic in x
(period Wc), so only x and z matter for hits; each 3D cosine direction's y-component free-streams and
is ignored. Verified: at s=1 (no re-emission) this reproduces the direct-flight geometric reference
(reference_floor_transmission), which matches the opposed-strip view factor sqrt(1+A_eff^2)-A_eff;
across s=0.1..0.5 the engine radiosity agrees with this MC to ~1-3%.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import qmc

from petch.boundary_state import MaxwellianFluxVelocityDensity

EPS = 1e-9


def _cosine_inplane(normal_axis, sign, n, rng):
    """(dirx, dirz) projection of a 3D cosine distribution about a surface normal.

    ``normal_axis='z'`` -> normal (0,0,sign); ``'x'`` -> normal (sign,0,0). The dropped component is
    the free-streaming y direction."""
    u1 = rng.random(n); u2 = rng.random(n)
    ct = np.sqrt(1.0 - u1)
    st = np.sqrt(u1)
    phi = 2.0 * np.pi * u2
    a = st * np.cos(phi)
    if normal_axis == "z":
        return a, sign * ct
    return sign * ct, a


def _t_to_wall_value(x, dirx, p, Wc):
    """Smallest t>EPS with (x + dirx t) == p (mod Wc); inf where dirx==0."""
    t = np.full_like(x, np.inf)
    pos = dirx > 1e-15
    neg = dirx < -1e-15
    rem_p = np.mod(p - x, Wc); rem_p = np.where(rem_p < EPS, rem_p + Wc, rem_p)
    t[pos] = rem_p[pos] / dirx[pos]
    rem_n = np.mod(x - p, Wc); rem_n = np.where(rem_n < EPS, rem_n + Wc, rem_n)
    t[neg] = rem_n[neg] / (-dirx[neg])
    return t


def mc_reactive_transmission(ar, sticking, *, opening=0.10, mask=0.05, dx=0.02,
                             log2n=19, max_bounce=200, seed=7):
    """Floor incident flux / source flux for one (AR, sticking) point, by particle Monte Carlo.

    Geometry constants match ``scripts/deboer_arde_static`` exactly. Per-area normalization: the source
    spans the full cell (width Wc = 2*opening) and the floor spans the opening, so the incident-flux
    ratio is 2 * (floor arrivals / launched particles)."""
    etched = ar * opening
    substrate_top = etched + max(4.0 * dx, 0.05)
    floor_z = substrate_top - etched
    mask_top = substrate_top + mask
    Wc = 2.0 * opening
    source_z = substrate_top + mask + max(6.0 * dx, 0.06)
    ox0 = (Wc - opening) / 2.0
    ox1 = ox0 + opening
    s = float(sticking)

    N = 2 ** int(log2n)
    u = qmc.Sobol(3, scramble=True, seed=seed).random_base2(int(log2n))
    vel = MaxwellianFluxVelocityDensity(0.05).sample_flux_velocity(u)
    d = vel / np.linalg.norm(vel, axis=1, keepdims=True)
    dirx = d[:, 0].copy()
    dirz = -d[:, 2].copy()                 # downward (toward feature = -z)
    rng = np.random.default_rng(seed + 1)
    x = rng.random(N) * Wc
    z = np.full(N, source_z)
    active = np.ones(N, dtype=bool)
    floor_incident = 0.0

    def opening_mask(xh):
        f = np.mod(xh, Wc)
        return (f >= ox0) & (f <= ox1)

    for _ in range(int(max_bounce)):
        if not active.any():
            break
        ax = np.where(active)[0]
        cx, cz, cdx, cdz = x[ax], z[ax], dirx[ax], dirz[ax]
        n = cx.shape[0]
        inf = np.full(n, np.inf)
        tmin = inf.copy()
        etype = np.full(n, -1, dtype=int)          # 0 floor,1 mask,2 wall,3 escape,4 transition
        wall_is_ox0 = np.zeros(n, dtype=bool)
        cdz_safe = np.where(cdz == 0.0, 1.0, cdz)

        open_reg = cz >= mask_top - EPS
        tube_reg = ~open_reg

        # OPEN region: reach mask_top (down) or source_z (up)
        o_down = open_reg & (cdz < -1e-15)
        td = np.where(o_down, (mask_top - cz) / cdz_safe, inf)
        xhd = cx + cdx * td
        into = opening_mask(xhd)
        set_trans = o_down & (td > EPS) & into
        set_mask = o_down & (td > EPS) & ~into
        tmin = np.where(set_trans, td, tmin); etype = np.where(set_trans, 4, etype)
        upd = set_mask & (td < tmin); tmin = np.where(upd, td, tmin); etype = np.where(upd, 1, etype)
        o_up = open_reg & (cdz > 1e-15)
        tu = np.where(o_up, (source_z - cz) / cdz_safe, inf)
        upd = o_up & (tu > EPS) & (tu < tmin); tmin = np.where(upd, tu, tmin); etype = np.where(upd, 3, etype)

        # TUBE region: floor (down), walls, tube-top (up -> transition)
        tf = np.where(tube_reg & (cdz < -1e-15), (floor_z - cz) / cdz_safe, inf)
        upd = tube_reg & (tf > EPS) & (tf < tmin); tmin = np.where(upd, tf, tmin); etype = np.where(upd, 0, etype)
        t0 = _t_to_wall_value(cx, cdx, ox0, Wc)
        t1 = _t_to_wall_value(cx, cdx, ox1, Wc)
        tw = np.minimum(t0, t1); is0 = t0 <= t1
        upd = tube_reg & (tw > EPS) & (tw < tmin)
        tmin = np.where(upd, tw, tmin); etype = np.where(upd, 2, etype)
        wall_is_ox0 = np.where(upd, is0, wall_is_ox0)
        tt = np.where(tube_reg & (cdz > 1e-15), (mask_top - cz) / cdz_safe, inf)
        upd = tube_reg & (tt > EPS) & (tt < tmin); tmin = np.where(upd, tt, tmin); etype = np.where(upd, 4, etype)

        dead = ~np.isfinite(tmin)
        step = np.where(np.isfinite(tmin), tmin, 0.0)
        cx = cx + cdx * step
        cz = cz + cdz * step
        trans = (etype == 4) & ~dead
        cz = np.where(trans, cz + np.sign(cdz) * 1e-7, cz)

        is_floor = (etype == 0) & ~dead
        is_mask = (etype == 1) & ~dead
        is_wall = (etype == 2) & ~dead
        is_esc = (etype == 3) & ~dead
        floor_incident += float(is_floor.sum())

        react = rng.random(n) < s
        stop = is_esc | dead | (react & (is_floor | is_mask | is_wall))
        ndx = cdx.copy(); ndz = cdz.copy()
        for surf in (is_floor & ~react, is_mask & ~react):
            if surf.any():
                a, b = _cosine_inplane("z", +1.0, int(surf.sum()), rng)
                ndx[surf] = a; ndz[surf] = b
        rw = is_wall & ~react
        if rw.any():
            sign = np.where(wall_is_ox0[rw], +1.0, -1.0)
            a, b = _cosine_inplane("x", 1.0, int(rw.sum()), rng)
            ndx[rw] = np.abs(a) * sign; ndz[rw] = b

        x[ax] = cx; z[ax] = cz; dirx[ax] = ndx; dirz[ax] = ndz
        active[ax] = ~stop

    return 2.0 * floor_incident / N
