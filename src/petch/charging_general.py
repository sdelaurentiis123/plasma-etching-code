"""General 2-D feature charging: charged particles in a self-consistent field, ANY geometry.

This is the geometry-agnostic engine the edge-array/trench solvers should have been. There are NO
named features ("edge line", "neighbour"): the input is just a material-tagged grid

    mat[ix, iz] in {GAS, INSULATOR, CONDUCTOR}   plus  cid[ix, iz] = conductor component id

and the physics is the same for a trench, a hole cross-section, a line-and-space edge array, or a
real device: launch ions (directional, from the sheath) and electrons (ISOTROPIC, Lambert flux
through the plasma boundary), push both through the Laplace field, deposit charge where they land,
let each insulator cell float to local current balance and each connected conductor float as one
equipotential, re-solve the field, iterate to steady state.

Electrons are TRACED (not a view-factor shortcut) so the electrostatic focusing HG requires --
electrons pulled into the positive trench, reaching the floor more than pure geometry allows -- is
captured self-consistently. The only thing that made the old electron source wrong was launching a
DOWN-going hemisphere at the mouth; here electrons enter isotropically from the top boundary.

Kushner MCFPM lineage; this is the W5 shared-interface direction from ROBUST_PHYSICS_MODEL_PLAN.
"""
from __future__ import annotations

import math

import numpy as np

try:
    from numba import njit, prange
except Exception:  # pragma: no cover
    njit = None
    prange = range

GAS = 0
INSULATOR = 1
CONDUCTOR = 2


def _trace_general_py(Ex, Ez, solid, x0, z0, vx0, vz0, q, nx, nz, max_steps):
    """Push particles through the field until they enter a solid cell; return the hit CELL.

    No feature labels -- just (hit_ix, hit_iz) and the impact kinematics. Leapfrog with adaptive dt
    (<=0.45 cell/step), periodic x. z<0.5 = escaped back to plasma (survivor)."""
    n = x0.shape[0]
    hit_ix = np.full(n, -1, np.int64)
    hit_iz = np.full(n, -1, np.int64)
    impact_E = np.zeros(n)
    hit_vx = np.zeros(n)
    survivor = np.zeros(n, np.uint8)
    xmax = float(nx)
    for p in prange(n):
        x = x0[p]; z = z0[p]; vx = vx0[p]; vz = vz0[p]
        alive = True
        for _ in range(max_steps):
            ix = int(x); iz = int(z)
            if ix < 0:
                ix = 0
            elif ix > nx - 2:
                ix = nx - 2
            if iz < 0:
                iz = 0
            elif iz > nz - 2:
                iz = nz - 2
            fx = Ex[ix, iz]; fz = Ez[ix, iz]
            ax = q * fx * 0.5; az = q * fz * 0.5
            avx = vx if vx >= 0.0 else -vx
            avz = vz if vz >= 0.0 else -vz
            vmax = avx if avx >= avz else avz
            if vmax < 0.8:
                vmax = 0.8
            dt_v = 0.45 / vmax
            field = (fx * fx + fz * fz) ** 0.5
            if field < 1.0e-9:
                field = 1.0e-9
            dt_e = 0.3 / (field ** 0.5)
            dt = dt_v if dt_v <= dt_e else dt_e
            vx_half = vx + 0.5 * ax * dt
            vz_half = vz + 0.5 * az * dt
            xa = x + vx_half * dt
            za = z + vz_half * dt
            ix2 = int(xa); iz2 = int(za)
            if ix2 < 0:
                ix2 = 0
            elif ix2 > nx - 2:
                ix2 = nx - 2
            if iz2 < 0:
                iz2 = 0
            elif iz2 > nz - 2:
                iz2 = nz - 2
            vx = vx_half + 0.25 * q * Ex[ix2, iz2] * dt
            vz = vz_half + 0.25 * q * Ez[ix2, iz2] * dt
            x = xa % xmax
            z = za
            if z < 0.5:
                alive = False
                break
            ixh = int(x); izh = int(z)
            if ixh < 0:
                ixh = 0
            elif ixh > nx - 1:
                ixh = nx - 1
            if izh < 0:
                izh = 0
            elif izh > nz - 1:
                izh = nz - 1
            if solid[ixh, izh]:
                hit_ix[p] = ixh
                hit_iz[p] = izh
                impact_E[p] = vx * vx + vz * vz
                hit_vx[p] = vx
                alive = False
                break
        if alive:
            survivor[p] = 1
    return hit_ix, hit_iz, impact_E, hit_vx, survivor


_trace_general = (njit(cache=True, parallel=True, fastmath=True)(_trace_general_py)
                  if njit is not None else _trace_general_py)


def _connected_conductor_ids(mat):
    """Label connected components of CONDUCTOR cells (4-connectivity). Each component is one
    floating equipotential. Returns an int grid (0 = not a conductor, >=1 = component id)."""
    try:
        from scipy.ndimage import label
        cid, ncomp = label(mat == CONDUCTOR)
        return cid.astype(np.int64), int(ncomp)
    except Exception:
        # tiny fallback: every conductor cell its own id (correct if conductors are pre-separated)
        cid = np.zeros(mat.shape, np.int64)
        idx = np.where(mat == CONDUCTOR)
        cid[idx] = np.arange(1, idx[0].size + 1)
        return cid, int(idx[0].size)


def sample_ions(n, rng, mouth, nx, V_dc, V_rf, iadf_hwhm_deg):
    """Directional ions from the sheath: energy over the RF cycle, narrow angular cone, launched
    at the mouth plane moving down. (Ion optics unchanged from the validated model.)"""
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    E0 = np.maximum(V_dc + V_rf * np.sin(phi), 0.5)
    sig = np.deg2rad(iadf_hwhm_deg) / 1.1774 * np.sqrt(V_dc / E0)
    th = rng.normal(0.0, sig, n)
    vx = np.sqrt(E0) * np.sin(th)
    vz = np.sqrt(E0) * np.abs(np.cos(th))
    x = rng.uniform(0.0, float(nx - 1), n)
    z = np.full(n, max(1.0, float(mouth) - 0.5))
    return x, z, vx, vz


def sample_electrons(n, rng, nx, Te, cos_power=1.0):
    """Isotropic thermal electrons entering through the plasma boundary (top, z~=1): thermal energy,
    Lambert (cos^p) flux distribution so they arrive at ALL angles including near-horizontal. The
    field then focuses them into positive features (the piece a down-going mouth-launch threw away).
    Launched at the top so they traverse the field region and reach outward-facing walls."""
    E0 = rng.gamma(2.0, Te, n)
    u = rng.uniform(0.0, 1.0, n)
    ct = (1.0 - u) ** (1.0 / (cos_power + 2.0))
    st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
    az = rng.uniform(0.0, 2.0 * np.pi, n)
    th = np.arctan2(st * np.cos(az), ct)
    vx = np.sqrt(E0) * np.sin(th)
    vz = np.sqrt(E0) * np.abs(np.cos(th))
    x = rng.uniform(0.0, float(nx - 1), n)
    z = np.full(n, 1.0)
    return x, z, vx, vz


def solve_charging(mat, mouth, Te=4.0, V_dc=37.0, V_rf=30.0, iadf_hwhm_deg=4.3,
                   n_per_iter=6000, n_iter=200, relax=None, seed=0,
                   insul_vguard=None, verbose=False):
    """Steady-state feature charging for ANY material grid `mat` (GAS/INSULATOR/CONDUCTOR).

    mat: (nx, nz) int grid. z=0 is the plasma boundary (Dirichlet 0), z increases into the wafer.
    Returns V (potential), per-cell insulator potential, and per-conductor equipotentials."""
    if _trace_general is None:
        raise RuntimeError("numba unavailable")
    rng = np.random.default_rng(seed)
    nx, nz = mat.shape
    solid = mat != GAS
    insul = mat == INSULATOR
    cid, ncomp = _connected_conductor_ids(mat)
    if relax is None:
        relax = 2.0 * Te
    if insul_vguard is None:
        insul_vguard = V_dc + V_rf

    Vs = np.zeros((nx, nz))          # per-cell insulator potential
    Vc = np.zeros(ncomp + 1)         # per-conductor-component equipotential (index 0 unused)

    ii, jj = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    red = ((ii + jj) % 2 == 0)
    inside = ~solid
    inside[:, 0] = False

    def apply_bc(V):
        V[:, 0] = 0.0
        V[insul] = Vs[insul]
        if ncomp > 0:
            V[cid > 0] = Vc[cid[cid > 0]]
        return V

    def laplace(V, sweeps=200, omega=1.88):
        for _ in range(sweeps):
            V = apply_bc(V)
            for color in (red, ~red):
                xm = np.empty_like(V); xp = np.empty_like(V)
                xm[1:, :] = V[:-1, :]; xm[0, :] = V[0, :]
                xp[:-1, :] = V[1:, :]; xp[-1, :] = V[1:, :][-1:]
                avg = np.zeros_like(V)
                avg[:, 1:-1] = 0.25 * (xm[:, 1:-1] + xp[:, 1:-1] + V[:, 2:] + V[:, :-2])
                m = inside & color
                V[m] = (1.0 - omega) * V[m] + omega * avg[m]
        return apply_bc(V)

    def trace(kind, n, Ex, Ez):
        if kind == "ion":
            x, z, vx, vz = sample_ions(n, rng, mouth, nx, V_dc, V_rf, iadf_hwhm_deg)
            q = 1.0
        else:
            x, z, vx, vz = sample_electrons(n, rng, nx, Te)
            q = -1.0
        hix, hiz, E, _, surv = _trace_general(Ex, Ez, solid, x, z, vx, vz, q, nx, nz, 40 * nz)
        counts = np.zeros((nx, nz))
        m = hix >= 0
        if m.any():
            np.add.at(counts, (hix[m], hiz[m]), 1.0)
        return counts, float(surv.mean())

    V = np.zeros((nx, nz))
    hist = []
    for it in range(n_iter):
        V = laplace(V)
        Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
        ci, si = trace("ion", n_per_iter, Ex, Ez)
        ce, se = trace("electron", n_per_iter, Ex, Ez)
        net = ci - ce
        anneal = max(1.0 / (1.0 + it / 25.0), 0.25)
        scale = anneal * relax / n_per_iter * nx
        Vs[insul] += scale * net[insul]
        Vs[insul] = np.clip(Vs[insul], -insul_vguard, V_dc + V_rf)
        for c in range(1, ncomp + 1):
            m = cid == c
            area = max(int(m.sum()), 1)
            Vc[c] += scale * float(net[m].sum()) / area
            Vc[c] = float(np.clip(Vc[c], -3.0 * Te, V_dc + V_rf))
        hist.append((float(si), float(se)))
        if verbose and it % 20 == 0:
            print(f"  it{it}: surv_i/e={si:.3f}/{se:.3f} "
                  f"Vc={np.round(Vc[1:], 1) if ncomp else '-'} "
                  f"Vs[max]={Vs.max():.1f}", flush=True)

    V = laplace(V, sweeps=180)
    return dict(V=V, Vs=Vs, Vc=Vc[1:], ncomp=ncomp, cid=cid,
                surv_ion=hist[-1][0], surv_electron=hist[-1][1])
