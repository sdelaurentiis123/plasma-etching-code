"""Minimal 3D feature-scale etch loop (Phase 1).

Reuses the dimension-agnostic chemistry (`chemistry.surface_rate`) and mirrors the 2D pipeline
in 3D: 3D level set -> marching-cubes triangle mesh -> Warp ray-traced flux (ions + neutrals,
with diffuse re-emission) -> chemistry -> 3D upwind advection -> reinit.

The flux kernel is a `wp.kernel` using wp.Mesh + wp.mesh_query_ray (BVH; RT cores on a GPU).
CPU-first; set DEVICE='cuda' on an NVIDIA box and the identical kernel runs on RT cores.
"""
import os
import warnings
import numpy as np
import skfmm
from skimage import measure
from scipy.spatial import cKDTree
from scipy.stats import qmc, norm
import warp as wp

from .params import PAR, DEFAULT_FLAGS
from .chemistry import surface_rate

wp.init()
DEVICE = os.environ.get("PETCH_DEVICE", "cpu")   # set PETCH_DEVICE=cuda on a GPU box

# ----------------------------- 3D geometry / level set -----------------------------
def make_trench_3d(Lx, Ly, Lz, dx, trench_width, mask_th, sub_top, hole=False):
    """3D level set (phi>0 solid). Trench (invariant in y) or circular hole (hole=True)."""
    nx, ny, nz = int(round(Lx/dx)), int(round(Ly/dx)), int(round(Lz/dx))
    xs, ys, zs = np.arange(nx)*dx, np.arange(ny)*dx, np.arange(nz)*dx
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    solid = Z < sub_top
    mask_band = (Z >= sub_top) & (Z < sub_top + mask_th)
    if hole:
        opening = (X - Lx/2)**2 + (Y - Ly/2)**2 < (trench_width/2)**2
    else:
        opening = np.abs(X - Lx/2) < trench_width/2
    mask = mask_band & (~opening)
    solid = solid | mask
    phi = skfmm.distance(np.where(solid, 1.0, -1.0), dx=dx)
    return dict(xs=xs, ys=ys, zs=zs, dx=dx, phi=phi, mask=mask,
                Lx=Lx, Ly=Ly, Lz=Lz, sub_top=sub_top)


def extract_mesh_3d(phi, dx):
    """Marching cubes -> triangle mesh. Returns verts (physical), faces, centroids, areas."""
    verts, faces, _, _ = measure.marching_cubes(phi, level=0.0, spacing=(dx, dx, dx))
    v = verts[faces]                                   # (F,3,3)
    centroids = v.mean(axis=1)
    cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    return verts.astype(np.float32), faces.astype(np.int32), centroids, areas


# ----------------------------- Warp ray-traced flux -----------------------------
@wp.kernel
def _trace3d(mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), dir0: wp.array(dtype=wp.vec3),
             sticking: float, n_reemit: int, seed: int,
             flux: wp.array(dtype=float), angacc: wp.array(dtype=float)):
    p = wp.tid()
    state = wp.rand_init(seed, p)
    o = origin[p]
    d = dir0[p]
    w = float(1.0)
    for bounce in range(n_reemit + 1):
        q = wp.mesh_query_ray(mesh, o, d, 1.0e6)
        if not q.result:
            break
        n = q.normal
        if wp.dot(n, d) > 0.0:                          # orient into the gas (toward the ray)
            n = -n
        cosang = -wp.dot(d, n)
        if cosang < 0.0:
            cosang = 0.0
        u = wp.randf(state)
        if u < sticking or bounce == n_reemit:
            wp.atomic_add(flux, q.face, w)
            wp.atomic_add(angacc, q.face, w * cosang)
            break
        else:                                           # diffuse 3D-cosine re-emission about n
            hit = o + q.t * d
            a = wp.vec3(1.0, 0.0, 0.0)
            if wp.abs(n[0]) > 0.9:
                a = wp.vec3(0.0, 1.0, 0.0)
            t = wp.normalize(wp.cross(a, n))
            b = wp.cross(n, t)
            r1 = wp.randf(state)
            r2 = wp.randf(state)
            ct = wp.sqrt(r1)
            st = wp.sqrt(1.0 - r1)
            phi = 6.2831853 * r2
            d = st * wp.cos(phi) * t + st * wp.sin(phi) * b + ct * n
            o = hit + 1.0e-4 * n


# (the old _launch_dirs helper was replaced by _source3d, which also does QMC)
def _source3d(kind, n, Lx, Ly, z_src, sigma, sampling, rng, sd):
    """Source launch (position + direction). 'pseudo' (PoC) or 'sobol' (QMC over the 4D source).

    Each particle uses (x, y, angle1, angle2); QMC over those four dims drives variance ~1/N.
    """
    if sampling == "sobol":
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*power of 2.*")
            u = qmc.Sobol(d=4, scramble=True, seed=sd).random(n)
        ox, oy = u[:, 0] * Lx, u[:, 1] * Ly
        u2 = np.clip(u[:, 2], 1e-9, 1 - 1e-9); u3 = u[:, 3]
        if kind == 'ion':
            ax = norm.ppf(u2) * sigma; ay = norm.ppf(u3) * sigma
            d = np.stack([np.sin(ax), np.sin(ay), -np.cos(ax) * np.cos(ay)], axis=1)
        else:
            ct = np.sqrt(u2); st = np.sqrt(1 - u2); ph = 2 * np.pi * u3
            d = np.stack([st * np.cos(ph), st * np.sin(ph), -ct], axis=1)
    else:
        ox, oy = rng.uniform(0, Lx, n), rng.uniform(0, Ly, n)
        if kind == 'ion':
            ax = rng.normal(0, sigma, n); ay = rng.normal(0, sigma, n)
            d = np.stack([np.sin(ax), np.sin(ay), -np.cos(ax) * np.cos(ay)], axis=1)
        else:
            ct = np.sqrt(rng.uniform(0, 1, n)); st = np.sqrt(1 - ct**2)
            ph = rng.uniform(0, 2 * np.pi, n)
            d = np.stack([st * np.cos(ph), st * np.sin(ph), -ct], axis=1)
    origin = np.stack([ox, oy, np.full(n, z_src)], axis=1).astype(np.float32)
    dirs = (d / np.linalg.norm(d, axis=1, keepdims=True)).astype(np.float32)
    return origin, dirs


def mc_flux_3d(mesh, verts, faces, areas, geo, par, n_ion=20000, n_neu=20000, seed=0,
               sampling="pseudo"):
    """Per-face normalized flux multipliers (m_i,m_F,m_O) + mean ion cos, via Warp ray tracing."""
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']
    A_src = Lx * Ly

    def run(kind, n, sticking, n_re, sd):
        origin, dirs = _source3d(kind, n, Lx, Ly, z_src, par['ion_ang_sigma'], sampling, rng, sd)
        flux = wp.zeros(F, dtype=float, device=DEVICE)
        ang = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d, dim=n, device=DEVICE,
                  inputs=[mesh.id, wp.array(origin, dtype=wp.vec3, device=DEVICE),
                          wp.array(dirs, dtype=wp.vec3, device=DEVICE),
                          float(sticking), int(n_re), int(sd), flux, ang])
        return flux.numpy(), ang.numpy()

    fi, ai = run('ion', n_ion, 1.0, 0, seed * 9 + 1)
    fF, _ = run('neutral', n_neu, par['s_F'], 12, seed * 9 + 2)
    fO, _ = run('neutral', n_neu, par['s_O'], 12, seed * 9 + 3)

    A = np.maximum(areas, 0.3 * np.median(areas))
    base_i = n_ion / A_src
    base_n = n_neu / A_src
    m_i = np.clip((fi / A) / base_i, 0.0, 1.5)
    m_F = np.clip((fF / A) / base_n, 0.0, 4.0)
    m_O = np.clip((fO / A) / base_n, 0.0, 4.0)
    cos_i = np.where(fi > 0, ai / np.maximum(fi, 1e-9), 0.0)
    return m_i, m_F, m_O, cos_i


# ----------------------------- 3D advection -----------------------------
def advect_3d(phi, Fspeed, dx, dt):
    """phi_t + F|grad phi| = 0, first-order upwind Godunov in 3D."""
    g = np.zeros_like(phi)
    for ax in range(3):
        dm = np.zeros_like(phi); dp = np.zeros_like(phi)
        sl_m = [slice(None)] * 3; sl_mm = [slice(None)] * 3
        sl_m[ax] = slice(1, None); sl_mm[ax] = slice(0, -1)
        diff = (phi[tuple(sl_m)] - phi[tuple(sl_mm)]) / dx
        dm[tuple(sl_m)] = diff
        dp[tuple(sl_mm)] = diff
        g += np.maximum(dm, 0)**2 + np.minimum(dp, 0)**2
    return phi - dt * Fspeed * np.sqrt(g)


def extend_velocity_3d(V, centroids, geo, band):
    """Nearest-face velocity extension into the 3D narrow band (CPU KDTree fallback)."""
    phi, dx = geo['phi'], geo['dx']
    Fs = np.zeros_like(phi)
    bm = np.abs(phi) < band
    ii, jj, kk = np.where(bm)
    pts = np.stack([geo['xs'][ii], geo['ys'][jj], geo['zs'][kk]], axis=1)
    if len(centroids) == 0:
        return Fs
    _, idx = cKDTree(centroids).query(pts)
    Fs[ii, jj, kk] = V[idx]
    return Fs


@wp.kernel
def _extend_kernel(mesh: wp.uint64, pts: wp.array(dtype=wp.vec3),
                   Vface: wp.array(dtype=float), out: wp.array(dtype=float)):
    i = wp.tid()
    q = wp.mesh_query_point_no_sign(mesh, pts[i], 1.0e6)
    if q.result:
        out[i] = Vface[q.face]


def extend_velocity_gpu(mesh, V, geo, band):
    """Nearest-face velocity extension via wp.mesh_query_point (BVH; GPU-resident)."""
    phi = geo['phi']
    Fs = np.zeros_like(phi)
    ii, jj, kk = np.where(np.abs(phi) < band)
    if len(ii) == 0:
        return Fs
    pts = np.stack([geo['xs'][ii], geo['ys'][jj], geo['zs'][kk]], axis=1).astype(np.float32)
    out = wp.zeros(len(pts), dtype=float, device=DEVICE)
    wp.launch(_extend_kernel, dim=len(pts), device=DEVICE,
              inputs=[mesh.id, wp.array(pts, dtype=wp.vec3, device=DEVICE),
                      wp.array(V.astype(np.float32), dtype=float, device=DEVICE), out])
    Fs[ii, jj, kk] = out.numpy()
    return Fs


@wp.kernel
def _reinit_iter(phi: wp.array3d(dtype=float), phi0: wp.array3d(dtype=float),
                 out: wp.array3d(dtype=float), inv_dx: float, dtau: float):
    """One Godunov reinitialization sweep of phi_tau + sgn(phi0)(|grad phi|-1)=0."""
    i, j, k = wp.tid()
    nx = phi.shape[0]; ny = phi.shape[1]; nz = phi.shape[2]
    c = phi[i, j, k]
    s0 = phi0[i, j, k]
    if wp.abs(s0) < 1.2 / inv_dx:        # interface band -> freeze (preserve phi=0 contour)
        out[i, j, k] = c
        return
    Dxm = (c - phi[wp.max(i - 1, 0), j, k]) * inv_dx
    Dxp = (phi[wp.min(i + 1, nx - 1), j, k] - c) * inv_dx
    Dym = (c - phi[i, wp.max(j - 1, 0), k]) * inv_dx
    Dyp = (phi[i, wp.min(j + 1, ny - 1), k] - c) * inv_dx
    Dzm = (c - phi[i, j, wp.max(k - 1, 0)]) * inv_dx
    Dzp = (phi[i, j, wp.min(k + 1, nz - 1)] - c) * inv_dx
    sgn = s0 / wp.sqrt(s0 * s0 + 0.25 / (inv_dx * inv_dx))     # smoothed sign (eps = dx/2)
    if sgn > 0.0:
        gx = wp.max(wp.max(Dxm, 0.0), -wp.min(Dxp, 0.0))
        gy = wp.max(wp.max(Dym, 0.0), -wp.min(Dyp, 0.0))
        gz = wp.max(wp.max(Dzm, 0.0), -wp.min(Dzp, 0.0))
    else:
        gx = wp.max(-wp.min(Dxm, 0.0), wp.max(Dxp, 0.0))
        gy = wp.max(-wp.min(Dym, 0.0), wp.max(Dyp, 0.0))
        gz = wp.max(-wp.min(Dzm, 0.0), wp.max(Dzp, 0.0))
    G = wp.sqrt(gx * gx + gy * gy + gz * gz)
    out[i, j, k] = c - dtau * sgn * (G - 1.0)


def reinit_gpu(phi_np, dx, n_iter=20):
    """GPU iterative SDF reinit (Warp), interface-band frozen to pin the phi=0 contour.

    EXPERIMENTAL / approximate: this first-order Godunov scheme plateaus at |grad phi| ~ 0.96
    (~4% error, ~0.2*dx SDF error), which compounds over an etch -> depth differs from skfmm.
    Use only where max throughput matters and small SDF error is acceptable. Production accuracy
    needs a WENO reinit + Russo-Smereka subcell interface fix. Default reinit_method stays 'skfmm'.
    """
    a = wp.array(phi_np.astype(np.float32), dtype=float, device=DEVICE)
    phi0 = wp.array(phi_np.astype(np.float32), dtype=float, device=DEVICE)
    b = wp.zeros_like(a)
    dtau = 0.5 * dx
    inv = 1.0 / dx
    for _ in range(n_iter):
        wp.launch(_reinit_iter, dim=phi_np.shape, device=DEVICE, inputs=[a, phi0, b, inv, dtau])
        a, b = b, a
    return a.numpy().astype(np.float64)


def faces_in_mask(centroids, geo, mask_th, trench_width, hole=False):
    """Mark faces whose centroid lies in the (un-etched) mask material."""
    x, y, z = centroids[:, 0], centroids[:, 1], centroids[:, 2]
    in_band = (z >= geo['sub_top']) & (z < geo['sub_top'] + mask_th)
    if hole:
        opening = (x - geo['Lx']/2)**2 + (y - geo['Ly']/2)**2 < (trench_width/2)**2
    else:
        opening = np.abs(x - geo['Lx']/2) < trench_width/2
    return in_band & (~opening)


# ----------------------------- driver -----------------------------
def run_etch_3d(Lx=10.0, Ly=4.0, Lz=14.0, dx=0.4, trench_width=4.0, mask_th=2.0,
                sub_top=10.0, t_end=2.0, n_steps=20, hole=False, par=None, flags=None,
                n_ion=20000, n_neu=20000, reinit_every=1, extend="gpu",
                reinit_method="skfmm", verbose=True):
    if par is None:
        par = PAR
    if flags is None:
        flags = DEFAULT_FLAGS
    geo = make_trench_3d(Lx, Ly, Lz, dx, trench_width, mask_th, sub_top, hole=hole)
    mask_phi = geo['phi'].copy()
    dt = t_end / n_steps
    band = 4 * dx
    import time
    timings = dict(flux=0.0, extend=0.0, reinit=0.0, total=0.0)
    t0 = time.time()
    for step in range(n_steps):
        verts, faces, centroids, areas = extract_mesh_3d(geo['phi'], dx)
        mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEVICE),
                       indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEVICE))
        tf = time.time()
        m_i, m_F, m_O, cos_i = mc_flux_3d(mesh, verts, faces, areas, geo, par,
                                          n_ion=n_ion, n_neu=n_neu, seed=step,
                                          sampling=getattr(flags, "sampling", "pseudo"))
        timings['flux'] += time.time() - tf
        is_mask = faces_in_mask(centroids, geo, mask_th, trench_width, hole=hole)
        V = surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=flags)
        te = time.time()
        Fs = (extend_velocity_gpu(mesh, V, geo, band) if extend == "gpu"
              else extend_velocity_3d(V, centroids, geo, band))
        timings['extend'] += time.time() - te
        Vmax = max(V.max(), 1e-6)
        nsub = max(1, min(int(np.ceil(Vmax * dt / (0.4 * dx))), 40))
        for _ in range(nsub):
            geo['phi'] = advect_3d(geo['phi'], Fs, dx, dt / nsub)
            geo['phi'][geo['mask']] = mask_phi[geo['mask']]
        # reinit_every>1 (lazy reinit) is faster but DRIFTS the result: |grad phi| deviates from
        # 1 between reinits and advect multiplies F*|grad phi|. Safe only with a proper extension
        # velocity (grad F . grad phi = 0). Keep =1 for fidelity.
        if (step + 1) % reinit_every == 0 or step == n_steps - 1:
            tr = time.time()
            geo['phi'] = (reinit_gpu(geo['phi'], dx) if reinit_method == "gpu"
                          else skfmm.distance(geo['phi'], dx=dx))
            timings['reinit'] += time.time() - tr
        if verbose and step % 5 == 0:
            depth = _depth3d(geo)
            print(f"  step {step:3d}/{n_steps}  faces {len(faces):5d}  depth ~ {depth:5.2f}  Vmax {Vmax:.3f}")
    timings['total'] = time.time() - t0
    geo['timings'] = timings
    return geo


def _depth3d(geo, half=1.0):
    """Center etch depth: robust over a small central region (median floor of central columns).

    For each central column, the floor is the deepest grid cell that is gas AND connected to the
    open top (no isolated-pocket spikes). Returns sub_top - median(floor_z) over the region.
    """
    phi, xs, ys, zs = geo['phi'], geo['xs'], geo['ys'], geo['zs']
    ic = np.where(np.abs(xs - geo['Lx']/2) <= half)[0]
    jc = np.where(np.abs(ys - geo['Ly']/2) <= half)[0]
    floors = []
    for i in ic:
        for j in jc:
            col = phi[i, j, :] < 0                       # gas mask, bottom..top
            # floor = lowest z index of the gas run that reaches the top
            k = len(col) - 1
            if not col[k]:
                continue
            while k > 0 and col[k - 1]:
                k -= 1
            floors.append(zs[k])
    if not floors:
        return 0.0
    return float(geo['sub_top'] - np.median(floors))


def center_depth_3d(geo, half=1.0):
    return _depth3d(geo, half)
