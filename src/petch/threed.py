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
             specular: int, cos_thr: float, eta: float,
             flux: wp.array(dtype=float), angacc: wp.array(dtype=float)):
    """Trace a species. specular=0: diffuse re-emission (neutrals). specular=1: ions deposit
    yield on hit and SPECULAR-REFLECT at grazing incidence (cosang < cos_thr) with energy
    retention eta -- feeds the bottom corners of deep HARC features (contributor #4)."""
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
        if specular == 1:                               # ION: deposit + grazing specular reflect
            if cosang >= cos_thr or bounce == n_reemit:
                wp.atomic_add(flux, q.face, w)
                wp.atomic_add(angacc, q.face, w * cosang)
                break
            R = (cos_thr - cosang) / cos_thr            # reflected fraction (0 at thr, 1 grazing)
            dep = w * (1.0 - R)
            wp.atomic_add(flux, q.face, dep)
            wp.atomic_add(angacc, q.face, dep * cosang)
            w = w * R * eta                             # reflected weight (energy loss)
            hit = o + q.t * d
            d = d - 2.0 * wp.dot(d, n) * n              # specular reflection
            o = hit + 1.0e-4 * n
        else:                                           # NEUTRAL: stick or diffuse re-emit
            u = wp.randf(state)
            if u < sticking or bounce == n_reemit:
                wp.atomic_add(flux, q.face, w)
                wp.atomic_add(angacc, q.face, w * cosang)
                break
            hit = o + q.t * d                           # diffuse 3D-cosine re-emission about n
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
               sampling="pseudo", ion_reflection=False):
    """Per-face normalized flux multipliers (m_i,m_F,m_O) + mean ion cos, via Warp ray tracing."""
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']
    A_src = Lx * Ly
    COS_THR, ETA = 0.34, 0.8                            # ion reflect onset (~70 deg), energy keep

    def run(kind, n, sticking, n_re, sd, specular=0):
        origin, dirs = _source3d(kind, n, Lx, Ly, z_src, par['ion_ang_sigma'], sampling, rng, sd)
        flux = wp.zeros(F, dtype=float, device=DEVICE)
        ang = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d, dim=n, device=DEVICE,
                  inputs=[mesh.id, wp.array(origin, dtype=wp.vec3, device=DEVICE),
                          wp.array(dirs, dtype=wp.vec3, device=DEVICE),
                          float(sticking), int(n_re), int(sd),
                          int(specular), float(COS_THR), float(ETA), flux, ang])
        return flux.numpy(), ang.numpy()

    if ion_reflection:
        fi, ai = run('ion', n_ion, 1.0, 3, seed * 9 + 1, specular=1)
    else:
        fi, ai = run('ion', n_ion, 1.0, 0, seed * 9 + 1, specular=0)
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


# ------------------- Langmuir coverage-dependent neutral transport (3D ARDE fix) -------------------
@wp.kernel
def _trace3d_cov(mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), dir0: wp.array(dtype=wp.vec3),
                 bare: wp.array(dtype=float), beta: float, n_reemit: int, seed: int,
                 flux: wp.array(dtype=float)):
    """Neutral trace with COVERAGE-DEPENDENT sticking S_eff = bare*beta (Langmuir: stick only on
    bare sites). Records ARRIVING flux on every hit (like ViennaPS). On saturated walls (low bare)
    radicals reflect -> penetrate to the under-fed floor -> keeps deep-hole floors fed (flat ARDE).
    """
    p = wp.tid()
    state = wp.rand_init(seed, p)
    o = origin[p]
    d = dir0[p]
    for bounce in range(n_reemit + 1):
        q = wp.mesh_query_ray(mesh, o, d, 1.0e6)
        if not q.result:
            break
        n = q.normal
        if wp.dot(n, d) > 0.0:
            n = -n
        wp.atomic_add(flux, q.face, 1.0)                # arriving flux (every hit)
        S = bare[q.face] * beta                         # sticks on bare sites only
        if wp.randf(state) < S or bounce == n_reemit:
            break
        hit = o + q.t * d                               # diffuse 3D-cosine re-emission
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


def _belen_coverages(m_i, m_F, m_O, cos_i, par, flags):
    """Belen steady-state coverages theta_F, theta_O from per-face fluxes (for the fixed point)."""
    from .chemistry import _yields, angular_factors
    Yie, Ysp, Yp = _yields(par)
    mode = "cosine" if flags is None else getattr(flags, "yield_angular", "cosine")
    f_ie, _ = angular_factors(cos_i, par, mode)
    Fi = par['ionFlux'] * m_i
    eps = 1e-9
    GY_ie = Yie * f_ie * Fi
    GY_p = Yp * f_ie * Fi
    Gb_E = par['Fflux'] * m_F * par.get('cal_F', 1.0) + eps
    Gb_P = par['Oflux'] * m_O + eps
    a = (par['k_sigma'] + 2.0 * GY_ie) / Gb_E
    b = (par['beta_sigma'] + GY_p) / Gb_P
    thF = 1.0 / (1.0 + a * (1.0 + 1.0 / (b + eps)))
    thO = 1.0 / (1.0 + b * (1.0 + 1.0 / (a + eps)))
    return thF, thO


def mc_flux_3d_coupled(mesh, verts, faces, areas, geo, par, n_ion=20000, n_neu=20000, seed=0,
                       sampling="pseudo", flags=None, n_fp=4):
    """Coverage-coupled flux: ions once, then a flux<->coverage fixed point with coverage-dependent
    neutral sticking. Returns per-face (m_i, m_F, m_O, cos_i) normalized to arriving open-field=1."""
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']
    A_src = Lx * Ly
    A = np.maximum(areas, 0.3 * np.median(areas))
    betaE = par.get('betaE', 0.7); betaO = par.get('betaO', 1.0)

    # ions (unchanged; not coverage-dependent)
    oi, di = _source3d('ion', n_ion, Lx, Ly, z_src, par['ion_ang_sigma'], sampling, rng, seed * 9 + 1)
    fi = wp.zeros(F, dtype=float, device=DEVICE); ai = wp.zeros(F, dtype=float, device=DEVICE)
    spec = 1 if (flags is not None and getattr(flags, "ion_reflection", False)) else 0
    nre_i = 3 if spec == 1 else 0
    wp.launch(_trace3d, dim=n_ion, device=DEVICE,
              inputs=[mesh.id, wp.array(oi, dtype=wp.vec3, device=DEVICE),
                      wp.array(di, dtype=wp.vec3, device=DEVICE),
                      1.0, int(nre_i), int(seed * 9 + 1), int(spec), 0.34, 0.8, fi, ai])
    fi = fi.numpy(); ai = ai.numpy()
    m_i = np.clip((fi / A) / (n_ion / A_src), 0.0, 1.5)
    cos_i = np.where(fi > 0, ai / np.maximum(fi, 1e-9), 0.0)

    def neutral(beta, bare, sd):
        o, d = _source3d('neutral', n_neu, Lx, Ly, z_src, par['ion_ang_sigma'], sampling, rng, sd)
        fl = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d_cov, dim=n_neu, device=DEVICE,
                  inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEVICE),
                          wp.array(d, dtype=wp.vec3, device=DEVICE),
                          wp.array(bare.astype(np.float32), dtype=float, device=DEVICE),
                          float(beta), 24, int(sd), fl])
        return np.clip((fl.numpy() / A) / (n_neu / A_src), 0.0, 8.0)

    bare = np.ones(F)
    m_F = m_O = np.zeros(F)
    for it in range(n_fp):
        m_F = neutral(betaE, bare, seed * 9 + 2 + 2 * it)
        m_O = neutral(betaO, bare, seed * 9 + 3 + 2 * it)
        thF, thO = _belen_coverages(m_i, m_F, m_O, cos_i, par, flags)
        bare = np.clip(1.0 - thF - thO, 0.0, 1.0)
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
    # phi0 face neighbors (clamped)
    o_xm = phi0[wp.max(i - 1, 0), j, k]; o_xp = phi0[wp.min(i + 1, nx - 1), j, k]
    o_ym = phi0[i, wp.max(j - 1, 0), k]; o_yp = phi0[i, wp.min(j + 1, ny - 1), k]
    o_zm = phi0[i, j, wp.max(k - 1, 0)]; o_zp = phi0[i, j, wp.min(k + 1, nz - 1)]
    # interface cell? (phi0 sign change with a face neighbor)
    interface = (s0 * o_xm < 0.0) or (s0 * o_xp < 0.0) or (s0 * o_ym < 0.0) \
        or (s0 * o_yp < 0.0) or (s0 * o_zm < 0.0) or (s0 * o_zp < 0.0)
    if interface:
        # Russo-Smereka subcell: drive phi toward the true distance D = phi0/|grad phi0|, which
        # pins the phi=0 contour at its exact sub-cell position (no interface drift).
        gx0 = (o_xp - o_xm) * 0.5 * inv_dx
        gy0 = (o_yp - o_ym) * 0.5 * inv_dx
        gz0 = (o_zp - o_zm) * 0.5 * inv_dx
        grad0 = wp.sqrt(gx0 * gx0 + gy0 * gy0 + gz0 * gz0) + 1.0e-9
        D = s0 / grad0
        sgn0 = 1.0
        if s0 < 0.0:
            sgn0 = -1.0
        out[i, j, k] = c - dtau * inv_dx * (sgn0 * wp.abs(c) - D)
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


def reinit_gpu(phi_np, dx, n_iter=24):
    """GPU SDF reinit (Warp) with the Russo-Smereka subcell interface fix.

    Interface cells are driven toward the exact sub-cell distance phi0/|grad phi0|, which PINS the
    phi=0 contour (no drift); other cells use the Godunov |grad phi|=1 sweep. Result: |grad phi| ~
    1.00 near the front and depth parity with skfmm is exact. The far field need not fully converge
    (the etch only uses the near-front band). This makes a fully GPU-resident loop accurate; on a
    GPU it is far faster than skfmm-on-CPU. (CPU default stays 'skfmm' since skfmm is fast there.)
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
    # belen chemistry uses ViennaPS sticking; keep transport re-emission consistent with it
    mc_par = par
    if getattr(flags, "chemistry", "langmuir") == "belen":
        mc_par = dict(par); mc_par['s_F'] = par['betaE']; mc_par['s_O'] = par['betaO']
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
        if getattr(flags, "coverage_sticking", False):   # Langmuir coverage-dependent sticking
            m_i, m_F, m_O, cos_i = mc_flux_3d_coupled(mesh, verts, faces, areas, geo, par,
                                                      n_ion=n_ion, n_neu=n_neu, seed=step,
                                                      sampling=getattr(flags, "sampling", "pseudo"),
                                                      flags=flags)
        else:
            m_i, m_F, m_O, cos_i = mc_flux_3d(mesh, verts, faces, areas, geo, mc_par,
                                              n_ion=n_ion, n_neu=n_neu, seed=step,
                                              sampling=getattr(flags, "sampling", "pseudo"),
                                              ion_reflection=getattr(flags, "ion_reflection", False))
        timings['flux'] += time.time() - tf
        is_mask = faces_in_mask(centroids, geo, mask_th, trench_width, hole=hole)
        V = surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=flags)
        V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)   # guard against blowup
        te = time.time()
        Fs = (extend_velocity_gpu(mesh, V, geo, band) if extend == "gpu"
              else extend_velocity_3d(V, centroids, geo, band))
        timings['extend'] += time.time() - te
        vmx = float(np.max(V)) if V.size else 0.0
        Vmax = max(vmx if np.isfinite(vmx) else 0.0, 1e-6)
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
