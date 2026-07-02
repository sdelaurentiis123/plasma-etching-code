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
try:                                             # optional: cupy for the GPU edge-adjacency sort
    import cupy as _cp
    _HAS_CUPY = True
except Exception:
    _cp = None
    _HAS_CUPY = False

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
    # Smooth analytic SDF seed (phi>0 solid), sub-cell accurate -- NOT a binary +-1 carve. The binary
    # carve quantizes the slot wall / substrate top to the nearest cell -> the effective trench width
    # shifts with dx -> grid-sensitive ARDE (the "feels unstable" drift). Here we place the interfaces at
    # their TRUE geometric positions via a CSG min/max field (substrate half-space UNION masked slab),
    # then skfmm cleans the distance while preserving the sub-cell zero contour (matches ViennaLS's smooth
    # geometry, grid-independent). r = lateral distance from the feature axis.
    r = np.sqrt((X - Lx/2)**2 + (Y - Ly/2)**2) if hole else np.abs(X - Lx/2)
    d_sub = sub_top - Z                                   # >0 inside substrate (z < sub_top)
    d_slab = np.minimum(Z - sub_top, sub_top + mask_th - Z)   # >0 inside the mask z-band
    d_mask = np.minimum(d_slab, r - trench_width / 2.0)   # mask solid = in z-band AND outside the slot
    phi_analytic = np.maximum(d_sub, d_mask)              # total solid = substrate UNION masked slab
    phi = skfmm.distance(phi_analytic, dx=dx)             # clean SDF, sub-cell zero contour preserved
    return dict(xs=xs, ys=ys, zs=zs, dx=dx, phi=phi, mask=mask,
                Lx=Lx, Ly=Ly, Lz=Lz, sub_top=sub_top,
                trench_width=trench_width, hole=hole)


def extract_mesh_3d(phi, dx):
    """Marching cubes -> triangle mesh. Returns verts (physical), faces, centroids, areas."""
    verts, faces, _, _ = measure.marching_cubes(phi, level=0.0, spacing=(dx, dx, dx))
    v = verts[faces]                                   # (F,3,3)
    centroids = v.mean(axis=1)
    cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    return verts.astype(np.float32), faces.astype(np.int32), centroids, areas


_mc_cache = {}


def extract_mesh_3d_gpu(phi, dx):
    """GPU marching cubes via Warp's built-in MarchingCubes (~24x faster than skimage incl. readback;
    mesh was the #2 loop cost at ~26%). Slightly different triangulation in ambiguous cubes than
    skimage's Lewiner (~96% of the faces) -> validate WITHIN-NOISE on depth before trusting. The
    MarchingCubes context is cached per grid shape (its kernels compile once)."""
    nx, ny, nz = phi.shape
    try:
        mc = _mc_cache.get((nx, ny, nz))
        if mc is None:
            mv = mt = nx * ny * nz // 2 + 1000
            mc = wp.MarchingCubes(nx, ny, nz, max_verts=mv, max_tris=mt, device=DEVICE)
            _mc_cache[(nx, ny, nz)] = mc
        fwp = wp.array(phi.astype(np.float32), dtype=float, device=DEVICE)
        mc.surface(fwp, 0.0)
        verts = mc.verts.numpy().astype(np.float64) * dx      # node-index coords -> physical
        faces = mc.indices.numpy().reshape(-1, 3)
        if len(faces) == 0:
            raise RuntimeError("GPU MC produced no faces")
    except Exception:
        # Warp MarchingCubes can fail on thin dims / certain Warp versions -> CPU skimage (reliable).
        return extract_mesh_3d(phi, dx)
    v = verts[faces]
    centroids = v.mean(axis=1)
    cross = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    return verts.astype(np.float32), faces.astype(np.int32), centroids, areas


def _edge_adjacency(faces):
    """Edge-neighbor face pairs (faces sharing an edge). The edge sort was the #1 flux host cost
    (~66ms/step on 180k faces). On GPU we offload the whole thing to cupy (argsort ~0.26ms); else a
    single-key int64 argsort on host. Pair order is irrelevant (the smooth scatters symmetrically)."""
    F = len(faces)
    if _HAS_CUPY and DEVICE == 'cuda':
        fg = _cp.asarray(faces)
        e = _cp.sort(fg[:, [[0, 1], [1, 2], [0, 2]]].reshape(-1, 2), axis=1)
        nv = int(faces.max()) + 1
        key = e[:, 0].astype(_cp.int64) * nv + e[:, 1]
        fid = _cp.repeat(_cp.arange(F), 3)
        order = _cp.argsort(key)
        key_s = key[order]; fid_s = fid[order]
        idx = _cp.where(key_s[:-1] == key_s[1:])[0]
        return _cp.asnumpy(_cp.stack([fid_s[idx], fid_s[idx + 1]], axis=1))
    e = np.sort(faces[:, [[0, 1], [1, 2], [0, 2]]].reshape(-1, 2), axis=1)   # (3F,2) sorted edges
    nv = int(faces.max()) + 1
    key = e[:, 0].astype(np.int64) * nv + e[:, 1]                            # encode each edge as one int64
    fid = np.repeat(np.arange(F), 3)
    order = np.argsort(key, kind='stable')
    key_s = key[order]; fid_s = fid[order]
    idx = np.where(key_s[:-1] == key_s[1:])[0]
    return np.stack([fid_s[idx], fid_s[idx + 1]], axis=1)



def smooth_flux(flux, normals, pairs, n_iter=1, alpha=1.0):
    """ViennaPS 1-neighbor normal-weighted flux smoothing (rayTraceDisk.hpp::smoothFlux): each face
    averaged with edge neighbors weighted by max(0, dot(n_i, n_j)), self-weight 1. Laterally diffuses
    flux into narrow HARC floors -> flattens small-feature ARDE.

    `alpha` (0..1) is a STRENGTH knob to calibrate our (dense edge-mesh) smoothing to ViennaPS's
    (spatial disk-radius) neighborhood: result = (1-alpha)*raw + alpha*fully_smoothed. alpha=1 is the
    default full smoothing; lower alpha = milder. (alpha=1 + n_iter=1 = exact ViennaPS default.)"""
    if len(pairs) == 0 or n_iter <= 0 or alpha <= 0.0:
        return flux
    raw = flux.astype(np.float64)
    out = raw.copy()
    i, j = pairs[:, 0], pairs[:, 1]
    w = np.maximum(0.0, np.sum(normals[i] * normals[j], axis=1))
    for _ in range(n_iter):
        num = out.copy(); den = np.ones(len(out))
        np.add.at(num, i, out[j] * w); np.add.at(den, i, w)
        np.add.at(num, j, out[i] * w); np.add.at(den, j, w)
        out = num / den
    return out if alpha >= 1.0 else (1.0 - alpha) * raw + alpha * out


@wp.kernel
def _smooth_scatter(src: wp.array(dtype=float), pi: wp.array(dtype=int), pj: wp.array(dtype=int),
                    w: wp.array(dtype=float), num: wp.array(dtype=float), den: wp.array(dtype=float)):
    """One edge-pair scatter for smooth_flux: num[i]+=src[j]*w, den[i]+=w (and symmetric j<-i).
    src is the PREVIOUS iterate (read-only); num is pre-seeded with src (self-weight 1), den with 1."""
    p = wp.tid()
    i = pi[p]; j = pj[p]; wij = w[p]
    wp.atomic_add(num, i, src[j] * wij)
    wp.atomic_add(den, i, wij)
    wp.atomic_add(num, j, src[i] * wij)
    wp.atomic_add(den, j, wij)


# ---- device-resident flux kernels (keep flux on the GPU through normalize+smooth -> 1 readback/step) ----
@wp.kernel
def _norm_clip(fl: wp.array(dtype=float), A: wp.array(dtype=float), scale: float,
               lo: float, hi: float, out: wp.array(dtype=float)):
    """out = clamp((fl/A)*scale, lo, hi) -- the host np.clip((fl/A)/(n/A_src),...) on device."""
    i = wp.tid()
    out[i] = wp.clamp((fl[i] / A[i]) * scale, lo, hi)


@wp.kernel
def _div_k(num: wp.array(dtype=float), den: wp.array(dtype=float), out: wp.array(dtype=float)):
    i = wp.tid()
    out[i] = num[i] / wp.max(den[i], 1.0e-12)


@wp.kernel
def _cos_k(fl: wp.array(dtype=float), ang: wp.array(dtype=float), out: wp.array(dtype=float)):
    i = wp.tid()
    if fl[i] > 0.0:
        out[i] = ang[i] / fl[i]
    else:
        out[i] = 0.0


@wp.kernel
def _lerp_k(raw: wp.array(dtype=float), sm: wp.array(dtype=float), a: float, out: wp.array(dtype=float)):
    """Device blend for the smoothing strength knob: out = (1-a)*raw + a*smoothed."""
    i = wp.tid()
    out[i] = (1.0 - a) * raw[i] + a * sm[i]


def smooth_flux_dev(src_wp, prep, alpha=1.0):
    """Device-resident smooth: src in/out are Warp arrays (no host round-trip). num seeded = src
    (self-weight), den = 1; one scatter, one divide -- all on device. `alpha` (0..1) is the strength
    knob: out = (1-alpha)*raw + alpha*fully_smoothed (alpha=1 = full smoothing; lower = milder). This
    matches the host smooth_flux/smooth_flux_gpu blend so flux_smooth_alpha works on the GPU path too
    (it was silently ignored here before -> the knob had no effect on the device-resident neutral flux)."""
    F = src_wp.shape[0]
    if alpha <= 0.0:
        return wp.clone(src_wp)
    num = wp.clone(src_wp)
    den = wp.full(F, 1.0, dtype=float, device=DEVICE)
    wp.launch(_smooth_scatter, dim=prep['pi'].shape[0], device=DEVICE,
              inputs=[src_wp, prep['pi'], prep['pj'], prep['ww'], num, den])
    out = wp.zeros(F, dtype=float, device=DEVICE)
    wp.launch(_div_k, dim=F, device=DEVICE, inputs=[num, den, out])
    if alpha < 1.0:                                    # blend raw <-> fully-smoothed (strength knob)
        res = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_lerp_k, dim=F, device=DEVICE, inputs=[src_wp, out, float(alpha), res])
        return res
    return out


def _smooth_prep_gpu(normals, pairs):
    """Build the device edge-pair index/weight arrays for smooth_flux_gpu ONCE per step. The 3 smooth
    calls/step (m_i, m_F, m_O) share the same mesh -> same pairs/weights; caching avoids rebuilding the
    weights (host np.sum over ~100k pairs) and re-uploading pi/pj/ww on every call."""
    i = np.ascontiguousarray(pairs[:, 0]).astype(np.int32)
    j = np.ascontiguousarray(pairs[:, 1]).astype(np.int32)
    w = np.maximum(0.0, np.sum(normals[i] * normals[j], axis=1)).astype(np.float32)
    return dict(pi=wp.array(i, dtype=int, device=DEVICE), pj=wp.array(j, dtype=int, device=DEVICE),
                ww=wp.array(w, dtype=float, device=DEVICE))


@wp.kernel
def _face_normals_k(verts: wp.array(dtype=wp.vec3), faces: wp.array(dtype=int), out: wp.array(dtype=wp.vec3)):
    """Unit face normal (cross of two edges) on the GPU -- the host np.cross over the face array was
    ~41ms/step on 180k faces."""
    i = wp.tid()
    a = verts[faces[3 * i + 0]]; b = verts[faces[3 * i + 1]]; c = verts[faces[3 * i + 2]]
    out[i] = wp.normalize(wp.cross(b - a, c - a))


@wp.kernel
def _edge_weights_k(fn: wp.array(dtype=wp.vec3), pi: wp.array(dtype=int), pj: wp.array(dtype=int),
                    w: wp.array(dtype=float)):
    """Smooth weight w[p] = max(0, dot(n_i, n_j)) on the GPU."""
    p = wp.tid()
    w[p] = wp.max(0.0, wp.dot(fn[pi[p]], fn[pj[p]]))


def _smooth_prep_dev(verts, faces, pairs):
    """Device smooth prep: face normals + edge weights computed ON THE GPU (host cross + np.sum were
    ~55ms/step on 180k faces). Only the edge-pair indices (from the host adjacency sort) are uploaded."""
    F = len(faces)
    vw = wp.array(verts.astype(np.float32), dtype=wp.vec3, device=DEVICE)
    fwp = wp.array(np.ascontiguousarray(faces).astype(np.int32).reshape(-1), dtype=int, device=DEVICE)
    fn = wp.zeros(F, dtype=wp.vec3, device=DEVICE)
    wp.launch(_face_normals_k, dim=F, device=DEVICE, inputs=[vw, fwp, fn])
    pi = wp.array(np.ascontiguousarray(pairs[:, 0]).astype(np.int32), dtype=int, device=DEVICE)
    pj = wp.array(np.ascontiguousarray(pairs[:, 1]).astype(np.int32), dtype=int, device=DEVICE)
    ww = wp.zeros(len(pairs), dtype=float, device=DEVICE)
    wp.launch(_edge_weights_k, dim=len(pairs), device=DEVICE, inputs=[fn, pi, pj, ww])
    return dict(pi=pi, pj=pj, ww=ww, n=len(pairs))


def smooth_flux_gpu(flux, normals, pairs, n_iter=1, alpha=1.0, prep=None):
    """GPU smooth_flux: the np.add.at edge scatter (the dominant flux-line host cost after FSM+warm) on a
    Warp atomic-add kernel. Matches numpy smooth_flux to ~1e-5 (float32 atomics). `prep` reuses cached
    device pairs/weights from _smooth_prep_gpu (built once/step) -> no per-call rebuild + re-upload."""
    if len(pairs) == 0 or n_iter <= 0 or alpha <= 0.0:
        return flux
    raw = flux.astype(np.float64)
    F = len(flux)
    if prep is None:
        prep = _smooth_prep_gpu(normals, pairs)
    pi, pj, ww = prep['pi'], prep['pj'], prep['ww']
    out = raw.copy()
    for _ in range(n_iter):
        of = out.astype(np.float32)
        src = wp.array(of, dtype=float, device=DEVICE)
        num = wp.array(of.copy(), dtype=float, device=DEVICE)               # seed num=src (self-weight)
        den = wp.array(np.ones(F, np.float32), dtype=float, device=DEVICE)  # seed den=1
        wp.launch(_smooth_scatter, dim=pi.shape[0], device=DEVICE, inputs=[src, pi, pj, ww, num, den])
        out = (num.numpy() / np.maximum(den.numpy(), 1e-12)).astype(np.float64)
    return out if alpha >= 1.0 else (1.0 - alpha) * raw + alpha * out


# ----------------------------- Warp ray-traced flux -----------------------------
@wp.struct
class _BCRay:
    o: wp.vec3
    d: wp.vec3


@wp.func
def _apply_bc(mesh: wp.uint64, o: wp.vec3, d: wp.vec3, lx: float, ly: float, lz: float, active: int):
    """Lateral boundary conditions for TRENCHES, matching ViennaPS MakeTrench(periodicBoundary=True):
    PERIODIC in x & y, infinite bottom, open top=sky. When a ray misses the surface and would exit a
    lateral face before the open top, it WRAPS to the opposite face with the same direction (it sees the
    periodic image of the feature). This keeps wide-angle cosine neutrals in the domain (open-field flux
    ~1.0, not the ~0.78 of an open boundary that lets them escape) WITHOUT the over-feeding of reflective
    BCs (which mirror flux straight back down a thin trench -> floor over-fed -> too-gentle ARDE; that was
    a regression). Holes pass active=0 (no lateral wrap; they use the full 3D reflective domain)."""
    r = _BCRay()
    r.o = o
    r.d = d
    if active == 0:
        return r
    oo = o
    dd = d
    for _w in range(256):
        q = wp.mesh_query_ray(mesh, oo, dd, 1.0e6)
        if q.result:
            r.o = oo
            r.d = dd
            return r                                     # hits a face from here
        tx = 1.0e30
        if dd[0] > 1.0e-9:
            tx = (lx - oo[0]) / dd[0]
        if dd[0] < -1.0e-9:
            tx = (0.0 - oo[0]) / dd[0]
        ty = 1.0e30
        if dd[1] > 1.0e-9:
            ty = (ly - oo[1]) / dd[1]
        if dd[1] < -1.0e-9:
            ty = (0.0 - oo[1]) / dd[1]
        tz = 1.0e30
        if dd[2] > 1.0e-9:
            tz = (lz - oo[2]) / dd[2]                     # reaches the open top first -> real escape
        # earliest lateral exit before the top? PERIODIC wrap (ViennaPS MakeTrench periodicBoundary=True):
        # the ray re-enters the opposite face with the SAME direction (it sees the periodic image of the
        # feature). NOT reflect -- reflective lateral BCs over-feed a thin trench floor -> too-gentle ARDE.
        if tx < tz and tx < ty and tx > 1.0e-7:          # exit x -> wrap to opposite x face
            hy = oo[1] + tx * dd[1]
            hz = oo[2] + tx * dd[2]
            nx = 1.0e-4                                  # exited +x -> re-enter at x=0
            if dd[0] < 0.0:
                nx = lx - 1.0e-4                          # exited -x -> re-enter at x=lx
            oo = wp.vec3(nx, hy, hz)
            dd = wp.vec3(dd[0], dd[1], dd[2])            # same direction (periodic, not mirrored)
        elif ty < tz and ty < 1.0e29 and ty > 1.0e-7:    # exit y -> wrap to opposite y face
            hx = oo[0] + ty * dd[0]
            hz = oo[2] + ty * dd[2]
            ny = 1.0e-4                                  # exited +y -> re-enter at y=0
            if dd[1] < 0.0:
                ny = ly - 1.0e-4                          # exited -y -> re-enter at y=ly
            oo = wp.vec3(hx, ny, hz)
            dd = wp.vec3(dd[0], dd[1], dd[2])            # same direction (periodic)
        else:
            r.o = oo
            r.d = dd
            return r                                     # escapes via top -> sky
    r.o = oo
    r.d = dd
    return r


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


@wp.kernel
def _gen_source(kind: int, Lx: float, Ly: float, z_src: float, sigma: float, seed: int,
                origin: wp.array(dtype=wp.vec3), dirs: wp.array(dtype=wp.vec3)):
    """Generate one source ray (origin on the source plane + launch direction) ON THE GPU -- kills the
    host Sobol gen + H2D upload (the ~4.3ms/call source-gen cost). PSEUDOrandom (not QMC); matches the
    `sampling='pseudo'` transform exactly: ion = near-vertical Gaussian tilt, neutral = cosine hemisphere."""
    p = wp.tid()
    st8 = wp.rand_init(seed, p)
    origin[p] = wp.vec3(wp.randf(st8) * Lx, wp.randf(st8) * Ly, z_src)
    if kind == 1:                                  # ion: near-vertical Gaussian
        ax = wp.randn(st8) * sigma
        ay = wp.randn(st8) * sigma
        d = wp.vec3(wp.sin(ax), wp.sin(ay), -wp.cos(ax) * wp.cos(ay))
    else:                                          # neutral: cosine-weighted hemisphere about -z
        r1 = wp.randf(st8); r2 = wp.randf(st8)
        ct = wp.sqrt(r1); s = wp.sqrt(1.0 - r1); ph = 6.2831853 * r2
        d = wp.vec3(s * wp.cos(ph), s * wp.sin(ph), -ct)
    dirs[p] = wp.normalize(d)


def gen_source_gpu(kind, n, Lx, Ly, z_src, sigma, seed):
    """On-device source rays -> (origin, dirs) wp.vec3 arrays, no host gen / no upload. Returns Warp
    arrays directly (pass straight to the trace kernels)."""
    o = wp.zeros(n, dtype=wp.vec3, device=DEVICE)
    d = wp.zeros(n, dtype=wp.vec3, device=DEVICE)
    wp.launch(_gen_source, dim=n, device=DEVICE,
              inputs=[1 if kind == 'ion' else 0, float(Lx), float(Ly), float(z_src), float(sigma),
                      int(seed), o, d])
    return o, d


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

    A = np.maximum(areas, 1.0e-9)              # actual triangle area (ViennaPS uses area; only a tiny degeneracy floor)
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


@wp.kernel
def _trace3d_cov_rr(mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), dir0: wp.array(dtype=wp.vec3),
                    bare: wp.array(dtype=float), beta: float, seed: int, flux: wp.array(dtype=float),
                    lx: float, ly: float, lz: float, periodic: int):
    """EXACT ViennaPS neutral transport: weighted ray + Russian roulette, NO fixed bounce cap.
    Verbatim from ViennaRay rayTraceKernel.hpp + psPlasmaEtching.hpp:
      - ray weight W starts at 1.0; DEPOSIT W into arriving flux on every hit (not +1.0);
      - S_eff = bare*beta (bare = 1 - eCov - pCov);  W *= (1 - S_eff);  stop if W <= 0;
      - rejection control: if W >= 0.1 always continue; else kill w.p. (1 - W/0.3), survivors W = 0.3;
      - diffuse cosine re-emission. The fixed cap (512) is only a safety bound (roulette/escape ends it).
    This is the UNBIASED estimator ViennaPS uses; a fixed bounce cap biases the deep-floor flux low."""
    p = wp.tid()
    state = wp.rand_init(seed, p)
    o = origin[p]
    d = dir0[p]
    w = float(1.0)
    for bounce in range(1024):
        bc = _apply_bc(mesh, o, d, lx, ly, lz, periodic)
        o = bc.o
        d = bc.d
        q = wp.mesh_query_ray(mesh, o, d, 1.0e6)
        if not q.result:
            break
        n = q.normal
        if wp.dot(n, d) > 0.0:
            n = -n
        wp.atomic_add(flux, q.face, w)                  # deposit ray WEIGHT (arriving flux)
        S = bare[q.face] * beta
        w = w - w * S                                   # W *= (1 - S_eff)
        if w <= 1.0e-6:
            break
        if w < 0.1:                                     # Russian roulette (rejection control)
            if wp.randf(state) < (1.0 - w / 0.3):
                break
            w = 0.3
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


@wp.kernel
def _trace3d_ion_yield(mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), dir0: wp.array(dtype=wp.vec3),
                       seed: int, meanE: float, sigmaE: float,
                       A_sp: float, Eth_sp: float, B_sp: float,
                       A_ie: float, Eth_ie: float, A_p: float, Eth_p: float,
                       inflect: float, minang: float, n_l: float, A_en: float, minE: float,
                       thrmin: float, thrmax: float,
                       lx: float, ly: float, lz: float, periodic: int,
                       sp: wp.array(dtype=float), ie: wp.array(dtype=float), ip: wp.array(dtype=float),
                       arr: wp.array(dtype=float), ang: wp.array(dtype=float)):
    """EXACT ViennaPS SF6O2 ion (psPlasmaEtching.hpp PlasmaEtchingIon): near-vertical ions have
    sticking=0 and REFLECT with ViennaRay's coned-cosine law, losing energy each bounce until the
    energy drops below min(Eth). On EVERY hit they deposit the yield-WEIGHTED sputter, ion-enhanced
    (Si) and passivation (O) rates (surfaceCollision) -- not just an arrival count. Grazing sidewall
    hits keep their energy and glance forward, so the ion funnels down a deep trench and feeds the
    floor: this is the deep-AR term petch's first-hit-stick ion model was missing (steep ARDE)."""
    p = wp.tid()
    state = wp.rand_init(seed, p)
    o = origin[p]
    d = dir0[p]
    E = meanE + sigmaE * wp.randn(state)
    if E < 0.0:
        E = 0.0
    w = float(1.0)                                                 # ray WEIGHT (ViennaPS: deposited yield *= w)
    PI_2 = 1.5707963
    sEsp = wp.sqrt(Eth_sp); sEie = wp.sqrt(Eth_ie); sEp = wp.sqrt(Eth_p)
    for bounce in range(1024):
        bc = _apply_bc(mesh, o, d, lx, ly, lz, periodic)
        o = bc.o
        d = bc.d
        q = wp.mesh_query_ray(mesh, o, d, 1.0e6)
        if not q.result:
            break
        n = q.normal
        if wp.dot(n, d) > 0.0:
            n = -n
        cosT = -wp.dot(d, n)
        if cosT < 0.0:
            cosT = 0.0
        if cosT > 1.0:
            cosT = 1.0
        incAng = wp.acos(cosT)
        sqrtE = wp.sqrt(E)
        f_sp = (1.0 + B_sp * (1.0 - cosT * cosT)) * cosT           # ViennaPS sputter angular factor
        if f_sp < 0.0:
            f_sp = 0.0
        f_ie = float(1.0)                                          # ion-enhanced: flat to 60deg then -> 0
        if cosT < 0.5:
            f_ie = 3.0 - 6.0 * incAng / 3.14159265
            if f_ie < 0.0:
                f_ie = 0.0
        wp.atomic_add(sp, q.face, A_sp * wp.max(sqrtE - sEsp, 0.0) * f_sp * w)   # deposit yield * rayWeight
        wp.atomic_add(ie, q.face, A_ie * wp.max(sqrtE - sEie, 0.0) * f_ie * w)
        wp.atomic_add(ip, q.face, A_p * wp.max(sqrtE - sEp, 0.0) * f_ie * w)
        wp.atomic_add(arr, q.face, w)                             # weighted arrival (m_i, diagnostics)
        wp.atomic_add(ang, q.face, cosT * w)
        # ViennaPS angle-dependent ion sticking: ABSORB near-normal (incAng<thetaRMin, e.g. the floor),
        # REFLECT at grazing (incAng>thetaRMin, e.g. sidewalls). This stops ions bouncing off the floor
        # onto the walls -> keeps wall passivation right -> deep neutral floor not starved.
        stick = float(1.0)
        if incAng > thrmin:
            rr = (incAng - thrmin) / (thrmax - thrmin)
            if rr < 0.0:
                rr = 0.0
            if rr > 1.0:
                rr = 1.0
            stick = 1.0 - rr
        if stick >= 1.0:                                          # near-normal -> fully absorbed
            break
        if incAng >= inflect:                                      # ViennaPS updateEnergy (reflected energy)
            Eref = 1.0 - (1.0 - A_en) * (PI_2 - incAng) / (PI_2 - inflect)
        else:
            Eref = A_en * wp.pow(incAng / inflect, n_l)
        newE = Eref * E
        for _e in range(4):                                       # N(Eref*E, 0.1E) clamped to [0,E]
            cand = Eref * E + 0.1 * E * wp.randn(state)
            if cand >= 0.0 and cand <= E:
                newE = cand
                break
        if newE < 0.0:
            newE = 0.0
        if newE > E:
            newE = E
        if newE <= minE:                                          # ion thermalized -> stop (ViennaPS)
            break
        E = newE
        w = w - w * stick                                         # reduce weight by the angle sticking
        if w <= 1.0e-4:
            break
        hit = o + q.t * d
        axis = wp.normalize(d - 2.0 * wp.dot(d, n) * n)           # specular direction (cone axis)
        if axis[2] < -0.9999999:                                  # Frisvad orthonormal basis about axis
            tt = wp.vec3(0.0, -1.0, 0.0)
            bb = wp.vec3(-1.0, 0.0, 0.0)
        else:
            aa = 1.0 / (1.0 + axis[2])
            cc = -axis[0] * axis[1] * aa
            tt = wp.vec3(1.0 - axis[0] * axis[0] * aa, cc, -axis[0])
            bb = wp.vec3(cc, 1.0 - axis[1] * axis[1] * aa, -axis[1])
        maxcone = PI_2 - wp.min(incAng, minang)                  # cone narrows toward specular at grazing
        theta = float(0.0)
        for _t in range(64):                                     # ViennaRay ConedCosine accept-reject
            uu = wp.sqrt(wp.randf(state))
            ss = wp.sqrt(wp.max(1.0 - uu, 0.0))
            theta = maxcone * ss
            if wp.randf(state) * theta * uu <= wp.cos(PI_2 * ss) * wp.sin(theta):
                break
        phi = 6.2831853 * wp.randf(state)
        nd = wp.sin(theta) * (wp.cos(phi) * tt + wp.sin(phi) * bb) + wp.cos(theta) * axis
        nd = wp.normalize(nd)
        if wp.dot(nd, n) <= 0.0:                                  # keep it in the upper hemisphere
            nd = nd - 2.0 * wp.dot(nd, n) * n
        d = nd
        o = hit + 1.0e-4 * n


def _belen_coverages(m_i, m_F, m_O, cos_i, par, flags):
    """Belen steady-state coverages theta_F, theta_O from per-face fluxes (for the fixed point)."""
    from .chemistry import _yields, angular_factors
    Yie, Ysp, Yp = _yields(par)
    mode = "cosine" if flags is None else getattr(flags, "yield_angular", "cosine")
    f_ie, _ = angular_factors(cos_i, par, mode)
    Fi = par['ionFlux'] * m_i
    eps = 1e-9
    iy = par.get('_ion_yield')
    if iy is not None:                                # faithful ViennaPS ion model: yields deposited
        _sp, _ie, _ip = iy                            # per-ion energy + reflection, during ray tracing
        GY_ie = par['ionFlux'] * _ie                  # ion-enhanced etchant removal (Si)
        GY_p = par['ionFlux'] * _ip                   # ion-enhanced passivation removal (O)
    else:
        GY_ie = Yie * f_ie * Fi
        GY_p = Yp * f_ie * Fi
    Gb_E = par['Fflux'] * m_F + eps                   # arriving F flux (ViennaPS convention, no fudge)
    Gb_P = par['Oflux'] * m_O + eps
    a = (par['k_sigma'] + 2.0 * GY_ie) / Gb_E
    b = (par['beta_sigma'] + GY_p) / Gb_P
    thF = 1.0 / (1.0 + a * (1.0 + 1.0 / (b + eps)))
    thO = 1.0 / (1.0 + b * (1.0 + 1.0 / (a + eps)))
    return thF, thO


def mc_flux_3d_coupled(mesh, verts, faces, areas, geo, par, n_ion=20000, n_neu=20000, seed=0,
                       sampling="pseudo", flags=None, n_fp=4, bare_init=None):
    """Coverage-coupled flux: ions once, then a flux<->coverage fixed point with coverage-dependent
    neutral sticking. Returns per-face (m_i, m_F, m_O, cos_i, bare); `bare` (1-thetaF-thetaO) is the
    converged coverage so the caller can WARM-START the next step (same fixed point in 1-2 iters, not 4
    from cold -- accuracy-neutral, the front moves <1 cell/step). `bare_init` seeds it when warm."""
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']
    A_src = Lx * Ly
    A = np.maximum(areas, 1.0e-9)              # actual triangle area (ViennaPS uses area; only a tiny degeneracy floor)
    betaE = par.get('betaE', 0.7); betaO = par.get('betaO', 1.0)
    # ViennaPS 1-neighbor normal-weighted flux smoothing (default on; par['flux_smooth']=0 to disable).
    n_smooth = int(par.get('flux_smooth', 1))
    sm_alpha = float(par.get('flux_smooth_alpha', 1.0))   # strength knob (calibrate to ViennaPS)
    _cuda = DEVICE == 'cuda'                                       # GPU speedups auto-enable on CUDA
    use_gpu_smooth = par.get('flux_smooth_gpu', _cuda)
    use_gpu_src = par.get('gpu_source', _cuda)     # on-device ray gen (pseudorandom, kills host source-gen+upload)
    dev = bool(par.get('device_flux', _cuda)) and use_gpu_smooth   # device-resident normalize+smooth (kernels)

    def _src(kind, n, sd, sig):
        if use_gpu_src:
            return gen_source_gpu(kind, n, Lx, Ly, z_src, sig, sd)
        o, d = _source3d(kind, n, Lx, Ly, z_src, sig, sampling, rng, sd)
        return wp.array(o, dtype=wp.vec3, device=DEVICE), wp.array(d, dtype=wp.vec3, device=DEVICE)

    if n_smooth > 0:
        pairs = _edge_adjacency(faces)
        if use_gpu_smooth:                       # normals + weights ON GPU (host cross was ~41ms/step on big meshes)
            _prep = _smooth_prep_dev(verts, faces, pairs)
            _smooth = lambda x: smooth_flux_gpu(x, None, pairs, n_smooth, sm_alpha, prep=_prep)
        else:
            vv = verts[faces]
            fn = np.cross(vv[:, 1] - vv[:, 0], vv[:, 2] - vv[:, 0])
            fn = fn / (np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12)
            _smooth = lambda x: smooth_flux(x, fn, pairs, n_smooth, sm_alpha)
    else:
        pairs = None
        _smooth = lambda x: x

    A_wp = wp.array(A.astype(np.float32), dtype=float, device=DEVICE) if dev else None

    def _finish(fl_wp, scale, hi):
        """Device-resident: normalize (clamp((fl/A)*scale,0,hi)) + smooth, all on GPU -> host (1 readback)."""
        m_wp = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_norm_clip, dim=F, device=DEVICE, inputs=[fl_wp, A_wp, float(scale), 0.0, float(hi), m_wp])
        if n_smooth > 0:
            for _ in range(n_smooth):
                m_wp = smooth_flux_dev(m_wp, _prep, alpha=sm_alpha)
        return m_wp.numpy()

    # ions (unchanged; not coverage-dependent)
    oi_w, di_w = _src('ion', n_ion, seed * 9 + 1, par['ion_ang_sigma'])
    # Faithful ViennaPS ion: sticking=0 + coned-cosine reflection + per-ion energy, depositing the
    # yield-WEIGHTED sputter/ion-enhanced/passivation rates directly (the deep-AR funneling term).
    # Active for belen+ion_reflection; otherwise the legacy first-hit ion (with optional crude specular).
    faithful_ion = (flags is not None and getattr(flags, "ion_reflection", False)
                    and getattr(flags, "chemistry", "langmuir") == "belen")
    if faithful_ion:
        inflect = float(par.get('inflectAngle', 1.55334303)); minang = float(par.get('minAngle', 1.3962634))
        n_l = float(par.get('n_l', 10.0)); A_en = 1.0 / (1.0 + n_l * (1.5707963 / inflect - 1.0))
        minE = float(min(par['Eth_ie'], par['Eth_sp']))
        sp_w = wp.zeros(F, dtype=float, device=DEVICE); ie_w = wp.zeros(F, dtype=float, device=DEVICE)
        ip_w = wp.zeros(F, dtype=float, device=DEVICE)
        fi = wp.zeros(F, dtype=float, device=DEVICE); ai = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d_ion_yield, dim=n_ion, device=DEVICE,
                  inputs=[mesh.id, oi_w, di_w, int(seed * 9 + 1), float(par['Emean']), float(par['Esig']),
                          float(par['A_sp']), float(par['Eth_sp']), float(par['B_sp']),
                          float(par['A_ie']), float(par['Eth_ie']), float(par['A_p']), float(par['Eth_p']),
                          inflect, minang, n_l, float(A_en), minE,
                          float(par.get('thetaRMin', 1.2217305)), float(par.get('thetaRMax', 1.5707963)),
                          float(Lx), float(Ly), float(Lz), int(par.get('periodic_y', 0)),
                          sp_w, ie_w, ip_w, fi, ai])
        norm = (A_src / n_ion)
        def _ynorm(a):                                  # normalize like m_i then smooth (no clip on yields)
            m = (a.numpy() / A) * norm
            return _smooth(m) if n_smooth > 0 else m
        par['_ion_yield'] = (_ynorm(sp_w), _ynorm(ie_w), _ynorm(ip_w))
        fi = fi.numpy(); ai = ai.numpy()
        m_i = (fi / A) * norm                            # arrival flux (diagnostics; rate uses _ion_yield)
        cos_i = np.where(fi > 0, ai / np.maximum(fi, 1e-9), 0.0)
        if n_smooth > 0:
            m_i = _smooth(m_i)
    else:
        par.pop('_ion_yield', None)                      # belen falls back to analytic yields
        fi = wp.zeros(F, dtype=float, device=DEVICE); ai = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d, dim=n_ion, device=DEVICE,
                  inputs=[mesh.id, oi_w, di_w,
                          1.0, int(0), int(seed * 9 + 1), int(0), 0.34, 0.8, fi, ai])
        if dev:
            cw = wp.zeros(F, dtype=float, device=DEVICE)
            wp.launch(_cos_k, dim=F, device=DEVICE, inputs=[fi, ai, cw])
            cos_i = cw.numpy()
            m_i = _finish(fi, A_src / n_ion, 1.5)        # cos from RAW fi/ai; m_i normalized+smoothed
        else:
            fi = fi.numpy(); ai = ai.numpy()
            m_i = np.clip((fi / A) / (n_ion / A_src), 0.0, 1.5)
            cos_i = np.where(fi > 0, ai / np.maximum(fi, 1e-9), 0.0)
            if n_smooth > 0:
                m_i = _smooth(m_i)

    # SURFACE CHARGING (beyond ViennaPS): electrons arrive DIFFUSELY (cosine, unity sticking) so they
    # are more geometrically shadowed in HARC than the directional ions. On a floating/insulating
    # surface the net current -> 0, so the floor's effective ion flux is throttled toward the electron
    # arrival rate: f = 1 - alpha*(1 - Gamma_e/Gamma_i). Applied to m_i BEFORE the coverage fixed point
    # so charging reduces both the physical etch and the ion-enhanced coverage coupling. Hwang-Giapis 1997.
    alpha = float(par.get('charge_alpha', 0.0))
    if flags is not None and getattr(flags, "surface_charging", False) and alpha > 0.0:
        # electrons: PARTIALLY collimated by the sheath -> a moderate-spread Gaussian source (the 'ion'
        # source with e_ang_sigma), NOT a full cosine. Wider e_ang_sigma -> more HARC shadowing -> steeper
        # charging rolloff. Calibrated to Hwang-Giapis (cosine over-predicted; ions are near-vertical).
        oe, de = _src('ion', n_neu, seed * 9 + 7, par.get('e_ang_sigma', 0.5))
        fe = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d_cov_rr, dim=n_neu, device=DEVICE,                    # unity sticking: bare=1, beta=1
                  inputs=[mesh.id, oe, de, wp.array(np.ones(F, np.float32), dtype=float, device=DEVICE),
                          1.0, int(seed * 9 + 7), fe, float(Lx), float(Ly), float(Lz), int(par.get('periodic_y', 0))])
        m_e = np.clip((fe.numpy() / A) / (n_neu / A_src) * par.get('eFlux', 1.0), 0.0, 8.0)
        if n_smooth > 0:
            m_e = _smooth(m_e)
        # Anchor each flux to its OWN open-field exposure (the most-exposed faces, top/field, ~unity
        # shadow) so the charging factor is 1 at the field and rolls off only where electrons are MORE
        # shadowed than ions. shadow = flux / open-field flux; f = 1 - alpha*(1 - shadow_e/shadow_i).
        ref_i = np.percentile(m_i[m_i > 1e-6], 90) if (m_i > 1e-6).any() else 1.0
        ref_e = np.percentile(m_e[m_e > 1e-6], 90) if (m_e > 1e-6).any() else 1.0
        sh_e = np.clip(m_e / max(ref_e, 1e-6), 0.0, 1.0)
        sh_i = np.clip(m_i / max(ref_i, 1e-6), 0.0, 1.0)
        f_charge = 1.0 - alpha * np.clip(1.0 - sh_e / np.maximum(sh_i, 1e-3), 0.0, 1.0)
        m_i = m_i * f_charge

    def neutral(beta, bare, sd):
        o_w, d_w = _src('neutral', n_neu, sd, par['ion_ang_sigma'])
        fl = wp.zeros(F, dtype=float, device=DEVICE)
        # EXACT ViennaPS transport: weighted ray + Russian roulette, no fixed bounce cap. The old
        # fixed/adaptive cap biased the deep-floor flux low (truncation); roulette is the unbiased
        # estimator ViennaPS uses (set par['cov_transport']='cap' to fall back to the capped kernel).
        if par.get('cov_transport', 'rr') == 'rr':
            wp.launch(_trace3d_cov_rr, dim=n_neu, device=DEVICE,
                      inputs=[mesh.id, o_w, d_w,
                              wp.array(bare.astype(np.float32), dtype=float, device=DEVICE),
                              float(beta), int(sd), fl, float(Lx), float(Ly), float(Lz), int(par.get('periodic_y', 0))])
        else:
            nre = int(par.get('n_reemit_cov', np.clip(8.0 / max(beta, 0.02), 24, 200)))
            wp.launch(_trace3d_cov, dim=n_neu, device=DEVICE,
                      inputs=[mesh.id, o_w, d_w,
                              wp.array(bare.astype(np.float32), dtype=float, device=DEVICE),
                              float(beta), int(nre), int(sd), fl])
        if dev:
            return _finish(fl, A_src / n_neu, 8.0)     # device normalize+smooth, 1 readback
        m = np.clip((fl.numpy() / A) / (n_neu / A_src), 0.0, 8.0)
        return _smooth(m) if n_smooth > 0 else m

    # coverage fixed-point iters: warm-started steps (bare_init given) converge in 1; cold needs ~4.
    n_fp = int(par.get('n_fp', 1 if bare_init is not None else 4))
    bare = np.ones(F) if bare_init is None else np.clip(np.asarray(bare_init, float), 0.0, 1.0)
    m_F = m_O = np.zeros(F)
    for it in range(n_fp):
        m_F = neutral(betaE, bare, seed * 9 + 2 + 2 * it)
        m_O = neutral(betaO, bare, seed * 9 + 3 + 2 * it)
        thF, thO = _belen_coverages(m_i, m_F, m_O, cos_i, par, flags)
        bare = np.clip(1.0 - thF - thO, 0.0, 1.0)
    return m_i, m_F, m_O, cos_i, bare


@wp.kernel
def _trace_ff(mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), dir0: wp.array(dtype=wp.vec3),
              hitface: wp.array(dtype=int), lx: float, ly: float, lz: float, periodic: int):
    """Single-bounce form-factor ray: cast one diffuse ray from a face, record the FIRST face it hits
    (or -1 = escapes to the open source/sky). Reflective-x / periodic-y lateral BCs (ViennaPS default):
    rays exiting the lateral faces re-enter, so the deep floor sees the trench conductance not 'sky'."""
    p = wp.tid()
    bc = _apply_bc(mesh, origin[p], dir0[p], lx, ly, lz, periodic)
    q = wp.mesh_query_ray(mesh, bc.o, bc.d, 1.0e6)
    if q.result:
        hitface[p] = q.face
    else:
        hitface[p] = -1


def mc_flux_3d_radiosity(mesh, verts, faces, centroids, areas, geo, par, n_ion=20000,
                         n_ff=64, seed=0, flags=None, n_fp=4):
    """DETERMINISTIC RADIOSITY neutral flux (vs many-bounce MC). Build the face-to-face form-factor
    matrix once by single-bounce ray casting (cheap, MC), then solve the multi-bounce equilibrium
    arriving flux EXACTLY by a sparse linear iteration -- so the deep HARC floor gets its true Knudsen
    conductance flux with NO under-sampling (the MC's deep-floor starvation), AND it is far cheaper than
    deep multi-bounce MC. Ions stay MC (directional, single-bounce). Solves the de-Boer real-wafer ARDE
    gap (MC under-sampling) and is the >14x neutral-flux speedup that pushes past ViennaPS.

      Gamma_i = D_i + sum_j (1-s_j) A[i,j] Gamma_j,   A[i,j]=F_{j->i},  D_i = sky view of i,  m_F=Gamma."""
    import scipy.sparse as sp
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']; A_src = Lx * Ly
    A = np.maximum(areas, 1.0e-9)              # actual triangle area (ViennaPS uses area; only a tiny degeneracy floor)
    betaE = par.get('betaE', 0.7); betaO = par.get('betaO', 1.0)
    fn = _gas_normals(verts, faces, centroids, geo)     # into-gas (essential for face emission)

    # ions: shared deterministic-path launch (faithful ViennaPS reflection when the config asks)
    oi, di = _source3d('ion', n_ion, Lx, Ly, z_src, par['ion_ang_sigma'], 'sobol', rng, seed * 9 + 1)
    oi_w = wp.array(oi, dtype=wp.vec3, device=DEVICE); di_w = wp.array(di, dtype=wp.vec3, device=DEVICE)
    m_i, cos_i = _ions_deterministic(mesh, F, A, A_src, Lx, Ly, Lz, n_ion, par, flags,
                                     seed * 9 + 1, oi_w, di_w)

    # form-factor matrix: n_ff cosine rays per face -> first-hit face (or escape=sky)
    src = np.repeat(np.arange(F), n_ff)
    o = (centroids[src] + 1e-3 * fn[src]).astype(np.float32)
    d = _cosine_dirs(fn[src], rng).astype(np.float32)
    hf = wp.zeros(F * n_ff, dtype=int, device=DEVICE)
    wp.launch(_trace_ff, dim=F * n_ff, device=DEVICE,
              inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEVICE), wp.array(d, dtype=wp.vec3, device=DEVICE),
                      hf, float(Lx), float(Ly), float(Lz), int(par.get('periodic_y', 0))])
    hit = hf.numpy()
    esc = hit < 0
    D = np.bincount(src[esc], minlength=F).astype(np.float64) / n_ff          # D_i = sky view of i
    j = src[~esc]; i = hit[~esc]                                              # ray from j hit i -> A[i,j]=F_{j->i}
    Aff = sp.coo_matrix((np.full(len(i), 1.0 / n_ff), (i, j)), shape=(F, F)).tocsr()

    # coverage <-> radiosity fixed point (re-emission factor (1-s_j) depends on bare coverage)
    bare = np.ones(F)
    rsolver = par.get('radiosity_solver', 'jacobi')       # 'jacobi' | 'gmres' (Craig's matrix-free GMRES)
    def _jacobi(D, M):
        G = D.copy()
        for _ in range(40):                               # Jacobi: converges in ~1/s iters, noise-free
            G = D + M.dot(G)
        return G
    def solve(beta):
        s = np.clip(bare * beta, 0.0, 1.0)
        M = Aff.multiply((1.0 - s)[None, :]).tocsr()      # M[i,j] = A[i,j]*(1-s_j)
        if rsolver == 'gmres':                            # solve (I-M)G=D matrix-free; better-conditioned
            from scipy.sparse.linalg import LinearOperator, gmres   # at high albedo (low s) than Jacobi
            n = len(D)
            op = LinearOperator((n, n), matvec=lambda x: x - M.dot(x), dtype=float)
            try:
                G, info = gmres(op, D, rtol=1e-3, atol=0.0, restart=12, maxiter=40)
            except TypeError:                             # older scipy: rtol -> tol
                G, info = gmres(op, D, tol=1e-3, restart=12, maxiter=40)
            if info != 0:                                 # non-convergence -> safe Jacobi fallback
                G = _jacobi(D, M)
        else:
            G = _jacobi(D, M)
        return np.clip(G, 0.0, 8.0)
    m_F = m_O = np.zeros(F)
    for _ in range(n_fp):
        m_F = solve(betaE); m_O = solve(betaO)
        thF, thO = _belen_coverages(m_i, m_F, m_O, cos_i, par, flags)
        bare = np.clip(1.0 - thF - thO, 0.0, 1.0)
    return m_i, m_F, m_O, cos_i


def _ions_deterministic(mesh, F, A, A_src, Lx, Ly, Lz, n_ion, par, flags, seed, oi_w, di_w):
    """Shared ion launch for the deterministic-neutral paths (knudsen / dda / radiosity).
    With flags.ion_reflection + belen (the validated accuracy config), runs the FAITHFUL ViennaPS
    coned-cosine reflection kernel and deposits par['_ion_yield'] -- ions funnel to the deep floor,
    keeping the ion-limited channel ~AR-independent (Blauw 2000). Previously these paths silently
    ignored ion_reflection and used the legacy first-hit ion, whose floor flux decays ~1/AR -- the
    dominant deep-tail steepener vs the de Boer wafer. Returns (m_i, cos_i)."""
    faithful = (flags is not None and getattr(flags, "ion_reflection", False)
                and getattr(flags, "chemistry", "langmuir") == "belen")
    fi = wp.zeros(F, dtype=float, device=DEVICE); ai = wp.zeros(F, dtype=float, device=DEVICE)
    if faithful:
        inflect = float(par.get('inflectAngle', 1.55334303)); minang = float(par.get('minAngle', 1.3962634))
        n_l = float(par.get('n_l', 10.0)); A_en = 1.0 / (1.0 + n_l * (1.5707963 / inflect - 1.0))
        minE = float(min(par['Eth_ie'], par['Eth_sp']))
        sp_w = wp.zeros(F, dtype=float, device=DEVICE); ie_w = wp.zeros(F, dtype=float, device=DEVICE)
        ip_w = wp.zeros(F, dtype=float, device=DEVICE)
        wp.launch(_trace3d_ion_yield, dim=n_ion, device=DEVICE,
                  inputs=[mesh.id, oi_w, di_w, int(seed), float(par['Emean']), float(par['Esig']),
                          float(par['A_sp']), float(par['Eth_sp']), float(par['B_sp']),
                          float(par['A_ie']), float(par['Eth_ie']), float(par['A_p']), float(par['Eth_p']),
                          inflect, minang, n_l, float(A_en), minE,
                          float(par.get('thetaRMin', 1.2217305)), float(par.get('thetaRMax', 1.5707963)),
                          Lx, Ly, Lz, int(par.get('periodic_y', 0)),
                          sp_w, ie_w, ip_w, fi, ai])
        norm = A_src / n_ion
        par['_ion_yield'] = ((sp_w.numpy() / A) * norm, (ie_w.numpy() / A) * norm, (ip_w.numpy() / A) * norm)
        fi = fi.numpy(); ai = ai.numpy()
        m_i = (fi / A) * norm
    else:
        par.pop('_ion_yield', None)
        wp.launch(_trace3d, dim=n_ion, device=DEVICE,
                  inputs=[mesh.id, oi_w, di_w, 1.0, 0, int(seed), 0, 0.34, 0.8, fi, ai])
        fi = fi.numpy(); ai = ai.numpy()
        m_i = np.clip((fi / A) / (n_ion / A_src), 0.0, 1.5)
    cos_i = np.where(fi > 0, ai / np.maximum(fi, 1e-9), 0.0)
    return m_i, cos_i


def mc_flux_3d_knudsen(mesh, verts, faces, centroids, areas, geo, par, n_ion=20000,
                       seed=0, flags=None, n_fp=4):
    """KNUDSEN deep-neutral tail (vs many-bounce MC). Ions stay MC (directional, single-bounce);
    neutrals come from a grid-native 1-D molecular-flow conductance solve down the feature
    (knudsen.knudsen_face_flux), coupled to the Belen coverage fixed point. A cheap high-AR
    alternative to deep MC re-emission, ported from Craig Xu Chen's plasma_sim/solver3d.py.
    By design accuracy-neutral vs petch's MC (which already samples the free-molecular regime):
    this is a SPEED/robustness option -- benchmark before claiming any accuracy delta."""
    from .knudsen import knudsen_face_flux
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']; A_src = Lx * Ly
    A = np.maximum(areas, 1.0e-9)
    betaE = par.get('betaE', 0.7); betaO = par.get('betaO', 1.0)
    fn = _gas_normals(verts, faces, centroids, geo)
    gas_nz = fn[:, 2]
    shape = "via" if geo.get('hole', False) else "trench"
    wls = float(par.get('knudsen_wall_loss_scale', 1.85))
    phi, zs, dx, sub_top = geo['phi'], geo['zs'], geo['dx'], geo['sub_top']

    # ions: shared deterministic-path launch (faithful ViennaPS reflection when the config asks).
    oi, di = _source3d('ion', n_ion, Lx, Ly, z_src, par['ion_ang_sigma'], 'sobol', rng, seed * 9 + 1)
    oi_w = wp.array(oi, dtype=wp.vec3, device=DEVICE); di_w = wp.array(di, dtype=wp.vec3, device=DEVICE)
    m_i, cos_i = _ions_deterministic(mesh, F, A, A_src, Lx, Ly, Lz, n_ion, par, flags,
                                     seed * 9 + 1, oi_w, di_w)

    # neutrals: 1-D Knudsen conductance profile, coupled to the Belen coverage fixed point.
    # knudsen_sink='local' (default): sink sticking = betaE*bare with the LOCAL coverage -- as the
    # floor starves, bare->1 and the sink GROWS, steepening the tail toward ~conc^2 decay.
    # knudsen_sink='field': clamp the sink at the FIELD (well-fed) coverage's value, so the sink
    # coefficient stays at its saturation level -- conc = 1/(1+k*AR), the flattening the real
    # de Boer tail shows. Physically: the starved floor's F consumption is capped by the ion-driven
    # demand it had when saturated, not boosted by emptier sites.
    sink_mode = par.get('knudsen_sink', 'local')
    field = centroids[:, 2] > sub_top - 0.5 * dx           # open-field faces (top surface)
    bare = np.ones(F)
    m_F = m_O = np.zeros(F)
    for _ in range(n_fp):
        s_bare = bare
        if sink_mode == 'field' and field.any():
            s_bare = np.minimum(bare, float(np.mean(bare[field])))
        m_F = knudsen_face_flux(phi, zs, dx, sub_top, shape, centroids, gas_nz,
                                np.clip(s_bare * betaE, 0.0, 1.0), wls)
        m_O = knudsen_face_flux(phi, zs, dx, sub_top, shape, centroids, gas_nz,
                                np.clip(s_bare * betaO, 0.0, 1.0), wls)
        thF, thO = _belen_coverages(m_i, m_F, m_O, cos_i, par, flags)
        bare = np.clip(1.0 - thF - thO, 0.0, 1.0)
    return m_i, m_F, m_O, cos_i


def mc_flux_3d_dda(mesh, verts, faces, centroids, areas, geo, par, n_ion=20000,
                   seed=0, flags=None, n_fp=4):
    """DETERMINISTIC DDA neutral transport (discrete-ordinates grid-march, ported from Craig Xu
    Chen's plasma_sim). Ions stay MC (directional, single-bounce); neutrals via dda.dda_neutral_flux
    -- a fixed-quadrature grid-march with diffuse re-emission, coupled to the Belen coverage fixed
    point. Noise-free deep-AR rolloff (no MC floor-starvation). neutral_transport='dda'."""
    from .dda import dda_neutral_flux
    Lx, Ly, Lz = geo['Lx'], geo['Ly'], geo['Lz']
    F = len(faces)
    rng = np.random.default_rng(seed)
    z_src = Lz - geo['dx']; A_src = Lx * Ly
    A = np.maximum(areas, 1.0e-9)
    betaE = par.get('betaE', 0.7); betaO = par.get('betaO', 1.0)
    fn = _gas_normals(verts, faces, centroids, geo)
    phi, zs, dx = geo['phi'], geo['zs'], geo['dx']
    n_dir = int(par.get('dda_n_dir', 64)); n_re = int(par.get('dda_n_reemit', 12))

    # ions: shared deterministic-path launch (faithful ViennaPS reflection when the config asks)
    oi, di = _source3d('ion', n_ion, Lx, Ly, z_src, par['ion_ang_sigma'], 'sobol', rng, seed * 9 + 1)
    oi_w = wp.array(oi, dtype=wp.vec3, device=DEVICE); di_w = wp.array(di, dtype=wp.vec3, device=DEVICE)
    m_i, cos_i = _ions_deterministic(mesh, F, A, A_src, Lx, Ly, Lz, n_ion, par, flags,
                                     seed * 9 + 1, oi_w, di_w)

    # neutrals: deterministic DDA, coupled to the Belen coverage fixed point
    bare = np.ones(F)
    m_F = m_O = np.zeros(F)
    for _ in range(n_fp):
        m_F = dda_neutral_flux(phi, dx, zs, centroids, fn, np.clip(bare * betaE, 0.0, 1.0), n_dir, n_re)
        m_O = dda_neutral_flux(phi, dx, zs, centroids, fn, np.clip(bare * betaO, 0.0, 1.0), n_dir, n_re)
        thF, thO = _belen_coverages(m_i, m_F, m_O, cos_i, par, flags)
        bare = np.clip(1.0 - thF - thO, 0.0, 1.0)
    return m_i, m_F, m_O, cos_i


# ----------------------------- 3D advection -----------------------------
def advect_3d(phi, Fspeed, dx, dt):
    """phi_t + F|grad phi| = 0, first-order upwind Godunov in 3D (F>=0 etch). NumPy reference."""
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


@wp.kernel
def _advect_iter(phi: wp.array3d(dtype=float), F: wp.array3d(dtype=float),
                 mask: wp.array3d(dtype=float), mask_phi: wp.array3d(dtype=float),
                 out: wp.array3d(dtype=float), inv_dx: float, dt: float):
    """One upwind-Godunov advection substep (matches advect_3d) + mask re-pin, fused. On GPU."""
    i, j, k = wp.tid()
    if mask[i, j, k] > 0.5:
        out[i, j, k] = mask_phi[i, j, k]                # re-pin mask material
        return
    nx = phi.shape[0]; ny = phi.shape[1]; nz = phi.shape[2]
    c = phi[i, j, k]
    g = float(0.0)
    dm = float(0.0); dp = float(0.0)
    if i >= 1:
        dm = (c - phi[i - 1, j, k]) * inv_dx
    if i <= nx - 2:
        dp = (phi[i + 1, j, k] - c) * inv_dx
    g += wp.max(dm, 0.0) * wp.max(dm, 0.0) + wp.min(dp, 0.0) * wp.min(dp, 0.0)
    dm = 0.0; dp = 0.0
    if j >= 1:
        dm = (c - phi[i, j - 1, k]) * inv_dx
    if j <= ny - 2:
        dp = (phi[i, j + 1, k] - c) * inv_dx
    g += wp.max(dm, 0.0) * wp.max(dm, 0.0) + wp.min(dp, 0.0) * wp.min(dp, 0.0)
    dm = 0.0; dp = 0.0
    if k >= 1:
        dm = (c - phi[i, j, k - 1]) * inv_dx
    if k <= nz - 2:
        dp = (phi[i, j, k + 1] - c) * inv_dx
    g += wp.max(dm, 0.0) * wp.max(dm, 0.0) + wp.min(dp, 0.0) * wp.min(dp, 0.0)
    out[i, j, k] = c - dt * F[i, j, k] * wp.sqrt(g)


def advect_3d_gpu(phi, Fspeed, mask, mask_phi, dx, dt, nsub):
    """All CFL substeps on the GPU: phi/F/mask stay on-device, ping-pong buffers, no per-substep CPU
    round-trip. The advect substep loop was the dominant CPU host op after narrow-band reinit."""
    a = wp.array(phi.astype(np.float32), dtype=float, device=DEVICE)
    Fw = wp.array(Fspeed.astype(np.float32), dtype=float, device=DEVICE)
    mk = wp.array(mask.astype(np.float32), dtype=float, device=DEVICE)
    mp = wp.array(mask_phi.astype(np.float32), dtype=float, device=DEVICE)
    b = wp.zeros_like(a)
    inv = 1.0 / dx
    for _ in range(nsub):
        wp.launch(_advect_iter, dim=phi.shape, device=DEVICE, inputs=[a, Fw, mk, mp, b, inv, dt])
        a, b = b, a
    return a.numpy().astype(np.float64)


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
        # Floor |grad phi0| at 0.5 (not 1e-9): a near-flat cell mis-flagged as interface (deep-hole
        # corner / re-pinned mask edge) otherwise makes D = s0/grad0 explode -> inf -> NaN spreads.
        grad0 = wp.max(wp.sqrt(gx0 * gx0 + gy0 * gy0 + gz0 * gz0), 0.5)
        dxl = 1.0 / inv_dx
        D = wp.clamp(s0 / grad0, -1.8 * dxl, 1.8 * dxl)   # |dist| <= cell diagonal for an interface cell
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
    # CFL: the Godunov |grad phi|=1 sweep is forward-Euler; in 3D the stable step is ~dx/sqrt(3).
    # 0.5*dx is borderline (works shallow, grows under stress -> blowup); 0.3*dx has clear margin.
    dtau = 0.3 * dx
    inv = 1.0 / dx
    for _ in range(n_iter):
        wp.launch(_reinit_iter, dim=phi_np.shape, device=DEVICE, inputs=[a, phi0, b, inv, dtau])
        a, b = b, a
    return a.numpy().astype(np.float64)


@wp.func
def _eik_solve(b1: float, b2: float, b3: float, dx: float) -> float:
    """Godunov solution of the Eikonal |grad U|=1 at one node: solve sum_d (x-a_d)_+^2 = dx^2 for the
    three upwind neighbor minima a_d. Staged 1->2->3 term form (sort ascending, add terms while the
    candidate exceeds the next neighbor). This is the EXACT Godunov update -> the Jacobi fixed point
    has |grad U|=1 by construction (no PDE/forward-Euler bias). Boundaries pass `big` for the absent
    neighbor so that term drops out (an absent axis never enters the active set)."""
    a1 = b1; a2 = b2; a3 = b3
    t = float(0.0)
    if a1 > a2:
        t = a1; a1 = a2; a2 = t
    if a2 > a3:
        t = a2; a2 = a3; a3 = t
    if a1 > a2:
        t = a1; a1 = a2; a2 = t
    x = a1 + dx
    if x <= a2:
        return x
    # two-term
    d2 = 2.0 * dx * dx - (a1 - a2) * (a1 - a2)
    if d2 < 0.0:
        d2 = 0.0
    x = 0.5 * (a1 + a2 + wp.sqrt(d2))
    if x <= a3:
        return x
    # three-term
    s1 = a1 + a2 + a3
    s2 = a1 * a1 + a2 * a2 + a3 * a3
    d3 = s1 * s1 - 3.0 * (s2 - dx * dx)
    if d3 < 0.0:
        d3 = 0.0
    return (s1 + wp.sqrt(d3)) / 3.0


@wp.kernel
def _eik_init(phi0: wp.array3d(dtype=float), U: wp.array3d(dtype=float),
              frozen: wp.array3d(dtype=float), inv_dx: float, big: float):
    """Seed the unsigned distance U: interface-adjacent cells get the Russo-Smereka sub-cell distance
    |phi0|/|grad phi0| (FROZEN = the boundary condition that pins the phi=0 contour); all others = big."""
    i, j, k = wp.tid()
    nx = phi0.shape[0]; ny = phi0.shape[1]; nz = phi0.shape[2]
    c = phi0[i, j, k]
    o_xm = phi0[wp.max(i - 1, 0), j, k]; o_xp = phi0[wp.min(i + 1, nx - 1), j, k]
    o_ym = phi0[i, wp.max(j - 1, 0), k]; o_yp = phi0[i, wp.min(j + 1, ny - 1), k]
    o_zm = phi0[i, j, wp.max(k - 1, 0)]; o_zp = phi0[i, j, wp.min(k + 1, nz - 1)]
    interface = (c * o_xm < 0.0) or (c * o_xp < 0.0) or (c * o_ym < 0.0) \
        or (c * o_yp < 0.0) or (c * o_zm < 0.0) or (c * o_zp < 0.0)
    if interface:
        gx0 = (o_xp - o_xm) * 0.5 * inv_dx
        gy0 = (o_yp - o_ym) * 0.5 * inv_dx
        gz0 = (o_zp - o_zm) * 0.5 * inv_dx
        grad0 = wp.max(wp.sqrt(gx0 * gx0 + gy0 * gy0 + gz0 * gz0), 0.5)
        dxl = 1.0 / inv_dx
        D = wp.min(wp.abs(c) / grad0, 1.8 * dxl)     # unsigned sub-cell distance, capped at the diagonal
        U[i, j, k] = D
        frozen[i, j, k] = 1.0
    else:
        U[i, j, k] = big
        frozen[i, j, k] = 0.0


@wp.kernel
def _eik_jacobi(U: wp.array3d(dtype=float), frozen: wp.array3d(dtype=float),
                out: wp.array3d(dtype=float), dx: float, big: float):
    """One Jacobi sweep of the Godunov Eikonal solver: out = min(U, godunov(neighbor minima)). Frozen
    cells are held. Monotone (U only decreases) and unconditionally stable -- no CFL, no blowup. The
    parallel-fast-sweeping fixed point (Detrixhe-Gibou-Min) reached by all-node-parallel min-updates."""
    i, j, k = wp.tid()
    if frozen[i, j, k] > 0.5:
        out[i, j, k] = U[i, j, k]
        return
    nx = U.shape[0]; ny = U.shape[1]; nz = U.shape[2]
    um = big; up = big
    if i >= 1:
        um = U[i - 1, j, k]
    if i <= nx - 2:
        up = U[i + 1, j, k]
    ax = wp.min(um, up)
    um = big; up = big
    if j >= 1:
        um = U[i, j - 1, k]
    if j <= ny - 2:
        up = U[i, j + 1, k]
    ay = wp.min(um, up)
    um = big; up = big
    if k >= 1:
        um = U[i, j, k - 1]
    if k <= nz - 2:
        up = U[i, j, k + 1]
    az = wp.min(um, up)
    cand = _eik_solve(ax, ay, az, dx)
    out[i, j, k] = wp.min(U[i, j, k], cand)


def reinit_fsm(phi_np, dx, band, n_iter=None):
    """GPU narrow-band reinit via Jacobi Godunov-Eikonal fast sweeping (Warp, fully on-device).

    Replaces CPU skfmm AND the old PDE `reinit_gpu`: the Godunov solve enforces |grad phi|=1 EXACTLY at
    the fixed point (no Russo-Smereka-PDE masked-front |grad|=1.32 bias) and the min-update is monotone
    (no forward-Euler CFL blowup). Distance info propagates 1 cell/sweep, so ~band/dx + a few sweeps
    cover the band -- the rest of the grid is left at the far-field sign*(band+dx) like reinit_narrow.
    On a GPU this removes the per-step CPU round-trip (the ~40% reinit bottleneck) and is graph-capturable."""
    nx = phi_np.shape
    big = float(band + 4.0 * dx)
    if n_iter is None:
        n_iter = int(np.ceil(band / dx)) + 4
    phi0 = wp.array(phi_np.astype(np.float32), dtype=float, device=DEVICE)
    U = wp.zeros(nx, dtype=float, device=DEVICE)
    frozen = wp.zeros(nx, dtype=float, device=DEVICE)
    wp.launch(_eik_init, dim=nx, device=DEVICE, inputs=[phi0, U, frozen, 1.0 / dx, big])
    b = wp.zeros(nx, dtype=float, device=DEVICE)
    for _ in range(n_iter):
        wp.launch(_eik_jacobi, dim=nx, device=DEVICE, inputs=[U, frozen, b, dx, big])
        U, b = b, U
    Un = U.numpy().astype(np.float64)
    sgn = np.sign(phi_np); sgn = np.where(sgn == 0, 1.0, sgn)
    return sgn * np.minimum(Un, band + dx)


def reinit_narrow(phi, dx, band):
    """NARROW-BAND reinit: skfmm fast-marching ONLY within `band` of the front (skfmm's `narrow=`),
    not the dense grid. The front only moves <1 cell/step (CFL) so the band is all that's needed; the
    far field keeps its sign at a large magnitude. ~5x faster than full reinit with EXACT agreement in
    the band -- this is the SOTA narrow-band approach (ViennaLS/HRLE do the same: no global reinit).
    Was the self-inflicted ~42% bottleneck vs ViennaPS, which never does a full per-step reinit."""
    d = skfmm.distance(phi, dx=dx, narrow=band)
    masked = np.ma.getmaskarray(d)
    out = np.ma.filled(d, 0.0).astype(np.float64)
    out[masked] = np.sign(phi[masked]) * (band + dx)   # far field: correct sign, |phi|>band
    return out


def _gas_normals(verts, faces, centroids, geo):
    """Face normals oriented INTO THE GAS (phi<0 side) -- essential for face-emission ray casting
    (radiosity form factors, redeposition). The marching-cubes cross-product normal has arbitrary sign;
    we flip it to align with -grad(phi) (the SDF normal points toward solid, so into-gas = -grad)."""
    vv = verts[faces]
    fn = np.cross(vv[:, 1] - vv[:, 0], vv[:, 2] - vv[:, 0])
    fn = fn / (np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12)
    dx = geo['dx']
    gx, gy, gz = np.gradient(geo['phi'], dx)
    ix = np.clip((centroids[:, 0] / dx).round().astype(int), 0, gx.shape[0] - 1)
    iy = np.clip((centroids[:, 1] / dx).round().astype(int), 0, gx.shape[1] - 1)
    iz = np.clip((centroids[:, 2] / dx).round().astype(int), 0, gx.shape[2] - 1)
    into_gas = -np.stack([gx[ix, iy, iz], gy[ix, iy, iz], gz[ix, iy, iz]], axis=1)
    flip = np.sum(fn * into_gas, axis=1) < 0.0
    fn[flip] *= -1.0
    return fn


def _cosine_dirs(normals, rng):
    """Cosine-weighted hemisphere launch directions about each unit normal (diffuse emission)."""
    n = normals
    a = np.where(np.abs(n[:, :1]) > 0.9, np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    t = np.cross(a, n); t /= (np.linalg.norm(t, axis=1, keepdims=True) + 1e-12)
    b = np.cross(n, t)
    u1 = rng.random(len(n)); u2 = rng.random(len(n))
    ct = np.sqrt(u1); st = np.sqrt(1.0 - u1); ph = 2.0 * np.pi * u2
    return (st * np.cos(ph))[:, None] * t + (st * np.sin(ph))[:, None] * b + ct[:, None] * n


def mc_redep_3d(mesh, centroids, areas, normals, src_flux, n_redep=20000, s_redep=0.5, seed=0):
    """Etch-product REDEPOSITION (beyond ViennaPS, which omits it). Emit product from each face in
    proportion to its etch/sputter rate (src_flux), trace ballistically with diffuse re-emission, and
    re-stick with probability s_redep -> redeposited flux per face. Because it re-emits and penetrates,
    it deposits preferentially on the lower sidewalls in line-of-sight of the etching floor -> sidewall
    passivation -> taper / top-narrowing (the Gomez SF6:O2-driven profile change). Reuses the exact RR
    trace kernel. Returns redeposited flux density per face (subtract k_redep*this from the etch rate)."""
    F = len(centroids)
    w = np.maximum(src_flux * areas, 0.0)
    tot = float(w.sum())
    if tot <= 0.0 or F == 0:
        return np.zeros(F)
    rng = np.random.default_rng(seed * 7 + 5)
    idx = rng.choice(F, size=n_redep, p=w / tot)
    o = (centroids[idx] + 1.0e-3 * normals[idx]).astype(np.float32)
    d = _cosine_dirs(normals[idx], rng).astype(np.float32)
    fl = wp.zeros(F, dtype=float, device=DEVICE)
    bare = np.ones(F, dtype=np.float32)            # product sticks anywhere it lands (S_eff = s_redep)
    wp.launch(_trace3d_cov_rr, dim=n_redep, device=DEVICE,
              inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEVICE),
                      wp.array(d, dtype=wp.vec3, device=DEVICE),
                      wp.array(bare, dtype=float, device=DEVICE), float(s_redep), int(seed), fl,
                      1.0, 1.0, 0])                     # redep: no periodic-y
    return (fl.numpy() * (tot / n_redep)) / np.maximum(areas, 1.0e-9)


def faces_in_mask(centroids, geo, mask_th, trench_width, hole=False):
    """Mark faces whose centroid lies in the (un-etched) mask material."""
    x, y, z = centroids[:, 0], centroids[:, 1], centroids[:, 2]
    in_band = (z >= geo['sub_top']) & (z < geo['sub_top'] + mask_th)
    if hole:
        opening = (x - geo['Lx']/2)**2 + (y - geo['Ly']/2)**2 < (trench_width/2)**2
    else:
        opening = np.abs(x - geo['Lx']/2) < trench_width/2
    return in_band & (~opening)


@wp.kernel
def _nearest_face_cov(mesh: wp.uint64, pts: wp.array(dtype=wp.vec3),
                      prev_bare: wp.array(dtype=float), out: wp.array(dtype=float)):
    """For each new face centroid, find the nearest face on the PREVIOUS mesh and copy its coverage.
    The GPU (BVH mesh_query_point) replacement for the scipy warm-start KDTree -- same nearest-face
    seed, no host round-trip (the KDTree was ~15ms/step on deep meshes, the top flux host cost)."""
    i = wp.tid()
    q = wp.mesh_query_point_no_sign(mesh, pts[i], 1.0e6)
    if q.result:
        out[i] = prev_bare[q.face]
    else:
        out[i] = 1.0


def gpu_warmstart_bare(prev_mesh, prev_bare_wp, centroids):
    """Warm-start coverage seed via GPU nearest-face query on the previous mesh (replaces scipy KDTree)."""
    n = len(centroids)
    pts = wp.array(centroids.astype(np.float32), dtype=wp.vec3, device=DEVICE)
    out = wp.zeros(n, dtype=float, device=DEVICE)
    wp.launch(_nearest_face_cov, dim=n, device=DEVICE, inputs=[prev_mesh.id, pts, prev_bare_wp, out])
    return out.numpy()


# ----------------------------- driver -----------------------------
def run_etch_3d(Lx=10.0, Ly=4.0, Lz=14.0, dx=0.4, trench_width=4.0, mask_th=2.0,
                sub_top=10.0, t_end=2.0, n_steps=20, hole=False, par=None, flags=None,
                n_ion=20000, n_neu=20000, reinit_every=1, extend="gpu",
                reinit_method="skfmm", verbose=True, record_depth_every=0, seed_offset=0,
                rays_per_point=None, record_frames=False, surf_smooth=0.0):
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
    timings = dict(flux=0.0, extend=0.0, reinit=0.0, mesh=0.0, advect=0.0, total=0.0, nsub_max=0)
    warm = getattr(flags, "warm_start_coverage", False)
    _cuda = DEVICE == 'cuda'                     # GPU speedups auto-enable on CUDA, fall back on CPU (portable)
    gpu_ws = par.get('gpu_warmstart', _cuda) and _cuda          # GPU nearest-face warm-start vs scipy KDTree
    use_gpu_mesh = par.get('gpu_mesh', _cuda) and _cuda         # Warp MarchingCubes is CUDA-only -> guard
    _extract = extract_mesh_3d_gpu if use_gpu_mesh else extract_mesh_3d
    cov_centroids = None; cov_bare = None       # previous-step coverage, for warm-starting the fixed point
    prev_mesh = None; prev_bare_wp = None       # previous wp.Mesh + coverage (device) for GPU warm-start
    t0 = time.time()
    for step in range(n_steps):
        tm = time.time()
        verts, faces, centroids, areas = _extract(geo['phi'], dx)
        mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEVICE),
                       indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEVICE))
        timings['mesh'] += time.time() - tm
        # ViennaPS-style ray budget: total rays = rays_per_point * #facets, so the per-facet
        # sampling stays CONSTANT under grid refinement (a fixed total budget starves the deep
        # floor at fine dx -> grid-sensitive ARDE). Sobol QMC lets petch use a far smaller
        # rays_per_point than ViennaPS's 1000 for the same noise floor (keeps it fast).
        if rays_per_point is not None:
            n_eff = max(int(rays_per_point) * len(faces), 1000)
            n_ion_s, n_neu_s = n_eff, n_eff
        else:
            n_ion_s, n_neu_s = n_ion, n_neu
        tf = time.time()
        _nt = getattr(flags, "neutral_transport", "mc")
        if _nt == "radiosity":   # deterministic radiosity neutrals
            m_i, m_F, m_O, cos_i = mc_flux_3d_radiosity(mesh, verts, faces, centroids, areas, geo, par,
                                                        n_ion=n_ion_s, seed=step + seed_offset, flags=flags)
        elif _nt == "knudsen":   # 1-D Knudsen conductance deep-neutral tail (ported from plasma_sim)
            m_i, m_F, m_O, cos_i = mc_flux_3d_knudsen(mesh, verts, faces, centroids, areas, geo, par,
                                                      n_ion=n_ion_s, seed=step + seed_offset, flags=flags)
        elif _nt == "dda":       # deterministic discrete-ordinates DDA neutrals (ported from plasma_sim)
            m_i, m_F, m_O, cos_i = mc_flux_3d_dda(mesh, verts, faces, centroids, areas, geo, par,
                                                  n_ion=n_ion_s, seed=step + seed_offset, flags=flags)
        elif getattr(flags, "coverage_sticking", False):   # Langmuir coverage-dependent sticking
            bi = None
            if warm and cov_bare is not None and np.all(np.isfinite(centroids)):   # seed from prev-step coverage
                if gpu_ws and prev_mesh is not None:           # GPU nearest-face query (no host KDTree)
                    bi = gpu_warmstart_bare(prev_mesh, prev_bare_wp, centroids)
                elif np.all(np.isfinite(cov_centroids)):       # scipy KDTree (fast build + parallel query)
                    tree = cKDTree(cov_centroids, balanced_tree=False, compact_nodes=False)
                    _, ix = tree.query(centroids, workers=-1)
                    bi = cov_bare[ix]
            m_i, m_F, m_O, cos_i, cov_bare = mc_flux_3d_coupled(mesh, verts, faces, areas, geo, par,
                                                      n_ion=n_ion_s, n_neu=n_neu_s, seed=step + seed_offset,
                                                      sampling=getattr(flags, "sampling", "pseudo"),
                                                      flags=flags, bare_init=bi)
            cov_centroids = centroids
            if gpu_ws:                                         # keep prev mesh + coverage on device
                prev_mesh = mesh
                prev_bare_wp = wp.array(cov_bare.astype(np.float32), dtype=float, device=DEVICE)
        else:
            m_i, m_F, m_O, cos_i = mc_flux_3d(mesh, verts, faces, areas, geo, mc_par,
                                              n_ion=n_ion_s, n_neu=n_neu_s, seed=step + seed_offset,
                                              sampling=getattr(flags, "sampling", "pseudo"),
                                              ion_reflection=getattr(flags, "ion_reflection", False))
        timings['flux'] += time.time() - tf
        trt = time.time()
        is_mask = faces_in_mask(centroids, geo, mask_th, trench_width, hole=hole)
        V = surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=flags)
        _nbad = int((~np.isfinite(V)).sum())                   # non-finite velocity = a real problem
        if _nbad:                                              # warn (don't silently hide) then guard
            import warnings
            warnings.warn(f"step {step}: {_nbad} non-finite surface velocities (NaN/inf) zeroed; "
                          "check flux normalization / coverage if this persists.", RuntimeWarning)
        V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)  # guard so the run doesn't crash
        timings['rate'] = timings.get('rate', 0.0) + time.time() - trt
        if getattr(flags, "redeposition", False):    # etch-product redeposition -> sidewall passivation
            fn = _gas_normals(verts, faces, centroids, geo)
            Rf = mc_redep_3d(mesh, centroids, areas, fn, V, n_redep=n_neu,
                             s_redep=par.get('s_redep', 0.5), seed=step)
            V = np.maximum(V - par.get('k_redep', 1.0) * Rf, 0.0)   # redeposited material slows etch
        te = time.time()
        Fs = (extend_velocity_gpu(mesh, V, geo, band) if extend == "gpu"
              else extend_velocity_3d(V, centroids, geo, band))
        timings['extend'] += time.time() - te
        vmx = float(np.max(V)) if V.size else 0.0
        Vmax = max(vmx if np.isfinite(vmx) else 0.0, 1e-6)
        # CFL substepping: nsub so each advect moves < 0.4*dx. Cap raised 40->160 because the exact
        # (Russian-roulette) neutral transport feeds the floor strongly -> high surface velocity; a
        # cap of 40 (max stable Vmax~53) blew up at moderate rates -> spurious deep depth. 160 -> ~190.
        nsub = max(1, min(int(np.ceil(Vmax * dt / (0.4 * dx))), 160))
        timings['nsub_max'] = max(timings['nsub_max'], nsub)
        ta = time.time()
        if extend == "gpu":            # all substeps on GPU (no per-substep CPU round-trip)
            geo['phi'] = advect_3d_gpu(geo['phi'], Fs, geo['mask'], mask_phi, dx, dt / nsub, nsub)
        else:
            for _ in range(nsub):
                geo['phi'] = advect_3d(geo['phi'], Fs, dx, dt / nsub)
                geo['phi'][geo['mask']] = mask_phi[geo['mask']]
        timings['advect'] += time.time() - ta
        # reinit_every>1 (lazy reinit) is faster but DRIFTS the result: |grad phi| deviates from
        # 1 between reinits and advect multiplies F*|grad phi|. Safe only with a proper extension
        # velocity (grad F . grad phi = 0). Keep =1 for fidelity.
        if (step + 1) % reinit_every == 0 or step == n_steps - 1:
            tr = time.time()
            if reinit_method == "fsm":          # GPU Jacobi Godunov-Eikonal fast sweep (no CPU round-trip)
                geo['phi'] = reinit_fsm(geo['phi'], dx, band + 2.0 * dx)
            elif reinit_method == "gpu":
                geo['phi'] = reinit_gpu(geo['phi'], dx)
            elif reinit_method == "skfmm_full":
                geo['phi'] = skfmm.distance(geo['phi'], dx=dx)
            else:                                  # default: SOTA narrow-band (skfmm 'narrow'), ~5x faster
                geo['phi'] = reinit_narrow(geo['phi'], dx, band + 2.0 * dx)
            timings['reinit'] += time.time() - tr
        if surf_smooth and surf_smooth > 0.0:      # light surface regularization (curvature/diffusion-like):
            from scipy.ndimage import gaussian_filter   # suppresses noise-seeded mask-edge fingering, like
            geo['phi'] = gaussian_filter(geo['phi'], surf_smooth)   # ViennaLS's implicit front regularization
            geo['phi'][geo['mask']] = mask_phi[geo['mask']]         # keep the mask pinned
        if record_depth_every and (step % record_depth_every == 0 or step == n_steps - 1):
            geo.setdefault('depth_history', []).append((step + 1, _depth3d(geo)))
            if record_frames:                               # stash the centre x-z slice for animation
                jc = geo['phi'].shape[1] // 2
                geo.setdefault('frames', []).append(
                    dict(step=step + 1, t=(step + 1) / n_steps * t_end,
                         depth=_depth3d(geo), phi_xz=geo['phi'][:, jc, :].copy()))
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
    # Clamp the sampling region to the FEATURE interior. Otherwise, for a trench/hole narrower than
    # `half`, the region also covers UNETCHED substrate columns (floor at sub_top) -> the median is
    # dragged to the surface and the reported depth collapses toward 0 (a metric artifact that looks
    # like grid 'instability' at fine dx). Sample only columns inside the opening.
    cx, cy = geo['Lx']/2, geo['Ly']/2
    r_feat = geo.get('trench_width', 2.0) / 2.0 - geo['dx']      # stay just inside the wall
    hx = max(min(half, r_feat), geo['dx'])
    ic = np.where(np.abs(xs - cx) <= hx)[0]
    if geo.get('hole', False):                                   # hole: radial footprint
        jc = np.where(np.abs(ys - cy) <= hx)[0]
        inside = lambda i, j: (xs[i]-cx)**2 + (ys[j]-cy)**2 <= r_feat**2
    else:                                                        # trench: invariant in y, full y span ok
        jc = np.where(np.abs(ys - cy) <= half)[0]
        inside = lambda i, j: True
    floors = []
    for i in ic:
        for j in jc:
            if not inside(i, j):
                continue
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


def max_depth_3d(geo):
    """Max etch depth = deepest etched point over the FEATURE footprint -- matches ViennaPS's
    `-surfaceNodes[:,2].min()` (global deepest surface node), NOT the median-center of `center_depth_3d`.
    For the smallest holes the two differ most (global-min picks the single deepest cell; median-center
    averages a floor that fills most of the hole), which skews the normalized ARDE ratio. Use this for a
    like-for-like comparison to ViennaPS hole depths."""
    phi, xs, ys, zs = geo['phi'], geo['xs'], geo['ys'], geo['zs']
    cx, cy = geo['Lx'] / 2.0, geo['Ly'] / 2.0
    r = geo.get('trench_width', 4.0) / 2.0
    if geo.get('hole', False):
        ii, jj = np.where((xs[:, None] - cx) ** 2 + (ys[None, :] - cy) ** 2 <= (r + geo['dx']) ** 2)
        cols = list(zip(ii, jj))
    else:  # trench: footprint is the opening strip in x, all y
        ix = np.where(np.abs(xs - cx) <= r + geo['dx'])[0]
        cols = [(i, j) for i in ix for j in range(len(ys))]
    floors = []
    for i, j in cols:
        col = phi[i, j, :] < 0
        k = len(col) - 1
        if not col[k]:
            continue
        while k > 0 and col[k - 1]:
            k -= 1
        floors.append(zs[k])
    if not floors:
        return 0.0
    # robust "deepest": 5th-percentile floor z (reject thin level-set filaments at the footprint
    # edge that a raw global-min would catch), matching ViennaPS's clean-surface deepest node.
    return float(geo['sub_top'] - np.percentile(floors, 5))
