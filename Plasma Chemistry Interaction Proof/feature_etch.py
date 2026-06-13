"""
feature_etch.py
A from-scratch 2D feature-scale plasma-etch simulator.

Physics (same class as ViennaPS SF6O2Etching, Belen/Langmuir model):
  - Level-set surface description (phi>0 solid, phi<0 gas), signed-distance reinit via skfmm
  - Forward Monte Carlo ballistic flux transport (ions + F etchant + O passivation),
    with shadowing and diffuse re-emission for neutrals  <-- the GPU-relevant kernel
  - Competitive-Langmuir surface coverage (fluorine etch vs oxygen passivation),
    ion-enhanced etch + physical sputter with energy/angle-dependent yields

Plasma parameters taken to match ViennaPS SF6O2 defaults:
  ion:F:O flux = 12:1800:100, ion mean energy 100 eV (sigma 10),
  Eth_ie=15, Eth_sp=20, A_ie=7, A_sp=0.0337, rho=5.02
"""
import numpy as np
import skfmm
from scipy.spatial import cKDTree
from skimage import measure
from numba import njit
import time

# ----------------------------- plasma / surface parameters -----------------------------
PAR = dict(
    ionFlux=12.0, Fflux=1800.0, Oflux=100.0,   # relative fluxes (ViennaPS SF6O2 defaults)
    Emean=100.0, Esig=10.0,                     # ion energy distribution (eV)
    ion_ang_sigma=np.deg2rad(2.5),              # ion angular spread (near-vertical)
    Eth_ie=15.0, Eth_sp=20.0, Eth_p=10.0,       # yield thresholds (eV)
    A_ie=7.0, A_sp=0.0337, A_p=3.0,             # yield prefactors
    s_F=0.20, s_O=0.30,                         # neutral sticking coefficients
    rho=5.02,                                   # substrate density factor
    rate_scale=1.0,                             # global rate calibration (set later)
)

# ----------------------------- geometry / level set -----------------------------
def make_trench(W, H, dx, trench_width, mask_thickness, sub_top):
    """Build phi (signed distance, >0 solid) and a static mask occupancy field."""
    nx, ny = int(round(W/dx)), int(round(H/dx))
    xs = (np.arange(nx)+0.5)*dx
    ys = (np.arange(ny)+0.5)*dx
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    solid = Y < sub_top                                   # substrate fills below sub_top
    mask_band = (Y >= sub_top) & (Y < sub_top+mask_thickness)
    opening = np.abs(X - W/2) < trench_width/2
    mask = mask_band & (~opening)                         # mask everywhere except the opening
    solid = solid | mask
    phi0 = np.where(solid, 1.0, -1.0)
    phi = skfmm.distance(phi0, dx=dx)                     # signed distance, >0 solid
    return X, Y, xs, ys, phi, mask, nx, ny

def extract_surface(phi, xs, ys, dx):
    """Marching-squares contour of phi=0 -> ordered polyline segments + outward normals."""
    contours = measure.find_contours(phi, 0.0)
    segs = []
    for c in contours:
        # c are (i,j) in index space; convert to physical (x,y)
        px = xs[0] + c[:,0]*dx
        py = ys[0] + c[:,1]*dx
        pts = np.column_stack([px, py])
        for k in range(len(pts)-1):
            segs.append((pts[k,0], pts[k,1], pts[k+1,0], pts[k+1,1]))
    segs = np.array(segs)  # (M,4): x0,y0,x1,y1
    if len(segs)==0:
        return segs, segs, segs, segs
    mid = np.column_stack([(segs[:,0]+segs[:,2])/2, (segs[:,1]+segs[:,3])/2])
    tang = np.column_stack([segs[:,2]-segs[:,0], segs[:,3]-segs[:,1]])
    L = np.hypot(tang[:,0], tang[:,1]); L[L==0]=1e-12
    # outward normal = -grad(phi) direction; estimate via phi gradient at midpoints
    nrm = np.column_stack([-tang[:,1], tang[:,0]]) / L[:,None]   # one of the two normals
    return segs, mid, nrm, L

def orient_normals(mid, nrm, phi, xs, ys, dx):
    """Flip normals so they point into gas (phi<0)."""
    gx, gy = np.gradient(phi, dx)
    # sample grad at midpoints (nearest node)
    ix = np.clip(((mid[:,0]-xs[0])/dx).astype(int), 0, phi.shape[0]-1)
    iy = np.clip(((mid[:,1]-ys[0])/dx).astype(int), 0, phi.shape[1]-1)
    g = np.column_stack([gx[ix,iy], gy[ix,iy]])
    gl = np.hypot(g[:,0], g[:,1]); gl[gl==0]=1e-12
    outward = -g/gl[:,None]            # outward (into gas) = -grad phi
    dot = (nrm*outward).sum(1)
    nrm = np.where(dot[:,None]<0, -nrm, nrm)
    return nrm

# ----------------------------- Monte Carlo ballistic flux (numba) -----------------------------
@njit(cache=True, fastmath=True)
def _trace(seg, nrm, is_mask, x0_src, y_src, dirs_x, dirs_y, sticking, n_reemit, rng_seed):
    """Forward-trace one species' particles. Returns per-segment hit weight + per-seg ion-angle accumulation."""
    M = seg.shape[0]
    flux = np.zeros(M)
    ang_acc = np.zeros(M)   # accumulated cos(incidence) weight (for ions)
    N = x0_src.shape[0]
    np.random.seed(rng_seed)
    for p in range(N):
        px = x0_src[p]; py = y_src
        dx_ = dirs_x[p]; dy_ = dirs_y[p]
        w = 1.0
        for bounce in range(n_reemit+1):
            # find nearest segment intersection ahead
            best_t = 1e18; best_s = -1
            for s in range(M):
                ax=seg[s,0]; ay=seg[s,1]; bx=seg[s,2]; by=seg[s,3]
                ex = bx-ax; ey = by-ay
                den = dx_*ey - dy_*ex
                if abs(den) < 1e-14:
                    continue
                # solve px+ t dx = ax + u ex
                t = ((ax-px)*ey - (ay-py)*ex)/den
                u = ((ax-px)*dy_ - (ay-py)*dx_)/den
                if t > 1e-6 and u >= -1e-6 and u <= 1.0+1e-6:
                    if t < best_t:
                        best_t = t; best_s = s
            if best_s < 0:
                break  # escaped domain
            s = best_s
            # incidence: cos angle between -dir and outward normal
            nx_ = nrm[s,0]; ny_ = nrm[s,1]
            cosang = -(dx_*nx_ + dy_*ny_)
            if cosang < 0:  # hit backface; treat as grazing
                cosang = 0.0
            if np.random.random() < sticking or bounce == n_reemit:
                flux[s] += w
                ang_acc[s] += w*cosang
                break
            else:
                # diffuse re-emit (cosine about outward normal) and continue
                hx = px + best_t*dx_; hy = py + best_t*dy_
                # tangent
                tx = -ny_; ty = nx_
                ca = np.sqrt(np.random.random())          # cosine law
                sa = np.sqrt(1.0-ca*ca)
                sign = 1.0 if np.random.random()<0.5 else -1.0
                dx_ = ca*nx_ + sign*sa*tx
                dy_ = ca*ny_ + sign*sa*ty
                px = hx + 1e-4*nx_; py = hy + 1e-4*ny_
                w *= 1.0   # weight conserved; absorption handled by sticking branch
    return flux, ang_acc

def mc_flux(seg, mid, nrm, is_mask, L, y_src, W, par, n_part_ion=20000, n_part_neu=20000):
    """Compute per-segment normalized flux multipliers + mean ion incidence cos for 3 species."""
    rng = np.random.default_rng(0)
    # --- ions: near-vertical, small angular spread ---
    xs0 = rng.uniform(0, W, n_part_ion)
    a = rng.normal(0, par['ion_ang_sigma'], n_part_ion)
    dix = np.sin(a); diy = -np.cos(a)
    fi, ai = _trace(seg, nrm, is_mask, xs0, y_src, dix, diy, 1.0, 0, 1)
    # --- F etchant: cosine launch, sticking s_F, re-emission ---
    xs1 = rng.uniform(0, W, n_part_neu)
    aF = np.arcsin(rng.uniform(-1,1,n_part_neu))     # cosine-ish into lower hemisphere
    dfx = np.sin(aF); dfy = -np.cos(aF)
    fF, _ = _trace(seg, nrm, is_mask, xs1, y_src, dfx, dfy, par['s_F'], 12, 2)
    # --- O passivation ---
    xs2 = rng.uniform(0, W, n_part_neu)
    aO = np.arcsin(rng.uniform(-1,1,n_part_neu))
    dox = np.sin(aO); doy = -np.cos(aO)
    fO, _ = _trace(seg, nrm, is_mask, xs2, y_src, dox, doy, par['s_O'], 12, 3)
    # normalize to open-field flux density (particles per unit x length)
    base_ion = n_part_ion / W
    base_neu = n_part_neu / W
    m_i = (fi/np.maximum(L,0.3*np.median(L))) / base_ion
    m_F = (fF/np.maximum(L,0.3*np.median(L))) / base_neu
    m_O = (fO/np.maximum(L,0.3*np.median(L))) / base_neu
    m_i = np.clip(m_i, 0.0, 1.5)        # ions: no flux focusing, <= open field
    m_F = np.clip(m_F, 0.0, 4.0)        # neutrals: mild enhancement via re-emission
    m_O = np.clip(m_O, 0.0, 4.0)
    cos_i = np.where(fi>0, ai/np.maximum(fi,1e-9), 0.0)   # mean incidence cosine for ions
    return m_i, m_F, m_O, cos_i

# ----------------------------- surface chemistry -> normal velocity -----------------------------
def surface_rate(m_i, m_F, m_O, cos_i, is_mask, par):
    E = par['Emean']
    sqrtE = np.sqrt(max(E,0.0))
    Yie = par['A_ie']*max(sqrtE-np.sqrt(par['Eth_ie']),0.0)
    Ysp = par['A_sp']*max(sqrtE-np.sqrt(par['Eth_sp']),0.0)
    Yp  = par['A_p'] *max(sqrtE-np.sqrt(par['Eth_p']),0.0)
    # angular factor for ion yields (forward-peaked; ~cos incidence)
    fang = np.clip(cos_i, 0.0, 1.0)
    Fi = par['ionFlux']*m_i*fang
    Fev = par['Fflux']*m_F
    Fp  = par['Oflux']*m_O
    eps = 1e-9
    # competitive Langmuir steady-state coverages
    rF = par['s_F']*Fev / (Yie*Fi + eps)        # theta_F / bare
    rO = par['s_O']*Fp  / (Yp*Fi + eps)         # theta_O / bare
    bare = 1.0/(1.0 + rF + rO)
    thF = rF*bare
    V = (1.0/par['rho'])*(Yie*Fi*thF + Ysp*Fi*bare)   # Si removal -> normal velocity
    V = V*par['rate_scale']
    V[is_mask] = 0.0                                    # mask not etched
    return V

# ----------------------------- level-set advection (Godunov upwind) -----------------------------
def advect(phi, F, dx, dt):
    """phi_t + F|grad phi| = 0, F>=0 (etch shrinks solid)."""
    dxm = np.zeros_like(phi); dxp = np.zeros_like(phi)
    dym = np.zeros_like(phi); dyp = np.zeros_like(phi)
    dxm[1:,:]  = (phi[1:,:]-phi[:-1,:])/dx
    dxp[:-1,:] = (phi[1:,:]-phi[:-1,:])/dx
    dym[:,1:]  = (phi[:,1:]-phi[:,:-1])/dx
    dyp[:,:-1] = (phi[:,1:]-phi[:,:-1])/dx
    grad = np.sqrt(np.maximum(dxm,0)**2 + np.minimum(dxp,0)**2 +
                   np.maximum(dym,0)**2 + np.minimum(dyp,0)**2)
    return phi - dt*F*grad

def extend_velocity(V, mid, phi, xs, ys, dx, band):
    """Extend surface velocities to grid narrow-band via nearest-segment lookup."""
    nx, ny = phi.shape
    F = np.zeros_like(phi)
    bandmask = np.abs(phi) < band
    ii, jj = np.where(bandmask)
    gx = xs[0]+ii*dx; gy = ys[0]+jj*dx
    if len(mid)==0:
        return F
    tree = cKDTree(mid)
    _, idx = tree.query(np.column_stack([gx, gy]))
    F[ii, jj] = V[idx]
    return F

# ----------------------------- main driver -----------------------------
def run_etch(W=20.0, H=20.0, dx=0.25, trench_width=8.0, mask_thickness=2.0,
             sub_top=15.0, t_end=3.0, n_steps=60, par=None, snapshots=None, verbose=True):
    if par is None: par = PAR
    X,Y,xs,ys,phi,mask,nx,ny = make_trench(W,H,dx,trench_width,mask_thickness,sub_top)
    mask_phi = phi.copy()                  # to re-stamp mask each step
    y_src = H - dx
    dt = t_end/n_steps
    band = 5*dx
    timings = dict(raytrace=0.0, total=0.0, chem=0.0, advect=0.0)
    snaps = {}
    t_total0 = time.time()
    for step in range(n_steps):
        segs, mid, nrm, L = extract_surface(phi, xs, ys, dx)
        if len(segs)==0: break
        nrm = orient_normals(mid, nrm, phi, xs, ys, dx)
        is_mask = _seg_in_mask(mid, mask, xs, ys, dx)
        t0=time.time()
        m_i,m_F,m_O,cos_i = mc_flux(segs, mid, nrm, is_mask, L, y_src, W, par)
        timings['raytrace'] += time.time()-t0
        t0=time.time()
        V = surface_rate(m_i,m_F,m_O,cos_i,is_mask,par)
        timings['chem'] += time.time()-t0
        # CFL
        F = extend_velocity(V, mid, phi, xs, ys, dx, band)
        t0=time.time()
        Vmax = max(V.max(), 1e-6)
        nsub = int(np.ceil(Vmax*dt/(0.4*dx)))
        nsub = max(1, min(nsub, 40))
        for _ in range(nsub):
            phi = advect(phi, F, dx, dt/nsub)
            phi[mask] = mask_phi[mask]
        phi = skfmm.distance(phi, dx=dx)                   # reinit to signed distance
        timings['advect'] += time.time()-t0
        if snapshots and step in snapshots:
            snaps[step] = extract_surface(phi, xs, ys, dx)[0]
        if verbose and step%10==0:
            depth = sub_top - _profile_bottom(phi, xs, ys, dx, W)
            print(f"  step {step:3d}/{n_steps}  etch depth ~ {depth:5.2f} um  Vmax {Vmax:.3f}")
    timings['total'] = time.time()-t_total0
    final_segs = extract_surface(phi, xs, ys, dx)[0]
    return dict(phi=phi, xs=xs, ys=ys, dx=dx, segs=final_segs, snaps=snaps,
                timings=timings, sub_top=sub_top, mask=mask, X=X, Y=Y)

def _seg_in_mask(mid, mask, xs, ys, dx):
    ix = np.clip(((mid[:,0]-xs[0])/dx).astype(int), 0, mask.shape[0]-1)
    iy = np.clip(((mid[:,1]-ys[0])/dx).astype(int), 0, mask.shape[1]-1)
    return mask[ix, iy]

def _profile_bottom(phi, xs, ys, dx, W):
    """Lowest y of the phi=0 contour near the trench centre."""
    segs = extract_surface(phi, xs, ys, dx)[0]
    if len(segs)==0: return ys[-1]
    cx = np.r_[segs[:,0], segs[:,2]]; cy = np.r_[segs[:,1], segs[:,3]]
    centre = np.abs(cx - W/2) < W*0.15
    return cy[centre].min() if centre.any() else cy.min()

if __name__ == "__main__":
    print("warm-up JIT...")
    _ = run_etch(W=12, H=12, dx=0.5, trench_width=5, sub_top=9, t_end=0.5, n_steps=3, verbose=False)
    print("done warm-up")
