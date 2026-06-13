"""Monte Carlo ballistic flux transport. Ported verbatim from feature_etch.py.

`_trace` is the numba hot loop (the ~94%-of-runtime ray kernel that the GPU/Warp port
targets). `mc_flux` launches ions + F + O species against the frozen surface. The default
path is byte-identical to the original PoC; variant hooks (QMC sampling, clip removal,
2D re-emission, ion/neutral split) are wired in later steps.
"""
import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def _trace(seg, nrm, is_mask, x0_src, y_src, dirs_x, dirs_y, sticking, n_reemit, rng_seed):
    """Forward-trace one species' particles. Returns per-segment hit weight + ion-angle accum."""
    M = seg.shape[0]
    flux = np.zeros(M)
    ang_acc = np.zeros(M)   # accumulated cos(incidence) weight (for ions)
    N = x0_src.shape[0]
    np.random.seed(rng_seed)
    for p in range(N):
        px = x0_src[p]; py = y_src
        dx_ = dirs_x[p]; dy_ = dirs_y[p]
        w = 1.0
        for bounce in range(n_reemit + 1):
            best_t = 1e18; best_s = -1
            for s in range(M):
                ax = seg[s, 0]; ay = seg[s, 1]; bx = seg[s, 2]; by = seg[s, 3]
                ex = bx - ax; ey = by - ay
                den = dx_ * ey - dy_ * ex
                if abs(den) < 1e-14:
                    continue
                t = ((ax - px) * ey - (ay - py) * ex) / den
                u = ((ax - px) * dy_ - (ay - py) * dx_) / den
                if t > 1e-6 and u >= -1e-6 and u <= 1.0 + 1e-6:
                    if t < best_t:
                        best_t = t; best_s = s
            if best_s < 0:
                break  # escaped domain
            s = best_s
            nx_ = nrm[s, 0]; ny_ = nrm[s, 1]
            cosang = -(dx_ * nx_ + dy_ * ny_)
            if cosang < 0:  # hit backface; treat as grazing
                cosang = 0.0
            if np.random.random() < sticking or bounce == n_reemit:
                flux[s] += w
                ang_acc[s] += w * cosang
                break
            else:
                # diffuse re-emit (cosine about outward normal) and continue
                hx = px + best_t * dx_; hy = py + best_t * dy_
                tx = -ny_; ty = nx_
                ca = np.sqrt(np.random.random())          # cosine law
                sa = np.sqrt(1.0 - ca * ca)
                sign = 1.0 if np.random.random() < 0.5 else -1.0
                dx_ = ca * nx_ + sign * sa * tx
                dy_ = ca * ny_ + sign * sa * ty
                px = hx + 1e-4 * nx_; py = hy + 1e-4 * ny_
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
    aF = np.arcsin(rng.uniform(-1, 1, n_part_neu))     # cosine-ish into lower hemisphere
    dfx = np.sin(aF); dfy = -np.cos(aF)
    fF, _ = _trace(seg, nrm, is_mask, xs1, y_src, dfx, dfy, par['s_F'], 12, 2)
    # --- O passivation ---
    xs2 = rng.uniform(0, W, n_part_neu)
    aO = np.arcsin(rng.uniform(-1, 1, n_part_neu))
    dox = np.sin(aO); doy = -np.cos(aO)
    fO, _ = _trace(seg, nrm, is_mask, xs2, y_src, dox, doy, par['s_O'], 12, 3)
    # normalize to open-field flux density (particles per unit x length)
    base_ion = n_part_ion / W
    base_neu = n_part_neu / W
    m_i = (fi / np.maximum(L, 0.3 * np.median(L))) / base_ion
    m_F = (fF / np.maximum(L, 0.3 * np.median(L))) / base_neu
    m_O = (fO / np.maximum(L, 0.3 * np.median(L))) / base_neu
    m_i = np.clip(m_i, 0.0, 1.5)        # ions: no flux focusing, <= open field
    m_F = np.clip(m_F, 0.0, 4.0)        # neutrals: mild enhancement via re-emission
    m_O = np.clip(m_O, 0.0, 4.0)
    cos_i = np.where(fi > 0, ai / np.maximum(fi, 1e-9), 0.0)   # mean incidence cosine for ions
    return m_i, m_F, m_O, cos_i
