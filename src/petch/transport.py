"""Monte Carlo ballistic flux transport. Ported verbatim from feature_etch.py.

`_trace` is the numba hot loop (the ~94%-of-runtime ray kernel that the GPU/Warp port
targets). `mc_flux` launches ions + F + O species against the frozen surface. The default
path is byte-identical to the original PoC; variant hooks (QMC sampling, clip removal,
2D re-emission, ion/neutral split) are wired in later steps.
"""
import warnings
from contextlib import contextmanager
import numpy as np
from numba import njit
from scipy.stats import qmc, norm


@contextmanager
def _suppress_qmc_warning():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*power of 2.*")
        yield


def _sobol_source(n, W, sigma, kind, seed):
    """Low-discrepancy (Sobol) source launch: stratified (x-position, angle) per particle.

    QMC over the *source* dimensions — the dominant variance source — drives variance toward
    ~1/N instead of ~1/sqrt(N). Bounce decisions inside _trace stay pseudorandom.
    kind: 'ion' (Gaussian angle, std=sigma) or 'neutral' (2D-cosine via arcsin).
    """
    eng = qmc.Sobol(d=2, scramble=True, seed=seed)
    with np.errstate(all='ignore'), _suppress_qmc_warning():
        u = eng.random(n)                   # (n, 2) in [0,1); non-pow2 n is fine (scrambled)
    u0 = u[:, 0]; u1 = np.clip(u[:, 1], 1e-9, 1 - 1e-9)
    xs = u0 * W
    if kind == 'ion':
        ang = norm.ppf(u1) * sigma
    else:
        ang = np.arcsin(2.0 * u1 - 1.0)     # inverse-CDF of the 2D cosine launch
    return xs, np.sin(ang), -np.cos(ang)


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


def mc_flux(seg, mid, nrm, is_mask, L, y_src, W, par, n_part_ion=20000, n_part_neu=20000,
            seed=0, n_reemit=12, sampling="pseudo"):
    """Compute per-segment normalized flux multipliers + mean ion incidence cos for 3 species.

    `seed` selects an independent Monte-Carlo realization (seed=0, sampling='pseudo' reproduces
    the PoC exactly). `n_reemit` is the neutral re-emission bounce cap. `sampling`: 'pseudo'
    (PoC) or 'sobol' (QMC source launch — variance ~1/N instead of ~1/sqrt(N)).
    """
    s_off = 3 * seed
    if sampling == "sobol":
        xs0, dix, diy = _sobol_source(n_part_ion, W, par['ion_ang_sigma'], 'ion', 10 + s_off)
        xs1, dfx, dfy = _sobol_source(n_part_neu, W, None, 'neutral', 11 + s_off)
        xs2, dox, doy = _sobol_source(n_part_neu, W, None, 'neutral', 12 + s_off)
    else:  # pseudo — exact PoC behavior (one shared rng, in species order)
        rng = np.random.default_rng(seed)
        xs0 = rng.uniform(0, W, n_part_ion)
        a = rng.normal(0, par['ion_ang_sigma'], n_part_ion)
        dix = np.sin(a); diy = -np.cos(a)
        xs1 = rng.uniform(0, W, n_part_neu)
        aF = np.arcsin(rng.uniform(-1, 1, n_part_neu))
        dfx = np.sin(aF); dfy = -np.cos(aF)
        xs2 = rng.uniform(0, W, n_part_neu)
        aO = np.arcsin(rng.uniform(-1, 1, n_part_neu))
        dox = np.sin(aO); doy = -np.cos(aO)
    fi, ai = _trace(seg, nrm, is_mask, xs0, y_src, dix, diy, 1.0, 0, 1 + s_off)
    fF, _ = _trace(seg, nrm, is_mask, xs1, y_src, dfx, dfy, par['s_F'], n_reemit, 2 + s_off)
    fO, _ = _trace(seg, nrm, is_mask, xs2, y_src, dox, doy, par['s_O'], n_reemit, 3 + s_off)
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
