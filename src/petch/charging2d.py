"""2-D feature-charging solver — the Hwang-Giapis method, minimally implemented.

The 0-D current-balance closure (charging.py) gets the rolloff SHAPE but not the magnitude:
it only decelerates ions vertically, while HG's 2-D field also DEFLECTS ions laterally into the
sidewalls before they reach the floor, and the floor potential is strongly non-uniform (center
~33 V vs ~59 V at the sidewall foot at AR 4). This module does the published thing:

  loop:  trace ions + electrons ballistically through E = -grad(V)
         accumulate net current per surface segment
         relax segment potentials  V_j += k * (I_ion - I_e)_j   (insulating surfaces)
         solve Laplace in the vacuum with Dirichlet surfaces + V=0 at the sheath boundary
  until net segment currents vanish  ->  steady state (HG: "equal positive and negative
  currents to all surface segments").

Simplifications vs HG (documented): relaxation on segment potentials instead of explicit charge
accounting; no RF-phase-resolved electron bursts (electrons sampled from the time-averaged
Maxwellian flux); no poly-Si equipotential (all surfaces insulating); 2-D trench cross-section.
Particle sources per HG JAP 82,566: ions bimodal IEDF (V_dc + V_rf cos phi), Gaussian IADF
(HWHM 4.3 deg); electrons Maxwellian T_e, EADF ~ cos^0.6.

Units: potentials in volts; particle energies in eV; masses folded so v = sqrt(E) per component.
Grid: x in [0, W_cells + 2*pad], z in [0, depth + mouth]; trench walls/floor = surface segments.
"""
from __future__ import annotations

import numpy as np


def solve_trench_charging(AR, W=32, pad=24, mouth=24, Te=4.0, V_dc=37.0, V_rf=30.0,
                          iadf_hwhm_deg=4.3, cos_power=0.6, n_per_iter=6000, n_iter=120,
                          relax=None, seed=0, verbose=False, smooth=False):
    """Steady-state charging of a 2-D insulating trench. Returns dict with:
    floor_flux (normalized ion flux reaching the floor), V_floor_center, V_foot_peak, V (grid).

    smooth=False (default): raw tail-averaged segment potentials. smooth=True applies a
    PRESENTATION-grade uniform filter to the converged segment potentials before the final
    Laplace solve of the returned map -- purely cosmetic (removes per-segment MC shot-noise
    raggedness from the walls); the gate quantities (floor_flux, V_floor_center, V_foot_peak)
    are always computed from the UNsmoothed tail averages, and scripts/charging_gate.py must
    always run with smooth=False."""
    rng = np.random.default_rng(seed)
    D = int(round(AR * W))
    nx = W + 2 * pad
    nz = D + mouth
    # geometry masks: vacuum interior; solid = everything below z=mouth outside the trench slot
    solid = np.zeros((nx, nz), bool)
    solid[:pad, mouth:] = True
    solid[pad + W:, mouth:] = True
    solid[:, nz - 1] = True                            # floor substrate line
    V = np.zeros((nx, nz))
    # surface segments: floor cells + two sidewalls (store index lists)
    floor_ix = np.arange(pad, pad + W)
    side_z = np.arange(mouth, nz - 1)
    if relax is None:
        relax = 2.0 * Te                               # volts per unit normalized net current

    Vfloor = np.zeros(W)                               # floor segment potentials
    Vleft = np.zeros(len(side_z)); Vright = np.zeros(len(side_z))
    Vtop_l = 0.0; Vtop_r = 0.0                         # top (mask) surfaces, one segment each

    def apply_dirichlet(V):
        V[:, 0] = 0.0                                  # sheath boundary
        V[floor_ix, nz - 1] = Vfloor
        V[pad - 1, mouth:nz - 1] = Vleft               # sidewall faces (solid boundary cells)
        V[pad + W, mouth:nz - 1] = Vright
        V[:pad, mouth] = Vtop_l                        # mask top surfaces
        V[pad + W:, mouth] = Vtop_r
        return V

    # red-black Gauss-Seidel SOR (Jacobi diverges for omega>1; np.roll would wrap boundaries)
    ii, jj = np.meshgrid(np.arange(nx), np.arange(nz), indexing='ij')
    _red = ((ii + jj) % 2 == 0)
    _inside = ~solid
    _inside[:, 0] = False
    _inside[0, :] = False; _inside[-1, :] = False; _inside[:, -1] = False

    def laplace(V, sweeps=120, omega=1.8):
        for _ in range(sweeps):
            V = apply_dirichlet(V)
            for color in (_red, ~_red):
                m = _inside & color
                avg = np.zeros_like(V)
                avg[1:-1, 1:-1] = 0.25 * (V[2:, 1:-1] + V[:-2, 1:-1] + V[1:-1, 2:] + V[1:-1, :-2])
                V[m] = (1 - omega) * V[m] + omega * avg[m]
        return apply_dirichlet(V)

    def trace(kind, n):
        """Ballistic leapfrog through E=-grad(V). Returns surface hit tallies."""
        # HG's IADF/EADF are quoted from THEIR 2-D simulation plane -> sample the in-plane angle
        # directly (no 3-D sin(theta) Jacobian, no azimuthal projection).
        if kind == 'ion':
            E0 = V_dc + V_rf * np.cos(rng.uniform(0, np.pi, n))
            sig = np.deg2rad(iadf_hwhm_deg) / 1.1774
            th = rng.normal(0.0, sig, n)               # in-plane angle (signed)
            q = +1.0
        else:
            Ee = rng.gamma(2.0, Te, n)
            # in-plane EADF ~ cos^p(theta) on (-pi/2, pi/2): numeric inverse CDF
            tg = np.linspace(-np.pi / 2 + 1e-4, np.pi / 2 - 1e-4, 512)
            cdf = np.cumsum(np.cos(tg) ** cos_power); cdf /= cdf[-1]
            th = np.interp(rng.uniform(0, 1, n), cdf, tg)
            E0 = Ee
            q = -1.0
        vx = np.sqrt(E0) * np.sin(th)
        vz = np.sqrt(E0) * np.abs(np.cos(th))
        x = rng.uniform(0, nx - 1.0, n)
        z = np.ones(n) * 1.0
        alive = np.ones(n, bool)
        hits_floor = np.zeros(W)
        hits_left = np.zeros(len(side_z)); hits_right = np.zeros(len(side_z))
        hit_top_l = 0.0; hit_top_r = 0.0
        foot_n = 0.0; foot_E = 0.0                     # sidewall-FOOT impacts (bottom 15% of depth)
        z_foot = mouth + 0.85 * D
        Ex = -(np.gradient(V, axis=0)); Ez = -(np.gradient(V, axis=1))
        for _ in range(int(14 * nz)):
            if not alive.any():
                break
            ix = np.clip(x[alive].astype(int), 0, nx - 2)
            iz = np.clip(z[alive].astype(int), 0, nz - 2)
            ax = q * Ex[ix, iz] * 0.5
            az = q * Ez[ix, iz] * 0.5
            vmax = np.maximum(np.abs(vx[alive]), np.abs(vz[alive]))
            dt = 0.45 / np.maximum(vmax, 0.8)          # per-particle: <=0.45 cell per step
            vx[alive] += ax * dt
            vz[alive] += az * dt
            x[alive] += vx[alive] * dt
            z[alive] += vz[alive] * dt
            ia = np.flatnonzero(alive)
            xi, zi = x[ia], z[ia]
            out = (zi < 0.5) | (xi < 0.5) | (xi > nx - 1.5)
            gone = ia[out]
            alive[gone] = False
            zi_i = zi.astype(int)
            # floor hit
            fh = (zi >= nz - 1.5)
            for j in ia[fh]:
                c = int(np.clip(x[j] - pad, 0, W - 1))
                if 0 <= x[j] - pad < W:
                    hits_floor[c] += 1
                alive[j] = False
            # mask top hits (z crossing mouth outside slot)
            th_hit = (zi >= mouth) & (zi < nz - 1.5) & ((xi < pad) | (xi >= pad + W))
            for j in ia[th_hit]:
                zj = int(np.clip(z[j] - mouth, 0, len(side_z) - 1))
                if z[j] < mouth + 1.5:                 # top surface
                    if x[j] < pad:
                        hit_top_l += 1
                    else:
                        hit_top_r += 1
                else:                                   # sidewall (entered solid laterally)
                    if x[j] < pad + W / 2:
                        hits_left[zj] += 1
                    else:
                        hits_right[zj] += 1
                    if z[j] >= z_foot:                  # deflected into the foot (notch driver)
                        foot_n += 1.0
                        foot_E += vx[j] ** 2 + vz[j] ** 2   # impact energy (v^2 = E in eV units)
                alive[j] = False
        return hits_floor, hits_left, hits_right, hit_top_l, hit_top_r, foot_n, foot_E

    hist = []
    vfloor_hist = []
    vleft_hist = []; vright_hist = []
    for it in range(n_iter):
        V = laplace(V)
        fi, li, ri, tli, tri, _, _ = trace('ion', n_per_iter)
        fe, le, re, tle, tre, _, _ = trace('electron', n_per_iter)
        # normalized net current per segment (per source particle); anneal the step so early
        # transients move fast and the steady state stops random-walking on shot noise
        anneal = max(1.0 / (1.0 + it / 25.0), 0.25)   # floor: late iters can still unpin clips
        scale = anneal * relax / n_per_iter * (nx)     # keep step size geometry-independent
        Vfloor += scale * (fi - fe)
        Vleft += scale * (li - le)
        Vright += scale * (ri - re)
        Vtop_l += scale * (tli - tle) / max(pad, 1)
        Vtop_r += scale * (tri - tre) / max(pad, 1)
        np.clip(Vfloor, 0.0, V_dc + V_rf, out=Vfloor)
        np.clip(Vleft, -3 * Te, V_dc + V_rf, out=Vleft)
        np.clip(Vright, -3 * Te, V_dc + V_rf, out=Vright)
        Vtop_l = float(np.clip(Vtop_l, -3 * Te, 0.0))
        Vtop_r = float(np.clip(Vtop_r, -3 * Te, 0.0))
        hist.append(fi.sum() / n_per_iter)
        vfloor_hist.append(Vfloor.copy())
        vleft_hist.append(Vleft.copy()); vright_hist.append(Vright.copy())
        if verbose and it % 20 == 0:
            print(f"  it{it}: floor_i={fi.sum()/n_per_iter:.3f} floor_e={fe.sum()/n_per_iter:.3f} "
                  f"Vc={Vfloor[W//2]:.1f} Vfoot={max(Vfloor[0], Vfloor[-1]):.1f}", flush=True)
    # steady state: average flux AND segment potentials over the last third (shot-noise suppression)
    k = max(n_iter // 3, 5)
    tail = np.mean(hist[-k:])
    Vf_avg = np.mean(np.array(vfloor_hist[-k:]), axis=0)
    Vl_avg = np.mean(np.array(vleft_hist[-k:]), axis=0)
    Vr_avg = np.mean(np.array(vright_hist[-k:]), axis=0)
    # final map: re-solve Laplace with the tail-averaged (optionally smoothed) segment potentials
    Vf_map, Vl_map, Vr_map = Vf_avg, Vl_avg, Vr_avg
    if smooth:                                         # cosmetic only -- see docstring
        from scipy.ndimage import uniform_filter1d
        Vf_map = uniform_filter1d(Vf_avg, 5, mode="nearest")
        Vl_map = uniform_filter1d(Vl_avg, 9, mode="nearest")
        Vr_map = uniform_filter1d(Vr_avg, 9, mode="nearest")
    Vfloor[:] = Vf_map; Vleft[:] = Vl_map; Vright[:] = Vr_map
    V = laplace(V, sweeps=240)
    # deflected-ion foot statistics at the converged field (one clean high-stat ion trace):
    # the sub-threshold IEDF slice bent into the sidewall foot -- the notching driver
    fi2, _, _, _, _, fn2, fE2 = trace('ion', 4 * n_per_iter)
    open_frac = W / nx                                 # fraction of source over the slot
    return dict(floor_flux=float(tail / open_frac), V_floor_center=float(Vf_avg[W // 2]),
                V_foot_peak=float(Vf_avg.max()), Vfloor=Vf_avg, V=V,
                foot_ion_flux=float(fn2 / (4 * n_per_iter) / open_frac),
                foot_ion_Emean=float(fE2 / max(fn2, 1.0)))


# Gate-validated reference curve (scripts/charging_gate.py, RMSE 0.039 vs Hwang-Giapis JAP 82,566
# Fig. 4; HG conditions: Cl2 HDP, V_s = 37+30 sin(wt), T_e = 4 eV). Model values at the HG ARs.
_GATE_AR = np.array([0.0, 1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
_GATE_FLUX = np.array([1.0, 0.648, 0.599, 0.504, 0.433, 0.334, 0.283, 0.213, 0.177])
_GATE_VFLOOR = np.array([0.0, 13.4, 16.8, 24.4, 30.9, 39.2, 44.9, 48.5, 52.9])
# Deflected-ion mean impact energy on the sidewall foot vs AR (eV). HG JAP 82,566 tabulated
# values (15 -> 27.5 eV over AR 1 -> 4); the notching-mechanism gate (scripts/notching_gate.py)
# validates our solver against these. AR=0 row: no charging, no deflection.
_FOOT_E = np.array([0.0, 15.0, 16.5, 17.5, 20.0, 23.0, 25.0, 26.5, 27.5])


def charging_floor_profile(AR):
    """Production hook (NOT yet wired into the flux pipeline): normalized floor ion flux and
    floor potential vs aspect ratio, from the gate-validated 2-D solver at the HG reference
    conditions. Returns (floor_flux_factor, V_floor, E_deflected) vs AR. Interpolates the gate
    curve; beyond AR 4 clamps (re-run solve_trench_charging for other plasma conditions).
    Applies to INSULATING floors (SiO2/SOI overetch, dielectric etch); a conductive grounded
    Si floor drains and should NOT be throttled by this (the de Boer deep-Si case)."""
    AR = np.asarray(AR, float)
    return (np.interp(AR, _GATE_AR, _GATE_FLUX),
            np.interp(AR, _GATE_AR, _GATE_VFLOOR),
            np.interp(AR, _GATE_AR, _FOOT_E))
