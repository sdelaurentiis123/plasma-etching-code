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
                          relax=None, seed=0, verbose=False):
    """Steady-state charging of a 2-D insulating trench. Returns dict with:
    floor_flux (normalized ion flux reaching the floor), V_floor_center, V_foot_peak, V (grid)."""
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
        if kind == 'ion':
            E0 = V_dc + V_rf * np.cos(rng.uniform(0, np.pi, n))
            sig = np.deg2rad(iadf_hwhm_deg) / 1.1774
            th = rng.normal(0.0, sig, n)
            q = +1.0
        else:
            Ee = rng.gamma(2.0, Te, n)
            u = rng.uniform(0, 1, n)
            ct = (1 - u) ** (1.0 / (cos_power + 1.0))
            th = np.arccos(ct) * np.where(rng.uniform(0, 1, n) < 0.5, 1, -1)
            E0 = Ee
            q = -1.0
        vx = np.sqrt(E0) * np.sin(th)
        vz = np.sqrt(E0) * np.cos(th)
        x = rng.uniform(0, nx - 1.0, n)
        z = np.ones(n) * 1.0
        alive = np.ones(n, bool)
        hits_floor = np.zeros(W)
        hits_left = np.zeros(len(side_z)); hits_right = np.zeros(len(side_z))
        hit_top_l = 0.0; hit_top_r = 0.0
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
                alive[j] = False
        return hits_floor, hits_left, hits_right, hit_top_l, hit_top_r

    hist = []
    for it in range(n_iter):
        V = laplace(V)
        fi, li, ri, tli, tri = trace('ion', n_per_iter)
        fe, le, re, tle, tre = trace('electron', n_per_iter)
        # normalized net current per segment (per source particle)
        scale = relax / n_per_iter * (nx)              # keep step size geometry-independent
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
        if verbose and it % 20 == 0:
            print(f"  it{it}: floor_i={fi.sum()/n_per_iter:.3f} floor_e={fe.sum()/n_per_iter:.3f} "
                  f"Vc={Vfloor[W//2]:.1f} Vfoot={max(Vfloor[0], Vfloor[-1]):.1f}", flush=True)
    # steady-state floor ion flux: average of the last third, normalized by the open fraction
    tail = np.mean(hist[-max(n_iter // 3, 5):])
    open_frac = W / nx                                 # fraction of source over the slot
    return dict(floor_flux=float(tail / open_frac), V_floor_center=float(Vfloor[W // 2]),
                V_foot_peak=float(max(Vfloor[0], Vfloor[-1])), Vfloor=Vfloor.copy(), V=V)
