"""2-D feature-charging solver — the Hwang-Giapis method with the poly-Si conductor.

The 0-D current-balance closure (charging.py) gets the rolloff SHAPE but not the magnitude:
it only decelerates ions vertically, while HG's 2-D field also DEFLECTS ions laterally into the
sidewalls before they reach the floor. This module does the published thing (HG JVST B 15, 70;
JAP 82, 566):

  loop:  trace ions + electrons ballistically through the in-feature field
         accumulate net current per surface segment (insulators) / per conductor (poly lines)
         relax potentials  V_j += k * (I_ion - I_e)_j ; conductors updated as ONE equipotential
         solve Laplace in the vacuum with Dirichlet surfaces + the sheath-edge boundary
  until net currents vanish  ->  steady state (HG: "equal positive and negative currents to
  all surface segments").

Geometry (HG JAP 82,566): photoresist mask (insulating, upper sidewall) over a 0.3 µm n+
poly-Si line (CONDUCTOR — one floating equipotential per line, explicit charge redistribution)
on gate oxide (insulating floor). Line/space 0.5/0.5 µm; PR thickness sets AR = 1 -> 4. In the
periodic array both walls' poly lines are equivalent -> tied to one V_poly.

Potential referencing (all potentials vs the GROUNDED SUBSTRATE, like HG):
  V_total(x,z; phi) = V_A(x,z) + V_s(phi) * V_B(x,z)      (Laplace is linear)
  V_A : surfaces at their ground-referenced potentials, sheath-edge boundary at 0
  V_B : surfaces at 0, sheath-edge boundary at 1;  V_s(phi) = V_dc + V_rf*sin(phi)
Species treatment (documented approximation):
  IONS enter with the measured wafer-arrival bimodal IEDF (E = V_s(phi), phi ~ U) and the
  arrival IADF (HWHM 4.3 deg) — their sheath fall is already folded into that arrival energy,
  so they see only the charging perturbation E_A (adding V_s*E_B would double-count the
  acceleration that produced their arrival KE).
  ELECTRONS must face the instantaneous sheath barrier explicitly: they arrive burst-weighted
  at sheath-voltage minima (flux ~ exp(-V_s(phi)/Te)) with Maxwellian residual energy, and see
  E_A + V_s(phi_e)*E_B — the V_B climb is the barrier; positively charged surfaces offset it
  locally (HG's electrostatic softening of electron shading).

Units: potentials in volts (vs ground); particle energies in eV; v = sqrt(E) per component.
"""
from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange
except Exception:  # pragma: no cover - optional acceleration
    njit = None
    prange = range


if njit is not None:
    @njit(cache=True, parallel=True, fastmath=True)
    def _trace_particles_adaptive(Ex, Ez, ExB, EzB, x0, vx0, vz0, sB, q,
                                  nx, nz, W, pad, mouth, n_side, z_poly0,
                                  max_steps):
        n = x0.shape[0]
        hit_type = np.zeros(n, np.int8)      # 0 escape/survivor, 1 floor, 2 left, 3 right, 4 top_l, 5 top_r
        hit_idx = np.full(n, -1, np.int64)
        foot_E = np.zeros(n)
        survivor = np.zeros(n, np.uint8)
        steps = np.zeros(n, np.int64)
        xmax = float(nx)

        for p in prange(n):
            x = x0[p]
            z = 1.0
            vx = vx0[p]
            vz = vz0[p]
            sb = sB[p]
            alive = True
            last_step = 0
            for st in range(max_steps):
                last_step = st + 1
                ix = int(x)
                if ix < 0:
                    ix = 0
                elif ix > nx - 2:
                    ix = nx - 2
                iz = int(z)
                if iz < 0:
                    iz = 0
                elif iz > nz - 2:
                    iz = nz - 2

                fx = Ex[ix, iz] + sb * ExB[ix, iz]
                fz = Ez[ix, iz] + sb * EzB[ix, iz]
                ax = q * fx * 0.5
                az = q * fz * 0.5
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

                ix2 = int(xa)
                if ix2 < 0:
                    ix2 = 0
                elif ix2 > nx - 2:
                    ix2 = nx - 2
                iz2 = int(za)
                if iz2 < 0:
                    iz2 = 0
                elif iz2 > nz - 2:
                    iz2 = nz - 2
                fx2 = Ex[ix2, iz2] + sb * ExB[ix2, iz2]
                fz2 = Ez[ix2, iz2] + sb * EzB[ix2, iz2]
                vx = vx_half + 0.25 * q * fx2 * dt
                vz = vz_half + 0.25 * q * fz2 * dt
                x = xa
                z = za

                # periodic pitch
                x = x % xmax

                if z < 0.5:
                    alive = False
                    break
                if z >= nz - 1.5:
                    c = int(x - pad)
                    if c >= 0 and c < W:
                        hit_type[p] = 1
                        hit_idx[p] = c
                    alive = False
                    break
                if z >= mouth and z < nz - 1.5 and (x < pad or x >= pad + W):
                    zj = int(z - mouth)
                    if zj < 0:
                        zj = 0
                    elif zj > n_side - 1:
                        zj = n_side - 1
                    if z < mouth + 1.5:
                        hit_type[p] = 4 if x < pad else 5
                    else:
                        if x < pad + 0.5 * W:
                            hit_type[p] = 2
                        else:
                            hit_type[p] = 3
                        hit_idx[p] = zj
                        if q > 0.0 and z >= z_poly0:
                            foot_E[p] = vx * vx + vz * vz
                    alive = False
                    break
            if alive:
                survivor[p] = 1
            steps[p] = last_step
        return hit_type, hit_idx, foot_E, survivor, steps
else:
    _trace_particles_adaptive = None


def solve_trench_charging(AR, W=32, pad=16, mouth=237, Te=4.0, V_dc=37.0, V_rf=30.0,
                          iadf_hwhm_deg=4.3, cos_power=0.6, n_per_iter=6000, n_iter=120,
                          relax=None, seed=0, verbose=False, smooth=False,
                          poly_um=0.3, feature_w_um=0.5, rf_bursts=True,
                          sheath_um=89.0, boundary_um=3.7, insul_vmin_Te=1.0,
                          trace_integrator="adaptive_numba", trace_step_cap_factor=40.0):
    """Steady-state charging of the HG poly-on-oxide trench. Returns dict with:
    floor_flux (normalized ion flux to the oxide floor), V_floor_center, V_foot_peak,
    V_poly (the poly-line equipotential), foot_ion_flux / foot_ion_Emean (ions striking the
    poly sidewall — the notch driver), V (grid map, ground-referenced).

    poly_um=0 recovers the all-insulator solver. rf_bursts=False recovers time-averaged
    electrons with no sheath barrier (the pre-burst model). smooth=True is PRESENTATION-grade
    filtering of the returned map only — gate quantities always computed raw
    (scripts/charging_gate.py and notching_gate.py must run smooth=False)."""
    rng = np.random.default_rng(seed)
    D = int(round(AR * W))
    nx = W + 2 * pad
    nz = D + mouth
    poly_cells = int(round(poly_um / feature_w_um * W)) if poly_um > 0 else 0
    poly_cells = min(poly_cells, max(D - 2, 0))
    solid = np.zeros((nx, nz), bool)
    solid[:pad, mouth:] = True
    solid[pad + W:, mouth:] = True
    solid[:, nz - 1] = True                            # floor substrate line
    V = np.zeros((nx, nz))
    floor_ix = np.arange(pad, pad + W)
    side_z = np.arange(mouth, nz - 1)
    n_side = len(side_z)
    z_poly0 = nz - 1 - poly_cells                      # poly band: bottom poly_cells of the wall
    is_poly = side_z >= z_poly0
    if relax is None:
        relax = 2.0 * Te

    Vfloor = np.zeros(W)
    Vleft = np.zeros(n_side); Vright = np.zeros(n_side)
    Vpoly = 0.0                                        # one floating equipotential (periodic array)
    Vtop_l = 0.0; Vtop_r = 0.0

    def apply_dirichlet(V, boundary=0.0, vf=None, vl=None, vr=None, vp=None, vtl=None, vtr=None):
        V[:, 0] = boundary
        V[floor_ix, nz - 1] = Vfloor if vf is None else vf
        vl_ = (Vleft if vl is None else vl).copy()
        vr_ = (Vright if vr is None else vr).copy()
        vp_ = Vpoly if vp is None else vp
        if poly_cells > 0:
            vl_[is_poly] = vp_
            vr_[is_poly] = vp_
        V[pad - 1, mouth:nz - 1] = vl_
        V[pad + W, mouth:nz - 1] = vr_
        V[:pad, mouth] = Vtop_l if vtl is None else vtl
        V[pad + W:, mouth] = Vtop_r if vtr is None else vtr
        return V

    ii, jj = np.meshgrid(np.arange(nx), np.arange(nz), indexing='ij')
    _red = ((ii + jj) % 2 == 0)
    _inside = ~solid
    _inside[:, 0] = False
    _inside[:, -1] = False

    def laplace(V, sweeps=340, omega=1.88, **bc):
        for _ in range(sweeps):
            V = apply_dirichlet(V, **bc)
            for color in (_red, ~_red):
                m = _inside & color
                avg = np.zeros_like(V)
                xm = np.roll(V, 1, axis=0); xp = np.roll(V, -1, axis=0)   # x periodic (pitch)
                avg[:, 1:-1] = 0.25 * (xm[:, 1:-1] + xp[:, 1:-1] + V[:, 2:] + V[:, :-2])
                V[m] = (1 - omega) * V[m] + omega * avg[m]
        return apply_dirichlet(V, **bc)

    # V_B: boundary=1, all surfaces 0 (solved once; the sheath-barrier shape for electrons)
    VB = np.zeros((nx, nz))
    if rf_bursts:
        z0 = np.zeros(W); zs_ = np.zeros(n_side)
        VB = laplace(VB, sweeps=300, boundary=1.0, vf=z0, vl=zs_, vr=zs_, vp=0.0, vtl=0.0, vtr=0.0)
    ExB = -(np.gradient(VB, axis=0)); EzB = -(np.gradient(VB, axis=1))

    # electron burst phase CDF: p(phi) ~ exp(-V_s(phi)/Te)
    _phi = np.linspace(0.0, 2.0 * np.pi, 720)
    _pw = np.exp(-(V_dc + V_rf * np.sin(_phi)) / Te)
    _pcdf = np.cumsum(_pw); _pcdf /= _pcdf[-1]

    trace_stats = []

    def trace(kind, n, Ex, Ez):
        """Ballistic trace. Ions: field E_A. Electrons: E_A + V_s(phi_e)*E_B (rf_bursts)."""
        if kind == 'ion':
            phi_p = rng.uniform(0.0, 2.0 * np.pi, n)
            E0 = np.maximum(V_dc + V_rf * np.sin(phi_p), 0.5)
            # HG energy-angle anticorrelation: transverse T_i is fixed by the presheath while the
            # sheath sets the vertical energy -> theta(E) ~ 1/sqrt(E) ("largest-angle ions have
            # least energy"). Normalized so the flux-mean HWHM = iadf_hwhm_deg at <E> = V_dc.
            sig = np.deg2rad(iadf_hwhm_deg) / 1.1774 * np.sqrt(V_dc / E0)
            th = rng.normal(0.0, sig, n)
            q = +1.0
            sB = np.zeros(n)                           # ions: no V_B (arrival KE includes the fall)
        else:
            Ee = rng.gamma(2.0, Te, n)
            # 3-D cos^p FLUX distribution reduced to the 2-D plane: sample the 3-D polar angle
            # (pdf ~ cos^p * sin, flux-weighted x cos -> cos^(p+1) sin), then project the lateral
            # component through a uniform azimuth: tan(th_2d) = sin(th3)cos(az)/cos(th3).
            # (Sampling cos^p directly in-plane over-widens the distribution -- a 2-D/3-D
            # reduction error, not physics.)
            u = rng.uniform(0.0, 1.0, n)
            ct3 = (1.0 - u) ** (1.0 / (cos_power + 2.0))
            st3 = np.sqrt(np.maximum(1.0 - ct3 * ct3, 0.0))
            az = rng.uniform(0.0, 2.0 * np.pi, n)
            th = np.arctan2(st3 * np.cos(az), ct3)
            E0 = Ee
            q = -1.0
            if rf_bursts:
                phi_e = np.interp(rng.uniform(0, 1, n), _pcdf, _phi)
                # residual unperturbed barrier below the sheath LOWER boundary (Child law
                # V(z)~z^(4/3)): fraction (boundary/sheath)^(4/3) of the instantaneous V_s --
                # ~1 V here, NOT the full sheath drop (electrons at the boundary already
                # climbed the rest; full V_s double-counts, over-suppressing the floor supply).
                frac = (boundary_um / sheath_um) ** (4.0 / 3.0)
                sB = frac * (V_dc + V_rf * np.sin(phi_e))
            else:
                sB = np.zeros(n)
        vx = np.sqrt(E0) * np.sin(th)
        vz = np.sqrt(E0) * np.abs(np.cos(th))
        x = rng.uniform(0, nx - 1.0, n)
        z = np.ones(n) * 1.0
        alive = np.ones(n, bool)
        hits_floor = np.zeros(W)
        hits_left = np.zeros(n_side); hits_right = np.zeros(n_side)
        hit_top_l = 0.0; hit_top_r = 0.0
        foot_n = 0.0; foot_E = 0.0                     # ion impacts on the POLY sidewall band
        max_steps = int(float(trace_step_cap_factor) * nz)
        if trace_integrator == "fixed":
            max_steps = int(14 * nz)
        steps_used = 0

        if trace_integrator == "adaptive_numba" and _trace_particles_adaptive is not None:
            ht, hi, fep, survivor, steps = _trace_particles_adaptive(
                Ex, Ez, ExB, EzB, x, vx, vz, sB, q, nx, nz, W, pad, mouth, n_side, z_poly0, max_steps
            )
            floor_sel = ht == 1
            if floor_sel.any():
                hits_floor += np.bincount(hi[floor_sel], minlength=W)[:W]
            left_sel = ht == 2
            if left_sel.any():
                hits_left += np.bincount(hi[left_sel], minlength=n_side)[:n_side]
            right_sel = ht == 3
            if right_sel.any():
                hits_right += np.bincount(hi[right_sel], minlength=n_side)[:n_side]
            hit_top_l = float((ht == 4).sum())
            hit_top_r = float((ht == 5).sum())
            foot_hit = fep > 0.0
            foot_n = float(foot_hit.sum())
            foot_E = float(fep[foot_hit].sum())
            survivors = int(survivor.sum())
            steps_used = int(steps.max()) if steps.size else 0
            trace_stats.append(dict(kind=kind, n=int(n), survivors=survivors,
                                    survivor_frac=float(survivors / max(n, 1)),
                                    steps=int(steps_used), cap=int(max_steps),
                                    integrator=trace_integrator))
            return hits_floor, hits_left, hits_right, hit_top_l, hit_top_r, foot_n, foot_E

        for _ in range(max_steps):
            steps_used += 1
            if not alive.any():
                break
            ix = np.clip(x[alive].astype(int), 0, nx - 2)
            iz = np.clip(z[alive].astype(int), 0, nz - 2)
            fx = Ex[ix, iz] + sB[alive] * ExB[ix, iz]
            fz = Ez[ix, iz] + sB[alive] * EzB[ix, iz]
            ax = q * fx * 0.5
            az = q * fz * 0.5
            vmax = np.maximum(np.abs(vx[alive]), np.abs(vz[alive]))
            if trace_integrator in ("adaptive", "adaptive_numba"):
                field = np.hypot(fx, fz)
                dt_v = 0.45 / np.maximum(vmax, 0.8)
                dt_e = 0.3 / np.sqrt(np.maximum(field, 1.0e-9))
                dt = np.minimum(dt_v, dt_e)
                vx_half = vx[alive] + 0.5 * ax * dt
                vz_half = vz[alive] + 0.5 * az * dt
                xa = x[alive] + vx_half * dt
                za = z[alive] + vz_half * dt
                ix2 = np.clip(xa.astype(int), 0, nx - 2)
                iz2 = np.clip(za.astype(int), 0, nz - 2)
                fx2 = Ex[ix2, iz2] + sB[alive] * ExB[ix2, iz2]
                fz2 = Ez[ix2, iz2] + sB[alive] * EzB[ix2, iz2]
                ax2 = q * fx2 * 0.5
                az2 = q * fz2 * 0.5
                vx[alive] = vx_half + 0.5 * ax2 * dt
                vz[alive] = vz_half + 0.5 * az2 * dt
                x[alive] = xa
                z[alive] = za
            elif trace_integrator == "fixed":
                dt = 0.45 / np.maximum(vmax, 0.8)
                vx[alive] += ax * dt
                vz[alive] += az * dt
                x[alive] += vx[alive] * dt
                z[alive] += vz[alive] * dt
            else:
                raise ValueError(f"unknown trace_integrator: {trace_integrator}")
            ia = np.flatnonzero(alive)
            xi, zi = x[ia], z[ia]
            x[ia] = np.mod(xi, float(nx))              # periodic pitch
            xi = x[ia]
            out = (zi < 0.5)
            alive[ia[out]] = False
            fh = (zi >= nz - 1.5)
            for j in ia[fh]:
                c = int(np.clip(x[j] - pad, 0, W - 1))
                if 0 <= x[j] - pad < W:
                    hits_floor[c] += 1
                alive[j] = False
            th_hit = (zi >= mouth) & (zi < nz - 1.5) & ((xi < pad) | (xi >= pad + W))
            for j in ia[th_hit]:
                zj = int(np.clip(z[j] - mouth, 0, n_side - 1))
                if z[j] < mouth + 1.5:
                    if x[j] < pad:
                        hit_top_l += 1
                    else:
                        hit_top_r += 1
                else:
                    if x[j] < pad + W / 2:
                        hits_left[zj] += 1
                    else:
                        hits_right[zj] += 1
                    if kind == 'ion' and z[j] >= z_poly0:      # ion striking the poly line
                        foot_n += 1.0
                        foot_E += vx[j] ** 2 + vz[j] ** 2
                alive[j] = False
        survivors = int(alive.sum())
        trace_stats.append(dict(kind=kind, n=int(n), survivors=survivors,
                                survivor_frac=float(survivors / max(n, 1)),
                                steps=int(steps_used), cap=int(max_steps),
                                integrator=trace_integrator))
        return hits_floor, hits_left, hits_right, hit_top_l, hit_top_r, foot_n, foot_E

    hist = []
    vfloor_hist, vleft_hist, vright_hist, vpoly_hist = [], [], [], []
    for it in range(n_iter):
        V = laplace(V)
        Ex = -(np.gradient(V, axis=0)); Ez = -(np.gradient(V, axis=1))
        fi, li, ri, tli, tri, _, _ = trace('ion', n_per_iter, Ex, Ez)
        fe, le, re, tle, tre, _, _ = trace('electron', n_per_iter, Ex, Ez)
        anneal = max(1.0 / (1.0 + it / 25.0), 0.25)
        scale = anneal * relax / n_per_iter * (nx)
        Vfloor += scale * (fi - fe)
        ins = ~is_poly                                 # insulating (PR) wall segments only
        Vleft[ins] += scale * (li - le)[ins]
        Vright[ins] += scale * (ri - re)[ins]
        if poly_cells > 0:                             # conductor: ONE equipotential from TOTAL net
            net_poly = (li - le)[is_poly].sum() + (ri - re)[is_poly].sum()
            Vpoly += scale * net_poly / max(2 * poly_cells, 1)
            Vpoly = float(np.clip(Vpoly, -3 * Te, V_dc + V_rf))
        Vtop_l += scale * (tli - tle) / max(pad, 1)
        Vtop_r += scale * (tri - tre) / max(pad, 1)
        vmin = -insul_vmin_Te * Te                     # interior insulator floating bound (HG: -0.5..-4 V)
        np.clip(Vfloor, 0.0, V_dc + V_rf, out=Vfloor)
        np.clip(Vleft, vmin, V_dc + V_rf, out=Vleft)
        np.clip(Vright, vmin, V_dc + V_rf, out=Vright)
        Vtop_l = float(np.clip(Vtop_l, vmin, 0.0))
        Vtop_r = float(np.clip(Vtop_r, vmin, 0.0))
        hist.append(fi.sum() / n_per_iter)
        vfloor_hist.append(Vfloor.copy())
        vleft_hist.append(Vleft.copy()); vright_hist.append(Vright.copy())
        vpoly_hist.append(Vpoly)
        if verbose and it % 20 == 0:
            tsi = next((s for s in reversed(trace_stats) if s["kind"] == "ion"), None)
            tse = next((s for s in reversed(trace_stats) if s["kind"] == "electron"), None)
            surv = ""
            if tsi is not None and tse is not None:
                surv = f" surv_i/e={tsi['survivor_frac']:.4f}/{tse['survivor_frac']:.4f}"
            print(f"  it{it}: floor_i={fi.sum()/n_per_iter:.3f} floor_e={fe.sum()/n_per_iter:.3f} "
                  f"Vc={Vfloor[W//2]:.1f} Vpoly={Vpoly:.1f}{surv}", flush=True)
    k = max(n_iter // 3, 5)
    tail = np.mean(hist[-k:])
    Vf_avg = np.mean(np.array(vfloor_hist[-k:]), axis=0)
    Vl_avg = np.mean(np.array(vleft_hist[-k:]), axis=0)
    Vr_avg = np.mean(np.array(vright_hist[-k:]), axis=0)
    Vp_avg = float(np.mean(vpoly_hist[-k:]))
    Vf_map, Vl_map, Vr_map = Vf_avg, Vl_avg, Vr_avg
    if smooth:
        from scipy.ndimage import uniform_filter1d
        Vf_map = uniform_filter1d(Vf_avg, 5, mode="nearest")
        Vl_map = uniform_filter1d(Vl_avg, 9, mode="nearest")
        Vr_map = uniform_filter1d(Vr_avg, 9, mode="nearest")
    Vfloor[:] = Vf_map; Vleft[:] = Vl_map; Vright[:] = Vr_map; Vpoly = Vp_avg
    V = laplace(V, sweeps=240)
    Ex = -(np.gradient(V, axis=0)); Ez = -(np.gradient(V, axis=1))
    ifl, ill, irr, itl, itr, fn2, fE2 = trace('ion', 4 * n_per_iter, Ex, Ez)
    efl, ell, err, etl, etr, _, _ = trace('electron', 4 * n_per_iter, Ex, Ez)
    open_frac = W / nx
    ntot = 4 * n_per_iter
    # per-species landing budget over the trench mouth (fractions of mouth-entering particles)
    ie_poly = (ill[is_poly].sum() + irr[is_poly].sum())    # ions on the poly sidewall band = foot
    ie_pr = (ill[~is_poly].sum() + irr[~is_poly].sum())     # ions on the PR (insulating) sidewalls
    ee_poly = (ell[is_poly].sum() + err[is_poly].sum())
    ee_pr = (ell[~is_poly].sum() + err[~is_poly].sum())
    diag = dict(
        ion=dict(floor=float(ifl.sum() / ntot / open_frac), poly=float(ie_poly / ntot / open_frac),
                 pr=float(ie_pr / ntot / open_frac), top=float((itl + itr) / ntot)),
        electron=dict(floor=float(efl.sum() / ntot / open_frac), poly=float(ee_poly / ntot / open_frac),
                      pr=float(ee_pr / ntot / open_frac), top=float((etl + etr) / ntot)),
        trace=dict(last_ion=next((s for s in reversed(trace_stats) if s["kind"] == "ion"), None),
                   last_electron=next((s for s in reversed(trace_stats) if s["kind"] == "electron"), None),
                   cap_factor=float(trace_step_cap_factor),
                   integrator=trace_integrator))
    _prwall = ~is_poly
    return dict(floor_flux=float(tail / open_frac), V_floor_center=float(Vf_avg[W // 2]),
                V_foot_peak=float(Vf_avg.max()), V_poly=Vp_avg, Vfloor=Vf_avg, V=V,
                foot_ion_flux=float(fn2 / (4 * n_per_iter) / open_frac),
                foot_ion_Emean=float(fE2 / max(fn2, 1.0)), diag=diag,
                geom=dict(pad=int(pad), W=int(W), mouth=int(mouth), D=int(D), nx=int(nx), nz=int(nz),
                          poly_cells=int(poly_cells)),
                Vprwall_mean=float(0.5 * (Vl_avg[_prwall].mean() + Vr_avg[_prwall].mean())) if _prwall.any() else 0.0,
                Vprwall_min=float(min(Vl_avg[_prwall].min(), Vr_avg[_prwall].min())) if _prwall.any() else 0.0)


# PRODUCTION TABLES = Hwang-Giapis PUBLISHED values (JAP 82,566): the 2-D solver is the
# validation instrument (it reproduces these mechanisms with nothing tuned -- flux RMSE 0.060,
# V_poly 6->39 V within 11%, Matsui asymptote pass; deep-AR foot energies under-predicted, the
# documented residual); the PUBLISHED data is what production interpolates. Solver-measured
# values for reference (2026-07-02 faithful config): flux 0.681/0.629/0.538/0.465/0.359/0.284/
# 0.209/0.174, V_floor 11.1/12.5/18.4/23.1/35.1/42.7/51.0/53.6, V_poly ~6->43.
_GATE_AR = np.array([0.0, 1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
_GATE_FLUX = np.array([1.0, 0.59, 0.55, 0.47, 0.40, 0.34, 0.30, 0.26, 0.22])
_GATE_VFLOOR = np.array([0.0, 8.0, 10.0, 15.0, 19.0, 25.0, 28.0, 31.0, 33.0])
_FOOT_E = np.array([0.0, 15.0, 16.5, 17.5, 20.0, 23.0, 25.0, 26.5, 27.5])


def charging_floor_profile(AR):
    """Production hook (NOT yet wired into the flux pipeline): normalized floor ion flux, floor
    potential, and deflected-ion foot energy vs aspect ratio, from the gate-validated 2-D solver
    at the HG reference conditions. INSULATING floors only (SiO2/SOI overetch, dielectric etch);
    a conductive grounded Si floor drains and should NOT be throttled (the de Boer deep-Si case)."""
    AR = np.asarray(AR, float)
    return (np.interp(AR, _GATE_AR, _GATE_FLUX),
            np.interp(AR, _GATE_AR, _GATE_VFLOOR),
            np.interp(AR, _GATE_AR, _FOOT_E))
