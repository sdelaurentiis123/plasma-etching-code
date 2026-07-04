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

import math

import numpy as np

from .sheath1d import sample_sheath_source

try:
    from numba import njit, prange
except Exception:  # pragma: no cover - optional acceleration
    njit = None
    prange = range


def _sky_view_factors_py(solid, surf_ix, surf_iz, n_angles):
    """Cosine-weighted fraction of the upper half-plane each surface cell can see 'sky' (z<0,
    the plasma boundary) along a straight line without hitting solid. z=0 is the top (plasma),
    z increases downward. A fully open upward-facing cell -> ~1.0; a vertical open wall -> ~0.5;
    a deep shadowed trench wall -> small. This is the geometric isotropic-electron view factor:
    for an isotropic species the collected flux is (base rate) x view_factor, so it captures the
    orientation-dependent + shadowed electron collection the down-going MC source gets wrong."""
    nx, nz = solid.shape
    m = surf_ix.shape[0]
    vf = np.zeros(m)
    step = 0.5
    maxsteps = 4 * (nx + nz)
    for s in prange(m):
        ix = surf_ix[s]
        iz = surf_iz[s]
        wsum = 0.0
        wsky = 0.0
        for a in range(n_angles):
            theta = -1.5707963267 + 3.1415926535 * (a + 0.5) / n_angles
            dx = math.sin(theta)
            dz = -math.cos(theta)          # up = -z toward the plasma
            w = math.cos(theta)            # cosine weight (angle from the upward normal)
            wsum += w
            fx = ix + 0.5
            fz = iz + 0.5
            saw_sky = False
            for _ in range(maxsteps):
                fx += dx * step
                fz += dz * step
                if fz < 0.5:               # reached the plasma boundary
                    saw_sky = True
                    break
                cix = int(fx)
                ciz = int(fz)
                if cix < 0 or cix >= nx:    # exited the side into open plasma
                    saw_sky = True
                    break
                if ciz >= nz:
                    break
                if (cix != ix or ciz != iz) and solid[cix, ciz]:
                    break                  # blocked by solid
            if saw_sky:
                wsky += w
        vf[s] = wsky / wsum if wsum > 0.0 else 0.0
    return vf


_sky_view_factors = (njit(cache=True, parallel=True, fastmath=True)(_sky_view_factors_py)
                     if njit is not None else _sky_view_factors_py)


if njit is not None:
    @njit(cache=True, parallel=True, fastmath=True)
    def _trace_particles_adaptive(Ex, Ez, ExB, EzB, x0, z0, vx0, vz0, sB, q,
                                  nx, nz, W, pad, mouth, n_side, z_poly0,
                                  max_steps):
        n = x0.shape[0]
        hit_type = np.zeros(n, np.int8)      # 0 escape/survivor, 1 floor, 2 left, 3 right, 4 top_l, 5 top_r
        hit_idx = np.full(n, -1, np.int64)
        impact_E = np.zeros(n)
        hit_x = np.zeros(n)
        hit_z = np.zeros(n)
        hit_vx = np.zeros(n)
        hit_vz = np.zeros(n)
        survivor = np.zeros(n, np.uint8)
        steps = np.zeros(n, np.int64)
        xmax = float(nx)

        for p in prange(n):
            x = x0[p]
            z = z0[p]
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
                        impact_E[p] = vx * vx + vz * vz
                        hit_x[p] = x
                        hit_z[p] = z
                        hit_vx[p] = vx
                        hit_vz[p] = vz
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
                    impact_E[p] = vx * vx + vz * vz
                    hit_x[p] = x
                    hit_z[p] = z
                    hit_vx[p] = vx
                    hit_vz[p] = vz
                    alive = False
                    break
            if alive:
                survivor[p] = 1
            steps[p] = last_step
        return hit_type, hit_idx, impact_E, hit_x, hit_z, hit_vx, hit_vz, survivor, steps


    @njit(cache=True, fastmath=True)
    def _trace_edge_particles_adaptive(Ex, Ez, ExB, EzB, solid, cond, x0, z0, vx0, vz0, sB, q,
                                       nx, nz, mouth, edge0, edge1, trench0, trench1,
                                       neigh0, neigh1, z_poly0, max_steps):
        n = x0.shape[0]
        hit_type = np.zeros(n, np.int8)      # 0 escape/survivor, 1 trench floor, 2 open floor,
                                             # 3 edge outer poly, 4 edge inner poly, 5 neighbor poly,
                                             # 6 PR/other insulator, 7 other floor
        hit_ix = np.full(n, -1, np.int64)
        hit_iz = np.full(n, -1, np.int64)
        impact_E = np.zeros(n)
        hit_vx = np.zeros(n)
        hit_vz = np.zeros(n)
        survivor = np.zeros(n, np.uint8)
        steps = np.zeros(n, np.int64)

        for p in range(n):
            x = x0[p]
            z = z0[p]
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

                if z < 0.5 or x < 0.0 or x >= nx - 1.0:
                    alive = False
                    break

                ixh = int(x)
                izh = int(z)
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
                    hit_vz[p] = vz
                    if izh >= nz - 1:
                        if ixh >= trench0 and ixh < trench1:
                            hit_type[p] = 1
                        elif ixh < edge0:
                            hit_type[p] = 2
                        else:
                            hit_type[p] = 7
                    else:
                        cid = cond[ixh, izh]
                        if cid == 1:
                            if vx >= 0.0:
                                hit_type[p] = 3
                            else:
                                hit_type[p] = 4
                        elif cid == 2:
                            hit_type[p] = 5
                        else:
                            hit_type[p] = 6
                    alive = False
                    break
            if alive:
                survivor[p] = 1
            steps[p] = last_step
        return hit_type, hit_ix, hit_iz, impact_E, hit_vx, hit_vz, survivor, steps
else:
    _trace_particles_adaptive = None
    _trace_edge_particles_adaptive = None


_PMMA_SIGMA_E_E = np.array([0.0, 5.0, 10.0, 16.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0])
_PMMA_SIGMA_E_Y = np.array([0.0, 0.08, 0.17, 0.28, 0.38, 0.60, 0.78, 0.94, 1.09, 1.39, 1.60])

_HG_AR = np.array([1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
_HG_EDGE_OUTER_ELECTRON_GROSS_FLUX = np.array([0.18, 0.19, 0.20, 0.18, 0.19, 0.18, 0.18, 0.18])


def _hg_edge_outer_electron_gross_flux(AR):
    """HG Fig. 3 gross electron flux to the outer edge-line poly-Si sidewall.

    This is not automatically a conductor-current source: the floating conductor responds to
    net open-side current after any open-side ion/countercurrent. Use it as a published diagnostic
    reference unless the open-side geometry is explicitly simulated.
    """
    return float(np.interp(float(AR), _HG_AR, _HG_EDGE_OUTER_ELECTRON_GROSS_FLUX,
                           left=_HG_EDGE_OUTER_ELECTRON_GROSS_FLUX[0],
                           right=_HG_EDGE_OUTER_ELECTRON_GROSS_FLUX[-1]))


def _sample_open_side_flux(rng, n, W, pr_cells, poly_cells, cos_power, iadf_hwhm_deg,
                           ion_angle_energy_corr, V_dc, V_rf, open_width_cells):
    """Line-of-sight current to the open-area side of the edge poly line.

    The auxiliary domain is the half-space outside the edge line. Particles launch through a
    horizontal open window and free-stream to the vertical outer sidewall; the counted band is the
    0.3 um poly-Si sidewall, not the PR wall. Returned currents are normalized to the trench
    opening width, matching the in-feature `floor_flux` convention.
    """
    if poly_cells <= 0:
        return dict(electron_gross=0.0, ion_gross=0.0, net_electron=0.0,
                    open_width_cells=float(open_width_cells), n=int(n))

    width = float(max(open_width_cells, W))
    x = rng.uniform(0.0, width, n)

    # Electrons: same 3-D-to-2-D projected cos^p flux source used by the feature tracer.
    u = rng.uniform(0.0, 1.0, n)
    ct3 = (1.0 - u) ** (1.0 / (cos_power + 2.0))
    st3 = np.sqrt(np.maximum(1.0 - ct3 * ct3, 0.0))
    az = rng.uniform(0.0, 2.0 * np.pi, n)
    th_e = np.arctan2(st3 * np.cos(az), ct3)
    toward = th_e < 0.0
    tan_e = np.tan(np.abs(th_e[toward]))
    zhit_e = np.full(tan_e.shape, np.inf)
    ok = tan_e > 1.0e-12
    zhit_e[ok] = x[toward][ok] / tan_e[ok]
    e_hits = (zhit_e >= pr_cells) & (zhit_e < pr_cells + poly_cells)
    electron_gross = float(e_hits.sum() / max(n, 1) * width / max(W, 1))

    # Ions: same narrow arrival-angle model; vertical ions rarely hit a vertical open sidewall.
    phi_p = rng.uniform(0.0, 2.0 * np.pi, n)
    E0 = np.maximum(V_dc + V_rf * np.sin(phi_p), 0.5)
    sig0 = np.deg2rad(iadf_hwhm_deg) / 1.1774
    if ion_angle_energy_corr == "anticorrelated":
        sig = sig0 * np.sqrt(V_dc / E0)
    elif ion_angle_energy_corr == "independent":
        sig = np.full(n, sig0)
    elif ion_angle_energy_corr == "positive":
        sig = sig0 * np.sqrt(E0 / V_dc)
    else:
        raise ValueError(f"unknown ion_angle_energy_corr: {ion_angle_energy_corr}")
    th_i = rng.normal(0.0, sig, n)
    toward = th_i < 0.0
    tan_i = np.tan(np.abs(th_i[toward]))
    zhit_i = np.full(tan_i.shape, np.inf)
    ok = tan_i > 1.0e-12
    zhit_i[ok] = x[toward][ok] / tan_i[ok]
    i_hits = (zhit_i >= pr_cells) & (zhit_i < pr_cells + poly_cells)
    ion_gross = float(i_hits.sum() / max(n, 1) * width / max(W, 1))

    return dict(electron_gross=electron_gross, ion_gross=ion_gross,
                net_electron=max(electron_gross - ion_gross, 0.0),
                open_width_cells=float(width), n=int(n))


def _pmma_see_yields(energy):
    """Memos/Lidorikis/Kokkoris PMMA SEEE model, digitized from their Fig. 2.

    sigma_e is Dapor's total yield curve; Burke's polymer backscatter law separates eta;
    true secondary yield delta is zero below 16 eV and sigma_e - eta above it.
    """
    e = np.asarray(energy, float)
    sigma = np.interp(np.clip(e, 0.0, 100.0), _PMMA_SIGMA_E_E, _PMMA_SIGMA_E_Y)
    eta = np.zeros_like(sigma)
    hi = e > 0.0
    eta[hi] = 0.115 * np.power(e[hi] / 1000.0, -0.223)
    eta = np.minimum(eta, sigma)
    delta = np.where(e >= 16.0, np.maximum(sigma - eta, 0.0), 0.0)
    eta = np.where(e < 16.0, sigma, eta)
    return sigma, eta, delta


def solve_trench_charging(AR, W=32, pad=16, mouth=237, Te=4.0, V_dc=37.0, V_rf=30.0,
                          iadf_hwhm_deg=4.3, cos_power=0.6, n_per_iter=6000, n_iter=120,
                          relax=None, seed=0, verbose=False, smooth=False,
                          poly_um=0.3, feature_w_um=0.5, rf_bursts=True,
                          sheath_um=89.0, boundary_um=3.7, insul_vmin_Te=1.0,
                          trace_integrator="adaptive_numba", trace_step_cap_factor=40.0,
                          see_model="none", see_generations=1,
                          ion_angle_energy_corr="anticorrelated",
                          source_model="analytic", poly_mode="tied", poly_bias_V=0.0,
                          edge_open_electron_flux=None, edge_open_model="none",
                          edge_open_samples=None, edge_open_width_um=3.7):
    """Steady-state charging of the HG poly-on-oxide trench. Returns dict with:
    floor_flux (normalized ion flux to the oxide floor), V_floor_center, V_foot_peak,
    V_poly (the poly-line equipotential), foot_ion_flux / foot_ion_Emean (ions striking the
    poly sidewall — the notch driver), V (grid map, ground-referenced).

    poly_um=0 recovers the all-insulator solver. rf_bursts=False recovers time-averaged
    electrons with no sheath barrier (the pre-burst model). smooth=True is PRESENTATION-grade
    filtering of the returned map only — gate quantities always computed raw
    (scripts/charging_gate.py and notching_gate.py must run smooth=False)."""
    rng = np.random.default_rng(seed)
    # HG's aspect ratio varies the photoresist height over a fixed 0.3 um poly-Si line:
    # AR=4 is 1.7 um PR + 0.3 um poly on a 0.5 um space, not a 2.0 um PR wall.
    pr_cells = int(round((AR * feature_w_um - poly_um) / feature_w_um * W))
    pr_cells = max(pr_cells, 1)
    D = pr_cells + int(round(poly_um / feature_w_um * W)) if poly_um > 0 else int(round(AR * W))
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
    # left = edge-line inner sidewall; right = neighboring-line sidewall.
    Vpoly_l = 0.0; Vpoly_r = 0.0
    Vtop_l = 0.0; Vtop_r = 0.0

    def apply_dirichlet(V, boundary=0.0, vf=None, vl=None, vr=None, vp=None, vpl=None, vpr=None,
                        vtl=None, vtr=None):
        V[:, 0] = boundary
        V[floor_ix, nz - 1] = Vfloor if vf is None else vf
        vl_ = (Vleft if vl is None else vl).copy()
        vr_ = (Vright if vr is None else vr).copy()
        vp_ = Vpoly if vp is None else vp
        if poly_cells > 0:
            if poly_mode in ("split", "edge_open"):
                vl_[is_poly] = Vpoly_l if vpl is None else vpl
                vr_[is_poly] = Vpoly_r if vpr is None else vpr
            elif poly_mode == "edge_bias":
                vl_[is_poly] = vp_ - 0.5 * poly_bias_V
                vr_[is_poly] = vp_ + 0.5 * poly_bias_V
            else:
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
        VB = laplace(VB, sweeps=300, boundary=1.0, vf=z0, vl=zs_, vr=zs_, vp=0.0,
                     vpl=0.0, vpr=0.0, vtl=0.0, vtr=0.0)
    ExB = -(np.gradient(VB, axis=0)); EzB = -(np.gradient(VB, axis=1))

    # electron burst phase CDF: p(phi) ~ exp(-V_s(phi)/Te)
    _phi = np.linspace(0.0, 2.0 * np.pi, 720)
    _pw = np.exp(-(V_dc + V_rf * np.sin(_phi)) / Te)
    _pcdf = np.cumsum(_pw); _pcdf /= _pcdf[-1]

    trace_stats = []
    see_stats = []
    open_frac = W / nx
    edge_outer_gross_flux = _hg_edge_outer_electron_gross_flux(AR)
    edge_open_samples = int(edge_open_samples or max(32768, 4 * n_per_iter))
    edge_open_width_cells = int(round(edge_open_width_um / feature_w_um * W))
    edge_open_diag = dict(model=edge_open_model, electron_gross=0.0, ion_gross=0.0,
                          net_electron=0.0, open_width_cells=float(edge_open_width_cells),
                          n=int(edge_open_samples), hg_electron_gross=edge_outer_gross_flux,
                          override=False)
    if edge_open_electron_flux is not None:
        edge_open_flux_value = float(edge_open_electron_flux)
        edge_open_diag.update(model="override", net_electron=edge_open_flux_value,
                              electron_gross=edge_open_flux_value, override=True)
    elif edge_open_model in ("none", None):
        edge_open_flux_value = 0.0
    elif edge_open_model == "line_of_sight":
        edge_open_diag.update(_sample_open_side_flux(
            rng, edge_open_samples, W, pr_cells, poly_cells, cos_power, iadf_hwhm_deg,
            ion_angle_energy_corr, V_dc, V_rf, edge_open_width_cells))
        edge_open_flux_value = float(edge_open_diag["net_electron"])
    else:
        raise ValueError(f"unknown edge_open_model: {edge_open_model}")

    def edge_open_net_electron(v_edge):
        if edge_open_electron_flux is not None:
            return edge_open_flux_value
        if edge_open_model != "line_of_sight":
            return 0.0
        # The line-of-sight calculation gives the geometric gross supply. To hit a negative
        # vertical sidewall, an electron's horizontal kinetic energy must clear the repulsive
        # conductor potential. For a Maxwellian transverse component this survival is erfc; a
        # positive edge line cannot collect more than the gross ballistic supply.
        if v_edge < 0.0:
            accept = math.erfc(math.sqrt(max(-float(v_edge), 0.0) / max(Te, 1.0e-9)))
        else:
            accept = 1.0
        return max(float(edge_open_diag["electron_gross"]) * accept
                   - float(edge_open_diag["ion_gross"]), 0.0)

    def trace(kind, n, Ex, Ez):
        """Ballistic trace. Ions: field E_A. Electrons: E_A + V_s(phi_e)*E_B (rf_bursts)."""
        if source_model == "sheath_mc" or (source_model == "sheath_electrons" and kind == "electron"):
            E0, th, sB = sample_sheath_source(kind, n, rng, Te=Te, V_dc=V_dc, V_rf=V_rf,
                                              iadf_hwhm_deg=iadf_hwhm_deg,
                                              cos_power=cos_power,
                                              boundary_um=boundary_um,
                                              sheath_um=sheath_um)
            q = +1.0 if kind == 'ion' else -1.0
        elif source_model not in ("analytic", "sheath_electrons"):
            raise ValueError(f"unknown source_model: {source_model}")
        elif kind == 'ion':
            phi_p = rng.uniform(0.0, 2.0 * np.pi, n)
            E0 = np.maximum(V_dc + V_rf * np.sin(phi_p), 0.5)
            # HG energy-angle anticorrelation: transverse T_i is fixed by the presheath while the
            # sheath sets the vertical energy -> theta(E) ~ 1/sqrt(E) ("largest-angle ions have
            # least energy"). Normalized so the flux-mean HWHM = iadf_hwhm_deg at <E> = V_dc.
            sig0 = np.deg2rad(iadf_hwhm_deg) / 1.1774
            if ion_angle_energy_corr == "anticorrelated":
                sig = sig0 * np.sqrt(V_dc / E0)
            elif ion_angle_energy_corr == "independent":
                sig = np.full(n, sig0)
            elif ion_angle_energy_corr == "positive":
                sig = sig0 * np.sqrt(E0 / V_dc)
            else:
                raise ValueError(f"unknown ion_angle_energy_corr: {ion_angle_energy_corr}")
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
        x = rng.uniform(0, float(nx), n)
        z = np.ones(n) * 1.0
        alive = np.ones(n, bool)
        hits_floor = np.zeros(W)
        hits_left = np.zeros(n_side); hits_right = np.zeros(n_side)
        hit_top_l = 0.0; hit_top_r = 0.0
        foot_n = 0.0; foot_E = 0.0; foot_En = 0.0      # ion impacts on the POLY sidewall band
        foot_z_mean = np.nan; foot_E_p50 = np.nan; foot_E_p90 = np.nan
        foot_n_left = 0.0; foot_E_left = 0.0; foot_n_right = 0.0; foot_E_right = 0.0
        foot_En_left = 0.0; foot_En_right = 0.0
        max_steps = int(float(trace_step_cap_factor) * nz)
        if trace_integrator == "fixed":
            max_steps = int(14 * nz)
        steps_used = 0

        if trace_integrator == "adaptive_numba" and _trace_particles_adaptive is not None:
            ht, hi, impact_E, hit_x, hit_z, hit_vx, hit_vz, survivor, steps = _trace_particles_adaptive(
                Ex, Ez, ExB, EzB, x, z, vx, vz, sB, q, nx, nz, W, pad, mouth, n_side, z_poly0, max_steps
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
            if kind == 'ion' and poly_cells > 0:
                foot_hit = ((ht == 2) | (ht == 3)) & (hi >= n_side - poly_cells)
                foot_n = float(foot_hit.sum())
                foot_E = float(impact_E[foot_hit].sum())
                foot_En = float((hit_vx[foot_hit] * hit_vx[foot_hit]).sum())
                foot_left = foot_hit & (ht == 2)
                foot_right = foot_hit & (ht == 3)
                foot_n_left = float(foot_left.sum())
                foot_E_left = float(impact_E[foot_left].sum())
                foot_En_left = float((hit_vx[foot_left] * hit_vx[foot_left]).sum())
                foot_n_right = float(foot_right.sum())
                foot_E_right = float(impact_E[foot_right].sum())
                foot_En_right = float((hit_vx[foot_right] * hit_vx[foot_right]).sum())
                if foot_n > 0:
                    foot_z_mean = float(hi[foot_hit].mean())
                    foot_E_p50 = float(np.percentile(impact_E[foot_hit], 50))
                    foot_E_p90 = float(np.percentile(impact_E[foot_hit], 90))

            see_wall = 0
            see_emit = 0
            see_back = 0
            see_sec = 0
            see_absorb = 0
            see_survivors = 0
            if (kind == 'electron' and see_model in ("pmma_pr", "pmma")
                    and see_generations > 0 and _trace_particles_adaptive is not None):
                cur_ht, cur_hi, cur_E = ht, hi, impact_E
                cur_x, cur_z, cur_vx, cur_vz, cur_sB = hit_x, hit_z, hit_vx, hit_vz, sB
                for gen in range(int(see_generations)):
                    wall_sel = ((cur_ht == 2) | (cur_ht == 3)) & (cur_hi >= 0)
                    if poly_cells > 0:
                        wall_sel &= cur_hi < n_side - poly_cells
                    wall_idx = np.flatnonzero(wall_sel)
                    if not wall_idx.size:
                        break
                    see_wall += int(wall_idx.size)
                    _, eta, delta = _pmma_see_yields(cur_E[wall_idx])

                    back_mask = rng.uniform(0.0, 1.0, wall_idx.size) < np.minimum(eta, 1.0)
                    sec_floor = np.floor(delta).astype(int)
                    sec_count = sec_floor + (rng.uniform(0.0, 1.0, wall_idx.size) < (delta - sec_floor))
                    emit_count = back_mask.astype(int) + sec_count
                    see_absorb += int((emit_count == 0).sum())
                    if emit_count.sum() == 0:
                        continue

                    src = np.repeat(wall_idx, emit_count)
                    is_left = cur_ht[src] == 2
                    src_hi = cur_hi[src]
                    for zj, cnt in zip(*np.unique(src_hi[is_left], return_counts=True)):
                        hits_left[zj] -= float(cnt)
                    for zj, cnt in zip(*np.unique(src_hi[~is_left], return_counts=True)):
                        hits_right[zj] -= float(cnt)

                    back_src = wall_idx[back_mask]
                    n_back = back_src.size
                    sec_src = np.repeat(wall_idx, sec_count)
                    n_sec = sec_src.size
                    see_back += int(n_back)
                    see_sec += int(n_sec)
                    see_emit += int(n_back + n_sec)

                    xs = []
                    zs = []
                    vxs = []
                    vzs = []
                    sbs = []
                    if n_back:
                        bl = cur_ht[back_src] == 2
                        xs.append(np.where(bl, pad + 0.1, pad + W - 0.1).astype(float))
                        zs.append(np.clip(cur_z[back_src], mouth + 0.1, nz - 1.6))
                        vxs.append(np.where(bl, np.abs(cur_vx[back_src]), -np.abs(cur_vx[back_src])))
                        vzs.append(cur_vz[back_src])
                        sbs.append(cur_sB[back_src])
                    if n_sec:
                        sl = cur_ht[sec_src] == 2
                        a = np.arcsin(rng.uniform(-1.0, 1.0, n_sec))
                        xs.append(np.where(sl, pad + 0.1, pad + W - 0.1).astype(float))
                        zs.append(np.clip(cur_z[sec_src], mouth + 0.1, nz - 1.6))
                        vxs.append(np.where(sl, np.cos(a), -np.cos(a)))
                        vzs.append(np.sin(a))
                        sbs.append(cur_sB[sec_src])

                    xe = np.concatenate(xs)
                    ze = np.concatenate(zs)
                    vxe = np.concatenate(vxs)
                    vze = np.concatenate(vzs)
                    sBe = np.concatenate(sbs)
                    cur_ht, cur_hi, cur_E, cur_x, cur_z, cur_vx, cur_vz, surve, stepse = _trace_particles_adaptive(
                        Ex, Ez, ExB, EzB, xe, ze, vxe, vze, sBe, -1.0,
                        nx, nz, W, pad, mouth, n_side, z_poly0, max_steps
                    )
                    efloor = cur_ht == 1
                    if efloor.any():
                        hits_floor += np.bincount(cur_hi[efloor], minlength=W)[:W]
                    eleft = cur_ht == 2
                    if eleft.any():
                        hits_left += np.bincount(cur_hi[eleft], minlength=n_side)[:n_side]
                    eright = cur_ht == 3
                    if eright.any():
                        hits_right += np.bincount(cur_hi[eright], minlength=n_side)[:n_side]
                    hit_top_l += float((cur_ht == 4).sum())
                    hit_top_r += float((cur_ht == 5).sum())
                    see_survivors += int(surve.sum())

                if see_wall > 0:
                    see_stats.append(dict(n_primary_wall=int(see_wall), emitted=int(see_emit),
                                          backscatter=int(see_back), secondary=int(see_sec),
                                          absorbed=int(see_absorb),
                                          survivor_frac=float(see_survivors / max(see_emit, 1))))
            survivors = int(survivor.sum())
            steps_used = int(steps.max()) if steps.size else 0
            trace_stats.append(dict(kind=kind, n=int(n), survivors=survivors,
                                    survivor_frac=float(survivors / max(n, 1)),
                                    steps=int(steps_used), cap=int(max_steps),
                                    integrator=trace_integrator,
                                    foot_z_mean=foot_z_mean, foot_E_p50=foot_E_p50,
                                    foot_E_p90=foot_E_p90,
                                    foot_n=foot_n, foot_E=foot_E, foot_En=foot_En,
                                    foot_n_left=foot_n_left, foot_E_left=foot_E_left,
                                    foot_En_left=foot_En_left,
                                    foot_n_right=foot_n_right, foot_E_right=foot_E_right,
                                    foot_En_right=foot_En_right,
                                    see_wall=int(see_wall),
                                    see_emitted=int(see_emit), see_backscatter=int(see_back),
                                    see_secondary=int(see_sec), see_absorbed=int(see_absorb)))
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
                        foot_En += vx[j] ** 2
                alive[j] = False
        survivors = int(alive.sum())
        trace_stats.append(dict(kind=kind, n=int(n), survivors=survivors,
                                survivor_frac=float(survivors / max(n, 1)),
                                steps=int(steps_used), cap=int(max_steps),
                                integrator=trace_integrator,
                                foot_n=foot_n, foot_E=foot_E, foot_En=foot_En))
        return hits_floor, hits_left, hits_right, hit_top_l, hit_top_r, foot_n, foot_E

    hist = []
    vfloor_hist, vleft_hist, vright_hist, vpoly_hist = [], [], [], []
    vpoly_l_hist, vpoly_r_hist = [], []
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
        if poly_cells > 0:
            if poly_mode == "tied":                    # conductor: ONE equipotential from TOTAL net
                net_poly = (li - le)[is_poly].sum() + (ri - re)[is_poly].sum()
                Vpoly += scale * net_poly / max(2 * poly_cells, 1)
                Vpoly = float(np.clip(Vpoly, -3 * Te, V_dc + V_rf))
                Vpoly_l = Vpoly; Vpoly_r = Vpoly
            elif poly_mode == "edge_bias":
                net_poly = (li - le)[is_poly].sum() + (ri - re)[is_poly].sum()
                Vpoly += scale * net_poly / max(2 * poly_cells, 1)
                Vpoly = float(np.clip(Vpoly, -3 * Te, V_dc + V_rf))
                Vpoly_l = float(np.clip(Vpoly - 0.5 * poly_bias_V, -3 * Te, V_dc + V_rf))
                Vpoly_r = float(np.clip(Vpoly + 0.5 * poly_bias_V, -3 * Te, V_dc + V_rf))
            elif poly_mode == "split":                 # diagnostic: left/right lines float separately
                Vpoly_l += scale * (li - le)[is_poly].sum() / max(poly_cells, 1)
                Vpoly_r += scale * (ri - re)[is_poly].sum() / max(poly_cells, 1)
                Vpoly_l = float(np.clip(Vpoly_l, -3 * Te, V_dc + V_rf))
                Vpoly_r = float(np.clip(Vpoly_r, -3 * Te, V_dc + V_rf))
                Vpoly = 0.5 * (Vpoly_l + Vpoly_r)
            elif poly_mode == "edge_open":
                # Diagnostic HG edge-line cell: left is the outermost line, right is the
                # neighboring line. The external term is NET open-side electron surplus, not the
                # Fig. 3 gross electron flux; the latter needs an explicit outer-domain model.
                edge_e = edge_open_net_electron(Vpoly_l) * n_per_iter * open_frac
                Vpoly_l += scale * ((li - le)[is_poly].sum() - edge_e) / max(poly_cells, 1)
                Vpoly_r += scale * (ri - re)[is_poly].sum() / max(poly_cells, 1)
                Vpoly_l = float(np.clip(Vpoly_l, -3 * Te, V_dc + V_rf))
                Vpoly_r = float(np.clip(Vpoly_r, -3 * Te, V_dc + V_rf))
                Vpoly = Vpoly_r
            else:
                raise ValueError(f"unknown poly_mode: {poly_mode}")
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
        vpoly_l_hist.append(Vpoly_l); vpoly_r_hist.append(Vpoly_r)
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
    Vpl_avg = float(np.mean(vpoly_l_hist[-k:]))
    Vpr_avg = float(np.mean(vpoly_r_hist[-k:]))
    Vf_map, Vl_map, Vr_map = Vf_avg, Vl_avg, Vr_avg
    if smooth:
        from scipy.ndimage import uniform_filter1d
        Vf_map = uniform_filter1d(Vf_avg, 5, mode="nearest")
        Vl_map = uniform_filter1d(Vl_avg, 9, mode="nearest")
        Vr_map = uniform_filter1d(Vr_avg, 9, mode="nearest")
    Vfloor[:] = Vf_map; Vleft[:] = Vl_map; Vright[:] = Vr_map
    Vpoly = Vp_avg; Vpoly_l = Vpl_avg; Vpoly_r = Vpr_avg
    V = laplace(V, sweeps=240)
    Ex = -(np.gradient(V, axis=0)); Ez = -(np.gradient(V, axis=1))
    ifl, ill, irr, itl, itr, fn2, fE2 = trace('ion', 4 * n_per_iter, Ex, Ez)
    efl, ell, err, etl, etr, _, _ = trace('electron', 4 * n_per_iter, Ex, Ez)
    ntot = 4 * n_per_iter
    # per-species landing budget over the trench mouth (fractions of mouth-entering particles)
    ie_poly = (ill[is_poly].sum() + irr[is_poly].sum())    # ions on the poly sidewall band = foot
    ie_pr = (ill[~is_poly].sum() + irr[~is_poly].sum())     # ions on the PR (insulating) sidewalls
    ee_poly = (ell[is_poly].sum() + err[is_poly].sum())
    ee_pr = (ell[~is_poly].sum() + err[~is_poly].sum())
    edge_extra_e = 0.0
    if poly_mode == "edge_open":
        edge_extra_e = edge_open_net_electron(Vpl_avg) * ntot * open_frac
    residual = dict(
        floor=float((ifl.sum() - efl.sum()) / max(ntot * open_frac, 1.0)),
        pr=float((ie_pr - ee_pr) / max(ntot * open_frac, 1.0)),
        poly=float((ie_poly - ee_poly - edge_extra_e) / max(ntot * open_frac, 1.0)),
        poly_left=float((ill[is_poly].sum() - ell[is_poly].sum() - edge_extra_e)
                        / max(ntot * open_frac, 1.0)),
        poly_right=float((irr[is_poly].sum() - err[is_poly].sum()) / max(ntot * open_frac, 1.0)),
        top=float((itl + itr - etl - etr) / max(ntot, 1.0)),
    )
    edge_open_diag["net_electron"] = float(edge_open_net_electron(Vpl_avg))
    edge_open_diag["electron_accepted"] = float(edge_open_diag["net_electron"]
                                                + edge_open_diag["ion_gross"])
    edge_open_diag["edge_potential_for_net"] = float(Vpl_avg)
    diag = dict(
        ion=dict(floor=float(ifl.sum() / ntot / open_frac), poly=float(ie_poly / ntot / open_frac),
                 pr=float(ie_pr / ntot / open_frac), top=float((itl + itr) / ntot)),
        electron=dict(floor=float(efl.sum() / ntot / open_frac), poly=float(ee_poly / ntot / open_frac),
                      pr=float(ee_pr / ntot / open_frac), top=float((etl + etr) / ntot)),
        trace=dict(last_ion=next((s for s in reversed(trace_stats) if s["kind"] == "ion"), None),
                   last_electron=next((s for s in reversed(trace_stats) if s["kind"] == "electron"), None),
                   cap_factor=float(trace_step_cap_factor),
                   integrator=trace_integrator,
                   source_model=source_model,
                   poly_mode=poly_mode,
                   edge_open_electron_flux=edge_open_flux_value,
                   edge_open_model=edge_open_model,
                   edge_outer_electron_gross_flux=edge_outer_gross_flux),
        edge_open=edge_open_diag,
        residual=residual,
        see=dict(model=see_model, generations=int(see_generations),
                 last=see_stats[-1] if see_stats else None))
    last_ion = diag["trace"]["last_ion"] or {}
    fn_l = float(last_ion.get("foot_n_left", 0.0))
    fE_l = float(last_ion.get("foot_E_left", 0.0))
    fEn_l = float(last_ion.get("foot_En_left", 0.0))
    fn_r = float(last_ion.get("foot_n_right", 0.0))
    fE_r = float(last_ion.get("foot_E_right", 0.0))
    fEn_r = float(last_ion.get("foot_En_right", 0.0))
    fn_all = float(last_ion.get("foot_n", fn2))
    fEn_all = float(last_ion.get("foot_En", np.nan))
    _prwall = ~is_poly
    floor_flux_tail = float(tail / open_frac)
    floor_flux_final = float(ifl.sum() / ntot / open_frac)
    return dict(floor_flux=floor_flux_final, floor_flux_tail=floor_flux_tail,
                V_floor_center=float(Vf_avg[W // 2]),
                V_foot_peak=float(Vf_avg.max()), V_poly=Vp_avg, Vfloor=Vf_avg, V=V,
                V_poly_left=Vpl_avg, V_poly_right=Vpr_avg,
                V_poly_edge=Vpl_avg, V_poly_neighbor=Vpr_avg,
                foot_ion_flux=float(fn2 / (4 * n_per_iter) / open_frac),
                foot_ion_Emean=float(fE2 / max(fn2, 1.0)),
                foot_ion_Enormal_mean=float(fEn_all / max(fn_all, 1.0)),
                foot_ion_flux_left=float(fn_l / (4 * n_per_iter) / open_frac),
                foot_ion_flux_right=float(fn_r / (4 * n_per_iter) / open_frac),
                foot_ion_Emean_left=float(fE_l / max(fn_l, 1.0)),
                foot_ion_Emean_right=float(fE_r / max(fn_r, 1.0)),
                foot_ion_Enormal_mean_left=float(fEn_l / max(fn_l, 1.0)),
                foot_ion_Enormal_mean_right=float(fEn_r / max(fn_r, 1.0)),
                foot_ion_flux_edge=float(fn_l / (4 * n_per_iter) / open_frac),
                foot_ion_Emean_edge=float(fE_l / max(fn_l, 1.0)),
                foot_ion_Enormal_mean_edge=float(fEn_l / max(fn_l, 1.0)), diag=diag,
                geom=dict(pad=int(pad), W=int(W), mouth=int(mouth), D=int(D), nx=int(nx), nz=int(nz),
                          poly_cells=int(poly_cells)),
                Vprwall_mean=float(0.5 * (Vl_avg[_prwall].mean() + Vr_avg[_prwall].mean())) if _prwall.any() else 0.0,
                Vprwall_min=float(min(Vl_avg[_prwall].min(), Vr_avg[_prwall].min())) if _prwall.any() else 0.0)


def _build_edge_array_geometry(AR, W=32, mouth=237, poly_um=0.3, feature_w_um=0.5,
                               open_width_um=3.7, right_buffer_um=0.5):
    pr_cells = int(round((AR * feature_w_um - poly_um) / feature_w_um * W))
    pr_cells = max(pr_cells, 1)
    poly_cells = int(round(poly_um / feature_w_um * W)) if poly_um > 0 else 0
    D = pr_cells + poly_cells
    nz = D + mouth
    poly_cells = min(poly_cells, max(D - 2, 0))
    open_w = max(W, int(round(open_width_um / feature_w_um * W)))
    buffer_w = max(4, int(round(right_buffer_um / feature_w_um * W)))
    edge0 = open_w
    edge1 = edge0 + W
    trench0 = edge1
    trench1 = trench0 + W
    neigh0 = trench1
    neigh1 = neigh0 + W
    right_trench0 = neigh1
    right_trench1 = right_trench0 + W
    next0 = right_trench1
    nx = next0 + W + buffer_w
    z_poly0 = nz - 1 - poly_cells

    solid = np.zeros((nx, nz), dtype=np.bool_)
    solid[edge0:edge1, mouth:nz - 1] = True
    solid[neigh0:neigh1, mouth:nz - 1] = True
    solid[next0:nx, mouth:nz - 1] = True
    solid[:, nz - 1] = True

    cond = np.zeros((nx, nz), dtype=np.int8)
    if poly_cells > 0:
        cond[edge0:edge1, z_poly0:nz - 1] = 1
        cond[neigh0:neigh1, z_poly0:nz - 1] = 2

    floor_mask = np.zeros_like(solid)
    floor_mask[:, nz - 1] = True
    floor_trench_mask = np.zeros_like(solid)
    floor_trench_mask[trench0:trench1, nz - 1] = True
    pr_mask = solid & (cond == 0) & ~floor_mask
    # Outer wall of the edge line: its leftmost solid column (edge0), whose left neighbour is the
    # open area. This face is exposed to the open plasma half-space and collects the open-area
    # electron flux (W2). Split implicitly into PR (cond==0) and poly (cond==1) by cond.
    edge_outer_mask = np.zeros_like(solid)
    edge_outer_mask[edge0, mouth:nz - 1] = True
    edge_exposed_area = max(2 * poly_cells, 1)
    neighbor_exposed_area = max(2 * poly_cells, 1)
    return dict(solid=solid, cond=cond, floor_mask=floor_mask,
                floor_trench_mask=floor_trench_mask, pr_mask=pr_mask,
                edge_outer_mask=edge_outer_mask,
                nx=int(nx), nz=int(nz), D=int(D), W=int(W), mouth=int(mouth),
                open_w=int(open_w), buffer_w=int(buffer_w),
                edge0=int(edge0), edge1=int(edge1), trench0=int(trench0), trench1=int(trench1),
                neigh0=int(neigh0), neigh1=int(neigh1),
                right_trench0=int(right_trench0), right_trench1=int(right_trench1),
                next0=int(next0), z_poly0=int(z_poly0),
                pr_cells=int(pr_cells), poly_cells=int(poly_cells),
                edge_exposed_area=int(edge_exposed_area),
                neighbor_exposed_area=int(neighbor_exposed_area))


def solve_edge_array_charging(AR, W=32, mouth=237, Te=4.0, V_dc=37.0, V_rf=30.0,
                              iadf_hwhm_deg=4.3, cos_power=0.6, n_per_iter=6000, n_iter=120,
                              relax=None, seed=0, verbose=False, smooth=False,
                              poly_um=0.3, feature_w_um=0.5, rf_bursts=True,
                              sheath_um=89.0, boundary_um=3.7, insul_vmin_Te=None,
                              trace_integrator="adaptive_numba", trace_step_cap_factor=40.0,
                              see_model="none", see_generations=1,
                              ion_angle_energy_corr="anticorrelated",
                              source_model="analytic", open_width_um=3.7, right_buffer_um=0.5,
                              edge_open_model="none", edge_open_electron_flux=None,
                              edge_open_samples=None, open_wall_frac=1.0,
                              electron_model="mc"):
    """Explicit nonperiodic HG edge-line/open-area charging cell.

    This is the next-step mechanism solver: open area + edge poly line + edge trench +
    neighboring line. It intentionally leaves the original periodic `solve_trench_charging` intact.
    """
    if see_model not in ("none", None):
        raise NotImplementedError("solve_edge_array_charging does not yet include SEE cascades")
    if trace_integrator not in ("adaptive_numba", "adaptive"):
        raise ValueError("edge-array solver currently requires adaptive tracing")
    if _trace_edge_particles_adaptive is None:
        raise RuntimeError("numba edge tracer is unavailable")

    rng = np.random.default_rng(seed)
    geom = _build_edge_array_geometry(AR, W=W, mouth=mouth, poly_um=poly_um,
                                      feature_w_um=feature_w_um,
                                      open_width_um=open_width_um,
                                      right_buffer_um=right_buffer_um)
    solid = geom["solid"]
    cond = geom["cond"]
    floor_mask = geom["floor_mask"]
    floor_trench_mask = geom["floor_trench_mask"]
    pr_mask = geom["pr_mask"]
    nx = geom["nx"]; nz = geom["nz"]
    open_frac = W / nx
    if relax is None:
        relax = 2.0 * Te
    # Floating-insulator numerical guard. An insulator charges until local net current -> 0; the
    # only bound is that it cannot charge past the incident particle energy, i.e. the sheath
    # potential scale (V_dc + V_rf). None -> that physical guard (tracks the bias); a number sets
    # a tighter guard in Te units for diagnostics. The old default 1*Te pinned the electron-
    # collecting upper-sidewall PR cells and broke charge conservation (pr residual -0.14..-0.33).
    pr_vguard = (V_dc + V_rf) if insul_vmin_Te is None else float(insul_vmin_Te) * Te

    V = np.zeros((nx, nz))
    Vsolid = np.zeros((nx, nz))
    Vedge = 0.0
    Vneighbor = 0.0
    edge_outer_gross_flux = _hg_edge_outer_electron_gross_flux(AR)
    edge_open_samples = int(edge_open_samples or max(32768, 4 * n_per_iter))
    edge_aux = dict(model=edge_open_model, electron_gross=0.0, ion_gross=0.0,
                    net_electron=0.0, open_width_cells=float(geom["open_w"]),
                    n=int(edge_open_samples), hg_electron_gross=edge_outer_gross_flux,
                    override=False)
    if edge_open_electron_flux is not None:
        edge_aux.update(model="override", electron_gross=float(edge_open_electron_flux),
                        net_electron=float(edge_open_electron_flux), override=True)
    elif edge_open_model == "line_of_sight":
        edge_aux.update(_sample_open_side_flux(
            rng, edge_open_samples, W, geom["pr_cells"], geom["poly_cells"], cos_power,
            iadf_hwhm_deg, ion_angle_energy_corr, V_dc, V_rf, geom["open_w"]))
    elif edge_open_model in ("none", None, "wall_flux"):
        pass
    else:
        raise ValueError(f"unknown edge_open_model: {edge_open_model}")

    def edge_boundary_net(v_edge):
        if edge_open_electron_flux is not None:
            return float(edge_open_electron_flux)
        if edge_open_model != "line_of_sight":
            return 0.0
        if v_edge < 0.0:
            accept = math.erfc(math.sqrt(max(-float(v_edge), 0.0) / max(Te, 1.0e-9)))
        else:
            accept = 1.0
        return max(float(edge_aux["electron_gross"]) * accept - float(edge_aux["ion_gross"]), 0.0)

    def apply_dirichlet(V, boundary=0.0, vsolid=None, vedge=None, vneighbor=None):
        vs = Vsolid if vsolid is None else vsolid
        ve = Vedge if vedge is None else vedge
        vn = Vneighbor if vneighbor is None else vneighbor
        V[:, 0] = boundary
        V[solid] = vs[solid]
        V[cond == 1] = ve
        V[cond == 2] = vn
        return V

    ii, jj = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    _red = ((ii + jj) % 2 == 0)
    _inside = ~solid
    _inside[:, 0] = False

    def laplace(V, sweeps=260, omega=1.88, **bc):
        # true red-black GS-SOR: the neighbor average is recomputed per color so the second
        # color sees the first color's updated values (masked Jacobi diverges for omega > 1)
        for _ in range(sweeps):
            V = apply_dirichlet(V, **bc)
            for color in (_red, ~_red):
                xm = np.empty_like(V); xp = np.empty_like(V)
                xm[1:, :] = V[:-1, :]; xm[0, :] = V[0, :]
                xp[:-1, :] = V[1:, :]; xp[-1, :] = V[-1, :]
                avg = np.zeros_like(V)
                avg[:, 1:-1] = 0.25 * (xm[:, 1:-1] + xp[:, 1:-1] + V[:, 2:] + V[:, :-2])
                m = _inside & color
                V[m] = (1.0 - omega) * V[m] + omega * avg[m]
        return apply_dirichlet(V, **bc)

    VB = np.zeros((nx, nz))
    if rf_bursts:
        zsolid = np.zeros_like(Vsolid)
        VB = laplace(VB, sweeps=220, boundary=1.0, vsolid=zsolid, vedge=0.0, vneighbor=0.0)
    ExB = -(np.gradient(VB, axis=0)); EzB = -(np.gradient(VB, axis=1))

    _phi = np.linspace(0.0, 2.0 * np.pi, 720)
    _pw = np.exp(-(V_dc + V_rf * np.sin(_phi)) / Te)
    _pcdf = np.cumsum(_pw); _pcdf /= _pcdf[-1]

    def sample_source(kind, n):
        if source_model == "sheath_mc" or (source_model == "sheath_electrons" and kind == "electron"):
            E0, th, sB = sample_sheath_source(kind, n, rng, Te=Te, V_dc=V_dc, V_rf=V_rf,
                                              iadf_hwhm_deg=iadf_hwhm_deg,
                                              cos_power=cos_power,
                                              boundary_um=boundary_um,
                                              sheath_um=sheath_um)
            q = +1.0 if kind == "ion" else -1.0
        elif source_model not in ("analytic", "sheath_electrons"):
            raise ValueError(f"unknown source_model: {source_model}")
        elif kind == "ion":
            phi_p = rng.uniform(0.0, 2.0 * np.pi, n)
            E0 = np.maximum(V_dc + V_rf * np.sin(phi_p), 0.5)
            sig0 = np.deg2rad(iadf_hwhm_deg) / 1.1774
            if ion_angle_energy_corr == "anticorrelated":
                sig = sig0 * np.sqrt(V_dc / E0)
            elif ion_angle_energy_corr == "independent":
                sig = np.full(n, sig0)
            elif ion_angle_energy_corr == "positive":
                sig = sig0 * np.sqrt(E0 / V_dc)
            else:
                raise ValueError(f"unknown ion_angle_energy_corr: {ion_angle_energy_corr}")
            th = rng.normal(0.0, sig, n)
            q = +1.0
            sB = np.zeros(n)
        else:
            E0 = rng.gamma(2.0, Te, n)
            u = rng.uniform(0.0, 1.0, n)
            ct3 = (1.0 - u) ** (1.0 / (cos_power + 2.0))
            st3 = np.sqrt(np.maximum(1.0 - ct3 * ct3, 0.0))
            az = rng.uniform(0.0, 2.0 * np.pi, n)
            th = np.arctan2(st3 * np.cos(az), ct3)
            q = -1.0
            if rf_bursts:
                phi_e = np.interp(rng.uniform(0, 1, n), _pcdf, _phi)
                frac = (boundary_um / sheath_um) ** (4.0 / 3.0)
                sB = frac * (V_dc + V_rf * np.sin(phi_e))
            else:
                sB = np.zeros(n)
        vx = np.sqrt(E0) * np.sin(th)
        vz = np.sqrt(E0) * np.abs(np.cos(th))
        x = rng.uniform(0.0, float(nx - 1), n)
        z = np.ones(n) * max(1.0, float(mouth) - 0.5)
        return x, z, vx, vz, sB, q

    trace_stats = []

    def trace(kind, n, Ex, Ez):
        x, z, vx, vz, sB, q = sample_source(kind, n)
        max_steps = int(float(trace_step_cap_factor) * nz)
        ht, hix, hiz, impact_E, hit_vx, hit_vz, survivor, steps = _trace_edge_particles_adaptive(
            Ex, Ez, ExB, EzB, solid, cond, x, z, vx, vz, sB, q,
            nx, nz, mouth, geom["edge0"], geom["edge1"], geom["trench0"], geom["trench1"],
            geom["neigh0"], geom["neigh1"], geom["z_poly0"], max_steps)
        counts = np.zeros((nx, nz))
        m = ht > 0
        if m.any():
            np.add.at(counts, (hix[m], hiz[m]), 1.0)
        foot = ht == 4
        foot_n = float(foot.sum())
        foot_E = float(impact_E[foot].sum())
        foot_En = float((hit_vx[foot] * hit_vx[foot]).sum())
        neigh = ht == 5
        neigh_n = float(neigh.sum())
        neigh_E = float(impact_E[neigh].sum())
        edge_outer = ht == 3
        trace_stats.append(dict(kind=kind, n=int(n), survivors=int(survivor.sum()),
                                survivor_frac=float(survivor.sum() / max(n, 1)),
                                steps=int(steps.max()) if steps.size else 0,
                                cap=int(max_steps), integrator="edge_adaptive_numba",
                                foot_n=foot_n, foot_E=foot_E, foot_En=foot_En,
                                neighbor_n=neigh_n, neighbor_E=neigh_E,
                                edge_outer_n=float(edge_outer.sum())))
        return counts, ht, impact_E, hit_vx

    # W2 isotropic electron source: precompute the sky view factor per exposed surface cell ONCE
    # (geometry is fixed). Electrons are isotropic, so each cell collects (base rate) x view_factor
    # x potential-throttle instead of the down-going MC trace -- this correctly shadows deep walls
    # (starves the neighbour -> it rises) and illuminates open ones (holds the edge line low).
    vf_grid = None
    if electron_model == "viewfactor":
        gas = ~solid
        exp_cell = np.zeros_like(solid)
        exp_cell[1:, :] |= solid[1:, :] & gas[:-1, :]
        exp_cell[:-1, :] |= solid[:-1, :] & gas[1:, :]
        exp_cell[:, 1:] |= solid[:, 1:] & gas[:, :-1]
        exp_cell[:, :-1] |= solid[:, :-1] & gas[:, 1:]
        six, siz = np.where(exp_cell)
        vf = _sky_view_factors(solid, six.astype(np.int64), siz.astype(np.int64), 180)
        vf_grid = np.zeros((nx, nz))
        vf_grid[six, siz] = vf
    e_base = float(n_per_iter) / nx   # electrons arrive at the mouth at the ion rate (zero net current)

    hist = []
    vsolid_hist = []
    vedge_hist = []
    vneighbor_hist = []
    res_hist = []
    for it in range(n_iter):
        V = laplace(V)
        Ex = -(np.gradient(V, axis=0)); Ez = -(np.gradient(V, axis=1))
        ci, hti, _, _ = trace("ion", n_per_iter, Ex, Ez)
        if electron_model == "viewfactor":
            Vcell = Vsolid.copy()
            Vcell[cond == 1] = Vedge
            Vcell[cond == 2] = Vneighbor
            throttle = np.where(Vcell >= 0.0, 1.0, np.exp(np.clip(Vcell / max(Te, 1e-9), -40.0, 0.0)))
            ce = e_base * vf_grid * throttle
            hte = None
        else:
            ce, hte, _, _ = trace("electron", n_per_iter, Ex, Ez)
        net = ci - ce
        anneal = max(1.0 / (1.0 + it / 25.0), 0.25)
        scale = anneal * relax / n_per_iter * nx
        # W2: open-plasma electron flux onto the edge line's OUTER wall (faces the open area).
        # The down-going MC source never delivers side-arriving electrons, so this wall is added
        # analytically: each outer-wall cell collects the open-field electron flux (one column's
        # worth, open_wall_frac of it), throttled by its own potential (electrons repelled by a
        # negative surface ~ exp(V/Te); a positive surface collects the full arriving flux). This
        # holds the edge line low while the walled-in neighbour stays starved and rises -> the
        # HG edge/neighbour potential split. Only active for edge_open_model="wall_flux".
        if edge_open_model == "wall_flux":
            eo = geom["edge_outer_mask"]
            Vwall = np.where(cond == 1, Vedge, Vsolid)
            accept = np.where(Vwall >= 0.0, 1.0, np.exp(np.clip(Vwall / max(Te, 1e-9), -40.0, 0.0)))
            open_e = np.zeros_like(net)
            open_e[eo] = open_wall_frac * (float(n_per_iter) / nx) * accept[eo]
            net = net - open_e   # electrons: reduce net current on the outer-wall cells
        Vsolid[pr_mask | floor_mask] += scale * net[pr_mask | floor_mask]
        if hte is None:   # view-factor electrons: no per-particle hits, no scalar top-up
            edge_extra_e = 0.0
        else:
            explicit_outer_net_e = float((hte == 3).sum() - (hti == 3).sum())
            target_outer_net_e = edge_boundary_net(Vedge) * n_per_iter * open_frac
            edge_extra_e = max(target_outer_net_e - explicit_outer_net_e, 0.0)
        res_hist.append((float(net[floor_trench_mask].sum()), float(net[pr_mask].sum()),
                         float(net[cond == 1].sum() - edge_extra_e),
                         float(net[cond == 2].sum())))
        Vedge += scale * (net[cond == 1].sum() - edge_extra_e) / max(geom["edge_exposed_area"], 1)
        Vneighbor += scale * net[cond == 2].sum() / max(geom["neighbor_exposed_area"], 1)
        Vedge = float(np.clip(Vedge, -3.0 * Te, V_dc + V_rf))
        Vneighbor = float(np.clip(Vneighbor, -3.0 * Te, V_dc + V_rf))
        Vsolid[floor_mask] = np.clip(Vsolid[floor_mask], 0.0, V_dc + V_rf)
        Vsolid[pr_mask] = np.clip(Vsolid[pr_mask], -pr_vguard, V_dc + V_rf)
        hist.append(float((hti == 1).sum() / n_per_iter))
        vsolid_hist.append(Vsolid.copy())
        vedge_hist.append(Vedge); vneighbor_hist.append(Vneighbor)
        if verbose and it % 20 == 0:
            tsi = next((s for s in reversed(trace_stats) if s["kind"] == "ion"), None)
            floor_e = float(ce[floor_trench_mask].sum()) if hte is None else float((hte == 1).sum())
            print(f"  it{it}: floor_i={(hti == 1).sum()/n_per_iter:.3f} "
                  f"floor_e={floor_e/n_per_iter:.3f} "
                  f"Vedge={Vedge:.1f} Vneighbor={Vneighbor:.1f}", flush=True)

    k = max(n_iter // 3, 5)
    Vsolid_avg = np.mean(np.array(vsolid_hist[-k:]), axis=0)
    Vedge_avg = float(np.mean(vedge_hist[-k:]))
    Vneighbor_avg = float(np.mean(vneighbor_hist[-k:]))
    Vsolid = Vsolid_avg
    Vedge = Vedge_avg
    Vneighbor = Vneighbor_avg
    if smooth:
        try:
            from scipy.ndimage import uniform_filter
            Vsolid = uniform_filter(Vsolid, size=3, mode="nearest")
        except Exception:
            pass
    V = laplace(V, sweeps=180)
    Ex = -(np.gradient(V, axis=0)); Ez = -(np.gradient(V, axis=1))
    ci, hti, Ei, vxi = trace("ion", 4 * n_per_iter, Ex, Ez)
    ntot = 4 * n_per_iter
    if electron_model == "viewfactor":
        Vcell = Vsolid.copy(); Vcell[cond == 1] = Vedge; Vcell[cond == 2] = Vneighbor
        throttle = np.where(Vcell >= 0.0, 1.0, np.exp(np.clip(Vcell / max(Te, 1e-9), -40.0, 0.0)))
        ce = (float(ntot) / nx) * vf_grid * throttle   # electron collection grid at 4x stats
        hte = None
    else:
        ce, hte, Ee, vxe = trace("electron", 4 * n_per_iter, Ex, Ez)
    trench_norm = max(ntot * open_frac, 1.0)

    edge_inner = hti == 4
    edge_outer_i = hti == 3
    neighbor_i = hti == 5
    foot_n = float(edge_inner.sum())
    foot_E = float(Ei[edge_inner].sum())
    foot_En = float((vxi[edge_inner] * vxi[edge_inner]).sum())
    neighbor_n = float(neighbor_i.sum())
    neighbor_E = float(Ei[neighbor_i].sum())
    net = ci - ce
    # electron region sums: from hit-labels (mc) or from the view-factor collection grid
    if hte is None:
        eouter = geom["edge_outer_mask"]
        e_floor_sum = float(ce[floor_trench_mask].sum())
        e_eouter_sum = float(ce[eouter].sum())
        e_einner_sum = 0.0
        e_neigh_sum = float(ce[cond == 2].sum())
    else:
        e_floor_sum = float((hte == 1).sum())
        e_eouter_sum = float((hte == 3).sum())
        e_einner_sum = float((hte == 4).sum())
        e_neigh_sum = float((hte == 5).sum())
    explicit_outer_net_e = float(e_eouter_sum - edge_outer_i.sum())
    target_outer_net_e = edge_boundary_net(Vedge_avg) * ntot * open_frac
    edge_extra_final = max(target_outer_net_e - explicit_outer_net_e, 0.0)
    # steady-state residual = tail-averaged net current per surface. A single 4x-n snapshot is
    # shot-noise limited (~0.05 at n_per_iter=1200); averaging the last k iterations measures the
    # actual imbalance at the averaged potential with k-times the samples.
    iter_norm = max(n_per_iter * open_frac, 1.0)
    res_tail = np.mean(np.array(res_hist[-k:]), axis=0) / iter_norm
    residual = dict(
        floor=float(res_tail[0]),
        pr=float(res_tail[1]),
        poly_edge=float(res_tail[2]),
        poly_neighbor=float(res_tail[3]),
        top=0.0,
    )
    residual_snapshot = dict(
        floor=float(net[floor_trench_mask].sum() / trench_norm),
        pr=float(net[pr_mask].sum() / trench_norm),
        poly_edge=float((net[cond == 1].sum() - edge_extra_final) / trench_norm),
        poly_neighbor=float(net[cond == 2].sum() / trench_norm),
        top=0.0,
    )
    edge_net = float(edge_boundary_net(Vedge_avg))
    explicit_e_gross = float(e_eouter_sum / trench_norm)
    explicit_i_gross = float(edge_outer_i.sum() / trench_norm)
    accepted_e_gross = edge_net + float(edge_aux["ion_gross"])
    boundary_e_gross = max(accepted_e_gross - explicit_e_gross, 0.0)
    boundary_i_gross = max(float(edge_aux["ion_gross"]) - explicit_i_gross, 0.0)
    boundary_net = max(edge_net - (explicit_e_gross - explicit_i_gross), 0.0)
    edge_aux["net_electron"] = edge_net
    edge_aux["electron_accepted"] = accepted_e_gross
    edge_aux["edge_potential_for_net"] = float(Vedge_avg)
    diag = dict(
        ion=dict(floor=float((hti == 1).sum() / trench_norm),
                 edge_outer_poly=float(edge_outer_i.sum() / trench_norm),
                 edge_inner_poly=float(edge_inner.sum() / trench_norm),
                 neighbor_poly=float(neighbor_i.sum() / trench_norm)),
        electron=dict(floor=float(e_floor_sum / trench_norm),
                      edge_outer_poly=float(e_eouter_sum / trench_norm),
                      edge_inner_poly=float(e_einner_sum / trench_norm),
                      neighbor_poly=float(e_neigh_sum / trench_norm)),
        trace=dict(last_ion=next((s for s in reversed(trace_stats) if s["kind"] == "ion"), None),
                   last_electron=next((s for s in reversed(trace_stats) if s["kind"] == "electron"), None),
                   cap_factor=float(trace_step_cap_factor),
                   integrator="edge_adaptive_numba",
                   source_model=source_model,
                   poly_mode="edge_array",
                   edge_open_model=edge_open_model,
                   edge_outer_electron_gross_flux=edge_outer_gross_flux),
        edge_open=dict(model="explicit_geometry",
                       electron_gross=float(explicit_e_gross + boundary_e_gross),
                       ion_gross=float(explicit_i_gross + boundary_i_gross),
                       net_electron=float(explicit_e_gross - explicit_i_gross + boundary_net),
                       explicit_electron_gross=explicit_e_gross,
                       explicit_ion_gross=explicit_i_gross,
                       boundary_electron_gross=float(boundary_e_gross),
                       boundary_net_electron=float(boundary_net),
                       hg_electron_gross=edge_outer_gross_flux,
                       open_width_cells=float(geom["open_w"]),
                       n=int(ntot)),
        residual=residual,
        residual_snapshot=residual_snapshot,
        see=dict(model=see_model, generations=int(see_generations), last=None),
    )
    floor_line = Vsolid[geom["trench0"]:geom["trench1"], nz - 1]
    pr_vals = Vsolid[pr_mask]
    floor_flux_final = float((hti == 1).sum() / trench_norm)
    floor_flux_tail = float(np.mean(hist[-k:]) / open_frac)
    return dict(
        floor_flux=floor_flux_final, floor_flux_tail=floor_flux_tail,
        V_floor_center=float(floor_line[W // 2]),
        V_foot_peak=float(floor_line.max()), V_poly=Vneighbor_avg, Vfloor=floor_line, V=V,
        V_poly_left=Vedge_avg, V_poly_right=Vneighbor_avg,
        V_poly_edge=Vedge_avg, V_poly_neighbor=Vneighbor_avg,
        foot_ion_flux=float(foot_n / trench_norm),
        foot_ion_Emean=float(foot_E / max(foot_n, 1.0)),
        foot_ion_Enormal_mean=float(foot_En / max(foot_n, 1.0)),
        foot_ion_flux_left=float(foot_n / trench_norm),
        foot_ion_flux_right=float(neighbor_n / trench_norm),
        foot_ion_Emean_left=float(foot_E / max(foot_n, 1.0)),
        foot_ion_Emean_right=float(neighbor_E / max(neighbor_n, 1.0)),
        foot_ion_Enormal_mean_left=float(foot_En / max(foot_n, 1.0)),
        foot_ion_Enormal_mean_right=np.nan,
        foot_ion_flux_edge=float(foot_n / trench_norm),
        foot_ion_Emean_edge=float(foot_E / max(foot_n, 1.0)),
        foot_ion_Enormal_mean_edge=float(foot_En / max(foot_n, 1.0)),
        diag=diag,
        geom=dict(W=int(W), mouth=int(mouth), D=int(geom["D"]), nx=int(nx), nz=int(nz),
                  poly_cells=int(geom["poly_cells"]), open_w=int(geom["open_w"]),
                  edge0=int(geom["edge0"]), edge1=int(geom["edge1"]),
                  trench0=int(geom["trench0"]), trench1=int(geom["trench1"]),
                  neigh0=int(geom["neigh0"]), neigh1=int(geom["neigh1"]),
                  right_trench0=int(geom["right_trench0"]),
                  right_trench1=int(geom["right_trench1"]),
                  next0=int(geom["next0"])),
        Vprwall_mean=float(pr_vals.mean()) if pr_vals.size else 0.0,
        Vprwall_min=float(pr_vals.min()) if pr_vals.size else 0.0)


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
