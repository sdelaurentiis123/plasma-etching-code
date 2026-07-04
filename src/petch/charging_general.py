"""General 2-D feature charging: charged particles in a self-consistent field, ANY geometry.

This is the geometry-agnostic engine the edge-array/trench solvers should have been. There are NO
named features ("edge line", "neighbour"): the input is just a material-tagged grid

    mat[ix, iz] in {GAS, INSULATOR, CONDUCTOR}   plus  cid[ix, iz] = conductor component id

and the physics is the same for a trench, a hole cross-section, a line-and-space edge array, or a
real device: launch ions (directional, from the sheath) and electrons (ISOTROPIC, Lambert flux
through the plasma boundary), push both through the Laplace field, deposit charge where they land,
let each insulator cell float to local current balance and each connected conductor float as one
equipotential, re-solve the field, iterate to steady state.

Electrons are TRACED (not a view-factor shortcut) so the electrostatic focusing HG requires --
electrons pulled into the positive trench, reaching the floor more than pure geometry allows -- is
captured self-consistently. The only thing that made the old electron source wrong was launching a
DOWN-going hemisphere at the mouth; here electrons enter isotropically from the top boundary.

Kushner MCFPM lineage; this is the W5 shared-interface direction from ROBUST_PHYSICS_MODEL_PLAN.
"""
from __future__ import annotations

import math

import numpy as np

try:
    from numba import njit, prange
except Exception:  # pragma: no cover
    njit = None
    prange = range

from .charging2d import _sky_view_factors

GAS = 0
INSULATOR = 1
CONDUCTOR = 2


def _cryo_conductivity(T_C):
    """Temperature-gated surface conductivity 0..1. Physical anchor: below ~-60 C an HF/H2O layer
    condenses and raises insulator surface conductivity 3-6 orders of magnitude (Appl. Phys. Lett.
    123, 212106, 2023), which dissipates feature charge. Sigmoid switch-on around the condensation
    onset (~-60 C), saturating cold. Warm (>0 C) -> ~0 (insulator holds charge)."""
    T_onset = -60.0
    width = 15.0
    return float(1.0 / (1.0 + np.exp((T_C - T_onset) / width)))


def _trace_general_py(Ex, Ez, solid, x0, z0, vx0, vz0, q, nx, nz, max_steps):
    """Push particles through the field until they enter a solid cell; return the hit CELL.

    No feature labels -- just (hit_ix, hit_iz) and the impact kinematics. Leapfrog with adaptive dt
    (<=0.45 cell/step), periodic x. z<0.5 = escaped back to plasma (survivor)."""
    n = x0.shape[0]
    hit_ix = np.full(n, -1, np.int64)
    hit_iz = np.full(n, -1, np.int64)
    impact_E = np.zeros(n)
    hit_vx = np.zeros(n)
    survivor = np.zeros(n, np.uint8)
    xmax = float(nx)
    for p in prange(n):
        x = x0[p]; z = z0[p]; vx = vx0[p]; vz = vz0[p]
        alive = True
        for _ in range(max_steps):
            ix = int(x); iz = int(z)
            if ix < 0:
                ix = 0
            elif ix > nx - 2:
                ix = nx - 2
            if iz < 0:
                iz = 0
            elif iz > nz - 2:
                iz = nz - 2
            fx = Ex[ix, iz]; fz = Ez[ix, iz]
            ax = q * fx * 0.5; az = q * fz * 0.5
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
            ix2 = int(xa); iz2 = int(za)
            if ix2 < 0:
                ix2 = 0
            elif ix2 > nx - 2:
                ix2 = nx - 2
            if iz2 < 0:
                iz2 = 0
            elif iz2 > nz - 2:
                iz2 = nz - 2
            vx = vx_half + 0.25 * q * Ex[ix2, iz2] * dt
            vz = vz_half + 0.25 * q * Ez[ix2, iz2] * dt
            x = xa % xmax
            z = za
            if z < 0.5:
                alive = False
                break
            ixh = int(x); izh = int(z)
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
                alive = False
                break
        if alive:
            survivor[p] = 1
    return hit_ix, hit_iz, impact_E, hit_vx, survivor


_trace_general = (njit(cache=True, parallel=True, fastmath=True)(_trace_general_py)
                  if njit is not None else _trace_general_py)


def _connected_conductor_ids(mat):
    """Label connected components of CONDUCTOR cells (4-connectivity). Each component is one
    floating equipotential. Returns an int grid (0 = not a conductor, >=1 = component id)."""
    try:
        from scipy.ndimage import label
        cid, ncomp = label(mat == CONDUCTOR)
        return cid.astype(np.int64), int(ncomp)
    except Exception:
        # tiny fallback: every conductor cell its own id (correct if conductors are pre-separated)
        cid = np.zeros(mat.shape, np.int64)
        idx = np.where(mat == CONDUCTOR)
        cid[idx] = np.arange(1, idx[0].size + 1)
        return cid, int(idx[0].size)


def sample_ions(n, rng, mouth, nx, V_dc, V_rf, iadf_hwhm_deg):
    """Directional ions from the sheath: energy over the RF cycle, narrow angular cone, launched
    at the mouth plane moving down. (Ion optics unchanged from the validated model.)"""
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    E0 = np.maximum(V_dc + V_rf * np.sin(phi), 0.5)
    sig = np.deg2rad(iadf_hwhm_deg) / 1.1774 * np.sqrt(V_dc / E0)
    th = rng.normal(0.0, sig, n)
    vx = np.sqrt(E0) * np.sin(th)
    vz = np.sqrt(E0) * np.abs(np.cos(th))
    x = rng.uniform(0.0, float(nx - 1), n)
    z = np.full(n, max(1.0, float(mouth) - 0.5))
    return x, z, vx, vz


def sample_electrons(n, rng, nx, Te, cos_power=1.0):
    """Isotropic thermal electrons entering through the plasma boundary (top, z~=1): thermal energy,
    Lambert (cos^p) flux distribution so they arrive at ALL angles including near-horizontal. The
    field then focuses them into positive features (the piece a down-going mouth-launch threw away).
    Launched at the top so they traverse the field region and reach outward-facing walls."""
    E0 = rng.gamma(2.0, Te, n)
    u = rng.uniform(0.0, 1.0, n)
    ct = (1.0 - u) ** (1.0 / (cos_power + 2.0))
    st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
    az = rng.uniform(0.0, 2.0 * np.pi, n)
    th = np.arctan2(st * np.cos(az), ct)
    vx = np.sqrt(E0) * np.sin(th)
    vz = np.sqrt(E0) * np.abs(np.cos(th))
    x = rng.uniform(0.0, float(nx - 1), n)
    z = np.full(n, 1.0)
    return x, z, vx, vz


def solve_charging(mat, mouth, Te=4.0, V_dc=37.0, V_rf=30.0, iadf_hwhm_deg=4.3,
                   n_per_iter=6000, n_iter=200, relax=None, seed=0,
                   insul_vguard=None, verbose=False,
                   field_model="laplace", eps_insulator=3.9, rho_coupling=1.0,
                   electron_open_vf=True, frame_every=0,
                   electron_model="trace", vf_focus=1.8, vf_focus_pot=0.0,
                   surface_conductivity=0.0, temperature_C=20.0):
    """Steady-state feature charging for ANY material grid `mat` (GAS/INSULATOR/CONDUCTOR).

    mat: (nx, nz) int grid. z=0 is the plasma boundary (Dirichlet 0), z increases into the wafer.
    Returns V (potential), per-cell insulator potential, and per-conductor equipotentials."""
    if _trace_general is None:
        raise RuntimeError("numba unavailable")
    rng = np.random.default_rng(seed)
    nx, nz = mat.shape
    solid = mat != GAS
    insul = mat == INSULATOR
    cid, ncomp = _connected_conductor_ids(mat)
    if relax is None:
        relax = 2.0 * Te
    if insul_vguard is None:
        insul_vguard = V_dc + V_rf

    Vs = np.zeros((nx, nz))          # per-cell insulator potential (laplace mode)
    Vc = np.zeros(ncomp + 1)         # per-conductor-component equipotential (index 0 unused)
    rho = np.zeros((nx, nz))         # per-cell free charge (poisson mode: on insulators)

    ii, jj = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    red = ((ii + jj) % 2 == 0)
    inside = ~solid
    inside[:, 0] = False

    # --- MCFPM-style variable-eps Poisson setup (field_model="poisson") ---
    # WARNING: this solves the field THROUGH the solid dielectric interiors. That is physical only
    # with a properly grounded substrate under the oxide (as MCFPM sets up). In this reduced geometry
    # (floating insulator blocks, no substrate ground, Neumann side boundary) the block interiors run
    # away to -100..-245 V -- a non-physical artifact. Default field_model="laplace" solves the field
    # in the GAS ONLY with solids as charged surface boundaries (HG's method); nothing runs through
    # the silicon and there is no artifact. Use "poisson" only once a grounded substrate layer is added.
    # Per-cell dielectric constant: gas 1, insulator eps_insulator, conductor Dirichlet (equipot,
    # the high-mobility limit). The Poisson domain SOLVES gas AND insulator cells (the field
    # penetrates the dielectric, carrying inter-feature coupling); conductors and the grounded
    # top/bottom are Dirichlet. This is Kushner MCFPM's div(eps grad phi) = -rho with per-cell eps.
    eps = np.ones((nx, nz))
    eps[insul] = float(eps_insulator)
    is_cond = cid > 0
    # face conductances (arithmetic-mean eps on each face)
    e_xm = np.empty_like(eps); e_xp = np.empty_like(eps)
    e_zm = np.empty_like(eps); e_zp = np.empty_like(eps)
    e_xm[1:, :] = 0.5 * (eps[1:, :] + eps[:-1, :]); e_xm[0, :] = eps[0, :]
    e_xp[:-1, :] = 0.5 * (eps[:-1, :] + eps[1:, :]); e_xp[-1, :] = eps[-1, :]
    e_zm[:, 1:] = 0.5 * (eps[:, 1:] + eps[:, :-1]); e_zm[:, 0] = eps[:, 0]
    e_zp[:, :-1] = 0.5 * (eps[:, :-1] + eps[:, 1:]); e_zp[:, -1] = eps[:, -1]
    e_diag = e_xm + e_xp + e_zm + e_zp
    e_diag[e_diag < 1e-9] = 1e-9
    poisson_inside = (~is_cond)
    poisson_inside[:, 0] = False          # grounded plasma boundary (top) is the only Dirichlet

    # --- open-plasma electron floor via sky view factor (closes the edge/neighbour split) ---
    # The down-going electron trace under-samples OPEN-FACING walls (the edge line's outer wall
    # facing the open area). For every surface cell, the electron collection should be at LEAST the
    # isotropic open-plasma flux it can see = e_base * view_factor. Interior cells (deep walls, floor)
    # get MORE from the field-focused trace, so max() keeps the trace there; open-facing cells (high
    # VF) get their view-factor flux, which holds the edge line low while the walled-in neighbour
    # stays starved and rises. General (VF-based), no hardcoded geometry.
    vf_grid = None
    if electron_open_vf:
        gas_c = ~solid
        exp_cell = np.zeros_like(solid)
        exp_cell[1:, :] |= solid[1:, :] & gas_c[:-1, :]
        exp_cell[:-1, :] |= solid[:-1, :] & gas_c[1:, :]
        exp_cell[:, 1:] |= solid[:, 1:] & gas_c[:, :-1]
        exp_cell[:, :-1] |= solid[:, :-1] & gas_c[:, 1:]
        six, siz = np.where(exp_cell)
        vf = _sky_view_factors(solid, six.astype(np.int64), siz.astype(np.int64), 180)
        vf_grid = np.zeros((nx, nz))
        vf_grid[six, siz] = vf
    # NB: the bottom row is the SiO2 FLOOR (an insulator) and must FLOAT with its charge, not be
    # grounded -- it charges positive relative to the plasma, which is the whole point. Zero-gradient
    # (insulating) closure below it. (HG references to the substrate a few cells lower; the relative
    # floor/edge/neighbour structure is reference-independent.)

    def apply_bc(V):
        V[:, 0] = 0.0
        V[insul] = Vs[insul]
        if ncomp > 0:
            V[cid > 0] = Vc[cid[cid > 0]]
        return V

    def laplace(V, sweeps=200, omega=1.88):
        for _ in range(sweeps):
            V = apply_bc(V)
            for color in (red, ~red):
                xm = np.empty_like(V); xp = np.empty_like(V)
                xm[1:, :] = V[:-1, :]; xm[0, :] = V[0, :]
                xp[:-1, :] = V[1:, :]; xp[-1, :] = V[1:, :][-1:]
                avg = np.zeros_like(V)
                avg[:, 1:-1] = 0.25 * (xm[:, 1:-1] + xp[:, 1:-1] + V[:, 2:] + V[:, :-2])
                m = inside & color
                V[m] = (1.0 - omega) * V[m] + omega * avg[m]
        return apply_bc(V)

    def poisson(V, sweeps=200, omega=1.7):
        """Variable-eps Poisson div(eps grad phi) = -rho*k, red-black SOR. Conductors are
        equipotential Dirichlet (set outside); top/bottom grounded; gas+insulator solved."""
        k = float(rho_coupling)
        for _ in range(sweeps):
            V[:, 0] = 0.0
            if ncomp > 0:
                V[is_cond] = Vc[cid[is_cond]]
            xm = np.empty_like(V); xp = np.empty_like(V)
            zm = np.empty_like(V); zp = np.empty_like(V)
            for color in (red, ~red):
                xm[1:, :] = V[:-1, :]; xm[0, :] = V[0, :]
                xp[:-1, :] = V[1:, :]; xp[-1, :] = V[-1, :]
                zm[:, 1:] = V[:, :-1]; zm[:, 0] = V[:, 0]
                zp[:, :-1] = V[:, 1:]; zp[:, -1] = V[:, -1]
                num = e_xm * xm + e_xp * xp + e_zm * zm + e_zp * zp + k * rho
                upd = num / e_diag
                m = poisson_inside & color
                V[m] = (1.0 - omega) * V[m] + omega * upd[m]
        V[:, 0] = 0.0
        if ncomp > 0:
            V[is_cond] = Vc[cid[is_cond]]
        return V

    def trace(kind, n, Ex, Ez, want_energy=False):
        if kind == "ion":
            x, z, vx, vz = sample_ions(n, rng, mouth, nx, V_dc, V_rf, iadf_hwhm_deg)
            q = 1.0
        else:
            x, z, vx, vz = sample_electrons(n, rng, nx, Te)
            q = -1.0
        hix, hiz, E, _, surv = _trace_general(Ex, Ez, solid, x, z, vx, vz, q, nx, nz, 40 * nz)
        counts = np.zeros((nx, nz))
        energy = np.zeros((nx, nz)) if want_energy else None
        m = hix >= 0
        if m.any():
            np.add.at(counts, (hix[m], hiz[m]), 1.0)
            if want_energy:
                np.add.at(energy, (hix[m], hiz[m]), E[m])
        if want_energy:
            return counts, float(surv.mean()), energy
        return counts, float(surv.mean())

    use_poisson = field_model == "poisson"
    V = np.zeros((nx, nz))
    hist = []
    frames = []          # (iter, V snapshot, Vc snapshot) for the dynamics movie
    vc_tail_sum = np.zeros(ncomp); vs_tail_sum = np.zeros((nx, nz)); tail_n = 0
    for it in range(n_iter):
        V = poisson(V) if use_poisson else laplace(V)
        Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
        ci, si = trace("ion", n_per_iter, Ex, Ez)
        Vcell = (Vs if not use_poisson else V).copy()
        if ncomp > 0:
            Vcell[is_cond] = Vc[cid[is_cond]]
        thr = np.where(Vcell >= 0.0, 1.0, np.exp(np.clip(Vcell / max(Te, 1e-9), -40.0, 0.0)))
        if electron_model == "vf" and vf_grid is not None:
            # view-factor electrons x POTENTIAL-DEPENDENT electrostatic focusing. Electrons are pulled
            # harder toward the MOST-POSITIVE surface (the floor), so the floor's collection is enhanced
            # beyond geometry (HG: "electrostatics decreases the geometric shadowing") while the low-V
            # sidewalls stay starved. focus = vf_focus + vf_focus_pot*max(V,0). This is the physical
            # form of the focusing the constant vf_focus faked; it lets the floor reach 0.22 (-> floor V
            # ~33) without over-feeding the walls, so the neighbour can rise ABOVE the floor to ~39.
            se = 0.0
            focus = float(vf_focus) + float(vf_focus_pot) * np.maximum(Vcell, 0.0)
            ce = (float(n_per_iter) / nx) * vf_grid * thr * focus
        else:
            ce, se = trace("electron", n_per_iter, Ex, Ez)
            if vf_grid is not None:
                ce = np.maximum(ce, (float(n_per_iter) / nx) * vf_grid * thr)
        net = ci - ce
        # Robbins-Monro decaying step (NO floor) so the stochastic relaxation reaches a true fixed
        # point instead of drifting; tail-average the potentials (Polyak) for the steady-state value.
        anneal = 1.0 / (1.0 + it / (0.15 * n_iter))
        scale = anneal * relax / n_per_iter * nx
        if use_poisson:
            # insulators accumulate free CHARGE; the Poisson solve maps it to potential via eps
            rho[insul] += scale * net[insul]
        else:
            Vs[insul] += scale * net[insul]
            Vs[insul] = np.clip(Vs[insul], -insul_vguard, V_dc + V_rf)
        # CRYO MODULE (frontier, nobody has it open-source): a temperature-gated condensed HF/H2O
        # layer raises insulator SURFACE conductivity 3-6 orders below ~-60C (APL 123,212106 2023),
        # letting accumulated feature charge FLOW LATERALLY along surfaces and dissipate toward the
        # grounded regions -- un-bending ion trajectories and relieving deep-AR over-charge. Modeled
        # as a T-gated lateral diffusion of the insulator potential. sigma(T): ~0 warm, ~1 cryo.
        sigma = _cryo_conductivity(temperature_C) if surface_conductivity == "auto" else float(surface_conductivity)
        if sigma > 0.0:
            from scipy.ndimage import uniform_filter
            Vloc = uniform_filter(Vs, size=3, mode="nearest")
            Vs[insul] += sigma * (Vloc[insul] - Vs[insul])   # charge spreads along the surface
            Vs[insul] *= (1.0 - 0.15 * sigma)                # + leaks toward ground (net dissipation)
        for c in range(1, ncomp + 1):
            m = cid == c
            area = max(int(m.sum()), 1)
            Vc[c] += scale * float(net[m].sum()) / area
            Vc[c] = float(np.clip(Vc[c], -3.0 * Te, V_dc + V_rf))
        hist.append((float(si), float(se)))
        if it >= int(0.6 * n_iter):        # accumulate tail for Polyak steady-state average
            vc_tail_sum += Vc[1:]; vs_tail_sum += Vs; tail_n += 1
        if frame_every and (it % frame_every == 0 or it == n_iter - 1):
            frames.append((it, V.copy(), Vc[1:].copy()))
        if verbose and it % 20 == 0:
            vmax = V.max() if use_poisson else Vs.max()
            print(f"  it{it}: surv_i/e={si:.3f}/{se:.3f} "
                  f"Vc={np.round(Vc[1:], 1) if ncomp else '-'} "
                  f"Vmax={vmax:.1f}", flush=True)

    # Polyak steady-state average over the tail (kills residual MC noise + any slow drift)
    if tail_n > 0:
        if ncomp > 0:
            Vc[1:] = vc_tail_sum / tail_n
        Vs = vs_tail_sum / tail_n
    V = poisson(V, sweeps=180) if use_poisson else laplace(V, sweeps=180)
    Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
    # final high-stat pass for diagnostics (ion counts + impact energy grid, electron collection)
    ci_f, _, Ei_f = trace("ion", 4 * n_per_iter, Ex, Ez, want_energy=True)
    ntot = 4 * n_per_iter
    Vcell = (Vs if not use_poisson else V).copy()
    if ncomp > 0:
        Vcell[is_cond] = Vc[cid[is_cond]]
    thr = np.where(Vcell >= 0.0, 1.0, np.exp(np.clip(Vcell / max(Te, 1e-9), -40.0, 0.0)))
    if electron_model == "vf" and vf_grid is not None:
        focus = float(vf_focus) + float(vf_focus_pot) * np.maximum(Vcell, 0.0)
        ce_f = (float(ntot) / nx) * vf_grid * thr * focus
    else:
        ce_f, _ = trace("electron", 4 * n_per_iter, Ex, Ez)
        if vf_grid is not None:
            ce_f = np.maximum(ce_f, (float(ntot) / nx) * vf_grid * thr)
    # per-cell insulator potential is the solved field at insulator cells (poisson) or Vs (laplace)
    Vs_out = np.where(insul, V, 0.0) if use_poisson else Vs
    return dict(V=V, Vs=Vs_out, Vc=Vc[1:], ncomp=ncomp, cid=cid, rho=rho,
                surv_ion=hist[-1][0], surv_electron=hist[-1][1],
                field_model=field_model,
                ion_counts=ci_f, ion_energy=Ei_f, electron_counts=ce_f, ntot=ntot,
                solid=solid, insul=insul, frames=frames)
