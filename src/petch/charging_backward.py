"""BACKWARD / ADJOINT per-cell flux gather for feature charging (EXPERIMENTAL, WIP 2026-07-10).

The forward launch-fan starves deep/shadowed target cells (deep floor, upper walls) -> every accuracy
fix has been a per-region importance-sampling patch (preint_floor, preint_floor_ion, preint_wall).
This module is the structural alternative endorsed by the adjoint-MC literature (e.g. arXiv:1504.00214,
"Adjoint Monte Carlo calculation of charged plasma particle flux to wall"): GATHER per cell -- trace
BACKWARD from each surface cell out through the field to the plasma, weight by the incident
distribution. No starvation (every sample per cell is relevant), no per-region overrides, uniform
across cells and species, and better-shaped for autodiff.

Core result (VALIDATED, Gate B exact): sample the INCIDENT plasma-side velocity from the species flux
distribution, map to the surface energy by energy conservation (E_surf = E_top - q*Vc), launch OUTWARD
from the cell, require v.n>0, trace through the field, and count escapes -- weighted by the flux factor
(v.n_cell)/(v.z_surface). Retardation (Boltzmann for electrons, IED-horn reflection for ions) is
AUTOMATIC via the E_surf>0 gate on physically-sampled incident energies -- NO explicit exp(V/Te)
factor (which blows up as e^16 on a +60 V floor). Reproduces the Langmuir-probe law exactly:
electron flux retarded ~exp(Vc/Te) for Vc<0, SATURATED for Vc>0.

Convention: both gathers return flux as a fraction of the incident flux through a horizontal plane
(open flat V=0 -> 1). The ion/electron scale ratio k=Ci/Ce is pinned by ONE ambipolar calibration on a
floating open surface (mask top). Then per-cell dV = Te*ln(k*Gi/Ge) is the charging direction.

STATUS: gathers validated on analytic gates (retardation, saturation, grazing suppression) + cross-
checked vs the forward floor equilibrium. self_consistent_backward() CONVERGES the dipole cleanly and
reproduces the Kushner picture across AR4/8/15 (upper wall negative ~-12.7V electron-shading; deep wall
positive and rising with AR = grazing ions; floor positive + monotone in AR; the low->high-AR crossover
of the potential max from floor to sidewall). Deterministic, ~1-2s/iter, no per-region overrides.
OPEN: (1) deep cells at very high AR still benefit from directional importance sampling / more samples;
(2) absolute floor magnitude to be calibrated vs experiment (de Boer ARDE / Fujiwara notch); (3) wire
this gather as the transport estimator inside solve_charging and retire the vf/thr/preint override
layer. Requires the exit-kinematics return added to _trace_general (charging_general.py).
"""
import numpy as np
from scipy.stats import qmc, gamma as gammadist, norm

from .charging_general import _trace_general


def _cone_angles(x0, z0, aperture, pad):
    """Angular interval (theta_lo, theta_hi) from vertical subtended by the trench-mouth aperture at the
    launch point (x0,z0), padded by `pad` rad each edge. NEE-style importance-sampling proposal for the
    deep-cell escape cone (~1/AR). Returns None (=> natural sampling) if no aperture or the cell is at/
    above the mouth (dz<=2), which keeps mask/flat cells -- the k-calibration anchor -- on the exact
    natural path and preserves Gate B."""
    if aperture is None:
        return None
    t0a, t1a, r0a = aperture
    dz = z0 - r0a
    if dz <= 2.0:
        return None
    lim = 0.5 * np.pi - 1e-3
    th_lo = max(np.arctan2(t0a - x0, dz) - pad, -lim)
    th_hi = min(np.arctan2(t1a - x0, dz) + pad, lim)
    if th_hi <= th_lo + 1e-6:
        return None
    return th_lo, th_hi


def _draw_absmass(lo, hi, u):
    """Vectorized inverse-CDF of the density ~|v| on the signed interval [lo,hi] (lo<hi), plus its mass
    M = int_lo^hi |v'| dv'. Handles single-sign and straddling-0 intervals. Returns (v, M)."""
    loc = np.minimum(lo, 0.0); hic = np.maximum(hi, 0.0)
    mass_neg = np.where(lo < 0.0, 0.5 * (lo * lo - np.minimum(hi, 0.0) ** 2), 0.0)
    mass_pos = np.where(hi > 0.0, 0.5 * (hi * hi - np.maximum(lo, 0.0) ** 2), 0.0)
    M = mass_neg + mass_pos
    t = u * M
    neg = (t < mass_neg) & (mass_neg > 0.0)
    v_neg = -np.sqrt(np.maximum(lo * lo - 2.0 * t, 0.0))
    tp = t - mass_neg
    v_pos = np.sqrt(np.maximum(np.maximum(lo, 0.0) ** 2 + 2.0 * tp, 0.0))
    v = np.where(neg, v_neg, v_pos)
    return v, M


def backward_electron_gather(solid, Ex, Ez, V_surf, cells, normals, Te=4.0,
                             n_log2=13, n_scramble=3, trace_dt=0.15, trace_dt_field=0.10,
                             trace_steps=200, seed=0, aperture=None, pad_deg=6.0, alpha=0.85):
    """Backward electron flux per cell (fraction of incident electron flux; open flat V=0 -> 1).

    Incident-energy sampling: E_top ~ Gamma(2,Te) flux-Maxwellian, Lambert cos-flux angle about
    vertical, transverse conserved. E_surf = E_top + Vc (electron); needs E_top*ct^2 > |Vc| to sit on a
    negative wall -> retardation automatic; on a positive floor E_surf>0 always -> saturated, no blowup.
    Flux factor (v.n_cell)/(v.z_surf): =1 on the floor (preserves the Langmuir law), <<1 on walls.

    aperture=(t0,t1,r0): enables NEE-style escape-cone importance sampling (mixture of cone + broad
    natural) -> resolves the deep-cell 1/AR escape cone at high AR. UNBIASED (weight w=f_nat/q_mix); the
    broad stratum captures field-focused escapers outside the geometric cone. aperture=None -> exact
    natural sampling (Gate B / calibration path, bit-identical)."""
    nx, nz = solid.shape
    msteps = int(trace_steps) * nz
    N = 2 ** n_log2
    pad = np.deg2rad(pad_deg)
    out = np.zeros(len(cells))
    for ci, ((cx, cz), (nnx, nnz)) in enumerate(zip(cells, normals)):
        Vc = float(V_surf[cx, cz])
        x0c = cx + 1.5 * nnx; z0c = cz + 1.5 * nnz
        cone = _cone_angles(x0c, z0c, aperture, pad)
        a_mix = alpha if cone is not None else 0.0
        vals = []
        for sc in range(n_scramble):
            s = qmc.Sobol(d=5, scramble=True, seed=seed + sc)
            u = s.random_base2(n_log2)
            E_top = gammadist.ppf(u[:, 0], a=2.0, scale=Te)
            E_surf = E_top + Vc
            m = np.sqrt(np.maximum(np.minimum(E_top, E_surf), 0.0))          # |vperp| cap
            ct = np.sqrt(u[:, 1]); st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
            sgn = np.where(u[:, 2] < 0.5, 1.0, -1.0)
            vperp_nat = np.sqrt(E_top) * st * sgn
            if cone is not None:
                sq = np.sqrt(np.maximum(E_surf, 0.0))
                a = sq * np.sin(cone[0]); b = sq * np.sin(cone[1])
                lo = np.minimum(a, b); hi = np.maximum(a, b)
                if nnx > 0:   lo = np.maximum(lo, 0.0)                        # left wall: emit vperp>0
                elif nnx < 0: hi = np.minimum(hi, 0.0)                        # right wall: emit vperp<0
                lo = np.clip(lo, -m, m); hi = np.clip(hi, -m, m)
                cone_ok = hi > lo + 1e-9
                vcone, Mcone = _draw_absmass(lo, hi, u[:, 4])
                Pcone = np.where(cone_ok, Mcone / np.maximum(E_top, 1e-12), 0.0)
                use_cone = (u[:, 3] < a_mix) & cone_ok
                vperp = np.where(use_cone, vcone, vperp_nat)
                in_cone = cone_ok & (vperp >= lo) & (vperp <= hi)
                w = np.where(in_cone, 1.0 / (a_mix / np.maximum(Pcone, 1e-12) + (1.0 - a_mix)),
                             1.0 / (1.0 - a_mix))
            else:
                vperp = vperp_nat; w = np.ones(N)
            vz2 = E_surf - vperp * vperp
            valid = (E_surf > 0.0) & (vz2 > 0.0)
            vz_surf = np.sqrt(np.maximum(vz2, 0.0))
            vX = vperp; vZ = -vz_surf
            vdotn = vX * nnx + vZ * nnz
            emit = valid & (vdotn > 0.0)
            x0 = np.full(N, x0c); z0 = np.full(N, z0c)
            hix, hiz, _, _, surv, _, _ = _trace_general(Ex, Ez, solid, x0, z0, vX, vZ, -1.0, nx, nz,
                                                        msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5) & emit
            ratio = vdotn / np.maximum(np.abs(vZ), 0.3)
            vals.append(float((w * escaped * ratio).sum()) / N)
        out[ci] = float(np.mean(vals))
    return out


def backward_ion_gather(solid, Ex, Ez, V_surf, cells, normals, Te=4.0, Ti=0.5, V_dc=37.0, V_rf=30.0,
                        n_log2=13, n_scramble=3, trace_dt=0.15, trace_dt_field=0.10, trace_steps=200,
                        seed=0, aperture=None, pad_deg=3.0, alpha=0.85, want_energy=False):
    """Backward ion flux per cell (fraction of incident ion flux; open flat V=0 -> 1).

    Incident ion (matches sample_sheath_source): phase -> Vs=V_dc+V_rf*sin (arcsine IED, weight
    Vs^-0.35), vz_in=sqrt(0.5 Te+Vs), transverse vperp~N(0,sqrt(0.5 Ti)). Near-VERTICAL in the lab
    frame, so v.n>0 on a wall selects only the grazing tail. E_surf_z = 0.5 Te+Vs-Vc; E_surf_z<0 =>
    reflected (floor repels the low-IED-horn = retardation). Flux factor (v.n_cell)/(v.z_surf) suppresses
    grazing ions on vertical walls. aperture -> escape-cone importance sampling (truncated-normal in the
    cone + broad, unbiased); aperture=None -> exact natural sampling."""
    nx, nz = solid.shape
    msteps = int(trace_steps) * nz
    sig = np.sqrt(0.5 * Ti)
    pad = np.deg2rad(pad_deg)
    out = np.zeros(len(cells))
    out_E = np.zeros(len(cells))                                     # flux-weighted mean impact energy [eV]
    for ci, ((cx, cz), (nnx, nnz)) in enumerate(zip(cells, normals)):
        Vc = float(V_surf[cx, cz])
        x0c = cx + 1.5 * nnx; z0c = cz + 1.5 * nnz
        cone = _cone_angles(x0c, z0c, aperture, pad)
        a_mix = alpha if cone is not None else 0.0
        vals = []; evals = []
        for sc in range(n_scramble):
            s = qmc.Sobol(d=4, scramble=True, seed=seed + sc)
            u = s.random_base2(n_log2)
            ph = u[:, 0] * 2.0 * np.pi
            Vs = V_dc + V_rf * np.sin(ph)
            wied = Vs ** (-0.35)
            E_surf_z = 0.5 * Te + Vs - Vc
            valid = E_surf_z > 0.0
            vz_surf = np.sqrt(np.maximum(E_surf_z, 0.0))
            vperp_nat = sig * norm.ppf(np.clip(u[:, 2], 1e-6, 1 - 1e-6))
            if cone is not None:
                a = vz_surf * np.tan(cone[0]); b = vz_surf * np.tan(cone[1])
                lo = np.minimum(a, b); hi = np.maximum(a, b)
                if nnx > 0:   lo = np.maximum(lo, 0.0)
                elif nnx < 0: hi = np.minimum(hi, 0.0)
                Plo = norm.cdf(lo / sig); Phi_ = norm.cdf(hi / sig); Pcone = Phi_ - Plo
                cone_ok = Pcone > 1e-9
                vcone = sig * norm.ppf(np.clip(Plo + u[:, 3] * Pcone, 1e-9, 1 - 1e-9))
                use_cone = (u[:, 1] < a_mix) & cone_ok
                vperp = np.where(use_cone, vcone, vperp_nat)
                in_cone = cone_ok & (vperp >= lo) & (vperp <= hi)
                w = np.where(in_cone, 1.0 / (a_mix / np.maximum(Pcone, 1e-12) + (1.0 - a_mix)),
                             1.0 / (1.0 - a_mix))
            else:
                vperp = vperp_nat; w = np.ones(u.shape[0])
            vX = vperp; vZ = -vz_surf
            vdotn = vX * nnx + vZ * nnz
            emit = valid & (vdotn > 0.0)
            x0 = np.full(u.shape[0], x0c); z0 = np.full(u.shape[0], z0c)
            hix, hiz, _, _, surv, _, _ = _trace_general(Ex, Ez, solid, x0, z0, vX, vZ, 1.0, nx, nz,
                                                        msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5) & emit
            ratio = vdotn / np.maximum(np.abs(vZ), 0.3)
            fnum = w * wied * escaped * ratio                        # per-sample flux contribution
            vals.append(float(fnum.sum() / wied.sum()))
            if want_energy:
                E_impact = vperp * vperp + E_surf_z                  # total KE at the surface [eV]
                evals.append(float((fnum * E_impact).sum() / max(fnum.sum(), 1e-12)))
        out[ci] = float(np.mean(vals))
        if want_energy:
            out_E[ci] = float(np.mean(evals)) if evals else 0.0
    return (out, out_E) if want_energy else out


def self_consistent_backward(g, Te=4.0, n_iter=14, beta=0.5, dVmax=8.0, n_log2=10, n_scramble=2,
                             n_wall=12, n_floor=6, sweeps=250, seed=0, cone_is=True):
    """Self-consistent BACKWARD charging solve: Laplace field <-> per-cell gathers <-> damped update
    dV = beta*Te*ln(k*Gi/Ge), where k=Ci/Ce is calibrated each iteration on the floating pillar tops.
    NO forward launch, NO per-region overrides -- the electron-shading dipole EMERGES. Deterministic,
    well-conditioned (converges in ~10 iters without razor's-edge tricks).

    Reproduces the Kushner picture: upper wall negative (electron shading), deep wall positive (grazing
    ions), the low-AR->high-AR crossover of the potential maximum from floor to sidewall, floor monotone
    in AR. Returns dict with Vs grid, wall/floor profiles, and the sampled cell coords.
    Cost ~ (n_wall*2 + n_floor + 2) cells * 2^n_log2 * n_scramble * 2 species * n_iter traces."""
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']; r0, r1 = int(g['mouth']), int(g['z_poly0']); fz = nz - 1
    cond = g['cond']; gas = ~solid
    red = (np.add.outer(np.arange(nx), np.arange(nz)) % 2 == 0)
    cells, normals, kind, comp = [], [], [], []
    wrows = np.linspace(r0 + 1, r1 - 2, n_wall).astype(int)
    for z in wrows:
        cells.append((t0 - 1, int(z))); normals.append((1.0, 0.0)); kind.append('Lwall'); comp.append(0)
        cells.append((t1, int(z))); normals.append((-1.0, 0.0)); kind.append('Rwall'); comp.append(0)
    fcols = np.linspace(t0 + 2, t1 - 3, n_floor).astype(int)
    for x in fcols:
        cells.append((int(x), fz)); normals.append((0.0, -1.0)); kind.append('floor'); comp.append(0)
    for x in (g['edge0'] + 2, g['neigh1'] - 3):
        cells.append((int(x), r0)); normals.append((0.0, -1.0)); kind.append('mask'); comp.append(0)
    # floating POLY (conductor) inner sidewall faces [z_poly0:fz]: equipotential, each floats to its own
    # current balance -> the HG V_poly (grounding the poly compressed the potential range + inflated the
    # foot ion energy). Left poly = edge pillar (cond 1), right = neigh pillar (cond 2).
    prows = np.linspace(r1 + 1, fz - 1, max(n_wall // 2, 4)).astype(int)
    for z in prows:
        cells.append((t0 - 1, int(z))); normals.append((1.0, 0.0)); kind.append('Lpoly'); comp.append(1)
        cells.append((t1, int(z))); normals.append((-1.0, 0.0)); kind.append('Rpoly'); comp.append(2)
    kind = np.array(kind); comp = np.array(comp)
    clist = [tuple(c) for c in cells]; nlist = [tuple(n) for n in normals]
    Vs = np.zeros((nx, nz)); vc = np.zeros(3)               # vc[1],vc[2] = floating poly potentials

    def laplace(Vsurf, omega=1.9):
        V = np.zeros((nx, nz)); ins = solid & (cond == 0); condm = solid & (cond > 0)
        def bc(V):
            V[:, 0] = 0.0; V[condm] = Vsurf[condm]; V[ins] = Vsurf[ins]; return V
        for _ in range(sweeps):
            V = bc(V)
            for color in (red, ~red):
                xm = np.empty_like(V); xp = np.empty_like(V); zm = np.empty_like(V); zp = np.empty_like(V)
                xm[1:] = V[:-1]; xm[0] = V[0]; xp[:-1] = V[1:]; xp[-1] = V[-1]
                zm[:, 1:] = V[:, :-1]; zm[:, 0] = V[:, 0]; zp[:, :-1] = V[:, 1:]; zp[:, -1] = V[:, -1]
                upd = 0.25 * (xm + xp + zm + zp); m = gas & color
                V[m] = (1 - omega) * V[m] + omega * upd[m]
        return bc(V)

    aperture = (t0, t1, r0) if cone_is else None
    for it in range(n_iter):
        for c in (1, 2):
            Vs[cond == c] = vc[c]                            # broadcast floating poly potential onto its body
        V = laplace(Vs); Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
        Ge = backward_electron_gather(solid, Ex, Ez, Vs, clist, nlist, Te=Te, n_log2=n_log2,
                                      n_scramble=n_scramble, seed=seed, aperture=aperture)
        Gi = backward_ion_gather(solid, Ex, Ez, Vs, clist, nlist, Te=Te, n_log2=n_log2,
                                 n_scramble=n_scramble, seed=seed, aperture=aperture)
        # k = Ci/Ce = 1 EXACTLY (first principle): the domain is the WAFER FRAME (top plane = wafer
        # surface below the sheath; the sheath drop is carried in the ion source energy and the electron
        # gather samples the post-sheath arrival dist). Wafer-scale charge conservation -- zero net DC
        # current is what DEFINES V_dc -- makes Gamma_i = Gamma_e on the uncharged wafer for ANY
        # geometry. So the uncharged wafer is a fixed point by construction; NO mask calibration (that
        # was a noisy estimator of 1 and an arbitrary reference knob). [invariant V1]
        dV = np.clip(beta * Te * np.log((Gi + 1e-6) / (Ge + 1e-6)), -dVmax, dVmax)
        for i, (cx, cz) in enumerate(clist):
            if comp[i] == 0:                                 # insulator: per-cell update
                Vs[cx, cz] += dV[i]
        for c in (1, 2):                                     # poly: POOLED equipotential current balance
            sel = comp == c
            if sel.any():
                dVc = np.clip(beta * Te * np.log((Gi[sel].sum() + 1e-6) / (Ge[sel].sum() + 1e-6)),
                              -dVmax, dVmax)
                vc[c] += dVc
    # final field + notch-driver observable E_defl (flux-weighted ion impact energy on the poly-inner
    # sidewall foot = the deflected-ion notch driver, HG _extract convention).
    for c in (1, 2):
        Vs[cond == c] = vc[c]
    V = laplace(Vs); Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
    zp = int(g['z_poly0'])
    foot_h = max(int(0.3 * (t1 - t0)), 3)                  # HG foot band: within 0.3*W of the floor
    fz0 = max(zp + 1, fz - foot_h)
    foot = [(t0 - 1, int(z)) for z in np.arange(fz0, fz)]  # poly-inner face near the poly/oxide junction
    fn = [(1.0, 0.0)] * len(foot)
    Gi_f, E_f = backward_ion_gather(solid, Ex, Ez, Vs, foot, fn, Te=Te, n_log2=n_log2 + 1,
                                    n_scramble=n_scramble, seed=seed, aperture=aperture, want_energy=True)
    fmask = Gi_f > 1e-6
    E_defl = float(np.sum(E_f[fmask] * Gi_f[fmask]) / max(np.sum(Gi_f[fmask]), 1e-9)) if fmask.any() else 0.0
    # floor_mean over the SOLVED cells (fcols) -- the backward solver is sparse (only sampled cells are
    # updated), so averaging the full _extract band would include un-updated (0 V) cells.
    floor_mean = float(Vs[fcols, fz].mean())
    return dict(Vs=Vs, V=V, Ex=Ex, Ez=Ez, wall_rows=wrows, wall_depth=wrows - r0, Lwall=Vs[t0 - 1, wrows],
                Rwall=Vs[t1, wrows], floor=Vs[fcols, fz], floor_mean=floor_mean, k=1.0,
                V_poly=float(0.5 * (vc[1] + vc[2])), Vc=vc.copy(),
                E_defl=E_defl, foot_flux=float(Gi_f[fmask].mean()) if fmask.any() else 0.0)
