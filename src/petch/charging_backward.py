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


def backward_electron_gather(solid, Ex, Ez, V_surf, cells, normals, Te=4.0,
                             n_log2=13, n_scramble=3, trace_dt=0.15, trace_dt_field=0.10,
                             trace_steps=200, seed=0):
    """Backward electron flux per cell (fraction of incident electron flux; open flat V=0 -> 1).

    Incident-energy sampling: E_top ~ Gamma(2,Te) flux-Maxwellian, Lambert cos-flux angle about
    vertical, transverse conserved. E_surf = E_top + Vc (electron); needs E_top*ct^2 > |Vc| to sit on a
    negative wall -> retardation automatic; on a positive floor E_surf>0 always -> saturated, no blowup.
    Flux factor (v.n_cell)/(v.z_surf): =1 on the floor (preserves the Langmuir law), <<1 on walls."""
    nx, nz = solid.shape
    msteps = int(trace_steps) * nz
    N = 2 ** n_log2
    out = np.zeros(len(cells))
    for ci, ((cx, cz), (nnx, nnz)) in enumerate(zip(cells, normals)):
        Vc = float(V_surf[cx, cz])
        vals = []
        for sc in range(n_scramble):
            s = qmc.Sobol(d=3, scramble=True, seed=seed + sc)
            u = s.random_base2(n_log2)
            E_top = gammadist.ppf(u[:, 0], a=2.0, scale=Te)
            ct = np.sqrt(u[:, 1]); st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
            sgn = np.where(u[:, 2] < 0.5, 1.0, -1.0)
            spd = np.sqrt(E_top)
            vperp = spd * st * sgn
            E_surf = E_top + Vc
            vz2 = E_surf - vperp * vperp
            valid = (E_surf > 0.0) & (vz2 > 0.0)
            vz_surf = np.sqrt(np.maximum(vz2, 0.0))
            vX = vperp; vZ = -vz_surf
            vdotn = vX * nnx + vZ * nnz
            emit = valid & (vdotn > 0.0)
            x0 = np.full(N, cx + 1.5 * nnx); z0 = np.full(N, cz + 1.5 * nnz)
            hix, hiz, _, _, surv, _, _ = _trace_general(Ex, Ez, solid, x0, z0, vX, vZ, -1.0, nx, nz,
                                                        msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5) & emit
            ratio = vdotn / np.maximum(np.abs(vZ), 0.3)
            vals.append(float((escaped * ratio).sum()) / N)
        out[ci] = float(np.mean(vals))
    return out


def backward_ion_gather(solid, Ex, Ez, V_surf, cells, normals, Te=4.0, Ti=0.5, V_dc=37.0, V_rf=30.0,
                        n_log2=13, n_scramble=3, trace_dt=0.15, trace_dt_field=0.10, trace_steps=200,
                        seed=0):
    """Backward ion flux per cell (fraction of incident ion flux; open flat V=0 -> 1).

    Incident ion (matches sample_sheath_source): phase -> Vs=V_dc+V_rf*sin (arcsine IED, weight
    Vs^-0.35), vz_in=sqrt(0.5 Te+Vs), transverse vperp~N(0,sqrt(0.5 Ti)). Near-VERTICAL in the lab
    frame, so v.n>0 on a wall selects only the grazing tail. E_surf_z = 0.5 Te+Vs-Vc; E_surf_z<0 =>
    reflected (floor repels the low-IED-horn = retardation). Flux factor (v.n_cell)/(v.z_surf) suppresses
    grazing ions on vertical walls."""
    nx, nz = solid.shape
    msteps = int(trace_steps) * nz
    sig = np.sqrt(0.5 * Ti)
    out = np.zeros(len(cells))
    for ci, ((cx, cz), (nnx, nnz)) in enumerate(zip(cells, normals)):
        Vc = float(V_surf[cx, cz])
        vals = []
        for sc in range(n_scramble):
            s = qmc.Sobol(d=3, scramble=True, seed=seed + sc)
            u = s.random_base2(n_log2)
            ph = u[:, 0] * 2.0 * np.pi
            Vs = V_dc + V_rf * np.sin(ph)
            wied = Vs ** (-0.35)
            vperp = sig * norm.ppf(np.clip(u[:, 2], 1e-6, 1 - 1e-6))
            E_surf_z = 0.5 * Te + Vs - Vc
            valid = E_surf_z > 0.0
            vz_surf = np.sqrt(np.maximum(E_surf_z, 0.0))
            vX = vperp.copy(); vZ = -vz_surf
            vdotn = vX * nnx + vZ * nnz
            emit = valid & (vdotn > 0.0)
            x0 = np.full(u.shape[0], cx + 1.5 * nnx); z0 = np.full(u.shape[0], cz + 1.5 * nnz)
            hix, hiz, _, _, surv, _, _ = _trace_general(Ex, Ez, solid, x0, z0, vX, vZ, 1.0, nx, nz,
                                                        msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5) & emit
            ratio = vdotn / np.maximum(np.abs(vZ), 0.3)
            vals.append(float((wied * escaped * ratio).sum() / wied.sum()))
        out[ci] = float(np.mean(vals))
    return out


def self_consistent_backward(g, Te=4.0, n_iter=14, beta=0.5, dVmax=8.0, n_log2=11, n_scramble=2,
                             n_wall=12, n_floor=6, sweeps=250, seed=0):
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
    cells, normals, kind = [], [], []
    wrows = np.linspace(r0 + 1, r1 - 2, n_wall).astype(int)
    for z in wrows:
        cells.append((t0 - 1, int(z))); normals.append((1.0, 0.0)); kind.append('Lwall')
        cells.append((t1, int(z))); normals.append((-1.0, 0.0)); kind.append('Rwall')
    fcols = np.linspace(t0 + 2, t1 - 3, n_floor).astype(int)
    for x in fcols:
        cells.append((int(x), fz)); normals.append((0.0, -1.0)); kind.append('floor')
    for x in (g['edge0'] + 2, g['neigh1'] - 3):
        cells.append((int(x), r0)); normals.append((0.0, -1.0)); kind.append('mask')
    kind = np.array(kind)
    clist = [tuple(c) for c in cells]; nlist = [tuple(n) for n in normals]
    Vs = np.zeros((nx, nz))

    def laplace(Vsurf, omega=1.9):
        V = np.zeros((nx, nz)); ins = solid & (cond == 0); condm = solid & (cond > 0)
        def bc(V):
            V[:, 0] = 0.0; V[condm] = 0.0; V[ins] = Vsurf[ins]; return V
        for _ in range(sweeps):
            V = bc(V)
            for color in (red, ~red):
                xm = np.empty_like(V); xp = np.empty_like(V); zm = np.empty_like(V); zp = np.empty_like(V)
                xm[1:] = V[:-1]; xm[0] = V[0]; xp[:-1] = V[1:]; xp[-1] = V[-1]
                zm[:, 1:] = V[:, :-1]; zm[:, 0] = V[:, 0]; zp[:, :-1] = V[:, 1:]; zp[:, -1] = V[:, -1]
                upd = 0.25 * (xm + xp + zm + zp); m = gas & color
                V[m] = (1 - omega) * V[m] + omega * upd[m]
        return bc(V)

    for it in range(n_iter):
        V = laplace(Vs); Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
        Ge = backward_electron_gather(solid, Ex, Ez, Vs, clist, nlist, Te=Te, n_log2=n_log2,
                                      n_scramble=n_scramble, seed=seed)
        Gi = backward_ion_gather(solid, Ex, Ez, Vs, clist, nlist, Te=Te, n_log2=n_log2,
                                 n_scramble=n_scramble, seed=seed)
        m = kind == 'mask'
        k = float(np.mean(Ge[m]) / max(np.mean(Gi[m]), 1e-6))
        dV = np.clip(beta * Te * np.log((k * Gi + 1e-6) / (Ge + 1e-6)), -dVmax, dVmax)
        for (cx, cz), d in zip(clist, dV):
            Vs[cx, cz] += d
    return dict(Vs=Vs, wall_rows=wrows, wall_depth=wrows - r0, Lwall=Vs[t0 - 1, wrows],
                Rwall=Vs[t1, wrows], floor=Vs[fcols, fz], floor_mean=float(Vs[fcols, fz].mean()), k=k)
