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
from .charging_nodal import solve_nodal_laplace, trace_nodal
from .boundary_transport import (
    adaptive_adjoint_boundary_state_face_flux,
    adjoint_boundary_state_face_flux,
    bidirectional_boundary_state_cell_flux,
)
from .adaptive_quadrature import adaptive_surface_quadrature


class AdaptiveQuadratureConvergenceError(RuntimeError):
    """Carries the rejected fixed-point state for estimator diagnosis without accepting it."""

    def __init__(self, message, *, iteration, species, quadrature, surface_voltage, potential,
                 cells, normals):
        super().__init__(message)
        self.iteration = int(iteration)
        self.species = species
        self.quadrature = quadrature
        self.surface_voltage = np.asarray(surface_voltage).copy()
        self.potential = np.asarray(potential).copy()
        self.cells = tuple(cells)
        self.normals = tuple(normals)


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


def _gas_faces(solid, target_mask, subsample=None):
    """Enumerate every gas-facing face of the target_mask solid cells: returns (cells, normals) where the
    outward normal points into the adjacent gas. A corner cell contributes one face per exposed side.
    General (any geometry); used to pool a floating conductor's current over its FULL surface. Optional
    even subsample keeps the pooled-current estimate cheap on long faces."""
    nx, nz = solid.shape; gas = ~solid
    cells, normals = [], []
    for dx, dz in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        sh = np.zeros_like(gas)
        i0 = max(dx, 0); i1 = nx + min(dx, 0); j0 = max(dz, 0); j1 = nz + min(dz, 0)
        sh[i0 - dx:i1 - dx, j0 - dz:j1 - dz] = gas[i0:i1, j0:j1]        # gas shifted opposite the normal
        face = target_mask & sh
        ii, jj = np.where(face)
        for i, j in zip(ii, jj):
            cells.append((int(i), int(j))); normals.append((float(dx), float(dz)))
    if subsample is not None and len(cells) > subsample:
        idx = np.linspace(0, len(cells) - 1, subsample).astype(int)
        cells = [cells[i] for i in idx]; normals = [normals[i] for i in idx]
    return cells, normals


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
                             trace_steps=200, seed=0, aperture=None, pad_deg=6.0, alpha=0.85,
                             nodal_potential=None):
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
        # Grid cell (cx,cz) occupies [cx,cx+1)x[cz,cz+1). Launch just OUTSIDE its exposed face.
        # The old cx+1.5*n formula was asymmetric and placed floor rays 1.5 cells above the actual
        # interface, producing a 22-44% electron reciprocity bias that only slowly vanished with dx.
        face_eps = 1e-3
        x0c = cx + 0.5 + (0.5 + face_eps) * nnx
        z0c = cz + 0.5 + (0.5 + face_eps) * nnz
        cone = _cone_angles(x0c, z0c, aperture, pad)
        a_mix = alpha if cone is not None else 0.0
        vals = []
        for sc in range(n_scramble):
            s = qmc.Sobol(d=6, scramble=True, seed=seed + sc)
            u = s.random_base2(n_log2)
            E_top = gammadist.ppf(u[:, 0], a=2.0, scale=Te)
            E_surf = E_top + Vc
            m = np.sqrt(np.maximum(np.minimum(E_top, E_surf), 0.0))          # |vperp| cap
            ct = np.sqrt(u[:, 1]); st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
            az = u[:, 2] * (2.0 * np.pi)                          # azimuth: cos(az) projects v_perp into
            vperp_nat = np.sqrt(E_top) * st * np.cos(az)          # the 2D x-z plane (matches the forward
            # source, charging_general.py:353-356). Omitting it (full |v_perp| in-plane) over-broadens the
            # electron angular dist -> over-shadows the floor -> floor over-charges, worse with AR.
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
            # normal-energy retardation from the POLAR angle (vz_surf^2 = E_top*ct^2 + Vc), decoupled
            # from the transverse cos(az) projection -- the two must be separate (Langmuir-exact + correct
            # angular width). Old code conflated them via vz2=E_surf-vperp^2, which only worked because it
            # also (wrongly) used the full |v_perp| in-plane.
            vz2 = E_top * ct * ct + Vc
            valid = (E_surf > 0.0) & (vz2 > 0.0)
            vz_surf = np.sqrt(np.maximum(vz2, 0.0))
            vX = vperp; vZ = -vz_surf
            vdotn = vX * nnx + vZ * nnz
            emit = valid & (vdotn > 0.0)
            # Integrate uniformly over the finite face, not only its center. Face-center collocation
            # creates a deterministic 6-16% reciprocity bias on coarse HAR grids that more rays cannot
            # remove. The tangent (-nz,nx) spans one cell-face length.
            face_s = u[:, 5] - 0.5
            x0 = x0c - nnz * face_s
            z0 = z0c + nnx * face_s
            if nodal_potential is None:
                hix, hiz, _, _, surv, _, _ = _trace_general(
                    Ex, Ez, solid, x0, z0, vX, vZ, -1.0, nx, nz,
                    msteps, trace_dt, trace_dt_field)
            else:
                hix, hiz, _, _, surv, _, _ = trace_nodal(
                    nodal_potential, solid, x0, z0, vX, vZ, -1.0, nx, nz,
                    msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5) & emit
            ratio = vdotn / np.maximum(np.abs(vZ), 0.3)
            # NOTE: a naive score-at-exit reweight |vz_exit|/|vz_implied| amplifies the tracer's dt energy
            # error for barely-escaping electrons on POSITIVE surfaces (the floor) -> worsens it. The
            # unbiased adjoint reweight needs both species + energy-robust exit velocity; deferred.
            vals.append(float((w * escaped * ratio).sum()) / N)
        out[ci] = float(np.mean(vals))
    return out


def backward_electron_floor_liouville(solid, nodal_potential, V_surf, cells, Te=4.0,
                                       n_log2=13, n_scramble=3, trace_dt=0.15,
                                       trace_dt_field=0.10, trace_steps=200, seed=0,
                                       shifted_fraction=0.8):
    """Low-variance, support-complete Liouville gather for horizontal floor faces.

    Surface ``vx`` is Maxwellian. Surface normal energy ``r=w^2`` is sampled from a mixture of the
    natural exponential and that exponential shifted by the local positive barrier ``B=max(Vc,0)``.
    The natural stratum preserves full support; the shifted stratum resolves the exponentially rare
    escaping population. The exact mixture density is divided out, so ``shifted_fraction`` controls
    variance only and is not physical input.
    """
    if not 0.0 <= shifted_fraction < 1.0:
        raise ValueError("shifted_fraction must lie in [0, 1)")
    solid = np.asarray(solid, dtype=bool); V_surf = np.asarray(V_surf, dtype=float)
    nx, nz = solid.shape; N = 2 ** int(n_log2); msteps = int(trace_steps) * nz
    sig = np.sqrt(0.5 * Te); out = np.zeros(len(cells))
    for ci, (cx, cz) in enumerate(cells):
        barrier = max(float(V_surf[cx, cz]), 0.0)
        vals = []
        for sc in range(n_scramble):
            u = qmc.Sobol(d=4, scramble=True, seed=seed + sc).random_base2(n_log2)
            shifted = u[:, 0] < shifted_fraction
            base_r = -Te * np.log(np.clip(1.0 - u[:, 1], 1e-12, 1.0))
            r = base_r + shifted * barrier
            vx = sig * norm.ppf(np.clip(u[:, 2], 1e-9, 1.0 - 1e-9))
            x0 = cx + u[:, 3]; z0 = np.full(N, cz - 1e-3)
            hix, _, _, _, surv, exit_vx, exit_vz = trace_nodal(
                nodal_potential, solid, x0, z0, vx, -np.sqrt(r), -1.0, nx, nz,
                msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5)
            natural = (1.0 - shifted_fraction) * np.exp(-r / Te)
            shifted_density = shifted_fraction * np.where(
                r >= barrier, np.exp(-(r - barrier) / Te), 0.0)
            density_scaled = natural + shifted_density  # common 1/Te cancels analytically
            exit_K = exit_vx * exit_vx + exit_vz * exit_vz
            weight = np.exp((vx * vx - exit_K) / Te) / np.maximum(density_scaled, 1e-300)
            vals.append(float(np.sum(escaped * weight)) / N)
        out[ci] = float(np.mean(vals))
    return out


def backward_ion_gather(solid, Ex, Ez, V_surf, cells, normals, Te=4.0, Ti=0.5, V_dc=37.0, V_rf=30.0,
                        n_log2=13, n_scramble=3, trace_dt=0.15, trace_dt_field=0.10, trace_steps=200,
                        seed=0, aperture=None, pad_deg=3.0, alpha=0.85, want_energy=False,
                        ied_phase_exponent=0.0, exit_state_weight=False,
                        exit_energy_mixture=0.0, nodal_potential=None):
    """Backward ion flux per cell (fraction of incident ion flux; open flat V=0 -> 1).

    Incident ion: uniform RF phase -> Vs=V_dc+V_rf*sin (analytic instantaneous-sheath arcsine IED),
    vz_in=sqrt(0.5 Te+Vs), transverse vperp~N(0,sqrt(0.5 Ti)). ``ied_phase_exponent`` optionally
    applies a Vs^-p phase weight; p=0.35 reproduces the Hwang-Giapis simulated horn ratio but is a
    named BENCHMARK convention, not first-principles input. Near-VERTICAL in the lab
    frame, so v.n>0 on a wall selects only the grazing tail. E_surf_z = 0.5 Te+Vs-Vc; E_surf_z<0 =>
    reflected (floor repels the low-IED-horn = retardation). Flux factor (v.n_cell)/(v.z_surf) suppresses
    grazing ions on vertical walls. aperture -> escape-cone importance sampling (truncated-normal in the
    cone + broad, unbiased); aperture=None -> exact natural sampling.

    ``exit_state_weight`` evaluates the incident phase-space density at the ACTUAL traced plasma-exit
    velocity.  This is the Liouville-consistent adjoint weight required in a nonuniform 2-D field; the
    older 1-D energy map is recovered exactly when vx is conserved and vz_exit^2=E_top_z.  It currently
    requires the first-principles uniform-RF-phase source (``ied_phase_exponent=0``)."""
    if exit_state_weight and ied_phase_exponent != 0.0:
        raise ValueError("exit-state weighting currently requires analytic uniform RF phase")
    if not 0.0 <= exit_energy_mixture < 1.0:
        raise ValueError("exit_energy_mixture must lie in [0, 1)")
    if exit_energy_mixture and not exit_state_weight:
        raise ValueError("broad exit-energy proposal requires exit-state weighting")
    nx, nz = solid.shape
    msteps = int(trace_steps) * nz
    sig = np.sqrt(0.5 * Ti)
    pad = np.deg2rad(pad_deg)
    out = np.zeros(len(cells))
    out_E = np.zeros(len(cells))                                     # flux-weighted mean impact energy [eV]
    for ci, ((cx, cz), (nnx, nnz)) in enumerate(zip(cells, normals)):
        Vc = float(V_surf[cx, cz])
        face_eps = 1e-3
        x0c = cx + 0.5 + (0.5 + face_eps) * nnx
        z0c = cz + 0.5 + (0.5 + face_eps) * nnz
        cone = _cone_angles(x0c, z0c, aperture, pad)
        a_mix = alpha if cone is not None else 0.0
        vals = []; evals = []
        for sc in range(n_scramble):
            s = qmc.Sobol(d=6 if exit_energy_mixture else 5, scramble=True, seed=seed + sc)
            u = s.random_base2(n_log2)
            ph = u[:, 0] * 2.0 * np.pi
            Vs = V_dc + V_rf * np.sin(ph)
            wied = Vs ** (-float(ied_phase_exponent))
            E_surf_z = 0.5 * Te + Vs - Vc
            if exit_energy_mixture:
                # Complete the proposal support in a genuinely 2-D field. Electrostatic work can
                # exchange vx and vz, so surface Ez is not confined to the 1-D shifted RF interval.
                # The upper bound covers the maximum source Ez, the sampled transverse tail, and
                # acceleration from a negative surface. It is a proposal bound, not physical data.
                center = 0.5 * Te + V_dc
                E_cap = center + V_rf + 0.5 * Ti * norm.ppf(1.0 - 1e-6) ** 2 + max(-Vc, 0.0)
                broad = u[:, 5] < exit_energy_mixture
                # Chebyshev/arcsine broad stratum. Uniform energy is a poor proposal for the
                # endpoint-singular RF source density and dominated the estimator variance.
                broad_E = 0.5 * E_cap * (1.0 + np.sin(ph))
                E_surf_z = np.where(broad, broad_E, E_surf_z)
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
            face_s = u[:, 4] - 0.5
            x0 = x0c - nnz * face_s
            z0 = z0c + nnx * face_s
            if nodal_potential is None:
                hix, hiz, _, _, surv, exit_vx, exit_vz = _trace_general(
                    Ex, Ez, solid, x0, z0, vX, vZ, 1.0, nx, nz,
                    msteps, trace_dt, trace_dt_field)
            else:
                hix, hiz, _, _, surv, exit_vx, exit_vz = trace_nodal(
                    nodal_potential, solid, x0, z0, vX, vZ, 1.0, nx, nz,
                    msteps, trace_dt, trace_dt_field)
            escaped = (hix < 0) & (surv < 0.5) & emit
            ratio = vdotn / np.maximum(np.abs(vZ), 0.3)
            if exit_state_weight:
                # At the source, E_z=0.5*Te+Vdc+Vrf*sin(phi) has the arcsine density
                # p(E)=1/[pi*sqrt(Vrf^2-(E-Ec)^2)].  Flux PDF = p_x(vx)*p_vz(vz), while
                # f/incident_flux = PDF/vz = 2*p_x*p(E).  Dividing by the surface proposal
                # p_x(vX)*2*|vZ|*p(E_top) leaves the existing v.n/|vZ| times this density ratio.
                center = 0.5 * Te + V_dc
                exit_Ez = exit_vz * exit_vz
                support_exit = np.abs(exit_Ez - center) < V_rf
                root_exit = np.sqrt(np.maximum(V_rf * V_rf - (exit_Ez - center) ** 2, 1e-24))
                # No weight clipping: vX is already drawn from the explicitly truncated normal above,
                # so the exponent is bounded on the positive side; clipping would bias the adjoint.
                log_px_ratio = (vX * vX - exit_vx * exit_vx) / (2.0 * sig * sig)
                p_exit = np.where(support_exit, 1.0 / (np.pi * root_exit), 0.0)
                mapped_delta = E_surf_z + Vc - center
                mapped_support = np.abs(mapped_delta) < V_rf
                mapped_root = np.sqrt(np.maximum(V_rf * V_rf - mapped_delta * mapped_delta, 1e-24))
                p_mapped = np.where(mapped_support, 1.0 / (np.pi * mapped_root), 0.0)
                proposal_E = (1.0 - exit_energy_mixture) * p_mapped
                if exit_energy_mixture:
                    broad_root = np.sqrt(np.maximum(E_surf_z * (E_cap - E_surf_z), 1e-24))
                    p_broad = np.where(
                        (E_surf_z > 0.0) & (E_surf_z < E_cap),
                        1.0 / (np.pi * broad_root), 0.0)
                    proposal_E = proposal_E + exit_energy_mixture * p_broad
                source_ratio = np.where(proposal_E > 0.0, np.exp(log_px_ratio) * p_exit / proposal_E, 0.0)
                ratio = ratio * source_ratio
            fnum = w * wied * escaped * ratio                        # per-sample flux contribution
            vals.append(float(fnum.sum() / wied.sum()))
            if want_energy:
                E_impact = vperp * vperp + E_surf_z                  # total KE at the surface [eV]
                evals.append(float((fnum * E_impact).sum() / max(fnum.sum(), 1e-12)))
        out[ci] = float(np.mean(vals))
        if want_energy:
            out_E[ci] = float(np.mean(evals)) if evals else 0.0
    return (out, out_E) if want_energy else out


def adaptive_backward_ion_gather(
        solid, Ex, Ez, V_surf, cells, normals, *, base_log2=10, max_log2=16,
        n_replicates=4, seed=0, absolute_tolerance=1e-3, relative_tolerance=5e-3,
        element_absolute_tolerance=None, refine_fraction=0.5, **gather_kwargs):
    """Apply the universal error-controlled surface quadrature to the ion adjoint gather.

    Refinement is selected only from replicate uncertainty and surface measure. No geometry labels or
    region-specific sample counts enter this adapter. Returns ``AdaptiveQuadratureResult``.
    """
    if gather_kwargs.get("want_energy", False):
        raise ValueError("adaptive flux adapter does not yet aggregate the energy observable")
    for reserved in ("n_log2", "n_scramble", "seed"):
        if reserved in gather_kwargs:
            raise ValueError(f"{reserved} is controlled by adaptive_backward_ion_gather")
    cells = [tuple(c) for c in cells]; normals = [tuple(n) for n in normals]

    def evaluator(indices, log2_samples, replicate_seed):
        subset_cells = [cells[int(i)] for i in indices]
        subset_normals = [normals[int(i)] for i in indices]
        return backward_ion_gather(
            solid, Ex, Ez, V_surf, subset_cells, subset_normals,
            n_log2=log2_samples, n_scramble=1, seed=replicate_seed, **gather_kwargs)

    return adaptive_surface_quadrature(
        evaluator, len(cells), base_log2=base_log2, max_log2=max_log2,
        n_replicates=n_replicates, seed=seed, absolute_tolerance=absolute_tolerance,
        relative_tolerance=relative_tolerance,
        element_absolute_tolerance=element_absolute_tolerance,
        refine_fraction=refine_fraction,
    )


def _current_balance_diagnostics(Gi, Ge, comp, cells=None, active_flux=1e-4):
    """Dimensionless floating-current residuals for insulators and pooled conductors.

    ``log(Gi/Ge)`` is the undamped fixed-point voltage residual in units of ``Te``. Cells for which
    both normalized fluxes are below ``active_flux`` are reported but excluded from the active maximum:
    their floating potential is physically underdetermined at the estimator's resolution. Conductors
    are evaluated from their pooled current, matching the voltage update used by the solver.
    """
    Gi = np.asarray(Gi, dtype=float); Ge = np.asarray(Ge, dtype=float); comp = np.asarray(comp)
    eps = 1e-12
    residual = np.full(Gi.shape, np.nan)
    active = np.zeros(Gi.shape, dtype=bool)
    ins = comp == 0
    if cells is None:
        ins_groups = [[i] for i in np.where(ins)[0]]
    else:
        by_cell = {}
        for i in np.where(ins)[0]:
            by_cell.setdefault(tuple(cells[i]), []).append(int(i))
        ins_groups = list(by_cell.values())
    for idx in ins_groups:
        gi = float(Gi[idx].sum()); ge = float(Ge[idx].sum())
        value = float(np.log((gi + eps) / (ge + eps)))
        residual[idx] = value
        active[idx] = (gi + ge) >= active_flux
    pooled = {}
    for c in np.unique(comp[comp > 0]):
        sel = comp == c
        gi = float(Gi[sel].sum()); ge = float(Ge[sel].sum())
        value = float(np.log((gi + eps) / (ge + eps)))
        residual[sel] = value
        active[sel] = (gi + ge) >= active_flux
        pooled[int(c)] = dict(Gi=gi, Ge=ge, log_ratio=value)
    active_values = np.abs(residual[active])
    return dict(
        log_ratio=residual,
        active=active,
        active_count=int(active.sum()),
        inactive_count=int((~active).sum()),
        max_abs_log_ratio=float(active_values.max()) if active_values.size else 0.0,
        rms_log_ratio=float(np.sqrt(np.mean(active_values ** 2))) if active_values.size else 0.0,
        pooled=pooled,
    )


def _laplace_residual(V, gas):
    """Five-point finite-difference Laplace residual on gas cells."""
    xm = np.empty_like(V); xp = np.empty_like(V); zm = np.empty_like(V); zp = np.empty_like(V)
    xm[1:] = V[:-1]; xm[0] = V[0]
    xp[:-1] = V[1:]; xp[-1] = V[-1]
    zm[:, 1:] = V[:, :-1]; zm[:, 0] = V[:, 0]
    zp[:, :-1] = V[:, 1:]; zp[:, -1] = V[:, -1]
    residual = V - 0.25 * (xm + xp + zm + zp)
    values = residual[gas]
    return dict(
        max_abs=float(np.max(np.abs(values))) if values.size else 0.0,
        rms=float(np.sqrt(np.mean(values ** 2))) if values.size else 0.0,
    )


def solve_boundary_state_charging(
        solid, conductor_ids, boundary_state, *, ion_species=None, electron_species=None,
        initial_surface_voltage=None, n_iter=40, beta=0.5, response_energy_eV=4.0, dVmax=8.0,
        balance_tol=1e-3, min_iter=2, field_sweeps=500, field_tolerance=1e-9,
        boundary_proposals=None, n_face_position=8, adaptive_quadrature=None,
        active_flux=1e-4):
    """Geometry- and chemistry-agnostic deterministic charging fixed point.

    ``solid`` is the 2-D material occupancy grid and ``conductor_ids`` assigns zero to locally floating
    dielectric cells and a positive connected-component id to each floating equipotential conductor.
    Every gas-facing material face participates; there are no named floor/wall/mask regions. Species
    physics comes only from ``boundary_state``. The log-current update is a nonlinear root iteration,
    not a physical charging-time integrator; ``response_energy_eV`` and ``beta`` set convergence speed
    but do not change a converged current-balance solution.
    """
    solid = np.asarray(solid, dtype=bool)
    conductor_ids = np.asarray(conductor_ids, dtype=int)
    if int(n_iter) <= 0 or int(min_iter) <= 0:
        raise ValueError("n_iter and min_iter must be positive")
    if conductor_ids.shape != solid.shape or np.any(conductor_ids < 0):
        raise ValueError("conductor_ids must be a nonnegative integer grid matching solid")
    if np.any((conductor_ids > 0) & ~solid):
        raise ValueError("conductor ids may only label solid cells")

    def selected_species(selection, positive):
        if selection is None:
            items = [item for item in boundary_state.species
                     if (item.charge_number > 0) == positive and item.charge_number != 0]
        else:
            names = [selection] if isinstance(selection, str) else list(selection)
            items = [boundary_state.get(name) for name in names]
        expected = 1 if positive else -1
        if not items or any(np.sign(item.charge_number) != expected for item in items):
            label = "positive" if positive else "negative"
            raise ValueError(f"charging requires at least one {label} species")
        if any(item.density_model is None for item in items):
            raise ValueError("charging species require continuous boundary density models")
        return items

    positive_species = selected_species(ion_species, True)
    negative_species = selected_species(electron_species, False)
    positive_incident_current = sum(
        item.flux_m2_s * abs(item.charge_number) for item in positive_species)
    negative_incident_current = sum(
        item.flux_m2_s * abs(item.charge_number) for item in negative_species)
    current_scale = max(positive_incident_current, negative_incident_current, 1e-300)
    proposals = {} if boundary_proposals is None else dict(boundary_proposals)
    cells, normals = _gas_faces(solid, solid)
    if not cells:
        raise ValueError("solid grid has no gas-facing material surface")
    components = np.asarray([conductor_ids[cell] for cell in cells], dtype=int)
    surface_voltage = (np.zeros(solid.shape) if initial_surface_voltage is None
                       else np.asarray(initial_surface_voltage, dtype=float).copy())
    if surface_voltage.shape != solid.shape or not np.all(np.isfinite(surface_voltage)):
        raise ValueError("initial_surface_voltage must be a finite grid matching solid")
    conductor_voltage = np.zeros(int(conductor_ids.max()) + 1)
    for component in range(1, conductor_voltage.size):
        values = surface_voltage[conductor_ids == component]
        conductor_voltage[component] = float(values.mean()) if values.size else 0.0
    history = []; field_history = []; quadrature_history = []; adaptive_levels = {}
    for iteration in range(int(n_iter)):
        for component in range(1, conductor_voltage.size):
            surface_voltage[conductor_ids == component] = conductor_voltage[component]
        potential, field_diag = solve_nodal_laplace(
            solid, surface_voltage, sweeps=field_sweeps, omega=1.7,
            tolerance=field_tolerance)
        field_history.append(field_diag)
        species_current = {}
        species_quadrature = {}
        for species in positive_species + negative_species:
            if adaptive_quadrature is None:
                result = adjoint_boundary_state_face_flux(
                    boundary_state, species.name, potential, solid, cells, normals,
                    proposal_species=proposals.get(species.name), n_face_position=n_face_position)
                normalized_face_flux = result["per_face"]
            else:
                options = dict(adaptive_quadrature)
                bidirectional = bool(options.pop("bidirectional", False))
                forward_options = dict(options.pop("forward_options", {}))
                options.setdefault("n_face_position", n_face_position)
                warm_start_backoff = int(options.pop("warm_start_backoff", 0))
                if warm_start_backoff < 0:
                    raise ValueError("warm_start_backoff must be nonnegative")
                if species.name in adaptive_levels:
                    base_level = int(options.get("base_log2", 6))
                    initial_level = np.maximum(
                        adaptive_levels[species.name] - warm_start_backoff, base_level)
                    options.setdefault("initial_log2_samples", initial_level)
                if bidirectional:
                    # The bidirectional wrapper applies its own cell-level uncertainty gate. It uses
                    # the same tolerances as the adjoint options unless explicitly overridden.
                    element_abs = options.get("element_absolute_tolerance", 1e-3)
                    element_rel = options.get("element_relative_tolerance", 0.05)
                    hybrid = bidirectional_boundary_state_cell_flux(
                        boundary_state, species.name, potential, solid, cells, normals,
                        proposal_species=proposals.get(species.name), adjoint_options=options,
                        forward_options=forward_options,
                        element_absolute_tolerance=element_abs,
                        element_relative_tolerance=element_rel)
                    species_quadrature[species.name] = hybrid
                    if not hybrid["converged"]:
                        raise AdaptiveQuadratureConvergenceError(
                            f"bidirectional phase-space quadrature did not converge for "
                            f"{species.name!r} at fixed-point iteration {iteration + 1}",
                            iteration=iteration + 1, species=species.name,
                            quadrature=hybrid, surface_voltage=surface_voltage,
                            potential=potential, cells=cells, normals=normals)
                    normalized_face_flux = hybrid["per_face"]
                    species_current[species.name] = (
                        normalized_face_flux * species.flux_m2_s * abs(species.charge_number))
                    continue
                adaptive = adaptive_adjoint_boundary_state_face_flux(
                    boundary_state, species.name, potential, solid, cells, normals,
                    proposal_species=proposals.get(species.name), **options)
                species_quadrature[species.name] = adaptive
                adaptive_levels[species.name] = adaptive.log2_samples
                if not adaptive.converged:
                    element_abs = options.get("element_absolute_tolerance")
                    element_rel = float(options.get("element_relative_tolerance", 0.0))
                    if element_abs is None:
                        severity = adaptive.element_stderr
                        allowed = np.full_like(severity, np.nan)
                    else:
                        allowed = float(element_abs) + element_rel * np.abs(adaptive.element_mean)
                        severity = adaptive.element_stderr / np.maximum(allowed, 1e-300)
                    worst = int(np.argmax(severity))
                    message = (
                        f"adaptive phase-space quadrature did not converge for {species.name!r} "
                        f"at fixed-point iteration {iteration + 1}: total stderr="
                        f"{adaptive.total_stderr:.3g}, max face stderr="
                        f"{adaptive.element_stderr[worst]:.3g} at cell={cells[worst]}, "
                        f"normal={normals[worst]}, mean={adaptive.element_mean[worst]:.3g}, "
                        f"allowed face stderr={allowed[worst]:.3g}, "
                        f"surface voltage={surface_voltage[cells[worst]]:.3g} V, "
                        f"level={adaptive.log2_samples[worst]}, max level="
                        f"{adaptive.log2_samples.max()}")
                    raise AdaptiveQuadratureConvergenceError(
                        message, iteration=iteration + 1, species=species.name,
                        quadrature=adaptive, surface_voltage=surface_voltage,
                        potential=potential, cells=cells, normals=normals)
                normalized_face_flux = adaptive.element_mean
            species_current[species.name] = (normalized_face_flux * species.flux_m2_s
                                             * abs(species.charge_number))
        quadrature_history.append(species_quadrature)
        ion_current = np.sum([species_current[item.name] for item in positive_species], axis=0)
        electron_current = np.sum([species_current[item.name] for item in negative_species], axis=0)
        balance = _current_balance_diagnostics(
            ion_current / current_scale, electron_current / current_scale,
            components, cells, active_flux=active_flux)
        history.append(balance)
        if (balance_tol is not None and len(history) >= int(min_iter)
                and balance["max_abs_log_ratio"] <= balance_tol):
            break
        by_cell = {}
        for face_index, cell in enumerate(cells):
            if components[face_index] == 0:
                by_cell.setdefault(cell, []).append(face_index)
        for cell, face_index in by_cell.items():
            gi = float(ion_current[face_index].sum())
            ge = float(electron_current[face_index].sum())
            if (gi + ge) / current_scale < active_flux:
                continue
            surface_voltage[cell] += np.clip(
                beta * response_energy_eV * np.log((gi + 1e-300) / (ge + 1e-300)),
                -dVmax, dVmax)
        for component in range(1, conductor_voltage.size):
            selected = components == component
            if selected.any():
                gi = float(ion_current[selected].sum()); ge = float(electron_current[selected].sum())
                if (gi + ge) / current_scale < active_flux:
                    continue
                conductor_voltage[component] += np.clip(
                    beta * response_energy_eV * np.log((gi + 1e-300) / (ge + 1e-300)),
                    -dVmax, dVmax)
    for component in range(1, conductor_voltage.size):
        surface_voltage[conductor_ids == component] = conductor_voltage[component]
    potential, field_final = solve_nodal_laplace(
        solid, surface_voltage, sweeps=field_sweeps, omega=1.7, tolerance=field_tolerance)
    return dict(
        surface_voltage=surface_voltage, potential=potential,
        cells=np.asarray(cells, dtype=int), normals=np.asarray(normals, dtype=float),
        components=components, ion_current=ion_current, electron_current=electron_current,
        species_current=species_current,
        iterations=len(history), balance_history=history, balance_final=history[-1],
        field_history=field_history, field_final=field_final,
        quadrature_history=quadrature_history,
        adaptive_levels={name: level.copy() for name, level in adaptive_levels.items()},
        conductor_voltage=conductor_voltage,
        current_scale_m2_s=current_scale, active_flux=active_flux,
    )


def self_consistent_backward(g, Te=4.0, n_iter=14, beta=0.5, dVmax=8.0, n_log2=10, n_scramble=2,
                             n_wall=12, n_floor=6, sweeps=250, seed=0, cone_is=False,
                             balance_tol=None, min_iter=6, ion_ied_phase_exponent=0.0,
                             ion_exit_state_weight=False, ion_exit_energy_mixture=0.0,
                             boundary_state=None, ion_species="ion", electron_species="electron",
                             boundary_proposals=None, n_face_position=8):
    """Self-consistent BACKWARD charging solve: Laplace field <-> per-cell gathers <-> damped update
    dV = beta*Te*ln(k*Gi/Ge), where k=Ci/Ce is calibrated each iteration on the floating pillar tops.
    NO forward launch, NO per-region overrides -- the electron-shading dipole EMERGES. Deterministic,
    The fixed-point residual is returned explicitly; a fixed iteration count must not be interpreted as
    convergence. ``balance_tol`` enables optional stopping on the maximum active ``|log(Gi/Ge)|``.

    Reproduces the Kushner picture: upper wall negative (electron shading), deep wall positive (grazing
    ions), the low-AR->high-AR crossover of the potential maximum from floor to sidewall, floor monotone
    in AR. Returns dict with Vs grid, wall/floor profiles, and the sampled cell coords.
    When ``boundary_state`` is supplied, both charged species are gathered from its joint velocity
    density through the boundary-fitted nodal field. This is the unified reactor/sheath/diagnostic path;
    the legacy analytic source remains the default until experimental gates close. ``boundary_proposals``
    may map species names to numerical quadratures and changes variance only, never source physics.

    Cost ~ (n_wall*2 + n_floor + 2) cells * 2^n_log2 * n_scramble * 2 species * n_iter traces."""
    solid = g['solid']; nx, nz = g['nx'], g['nz']
    t0, t1 = g['trench0'], g['trench1']; r0, r1 = int(g['mouth']), int(g['z_poly0']); fz = nz - 1
    cond = g['cond']; gas = ~solid
    red = (np.add.outer(np.arange(nx), np.arange(nz)) % 2 == 0)
    # DENSE surface state: solve EVERY gas-facing insulator cell (per-cell float), so the field is never
    # pinned at 0 V between collocation points (the sparse-solve field-corruption artifact). Observables
    # then average the real solved field, convention-independently.
    ic0, inn0 = _gas_faces(solid, solid & (cond == 0))
    # keep only PHYSICAL trench-facing faces: the gas neighbor must be inside the trench gap [t0,t1) or
    # above the mouth. This drops the finite-array pillars' OUTER faces (facing the open-area simulation
    # boundary -- an array-edge artifact, not a real neighboring feature).
    cells = []; normals = []; kind = []; comp = []
    for (cx, cz), (nnx, nnz) in zip(ic0, inn0):
        gx = cx + int(nnx); gz = cz + int(nnz)
        if not ((t0 <= gx < t1) or gz <= r0):
            continue
        cells.append((cx, cz)); normals.append((float(nnx), float(nnz))); comp.append(0)
        kind.append('floor' if cz == fz else ('mask' if cz <= r0 else 'wall'))
    wrows = np.linspace(r0 + 1, r1 - 2, n_wall).astype(int)      # display sampling of the wall profile
    fcols = np.linspace(t0 + 2, t1 - 3, n_floor).astype(int)
    # floating POLY conductors: pool the current over the TRENCH-FACING inner faces (the physical notch
    # surfaces). NOT the pillar's outer face -- in this finite edge-array that face borders the open-area
    # simulation boundary (an array-edge artifact), and pooling it floats the poly the wrong way
    # (empirically worse). In a real dense line/space array every poly face is trench-facing; _gas_faces
    # restricted to trench-side normals is the general rule. cond 1 = edge pillar, cond 2 = neigh pillar.
    prows = np.linspace(r1 + 1, fz - 1, max(n_wall // 2, 4)).astype(int)
    for z in prows:
        cells.append((t0 - 1, int(z))); normals.append((1.0, 0.0)); kind.append('Lpoly'); comp.append(1)
        cells.append((t1, int(z))); normals.append((-1.0, 0.0)); kind.append('Rpoly'); comp.append(2)
    kind = np.array(kind); comp = np.array(comp)
    clist = [tuple(c) for c in cells]; nlist = [tuple(n) for n in normals]
    Vs = np.zeros((nx, nz)); vc = np.zeros(3)               # vc[1],vc[2] = floating poly potentials
    if boundary_state is not None:
        ion_state = boundary_state.get(ion_species)
        electron_state = boundary_state.get(electron_species)
        if ion_state.charge_number <= 0 or electron_state.charge_number >= 0:
            raise ValueError("unified charging requires positive-ion and negative-electron species")
        if ion_state.density_model is None or electron_state.density_model is None:
            raise ValueError("unified charging species require continuous boundary density models")
        proposals = {} if boundary_proposals is None else dict(boundary_proposals)

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
    balance_history = []
    field_history = []
    for it in range(n_iter):
        for c in (1, 2):
            Vs[cond == c] = vc[c]                            # broadcast floating poly potential onto its body
        if boundary_state is None:
            V = laplace(Vs)
            field_history.append(_laplace_residual(V, gas))
            Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
            Ge = backward_electron_gather(solid, Ex, Ez, Vs, clist, nlist, Te=Te, n_log2=n_log2,
                                          n_scramble=n_scramble, seed=seed, aperture=aperture)
            Gi = backward_ion_gather(solid, Ex, Ez, Vs, clist, nlist, Te=Te, n_log2=n_log2,
                                     n_scramble=n_scramble, seed=seed, aperture=aperture,
                                     ied_phase_exponent=ion_ied_phase_exponent,
                                     exit_state_weight=ion_exit_state_weight,
                                     exit_energy_mixture=ion_exit_energy_mixture)
        else:
            V, field_diag = solve_nodal_laplace(
                solid, Vs, sweeps=sweeps, omega=1.7, tolerance=1e-9)
            field_history.append(field_diag)
            ion_result = adjoint_boundary_state_face_flux(
                boundary_state, ion_species, V, solid, clist, nlist,
                proposal_species=proposals.get(ion_species), n_face_position=n_face_position)
            electron_result = adjoint_boundary_state_face_flux(
                boundary_state, electron_species, V, solid, clist, nlist,
                proposal_species=proposals.get(electron_species), n_face_position=n_face_position)
            Gi = (ion_result["per_face"] * ion_state.flux_m2_s
                  * abs(ion_state.charge_number))
            Ge = (electron_result["per_face"] * electron_state.flux_m2_s
                  * abs(electron_state.charge_number))
        balance = _current_balance_diagnostics(Gi, Ge, comp, clist)
        balance_history.append(balance)
        if balance_tol is not None and it + 1 >= min_iter and balance['max_abs_log_ratio'] <= balance_tol:
            break
        # k = Ci/Ce = 1 EXACTLY (first principle): the domain is the WAFER FRAME (top plane = wafer
        # surface below the sheath; the sheath drop is carried in the ion source energy and the electron
        # gather samples the post-sheath arrival dist). Wafer-scale charge conservation -- zero net DC
        # current is what DEFINES V_dc -- makes Gamma_i = Gamma_e on the uncharged wafer for ANY
        # geometry. So the uncharged wafer is a fixed point by construction; NO mask calibration (that
        # was a noisy estimator of 1 and an arbitrary reference knob). [invariant V1]
        # A corner cell has multiple exposed faces, but only ONE floating potential. Pool the current
        # over all faces of that physical cell and apply one update; treating faces independently gives
        # contradictory updates to the same capacitor and cannot converge.
        ins_by_cell = {}
        for i, cell in enumerate(clist):
            if comp[i] == 0:
                ins_by_cell.setdefault(cell, []).append(i)
        for (cx, cz), idx in ins_by_cell.items():
            gi = Gi[idx].sum(); ge = Ge[idx].sum()
            dVcell = np.clip(beta * Te * np.log((gi + 1e-6) / (ge + 1e-6)), -dVmax, dVmax)
            Vs[cx, cz] += dVcell
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
    if boundary_state is None:
        V = laplace(Vs); final_field = _laplace_residual(V, gas)
    else:
        V, final_field = solve_nodal_laplace(solid, Vs, sweeps=sweeps, omega=1.7, tolerance=1e-9)
    Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
    zp = int(g['z_poly0'])
    foot_h = max(int(0.3 * (t1 - t0)), 3)                  # HG foot band: within 0.3*W of the floor
    fz0 = max(zp + 1, fz - foot_h)
    foot = [(t0 - 1, int(z)) for z in np.arange(fz0, fz)]  # poly-inner face near the poly/oxide junction
    fn = [(1.0, 0.0)] * len(foot)
    if boundary_state is None:
        Gi_f, E_f = backward_ion_gather(
            solid, Ex, Ez, Vs, foot, fn, Te=Te, n_log2=n_log2 + 1,
            n_scramble=n_scramble, seed=seed, aperture=aperture, want_energy=True,
            ied_phase_exponent=ion_ied_phase_exponent, exit_state_weight=ion_exit_state_weight,
            exit_energy_mixture=ion_exit_energy_mixture)
    else:
        foot_result = adjoint_boundary_state_face_flux(
            boundary_state, ion_species, V, solid, foot, fn,
            proposal_species=proposals.get(ion_species), n_face_position=n_face_position,
            want_energy=True)
        Gi_f = foot_result["per_face"]
        E_f = foot_result["mean_impact_energy_eV_per_face"]
    fmask = Gi_f > 1e-6
    E_defl = float(np.sum(E_f[fmask] * Gi_f[fmask]) / max(np.sum(Gi_f[fmask]), 1e-9)) if fmask.any() else 0.0
    # floor_mean over the full _extract band (t0+4:t1-4) -- valid now that the state is DENSE (every
    # floor cell solved), so this is convention-independent and matches the forward _extract exactly.
    floor_mean = float(Vs[t0 + 4:t1 - 4, fz].mean())
    return dict(Vs=Vs, V=V, Ex=Ex, Ez=Ez, wall_rows=wrows, wall_depth=wrows - r0, Lwall=Vs[t0 - 1, wrows],
                Rwall=Vs[t1, wrows], floor=Vs[fcols, fz], floor_mean=floor_mean, k=1.0,
                V_poly=float(0.5 * (vc[1] + vc[2])), Vc=vc.copy(),
                E_defl=E_defl, foot_flux=float(Gi_f[fmask].mean()) if fmask.any() else 0.0,
                iterations=len(balance_history), balance_history=balance_history,
                field_history=field_history, field_final=final_field,
                balance_preupdate=balance_history[-1], sampled_cells=np.asarray(clist, dtype=int),
                sampled_normals=np.asarray(nlist, dtype=float), sampled_kind=kind.copy(),
                sampled_component=comp.copy(), sampled_Gi=Gi.copy(), sampled_Ge=Ge.copy())
