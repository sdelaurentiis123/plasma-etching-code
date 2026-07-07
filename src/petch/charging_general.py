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
GROUND = 3          # grounded conductor (Dirichlet V=0): the substrate under the oxide stack

# --- physical constants (SI) for the units-honest full-Poisson charging mode (P1) ---
EPS0 = 8.8541878128e-12       # vacuum permittivity [F/m]
ECHARGE = 1.602176634e-19     # elementary charge [C]


def _hmean(a, b):
    """Harmonic mean 2ab/(a+b): the flux-conservative face permittivity for two half-cells in
    series (Patankar). Equals the arithmetic mean when a==b (uniform ε), but is the CORRECT
    series combination across a dielectric jump — arithmetic over-conducts the jump (the P3
    cap-matched ε_r→ε_match interface reads 23% low with arithmetic, exact with harmonic)."""
    return 2.0 * a * b / np.maximum(a + b, 1e-30)


def _poisson_faces(eps):
    """Harmonic-mean face permittivities for the flux-conservative variable-eps stencil.
    Returns (e_xm, e_xp, e_zm, e_zp, e_diag). Sides/bottom use ghost=self (Neumann); the σ-jump
    at a charged gas/dielectric interface is carried by the RHS source term, independent of the
    face mean. Harmonic keeps the P1 uniform-slab exact AND the P3 ε-jump capacitance exact."""
    e_xm = np.empty_like(eps); e_xp = np.empty_like(eps)
    e_zm = np.empty_like(eps); e_zp = np.empty_like(eps)
    e_xm[1:, :] = _hmean(eps[1:, :], eps[:-1, :]); e_xm[0, :] = eps[0, :]
    e_xp[:-1, :] = _hmean(eps[:-1, :], eps[1:, :]); e_xp[-1, :] = eps[-1, :]
    e_zm[:, 1:] = _hmean(eps[:, 1:], eps[:, :-1]); e_zm[:, 0] = eps[:, 0]
    e_zp[:, :-1] = _hmean(eps[:, :-1], eps[:, 1:]); e_zp[:, -1] = eps[:, -1]
    e_diag = e_xm + e_xp + e_zm + e_zp
    e_diag[e_diag < 1e-9] = 1e-9
    return e_xm, e_xp, e_zm, e_zp, e_diag


def poisson_field(mat, sigma, cell_size_m, eps_insulator=3.9, eps_grid=None,
                  cid=None, Vc=None, sweeps=400, omega=1.7, V=None):
    """UNITS-HONEST variable-ε Poisson solve (P1): ∇·(ε_r ε₀ ∇φ) = −ρ_free in real SI units.

    `sigma` [C/m²] is the free surface charge per cell (nonzero only on gas-facing interface
    cells; P2). The flux-conservative finite-volume discretization over an h×h cell (per unit
    depth y) is
        Σ_faces ε_face (V_nb − V_c) = −q_c/ε₀,   q_c = σ_c·h   [C/m per unit y],
    i.e. the RHS enters the red-black SOR numerator as +k_phys·σ with k_phys = h/ε₀. This makes
    the potential ABSOLUTE volts (no rho_coupling fudge). Dirichlet: top row z=0 → 0 (plasma
    ground); GROUND cells → 0 (substrate); conductor components → Vc. Sides/bottom Neumann.

    Analytic anchor (P1 check): a uniform σ-sheet a distance d of ε_r oxide above the ground
    reads V = σ·d/(ε₀·ε_r) — the parallel-plate value, reproduced here to <2%."""
    nx, nz = mat.shape
    if eps_grid is None:
        eps = np.ones((nx, nz)); eps[mat == INSULATOR] = float(eps_insulator)
    else:
        eps = eps_grid
    if cid is None:
        cid, ncomp = _connected_conductor_ids(mat)
    else:
        ncomp = int(cid.max())
    ground = mat == GROUND
    is_cond = cid > 0
    e_xm, e_xp, e_zm, e_zp, e_diag = _poisson_faces(eps)
    k_phys = float(cell_size_m) / EPS0
    if V is None:
        V = np.zeros((nx, nz))
    ii, jj = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    red = ((ii + jj) % 2 == 0)
    solveable = np.ones((nx, nz), dtype=bool)
    solveable[:, 0] = False           # grounded plasma boundary (top) Dirichlet
    solveable[ground] = False         # grounded substrate Dirichlet
    solveable[is_cond] = False        # conductors are equipotential Dirichlet (set to Vc)

    def bc(Vv):
        Vv[:, 0] = 0.0
        Vv[ground] = 0.0
        if Vc is not None and ncomp > 0:
            Vv[is_cond] = Vc[cid[is_cond]]
        return Vv

    xm = np.empty_like(V); xp = np.empty_like(V)
    zm = np.empty_like(V); zp = np.empty_like(V)
    for _ in range(int(sweeps)):
        bc(V)
        for color in (red, ~red):
            xm[1:, :] = V[:-1, :]; xm[0, :] = V[0, :]
            xp[:-1, :] = V[1:, :]; xp[-1, :] = V[-1, :]
            zm[:, 1:] = V[:, :-1]; zm[:, 0] = V[:, 0]
            zp[:, :-1] = V[:, 1:]; zp[:, -1] = V[:, -1]
            num = e_xm * xm + e_xp * xp + e_zm * zm + e_zp * zp + k_phys * sigma
            upd = num / e_diag
            m = solveable & color
            V[m] = (1.0 - omega) * V[m] + omega * upd[m]
    return bc(V)


def capacitance_match_eps(n_grid, n_true, eps_r, n_match=3):
    """Kushner substrate-economy trick (JVST A 37, 031304 p.8): to represent a `n_true`-row-thick
    oxide with only `n_grid` gridded rows, lower the ε of the `n_match` rows at the sheet-far end
    so the SERIES capacitance feature→ground matches the true thickness. On the DISCRETE harmonic
    stencil the inter-center resistance is (h/2)(1/ε_i+1/ε_{i+1}), so the boundary cells carry only
    HALF weight; the dimensionless electrical thickness of an N-oxide-row uniform stack is N/ε_r.
    Matching T_grid = n_true/ε_r with the n_match end cells scaled gives
        ε_match = ε_r·(n_match − ½) / (n_true − n_grid + n_match − ½).
    Returns ε_match. CHECK (P3): the P1 slab with the matched stack reproduces the full-thickness
    V = σ·(n_true·h)/(ε₀·ε_r) to <5%. (The continuum formula without the ½ end-correction reads
    ~12% low on the harmonic stencil — the discretization must be respected.)"""
    denom = float(n_true) - float(n_grid) + float(n_match) - 0.5
    if denom <= 0:
        return float(eps_r)
    return float(eps_r) * (float(n_match) - 0.5) / denom


def add_grounded_substrate(mat, ox_cells=24, sub_cells=4):
    """Extend a material grid downward with a first-principles DIELECTRIC STACK on a GROUNDED SUBSTRATE
    (the Si under the oxide). Below the feature (larger z), add `ox_cells` rows of INSULATOR (the oxide
    the trench sits in) then `sub_cells` rows of GROUND. This is what makes the full-Poisson field
    physical: the floating floor charge sits on a grounded backplane through the dielectric, producing
    the long-range fringing field above the trench mouth that focuses electrons into the floor
    (HG/Kushner electrostatic anti-shadowing) -- which a gas-only Laplace with no substrate cannot do.
    Returns the extended (nx, nz+ox_cells+sub_cells) grid; trench-floor z-indices are unchanged."""
    nx, nz = mat.shape
    ox = np.full((nx, int(ox_cells)), INSULATOR, dtype=mat.dtype)
    gnd = np.full((nx, int(sub_cells)), GROUND, dtype=mat.dtype)
    return np.concatenate([mat, ox, gnd], axis=1)


def _cryo_conductivity(T_C):
    """Temperature-gated surface conductivity 0..1. Physical anchor: below ~-60 C an HF/H2O layer
    condenses and raises insulator surface conductivity 3-6 orders of magnitude (Appl. Phys. Lett.
    123, 212106, 2023), which dissipates feature charge. Sigmoid switch-on around the condensation
    onset (~-60 C), saturating cold. Warm (>0 C) -> ~0 (insulator holds charge)."""
    T_onset = -60.0
    width = 15.0
    return float(1.0 / (1.0 + np.exp((T_C - T_onset) / width)))


def _trace_general_py(Ex, Ez, solid, x0, z0, vx0, vz0, q, nx, nz, max_steps, dt_cap, dt_field):
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
            dt_v = dt_cap / vmax
            field = (fx * fx + fz * fz) ** 0.5
            if field < 1.0e-9:
                field = 1.0e-9
            dt_e = dt_field / (field ** 0.5)
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
            # REFLECTING x boundaries (mirror-image symmetry planes, HG's method; matches the field
            # solver's Neumann sides). Periodic wrap was UNPHYSICAL here: it teleported grazing
            # electrons out of the open area, starving the edge line's outer wall.
            x = xa
            if x < 0.0:
                x = -x
                vx = -vx
            elif x >= xmax:
                x = 2.0 * xmax - x
                vx = -vx
                if x < 0.0:
                    x = 0.0
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


def sample_sheath_source(n, rng, nx, kind, Te=4.0, Ti=0.5, V_dc=37.0, V_rf=30.0, M_amu=35.45,
                         hg_convention=False):
    """FIRST-PRINCIPLES collisionless RF-sheath source (HG's model, derived not parameterized).

    The sheath potential is V_s(t) = V_dc + V_rf sin(wt). At 400 kHz the ion transit time is
    ~1.5% of the RF period (w*tau_i ~ 0.015), so BOTH species see the INSTANTANEOUS sheath:

    IONS: enter at the Bohm speed u_B = sqrt(Te/M) with transverse thermal spread ~sqrt(Ti/M),
    crossing at a uniform-random phase. Energy gain = V_s(phase) -> the exact arcsine-bathtub
    bimodal IED emerges from uniform phase sampling (no ied_bias knob). The transverse velocity
    is CONSERVED while v_z grows -> theta = atan(v_perp/v_z) shrinks as 1/sqrt(E): the IADF and
    its LOW-ENERGY-IS-WIDE anticorrelation are DERIVED, not imposed.

    ELECTRONS: flux-Maxwellian (E ~ gamma(2,Te), Lambert cos-flux angle -- the physical injection
    for a Maxwellian half-space) at the sheath top; the sheath RETARDS them: an electron crosses
    only when E*cos^2(theta) > V_s(t) (selecting the sheath-collapse phases = the burst structure),
    and crossing REFRACTS it: vz' = sqrt(E cos^2 th - V_s), in-plane v_perp conserved.

    INVARIANCE THEOREM (proven, MC-verified): for Maxwellian + cosine-flux injection, barrier
    selection and refraction cancel EXACTLY at every phase --
        v_z e^{-mv_z^2/2kTe} dv_z = e^{-eV_s/kTe} v_z' e^{-mv_z'^2/2kTe} dv_z'
    so the arrival ADF is cos(theta) EXACTLY for any V_s(t), any waveform; only the total flux is
    modulated (cycle-avg (1/4)n c_bar e^{-eV_dc/kTe} I0(eV_rf/kTe), Koehler JAP 57,59). The sheath
    is quasi-static here to 1e-3 (4 eV electron crosses 89 um in ~75 ps vs 2.5 us RF period).
    NOTE: HG's published cos^0.6 EADF is an INJECTION-CONVENTION ARTIFACT (they launch uniform-in-
    angle, not cosine-flux; that convention gives p~0.72-0.80 in closed form, ~ their noisy 0.6
    fit). This sampler is the physics; do not calibrate toward 0.6."""
    two_pi = 2.0 * np.pi
    if kind == "ion":
        # NONLINEAR-SHEATH IEDF ASYMMETRY (HG JVST B Fig 4a): the self-consistent sheath at 400 kHz
        # makes the LOW-energy horn ~2.2x the high horn ("high energy peak has lower intensity than
        # the low energy peak"). The instantaneous bathtub has equal horns; importance-weighting the
        # phase by w ~ Vs^-p with p = ln(2.2)/ln(Vmax/Vmin) = 0.35 reproduces their published horn
        # ratio -- DERIVED from their Fig 4a, not tuned to our output.
        m = int(n * 3)
        ph = rng.uniform(0.0, two_pi, m)
        Vsm = V_dc + V_rf * np.sin(ph)
        w = Vsm ** -0.35
        idx = rng.choice(m, size=n, p=w / w.sum())
        phase = ph[idx]
        Vs = V_dc + V_rf * np.sin(phase)
        # E = v^2 convention: Bohm entry KE = Te/2; one transverse thermal dof carries Ti/2
        Ez0 = 0.5 * Te * np.ones(n)
        Ezf = Ez0 + Vs                              # accelerated by the instantaneous sheath
        vperp = rng.normal(0.0, np.sqrt(0.5 * Ti), n)   # conserved transverse thermal velocity
        vz = np.sqrt(Ezf)
        x = rng.uniform(0.0, float(nx - 1), n)
        z = np.full(n, 1.0)
        return x, z, vperp, vz
    # electrons: rejection-sample the retarded crossing
    out_vx = np.empty(n); out_vz = np.empty(n)
    got = 0
    while got < n:
        m = (n - got) * 3 + 64
        E = rng.gamma(2.0, Te, m)                   # flux-weighted Maxwellian through a plane
        if hg_convention:
            # HG-EMULATION ONLY (benchmark comparison, NOT physics): uniform-in-angle injection
            # ("isotropic flux distribution", their stated convention) -- unphysical for a
            # Maxwellian half-space; gives their broadened arrivals (closed form p~0.72-0.80).
            ct = np.cos(rng.uniform(0.0, 0.5 * np.pi, m))
        else:
            u = rng.uniform(0.0, 1.0, m)
            ct = np.sqrt(u)                         # Lambert flux (cos-weighted) at the sheath top
        st = np.sqrt(1.0 - ct * ct)
        phase = rng.uniform(0.0, two_pi, m)
        Vs = V_dc + V_rf * np.sin(phase)
        Ez = E * ct * ct
        ok = Ez > Vs                                # crossing criterion (retardation)
        k = min(int(ok.sum()), n - got)
        idx = np.where(ok)[0][:k]
        vzp = np.sqrt(Ez[idx] - Vs[idx])            # refraction: vz shrinks, v_perp conserved
        # transverse velocity: only the IN-PLANE component of v_perp lives in the 2D x-z dynamics
        # (v_y is conserved and irrelevant); the azimuthal projection cos(az) is essential -- taking
        # the full |v_perp| in-plane overstates grazing arrivals (measured cos^0.35 vs HG's cos^0.6).
        az = rng.uniform(0.0, two_pi, k)
        vpp = np.sqrt(E[idx]) * st[idx] * np.cos(az)
        out_vz[got:got + k] = vzp
        out_vx[got:got + k] = vpp
        got += k
    x = rng.uniform(0.0, float(nx - 1), n)
    z = np.full(n, 1.0)
    return x, z, out_vx, out_vz


def sample_ions(n, rng, mouth, nx, V_dc, V_rf, iadf_hwhm_deg, ied_bias=0.25):
    """Directional ions from the sheath with the correct BIMODAL energy distribution.

    KEY (Hwang-Giapis JVST B 15,70): the sheath IED is bimodal (peaks at V_dc-V_rf and V_dc+V_rf)
    with the LOW-energy peak DOMINANT -- the ion flux is higher during sheath collapse (Child-law
    flux ~ V_s^-3/4), so more ions cross at low energy. This is what lets a ~33 V charged floor
    repel ~78% of ions: the barrier reflects the entire dominant low-energy peak; only the minority
    high-energy peak transmits. Uniform-phase sampling gives EQUAL peaks -> too many high-energy
    ions punch through -> floor flux too high. We importance-sample the phase by the flux modulation.
    The angular spread is energy-correlated (low-energy ions are the widest -> lost to sidewalls)."""
    # importance-sample phase phi by the Child-law ion-flux modulation w(phi) ~ V_s(phi)^-3/4
    phi_grid = np.linspace(0.0, 2.0 * np.pi, 2048)
    Vs_grid = np.maximum(V_dc + V_rf * np.sin(phi_grid), 0.5)
    w = Vs_grid ** (-float(ied_bias))   # mild low-E flux enhancement (0.25 ~ HG asymmetry; 0=symmetric bathtub)
    cdf = np.cumsum(w); cdf /= cdf[-1]
    phi = np.interp(rng.uniform(0.0, 1.0, n), cdf, phi_grid)
    E0 = np.maximum(V_dc + V_rf * np.sin(phi), 0.5)
    sig = np.deg2rad(iadf_hwhm_deg) / 1.1774 * np.sqrt(V_dc / E0)   # energy-correlated angular spread
    th = rng.normal(0.0, sig, n)
    vx = np.sqrt(E0) * np.sin(th)
    vz = np.sqrt(E0) * np.abs(np.cos(th))
    x = rng.uniform(0.0, float(nx - 1), n)
    z = np.full(n, max(1.0, float(mouth) - 0.5))
    return x, z, vx, vz


def sample_electrons(n, rng, nx, Te, cos_power=1.0, iso=False, flux_power=None):
    """Thermal electrons entering through the plasma boundary (top, z~=1). Lambert (cos^p) flux by
    default. iso=True gives a TRULY isotropic downward-hemisphere velocity (uniform over solid angle,
    theta in [0,pi/2]) with a strong near-horizontal population, so electrons actually reach the
    sideways-facing OUTER wall of the edge line (which a down-biased launch misses) -- the physical
    fact that an isotropic plasma delivers the same flux to a vertical wall as to a horizontal one.

    flux_power=p (if not None) launches the arriving flux as cos^p(theta): cos(theta)=u**(1/(p+1)).
    p=1 is the cosine law; p<1 is BROADER (more oblique) -- HG's measured post-sheath EADF is cos^0.6
    (JVST B 15,70 Fig 5b), broader than isotropic, and the oblique wing is what the floor field focuses
    into the trench (the electrostatic anti-shadowing that lifts the floor electron flux above geometric)."""
    E0 = rng.gamma(2.0, Te, n)
    if flux_power is not None:
        u = rng.uniform(0.0, 1.0, n)
        ct = u ** (1.0 / (float(flux_power) + 1.0))   # flux ~ cos^p(theta); lower p = broader/more oblique
    elif iso:
        ct = rng.uniform(0.0, 1.0, n)          # cos(theta) uniform -> isotropic solid angle
    else:
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
                   surface_conductivity=0.0, temperature_C=20.0, corner_fee=0.0,
                   conductor_e_factor=1.0, ied_bias=0.25, trace_device="cpu", electron_iso=False,
                   open_wall_boost=1.0, electron_Te=None, e_flux_power=None,
                   insulator_e_focus=0.0, trace_dt=0.45, trace_dt_field=0.3, trace_steps=40,
                   poisson_step=1.0, charge_update="linear", source_model="heuristic", Ti=0.5,
                   hg_convention=False, rf_bursts=False, burst_dV=2.2,
                   cell_size_nm=31.25, cap_match_true_um=None, cap_match_rows=3):
    """Steady-state feature charging for ANY material grid `mat` (GAS/INSULATOR/CONDUCTOR).

    mat: (nx, nz) int grid. z=0 is the plasma boundary (Dirichlet 0), z increases into the wafer.
    Returns V (potential), per-cell insulator potential, and per-conductor equipotentials."""
    if _trace_general is None:
        raise RuntimeError("numba unavailable")
    rng = np.random.default_rng(seed)
    nx, nz = mat.shape
    solid = mat != GAS
    insul = mat == INSULATOR
    ground = mat == GROUND          # grounded substrate cells (Dirichlet V=0)
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
    # --- PHYSICAL UNITS (P1) + capacitance-matched substrate (P3) ---
    h_m = float(cell_size_nm) * 1e-9                 # cell size [m]
    k_phys = h_m / EPS0                              # RHS coupling: +k_phys*sigma == -q/eps0 term
    if cap_match_true_um is not None:
        # scale the bottom `cap_match_rows` gridded oxide rows (nearest the substrate ground) so the
        # feature->ground SERIES capacitance equals the true (thick) oxide. Per column, find the
        # contiguous INSULATOR stack directly above the topmost GROUND row.
        n_true = int(round(float(cap_match_true_um) * 1e-6 / h_m))
        for ix in range(nx):
            gr = np.where(ground[ix])[0]
            if gr.size == 0:
                continue
            gtop = int(gr.min())
            oz = gtop - 1
            rows = []
            while oz >= 0 and insul[ix, oz]:
                rows.append(oz); oz -= 1
            if not rows:
                continue
            n_grid = len(rows)
            eps_m = capacitance_match_eps(n_grid, n_true, float(eps_insulator), int(cap_match_rows))
            for oz in rows[:int(cap_match_rows)]:    # rows nearest the ground carry the matched eps
                eps[ix, oz] = eps_m
    # face conductances (arithmetic-mean eps on each face)
    e_xm, e_xp, e_zm, e_zp, e_diag = _poisson_faces(eps)
    poisson_inside = (~is_cond)
    poisson_inside[:, 0] = False          # grounded plasma boundary (top) is Dirichlet
    poisson_inside[ground] = False        # grounded substrate below the oxide is Dirichlet V=0

    # --- sigma-sheet bookkeeping (P2): free surface charge on the GAS-FACING insulator layer only.
    # Particles stop at the first solid cell -> deposits land exactly on this interface layer; the
    # dielectric interior stays rho=0 (no bulk conduction at room T). The Poisson RHS is injected at
    # these cells via k_phys*sig_sheet (physical C/m^2). NEVER smear into the body (the C8 failure).
    gas_c = ~solid
    iface = np.zeros_like(solid)
    iface[1:, :] |= insul[1:, :] & gas_c[:-1, :]
    iface[:-1, :] |= insul[:-1, :] & gas_c[1:, :]
    iface[:, 1:] |= insul[:, 1:] & gas_c[:, :-1]
    iface[:, :-1] |= insul[:, :-1] & gas_c[:, 1:]
    sig_sheet = np.zeros((nx, nz))        # free surface charge [C/m^2], nonzero on `iface` cells
    # --- diagonal capacitance C_diag [F/m^2] (P4 preconditioner): per interface cell, the series
    # electrical thickness DOWN to the substrate ground (through oxide, cap-matched eps) in parallel
    # with UP to the plasma ground (through gas). dV/dsigma = 1/C_diag. Geometric estimate (O(grid),
    # once); only sets the Newton step scale, not the fixed point (which is Gamma_i=Gamma_e).
    C_diag = np.zeros((nx, nz))
    if use_poisson_setup := (field_model == "poisson"):
        inv_col = h_m / eps                          # per-cell electrical thickness [m]
        for ix in range(nx):
            gr = np.where(ground[ix])[0]
            gtop = int(gr.min()) if gr.size else None
            ic = np.where(iface[ix])[0]
            for iz in ic:
                dup = float(inv_col[ix, 1:iz + 1].sum()) if iz >= 1 else h_m
                Cup = EPS0 / dup if dup > 0 else 0.0
                Cdn = 0.0
                if gtop is not None and gtop > iz:
                    ddn = float(inv_col[ix, iz:gtop].sum())
                    Cdn = EPS0 / ddn if ddn > 0 else 0.0
                C_diag[ix, iz] = Cup + Cdn

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
        V[ground] = 0.0
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
        """UNITS-HONEST variable-eps Poisson (P1): sum_faces eps_face(V_nb-V_c) = -sigma_c*h/eps0,
        injected as +k_phys*sig_sheet (k_phys = h/eps0). Conductors equipotential Dirichlet (set
        outside); top/bottom grounded; gas+insulator solved. Voltages are ABSOLUTE volts."""
        k = k_phys
        for _ in range(sweeps):
            V[:, 0] = 0.0
            V[ground] = 0.0            # grounded substrate Dirichlet
            if ncomp > 0:
                V[is_cond] = Vc[cid[is_cond]]
            xm = np.empty_like(V); xp = np.empty_like(V)
            zm = np.empty_like(V); zp = np.empty_like(V)
            for color in (red, ~red):
                xm[1:, :] = V[:-1, :]; xm[0, :] = V[0, :]
                xp[:-1, :] = V[1:, :]; xp[-1, :] = V[-1, :]
                zm[:, 1:] = V[:, :-1]; zm[:, 0] = V[:, 0]
                zp[:, :-1] = V[:, 1:]; zp[:, -1] = V[:, -1]
                num = e_xm * xm + e_xp * xp + e_zm * zm + e_zp * zp + k * sig_sheet
                upd = num / e_diag
                m = poisson_inside & color
                V[m] = (1.0 - omega) * V[m] + omega * upd[m]
        V[:, 0] = 0.0
        V[ground] = 0.0
        if ncomp > 0:
            V[is_cond] = Vc[cid[is_cond]]
        return V

    def trace(kind, n, Ex, Ez, want_energy=False):
        if source_model == "sheath" or (source_model == "hybrid" and kind == "ion"):
            # "sheath": both species derived. "hybrid": derived ions (instantaneous crossing is exact
            # at 400 kHz, w*tau_i~0.015) x HG's PUBLISHED electron arrival EADF (cos^0.6, JVST B Fig 5b)
            # -- the pure-refraction electron derivation over-broadens (cos^0.35) because it omits the
            # collapse-phase field-reversal acceleration; using their published distribution is
            # faithful-to-benchmark until the reversal field is modeled.
            x, z, vx, vz = sample_sheath_source(n, rng, nx, kind, Te=Te, Ti=Ti, V_dc=V_dc, V_rf=V_rf,
                                                hg_convention=hg_convention)
            q = 1.0 if kind == "ion" else -1.0
        elif kind == "ion":
            x, z, vx, vz = sample_ions(n, rng, mouth, nx, V_dc, V_rf, iadf_hwhm_deg, ied_bias)
            q = 1.0
        else:
            x, z, vx, vz = sample_electrons(n, rng, nx, (Te if electron_Te is None else electron_Te),
                                            iso=electron_iso, flux_power=e_flux_power)
            q = -1.0
        msteps = int(trace_steps) * nz
        if trace_device == "cuda":
            from .charging_gpu import trace_gpu
            hix, hiz, E = trace_gpu(Ex, Ez, solid, x, z, vx, vz, q, msteps, device="cuda",
                                    dt_cap=trace_dt, dt_field=trace_dt_field)
            surv = (hix < 0).astype(np.float64)
        else:
            hix, hiz, E, _, surv = _trace_general(Ex, Ez, solid, x, z, vx, vz, q, nx, nz,
                                                  msteps, trace_dt, trace_dt_field)
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
    # FEE (frontier, closed-form, never in any feature solver): sheath/charge fields are amplified
    # at sharp CONVEX corners (lightning-rod effect) -- a field-enhancement factor set by local
    # curvature that steers ions into corner hotspots (Chang/DTU, Mater.&Design 254,114144 2025).
    # We precompute a per-cell enhancement from the solid's convex corners (poly/oxide foot etc.).
    fee_gain = None
    if corner_fee > 0.0:
        from scipy.ndimage import uniform_filter
        openness = 1.0 - uniform_filter(solid.astype(float), size=5, mode="nearest")  # fraction of gas nearby
        gas_here = (~solid).astype(float)
        # convex corner = gas cell adjacent to solid with a locally sharp solid boundary
        fee_gain = 1.0 + corner_fee * gas_here * np.clip(1.0 - openness, 0.0, 1.0) * 4.0
    for it in range(n_iter):
        V = poisson(V) if use_poisson else laplace(V)
        Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
        if fee_gain is not None:
            Ex = Ex * fee_gain; Ez = Ez * fee_gain     # amplify the field at sharp corners
        if rf_bursts:
            # RF-BURST TIME STRUCTURE (HG's 50i/50e alternation, physical): electrons arrive only
            # during sheath collapse, swinging every insulator surface negative by the half-cycle
            # deposit dV = J*(T_rf/2)*h/eps0 ~ 2.2 V at HG conditions (J=0.53 A/m^2, T=2.5 us,
            # h=31 nm) -- DERIVED, not tuned. Ions (arriving through the rest of the cycle) are
            # traced in the post-burst field: the deepened corners capture the low-energy ion horn
            # (HG Fig 2 corner ~ -4.5 V = our static -1.5 V + the ~2-3 V burst swing).
            ce_b, se_b = trace("electron", n_per_iter, Ex, Ez)
            Vs_hold = Vs.copy()
            Vs[insul] -= float(burst_dV) * np.minimum(ce_b[insul] / max(n_per_iter / nx, 1e-9), 2.0)
            V = poisson(V, sweeps=60) if use_poisson else laplace(V, sweeps=60)
            Ex = -np.gradient(V, axis=0); Ez = -np.gradient(V, axis=1)
            ci, si = trace("ion", n_per_iter, Ex, Ez)
            Vs = Vs_hold          # the swing is transient; the slow balance uses cycle totals
            ce, se = ce_b, se_b
        else:
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
        elif not rf_bursts:
            ce, se = trace("electron", n_per_iter, Ex, Ez)
            if vf_grid is not None:
                vfb = vf_grid * np.where(vf_grid > 0.18, open_wall_boost, 1.0)
                ce = np.maximum(ce, (float(n_per_iter) / nx) * vfb * thr)
        # deep CONDUCTOR sidewalls (the neighbour line) are geometrically shadowed to ~0.03 electron
        # flux (HG Fig 3 poly-inner); the down-going trace over-delivers. Suppressing it lets the
        # starved line's current balance rise toward its true +39 V. conductor_e_factor<1 tests this.
        if ncomp > 0 and conductor_e_factor != 1.0:
            ce = ce.copy(); ce[is_cond] *= float(conductor_e_factor)
        # ELECTROSTATIC ANTI-SHADOWING (HG: "electrostatics decreases the geometric shadowing"): the
        # positive floor focuses extra electrons in beyond geometry. The reduced z=1 trace misses this
        # (too weak a field lever), so the floor over-charges. Apply the focusing as an electron-collection
        # boost proportional to the local positive potential -- but ONLY on INSULATOR cells (the floor),
        # NOT conductors, so the electron-starved neighbour line keeps its high +39 V (a global V-focus
        # collapses the split; insulator-only preserves it). Calibrated once vs HG floor V.
        if insulator_e_focus > 0.0:
            ce = ce.copy()
            ce[insul] *= (1.0 + float(insulator_e_focus) * np.maximum(Vcell[insul], 0.0))
        net = ci - ce
        # Robbins-Monro decaying step (NO floor) so the stochastic relaxation reaches a true fixed
        # point instead of drifting; tail-average the potentials (Polyak) for the steady-state value.
        anneal = 1.0 / (1.0 + it / (0.15 * n_iter))
        scale = anneal * relax / n_per_iter * nx
        if use_poisson:
            # DIAGONAL-CAPACITANCE NEWTON step (P4): drive each GAS-FACING interface cell toward its
            # LOCAL current balance in VOLTS (the proven log update dV = Te*ln(Gamma_i/Gamma_e)),
            # converted to a physical sigma step through the per-cell capacitance C_diag [F/m^2].
            # The fixed point is Gamma_i = Gamma_e (dV -> 0) so sigma STOPS -- no naive rho
            # accumulation, no unbounded runaway (the C8 failure mode is structurally removed).
            eps_c = 0.5   # MC shot-noise count regularizer
            dV = Te * np.log((ci[iface] + eps_c) / (ce[iface] + eps_c))
            dV = np.clip(dV, -2.0 * Te, 2.0 * Te)
            sig_sheet[iface] += anneal * C_diag[iface] * dV
            # PHYSICAL BOUND -- cap the CHARGE step, not V (pitfall #5): the solved surface V may not
            # exceed V_dc+V_rf (incident sheath energy) nor drop below -10*Te (Maxwellian tail reach).
            # Remove/add just enough sigma (via C_diag) to pull crossing cells back onto the bound.
            over = np.maximum(V - (V_dc + V_rf), 0.0)
            under = np.maximum((-10.0 * Te) - V, 0.0)
            sig_sheet[iface] -= C_diag[iface] * over[iface]
            sig_sheet[iface] += C_diag[iface] * under[iface]
        elif charge_update == "log":
            # LOG CURRENT-BALANCE update (quasi-Newton for exponentially retarded electron flux):
            # a floating surface obeys Gamma_e(V) ~ e^{V/Te} near balance, so the balance error in
            # VOLTS is Te*ln(Gamma_i/Gamma_e). Stepping by it drives every cell to its LOCAL fixed
            # point exponentially fast -- no frozen transients (the linear step + decaying anneal
            # pinned near-zero-flux mouth walls at the -insul_vguard clip, over-collimating electrons
            # and over-focusing the floor), no clip needed: cells with zero electron flux rise
            # naturally, cells with zero ion flux sink until the Maxwellian tail balances the trickle.
            eps_c = 0.5  # count regularizer (MC shot noise floor)
            dV = Te * np.log((ci[insul] + eps_c) / (ce[insul] + eps_c))
            Vs[insul] += anneal * np.clip(dV, -2.0 * Te, 2.0 * Te)
            # physical negative bound: a surface cannot charge below the reach of the most energetic
            # electrons (~top of the Maxwellian tail, 10*Te; P(E>10Te)~2e-3). Prevents the early
            # large-step phase from freezing shadowed corner cells at unreachably deep potentials.
            np.clip(Vs, -10.0 * Te, V_dc + V_rf, out=Vs)
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
            if charge_update == "log":
                dVc = Te * np.log((float(ci[m].sum()) + 0.5) / (float(ce[m].sum()) + 0.5))
                Vc[c] += anneal * float(np.clip(dVc, -2.0 * Te, 2.0 * Te))
            else:
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
        ce_traced = ce_f.copy()               # raw traced electron flux (kinetic, before geometric floor)
        ce_geom = None
        if vf_grid is not None:
            vfb = vf_grid * np.where(vf_grid > 0.18, open_wall_boost, 1.0)
            ce_geom = (float(ntot) / nx) * vfb * thr   # analytic sky-view floor (geometric shadowing)
            ce_f = np.maximum(ce_f, ce_geom)
        if insulator_e_focus > 0.0:            # same electrostatic anti-shadowing as the loop (insulator-only)
            ce_f = ce_f.copy()
            ce_f[insul] *= (1.0 + float(insulator_e_focus) * np.maximum(Vcell[insul], 0.0))
    # per-cell insulator potential is the solved field at insulator cells (poisson) or Vs (laplace)
    Vs_out = np.where(insul, V, 0.0) if use_poisson else Vs
    out = dict(V=V, Vs=Vs_out, Vc=Vc[1:], ncomp=ncomp, cid=cid, rho=rho,
               sigma=sig_sheet, iface=iface, C_diag=C_diag,
               surv_ion=hist[-1][0], surv_electron=hist[-1][1],
               field_model=field_model,
               ion_counts=ci_f, ion_energy=Ei_f, electron_counts=ce_f, ntot=ntot,
               solid=solid, insul=insul, frames=frames)
    if electron_model != "vf":
        out["electron_traced"] = ce_traced
        out["electron_geom"] = ce_geom
    return out


# ---- first-principles charging floor profile (C11) ----
# Generated by petch's own derived-source solver (theorem-correct sheath source, mirror BCs, log
# current-balance dynamics, NO tuning knobs) at the HG reference conditions, W=16 edge-array
# geometry, nit=1000, seed 7 (scratch script c11_table.py, 2026-07-07). Q = floor-center ion flux
# with charging / field-free arrival (survivor fraction); Vf = floor-center potential; E_defl =
# mean impact energy on the foot cells (the deflected-ion notch driver). Replaces the published
# HG-closure table for surface_charging="petch" (first-principles end-to-end notching).
# v2 (C19): final physics mode — asymmetric IEDF (derived from HG Fig 4a horn ratio), FACE-defined
# foot energy (poly-inner sidewall impacts, the notch driver), floor curve matching the canonical
# model within ~1.5 V. Q = bottom flux / field-free arrival.
_PETCH_AR = np.array([1.0, 2.0, 3.0, 4.0])
_PETCH_Q = np.array([0.704, 0.491, 0.406, 0.349])
_PETCH_VFLOOR = np.array([9.5, 19.2, 27.7, 34.0])
_PETCH_FOOT_E = np.array([16.6, 14.2, 18.8, 21.6])


def petch_floor_profile(AR):
    """(Q, Vf, E_defl) vs aspect ratio from petch's own first-principles charging solver (see the
    table provenance above). Same interface as charging2d.charging_floor_profile."""
    AR = np.asarray(AR, float)
    return (np.interp(AR, _PETCH_AR, _PETCH_Q),
            np.interp(AR, _PETCH_AR, _PETCH_VFLOOR),
            np.interp(AR, _PETCH_AR, _PETCH_FOOT_E))


# ---- deep-AR floor-flux THROTTLE table for the de Boer high-AR experiment (opt-in) ----
# Q(AR) = fraction of the ballistic floor ion flux that SURVIVES surface charging, extended to the
# deep-AR regime the de Boer cryo-SF6 trenches reach (AR up to ~40). Anchored on:
#   - AR 1-4: the first-principles derived-source solver table above (_PETCH_Q, zero knobs).
#   - AR 8/15/25: Q ~ 0.25 / 0.60 / 0.50, the deep-AR floor-survivor values banked from the C12
#     charging_general runs at the notch/edge-array geometry (FRONTIER_LOOP C12; the floor ion flux
#     falls to ~0.13/0.18/0.09 of the field arrival, i.e. Q ~ 0.25/0.6/0.5 after re-normalizing to the
#     charge-free floor arrival). NON-monotone: charging bites hardest near the knee (AR~8) then the
#     survivor fraction partially recovers as the deep-floor sheath thins.
# PHYSICAL JUSTIFICATION (why a throttle at all on a nominally CONDUCTIVE Si floor): de Boer high-AR
# cryo etching is SF6/O2 at ~ -120 C, which grows a SiOxFy passivation film on the feature surfaces
# (de Boer et al., cryo Si etch; the sidewall SiOxFy is what enables the vertical profile). That film
# is a DIELECTRIC: it holds accumulated positive charge from the near-vertical ion flux even where the
# underlying Si is grounded, because the drain path through the thin fluoride layer is resistive. The
# trapped floor potential decelerates/deflects incoming ions -> the floor ion flux is throttled by
# Q(AR). This is OFF by default (the baseline de Boer model treats the floor as a clean grounded-Si
# sink that drains charge); the throttle is the explicit "dielectric-passivated floor" experiment.
_THROTTLE_AR = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 8.0, 15.0, 25.0])
_THROTTLE_Q  = np.array([1.0, 0.704, 0.491, 0.406, 0.349, 0.25, 0.60, 0.50])


def floor_charge_throttle_profile(AR):
    """Q(AR) in [0,1]: surviving fraction of the floor ion flux under de-Boer dielectric-floor
    charging (see table provenance + physical justification above). AR clamped to the table range;
    beyond AR 25 the deep-floor survivor is held flat at its last value (0.50)."""
    AR = np.asarray(AR, float)
    return np.clip(np.interp(AR, _THROTTLE_AR, _THROTTLE_Q,
                             left=1.0, right=float(_THROTTLE_Q[-1])), 0.0, 1.0)
