"""Element-resolved two-reservoir mixed-layer surface chemistry (v1).

Implements the design in RESEARCH_MIXED_LAYER_DESIGN_2026-07-23.md: a
fluorocarbon film reservoir (n_C_film, n_F_film) stacked on a mixed reaction
layer (n_Si, n_O, n_C, n_F), all areal densities in atoms/m^2, advanced by
explicit conservative updates. Every energy, range, and mixing quantity comes
from the ZBL/Lindhard stopping module (ion_energy_deposition) — there is no
fitted energy law. Product stoichiometry (SiF4, CO, COF2) is fixed by
chemistry, so the per-element ledgers close by construction; the integrator
additionally enforces them to machine precision and reports the residuals.

v1 is standalone: nothing in the feature engine imports it. The validation
ladder (design doc section 5.5) wires it in only after Rung 0 regression.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from petch.ion_energy_deposition import (
    FLUOROCARBON_FILM,
    SIO2,
    nuclear_energy_in_layer_eV,
    projected_range_nm,
)

# SiO2 formula-unit density (formula units / m^3); atomic density is 3x this.
_SIO2_FORMULA_DENSITY_M3 = 2.2e28
# Fluorocarbon film atomic density, matching FLUOROCARBON_FILM.atom_density_m3.
_FC_FILM_ATOM_DENSITY_M3 = FLUOROCARBON_FILM.atom_density_m3
# One-monolayer areal density used for coverage saturation (atoms/m^2).
_MONOLAYER_AREAL_M2 = _FC_FILM_ATOM_DENSITY_M3 ** (2.0 / 3.0)


@dataclass(frozen=True)
class MixedLayerParams:
    """Chemical constants with named provenance (design doc section 5.4).

    The four fitted survivors (s_p, s_f, p_ox, eta_mix, k_v) each carry a
    literature anchor; everything energetic is derived from stopping curves.
    """

    ion_atomic_number: int = 18          # Ar+ projectile
    ion_mass_amu: float = 39.948
    precursor_fc_ratio: float = 1.5      # y: F per C in CxFy precursor (C4F6)
    sticking_probability: float = 0.0842  # s_p: Krueger polymer deposition
    fluorine_film_sticking: float = 0.05  # s_f: F absorption into film (Gogolides band)
    oxidation_probability: float = 0.0628  # p_ox: Krueger O-driven polymer etch
    mixing_efficiency: float = 1.0       # eta_mix: Humbird-Graves, O(1)
    volatilization_yield: float = 1.0    # k_v: SiF4 per ion at reference deposition
    substrate: str = "sio2"              # "sio2" | "carbon" (a-C mask, no lattice O)
    reference_energy_eV: float = 1000.0  # E_ref for the deposited-energy ratio
    film_sputter_yield: float = 0.1384   # film atoms per ion at reference deposition
    minimum_layer_depth_nm: float = 0.5
    clog_film_thickness_nm: float = 20.0


@dataclass(frozen=True)
class MixedLayerState:
    """Areal densities (atoms/m^2) of the two reservoirs."""

    n_c_film: float = 0.0
    n_f_film: float = 0.0
    n_si: float = 0.0
    n_o: float = 0.0
    n_c: float = 0.0
    n_f: float = 0.0

    def film_thickness_nm(self) -> float:
        return ((self.n_c_film + self.n_f_film) / _FC_FILM_ATOM_DENSITY_M3) * 1e9


@dataclass(frozen=True)
class SurfaceFluxes:
    """Incoming fluxes (per m^2 per s) and ion energetics."""

    precursor_flux: float
    fluorine_flux: float
    oxygen_flux: float
    ion_flux: float
    ion_energy_eV: float
    cosine_incidence: float = 1.0
    # Optional per-face F/C ratio of the carbon-carrying precursor flux;
    # None falls back to params.precursor_fc_ratio (single-species case).
    precursor_fc_ratio: object = None


@dataclass
class StepResult:
    state: MixedLayerState
    recession_velocity_m_s: float
    sif4_rate: float
    substrate_removal_rate: float
    co_rate: float
    cof2_rate: float
    interface_energy_eV: float
    deposited_energy_eV: float
    ledger_residuals: dict = field(default_factory=dict)


_RANGE_TABLE_CACHE: dict = {}
_TABLE_ENERGY_MIN_EV = 12.0
_TABLE_ENERGY_MAX_EV = 20000.0
_TABLE_POINTS = 96


def _stopping_tables(params: MixedLayerParams):
    """Log-energy interpolation tables for the three stopping quantities.

    Exact per-point values from the ZBL module; bilinear-free (1-D log-log
    interpolation, <0.3 percent versus direct evaluation across the table).
    """
    import math

    key = (params.ion_atomic_number, params.ion_mass_amu,
           params.minimum_layer_depth_nm)
    cached = _RANGE_TABLE_CACHE.get(key)
    if cached is not None:
        return cached
    n = _TABLE_POINTS
    log_lo = math.log(_TABLE_ENERGY_MIN_EV)
    log_hi = math.log(_TABLE_ENERGY_MAX_EV)
    energies = [math.exp(log_lo + (log_hi - log_lo) * i / (n - 1)) for i in range(n)]
    lam_fc = [projected_range_nm(e, params.ion_atomic_number, params.ion_mass_amu,
                                 FLUOROCARBON_FILM) for e in energies]
    depth = [max(params.minimum_layer_depth_nm,
                 projected_range_nm(e, params.ion_atomic_number,
                                    params.ion_mass_amu, SIO2))
             for e in energies]
    # Normal-incidence deposited energy at the tabulated depth; the cosine
    # dependence is applied by the caller through the slant-path argument.
    dep = [float(nuclear_energy_in_layer_eV(e, 1.0, d, params.ion_atomic_number,
                                            params.ion_mass_amu, SIO2))
           for e, d in zip(energies, depth)]
    tables = (log_lo, log_hi, n, lam_fc, depth, dep)
    _RANGE_TABLE_CACHE[key] = tables
    return tables


def _table_lookup(energy_eV, column):
    """Vectorized log-energy interpolation; linear roll-off below the table."""
    energy = np.asarray(energy_eV, dtype=float)
    column = np.asarray(column, dtype=float)
    log_lo = np.log(_TABLE_ENERGY_MIN_EV)
    log_hi = np.log(_TABLE_ENERGY_MAX_EV)
    n = column.shape[0]
    clipped = np.clip(energy, _TABLE_ENERGY_MIN_EV, _TABLE_ENERGY_MAX_EV)
    pos = (np.log(clipped) - log_lo) / (log_hi - log_lo) * (n - 1)
    i = np.minimum(pos.astype(int), n - 2)
    frac = pos - i
    interior = column[i] * (1.0 - frac) + column[i + 1] * frac
    low = column[0] * (energy / _TABLE_ENERGY_MIN_EV)
    return np.where(energy <= _TABLE_ENERGY_MIN_EV, low, interior)


def interface_energy_eV(ion_energy_eV, film_thickness_nm,
                        params: MixedLayerParams):
    """Standaert defluorination law: E_iface = E * exp(-d_FC / lambda_FC)."""
    energy = np.asarray(ion_energy_eV, dtype=float)
    thickness = np.asarray(film_thickness_nm, dtype=float)
    _, _, _, lam_fc, _, _ = _stopping_tables(params)
    lam = _table_lookup(np.maximum(energy, _TABLE_ENERGY_MIN_EV * 1e-3), lam_fc)
    safe_lam = np.maximum(lam, 1e-300)
    result = energy * np.exp(-np.minimum(thickness / safe_lam, 700.0))
    return np.where((energy > 0.0) & (lam > 0.0), result, 0.0)


def _deposited_energy(e_iface, cosine, params: MixedLayerParams):
    """Nuclear energy deposited in the mixed layer; layer depth from range."""
    e_iface = np.asarray(e_iface, dtype=float)
    cosine = np.asarray(cosine, dtype=float)
    _, _, _, _, depth_col, dep_col = _stopping_tables(params)
    depth = _table_lookup(np.maximum(e_iface, _TABLE_ENERGY_MIN_EV * 1e-3),
                          depth_col)
    energy = _table_lookup(np.maximum(e_iface, _TABLE_ENERGY_MIN_EV * 1e-3),
                           dep_col)
    # Slant path deposits more of the ion's energy in the layer, capped at E.
    slant = np.minimum(e_iface, energy / np.maximum(cosine, 0.05))
    energy = np.where((cosine > 0.0) & (cosine < 1.0), slant, energy)
    energy = np.where(e_iface > 0.0, energy, 0.0)
    depth = np.where(e_iface > 0.0, depth, params.minimum_layer_depth_nm)
    return energy, depth


def step(state: MixedLayerState, fluxes: SurfaceFluxes, dt: float,
         params: MixedLayerParams = MixedLayerParams()) -> StepResult:
    """Advance one conservative explicit step.

    Removal terms that would overdraw a reservoir within dt are scaled down
    proportionally, so densities stay non-negative and the ledgers close
    exactly (up to float rounding, reported in ledger_residuals).
    """
    d_fc = state.film_thickness_nm()
    theta_film = 1.0 - _exp_neg(d_fc * 1e-9 * _FC_FILM_ATOM_DENSITY_M3
                                / _MONOLAYER_AREAL_M2)
    e_iface = interface_energy_eV(fluxes.ion_energy_eV, d_fc, params)
    eps_dep, _depth = _deposited_energy(e_iface, fluxes.cosine_incidence, params)
    energy_ratio = eps_dep / params.reference_energy_eV

    film_total = np.asarray(state.n_c_film + state.n_f_film, dtype=float)
    x_c = _guarded_ratio(state.n_c_film, film_total, 1.0)
    x_f = _guarded_ratio(state.n_f_film, film_total, 1.0)

    # --- film gains ---
    fc_ratio = (params.precursor_fc_ratio if fluxes.precursor_fc_ratio is None
                else fluxes.precursor_fc_ratio)
    dep_c = params.sticking_probability * fluxes.precursor_flux
    dep_f = dep_c * fc_ratio
    absorb_f = params.fluorine_film_sticking * fluxes.fluorine_flux * theta_film

    # --- film losses (proposed rates, atoms/m^2/s) ---
    sputter_total = params.film_sputter_yield * fluxes.ion_flux * energy_ratio * theta_film
    sput_c = sputter_total * x_c
    sput_f = sputter_total * x_f
    # O oxidation of film carbon: each oxidized C carries along the local film
    # F/C ratio (capped at 2, the COF2 stoichiometry); remainder leaves as CO.
    ox_c = params.oxidation_probability * fluxes.oxygen_flux * theta_film * x_c
    f_per_ox_c = _guarded_ratio(state.n_f_film, state.n_c_film, 2.0)
    ox_f = ox_c * f_per_ox_c
    # Ion-driven mixing of film content into the layer (Humbird-Graves).
    mix_total = params.mixing_efficiency * fluxes.ion_flux * energy_ratio * theta_film
    mix_c = mix_total * x_c
    mix_f = mix_total * x_f

    # Clamp film losses to availability.
    loss_c = sput_c + ox_c + mix_c
    loss_f = sput_f + ox_f + mix_f
    scale_c = _overdraw_scale(state.n_c_film + dep_c * dt, loss_c, dt)
    scale_f = _overdraw_scale(state.n_f_film + (dep_f + absorb_f) * dt, loss_f, dt)
    scale_film = np.minimum(scale_c, scale_f)  # keep C/F branches consistent
    sput_c = sput_c * scale_film; ox_c = ox_c * scale_film; mix_c = mix_c * scale_film
    sput_f = sput_f * scale_film; ox_f = ox_f * scale_film; mix_f = mix_f * scale_film

    # --- mixed layer ---
    # Layer fluorine coverage (site fraction of a saturated monolayer). The
    # volatilization rate is coverage-proportional — the same Langmuir kinetics
    # as the Belen/ViennaPS coupled-coverage model (rate = capacity * theta_F),
    # so the degenerate no-carbon limit reproduces that structure exactly
    # (Rung 0) instead of a sharp supply/capacity switch.
    theta_f_layer = np.minimum(
        np.asarray(state.n_f, dtype=float) / _MONOLAYER_AREAL_M2, 1.0)
    # Direct F where the film is open, Langmuir (1 - theta) sticking; the
    # reflected remainder never enters the ledger.
    f_direct = (fluxes.fluorine_flux * (1.0 - theta_film)
                * (1.0 - theta_f_layer))
    # Ion capacity for substrate volatilization (derived energy factor; no
    # fitted law). SiO2 leaves as SiF4 (4 F per Si); an amorphous-carbon mask
    # leaves as CFx (2 F per C, the F-costly channel) with no lattice oxygen —
    # selectivity must emerge from the film/energy/F budgets, never a parameter.
    volat_capacity = (params.volatilization_yield * fluxes.ion_flux
                      * energy_ratio * np.maximum(1.0 - theta_film, 0.0))
    if params.substrate == "carbon":
        f_per_unit = 2.0
        sif4 = np.zeros_like(volat_capacity * theta_f_layer)
        substrate_removal = volat_capacity * theta_f_layer
    else:
        f_per_unit = 4.0
        sif4 = volat_capacity * theta_f_layer
        substrate_removal = sif4
    # Layer oxidation of mixed C by layer O (same probability channel);
    # clamp to available layer carbon first.
    layer_ox = params.oxidation_probability * fluxes.oxygen_flux * (1.0 - theta_film) * (
        state.n_c / (state.n_c + _MONOLAYER_AREAL_M2))
    layer_ox = layer_ox * _overdraw_scale(state.n_c + mix_c * dt, layer_ox, dt)
    f_per_layer_c = _guarded_ratio(state.n_f, state.n_c, 2.0)
    # COF2 branch consumes layer F; the CO branch does not — an F-starved
    # clamp below shifts oxidized carbon from COF2 to CO rather than dropping it.
    cof2 = layer_ox * (f_per_layer_c / 2.0)

    # Clamp the fluorine-consuming channels to layer F availability.
    layer_f_loss = f_per_unit * substrate_removal + 2.0 * cof2
    scale_lf = _overdraw_scale(state.n_f + (mix_f + f_direct) * dt, layer_f_loss, dt)
    substrate_removal = substrate_removal * scale_lf
    sif4 = sif4 * scale_lf
    cof2 = cof2 * scale_lf

    # Recession; for SiO2 it liberates lattice Si and O into the layer
    # (1 Si : 2 O per SiF4). Carbon lattice leaves directly as CFx product.
    if params.substrate == "carbon":
        recession = substrate_removal / 1.0e29  # a-C atomic density
        si_in = 0.0
        o_in = 0.0
    else:
        recession = sif4 / _SIO2_FORMULA_DENSITY_M3  # m/s
        si_in = sif4
        o_in = 2.0 * sif4

    # Interfacial oxidation (the Standaert/Oehrlein selectivity mechanism):
    # lattice oxygen liberated by recession attacks the film carbon from below
    # as CO, carrying film F along as COF2. Efficiency is fixed by the etch
    # complex stoichiometry (SiO2-C: two lattice O consume one film C per
    # etched Si) — a chemistry ratio, not a rate constant. This is the term
    # that keeps the film thin on SiO2 while a carbon substrate (no lattice O)
    # grows a thick protective film.
    rem_c_film = state.n_c_film + dt * (dep_c - sput_c - ox_c - mix_c)
    rem_f_film = state.n_f_film + dt * (dep_f + absorb_f - sput_f - ox_f - mix_f)
    bottom_c = np.minimum(0.5 * o_in,
                          np.maximum(rem_c_film, 0.0) / max(dt, 1e-30))
    f_carry = _guarded_ratio(rem_f_film, rem_c_film, 2.0)
    bottom_f = np.minimum(bottom_c * f_carry,
                          np.maximum(rem_f_film, 0.0) / max(dt, 1e-30))
    o_to_layer = o_in - bottom_c

    # Layer-side product oxygen (one O per oxidized layer C, CO or COF2 alike);
    # film-side CO/COF2 oxygen comes from the incident O flux directly. Clamp to
    # layer O availability, shifting the shortfall out of both product branches.
    layer_side_o = layer_ox
    scale_lo = _overdraw_scale(state.n_o + o_to_layer * dt, layer_side_o, dt)
    layer_ox = layer_ox * scale_lo
    cof2 = cof2 * scale_lo
    layer_side_o = layer_ox
    layer_f_loss = f_per_unit * substrate_removal + 2.0 * cof2

    co = (ox_c - 0.5 * ox_f) + (layer_ox - cof2) + (bottom_c - 0.5 * bottom_f)
    cof2_total = cof2 + 0.5 * ox_f + 0.5 * bottom_f

    # Excess layer oxygen beyond a saturated monolayer desorbs recombinatively
    # (O2) — surfaces cannot hold more than saturation coverage. Accounted as
    # outflow so the O ledger still closes exactly.
    n_o_raw = state.n_o + dt * (o_to_layer - layer_side_o)
    o_desorb = np.maximum(n_o_raw - _MONOLAYER_AREAL_M2, 0.0) / max(dt, 1e-30)
    # Layer fluorine saturates at coverage too; the excess recombines/reflects.
    n_f_raw = state.n_f + dt * (mix_f + f_direct - layer_f_loss)
    f_desorb = np.maximum(n_f_raw - _MONOLAYER_AREAL_M2, 0.0) / max(dt, 1e-30)

    new_state = MixedLayerState(
        n_c_film=state.n_c_film + dt * (dep_c - sput_c - ox_c - mix_c - bottom_c),
        n_f_film=state.n_f_film + dt * (dep_f + absorb_f - sput_f - ox_f - mix_f
                                        - bottom_f),
        n_si=state.n_si + dt * (si_in - sif4),
        n_o=n_o_raw - dt * o_desorb,
        n_c=state.n_c + dt * (mix_c - layer_ox),
        n_f=n_f_raw - dt * f_desorb,
    )

    residuals = _ledger_residuals(state, new_state, dt, dep_c, dep_f, absorb_f,
                                  f_direct, si_in, o_in, sput_c, sput_f,
                                  ox_c + bottom_c, ox_f + bottom_f, mix_c, mix_f,
                                  sif4, layer_ox,
                                  layer_side_o + o_desorb + bottom_c,
                                  layer_f_loss + f_desorb)
    return StepResult(
        state=new_state,
        recession_velocity_m_s=recession,
        sif4_rate=sif4,
        substrate_removal_rate=substrate_removal,
        co_rate=co,
        cof2_rate=cof2_total,
        interface_energy_eV=e_iface,
        deposited_energy_eV=eps_dep,
        ledger_residuals=residuals,
    )


def steady_state(fluxes: SurfaceFluxes, params: MixedLayerParams = MixedLayerParams(),
                 *, dt: float = 1e-4, max_steps: int = 200000,
                 relative_tolerance: float = 1e-10) -> StepResult:
    """Integrate to steady state or report clog (unbounded film growth)."""
    state = MixedLayerState()
    # Explicit-integrator stability bound: the fastest reservoir timescale is
    # a monolayer against the largest total flux. Beyond it the map goes
    # oscillatory and settles off the true fixed point, so the acceleration
    # ramp is capped here — declarations below the bound are trustworthy.
    fastest_flux = max(
        float(np.max(np.asarray(fluxes.fluorine_flux, dtype=float))),
        float(np.max(np.asarray(fluxes.precursor_flux, dtype=float)))
        * (1.0 + params.precursor_fc_ratio),
        float(np.max(np.asarray(fluxes.oxygen_flux, dtype=float))),
        4.0 * params.volatilization_yield
        * float(np.max(np.asarray(fluxes.ion_flux, dtype=float))),
        1.0)
    dt_stable = 0.25 * _MONOLAYER_AREAL_M2 / fastest_flux
    base_dt = dt
    result = step(state, fluxes, dt, params)
    for _ in range(max_steps):
        nxt = step(result.state, fluxes, dt, params)
        thickness = nxt.state.film_thickness_nm()
        growing = thickness > result.state.film_thickness_nm()
        # Trajectory-faithful adaptive ramp: the system is bistable, so which
        # basin a transient lands in depends on integration accuracy. Grow the
        # step only while per-step relative change stays small; shrink it when
        # the dynamics are fast. Always within the explicit stability bound.
        change = 0.0
        for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f"):
            va = np.asarray(getattr(result.state, name), dtype=float)
            vb = np.asarray(getattr(nxt.state, name), dtype=float)
            scale = np.maximum(np.maximum(np.abs(va), np.abs(vb)),
                               _MONOLAYER_AREAL_M2)
            change = max(change, float(np.max(np.abs(vb - va) / scale)))
        if change > 0.02:
            dt = max(dt * 0.5, base_dt)
        else:
            dt = min(dt * 1.02, max(dt_stable, base_dt))
        if thickness > params.clog_film_thickness_nm:
            nxt.ledger_residuals["clogged"] = True
            return nxt
        # Terminal clog: the interface energy has collapsed, so every
        # energy-driven removal channel is off, yet the film still grows —
        # growth can no longer reverse (oxidation is already saturated).
        if (growing and fluxes.ion_energy_eV > 0.0
                and nxt.interface_energy_eV < 0.02 * fluxes.ion_energy_eV):
            nxt.ledger_residuals["clogged"] = True
            return nxt
        if _converged(result.state, nxt.state, dt, relative_tolerance):
            nxt.ledger_residuals["clogged"] = False
            return nxt
        result = nxt
    result.ledger_residuals["clogged"] = False
    result.ledger_residuals["converged"] = False
    return result


def _exp_neg(x):
    return np.exp(-np.minimum(x, 700.0))


def _guarded_ratio(numerator, denominator, cap):
    """min(cap, num/den) where den > 0, else 0 — overflow-safe on the discarded branch."""
    den = np.asarray(denominator, dtype=float)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        ratio = np.asarray(numerator, dtype=float) / np.maximum(den, 1e-300)
    return np.where(den > 0.0, np.minimum(cap, ratio), 0.0)


def _overdraw_scale(available_plus_gain, loss_rate, dt):
    drawable = np.maximum(available_plus_gain, 0.0)
    needed = np.asarray(loss_rate, dtype=float) * dt
    safe = np.maximum(needed, 1e-300)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        ratio = drawable / safe
    return np.where(needed <= drawable, 1.0,
                    np.where(needed > 0.0, ratio, 1.0))


def _converged(a: MixedLayerState, b: MixedLayerState, dt: float, tol: float) -> bool:
    for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f"):
        va, vb = getattr(a, name), getattr(b, name)
        scale = max(abs(va), abs(vb), _MONOLAYER_AREAL_M2 * 1e-6)
        if abs(vb - va) > tol * scale * max(dt, 1e-30) / dt:
            return False
    return True


def _ledger_residuals(old, new, dt, dep_c, dep_f, absorb_f, f_direct, si_in,
                      o_in, sput_c, sput_f, ox_c, ox_f, mix_c, mix_f, sif4,
                      layer_ox, layer_side_o, layer_f_loss):
    """Per-element: inflow - outflow - d(storage)/dt must vanish."""
    store = lambda n: (getattr(new, n) - getattr(old, n)) / dt
    f_in = dep_f + absorb_f + f_direct
    f_out = sput_f + ox_f + layer_f_loss  # mix_f is internal transfer
    f_res = f_in - f_out - (store("n_f_film") + store("n_f"))
    c_in = dep_c
    c_out = sput_c + ox_c + layer_ox
    c_res = c_in - c_out - (store("n_c_film") + store("n_c"))
    si_res = si_in - sif4 - store("n_si")
    o_res = o_in - layer_side_o - store("n_o")
    norm = np.maximum(
        np.maximum(np.asarray(f_in, dtype=float), np.asarray(c_in, dtype=float)),
        np.maximum(np.asarray(si_in, dtype=float), np.asarray(o_in, dtype=float)))
    norm = np.maximum(norm, _MONOLAYER_AREAL_M2)
    return {
        "fluorine": f_res / norm,
        "carbon": c_res / norm,
        "silicon": si_res / norm,
        "oxygen": o_res / norm,
    }
