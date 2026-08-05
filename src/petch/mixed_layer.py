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
# Krueger oxide site density (site_density_m2 in the reduced projection).
_SITE_DENSITY_M2 = 1.0e19


@dataclass(frozen=True)
class MixedLayerParams:
    """Chemical constants with named provenance (design doc section 5.4).

    The four fitted survivors (s_p, s_f, p_ox, eta_mix, k_v) each carry a
    literature anchor; everything energetic is derived from stopping curves.
    """

    displacement_energy_eV: float = 25.0  # polymer displacement energy (literature band 10-80)
    ion_atomic_number: int = 18          # Ar+ projectile
    ion_mass_amu: float = 39.948
    precursor_fc_ratio: float = 1.5      # y: F per C in CxFy precursor (C4F6)
    sticking_probability: float = 0.0842  # s_p: Krueger polymer deposition
    fluorine_film_sticking: float = 0.05  # s_f: F absorption into film (Gogolides band)
    oxidation_probability: float = 0.0423  # p_ox: Table-6.5 converged O polymer etch
    mixing_efficiency: float = 1.0       # eta_mix: Humbird-Graves, O(1)
    volatilization_yield: float = 1.0    # multiplier on published yields (1.0 = no anchor)
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
    # Crosslinked film atoms (subset of n_c_film + n_f_film): ion-processed
    # skin with reduced radical attachment (Krueger Appendix B: 0.02).
    n_xl_film: float = 0.0
    # Ion-activated oxide sites (Appendix A1: SiO2+ion->SiO2* at 0.9/ion);
    # chemisorption on activated sites runs at 0.8-0.9 vs 0.278 bare (CH2).
    # Mass-neutral surface state, capped at the site density.
    n_act: float = 0.0

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
    # Direct chemisorption into the mixed layer where the film is open
    # (Krueger complex-formation channel; per-species published
    # probabilities are applied by the adapter). Atoms/m^2/s of C and F.
    chemisorption_carbon_flux: object = 0.0
    chemisorption_fluorine_flux: object = 0.0
    # Substrate-dependent deposition (Krueger: sticking on polymer ~0.1 vs
    # on bare substrate ~0.001-0.002). When provided, these override the
    # single sticking_probability path: deposition blends with theta_film.
    film_deposition_carbon_flux: object = None
    film_deposition_fluorine_flux: object = None
    substrate_deposition_carbon_flux: object = None
    substrate_deposition_fluorine_flux: object = None
    crosslinked_deposition_carbon_flux: object = None
    crosslinked_deposition_fluorine_flux: object = None
    # Chemisorption aggregates at the ACTIVATED-site probabilities (CH2:
    # 0.8/0.85/0.9); None disables the two-state blend.
    chemisorption_activated_carbon_flux: object = None
    chemisorption_activated_fluorine_flux: object = None
    # Deposition-weighted maximum crosslink partners per arriving radical
    # (Krueger 2024: CF 3, CF2 2, CF3 1).  None falls back to the published
    # row stoichiometry alone (one partner).
    deposition_available_bonds: object = None
    # Atom-resolved ion spectrum (per-event): sparse (face, flux, E, cos)
    # arrays. When present, EVERY ion-driven term is evaluated per atom
    # against the live face state and segment-summed (research doc
    # RESEARCH_EVENT_RESOLVED_CHEMISTRY: clamps stay at face level).
    ion_atom_face: object = None
    ion_atom_flux: object = None
    ion_atom_energy_eV: object = None
    ion_atom_cosine: object = None


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
    film_deposition_rate: float = 0.0
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


def _threshold_power_yield(energy_eV, p0, threshold_eV, reference_eV, exponent):
    """Krueger Eq-2 threshold power law: p0*(E^q - Eth^q)/(E0^q - Eth^q), >= 0."""
    energy = np.asarray(energy_eV, dtype=float)
    numerator = np.maximum(energy, 0.0) ** exponent - threshold_eV ** exponent
    denominator = reference_eV ** exponent - threshold_eV ** exponent
    return p0 * np.maximum(numerator, 0.0) / denominator


# --- Appendix-B angular classes (Krueger Table B.0.1 column 5) -----------------
# The legend (thesis L5332) reads: "p0 is modified according to Eq. (2.40) if
# angular or energy dependence of the reaction is present ... and (angle) defines
# the nature of the angular dependence, with 1 corresponding to the results
# obtained by [1] and 2 corresponding to the results obtained by [2]", with
# [1] = Kress et al., JVST A 17, 2819 (1999) and
# [2] = Chang & Sawin, JVST A 15, 610 (1997).
#
# The two shapes are stated verbatim and identically across three theses of the
# same MCFPM lineage:
#   Huang (Eq. 2.32 discussion): "For physical sputtering, f(theta) is an
#     empirical function with a maximum at 60 deg, reduced probability at normal
#     incidence and zero probability at grazing incidence. For chemically
#     enhanced etching, f(theta) is unity for normal incidence and angles up to
#     45 deg, with a monotonic roll-off to zero probability at grazing
#     incidence."
#   Huard (Eq. 2.24): "the total yield of a sputtering reaction is given by
#     P(e,q) = P(e) P(q)" -- energy and angle factors are separable, which is
#     why the angular class multiplies our energy term rather than replacing it.
#   Qu (Eq. 2.40 discussion): same two classes.
#
# Normalisation convention: f(0) = 1 for BOTH classes, so p0 in the table is the
# normal-incidence probability at the reference energy.  This is the convention
# the already-validated polymer row (CF(s)+Ar+, class 1) has used since the
# audit-corrected campaign, and mixing conventions across rows of one table
# would be incoherent.  (The sources describe class 1 as "less than unity at
# normal incidence" relative to its 60-deg peak; that is a statement about the
# peak/normal ratio, which this form reproduces as 4.17, not about p0.)


# Class-1 shape parameter B in (1 + B sin^2 t) cos t.
#
# KRUEGER_B is the value his cited source implies and the one the polymer row
# (CF(s)+Ar+) has carried through every validated lip/mouth result.
#
# OXIDE_B is the value bounded by the only angular sputter measurements taken
# in THIS chemistry on THIS material:
#   Cho et al., JVST A 18, 2705 (2000): SiO2 in CF4, peak/normal ~1.3 at all
#     biases;  Schaepkens et al., JVST A 16, 3281 (1998): 54.7-deg V-groove,
#     peak/normal ~1.33.
# Kress et al. (1999), the citation behind class 1, is a molecular-dynamics
# study of "Cu and Ar ion sputtering of Cu(111) surfaces" at 50-250 eV --
# a different material system, and its shape gives peak/normal 4.17, i.e.
# 3.2x above every in-chemistry measurement (RESEARCH_VERIFY_HUNT_2026-08-05).
#
# Landing OXIDE_B on the oxide/mask rows is a declared model choice between two
# sourced options, made on chemistry match.  It is not fitted: f(0) = 1 for any
# B, so every normal-incidence and blanket result is bitwise unchanged, and the
# resulting peak (1.31) sits inside the measured band [1.30, 1.33].
# The polymer row keeps KRUEGER_B; the measured FC-film counterpart
# (Barklund & Blom, JVST A 10, 1212 (1992), Ar+ peak 1.448 at 65 deg) is a
# flagged follow-up that needs its own graded run against the lip results.
_KRUEGER_CLASS1_B = 9.3
_OXIDE_CLASS1_B = 1.7


def _class1_shape(cosine, b):
    cos_t = np.clip(np.asarray(cosine, dtype=float), 0.0, 1.0)
    return np.maximum((1.0 + b * (1.0 - cos_t ** 2)) * cos_t, 0.0)


def _angular_physical_sputter(cosine):
    """Class 1 on the polymer row (Krueger's cited Kress form), f(0)=1."""
    return _class1_shape(cosine, _KRUEGER_CLASS1_B)


def _angular_oxide_sputter(cosine):
    """Class 1 on the oxide/mask rows, bounded by in-chemistry measurement."""
    return _class1_shape(cosine, _OXIDE_CLASS1_B)


def _angular_chemical_sputter(cosine):
    """Class 2 (Chang & Sawin 1997): unity to 45 deg, monotone roll-off to zero.

    The plateau edge (45 deg), both endpoints (unity at normal, zero at grazing)
    and monotonicity between them are verbatim from the source lineage.  The
    interpolation across the roll-off is the minimal projected-flux choice
    cos(theta)/cos(45 deg), which introduces no constant beyond the stated
    plateau edge.  The published curve itself is paywalled; the specific
    roll-off shape is therefore recorded [VERIFY] while its stated properties
    are gated.
    """
    cos_t = np.clip(np.asarray(cosine, dtype=float), 0.0, 1.0)
    return np.minimum(cos_t / _COS45, 1.0)


_COS45 = float(np.cos(np.pi / 4.0))


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
    ref_dep_140, _ = _deposited_energy(np.asarray(140.0), np.asarray(1.0), params)
    ref_dep_140 = np.maximum(float(ref_dep_140), 1e-300)
    if fluxes.ion_atom_face is not None:
        atom_face = np.asarray(fluxes.ion_atom_face, dtype=int)
        atom_flux = np.asarray(fluxes.ion_atom_flux, dtype=float)
        atom_energy = np.asarray(fluxes.ion_atom_energy_eV, dtype=float)
        atom_cos = np.clip(np.asarray(fluxes.ion_atom_cosine, dtype=float), 0.0, 1.0)
        shape = np.broadcast(np.asarray(d_fc), np.asarray(state.n_f)).shape
        d_fc_arr = np.broadcast_to(np.asarray(d_fc, dtype=float), shape)
        atom_e_iface = np.asarray(interface_energy_eV(
            atom_energy, d_fc_arr.ravel()[atom_face] if shape else d_fc_arr,
            params))
        atom_eps, _ = _deposited_energy(atom_e_iface, atom_cos, params)
        atom_kress = _angular_physical_sputter(atom_cos)
        atom_chem_ang = _angular_chemical_sputter(atom_cos)

        def _segment(values):
            # bincount == add.at bitwise (both accumulate in input order)
            # and ~28x faster on the segment sums that dominate step cost.
            size = int(np.prod(shape)) if shape else 1
            out = np.bincount(
                atom_face, weights=np.broadcast_to(values, atom_face.shape),
                minlength=size)
            return out.reshape(shape) if shape else out[0]

        kernel_sputter = _segment(
            atom_flux * np.asarray(_threshold_power_yield(
                atom_energy, 0.9, 20.0, 500.0, 0.5)) * atom_kress)
        kernel_mix = _segment(atom_flux * atom_eps
                              / params.reference_energy_eV)
        kernel_xl = _segment(
            atom_flux * np.maximum(atom_energy - atom_e_iface, 0.0))
        # Appendix-B angular classes: SiO2CF(s)+Ar+ is class 2, SiO2(s)+Ar+ and
        # AC(s)+Ar+ are class 1 (RESULTS_LIP_REMOVAL_AUDIT_2026-08-04).
        kernel_complex = _segment(
            0.1471 * atom_flux * atom_eps / ref_dep_140 * atom_chem_ang)
        atom_kress_ox = _angular_oxide_sputter(atom_cos)
        kernel_bare = _segment(atom_flux * np.asarray(_threshold_power_yield(
            atom_e_iface, 0.0852, 70.0, 140.0, 1.0)) * atom_kress_ox)
        kernel_ac = _segment(atom_flux * np.asarray(_threshold_power_yield(
            atom_e_iface, 0.001, 200.0, 250.0, 0.4)) * atom_kress_ox)
        kernel_dexl = _segment(atom_flux * np.asarray(_threshold_power_yield(
            atom_energy, 0.3, 8.0, 500.0, 0.5)))
        kernel_act = _segment(atom_flux) * 0.9
        # Flux-weighted diagnostics for receipts and clog logic.
        total_atom_flux = _segment(atom_flux)
        safe_total = np.maximum(total_atom_flux, 1e-300)
        e_iface = _segment(atom_flux * atom_e_iface) / safe_total
        eps_dep = _segment(atom_flux * atom_eps) / safe_total
        energy_ratio = eps_dep / params.reference_energy_eV
    else:
        e_iface = interface_energy_eV(fluxes.ion_energy_eV, d_fc, params)
        eps_dep, _depth = _deposited_energy(
            e_iface, fluxes.cosine_incidence, params)
        energy_ratio = eps_dep / params.reference_energy_eV
        cos_scalar = np.clip(
            np.asarray(fluxes.cosine_incidence, dtype=float), 0.0, 1.0)
        kress_scalar = _angular_physical_sputter(cos_scalar)
        chem_ang_scalar = _angular_chemical_sputter(cos_scalar)
        kernel_sputter = (np.asarray(_threshold_power_yield(
            fluxes.ion_energy_eV, 0.9, 20.0, 500.0, 0.5))
            * kress_scalar * fluxes.ion_flux)
        kernel_mix = fluxes.ion_flux * energy_ratio
        kernel_xl = fluxes.ion_flux * np.maximum(
            np.asarray(fluxes.ion_energy_eV, dtype=float) - e_iface, 0.0)
        kernel_complex = (0.1471 * fluxes.ion_flux * np.asarray(
            eps_dep, dtype=float) / ref_dep_140 * chem_ang_scalar)
        kress_ox_scalar = _angular_oxide_sputter(cos_scalar)
        kernel_bare = fluxes.ion_flux * np.asarray(_threshold_power_yield(
            e_iface, 0.0852, 70.0, 140.0, 1.0)) * kress_ox_scalar
        kernel_ac = fluxes.ion_flux * np.asarray(_threshold_power_yield(
            e_iface, 0.001, 200.0, 250.0, 0.4)) * kress_ox_scalar
        kernel_dexl = fluxes.ion_flux * np.asarray(_threshold_power_yield(
            fluxes.ion_energy_eV, 0.3, 8.0, 500.0, 0.5))
        kernel_act = fluxes.ion_flux * 0.9

    film_total = np.asarray(state.n_c_film + state.n_f_film, dtype=float)
    x_c = _guarded_ratio(state.n_c_film, film_total, 1.0)
    x_f = _guarded_ratio(state.n_f_film, film_total, 1.0)

    # --- film gains ---
    if fluxes.film_deposition_carbon_flux is not None:
        # Substrate-dependent sticking (published on-polymer vs on-substrate
        # probabilities): deposition rate blends with film coverage, and the
        # film-covered share further blends fresh vs crosslinked attachment
        # by the crosslinked fraction of the film (ion-processed skin).
        x_frac = _guarded_ratio(state.n_xl_film, film_total, 1.0)
        film_dep_c = np.asarray(fluxes.film_deposition_carbon_flux, dtype=float)
        film_dep_f = np.asarray(fluxes.film_deposition_fluorine_flux, dtype=float)
        if fluxes.crosslinked_deposition_carbon_flux is not None:
            film_dep_c = ((1.0 - x_frac) * film_dep_c
                          + x_frac * np.asarray(
                              fluxes.crosslinked_deposition_carbon_flux,
                              dtype=float))
            film_dep_f = ((1.0 - x_frac) * film_dep_f
                          + x_frac * np.asarray(
                              fluxes.crosslinked_deposition_fluorine_flux,
                              dtype=float))
        dep_c = (theta_film * film_dep_c
                 + (1.0 - theta_film)
                 * np.asarray(fluxes.substrate_deposition_carbon_flux,
                              dtype=float))
        dep_f = (theta_film * film_dep_f
                 + (1.0 - theta_film)
                 * np.asarray(fluxes.substrate_deposition_fluorine_flux,
                              dtype=float))
    else:
        fc_ratio = (params.precursor_fc_ratio
                    if fluxes.precursor_fc_ratio is None
                    else fluxes.precursor_fc_ratio)
        dep_c = params.sticking_probability * fluxes.precursor_flux
        dep_f = dep_c * fc_ratio
    absorb_f = params.fluorine_film_sticking * fluxes.fluorine_flux * theta_film

    # --- film losses (proposed rates, atoms/m^2/s) ---
    # Film sputter: Krueger's PUBLISHED polymer sputter law (p0=0.9,
    # eth=20 eV, q=0.5, e0=500 eV; Kress B=9.3) at the INCIDENT ion energy,
    # evaluated per atom when the spectrum is atom-resolved.
    sputter_total = kernel_sputter * theta_film
    sput_c = sputter_total * x_c
    sput_f = sputter_total * x_f
    # O oxidation of the film. Krueger's published row is per collision with
    # EXPOSED POLYMER -- O(g) + P(s) -> products at p_ox -- so the reacting
    # fraction is the film coverage theta_film alone, and each reaction removes
    # one polymer unit: one carbon plus the local film F/C ratio (capped at 2,
    # the COF2 stoichiometry; the remainder leaves as CO). Scaling the carbon by
    # the film composition x_c as well double-counts the composition: p_ox is a
    # per-cell probability that already subsumes which atom the O lands on. That
    # spurious factor throttled the channel by 1/x_c = 2.7x at the film
    # composition the evolution runs actually reach (F/C = 1.69), which is why
    # a +48% change in p_ox moved the neck only +11%.
    # See RESULTS_O_CHANNEL_2026-08-04.md.
    ox_c = params.oxidation_probability * fluxes.oxygen_flux * theta_film
    f_per_ox_c = _guarded_ratio(state.n_f_film, state.n_c_film, 2.0)
    ox_f = ox_c * f_per_ox_c
    # Ion-driven mixing of film content into the layer (Humbird-Graves).
    mix_total = params.mixing_efficiency * kernel_mix * theta_film
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
    # Direct precursor chemisorption into the layer where the film is open
    # (the Krueger complex-formation channel): fluorine arrives bound to
    # carbon and both enter the layer ledgers, site-limited like f_direct.
    # Two-state oxide (Appendix A1/CH2): activated sites chemisorb ~3x the
    # bare rate; the blend follows the activated fraction theta_act.
    theta_act = np.minimum(
        np.asarray(state.n_act, dtype=float) / _SITE_DENSITY_M2, 1.0)
    site_open = (1.0 - theta_film) * (1.0 - theta_f_layer)
    if fluxes.chemisorption_activated_carbon_flux is not None:
        chem_c = ((1.0 - theta_act)
                  * np.asarray(fluxes.chemisorption_carbon_flux, dtype=float)
                  + theta_act * np.asarray(
                      fluxes.chemisorption_activated_carbon_flux, dtype=float)
                  ) * site_open
        chem_f = ((1.0 - theta_act)
                  * np.asarray(fluxes.chemisorption_fluorine_flux, dtype=float)
                  + theta_act * np.asarray(
                      fluxes.chemisorption_activated_fluorine_flux, dtype=float)
                  ) * site_open
    else:
        chem_c = (np.asarray(fluxes.chemisorption_carbon_flux, dtype=float)
                  * site_open)
        chem_f = (np.asarray(fluxes.chemisorption_fluorine_flux, dtype=float)
                  * site_open)
    # Two-state oxide removal at PUBLISHED magnitudes, anchored to the
    # DEKNOB-validated ZBL deposited-energy shape at Krueger's reference
    # energy (140 eV): complex channel 0.1384 (F-covered sites, costs 4F as
    # SiF4) + bare-SiO2 physical sputter 0.0909 @ 70 eV threshold (no F
    # cost — formula units leave as ejecta). volatilization_yield is a
    # multiplier defaulting to 1.0: published magnitude, no anchor constant.
    open_area = np.maximum(1.0 - theta_film, 0.0)
    if params.substrate == "carbon":
        # Krueger AC mask: pure physical sputter, 0.001 @ 200 eV threshold
        # (q=0.4, e0=250), NO chemical F channel, O-inert (1e-5, declared
        # omitted). The mask erodes only where the film exposes it.
        f_per_unit = 2.0
        sif4 = np.zeros_like(kernel_ac * open_area)
        substrate_removal = (params.volatilization_yield * kernel_ac
                             * open_area)
        bare_removal = np.zeros_like(substrate_removal)
        f_costed_removal = np.zeros_like(substrate_removal)
    else:
        f_per_unit = 4.0
        capacity = (params.volatilization_yield * kernel_complex * open_area)
        sif4 = capacity * theta_f_layer
        bare_removal = (params.volatilization_yield * kernel_bare * open_area
                        * np.maximum(1.0 - theta_f_layer, 0.0))
        substrate_removal = sif4 + bare_removal
        f_costed_removal = sif4
    # Layer oxidation of mixed C by layer O (same probability channel);
    # clamp to available layer carbon first.
    layer_ox = params.oxidation_probability * fluxes.oxygen_flux * (1.0 - theta_film) * (
        state.n_c / (state.n_c + _MONOLAYER_AREAL_M2))
    layer_ox = layer_ox * _overdraw_scale(
        state.n_c + (mix_c + chem_c) * dt, layer_ox, dt)
    f_per_layer_c = _guarded_ratio(state.n_f, state.n_c, 2.0)
    # COF2 branch consumes layer F; the CO branch does not — an F-starved
    # clamp below shifts oxidized carbon from COF2 to CO rather than dropping it.
    cof2 = layer_ox * (f_per_layer_c / 2.0)

    # Clamp the fluorine-consuming channels to layer F availability.
    layer_f_loss = f_per_unit * f_costed_removal + 2.0 * cof2
    scale_lf = _overdraw_scale(
        state.n_f + (mix_f + f_direct + chem_f) * dt, layer_f_loss, dt)
    f_costed_removal = f_costed_removal * scale_lf
    sif4 = sif4 * scale_lf
    cof2 = cof2 * scale_lf
    substrate_removal = f_costed_removal + (
        bare_removal if params.substrate != "carbon" else 0.0)

    # Recession; for SiO2 it liberates lattice Si and O into the layer
    # (1 Si : 2 O per SiF4). Carbon lattice leaves directly as CFx product.
    if params.substrate == "carbon":
        recession = substrate_removal / 1.0e29  # a-C atomic density
        si_in = 0.0
        o_in = 0.0
    else:
        recession = (sif4 + bare_removal) / _SIO2_FORMULA_DENSITY_M3  # m/s
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
    layer_f_loss = f_per_unit * f_costed_removal + 2.0 * cof2

    co = (ox_c - 0.5 * ox_f) + (layer_ox - cof2) + (bottom_c - 0.5 * bottom_f)
    cof2_total = cof2 + 0.5 * ox_f + 0.5 * bottom_f

    # Excess layer oxygen beyond a saturated monolayer desorbs recombinatively
    # (O2) — surfaces cannot hold more than saturation coverage. Accounted as
    # outflow so the O ledger still closes exactly.
    n_o_raw = state.n_o + dt * (o_to_layer - layer_side_o)
    o_desorb = np.maximum(n_o_raw - _MONOLAYER_AREAL_M2, 0.0) / max(dt, 1e-30)
    # Layer fluorine saturates at coverage too; the excess recombines/reflects.
    n_f_raw = state.n_f + dt * (mix_f + f_direct + chem_f - layer_f_loss)
    f_desorb = np.maximum(n_f_raw - _MONOLAYER_AREAL_M2, 0.0) / max(dt, 1e-30)
    n_c_raw = state.n_c + dt * (mix_c + chem_c - layer_ox)
    c_desorb = np.maximum(n_c_raw - _MONOLAYER_AREAL_M2, 0.0) / max(dt, 1e-30)

    # Ion-processed-skin crosslinking (zero-knob): every ion converts fresh
    # film atoms at (energy absorbed in the film)/E_displacement per ion;
    # removal channels draw down crosslinked atoms in proportion to their
    # share, so n_xl stays a subset of the film inventory exactly.
    fresh_fraction = np.maximum(
        1.0 - _guarded_ratio(state.n_xl_film, film_total, 1.0), 0.0)
    xl_rate = kernel_xl / params.displacement_energy_eV * fresh_fraction
    # Deposition-driven crosslinking (Krueger Table 6.2 row
    # `P(s) + P(s) -> PC(s) + PC(s)`, module described verbatim in thesis
    # sec. 2.2.3: "Crosslinking occurs during the deposition of eligible
    # materials... During deposition bonds to random eligible cell neighbors
    # can be formed"). Creation is a property of the DEPOSITION event, not of
    # ion dose; ions only break crosslinks (`CF(xs)+M -> CF(s)+M`, 0.3 @ 8 eV,
    # already implemented as dexl). The published row converts BOTH partners,
    # so each deposited unit converts itself plus one eligible fresh neighbour:
    # rate = 2 x (deposited atoms) x (fresh fraction) x (film present).
    # No constant is introduced: the 2 is the row's own stoichiometry and the
    # gates are the eligibility conditions the module states.
    # The partner COUNT is published per species -- Krueger et al., JVST A 42,
    # 043008 (2024), sec. III (tmp/pdfs/krueger-2024.txt L388-390): "based on
    # the number of available bonds (three in the example in Fig. 5).  For
    # example, CF2 would have a maximum of two crosslinks and CF3 would have a
    # maximum of a single crosslink."  Those worked examples pin the rule
    # uniquely (available = 4n - m - 2(n-1) for C_nF_m; CF->3, CF2->2,
    # CF3->1), and the adapter passes the deposition-weighted mean here.  When
    # the caller supplies no bond count the published row stoichiometry alone
    # applies (one partner).
    # This channel does not collapse on near-vertical walls (deposition is
    # isotropic) whereas ion-driven creation does by ~200x through the double
    # cosine -- which is why the lip film crosslinks in his model and did not
    # in ours (x_xl = 0.163 measured against ~0.9 required).
    # See RESULTS_LIP_CROSSLINK_2026-08-04.md.
    if fluxes.deposition_available_bonds is None:
        partners = 1.0
    else:
        partners = np.maximum(
            np.asarray(fluxes.deposition_available_bonds, dtype=float), 0.0)
    xl_rate = xl_rate + ((1.0 + partners)
                         * (dep_c + dep_f) * fresh_fraction * theta_film)
    # Ion de-crosslinking (Appendix S3d: 0.3 @ 8 eV, q=0.5, e0=500): converts
    # crosslinked film back to fresh — mass-neutral within the film.
    dexl_rate = kernel_dexl * _guarded_ratio(state.n_xl_film, film_total, 1.0)
    fresh_available = np.maximum(
        np.asarray(film_total, dtype=float) - state.n_xl_film, 0.0)
    xl_rate = xl_rate * _overdraw_scale(fresh_available, xl_rate, dt)
    film_loss_total = sput_c + sput_f + ox_c + ox_f + mix_c + mix_f + bottom_c + bottom_f
    xl_share = _guarded_ratio(state.n_xl_film, film_total, 1.0)
    n_xl_new = np.maximum(
        state.n_xl_film + dt * (xl_rate - dexl_rate
                                - film_loss_total * xl_share), 0.0)

    if params.substrate == "carbon":
        n_act_new = np.zeros_like(np.asarray(state.n_act, dtype=float))
    else:
        act_formation = kernel_act * open_area * (1.0 - theta_act)
        # Exact site balance: activated sites are consumed by chemisorption
        # THROUGH THE ACTIVATED CHANNEL (theta_act x activated-probability
        # flux), not by the blended total, plus removal of activated surface.
        if fluxes.chemisorption_activated_carbon_flux is not None:
            activated_chem_c = (theta_act * np.asarray(
                fluxes.chemisorption_activated_carbon_flux, dtype=float)
                * site_open)
        else:
            activated_chem_c = theta_act * chem_c
        act_consumption = activated_chem_c + substrate_removal * theta_act
        n_act_new = np.clip(
            state.n_act + dt * (act_formation - act_consumption),
            0.0, _SITE_DENSITY_M2)
    new_state = MixedLayerState(
        n_c_film=state.n_c_film + dt * (dep_c - sput_c - ox_c - mix_c - bottom_c),
        n_f_film=state.n_f_film + dt * (dep_f + absorb_f - sput_f - ox_f - mix_f
                                        - bottom_f),
        n_xl_film=np.minimum(
            n_xl_new,
            np.maximum(
                state.n_c_film + dt * (dep_c - sput_c - ox_c - mix_c - bottom_c)
                + state.n_f_film + dt * (dep_f + absorb_f - sput_f - ox_f - mix_f
                                         - bottom_f), 0.0)),
        n_si=state.n_si + dt * (si_in - sif4),
        n_o=n_o_raw - dt * o_desorb,
        n_c=n_c_raw - dt * c_desorb,
        n_f=n_f_raw - dt * f_desorb,
        n_act=n_act_new,
    )

    residuals = _ledger_residuals(state, new_state, dt,
                                  dep_c + chem_c, dep_f + chem_f, absorb_f,
                                  f_direct, si_in, o_in, sput_c, sput_f,
                                  ox_c + bottom_c, ox_f + bottom_f, mix_c, mix_f,
                                  sif4, layer_ox + c_desorb,
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
        film_deposition_rate=dep_c + dep_f,
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
        for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f",
                     "n_xl_film", "n_act"):
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
    for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f",
                 "n_xl_film", "n_act"):
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
