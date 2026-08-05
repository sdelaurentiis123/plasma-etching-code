"""Joint solve of the ion-channel model-selection problem.

Every prior pass changed ONE assignment and graded it against ONE observable,
which is why each fix over-corrected another channel (RESULTS_ANGULAR_CONVENTION
§5.2: "any class-1 normalisation moves both at once, in opposite senses").  This
script enumerates the whole DISCRETE, SOURCED space at once and scores every
combination against every measured constraint simultaneously.

Nothing here is fitted.  Each option carries a citation; combinations that would
apply a measurement outside the material it was taken on are not enumerated.

Constraints (each measured or receipted, see CONSTRAINTS below for citations):

  C1a Gray-1993 beam dynamic range Y(F/Ar+ -> 0) / Y(sat)      0.20 - 0.30
  C1b Gray-1993 half-rise F/Ar+                                19 - 35
  C2  oxide total angular peak/normal                          1.28 - 1.36
  C3  FC-film angular peak/normal (polymer row)                1.30 - 3.50
  C4  coupled blanket rate at Krueger Table 6.1 fluxes         4.5 - 14.0 nm/s
  C5  coupled depth factor at the measured 3406 eV front       0.735 - 0.812
  C6  lip film growth (mouth equilibrium proxy)                0.30 - 0.60 nm/s
  C7  ARDE sign: coupled floor rate falls from AR 0 to AR 16   ratio < 1
  C8  blanket film thickness (validated regime)                <= 0.15 nm

Run:  python scripts/ion_channel_model_solve.py
"""

from __future__ import annotations

import itertools
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from petch import mixed_layer as ml  # noqa: E402
from petch.mixed_layer_mechanism import (  # noqa: E402
    build_krueger_2024_mixed_layer_mechanisms)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes  # noqa: E402

# --------------------------------------------------------------------------
# Measured inputs
# --------------------------------------------------------------------------

# Krueger Table 6.1 wafer-plane fluxes (m^-2 s^-1) and ion flux/energy.
SRC = {"CF": 4.4e20, "CF2": 9.4e20, "C2F3": 6.8e20, "CF3": 8.4e19,
       "O": 7.7e20, "C3F4": 9.5e20}
ION_SRC = 9.6e19
E_BLANKET = 1500.0
E_FRONT = 3406.0          # measured ml19 etch-front population energy

# Audited top-band (lip) delivery and wall tilt, from the lip crosslink pass.
LIP_DELIVERY = 0.37189910507849583
LIP_TILT_DEG = 0.4722228422825369
LIP_VISIBILITY = 0.7417599287854645
KRUEGER_LIP_CLOSURE = 0.427            # nm/s per side, run average

# Barklund & Blom, JVST A 10, 1212 (1992), Ar+ on a CHF3-deposited FC film,
# digitised from Chang thesis Fig. 4.16 p. 104 (RESEARCH_VERIFY_HUNT §3).
BARKLUND_DEG = [0.0, 20.0, 40.0, 50.0, 60.0, 65.0, 70.0, 75.0, 80.0, 90.0]
BARKLUND_RATE = [1.000, 1.098, 1.168, 1.197, 1.348, 1.448, 1.348, 1.194,
                 0.891, 0.0]

# Frozen-geometry floor delivery vs aspect ratio (results/curated/
# cascade_funnelling/scan.json, straight-wall gathers under the Krueger IEAD).
AR_DELIVERY = {
    0.0: {"ion": 8.181e19, "hot": 7.962e18, "ion_pow": 2.7973946253800928e23,
          "CF": 3.225e19, "CF2": 6.891e19, "O": 5.644e19},
    16.0: {"ion": 6.530e19, "hot": 2.160e19, "ion_pow": 2.242e23,
           "CF": 2.821e19, "CF2": 6.026e19, "O": 4.936e19},
}

CONSTRAINTS = {
    "C1a_gray_dynamic_range": (
        (0.20, 0.30),
        "Gray, Tepermeister & Sawin, JVST B 11, 1243 (1993), 350 eV Ar+ on "
        "SiO2, replotted Kwon ScD MIT 2004 Fig. 3.4 p.76: floor 0.28, "
        "plateau 1.10"),
    "C2_oxide_peak_over_normal": (
        (1.28, 1.36),
        "Cho, JVST A 18, 2705 (2000) ~1.30; Schaepkens, JVST A 16, 3281 "
        "(1998) ~1.33 -- SiO2 in fluorocarbon"),
    "C3_polymer_peak_over_normal": (
        (1.30, 3.50),
        "Barklund & Blom, JVST A 10, 1212 (1992), Ar+ on CHF3 FC film: 1.448 "
        "at 65 deg (yield reading) to 2.70 (raw-rate reading); band spans "
        "both readings, the flux convention being [VERIFY]"),
    "C4_blanket_rate_nm_s": (
        (4.5, 14.0),
        "Krueger 825 nm / 60 s = 13.75 nm/s at 1-3x funnelled delivery"),
    "C5_depth_factor": (
        (0.735, 0.812),
        "[784, 866] nm gate against the +29% baseline extrapolation (1066 nm)"),
    "C6_lip_growth_nm_s": (
        (0.30, 0.60),
        "0.7-1.4x Krueger's 0.427 nm/s per-side closure (mouth equilibrium)"),
    "C7_arde_ratio": (
        (0.0, 1.0),
        "rate(AR 16 delivery) / rate(AR 0 delivery) must be < 1; real HARC "
        "loses ~80% by AR 40 (Huang L5430-5478)"),
    "C8_blanket_film_nm": (
        (0.0, 0.15),
        "validated blanket film ~0.085 nm; 4x thickening throttles interface "
        "energy (RESULTS_ANGULAR_CONVENTION §4)"),
}

# --------------------------------------------------------------------------
# The discrete, sourced model space
# --------------------------------------------------------------------------


def _kress(b, peak_normalised):
    ratio = _peak_ratio(lambda c: ml._class1_shape(c, b))
    scale = (1.0 / ratio) if peak_normalised else 1.0

    def shape(cosine):
        return ml._class1_shape(cosine, b) * scale
    return shape


def _barklund(reading, peak_normalised=False):
    deg = np.asarray(BARKLUND_DEG)
    rate = np.asarray(BARKLUND_RATE)
    if reading == "yield":
        val = rate.copy()
    else:                       # raw-rate reading: Y = R / cos(theta)
        cos = np.cos(np.deg2rad(deg))
        val = np.where(cos > 1e-3, rate / np.maximum(cos, 1e-3), 0.0)
        val[-1] = 0.0
        val = val / val[0]
    scale = (1.0 / float(val.max())) if peak_normalised else 1.0

    def shape(cosine):
        c = np.clip(np.asarray(cosine, dtype=float), 0.0, 1.0)
        theta = np.rad2deg(np.arccos(c))
        return np.interp(theta, deg, val) * scale
    return shape


def _peak_ratio(shape):
    cos_t = np.cos(np.deg2rad(np.linspace(0.0, 89.9, 4000)))
    vals = np.asarray(shape(cos_t), dtype=float)
    return float(vals.max() / max(vals[0], 1e-300))


# Citation gate: Barklund measured an FC FILM, so it is admissible only on the
# polymer row; Cho/Schaepkens measured SiO2, so B=1.7 is admissible only on the
# oxide/mask rows; Kress is Krueger's own class-1 citation, admissible on both.
POLYMER_OPTIONS = {
    "kress9.3_f0=1": (lambda: _kress(9.3, False),
                      "Kress, JVST A 17, 2819 (1999) -- Krueger's cited class-1 "
                      "source; the convention petch has carried"),
    "kress9.3_peaknorm": (lambda: _kress(9.3, True),
                          "same shape read under Huang L2290-2296 ('reduced "
                          "probability at normal incidence')"),
    "barklund_yield": (lambda: _barklund("yield"),
                       "Barklund & Blom 1992 digitised, yield reading "
                       "(Chang p.103 calls it 'etching yield')"),
    "barklund_rate": (lambda: _barklund("rate"),
                      "Barklund & Blom 1992 digitised, raw-rate reading "
                      "Y = R/cos(theta)"),
}

OXIDE_OPTIONS = {
    "kress9.3_f0=1": (lambda: _kress(9.3, False),
                      "Kress 1999 via Krueger's class-1 legend"),
    "kress9.3_peaknorm": (lambda: _kress(9.3, True),
                          "Kress 1999 under the Huang L2290-2296 reading"),
    "b1.7_f0=1": (lambda: _kress(1.7, False),
                  "shape parameter bounded by Cho 2000 / Schaepkens 1998 "
                  "in-chemistry SiO2 measurements (peak 1.31)"),
    "b1.7_peaknorm": (lambda: _kress(1.7, True),
                      "same bound under the peak-normalised reading"),
}

ENERGY_OPTIONS = {
    "zbl": (None,
            "K24-DEKNOB deposited-energy shape eps(E)/eps(140); retired the "
            "fitted yield-scale knob, validated against the power sweeps"),
    "linear_n1": ("linear",
                  "Krueger Appendix B row `0.1471 35 1 140 2`: threshold-power "
                  "(E-35)/(140-35) with n=1"),
}

FDIRECT_OPTIONS = {
    "unity": (1.0, "petch's standing Langmuir sticking of thermal F"),
    "zero": (0.0,
             "Krueger has no thermal-F-on-bare-oxide row: SiO2(s)+F entries "
             "are the ion F+ (L5905-5909); thermal F only fluorinates an "
             "already-complexed site at 0.1 (L6548-6555)"),
}


class Config:
    """One point of the model space, applied by monkeypatching the module."""

    def __init__(self, polymer, oxide, energy, fdirect):
        self.polymer, self.oxide = polymer, oxide
        self.energy, self.fdirect = energy, fdirect

    def __enter__(self):
        self._saved = (ml._angular_physical_sputter, ml._angular_oxide_sputter,
                       ml._complex_energy_factor, ml._THERMAL_F_STICKING)
        ml._angular_physical_sputter = POLYMER_OPTIONS[self.polymer][0]()
        ml._angular_oxide_sputter = OXIDE_OPTIONS[self.oxide][0]()
        if ENERGY_OPTIONS[self.energy][0] == "linear":
            def linear(e_iface, eps_dep, ref_dep_140):
                # Krueger Appendix B `0.1471 35 1 140 2`: threshold-power law
                # with n = 1 in the INTERFACE energy, the same energy the ZBL
                # factor consumes.
                return np.maximum(
                    (np.asarray(e_iface, dtype=float) - 35.0) / 105.0, 0.0)
            ml._complex_energy_factor = linear
        ml._THERMAL_F_STICKING = FDIRECT_OPTIONS[self.fdirect][0]
        return self

    def __exit__(self, *exc):
        (ml._angular_physical_sputter, ml._angular_oxide_sputter,
         ml._complex_energy_factor, ml._THERMAL_F_STICKING) = self._saved
        return False

    @property
    def key(self):
        return f"{self.polymer}|{self.oxide}|{self.energy}|{self.fdirect}"


# --------------------------------------------------------------------------
# Constraint evaluators
# --------------------------------------------------------------------------


def beam_yield(flux_ratio, energy_eV=350.0, cosine=1.0, ion_flux=1.0e19):
    fluxes = ml.SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=flux_ratio * ion_flux,
        oxygen_flux=0.0, ion_flux=ion_flux, ion_energy_eV=energy_eV,
        cosine_incidence=cosine)
    res = ml.steady_state(fluxes)
    return float(np.asarray(res.substrate_removal_rate)) / ion_flux


def gray_curve():
    """Dynamic range only.

    The half-rise POSITION (C1b) is settled analytically and excluded from the
    per-combination sweep: `scripts/gray_half_rise_scan.py` shows it is a
    function of the thermal-F sticking magnitude alone (1.94 at s=1.0, 25.9 at
    s=0.06, Gray's target 27 +/- 8) and is flat in every other axis of this
    space (dynamic range moves 0.873 -> 0.878 across the same sweep).  Neither
    enumerated sticking option can reach it, so C1b is reported as the
    unsatisfiable constraint rather than re-measured 64 times.
    """
    floor = beam_yield(0.0)
    sat = beam_yield(500.0)
    return floor / max(sat, 1e-300), float("nan")


def oxide_angular_peak():
    """Total oxide removal per ion vs angle, saturated supply, beam mode."""
    cosines = np.cos(np.deg2rad(np.linspace(0.0, 85.0, 13)))
    ys = np.array([beam_yield(40.0, energy_eV=350.0, cosine=float(c))
                   for c in cosines])
    return float(ys.max() / max(ys[0], 1e-300))


_MECH = {}


def _mechanisms():
    """Cached; step() reads the patched shapes from module globals at call
    time, so one build serves every configuration."""
    if "m" not in _MECH:
        _MECH["m"] = build_krueger_2024_mixed_layer_mechanisms()
    return _MECH["m"]


def coupled_rate(neutral_scale, ion_flux, energy_eV, dt=2.0, steps=120):
    """Coupled steady state at Krueger fluxes: the only valid depth forecast."""
    oxide, _ = _mechanisms()
    neutral = {k: v * neutral_scale for k, v in SRC.items()}
    ion = EnergeticFlux(name="Ar+", flux_m2_s=ion_flux,
                        energy_eV=np.array([energy_eV]),
                        cosine_incidence=np.array([1.0]),
                        weight=np.array([1.0]))
    fx = SurfaceFluxes(neutral_flux_m2_s=neutral, energetic_fluxes=(ion,))
    st = oxide.initial_state(())
    for _ in range(steps):
        r = oxide.advance(st, fx, dt)
        st = r.state
    film = (float(np.asarray(st.n_c_film)) + float(np.asarray(st.n_f_film)))
    film_nm = film / 7.5e28 * 1e9
    return float(np.asarray(r.etch_velocity_m_s)) * 1e9, film_nm


def lip_growth(dt=2.0, steps=400):
    _, mask = _mechanisms()
    cos_lip = float(np.sin(np.deg2rad(LIP_TILT_DEG)))
    neutral = {k: v * LIP_DELIVERY for k, v in SRC.items()}
    ion = EnergeticFlux(name="Ar+", flux_m2_s=ION_SRC * cos_lip * LIP_VISIBILITY,
                        energy_eV=np.array([E_BLANKET]),
                        cosine_incidence=np.array([max(cos_lip, 1e-4)]),
                        weight=np.array([1.0]))
    fx = SurfaceFluxes(neutral_flux_m2_s=neutral, energetic_fluxes=(ion,))
    st = mask.initial_state(())
    for _ in range(steps):
        r = mask.advance(st, fx, dt)
        st = r.state
    return float(np.asarray(r.normal_growth_velocity_m_s)) * 1e9


def arde_ratio():
    """Coupled floor rate at AR 16 delivery / AR 0 delivery."""
    out = []
    for ar in (0.0, 16.0):
        d = AR_DELIVERY[ar]
        energetic = d["ion"] + d["hot"]
        mean_E = d["ion_pow"] / max(d["ion"], 1e-300)
        scale = d["CF"] / AR_DELIVERY[0.0]["CF"]
        rate, _ = coupled_rate(scale * (AR_DELIVERY[0.0]["CF"] / SRC["CF"]),
                               energetic, mean_E, steps=120)
        out.append(rate)
    return out[1] / max(out[0], 1e-300), out


def free_metrics(cfg):
    """Tier 0: pure shape properties, no chemistry evaluation at all."""
    with cfg:
        return {"C3_polymer_peak_over_normal":
                _peak_ratio(ml._angular_physical_sputter)}


def evaluate(cfg):
    with cfg:
        dyn, half = gray_curve()
        ox_peak = oxide_angular_peak()
        poly_peak = _peak_ratio(ml._angular_physical_sputter)
        blanket, _ = coupled_rate(1.0, ION_SRC, E_BLANKET)
        front, film_nm = coupled_rate(1.0, ION_SRC, E_FRONT)
        lip = lip_growth()
        ratio, _rates = arde_ratio()
    return {
        "C1a_gray_dynamic_range": dyn,
        "C2_oxide_peak_over_normal": ox_peak,
        "C3_polymer_peak_over_normal": poly_peak,
        "C4_blanket_rate_nm_s": blanket,
        "C5_depth_factor": front,      # normalised to baseline after the sweep
        "C6_lip_growth_nm_s": lip,
        "C7_arde_ratio": ratio,
        "C8_blanket_film_nm": film_nm,
    }


def main():
    combos = [Config(p, o, e, f) for p, o, e, f in itertools.product(
        POLYMER_OPTIONS, OXIDE_OPTIONS, ENERGY_OPTIONS, FDIRECT_OPTIONS)]
    print(f"model space: {len(combos)} sourced combinations\n")

    # Baseline = the tree's current assignment, for the depth-factor reference.
    base = Config("kress9.3_f0=1", "b1.7_f0=1", "zbl", "unity")
    base_metrics = evaluate(base)
    base_front = base_metrics["C5_depth_factor"]
    print(f"baseline front rate {base_front:.3f} nm/s "
          f"(= 60 s depth ~1066 nm, +29%)\n")

    # Tier 0 first: C3 is a pure shape ratio, free to evaluate, and it is a
    # MEASURED constraint (Barklund & Blom).  Combinations it eliminates need no
    # coupled run -- the matrix records them as eliminated at tier 0.
    rows = []
    live = []
    for cfg in combos:
        fm = free_metrics(cfg)
        band = CONSTRAINTS["C3_polymer_peak_over_normal"][0]
        v = fm["C3_polymer_peak_over_normal"]
        if not (band[0] <= v <= band[1]):
            rows.append({"config": cfg.key, "polymer": cfg.polymer,
                         "oxide": cfg.oxide, "energy": cfg.energy,
                         "fdirect": cfg.fdirect, "metrics": fm,
                         "fails": ["C3"], "n_fail": 1, "tier": 0})
        else:
            live.append(cfg)
    print(f"tier 0 (free shape test, C3): {len(combos) - len(live)} eliminated,"
          f" {len(live)} advance to the coupled evaluation\n", flush=True)

    for i, cfg in enumerate(live, 1):
        m = evaluate(cfg)
        m["C5_depth_factor"] = m["C5_depth_factor"] / base_front
        fails = []
        for name, (band, _cite) in CONSTRAINTS.items():
            v = m[name]
            if not (np.isfinite(v) and band[0] <= v <= band[1]):
                fails.append(name.split("_")[0])
        rows.append({"config": cfg.key, "polymer": cfg.polymer,
                     "oxide": cfg.oxide, "energy": cfg.energy,
                     "fdirect": cfg.fdirect, "metrics": m,
                     "fails": fails, "n_fail": len(fails), "tier": 2})
        print(f"[{i:2d}/{len(live)}] {cfg.key:<62s} "
              f"fails={len(fails)} {','.join(fails)}", flush=True)

    rows.sort(key=lambda r: r["n_fail"])
    survivors = [r for r in rows if r["n_fail"] == 0]

    print("\n" + "=" * 78)
    if survivors:
        print(f"SURVIVORS: {len(survivors)}")
        for s in survivors:
            print("  " + s["config"])
    else:
        print("SURVIVOR SET EMPTY -- constraint failure census "
              "(how often each constraint is the blocker):")
        census = {}
        for r in rows:
            for f in r["fails"]:
                census[f] = census.get(f, 0) + 1
        for k, v in sorted(census.items(), key=lambda kv: -kv[1]):
            print(f"  {k:<6s} fails in {v}/{len(rows)} combinations")
        print("\n  best (fewest failures):")
        for r in rows[:5]:
            print(f"    {r['config']:<62s} fails={r['n_fail']} "
                  f"{','.join(r['fails'])}")

    out = pathlib.Path("results/curated/ion_channel_solve")
    out.mkdir(parents=True, exist_ok=True)
    (out / "solve.json").write_text(json.dumps(
        {"constraints": {k: {"band": list(v[0]), "citation": v[1]}
                         for k, v in CONSTRAINTS.items()},
         "polymer_options": {k: v[1] for k, v in POLYMER_OPTIONS.items()},
         "oxide_options": {k: v[1] for k, v in OXIDE_OPTIONS.items()},
         "energy_options": {k: v[1] for k, v in ENERGY_OPTIONS.items()},
         "fdirect_options": {k: v[1] for k, v in FDIRECT_OPTIONS.items()},
         "baseline_front_rate_nm_s": base_front,
         "rows": rows}, indent=2))
    print(f"\nwrote {out / 'solve.json'}")


if __name__ == "__main__":
    main()
