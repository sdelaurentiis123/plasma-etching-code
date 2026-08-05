"""Discriminate the class-1 angular NORMALISATION convention against data.

The question: in Krueger Eq. (2.40) `p(E,t) = p0 * energy(E) * f(t)`, is the
tabulated `p0` the yield at NORMAL incidence (f(0) = 1, petch's current
reading) or the yield at the shape's PEAK (f(peak) = 1, so f(0) < 1)?

No source states the normalisation directly.  The discriminating verbatim is
Huang L2290-2296, which describes the two classes in one breath:

  "For physical sputtering, f(theta) is an empirical function with a maximum
   at 60 deg, REDUCED PROBABILITY AT NORMAL INCIDENCE and zero probability at
   grazing incidence.  For chemically enhanced etching, f(theta) is UNITY FOR
   NORMAL INCIDENCE and angles up to 45 deg, with a monotonic roll-off to
   zero probability at grazing incidence."

The sentence only carries information if the two classes share one scale:
class 2 is stated to be unity at normal, class 1 is stated to be reduced
there.  Under petch's current f(0) = 1 for both, "reduced probability at
normal incidence" would be false for class 1 -- it would equal class 2 there.

This script does not assert that reading; it measures what each convention
predicts against published data, at the HARC floor condition (near-normal
incidence) where the conventions differ most:

  GATE N1 dynamic range  Y(F/Ar+ -> 0) / Y(saturated), target 0.25 +/- 0.05
  (Gray 1993 via Kwon Fig. 3.4, 350 eV);
  floor removal rate at the measured ml19 etch-front energy, whose current
  value is +29% against the 825 nm depth gate.

Run:  python scripts/angular_convention_discriminator.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from petch import mixed_layer as ml  # noqa: E402

GRAY_DYNAMIC_RANGE = 0.28 / 1.10
FRONT_ENERGY_EV = 3406.0            # measured ml19 etch-front population
ION_FLUX = 1.0e19


def class1_peak_ratio(b):
    """max f / f(0) for (1 + B sin^2 t) cos t, and the peak angle."""
    cos_t = np.cos(np.deg2rad(np.linspace(0.0, 89.9, 20000)))
    shape = (1.0 + b * (1.0 - cos_t ** 2)) * cos_t
    idx = int(np.argmax(shape))
    peak_deg = float(np.rad2deg(np.arccos(cos_t[idx])))
    return float(shape[idx] / shape[0]), peak_deg


def with_convention(b_oxide, peak_normalised, fn):
    """Run fn() with the oxide class-1 shape under a given convention."""
    original = ml._angular_oxide_sputter
    ratio, _ = class1_peak_ratio(b_oxide)
    scale = (1.0 / ratio) if peak_normalised else 1.0

    def patched(cosine):
        return ml._class1_shape(cosine, b_oxide) * scale

    ml._angular_oxide_sputter = patched
    try:
        return fn()
    finally:
        ml._angular_oxide_sputter = original


def beam_yield(flux_ratio, energy_eV=350.0, cosine=1.0):
    fluxes = ml.SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=flux_ratio * ION_FLUX,
        oxygen_flux=0.0, ion_flux=ION_FLUX, ion_energy_eV=energy_eV,
        cosine_incidence=cosine)
    res = ml.steady_state(fluxes)
    return float(np.asarray(res.substrate_removal_rate)) / ION_FLUX


def floor_rate(cosine=1.0):
    """Removal per ion at the measured etch-front energy, saturated supply."""
    return beam_yield(40.0, energy_eV=FRONT_ENERGY_EV, cosine=cosine)


def main():
    cases = [("f(0)=1 (current)", 1.7, False),
             ("peak-normalised", 1.7, True),
             ("f(0)=1, Kress B", 9.3, False),
             ("peak-normalised, Kress B", 9.3, True)]

    rows = []
    for label, b, peak_norm in cases:
        ratio, peak_deg = class1_peak_ratio(b)

        def measure():
            floor = beam_yield(0.0)
            sat = beam_yield(500.0)
            return floor, sat, floor_rate(1.0)

        floor, sat, front = with_convention(b, peak_norm, measure)
        rows.append({
            "convention": label, "B": b, "peak_normalised": peak_norm,
            "peak_over_normal": ratio, "peak_deg": peak_deg,
            "f_at_normal": (1.0 / ratio) if peak_norm else 1.0,
            "n1_floor": floor, "n1_saturated": sat,
            "n1_dynamic_range": floor / sat,
            "front_removal_per_ion": front,
        })

    base_front = rows[0]["front_removal_per_ion"]
    print(f"{'convention':>26} {'f(0)':>7} {'peak/n':>7} {'Y0/Ysat':>9}"
          f" {'vs Gray':>9} {'front':>8} {'rel':>7}")
    for r in rows:
        err = r["n1_dynamic_range"] - GRAY_DYNAMIC_RANGE
        r["front_relative"] = r["front_removal_per_ion"] / base_front
        print(f"{r['convention']:>26} {r['f_at_normal']:7.3f}"
              f" {r['peak_over_normal']:7.2f} {r['n1_dynamic_range']:9.3f}"
              f" {err:+9.3f} {r['front_removal_per_ion']:8.3f}"
              f" {r['front_relative']:7.3f}")
    print(f"\n  Gray 1993 target dynamic range = {GRAY_DYNAMIC_RANGE:.3f}"
          f"  (band 0.20-0.30)")
    print("  'front' = SiO2 units removed per ion at the measured 3406 eV")
    print("  etch-front energy, saturated supply; 'rel' is vs the current")
    print("  convention.  Depth is currently +29% against the 825 nm gate,")
    print("  so a convention landing rel ~ 0.78 would close it.")

    out = pathlib.Path("results/curated/neutral_limited_gates")
    out.mkdir(parents=True, exist_ok=True)
    (out / "angular_convention.json").write_text(json.dumps(
        {"gray_dynamic_range": GRAY_DYNAMIC_RANGE,
         "front_energy_eV": FRONT_ENERGY_EV, "cases": rows}, indent=2))
    print(f"\nwrote {out / 'angular_convention.json'}")


if __name__ == "__main__":
    main()
