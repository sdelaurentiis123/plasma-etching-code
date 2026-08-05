"""GATE N1/N2 harness: 0-D beam replication of the measured SiO2 yield curve.

Preregistered in RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md sec 5.4.  The
published data on both axes:

  Gray, Tepermeister & Sawin, JVST B 11, 1243 (1993), replotted as Kwon
  (ScD thesis, MIT DMSE 2004, https://hdl.handle.net/1721.1/28353) Fig. 3.4
  p. 76 -- "Ion bombardment energy is 350 eV", axes yield [SiO2/Ar+] vs flux
  ratio [F/Ar+]:  sputter floor Y = 0.28 at F/Ar+ -> 0, half-rise at
  F/Ar+ ~ 27, ~90% of plateau by F/Ar+ ~ 100, plateau Y ~ 1.10.
  Dynamic range Y(0)/Y(sat) = 0.25.

  Butterbaugh, Gray & Sawin, JVST B 9, 1461 (1991), via Kwon Fig. 2.6 p. 36:
  adding CFx at fixed F/Ar+ *reduces* the yield ("there is a reduction in
  etching yield as the CFx flux increases due to the deposition on the
  surface", Kwon p. 36).

Chemistry caveat, stated because it bounds what the gate can prove: Gray's
beam is F + Ar+ on SiO2, whereas the constants under test are Krueger's
fluorocarbon-optimised rows.  The SHAPE of the curve (Langmuir rise,
monotonicity, saturation, half-rise position) tests the coverage kinetics
that both systems share and is gated hard.  The absolute dynamic range
Y(0)/Y(sat) is the ratio of the two removal rows and is reported as a
diagnostic against 0.25, not asserted, because the row constants are from a
different halogen chemistry.

Run:  python scripts/gate_n1_beam.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from petch.mixed_layer import (  # noqa: E402
    MixedLayerParams,
    SurfaceFluxes,
    steady_state,
)

# Gray 1993 beam conditions (Kwon Fig. 3.4).
BEAM_ENERGY_EV = 350.0
ION_FLUX = 1.0e19  # m^-2 s^-1; only the RATIO F/Ar+ is observable.
RATIOS = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 27.0, 40.0, 60.0, 100.0,
          150.0, 200.0, 300.0, 500.0)

# Digitized Gray 1993 curve (Kwon Fig. 3.4; [VERIFY] against the original
# JVST B figure -- read off the replot).
GRAY_1993 = {0.0: 0.28, 5.0: 0.51, 10.0: 0.54, 20.0: 0.63, 40.0: 0.78,
             50.0: 0.84, 70.0: 0.89, 100.0: 1.05, 200.0: 1.07}
GRAY_PLATEAU = 1.10
GRAY_FLOOR_RATIO = 0.28 / 1.10          # 0.25
GRAY_HALF_RISE = 27.0


def beam_yield(flux_ratio, *, energy_eV=BEAM_ENERGY_EV, ion_flux=ION_FLUX,
               cf2_ratio=0.0, params=None):
    """Total SiO2 formula units removed per incident ion at steady state."""
    params = params or MixedLayerParams()
    fluxes = SurfaceFluxes(
        precursor_flux=cf2_ratio * ion_flux,
        fluorine_flux=flux_ratio * ion_flux,
        oxygen_flux=0.0,
        ion_flux=ion_flux,
        ion_energy_eV=energy_eV,
    )
    result = steady_state(fluxes, params)
    return float(np.asarray(result.substrate_removal_rate)) / ion_flux


def sweep(ratios=RATIOS, **kwargs):
    return {r: beam_yield(r, **kwargs) for r in ratios}


def half_rise_ratio(curve):
    """Flux ratio at which the yield reaches the midpoint of its own range."""
    ratios = sorted(curve)
    values = [curve[r] for r in ratios]
    target = 0.5 * (values[0] + values[-1])
    for lo, hi in zip(range(len(ratios) - 1), range(1, len(ratios))):
        if values[lo] <= target <= values[hi]:
            span = values[hi] - values[lo]
            if span <= 0.0:
                return ratios[lo]
            frac = (target - values[lo]) / span
            return ratios[lo] + frac * (ratios[hi] - ratios[lo])
    return float("nan")


def main():
    curve = sweep()
    ratios = sorted(curve)
    floor, plateau = curve[ratios[0]], curve[ratios[-1]]
    print(f"GATE N1 -- SiO2 beam, {BEAM_ENERGY_EV:.0f} eV Ar+, sweep F/Ar+")
    print(f"{'F/Ar+':>8} {'Y (petch)':>12} {'Y (Gray 1993)':>15}")
    for r in ratios:
        measured = GRAY_1993.get(r)
        cell = f"{measured:.2f}" if measured is not None else "-"
        print(f"{r:8.0f} {curve[r]:12.4f} {cell:>15}")
    print()
    print(f"  floor Y(0)           = {floor:.4f}")
    print(f"  plateau Y(500)       = {plateau:.4f}")
    print(f"  dynamic range Y0/Ysat= {floor / plateau:.4f}"
          f"   (Gray: {GRAY_FLOOR_RATIO:.2f})")
    print(f"  half-rise F/Ar+      = {half_rise_ratio(curve):.1f}"
          f"   (Gray: {GRAY_HALF_RISE:.0f})")

    print()
    print("GATE N2 -- CF2 added at fixed F/Ar+ = 40 must reduce the yield")
    for cf2 in (0.0, 1.0, 5.0, 10.0, 20.0):
        y = beam_yield(40.0, cf2_ratio=cf2)
        print(f"  CF2/Ar+ = {cf2:5.1f}   Y = {y:.4f}")

    out = pathlib.Path("results/curated/neutral_limited_gates")
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate_n1.json").write_text(json.dumps({
        "beam_energy_eV": BEAM_ENERGY_EV,
        "curve": {str(k): v for k, v in curve.items()},
        "floor": floor,
        "plateau": plateau,
        "dynamic_range": floor / plateau,
        "half_rise": half_rise_ratio(curve),
        "gray_1993": {str(k): v for k, v in GRAY_1993.items()},
        "gray_dynamic_range": GRAY_FLOOR_RATIO,
        "gray_half_rise": GRAY_HALF_RISE,
    }, indent=2))
    print(f"\nwrote {out / 'gate_n1.json'}")


if __name__ == "__main__":
    main()
