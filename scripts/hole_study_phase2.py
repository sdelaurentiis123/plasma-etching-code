"""Phase-2 of the HAR hole study: coupled profile evolution.

Phase 1 (``scripts/hole_study_phase1.py``) characterised transport through a
frozen hole and had to stop there: converting delivered flux into an etch rate
needs the surface model, and the axisymmetric path had no evolution driver.
:mod:`petch.axisymmetric_evolution` supplies one, so this script reports what
phase 1 could not.

Series A  coupled rate-ARDE: floor etch rate vs aspect ratio, swept over the
          declared tail-fraction band.  This is the quantity phase 1's honesty
          appendix marked blocked -- its section 5 reports *energetic delivery*
          vs AR and explicitly forbids reading it as an etch-rate prediction.
Series B  thermalised-return (E8) sweep at high aspect ratio: does the reborn
          radical channel relieve the neutral starvation that sets the rate?
          The fluorocarbon share is unpublished for every reactor in the
          corpus, so it is swept over its physical band, never fitted.
Series C  evolution to the declared straight-wall envelope: how far a hole
          etches before wall passivation makes the profile non-straight, which
          is where the exact cylinder operators stop being the right tool.

Everything here is deterministic and runs on a laptop in minutes -- no Monte
Carlo, no box, no fitted constant.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from petch.axisymmetric_evolution import (  # noqa: E402
    HoleEvolutionState,
    HoleGeometry,
    advance_hole_step,
    evolve_hole,
)
from petch.iadf_two_component import kim_2025_reference_iadf  # noqa: E402
from petch.mixed_layer_mechanism import (  # noqa: E402
    MixedLayerSurfaceState,
    build_krueger_2024_mixed_layer_mechanisms,
)
from petch.reactor_boundary import (  # noqa: E402
    build_krueger_2024_development_boundary,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "krueger_2024"

#: Hole radius: the plan's 90 nm opening, so AR = depth / 90 nm.
RADIUS_M = 45e-9
#: Reference ion energy, phase-1 convention (inside the Kim 2025 1.4-2.0 keV band).
REFERENCE_ENERGY_EV = 1500.0
TAIL_FRACTIONS = (0.0, 0.35, 0.50, 0.65)
ASPECT_RATIOS = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 50.0, 100.0, 150.0, 200.0)
#: Cascade quadrature.  Phase-1 production orders; the regression gate in
#: tests/test_axisymmetric_evolution.py pins this path to phase 1 bitwise.
CASCADE = dict(n_polar=192, n_azimuth=64, n_radial=24)
#: Relaxation time before reading a rate: the mixed layer must reach its own
#: steady state at the delivered fluxes, otherwise the reported rate is the
#: bare-surface transient.
RELAX_S = 0.2
#: Stoichiometry of the thermalised return, by carbon/fluorine content of the
#: reborn species.  Huang sec. 6.4.3: thermalised CFx+/CxFy+ become CFx/CxFy
#: radicals; the partners of other ions are non-reactive and diffuse out.
THERMALIZED_STOICHIOMETRY = {"CF2": 1.0}


def boundary_fluxes():
    boundary = build_krueger_2024_development_boundary(
        DATA, reference_plane_m=2.0e-6)
    mouth = {item.name: float(item.flux_m2_s)
             for item in boundary.species if item.charge_number == 0}
    return mouth, float(boundary.get("ions").flux_m2_s)


def _fresh(aspect, bands=40):
    depth = float(aspect) * 2.0 * RADIUS_M
    height = max(depth / float(bands), 2.0 * RADIUS_M / 4.0)
    geometry = HoleGeometry.straight(RADIUS_M, depth, height)
    return HoleEvolutionState(
        geometry, MixedLayerSurfaceState.bare((geometry.band_count + 1,)))


def frozen_rate(mechanism, mouth, ion, iadf, aspect, **kwargs):
    """Relaxed floor etch rate at a frozen hole of the given aspect ratio."""
    state = _fresh(aspect)
    _next, record = advance_hole_step(
        state, mechanism, mouth_flux_m2_s=mouth, ion_flux_m2_s=ion, iadf=iadf,
        energy_eV=REFERENCE_ENERGY_EV, dt_s=RELAX_S, cascade_kwargs=CASCADE,
        **kwargs)
    return record


def series_a_rate_arde(mechanism, mouth, ion):
    rows = []
    for tail in TAIL_FRACTIONS:
        iadf = kim_2025_reference_iadf(tail_fraction=tail)
        for aspect in ASPECT_RATIOS:
            record = frozen_rate(mechanism, mouth, ion, iadf, aspect)
            rows.append({
                "tail_fraction": float(tail),
                "aspect_ratio": float(aspect),
                "floor_etch_nm_s": float(record["floor_etch_velocity_m_s"] * 1e9),
                "energetic_delivery": float(record["cascade_bottom_delivery"]),
                "max_wall_growth_nm_s": float(
                    record["max_wall_growth_velocity_m_s"] * 1e9),
            })
    return rows


def series_b_thermalized_sweep(mechanism, mouth, ion, fractions, aspects):
    iadf = kim_2025_reference_iadf(tail_fraction=0.65)
    rows = []
    for aspect in aspects:
        for fraction in fractions:
            record = frozen_rate(
                mechanism, mouth, ion, iadf, aspect,
                thermalized_return_fraction=float(fraction),
                thermalized_stoichiometry=THERMALIZED_STOICHIOMETRY)
            rows.append({
                "aspect_ratio": float(aspect),
                "thermalized_return_fraction": float(fraction),
                "floor_etch_nm_s": float(record["floor_etch_velocity_m_s"] * 1e9),
                "thermalised_born_rate_s": float(record["thermalised_born_rate_s"]),
            })
    return rows


def series_c_envelope(mechanism, mouth, ion, aspect, *, duration_s, dt_s,
                      tolerance):
    iadf = kim_2025_reference_iadf(tail_fraction=0.65)
    state = _fresh(aspect)
    initial = state.geometry
    final, records, reason = evolve_hole(
        state, mechanism, mouth_flux_m2_s=mouth, ion_flux_m2_s=ion, iadf=iadf,
        energy_eV=REFERENCE_ENERGY_EV, duration_s=duration_s, dt_s=dt_s,
        straightness_tolerance=tolerance, cascade_kwargs=CASCADE)
    return {
        "aspect_ratio_initial": float(initial.aspect_ratio),
        "straightness_tolerance": float(tolerance),
        "stop_reason": reason,
        "steps": len(records),
        "process_time_s": float(records[-1]["time_s"]) if records else 0.0,
        "depth_advance_nm": float(
            (final.geometry.floor_depth - initial.floor_depth) * 1e9),
        "mouth_radius_loss_nm": float(
            (initial.radius[0] - final.geometry.radius[0]) * 1e9),
        "final_straightness_deviation": float(
            records[-1]["straightness_deviation"]) if records else 0.0,
        "trace": [
            {"time_s": row["time_s"],
             "floor_depth_nm": row["floor_depth_m"] * 1e9,
             "mouth_radius_nm": row["mouth_radius_m"] * 1e9,
             "straightness_deviation": row["straightness_deviation"],
             "floor_etch_nm_s": row["floor_etch_velocity_m_s"] * 1e9}
            for row in records],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(
        ROOT / "results" / "curated" / "hole_study" / "phase2.json"))
    parser.add_argument("--envelope-duration-s", type=float, default=20.0)
    parser.add_argument("--envelope-dt-s", type=float, default=0.25)
    parser.add_argument("--straightness-tolerance", type=float, default=0.02)
    args = parser.parse_args()

    mouth, ion = boundary_fluxes()
    oxide, _mask = build_krueger_2024_mixed_layer_mechanisms()

    payload = {
        "provenance": {
            "radius_m": RADIUS_M,
            "reference_energy_eV": REFERENCE_ENERGY_EV,
            "relaxation_s": RELAX_S,
            "cascade_quadrature": CASCADE,
            "mouth_flux_m2_s": mouth,
            "ion_flux_m2_s": ion,
            "beam": "kim_2025_reference_iadf",
            "chemistry": "krueger_2024 deck 1 (oxide arm), Gray-anchored ion laws",
        },
        "series_a_rate_arde": series_a_rate_arde(oxide, mouth, ion),
        "series_b_thermalized": series_b_thermalized_sweep(
            oxide, mouth, ion, (0.0, 0.3, 0.65, 1.0), (50.0, 100.0, 200.0)),
        "series_c_envelope": [
            series_c_envelope(oxide, mouth, ion, aspect,
                              duration_s=args.envelope_duration_s,
                              dt_s=args.envelope_dt_s,
                              tolerance=args.straightness_tolerance)
            for aspect in (10.0, 50.0)],
    }
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))

    rows = payload["series_a_rate_arde"]
    print("series A -- coupled floor etch rate vs aspect ratio (nm/s)")
    header = "  AR   " + "".join(f"  tail {t:.2f}" for t in TAIL_FRACTIONS)
    print(header)
    for aspect in ASPECT_RATIOS:
        cells = []
        for tail in TAIL_FRACTIONS:
            value = next(r["floor_etch_nm_s"] for r in rows
                         if r["aspect_ratio"] == aspect and r["tail_fraction"] == tail)
            cells.append(f"   {value:7.3f}")
        print(f"  {aspect:5.0f}" + "".join(cells))
    for tail in TAIL_FRACTIONS:
        values = [r["floor_etch_nm_s"] for r in rows if r["tail_fraction"] == tail]
        print(f"  tail {tail:.2f}: rate(AR200)/rate(AR1) = {values[-1] / values[0]:.4f}")
    print("\nseries B -- thermalised return sweep")
    for row in payload["series_b_thermalized"]:
        print(f"  AR {row['aspect_ratio']:5.0f}  fraction {row['thermalized_return_fraction']:.2f}"
              f"  {row['floor_etch_nm_s']:7.4f} nm/s")
    print("\nseries C -- straight-wall envelope")
    for row in payload["series_c_envelope"]:
        print(f"  AR {row['aspect_ratio_initial']:5.1f}: {row['stop_reason']} after "
              f"{row['process_time_s']:.2f} s, depth +{row['depth_advance_nm']:.2f} nm, "
              f"mouth radius -{row['mouth_radius_loss_nm']:.2f} nm")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
