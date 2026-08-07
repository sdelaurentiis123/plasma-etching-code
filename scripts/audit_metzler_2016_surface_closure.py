#!/usr/bin/env python3
"""Audit the mixed-layer closure against Metzler's ion-normalized cyclic data.

The useful feature of Metzler Figure 6.9 is that its two axes have the same
incident-Ar-ion denominator.  For a declared initial film inventory,

    ion fluence = initial film F inventory / plotted (film F / incident ion).

That identity permits a surface-only replay without manufacturing a reactor
flux.  It does *not* reconstruct Metzler's wafer flux: the initial film atom
density/composition and the ion-energy distribution were not published.
Those missing quantities are therefore bracketed, never fitted.

The audit intentionally records the pre-correction mixed-layer model's
response before its film transport and surface/subsurface reservoirs are
changed.  It is a falsification receipt, not a calibration routine.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from petch.mixed_layer import (
    MixedLayerParams,
    MixedLayerState,
    SurfaceFluxes,
    _FC_FILM_ATOM_DENSITY_M3,
    _MONOLAYER_AREAL_M2,
    step,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "data" / "experimental" / "metzler_2016"
    / "figure6_9_cycle_averaged_yield.csv")
OUTPUT = (
    ROOT / "results" / "curated" / "metzler_2016_surface_closure"
    / "legacy_mixed_layer_audit.json")
FILM_THICKNESS_NM = 0.5
FC_RATIOS = (0.4, 1.0, 1.5, 2.0)
ARCSINE_NODES = 32


def _source_rows():
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return [
        row for row in rows
        if row["panel"] == "6.9b" and row["substrate"] == "SiO2"
    ]


def _arcsine_nodes(lower_eV, upper_eV, count=ARCSINE_NODES):
    """Equal-probability midpoint quadrature for an arcsine distribution."""
    phase = (np.arange(count, dtype=float) + 0.5) * np.pi / count
    middle = 0.5 * (lower_eV + upper_eV)
    half_width = 0.5 * (upper_eV - lower_eV)
    return middle + half_width * np.cos(phase)


def _surface_flux(ion_flux, maximum_energy_eV, spectrum):
    if spectrum == "delta_at_reported_maximum":
        return SurfaceFluxes(
            precursor_flux=0.0,
            fluorine_flux=0.0,
            oxygen_flux=0.0,
            ion_flux=ion_flux,
            ion_energy_eV=maximum_energy_eV,
        )
    lower = {
        "arcsine_0_to_reported_maximum": 0.0,
        "arcsine_14eV_to_reported_maximum": 14.0,
    }[spectrum]
    energies = _arcsine_nodes(lower, maximum_energy_eV)
    event_flux = np.full(energies.shape, ion_flux / len(energies))
    return SurfaceFluxes(
        precursor_flux=0.0,
        fluorine_flux=0.0,
        oxygen_flux=0.0,
        ion_flux=ion_flux,
        ion_energy_eV=float(np.mean(energies)),
        ion_atom_face=np.zeros(energies.shape, dtype=int),
        ion_atom_flux=event_flux,
        ion_atom_energy_eV=energies,
        ion_atom_cosine=np.ones(energies.shape),
    )


def _replay(
        row, film_f_over_c, spectrum, maximum_step_s=0.01,
        film_energy_transport="legacy_exponential"):
    energy_eV = float(row["energy_eV"])
    duration_s = float(row["etch_step_s"])
    source_f_per_ion = float(row["fluorine_per_incident_ion"])
    measured_yield = float(row["substrate_units_per_incident_ion"])

    film_atoms_m2 = FILM_THICKNESS_NM * 1e-9 * _FC_FILM_ATOM_DENSITY_M3
    film_c_m2 = film_atoms_m2 / (1.0 + film_f_over_c)
    film_f_m2 = film_atoms_m2 - film_c_m2
    ion_fluence_m2 = film_f_m2 / source_f_per_ion
    ion_flux_m2_s = ion_fluence_m2 / duration_s

    state = MixedLayerState(
        n_c_film=film_c_m2,
        n_f_film=film_f_m2,
    )
    flux = _surface_flux(ion_flux_m2_s, energy_eV, spectrum)
    # Bound each update to at most one tenth of a monolayer per incident-ion
    # fluence as well as the declared wall-clock ceiling.
    dt = min(maximum_step_s, 0.1 * _MONOLAYER_AREAL_M2 / ion_flux_m2_s)
    elapsed_s = 0.0
    removed_formula_m2 = 0.0
    largest_ledger_residual = 0.0
    params = MixedLayerParams(
        substrate="sio2",
        film_energy_transport=film_energy_transport,
    )
    while elapsed_s < duration_s:
        step_s = min(dt, duration_s - elapsed_s)
        result = step(state, flux, step_s, params)
        removed_formula_m2 += (
            float(np.asarray(result.substrate_removal_rate)) * step_s)
        largest_ledger_residual = max(
            largest_ledger_residual,
            max(abs(float(np.asarray(value)))
                for value in result.ledger_residuals.values()),
        )
        state = result.state
        elapsed_s += step_s

    predicted_yield = removed_formula_m2 / ion_fluence_m2
    return {
        "observation_id": row["observation_id"],
        "reported_maximum_energy_eV": energy_eV,
        "etch_step_s": duration_s,
        "film_F_per_incident_ion": source_f_per_ion,
        "measured_SiO2_per_incident_ion": measured_yield,
        "assumed_initial_film_F_over_C": film_f_over_c,
        "assumed_initial_film_atoms_m2": film_atoms_m2,
        "algebraic_ion_fluence_m2": ion_fluence_m2,
        "algebraic_average_ion_flux_m2_s": ion_flux_m2_s,
        "spectrum_sensitivity_case": spectrum,
        "predicted_SiO2_per_incident_ion": predicted_yield,
        "prediction_over_measurement": predicted_yield / measured_yield,
        "final_film_thickness_nm": float(state.film_thickness_nm()),
        "maximum_absolute_element_ledger_residual_atoms_m2":
            largest_ledger_residual,
        "integration_step_s": dt,
    }


def _metrics(rows):
    ratios = np.asarray(
        [row["prediction_over_measurement"] for row in rows], dtype=float)
    predicted = np.asarray(
        [row["predicted_SiO2_per_incident_ion"] for row in rows], dtype=float)
    measured = np.asarray(
        [row["measured_SiO2_per_incident_ion"] for row in rows], dtype=float)
    return {
        "points": len(rows),
        "minimum_prediction_over_measurement": float(np.min(ratios)),
        "median_prediction_over_measurement": float(np.median(ratios)),
        "maximum_prediction_over_measurement": float(np.max(ratios)),
        "mean_absolute_percentage_error": float(
            np.mean(np.abs(predicted - measured) / measured)),
        "root_mean_square_log_error": float(
            np.sqrt(np.mean(np.log(predicted / measured) ** 2))),
    }


def build_report():
    source_rows = _source_rows()
    replays = []
    for film_f_over_c in FC_RATIOS:
        for row in source_rows:
            replays.append(_replay(
                row, film_f_over_c, "delta_at_reported_maximum"))
    # IEDF support was not published.  These two spectra are deliberately
    # sensitivity extremes at one composition, not candidate reconstructions.
    for spectrum in (
            "arcsine_0_to_reported_maximum",
            "arcsine_14eV_to_reported_maximum"):
        for row in source_rows:
            replays.append(_replay(row, 0.4, spectrum))

    grouped = {}
    for row in replays:
        key = (
            row["spectrum_sensitivity_case"],
            row["assumed_initial_film_F_over_C"],
        )
        grouped.setdefault(key, []).append(row)
    metrics = [
        {
            "spectrum_sensitivity_case": key[0],
            "assumed_initial_film_F_over_C": key[1],
            **_metrics(rows),
        }
        for key, rows in sorted(grouped.items())
    ]

    # Verify that the stored integration ceiling is converged rather than
    # silently setting the result.
    coarse = [
        _replay(row, 0.4, "delta_at_reported_maximum", maximum_step_s=0.01)
        for row in source_rows
    ]
    fine = [
        _replay(row, 0.4, "delta_at_reported_maximum", maximum_step_s=0.005)
        for row in source_rows
    ]
    integration_change = [
        abs(a["predicted_SiO2_per_incident_ion"]
            - b["predicted_SiO2_per_incident_ion"])
        / b["predicted_SiO2_per_incident_ion"]
        for a, b in zip(coarse, fine)
    ]

    return {
        "audit_id": "METZLER-2016-SURFACE-CLOSURE-LEGACY-R1",
        "operation": (
            "surface-only falsification replay; no model parameter, film "
            "composition, ion spectrum, or measured point was fitted"
        ),
        "source": {
            "citation": (
                "D. Metzler, High Precision Plasma Etch for Pattern "
                "Transfer, PhD thesis, University of Maryland (2016), "
                "Figure 6.9(b)"
            ),
            "source_pdf_sha256": (
                "ea5701d0bcf67b56403625253f7bb619c0e3e3a5e0a9cfd2ff6f1e435fa90f62"
            ),
            "data_path": str(SOURCE.relative_to(ROOT)),
            "points_used": len(source_rows),
            "observable": (
                "cycle-averaged SiO2 formula units removed per incident Ar "
                "ion versus deposited-film F atoms per the same incident ion"
            ),
        },
        "inference_chain": [
            (
                "initialize the reported 5 A optical film using the declared "
                "model film atom density and one bracketed initial F/C ratio"
            ),
            (
                "divide that assumed initial F inventory by the author's "
                "plotted F-per-ion ratio to obtain an algebraic fluence"
            ),
            (
                "divide fluence by the reported 20, 40, or 60 s ion-step "
                "duration only to integrate the transient; renormalize "
                "removed SiO2 by the same fluence"
            ),
            (
                "set precursor, free-F, and O influxes to zero during the "
                "Ar ion step, as required by the separated cyclic sequence"
            ),
        ],
        "non_inferences": [
            (
                "the algebraic average ion flux is not an independently "
                "measured or reconstructed reactor boundary"
            ),
            (
                "ellipsometric thickness is not promoted to an experimental "
                "absolute C/F inventory"
            ),
            (
                "neither arcsine sensitivity spectrum is asserted to be "
                "Metzler's unreported IEDF"
            ),
        ],
        "declared_assumptions": {
            "film_thickness_nm": FILM_THICKNESS_NM,
            "film_atom_density_m3": _FC_FILM_ATOM_DENSITY_M3,
            "film_F_over_C_sensitivity": list(FC_RATIOS),
            "ion_incidence_cosine": 1.0,
            "arcsine_quadrature_nodes": ARCSINE_NODES,
            "unreported_IEDF_sensitivity_support_eV": [
                "[0, reported maximum]",
                "[14, reported maximum]",
            ],
        },
        "integration_convergence": {
            "coarse_maximum_step_s": 0.01,
            "fine_maximum_step_s": 0.005,
            "maximum_relative_prediction_change": float(
                max(integration_change)),
        },
        "metrics": metrics,
        "replays": replays,
        "verdict": {
            "legacy_delta_energy_closure_passes": False,
            "reason": (
                "for every bracketed initial film F/C ratio, every "
                "monoenergetic replay overpredicts the direct SiO2 response "
                "by multiple factors; the failure is far larger than "
                "digitization placement and integration error"
            ),
            "unreported_IEDF_closes_model": False,
            "reason_IEDF": (
                "unfitted broad-spectrum sensitivities span portions of the "
                "magnitude but do not identify the source IEDF or repair the "
                "duration/energy response as one declared boundary"
            ),
            "model_form_defects_to_test_next": [
                (
                    "exponential interface-energy attenuation is not a "
                    "published Standaert law and violates finite ion range"
                ),
                (
                    "chemical substrate removal is multiplied by lateral "
                    "bare area even though Metzler directly observes etching "
                    "while an ultrathin fluorocarbon film remains"
                ),
                (
                    "one-monolayer reactive coverage is conflated with the "
                    "measured 1.5-3 nm ion-mixed subsurface reservoir"
                ),
            ],
        },
    }


def canonical_payload(report):
    return json.dumps(report, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = canonical_payload(build_report())
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(payload, encoding="utf-8")
        print(f"wrote {OUTPUT.relative_to(ROOT)}")
    elif args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"stale Metzler surface audit: {OUTPUT}")
        print(f"verified {OUTPUT.relative_to(ROOT)}")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
