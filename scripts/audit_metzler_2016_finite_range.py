#!/usr/bin/env python3
"""Grade finite-range FC-film transport on Metzler's cyclic SiO2 response.

This is a model-form replacement audit, not a calibration.  It repeats the
frozen Figure 6.9(b) surface-only replay with the ZBL/Lindhard CSDA
finite-range transport added to ``mixed_layer.py``.  All legacy surface
kinetics, initial-film brackets, and deliberately non-identified IEDF
sensitivities are held fixed.

The result asks one narrow question: does removing the exponential
transmission tail repair the measured cyclic response?  A failure is
informative because it localizes the remaining problem to the low-energy
boundary and/or mixed reaction-volume kinetics rather than permitting an
attenuation-length fit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from audit_metzler_2016_surface_closure import (
    FC_RATIOS,
    _arcsine_nodes,
    _replay,
    _source_rows,
)
from petch.mixed_layer import _bare_sputter_yield


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "results" / "curated" / "metzler_2016_surface_closure"
    / "finite_range_audit.json"
)
TRANSPORT = "csda_finite_range"


def _metrics(rows):
    measured = np.asarray(
        [row["measured_SiO2_per_incident_ion"] for row in rows],
        dtype=float,
    )
    predicted = np.asarray(
        [row["predicted_SiO2_per_incident_ion"] for row in rows],
        dtype=float,
    )
    ratio = predicted / measured
    return {
        "points": len(rows),
        "minimum_prediction_over_measurement": float(np.min(ratio)),
        "median_prediction_over_measurement": float(np.median(ratio)),
        "maximum_prediction_over_measurement": float(np.max(ratio)),
        "mean_absolute_percentage_error": float(
            np.mean(np.abs(predicted - measured) / measured)
        ),
        "zero_predictions": int(np.count_nonzero(predicted == 0.0)),
        "root_mean_square_log_error_on_positive_predictions": (
            float(np.sqrt(np.mean(np.log(ratio[ratio > 0.0]) ** 2)))
            if np.any(ratio > 0.0) else None
        ),
    }


def _replays(maximum_step_s=0.01):
    rows = _source_rows()
    output = []
    for film_f_over_c in FC_RATIOS:
        for row in rows:
            result = _replay(
                row,
                film_f_over_c,
                "delta_at_reported_maximum",
                maximum_step_s=maximum_step_s,
                film_energy_transport=TRANSPORT,
            )
            result["film_energy_transport"] = TRANSPORT
            output.append(result)
    for spectrum in (
        "arcsine_0_to_reported_maximum",
        "arcsine_14eV_to_reported_maximum",
    ):
        for row in rows:
            result = _replay(
                row,
                0.4,
                spectrum,
                maximum_step_s=maximum_step_s,
                film_energy_transport=TRANSPORT,
            )
            result["film_energy_transport"] = TRANSPORT
            output.append(result)
    return output


def _group_metrics(rows):
    groups = {}
    for row in rows:
        key = (
            row["spectrum_sensitivity_case"],
            row["assumed_initial_film_F_over_C"],
        )
        groups.setdefault(key, []).append(row)
    return [
        {
            "spectrum_sensitivity_case": key[0],
            "assumed_initial_film_F_over_C": key[1],
            **_metrics(members),
        }
        for key, members in sorted(groups.items())
    ]


def _bare_channel_receipt():
    rows = _source_rows()
    receipt = []
    for energy in (25.0, 30.0):
        measured = [
            float(row["substrate_units_per_incident_ion"])
            for row in rows if float(row["energy_eV"]) == energy
        ]
        cases = {
            "delta_at_reported_maximum": np.asarray([energy]),
            "arcsine_0_to_reported_maximum": _arcsine_nodes(
                0.0, energy, count=4096
            ),
            "arcsine_14eV_to_reported_maximum": _arcsine_nodes(
                14.0, energy, count=4096
            ),
        }
        for name, energies in cases.items():
            yield_value = float(np.mean(_bare_sputter_yield(energies)))
            receipt.append({
                "reported_energy_eV": energy,
                "spectrum_sensitivity_case": name,
                "unshielded_Gray_bare_SiO2_yield_per_ion": yield_value,
                "minimum_measured_total_cyclic_yield_per_ion": min(measured),
                "maximum_measured_total_cyclic_yield_per_ion": max(measured),
                "bare_yield_over_maximum_measured_total":
                    yield_value / max(measured),
            })
    return receipt


def build_report():
    replays = _replays()
    fine = _replays(maximum_step_s=0.005)
    relative_changes = []
    for coarse, refined in zip(replays, fine):
        denominator = max(
            refined["predicted_SiO2_per_incident_ion"], 1e-300
        )
        relative_changes.append(
            abs(
                coarse["predicted_SiO2_per_incident_ion"]
                - refined["predicted_SiO2_per_incident_ion"]
            ) / denominator
        )
    metrics = _group_metrics(replays)
    delta = [
        row for row in metrics
        if row["spectrum_sensitivity_case"]
        == "delta_at_reported_maximum"
    ]
    return {
        "audit_id": "METZLER-2016-FINITE-RANGE-SURFACE-CLOSURE-R1",
        "operation": (
            "one-at-a-time model-form replacement; exponential film-energy "
            "transmission is replaced by ZBL/Lindhard CSDA path inversion, "
            "with no parameter or target value fitted"
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
            "surface_points": len(_source_rows()),
        },
        "transport_model": {
            "name": TRANSPORT,
            "physics": (
                "subtract depth/cos(theta) from the CSDA path-to-rest and "
                "invert the same Bragg-additive ZBL nuclear plus "
                "Lindhard-Scharff electronic stopping integral"
            ),
            "free_attenuation_parameters": [],
            "finite_range": True,
            "declared_omissions": [
                "collision-cascade straggling",
                "recoil transport across the film/substrate interface",
                "dynamic film density/composition in the stopping target",
            ],
        },
        "metrics": metrics,
        "replays": replays,
        "bare_channel_spectrum_receipt": _bare_channel_receipt(),
        "integration_convergence": {
            "coarse_maximum_step_s": 0.01,
            "fine_maximum_step_s": 0.005,
            "maximum_relative_prediction_change": float(
                max(relative_changes)
            ),
        },
        "verdict": {
            "finite_range_transport_repairs_surface_response": False,
            "all_maximum_energy_replays_still_overpredict": (
                min(
                    row["minimum_prediction_over_measurement"]
                    for row in delta
                ) > 1.0
            ),
            "smallest_maximum_energy_overprediction_factor": min(
                row["minimum_prediction_over_measurement"] for row in delta
            ),
            "reason": (
                "finite primary-ion range removes the unphysical "
                "transmission tail, but once the modeled film is cleared the "
                "legacy bare/complex kinetics still remove too much oxide; "
                "the direct maximum-energy assumption misses every point"
            ),
            "iedf_sensitivity_identifies_one_boundary": False,
            "reason_iedf": (
                "two unfitted arcsine supports move individual conditions "
                "from zero response to multi-fold overprediction; neither is "
                "a measured IEDF and neither closes all six conditions"
            ),
            "next_required_physics": [
                (
                    "measured or independently derived low-energy IEDF, "
                    "because the nominal 25/30 eV values are maxima in the "
                    "chapter-6 methods while threshold averaging is dominant"
                ),
                (
                    "a stratified mixed reaction volume carrying recoil "
                    "transfer and finite-residence Si/O/C/F inventories, "
                    "rather than lateral bare area times a bulk yield"
                ),
                (
                    "a complete incoming fluorine inventory including "
                    "chamber-wall/redeposition supply, which the thesis says "
                    "can persist after the precursor pulse"
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
            raise SystemExit(f"stale finite-range audit: {OUTPUT}")
        print(f"verified {OUTPUT.relative_to(ROOT)}")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
