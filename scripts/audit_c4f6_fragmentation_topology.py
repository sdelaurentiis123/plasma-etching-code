#!/usr/bin/env python3
"""Test whether direct C4F6 EI fragmentation can explain reactor ions.

This is a topology audit, not a plasma fit.  The NIST 70 eV spectrum fixes a
single-collision reference ratio, while Benck et al. independently measured
mass-resolved positive-ion currents in a C4F6/Ar ICP.  Their absolute scales
are deliberately not equated.  Only within-board ratios and trends are used
to decide whether a direct-fragmentation-only chemistry is structurally
adequate.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NIST_DIRECTORY = ROOT / "data" / "experimental" / "nist_c4f6_mass_spectrum"
NIST = NIST_DIRECTORY / "electron_ionization_sticks.csv"
NIST_MANIFEST = NIST_DIRECTORY / "digitization_manifest.json"
BENCK_DIRECTORY = ROOT / "data" / "experimental" / "benck_2003_c4f6"
BENCK_FEED = BENCK_DIRECTORY / "figure9_mass_resolved_ion_current.csv"
BENCK_FEED_MANIFEST = BENCK_DIRECTORY / "digitization_manifest.json"
BENCK_PRESSURE = BENCK_DIRECTORY / "figure10_pressure_mass_resolved_ion_current.csv"
BENCK_PRESSURE_MANIFEST = BENCK_DIRECTORY / "figure10_digitization_manifest.json"
DEFAULT_OUTPUT = ROOT / "results" / "curated" / "c4f6_fragmentation_topology_v1"


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open(newline="", encoding="utf-8") as stream:
        return tuple(csv.DictReader(stream))


def _direct_ratios() -> dict[str, float]:
    intensity = {
        row["assignment"]: float(row["relative_intensity_percent"])
        for row in _rows(NIST)
        if row["assignment_class"].startswith("monoisotopic")
    }
    return {
        "CF2+/CF+": intensity["CF2+"] / intensity["CF+"],
        "CF3+/CF+": intensity["CF3+"] / intensity["CF+"],
        "heavy_C3F3+_plus_parent_fraction": (
            (intensity["C3F3+"] + intensity["C4F6+"])
            / sum(intensity.values())
        ),
    }


def _condition(current: dict[str, float], condition: dict[str, float]) -> dict:
    # A tuple, not a set: the resolved-current sum below must run in one
    # declared order.  Iterating a set literal made the float summation order
    # follow PYTHONHASHSEED, so the committed exact-replay board differed by
    # one ulp on some interpreter launches.
    required = ("Ar+", "CF+", "CF2+", "CF3+", "total_positive_ion_current")
    if set(current) != set(required):
        raise RuntimeError(f"incomplete Benck board at {condition}")
    resolved = sum(current[name] for name in required if name != "total_positive_ion_current")
    total = current["total_positive_ion_current"]
    if resolved > total:
        raise RuntimeError("digitized resolved current exceeds Benck total")
    return {
        **condition,
        "CF2+/CF+": current["CF2+"] / current["CF+"],
        "CF3+/CF+": current["CF3+"] / current["CF+"],
        "resolved_fraction_of_total_positive_current": resolved / total,
        "unresolved_fraction_of_total_positive_current": 1.0 - resolved / total,
    }


def _feed_conditions() -> list[dict]:
    grouped: dict[int, dict[str, float]] = {}
    for row in _rows(BENCK_FEED):
        feed = int(row["c4f6_feed_percent"])
        grouped.setdefault(feed, {})[row["species"]] = float(
            row["ion_current_density_mA_cm2"])
    results = []
    for feed, current in sorted(grouped.items()):
        # Figure 9 does not draw an Ar+ marker at the zero-Ar endpoint.  Its
        # current is physically fixed to zero by the declared feed, not
        # inferred from missing pixels.
        if feed == 100:
            current["Ar+"] = 0.0
        results.append(_condition(current, {
            "board": "Benck_Figure_9_feed_sweep",
            "c4f6_feed_percent": feed,
            "pressure_Pa": 1.33,
        }))
    return results


def _pressure_conditions() -> list[dict]:
    grouped: dict[float, dict[str, float]] = {}
    for row in _rows(BENCK_PRESSURE):
        pressure = float(row["pressure_Pa"])
        grouped.setdefault(pressure, {})[row["species"]] = float(
            row["ion_current_density_mA_cm2"])
    return [
        _condition(current, {
            "board": "Benck_Figure_10_pressure_sweep",
            "c4f6_feed_percent": 50,
            "pressure_Pa": pressure,
        })
        for pressure, current in sorted(grouped.items())
    ]


def audit() -> dict:
    direct = _direct_ratios()
    conditions = _feed_conditions() + _pressure_conditions()
    for condition in conditions:
        condition["CF2+/CF+_enhancement_over_direct_EI"] = (
            condition["CF2+/CF+"] / direct["CF2+/CF+"]
        )
        condition["CF3+/CF+_ratio_to_direct_EI"] = (
            condition["CF3+/CF+"] / direct["CF3+/CF+"]
        )

    cf2_enhancement = [
        row["CF2+/CF+_enhancement_over_direct_EI"] for row in conditions
    ]
    cf3_ratio = [row["CF3+/CF+"] for row in conditions]
    unresolved = [
        row["unresolved_fraction_of_total_positive_current"]
        for row in conditions
    ]
    return {
        "schema": "petch.c4f6_fragmentation_topology_audit.v1",
        "claim_class": "independent_shape_comparison_no_absolute_scale_transfer",
        "question": (
            "Can one-step direct electron-impact fragmentation of C4F6 explain "
            "the light positive-ion composition measured in a C4F6/Ar reactor?"
        ),
        "calibration_firewall": {
            "absolute_NIST_and_Benck_scales_equated": False,
            "reaction_parameters_fitted": False,
            "feature_depth_used": False,
            "krueger_825_nm_used": False,
        },
        "direct_70_eV_EI_reference": direct,
        "reactor_conditions": conditions,
        "diagnostics": {
            "condition_count_including_repeated_cross_figure_condition": len(conditions),
            "minimum_CF2+/CF+_enhancement_over_direct_EI": min(cf2_enhancement),
            "maximum_CF2+/CF+_enhancement_over_direct_EI": max(cf2_enhancement),
            "minimum_reactor_CF3+/CF+": min(cf3_ratio),
            "maximum_reactor_CF3+/CF+": max(cf3_ratio),
            "maximum_to_minimum_reactor_CF3+/CF+": max(cf3_ratio) / min(cf3_ratio),
            "minimum_unresolved_positive_ion_fraction": min(unresolved),
            "maximum_unresolved_positive_ion_fraction": max(unresolved),
        },
        "physics_decision": {
            "direct_fragmentation_only_is_adequate": False,
            "heavy_primary_products_must_be_retained": True,
            "secondary_fragment_ionization_required": True,
            "ion_neutral_conversion_required": True,
            "pressure_dependent_residence_and_wall_loss_required": True,
            "smallest_authorized_next_model": (
                "An atom/charge-conserving C4F6/Ar global model retaining C4F6, "
                "C3F3, CF, CF2, CF3 and their measured positive ions, with "
                "secondary CFx electron ionization, source-backed ion-neutral "
                "conversion, Bohm/wall losses, and gas residence time."
            ),
        },
        "certification": {
            "supports_required_reaction_topology": True,
            "supports_absolute_branching_fraction": False,
            "supports_absolute_reactor_flux": False,
            "supports_krueger_boundary": False,
            "supports_feature_depth": False,
        },
        "sources": {
            "nist_spectrum_csv_sha256": _sha(NIST),
            "nist_manifest_sha256": _sha(NIST_MANIFEST),
            "benck_feed_csv_sha256": _sha(BENCK_FEED),
            "benck_feed_manifest_sha256": _sha(BENCK_FEED_MANIFEST),
            "benck_pressure_csv_sha256": _sha(BENCK_PRESSURE),
            "benck_pressure_manifest_sha256": _sha(BENCK_PRESSURE_MANIFEST),
        },
    }


def _write(result: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    direct = result["direct_70_eV_EI_reference"]
    diagnostics = result["diagnostics"]
    report = f"""# C4F6 fragmentation-topology audit

This audit compares shapes, not absolute scales. The NIST 70 eV electron-
ionization spectrum and Benck's mass-resolved reactor currents are independent
observables from different apparatus; no intensity is transplanted as a flux.

Direct C4F6 EI has `CF2+/CF+ = {direct['CF2+/CF+']:.4f}`. Across Benck's
feed and pressure boards, the reactor value is enhanced by
`{diagnostics['minimum_CF2+/CF+_enhancement_over_direct_EI']:.2f}--{diagnostics['maximum_CF2+/CF+_enhancement_over_direct_EI']:.2f}x`.
The measured `CF3+/CF+` ratio also moves by
`{diagnostics['maximum_to_minimum_reactor_CF3+/CF+']:.2f}x` across conditions,
while as much as `{100.0 * diagnostics['maximum_unresolved_positive_ion_fraction']:.1f}%`
of total positive current lies outside Ar+, CF+, CF2+, and CF3+.

Therefore a one-collision light-fragment map is structurally rejected. The
next reactor must retain the heavy direct products and resolve secondary CFx
ionization, ion-neutral conversion, pressure-dependent residence time, and
wall/Bohm losses. This result authorizes that topology only. It does not
provide an absolute reactor flux, a Krueger boundary, or feature depth, and no
depth target was used.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit()
    _write(result, args.output)
    print(json.dumps({
        "diagnostics": result["diagnostics"],
        "physics_decision": result["physics_decision"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
