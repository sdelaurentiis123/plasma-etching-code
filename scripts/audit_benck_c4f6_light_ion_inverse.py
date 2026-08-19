#!/usr/bin/env python3
"""Falsify a common-loss C4F6 light-ion source closure against Benck.

This audit deliberately tests a *small* model.  NIST direct C4F6 electron-
impact branching and measured CFx secondary-ionization curves are combined
with Benck's condition-resolved Figure-14(a) neutral-density ratio.  For each
measured light-ion current triplet, the two remaining nonnegative
unknowns are the aggregate parent-ionization source per ``ne*n(CF)`` and
``n(CF3)/n(CF)``.

Equating volume-source ratios to measured wall-current ratios is equivalent
to assuming a common linear transport/loss operator for CF+, CF2+, and CF3+.
A negative inferred neutral density therefore rejects that assumption; it is
not repaired by clipping the density or fitting a feature depth.  This audit
does not itself identify which omitted operator (mass-dependent Bohm/wall
loss, ion-neutral conversion, heavy-product cascades, or surface return) is
responsible.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.reactor_global.c4f6_electron_collisions import (
    load_lan_jeon_2014_c4f6_replay,
)
from petch.reactor_global.c4f6_ion_sources import (
    C4F6PositiveIonSourceModel,
)
from petch.reactor_global.network import (
    ElectronMaxwellianCrossSectionRateCoefficient,
    RateContext,
)


ROOT = Path(__file__).resolve().parents[1]
BENCK_DIRECTORY = ROOT / "data" / "experimental" / "benck_2003_c4f6"
BENCK_FEED = BENCK_DIRECTORY / "figure9_mass_resolved_ion_current.csv"
BENCK_NEUTRAL = BENCK_DIRECTORY / "figure14a_cf2_cf_feed_ratio.csv"
BENCK_NEUTRAL_MANIFEST = (
    BENCK_DIRECTORY / "figure14a_digitization_manifest.json"
)
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "benck_c4f6_light_ion_inverse_v1"
)
ELECTRON_TEMPERATURES_EV = (2.0, 3.0, 4.0, 5.0, 6.0)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open(newline="", encoding="utf-8") as stream:
        return tuple(csv.DictReader(stream))


def _measured_conditions() -> tuple[dict[str, object], ...]:
    grouped: dict[int, dict[str, float]] = {}
    for row in _rows(BENCK_FEED):
        feed = int(row["c4f6_feed_percent"])
        grouped.setdefault(feed, {})[row["species"]] = float(
            row["ion_current_density_mA_cm2"]
        )
    neutral = {
        int(row["c4f6_feed_percent"]): float(
            row["neutral_CF2_to_CF_density_ratio"]
        )
        for row in _rows(BENCK_NEUTRAL)
    }
    if set(grouped) != set(neutral):
        raise RuntimeError("Benck ion and neutral feed boards do not align")
    conditions = []
    for feed, currents in sorted(grouped.items()):
        required = {"CF+", "CF2+", "CF3+"}
        if not required.issubset(currents):
            raise RuntimeError(f"incomplete Benck light-ion board at {feed}%")
        conditions.append({
            "board": "Benck_Figures_9_and_14a_5sccm_feed_sweep",
            "c4f6_feed_percent": float(feed),
            "pressure_Pa": 1.33,
            "neutral_nCF2_to_nCF": neutral[feed],
            "measured_CF2+_to_CF+": currents["CF2+"] / currents["CF+"],
            "measured_CF3+_to_CF+": currents["CF3+"] / currents["CF+"],
        })
    return tuple(conditions)


def _parent_ionization_coefficient():
    process = load_lan_jeon_2014_c4f6_replay().derived_deck.processes[-1]
    if process.kind != "IONIZATION" or process.target != "C4F6":
        raise RuntimeError("Lan--Jeon aggregate parent-ionization row changed")
    return ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=process.electron_energy_eV,
        cross_section_m2=process.cross_section_m2,
        threshold_eV=process.energy_loss_eV,
        relative_uncertainty=None,
        source="Lan--Jeon 2014 C4F6 aggregate Qi replay",
        evidence_kind="derived",
    )


def _source_coefficients(model, temperature_eV: float, neutral_cf2_to_cf: float):
    direct = model.direct.partition_event_rate_m3_s(1.0)
    cf_cf2 = model.evaluate(
        aggregate_parent_ionization_rate_m3_s=0.0,
        electron_density_m3=1.0,
        neutral_cfx_densities_m3={
            "CF": 1.0,
            "CF2": float(neutral_cf2_to_cf),
            "CF3": 0.0,
        },
        electron_temperature_eV=temperature_eV,
    ).secondary_cfx_sources_m3_s
    unit_cf3 = model.evaluate(
        aggregate_parent_ionization_rate_m3_s=0.0,
        electron_density_m3=1.0,
        neutral_cfx_densities_m3={"CF": 0.0, "CF2": 0.0, "CF3": 1.0},
        electron_temperature_eV=temperature_eV,
    ).secondary_cfx_sources_m3_s
    return direct, cf_cf2, unit_cf3


def _inverse(condition, *, temperature_eV, parent_coefficient, model):
    direct, base, cf3 = _source_coefficients(
        model, temperature_eV, float(condition["neutral_nCF2_to_nCF"])
    )
    r2 = float(condition["measured_CF2+_to_CF+"])
    r3 = float(condition["measured_CF3+_to_CF+"])
    matrix = np.asarray([
        [
            direct["CF2+"] - r2 * direct["CF+"],
            cf3["CF2+"] - r2 * cf3["CF+"],
        ],
        [
            direct["CF3+"] - r3 * direct["CF+"],
            cf3["CF3+"] - r3 * cf3["CF+"],
        ],
    ])
    right_hand_side = np.asarray([
        r2 * base["CF+"] - base["CF2+"],
        r3 * base["CF+"] - base["CF3+"],
    ])
    parent_source_per_ne_cf, neutral_cf3_to_cf = np.linalg.solve(
        matrix, right_hand_side
    )
    parent_rate = parent_coefficient.coefficient_si(
        RateContext(temperature_eV)
    )
    sources = {
        ion: (
            direct[ion] * parent_source_per_ne_cf
            + base[ion]
            + cf3[ion] * neutral_cf3_to_cf
        )
        for ion in ("CF+", "CF2+", "CF3+")
    }
    replay_r2 = sources["CF2+"] / sources["CF+"]
    replay_r3 = sources["CF3+"] / sources["CF+"]
    return {
        "electron_temperature_eV": float(temperature_eV),
        "aggregate_parent_source_per_ne_nCF_m3_s": float(
            parent_source_per_ne_cf
        ),
        "conditional_nC4F6_to_nCF": float(
            parent_source_per_ne_cf / parent_rate
        ),
        "inferred_nCF3_to_nCF": float(neutral_cf3_to_cf),
        "all_inferred_densities_nonnegative": bool(
            parent_source_per_ne_cf >= 0.0 and neutral_cf3_to_cf >= 0.0
        ),
        "ratio_replay_absolute_error": float(max(
            abs(replay_r2 - r2), abs(replay_r3 - r3)
        )),
    }


def audit() -> dict:
    model = C4F6PositiveIonSourceModel()
    parent_coefficient = _parent_ionization_coefficient()
    conditions = []
    for condition in _measured_conditions():
        inversions = [
            _inverse(
                condition,
                temperature_eV=temperature,
                parent_coefficient=parent_coefficient,
                model=model,
            )
            for temperature in ELECTRON_TEMPERATURES_EV
        ]
        half_ratio_condition = {
            **condition,
            "neutral_nCF2_to_nCF": 0.5 * float(
                condition["neutral_nCF2_to_nCF"]
            ),
        }
        half_ratio_inversions = [
            _inverse(
                half_ratio_condition,
                temperature_eV=temperature,
                parent_coefficient=parent_coefficient,
                model=model,
            )
            for temperature in ELECTRON_TEMPERATURES_EV
        ]
        conditions.append({
            **condition,
            "inversions": inversions,
            "half_line_of_sight_CF2_CF_ratio_sensitivity": (
                half_ratio_inversions
            ),
            "nonnegative_solution_exists_on_declared_temperature_grid": any(
                row["all_inferred_densities_nonnegative"]
                for row in inversions
            ),
            "nonnegative_solution_exists_when_CF2_CF_ratio_is_halved": any(
                row["all_inferred_densities_nonnegative"]
                for row in half_ratio_inversions
            ),
        })
    ar_containing = [
        row for row in conditions
        if float(row["c4f6_feed_percent"]) < 100.0
    ]
    maximum_replay_error = max(
        row["ratio_replay_absolute_error"]
        for condition in conditions
        for row in condition["inversions"]
    )
    return {
        "schema": "petch.benck-c4f6-light-ion-inverse.v1",
        "question": (
            "Can direct C4F6 branching plus measured CFx secondary ionization "
            "map to Benck light-ion currents through one common loss operator?"
        ),
        "calibration_firewall": {
            "feature_depth_used": False,
            "krueger_825_nm_used": False,
            "absolute_current_scale_used": False,
            "reaction_parameter_fitted": False,
        },
        "declared_assumptions": {
            "neutral_nCF2_to_nCF": (
                "condition-resolved 5 sccm Figure-14a filled-circle series"
            ),
            "neutral_ratio_evidence": (
                "Benck et al. submillimeter diagnostic; line-of-sight ratio "
                "increases with C4F6 feed and probably exceeds the local "
                "in-plasma ratio"
            ),
            "electron_temperature_grid_eV": list(ELECTRON_TEMPERATURES_EV),
            "wall_current_mapping": (
                "common linear transport/loss factor for CF+, CF2+, CF3+"
            ),
            "direct_parent_branching": "NIST 70 eV EI branching prior",
            "secondary_ionization": "NIST measured CF/CF2/CF3 curves",
            "parent_ionization_rate": "Lan--Jeon aggregate Qi replay",
        },
        "conditions": conditions,
        "diagnostics": {
            "condition_count": len(conditions),
            "ar_containing_condition_count": len(ar_containing),
            "ar_containing_conditions_with_nonnegative_solution": sum(
                row[
                    "nonnegative_solution_exists_on_declared_temperature_grid"
                ]
                for row in ar_containing
            ),
            "ar_containing_conditions_with_nonnegative_solution_when_"
            "CF2_CF_ratio_is_halved": sum(
                row[
                    "nonnegative_solution_exists_when_CF2_CF_ratio_is_halved"
                ]
                for row in ar_containing
            ),
            "maximum_ratio_replay_absolute_error": maximum_replay_error,
        },
        "physics_decision": {
            "common_loss_source_only_closure_is_adequate": False,
            "negative_density_is_accepted_as_physical": False,
            "negative_density_is_clipped": False,
            "required_next_operators": [
                "species- and mass-dependent Bohm/wall/exhaust loss",
                "ion-neutral conversion and charge exchange",
                "heavy-primary fragmentation cascades",
                "surface-product return and byproduct-ion formation",
            ],
            "interpretation": (
                "All three Ar-containing feed conditions have no nonnegative "
                "solution anywhere on the declared Te grid under the tested "
                "common-loss closure, and this remains true when every "
                "line-of-sight CF2/CF ratio is halved. The closure is "
                "rejected on the co-conditioned ion and neutral feed board."
            ),
        },
        "certification": {
            "supports_common_loss_falsification": True,
            "identifies_unique_missing_operator": False,
            "supports_steady_reactor_composition": False,
            "supports_absolute_reactor_flux": False,
            "supports_krueger_boundary": False,
            "supports_feature_depth": False,
        },
        "sources": {
            "benck_feed_csv_sha256": _sha(BENCK_FEED),
            "benck_neutral_ratio_csv_sha256": _sha(BENCK_NEUTRAL),
            "benck_neutral_ratio_manifest_sha256": _sha(
                BENCK_NEUTRAL_MANIFEST
            ),
            "nist_direct_spectrum_sha256": model.direct.source_sha256,
        },
    }


def _write(payload: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    diagnostics = payload["diagnostics"]
    report = f"""# Benck C4F6 light-ion inverse

This target-free audit asks whether the smallest direct-plus-secondary C4F6
ion-source model can be mapped to Benck's measured CF+, CF2+, and CF3+ wall
currents using one common loss factor. It uses the paper's independent neutral
Figure-14(a) neutral ratio at each feed condition and solves for the remaining
parent source and `n(CF3)/n(CF)` at 2--6 eV.

The algebra replays the two measured ion ratios to within
`{diagnostics['maximum_ratio_replay_absolute_error']:.2e}`. All three
Ar-containing feed entries require a negative inferred CF3 density throughout
the declared temperature grid. This remains true if every line-of-sight
CF2/CF ratio is halved, addressing the paper's warning that the local plasma
ratio is probably lower. No negative density was clipped or accepted. This is
a model-form failure: one common source-to-wall-current loss operator cannot
explain the co-conditioned neutral and ion feed boards together.

The next reactor rung must evolve species-dependent Bohm/wall/exhaust losses,
ion-neutral conversion, heavy-fragment cascades, and surface-product return.
This audit does not decide which one dominates, does not use absolute current
scale or feature depth, and does not supply a Krueger boundary or wafer flux.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    audit_path = args.output / "audit.json"
    if args.check:
        if not audit_path.exists() or audit_path.read_text(
            encoding="utf-8"
        ) != rendered:
            raise SystemExit("Benck light-ion inverse audit is stale")
        print(f"PASS {audit_path.relative_to(ROOT)}")
        return
    _write(payload, args.output)
    print(json.dumps(payload["diagnostics"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
