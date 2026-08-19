#!/usr/bin/env python3
"""Bound the differential light-ion loss required by Benck's C4F6 board.

The preceding common-loss inverse is deliberately overconstrained and returns
an unphysical negative CF3 neutral density for every Ar-containing condition.
This audit removes that nonphysical degree of freedom: ``n(CF3)`` is fixed to
zero, the nonnegative parent source is determined only from CF2+/CF+, and the
remaining CF3+/CF+ discrepancy is expressed as a required *differential*
CF3+ loss relative to the common wall/exhaust loss.

No absolute current, feature depth, or Krueger result enters this calculation.
The result is a necessary effective-operator bound, not an identification of a
specific ion-neutral reaction or a transferable reactor calibration.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.reactor_global.c4f6_ion_sources import C4F6PositiveIonSourceModel
from scripts.audit_benck_c4f6_light_ion_inverse import (
    BENCK_FEED,
    BENCK_NEUTRAL,
    BENCK_NEUTRAL_MANIFEST,
    ELECTRON_TEMPERATURES_EV,
    _measured_conditions,
    _source_coefficients,
)

DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "benck_c4f6_differential_loss_v1"
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _condition_at_temperature(model, condition, temperature_eV: float) -> dict:
    direct, base, _ = _source_coefficients(
        model,
        temperature_eV,
        float(condition["neutral_nCF2_to_nCF"]),
    )
    measured_cf2_ratio = float(condition["measured_CF2+_to_CF+"])
    measured_cf3_ratio = float(condition["measured_CF3+_to_CF+"])
    denominator = (
        direct["CF2+"] - measured_cf2_ratio * direct["CF+"]
    )
    if denominator == 0.0:
        raise RuntimeError("CF2 ratio does not identify the parent source")
    parent_source_per_ne_ncf = (
        measured_cf2_ratio * base["CF+"] - base["CF2+"]
    ) / denominator
    sources = {
        ion: direct[ion] * parent_source_per_ne_ncf + base[ion]
        for ion in ("CF+", "CF2+", "CF3+")
    }
    if parent_source_per_ne_ncf < 0.0 or any(value <= 0.0 for value in sources.values()):
        raise RuntimeError("CF2-conditioned source state is not physical")
    source_cf2_ratio = sources["CF2+"] / sources["CF+"]
    source_cf3_ratio = sources["CF3+"] / sources["CF+"]
    required_differential_loss = (
        source_cf3_ratio / measured_cf3_ratio - 1.0
    )
    return {
        "electron_temperature_eV": float(temperature_eV),
        "inferred_parent_source_per_ne_nCF_m3_s": float(
            parent_source_per_ne_ncf
        ),
        "assumed_nCF3_to_nCF": 0.0,
        "all_source_terms_nonnegative": True,
        "source_CF2plus_to_CFplus": float(source_cf2_ratio),
        "measured_CF2plus_to_CFplus": measured_cf2_ratio,
        "CF2_ratio_replay_absolute_error": float(
            abs(source_cf2_ratio - measured_cf2_ratio)
        ),
        "source_CF3plus_to_CFplus_before_differential_loss": float(
            source_cf3_ratio
        ),
        "measured_CF3plus_to_CFplus": measured_cf3_ratio,
        "required_CF3plus_selective_loss_over_common_loss": float(
            required_differential_loss
        ),
        "positive_selective_CF3plus_loss_is_sufficient": bool(
            required_differential_loss >= 0.0
        ),
    }


def audit() -> dict:
    model = C4F6PositiveIonSourceModel()
    conditions = []
    for condition in _measured_conditions():
        rows = [
            _condition_at_temperature(model, condition, temperature)
            for temperature in ELECTRON_TEMPERATURES_EV
        ]
        conditions.append({
            **condition,
            "temperature_sweep": rows,
            "positive_selective_CF3plus_loss_suffices_at_every_temperature": all(
                row["positive_selective_CF3plus_loss_is_sufficient"]
                for row in rows
            ),
        })
    ar_rows = [
        row
        for condition in conditions
        if float(condition["c4f6_feed_percent"]) < 100.0
        for row in condition["temperature_sweep"]
    ]
    pure_rows = conditions[-1]["temperature_sweep"]
    required_ar_loss = [
        float(row["required_CF3plus_selective_loss_over_common_loss"])
        for row in ar_rows
    ]
    maximum_cf2_error = max(
        float(row["CF2_ratio_replay_absolute_error"])
        for condition in conditions
        for row in condition["temperature_sweep"]
    )
    return {
        "schema": "petch.benck-c4f6-differential-light-ion-loss.v1",
        "question": (
            "What species-selective CF3+ loss is required after a nonnegative "
            "parent source is fixed independently by Benck CF2+/CF+?"
        ),
        "calibration_firewall": {
            "feature_depth_used": False,
            "krueger_825_nm_used": False,
            "absolute_current_scale_used": False,
            "reaction_rate_fitted": False,
            "negative_density_used": False,
        },
        "declared_model": {
            "direct_parent_branching": "NIST 70 eV EI branching prior",
            "secondary_ionization": "NIST measured CF/CF2 curves",
            "neutral_CF2_to_CF": (
                "condition-resolved Benck Figure-14a filled-circle series"
            ),
            "neutral_CF3_to_CF": 0.0,
            "parent_source_identification": (
                "solve the measured CF2+/CF+ ratio exactly"
            ),
            "differential_operator": (
                "additional first-order CF3+ removal relative to the common "
                "light-ion wall/exhaust loss"
            ),
        },
        "conditions": conditions,
        "diagnostics": {
            "condition_count": len(conditions),
            "temperature_count_per_condition": len(ELECTRON_TEMPERATURES_EV),
            "maximum_CF2_ratio_replay_absolute_error": maximum_cf2_error,
            "minimum_required_Ar_mixture_CF3plus_selective_loss_over_common_loss": min(
                required_ar_loss
            ),
            "maximum_required_Ar_mixture_CF3plus_selective_loss_over_common_loss": max(
                required_ar_loss
            ),
            "all_Ar_mixture_rows_have_nonnegative_sources": all(
                row["all_source_terms_nonnegative"] for row in ar_rows
            ),
            "all_Ar_mixture_rows_admit_positive_CF3plus_selective_loss": all(
                row["positive_selective_CF3plus_loss_is_sufficient"]
                for row in ar_rows
            ),
            "pure_C4F6_rows_requiring_positive_loss": sum(
                row["positive_selective_CF3plus_loss_is_sufficient"]
                for row in pure_rows
            ),
            "pure_C4F6_rows_requiring_missing_CF3plus_source_or_changed_branching": sum(
                not row["positive_selective_CF3plus_loss_is_sufficient"]
                for row in pure_rows
            ),
        },
        "physics_decision": {
            "negative_CF3_density_is_required": False,
            "common_loss_model_is_reinstated": False,
            "Ar_mixture_model_form_can_be_repaired_by_positive_differential_loss": True,
            "single_condition_independent_CF3plus_loss_is_certified": False,
            "narrowed_missing_operator": (
                "condition-dependent CF3+ removal/source balance; candidate "
                "mechanisms include ion-neutral conversion, heavy-fragment "
                "cascades, and surface-product return, but this audit does not "
                "select among them"
            ),
            "interpretation": (
                "The nonphysical negative-density result was caused by forcing "
                "all light ions through one loss operator. Every Ar-containing "
                "condition is compatible with nonnegative sources and a modest "
                "positive CF3+-selective removal term. Pure C4F6 changes sign "
                "across the Te grid, proving that one fixed loss factor is still "
                "not a complete reactor model."
            ),
        },
        "certification": {
            "supports_required_effective_loss_bound": True,
            "identifies_unique_reaction": False,
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
    diagnostic = payload["diagnostics"]
    report = f"""# Benck C4F6 differential light-ion loss bound

The previous common-loss inverse produced an unphysical negative CF3 neutral
density in all three Ar-containing feed conditions. This audit fixes
`n(CF3)/n(CF) = 0`, identifies a nonnegative parent source only from the
independent measured `CF2+/CF+` ratio, and asks how much additional CF3+ loss
is then required by `CF3+/CF+`.

All 15 Ar-mixture temperature rows close with nonnegative sources. The required
CF3+-selective first-order loss is
`{diagnostic['minimum_required_Ar_mixture_CF3plus_selective_loss_over_common_loss']:.3f}`--
`{diagnostic['maximum_required_Ar_mixture_CF3plus_selective_loss_over_common_loss']:.3f}`
times the common light-ion wall/exhaust loss over 2--6 eV. The CF2 ratio replay
error is below `{diagnostic['maximum_CF2_ratio_replay_absolute_error']:.2e}`.

This repairs the *sign* of the Ar-mixture inverse without clipping a density or
using depth. It does not identify a unique reaction. The pure-C4F6 condition
changes the required operator sign across the temperature grid, so one fixed
CF3+ loss coefficient is still rejected. The next forward reactor must evolve
condition-dependent ion-neutral conversion, heavy-fragment cascades, and
surface-product return, then pass both Benck feed and pressure boards.

No absolute current, Krueger depth, or feature result was used. This receipt
does not provide a Krueger boundary or wafer flux.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path = args.output / "audit.json"
    if args.check:
        if not path.exists() or path.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed differential-loss audit is stale")
        print(f"PASS {path.relative_to(ROOT)}")
        return
    _write(payload, args.output)
    print(json.dumps(payload["diagnostics"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
