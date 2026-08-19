#!/usr/bin/env python3
"""Freeze the Oxford TiO2 surface-deck readiness ledger."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from petch.tio2_surface_deck import (
    TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE,
    TIO2_TARGET_MODEL_FORM_GAPS,
)


ROOT = Path(__file__).resolve().parents[1]
TOPOLOGY = (
    ROOT / "results" / "curated" / "zhu_npg80_tio2_surface_topology_v1"
    / "audit.json"
)
CHOI = (
    ROOT / "data" / "experimental" / "choi_2013_tio2_cf4"
    / "bias_response_audit.json"
)
DEPLA = (
    ROOT / "data" / "experimental" / "depla_2024_tio2_sputter"
    / "digitization_manifest.json"
)
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_tio2_surface_deck_v1"
    / "audit.json"
)


def _receipt(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }


def build() -> dict[str, object]:
    topology = json.loads(TOPOLOGY.read_text(encoding="utf-8"))
    if topology["supports_absolute_oxford_profile_prediction"] is not False:
        raise RuntimeError("surface topology unexpectedly became predictive")

    slots = {
        "site_density_m2": {
            "numerical_target_value_identified": False,
            "available_evidence": "cross-material model assumptions only",
        },
        "bulk_formula_density_m3": {
            "numerical_target_value_identified": False,
            "available_evidence": (
                "published cross-process ALD TiO2 mass-density bracket "
                "3.25--4.15 g/cm3"
            ),
            "valid_current_use": "clearance and sensitivity interval",
        },
        "polymer_monolayer_density_m2": {
            "numerical_target_value_identified": False,
            "available_evidence": "none for the Oxford retained layer",
        },
        "polymer_bulk_unit_density_m3": {
            "numerical_target_value_identified": False,
            "available_evidence": "none for the Oxford retained layer",
        },
        "complex_formation_probability": {
            "numerical_target_value_identified": False,
            "available_evidence": "Choi XPS fixes Ti-F topology and response sign only",
        },
        "polymer_deposition_probability_on_substrate": {
            "numerical_target_value_identified": False,
            "available_evidence": "Ji morphology fixes positive retained-volume topology only",
        },
        "polymer_deposition_probability_on_polymer": {
            "numerical_target_value_identified": False,
            "available_evidence": "none for the Oxford feed and film",
        },
        "oxygen_polymer_etch_probability": {
            "numerical_target_value_identified": False,
            "available_evidence": "Choi nonmonotonic O2 response; coefficient not isolated",
        },
        "oxygen_blocking_probability": {
            "numerical_target_value_identified": False,
            "available_evidence": "Choi excess-O blocking interpretation; coefficient absent",
        },
        "oxygen_blocker_ion_removal_yield": {
            "numerical_target_value_identified": False,
            "available_evidence": "no state-resolved beam or same-tool measurement",
        },
        "bare_sio2_yield": {
            "meaning_in_tio2_deck": "bare TiO2 energetic removal yield",
            "numerical_target_value_identified": False,
            "available_evidence": (
                "Depla/Van Bever semi-empirical bare Ar->TiO2 reference "
                "curve plus Choi reactive bias response; neither identifies "
                "the Oxford mixed-ion fluorinated-surface target law"
            ),
            "reference_formula_unit_yield_at_276_eV": 0.192143,
            "reference_evidence_class": "digitized_semiempirical_fit_curve",
            "reference_transferable_as_target_coefficient": False,
        },
        "complex_sio2_yield": {
            "meaning_in_tio2_deck": "fluorinated TiO2 energetic removal yield",
            "numerical_target_value_identified": False,
            "available_evidence": "Choi bias response without state-resolved yield",
        },
        "polymer_sputter_yield": {
            "numerical_target_value_identified": False,
            "available_evidence": "Ji RF morphology response without IEAD or film inventory",
        },
    }
    if set(slots) != set(TIO2_REDUCED_SURFACE_REQUIRED_EVIDENCE):
        raise RuntimeError("TiO2 surface-deck ledger drifted from executable contract")

    missing = sorted(
        name for name, status in slots.items()
        if not status["numerical_target_value_identified"]
    )
    return {
        "schema": "petch.zhu-npg80-tio2-surface-deck.v1",
        "condition_id": topology["condition_id"],
        "target_sem_used": False,
        "target_depth_used": False,
        "inputs": {
            "surface_topology": _receipt(TOPOLOGY),
            "choi_multiaxis_response": _receipt(CHOI),
            "depla_bare_tio2_ar_sputter_reference": _receipt(DEPLA),
        },
        "executable_contract": {
            "class": "petch.tio2_surface_deck.Tio2ReducedSurfaceDeck",
            "material_name": "ALD TiO2",
            "inventory_name": "TiO2_formula_unit",
            "silent_sio2_coefficient_transfer_allowed": False,
            "numerical_defaults_supplied": False,
            "reduced_sensitivity_execution_available": True,
            "competitive_oxygen_state_implemented": True,
        },
        "parameter_slots": slots,
        "target_parameters_not_identified": missing,
        "unresolved_model_form": list(TIO2_TARGET_MODEL_FORM_GAPS),
        "current_best_use": [
            "propagate sourced density and reactor-boundary intervals",
            "falsify candidate surface laws against all Choi and Ji response axes",
            "freeze a blind Oxford sensitivity ensemble without fitting the target SEM",
        ],
        "supports_absolute_oxford_profile_prediction": False,
        "supports_atomic_accuracy": False,
        "next_closure_observables": [
            "same-run achieved electrode/DC-bias waveform",
            "same-run blanket TiO2 thickness loss",
            "same-run Cr thickness loss",
            "one O2 or bias perturbation around the target recipe",
            "cross-section/top-down SEM answer key after prediction freeze",
        ],
    }


def _render(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build())
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
        return
    if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("TiO2 surface-deck audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
