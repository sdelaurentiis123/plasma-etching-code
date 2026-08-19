#!/usr/bin/env python3
"""Bind independent TiO2 response boards into a fail-closed model contract."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHOI = (
    ROOT / "data" / "experimental" / "choi_2013_tio2_cf4"
    / "bias_response_audit.json"
)
JI3 = (
    ROOT / "data" / "experimental" / "ji_2024_tio2_hierarchical"
    / "figure3_digitization_manifest.json"
)
JI4 = (
    ROOT / "data" / "experimental" / "ji_2024_tio2_hierarchical"
    / "figure4_digitization_manifest.json"
)
CONDITIONAL = (
    ROOT / "results" / "curated" / "zhu_npg80_conditional_profiles_v1"
    / "audit.json"
)
OUTPUT = (
    ROOT / "results" / "curated"
    / "zhu_npg80_tio2_surface_topology_v1" / "audit.json"
)


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _receipt(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }


def build() -> dict[str, object]:
    choi = _load(CHOI)
    ji3 = _load(JI3)
    ji4 = _load(JI4)
    conditional = _load(CONDITIONAL)
    if (
        conditional["target_sem_used"] is not False
        or conditional["target_depth_used"] is not False
        or ji3["physics_use"][
            "removal_only_model_can_reproduce_strict_gap_narrowing"
        ] is not False
        or ji4["derived_checks"]["strict_threshold_identified"] is not False
        or choi["freddie_boundary"]["coefficient_transfer_allowed"] is not False
    ):
        raise RuntimeError("surface-topology evidence boundary changed")
    current = {
        "local_geometric_ion_dose": True,
        "removal_velocity_nonnegative": True,
        "ion_energy_dependent_surface_yield": False,
        "neutral_radical_surface_reactions": False,
        "fluorinated_tio2_surface_inventory": False,
        "fluorocarbon_or_tifx_passivation_inventory": False,
        "positive_surface_volume_deposition_or_retention": False,
        "cr_mask_geometry_evolution": False,
        "surface_charging_feedback": False,
        "removed_product_redeposition": False,
    }
    required = {
        "local_geometric_ion_dose": {
            "required": True,
            "evidence": "feature shadowing and AR-dependent transport",
        },
        "ion_energy_dependent_surface_yield": {
            "required": True,
            "evidence": "Choi DC-bias response",
        },
        "neutral_radical_surface_reactions": {
            "required": True,
            "evidence": "Choi Ti-F/O blocking plus Ji feed chemistry",
        },
        "fluorinated_tio2_surface_inventory": {
            "required": True,
            "evidence": "Choi XPS Ti-F formation",
        },
        "passivation_inventory_with_physical_thickness": {
            "required": True,
            "evidence": "Ji gap narrowing and reported lower-feature passivation",
        },
        "ion_assisted_passivation_removal": {
            "required": True,
            "evidence": "Ji RF response and Choi product-desorption response",
        },
        "cr_mask_geometry_evolution": {
            "required": True,
            "evidence": "Ji observed lateral Cr shrink",
        },
        "pattern_dependent_neutral_and_ion_transport": {
            "required": True,
            "evidence": "Ji spacing response",
        },
    }
    missing = [
        name for name, value in required.items()
        if value["required"] and not current.get(name, False)
    ]
    return {
        "schema": "petch.zhu-npg80-tio2-surface-topology.v1",
        "condition_id": conditional["condition_id"],
        "target_sem_used": False,
        "target_depth_used": False,
        "inputs": {
            "choi_bias_response": _receipt(CHOI),
            "ji_rf_morphology_response": _receipt(JI3),
            "ji_spacing_morphology_response": _receipt(JI4),
            "conditional_profile_board": _receipt(CONDITIONAL),
        },
        "experimental_discriminants": {
            "energy": {
                "observable": "TiO2 etch rate",
                "response": "130.9 to 197.2 nm/min as bias magnitude rises 50 to 250 V",
                "forces": [
                    "ion-energy-dependent bond breaking",
                    "ion-assisted reaction-product desorption",
                ],
            },
            "rf_morphology": {
                "observable": "upper height, corner radius, and interfeature gap",
                "response": "gap narrows 95.96 to 18.02 nm over 90 to 210 W",
                "forces": [
                    "positive retained or deposited surface volume",
                    "ion-assisted removal of that passivation",
                    "evolving mask or equivalent measured boundary",
                ],
            },
            "spacing": {
                "observable": "two-zone TiO2 morphology",
                "response": "100 and 70 nm gap points separate from 350 to 750 nm cluster",
                "forces": [
                    "pattern-dependent transport and surface-state feedback",
                ],
                "sharp_threshold_identified": False,
                "sampled_transition_bracket_nm": [100.0, 350.0],
            },
        },
        "current_conditional_mechanism": {
            "name": "RateNormalizedRemovalMechanism",
            "capabilities": current,
            "physical_within_declared_scope": True,
            "sufficient_for_target_surface_prediction": False,
            "interpretation": (
                "a converged wider lower section can be a physical consequence "
                "of slower shadowed-wall recession, but this mechanism cannot "
                "attribute the real Oxford profile"
            ),
        },
        "minimum_physical_model_contract": required,
        "missing_from_current_conditional_mechanism": missing,
        "deterministic_differentiable_state_design": {
            "per_surface_element_state": [
                "bounded fluorinated-site fraction",
                "nonnegative passivation units per square metre",
                "bounded activated-passivation fraction",
                "cumulative removed TiO2 formula units per square metre",
            ],
            "moving_materials": ["TiO2", "Cr mask", "passivation overlayer"],
            "operators": [
                "analytic bounded adsorption/activation update",
                "deterministic angular quadrature for ion and neutral delivery",
                "energy- and angle-resolved ion-assisted removal",
                "conservative surface-state remap during level-set motion",
            ],
            "monte_carlo_required": False,
            "compatible_with_automatic_differentiation": True,
        },
        "target_parameter_status": {
            "identified_from_existing_cross_process_boards": False,
            "reason": (
                "the boards identify topology and response signs but change tool, "
                "feed, film state, mask, and unreported wafer fluxes"
            ),
            "required_same_condition_observables": [
                "achieved DC self-bias or electrode voltage waveform",
                "blanket TiO2 loss over the 1200 s recipe",
                "remaining Cr thickness or Cr loss",
                "target GDS dimensions and sample radius",
                "cross-section/top-down SEM answer key",
            ],
        },
        "supports_absolute_oxford_profile_prediction": False,
        "supports_atomic_accuracy": False,
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
        raise SystemExit("TiO2 surface-topology audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
