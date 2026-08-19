#!/usr/bin/env python3
"""Freeze Oxford reactor composition and expose daughter-EEDF feedback.

The accepted 60--120 W states were solved with energy-resolved collisions for
CHF3, SF6, and O2 only.  This audit holds each conserved heavy-particle state
fixed, then repeats the electron Boltzmann solve after adding the dominant HF
and F2 collision targets.  Electric field and RF frequency are preserved in
dimensional units as the represented neutral basis changes.

This is the correct diagnostic before modifying the nonlinear reactor.  It is
not itself a reclosed reactor state, wafer flux, TiO2 depth, or SEM prediction.
The local LXCat exports are hash locked and never packaged in the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.reactor_global.electron_collision_deck import (
    parse_bolsig_lxcat_bytes,
)
from petch.reactor_global.electron_collision_mixture import (
    compose_electron_collision_decks,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    TwoTermBoltzmannCondition,
)
from petch.reactor_global.zhu_daughter_electron_collisions import (
    deconvolve_siglo_f2_effective_momentum,
    derive_huang_2020_partial_hf_replay,
)
from petch.reactor_global.zhu_parent_collision_chemistry import (
    build_zhu_parent_collision_chemistry,
)
from petch.reactor_global.zhu_supplemental_chemistry import (
    build_zhu_supplemental_chemistry,
)
from scripts.run_zhu_open_reactor import _grid


ENSEMBLE_DIR = (
    ROOT / "results" / "curated"
    / "zhu_npg80_absorbed_power_ensemble_v1"
)
STATE_PATHS = {
    60: ENSEMBLE_DIR / "power_60W.json",
    90: (
        ROOT / "results" / "curated" / "zhu_npg80_sheath_coupled_v1"
        / "central_276V.json"
    ),
    105: ENSEMBLE_DIR / "power_105W.json",
    120: ENSEMBLE_DIR / "power_120W.json",
}
OUTPUT_DIR = (
    ROOT / "results" / "curated"
    / "zhu_npg80_daughter_eedf_sensitivity_v1"
)
OUTPUT = OUTPUT_DIR / "audit.json"
RF_FREQUENCY_HZ = 13.56e6
PARENT_TARGETS = ("CHF3", "SF6", "O2")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_source_decks(hcl_path: Path, f2_path: Path):
    hcl = parse_bolsig_lxcat_bytes(
        hcl_path.read_bytes(),
        source_database="Hayashi database",
        retrieved_at="2026-08-18",
        source_reference=(
            "local user-supplied/secondary-mirror LXCat export; official "
            "database reference https://www.lxcat.net/Hayashi"
        ),
        target="HCl",
        database_filter="Hayashi database",
    )
    f2 = parse_bolsig_lxcat_bytes(
        f2_path.read_bytes(),
        source_database="SIGLO database",
        retrieved_at="2026-08-18",
        source_reference=(
            "local user-supplied/secondary-mirror LXCat export; official "
            "database reference https://www.lxcat.net/SIGLO"
        ),
        target="F2",
        database_filter="SIGLO database",
    )
    return (
        hcl,
        f2,
        derive_huang_2020_partial_hf_replay(hcl),
        deconvolve_siglo_f2_effective_momentum(f2),
    )


def _neutral_names() -> tuple[str, ...]:
    network = build_zhu_supplemental_chemistry().network
    return tuple(sorted({
        item.name for item in network.species
        if item.role in {"neutral", "excited_neutral"}
    }))


def _collision_source_rows(solution, fractions: dict[str, float]) -> list[dict]:
    rows = []
    for moment in solution.collision_moments:
        if moment.process_kind not in {"IONIZATION", "ATTACHMENT"}:
            continue
        electron_change = 1 if moment.process_kind == "IONIZATION" else -1
        coefficient = float(moment.rate_coefficient_m3_s)
        rows.append({
            "target": moment.target,
            "kind": moment.process_kind,
            "product": moment.product,
            "rate_coefficient_m3_s": coefficient,
            "target_fraction_weighted_electron_source_m3_s": (
                electron_change * fractions[moment.target] * coefficient
            ),
        })
    return rows


def _solve_variant(
    *,
    deck,
    target_names: tuple[str, ...],
    payload: dict,
    total_neutral_density_m3: float,
) -> dict:
    densities = payload["state"]["densities_m3"]
    represented_density = sum(float(densities[name]) for name in target_names)
    basis_fraction = represented_density / total_neutral_density_m3
    fractions = {
        name: float(densities[name]) / represented_density
        for name in target_names
    }
    total_field = float(payload["state"][
        "implied_total_neutral_reduced_electric_field_Td"
    ])
    represented_field = total_field / basis_fraction
    reduced_frequency = (
        2.0 * math.pi * RF_FREQUENCY_HZ / represented_density
    )
    solution = DeterministicTwoTermBoltzmannSolver(
        _grid(deck), deck
    ).solve(
        TwoTermBoltzmannCondition(
            reduced_electric_field_Td=represented_field,
            gas_temperature_K=float(payload["input"]["gas_temperature_K"]),
            target_mole_fractions=fractions,
            growth_model="temporal_growth",
            initial_electron_temperature_eV=max(
                0.3,
                2.0 / 3.0 * float(
                    payload["state"]["mean_electron_energy_eV"]
                ),
            ),
            angular_field_frequency_over_density_m3_s=reduced_frequency,
        ),
        relative_tolerance=1.0e-8,
        maximum_iterations=220,
        maximum_tail_population_fraction=2.0e-6,
    )
    source_rows = _collision_source_rows(solution, fractions)
    reconstructed_growth = sum(
        row["target_fraction_weighted_electron_source_m3_s"]
        for row in source_rows
    )
    field_V_m = total_field * 1.0e-21 * total_neutral_density_m3
    represented_field_V_m = (
        represented_field * 1.0e-21 * represented_density
    )
    if abs(represented_field_V_m / field_V_m - 1.0) > 2.0e-14:
        raise ValueError("dimensional electric field changed across basis")
    return {
        "collision_targets": list(target_names),
        "collision_deck_sha256": deck.payload_sha256,
        "represented_neutral_density_m3": represented_density,
        "represented_neutral_fraction": basis_fraction,
        "target_mole_fractions_within_basis": fractions,
        "represented_neutral_reduced_field_Td": represented_field,
        "total_neutral_reduced_field_Td": total_field,
        "dimensional_bulk_field_V_m": field_V_m,
        "angular_frequency_over_represented_density_m3_s": reduced_frequency,
        "mean_electron_energy_eV": float(
            solution.distribution.mean_energy_eV
        ),
        "net_growth_rate_coefficient_m3_s": float(
            solution.net_growth_rate_coefficient_m3_s
        ),
        "reconstructed_net_growth_rate_coefficient_m3_s": (
            reconstructed_growth
        ),
        "growth_reconstruction_error_m3_s": (
            reconstructed_growth
            - solution.net_growth_rate_coefficient_m3_s
        ),
        "flux_reduced_mobility_m_inv_V_inv_s_inv": float(
            solution.transport_moments
            .flux_reduced_mobility_m_inv_V_inv_s_inv
        ),
        "reduced_field_power_gain_eV_m3_s": float(
            solution.transport_moments.reduced_field_power_gain_eV_m3_s
        ),
        "iteration_count": solution.iteration_count,
        "weighted_iteration_residual": float(
            solution.weighted_iteration_residual
        ),
        "electron_creation_and_attachment": source_rows,
    }


def _relative(candidate: dict, reference: dict, key: str) -> float:
    return float(candidate[key] / reference[key] - 1.0)


def audit(
    source_workbook: Path,
    hcl_path: Path,
    f2_path: Path,
) -> dict:
    hcl, f2, hf_replay, f2_replay = _load_source_decks(
        hcl_path, f2_path
    )
    parent = build_zhu_parent_collision_chemistry(source_workbook)
    parent_hf = compose_electron_collision_decks(
        (parent.mixed_deck, hf_replay.derived_deck),
        retrieved_at="2026-08-18",
        mixture_name="Zhu frozen state parents plus partial HF",
    )
    parent_hf_f2 = compose_electron_collision_decks(
        (parent.mixed_deck, hf_replay.derived_deck, f2_replay.derived_deck),
        retrieved_at="2026-08-18",
        mixture_name="Zhu frozen state parents plus partial HF and F2",
    )
    variants = (
        ("parent_only_replay", parent.mixed_deck, PARENT_TARGETS),
        ("partial_hf", parent_hf, (*PARENT_TARGETS, "HF")),
        ("partial_hf_plus_f2", parent_hf_f2, (*PARENT_TARGETS, "HF", "F2")),
    )
    neutral_names = _neutral_names()
    board = []
    for power_W, path in sorted(STATE_PATHS.items()):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload["input"]["feature_or_sem_target_used"]:
            raise ValueError("feature target entered daughter-EEDF audit")
        densities = payload["state"]["densities_m3"]
        total_neutral_density = sum(
            float(densities[name]) for name in neutral_names
        )
        predictions = {
            name: _solve_variant(
                deck=deck,
                target_names=targets,
                payload=payload,
                total_neutral_density_m3=total_neutral_density,
            )
            for name, deck, targets in variants
        }
        baseline = predictions["parent_only_replay"]
        stored_mean = float(payload["state"]["mean_electron_energy_eV"])
        if abs(baseline["mean_electron_energy_eV"] / stored_mean - 1) > 2e-13:
            raise ValueError("parent-only frozen EEDF does not replay state")
        for name, prediction in predictions.items():
            prediction["relative_to_parent_only"] = {
                "mean_electron_energy": _relative(
                    prediction, baseline, "mean_electron_energy_eV"
                ),
                "flux_reduced_mobility": _relative(
                    prediction,
                    baseline,
                    "flux_reduced_mobility_m_inv_V_inv_s_inv",
                ),
                "reduced_field_power_gain": _relative(
                    prediction,
                    baseline,
                    "reduced_field_power_gain_eV_m3_s",
                ),
            }
        board.append({
            "absorbed_power_W": power_W,
            "frozen_state_path": str(path.relative_to(ROOT)),
            "frozen_state_sha256": _sha256(path),
            "total_neutral_density_m3": total_neutral_density,
            "stored_parent_only_mean_energy_eV": stored_mean,
            "variants": predictions,
        })

    hf_changes = [
        row["variants"]["partial_hf"]["relative_to_parent_only"]
        ["mean_electron_energy"]
        for row in board
    ]
    f2_increment = [
        row["variants"]["partial_hf_plus_f2"]["mean_electron_energy_eV"]
        / row["variants"]["partial_hf"]["mean_electron_energy_eV"] - 1.0
        for row in board
    ]
    growth_signs = {
        name: [
            math.copysign(1.0, row["variants"][name][
                "net_growth_rate_coefficient_m3_s"
            ])
            for row in board
        ]
        for name in ("parent_only_replay", "partial_hf", "partial_hf_plus_f2")
    }
    return {
        "schema": "petch.zhu-npg80-daughter-eedf-sensitivity.v1",
        "condition_id": "zhu-2026-npg80-tio2-chf3-sf6-o2-20min",
        "sem_or_depth_target_used": False,
        "source_inputs": {
            "o2_workbook_sha256": _sha256(source_workbook),
            "hcl_lxcat_payload_sha256": hcl.payload_sha256,
            "f2_lxcat_payload_sha256": f2.payload_sha256,
            "hcl_source_database": hcl.source_database,
            "f2_source_database": f2.source_database,
            "raw_lxcat_bytes_committed": False,
            "local_lxcat_acquisition_class": (
                "secondary mirror inspected locally; official LXCat database "
                "references retained; replace with an official export and "
                "hash-compare before certification"
            ),
        },
        "derived_inputs": {
            "partial_hf_deck_sha256": hf_replay.derived_deck.payload_sha256,
            "partial_hf_evidence_class": hf_replay.evidence_class,
            "omitted_hf_channels": list(hf_replay.omitted_hf_channels),
            "complete_hf_eedf": hf_replay.supports_complete_hf_eedf,
            "f2_deconvolved_deck_sha256": (
                f2_replay.derived_deck.payload_sha256
            ),
            "f2_minimum_elastic_cross_section_m2": (
                f2_replay.minimum_elastic_cross_section_m2
            ),
            "f2_working_domain_eV": [0.0, f2_replay.maximum_energy_eV],
        },
        "method": {
            "heavy_particle_composition_frozen": True,
            "dimensional_electric_field_preserved": True,
            "rf_angular_frequency_preserved": True,
            "represented_reduced_field_rescaled_with_collision_basis": True,
            "nonlinear_reactor_reclosed": False,
        },
        "power_board": board,
        "finding": {
            "partial_hf_mean_energy_relative_change_range": [
                min(hf_changes), max(hf_changes)
            ],
            "f2_incremental_mean_energy_relative_change_range": [
                min(f2_increment), max(f2_increment)
            ],
            "growth_sign_by_variant": growth_signs,
            "daughter_collisions_material_to_eedf": (
                max(abs(value) for value in hf_changes) > 0.1
            ),
            "nonlinear_reclose_required": True,
            "interpretation": (
                "HF momentum and threshold-shifted losses materially cool "
                "the frozen EEDF. F2 changes mean energy weakly at the frozen "
                "composition but its attachment changes electron balance. "
                "The parent-only heavy state is therefore not a self-"
                "consistent daughter-collision solution."
            ),
        },
        "certification": {
            "supports_hf_f2_eedf_sensitivity": True,
            "supports_complete_hf_collision_basis": False,
            "supports_reclosed_reactor_state": False,
            "supports_wafer_flux_prediction": False,
            "supports_feature_depth_prediction": False,
        },
    }


def _check(path: Path) -> None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    board = payload.get("power_board", [])
    finding = payload.get("finding", {})
    if (
        payload.get("schema")
        != "petch.zhu-npg80-daughter-eedf-sensitivity.v1"
        or payload.get("sem_or_depth_target_used") is not False
        or len(board) != 4
        or [row.get("absorbed_power_W") for row in board]
        != [60, 90, 105, 120]
        or finding.get("daughter_collisions_material_to_eedf") is not True
        or finding.get("nonlinear_reclose_required") is not True
        or payload.get("certification", {}).get(
            "supports_feature_depth_prediction"
        ) is not False
    ):
        raise RuntimeError("committed daughter-EEDF sensitivity is stale")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-workbook", type=Path)
    parser.add_argument("--hcl-lxcat", type=Path)
    parser.add_argument("--f2-lxcat", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        _check(args.output)
        return
    if any(value is None for value in (
        args.source_workbook, args.hcl_lxcat, args.f2_lxcat
    )):
        parser.error(
            "--source-workbook, --hcl-lxcat, and --f2-lxcat are required"
        )
    result = audit(
        args.source_workbook, args.hcl_lxcat, args.f2_lxcat
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
