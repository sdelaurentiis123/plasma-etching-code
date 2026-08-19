#!/usr/bin/env python3
"""Rank the missing electron-collision basis in Oxford reactor states.

The conserved Zhu reactor includes daughter-particle balances and a subset of
published daughter electron-rate regressions.  Its two-term Boltzmann operator,
however, presently contains energy-dependent cross sections only for the three
feed parents.  This audit separates those two levels and ranks the daughter
targets whose cross-section decks would close the largest fraction of neutral
pressure.  No SEM, etch depth, or surface coefficient enters the ranking.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.reactor_global.zhu_supplemental_chemistry import (
    build_zhu_supplemental_chemistry,
)


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
    / "zhu_npg80_daughter_collision_gap_v1"
)
OUTPUT = OUTPUT_DIR / "audit.json"
FULL_EEDF_TARGETS = frozenset({"CHF3", "SF6", "O2"})


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _electron_rate_targets(network) -> dict[str, dict[str, object]]:
    neutral_names = {
        item.name
        for item in network.species
        if item.role in {"neutral", "excited_neutral"}
    }
    reactions: dict[str, list[str]] = defaultdict(list)
    sources: dict[str, set[str]] = defaultdict(set)
    for reaction in network.reactions:
        if "e" not in reaction.kinetic_orders:
            continue
        for name in reaction.kinetic_orders:
            if name in neutral_names:
                reactions[name].append(reaction.name)
                sources[name].add(reaction.source)
    return {
        name: {
            "reaction_count": len(reactions[name]),
            "reaction_names": sorted(reactions[name]),
            "sources": sorted(sources[name]),
        }
        for name in reactions
    }


def _kinetic_tier(
    species: str,
    reduced_targets: dict[str, dict[str, object]],
) -> str:
    if species in FULL_EEDF_TARGETS:
        return "energy_resolved_two_term_eedf"
    if species in reduced_targets:
        return "electron_rate_regression_only"
    return "electron_transparent"


def _neutral_rows(
    payload: dict,
    *,
    neutral_names: tuple[str, ...],
    reduced_targets: dict[str, dict[str, object]],
) -> tuple[list[dict], dict[str, float]]:
    densities = payload["state"]["densities_m3"]
    total = sum(float(densities[name]) for name in neutral_names)
    rows = []
    tier_fractions = defaultdict(float)
    for name in neutral_names:
        density = float(densities[name])
        fraction = density / total
        tier = _kinetic_tier(name, reduced_targets)
        tier_fractions[tier] += fraction
        rows.append({
            "species": name,
            "density_m3": density,
            "neutral_fraction": fraction,
            "electron_kinetic_tier": tier,
            "reduced_electron_reaction_count": int(
                reduced_targets.get(name, {}).get("reaction_count", 0)
            ),
        })
    rows.sort(key=lambda row: (-row["neutral_fraction"], row["species"]))
    return rows, dict(tier_fractions)


def _closure_milestones(
    rows: list[dict], *, current_fraction: float,
) -> list[dict]:
    missing = [
        row for row in rows
        if row["electron_kinetic_tier"] != "energy_resolved_two_term_eedf"
    ]
    milestones = []
    cumulative = current_fraction
    selected = []
    thresholds = iter((0.50, 0.75, 0.90, 0.95, 0.99))
    target = next(thresholds)
    for row in missing:
        selected.append(row["species"])
        cumulative += row["neutral_fraction"]
        while cumulative >= target:
            milestones.append({
                "total_neutral_fraction_represented": cumulative,
                "minimum_ranked_daughter_targets": list(selected),
                "threshold": target,
            })
            try:
                target = next(thresholds)
            except StopIteration:
                return milestones
    return milestones


def build_receipt() -> dict:
    supplemental = build_zhu_supplemental_chemistry()
    network = supplemental.network
    neutral_names = tuple(sorted({
        item.name
        for item in network.species
        if item.role in {"neutral", "excited_neutral"}
    }))
    reduced_targets = _electron_rate_targets(network)
    power_board = []
    peak_fraction = defaultdict(float)
    rows_120 = None
    tier_120 = None
    for power_W, path in sorted(STATE_PATHS.items()):
        payload = _load(path)
        if payload["condition_id"] != (
            "zhu-2026-npg80-tio2-chf3-sf6-o2-20min"
        ):
            raise ValueError("collision-gap state belongs to another condition")
        if payload["input"]["feature_or_sem_target_used"]:
            raise ValueError("feature target entered collision-gap audit")
        rows, tier_fractions = _neutral_rows(
            payload,
            neutral_names=neutral_names,
            reduced_targets=reduced_targets,
        )
        represented = tier_fractions["energy_resolved_two_term_eedf"]
        reported = payload["state"][
            "electron_collision_basis_neutral_fraction"
        ]
        if abs(represented - reported) > 2.0e-14:
            raise ValueError("reported EEDF collision fraction does not replay")
        for row in rows:
            peak_fraction[row["species"]] = max(
                peak_fraction[row["species"]], row["neutral_fraction"]
            )
        power_board.append({
            "absorbed_power_W": power_W,
            "state_path": str(path.relative_to(ROOT)),
            "state_sha256": _sha256(path),
            "energy_resolved_two_term_eedf_fraction": represented,
            "electron_rate_regression_only_fraction": tier_fractions.get(
                "electron_rate_regression_only", 0.0
            ),
            "electron_transparent_fraction": tier_fractions.get(
                "electron_transparent", 0.0
            ),
        })
        if power_W == 120:
            rows_120 = rows
            tier_120 = tier_fractions
    if rows_120 is None or tier_120 is None:
        raise RuntimeError("120 W state missing from collision-gap audit")

    ranked = []
    for row in rows_120:
        if row["electron_kinetic_tier"] == "energy_resolved_two_term_eedf":
            continue
        enriched = dict(row)
        enriched["peak_neutral_fraction_over_power_board"] = peak_fraction[
            row["species"]
        ]
        reduced = reduced_targets.get(row["species"])
        enriched["reduced_rate_sources"] = (
            [] if reduced is None else reduced["sources"]
        )
        ranked.append(enriched)

    return {
        "schema": "petch.zhu-npg80-daughter-collision-gap.v1",
        "condition_id": "zhu-2026-npg80-tio2-chf3-sf6-o2-20min",
        "sem_or_depth_target_used": False,
        "scope": {
            "energy_resolved_two_term_eedf_targets": sorted(FULL_EEDF_TARGETS),
            "electron_rate_regression_only_definition": (
                "Published or regressed scalar electron-rate laws contribute "
                "chemistry and power but do not alter the solved EEPF momentum "
                "and energy operator as an energy-resolved collision target."
            ),
            "electron_transparent_definition": (
                "No direct electron-impact row is active for this neutral in "
                "either the two-term collision deck or reduced rate network."
            ),
        },
        "power_board": power_board,
        "accepted_120W_state": {
            "energy_resolved_two_term_eedf_fraction": tier_120[
                "energy_resolved_two_term_eedf"
            ],
            "electron_rate_regression_only_fraction": tier_120.get(
                "electron_rate_regression_only", 0.0
            ),
            "electron_transparent_fraction": tier_120.get(
                "electron_transparent", 0.0
            ),
            "ranked_missing_targets": ranked,
            "collision_basis_closure_milestones": _closure_milestones(
                rows_120,
                current_fraction=tier_120[
                    "energy_resolved_two_term_eedf"
                ],
            ),
        },
        "priority_verdict": {
            "first_target": ranked[0]["species"],
            "first_target_120W_neutral_fraction": ranked[0][
                "neutral_fraction"
            ],
            "interpretation": (
                "HF is the first kinetic-closure target by a wide margin. "
                "The ranking is a leverage audit of the current conditional "
                "state, not a converged composition: adding collisions feeds "
                "back on the EEPF and all daughter densities."
            ),
            "supports_unique_reactor_state": False,
            "supports_wafer_flux_prediction": False,
            "supports_feature_depth_prediction": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    receipt = build_receipt()
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != encoded:
            raise SystemExit("committed daughter-collision audit is stale")
        print(encoded, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded, encoding="utf-8")
    print(f"wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
