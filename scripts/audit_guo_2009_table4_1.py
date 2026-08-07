#!/usr/bin/env python3
"""Audit the visually transcribed Guo 2009 Table 4.1 source board."""
from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path
import re


ELEMENTS = ("Si", "C", "O", "F")
FORMULA_TOKEN = re.compile(r"(Si|C|O|F)([0-9]*)")
EXPECTED_THRESHOLDS = {
    5: 44.0,
    6: 60.0,
    7: 72.0,
    8: 55.0,
    9: 22.0,
    10: 20.0,
    11: 10.0,
    12: 20.0,
    13: 35.0,
    14: 14.0,
    15: 0.0,
    16: 40.0,
    17: 20.0,
}


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _formula_counts(formula: str) -> dict[str, int]:
    if not formula:
        return {}
    counts = {element: 0 for element in ELEMENTS}
    position = 0
    for match in FORMULA_TOKEN.finditer(formula):
        if match.start() != position:
            raise ValueError(f"unparsed formula suffix in {formula!r}")
        element, count = match.groups()
        counts[element] += int(count or "1")
        position = match.end()
    if position != len(formula):
        raise ValueError(f"unparsed formula suffix in {formula!r}")
    return {element: count for element, count in counts.items() if count}


def build_audit(root: Path | None = None) -> dict:
    root = root or Path(__file__).resolve().parents[1]
    source = root / "data" / "surface_interactions" / "guo_2009"
    csv_path = source / "table4_1_reaction_deck.csv"
    manifest_path = source / "source_manifest.json"
    rows = list(csv.DictReader(csv_path.open()))
    manifest = json.loads(manifest_path.read_text())
    csv_sha256 = _sha(csv_path)
    csv_hash_matches_manifest = (
        csv_sha256 == manifest["transcription"]["csv_sha256"])

    ids = [int(row["reaction_id"]) for row in rows]
    threshold_failures = []
    balance = []
    for row in rows:
        reaction_id = int(row["reaction_id"])
        if reaction_id in EXPECTED_THRESHOLDS:
            observed = float(row["threshold_eV"])
            expected = EXPECTED_THRESHOLDS[reaction_id]
            if observed != expected:
                threshold_failures.append({
                    "reaction_id": reaction_id,
                    "expected_eV": expected,
                    "observed_eV": observed,
                })
        if row["balance_class"] in {"exact", "exact_state_transition"}:
            reactants = _formula_counts(row["reactant_formula"])
            products = _formula_counts(row["product_formula"])
            balance.append({
                "reaction_id": reaction_id,
                "reactants": reactants,
                "products": products,
                "closed": reactants == products,
            })

    inconsistencies = manifest["transcription"]["source_inconsistencies"]
    payload = {
        "status": "passed",
        "claim": (
            "Source-table and elemental-accounting audit only; this is not "
            "independent validation, an atomistic potential, or a Krueger "
            "absolute-depth fit."
        ),
        "source_pdf_sha256": manifest["source"]["pdf_sha256"],
        "table_csv_sha256": csv_sha256,
        "table_csv_hash_matches_manifest": csv_hash_matches_manifest,
        "manifest_sha256": _sha(manifest_path),
        "row_count": len(rows),
        "sequential_reaction_ids": ids == list(range(1, 21)),
        "thresholds_match_visual_transcription": not threshold_failures,
        "threshold_failures": threshold_failures,
        "atom_counted_reactions": balance,
        "all_atom_counted_reactions_close": all(
            item["closed"] for item in balance),
        "coefficient_semantics": (
            "kinetic rate-expression coefficients, not bounded sticking "
            "probabilities"
        ),
        "coefficients_greater_than_one": {
            row["coefficient_symbol"]: float(row["coefficient_value"])
            for row in rows
            if float(row["coefficient_value"]) > 1.0
        },
        "source_inconsistencies_preserved": inconsistencies,
        "visual_audit": manifest["visual_audit"],
        "evidence_ceiling": manifest["transcription"]["evidence_ceiling"],
    }
    if (
        len(rows) != 20
        or not payload["sequential_reaction_ids"]
        or not csv_hash_matches_manifest
        or threshold_failures
        or not payload["all_atom_counted_reactions_close"]
        or len(inconsistencies) != 2
    ):
        payload["status"] = "failed"
    return payload


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "results" / "curated" / "guo_2009_table4_1"
    output.mkdir(parents=True, exist_ok=True)
    audit = build_audit(root)
    (output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if audit["status"] != "passed":
        raise SystemExit("Guo Table 4.1 audit failed")


if __name__ == "__main__":
    main()
