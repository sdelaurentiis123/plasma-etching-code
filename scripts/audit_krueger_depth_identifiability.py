#!/usr/bin/env python3
"""Reproduce the Krueger absolute-depth identifiability audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from petch.depth_identifiability import krueger_2024_depth_identifiability
from petch.experimental_data import (
    load_karahashi_2007_reactive_ion_yields,
    load_takada_2005_coincidence_yields,
)


ROOT = Path(__file__).resolve().parents[1]
KARAHASHI = (
    ROOT / "data" / "experimental" / "karahashi_2007"
    / "figure4_reactive_ion_yields.csv")
TAKADA = (
    ROOT / "data" / "experimental" / "takada_2005"
    / "figure3_sio2_coincidence_yields.csv")
SCORECARD = (
    ROOT / "results" / "curated" / "scorecard_endpoint" / "sc-base"
    / "audit.json")
OUTPUT = (
    ROOT / "results" / "curated" / "depth_identifiability"
    / "audit.json")


def build_audit():
    karahashi = load_karahashi_2007_reactive_ion_yields(KARAHASHI)
    takada = load_takada_2005_coincidence_yields(TAKADA)
    scorecard = json.loads(SCORECARD.read_text(encoding="utf-8"))
    return krueger_2024_depth_identifiability(
        simulated_depth_nm=scorecard["final_metrics"]["etch_depth_nm"],
        karahashi_pure_ion_peak_yield=max(
            row.yield_sio2_per_ion for row in karahashi
            if row.species == "CF3+"),
        takada_400eV_peak_yield=max(
            row.yield_sio2_per_ar_ion for row in takada
            if row.coincident_species == "C5F8"),
    )


def audit_text():
    return json.dumps(build_audit(), indent=2) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    payload = audit_text()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(payload, encoding="utf-8")
    elif OUTPUT.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed depth-identifiability audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
