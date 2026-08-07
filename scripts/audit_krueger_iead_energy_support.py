#!/usr/bin/env python3
"""Quantify source-domain support over the Krueger Figure 4(a) IEAD.

This audit does not score an etch model.  It answers the narrower provenance
question that must precede one: what probability mass of the published
reactor IEAD lies inside each independently established surface-physics
energy domain?

The IEAD is an HPEM result digitized from Krueger Figure 4(a), not a measured
wafer distribution.  Energy overlap alone also does not resolve its missing
positive-ion composition or transfer a C4F8 beam closure to C4F6 plasma.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IEAD = (
    ROOT
    / "data"
    / "experimental"
    / "krueger_2024"
    / "digitized_figure4_iead.csv"
)
METADATA = IEAD.with_name("digitized_figure4_iead_metadata.json")
OUTPUT = (
    ROOT
    / "results"
    / "curated"
    / "krueger_iead_energy_support"
    / "audit.json"
)

IEAD_SHA256 = (
    "913d31be623ec5d52d226c8cea499e7f014cf4f5a27e017b519633c96e5e3ee3"
)
METADATA_SHA256 = (
    "7904a700afdcd116c6f57ef35aeb5555661ffed0d004304d815d20f79840ca55"
)

DOMAINS = (
    {
        "domain_id": "guo_yin_regression",
        "maximum_energy_eV": 370.0,
        "evidence": (
            "Guo/Kwon translating-layer coefficients regressed on the "
            "Yin C4F8/Ar-SiO2 beam board"
        ),
        "authority": "beam_regressed_mechanism",
        "target_match": True,
        "species_match": False,
        "chemistry_match": False,
    },
    {
        "domain_id": "an_2026_released_nnp_outputs",
        "maximum_energy_eV": 1000.0,
        "evidence": (
            "released DFT-trained NNP plus ZBL molecular-dynamics outputs "
            "for mass-selected CFx+ bombardment of SiO2"
        ),
        "authority": "atomistic_model_tested_against_direct_beam_data",
        "target_match": True,
        "species_match": False,
        "chemistry_match": False,
    },
    {
        "domain_id": "karahashi_mass_selected_beam",
        "maximum_energy_eV": 2000.0,
        "evidence": (
            "direct normal-incidence F+, CF+, CF2+, and CF3+ SiO2 beam "
            "measurements"
        ),
        "authority": "direct_beam_measurement",
        "target_match": True,
        "species_match": False,
        "chemistry_match": False,
    },
    {
        "domain_id": "tachi_1982_si_target_lead",
        "maximum_energy_eV": 3000.0,
        "evidence": (
            "mass-selected F+, CF+, CF2+, CF3+, and C+ beam study on "
            "elemental Si; abstract/metadata lead only, not landed SiO2 data"
        ),
        "authority": "target_mismatched_lead_not_surface_support",
        "target_match": False,
        "species_match": False,
        "chemistry_match": False,
    },
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [
            {
                "energy_eV": float(row["energy_eV"]),
                "signed_angle_deg": float(row["signed_angle_deg"]),
                "probability_weight": float(row["probability_weight"]),
            }
            for row in csv.DictReader(stream)
        ]
    if not rows:
        raise ValueError("empty IEAD")
    if any(row["probability_weight"] < 0.0 for row in rows):
        raise ValueError("negative IEAD probability weight")
    return rows


def build_audit(
    iead_path: Path = IEAD,
    metadata_path: Path = METADATA,
) -> dict[str, object]:
    if _sha256(iead_path) != IEAD_SHA256:
        raise ValueError("Krueger IEAD checksum changed")
    if _sha256(metadata_path) != METADATA_SHA256:
        raise ValueError("Krueger IEAD metadata checksum changed")

    rows = _load_rows(iead_path)
    total_weight = sum(row["probability_weight"] for row in rows)
    if abs(total_weight - 1.0) > 1.0e-12:
        raise ValueError(f"IEAD probability is not normalized: {total_weight}")
    mean_energy = (
        sum(row["energy_eV"] * row["probability_weight"] for row in rows)
        / total_weight
    )

    domain_rows = []
    for domain in DOMAINS:
        maximum = float(domain["maximum_energy_eV"])
        supported = (
            sum(
                row["probability_weight"]
                for row in rows
                if row["energy_eV"] <= maximum
            )
            / total_weight
        )
        domain_rows.append(
            {
                **domain,
                "iead_probability_at_or_below_maximum": supported,
                "iead_probability_above_maximum": 1.0 - supported,
                "grants_krueger_prediction_authority": False,
                "authority_limit": (
                    "Energy overlap is necessary but insufficient: the "
                    "published IEAD is aggregate and reactor-computed, C4F6 "
                    "wafer radical flux is missing, and no domain identifies "
                    "the incident positive-ion mixture."
                ),
            }
        )

    return {
        "schema": "petch.krueger-iead.energy-support.v1",
        "audit_id": "KRUEGER-IEAD-ENERGY-SUPPORT-R1",
        "status": "completed_provenance_audit",
        "question": (
            "What fraction of the published Krueger IEAD lies inside each "
            "surface-physics energy domain?"
        ),
        "source": {
            "citation": (
                "Krueger et al., J. Vac. Sci. Technol. A 42, 043008 "
                "(2024), Figure 4(a), DOI 10.1116/6.0003554"
            ),
            "iead_kind": "HPEM reactor-model output_not_wafer_measurement",
            "species_resolution": "all_positive_ions_combined",
            "iead_sha256": IEAD_SHA256,
            "metadata_sha256": METADATA_SHA256,
            "node_count": len(rows),
            "probability_sum": total_weight,
            "minimum_energy_eV": min(row["energy_eV"] for row in rows),
            "maximum_energy_eV": max(row["energy_eV"] for row in rows),
            "probability_weighted_mean_energy_eV": mean_energy,
        },
        "domains": domain_rows,
        "atomic_accuracy_verdict": {
            "granted": False,
            "reason": (
                "The released atomistic output ends at 1000 eV, below "
                "94.898% of this IEAD. The direct species-resolved SiO2 "
                "beam board ends at 2000 eV, below 86.817% of this IEAD. "
                "A depth calculation can therefore be a declared "
                "high-energy/species sensitivity, but not an atomic-level "
                "prediction, until a validated high-energy reactive-event "
                "closure and species-resolved reactor boundary exist."
            ),
        },
        "forbidden_inference": (
            "Do not select an extrapolation, ion mixture, or flux scale by "
            "agreement with Krueger's 825 nm endpoint."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = build_audit()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"stale audit: {args.output}")
        print(f"verified {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
