#!/usr/bin/env python3
"""Build the checksum-bound Guo/Krueger finite-fluence prefix receipt.

The four source directories are full pilot artifacts and intentionally remain
outside git.  This script distills only the convergence and conservation
evidence needed to audit the 0.5 s numerical gate.  It never extrapolates the
prefix to the 60 s experimental endpoint.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


CASE_DIRECTORIES = {
    "temporal_dt_0p25_dx_10nm": "krueger_guo_transient_dx10",
    "temporal_dt_0p125_dx_10nm": "krueger_guo_transient_dt125_dx10",
    "temporal_dt_0p0625_dx_10nm": "krueger_guo_transient_dt0625_dx10",
    "spatial_dt_0p0625_dx_5nm": "krueger_guo_transient_dt0625_dx5",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_case(path: Path) -> dict:
    audit_path = path / "audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit["status"] != "complete":
        raise ValueError(f"{audit_path} is not complete")
    final = audit["history"][-1]
    if not math.isclose(
        float(final["physical_time_s"]), 0.5, rel_tol=0.0, abs_tol=1e-14
    ):
        raise ValueError(f"{audit_path} does not end at the 0.5 s gate")
    if len(audit["topology_events"]) != 0:
        raise ValueError(f"{audit_path} contains a topology event")

    profile_path = path / "profile.png"
    return {
        "artifact_directory": str(path),
        "audit_sha256": _sha256(audit_path),
        "profile_sha256": _sha256(profile_path) if profile_path.exists() else None,
        "config_hash": audit["config_hash"],
        "implementation_revision": audit["git"]["revision"],
        "implementation_worktree_was_dirty": bool(audit["git"]["dirty"]),
        "dx_nm": 1000.0 * float(audit["configuration"]["dx_um"]),
        "nominal_step_s": (
            float(audit["configuration"]["duration_s"])
            / int(audit["configuration"]["n_steps"])
        ),
        "accepted_steps": len(audit["history"]) - 1,
        "depth_nm": float(final["metrics"]["etch_depth_nm"]),
        "mean_depth_rate_nm_s": 2.0 * float(final["metrics"]["etch_depth_nm"]),
        "mask_opening_nm": float(final["metrics"]["mask_opening_nm"]),
        "maximum_material_ledger_residual_units_m2": max(
            float(row.get("maximum_material_ledger_residual_units_m2", 0.0))
            for row in audit["history"]
        ),
        "maximum_neutral_radiosity_relative_balance_error": max(
            float(row.get("maximum_neutral_radiosity_relative_balance_error", 0.0))
            for row in audit["history"]
        ),
        "wall_time_s": float(audit["wall_time_s"]),
    }


def build_receipt(input_root: Path) -> dict:
    cases = {
        name: _load_case(input_root / directory)
        for name, directory in CASE_DIRECTORIES.items()
    }

    coarse = cases["temporal_dt_0p25_dx_10nm"]["depth_nm"]
    medium = cases["temporal_dt_0p125_dx_10nm"]["depth_nm"]
    fine = cases["temporal_dt_0p0625_dx_10nm"]["depth_nm"]
    if not (coarse < medium < fine):
        raise ValueError("temporal depth sequence is not monotone")
    observed_order = math.log((medium - coarse) / (fine - medium), 2.0)
    richardson_limit = fine + (fine - medium) / (2.0**observed_order - 1.0)
    fine_to_limit = abs(fine - richardson_limit) / richardson_limit

    spatial_fine = cases["spatial_dt_0p0625_dx_5nm"]["depth_nm"]
    spatial_relative = (spatial_fine - fine) / fine
    maximum_ledger = max(
        row["maximum_material_ledger_residual_units_m2"]
        for row in cases.values()
    )
    maximum_radiosity = max(
        row["maximum_neutral_radiosity_relative_balance_error"]
        for row in cases.values()
    )

    gates = {
        "temporal_fine_to_richardson_limit_le_1pct": {
            "limit": 0.01,
            "value": fine_to_limit,
            "passed": fine_to_limit <= 0.01,
        },
        "spatial_5nm_vs_10nm_abs_relative_le_5pct": {
            "limit": 0.05,
            "value": abs(spatial_relative),
            "signed_value": spatial_relative,
            "passed": abs(spatial_relative) <= 0.05,
        },
        "material_ledger_exact": {
            "limit": 0.0,
            "value": maximum_ledger,
            "passed": maximum_ledger == 0.0,
        },
        "neutral_radiosity_balance_le_1e_9": {
            "limit": 1.0e-9,
            "value": maximum_radiosity,
            "passed": maximum_radiosity <= 1.0e-9,
        },
    }
    if not all(row["passed"] for row in gates.values()):
        raise ValueError("one or more finite-fluence prefix gates failed")

    return {
        "schema": "petch-guo-krueger-finite-fluence-prefix-v1",
        "claim": (
            "Numerically converged 0.5 s prefix for the source-fixed finite-fluence "
            "Guo/Kwon transfer; not an absolute-depth prediction and not a 60 s "
            "extrapolation."
        ),
        "cases": cases,
        "temporal_convergence": {
            "observed_order": observed_order,
            "richardson_limit_nm": richardson_limit,
            "fine_to_limit_relative": fine_to_limit,
        },
        "spatial_convergence": {
            "coarse_dx_nm": cases["temporal_dt_0p0625_dx_10nm"]["dx_nm"],
            "fine_dx_nm": cases["spatial_dt_0p0625_dx_5nm"]["dx_nm"],
            "fine_minus_coarse_relative": spatial_relative,
        },
        "gates": gates,
        "physical_evidence_boundary": {
            "supports_absolute_depth_prediction": False,
            "blockers": [
                "Krueger publishes an aggregate energetic-ion flux and IEAD but not the ion-species mixture.",
                "Krueger publishes no stable C4F6 wafer flux.",
                "C2F3 and C3F4 require declared topology transfers outside Guo's printed neutral list.",
                "Most Krueger IEAD support lies above the <=370 eV Guo/Yin regression board.",
                "The source's printed off-normal angular polynomial requires the declared Table 4.2 repair.",
                "The 2.5 nm translating-layer capacity is a cross-chemistry source transfer.",
                "The amorphous-carbon mask includes unmeasured density and reduced-film parameters.",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=Path("/private/tmp"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = build_receipt(args.input_root)
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
