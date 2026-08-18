#!/usr/bin/env python3
"""Audit the source-correct NPG80 diameter without touching the blind v1 call."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from petch.tio2_ion_dose import required_formula_units_per_incident_ion


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
V1 = ROOT / "results" / "curated" / "zhu_npg80_open_reactor_v1" / "central.json"
V2 = ROOT / "results" / "curated" / "zhu_npg80_open_reactor_v2"
CENTRAL = V2 / "source_geometry_central.json"
ALTERNATE = V2 / "source_geometry_chf3_rate_alternate.json"
EVIDENCE = DATA / "machine_geometry_evidence.json"
FROZEN_GATE = DATA / "pre_sem_depth_gate.json"
OUTPUT = V2 / "audit.json"

FILM_NM = 700.0
ETCH_TIME_S = 1200.0
ALD_TIO2_DENSITY_KG_M3 = (3250.0, 4150.0)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _state_summary(payload: dict) -> dict:
    state = payload["state"]
    return {
        "electrode_diameter_mm": payload["input"]["electrode_diameter_mm"],
        "chf3_f_rate_branch": payload["input"]["chf3_f_rate_branch"],
        "maximum_normalized_residual": payload["numerics"][
            "maximum_normalized_residual"
        ],
        "mean_electron_energy_eV": state["mean_electron_energy_eV"],
        "electron_density_m3": state["electron_density_m3"],
        "electronegativity": state["electronegativity"],
        "total_axial_positive_ion_flux_m2_s": state[
            "total_axial_positive_ion_flux_m2_s"
        ],
        "neutral_F_thermal_flux_m2_s": state["neutral_thermal_flux_m2_s"]["F"],
        "electron_collision_basis_neutral_fraction": state[
            "electron_collision_basis_neutral_fraction"
        ],
    }


def _required_yield(ion_flux_m2_s: float) -> list[float]:
    return [
        required_formula_units_per_incident_ion(
            FILM_NM, density, ion_flux_m2_s, ETCH_TIME_S
        )
        for density in ALD_TIO2_DENSITY_KG_M3
    ]


def build_receipt() -> dict:
    evidence = _load(EVIDENCE)
    old = _load(V1)
    central = _load(CENTRAL)
    alternate = _load(ALTERNATE)
    frozen = _load(FROZEN_GATE)
    if evidence["manufacturer_family_specification"]["electrode_diameter_mm"] != 240.0:
        raise ValueError("source-correct electrode diameter changed")
    if old["input"]["electrode_diameter_mm"] != 170.0:
        raise ValueError("immutable v1 state no longer carries its declared geometry")
    for payload in (central, alternate):
        if payload["input"]["electrode_diameter_mm"] != 240.0:
            raise ValueError("v2 state does not use the source-correct diameter")
        if payload["input"]["feature_or_sem_target_used"]:
            raise ValueError("target data entered a target-free reactor solve")
        if payload["numerics"]["maximum_normalized_residual"] > 2.0e-6:
            raise ValueError("v2 state failed its conservation residual gate")
    if central["input"]["chf3_f_rate_branch"] != "voloshin_350K":
        raise ValueError("unexpected central chemistry branch")
    if alternate["input"]["chf3_f_rate_branch"] != "lim_700K":
        raise ValueError("unexpected alternate chemistry branch")
    if not frozen["frozen_before_specific_condition_sem"]:
        raise ValueError("the preregistered depth gate is no longer frozen")

    old_flux = old["state"]["total_axial_positive_ion_flux_m2_s"]
    central_flux = central["state"]["total_axial_positive_ion_flux_m2_s"]
    alternate_flux = alternate["state"]["total_axial_positive_ion_flux_m2_s"]
    old_f = old["state"]["neutral_thermal_flux_m2_s"]["F"]
    central_f = central["state"]["neutral_thermal_flux_m2_s"]["F"]
    alternate_f = alternate["state"]["neutral_thermal_flux_m2_s"]["F"]
    return {
        "schema": "petch.zhu-npg80-geometry-correction-audit.v1",
        "condition_id": central["condition_id"],
        "source_evidence": {
            "path": str(EVIDENCE.relative_to(ROOT)),
            "sha256": _sha(EVIDENCE),
            "exact_tool": evidence["exact_tool"],
            "manufacturer_family_specification": evidence[
                "manufacturer_family_specification"
            ],
        },
        "immutable_preregistered_v1": {
            "path": str(V1.relative_to(ROOT)),
            "sha256": _sha(V1),
            "blind_gate_path": str(FROZEN_GATE.relative_to(ROOT)),
            "blind_gate_sha256": _sha(FROZEN_GATE),
            "forecast_changed_retroactively": False,
            "state": _state_summary(old),
        },
        "source_geometry_v2": {
            "central_path": str(CENTRAL.relative_to(ROOT)),
            "central_sha256": _sha(CENTRAL),
            "alternate_path": str(ALTERNATE.relative_to(ROOT)),
            "alternate_sha256": _sha(ALTERNATE),
            "central": _state_summary(central),
            "published_rate_alternate": _state_summary(alternate),
        },
        "geometry_effect_central_over_v1": {
            "positive_ion_flux_ratio": central_flux / old_flux,
            "neutral_F_flux_ratio": central_f / old_f,
            "electron_density_ratio": (
                central["state"]["electron_density_m3"]
                / old["state"]["electron_density_m3"]
            ),
        },
        "published_chf3_f_rate_branch_effect_at_240mm": {
            "positive_ion_flux_ratio": alternate_flux / central_flux,
            "neutral_F_flux_ratio": alternate_f / central_f,
        },
        "revised_conditional_clearance_ledger": {
            "film_thickness_nm": FILM_NM,
            "etch_time_s": ETCH_TIME_S,
            "ald_tio2_density_sensitivity_kg_m3": list(ALD_TIO2_DENSITY_KG_M3),
            "central_required_formula_units_per_positive_ion": _required_yield(
                central_flux
            ),
            "alternate_required_formula_units_per_positive_ion": _required_yield(
                alternate_flux
            ),
            "interpretation": (
                "This is an exact atom/dose requirement conditional on each global "
                "state, not a fitted TiO2 yield or a certified local wafer flux."
            ),
        },
        "verdict": {
            "frozen_binary_forecast": frozen["blind_forecast"]["primary_outcome"],
            "current_reactor_only_call": "clearance remains unresolved",
            "reason": (
                "The source-correct diameter raises the central required blanket "
                "yield from 0.62--0.80 to about 1.14--1.46 formula units per ion. "
                "That is physically reachable but no target-tool TiO2 yield, local "
                "ion flux, feature transmission, or exact self-bias has been measured."
            ),
            "spatial_profile_authorized": False,
        },
    }


def _render(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = _render(build_receipt())
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed NPG80 geometry audit is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
