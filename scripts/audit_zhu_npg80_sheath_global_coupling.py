#!/usr/bin/env python3
"""Audit the target-free Oxford self-bias/global-power fixed point."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from petch.reactor_global.argon_reactor_to_wafer import (
    maxwellian_floating_sheath_potential_eV,
)
from petch.reactor_global.zhu_axisymmetric_transport import (
    ZHU_POSITIVE_ION_IDENTITIES,
)
from petch.tio2_ion_dose import required_formula_units_per_incident_ion


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
BASELINE = (
    ROOT / "results" / "curated" / "zhu_npg80_open_reactor_v2"
    / "source_geometry_central.json"
)
COUPLED_DIR = (
    ROOT / "results" / "curated" / "zhu_npg80_sheath_coupled_v1"
)
COUPLED = COUPLED_DIR / "central_276V.json"
OUTPUT = COUPLED_DIR / "audit.json"
SELF_BIAS = DATA / "oxford80_self_bias_transfer.json"

FILM_NM = 700.0
ETCH_TIME_S = 1200.0
ALD_TIO2_DENSITY_KG_M3 = (3250.0, 4150.0)
SOURCE_WORKBOOK_SHA256 = (
    "6f98ac82e169d25d0a4328b1a3703f733668539adb8141d736209d199013c860"
)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _effective_bohm_mass_amu(payload: dict) -> float:
    state = payload["state"]
    densities = state["densities_m3"]
    names = state["axial_positive_ion_flux_m2_s"]
    total = sum(densities[name] for name in names)
    inverse_sqrt_mass = sum(
        densities[name]
        * ZHU_POSITIVE_ION_IDENTITIES[name].charge_number ** 0.5
        / ZHU_POSITIVE_ION_IDENTITIES[name].mass_amu ** 0.5
        for name in names
    ) / total
    return 1.0 / inverse_sqrt_mass ** 2


def _summary(payload: dict) -> dict:
    state = payload["state"]
    return {
        "mean_electron_energy_eV": state["mean_electron_energy_eV"],
        "electron_density_m3": state["electron_density_m3"],
        "electronegativity": state["electronegativity"],
        "total_axial_positive_ion_flux_m2_s": state[
            "total_axial_positive_ion_flux_m2_s"
        ],
        "neutral_F_thermal_flux_m2_s": state[
            "neutral_thermal_flux_m2_s"
        ]["F"],
        "maximum_normalized_residual": payload["numerics"][
            "maximum_normalized_residual"
        ],
        "power_density_W_m3": payload["power_density_W_m3"],
    }


def _required_yield(ion_flux_m2_s: float) -> list[float]:
    return [
        required_formula_units_per_incident_ion(
            FILM_NM, density, ion_flux_m2_s, ETCH_TIME_S
        )
        for density in ALD_TIO2_DENSITY_KG_M3
    ]


def build_receipt() -> dict:
    baseline = _load(BASELINE)
    coupled = _load(COUPLED)
    bias = _load(SELF_BIAS)
    if baseline["condition_id"] != coupled["condition_id"]:
        raise ValueError("baseline and coupled states differ in condition")
    if coupled["input"]["feature_or_sem_target_used"]:
        raise ValueError("target outcome entered the coupled solve")
    if not coupled["input"]["wall_resolved_sheath_power"]:
        raise ValueError("coupled state did not resolve powered/grounded walls")
    if coupled["input"]["o2_source_workbook_sha256"] != SOURCE_WORKBOOK_SHA256:
        raise ValueError("unexpected O2 collision workbook")
    if coupled["numerics"]["maximum_normalized_residual"] > 2.0e-6:
        raise ValueError("coupled state failed conservation")
    anchor = float(bias["mechanical_anchor_selection"]["anchor_V"])
    if bias["mechanical_anchor_selection"]["is_target_measurement"]:
        raise ValueError("family anchor was mislabeled as a target measurement")

    state = coupled["state"]
    mean_energy = float(state["mean_electron_energy_eV"])
    electron_temperature = (2.0 / 3.0) * mean_energy
    effective_mass = _effective_bohm_mass_amu(coupled)
    plasma_potential = maxwellian_floating_sheath_potential_eV(
        electron_temperature, effective_mass
    )
    powered_drop = float(coupled["input"][
        "powered_electrode_sheath_drop_V"
    ])
    grounded_drop = float(coupled["input"][
        "grounded_surface_sheath_drop_V"
    ])
    fixed_point_residual = max(
        abs(grounded_drop - plasma_potential),
        abs(powered_drop - anchor - plasma_potential),
    )
    if fixed_point_residual > 0.01:
        raise ValueError("plasma-potential/sheath fixed point did not converge")

    total_loss = sum(state["positive_ion_wall_loss_m3_s"].values())
    powered_loss = sum(
        state["powered_electrode_positive_ion_wall_loss_m3_s"].values()
    )
    grounded_loss = sum(
        state["grounded_positive_ion_wall_loss_m3_s"].values()
    )
    if abs((powered_loss + grounded_loss) / total_loss - 1.0) > 2.0e-14:
        raise ValueError("wall-resolved particle ledger does not close")

    old_flux = float(baseline["state"][
        "total_axial_positive_ion_flux_m2_s"
    ])
    new_flux = float(state["total_axial_positive_ion_flux_m2_s"])
    old_f = float(baseline["state"]["neutral_thermal_flux_m2_s"]["F"])
    new_f = float(state["neutral_thermal_flux_m2_s"]["F"])
    return {
        "schema": "petch.zhu-npg80-sheath-global-coupling-audit.v1",
        "condition_id": coupled["condition_id"],
        "target_outcome_used": False,
        "source_state": {
            "path": str(BASELINE.relative_to(ROOT)),
            "sha256": _sha(BASELINE),
            "summary": _summary(baseline),
        },
        "coupled_state": {
            "path": str(COUPLED.relative_to(ROOT)),
            "sha256": _sha(COUPLED),
            "summary": _summary(coupled),
        },
        "boundary_evidence": {
            "self_bias_transfer_path": str(SELF_BIAS.relative_to(ROOT)),
            "self_bias_transfer_sha256": _sha(SELF_BIAS),
            "family_anchor_bias_magnitude_V": anchor,
            "anchor_is_target_measurement": False,
            "o2_workbook": {
                "doi": "10.60893/figshare.jpr.30850013.v1",
                "sha256": SOURCE_WORKBOOK_SHA256,
                "license": "CC BY-NC 4.0",
                "source_file_committed": False,
            },
        },
        "fixed_point": {
            "electron_temperature_eV": electron_temperature,
            "effective_multi_ion_bohm_mass_amu": effective_mass,
            "maxwellian_floating_plasma_potential_V": plasma_potential,
            "powered_electrode_sheath_drop_V": powered_drop,
            "grounded_surface_sheath_drop_V": grounded_drop,
            "maximum_voltage_residual_V": fixed_point_residual,
            "converged_below_0p01_V": fixed_point_residual < 0.01,
        },
        "wall_ledger": {
            "powered_electrode_positive_ion_loss_fraction": (
                powered_loss / total_loss
            ),
            "grounded_positive_ion_loss_fraction": grounded_loss / total_loss,
            "particle_closure_relative_residual": (
                (powered_loss + grounded_loss) / total_loss - 1.0
            ),
            "charged_wall_power_density_W_m3": coupled["power_density_W_m3"][
                "charged_wall"
            ],
            "powered_electrode_charged_wall_power_density_W_m3": coupled[
                "power_density_W_m3"
            ]["powered_electrode_charged_wall"],
            "grounded_charged_wall_power_density_W_m3": coupled[
                "power_density_W_m3"
            ]["grounded_charged_wall"],
        },
        "effect_of_corrected_wall_power": {
            "positive_ion_flux_ratio": new_flux / old_flux,
            "neutral_F_flux_ratio": new_f / old_f,
            "electron_density_ratio": (
                state["electron_density_m3"]
                / baseline["state"]["electron_density_m3"]
            ),
            "charged_wall_power_ratio": (
                coupled["power_density_W_m3"]["charged_wall"]
                / baseline["power_density_W_m3"]["charged_wall"]
            ),
        },
        "conditional_tio2_dose": {
            "film_thickness_nm": FILM_NM,
            "etch_time_s": ETCH_TIME_S,
            "ald_tio2_density_sensitivity_kg_m3": list(
                ALD_TIO2_DENSITY_KG_M3
            ),
            "required_blanket_formula_units_per_positive_ion": (
                _required_yield(new_flux)
            ),
            "interpretation": (
                "Exact atom/dose requirement conditional on the family-anchor "
                "self-bias, 90 W absorbed-power, geometry, and wall sensitivities; "
                "not a fitted TiO2 yield or target-tool depth prediction."
            ),
        },
        "certification": {
            "global_particle_and_power_conservation_passed": True,
            "sheath_global_fixed_point_passed": True,
            "target_machine_self_bias_measured": False,
            "forward_to_absorbed_power_transfer_measured": False,
            "target_tio2_surface_response_measured": False,
            "supports_unique_absolute_depth_prediction": False,
            "result": (
                "The corrected wall physics materially strengthens the clearance "
                "case but does not identify a unique target-tool profile."
            ),
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
            raise SystemExit("committed sheath/global coupling audit is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
