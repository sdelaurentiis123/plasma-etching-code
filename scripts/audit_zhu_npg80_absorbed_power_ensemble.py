#!/usr/bin/env python3
"""Audit the Oxford absorbed-power/global-state sensitivity ensemble."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.reactor_global.argon_reactor_to_wafer import (
    maxwellian_floating_sheath_potential_eV,
)
from petch.tio2_ion_dose import required_formula_units_per_incident_ion
from scripts.audit_zhu_npg80_self_bias_global_ensemble import (
    _radial_transfer,
    _transport_state,
)
from scripts.audit_zhu_npg80_sheath_global_coupling import (
    _effective_bohm_mass_amu,
)


ENSEMBLE_DIR = (
    ROOT / "results" / "curated"
    / "zhu_npg80_absorbed_power_ensemble_v1"
)
CENTRAL_90W = (
    ROOT / "results" / "curated" / "zhu_npg80_sheath_coupled_v1"
    / "central_276V.json"
)
STATE_PATHS = {
    60: ENSEMBLE_DIR / "power_60W.json",
    90: CENTRAL_90W,
    105: ENSEMBLE_DIR / "power_105W.json",
    120: ENSEMBLE_DIR / "power_120W.json",
}
REJECTED_120W = ENSEMBLE_DIR / "power_120W_default_domain_failure.json"
HONG_CSV = (
    ROOT / "data" / "experimental" / "hong_2023_tio2"
    / "figure2_feature_response.csv"
)
OUTPUT = ENSEMBLE_DIR / "audit.json"

FILM_NM = 700.0
ETCH_TIME_S = 1200.0
FORWARD_POWER_W = 150.0
SELF_BIAS_MAGNITUDE_V = 276.0
ALD_TIO2_DENSITY_KG_M3 = (3250.0, 4150.0)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _continuation_check(payload: dict) -> dict:
    continuation = payload["input"]["continuation_state"]
    path = Path(continuation["path"])
    if not path.is_absolute():
        path = ROOT / path
    if not path.exists() or _sha(path) != continuation["sha256"]:
        raise ValueError("power-state continuation checksum mismatch")
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": continuation["sha256"],
    }


def _required_yield(flux_m2_s: float) -> list[float]:
    return [
        required_formula_units_per_incident_ion(
            FILM_NM, density, flux_m2_s, ETCH_TIME_S
        )
        for density in ALD_TIO2_DENSITY_KG_M3
    ]


def _state_row(power_W: int, path: Path) -> dict:
    payload = _load(path)
    inputs = payload["input"]
    state = payload["state"]
    if inputs["feature_or_sem_target_used"]:
        raise ValueError("target outcome entered a power state")
    if not inputs["wall_resolved_sheath_power"]:
        raise ValueError("power state is not wall resolved")
    if inputs["absorbed_power_sensitivity_W"] != power_W:
        raise ValueError("power-state label mismatch")
    if payload["numerics"]["maximum_normalized_residual"] > 2.0e-6:
        raise ValueError("power state failed conservation")
    powered = float(inputs["powered_electrode_sheath_drop_V"])
    grounded = float(inputs["grounded_surface_sheath_drop_V"])
    if abs(powered - grounded - SELF_BIAS_MAGNITUDE_V) > 1.0e-9:
        raise ValueError("power state does not preserve self-bias")

    electron_temperature = 2.0 * state["mean_electron_energy_eV"] / 3.0
    plasma_potential = maxwellian_floating_sheath_potential_eV(
        electron_temperature, _effective_bohm_mass_amu(payload)
    )
    voltage_residual = max(
        abs(grounded - plasma_potential),
        abs(powered - SELF_BIAS_MAGNITUDE_V - plasma_potential),
    )
    if voltage_residual > 0.01:
        raise ValueError("power-state sheath/global fixed point failed")

    radial = _radial_transfer(_transport_state(path.resolve()))
    optic_flux = float(radial.optic_average_flux_m2_s)
    collision_fraction = float(
        state["electron_collision_basis_neutral_fraction"]
    )
    represented_field = float(state["reduced_electric_field_Td"])
    implied_total_field = float(
        state["implied_total_neutral_reduced_electric_field_Td"]
    )
    if abs(represented_field * collision_fraction / implied_total_field - 1) > 1e-12:
        raise ValueError("reduced-field convention ledger failed")
    return {
        "absorbed_power_W": power_W,
        "absorbed_to_forward_power_fraction": power_W / FORWARD_POWER_W,
        "state_path": str(path.relative_to(ROOT)),
        "state_sha256": _sha(path),
        "continuation": _continuation_check(payload),
        "reduced_field_bounds_Td": inputs.get(
            "reduced_field_bounds_Td", [205.0, 600.0]
        ),
        "represented_neutral_reduced_field_Td": represented_field,
        "electron_collision_basis_neutral_fraction": collision_fraction,
        "implied_total_neutral_reduced_field_Td": implied_total_field,
        "mean_electron_energy_eV": state["mean_electron_energy_eV"],
        "electron_density_m3": state["electron_density_m3"],
        "electronegativity": state["electronegativity"],
        "plasma_potential_V": plasma_potential,
        "sheath_fixed_point_residual_V": voltage_residual,
        "full_electrode_global_ion_flux_m2_s": state[
            "total_axial_positive_ion_flux_m2_s"
        ],
        "central_3mm_optic_ion_flux_m2_s": optic_flux,
        "central_3mm_to_full_electrode_flux_ratio": (
            radial.optic_to_full_electrode_flux_ratio
        ),
        "axisymmetric_global_flux_relative_residual": (
            radial.global_to_spatial_relative_residual
        ),
        "neutral_F_thermal_flux_m2_s": state[
            "neutral_thermal_flux_m2_s"
        ]["F"],
        "required_central_3mm_formula_units_per_positive_ion": (
            _required_yield(optic_flux)
        ),
        "maximum_normalized_reactor_residual": payload["numerics"][
            "maximum_normalized_residual"
        ],
    }


def _hong_arde() -> dict[str, float]:
    with HONG_CSV.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        row["plasma_mode"]: float(row["value"])
        for row in rows
        if row["quantity"] == "arde_p1_depth_over_p3_depth"
        and row["chemistry"] == "C4F8/SF6/Ar"
    }


def build_receipt() -> dict:
    rows = [_state_row(power, path) for power, path in STATE_PATHS.items()]
    rows.sort(key=lambda row: row["absorbed_power_W"])
    ion_flux = [row["central_3mm_optic_ion_flux_m2_s"] for row in rows]
    f_flux = [row["neutral_F_thermal_flux_m2_s"] for row in rows]
    if not all(left < right for left, right in zip(ion_flux, ion_flux[1:])):
        raise ValueError("ion flux is not monotone with absorbed power")
    if not all(left < right for left, right in zip(f_flux, f_flux[1:])):
        raise ValueError("F flux is not monotone with absorbed power")
    rejected = _load(REJECTED_120W)
    if rejected["accepted"] or rejected["maximum_normalized_residual"] <= 2e-6:
        raise ValueError("rejected field-domain node is mislabeled")
    hong = _hong_arde()
    all_yields = [
        value
        for row in rows
        for value in row[
            "required_central_3mm_formula_units_per_positive_ion"
        ]
    ]
    return {
        "schema": "petch.zhu-npg80-absorbed-power-ensemble.v1",
        "condition_id": "zhu-2026-npg80-tio2-chf3-sf6-o2-20min",
        "sem_or_depth_target_used": False,
        "fixed_inputs": {
            "forward_power_setpoint_W": FORWARD_POWER_W,
            "self_bias_magnitude_V": SELF_BIAS_MAGNITUDE_V,
            "pressure_mTorr": 30.0,
            "etch_time_s": ETCH_TIME_S,
            "film_thickness_nm": FILM_NM,
            "ald_tio2_density_sensitivity_kg_m3": list(
                ALD_TIO2_DENSITY_KG_M3
            ),
        },
        "state_board": rows,
        "power_response": {
            "ion_flux_monotonically_increases_with_absorbed_power": True,
            "F_flux_monotonically_increases_with_absorbed_power": True,
            "central_ion_flux_120W_to_60W_ratio": ion_flux[-1] / ion_flux[0],
            "F_flux_120W_to_60W_ratio": f_flux[-1] / f_flux[0],
            "required_surface_yield_envelope": [
                min(all_yields), max(all_yields)
            ],
            "interpretation": (
                "The recipe's 150 W forward setpoint does not choose one row. "
                "These are conserved conditional states spanning 40--80% "
                "forward-to-absorbed transfer at fixed self-bias."
            ),
        },
        "field_domain_audit": {
            "rejected_default_domain_path": str(
                REJECTED_120W.relative_to(ROOT)
            ),
            "rejected_default_domain_sha256": _sha(REJECTED_120W),
            "default_upper_represented_field_Td": 600.0,
            "extended_upper_represented_field_Td": 900.0,
            "accepted_120W_represented_field_Td": rows[-1][
                "represented_neutral_reduced_field_Td"
            ],
            "accepted_120W_implied_total_neutral_field_Td": rows[-1][
                "implied_total_neutral_reduced_field_Td"
            ],
            "accepted_120W_collision_basis_fraction": rows[-1][
                "electron_collision_basis_neutral_fraction"
            ],
            "interpretation": (
                "The extended coordinate is required because the electron "
                "collision deck represents only a shrinking parent-neutral "
                "fraction. The total-neutral E/N remains below 100 Td. This "
                "exposes missing daughter-gas electron collisions and makes "
                "the 120 W row a domain sensitivity, not a stronger grade."
            ),
        },
        "independent_feature_response": {
            "source": "hong-2023-tio2",
            "same_tool_or_chemistry": False,
            "c4f8_sf6_ar_p1_over_p3_depth_ratio": hong,
            "cw_layout_depth_ratio": hong["CW"],
            "warning": (
                "Constrains the existence and scale of TiO2 pattern response; "
                "never transferred as a Zhu absolute coefficient."
            ),
        },
        "corrected_depth_verdict": {
            "original_frozen_binary_call_preserved": True,
            "current_unique_depth_supported": False,
            "current_clearance_supported_over_full_power_envelope": False,
            "reason": (
                "The atom-count requirement spans the independent plausible "
                "surface-response range as absorbed power changes; target "
                "TiO2/Cr surface kinetics and feature transport remain "
                "unmeasured."
            ),
        },
        "certification": {
            "all_accepted_states_conserve_below_2e_6": all(
                row["maximum_normalized_reactor_residual"] < 2e-6
                for row in rows
            ),
            "all_sheath_fixed_points_below_0p01_V": all(
                row["sheath_fixed_point_residual_V"] < 0.01
                for row in rows
            ),
            "all_axisymmetric_global_residuals_below_1_percent": all(
                abs(row["axisymmetric_global_flux_relative_residual"]) < 0.01
                for row in rows
            ),
            "absorbed_power_measured": False,
            "daughter_gas_electron_collision_basis_complete": False,
            "target_tio2_cr_surface_law_validated": False,
            "supports_unique_profile_depth": False,
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
            raise SystemExit("committed absorbed-power ensemble is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
