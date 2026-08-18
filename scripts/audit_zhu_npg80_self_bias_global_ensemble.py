#!/usr/bin/env python3
"""Audit the Oxford self-bias/global-state sensitivity ensemble."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from petch.reactor_global.argon_reactor_to_wafer import (
    maxwellian_floating_sheath_potential_eV,
)
from petch.reactor_global.axisymmetric_reaction_diffusion import (
    AxisymmetricFiniteVolumeGrid,
)
from petch.reactor_global.geometry import CylindricalReactor
from petch.reactor_global.zhu_axisymmetric_transport import (
    DeterministicZhuAxisymmetricCCPTransport,
    ZhuAxisymmetricTransportInput,
)
from petch.tio2_ion_dose import required_formula_units_per_incident_ion


ROOT = Path(__file__).resolve().parents[1]
ENSEMBLE_DIR = (
    ROOT / "results" / "curated"
    / "zhu_npg80_self_bias_global_ensemble_v1"
)
CENTRAL = (
    ROOT / "results" / "curated" / "zhu_npg80_sheath_coupled_v1"
    / "central_276V.json"
)
SELF_BIAS = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "oxford80_self_bias_transfer.json"
)
OUTPUT = ENSEMBLE_DIR / "audit.json"
BIAS_PATHS = {
    200: ENSEMBLE_DIR / "bias_200V.json",
    250: ENSEMBLE_DIR / "bias_250V.json",
    276: CENTRAL,
    300: ENSEMBLE_DIR / "bias_300V.json",
    360: ENSEMBLE_DIR / "bias_360V.json",
    387: ENSEMBLE_DIR / "bias_387V.json",
    400: ENSEMBLE_DIR / "bias_400V.json",
}

FILM_NM = 700.0
ETCH_TIME_S = 1200.0
ALD_TIO2_DENSITY_KG_M3 = (3250.0, 4150.0)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _required_yield(flux_m2_s: float) -> list[float]:
    return [
        required_formula_units_per_incident_ion(
            FILM_NM, density, flux_m2_s, ETCH_TIME_S
        )
        for density in ALD_TIO2_DENSITY_KG_M3
    ]


def _transport_state(path: Path) -> ZhuAxisymmetricTransportInput:
    payload = _load(path)
    state = payload["state"]
    model_input = payload["input"]
    densities = state["densities_m3"]
    positive_flux = state["axial_positive_ion_flux_m2_s"]
    positive_density = {name: densities[name] for name in positive_flux}
    total_neutral = sum(
        value for name, value in densities.items()
        if name != "e" and not name.endswith("+") and not name.endswith("-")
    )
    return ZhuAxisymmetricTransportInput(
        condition_id=payload["condition_id"],
        geometry=CylindricalReactor(
            radius_m=0.5e-3 * model_input["electrode_diameter_mm"],
            length_m=1.0e-3 * model_input["plasma_height_mm"],
        ),
        positive_ion_density_m3=positive_density,
        global_axial_positive_ion_flux_m2_s=positive_flux,
        electron_density_m3=state["electron_density_m3"],
        electronegativity=state["electronegativity"],
        mean_electron_energy_eV=state["mean_electron_energy_eV"],
        total_neutral_density_m3=total_neutral,
        ion_temperature_eV=0.03,
        ion_momentum_mean_free_path_m=1.0e-6 * model_input["ion_mfp_um"],
        source="committed self-bias/global ensemble state",
    )


def _radial_transfer(state: ZhuAxisymmetricTransportInput):
    grid = AxisymmetricFiniteVolumeGrid.uniform(
        state.geometry, radial_cell_count=48, axial_cell_count=16
    )
    return DeterministicZhuAxisymmetricCCPTransport(
        grid=grid, mobility_reduced_field_Td=50.0
    ).predict(state, optic_radius_m=1.5e-3)


def _state_row(bias_V: int, path: Path) -> dict:
    payload = _load(path)
    state = payload["state"]
    inputs = payload["input"]
    if inputs["feature_or_sem_target_used"]:
        raise ValueError("target outcome entered a voltage state")
    if not inputs["wall_resolved_sheath_power"]:
        raise ValueError("voltage state is not wall resolved")
    if payload["numerics"]["maximum_normalized_residual"] > 2.0e-6:
        raise ValueError("voltage state failed conservation")
    powered = float(inputs["powered_electrode_sheath_drop_V"])
    grounded = float(inputs["grounded_surface_sheath_drop_V"])
    if abs(powered - grounded - bias_V) > 1.0e-9:
        raise ValueError("powered/grounded drops do not recover self-bias")
    continuation = inputs["continuation_state"]
    continuation_path = ROOT / continuation["path"]
    if _sha(continuation_path) != continuation["sha256"]:
        raise ValueError("continuation-state checksum mismatch")

    transport_state = _transport_state(path.resolve())
    plasma_potential = maxwellian_floating_sheath_potential_eV(
        transport_state.electron_temperature_eV,
        transport_state.effective_bohm_mass_amu,
    )
    fixed_point_residual = max(
        abs(grounded - plasma_potential),
        abs(powered - bias_V - plasma_potential),
    )
    if fixed_point_residual > 0.01:
        raise ValueError("voltage-state plasma-potential fixed point failed")
    radial = _radial_transfer(transport_state)
    optic_flux = radial.optic_average_flux_m2_s
    return {
        "self_bias_magnitude_V": bias_V,
        "state_path": str(path.relative_to(ROOT)),
        "state_sha256": _sha(path),
        "continuation_path": continuation["path"],
        "continuation_sha256": continuation["sha256"],
        "plasma_potential_V": plasma_potential,
        "fixed_point_residual_V": fixed_point_residual,
        "mean_electron_energy_eV": state["mean_electron_energy_eV"],
        "electron_density_m3": state["electron_density_m3"],
        "electronegativity": state["electronegativity"],
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


def build_receipt() -> dict:
    transfer = _load(SELF_BIAS)
    rows = [_state_row(bias, path) for bias, path in BIAS_PATHS.items()]
    rows.sort(key=lambda row: row["self_bias_magnitude_V"])
    by_bias = {row["self_bias_magnitude_V"]: row for row in rows}
    ion_flux = [row["full_electrode_global_ion_flux_m2_s"] for row in rows]
    fluorine_flux = [row["neutral_F_thermal_flux_m2_s"] for row in rows]
    if not all(left > right for left, right in zip(ion_flux, ion_flux[1:])):
        raise ValueError("ion-flux voltage trend is not monotone")
    if not all(
        left > right for left, right in zip(fluorine_flux, fluorine_flux[1:])
    ):
        raise ValueError("F-flux voltage trend is not monotone")

    # The exact-tool source reports only censored 300 V -> 200 V thresholds.
    # Simpson integration of the three solved nodes is exact for a quadratic
    # response along the explicitly linear threshold witness; it is not a
    # claim that the censored endpoints are the true target waveform.
    threshold_ion_average = (
        by_bias[300]["full_electrode_global_ion_flux_m2_s"]
        + 4.0 * by_bias[250]["full_electrode_global_ion_flux_m2_s"]
        + by_bias[200]["full_electrode_global_ion_flux_m2_s"]
    ) / 6.0
    threshold_f_average = (
        by_bias[300]["neutral_F_thermal_flux_m2_s"]
        + 4.0 * by_bias[250]["neutral_F_thermal_flux_m2_s"]
        + by_bias[200]["neutral_F_thermal_flux_m2_s"]
    ) / 6.0
    histories = {item["name"]: item for item in transfer[
        "sensitivity_histories"
    ]}
    if "exact-NGP80 conditioning thresholds" not in histories:
        raise ValueError("missing exact-tool threshold witness")
    return {
        "schema": "petch.zhu-npg80-self-bias-global-ensemble.v1",
        "condition_id": transfer["condition_id"],
        "sem_or_depth_target_used": False,
        "self_bias_transfer": {
            "path": str(SELF_BIAS.relative_to(ROOT)),
            "sha256": _sha(SELF_BIAS),
            "supports_unique_target_bias": transfer["certification"][
                "supports_unique_target_bias"
            ],
        },
        "state_board": rows,
        "voltage_response": {
            "ion_flux_monotonically_decreases_with_bias_at_fixed_absorbed_power": (
                True
            ),
            "F_flux_monotonically_decreases_with_bias_at_fixed_absorbed_power": (
                True
            ),
            "ion_flux_200V_to_400V_ratio": ion_flux[0] / ion_flux[-1],
            "F_flux_200V_to_400V_ratio": fluorine_flux[0] / fluorine_flux[-1],
            "mechanism": (
                "At fixed absorbed power, a larger powered-electrode sheath "
                "drop consumes more charged-particle wall power and leaves less "
                "power for ionization/dissociation sustaining the bulk plasma."
            ),
        },
        "exact_ngp80_conditioning_threshold_history": {
            "endpoint_bias_magnitude_V": [300.0, 200.0],
            "endpoints_are_censor_thresholds": True,
            "linear_history_is_sensitivity_not_measurement": True,
            "simpson_node_bias_magnitude_V": [300.0, 250.0, 200.0],
            "time_averaged_full_electrode_ion_flux_m2_s": (
                threshold_ion_average
            ),
            "time_averaged_neutral_F_thermal_flux_m2_s": threshold_f_average,
            "required_blanket_formula_units_per_positive_ion": (
                _required_yield(threshold_ion_average)
            ),
        },
        "certification": {
            "all_reactor_conservation_gates_passed": True,
            "all_sheath_global_fixed_points_below_0p01_V": True,
            "all_axisymmetric_global_residuals_below_1_percent": all(
                abs(row["axisymmetric_global_flux_relative_residual"]) < 0.01
                for row in rows
            ),
            "target_machine_self_bias_measured": False,
            "absorbed_power_measured": False,
            "tio2_cr_surface_law_validated_for_target": False,
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
            raise SystemExit("committed self-bias/global ensemble is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
