#!/usr/bin/env python3
"""Reclose the Zhu NPG80 daughter collision and sheath power board.

This runner does not use a TiO2 depth, SEM, blanket rate, or wafer-flux target.
For every absorbed-power sensitivity node it alternates the conserved global
reactor with the Maxwellian floating-potential sheath closure while preserving
the independently declared 276 V family self-bias anchor.  A final replay from
a stable, retained continuation path removes temporary fixed-point iterates
from the provenance chain.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import sys
from types import SimpleNamespace
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.reactor_global.argon_reactor_to_wafer import (
    maxwellian_floating_sheath_potential_eV,
)
from scripts.audit_zhu_npg80_sheath_global_coupling import (
    _effective_bohm_mass_amu,
)
from scripts.run_zhu_open_reactor import run


OUTPUT_DIR = (
    ROOT / "results" / "curated"
    / "zhu_npg80_daughter_reclosed_v1"
)
INITIAL_60W = Path(
    "results/curated/zhu_npg80_absorbed_power_ensemble_v1/power_60W.json"
)
POWER_NODES_W = (60, 90, 105, 120)
SELF_BIAS_MAGNITUDE_V = 276.0
INITIAL_GROUNDED_SHEATH_V = {
    60: 21.517845741068886,
    90: 20.474843685564036,
    105: 20.22140392984418,
    120: 20.077811323036663,
}


def _render(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _plasma_potential_V(payload: dict) -> float:
    electron_temperature_eV = (
        2.0 / 3.0 * float(payload["state"]["mean_electron_energy_eV"])
    )
    return float(maxwellian_floating_sheath_potential_eV(
        electron_temperature_eV,
        _effective_bohm_mass_amu(payload),
    ))


def _run_args(
    *,
    source_workbook: Path,
    hcl_lxcat: Path,
    f2_lxcat: Path,
    initial_state: Path,
    absorbed_power_W: float,
    grounded_sheath_V: float,
    maximum_evaluations: int,
    nonlinear_verbose: int,
) -> SimpleNamespace:
    return SimpleNamespace(
        source_workbook=source_workbook,
        hcl_lxcat=hcl_lxcat,
        f2_lxcat=f2_lxcat,
        output=None,
        initial_state_json=initial_state,
        absorbed_power_W=float(absorbed_power_W),
        minimum_reduced_field_Td=40.0,
        maximum_reduced_field_Td=900.0,
        electrode_diameter_mm=240.0,
        plasma_height_mm=30.0,
        gas_temperature_K=350.0,
        ion_temperature_eV=0.03,
        ion_mfp_um=100.0,
        mean_wall_ion_energy_eV=250.0,
        powered_electrode_sheath_drop_V=(
            SELF_BIAS_MAGNITUDE_V + float(grounded_sheath_V)
        ),
        grounded_surface_sheath_drop_V=float(grounded_sheath_V),
        neutral_reduced_diffusivity_m_inv_s=6.0e20,
        f_wall_probability=0.05,
        h_wall_probability=0.05,
        o_wall_probability=0.10,
        excited_o_wall_probability=1.0,
        c_wall_probability=1.0,
        cf3_wall_probability=0.05,
        cf2_wall_probability=0.10,
        cf_wall_probability=0.10,
        kokkoris_eedf_shape="druyvesteyn",
        chf3_f_rate_branch="voloshin_350K",
        initial_field_Td=None,
        maximum_evaluations=int(maximum_evaluations),
        residual_tolerance=2.0e-6,
        nonlinear_verbose=int(nonlinear_verbose),
    )


def solve_sheath_fixed_point(
    *,
    source_workbook: Path,
    hcl_lxcat: Path,
    f2_lxcat: Path,
    stable_initial_state: Path,
    absorbed_power_W: int,
    initial_grounded_sheath_V: float,
    maximum_fixed_point_iterations: int,
    maximum_evaluations: int,
    voltage_tolerance_V: float,
    nonlinear_verbose: int,
) -> tuple[dict, list[dict]]:
    """Return a final stable-continuation replay and iteration receipt."""

    grounded = float(initial_grounded_sheath_V)
    history: list[dict] = []
    with TemporaryDirectory(prefix=f"petch-zhu-{absorbed_power_W}W-") as tmp:
        continuation = stable_initial_state
        for iteration in range(maximum_fixed_point_iterations):
            payload = run(_run_args(
                source_workbook=source_workbook,
                hcl_lxcat=hcl_lxcat,
                f2_lxcat=f2_lxcat,
                initial_state=continuation,
                absorbed_power_W=absorbed_power_W,
                grounded_sheath_V=grounded,
                maximum_evaluations=maximum_evaluations,
                nonlinear_verbose=nonlinear_verbose,
            ))
            plasma = _plasma_potential_V(payload)
            residual = grounded - plasma
            history.append({
                "iteration": iteration,
                "grounded_sheath_drop_V": grounded,
                "powered_electrode_sheath_drop_V": (
                    grounded + SELF_BIAS_MAGNITUDE_V
                ),
                "predicted_plasma_potential_V": plasma,
                "voltage_residual_V": residual,
                "maximum_normalized_reactor_residual": payload["numerics"][
                    "maximum_normalized_residual"
                ],
                "solver_evaluations": payload["numerics"][
                    "solver_evaluations"
                ],
            })
            if abs(residual) <= voltage_tolerance_V:
                break
            grounded = plasma
            temporary = Path(tmp) / f"iteration_{iteration}.json"
            temporary.write_text(_render(payload), encoding="utf-8")
            continuation = temporary
        else:
            raise RuntimeError(
                f"{absorbed_power_W} W sheath/global fixed point did not "
                f"close in {maximum_fixed_point_iterations} iterations"
            )

    # Recompute from the retained continuation so the final state never
    # depends on a temporary file.  The nonlinear root must be independent of
    # the continuation path; the check below makes that claim executable.
    final_payload = run(_run_args(
        source_workbook=source_workbook,
        hcl_lxcat=hcl_lxcat,
        f2_lxcat=f2_lxcat,
        initial_state=stable_initial_state,
        absorbed_power_W=absorbed_power_W,
        grounded_sheath_V=grounded,
        maximum_evaluations=maximum_evaluations,
        nonlinear_verbose=nonlinear_verbose,
    ))
    final_plasma = _plasma_potential_V(final_payload)
    final_residual = grounded - final_plasma
    if abs(final_residual) > voltage_tolerance_V:
        raise RuntimeError(
            f"{absorbed_power_W} W stable replay shifted sheath fixed point "
            f"by {final_residual:.12g} V"
        )
    if final_payload["numerics"]["maximum_normalized_residual"] > 2.0e-6:
        raise RuntimeError("stable replay failed reactor conservation")
    history.append({
        "iteration": "stable_replay",
        "grounded_sheath_drop_V": grounded,
        "powered_electrode_sheath_drop_V": (
            grounded + SELF_BIAS_MAGNITUDE_V
        ),
        "predicted_plasma_potential_V": final_plasma,
        "voltage_residual_V": final_residual,
        "maximum_normalized_reactor_residual": final_payload["numerics"][
            "maximum_normalized_residual"
        ],
        "solver_evaluations": final_payload["numerics"][
            "solver_evaluations"
        ],
    })
    return final_payload, history


def build_board(args) -> dict:
    os.chdir(ROOT)
    output_dir = (
        args.output_dir
        if args.output_dir.is_absolute()
        else ROOT / args.output_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    stable_initial = INITIAL_60W
    rows = []
    for power_W in POWER_NODES_W:
        state_path = output_dir / f"power_{power_W}W.json"
        if args.resume_existing and state_path.exists():
            payload = json.loads(state_path.read_text(encoding="utf-8"))
            daughter = payload["input"]["daughter_collision_basis"]
            plasma = _plasma_potential_V(payload)
            grounded = float(payload["input"][
                "grounded_surface_sheath_drop_V"
            ])
            if (
                payload["input"]["absorbed_power_sensitivity_W"] != power_W
                or payload["numerics"]["maximum_normalized_residual"]
                > 2.0e-6
                or abs(grounded - plasma) > args.voltage_tolerance_V
                or payload["input"]["o2_source_workbook_sha256"]
                != _hash(args.source_workbook)
                or daughter["hcl_lxcat_payload_sha256"]
                != _hash(args.hcl_lxcat)
                or daughter["f2_lxcat_payload_sha256"]
                != _hash(args.f2_lxcat)
            ):
                raise RuntimeError(
                    f"existing {power_W} W state is not a valid resume node"
                )
            history = [{
                "iteration": "resumed_existing",
                "grounded_sheath_drop_V": grounded,
                "powered_electrode_sheath_drop_V": (
                    grounded + SELF_BIAS_MAGNITUDE_V
                ),
                "predicted_plasma_potential_V": plasma,
                "voltage_residual_V": grounded - plasma,
                "maximum_normalized_reactor_residual": payload["numerics"][
                    "maximum_normalized_residual"
                ],
                "solver_evaluations": payload["numerics"][
                    "solver_evaluations"
                ],
            }]
        else:
            payload, history = solve_sheath_fixed_point(
                source_workbook=args.source_workbook,
                hcl_lxcat=args.hcl_lxcat,
                f2_lxcat=args.f2_lxcat,
                stable_initial_state=stable_initial,
                absorbed_power_W=power_W,
                initial_grounded_sheath_V=INITIAL_GROUNDED_SHEATH_V[power_W],
                maximum_fixed_point_iterations=(
                    args.maximum_fixed_point_iterations
                ),
                maximum_evaluations=args.maximum_evaluations,
                voltage_tolerance_V=args.voltage_tolerance_V,
                nonlinear_verbose=args.nonlinear_verbose,
            )
            state_path.write_text(_render(payload), encoding="utf-8")
        state = payload["state"]
        grounded = float(payload["input"][
            "grounded_surface_sheath_drop_V"
        ])
        plasma = _plasma_potential_V(payload)
        rows.append({
            "absorbed_power_W": power_W,
            "state_path": str(state_path.relative_to(ROOT)),
            "state_sha256": _hash(state_path),
            "stable_continuation_path": payload["input"][
                "continuation_state"
            ]["path"],
            "stable_continuation_sha256": payload["input"][
                "continuation_state"
            ]["sha256"],
            "fixed_point_history": history,
            "grounded_sheath_drop_V": grounded,
            "powered_electrode_sheath_drop_V": (
                grounded + SELF_BIAS_MAGNITUDE_V
            ),
            "plasma_potential_V": plasma,
            "sheath_fixed_point_residual_V": grounded - plasma,
            "maximum_normalized_reactor_residual": payload["numerics"][
                "maximum_normalized_residual"
            ],
            "solver_evaluations": payload["numerics"]["solver_evaluations"],
            "represented_reduced_field_Td": state[
                "reduced_electric_field_Td"
            ],
            "total_neutral_reduced_field_Td": state[
                "implied_total_neutral_reduced_electric_field_Td"
            ],
            "mean_electron_energy_eV": state["mean_electron_energy_eV"],
            "electron_density_m3": state["electron_density_m3"],
            "electronegativity": state["electronegativity"],
            "HF_density_m3": state["densities_m3"]["HF"],
            "F2_density_m3": state["densities_m3"]["F2"],
            "HF_positive_ion_density_m3": state["densities_m3"]["HF+"],
            "total_axial_positive_ion_flux_m2_s": state[
                "total_axial_positive_ion_flux_m2_s"
            ],
            "neutral_F_thermal_flux_m2_s": state[
                "neutral_thermal_flux_m2_s"
            ]["F"],
        })
        stable_initial = state_path.relative_to(ROOT)
    return {
        "schema": "petch.zhu-npg80-daughter-reclosed-board.v1",
        "condition_id": rows[0] and (
            "zhu-2026-npg80-tio2-chf3-sf6-o2-20min"
        ),
        "target_outcome_used": False,
        "fixed_inputs": {
            "absorbed_power_sensitivity_W": list(POWER_NODES_W),
            "self_bias_family_anchor_magnitude_V": SELF_BIAS_MAGNITUDE_V,
            "self_bias_is_target_machine_measurement": False,
            "feature_or_sem_target_used": False,
            "hf_collision_basis_complete": False,
            "raw_lxcat_bytes_committed": False,
        },
        "source_inputs": {
            "o2_workbook_sha256": _hash(args.source_workbook),
            "hcl_lxcat_payload_sha256": _hash(args.hcl_lxcat),
            "f2_lxcat_payload_sha256": _hash(args.f2_lxcat),
        },
        "state_board": rows,
        "certification": {
            "all_global_particle_power_residuals_below_2e_6": all(
                row["maximum_normalized_reactor_residual"] < 2.0e-6
                for row in rows
            ),
            "all_sheath_fixed_points_below_voltage_tolerance": all(
                abs(row["sheath_fixed_point_residual_V"])
                <= args.voltage_tolerance_V
                for row in rows
            ),
            "daughter_collision_nonlinear_reclose_completed": True,
            "target_machine_absorbed_power_measured": False,
            "target_machine_self_bias_measured": False,
            "hf_vibrational_and_attachment_basis_complete": False,
            "target_tio2_surface_law_validated": False,
            "supports_unique_sem_profile": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-workbook", type=Path, required=True)
    parser.add_argument("--hcl-lxcat", type=Path, required=True)
    parser.add_argument("--f2-lxcat", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--maximum-fixed-point-iterations", type=int, default=8)
    parser.add_argument("--maximum-evaluations", type=int, default=800)
    parser.add_argument("--voltage-tolerance-V", type=float, default=0.01)
    parser.add_argument("--resume-existing", action="store_true")
    parser.add_argument(
        "--nonlinear-verbose", type=int, choices=(0, 1, 2), default=0
    )
    args = parser.parse_args()
    if (
        args.maximum_fixed_point_iterations < 1
        or args.maximum_evaluations < 1
        or not math.isfinite(args.voltage_tolerance_V)
        or args.voltage_tolerance_V <= 0.0
    ):
        raise SystemExit("invalid iteration controls")
    board = build_board(args)
    audit_path = args.output_dir / "audit.json"
    audit_path.write_text(_render(board), encoding="utf-8")
    print(f"wrote {audit_path}")


if __name__ == "__main__":
    main()
