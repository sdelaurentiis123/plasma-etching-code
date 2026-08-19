#!/usr/bin/env python3
"""Lift the daughter-reclosed Zhu reactor board to wafer dose coordinates.

The receipt is prospective: it never reads the withheld SEM or a target TiO2
depth.  Each conserved 0-D power node is independently lifted through the
fixed-topology axisymmetric transport operator and converted to an atom-counted
TiO2 clearance requirement.  Surface yields and feature-floor transmission
remain explicit sensitivities rather than fitted coefficients.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.tio2_ion_dose import (
    build_clearance_gate,
    depth_nm_from_positive_ion_dose,
)
from scripts.audit_zhu_npg80_axisymmetric_ccp import (
    build_receipt as build_axisymmetric_receipt,
)

REACTOR_BOARD = (
    ROOT / "results" / "curated" / "zhu_npg80_daughter_reclosed_v1"
    / "audit.json"
)
TARGET_MANIFEST = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "recipe_manifest.json"
)
OUTPUT = (
    ROOT / "results" / "curated"
    / "zhu_npg80_daughter_wafer_dose_v1" / "audit.json"
)
ALD_TIO2_DENSITY_KG_M3 = (3250.0, 4150.0)
CANDIDATE_SURFACE_YIELDS = (0.5, 1.0, 1.5, 2.0)
CANDIDATE_FEATURE_TRANSMISSIONS = (0.25, 0.5, 0.75, 1.0)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _render(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _depth_sensitivity(
    ion_flux_m2_s: float,
    duration_s: float,
    film_thickness_nm: float,
) -> dict:
    rows = []
    for surface_yield in CANDIDATE_SURFACE_YIELDS:
        for transmission in CANDIDATE_FEATURE_TRANSMISSIONS:
            depths = [
                depth_nm_from_positive_ion_dose(
                    surface_yield,
                    density,
                    ion_flux_m2_s,
                    duration_s,
                    feature_transmission=transmission,
                    film_thickness_nm=film_thickness_nm,
                )
                for density in ALD_TIO2_DENSITY_KG_M3
            ]
            rows.append({
                "surface_yield_formula_units_per_incident_ion": (
                    surface_yield
                ),
                "run_averaged_feature_floor_transmission": transmission,
                "film_capped_depth_nm_by_density_endpoint": depths,
            })
    return {
        "surface_yield_candidates_are_fitted": False,
        "feature_transmission_candidates_are_fitted": False,
        "rows": rows,
    }


def build_receipt(
    reactor_board_path: Path = REACTOR_BOARD,
    target_manifest_path: Path = TARGET_MANIFEST,
) -> dict:
    reactor_board_path = Path(reactor_board_path).resolve()
    target_manifest_path = Path(target_manifest_path).resolve()
    reactor_board = _load(reactor_board_path)
    manifest = _load(target_manifest_path)
    if reactor_board["schema"] != (
        "petch.zhu-npg80-daughter-reclosed-board.v1"
    ):
        raise ValueError("unexpected daughter-reclosed reactor board")
    if reactor_board["target_outcome_used"]:
        raise ValueError("reactor board is not prospective")
    if (
        manifest["outcomes"]["specific_condition_sem_received"]
        or manifest["outcomes"]["post_etch_tio2_depth_nm"] is not None
    ):
        raise ValueError("wafer-dose board must remain pre-SEM")
    if reactor_board["condition_id"] != manifest["condition_id"]:
        raise ValueError("reactor and target condition IDs differ")

    duration_s = float(manifest["process"]["etch_time_s"])
    film_thickness_nm = float(
        manifest["stack"]["film_initial_thickness_nm"]
    )
    rows = []
    for reactor_row in reactor_board["state_board"]:
        state_path = ROOT / reactor_row["state_path"]
        if _hash(state_path) != reactor_row["state_sha256"]:
            raise ValueError(f"reactor state hash changed: {state_path}")
        state_payload = _load(state_path)
        if (
            state_payload["numerics"]["maximum_normalized_residual"]
            >= 2.0e-6
        ):
            raise ValueError(f"reactor state is not conserved: {state_path}")
        axisymmetric = build_axisymmetric_receipt(state_path)
        if not axisymmetric["certification"]["conserved_spatial_solve"]:
            raise ValueError(f"spatial lift is not conserved: {state_path}")
        if not axisymmetric["grid_convergence"]["passed_0p1_percent"]:
            raise ValueError(f"spatial lift is not grid converged: {state_path}")

        fine = axisymmetric["resolution_board"][-1]
        local_ion_flux = float(
            fine["central_3mm_optic_average_flux_m2_s"]
        )
        gates = [
            build_clearance_gate(
                film_thickness_nm,
                density,
                local_ion_flux,
                duration_s,
            )
            for density in ALD_TIO2_DENSITY_KG_M3
        ]
        state = state_payload["state"]
        power_density = state_payload["power_density_W_m3"]
        power_closure_residual = (
            float(power_density["absorbed"])
            - float(power_density["parent_collision"])
            - float(power_density["supplemental_collision"])
            - float(power_density["charged_wall"])
        )
        if abs(power_closure_residual) > (
            2.0e-6 * float(power_density["absorbed"])
        ):
            raise ValueError(f"reactor power is not closed: {state_path}")
        f_flux = float(state["neutral_thermal_flux_m2_s"]["F"])
        rows.append({
            "absorbed_power_sensitivity_W": reactor_row[
                "absorbed_power_W"
            ],
            "state_path": reactor_row["state_path"],
            "state_sha256": reactor_row["state_sha256"],
            "grounded_sheath_drop_V": reactor_row[
                "grounded_sheath_drop_V"
            ],
            "powered_electrode_sheath_drop_V": reactor_row[
                "powered_electrode_sheath_drop_V"
            ],
            "power_density_W_m3": power_density,
            "power_closure_residual_W_m3": power_closure_residual,
            "represented_reduced_field_Td": state[
                "reduced_electric_field_Td"
            ],
            "total_neutral_reduced_field_Td": state[
                "implied_total_neutral_reduced_electric_field_Td"
            ],
            "electron_collision_basis_neutral_fraction": state[
                "electron_collision_basis_neutral_fraction"
            ],
            "mean_electron_energy_eV": state[
                "mean_electron_energy_eV"
            ],
            "electron_density_m3": state["electron_density_m3"],
            "electronegativity": state["electronegativity"],
            "global_positive_ion_flux_m2_s": state[
                "total_axial_positive_ion_flux_m2_s"
            ],
            "central_3mm_positive_ion_flux_m2_s": local_ion_flux,
            "central_3mm_to_full_electrode_flux_ratio": fine[
                "central_3mm_to_full_electrode_flux_ratio"
            ],
            "axisymmetric_grid_convergence": axisymmetric[
                "grid_convergence"
            ],
            "axisymmetric_species_flux_fraction_at_center": (
                axisymmetric["central_48x16_result"][
                    "species_flux_fraction"
                ]
            ),
            "global_neutral_F_thermal_flux_m2_s": f_flux,
            "global_F_to_local_positive_ion_flux_ratio": (
                f_flux / local_ion_flux
            ),
            "neutral_F_radial_transfer_solved": False,
            "tio2_clearance_gate_by_density_endpoint": [
                asdict(gate) for gate in gates
            ],
            "required_blanket_formula_units_per_positive_ion": [
                min(
                    gate.required_formula_units_per_incident_ion
                    for gate in gates
                ),
                max(
                    gate.required_formula_units_per_incident_ion
                    for gate in gates
                ),
            ],
            "conditional_depth_sensitivity": _depth_sensitivity(
                local_ion_flux,
                duration_s,
                film_thickness_nm,
            ),
        })

    return {
        "schema": "petch.zhu-npg80-daughter-wafer-dose-board.v1",
        "condition_id": manifest["condition_id"],
        "target_outcome_used": False,
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "input": {
            "reactor_board_path": str(reactor_board_path.relative_to(ROOT)),
            "reactor_board_sha256": _hash(reactor_board_path),
            "target_manifest_path": str(target_manifest_path.relative_to(ROOT)),
            "target_manifest_sha256": _hash(target_manifest_path),
            "etch_duration_s": duration_s,
            "film_thickness_nm": film_thickness_nm,
            "ald_tio2_density_sensitivity_kg_m3": list(
                ALD_TIO2_DENSITY_KG_M3
            ),
            "surface_yield_candidates_formula_units_per_ion": list(
                CANDIDATE_SURFACE_YIELDS
            ),
            "feature_floor_transmission_candidates": list(
                CANDIDATE_FEATURE_TRANSMISSIONS
            ),
        },
        "power_board": rows,
        "certification": {
            "all_reactor_states_conserved": True,
            "all_axisymmetric_lifts_conserved_and_grid_converged": True,
            "target_outcome_used": False,
            "surface_coefficient_fitted_to_target": False,
            "target_machine_absorbed_power_measured": False,
            "target_machine_self_bias_measured": False,
            "species_resolved_tio2_surface_law_validated": False,
            "feature_floor_transmission_predicted": False,
            "supports_conditional_atom_counted_depths": True,
            "supports_unique_sem_profile": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactor-board", type=Path, default=REACTOR_BOARD)
    parser.add_argument("--target-manifest", type=Path, default=TARGET_MANIFEST)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = _render(build_receipt(
        args.reactor_board,
        args.target_manifest,
    ))
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing committed receipt: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed daughter wafer-dose board is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
