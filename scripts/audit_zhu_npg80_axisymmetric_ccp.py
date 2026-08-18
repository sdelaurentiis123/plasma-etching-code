#!/usr/bin/env python3
"""Build/check the target-free Oxford NPG80 axisymmetric ion-flux receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from petch.reactor_global.axisymmetric_reaction_diffusion import (
    AxisymmetricFiniteVolumeGrid,
)
from petch.reactor_global.geometry import CylindricalReactor
from petch.reactor_global.zhu_axisymmetric_transport import (
    DeterministicZhuAxisymmetricCCPTransport,
    ZhuAxisymmetricTransportInput,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATE = (
    ROOT / "results" / "curated" / "zhu_npg80_open_reactor_v2"
    / "source_geometry_central.json"
)
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_axisymmetric_ccp_v1"
    / "audit.json"
)


def _load_input(path: Path) -> ZhuAxisymmetricTransportInput:
    payload = json.loads(path.read_text(encoding="utf-8"))
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
        source=(
            "conserved Zhu open-reactor state "
            + str(path.relative_to(ROOT))
        ),
    )


def _run_resolution(state, radial_count: int, axial_count: int):
    grid = AxisymmetricFiniteVolumeGrid.uniform(
        state.geometry,
        radial_cell_count=radial_count,
        axial_cell_count=axial_count,
    )
    result = DeterministicZhuAxisymmetricCCPTransport(
        grid=grid,
        mobility_reduced_field_Td=50.0,
    ).predict(state, optic_radius_m=1.5e-3)
    profile = result.total_lower_endcap_flux_m2_s
    return result, {
        "radial_cell_count": radial_count,
        "axial_cell_count": axial_count,
        "full_electrode_average_flux_m2_s": (
            result.full_electrode_average_flux_m2_s),
        "global_model_average_flux_m2_s": (
            result.global_model_average_flux_m2_s),
        "global_to_spatial_relative_residual": (
            result.global_to_spatial_relative_residual),
        "central_3mm_optic_average_flux_m2_s": (
            result.optic_average_flux_m2_s),
        "central_3mm_to_full_electrode_flux_ratio": (
            result.optic_to_full_electrode_flux_ratio),
        "first_to_last_annulus_flux_ratio": float(profile[0] / profile[-1]),
        "maximum_species_ledger_relative_residual": (
            result.lift_result.solution.maximum_species_ledger_relative_residual
        ),
        "maximum_inventory_relative_residual": (
            result.lift_result.maximum_inventory_relative_residual),
    }


def build_receipt(state_path: Path = DEFAULT_STATE) -> dict:
    state_path = Path(state_path).resolve()
    state = _load_input(state_path)
    resolutions = []
    central_result = None
    for radial, axial in ((24, 8), (48, 16), (96, 32)):
        result, row = _run_resolution(state, radial, axial)
        resolutions.append(row)
        if (radial, axial) == (48, 16):
            central_result = result
    assert central_result is not None
    coarse, central, fine = resolutions
    full_grid_change = abs(
        fine["full_electrode_average_flux_m2_s"]
        / central["full_electrode_average_flux_m2_s"] - 1.0)
    optic_grid_change = abs(
        fine["central_3mm_optic_average_flux_m2_s"]
        / central["central_3mm_optic_average_flux_m2_s"] - 1.0)
    profile = central_result.total_lower_endcap_flux_m2_s
    return {
        "schema": "petch.zhu-npg80-axisymmetric-ccp-audit.v1",
        "condition_id": state.condition_id,
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "input": {
            "open_reactor_state": str(state_path.relative_to(ROOT)),
            "electrode_radius_m": state.geometry.radius_m,
            "plasma_height_m_sensitivity": state.geometry.length_m,
            "frequency_hz": 13.56e6,
            "optic_radius_m": central_result.optic_radius_m,
            "electron_density_m3": state.electron_density_m3,
            "electronegativity": state.electronegativity,
            "electron_temperature_eV": state.electron_temperature_eV,
            "total_positive_ion_density_m3": (
                state.total_positive_ion_density_m3),
            "effective_multi_ion_bohm_mass_amu": (
                state.effective_bohm_mass_amu),
            "ion_temperature_eV": state.ion_temperature_eV,
            "ion_momentum_mean_free_path_m": (
                state.ion_momentum_mean_free_path_m),
        },
        "physics": {
            "source_moment": central_result.source_moment,
            "source_moment_provenance": (
                central_result.source_moment_provenance),
            "source_reference": {
                "source_id": "zhao-2019-prl-ccp-standing-waves",
                "doi": "10.1103/PhysRevLett.122.185002",
                "relevant_measurement": (
                    "At 13.56 MHz and 3 Pa the fundamental current was "
                    "radially uniform and center/edge sheath motion nearly "
                    "identical in a 21 cm electrode, 3 cm gap CCP."
                ),
                "target_gas_is_same": False,
                "target_tool_is_same": False,
            },
            "reference_mobility": {
                "source": "basurto-2002-chf3-ion-mobility",
                "species_pair": "CHF2+ in CHF3",
                "reduced_field_Td": (
                    central_result.reference_mobility.reduced_field_Td),
                "actual_mobility_m2_V_s": (
                    central_result.reference_mobility.actual_mobility_m2_V_s),
                "ambipolar_diffusion_m2_s": (
                    central_result.reference_mobility.actual_mobility_m2_V_s
                    * state.electron_temperature_eV),
                "other_ion_mobilities_measured": False,
            },
            "edge_model": (
                "Lee-Lieberman electronegative axial/radial edge factors, "
                "same closure used by the conserved 0-D state"
            ),
        },
        "resolution_board": resolutions,
        "grid_convergence": {
            "central_to_fine_full_flux_relative_change": full_grid_change,
            "central_to_fine_optic_flux_relative_change": optic_grid_change,
            "passed_0p1_percent": max(
                full_grid_change, optic_grid_change) < 1.0e-3,
        },
        "central_48x16_result": {
            "radial_center_m": central_result.radial_center_m.tolist(),
            "total_lower_endcap_flux_m2_s": profile.tolist(),
            "full_electrode_average_flux_m2_s": (
                central_result.full_electrode_average_flux_m2_s),
            "global_model_average_flux_m2_s": (
                central_result.global_model_average_flux_m2_s),
            "global_to_spatial_relative_residual": (
                central_result.global_to_spatial_relative_residual),
            "central_3mm_optic_average_flux_m2_s": (
                central_result.optic_average_flux_m2_s),
            "central_3mm_to_full_electrode_flux_ratio": (
                central_result.optic_to_full_electrode_flux_ratio),
            "first_to_last_annulus_flux_ratio": float(
                profile[0] / profile[-1]),
            "species_flux_fraction": {
                name: float(
                    central_result.species_lower_endcap_flux_m2_s[name][0]
                    / profile[0])
                for name in sorted(
                    central_result.species_lower_endcap_flux_m2_s)
            },
        },
        "certification": {
            "conserved_spatial_solve": True,
            "differentiable_fixed_topology_operator": True,
            "zero_fit_to_sem_or_depth": True,
            "global_flux_reproduced_within_1_percent": abs(
                central_result.global_to_spatial_relative_residual) < 0.01,
            "central_3mm_flux_nonuniformity_vs_full_electrode_percent": (
                100.0 * (
                    central_result.optic_to_full_electrode_flux_ratio - 1.0
                )
            ),
            "conditioned_radial_partition_supported": True,
            "absolute_target_wafer_flux_supported": False,
            "reason_absolute_is_not_yet_supported": (
                "The 30 mm active height, 50 Td CHF2+ mobility transfer, "
                "uniform target-gas source moment, and 0-D absorbed power "
                "remain physics sensitivities rather than measurements on "
                "the target serial tool."
            ),
        },
    }


def _canonical(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = _canonical(build_receipt(args.state))
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing committed receipt: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed axisymmetric CCP receipt is stale")
        print(rendered, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
