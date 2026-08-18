#!/usr/bin/env python3
"""Build/check the blind Zhu TiO2 reactor-dose clearance gate."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

from petch.tio2_ion_dose import (
    build_clearance_gate,
    minimum_feature_transmission_for_depth,
    required_formula_units_per_incident_ion,
)


ROOT = Path(__file__).resolve().parents[1]
TARGET_DIR = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
MANIFEST = TARGET_DIR / "recipe_manifest.json"
ANALOG_BOARD = TARGET_DIR / "janissen_tio2_analog_board.json"
REACTOR_STATE = (
    ROOT / "results" / "curated" / "zhu_npg80_open_reactor_v1"
    / "central.json"
)
OUTPUT = TARGET_DIR / "pre_sem_depth_gate.json"

# Published ALD TiO2 density measurements span this bracket.  The low endpoint
# is a direct Piercy et al. XRR measurement.  The high endpoint is a direct
# result from the Go et al. process as reported with explicit scope in Saari
# et al.'s primary density/crystallization paper.  It is deliberately retained
# as a material sensitivity because Zhu's deposition recipe/phase is not yet
# reported.
ALD_TIO2_DENSITY_KG_M3 = (3250.0, 4150.0)
DENSITY_SOURCES = (
    {
        "bibkey": "piercy-2017-ald-tio2-density",
        "doi": "10.1116/1.4979047",
        "measurement": "XRR and ellipsometry",
        "reported_density_kg_m3": [3250.0, 3680.0],
        "scope": "TiCl4/H2O ALD, 38--125 C",
    },
    {
        "bibkey": "saari-2022-ald-tio2-density",
        "doi": "10.1021/acs.jpcc.2c04905",
        "measurement": (
            "XRR in the source paper plus its explicitly scoped report of "
            "the Go et al. oxygen-plasma-duration split"
        ),
        "reported_density_kg_m3": [3730.0, 4150.0],
        "scope": "TDMAT/O3 ALD; amorphous to anatase process split",
    },
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_receipt() -> dict:
    manifest = _load(MANIFEST)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_STATE)
    if (
        manifest["measurement_state"] != "pre_sem_specific_condition"
        or manifest["outcomes"]["specific_condition_sem_received"]
        or manifest["outcomes"]["post_etch_tio2_depth_nm"] is not None
    ):
        raise ValueError("depth gate must be frozen before the target SEM")
    condition_id = manifest["condition_id"]
    if analog["condition_id"] != condition_id:
        raise ValueError("analog board belongs to another condition")
    if reactor["condition_id"] != condition_id:
        raise ValueError("reactor state belongs to another condition")

    film_nm = float(manifest["stack"]["film_initial_thickness_nm"])
    duration_s = float(manifest["process"]["etch_time_s"])
    mask_nm = float(manifest["stack"]["mask_initial_thickness_nm"])
    ion_flux = float(
        reactor["state"]["total_axial_positive_ion_flux_m2_s"])
    f_flux = float(reactor["state"]["neutral_thermal_flux_m2_s"]["F"])
    density_gates = [
        build_clearance_gate(film_nm, density, ion_flux, duration_s)
        for density in ALD_TIO2_DENSITY_KG_M3
    ]
    required_yields = [
        gate.required_formula_units_per_incident_ion
        for gate in density_gates
    ]
    transmission_board = {
        str(candidate_yield): [
            minimum_feature_transmission_for_depth(
                film_nm,
                density,
                ion_flux,
                duration_s,
                candidate_yield,
            )
            for density in ALD_TIO2_DENSITY_KG_M3
        ]
        for candidate_yield in (0.5, 1.0, 1.5, 2.0)
    }

    source_feature = analog["source_feature_depth_board"]
    source_rate_range = (
        float(source_feature["minimum_implied_rate_nm_min"]),
        float(source_feature["maximum_implied_rate_nm_min"]),
    )
    source_depth_comparison = [
        min(film_nm, rate * duration_s / 60.0)
        for rate in source_rate_range
    ]
    source_effective_yield_range = [
        required_formula_units_per_incident_ion(
            rate * duration_s / 60.0,
            density,
            ion_flux,
            duration_s,
        )
        for rate, density in (
            (source_rate_range[0], ALD_TIO2_DENSITY_KG_M3[0]),
            (source_rate_range[1], ALD_TIO2_DENSITY_KG_M3[1]),
        )
    ]

    sweep = analog["source_power_sweep_interpolation"]
    closest = analog["closest_stack_witness"]
    mask_supported_closest_nm = (
        mask_nm * float(closest["reported_approximate_selectivity"])
    )
    return {
        "schema": "petch.zhu-tio2-pre-sem-depth-gate.v1",
        "condition_id": condition_id,
        "frozen_before_specific_condition_sem": True,
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "coefficient_selected_from_target": None,
        "inputs": {
            "manifest": {
                "path": str(MANIFEST.relative_to(ROOT)),
                "sha256": _sha256(MANIFEST),
            },
            "reactor_sensitivity_state": {
                "path": str(REACTOR_STATE.relative_to(ROOT)),
                "sha256": _sha256(REACTOR_STATE),
                "total_axial_positive_ion_flux_m2_s": ion_flux,
                "neutral_F_thermal_flux_m2_s": f_flux,
                "F_to_positive_ion_incident_flux_ratio": f_flux / ion_flux,
                "certified_as_local_wafer_flux": False,
            },
            "ald_tio2_density_sensitivity_kg_m3": list(
                ALD_TIO2_DENSITY_KG_M3),
            "density_sources": list(DENSITY_SOURCES),
            "film_thickness_nm": film_nm,
            "etch_duration_s": duration_s,
            "cr_mask_thickness_nm": mask_nm,
        },
        "central_reactor_dose_clearance_gate": {
            "density_endpoint_ledgers": [
                asdict(gate) for gate in density_gates
            ],
            "required_blanket_formula_units_per_positive_ion": [
                min(required_yields), max(required_yields)
            ],
            "minimum_run_averaged_feature_transmission_by_candidate_surface_yield": (
                transmission_board
            ),
            "interpretation": (
                "These are exact atom/dose requirements conditional on the "
                "central global axial flux. They are not fitted yields and the "
                "global axial flux is not certified as a local wafer flux."
            ),
        },
        "independent_tio2_process_comparison": {
            "source": "janissen-2016-tio2-rie",
            "same_recipe_or_machine": False,
            "feature_rate_range_nm_min": list(source_rate_range),
            "target_required_rate_nm_min": film_nm / (duration_s / 60.0),
            "twenty_minute_film_capped_depth_comparison_nm": (
                source_depth_comparison
            ),
            "effective_blanket_yield_range_if_central_flux_were_shared": (
                source_effective_yield_range
            ),
            "warning": (
                "The depth pair is a cross-machine comparison, not a "
                "probability or target uncertainty interval."
            ),
        },
        "independent_sf6_direction_evidence": {
            "source": "hegeman-2020-tio2-rie",
            "common_source_condition": {
                "pressure_mTorr": 10.0,
                "bias_power_W": 20.0,
                "icp_power_W": 1500.0,
                "temperature_C": 10.0,
                "single_gas_flow_sccm": 25.0,
            },
            "reported_tio2_rate_nm_min": {
                "SF6": 55.0,
                "CHF3": 15.0,
            },
            "supports_direction_not_magnitude": True,
            "warning": (
                "The source is an ICP tool with sputtered TiO2 and a "
                "single-gas comparison; neither rate is a Zhu coefficient."
            ),
        },
        "cr_mask_survival_straddle": {
            "closest_feature_witness_selectivity": float(
                closest["reported_approximate_selectivity"]),
            "closest_feature_witness_supported_tio2_nm": (
                mask_supported_closest_nm
            ),
            "source_power_sweep_interpolated_selectivity": float(
                sweep["source_system_selectivity"]),
            "source_power_sweep_supported_tio2_nm": (
                mask_nm * float(sweep["source_system_selectivity"])
            ),
            "full_clear_selectivity_requirement": film_nm / mask_nm,
            "interpretation": (
                "Independent source rows straddle mask survival; mask loss "
                "cannot be collapsed into a TiO2 yield."
            ),
        },
        "blind_forecast": {
            "primary_outcome": "TiO2 film clearance expected",
            "predicted_tio2_depth_nm": film_nm,
            "forecast_type": "film-capped binary clearance call",
            "basis": [
                "the central reactor dose requires 0.62--0.80 TiO2 formula units per incident positive ion before feature attenuation",
                "the audited adjacent TiO2 feature-rate board straddles but mostly exceeds the 35 nm/min clearance threshold",
                "an independent common-condition TiO2 experiment measures faster SF6 than CHF3 etching, fixing only the sign of the added-SF6 evidence",
            ],
            "highest_risk_failure_mode": (
                "Cr mask exhaustion or feature-bottom transport can prevent "
                "a usable 700 nm pillar even if open TiO2 clears."
            ),
            "not_claimed": [
                "atomic-accuracy profile",
                "unique continuous depth uncertainty interval",
                "validated local wafer ion flux",
                "validated TiO2 or Cr surface yield for this condition",
            ],
        },
        "decisive_post_forecast_discriminators": [
            "scale-bearing cross-section SEM for the frozen condition",
            "residual Cr thickness or direct evidence of mask exhaustion",
            "blanket TiO2 loss on an unpatterned witness",
            "achieved powered-electrode DC self-bias",
        ],
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
        if not args.output.exists():
            raise SystemExit(f"missing committed depth gate: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed pre-SEM depth gate is stale")
        print(rendered, end="")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
