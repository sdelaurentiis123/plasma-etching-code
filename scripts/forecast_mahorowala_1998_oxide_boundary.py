#!/usr/bin/env python3
"""Condition Mahorowala's ion boundary on oxide-mask erosion, then test Si.

No polysilicon rate, feature depth, or profile is used in the inversion.  For
each recipe, the published oxide-mask loss is matched with Chang's independent
Ar+/Cl/SiO2 beam law by varying only delivered bias power through the
deterministic RF-sheath operator.  The inferred IEAD is then passed unchanged
to the independently measured Cl+/Cl2+/Cl/poly-Si laws.

The target ions are not Ar+.  Two explicit projectile-transfer sensitivities
bound Cl2+ impact: intact energy transfer and two Cl fragments at half energy.
Neither grants a formal prediction; the comparison asks whether the oxide
observable is capable of identifying the missing boundary without touching
the silicon answer.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from petch.chang_sawin_chlorine_sio2 import (
    ChangSawinArClSiO2Mechanism,
)
from petch.chlorine_species_resolved_si import (
    SpeciesResolvedChlorineSiMechanism,
)
from petch.reactor_global import DiagnosticConditionedRFSheathTransfer
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECEIPT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "mahorowala_1998_diagnostic_conditioned_depth_projection.json"
)
SOURCE_CSV = (
    ROOT / "data" / "experimental" / "mahorowala_1998_cl2"
    / "table2_2_oxide_mask_fixed_time.csv"
)
OXIDE_MANIFEST = (
    ROOT / "data" / "experimental" / "chang_1998_figure4_14"
    / "digitization_manifest.json"
)
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "mahorowala_1998_oxide_conditioned_boundary.json"
)
DEFAULT_REPORT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "MAHOROWALA_1998_OXIDE_CONDITIONED_BOUNDARY.md"
)

ELECTRODE_AREA_M2 = 0.04
PLASMA_POTENTIAL_EV = 20.0
ETCH_DURATION_S = 75.0
BIAS_DC_SEARCH_MAXIMUM_V = 2000.0
PROJECTILE_MODES = ("cl2plus_intact", "cl2plus_two_half_energy_fragments")


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_rows():
    with SOURCE_CSV.open(newline="", encoding="utf-8") as stream:
        return {int(row["run"]): row for row in csv.DictReader(stream)}


def _sheath_transfer(row):
    return DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=ELECTRODE_AREA_M2,
        plasma_potential_eV=PLASMA_POTENTIAL_EV,
        frequency_hz=float(row["rf_sheath_frequency_hz"]),
        collapse_fraction=1.0,
        phase_count=96,
        steps_per_period=128,
        steps_per_transit=128,
        source=(
            "deterministic sinusoidal bias-component sheath; static plasma "
            "potential preserved in the zero-bias limit"
        ),
    )


def _predict_sheath(transfer, row, delivered_power_W):
    return transfer.predict(
        positive_ion_flux_m2_s={
            "Cl+": float(row["wafer_clplus_flux_m2_s"]),
            "Cl2+": float(row["wafer_cl2plus_flux_m2_s"]),
        },
        electron_temperature_eV=(
            (2.0 / 3.0) * float(row["reactor_mean_electron_energy_eV"])
        ),
        electron_density_m3=float(row["reactor_electron_density_m3"]),
        delivered_bias_power_W=float(delivered_power_W),
    )


def _project_bias(transfer, row, bias_dc_component_v):
    return transfer.project_from_bias_dc_component(
        positive_ion_flux_m2_s={
            "Cl+": float(row["wafer_clplus_flux_m2_s"]),
            "Cl2+": float(row["wafer_cl2plus_flux_m2_s"]),
        },
        electron_temperature_eV=(
            (2.0 / 3.0) * float(row["reactor_mean_electron_energy_eV"])
        ),
        electron_density_m3=float(row["reactor_electron_density_m3"]),
        bias_dc_component_v=float(bias_dc_component_v),
    )


def _oxide_fluxes(row, sheath, projectile_mode):
    cl = sheath.distributions["Cl+"]
    cl2 = sheath.distributions["Cl2+"]
    populations = [
        EnergeticFlux(
            "Ar+", cl.flux_m2_s, cl.energy_eV,
            np.ones(cl.energy_eV.shape), cl.weight,
        )
    ]
    if projectile_mode == "cl2plus_intact":
        populations.append(EnergeticFlux(
            "Ar+", cl2.flux_m2_s, cl2.energy_eV,
            np.ones(cl2.energy_eV.shape), cl2.weight,
        ))
    elif projectile_mode == "cl2plus_two_half_energy_fragments":
        populations.append(EnergeticFlux(
            "Ar+", 2.0 * cl2.flux_m2_s, 0.5 * cl2.energy_eV,
            np.ones(cl2.energy_eV.shape), cl2.weight,
        ))
    else:
        raise ValueError(f"unknown projectile mode: {projectile_mode}")
    return SurfaceFluxes(
        {"Cl": float(row["wafer_atomic_chlorine_flux_m2_s"])},
        tuple(populations),
    )


def _oxide_at_bias(mechanism, transfer, row, bias_dc_v, projectile_mode):
    sheath = _project_bias(transfer, row, bias_dc_v)
    surface = mechanism.advance(
        mechanism.initial_state(),
        _oxide_fluxes(row, sheath, projectile_mode),
        ETCH_DURATION_S,
        strict=False,
    )
    depth_nm = float(surface.etch_velocity_m_s * ETCH_DURATION_S * 1.0e9)
    return depth_nm, surface, sheath


def _invert_oxide_depth(mechanism, transfer, row, target_depth_nm, mode):
    lower_bias = 0.0
    lower_depth, lower_surface, lower_sheath = _oxide_at_bias(
        mechanism, transfer, row, lower_bias, mode
    )
    if target_depth_nm < lower_depth:
        return {
            "status": "below_zero_bias_floor",
            "delivered_bias_power_W": None,
            "bias_dc_component_v": None,
            "zero_bias_depth_nm": lower_depth,
            "maximum_search_depth_nm": None,
            "surface": lower_surface,
            "sheath": lower_sheath,
        }
    upper_bias = max(float(row["rf_sheath_bias_dc_component_v"]), 1.0)
    upper_depth, upper_surface, upper_sheath = _oxide_at_bias(
        mechanism, transfer, row, upper_bias, mode
    )
    while (
        upper_depth < target_depth_nm
        and upper_bias < BIAS_DC_SEARCH_MAXIMUM_V
    ):
        upper_bias = min(2.0 * upper_bias, BIAS_DC_SEARCH_MAXIMUM_V)
        upper_depth, upper_surface, upper_sheath = _oxide_at_bias(
            mechanism, transfer, row, upper_bias, mode
        )
    if upper_depth < target_depth_nm:
        return {
            "status": "above_search_ceiling",
            "delivered_bias_power_W": None,
            "bias_dc_component_v": None,
            "zero_bias_depth_nm": lower_depth,
            "maximum_search_depth_nm": upper_depth,
            "surface": upper_surface,
            "sheath": upper_sheath,
        }
    def residual(bias_dc_v):
        depth, _, _ = _oxide_at_bias(
            mechanism, transfer, row, bias_dc_v, mode
        )
        return depth - target_depth_nm

    inferred_bias = float(brentq(
        residual,
        lower_bias,
        upper_bias,
        xtol=1.0e-7,
        rtol=1.0e-10,
        maxiter=32,
    ))
    depth, surface, sheath = _oxide_at_bias(
        mechanism, transfer, row, inferred_bias, mode
    )
    return {
        "status": "solved",
        "delivered_bias_power_W": sheath.delivered_bias_power_W,
        "bias_dc_component_v": inferred_bias,
        "zero_bias_depth_nm": lower_depth,
        "maximum_search_depth_nm": upper_depth,
        "matched_depth_nm": depth,
        "surface": surface,
        "sheath": sheath,
    }


def _silicon_surface(row, sheath, *, sicl2_flux_m2_s=0.0):
    mechanism = SpeciesResolvedChlorineSiMechanism(
        strict_by_default=False
    )
    cl = sheath.distributions["Cl+"]
    cl2 = sheath.distributions["Cl2+"]
    result = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes(
            {
                "Cl": float(row["wafer_atomic_chlorine_flux_m2_s"]),
                "SiCl2": float(sicl2_flux_m2_s),
            },
            (
                EnergeticFlux(
                    "Cl+", cl.flux_m2_s, cl.energy_eV,
                    np.ones(cl.energy_eV.shape), cl.weight,
                ),
                EnergeticFlux(
                    "Cl2+", cl2.flux_m2_s, cl2.energy_eV,
                    np.ones(cl2.energy_eV.shape), cl2.weight,
                ),
            ),
        ),
        ETCH_DURATION_S,
        strict=False,
    )
    return result


def _mode_record(row, observed_oxide_depth_nm, inversion, mode):
    surface = inversion["surface"]
    sheath = inversion["sheath"]
    record = {
        "projectile_mode": mode,
        "oxide_inversion_status": inversion["status"],
        "observed_oxide_mask_loss_nm": observed_oxide_depth_nm,
        "zero_bias_predicted_oxide_loss_nm": inversion["zero_bias_depth_nm"],
        "inferred_delivered_bias_power_W": inversion["delivered_bias_power_W"],
        "inferred_bias_dc_component_v": inversion["bias_dc_component_v"],
        "applied_rf_bias_power_setpoint_W": float(row["rf_bias_power_W"]),
        "inferred_to_setpoint_power_ratio": None,
        "inferred_power_physically_available_from_setpoint": False,
        "oxide_energy_cards_formally_cover_iead": False,
        "oxide_ar_projectile_transfer_measured": False,
        "oxide_surface_chlorination_fraction": float(
            surface.chlorination_fraction
        ),
        "oxide_mean_yield_formula_per_projected_ar_ion": float(
            surface.mean_yield_sio2_formula_per_ion
        ),
    }
    if inversion["status"] != "solved":
        record["maximum_search_depth_nm"] = inversion[
            "maximum_search_depth_nm"
        ]
        return record
    inferred = float(inversion["delivered_bias_power_W"])
    setpoint = float(row["rf_bias_power_W"])
    cl = sheath.distributions["Cl+"]
    cl2 = sheath.distributions["Cl2+"]
    record.update({
        "matched_oxide_mask_loss_nm": inversion["matched_depth_nm"],
        "inferred_to_setpoint_power_ratio": inferred / setpoint,
        "inferred_power_physically_available_from_setpoint": (
            inferred <= setpoint * (1.0 + 1.0e-6)
        ),
        "rf_sheath_bias_dc_component_v": sheath.bias_dc_component_v,
        "rf_sheath_amplitude_v": sheath.sheath_rf_amplitude_v,
        "clplus_mean_energy_eV": cl.mean_energy_eV,
        "cl2plus_mean_energy_eV": cl2.mean_energy_eV,
        "clplus_probability_inside_oxide_70_100eV": cl.probability_inside(
            70.0, 100.0
        ),
        "cl2plus_probability_inside_oxide_70_100eV": cl2.probability_inside(
            70.0, 100.0
        ),
    })
    if mode == "cl2plus_two_half_energy_fragments":
        cl2_oxide_probability = float(np.sum(
            cl2.weight
            * (0.5 * cl2.energy_eV >= 70.0)
            * (0.5 * cl2.energy_eV <= 100.0)
        ))
        record["cl2plus_probability_inside_oxide_70_100eV"] = (
            cl2_oxide_probability
        )
    record["oxide_energy_cards_formally_cover_iead"] = bool(
        record["clplus_probability_inside_oxide_70_100eV"] >= 0.999
        and record["cl2plus_probability_inside_oxide_70_100eV"] >= 0.999
    )

    clean = _silicon_surface(row, sheath)
    clean_depth = float(clean.etch_velocity_m_s * ETCH_DURATION_S * 1.0e9)
    frozen_reflective_flux = float(row["sicl2_reflective_wall_flux_m2_s"])
    frozen_product = _silicon_surface(
        row, sheath, sicl2_flux_m2_s=frozen_reflective_flux
    )
    frozen_depth = float(
        frozen_product.etch_velocity_m_s * ETCH_DURATION_S * 1.0e9
    )
    observed_si = row["observed_feature_depth_nm"]
    record.update({
        "predicted_clean_surface_plane_si_depth_nm": clean_depth,
        "clean_surface_plane_si_signed_error_percent": (
            None if observed_si is None
            else 100.0 * (clean_depth / float(observed_si) - 1.0)
        ),
        "frozen_old_reflective_sicl2_flux_m2_s": frozen_reflective_flux,
        "frozen_product_surface_plane_si_depth_nm": frozen_depth,
        "frozen_product_surface_plane_si_signed_error_percent": (
            None if observed_si is None
            else 100.0 * (frozen_depth / float(observed_si) - 1.0)
        ),
        "frozen_product_flux_is_self_consistent_with_new_iead": False,
        "silicon_depth_used_for_inference": False,
    })
    return record


def run(receipt_path: Path):
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    source = _source_rows()
    oxide_mechanism = ChangSawinArClSiO2Mechanism()
    rows = []
    for row in receipt["rows"]:
        run_id = int(row["run"])
        source_row = source[run_id]
        observed_oxide = (
            None if not source_row["derived_oxide_removed_nm"]
            else float(source_row["derived_oxide_removed_nm"])
        )
        output = {
            "run": run_id,
            "quantitative_status": source_row["quantitative_status"],
            "inductive_power_W": float(row["inductive_power_W"]),
            "rf_bias_power_W": float(row["rf_bias_power_W"]),
            "cl2_flow_sccm": float(row["cl2_flow_sccm"]),
            "observed_feature_si_depth_nm": row["observed_feature_depth_nm"],
            "observed_oxide_mask_loss_nm": observed_oxide,
            "modes": {},
        }
        if observed_oxide is not None:
            transfer = _sheath_transfer(row)
            for mode in PROJECTILE_MODES:
                inversion = _invert_oxide_depth(
                    oxide_mechanism, transfer, row, observed_oxide, mode
                )
                output["modes"][mode] = _mode_record(
                    row, observed_oxide, inversion, mode
                )
        rows.append(output)

    summary = {}
    for mode in PROJECTILE_MODES:
        solved = [
            row["modes"][mode] for row in rows if mode in row["modes"]
            and row["modes"][mode]["oxide_inversion_status"] == "solved"
        ]
        clean_error = [
            value["clean_surface_plane_si_signed_error_percent"]
            for value in solved
            if value["clean_surface_plane_si_signed_error_percent"] is not None
        ]
        frozen_error = [
            value["frozen_product_surface_plane_si_signed_error_percent"]
            for value in solved
            if value["frozen_product_surface_plane_si_signed_error_percent"]
            is not None
        ]
        summary[mode] = {
            "solved_oxide_inversion_count": len(solved),
            "inferred_power_not_above_setpoint_count": sum(
                value["inferred_power_physically_available_from_setpoint"]
                for value in solved
            ),
            "iead_fully_inside_measured_oxide_energy_cards_count": sum(
                value["oxide_energy_cards_formally_cover_iead"]
                for value in solved
            ),
            "clean_surface_plane_si_mape_percent": (
                None if not clean_error else float(np.mean(np.abs(clean_error)))
            ),
            "frozen_product_surface_plane_si_mape_percent": (
                None if not frozen_error else float(np.mean(np.abs(frozen_error)))
            ),
        }
    return {
        "schema": "petch.mahorowala-1998-oxide-conditioned-boundary.v1",
        "claim_class": "cross-material boundary-identification diagnostic",
        "source_receipt": str(receipt_path),
        "source_receipt_sha256": _hash(receipt_path),
        "source_table": str(SOURCE_CSV),
        "source_table_sha256": _hash(SOURCE_CSV),
        "oxide_angular_digitization_manifest": str(OXIDE_MANIFEST),
        "oxide_angular_digitization_manifest_sha256": _hash(OXIDE_MANIFEST),
        "conditioning_observable": "published simultaneous oxide-mask loss",
        "held_out_observable": "published polysilicon feature depth",
        "polysilicon_depth_used_for_conditioning": False,
        "projectile_transfer_is_formally_validated": False,
        "rows": rows,
        "summary": summary,
        "formal_knobs_to_depth_pass": False,
        "evidence_blockers": [
            "measured target-tool species-resolved IEAD or bias waveform/self-bias",
            "mass-selected Cl+ and Cl2+ oxide yields including impact fragmentation",
            "independent uncertainty for the small oxide-mask erosion measurements",
            "self-consistent target-tool SiClx product-ion and neutral-product boundary",
        ],
    }


def _report(payload):
    lines = [
        "# Mahorowala 1998 oxide-conditioned wafer boundary",
        "",
        "The published oxide-mask loss conditions delivered bias power through "
        "an independent direct-beam SiO2 law. Polysilicon depth is held out.",
        "",
        "| transfer | oxide inversions | power <= setpoint | IEAD inside 70--100 eV | clean-Si MAPE | frozen-product Si MAPE |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for mode in PROJECTILE_MODES:
        value = payload["summary"][mode]
        lines.append(
            f"| {mode} | {value['solved_oxide_inversion_count']} | "
            f"{value['inferred_power_not_above_setpoint_count']} | "
            f"{value['iead_fully_inside_measured_oxide_energy_cards_count']} | "
            f"{value['clean_surface_plane_si_mape_percent']:.2f}% | "
            f"{value['frozen_product_surface_plane_si_mape_percent']:.2f}% |"
        )
    lines.extend([
        "",
        "This is an evidence-gated diagnostic, not a certified depth pass. The "
        "oxide law was measured with Ar+, while the reactor delivers Cl+ and "
        "Cl2+; both explicit transfer limits remain unvalidated.",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    arguments = parser.parse_args()
    payload = run(arguments.receipt)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    arguments.report.write_text(_report(payload), encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
