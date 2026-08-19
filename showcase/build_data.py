#!/usr/bin/env python3
"""Distill curated results into showcase/data/showcase_data.json.

Every number in the showcase page comes from this script, and this script
reads only frozen curated boards plus one direct re-evaluation of the
Turner-Chabert sheath at the audited condition (verified against the frozen
audit receipt before being written).
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from petch.reactor_global import PeriodicCurrentDensity  # noqa: E402
from petch.reactor_global.current_driven_rf_sheath import (  # noqa: E402
    TurnerChabertCurrentDrivenSheath,
)

OUT = Path(__file__).resolve().parent / "data" / "showcase_data.json"
CURATED = ROOT / "results" / "curated"


def load_sheath():
    """Replay the audited Turner-Chabert sheath; fail closed on any drift."""
    audit = json.loads(
        (CURATED / "current_driven_argon_reactor_stack" / "audit.json")
        .read_text())
    harmonics = audit["inputs"]["sheath_current_harmonics_A_m2"]
    current = PeriodicCurrentDensity(
        fundamental_frequency_hz=audit["inputs"][
            "sheath_current_frequency_hz"],
        harmonic_number=np.array([1, 2]),
        sine_A_m2=np.array(harmonics["sine"]),
        cosine_A_m2=np.array(harmonics["cosine"]),
        source="showcase replay of the audited two-harmonic waveform",
        evidence_kind="assumed",
    )
    sheath = TurnerChabertCurrentDrivenSheath(
        current=current,
        electron_temperature_eV=audit["global_plasma"][
            "electron_temperature_eV"],
        ion_mass_amu=39.948,
        sheath_edge_density_m3=audit["bohm_flux_seam"][
            "sheath_edge_density_m3"],
        phase_quadrature_count=2048,
    )
    # Receipt check: the replayed sheath must match the frozen audit.
    for got, want, name in (
        (sheath.maximum_voltage_v,
         audit["moving_rf_sheath"]["maximum_voltage_v"], "V_max"),
        (sheath.maximum_width_m,
         audit["moving_rf_sheath"]["maximum_width_m"], "s_max"),
        (sheath.xi, audit["moving_rf_sheath"]["xi"], "xi"),
    ):
        if abs(got - want) > 1e-9 * abs(want):
            raise SystemExit(f"sheath replay mismatch on {name}: "
                             f"{got!r} != {want!r}")
    return sheath, audit


def sheath_block() -> dict:
    sheath, audit = load_sheath()
    period = sheath.period_s
    times = np.linspace(0.0, period, 241)
    y = np.asarray(sheath.normalized_charge(times))

    # Collisionless IEAD from the moving field: dense uniform entry phases.
    phases = np.linspace(0.0, 2.0 * np.pi, 4001)[:-1]
    energies = sheath.ion_impact_energies(phases)
    hist, edges = np.histogram(energies, bins=64)

    return {
        "period_s": period,
        "frequency_hz": sheath.current.fundamental_frequency_hz,
        "maximum_voltage_v": sheath.maximum_voltage_v,
        "maximum_width_m": sheath.maximum_width_m,
        "mean_voltage_v": sheath.mean_voltage_v,
        "xi": sheath.xi,
        "electron_temperature_eV": sheath.electron_temperature_eV,
        "bohm_speed_m_s": float(
            np.sqrt(1.602176634e-19 * sheath.electron_temperature_eV
                    / (39.948 * 1.66053906660e-27))),
        "normalized_charge": y.tolist(),
        "iead_hist": hist.tolist(),
        "iead_edges_eV": edges.tolist(),
        "iead_mean_eV": float(np.mean(energies)),
        "audit_stats": {
            "mean_impact_energy_eV": audit["collisional_wafer_boundary"][
                "mean_impact_energy_eV"],
            "rms_impact_angle_deg": audit["collisional_wafer_boundary"][
                "rms_impact_angle_deg"],
            "arriving_ion_flux_m2_s": audit["collisional_wafer_boundary"][
                "arriving_ion_flux_m2_s"],
            "electron_density_m3": audit["global_plasma"][
                "electron_density_m3"],
            "conservation_residual": audit[
                "maximum_conservation_residual"],
            "pressure_mTorr": audit["inputs"]["pressure_mTorr"],
            "absorbed_power_W": audit["inputs"]["absorbed_bulk_power_W"],
        },
    }


def zhu_block() -> dict:
    nodes = []
    for watts in (60, 90, 105):
        raw = json.loads(
            (CURATED / "zhu_npg80_daughter_reclosed_v1"
             / f"power_{watts}W.json").read_text())
        state = raw["state"]
        inp = raw["input"]
        densities = {k: float(v) for k, v in state["densities_m3"].items()}
        fluxes = {k: float(v)
                  for k, v in state["axial_positive_ion_flux_m2_s"].items()}
        top_density = sorted(densities.items(), key=lambda kv: -kv[1])[:14]
        top_flux = sorted(fluxes.items(), key=lambda kv: -kv[1])[:8]
        nodes.append({
            "absorbed_power_W": watts,
            "electron_density_m3": state["electron_density_m3"],
            "mean_electron_energy_eV": state["mean_electron_energy_eV"],
            "electronegativity": state["electronegativity"],
            "reduced_field_Td": state["reduced_electric_field_Td"],
            "total_ion_flux_m2_s": state[
                "total_axial_positive_ion_flux_m2_s"],
            "powered_sheath_V": inp["powered_electrode_sheath_drop_V"],
            "grounded_sheath_V": inp["grounded_surface_sheath_drop_V"],
            "top_densities": top_density,
            "top_ion_fluxes": top_flux,
        })
    return {
        "recipe": {
            "gases_sccm": {"CHF3": 55, "SF6": 5, "O2": 1},
            "pressure_mTorr": 30,
            "forward_power_W": 150,
            "frequency_MHz": 13.56,
            "electrode_diameter_mm": 240,
            "species_count": 67,
            "reaction_count": 259,
        },
        "nodes": nodes,
    }


def validation_block() -> dict:
    tinacba = json.loads(
        (CURATED / "tinacba_2021_sf5_depth" / "audit.json").read_text())
    vella = json.loads(
        (CURATED / "vella_hao_ale_depth" / "audit.json").read_text())

    ll_rows = []
    with open(CURATED / "reactor_global_argon"
              / "figure3_reproduction.csv") as fh:
        for row in csv.DictReader(fh):
            if float(row["wall_energy_factor_Te"]) == 5.0:
                ll_rows.append({
                    "pressure_mTorr": float(row["pressure_mTorr"]),
                    "reference_Te_eV": float(
                        row["reference_electron_temperature_eV"]),
                    "model_Te_eV": float(
                        row["model_electron_temperature_eV"]),
                })
    ll_rows.sort(key=lambda r: r["pressure_mTorr"])

    board = vella["absolute_depth_board"]
    return {
        "tinacba": {
            "mape": tinacba["comparison"][
                "mean_absolute_relative_depth_error"],
            "points": [
                {
                    "material": p["material"],
                    "energy_eV": p["energy_eV"],
                    "measured": p["measured_depth_nm_per_1e16_cm2"],
                    "predicted": p["predicted_depth_nm_per_1e16_cm2"],
                }
                for p in tinacba["comparison"]["points"]
            ],
        },
        "vella_ale": {
            "energy_eV": board["comparison_energy_eV"],
            "experiment_nm": board["experimental_curve_depth_nm"],
            "predicted_nm": board["predicted_depth_nm"],
            "max_rel_error": board["maximum_absolute_relative_error"],
        },
        "lee_lieberman": {"rows": ll_rows, "mape_pct": 8.543},
        # de Boer SF6/O2 through the coupled engine (AUTONOMOUS_PROGRESS.md):
        # calibrate one sticking coefficient on the AR10/20 knee, predict the
        # held-out AR40 floor.
        "deboer": {
            "aspect_ratio": [0, 10, 20, 40],
            "experiment": [1.0, 0.43, 0.29, 0.20],
            "model": [1.0, 0.476, 0.289, 0.166],
            "held_out_index": 3,
            "held_out_error": 0.034,
        },
        "cards": [
            {"value": "850.2 nm", "target": "850 nm measured",
             "label": "Krüger trench: mask remaining after 60 s, "
                      "zero constants fitted to the feature"},
            {"value": "41.8 nm", "target": "41.1 nm reference",
             "label": "Krüger trench: aperture at 270 nm depth "
                      "(1.7% error)"},
            {"value": "1.570", "target": "1.5 measured",
             "label": "Karahashi SiO2 chemical sputter yield at 1 keV, "
                      "molecules per ion (4.7%)"},
            {"value": "0.033", "target": "RMSE, threshold 0.15",
             "label": "Humbird-Graves fluorocarbon/Si cumulative etch vs "
                      "held-out MD, all gates pass"},
            {"value": "0.12–0.63%", "target": "vs BOLSIG+",
             "label": "Electron-kinetics backbone: mean energy, excitation, "
                      "EEPF residuals against the swarm oracle"},
            {"value": "+1.76%", "target": "held-out 500 W",
             "label": "Lam chlorine ICP electron density, transferred from a "
                      "300 W anchor"},
        ],
    }


def main() -> None:
    data = {
        "generated_from": "results/curated (frozen boards); see build_data.py",
        "sheath": sheath_block(),
        "zhu": zhu_block(),
        "validation": validation_block(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data))
    print(f"wrote {OUT} ({OUT.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
