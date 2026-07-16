#!/usr/bin/env python3
"""Develop the common Belen silicon mechanism against the exposed de Boer Figure-9 data.

The original preregistered split has already been scored, so its twelve transfer markers are no
longer held out.  This script may diagnose them, but it never labels them validation data.  Model
selection remains restricted to the original single calibration marker; the remaining markers are
reported as development diagnostics.  A different experiment is required for the next validation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from deboer_2002_direct_validation import predict_depths  # noqa: E402
from deboer_feature3d import build_common_belen_si_mechanism, floor_rate  # noqa: E402
from petch.charged_surface_response_3d import (  # noqa: E402
    GrazingSpecularIonReflection3D,
)
from petch.experimental_data import load_deboer_2002_figure9_depths  # noqa: E402


EVIDENCE_SHA256 = "ed0b72235887df70552356838e376540b26234a9084ce5c886fab45ed40d7b1b"
DEVELOPMENT_STATUS_SHA256 = "d662b883330c21bdcf70be4a1f8aecde28bd2498bfef8b1e4eaeae425e66bb21"
ASPECT_RATIO_NODES = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0)


def _sha256(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _atomic_json(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _score(predictions, split):
    rows = [row for row in predictions if row["split"] == split]
    residual = np.asarray([row["residual_um"] for row in rows], dtype=float)
    return {
        "count": len(rows),
        "rmse_um": float(np.sqrt(np.mean(residual ** 2))),
        "mae_um": float(np.mean(np.abs(residual))),
        "mean_bias_um": float(np.mean(residual)),
        "within_digitization_count": int(sum(
            abs(row["residual_um"]) <= row["digitization_uncertainty_y_um"]
            for row in rows)),
    }


def _write_csv(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def _cache_key(*, s_f, dx_um, seed, aspect_ratio, tolerance, ion_energy_eV,
               iad_sigma_deg, ion_reflection):
    return (
        f"belen-v3|sF={s_f:.12g}|dx={dx_um:.12g}|seed={seed}|"
        f"AR={aspect_ratio:.12g}|tol={tolerance:.12g}|"
        f"Eion={ion_energy_eV:.12g}|iadSigma={iad_sigma_deg:.12g}|area|"
        f"ionReflection={int(bool(ion_reflection))}|response=v1-p095-a3-r090|"
        "dt=0.005|maxSteps=256|adaptiveSteps=4096|bounces=64|adaptiveBounces=512|"
        "tail=1e-8|periodic=1")


def _curve(*, s_f, dx_um, seed, tolerance, ion_energy_eV, iad_sigma_deg,
           ion_reflection, output, cache):
    mechanism = build_common_belen_si_mechanism(s_F=s_f)
    response = (GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
        1, ion_species_name="ion") if ion_reflection else None)
    response_options = (None if response is None else {
        "fixed_dt": 0.005,
        "max_steps": 256,
        "trajectory_adaptive_horizon": True,
        "trajectory_emergency_max_steps": 4096,
        "max_bounces": 64,
        "relative_tail_tolerance": 1e-8,
        "adaptive_bounce_extension": True,
        "emergency_max_bounces": 512,
        "periodic_lateral": True,
    })
    raw = []; diagnostics = []
    for aspect_ratio in ASPECT_RATIO_NODES:
        key = _cache_key(
            s_f=s_f, dx_um=dx_um, seed=seed, aspect_ratio=aspect_ratio,
            tolerance=tolerance, ion_energy_eV=ion_energy_eV,
            iad_sigma_deg=iad_sigma_deg, ion_reflection=ion_reflection)
        if key not in cache["evaluations"]:
            started = time.monotonic()
            rate, audit = floor_rate(
                aspect_ratio, mechanism, dx_um=dx_um, seed=seed,
                ion_energy_eV=ion_energy_eV, iad_sigma_deg=iad_sigma_deg,
                floor_average="area", surface_equilibration_steps=1,
                surface_fixed_point_tolerance=tolerance,
                surface_fixed_point_max_iterations=20,
                charged_surface_response=response,
                charged_surface_response_options=response_options,
                return_diagnostics=True)
            cache["evaluations"][key] = {
                "rate_m_s": rate, "diagnostics": audit,
                "elapsed_s": time.monotonic() - started,
                "completed_utc": datetime.now(timezone.utc).isoformat(),
            }
            _atomic_json(output / "rate_cache.json", cache)
            print(
                f"s_F={s_f:.3f} dx={dx_um:g} AR={aspect_ratio:g}: "
                f"{rate:.6e} m/s, fp={audit['surface_fixed_point_residual']:.3e} "
                f"in {audit['surface_fixed_point_iterations']} iterations",
                flush=True)
        item = cache["evaluations"][key]
        raw.append(float(item["rate_m_s"])); diagnostics.append(item["diagnostics"])
    raw = np.asarray(raw, dtype=float)
    return {
        "aspect_ratio": list(ASPECT_RATIO_NODES),
        "raw_rate_m_s": raw.tolist(),
        "normalized_rate": (raw / raw[0]).tolist(),
        "floor_average": "area",
        "neutral_surface_fixed_point_tolerance": tolerance,
        "ion_energy_eV": ion_energy_eV,
        "iad_component_sigma_deg": iad_sigma_deg,
        "approximate_polar_iad_fwhm_deg": 1.6025 * iad_sigma_deg,
        "ion_reflection": bool(ion_reflection),
        "charged_surface_response": (
            None if response is None else dict(response.provenance)),
        "charged_surface_response_options": response_options,
        "floor_diagnostics": diagnostics,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=(
        ROOT / "results" / "deboer_2002_belen_development"))
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=41)
    parser.add_argument("--fixed-point-tolerance", type=float, default=1e-3)
    parser.add_argument(
        "--ion-energy-e-v", type=float, default=40.0,
        help=("representative impact energy; 40 eV = Figure-9 30-V self bias plus the "
              "paper's approximately 10-V ICP plasma potential"))
    parser.add_argument(
        "--iad-component-sigma-deg", type=float, default=3.0,
        help=("transverse component sigma; 3 degrees gives approximately 4.8-degree polar "
              "FWHM, matching the primary RIE-lag paper's low-pressure scale"))
    parser.add_argument(
        "--candidate-s-f", type=float, nargs="+", default=[0.03, 0.05, 0.1, 0.15, 0.2])
    parser.add_argument(
        "--ion-reflection", action=argparse.BooleanOptionalAction, default=True,
        help=("apply the common-engine literature-bounded grazing-ion response and certified "
              "reimpact cascade; enabled by default for the repaired development path"))
    args = parser.parse_args()

    evidence = ROOT / "data" / "experimental" / "deboer_2002" / "digitized_figure9.csv"
    if _sha256(evidence) != EVIDENCE_SHA256:
        raise ValueError("de Boer evidence changed after its original score")
    development_status = (
        ROOT / "data" / "experimental" / "deboer_2002"
        / "development_status_2026-07-15.json")
    if _sha256(development_status) != DEVELOPMENT_STATUS_SHA256:
        raise ValueError("de Boer development-status declaration changed")
    rows = load_deboer_2002_figure9_depths(evidence)
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    cache_path = output / "rate_cache.json"
    cache = (json.loads(cache_path.read_text()) if cache_path.exists()
             else {"schema": "deboer-belen-development-rate-cache-v2", "evaluations": {}})

    candidates = []
    for s_f in args.candidate_s_f:
        curve = _curve(
            s_f=float(s_f), dx_um=args.dx_um, seed=args.seed,
            tolerance=args.fixed_point_tolerance, ion_energy_eV=args.ion_energy_e_v,
            iad_sigma_deg=args.iad_component_sigma_deg,
            ion_reflection=args.ion_reflection, output=output, cache=cache)
        predictions = predict_depths(rows, curve)
        candidates.append({
            "s_F": float(s_f), "curve": curve,
            "calibration_score": _score(predictions, "calibration")})
        print(
            f"candidate {s_f:.3f}: calibration error "
            f"{candidates[-1]['calibration_score']['rmse_um']:.3f} um",
            flush=True)
    selected = min(candidates, key=lambda item: (
        item["calibration_score"]["rmse_um"], item["s_F"]))
    predictions = predict_depths(rows, selected["curve"])
    development = _score(predictions, "held_out_transfer")
    result = {
        "campaign": "deboer_2002_belen_mechanism_development",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_sha256": EVIDENCE_SHA256,
        "development_status_sha256": DEVELOPMENT_STATUS_SHA256,
        "data_status": (
            "DEVELOPMENT_ONLY: all Figure-9 markers were exposed by the prior validation run; "
            "a different experiment is required for validation"),
        "selection_rule": "single original calibration marker only; smallest s_F breaks ties",
        "mechanism": "BelenSiliconSF6O2Mechanism through common feature-3d transport",
        "photon_channel": "absent; no benchmark photon flux/yield evidence",
        "numerics": {
            "dx_um": args.dx_um, "seed": args.seed,
            "aspect_ratio_nodes": list(ASPECT_RATIO_NODES),
            "floor_average": "area",
            "neutral_surface_fixed_point_tolerance": args.fixed_point_tolerance,
            "ion_energy_eV": args.ion_energy_e_v,
            "ion_energy_basis": (
                "de Boer Figure 9 declares -30 V self bias; de Boer text states an ICP plasma "
                "potential of order 10 V; representative singly charged ion energy is 40 eV"),
            "ion_energy_development_bounds_eV": [30.0, 40.0],
            "iad_component_sigma_deg": args.iad_component_sigma_deg,
            "approximate_polar_iad_fwhm_deg": 1.6025 * args.iad_component_sigma_deg,
            "iad_basis": (
                "Jansen et al. 1997 DOI 10.1016/S0167-9317(96)00142-6: low-pressure "
                "ion angular widths below about 0.1 rad and an illustrative 5-degree FWHM"),
            "iad_component_sigma_development_bounds_deg": [1.0, 3.5],
            "ion_reflection_enabled": bool(args.ion_reflection),
            "ion_reflection_basis": (
                "common-engine material-tagged grazing specular sensitivity; parameters and bounds "
                "are carried in each rate curve and are not calibrated to Figure 9"),
        },
        "candidate_results": candidates,
        "selected_s_F": selected["s_F"],
        "calibration_score": selected["calibration_score"],
        "development_transfer_score": development,
        "earned_claim": (
            "mechanism-development diagnostic only; not an experimental validation and not a "
            "predictive-chemistry claim"),
    }
    _write_csv(output / "predictions.csv", predictions)
    _atomic_json(output / "audit.json", result)
    print(json.dumps({
        "selected_s_F": result["selected_s_F"],
        "calibration_score": result["calibration_score"],
        "development_transfer_score": result["development_transfer_score"],
        "earned_claim": result["earned_claim"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
