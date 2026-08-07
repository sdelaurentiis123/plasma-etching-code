#!/usr/bin/env python3
"""Build the no-depth-fit Vella–Hao Si/Cl2/Ar+ ALE validation board."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image

from petch import ale
from petch.interaction_data import load_kounis_melas_2024_tables
from petch.si_cl_ale_depth import (
    DEEPMD_AR_IMPACTS_PER_ALE_CYCLE,
    DEEPMD_CL2_IMPACTS_PER_ALE_CYCLE,
    DEEPMD_SI_ATOMS_PER_MATERIAL_ML,
    VellaHaoAleBoundary,
    predict_vella_hao_ale_depth,
    printed_rom_chlorine_creation_per_ar,
)


PAPER_SHA256 = (
    "789bf50302fc2ed9175403c47f895e1fb8db5481be5fb980739cad97f33c2218")
PAPER_URL = "https://www.osti.gov/servlets/purl/2248044"
FIGURE_DATA_SHA256 = (
    "5e041f07c1423a55312940f8e497cce62ddf38b94a3be1be9c5be12b3fcb9d1b")
RENDER_DPI = 240
PAPER_PAGE = 17
SOURCE_FRAME = {"left": 657.0, "right": 1450.0, "top": 268.0, "bottom": 827.0}


def _load_experiment(path):
    payload = Path(path).read_bytes()
    digest = sha256(payload).hexdigest()
    if digest != FIGURE_DATA_SHA256:
        raise ValueError(f"Figure 8 data checksum mismatch: {digest}")
    with Path(path).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    expected = [
        "Point", "Marker center x (pixel)", "Marker center y (pixel)",
        "Mean ion energy (eV)", "Etch per cycle (nm)",
    ]
    if not rows or list(rows[0]) != expected:
        raise ValueError("unexpected Figure 8 digitization schema")
    return {
        "point": np.asarray([int(row["Point"]) for row in rows]),
        "x_pixel": np.asarray([
            float(row["Marker center x (pixel)"]) for row in rows]),
        "y_pixel": np.asarray([
            float(row["Marker center y (pixel)"]) for row in rows]),
        "energy_eV": np.asarray([
            float(row["Mean ion energy (eV)"]) for row in rows]),
        "depth_nm": np.asarray([
            float(row["Etch per cycle (nm)"]) for row in rows]),
        "sha256": digest,
    }


def _render_page(pdf_path, output_prefix):
    subprocess.run(
        [
            "pdftoppm",
            "-f", str(PAPER_PAGE),
            "-l", str(PAPER_PAGE),
            "-singlefile",
            "-png",
            "-r", str(RENDER_DPI),
            str(pdf_path),
            str(output_prefix),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return output_prefix.with_suffix(".png")


def _line_center(score, expected, half_width=5):
    start = int(round(expected)) - half_width
    stop = int(round(expected)) + half_width + 1
    local = np.asarray(score[start:stop])
    peak = np.max(local)
    # Anti-aliasing can change one border row by a single black pixel.  Treat
    # the contiguous three-pixel vector stroke as one line rather than choosing
    # its darkest edge.
    indices = np.flatnonzero(local >= peak - 1) + start
    return float(np.mean(indices))


def _vision_audit(pdf_path, experiment):
    digest = sha256(Path(pdf_path).read_bytes()).hexdigest()
    if digest != PAPER_SHA256:
        raise ValueError(f"Vella–Hao manuscript checksum mismatch: {digest}")
    with tempfile.TemporaryDirectory(prefix="petch-vella-vision-") as temporary:
        image_path = _render_page(
            Path(pdf_path), Path(temporary) / "figure8-page")
        rgb = np.asarray(Image.open(image_path).convert("RGB"))
    black = np.all(rgb < 60, axis=2)
    row_score = black[:, 650:1460].sum(axis=1)
    column_score = black[260:835, :].sum(axis=0)
    frame = {
        "left": _line_center(column_score, SOURCE_FRAME["left"]),
        "right": _line_center(column_score, SOURCE_FRAME["right"]),
        "top": _line_center(row_score, SOURCE_FRAME["top"]),
        "bottom": _line_center(row_score, SOURCE_FRAME["bottom"]),
    }
    if any(abs(frame[name] - SOURCE_FRAME[name]) > 0.25 for name in frame):
        raise RuntimeError(f"source Figure 8 frame moved: {frame}")

    recovered = []
    for expected_x, expected_y in zip(
            experiment["x_pixel"], experiment["y_pixel"]):
        x_center = int(round(expected_x))
        y_center = int(round(expected_y))
        x0, x1 = x_center - 9, x_center + 10
        y0, y1 = y_center - 9, y_center + 10
        y_local, x_local = np.nonzero(black[y0:y1, x0:x1])
        if x_local.size < 100:
            raise RuntimeError("failed to recover a Figure 8 black square")
        recovered.append((x0 + np.mean(x_local), y0 + np.mean(y_local)))
    recovered = np.asarray(recovered)
    pixel_error = np.sqrt(
        (recovered[:, 0] - experiment["x_pixel"]) ** 2
        + (recovered[:, 1] - experiment["y_pixel"]) ** 2)
    pixel_tolerance = 1.25
    x_scale = 250.0 / (frame["right"] - frame["left"])
    y_scale = 8.0 / (frame["bottom"] - frame["top"])
    recovered_energy = (
        (recovered[:, 0] - frame["left"]) * x_scale)
    recovered_depth = (
        (frame["bottom"] - recovered[:, 1]) * y_scale)
    return {
        "status": "passed",
        "method": (
            "Poppler 240-dpi render; plot borders recovered from black-pixel "
            "line scores; PIL black-mask centroids recovered independently "
            "inside each 19-by-19-pixel marker window"),
        "paper_url": PAPER_URL,
        "paper_sha256": digest,
        "paper_page": PAPER_PAGE,
        "render_dpi": RENDER_DPI,
        "plot_frame_pixels": frame,
        "marker_centers_pixels": recovered.tolist(),
        "marker_energies_eV": recovered_energy.tolist(),
        "marker_depths_nm": recovered_depth.tolist(),
        "maximum_marker_center_error_pixels": float(np.max(pixel_error)),
        "acceptance_tolerance_pixels": pixel_tolerance,
        "all_markers_match_within_1p25_pixels": bool(
            np.all(pixel_error <= pixel_tolerance)),
        "digitization_bound_eV": 1.25 * x_scale,
        "digitization_bound_nm": 1.25 * y_scale,
    }


def _rom_epc_angstrom(energy_eV, *, conservative_sicl2, dt):
    theta_top = 0.0
    theta_mixed = 0.0
    history = []
    for _ in range(8):
        theta_top, _ = ale.modification_step(theta_top, 2.0)
        yields = ale.yields(energy_eV)
        steps = int(round(3.0 / dt))
        step = 3.0 / steps
        removed_si_cm2 = 0.0
        for _ in range(steps):
            top_sink = (
                yields["Y_Cl"] * theta_top
                + yields["Y_SiCl"] * theta_top
                + 2.0 * yields["Y_SiCl2"] * theta_top ** 2
                + ale.K_CL_MIX * theta_top
            )
            mixed_sink = (
                yields["Y_Cl"] * theta_mixed
                + yields["Y_SiCl"] * theta_mixed
                + 2.0 * yields["Y_SiCl2"] * theta_mixed ** 2
                - ale.K_CL_MIX * theta_top
            )
            d_top_dt = ale.J_AR / ale.SIGMA1 * top_sink
            sigma_mixed = yields["s2_ratio"] * ale.SIGMA1
            d_mixed_dt = ale.J_AR / sigma_mixed * mixed_sink
            if conservative_sicl2:
                sicl2_coverage = theta_top ** 2 + theta_mixed ** 2
            else:
                sicl2_coverage = theta_top + theta_mixed
            product_yield = (
                yields["Y_SiCl"] * (theta_top + theta_mixed)
                + yields["Y_SiCl2"] * sicl2_coverage
                + yields["Y_Si"] * (1.0 - theta_mixed)
            )
            removed_si_cm2 += step * ale.J_AR * product_yield
            theta_top = max(0.0, theta_top - step * d_top_dt)
            theta_mixed = min(
                1.0, max(0.0, theta_mixed - step * d_mixed_dt))
        history.append(removed_si_cm2 / ale.N_SI * 1.0e8)
    return float(np.mean(history[-3:]))


def _rom_atom_balance_audit(first_experimental_energy_eV):
    energy = float(first_experimental_energy_eV)
    yield_sicl2 = float(ale.yields(energy)["Y_SiCl2"])
    printed_fine = _rom_epc_angstrom(
        energy, conservative_sicl2=False, dt=5.0e-4)
    printed_coarse = _rom_epc_angstrom(
        energy, conservative_sicl2=False, dt=1.0e-3)
    conservative_fine = _rom_epc_angstrom(
        energy, conservative_sicl2=True, dt=5.0e-4)
    conservative_coarse = _rom_epc_angstrom(
        energy, conservative_sicl2=True, dt=1.0e-3)
    return {
        "status": "failed_elemental_conservation",
        "source_equation_conflict": (
            "coverage equations consume the SiCl2 channel as 2*Y*theta^2 "
            "Cl atoms, while the printed product equation emits "
            "2*Y*theta Cl atoms"),
        "chlorine_created_per_ar_at_theta_top_0p5_theta_mixed_0p5": (
            printed_rom_chlorine_creation_per_ar(0.5, 0.5, yield_sicl2)),
        "first_experimental_energy_eV": energy,
        "experiment_schedule_s": {"cl2": 2.0, "ar_ion": 3.0},
        "printed_rom_depth_nm": printed_fine / 10.0,
        "atom_conservative_theta_squared_depth_nm": (
            conservative_fine / 10.0),
        "time_step_refinement_difference_nm": {
            "printed": abs(printed_fine - printed_coarse) / 10.0,
            "conservative": abs(
                conservative_fine - conservative_coarse) / 10.0,
        },
        "interpretation": (
            "the legacy ROM may replay the paper's empirical curve but cannot "
            "serve as an atom-conservative absolute-depth closure"),
    }


def build_audit(paper_pdf=None):
    root = Path(__file__).resolve().parents[1]
    source_directory = (
        root / "data" / "surface_interactions" / "kounis_melas_2024")
    experiment_path = (
        root / "data" / "experimental" / "vella_hao_2023"
        / "figure8_epc.csv")
    tables = load_kounis_melas_2024_tables(source_directory)
    experiment = _load_experiment(experiment_path)
    comparison_energies = np.asarray([60.0, 80.0, 100.0])
    predictions = [
        predict_vella_hao_ale_depth(
            energy, tables.ale_cycles, tables.sputtering)
        for energy in comparison_energies
    ]
    predicted_depth = np.asarray([
        item.total_depth_nm for item in predictions])
    observed_depth = np.interp(
        comparison_energies,
        experiment["energy_eV"],
        experiment["depth_nm"],
    )
    error = predicted_depth - observed_depth
    relative_error = error / observed_depth
    boundary = VellaHaoAleBoundary()
    vision = (
        _vision_audit(Path(paper_pdf), experiment)
        if paper_pdf is not None
        else {
            "status": "not_run",
            "paper_url": PAPER_URL,
            "required_paper_sha256": PAPER_SHA256,
        }
    )
    return {
        "schema_version": 1,
        "claim": (
            "no-depth-fit absolute Si/Cl2/Ar+ ALE cross-source validation; "
            "not a blind prediction and not a measured species-resolved IEAD"),
        "source": {
            "experiment_paper_doi": "10.1002/ppap.202200198",
            "experiment_osti_id": "2248044",
            "experiment_figure": "Figure 8",
            "digitization_sha256": experiment["sha256"],
            "deepmd_dataset_doi": tables.ale_cycles.provenance[
                "dataset_doi"],
            "deepmd_archive_sha256": tables.ale_cycles.provenance[
                "archive_sha256"],
            "deepmd_ale_trajectory_sha256": tables.ale_cycles.provenance[
                "conditions"]["parent_source_table_sha256"],
            "ale_cycle_table_fingerprint": tables.ale_cycles.fingerprint,
            "sputter_table_fingerprint": tables.sputtering.fingerprint,
        },
        "boundary_and_units": {
            "measured_positive_ion_flux_cm2_s": (
                boundary.positive_ion_flux_cm2_s),
            "ion_bombardment_duration_s": (
                boundary.ion_bombardment_duration_s),
            "experimental_ion_fluence_cm2": boundary.ion_fluence_cm2,
            "deepmd_ar_impacts_per_cycle": (
                DEEPMD_AR_IMPACTS_PER_ALE_CYCLE),
            "deepmd_cl2_impacts_per_cycle": (
                DEEPMD_CL2_IMPACTS_PER_ALE_CYCLE),
            "deepmd_si_atoms_per_material_ml": (
                DEEPMD_SI_ATOMS_PER_MATERIAL_ML),
            "cell_area_cm2": predictions[0].simulated_cell_area_cm2,
            "deepmd_ar_fluence_cm2": predictions[0].simulated_ar_fluence_cm2,
            "remaining_physical_sputter_fluence_cm2": (
                predictions[0].physical_sputter_tail_fluence_cm2),
            "monolayer_alias_used_for_fluence_conversion": False,
            "conversion": (
                "atomistic removal uses explicit atoms/cell and cell area; "
                "experiment uses measured flux times duration"),
        },
        "absolute_depth_board": {
            "comparison_energy_eV": comparison_energies.tolist(),
            "comparison_curve_operation": (
                "piecewise-linear interpolation between adjacent digitized "
                "experimental markers; no extrapolation"),
            "fit_to_depth_used": False,
            "predicted_depth_nm": predicted_depth.tolist(),
            "experimental_curve_depth_nm": observed_depth.tolist(),
            "signed_error_nm": error.tolist(),
            "relative_error": relative_error.tolist(),
            "mae_nm": float(np.mean(np.abs(error))),
            "rmse_nm": float(np.sqrt(np.mean(error ** 2))),
            "maximum_absolute_relative_error": float(
                np.max(np.abs(relative_error))),
            "nominal_gate_maximum_relative_error": 0.15,
            "nominal_gate_passed": bool(
                np.max(np.abs(relative_error)) <= 0.15),
            "point_ledgers": [
                {
                    "energy_eV": item.mean_ion_energy_eV,
                    "steady_cycle_si_material_ml": (
                        item.steady_cycle_si_material_ml),
                    "chlorinated_transient_si_atoms_cm2": (
                        item.chlorinated_transient_si_atoms_cm2),
                    "physical_tail_si_atoms_cm2": (
                        item.physical_sputter_tail_si_atoms_cm2),
                    "atom_balance_residual_cm2": (
                        item.dimensional_atom_balance_residual_cm2),
                    "chlorinated_transient_depth_nm": (
                        item.chlorinated_transient_depth_nm),
                    "physical_tail_depth_nm": (
                        item.physical_sputter_tail_depth_nm),
                    "source_reported_sputter_tail_uncertainty_nm": (
                        item.source_reported_tail_depth_uncertainty_nm),
                }
                for item in predictions
            ],
        },
        "uncertainty_and_validity": {
            "digitization_bound_nm": 0.01788908765652952,
            "source_reported_sputter_tail_uncertainty_nm": [
                item.source_reported_tail_depth_uncertainty_nm
                for item in predictions
            ],
            "combined_uncertainty_claimed": False,
            "incomplete_terms": [
                "experimental mean-energy uncertainty",
                "experimental IEDF and ion angular distribution",
                "experimental ion-species fractions",
                "finite-cell/statistical uncertainty of ALE cycle increments",
                "model-form error of pure-Si tail after the chlorinated transient",
            ],
            "known_model_form_assumptions": [
                "normal-incidence monoenergetic DeepMD impacts represent the inferred experimental mean energy",
                "the released 1000-impact ALE sequence exhausts the chlorinated transient",
                "the remaining experimental fluence removes Si through the released pure-Si Ar+ sputter law",
                "slow between-impact chemistry and facility wall conditioning are not represented",
            ],
            "evidence_tier": (
                "retrospective no-fit cross-source transfer with a measured "
                "fluence and inferred mean-energy boundary"),
        },
        "vision_digitization_audit": vision,
        "legacy_transient_rom_atom_balance": _rom_atom_balance_audit(
            experiment["energy_eV"][0]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--paper-pdf",
        help=(
            "optional pinned Vella–Hao DOE manuscript; checksum must match "
            "before the PIL/Poppler vision gate runs"),
    )
    parser.add_argument("--output", help="write JSON to this path instead of stdout")
    arguments = parser.parse_args()
    payload = json.dumps(
        build_audit(arguments.paper_pdf), indent=2, sort_keys=True) + "\n"
    if arguments.output:
        Path(arguments.output).write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
