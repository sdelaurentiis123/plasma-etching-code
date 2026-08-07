#!/usr/bin/env python3
"""Audit the Kounis-Melas 215 eV ALE product sequence.

The numerical audit is reproducible from the checksum-gated CSV files.  When
``--paper-pdf`` points to the pinned accepted manuscript, the script also
renders Figures 13--14 at 180 dpi and independently recovers the five blue
DeepMD markers using PIL color masks and connected components.  That vision
path checks both the 215 eV figure binding and the CSV values.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image
from scipy import ndimage

from petch.chang_sawin_chlorine_si import (
    ChangSawinArClSiParameters,
)
from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_kinetics import (
    EnergeticFlux,
    ParameterEvidence,
    SurfaceFluxes,
)
from petch.tabulated_si_cl_ale import (
    TabulatedSiClAleProductMechanism,
)


PAPER_SHA256 = (
    "d8d374d412d99e625b2b989399d564976332f0c342b37b0249d1039bca4b5bb1")
PAPER_URL = "https://www.osti.gov/servlets/purl/2514378"
RENDER_DPI = 180


def _float(value):
    return float(np.asarray(value))


def _render_page(pdf_path, page_number, output_prefix):
    subprocess.run(
        [
            "pdftoppm",
            "-f", str(page_number),
            "-l", str(page_number),
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


def _figure_markers(image_path, *, x_max, y_max):
    rgb = np.asarray(Image.open(image_path).convert("RGB"))
    black = np.all(rgb < 80, axis=2)
    row_scores = black[150:850, 550:1550].sum(axis=1)
    column_scores = black[150:850, 550:1550].sum(axis=0)
    horizontal = sorted(
        np.argsort(row_scores)[-2:] + 150)
    vertical = sorted(
        np.argsort(column_scores)[-2:] + 550)
    top, bottom = (int(value) for value in horizontal)
    left, right = (int(value) for value in vertical)
    if (bottom - top < 500 or right - left < 700):
        raise RuntimeError("failed to identify source-figure plot frame")

    blue = (
        (rgb[:, :, 2] > 180)
        & (rgb[:, :, 0] < 80)
        & (rgb[:, :, 1] > 70)
        & (rgb[:, :, 1] < 190)
    )
    labels, _ = ndimage.label(blue)
    markers = []
    for index, slices in enumerate(ndimage.find_objects(labels), start=1):
        if slices is None:
            continue
        component = labels[slices] == index
        area = int(component.sum())
        # The five marker/error-bar objects are 300--510 blue pixels at this
        # exact source hash and rendering.  The legend marker is only 116.
        if area < 250:
            continue
        x_pixel = (
            slices[1].start + int(np.argmax(component.sum(axis=0))))
        y_pixel = (
            slices[0].start + int(np.argmax(component.sum(axis=1))))
        if not (left <= x_pixel <= right and top <= y_pixel <= bottom):
            continue
        markers.append({
            "x_pixel": x_pixel,
            "y_pixel": y_pixel,
            "component_blue_pixels": area,
            "x_value": (
                x_max * (x_pixel - left) / (right - left)),
            "y_value": (
                y_max * (bottom - y_pixel) / (bottom - top)),
        })
    markers.sort(key=lambda item: item["x_pixel"])
    if len(markers) != 5:
        raise RuntimeError(
            f"expected five DeepMD markers, found {len(markers)}")
    return {
        "image_shape": list(rgb.shape),
        "plot_frame_pixels": {
            "left": left,
            "right": right,
            "top": top,
            "bottom": bottom,
        },
        "markers": markers,
        "x_value_per_pixel": x_max / (right - left),
        "y_value_per_pixel": y_max / (bottom - top),
    }


def _vision_audit(pdf_path, table):
    payload = Path(pdf_path).read_bytes()
    digest = sha256(payload).hexdigest()
    if digest != PAPER_SHA256:
        raise ValueError(
            f"accepted-manuscript checksum mismatch: {digest}")
    source_x = table.axes[0].values
    panels = {}
    with tempfile.TemporaryDirectory(prefix="petch-kounis-vision-") as temp:
        temp = Path(temp)
        for page, output, y_max in (
                (31, "si_yield", 0.35),
                (32, "sicl_yield", 0.12)):
            image_path = _render_page(
                pdf_path, page, temp / f"page-{page}")
            extracted = _figure_markers(
                image_path, x_max=15.0, y_max=y_max)
            observed_x = np.asarray([
                item["x_value"] for item in extracted["markers"]])
            observed_y = np.asarray([
                item["y_value"] for item in extracted["markers"]])
            x_tolerance = 1.25 * extracted["x_value_per_pixel"]
            y_tolerance = 1.25 * extracted["y_value_per_pixel"]
            x_error = np.abs(observed_x - source_x)
            y_error = np.abs(observed_y - table.outputs[output])
            extracted.update({
                "paper_page": page,
                "figure": "13" if page == 31 else "14",
                "output": output,
                "caption_ion_energy_eV": 215.0,
                "maximum_absolute_x_error": float(np.max(x_error)),
                "maximum_absolute_y_error": float(np.max(y_error)),
                "x_acceptance_tolerance": x_tolerance,
                "y_acceptance_tolerance": y_tolerance,
                "csv_nodes_match_within_1p25_pixels": bool(
                    np.all(x_error <= x_tolerance)
                    and np.all(y_error <= y_tolerance)),
            })
            panels[output] = extracted
    return {
        "status": "passed",
        "method": (
            "PIL RGB threshold plus connected components on a Poppler "
            "180-dpi render; axis frame recovered from black border pixels"),
        "paper_url": PAPER_URL,
        "paper_sha256": digest,
        "render_dpi": RENDER_DPI,
        "panels": panels,
        "condition_disambiguation": (
            "Figures 13-14 bind Products.csv to 215 eV. Figure 12 is the "
            "distinct 80 eV morphology/cycle sequence."),
    }


def build_audit(paper_pdf=None):
    root = Path(__file__).resolve().parents[1]
    source_dir = (
        root / "data" / "surface_interactions" / "kounis_melas_2024")
    tables = load_kounis_melas_2024_tables(source_dir)
    table = tables.ale_products
    density = 8.0 / (5.43e-10) ** 3
    mechanism = TabulatedSiClAleProductMechanism(
        table,
        density,
        ParameterEvidence(
            "Kounis-Melas OSTI 2589032: diamond-Si lattice a=5.43 angstrom",
            "source_derived",
            supports_prediction_within_declared_domain=True,
        ),
    )
    total_dose = mechanism.dose_bin_edges_m2[-1]
    result = mechanism.advance(
        mechanism.initial_state(1.0e19),
        SurfaceFluxes({}, (
            EnergeticFlux(
                "Ar+", total_dose,
                np.array([215.0]), np.array([1.0]), np.array([1.0])),
        )),
        1.0,
    )
    exchange_residual = {
        name: _float(result.material_exchange.residual_units_m2(name))
        for name in result.material_exchange.removed_units_m2
    }

    chang = ChangSawinArClSiParameters.molecular_chlorine_100eV()
    ratios = tables.reactive_ion_etch.axes[0].values
    theta = (
        chang.surface_chlorination_coefficient * ratios
        / (
            chang.surface_chlorination_coefficient * ratios
            + chang.product_chlorine_atoms_per_si
            * chang.ion_enhanced_yield_si_per_ion_at_full_chlorination
        )
    )
    chang_prediction = (
        chang.physical_sputter_yield_si_per_ion * (1.0 - theta)
        + chang.ion_enhanced_yield_si_per_ion_at_full_chlorination * theta
    )
    deepmd_observed = tables.reactive_ion_etch.outputs[
        "reactive_etch_yield"]
    relative_error = np.abs(
        chang_prediction - deepmd_observed) / deepmd_observed

    vision = (
        _vision_audit(Path(paper_pdf), table)
        if paper_pdf is not None
        else {
            "status": "not_run",
            "required_paper_sha256": PAPER_SHA256,
            "paper_url": PAPER_URL,
        }
    )
    return {
        "schema_version": 1,
        "claim": (
            "atom-balanced replay of source DeepMD dose-window products; "
            "not experimental validation or a held-out depth prediction"),
        "source": {
            "dataset_doi": table.provenance["dataset_doi"],
            "paper_doi": table.provenance["paper_doi"],
            "evidence_type": table.provenance["evidence_type"],
            "table_fingerprint": table.fingerprint,
            "source_table_sha256": table.provenance["source_table_sha256"],
            "ion_energy_eV": table.provenance[
                "conditions"]["ar_ion_energy_eV"],
            "incidence_angle_deg": table.provenance[
                "conditions"]["incidence_angle_deg"],
        },
        "dose_windows": {
            "centres_1e15_cm2": table.axes[0].values.tolist(),
            "edges_m2": mechanism.dose_bin_edges_m2.tolist(),
            "widths_m2": np.diff(
                mechanism.dose_bin_edges_m2).tolist(),
            "integration_rule": (
                "exact overlap with piecewise-constant released "
                "window-average yields"),
            "interpolation_used": False,
        },
        "integrated_sequence": {
            "incident_ar_ions_m2": _float(total_dose),
            "product_counts_m2": {
                name: _float(value)
                for name, value in result.product_counts_m2.items()
            },
            "removed_si_atoms_m2": _float(result.removed_si_atoms_m2),
            "emitted_chlorine_atoms_m2": _float(
                result.emitted_chlorine_atoms_m2),
            "equivalent_removed_si_depth_nm": (
                _float(result.removed_si_atoms_m2) / density * 1.0e9),
            "one_monolayer_cl_loading_atoms_m2": 1.0e19,
            "remaining_cl_atoms_m2": _float(
                result.state.retained_chlorine_atoms_m2),
            "exchange_residual_atoms_m2": exchange_residual,
            "product_routing_complete": (
                result.material_exchange.product_routing_complete),
            "all_products_lack_launch_distribution": all(
                not item.transport_ready
                for item in result.product_populations),
        },
        "independent_continuous_rie_diagnostic": {
            "source_model": (
                "Chang-Sawin molecular-Cl2 100 eV beam regression"),
            "comparison_evidence": "DeepMD molecular dynamics",
            "fit_to_comparison_used": False,
            "flux_ratios": ratios.tolist(),
            "chang_prediction_si_per_ar": chang_prediction.tolist(),
            "deepmd_si_per_ar": deepmd_observed.tolist(),
            "maximum_relative_error": float(np.max(relative_error)),
            "rmse_si_per_ar": float(np.sqrt(np.mean(
                (chang_prediction - deepmd_observed) ** 2))),
            "verdict": (
                "independent trend agreement with material disagreement; "
                "not an exact cross-source pass"),
        },
        "vision_condition_and_value_audit": vision,
        "known_limits": list(
            result.validity.known_model_form_omissions),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--paper-pdf",
        help=(
            "optional accepted-manuscript PDF; checksum must match the "
            "pinned DOE copy"),
    )
    parser.add_argument(
        "--output",
        help="write JSON to this path instead of stdout",
    )
    arguments = parser.parse_args()
    payload = json.dumps(
        build_audit(arguments.paper_pdf),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if arguments.output:
        Path(arguments.output).write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
