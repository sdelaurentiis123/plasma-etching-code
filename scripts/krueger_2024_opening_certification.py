#!/usr/bin/env python3
"""Certify the Krueger-2024 mask-opening observable without evolving a profile.

The paper's :math:`w_m` is the minimum horizontal width of the mask opening,
including deposited material.  It is a connectivity-qualified observable: the
paper explicitly assigns ``w_m = 0`` after clogging.  A sealed gas pocket inside
the mask is therefore not a nonzero opening.

This bounded diagnostic reads an existing checkpoint and audit only.  It does
not load held-out profiles, run transport, evolve geometry, smooth visibility,
or change the production observable.  It reports three deliberately separate
quantities:

* the minimum width reconstructed from the mask material level set;
* whether exterior gas reaches the lowest resolved plane inside the mask;
* the one-step temporal flicker and cross-extrusion spread that must be carried
  as numerical evidence rather than silently filtered.

The dense vertical scan evaluates the same bilinear level-set reconstruction at
more heights.  It is a measurement refinement, not a geometry mollifier.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label


PAPER = {
    "citation": (
        "Krueger et al., J. Vac. Sci. Technol. A 42, 043008 (2024), "
        "doi:10.1116/6.0003554"
    ),
    "url": "https://cpseg.eecs.umich.edu/pub/articles/JVSTA_42_043008_2024.pdf",
    "definition": (
        "minimum horizontal width of the mask opening including deposition; "
        "zero while the opening is clogged"
    ),
    "source_locations": [
        "Sec. V and Fig. 7 (definition and geometry)",
        "Table IV (45 nm base target)",
        "Sec. VIII (w_m = 0 for fully clogged features)",
    ],
}


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def _linear_zero(a, fa, b, fb):
    denominator = float(fa) - float(fb)
    if denominator == 0.0:
        return 0.5 * (float(a) + float(b))
    return float(a) + (float(b) - float(a)) * float(fa) / denominator


def _gas_interval(field, x, center_index):
    """Return the center-containing horizontal negative interval."""
    values = np.asarray(field, dtype=float)
    if values[center_index] >= 0.0:
        return None
    left = int(center_index)
    right = int(center_index)
    while left > 0 and values[left - 1] < 0.0:
        left -= 1
    while right + 1 < len(values) and values[right + 1] < 0.0:
        right += 1
    if left == 0 or right + 1 == len(values):
        return None
    x_left = _linear_zero(
        x[left - 1], values[left - 1], x[left], values[left])
    x_right = _linear_zero(
        x[right], values[right], x[right + 1], values[right + 1])
    return float(x_left), float(x_right)


def _load_checkpoint(path: Path, mask_material_id: int):
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"]))
        name = f"material_levelset_{int(mask_material_id)}"
        if name not in archive.files:
            raise ValueError(
                f"checkpoint lacks mask material level set {mask_material_id}")
        payload = {
            "phi": np.asarray(archive["phi"], dtype=float).copy(),
            "mask_levelset": np.asarray(archive[name], dtype=float).copy(),
            "material_id": np.asarray(archive["material_id"], dtype=int).copy(),
            "dx": float(metadata["dx"]),
            "mesh_length_unit_m": float(metadata["mesh_length_unit_m"]),
            "metadata": metadata,
        }
    if (payload["phi"].ndim != 3
            or payload["mask_levelset"].shape != payload["phi"].shape
            or payload["material_id"].shape != payload["phi"].shape
            or not np.all(np.isfinite(payload["phi"]))
            or not np.all(np.isfinite(payload["mask_levelset"]))):
        raise ValueError("checkpoint geometry arrays are invalid")
    return payload


def _mask_top_z(mask_levelset, x, z, y_rows, opening_center, dx):
    outside = np.flatnonzero(np.abs(x - float(opening_center)) >= 0.055)
    values = []
    for i in outside:
        for j in y_rows:
            field = mask_levelset[i, j]
            crossing = np.flatnonzero(
                (field[:-1] >= 0.0) & (field[1:] < 0.0))
            for k in crossing:
                top = _linear_zero(
                    z[k], field[k], z[k + 1], field[k + 1])
                values.append(float(top))
    if not values:
        raise ValueError("mask top is unresolved")
    top = float(np.median(values))
    if not np.isfinite(top) or top <= dx:
        raise ValueError("mask top is invalid")
    return top


def _width_at_height(field, x, z, y_rows, center_index, z_value):
    if z_value <= z[0]:
        lower, alpha = 0, 0.0
    elif z_value >= z[-1]:
        lower, alpha = len(z) - 2, 1.0
    else:
        lower = int(np.searchsorted(z, z_value) - 1)
        alpha = float((z_value - z[lower]) / (z[lower + 1] - z[lower]))
    intervals = []
    for j in y_rows:
        slice_field = (
            (1.0 - alpha) * field[:, j, lower]
            + alpha * field[:, j, lower + 1]
        )
        interval = _gas_interval(slice_field, x, center_index)
        if interval is not None:
            intervals.append(interval)
    if not intervals:
        return None
    widths = np.asarray([right - left for left, right in intervals], dtype=float)
    return {
        "mean_width": float(np.mean(widths)),
        "per_y_widths": widths,
        "mean_left": float(np.mean([left for left, _ in intervals])),
        "mean_right": float(np.mean([right for _, right in intervals])),
    }


def _minimum_width(
        field, *, dx, substrate_top, mask_top, opening_center,
        dense_samples_per_cell):
    shape = field.shape
    x, y, z = (np.arange(size, dtype=float) * dx for size in shape)
    center_index = int(np.argmin(np.abs(x - float(opening_center))))
    y_rows = tuple(range(1, len(y) - 1)) if len(y) > 2 else tuple(range(len(y)))
    lower = float(substrate_top) + 0.25 * dx
    upper = float(mask_top) - 0.25 * dx
    if lower >= upper:
        raise ValueError("mask has no resolved interior height")

    nodal = z[(z >= lower) & (z <= upper)]
    intervals = max(1, int(np.ceil((upper - lower) / dx)))
    # Always include the native rows.  Reinitialization can place a sharp but
    # valid piecewise-bilinear throat exactly on one row; an offset dense lattice
    # must not accidentally step over it and report a wider opening.
    dense = np.unique(np.concatenate((
        nodal,
        np.linspace(
            lower, upper, intervals * int(dense_samples_per_cell) + 1),
    )))

    def best(values):
        measured = []
        for z_value in values:
            item = _width_at_height(
                field, x, z, y_rows, center_index, float(z_value))
            if item is not None:
                measured.append((item["mean_width"], float(z_value), item))
        if not measured:
            raise ValueError("mask opening is unresolved at every sampled height")
        width, location, item = min(measured, key=lambda row: row[0])
        return {
            "width_mesh_units": float(width),
            "z_mesh_units": float(location),
            "per_y_widths_mesh_units": item["per_y_widths"].tolist(),
            "cross_y_span_mesh_units": float(np.ptp(item["per_y_widths"])),
            "resolved_y_rows": int(len(item["per_y_widths"])),
        }

    return {"nodal_rows": best(nodal), "dense_bilinear_scan": best(dense)}


def _periodic_component_roots(field):
    occupied = np.asarray(field, dtype=bool)
    component, count = label(occupied)
    parent = np.arange(int(count) + 1)

    def find(index):
        index = int(index)
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    def union(left, right):
        left = find(left)
        right = find(right)
        if left and right and left != right:
            parent[right] = left

    for axis in (0, 1):
        first = np.take(component, 0, axis=axis)
        last = np.take(component, -1, axis=axis)
        selected = (first > 0) & (last > 0)
        for left, right in zip(first[selected], last[selected]):
            union(left, right)
    roots = np.zeros(int(count) + 1, dtype=int)
    for index in range(1, int(count) + 1):
        roots[index] = find(index)
    return roots[component]


def _opening_connectivity(
        phi, *, dx, substrate_top, opening_center, opening_width):
    """Test exterior connection to the lowest resolved mask-interior plane.

    Requiring gas below the substrate falsely labels a shallow but open etch as
    closed until its depth exceeds one cell.  The paper observable concerns the
    mask throat, so its correct lower witness is the first resolved plane inside
    the mask, within the declared opening region of interest.
    """
    # Checkpoint x/y arrays include duplicate periodic endpoints, matching the
    # production topology audit.  Remove them before wrapping the remaining core.
    gas = np.asarray(phi[:-1, :-1, :], dtype=float) <= 0.0
    rooted = _periodic_component_roots(gas)
    upper = {int(value) for value in np.unique(rooted[:, :, -1]) if value > 0}
    x = np.arange(rooted.shape[0], dtype=float) * dx
    z = np.arange(rooted.shape[2], dtype=float) * dx
    exit_candidates = np.flatnonzero(z > float(substrate_top) + 0.25 * dx)
    if not exit_candidates.size:
        raise ValueError("no resolved mask-interior exit plane")
    exit_index = int(exit_candidates[0])
    x_roi = np.flatnonzero(
        np.abs(x - float(opening_center))
        <= 0.5 * float(opening_width) + dx)
    if not x_roi.size:
        raise ValueError("declared opening ROI has no resolved lateral node")
    exit_roots = {
        int(value)
        for value in np.unique(rooted[x_roi, :, exit_index])
        if value > 0
    }
    shared = upper & exit_roots
    return {
        "exterior_to_mask_exit_open": bool(shared),
        "upper_boundary_gas_component_count": len(upper),
        "shared_component_count": len(shared),
        "mask_exit_z_mesh_units": float(z[exit_index]),
    }


def _timeline(audit, *, dx_nm):
    history = list(audit.get("history", []))
    open_state = True
    rows = []
    for row in history:
        event = row.get("topology_event")
        if isinstance(event, dict) and event.get("accepted"):
            if event.get("kind") == "gas_cavity_enclosed":
                open_state = False
            elif event.get("kind") == "gas_cavity_opened":
                open_state = True
        legacy = float(row["metrics"]["mask_opening_nm"])
        rows.append({
            "step": int(row["step"]),
            "physical_time_s": float(row["physical_time_s"]),
            "connectivity_open": bool(open_state),
            "legacy_mask_opening_nm": legacy,
            "paper_connectivity_qualified_opening_nm": legacy if open_state else 0.0,
            "topology_event_kind": (
                event.get("kind") if isinstance(event, dict) else None),
        })

    spikes = []
    for left, middle, right in zip(rows, rows[1:], rows[2:]):
        if (not left["connectivity_open"]
                or not middle["connectivity_open"]
                or not right["connectivity_open"]
                or middle["topology_event_kind"] is not None
                or middle["step"] != left["step"] + 1
                or right["step"] != middle["step"] + 1):
            continue
        fraction = (
            (middle["physical_time_s"] - left["physical_time_s"])
            / (right["physical_time_s"] - left["physical_time_s"])
        )
        expected = (
            left["legacy_mask_opening_nm"]
            + fraction * (
                right["legacy_mask_opening_nm"]
                - left["legacy_mask_opening_nm"])
        )
        residual = middle["legacy_mask_opening_nm"] - expected
        spikes.append({
            "step": middle["step"],
            "physical_time_s": middle["physical_time_s"],
            "signed_residual_nm": float(residual),
            "absolute_residual_nm": float(abs(residual)),
            "absolute_residual_cells": float(abs(residual) / dx_nm),
        })
    event_times = [
        row["physical_time_s"] for row in rows
        if row["topology_event_kind"] is not None
    ]
    for item in spikes:
        item["distance_to_nearest_topology_event_s"] = (
            min(abs(item["physical_time_s"] - value) for value in event_times)
            if event_times else None
        )
    maximum = max(spikes, key=lambda item: item["absolute_residual_nm"], default=None)
    largest = sorted(
        spikes, key=lambda item: item["absolute_residual_nm"], reverse=True)[:10]
    closed = [row for row in rows if not row["connectivity_open"]]
    return {
        "connectivity_inference": (
            "history starts open and toggles only on accepted gas-cavity "
            "enclosure/opening events; final checkpoint is checked directly"
        ),
        "history_count": len(rows),
        "closed_history_count": len(closed),
        "legacy_nonzero_while_closed_count": sum(
            row["legacy_mask_opening_nm"] > 0.0 for row in closed),
        "maximum_open_state_one_step_flicker": maximum,
        "largest_open_state_one_step_flickers": largest,
        "final_connectivity_open": rows[-1]["connectivity_open"] if rows else None,
    }


def build_report(
        audit_path, checkpoint_path, *, mask_material_id=2,
        substrate_top_um=None, opening_center_um=None, opening_width_um=None,
        dense_samples_per_cell=64, paper_pdf=None):
    audit_path = Path(audit_path)
    checkpoint_path = Path(checkpoint_path)
    if dense_samples_per_cell < 2:
        raise ValueError("dense_samples_per_cell must be at least two")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    geometry_config = dict(
        audit.get("configuration", {}).get("geometry", {}) or {})
    declared_substrate_top = geometry_config.get("substrate_top_um")
    declared_opening_width = geometry_config.get("opening_width_um")
    declared_center = (
        0.5 * float(geometry_config["cell_width_um"])
        if geometry_config.get("cell_width_um") is not None else None)
    for name, supplied, declared in (
            ("substrate_top_um", substrate_top_um, declared_substrate_top),
            ("opening_center_um", opening_center_um, declared_center),
            ("opening_width_um", opening_width_um, declared_opening_width)):
        if (supplied is not None and declared is not None
                and not np.isclose(float(supplied), float(declared), rtol=0.0, atol=1.0e-12)):
            raise ValueError(
                f"{name}={float(supplied):.12g} disagrees with "
                f"audit geometry {float(declared):.12g}")
    if substrate_top_um is None:
        substrate_top_um = declared_substrate_top
    if opening_width_um is None:
        opening_width_um = declared_opening_width
    if opening_center_um is None:
        opening_center_um = declared_center
    if (substrate_top_um is None or opening_center_um is None
            or opening_width_um is None):
        raise ValueError(
            "substrate top, opening center, and opening width must be provided "
            "or recoverable from audit.configuration.geometry")
    if (not np.isfinite(substrate_top_um) or substrate_top_um <= 0.0
            or not np.isfinite(opening_center_um)
            or not np.isfinite(opening_width_um) or opening_width_um <= 0.0):
        raise ValueError("declared opening geometry is invalid")
    checkpoint = _load_checkpoint(checkpoint_path, mask_material_id)
    dx = checkpoint["dx"]
    dx_nm = dx * checkpoint["mesh_length_unit_m"] * 1.0e9
    x, y, z = (
        np.arange(size, dtype=float) * dx for size in checkpoint["phi"].shape)
    y_rows = tuple(range(1, len(y) - 1)) if len(y) > 2 else tuple(range(len(y)))
    mask_top = _mask_top_z(
        checkpoint["mask_levelset"], x, z, y_rows,
        float(opening_center_um), dx)
    mask_width = _minimum_width(
        checkpoint["mask_levelset"], dx=dx,
        substrate_top=float(substrate_top_um), mask_top=mask_top,
        opening_center=float(opening_center_um),
        dense_samples_per_cell=int(dense_samples_per_cell))
    combined_width = _minimum_width(
        checkpoint["phi"], dx=dx,
        substrate_top=float(substrate_top_um), mask_top=mask_top,
        opening_center=float(opening_center_um),
        dense_samples_per_cell=int(dense_samples_per_cell))
    connectivity = _opening_connectivity(
        checkpoint["phi"], dx=dx, substrate_top=float(substrate_top_um),
        opening_center=float(opening_center_um),
        opening_width=float(opening_width_um))
    timeline = _timeline(audit, dx_nm=dx_nm)

    scale = checkpoint["mesh_length_unit_m"] * 1.0e9
    for family in (mask_width, combined_width):
        for item in family.values():
            item["width_nm"] = item["width_mesh_units"] * scale
            item["z_um"] = item["z_mesh_units"] * (
                checkpoint["mesh_length_unit_m"] * 1.0e6)
            item["per_y_widths_nm"] = (
                np.asarray(item["per_y_widths_mesh_units"]) * scale).tolist()
            item["cross_y_span_nm"] = item["cross_y_span_mesh_units"] * scale

    final_reported = float(audit["final_metrics"]["mask_opening_nm"])
    mask_nodal = mask_width["nodal_rows"]["width_nm"]
    mask_dense = mask_width["dense_bilinear_scan"]["width_nm"]
    combined_dense = combined_width["dense_bilinear_scan"]["width_nm"]
    paper_width = mask_dense if connectivity["exterior_to_mask_exit_open"] else 0.0
    report_match_tolerance_nm = max(1.0e-6, 0.01 * dx_nm)
    flicker = timeline["maximum_open_state_one_step_flicker"]
    gates = {
        "checkpoint_connectivity_matches_timeline": (
            timeline["final_connectivity_open"] is None
            or bool(timeline["final_connectivity_open"])
            == connectivity["exterior_to_mask_exit_open"]
        ),
        "reported_final_matches_mask_layer_nodal_metric": (
            abs(final_reported - mask_nodal) <= report_match_tolerance_nm
        ),
        "reported_final_respects_connectivity_qualification": (
            abs(final_reported - (
                mask_nodal
                if connectivity["exterior_to_mask_exit_open"] else 0.0
            )) <= report_match_tolerance_nm
        ),
        "dense_reconstruction_change_below_one_cell": (
            abs(mask_dense - mask_nodal) <= dx_nm
        ),
        "maximum_open_state_one_step_flicker_below_one_cell": (
            flicker is None or flicker["absolute_residual_nm"] <= dx_nm
        ),
    }
    definition_certified = bool(all(gates.values()))
    paper_provenance = dict(PAPER)
    if paper_pdf is not None:
        paper_pdf = Path(paper_pdf)
        paper_provenance["local_pdf_name"] = paper_pdf.name
        paper_provenance["local_pdf_sha256"] = _digest(paper_pdf)

    return {
        "schema": "petch.krueger_2024.mask_opening_certification.v1",
        "implementation": {
            "name": Path(__file__).name,
            "sha256": _digest(Path(__file__)),
        },
        "scope": (
            "base development checkpoint only; no held-out profile was loaded; "
            "no profile evolution was run"
        ),
        "paper_definition": paper_provenance,
        "inputs": {
            "audit_name": audit_path.name,
            "audit_sha256": _digest(audit_path),
            "checkpoint_name": checkpoint_path.name,
            "checkpoint_sha256": _digest(checkpoint_path),
            "mask_material_id": int(mask_material_id),
            "substrate_top_um": float(substrate_top_um),
            "opening_center_um": float(opening_center_um),
            "opening_width_um": float(opening_width_um),
            "dense_samples_per_cell": int(dense_samples_per_cell),
        },
        "mesh": {
            "shape": list(checkpoint["phi"].shape),
            "dx_nm": float(dx_nm),
            "physical_time_s": float(
                checkpoint["metadata"].get("physical_time_s", np.nan)),
            "mask_top_um": float(mask_top),
        },
        "connectivity": connectivity,
        "mask_material_levelset_measurement": mask_width,
        "combined_union_contrast_measurement": combined_width,
        "timeline": timeline,
        "paper_connectivity_qualified_opening_nm": float(paper_width),
        "reported_final_mask_opening_nm": final_reported,
        "combined_union_minus_mask_layer_dense_nm": float(
            combined_dense - mask_dense),
        "gates": gates,
        "definition_certified": definition_certified,
        "grid_authority_certified": False,
        "authority_status": (
            "definition-certified; same-operator 10/5 nm refinement still required"
            if definition_certified else
            "definition certification failed; do not use opening for calibration"
        ),
        "interpretation": [
            "Mask-layer geometry, not the combined substrate-plus-mask union, defines w_m.",
            "An enclosed gas pocket is assigned w_m = 0 exactly; its pocket width "
            "is diagnostic only.",
            "Dense reconstruction and flicker diagnostics change no geometry and smooth no result.",
            "Cross-y span and one-step flicker remain reported numerical evidence.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mask-material-id", type=int, default=2)
    parser.add_argument("--substrate-top-um", type=float)
    parser.add_argument("--opening-center-um", type=float)
    parser.add_argument("--opening-width-um", type=float)
    parser.add_argument("--dense-samples-per-cell", type=int, default=64)
    parser.add_argument("--paper-pdf", type=Path)
    args = parser.parse_args()
    report = build_report(
        args.audit, args.checkpoint,
        mask_material_id=args.mask_material_id,
        substrate_top_um=args.substrate_top_um,
        opening_center_um=args.opening_center_um,
        opening_width_um=args.opening_width_um,
        dense_samples_per_cell=args.dense_samples_per_cell,
        paper_pdf=args.paper_pdf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "definition_certified": report["definition_certified"],
        "grid_authority_certified": report["grid_authority_certified"],
        "paper_connectivity_qualified_opening_nm": (
            report["paper_connectivity_qualified_opening_nm"]),
        "authority_status": report["authority_status"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
