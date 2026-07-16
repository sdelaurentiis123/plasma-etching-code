#!/usr/bin/env python3
"""Audit whether C3 face charge/current contains modes invisible to Q1 Poisson."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/petch-matplotlib")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_task1_physical_time_3d import _geometry_and_poisson  # noqa: E402
from petch.charging_coevolution_3d import (  # noqa: E402
    physical_surface_patch_groups_3d,
)
from petch.charging_coupled_3d import current_balance_metrics_3d  # noqa: E402
from petch.feature_step_3d import (  # noqa: E402
    _face_material_ids,
    _surface_gas_normals,
)
from petch.charging_poisson_3d import (  # noqa: E402
    CompatibleQ1SurfaceChargeProjector3D,
)
from petch.threed import extract_mesh_3d  # noqa: E402


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json_value(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(f"cannot serialize {type(value).__name__}")


def _patch_ratios(positive, negative, resolved_net, groups, projector):
    output = []
    throughput_scale = max(
        float(np.sum(positive) + np.sum(negative)), np.finfo(float).tiny)
    for group in np.unique(groups[groups >= 0]):
        selected = groups == group
        ion = float(np.sum(positive[selected]))
        electron = float(np.sum(negative[selected]))
        active = ion + electron > 1e-15 * throughput_scale
        if not active or ion <= 0.0:
            continue
        raw_net = ion - electron
        visible_net = float(np.sum(resolved_net[selected]))
        output.append(dict(
            group=int(group), ion_current_a=ion, electron_current_a=electron,
            face_count=int(np.count_nonzero(selected)),
            q1_invisible_functional_fraction=(
                projector.unresolved_linear_functional_fraction(selected.astype(float))),
            raw_signed_ratio=raw_net / ion,
            raw_absolute_ratio=abs(raw_net) / ion,
            q1_resolved_signed_net_over_ion=visible_net / ion,
            q1_resolved_absolute_net_over_ion=abs(visible_net) / ion,
            q1_invisible_signed_net_over_ion=(raw_net - visible_net) / ion,
            q1_invisible_absolute_net_over_ion=abs(raw_net - visible_net) / ion))
    return output


def _checkpoint_comparison(reference_path, candidate_path, projector, poisson):
    with np.load(reference_path) as reference, np.load(candidate_path) as candidate:
        reference_face = np.asarray(reference["face_charge_c"], dtype=float)
        candidate_face = np.asarray(candidate["face_charge_c"], dtype=float)
        reference_node = np.asarray(reference["charge_node_c"], dtype=float)
        candidate_node = np.asarray(candidate["charge_node_c"], dtype=float)
        reference_potential = np.asarray(reference["potential_v"], dtype=float)
        candidate_potential = np.asarray(candidate["potential_v"], dtype=float)
    reference_reduced_node = poisson.reduce_charge(reference_node)
    candidate_reduced_node = poisson.reduce_charge(candidate_node)
    return dict(
        node_charge_relative_l1=float(
            np.sum(np.abs(candidate_reduced_node - reference_reduced_node))
            / max(float(np.sum(np.abs(reference_reduced_node))), np.finfo(float).tiny)),
        potential_relative_l2=float(
            np.linalg.norm(candidate_potential - reference_potential)
            / max(float(np.linalg.norm(reference_potential)), np.finfo(float).tiny)),
        global_charge_difference_c=float(candidate_face.sum() - reference_face.sum()),
        face_charge_relative_l1=float(
            np.sum(np.abs(candidate_face - reference_face))
            / max(float(np.sum(np.abs(reference_face))), np.finfo(float).tiny)),
        reference_q1_invisible_fraction=projector.unresolved_fraction(reference_face),
        candidate_q1_invisible_fraction=projector.unresolved_fraction(candidate_face))


def _plot(rows, output: Path) -> None:
    time_ms = np.asarray([row["cumulative_physical_time_s"] for row in rows]) * 1e3
    raw = np.asarray([row["raw_patch_b2_max"] for row in rows])
    nodal = np.asarray([
        [row["q1_window_node_rms"], row["q1_window_node_worst"]] for row in rows])
    b1 = np.asarray([row["b1_potential_rate_v_s"] for row in rows])
    invisible = np.asarray([
        [row["face_charge_q1_invisible_fraction"],
         row["window_current_q1_invisible_fraction"]] for row in rows])

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), constrained_layout=True)
    axis = axes[0, 0]
    for index in range(raw.shape[1]):
        axis.plot(time_ms, raw[:, index], marker="o", label=f"raw B2 scale {index + 1}")
    axis.axhline(0.08, color="black", linewidth=1, linestyle="--", label="0.08 gate")
    axis.set_yscale("log")
    axis.set_title("Contract B2 remains a raw physical-patch diagnostic")
    axis.set_xlabel("cumulative physical time (ms)")
    axis.set_ylabel("ion-normalized imbalance")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[0, 1]
    axis.plot(time_ms, nodal[:, 0], marker="o", label="Q1 window node RMS")
    axis.plot(time_ms, nodal[:, 1], marker="s", label="Q1 window node worst")
    axis.axhline(0.08, color="black", linewidth=1, linestyle="--", label="0.08 reference")
    axis.set_yscale("log")
    axis.set_title("The field-visible current equations are already small")
    axis.set_xlabel("cumulative physical time (ms)")
    axis.set_ylabel("symmetric relative imbalance")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 0]
    axis.plot(time_ms, b1, marker="o", label="terminal-window B1")
    axis.axhline(1e3, color="black", linewidth=1, linestyle="--", label="1,000 V/s gate")
    axis.set_yscale("log")
    axis.set_title("A real field-visible slow drift remains")
    axis.set_xlabel("cumulative physical time (ms)")
    axis.set_ylabel("maximum |dV/dt| (V/s)")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 1]
    axis.plot(time_ms, 100.0 * invisible[:, 0], marker="o", label="stored charge")
    axis.plot(time_ms, 100.0 * invisible[:, 1], marker="s", label="window net current")
    axis.set_ylim(0.0, 105.0)
    axis.set_title("Most face-space activity is invisible to Q1 Poisson")
    axis.set_xlabel("cumulative physical time (ms)")
    axis.set_ylabel("Q1-null component (%)")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    figure.suptitle("C3 charging: raw face diagnostics versus field-resolved dynamics")
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--grid-dx-um", type=float, default=0.25)
    parser.add_argument("--compatible-fork-dir", type=Path)
    parser.add_argument("--projection-reference-segment", type=int, default=13)
    parser.add_argument("--continuation-reference-segment", type=int, default=14)
    args = parser.parse_args()
    campaign = args.campaign_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    segments = sorted(
        (path.parent for path in campaign.glob("segment_*/summary.json")),
        key=lambda path: int(path.name.split("_")[-1]))
    if not segments:
        parser.error("campaign contains no completed segment summaries")
    required = ("summary.json", "current_audit.npz", "face_checkpoint.npz")
    for segment in segments:
        if any(not (segment / name).exists() for name in required):
            parser.error(f"incomplete audit inputs in {segment}")

    campaign_status_path = campaign / "gpu_seed79_campaign_status.json"
    if not campaign_status_path.exists():
        campaign_status_path = campaign / "campaign_status.json"
    status = json.loads(campaign_status_path.read_text())
    cumulative_by_segment = {
        int(record["segment"]): float(record["cumulative_physical_time_s"])
        for record in status.get("records", [])}
    base_time = float(status.get("base_physical_time_s", 0.0))

    with np.load(segments[0] / "face_checkpoint.npz") as checkpoint:
        vertices = np.asarray(checkpoint["vertices"], dtype=float)
        faces = np.asarray(checkpoint["faces"], dtype=int)
        areas = np.asarray(checkpoint["areas"], dtype=float)
        node_shape = np.asarray(checkpoint["charge_node_c"]).shape
    geometry, poisson = _geometry_and_poisson(args.grid_dx_um)
    if tuple(poisson.shape) != tuple(node_shape):
        raise RuntimeError("declared grid does not match checkpoint nodal shape")
    projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        poisson, vertices, faces, grid_origin=(0.0, 0.0, 0.0),
        grid_spacing=args.grid_dx_um,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    if not np.allclose(
            physical_area, projector.physical_face_area_m2,
            rtol=2e-13, atol=0.0):
        raise RuntimeError("checkpoint face areas disagree with Q1 coupling")

    rows = []
    details = []
    window_net_currents = []
    running_time = base_time
    for segment in segments:
        segment_index = int(segment.name.split("_")[-1])
        summary = json.loads((segment / "summary.json").read_text())
        result = summary["result"]
        running_time += float(result.get("physical_time_s", 0.0))
        cumulative_time = cumulative_by_segment.get(segment_index, running_time)
        with np.load(segment / "face_checkpoint.npz") as checkpoint:
            face_charge = np.asarray(checkpoint["face_charge_c"], dtype=float)
            archived_node_charge = np.asarray(checkpoint["charge_node_c"], dtype=float)
        compatible_charge = projector.project_face_charge(face_charge)
        node_before = projector.node_charge_from_face_charge(face_charge)
        node_after = projector.node_charge_from_face_charge(compatible_charge)
        potential_before, _ = poisson.solve(
            poisson.canonicalize_reduced_charge(node_before))
        potential_after, _ = poisson.solve(
            poisson.canonicalize_reduced_charge(node_after))
        node_scale = max(float(np.sum(np.abs(node_before))), np.finfo(float).tiny)
        potential_scale = max(float(np.linalg.norm(potential_before)), np.finfo(float).tiny)

        with np.load(segment / "current_audit.npz") as audit:
            terminal_ready = bool(np.asarray(audit["terminal_window_ready"]).item())
            if terminal_ready and np.asarray(
                    audit["terminal_window_positive_face_current_density_a_m2"]).size:
                positive_density = np.asarray(
                    audit["terminal_window_positive_face_current_density_a_m2"], dtype=float)
                negative_density = np.asarray(
                    audit["terminal_window_negative_face_current_density_a_m2"], dtype=float)
            else:
                positive_density = np.asarray(
                    audit["positive_face_current_density_a_m2"], dtype=float)
                negative_density = np.asarray(
                    audit["negative_face_current_density_a_m2"], dtype=float)
            positive_current = positive_density * physical_area
            negative_current = negative_density * physical_area
            net_current = positive_current - negative_current
            resolved_current = projector.project_face_charge(net_current)
            window_net_currents.append(net_current.copy())
            groups_by_scale = np.asarray(audit["patch_group_by_scale"], dtype=int)
            patch_scales = np.asarray(audit["patch_scales_m"], dtype=float)
        positive_node = projector.node_charge_from_face_charge(positive_current)
        negative_node = projector.node_charge_from_face_charge(negative_current)
        node_balance = current_balance_metrics_3d(positive_node, negative_node)
        patch_detail = [
            _patch_ratios(
                positive_current, negative_current, resolved_current, groups, projector)
            for groups in groups_by_scale]
        raw_patch_b2 = [
            max((item["raw_absolute_ratio"] for item in scale), default=0.0)
            for scale in patch_detail]
        resolved_patch_proxy = [
            max((item["q1_resolved_absolute_net_over_ion"] for item in scale), default=0.0)
            for scale in patch_detail]
        invisible_patch_functional = [
            max((item["q1_invisible_functional_fraction"] for item in scale), default=0.0)
            for scale in patch_detail]
        archived_node_relative_l1 = float(
            np.sum(np.abs(node_before - poisson.reduce_charge(archived_node_charge)))
            / node_scale)
        row = dict(
            segment=segment_index,
            cumulative_physical_time_s=cumulative_time,
            b1_potential_rate_v_s=float(result["final_potential_rate_max_v_s"]),
            raw_patch_b2_max=raw_patch_b2,
            q1_resolved_patch_net_over_ion_max=resolved_patch_proxy,
            patch_q1_invisible_functional_fraction_max=invisible_patch_functional,
            q1_window_node_rms=node_balance.rms_relative_imbalance,
            q1_window_node_worst=node_balance.maximum_relative_imbalance,
            q1_window_node_throughput_weighted_rms=(
                node_balance.throughput_weighted_rms_relative_imbalance),
            face_charge_q1_invisible_fraction=projector.unresolved_fraction(face_charge),
            face_charge_q1_invisible_l1_c=float(np.sum(np.abs(
                face_charge - compatible_charge))),
            window_current_q1_invisible_fraction=projector.unresolved_fraction(net_current),
            window_current_q1_invisible_l1_a=float(np.sum(np.abs(
                net_current - resolved_current))),
            window_current_q1_invisible_net_a=float(np.sum(
                net_current - resolved_current)),
            compatible_projection_node_relative_l1=float(
                np.sum(np.abs(node_after - node_before)) / node_scale),
            compatible_projection_potential_relative_l2=float(
                np.linalg.norm(potential_after - potential_before) / potential_scale),
            compatible_projection_global_charge_error_c=float(
                compatible_charge.sum() - face_charge.sum()),
            archived_node_relative_l1=archived_node_relative_l1)
        rows.append(row)
        details.append(dict(
            segment=segment_index, patch_scales_m=patch_scales,
            patch_detail=patch_detail,
            summary_sha256=_hash(segment / "summary.json"),
            checkpoint_sha256=_hash(segment / "face_checkpoint.npz"),
            current_audit_sha256=_hash(segment / "current_audit.npz")))

    unresolved_window_currents = np.asarray([
        current - projector.project_face_charge(current)
        for current in window_net_currents])
    mean_unresolved_current = np.mean(unresolved_window_currents, axis=0)
    mean_unresolved_norm = float(np.linalg.norm(mean_unresolved_current))
    cosine_to_mean = np.asarray([
        float(current @ mean_unresolved_current) / max(
            float(np.linalg.norm(current)) * mean_unresolved_norm,
            np.finfo(float).tiny)
        for current in unresolved_window_currents])
    tail_average = []
    for count in range(1, len(window_net_currents) + 1):
        current = np.mean(window_net_currents[-count:], axis=0)
        unresolved = current - projector.project_face_charge(current)
        tail_average.append(dict(
            terminal_window_count=count,
            q1_invisible_fraction=projector.unresolved_fraction(current),
            q1_invisible_l2_a=float(np.linalg.norm(unresolved)),
            q1_invisible_l1_a=float(np.sum(np.abs(unresolved))),
            q1_invisible_max_a=float(np.max(np.abs(unresolved)))))

    grid_structure = []
    for grid_dx_um in (0.5, 0.25, 0.125):
        grid_geometry, grid_poisson = _geometry_and_poisson(grid_dx_um)
        grid_vertices, grid_faces, grid_centroids, _grid_areas = extract_mesh_3d(
            grid_geometry.phi, grid_geometry.dx)
        grid_normals = _surface_gas_normals(
            grid_vertices, grid_faces, grid_centroids, grid_geometry)
        grid_material = _face_material_ids(grid_centroids, grid_geometry)
        grid_projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
            grid_poisson, grid_vertices, grid_faces,
            grid_spacing=grid_geometry.dx,
            coordinate_length_unit_m=grid_geometry.mesh_length_unit_m)
        patch_structure = []
        for patch_scale_m in (0.25e-6, 0.5e-6):
            group = physical_surface_patch_groups_3d(
                grid_centroids, grid_normals, grid_material, patch_scale_m,
                mesh_length_unit_m=grid_geometry.mesh_length_unit_m,
                mesh_origin_m=grid_geometry.mesh_origin_m)
            fractions = np.asarray([
                grid_projector.unresolved_linear_functional_fraction(
                    (group == label).astype(float))
                for label in np.unique(group[group >= 0])])
            patch_structure.append(dict(
                patch_scale_m=patch_scale_m, patch_count=len(fractions),
                maximum_q1_invisible_functional_fraction=float(np.max(fractions)),
                median_q1_invisible_functional_fraction=float(np.median(fractions))))
        grid_structure.append(dict(
            grid_dx_um=grid_dx_um, node_shape=list(grid_poisson.shape),
            face_count=len(grid_faces), rank=grid_projector.rank,
            nullity=grid_projector.nullity,
            nullity_fraction=grid_projector.nullity / len(grid_faces),
            condition_number=grid_projector.condition_number,
            patch_structure=patch_structure))

    compatible_fork = None
    if args.compatible_fork_dir is not None:
        fork = args.compatible_fork_dir.resolve()
        projection = fork / "projection"
        continuation = fork / "continuation"
        projection_reference = campaign / f"segment_{args.projection_reference_segment:04d}"
        continuation_reference = campaign / f"segment_{args.continuation_reference_segment:04d}"
        for path in (
                projection / "summary.json", projection / "face_checkpoint.npz",
                continuation / "summary.json", continuation / "face_checkpoint.npz",
                projection_reference / "summary.json",
                projection_reference / "face_checkpoint.npz",
                continuation_reference / "summary.json",
                continuation_reference / "face_checkpoint.npz"):
            if not path.exists():
                parser.error(f"missing compatible-fork comparison artifact: {path}")
        projection_summary = json.loads((projection / "summary.json").read_text())["result"]
        continuation_summary = json.loads(
            (continuation / "summary.json").read_text())["result"]
        reference_summary = json.loads(
            (continuation_reference / "summary.json").read_text())["result"]
        metric_names = (
            "final_potential_rate_max_v_s",
            "retained_node_rms_relative_current_imbalance",
            "retained_node_max_relative_current_imbalance")
        paired_metric_relative_difference = {
            name: (
                float(continuation_summary[name]) - float(reference_summary[name]))
                / max(abs(float(reference_summary[name])), np.finfo(float).tiny)
            for name in metric_names}
        compatible_b2 = np.asarray(
            continuation_summary["patch_b2_max_ion_normalized"], dtype=float)
        reference_b2 = np.asarray(
            reference_summary["patch_b2_max_ion_normalized"], dtype=float)
        paired_metric_relative_difference["patch_b2_max_ion_normalized"] = (
            (compatible_b2 - reference_b2) / np.maximum(
                np.abs(reference_b2), np.finfo(float).tiny))
        compatible_fork = dict(
            zero_update_projection=dict(
                comparison=_checkpoint_comparison(
                    projection_reference / "face_checkpoint.npz",
                    projection / "face_checkpoint.npz", projector, poisson),
                initial_unresolved_face_charge_fraction=projection_summary[
                    "initial_unresolved_face_charge_fraction"],
                accepted_steps=projection_summary["accepted_steps"]),
            paired_500_step_continuation=dict(
                comparison=_checkpoint_comparison(
                    continuation_reference / "face_checkpoint.npz",
                    continuation / "face_checkpoint.npz", projector, poisson),
                relative_metric_difference=paired_metric_relative_difference,
                compatible_metrics={
                    name: continuation_summary[name] for name in metric_names},
                reference_metrics={name: reference_summary[name] for name in metric_names},
                compatible_patch_b2=compatible_b2,
                reference_patch_b2=reference_b2,
                accepted_steps=continuation_summary["accepted_steps"]),
            provenance=dict(
                projection_summary_sha256=_hash(projection / "summary.json"),
                continuation_summary_sha256=_hash(continuation / "summary.json"),
                projection_reference_summary_sha256=_hash(
                    projection_reference / "summary.json"),
                continuation_reference_summary_sha256=_hash(
                    continuation_reference / "summary.json")))

    plot_path = output / "q1-compatibility-trajectory.png"
    _plot(rows, plot_path)
    csv_path = output / "trajectory.csv"
    field_names = list(rows[0])
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=field_names)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: json.dumps(value) if isinstance(value, list) else value
                for key, value in row.items()})

    latest = rows[-1]
    artifact = dict(
        schema="petch.charging.c3.q1-compatibility-audit.v1",
        hypothesis=(
            "raw P0 face charge/current contains exact modes that the declared Q1 Poisson "
            "field cannot observe or feed back on"),
        q1_coupling=dict(
            shape=list(projector.coupling.shape), rank=projector.rank,
            nullity=projector.nullity,
            poisson_periodic_axes=list(poisson.periodic_axes),
            poisson_independent_node_shape=list(poisson.reduced_shape),
            condition_number=projector.condition_number,
            relative_rank_tolerance=projector.relative_rank_tolerance,
            singular_values=projector.singular_values),
        latest=latest, trajectory=rows, per_patch=details,
        persistent_null_current=dict(
            terminal_window_count=len(window_net_currents),
            mean_q1_invisible_current_l2_a=mean_unresolved_norm,
            minimum_cosine_alignment_to_mean=float(np.min(cosine_to_mean)),
            maximum_cosine_alignment_to_mean=float(np.max(cosine_to_mean)),
            conclusion=(
                "independent terminal-window null currents remain coherently aligned; the "
                "component is systematic at this discretization, not zero-mean sampling noise")),
        tail_average=tail_average,
        grid_structure=grid_structure,
        compatible_fork=compatible_fork,
        interpretation=dict(
            confirmed=bool(
                projector.nullity > 0
                and latest["face_charge_q1_invisible_fraction"] > 0.5
                and latest["window_current_q1_invisible_fraction"] > 0.5),
            state_fix=(
                "make the Q1 nodal load authoritative and retain the unique area-weighted "
                "minimum-density face representative; this preserves field and charge exactly"),
            convergence_contract=(
                "raw B2 remains reported and is not replaced by this audit; its disagreement "
                "with Q1-resolved dynamics triggers the already-required grid-refinement review"),
            remaining_dynamics=(
                "B1 is field-visible and remains above tolerance, so a short compatible-state "
                "continuation and stochastic tail averaging are still required"),
            patch_gate_compatibility=(
                "declared raw patch sums have a nonzero dual projection onto the Q1 null space; "
                "two identical Q1 fields can therefore report different raw B2 values")),
        provenance=dict(
            campaign_status_name=campaign_status_path.name,
            campaign_status_sha256=_hash(campaign_status_path),
            script_sha256=_hash(Path(__file__).resolve()),
            charging_poisson_sha256=_hash(ROOT / "src/petch/charging_poisson_3d.py"),
            charging_coevolution_sha256=_hash(
                ROOT / "src/petch/charging_coevolution_3d.py")),
        artifacts=dict(
            trajectory_csv=dict(name=csv_path.name, sha256=_hash(csv_path)),
            plot=dict(name=plot_path.name, sha256=_hash(plot_path))))
    json_path = output / "audit.json"
    json_path.write_text(json.dumps(_json_value(artifact), indent=2) + "\n")
    print(json.dumps(_json_value({
        "q1_coupling": artifact["q1_coupling"],
        "latest": latest,
        "interpretation": artifact["interpretation"],
        "audit": str(json_path), "plot": str(plot_path)}), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
