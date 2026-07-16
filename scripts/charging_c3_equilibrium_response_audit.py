#!/usr/bin/env python3
"""Test current-balance accessibility along a physical C3 face-charge mode.

The target patch is selected once from an independent ensemble-mean current audit.  Every voltage
shift in the sweep is then produced by adding charge to that patch in authoritative face-charge
space, projecting it with the engine's compatible Q1 map, and solving the unchanged Poisson system.
No local voltage is edited and no visibility or response physics is softened.  Independent exact
zero-step transport scrambles score each perturbed state.

A confidence-bracketed current crossing proves that the controlling patch can balance under the
declared physics along this admissible charge mode.  Failure to cross along one mode is explicitly
not treated as proof that no equilibrium exists.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_coevolution_c3_trench import (  # noqa: E402
    _geometry_and_poisson, atomic_npz,
)
from petch.charging_coevolution_3d import _patch_balances  # noqa: E402
from petch.charging_coupled_3d import current_balance_metrics_3d  # noqa: E402
from petch.charging_poisson_3d import lump_triangle_sheet_charge_3d  # noqa: E402


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _atomic_json(path: Path, payload) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_json_value(payload), indent=2) + "\n")
    temporary.replace(path)


def _mean_ci(values, confidence=0.95):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1))
    half = float(student_t.ppf(0.5 + confidence / 2.0, len(values) - 1))
    half *= sd / math.sqrt(len(values))
    return dict(mean=mean, sd=sd, ci_half_width=half, lower=mean - half, upper=mean + half)


def _face_potential(potential, points, spacing):
    """Trilinearly sample a nodal field at mesh-coordinate points."""
    potential = np.asarray(potential, dtype=float)
    coordinate = np.asarray(points, dtype=float) / float(spacing)
    lower = np.floor(coordinate).astype(int)
    lower = np.minimum(np.maximum(lower, 0), np.asarray(potential.shape) - 2)
    fraction = coordinate - lower
    value = np.zeros(len(coordinate), dtype=float)
    for x_bit in (0, 1):
        for y_bit in (0, 1):
            for z_bit in (0, 1):
                weight = (
                    (fraction[:, 0] if x_bit else 1.0 - fraction[:, 0])
                    * (fraction[:, 1] if y_bit else 1.0 - fraction[:, 1])
                    * (fraction[:, 2] if z_bit else 1.0 - fraction[:, 2]))
                value += weight * potential[
                    lower[:, 0] + x_bit,
                    lower[:, 1] + y_bit,
                    lower[:, 2] + z_bit]
    return value


def _runner_command(args, checkpoint, seed, output):
    command = [
        sys.executable, str(ROOT / "scripts/charging_coevolution_c3_trench.py"),
        "--output-dir", str(output),
        "--initial-face-state", str(checkpoint),
        "--method-map", str(args.method_map),
        "--method-key", args.method_key,
        "--maximum-steps", "0",
        "--timestep-s", "1.25e-7",
        "--terminal-window-s", "5e-5",
        "--forward-level", str(args.level),
        "--adjoint-level", str(args.level - 2),
        "--electron-estimator", "forward",
        "--n-position", "256",
        "--seed", str(seed),
        "--scramble-mode", "frozen",
        "--trajectory-dt", "0.000078125",
        "--trajectory-max-steps", "16384000",
        "--trajectory-adaptive-horizon",
        "--trajectory-emergency-max-steps", "32768000",
        "--transport-device", args.transport_device,
        "--response-max-bounces", "512",
        "--response-adaptive-bounce-extension",
        "--response-emergency-max-bounces", "1024",
        "--response-tail-tolerance", "1e-10",
        "--response-launch-offset", "5e-6",
    ]
    if args.compatible_q1_charge_state:
        command.append("--compatible-q1-charge-state")
    return command


def _compact_patch_currents(positive_density, negative_density, area, group, label):
    labels = np.unique(group[group >= 0])
    match = np.flatnonzero(labels == label)
    if len(match) != 1:
        raise RuntimeError("fixed target patch disappeared from the unchanged geometry")
    metrics = current_balance_metrics_3d(
        positive_density * area, negative_density * area, group=group)
    index = int(match[0])
    return (
        float(metrics.positive_current_a[index]),
        float(metrics.negative_current_a[index]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--method-map", type=Path, required=True)
    parser.add_argument("--method-key", default="refined_method_hint_Ar+")
    parser.add_argument("--sampling-campaign", type=Path, required=True)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--level", type=int, default=13)
    parser.add_argument(
        "--seeds", type=int, nargs="+",
        default=(10079, 10179, 10279, 10379, 10479, 10579, 10679, 10779))
    parser.add_argument("--target-patch-scale-um", type=float, default=0.5)
    parser.add_argument(
        "--target-potential-shifts-v", type=float, nargs="+",
        default=(-40.0, -30.0, -20.0, -15.0, -10.0, -5.0, 0.0, 10.0))
    parser.add_argument("--transport-device", choices=("cpu", "cuda", "cuda:0"), default="cpu")
    args = parser.parse_args()
    args.checkpoint = args.checkpoint.resolve()
    args.method_map = args.method_map.resolve()
    args.sampling_campaign = args.sampling_campaign.resolve()
    args.campaign_dir = args.campaign_dir.resolve()
    if len(args.seeds) < 2 or args.level < 3:
        parser.error("at least two seeds and a forward sample level of at least three are required")
    if len(set(args.seeds)) != len(args.seeds):
        parser.error("seeds must be unique")
    shifts = np.unique(np.asarray(args.target_potential_shifts_v, dtype=float))
    if len(shifts) < 3 or not np.any(shifts == 0.0) or np.any(~np.isfinite(shifts)):
        parser.error("the potential sweep must contain zero and at least two other finite values")
    args.campaign_dir.mkdir(parents=True, exist_ok=True)

    with np.load(args.checkpoint) as archived:
        sigma_base = np.asarray(archived["sigma_c_per_m2"], dtype=float)
        vertices = np.asarray(archived["vertices"])
        faces = np.asarray(archived["faces"])
        centroids = np.asarray(archived["centroids"], dtype=float)
        areas = np.asarray(archived["areas"], dtype=float)
        face_material = np.asarray(archived["face_material_id"], dtype=int)
        method_hint = np.asarray(archived["method_hint_Ar"])
        checkpoint_potential = np.asarray(archived["potential_v"], dtype=float)
        sampling_metadata = {
            key: np.asarray(archived[key]).copy() for key in (
                "resume_sampling_epoch", "scramble_mode", "scramble_base_seed",
                "sampling_seed_stride", "compatible_q1_charge_state",
                "poisson_periodic_axes", "poisson_independent_node_shape")
            if key in archived}
    args.compatible_q1_charge_state = bool(np.asarray(
        sampling_metadata.get("compatible_q1_charge_state", False)).item())

    sample_audits = []
    for seed in args.seeds:
        path = args.sampling_campaign / f"level{args.level}_seed{seed}" / "current_audit.npz"
        if not path.exists():
            raise RuntimeError(f"missing independent sampling audit: {path}")
        sample_audits.append(path)
    with np.load(sample_audits[0]) as current:
        scales = np.asarray(current["patch_scales_m"], dtype=float)
        groups = np.asarray(current["patch_group_by_scale"], dtype=int)
        physical_area = np.asarray(current["physical_face_area_m2"], dtype=float)
    scale_target = args.target_patch_scale_um * 1e-6
    scale_matches = np.flatnonzero(np.isclose(scales, scale_target, rtol=0.0, atol=1e-15))
    if len(scale_matches) != 1:
        raise RuntimeError("requested physical patch scale is absent or ambiguous")
    scale_index = int(scale_matches[0])
    positive_base = []
    negative_base = []
    for path in sample_audits:
        with np.load(path) as current:
            positive_base.append(np.asarray(
                current["positive_face_current_density_a_m2"], dtype=float))
            negative_base.append(np.asarray(
                current["negative_face_current_density_a_m2"], dtype=float))
    positive_mean = np.mean(positive_base, axis=0)
    negative_mean = np.mean(negative_base, axis=0)
    metrics = current_balance_metrics_3d(
        positive_mean * physical_area, negative_mean * physical_area,
        group=groups[scale_index])
    ratio = np.divide(
        np.abs(metrics.positive_current_a - metrics.negative_current_a),
        metrics.positive_current_a,
        out=np.full(metrics.positive_current_a.shape, np.inf),
        where=metrics.positive_current_a > 0.0)
    compact_target = int(np.argmax(np.where(metrics.active, ratio, -np.inf)))
    patch_labels = np.unique(groups[scale_index][groups[scale_index] >= 0])
    target_label = int(patch_labels[compact_target])
    target_faces = groups[scale_index] == target_label

    geometry, poisson = _geometry_and_poisson(0.25)
    coordinate_unit = geometry.mesh_length_unit_m
    base_charge = lump_triangle_sheet_charge_3d(
        poisson.shape, vertices, faces, sigma_base,
        grid_origin=(0.0, 0.0, 0.0), grid_spacing=geometry.dx,
        coordinate_length_unit_m=coordinate_unit)
    base_potential, _ = poisson.solve(base_charge)
    archive_relative_error = float(
        np.linalg.norm(base_potential - checkpoint_potential)
        / max(np.linalg.norm(checkpoint_potential), 1e-300))
    if archive_relative_error > 1e-12:
        raise RuntimeError("checkpoint potential is inconsistent with its authoritative face charge")
    unit_direction = np.zeros_like(sigma_base)
    unit_direction[target_faces] = 1.0
    response_charge = lump_triangle_sheet_charge_3d(
        poisson.shape, vertices, faces, unit_direction,
        grid_origin=(0.0, 0.0, 0.0), grid_spacing=geometry.dx,
        coordinate_length_unit_m=coordinate_unit)
    response_potential, _ = poisson.solve(response_charge)
    voltage_gain = float(np.mean(_face_potential(
        response_potential, centroids[target_faces], geometry.dx)))
    if not math.isfinite(voltage_gain) or abs(voltage_gain) < 1e-12:
        raise RuntimeError("selected face-charge mode has no resolvable target-patch voltage response")
    base_target_potential = float(np.mean(_face_potential(
        base_potential, centroids[target_faces], geometry.dx)))

    perturbation_records = []
    for shift in shifts:
        tag = f"dv_{shift:+09.3f}".replace("+", "p").replace("-", "m").replace(".", "p")
        perturbation_dir = args.campaign_dir / tag
        perturbation_dir.mkdir(parents=True, exist_ok=True)
        checkpoint = perturbation_dir / "face_checkpoint.npz"
        delta_sigma = float(shift / voltage_gain)
        sigma = sigma_base + delta_sigma * unit_direction
        charge = lump_triangle_sheet_charge_3d(
            poisson.shape, vertices, faces, sigma,
            grid_origin=(0.0, 0.0, 0.0), grid_spacing=geometry.dx,
            coordinate_length_unit_m=coordinate_unit)
        potential, _ = poisson.solve(charge)
        actual_target = float(np.mean(_face_potential(
            potential, centroids[target_faces], geometry.dx)))
        if not checkpoint.exists():
            atomic_npz(
                checkpoint,
                sigma_c_per_m2=sigma,
                face_charge_c=sigma * areas * coordinate_unit ** 2,
                charge_node_c=charge, potential_v=potential,
                vertices=vertices, faces=faces, centroids=centroids, areas=areas,
                face_material_id=face_material, method_hint_Ar=method_hint,
                **sampling_metadata)
        checkpoint_hash = _hash(checkpoint)
        run_paths = []
        for seed in args.seeds:
            output = perturbation_dir / f"level{args.level}_seed{seed}"
            summary_path = output / "summary.json"
            current_path = output / "current_audit.npz"
            if summary_path.exists() and current_path.exists():
                summary = json.loads(summary_path.read_text())
                config = summary["config"]
                compatible = (
                    config.get("initial_face_state_sha256") == checkpoint_hash
                    and config.get("forward_level") == args.level
                    and config.get("adjoint_level") == args.level - 2
                    and config.get("seed") == seed
                    and config.get("maximum_steps") == 0)
                if not compatible:
                    raise RuntimeError(f"incompatible existing response run: {output}")
            else:
                output.mkdir(parents=True, exist_ok=True)
                subprocess.run(_runner_command(args, checkpoint, seed, output), check=True)
            run_paths.append(current_path)

        positive_stack = []
        negative_stack = []
        positive_target = []
        negative_target = []
        for path in run_paths:
            with np.load(path) as current:
                positive = np.asarray(current["positive_face_current_density_a_m2"], dtype=float)
                negative = np.asarray(current["negative_face_current_density_a_m2"], dtype=float)
            positive_stack.append(positive); negative_stack.append(negative)
            target_positive, target_negative = _compact_patch_currents(
                positive, negative, physical_area, groups[scale_index], target_label)
            positive_target.append(target_positive); negative_target.append(target_negative)
        positive_stack = np.stack(positive_stack)
        negative_stack = np.stack(negative_stack)
        positive_ensemble = np.mean(positive_stack, axis=0)
        negative_ensemble = np.mean(negative_stack, axis=0)
        aggregate = _patch_balances(
            positive_ensemble, negative_ensemble, physical_area,
            tuple(groups[index] for index in range(len(scales))), tuple(scales))
        positive_target = np.asarray(positive_target)
        negative_target = np.asarray(negative_target)
        net_target = positive_target - negative_target
        positive_ci = _mean_ci(positive_target)
        negative_ci = _mean_ci(negative_target)
        net_ci = _mean_ci(net_target)
        target_b2 = abs(positive_ci["mean"] - negative_ci["mean"]) / positive_ci["mean"]
        perturbation_records.append(dict(
            requested_target_potential_shift_v=float(shift),
            actual_target_potential_v=actual_target,
            actual_target_potential_shift_v=actual_target - base_target_potential,
            delta_sigma_c_per_m2=delta_sigma,
            added_charge_c=float(delta_sigma * np.sum(physical_area[target_faces])),
            checkpoint_sha256=checkpoint_hash,
            target_positive_current_a=positive_ci,
            target_negative_current_a=negative_ci,
            target_signed_net_current_a=net_ci,
            target_b2=float(target_b2),
            all_patch_b2=[item.b2_maximum_ion_normalized_imbalance for item in aggregate],
            current_audit_sha256_by_seed={
                str(seed): _hash(path) for seed, path in zip(args.seeds, run_paths)}))

    crossings = []
    ordered = sorted(perturbation_records, key=lambda item: item["actual_target_potential_v"])
    for left, right in zip(ordered[:-1], ordered[1:]):
        left_net = left["target_signed_net_current_a"]
        right_net = right["target_signed_net_current_a"]
        opposite_mean = left_net["mean"] * right_net["mean"] <= 0.0
        confidence_bracketed = (
            (left_net["upper"] < 0.0 and right_net["lower"] > 0.0)
            or (left_net["lower"] > 0.0 and right_net["upper"] < 0.0))
        if opposite_mean:
            denominator = right_net["mean"] - left_net["mean"]
            crossing_voltage = (
                left["actual_target_potential_v"]
                if denominator == 0.0 else
                left["actual_target_potential_v"]
                - left_net["mean"] * (
                    right["actual_target_potential_v"]
                    - left["actual_target_potential_v"]) / denominator)
            crossings.append(dict(
                left_shift_v=left["actual_target_potential_shift_v"],
                right_shift_v=right["actual_target_potential_shift_v"],
                interpolated_target_potential_v=float(crossing_voltage),
                confidence_bracketed=bool(confidence_bracketed)))

    confidence_brackets = []
    for left_index, left in enumerate(ordered[:-1]):
        left_net = left["target_signed_net_current_a"]
        left_sign = 1 if left_net["lower"] > 0.0 else (-1 if left_net["upper"] < 0.0 else 0)
        if left_sign == 0:
            continue
        for right in ordered[left_index + 1:]:
            right_net = right["target_signed_net_current_a"]
            right_sign = (
                1 if right_net["lower"] > 0.0 else
                (-1 if right_net["upper"] < 0.0 else 0))
            if right_sign == -left_sign:
                confidence_brackets.append(dict(
                    left_shift_v=left["actual_target_potential_shift_v"],
                    right_shift_v=right["actual_target_potential_shift_v"],
                    voltage_span_v=(
                        right["actual_target_potential_v"]
                        - left["actual_target_potential_v"])))
    confidence_brackets.sort(key=lambda item: item["voltage_span_v"])
    direct_gate_hits = [
        dict(
            target_potential_shift_v=item["actual_target_potential_shift_v"],
            target_potential_v=item["actual_target_potential_v"],
            target_b2=item["target_b2"])
        for item in ordered if item["target_b2"] <= 0.08]
    conclusion = (
        "confidence-bracketed current crossing and a directly scored B2 gate hit exist along "
        "the physical face-charge mode; a missing-physics wall is rejected for target-patch "
        "balance accessibility"
        if confidence_brackets and direct_gate_hits else
        "current crossing evidence exists along the physical face-charge mode, but either its "
        "confidence bracket or a directly scored B2 gate hit remains unresolved"
        if confidence_brackets or direct_gate_hits else
        "no confidence-bracketed crossing was found along this one physical charge mode; "
        "this alone does not prove equilibrium nonexistence")
    audit = dict(
        schema="petch.charging.c3.equilibrium-response-audit.v1",
        status="exact zero-step globally consistent face-charge response sweep",
        engine_git_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        campaign_script_sha256=_hash(Path(__file__).resolve()),
        baseline_checkpoint=dict(name=args.checkpoint.name, sha256=_hash(args.checkpoint)),
        method_map=dict(name=args.method_map.name, sha256=_hash(args.method_map)),
        protocol=dict(
            sample_level=args.level, seeds=args.seeds,
            target_patch_scale_m=float(scales[scale_index]),
            target_patch_group_label=target_label,
            target_face_count=int(np.sum(target_faces)),
            perturbation=(
                "uniform surface-charge-density increment on the fixed controlling patch; "
                "compatible face-to-node projection and exact Poisson solve"),
            archive_potential_relative_reconstruction_error=archive_relative_error,
            target_voltage_gain_v_per_c_m2=voltage_gain,
            base_target_potential_v=base_target_potential,
            exact_operator="hard visibility and declared charged-surface response",
            transport_device=args.transport_device),
        perturbations=perturbation_records,
        mean_crossing_brackets=crossings,
        confidence_brackets=confidence_brackets,
        confidence_bracketed_crossing=bool(confidence_brackets),
        direct_target_b2_gate_hits=direct_gate_hits,
        conclusion=conclusion)
    _atomic_json(args.campaign_dir / "audit.json", audit)

    potential = np.asarray([item["actual_target_potential_v"] for item in ordered])
    ion = np.asarray([item["target_positive_current_a"]["mean"] for item in ordered])
    ion_error = np.asarray([item["target_positive_current_a"]["ci_half_width"] for item in ordered])
    electron = np.asarray([item["target_negative_current_a"]["mean"] for item in ordered])
    electron_error = np.asarray([
        item["target_negative_current_a"]["ci_half_width"] for item in ordered])
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    axes[0].errorbar(potential, ion * 1e12, yerr=ion_error * 1e12, marker="o", label="ion")
    axes[0].errorbar(
        potential, electron * 1e12, yerr=electron_error * 1e12,
        marker="o", label="electron")
    axes[0].set_xlabel("ensemble target-patch potential (V)")
    axes[0].set_ylabel("patch current (pA)")
    axes[0].set_title("Physical charge-mode current crossing")
    axes[0].grid(True, alpha=0.25); axes[0].legend()
    axes[1].plot(
        [item["actual_target_potential_shift_v"] for item in ordered],
        [item["target_b2"] for item in ordered], marker="o")
    axes[1].axhline(0.08, color="black", linestyle="--", linewidth=1.2, label="B2 gate")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("target-patch potential shift (V)")
    axes[1].set_ylabel("fixed target-patch B2")
    axes[1].set_title("Balance along admissible face charge")
    axes[1].grid(True, which="both", alpha=0.25); axes[1].legend()
    figure.tight_layout()
    figure.savefig(args.campaign_dir / "equilibrium_response.png", dpi=180)
    plt.close(figure)
    print(json.dumps(_json_value(audit), indent=2))


if __name__ == "__main__":
    main()
