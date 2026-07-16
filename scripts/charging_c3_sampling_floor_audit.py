#!/usr/bin/env python3
"""Audit whether a C3 checkpoint residual survives charged-transport refinement.

Each scoring run evaluates the exact hard-visibility operator without advancing surface charge.
Independent scrambled-Sobol replicates are paired across nested sample levels.  The campaign
averages signed face currents before computing nonlinear patch maxima; averaging per-replicate B2
values would answer a different and biased question.  The runner's separately certified estimator
method map remains frozen throughout.
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

from charging_task1_physical_time_3d import _geometry_and_poisson  # noqa: E402
from petch.charging_coevolution_3d import (  # noqa: E402
    _patch_balances,
    _q1_patch_balance_diagnostics,
)
from petch.charging_coupled_3d import current_balance_metrics_3d  # noqa: E402
from petch.charging_poisson_3d import (  # noqa: E402
    CompatibleQ1SurfaceChargeProjector3D,
)


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
    array = np.asarray(values, dtype=float)
    mean = float(np.mean(array))
    if len(array) < 2:
        return dict(mean=mean, sd=None, ci_half_width=None)
    sd = float(np.std(array, ddof=1))
    quantile = float(student_t.ppf(0.5 + confidence / 2.0, len(array) - 1))
    return dict(
        mean=mean, sd=sd,
        ci_half_width=quantile * sd / math.sqrt(len(array)))


def _runner_command(args, level, seed, output):
    runner = ROOT / "scripts/charging_coevolution_c3_trench.py"
    command = [
        sys.executable, str(runner),
        "--output-dir", str(output),
        "--initial-face-state", str(args.checkpoint),
        "--method-map", str(args.method_map),
        "--method-key", args.method_key,
        "--maximum-steps", "0",
        "--timestep-s", "1.25e-7",
        "--terminal-window-s", "5e-5",
        "--forward-level", str(level),
        "--adjoint-level", str(level - 2),
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


def _load_run(
        path: Path, checkpoint_hash: str, level: int, seed: int,
        compatible_q1_charge_state: bool):
    summary_path = path / "summary.json"
    current_path = path / "current_audit.npz"
    if not summary_path.exists() or not current_path.exists():
        return None
    summary = json.loads(summary_path.read_text())
    config = summary["config"]
    expected = {
        "maximum_steps": 0,
        "forward_level": level,
        "adjoint_level": level - 2,
        "seed": seed,
        "electron_estimator": "forward",
        "initial_face_state_sha256": checkpoint_hash,
        "compatible_q1_charge_state": bool(compatible_q1_charge_state),
    }
    if any(config.get(key) != value for key, value in expected.items()):
        raise RuntimeError(f"existing run has incompatible provenance: {path}")
    return summary, current_path


def _patch_record(positive_density, negative_density, area, group):
    positive = np.asarray(positive_density) * area
    negative = np.asarray(negative_density) * area
    metrics = current_balance_metrics_3d(positive, negative, group=group)
    ratio = np.divide(
        np.abs(metrics.positive_current_a - metrics.negative_current_a),
        metrics.positive_current_a,
        out=np.full(metrics.positive_current_a.shape, np.inf),
        where=metrics.positive_current_a > 0.0)
    masked = np.where(metrics.active, ratio, -np.inf)
    index = int(np.argmax(masked))
    labels = np.unique(np.asarray(group)[np.asarray(group) >= 0])
    return dict(
        compact_index=index,
        group_label=int(labels[index]),
        b2=float(ratio[index]),
        positive_current_a=float(metrics.positive_current_a[index]),
        negative_current_a=float(metrics.negative_current_a[index]),
        signed_net_current_a=float(
            metrics.positive_current_a[index] - metrics.negative_current_a[index]),
        active_count=int(metrics.active_count))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--method-map", type=Path, required=True)
    parser.add_argument("--method-key", default="refined_method_hint_Ar+")
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--levels", type=int, nargs="+", default=(11, 12, 13))
    parser.add_argument(
        "--seeds", type=int, nargs="+",
        default=(10079, 10179, 10279, 10379, 10479, 10579, 10679, 10779))
    parser.add_argument("--transport-device", choices=("cpu", "cuda", "cuda:0"), default="cpu")
    args = parser.parse_args()
    args.checkpoint = args.checkpoint.resolve()
    args.method_map = args.method_map.resolve()
    args.campaign_dir = args.campaign_dir.resolve()
    if len(args.seeds) < 2 or len(args.levels) < 2:
        parser.error("at least two seeds and two nested sample levels are required")
    if any(level < 3 for level in args.levels):
        parser.error("forward levels must be at least three")
    if len(set(args.seeds)) != len(args.seeds) or len(set(args.levels)) != len(args.levels):
        parser.error("levels and seeds must be unique")
    args.campaign_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_hash = _hash(args.checkpoint)
    method_hash = _hash(args.method_map)
    with np.load(args.checkpoint) as checkpoint:
        args.compatible_q1_charge_state = bool(np.asarray(
            checkpoint["compatible_q1_charge_state"]).item()) if (
                "compatible_q1_charge_state" in checkpoint) else False
        checkpoint_node_charge = np.asarray(checkpoint["charge_node_c"], dtype=float)
        checkpoint_potential = np.asarray(checkpoint["potential_v"], dtype=float)
        checkpoint_vertices = np.asarray(checkpoint["vertices"])
        checkpoint_faces = np.asarray(checkpoint["faces"], dtype=int)
        checkpoint_areas = np.asarray(checkpoint["areas"], dtype=float)
    geometry, poisson = _geometry_and_poisson(0.25)
    projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        poisson, checkpoint_vertices, checkpoint_faces, grid_spacing=geometry.dx,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)

    runs = {}
    for level in sorted(args.levels):
        for seed in args.seeds:
            output = args.campaign_dir / f"level{level}_seed{seed}"
            loaded = _load_run(
                output, checkpoint_hash, level, seed,
                args.compatible_q1_charge_state)
            if loaded is None:
                output.mkdir(parents=True, exist_ok=True)
                subprocess.run(_runner_command(args, level, seed, output), check=True)
                loaded = _load_run(
                    output, checkpoint_hash, level, seed,
                    args.compatible_q1_charge_state)
                if loaded is None:
                    raise RuntimeError(f"runner did not produce complete artifacts: {output}")
            runs[level, seed] = loaded

    level_records = []
    aggregate_b2 = []
    aggregate_q1_b2 = []
    aggregate_potential_rate = []
    first_groups = None
    first_scales = None
    first_area = None
    for level in sorted(args.levels):
        positive_stack = []
        negative_stack = []
        positive_node_stack = []
        negative_node_stack = []
        potential_rate_stack = []
        scalar = []
        for seed in args.seeds:
            summary, current_path = runs[level, seed]
            result = summary["result"]
            with np.load(current_path) as current:
                positive_stack.append(np.asarray(
                    current["positive_face_current_density_a_m2"], dtype=float))
                negative_stack.append(np.asarray(
                    current["negative_face_current_density_a_m2"], dtype=float))
                positive_node_stack.append(np.asarray(
                    current["positive_current_node_a"], dtype=float))
                negative_node_stack.append(np.asarray(
                    current["negative_current_node_a"], dtype=float))
                potential_rate_stack.append(np.asarray(
                    current["potential_rate_v_s"], dtype=float))
                scales = np.asarray(current["patch_scales_m"], dtype=float)
                groups = np.asarray(current["patch_group_by_scale"], dtype=int)
                area = np.asarray(current["physical_face_area_m2"], dtype=float)
            if first_groups is None:
                first_groups = groups.copy(); first_scales = scales.copy(); first_area = area.copy()
            elif (not np.array_equal(groups, first_groups)
                    or not np.array_equal(scales, first_scales)
                    or not np.array_equal(area, first_area)):
                raise RuntimeError("patch maps or face areas changed within the paired audit")
            scalar.append(dict(
                seed=seed,
                node_rms=result["retained_node_rms_relative_current_imbalance"],
                node_worst=result["retained_node_max_relative_current_imbalance"],
                b2=list(result["patch_b2_max_ion_normalized"]),
                potential_rate_max_v_s=result["final_instantaneous_potential_rate_max_v_s"],
                current_audit_sha256=_hash(current_path)))
        positive_stack = np.stack(positive_stack)
        negative_stack = np.stack(negative_stack)
        positive_node_stack = np.stack(positive_node_stack)
        negative_node_stack = np.stack(negative_node_stack)
        potential_rate_stack = np.stack(potential_rate_stack)
        positive_mean = np.mean(positive_stack, axis=0)
        negative_mean = np.mean(negative_stack, axis=0)
        positive_node_mean = np.mean(positive_node_stack, axis=0)
        negative_node_mean = np.mean(negative_node_stack, axis=0)
        aggregate = _patch_balances(
            positive_mean, negative_mean, first_area,
            tuple(first_groups[index] for index in range(len(first_scales))),
            tuple(first_scales))
        b2 = [item.b2_maximum_ion_normalized_imbalance for item in aggregate]
        aggregate_b2.append(b2)
        mean_net_face_a = (positive_mean - negative_mean) * first_area
        resolved_net_face_a = (
            projector.project_face_charge(mean_net_face_a)
            if args.compatible_q1_charge_state else mean_net_face_a)
        sensitivities = tuple(
            (max(projector.unresolved_linear_functional_fraction(
                (group == label).astype(float)) for label in np.unique(group))
             if args.compatible_q1_charge_state else 0.0)
            for group in first_groups)
        q1_patch = _q1_patch_balance_diagnostics(
            positive_mean, negative_mean, first_area,
            tuple(first_groups[index] for index in range(len(first_scales))),
            tuple(first_scales), resolved_net_face_a, sensitivities)
        q1_b2 = [
            item["q1_resolved_maximum_ion_normalized_imbalance"]
            for item in q1_patch]
        aggregate_q1_b2.append(q1_b2)
        node_balance = current_balance_metrics_3d(
            poisson.reduce_charge(positive_node_mean),
            poisson.reduce_charge(negative_node_mean))
        audit_dt = 1.25e-7
        candidate_charge = poisson.canonicalize_charge(
            checkpoint_node_charge
            + (positive_node_mean - negative_node_mean) * audit_dt)
        candidate_potential, _ = poisson.solve(candidate_charge)
        potential_rate_vector = (
            candidate_potential - checkpoint_potential) / audit_dt
        potential_rate = float(np.max(np.abs(potential_rate_vector)))
        aggregate_potential_rate.append(potential_rate)
        replicate_rate_flat = potential_rate_stack.reshape(len(args.seeds), -1)
        mean_rate_flat = potential_rate_vector.reshape(-1)
        rate_standard_error_flat = np.std(
            replicate_rate_flat, axis=0, ddof=1) / math.sqrt(len(args.seeds))
        mean_rate_l2 = float(np.linalg.norm(mean_rate_flat))
        rate_standard_error_l2 = float(np.linalg.norm(rate_standard_error_flat))
        largest_mean_component = int(np.argmax(np.abs(mean_rate_flat)))
        pointwise_quantile = float(student_t.ppf(0.975, len(args.seeds) - 1))
        largest_component_ci_half_width = float(
            pointwise_quantile * rate_standard_error_flat[largest_mean_component])
        controlling = []
        for scale_index, scale in enumerate(first_scales):
            mean_patch = _patch_record(
                positive_mean, negative_mean, first_area, first_groups[scale_index])
            label = mean_patch["group_label"]
            labels = np.unique(first_groups[scale_index][first_groups[scale_index] >= 0])
            compact = int(np.flatnonzero(labels == label)[0])
            positive_values = []
            negative_values = []
            net_values = []
            for positive, negative in zip(positive_stack, negative_stack):
                metrics = current_balance_metrics_3d(
                    positive * first_area, negative * first_area,
                    group=first_groups[scale_index])
                positive_values.append(metrics.positive_current_a[compact])
                negative_values.append(metrics.negative_current_a[compact])
                net_values.append(
                    metrics.positive_current_a[compact] - metrics.negative_current_a[compact])
            mean_patch.update(
                patch_scale_m=float(scale),
                positive_current=_mean_ci(positive_values),
                negative_current=_mean_ci(negative_values),
                signed_net_current=_mean_ci(net_values))
            controlling.append(mean_patch)
        level_records.append(dict(
            forward_level=level, adjoint_level=level - 2,
            forward_samples=2 ** level, adjoint_samples=2 ** (level - 2),
            replicate_count=len(args.seeds),
            aggregate_signed_current_b2=b2,
            aggregate_q1_resolved_signed_current_b2=q1_b2,
            ensemble_mean_node_rms=node_balance.rms_relative_imbalance,
            ensemble_mean_node_worst=node_balance.maximum_relative_imbalance,
            ensemble_mean_potential_rate_max_v_s=potential_rate,
            potential_rate_vector_diagnostic=dict(
                note=(
                    "Between-scramble standard errors are pointwise diagnostics. The largest "
                    "component was selected after seeing the mean and its interval is not a "
                    "simultaneous B1 acceptance interval."),
                mean_vector_l2_v_s=mean_rate_l2,
                standard_error_vector_l2_v_s=rate_standard_error_l2,
                global_signal_to_standard_error=(
                    mean_rate_l2 / rate_standard_error_l2
                    if rate_standard_error_l2 > 0.0 else math.inf),
                largest_mean_component_flat_index=largest_mean_component,
                largest_mean_component_v_s=float(
                    mean_rate_flat[largest_mean_component]),
                pointwise_95_ci_half_width_v_s=largest_component_ci_half_width,
                pointwise_interval_contains_zero=bool(
                    abs(mean_rate_flat[largest_mean_component])
                    <= largest_component_ci_half_width),
                mean_vector_v_s=potential_rate_vector,
                standard_error_vector_v_s=rate_standard_error_flat.reshape(
                    potential_rate_vector.shape)),
            controlling_patches=controlling,
            per_replicate=scalar))

    paired = []
    ordered = sorted(args.levels)
    for lower, upper in zip(ordered[:-1], ordered[1:]):
        low = level_records[ordered.index(lower)]
        high = level_records[ordered.index(upper)]
        paired.append(dict(
            lower_level=lower, upper_level=upper,
            aggregate_b2_relative_change=[
                (high["aggregate_signed_current_b2"][index]
                 / low["aggregate_signed_current_b2"][index] - 1.0)
                for index in range(len(first_scales))]))
        paired[-1].update(
            aggregate_q1_b2_relative_change=[
                (high["aggregate_q1_resolved_signed_current_b2"][index]
                 / low["aggregate_q1_resolved_signed_current_b2"][index] - 1.0)
                for index in range(len(first_scales))],
            ensemble_mean_potential_rate_relative_change=(
                high["ensemble_mean_potential_rate_max_v_s"]
                / low["ensemble_mean_potential_rate_max_v_s"] - 1.0))

    audit = dict(
        schema="petch.charging.c3.sampling-floor-audit.v1",
        status="exact zero-step paired nested-sample audit",
        engine_git_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        campaign_script_sha256=_hash(Path(__file__).resolve()),
        checkpoint=dict(name=args.checkpoint.name, sha256=checkpoint_hash),
        method_map=dict(name=args.method_map.name, sha256=method_hash),
        protocol=dict(
            operator="exact hard visibility and declared charged-surface response",
            maximum_steps=0,
            independent_scramble_seeds=args.seeds,
            nested_forward_levels=ordered,
            nonlinear_aggregation=(
                "average signed positive/negative face currents across independent scrambles, "
                "then compute raw and Q1-resolved physical-patch B2 plus the Q1 potential rate"),
            compatible_q1_charge_state=args.compatible_q1_charge_state,
            transport_device=args.transport_device),
        patch_scales_m=first_scales,
        levels=level_records,
        paired_level_changes=paired,
        gate=0.08)
    audit_path = args.campaign_dir / "audit.json"
    _atomic_json(audit_path, audit)

    figure, axis = plt.subplots(figsize=(7.2, 4.6))
    sample_counts = np.asarray([2 ** level for level in ordered])
    aggregate_b2 = np.asarray(aggregate_b2)
    aggregate_q1_b2 = np.asarray(aggregate_q1_b2)
    for scale_index, scale in enumerate(first_scales):
        axis.plot(
            sample_counts, aggregate_b2[:, scale_index], marker="o",
            label=f"raw {scale * 1e6:.2f} um")
        axis.plot(
            sample_counts, aggregate_q1_b2[:, scale_index], marker="s", linestyle="--",
            label=f"Q1-resolved {scale * 1e6:.2f} um")
    axis.axhline(0.08, color="black", linestyle="--", linewidth=1.2, label="B2 gate")
    axis.set_xscale("log", base=2)
    axis.set_yscale("log")
    axis.set_xlabel("forward scrambled-Sobol samples per replicate")
    axis.set_ylabel("B2 of ensemble-mean signed currents")
    axis.set_title("C3 checkpoint sampling-floor decision audit")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(args.campaign_dir / "sampling_floor.png", dpi=180)
    plt.close(figure)
    print(json.dumps(_json_value(audit), indent=2))


if __name__ == "__main__":
    main()
