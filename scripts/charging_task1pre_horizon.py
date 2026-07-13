"""Paired trajectory-horizon audit for the canonical stuck charging current map."""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import sys
import time

# Numba's workqueue backend is not thread-safe. Process workers must not start nested math threads.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import charging_task0a_ensemble_response as task0a  # noqa: E402


_MODEL = None
_POTENTIAL = None


def _initialize_worker(checkpoint):
    global _MODEL, _POTENTIAL
    _MODEL = task0a.load_model(checkpoint)
    _POTENTIAL, _ = _MODEL.poisson.solve(_MODEL.base_charge)


def _unresolved(detail):
    result = {}
    for method in ("forward", "adjoint"):
        if method not in detail:
            result[method] = np.nan
            continue
        status = np.asarray(detail[method]["trajectory_outcomes"])[..., 0]
        result[method] = float(np.mean(status == 2))
    return result


def _evaluate(job):
    label, fixed_dt, max_steps, level, seed, config_hash, raw_path = job
    raw_path = Path(raw_path)
    if raw_path.exists():
        with np.load(raw_path) as saved:
            if str(saved.get("config_hash", "")) == config_hash:
                return raw_path.name, True, float(saved["elapsed_seconds"])
    started = time.perf_counter()
    currents = {}; unresolved = {}
    for species_name in ("ion", "electron"):
        endpoint, detail = task0a._species_evaluation(
            _MODEL, species_name, _POTENTIAL, int(level), int(seed),
            max_steps=int(max_steps), fixed_dt=float(fixed_dt))
        currents[species_name] = task0a._assemble_nodes(
            _MODEL, species_name, endpoint)
        unresolved[species_name] = _unresolved(detail)
    ion = currents["ion"]; electron = currents["electron"]
    active = ion + electron >= 1e-4
    log_ratio = np.log((ion + 1e-12) / (electron + 1e-12))
    imbalance = np.divide(
        ion - electron, ion + electron,
        out=np.zeros_like(ion), where=(ion + electron) > 0.0)
    elapsed = time.perf_counter() - started
    np.savez(
        raw_path, schema="petch.charging.task1pre.horizon.raw.v1",
        config_hash=config_hash, label=label, fixed_dt=fixed_dt,
        max_steps=max_steps, level=level, seed=seed,
        ion_current=ion, electron_current=electron,
        active=active, log_ratio=log_ratio, imbalance=imbalance,
        ion_forward_unresolved=unresolved["ion"]["forward"],
        ion_adjoint_unresolved=unresolved["ion"]["adjoint"],
        electron_forward_unresolved=unresolved["electron"]["forward"],
        electron_adjoint_unresolved=unresolved["electron"]["adjoint"],
        elapsed_seconds=elapsed)
    return raw_path.name, False, elapsed


def _stderr(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return 0.0
    return float(finite.std(ddof=1) / np.sqrt(finite.size))


def _write_csv(path, rows):
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _load_raw(path):
    with np.load(path) as saved:
        return {name: saved[name] for name in saved.files}


def _analyze(output, config):
    raw = output / "raw"
    by_variant = {}
    for variant in config["variants"]:
        label = variant["label"]
        by_variant[label] = [
            _load_raw(raw / f"{label}_s{seed}.npz") for seed in config["seeds"]]
    reference_label = config["reference_label"]
    reference = by_variant[reference_label]
    rows = []
    for variant in config["variants"]:
        label = variant["label"]; items = by_variant[label]
        ion_mean = np.mean([item["ion_current"] for item in items], axis=0)
        electron_mean = np.mean([item["electron_current"] for item in items], axis=0)
        ref_ion = np.mean([item["ion_current"] for item in reference], axis=0)
        ref_electron = np.mean([item["electron_current"] for item in reference], axis=0)
        active = np.any([item["active"] for item in items + reference], axis=0)
        log_ratio = np.log((ion_mean + 1e-12) / (electron_mean + 1e-12))
        imbalance = np.divide(
            ion_mean - electron_mean, ion_mean + electron_mean,
            out=np.zeros_like(ion_mean), where=(ion_mean + electron_mean) > 0.0)
        ref_imbalance = np.divide(
            ref_ion - ref_electron, ref_ion + ref_electron,
            out=np.zeros_like(ref_ion), where=(ref_ion + ref_electron) > 0.0)
        ion_change = float(
            np.linalg.norm(ion_mean - ref_ion) / max(np.linalg.norm(ref_ion), 1e-300))
        electron_change = float(
            np.linalg.norm(electron_mean - ref_electron)
            / max(np.linalg.norm(ref_electron), 1e-300))
        paired_rms_changes = []
        for item, ref in zip(items, reference):
            selected = np.asarray(item["active"] | ref["active"], dtype=bool)
            paired_rms_changes.append(float(np.sqrt(np.mean(
                (item["imbalance"][selected] - ref["imbalance"][selected]) ** 2))))
        unresolved_names = (
            "ion_forward_unresolved", "ion_adjoint_unresolved",
            "electron_forward_unresolved", "electron_adjoint_unresolved")
        unresolved_values = {
            name: float(np.nanmean([float(item[name]) for item in items]))
            if np.any(np.isfinite([float(item[name]) for item in items])) else np.nan
            for name in unresolved_names
        }
        max_unresolved = float(np.nanmax(list(unresolved_values.values())))
        rms = float(np.sqrt(np.mean(imbalance[active] ** 2)))
        worst = float(np.max(np.abs(imbalance[active])))
        ref_rms = float(np.sqrt(np.mean(ref_imbalance[active] ** 2)))
        ref_worst = float(np.max(np.abs(ref_imbalance[active])))
        passes = bool(
            max_unresolved <= config["gates"]["max_unresolved_fraction"]
            and ion_change <= config["gates"]["max_relative_current_change"]
            and electron_change <= config["gates"]["max_relative_current_change"]
            and abs(rms - ref_rms) <= config["gates"]["max_rms_imbalance_change"]
            and abs(worst - ref_worst) <= config["gates"]["max_worst_imbalance_change"])
        rows.append(dict(
            label=label, fixed_dt=variant["fixed_dt"], max_steps=variant["max_steps"],
            physical_horizon=variant["fixed_dt"] * variant["max_steps"],
            level=config["level"], replicates=len(items),
            ion_forward_unresolved=unresolved_values["ion_forward_unresolved"],
            ion_adjoint_unresolved=unresolved_values["ion_adjoint_unresolved"],
            electron_forward_unresolved=unresolved_values["electron_forward_unresolved"],
            electron_adjoint_unresolved=unresolved_values["electron_adjoint_unresolved"],
            max_unresolved=max_unresolved,
            ion_relative_change_vs_reference=ion_change,
            electron_relative_change_vs_reference=electron_change,
            paired_rms_change_vs_reference=float(np.mean(paired_rms_changes)),
            paired_rms_change_stderr=_stderr(paired_rms_changes),
            log_residual_rms=float(np.sqrt(np.mean(log_ratio[active] ** 2))),
            log_residual_worst=float(np.max(np.abs(log_ratio[active]))),
            imbalance_rms=rms, imbalance_worst=worst,
            imbalance_rms_change_vs_reference=abs(rms - ref_rms),
            imbalance_worst_change_vs_reference=abs(worst - ref_worst),
            pass_reference_gate=passes))
    _write_csv(output / "horizon_summary.csv", rows)
    production = [row for row in rows if row["fixed_dt"] == 0.01]
    passing = [row for row in production if row["pass_reference_gate"]]
    selected = min(passing, key=lambda row: row["physical_horizon"]) if passing else None
    dt_row = next(row for row in rows if row["label"] == config["dt_comparison_label"])
    reference_row = next(row for row in rows if row["label"] == reference_label)
    dt_pass = bool(
        dt_row["ion_relative_change_vs_reference"]
        <= config["gates"]["max_relative_current_change"]
        and dt_row["electron_relative_change_vs_reference"]
        <= config["gates"]["max_relative_current_change"]
        and dt_row["imbalance_rms_change_vs_reference"]
        <= config["gates"]["max_rms_imbalance_change"]
        and dt_row["imbalance_worst_change_vs_reference"]
        <= config["gates"]["max_worst_imbalance_change"])
    summary = dict(
        schema="petch.charging.task1pre.horizon.summary.v1",
        config_hash=config["config_hash"], checkpoint=config["checkpoint"],
        checkpoint_sha256=config["checkpoint_sha256"], exact_operator="hard visibility",
        reference_label=reference_label,
        reference_imbalance_rms=reference_row["imbalance_rms"],
        reference_imbalance_worst=reference_row["imbalance_worst"],
        timestep_refinement_pass=dt_pass,
        selected_production_label=(None if selected is None else selected["label"]),
        selected_production_max_steps=(None if selected is None else selected["max_steps"]),
        horizon_gate_pass=bool(selected is not None and dt_pass),
        decision=("Task 1 horizon entry gate passed" if selected is not None and dt_pass else
                  "hold Task 1; extend horizon/timestep audit"))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _plot(output, rows, reference_label)
    return summary


def _plot(output, rows, reference_label):
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)
    dt01 = sorted(
        [row for row in rows if row["fixed_dt"] == 0.01],
        key=lambda row: row["physical_horizon"])
    horizon = [row["physical_horizon"] for row in dt01]
    axes[0].plot(horizon, [row["ion_forward_unresolved"] for row in dt01], "o-",
                 label="ion forward")
    axes[0].plot(horizon, [row["ion_adjoint_unresolved"] for row in dt01], "o-",
                 label="ion adjoint")
    axes[0].plot(horizon, [row["electron_adjoint_unresolved"] for row in dt01], "o-",
                 label="electron adjoint")
    axes[0].axhline(0.005, color="black", ls="--", lw=1, label="gate")
    axes[0].set_xscale("log"); axes[0].set_yscale("symlog", linthresh=1e-5)
    axes[0].set_xlabel("physical trace horizon"); axes[0].set_ylabel("unresolved fraction")
    axes[0].set_title("Trajectory completion"); axes[0].legend(fontsize=7)
    axes[1].plot(horizon, [row["imbalance_rms"] for row in dt01], "o-", label="RMS")
    axes[1].plot(horizon, [row["imbalance_worst"] for row in dt01], "o-", label="worst")
    axes[1].set_xscale("log"); axes[1].set_xlabel("physical trace horizon")
    axes[1].set_ylabel("|Ii-Ie| / (Ii+Ie)"); axes[1].set_title("Stuck-state residual")
    axes[1].legend()
    axes[2].plot(horizon, [row["ion_relative_change_vs_reference"] for row in dt01],
                 "o-", label="ion")
    axes[2].plot(horizon, [row["electron_relative_change_vs_reference"] for row in dt01],
                 "o-", label="electron")
    axes[2].axhline(0.01, color="black", ls="--", lw=1, label="gate")
    axes[2].set_xscale("log"); axes[2].set_yscale("symlog", linthresh=1e-5)
    axes[2].set_xlabel("physical trace horizon"); axes[2].set_ylabel("relative change vs reference")
    axes[2].set_title("Current-map sensitivity"); axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.suptitle(f"Task 1-pre horizon audit; reference={reference_label}")
    figure.savefig(output / "trajectory_horizon_audit.png", dpi=180)
    plt.close(figure)


def run(args):
    output = Path(args.output_dir); raw = output / "raw"
    output.mkdir(parents=True, exist_ok=True); raw.mkdir(parents=True, exist_ok=True)
    model = task0a.load_model(args.checkpoint)
    nz = model.solid.shape[1]
    variants = [
        dict(label="default", fixed_dt=0.01, max_steps=200 * nz),
        dict(label="horizon_x4", fixed_dt=0.01, max_steps=800 * nz),
        dict(label="horizon_x8", fixed_dt=0.01, max_steps=1600 * nz),
        dict(label="dt_half_reference", fixed_dt=0.005, max_steps=3200 * nz),
    ]
    config = dict(
        schema="petch.charging.task1pre.horizon.config.v1",
        checkpoint=Path(args.checkpoint).name,
        checkpoint_sha256=task0a._sha256(args.checkpoint),
        level=args.level, seeds=list(args.seeds), variants=variants,
        reference_label="dt_half_reference", dt_comparison_label="horizon_x8",
        common_random_numbers=True, hard_visibility=True,
        gates=dict(
            max_unresolved_fraction=0.005, max_relative_current_change=0.01,
            max_rms_imbalance_change=0.01, max_worst_imbalance_change=0.02))
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config["config_hash"] = hashlib.sha256(canonical).hexdigest()
    (output / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    jobs = []
    for variant in variants:
        for seed in args.seeds:
            path = raw / f"{variant['label']}_s{seed}.npz"
            jobs.append((
                variant["label"], variant["fixed_dt"], variant["max_steps"],
                args.level, seed, config["config_hash"], str(path)))
    with ProcessPoolExecutor(
            max_workers=args.workers, initializer=_initialize_worker,
            initargs=(args.checkpoint,)) as pool:
        futures = [pool.submit(_evaluate, job) for job in jobs]
        for index, future in enumerate(as_completed(futures), 1):
            name, reused, elapsed = future.result()
            print(f"[{index}/{len(jobs)}] {name} "
                  f"{'reused' if reused else f'{elapsed:.1f}s'}", flush=True)
    summary = _analyze(output, config)
    print(json.dumps(summary, indent=2))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--level", type=int, default=9)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(task0a.SEEDS))
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
