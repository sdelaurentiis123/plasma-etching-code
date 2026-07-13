"""Audit a planar Maxwellian electron barrier term before any IMEX charging work."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

# Reuse the exact checkpoint reconstruction and physical boundary declaration from Task 0A.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import charging_task0a_ensemble_response as task0a  # noqa: E402

from petch.boundary_state import (  # noqa: E402
    PlasmaBoundaryState,
    maxwellian_electron_boundary_state,
    qmc_boundary_proposal_with_auxiliary,
)
from petch.boundary_transport import adjoint_boundary_state_face_flux  # noqa: E402


FLAT_VOLTAGES = (-12.0, -8.0, -4.0, -2.0, 0.0, 2.0, 4.0, 8.0)
TRENCH_VOLTAGES = FLAT_VOLTAGES
SEEDS = task0a.SEEDS
TEMPERATURE_EV = 4.0
_WORKER_MODEL = None


def barrier_factor(voltage_v, temperature_eV=TEMPERATURE_EV):
    voltage = np.asarray(voltage_v, dtype=float)
    return np.exp(np.minimum(voltage, 0.0) / float(temperature_eV))


def _flat_flux(voltage, level, seed, max_steps):
    nx, nz = 8, 8
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    cells = np.asarray([(index, nz - 1) for index in range(nx)], dtype=int)
    normals = np.tile((0, -1), (nx, 1))
    potential = np.empty((nx + 1, nz + 1))
    for z in range(nz + 1):
        potential[:, z] = voltage * min(z / (nz - 1), 1.0)
    electron = maxwellian_electron_boundary_state(
        TEMPERATURE_EV, 1e19, n_transverse=3, n_normal=6).get("electron")
    boundary = PlasmaBoundaryState((electron,), reference_plane_m=0.0)
    proposal, auxiliary = qmc_boundary_proposal_with_auxiliary(
        electron, int(level), 1, int(seed), name="electron-task0b-flat")
    estimate = adjoint_boundary_state_face_flux(
        boundary, "electron", potential, solid, cells, normals,
        proposal_species=proposal, face_position_samples=auxiliary[:, 0],
        fixed_dt=0.01, face_offset=1e-6, max_steps=int(max_steps),
        return_trajectory_outcomes=True)
    outcomes = estimate["trajectory_outcomes"][:, :, 0]
    return float(np.mean(estimate["per_face"])), float(np.mean(outcomes == 2))


def _regions(model):
    horizontal = model.normals[:, 1] == -1
    wall = model.normals[:, 0] != 0
    floor_z = int(np.max(model.cells[horizontal, 1]))
    wall_z = model.cells[wall, 1]
    split = float(np.median(wall_z))
    return {
        "top": np.where(horizontal & (model.cells[:, 1] < floor_z))[0],
        "upper_wall": np.where(wall & (model.cells[:, 1] <= split))[0],
        "lower_wall": np.where(wall & (model.cells[:, 1] > split))[0],
        "floor": np.where(horizontal & (model.cells[:, 1] == floor_z))[0],
    }


def _region_nodes(model, faces):
    selected = {node for face in faces for node in model.face_nodes[face] if node[1] != 0}
    lookup = {tuple(node): index for index, node in enumerate(model.dielectric_nodes)}
    return np.asarray([lookup[tuple(node)] for node in sorted(selected)], dtype=int)


def _trench_point(model, faces, coordinate_shift, level, seed):
    charge = model.base_charge.copy()
    node_indices = _region_nodes(model, faces)
    for index in node_indices:
        node = tuple(model.dielectric_nodes[index])
        charge[node] += coordinate_shift * model.capacitance[index]
    potential, _ = model.poisson.solve(charge)
    endpoint, _ = task0a._species_evaluation(
        model, "electron", potential, int(level), int(seed))
    face_voltage = np.asarray([
        np.mean([potential[node] for node in model.face_nodes[face]]) for face in faces])
    return dict(
        local_voltage_v=float(np.mean(face_voltage)),
        barrier_factor=float(np.mean(barrier_factor(face_voltage))),
        kinetic_flux=float(np.mean(np.sum(endpoint[faces], axis=1))))


def _region_voltage_calibration(model, faces):
    """Exploit Poisson linearity to map a charge-coordinate shift to mean face voltage."""
    values = []
    for coordinate_shift in (0.0, 1.0):
        charge = model.base_charge.copy()
        for index in _region_nodes(model, faces):
            node = tuple(model.dielectric_nodes[index])
            charge[node] += coordinate_shift * model.capacitance[index]
        potential, _ = model.poisson.solve(charge)
        values.append(float(np.mean([
            np.mean([potential[node] for node in model.face_nodes[face]]) for face in faces])))
    slope = values[1] - values[0]
    if not np.isfinite(slope) or abs(slope) < 1e-12:
        raise RuntimeError("regional charge coordinate does not change local surface voltage")
    return values[0], slope


def _initialize_worker(checkpoint):
    global _WORKER_MODEL
    _WORKER_MODEL = task0a.load_model(checkpoint)


def _flat_job(job):
    voltage, level, seed, max_steps = job
    return voltage, seed, _flat_flux(voltage, level, seed, max_steps)


def _trench_job(job):
    region, faces, target_voltage, shift, level, seed = job
    return (region, target_voltage, seed,
            _trench_point(_WORKER_MODEL, faces, shift, level, seed))


def _derivative_metrics(voltage, kinetic, analytic):
    order = np.argsort(voltage)
    voltage = np.asarray(voltage)[order]
    kinetic = np.asarray(kinetic)[order]
    analytic = np.asarray(analytic)[order]
    dkinetic = np.gradient(kinetic, voltage)
    danalytic = np.gradient(analytic, voltage)
    denominator = float(np.linalg.norm(dkinetic))
    capture = 1.0 - float(np.linalg.norm(dkinetic - danalytic)) / max(denominator, 1e-300)
    correlation = (float(np.corrcoef(dkinetic, danalytic)[0, 1])
                   if np.std(dkinetic) > 0.0 and np.std(danalytic) > 0.0 else 0.0)
    return order, dkinetic, danalytic, capture, correlation


def run(args):
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    model = task0a.load_model(args.checkpoint)
    config = dict(
        schema="petch.charging.task0b.v1",
        checkpoint=Path(args.checkpoint).name,
        checkpoint_sha256=task0a._sha256(args.checkpoint), level=args.level,
        seeds=list(args.seeds), flat_voltages=list(args.flat_voltages),
        flat_max_steps=args.flat_max_steps,
        trench_target_voltages=list(args.trench_voltages),
        hard_visibility=True, electron_temperature_eV=TEMPERATURE_EV)
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config["config_hash"] = hashlib.sha256(canonical).hexdigest()
    (output / "config.json").write_text(json.dumps(config, indent=2) + "\n")

    executor = ProcessPoolExecutor(
        max_workers=args.workers, initializer=_initialize_worker,
        initargs=(args.checkpoint,))
    flat_values = {}
    for voltage, seed, value in executor.map(
            _flat_job, [(voltage, args.level, seed, args.flat_max_steps)
                        for voltage in args.flat_voltages for seed in args.seeds]):
        flat_values.setdefault(voltage, {})[seed] = value
    flat_rows = []
    for voltage in args.flat_voltages:
        replicate = np.asarray([flat_values[voltage][seed][0] for seed in args.seeds])
        unresolved = np.asarray([flat_values[voltage][seed][1] for seed in args.seeds])
        flat_rows.append(dict(
            voltage_v=voltage, kinetic_flux=float(replicate.mean()),
            kinetic_stderr=float(replicate.std(ddof=1) / np.sqrt(replicate.size)),
            unresolved_fraction=float(unresolved.mean()),
            analytic_flux=float(barrier_factor(voltage))))
        print("flat", voltage, flat_rows[-1], flush=True)
    flat_error = np.asarray([
        row["kinetic_flux"] - row["analytic_flux"] for row in flat_rows])
    flat_pass = bool(
        np.sqrt(np.mean(flat_error ** 2)) <= 0.02
        and np.max(np.abs(flat_error)) <= 0.04
        and max(row["unresolved_fraction"] for row in flat_rows) <= 0.005)

    if not flat_pass:
        summary = dict(
            schema="petch.charging.task0b.summary.v1", config_hash=config["config_hash"],
            flat_gate_pass=False, flat_rmse=float(np.sqrt(np.mean(flat_error ** 2))),
            flat_max_abs_error=float(np.max(np.abs(flat_error))),
            flat_max_unresolved_fraction=float(max(
                row["unresolved_fraction"] for row in flat_rows)),
            trench_regions_passing=0, trench_regions_total=0,
            boundary_preconditioner_promoted=False,
            decision="drop boundary-current preconditioner; trench audit not entered",
            exact_operator="hard visibility kinetic current; analytic curve is diagnostic only")
        _write_csv(output / "flat_barrier.csv", flat_rows)
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        _plot(output, flat_rows, [], [])
        print(json.dumps(summary, indent=2))
        executor.shutdown(wait=True, cancel_futures=True)
        return

    trench_rows = []; region_summary = []
    regions = _regions(model)
    calibrations = {
        region: _region_voltage_calibration(model, faces)
        for region, faces in regions.items()
    }
    region_shifts = {
        region: {
            float(target): (float(target) - intercept) / slope
            for target in args.trench_voltages
        }
        for region, (intercept, slope) in calibrations.items()
    }
    trench_values = {}
    jobs = [(region, faces, target, region_shifts[region][float(target)], args.level, seed)
            for region, faces in regions.items() for target in args.trench_voltages
            for seed in args.seeds]
    for region, target, seed, value in executor.map(_trench_job, jobs):
        trench_values.setdefault((region, target), {})[seed] = value
    executor.shutdown(wait=True, cancel_futures=True)
    for region, faces in regions.items():
        points = []
        for target in args.trench_voltages:
            shift = region_shifts[region][float(target)]
            replicate = [trench_values[(region, target)][seed] for seed in args.seeds]
            points.append(dict(
                region=region, target_local_voltage_v=target, coordinate_shift_v=shift,
                local_voltage_v=float(np.mean([item["local_voltage_v"] for item in replicate])),
                barrier_factor=float(np.mean([item["barrier_factor"] for item in replicate])),
                kinetic_flux=float(np.mean([item["kinetic_flux"] for item in replicate])),
                kinetic_stderr=float(np.std(
                    [item["kinetic_flux"] for item in replicate], ddof=1)
                    / np.sqrt(len(replicate)))))
            print(region, target, points[-1], flush=True)
        factor = np.asarray([item["barrier_factor"] for item in points])
        kinetic = np.asarray([item["kinetic_flux"] for item in points])
        amplitude = float(np.dot(factor, kinetic) / max(np.dot(factor, factor), 1e-300))
        analytic = amplitude * factor
        voltage = np.asarray([item["local_voltage_v"] for item in points])
        order, dkinetic, danalytic, capture, correlation = _derivative_metrics(
            voltage, kinetic, analytic)
        derivative_by_index = {int(source): (dkinetic[target], danalytic[target])
                               for target, source in enumerate(order)}
        for index, point in enumerate(points):
            point["analytic_flux"] = float(analytic[index])
            point["kinetic_dJ_dV"] = float(derivative_by_index[index][0])
            point["analytic_dJ_dV"] = float(derivative_by_index[index][1])
            trench_rows.append(point)
        region_pass = bool(capture >= 0.5 and correlation >= 0.8)
        region_summary.append(dict(
            region=region, fitted_visibility_amplitude=amplitude,
            derivative_capture_fraction=capture,
            derivative_correlation=correlation, pass_stiffness=region_pass,
            clipping_required=bool(np.any(voltage >= 0.0))))
    most_regions = sum(item["pass_stiffness"] for item in region_summary) >= 3
    promoted = bool(flat_pass and most_regions)
    summary = dict(
        schema="petch.charging.task0b.summary.v1", config_hash=config["config_hash"],
        flat_gate_pass=flat_pass, flat_rmse=float(np.sqrt(np.mean(flat_error ** 2))),
        flat_max_abs_error=float(np.max(np.abs(flat_error))),
        flat_max_unresolved_fraction=float(max(
            row["unresolved_fraction"] for row in flat_rows)),
        trench_regions_passing=int(sum(item["pass_stiffness"] for item in region_summary)),
        trench_regions_total=len(region_summary),
        boundary_preconditioner_promoted=promoted,
        decision=("promote Task 2b after Task 1/2a entry gates" if promoted else
                  "drop boundary-current preconditioner" if not flat_pass else
                  "do not promote: planar term does not capture most trench stiffness"),
        exact_operator="hard visibility kinetic current; analytic curve is diagnostic only")
    _write_csv(output / "flat_barrier.csv", flat_rows)
    _write_csv(output / "trench_regions.csv", trench_rows)
    _write_csv(output / "region_summary.csv", region_summary)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _plot(output, flat_rows, trench_rows, region_summary)
    print(json.dumps(summary, indent=2))


def _write_csv(path, rows):
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _plot(output, flat_rows, trench_rows, region_summary):
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    voltage = [row["voltage_v"] for row in flat_rows]
    axes[0].errorbar(
        voltage, [row["kinetic_flux"] for row in flat_rows],
        yerr=[row["kinetic_stderr"] for row in flat_rows], fmt="o-", label="kinetic")
    axes[0].plot(voltage, [row["analytic_flux"] for row in flat_rows], "--", label="barrier")
    axes[0].set_title("Flat Maxwellian gate"); axes[0].set_xlabel("surface voltage (V)")
    axes[0].set_ylabel("normalized electron flux"); axes[0].legend(); axes[0].grid(alpha=0.25)
    for summary in region_summary:
        rows = [row for row in trench_rows if row["region"] == summary["region"]]
        rows.sort(key=lambda row: row["local_voltage_v"])
        axes[1].plot(
            [row["local_voltage_v"] for row in rows],
            [row["kinetic_dJ_dV"] for row in rows], "o-", label=f"{summary['region']} kinetic")
        axes[1].plot(
            [row["local_voltage_v"] for row in rows],
            [row["analytic_dJ_dV"] for row in rows], "--", alpha=0.8)
    axes[1].axvline(0.0, color="black", lw=0.8)
    axes[1].set_title("Trench regional stiffness (dashed = barrier)")
    axes[1].set_xlabel("mean local surface voltage (V)"); axes[1].set_ylabel("dJ/dV")
    axes[1].grid(alpha=0.25); axes[1].legend(fontsize=7, ncol=2)
    figure.savefig(output / "electron_boundary_audit.png", dpi=180)
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--level", type=int, default=10)
    parser.add_argument("--flat-max-steps", type=int, default=25600)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    parser.add_argument("--flat-voltages", type=float, nargs="+", default=list(FLAT_VOLTAGES))
    parser.add_argument("--trench-voltages", type=float, nargs="+", default=list(TRENCH_VOLTAGES))
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
