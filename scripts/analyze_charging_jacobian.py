"""Assemble and condition a fixed-map charging Jacobian from campaign checkpoints."""
import argparse
from pathlib import Path

import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--baseline", required=True)
parser.add_argument("--plus", required=True)
parser.add_argument("--minus", required=True)
parser.add_argument("--step", type=float, required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--stderr-floor", type=float, default=0.05)
args = parser.parse_args()
if args.step <= 0.0 or args.stderr_floor <= 0.0:
    raise ValueError("step and stderr floor must be positive")


def load(path):
    with np.load(path) as saved:
        return {name: saved[name] for name in saved.files}


baseline = load(args.baseline)
residual = np.asarray(baseline["mean_log_ratio_history"][-1], dtype=float)
dof_count = residual.size
jacobian = np.empty((dof_count, dof_count))
schemes = []
plus_directory = Path(args.plus)
minus_directory = Path(args.minus)
for dof in range(dof_count):
    plus = load(plus_directory / f"dof_{dof:02d}_plus.npz")
    minus = load(minus_directory / f"dof_{dof:02d}_minus.npz")
    plus_ok = str(plus["status"]) != "quadrature_failure"
    minus_ok = str(minus["status"]) != "quadrature_failure"
    if plus_ok and minus_ok:
        jacobian[:, dof] = (
            plus["mean_log_ratio_history"][-1]
            - minus["mean_log_ratio_history"][-1]) / (2.0 * args.step)
        schemes.append("central")
    elif plus_ok:
        jacobian[:, dof] = (
            plus["mean_log_ratio_history"][-1] - residual) / args.step
        schemes.append("forward")
    elif minus_ok:
        jacobian[:, dof] = (
            residual - minus["mean_log_ratio_history"][-1]) / args.step
        schemes.append("backward")
    else:
        raise RuntimeError(f"both perturbations failed quadrature certification for dof {dof}")

ion = np.asarray(baseline["ion_current"], dtype=float)
electron = np.asarray(baseline["electron_current"], dtype=float)
ion_stderr = np.asarray(baseline["ion_current_stderr"], dtype=float)
electron_stderr = np.asarray(baseline["electron_current_stderr"], dtype=float)
log_current_stderr = np.sqrt(
    (ion_stderr / np.maximum(ion, 1e-300)) ** 2
    + (electron_stderr / np.maximum(electron, 1e-300)) ** 2)
row_weight = 1.0 / np.maximum(log_current_stderr, args.stderr_floor)
row_weight /= np.median(row_weight)
singular_values = np.linalg.svd(jacobian, compute_uv=False)
weighted_singular_values = np.linalg.svd(row_weight[:, None] * jacobian, compute_uv=False)

np.savez(
    args.output, residual=residual, jacobian=jacobian,
    finite_difference_step_volts=np.asarray(args.step),
    finite_difference_scheme=np.asarray(schemes),
    log_current_stderr=log_current_stderr, row_weight=row_weight,
    singular_values=singular_values,
    weighted_singular_values=weighted_singular_values)

unique, counts = np.unique(schemes, return_counts=True)
print("schemes", dict(zip(unique.tolist(), counts.tolist())))
print("residual_max", float(np.max(np.abs(residual))))
print("residual_rms", float(np.sqrt(np.mean(residual ** 2))))
print("condition", float(singular_values[0] / singular_values[-1]))
print("weighted_condition", float(
    weighted_singular_values[0] / weighted_singular_values[-1]))
