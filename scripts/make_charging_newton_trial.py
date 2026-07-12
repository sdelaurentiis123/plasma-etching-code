"""Build a Poisson-consistent AR charging trial from a fixed-map Jacobian artifact."""
import argparse

import numpy as np

from petch.charging_nodal import material_face_nodes
from petch.charging_poisson import NodalPoissonSystem


parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", required=True)
parser.add_argument("--jacobian", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--regularization", type=float, required=True)
parser.add_argument("--fraction", type=float, default=1.0)
parser.add_argument("--maximum-coordinate-step", type=float, default=1.0)
parser.add_argument("--maximum-row-weight", type=float, default=0.0)
parser.add_argument("--epsilon-solid", type=float, default=3.9)
args = parser.parse_args()
if (args.regularization <= 0.0 or not 0.0 < args.fraction <= 1.0
        or args.maximum_coordinate_step <= 0.0 or args.epsilon_solid <= 0.0):
    raise ValueError("regularization, fraction, step limit, and permittivity must be positive")
if args.maximum_row_weight < 0.0:
    raise ValueError("maximum row weight must be nonnegative")


def load(path):
    with np.load(path) as saved:
        return {name: saved[name] for name in saved.files}


checkpoint = load(args.checkpoint)
linearization = load(args.jacobian)
residual = np.asarray(checkpoint["mean_log_ratio_history"][-1], dtype=float)
jacobian = np.asarray(linearization["jacobian"], dtype=float)
if jacobian.shape != (residual.size, residual.size):
    raise ValueError("Jacobian and checkpoint charging dimensions do not match")

ion = np.asarray(checkpoint["ion_current"], dtype=float)
electron = np.asarray(checkpoint["electron_current"], dtype=float)
ion_stderr = np.asarray(checkpoint["ion_current_stderr"], dtype=float)
electron_stderr = np.asarray(checkpoint["electron_current_stderr"], dtype=float)
log_stderr = np.sqrt(
    (ion_stderr / np.maximum(ion, 1e-300)) ** 2
    + (electron_stderr / np.maximum(electron, 1e-300)) ** 2)
row_weight = 1.0 / np.maximum(log_stderr, 0.05)
row_weight /= np.median(row_weight)
row_weight *= 1.0 + args.maximum_row_weight * (
    np.abs(residual) / max(float(np.max(np.abs(residual))), 1e-300)) ** 4
row_weight /= np.median(row_weight)
identity = np.eye(residual.size)
system = np.vstack((row_weight[:, None] * jacobian, args.regularization * identity))
target = np.concatenate((-row_weight * residual, np.zeros(residual.size)))
step = args.fraction * np.linalg.lstsq(system, target, rcond=1e-10)[0]
largest = float(np.max(np.abs(step)))
step_scale = min(1.0, args.maximum_coordinate_step / max(largest, 1e-300))
step *= step_scale

solid = np.asarray(checkpoint["solid"], dtype=bool)
epsilon_r = np.ones_like(solid, dtype=float)
epsilon_r[solid] = args.epsilon_solid
fixed = np.zeros((solid.shape[0] + 1, solid.shape[1] + 1), dtype=bool)
fixed[:, 0] = True
fixed[:, -1] = True
poisson = NodalPoissonSystem(epsilon_r, fixed, np.zeros(fixed.shape))
nodes = np.asarray(checkpoint["dielectric_nodes"], dtype=int)
if nodes.shape != (residual.size, 2):
    raise ValueError("trial builder currently requires all charging degrees of freedom to be dielectric")
capacitance = poisson.diagonal_surface_capacitance(nodes)
charge = np.asarray(checkpoint["surface_charge_node_c_per_m"], dtype=float).copy()
charge[nodes[:, 0], nodes[:, 1]] += capacitance * step
potential, diagnostics = poisson.solve(charge)
surface_voltage = np.zeros_like(solid, dtype=float)
for cell, normal in zip(checkpoint["cells"], checkpoint["normals"]):
    cell_tuple = tuple(map(int, cell))
    endpoints = material_face_nodes(cell_tuple, tuple(map(int, normal)))
    surface_voltage[cell_tuple] = np.mean([potential[node] for node in endpoints])

predicted = residual + jacobian @ step
payload = dict(checkpoint)
payload.update(
    status=np.asarray("newton_trial"), surface_charge_node_c_per_m=charge,
    boundary_nodal_voltage=potential, potential=potential, surface_voltage=surface_voltage,
    anderson_x_history=np.empty((0, residual.size)),
    anderson_residual_history=np.empty((0, residual.size)),
    newton_step_coordinate_volts=step,
    newton_regularization=np.asarray(args.regularization),
    newton_fraction=np.asarray(args.fraction), newton_step_scale=np.asarray(step_scale),
    newton_maximum_row_weight=np.asarray(args.maximum_row_weight),
    newton_predicted_raw_max=np.asarray(np.max(np.abs(predicted))),
    newton_predicted_raw_rms=np.asarray(np.sqrt(np.mean(predicted ** 2))),
    newton_poisson_charge_balance_c_per_m=np.asarray(
        diagnostics.charge_balance_c_per_m))
np.savez(args.output, **payload)

print("step_scale", step_scale)
print("step_max", float(np.max(np.abs(step))))
print("step_rms", float(np.sqrt(np.mean(step ** 2))))
print("predicted_raw_max", float(np.max(np.abs(predicted))))
print("predicted_raw_rms", float(np.sqrt(np.mean(predicted ** 2))))
print("charge_closure", diagnostics.charge_balance_c_per_m)
