"""Reproducible AR charging, restart, fixed-map, and local-response campaign driver."""
import argparse
import time

import numpy as np

from petch.boundary_state import (
    MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    folded_normal_tangential_proposal,
    maxwellian_electron_boundary_state,
    mixture_boundary_proposal,
    qmc_boundary_proposal,
)
from petch.charging_backward import AdaptiveQuadratureConvergenceError, solve_boundary_state_charging
from petch.charging_backward import _gas_faces
from petch.charging_nodal import material_face_nodes
from petch.charging_nodal_fixed_point import solve_boundary_state_charging_nodal
from petch.charging_poisson import NodalPoissonSystem
from petch.sheath import CollisionlessRFSheath


parser = argparse.ArgumentParser()
parser.add_argument("--iterations", type=int, default=40)
parser.add_argument("--beta", type=float, default=0.05)
parser.add_argument("--gain-decay", type=float, default=0.0)
parser.add_argument("--gain-offset", type=float, default=5.0)
parser.add_argument("--dvmax", type=float, default=0.5)
parser.add_argument("--energy-bins", type=int, default=64)
parser.add_argument("--adjoint-base", type=int, default=6)
parser.add_argument("--adjoint-max", type=int, default=12)
parser.add_argument("--forward-base", type=int, default=8)
parser.add_argument("--forward-max", type=int, default=14)
parser.add_argument("--element-absolute", type=float, default=0.01)
parser.add_argument("--element-relative", type=float, default=0.15)
parser.add_argument("--fixed-dt", type=float, default=0.0)
parser.add_argument("--face-offset", type=float, default=1e-3)
parser.add_argument("--grazing", action="store_true")
parser.add_argument("--freeze-method", action="store_true")
parser.add_argument("--freeze-levels", action="store_true")
parser.add_argument("--update", choices=("picard", "anderson"), default="picard")
parser.add_argument("--anderson-depth", type=int, default=4)
parser.add_argument("--output", default="/tmp/petch_charging_diag_result.npz")
parser.add_argument("--initial", default=None)
parser.add_argument("--initial-state", choices=("accepted", "rejected"), default="accepted")
parser.add_argument("--clear-acceleration-history", action="store_true")
parser.add_argument("--override-restart-beta", type=float)
parser.add_argument("--geometry", choices=("bulk", "sheet"), default="bulk")
parser.add_argument("--trench-width", type=int, default=10)
parser.add_argument("--trench-depth", type=int, default=14)
parser.add_argument("--side-thickness", type=int, default=5)
parser.add_argument("--poisson", action="store_true")
parser.add_argument("--cell-size-nm", type=float, default=31.25)
parser.add_argument("--epsilon-solid", type=float, default=3.9)
parser.add_argument("--perturb-node-i", type=int, default=None)
parser.add_argument("--perturb-node-j", type=int, default=None)
parser.add_argument("--perturb-node-volts", type=float, default=0.0)
parser.add_argument("--perturb-charge-dof", type=int, default=None)
parser.add_argument("--perturb-charge-coordinate-volts", type=float, default=0.0)
parser.add_argument("--trust", action="store_true")
parser.add_argument("--trust-merit", choices=("rms", "max", "pareto"), default="rms")
parser.add_argument("--nodal", action="store_true")
args = parser.parse_args()

if args.trench_width <= 0 or args.trench_depth <= 0 or args.side_thickness <= 0:
    raise ValueError("trench dimensions must be positive")
left = args.side_thickness
right = left + args.trench_width
floor_z = 1 + args.trench_depth
nx, nz = right + args.side_thickness, floor_z + 3
solid = np.zeros((nx, nz), dtype=bool)
if args.geometry == "sheet":
    solid[left, 1:] = True
    solid[right - 1, 1:] = True
    solid[left:right, floor_z:] = True
else:
    solid[:left, 1:] = True
    solid[right:, 1:] = True
    solid[left:right, floor_z:] = True
conductors = np.zeros_like(solid, dtype=int)

sheath = CollisionlessRFSheath(
    40.0, 10.0, 2e6, 4.0, 40.0, thickness_m=5e-4)
ion = collisionless_sheath_boundary_state(
    sheath, 1e19, n_phase=16, tangential_temperature_eV=0.2,
    n_transverse=3, normal_energy_bins=args.energy_bins).get("ion")
electron = maxwellian_electron_boundary_state(
    4.0, 1e19, n_transverse=3, n_normal=6).get("electron")
boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)

ion_tail = SpeciesBoundaryState(
    "ion_tail", 1, 40.0, 1.0, [[0.0, 0.0, np.sqrt(40.0)]], [1.0],
    density_model=MaxwellianFluxVelocityDensity(40.0))
electron_tail = SpeciesBoundaryState(
    "electron_tail", -1, electron.mass_amu, 1.0,
    [[0.0, 0.0, np.sqrt(16.0)]], [1.0],
    density_model=MaxwellianFluxVelocityDensity(16.0))
ion_proposal = mixture_boundary_proposal((
    qmc_boundary_proposal(ion, 8, seed=101),
    qmc_boundary_proposal(ion_tail, 8, seed=102),
), (0.8, 0.2), name="ion_proposal")
electron_proposal = mixture_boundary_proposal((
    qmc_boundary_proposal(electron, 8, seed=201),
    qmc_boundary_proposal(electron_tail, 8, seed=202),
), (0.8, 0.2), name="electron_proposal")
if args.grazing:
    ion_proposal = mixture_boundary_proposal((
        ion,
        folded_normal_tangential_proposal(ion, +1),
        folded_normal_tangential_proposal(ion, -1),
        ion_tail,
    ), (0.4, 0.25, 0.25, 0.1), name="ion_grazing_proposal")
    electron_proposal = mixture_boundary_proposal((
        electron,
        folded_normal_tangential_proposal(electron, +1),
        folded_normal_tangential_proposal(electron, -1),
        electron_tail,
    ), (0.4, 0.25, 0.25, 0.1), name="electron_grazing_proposal")

adaptive = dict(
    bidirectional=True,
    base_log2=args.adjoint_base,
    max_log2=args.adjoint_max,
    n_replicates=4,
    absolute_tolerance=0.01,
    relative_tolerance=0.05,
    element_absolute_tolerance=args.element_absolute,
    element_relative_tolerance=args.element_relative,
    method_switch_factor=2.0,
    fixed_dt=args.fixed_dt,
    face_offset=args.face_offset,
    freeze_method_hint=args.freeze_method,
    freeze_levels=args.freeze_levels,
    forward_options=dict(
        base_log2=args.forward_base,
        max_log2=args.forward_max,
        n_replicates=4,
        absolute_tolerance=0.01,
        relative_tolerance=0.05,
        element_absolute_tolerance=args.element_absolute,
        element_relative_tolerance=args.element_relative,
        fixed_dt=args.fixed_dt,
        freeze_levels=args.freeze_levels,
    ),
)

start = time.perf_counter()
initial_surface_voltage = None
initial_boundary_nodal_voltage = None
initial_surface_charge_node_c_per_m = None
initial_adaptive_levels = {}
initial_forward_adaptive_levels = {}
initial_method_hint = {}
initial_accepted_iterations = 0
initial_beta = None
initial_anderson_x = None
initial_anderson_residual = None
initial_trust_best_rms = None
initial_trust_best_max = None
if args.initial is not None:
    with np.load(args.initial) as saved:
        if "solid" in saved and not np.array_equal(saved["solid"], solid):
            raise ValueError("checkpoint material topology does not match requested geometry")
        prefix = "failed_" if args.initial_state == "rejected" else ""
        if args.initial_state == "rejected" and str(saved.get("status", "")) != "quadrature_failure":
            raise ValueError("a rejected initial state requires a quadrature-failure checkpoint")
        initial_surface_voltage = saved[f"{prefix}surface_voltage"]
        if f"{prefix}boundary_nodal_voltage" in saved:
            initial_boundary_nodal_voltage = saved[f"{prefix}boundary_nodal_voltage"]
        elif args.initial_state == "rejected" and "failed_potential" in saved:
            initial_boundary_nodal_voltage = saved["failed_potential"]
        elif args.nodal and "potential" in saved:
            initial_boundary_nodal_voltage = saved["potential"]
        if args.perturb_node_i is not None or args.perturb_node_j is not None:
            if args.perturb_node_i is None or args.perturb_node_j is None:
                raise ValueError("both perturbation node indices are required")
            initial_boundary_nodal_voltage = initial_boundary_nodal_voltage.copy()
            initial_boundary_nodal_voltage[
                args.perturb_node_i, args.perturb_node_j] += args.perturb_node_volts
        if f"{prefix}surface_charge_node_c_per_m" in saved:
            initial_surface_charge_node_c_per_m = saved[
                f"{prefix}surface_charge_node_c_per_m"]
        for species_name in ("ion", "electron"):
            if f"{prefix}adaptive_{species_name}" in saved:
                initial_adaptive_levels[species_name] = saved[
                    f"{prefix}adaptive_{species_name}"]
            if f"{prefix}forward_adaptive_{species_name}" in saved:
                initial_forward_adaptive_levels[species_name] = saved[
                    f"{prefix}forward_adaptive_{species_name}"]
            if f"{prefix}method_hint_{species_name}" in saved:
                initial_method_hint[species_name] = saved[
                    f"{prefix}method_hint_{species_name}"]
        if "restart_accepted_iterations" in saved:
            initial_accepted_iterations = int(saved["restart_accepted_iterations"])
        elif "accepted_iterations_total" in saved:
            initial_accepted_iterations = int(saved["accepted_iterations_total"])
        if "restart_beta" in saved:
            initial_beta = float(saved["restart_beta"])
        if "anderson_x_history" in saved:
            initial_anderson_x = saved["anderson_x_history"]
            initial_anderson_residual = saved["anderson_residual_history"]
        if "trust_best_rms" in saved:
            initial_trust_best_rms = float(saved["trust_best_rms"])
        if "trust_best_max" in saved:
            initial_trust_best_max = float(saved["trust_best_max"])
    if args.initial_state == "rejected":
        if args.poisson and args.nodal and initial_surface_charge_node_c_per_m is None:
            raise ValueError("rejected Poisson state is missing its physical surface charge")
        if (set(initial_adaptive_levels) != {"ion", "electron"}
                or set(initial_forward_adaptive_levels) != {"ion", "electron"}
                or set(initial_method_hint) != {"ion", "electron"}):
            raise ValueError("rejected state is missing its frozen estimator rule")
        # Refining a rejected trial changes the deterministic sample-average map. Secant/Anderson
        # history from the preceding rule is no longer a Jacobian approximation to this map.
        initial_anderson_x = None
        initial_anderson_residual = None
if args.clear_acceleration_history:
    initial_anderson_x = None
    initial_anderson_residual = None
if args.override_restart_beta is not None:
    if not np.isfinite(args.override_restart_beta) or args.override_restart_beta <= 0.0:
        raise ValueError("override restart beta must be finite and positive")
    initial_beta = float(args.override_restart_beta)
if args.perturb_charge_dof is not None:
    if not args.poisson or not args.nodal or initial_surface_charge_node_c_per_m is None:
        raise ValueError("charge-coordinate perturbation requires a nodal Poisson checkpoint")
    cells_for_dof, normals_for_dof = _gas_faces(solid, solid)
    dielectric_nodes_for_dof = sorted({
        node
        for cell, normal in zip(cells_for_dof, normals_for_dof)
        for node in material_face_nodes(cell, normal)
        if node[1] != 0
    })
    if not 0 <= args.perturb_charge_dof < len(dielectric_nodes_for_dof):
        raise ValueError("perturb-charge-dof is outside the dielectric state vector")
    epsilon_for_dof = np.ones_like(solid, dtype=float)
    epsilon_for_dof[solid] = args.epsilon_solid
    fixed_for_dof = np.zeros((nx + 1, nz + 1), dtype=bool)
    fixed_for_dof[:, 0] = True
    fixed_for_dof[:, -1] = True
    poisson_for_dof = NodalPoissonSystem(
        epsilon_for_dof, fixed_for_dof, np.zeros_like(fixed_for_dof, dtype=float))
    node_for_dof = dielectric_nodes_for_dof[args.perturb_charge_dof]
    capacitance_for_dof = poisson_for_dof.diagonal_surface_capacitance(
        np.asarray([node_for_dof], dtype=int))[0]
    initial_surface_charge_node_c_per_m = initial_surface_charge_node_c_per_m.copy()
    initial_surface_charge_node_c_per_m[node_for_dof] += (
        capacitance_for_dof * args.perturb_charge_coordinate_volts)
try:
    solve_options = dict(
        initial_surface_voltage=initial_surface_voltage,
        n_iter=args.iterations, min_iter=2, balance_tol=0.15,
        beta=args.beta, response_energy_eV=4.0, dVmax=args.dvmax,
        field_sweeps=150, boundary_proposals={
            "ion": ion_proposal, "electron": electron_proposal},
        adaptive_quadrature=adaptive)
    if args.nodal:
        poisson_options = {}
        if args.poisson:
            epsilon_r = np.ones_like(solid, dtype=float)
            epsilon_r[solid] = args.epsilon_solid
            grounded = np.zeros((nx + 1, nz + 1), dtype=bool)
            grounded[:, -1] = True
            poisson_options = dict(
                epsilon_r=epsilon_r, cell_size_m=args.cell_size_nm * 1e-9,
                grounded_nodes=grounded,
                initial_surface_charge_node_c_per_m=initial_surface_charge_node_c_per_m)
        result = solve_boundary_state_charging_nodal(
            solid, conductors, boundary,
            initial_boundary_nodal_voltage=initial_boundary_nodal_voltage,
            trust_region=args.trust, gain_decay=args.gain_decay,
            gain_offset=args.gain_offset,
            initial_adaptive_levels=initial_adaptive_levels,
            initial_forward_adaptive_levels=initial_forward_adaptive_levels,
            initial_method_hint=initial_method_hint,
            initial_accepted_iterations=initial_accepted_iterations,
            initial_beta=initial_beta,
            initial_anderson_x=initial_anderson_x,
            initial_anderson_residual=initial_anderson_residual,
            initial_trust_best_rms=initial_trust_best_rms,
            initial_trust_best_max=initial_trust_best_max,
            nonlinear_update=args.update, anderson_depth=args.anderson_depth,
            trust_merit=args.trust_merit,
            **poisson_options,
            **solve_options)
    else:
        result = solve_boundary_state_charging(
            solid, conductors, boundary, trust_region=args.trust,
            nonlinear_update=args.update, **solve_options)
except AdaptiveQuadratureConvergenceError as error:
    print("quadrature_failure", str(error))
    print("iteration", error.iteration, "species", error.species)
    print("voltage_range", error.surface_voltage[solid].min(), error.surface_voltage[solid].max())
    q = error.quadrature
    for method in ("adjoint", "forward"):
        data = q[method]
        print(method, "max_stderr", data.element_stderr.max(), "max_level", data.log2_samples.max())
    accepted = error.accepted_state
    restart = {}
    if accepted is not None:
        restart = dict(
            surface_voltage=accepted["surface_voltage"],
            boundary_nodal_voltage=accepted["boundary_nodal_voltage"],
            surface_charge_node_c_per_m=accepted["surface_charge_node_c_per_m"],
            accepted_iterations_total=accepted["accepted_iterations_total"],
            restart_accepted_iterations=accepted["restart_accepted_iterations"],
            restart_beta=accepted["beta_current"],
            accepted_raw_max=accepted["raw_max_abs_log_ratio"],
            accepted_raw_rms=accepted["raw_rms_log_ratio"],
            accepted_confidence_max=(
                accepted["confidence_envelope_max_abs_log_ratio"]),
            accepted_confidence_rms=(
                accepted["confidence_envelope_rms_log_ratio"]),
            trust_best_rms=accepted["trust_best_rms"],
            trust_best_max=accepted["trust_best_max"],
            anderson_x_history=accepted["anderson_x"],
            anderson_residual_history=accepted["anderson_residual"],
            **{f"adaptive_{name}": value
               for name, value in accepted["adaptive_levels"].items()},
            **{f"forward_adaptive_{name}": value
               for name, value in accepted["forward_adaptive_levels"].items()},
            **{f"method_hint_{name}": value
               for name, value in accepted["method_hint"].items()})
    rejected = error.rejected_state
    failed = {}
    if rejected is not None:
        failed = dict(
            failed_boundary_nodal_voltage=rejected["boundary_nodal_voltage"],
            failed_surface_charge_node_c_per_m=(
                rejected["surface_charge_node_c_per_m"]),
            **{f"failed_adaptive_{name}": value
               for name, value in rejected["adaptive_levels"].items()},
            **{f"failed_forward_adaptive_{name}": value
               for name, value in rejected["forward_adaptive_levels"].items()},
            **{f"failed_method_hint_{name}": value
               for name, value in rejected["method_hint"].items()})
    np.savez(
        args.output, status="quadrature_failure", iteration=error.iteration,
        species=error.species, solid=solid,
        failed_surface_voltage=error.surface_voltage,
        failed_potential=error.potential, cells=np.asarray(error.cells), normals=np.asarray(error.normals),
        forward_cell_mean=q["forward_cell_mean"],
        forward_cell_stderr=q["forward_cell_stderr"],
        adjoint_cell_mean=q["adjoint_cell_mean"],
        adjoint_cell_stderr=q["adjoint_cell_stderr"],
        discrepancy_sigma=q["estimator_discrepancy_sigma"],
        method=q["method"], method_within_tolerance=q["method_within_tolerance"],
        estimator_consistent=q["estimator_consistent"], cell_converged=q["cell_converged"],
        consistency_sigma=q["consistency_sigma"], **restart, **failed)
    raise

print("seconds", time.perf_counter() - start)
print("iterations", result["iterations"], "rejected", result.get("rejected_steps", 0))
print("converged", result.get("converged", False), result.get("termination_reason"))
print("beta_final", result.get("beta_final", args.beta))
print("raw_max", [round(x["max_abs_log_ratio"], 4) for x in result["balance_history"]])
print("raw_rms", [round(x["rms_log_ratio"], 4) for x in result["balance_history"]])
print("certified_update_max", [
    round(x["max_abs_log_ratio"], 4) for x in result["interval_balance_history"]])
print("certified_update_rms", [
    round(x["rms_log_ratio"], 4) for x in result["interval_balance_history"]])
print("confidence_max", [
    round(x["confidence_envelope_max_abs_log_ratio"], 4)
    for x in result["interval_balance_history"]])
print("confidence_rms", [
    round(x["confidence_envelope_rms_log_ratio"], 4)
    for x in result["interval_balance_history"]])
final = result["interval_balance_final"]
active = final["active"]
residual = np.abs(final["log_ratio"])
worst = int(np.nanargmax(np.where(active, residual, np.nan)))
if args.nodal:
    if worst < len(result["dielectric_nodes"]):
        node = tuple(result["dielectric_nodes"][worst])
        print("worst_node", node, final["log_ratio"][worst],
              "voltage", result["boundary_nodal_voltage"][node])
    else:
        print("worst_conductor", worst, final["log_ratio"][worst])
    print("currents", result["ion_current"][worst], result["electron_current"][worst],
          result["ion_current_stderr"][worst], result["electron_current_stderr"][worst])
    if args.perturb_node_i is not None:
        probe_node = np.array([args.perturb_node_i, args.perturb_node_j])
        probe_matches = np.where(np.all(result["dielectric_nodes"] == probe_node, axis=1))[0]
        if probe_matches.size != 1:
            raise ValueError("perturbed node is not one dielectric degree of freedom")
        probe = int(probe_matches[0])
        print("probe_node", tuple(probe_node), "residual", final["log_ratio"][probe],
              "voltage", result["boundary_nodal_voltage"][tuple(probe_node)])
        print("probe_currents", result["ion_current"][probe], result["electron_current"][probe],
              result["ion_current_stderr"][probe], result["electron_current_stderr"][probe])
else:
    cell = tuple(result["cells"][worst])
    print("worst", cell, tuple(result["normals"][worst]), final["log_ratio"][worst])
    print("detail", final["detail"].get(("cell", cell)))
    print("voltage", result["surface_voltage"][cell])
print("field", result["field_final"])
for species, q in result["quadrature_history"][-1].items():
    methods, counts = np.unique(q["method"], return_counts=True)
    print(species, dict(zip(methods.tolist(), counts.tolist())))
    discrepancy_index = int(np.argmax(q["estimator_discrepancy_sigma"]))
    print("max_discrepancy", tuple(q["unique_cells"][discrepancy_index]),
          q["estimator_discrepancy_sigma"][discrepancy_index],
          "forward", q["forward_cell_mean"][discrepancy_index],
          q["forward_cell_stderr"][discrepancy_index],
          "adjoint", q["adjoint_cell_mean"][discrepancy_index],
          q["adjoint_cell_stderr"][discrepancy_index])
np.savez(
    args.output, status=("success" if result.get("converged", False) else "iteration_limit"),
    converged=result.get("converged", False), solid=solid,
    surface_voltage=result["surface_voltage"], potential=result["potential"],
    cells=result["cells"], normals=result["normals"],
    raw_max=np.asarray([x["max_abs_log_ratio"] for x in result["balance_history"]]),
    raw_rms=np.asarray([x["rms_log_ratio"] for x in result["balance_history"]]),
    interval_max=np.asarray([x["max_abs_log_ratio"] for x in result["interval_balance_history"]]),
    interval_rms=np.asarray([x["rms_log_ratio"] for x in result["interval_balance_history"]]),
    mean_log_ratio_history=np.stack([
        x["log_ratio"] for x in result["balance_history"]]),
    confidence_separated_log_ratio_history=np.stack([
        x["log_ratio"] for x in result["interval_balance_history"]]),
    active_history=np.stack([
        x["active"] for x in result["interval_balance_history"]]),
    confidence_envelope_max=np.asarray([
        x["confidence_envelope_max_abs_log_ratio"]
        for x in result["interval_balance_history"]]),
    confidence_envelope_rms=np.asarray([
        x["confidence_envelope_rms_log_ratio"]
        for x in result["interval_balance_history"]]),
    **({"boundary_nodal_voltage": result["boundary_nodal_voltage"]}
       if "boundary_nodal_voltage" in result else {}),
    **({"surface_charge_node_c_per_m": result["surface_charge_node_c_per_m"]}
       if "surface_charge_node_c_per_m" in result else {}),
    **({"dielectric_nodes": result["dielectric_nodes"],
        "ion_current": result["ion_current"],
        "electron_current": result["electron_current"],
        "ion_current_stderr": result["ion_current_stderr"],
        "electron_current_stderr": result["electron_current_stderr"]}
       if "dielectric_nodes" in result else {}),
    **{
        f"quadrature_{species_name}_{field_name}": quadrature[field_name]
        for species_name, quadrature in result.get("quadrature_final", {}).items()
        for field_name in (
            "unique_cells", "method", "forward_cell_mean", "forward_cell_stderr",
            "adjoint_cell_mean", "adjoint_cell_stderr", "estimator_discrepancy_sigma",
            "method_within_tolerance", "estimator_consistent", "cell_converged")
    },
    accepted_iterations_total=result.get("accepted_iterations_total", result["iterations"]),
    restart_accepted_iterations=result.get(
        "restart_accepted_iterations", max(result.get("accepted_iterations_total", 1) - 1, 0)),
    restart_beta=result.get("restart_beta", result.get("beta_final", args.beta)),
    trust_best_rms=result.get("trust_best_rms", np.inf),
    trust_best_max=result.get("trust_best_max", np.inf),
    anderson_x_history=result.get("anderson_x_history", np.empty((0, 0))),
    anderson_residual_history=result.get("anderson_residual_history", np.empty((0, 0))),
    **{f"adaptive_{name}": value for name, value in result.get("adaptive_levels", {}).items()},
    **{f"forward_adaptive_{name}": value
       for name, value in result.get("forward_adaptive_levels", {}).items()},
    **{f"method_hint_{name}": value for name, value in result.get("method_hint", {}).items()})
