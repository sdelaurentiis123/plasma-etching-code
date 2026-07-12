"""Compatible boundary-node dielectric charging fixed point.

This is the promotion candidate for the unified charging engine. It deliberately lives beside the
cell-voltage solver until convergence, reciprocity, and invariance gates close.
"""
from __future__ import annotations

import numpy as np

from .boundary_transport import (
    adaptive_adjoint_boundary_state_face_flux,
    adjoint_boundary_state_face_flux,
    bidirectional_boundary_state_cell_flux,
)
from .charging_backward import AdaptiveQuadratureConvergenceError, _gas_faces
from .charging_nodal import material_face_nodes, nodal_domain, solve_nodal_laplace
from .charging_poisson import NodalPoissonSystem


def _anderson_step(x, residual, x_history, residual_history, gain, depth):
    """Return a type-II Anderson step for one preconditioned fixed-point residual."""
    x = np.asarray(x, dtype=float); residual = np.asarray(residual, dtype=float)
    if x.shape != residual.shape or x.ndim != 1:
        raise ValueError("Anderson state and residual must be matching vectors")
    if int(depth) != depth or depth <= 0 or not np.isfinite(gain) or gain <= 0.0:
        raise ValueError("Anderson depth and gain must be positive")
    x_history.append(x.copy()); residual_history.append(residual.copy())
    if len(x_history) > int(depth) + 1:
        x_history.pop(0); residual_history.pop(0)
    step = float(gain) * residual
    if len(residual_history) >= 2:
        delta_residual = np.stack([
            residual_history[index + 1] - residual_history[index]
            for index in range(len(residual_history) - 1)], axis=1)
        delta_x = np.stack([
            x_history[index + 1] - x_history[index]
            for index in range(len(x_history) - 1)], axis=1)
        gamma, *_ = np.linalg.lstsq(delta_residual, residual, rcond=1e-8)
        step = step - (delta_x + float(gain) * delta_residual) @ gamma
    return step


def _confidence_separated_log_ratio(
        ion_current, electron_current, ion_stderr, electron_stderr, confidence_sigma):
    """Return the current imbalance resolved outside estimator confidence intervals.

    The exact-current limit is ``log(Gi/Ge)``. Overlapping intervals return zero because the sampled
    transport cannot determine the update direction. Tightening the current uncertainty shrinks that
    unresolved band toward the physical mean-current root.
    """
    ion = np.asarray(ion_current, dtype=float)
    electron = np.asarray(electron_current, dtype=float)
    ion_error = np.asarray(ion_stderr, dtype=float)
    electron_error = np.asarray(electron_stderr, dtype=float)
    if (ion.shape != electron.shape or ion.shape != ion_error.shape
            or ion.shape != electron_error.shape):
        raise ValueError("current means and standard errors must have identical shapes")
    if (not np.all(np.isfinite(ion)) or not np.all(np.isfinite(electron))
            or not np.all(np.isfinite(ion_error)) or not np.all(np.isfinite(electron_error))
            or np.any(ion < 0.0) or np.any(electron < 0.0)
            or np.any(ion_error < 0.0) or np.any(electron_error < 0.0)):
        raise ValueError("current means and standard errors must be finite and nonnegative")
    sigma = float(confidence_sigma)
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("confidence_sigma must be finite and positive")
    ion_low = np.maximum(ion - sigma * ion_error, 0.0)
    ion_high = ion + sigma * ion_error
    electron_low = np.maximum(electron - sigma * electron_error, 0.0)
    electron_high = electron + sigma * electron_error
    residual = np.zeros_like(ion)
    ion_dominant = ion_low > electron_high
    electron_dominant = electron_low > ion_high
    residual[ion_dominant] = np.log(
        ion_low[ion_dominant] / np.maximum(electron_high[ion_dominant], 1e-300))
    residual[electron_dominant] = np.log(
        np.maximum(ion_high[electron_dominant], 1e-300)
        / electron_low[electron_dominant])
    return residual, ion_low, ion_high, electron_low, electron_high


def _trust_merit_worsened(current_rms, current_max, previous_rms, previous_max, tolerance, mode):
    """Return whether a trial violates the selected physical mean-current trust merit."""
    rms_worse = current_rms > previous_rms * (1.0 + tolerance)
    max_worse = current_max > previous_max * (1.0 + tolerance)
    if mode == "rms":
        return rms_worse
    if mode == "max":
        return max_worse
    if mode == "pareto":
        return rms_worse or max_worse
    raise ValueError("trust_merit must be 'rms', 'max', or 'pareto'")


def _trust_merit_strongly_improved(
        current_rms, current_max, previous_rms, previous_max, mode, factor=0.8):
    if mode == "rms":
        return current_rms < factor * previous_rms
    if mode == "max":
        return current_max < factor * previous_max
    if mode == "pareto":
        return (current_rms < factor * previous_rms
                and current_max < factor * previous_max)
    raise ValueError("trust_merit must be 'rms', 'max', or 'pareto'")


def solve_boundary_state_charging_nodal(
        solid, conductor_ids, boundary_state, *, ion_species=None, electron_species=None,
        initial_surface_voltage=None, initial_boundary_nodal_voltage=None,
        n_iter=40, beta=0.5, response_energy_eV=4.0,
        dVmax=8.0, balance_tol=1e-3, min_iter=2, field_sweeps=500,
        field_tolerance=1e-9, boundary_proposals=None, n_face_position=8,
        adaptive_quadrature=None, active_flux=1e-4, current_confidence_sigma=2.0,
        trust_region=True, trust_growth_tolerance=0.02, minimum_beta=1e-4,
        trust_merit="rms",
        gain_decay=0.0, gain_offset=5.0, epsilon_r=None, cell_size_m=None,
        grounded_nodes=None, initial_surface_charge_node_c_per_m=None,
        initial_adaptive_levels=None, initial_forward_adaptive_levels=None,
        initial_method_hint=None, initial_accepted_iterations=0,
        initial_beta=None, initial_anderson_x=None, initial_anderson_residual=None,
        initial_trust_best_rms=None, initial_trust_best_max=None,
        nonlinear_update="picard", anderson_depth=4):
    """Solve steady floating-current balance on physical material-boundary vertices.

    Dielectric unknowns and the electrostatic field share the same nodal basis. Particle histories and
    adjoint face quadrature deposit current through the same linear edge shape functions; connected
    conductors retain one pooled equipotential. The iteration is a nonlinear root finder, not physical
    charging time.
    """
    solid = np.asarray(solid, dtype=bool)
    conductor_ids = np.asarray(conductor_ids, dtype=int)
    if (conductor_ids.shape != solid.shape or np.any(conductor_ids < 0)
            or np.any((conductor_ids > 0) & ~solid)):
        raise ValueError("conductor_ids must be nonnegative, match solid, and label only material")
    if np.any(solid[:, 0]):
        raise ValueError("plasma reference plane must be a gas-only top row above all material")
    if n_iter <= 0 or min_iter <= 0 or beta <= 0.0 or dVmax <= 0.0:
        raise ValueError("iteration counts, beta, and dVmax must be positive")
    if not 0.0 <= gain_decay <= 1.0 or gain_offset <= 0.0:
        raise ValueError("gain_decay must lie in [0,1] and gain_offset must be positive")
    if int(initial_accepted_iterations) != initial_accepted_iterations or initial_accepted_iterations < 0:
        raise ValueError("initial_accepted_iterations must be a nonnegative integer")
    if nonlinear_update not in {"picard", "anderson"}:
        raise ValueError("nonlinear_update must be 'picard' or 'anderson'")
    if int(anderson_depth) != anderson_depth or anderson_depth <= 0:
        raise ValueError("anderson_depth must be a positive integer")
    if trust_merit not in {"rms", "max", "pareto"}:
        raise ValueError("trust_merit must be 'rms', 'max', or 'pareto'")
    if initial_beta is not None and (
            not np.isfinite(initial_beta) or initial_beta <= 0.0 or initial_beta > beta):
        raise ValueError("initial_beta must be positive and no larger than beta")
    for name, value in (("initial_trust_best_rms", initial_trust_best_rms),
                        ("initial_trust_best_max", initial_trust_best_max)):
        if value is not None and (not np.isfinite(value) or value < 0.0):
            raise ValueError(f"{name} must be finite and nonnegative")
    charge_mode = epsilon_r is not None
    if charge_mode:
        epsilon_r = np.asarray(epsilon_r, dtype=float)
        if (epsilon_r.shape != solid.shape or not np.all(np.isfinite(epsilon_r))
                or np.any(epsilon_r <= 0.0)):
            raise ValueError("epsilon_r must be a finite positive cell grid matching solid")
        if cell_size_m is None or not np.isfinite(cell_size_m) or cell_size_m <= 0.0:
            raise ValueError("charge-space Poisson mode requires positive cell_size_m")
    elif grounded_nodes is not None or initial_surface_charge_node_c_per_m is not None:
        raise ValueError("grounded nodes and physical surface charge require epsilon_r")

    def select(selection, positive):
        if selection is None:
            items = [item for item in boundary_state.species
                     if item.charge_number != 0 and (item.charge_number > 0) == positive]
        else:
            names = [selection] if isinstance(selection, str) else list(selection)
            items = [boundary_state.get(name) for name in names]
        sign = 1 if positive else -1
        if (not items or any(np.sign(item.charge_number) != sign for item in items)
                or any(item.density_model is None for item in items)):
            raise ValueError("charging species require signed continuous boundary densities")
        return items

    positive_species = select(ion_species, True)
    negative_species = select(electron_species, False)
    current_scale = max(
        sum(item.flux_m2_s * abs(item.charge_number) for item in positive_species),
        sum(item.flux_m2_s * abs(item.charge_number) for item in negative_species), 1e-300)
    proposals = {} if boundary_proposals is None else dict(boundary_proposals)
    cells, normals = _gas_faces(solid, solid)
    if not cells:
        raise ValueError("solid grid has no gas-facing material surface")
    normals_array = np.asarray(normals, dtype=float)
    face_components = np.asarray([conductor_ids[cell] for cell in cells], dtype=int)
    face_nodes = [material_face_nodes(cell, normal) for cell, normal in zip(cells, normals)]

    initial_cell = (np.zeros(solid.shape) if initial_surface_voltage is None
                    else np.asarray(initial_surface_voltage, dtype=float).copy())
    if initial_cell.shape != solid.shape or not np.all(np.isfinite(initial_cell)):
        raise ValueError("initial_surface_voltage must be a finite grid matching solid")
    if initial_boundary_nodal_voltage is None:
        _, _, boundary_voltage = nodal_domain(solid, initial_cell)
    else:
        boundary_voltage = np.asarray(initial_boundary_nodal_voltage, dtype=float).copy()
        expected = (solid.shape[0] + 1, solid.shape[1] + 1)
        if boundary_voltage.shape != expected or not np.all(np.isfinite(boundary_voltage)):
            raise ValueError("initial_boundary_nodal_voltage must be a finite nodal grid")
        boundary_voltage[:, 0] = 0.0

    node_component = {}
    for component, endpoints in zip(face_components, face_nodes):
        if component <= 0:
            continue
        for node in endpoints:
            previous = node_component.get(node, int(component))
            if previous != int(component):
                raise ValueError("distinct conductor ids meet at one electrical boundary node")
            node_component[node] = int(component)
    boundary_nodes = sorted({node for endpoints in face_nodes for node in endpoints if node[1] != 0})
    dielectric_nodes = [node for node in boundary_nodes if node_component.get(node, 0) == 0]
    dielectric_index = {node: index for index, node in enumerate(dielectric_nodes)}
    components = sorted(set(int(value) for value in face_components if value > 0))
    component_index = {
        component: len(dielectric_nodes) + index for index, component in enumerate(components)}
    dof_count = len(dielectric_nodes) + len(components)
    if dof_count == 0:
        raise ValueError("charging surface has no floating electrical degrees of freedom")
    if (initial_anderson_x is None) != (initial_anderson_residual is None):
        raise ValueError("both Anderson restart histories are required together")
    if initial_anderson_x is None:
        initial_anderson_x_array = np.empty((0, dof_count))
        initial_anderson_residual_array = np.empty((0, dof_count))
    else:
        initial_anderson_x_array = np.asarray(initial_anderson_x, dtype=float)
        initial_anderson_residual_array = np.asarray(initial_anderson_residual, dtype=float)
        if (initial_anderson_x_array.ndim != 2
                or initial_anderson_x_array.shape[1:] != (dof_count,)
                or initial_anderson_residual_array.shape != initial_anderson_x_array.shape
                or initial_anderson_x_array.shape[0] > int(anderson_depth) + 1
                or not np.all(np.isfinite(initial_anderson_x_array))
                or not np.all(np.isfinite(initial_anderson_residual_array))):
            raise ValueError("invalid Anderson restart histories")
    conductor_voltage = np.zeros(int(conductor_ids.max()) + 1)
    for component in components:
        values = (np.asarray([boundary_voltage[node] for node, owner in node_component.items()
                             if owner == component]) if initial_boundary_nodal_voltage is not None
                  else initial_cell[conductor_ids == component])
        conductor_voltage[component] = float(values.mean()) if values.size else 0.0

    poisson_system = None
    surface_charge_node = np.zeros_like(boundary_voltage)
    surface_capacitance = np.zeros(dof_count)
    node_surface_length_m = np.zeros_like(boundary_voltage)
    if charge_mode:
        fixed = np.zeros_like(boundary_voltage, dtype=bool); fixed[:, 0] = True
        fixed_voltage = np.zeros_like(boundary_voltage)
        for node in node_component:
            fixed[node] = True
        if grounded_nodes is not None:
            grounded = np.asarray(grounded_nodes, dtype=bool)
            if grounded.shape != boundary_voltage.shape:
                raise ValueError("grounded_nodes must match the nodal grid")
            if any(grounded[node] for node in node_component):
                raise ValueError("a floating-conductor node cannot also be grounded")
            fixed |= grounded
        poisson_system = NodalPoissonSystem(epsilon_r, fixed, fixed_voltage)
        if initial_surface_charge_node_c_per_m is not None:
            initial_charge = np.asarray(initial_surface_charge_node_c_per_m, dtype=float)
            if initial_charge.shape != boundary_voltage.shape or not np.all(np.isfinite(initial_charge)):
                raise ValueError("initial surface charge must be a finite nodal grid")
            surface_charge_node[:] = initial_charge
        dielectric_array = np.asarray(dielectric_nodes, dtype=int)
        if dielectric_array.size:
            surface_capacitance[:len(dielectric_nodes)] = (
                poisson_system.diagonal_surface_capacitance(dielectric_array))
        for endpoints in face_nodes:
            for node in endpoints:
                if node in dielectric_index:
                    node_surface_length_m[node] += 0.5 * float(cell_size_m)

    def impose_conductors():
        for node, component in node_component.items():
            if poisson_system is None:
                boundary_voltage[node] = conductor_voltage[component]
            else:
                poisson_system.dirichlet_voltage[node] = conductor_voltage[component]
        if poisson_system is None:
            boundary_voltage[:, 0] = 0.0
        else:
            poisson_system.dirichlet_voltage[:, 0] = 0.0

    by_cell = {}
    for cell, endpoints in zip(cells, face_nodes):
        by_cell.setdefault(cell, []).extend(endpoints)

    def surface_readout():
        result = np.zeros(solid.shape)
        for cell, nodes in by_cell.items():
            result[cell] = float(np.mean([boundary_voltage[node] for node in nodes]))
        for component in components:
            result[conductor_ids == component] = conductor_voltage[component]
        return result

    hybrid_hint = ({} if initial_method_hint is None else {
        name: np.asarray(value).copy() for name, value in initial_method_hint.items()})
    adaptive_levels = ({} if initial_adaptive_levels is None else {
        name: np.asarray(value, dtype=int).copy()
        for name, value in initial_adaptive_levels.items()})
    forward_adaptive_levels = ({} if initial_forward_adaptive_levels is None else {
        name: np.asarray(value, dtype=int).copy()
        for name, value in initial_forward_adaptive_levels.items()})
    history = []; interval_history = []; field_history = []; quadrature_history = []
    species_face_current = {}; species_face_stderr = {}; species_face_replicates = {}
    species_endpoint_stderr = {}; species_endpoint_replicates = {}
    beta_current = float(beta if initial_beta is None else initial_beta)
    pending_step = None; rejected_steps = 0
    anderson_x = [row.copy() for row in initial_anderson_x_array]
    anderson_residual = [row.copy() for row in initial_anderson_residual_array]
    last_accepted_state = None
    trial_merit_history = []; trial_max_history = []
    accepted_beta_history = []; accepted_gain_history = []
    trust_best_rms = (float("inf") if initial_trust_best_rms is None
                      else float(initial_trust_best_rms))
    trust_best_max = (float("inf") if initial_trust_best_max is None
                      else float(initial_trust_best_max))

    def accepted_checkpoint():
        """Return only the serializable state needed to restart the last accepted iterate."""
        if last_accepted_state is None:
            return None
        return dict(
            solid=solid.copy(),
            surface_voltage=last_accepted_state["surface_voltage"].copy(),
            boundary_nodal_voltage=last_accepted_state["boundary_voltage"].copy(),
            surface_charge_node_c_per_m=last_accepted_state["surface_charge_node"].copy(),
            adaptive_levels={
                name: value.copy()
                for name, value in last_accepted_state["adaptive_levels"].items()},
            forward_adaptive_levels={
                name: value.copy()
                for name, value in last_accepted_state["forward_adaptive_levels"].items()},
            method_hint={
                name: value.copy() for name, value in last_accepted_state["method_hint"].items()},
            beta_current=float(last_accepted_state["beta_current"]),
            anderson_x=np.asarray(last_accepted_state["anderson_x"]).copy(),
            anderson_residual=np.asarray(
                last_accepted_state["anderson_residual"]).copy(),
            raw_max_abs_log_ratio=float(last_accepted_state["raw_max_abs_log_ratio"]),
            raw_rms_log_ratio=float(last_accepted_state["raw_rms_log_ratio"]),
            confidence_envelope_max_abs_log_ratio=float(
                last_accepted_state["confidence_envelope_max_abs_log_ratio"]),
            confidence_envelope_rms_log_ratio=float(
                last_accepted_state["confidence_envelope_rms_log_ratio"]),
            trust_best_rms=float(last_accepted_state["trust_best_rms"]),
            trust_best_max=float(last_accepted_state["trust_best_max"]),
            accepted_iterations_total=int(
                last_accepted_state["accepted_iterations_total"]),
            restart_accepted_iterations=max(
                int(last_accepted_state["accepted_iterations_total"]) - 1, 0))

    def rejected_quadrature_checkpoint(potential):
        """Capture the failed trial solely so a new quadrature epoch can certify its neighborhood."""
        return dict(
            solid=solid.copy(), surface_voltage=surface_readout(),
            boundary_nodal_voltage=np.asarray(potential).copy(),
            surface_charge_node_c_per_m=surface_charge_node.copy(),
            adaptive_levels={name: value.copy() for name, value in adaptive_levels.items()},
            forward_adaptive_levels={
                name: value.copy() for name, value in forward_adaptive_levels.items()},
            method_hint={name: value.copy() for name, value in hybrid_hint.items()})

    for iteration in range(int(n_iter)):
        impose_conductors()
        if poisson_system is None:
            potential, field_diag = solve_nodal_laplace(
                solid, boundary_nodal_voltage=boundary_voltage,
                sweeps=field_sweeps, omega=1.7, tolerance=field_tolerance)
        else:
            potential, poisson_diag = poisson_system.solve(surface_charge_node)
            boundary_voltage[:] = potential
            field_diag = dict(
                sweeps=1, max_abs=poisson_diag.max_abs_residual_v,
                rms=poisson_diag.rms_residual_v,
                active_nodes=int(np.prod(potential.shape)),
                free_nodes=poisson_diag.free_nodes,
                electrostatic_energy_j_per_m=poisson_diag.electrostatic_energy_j_per_m,
                specified_charge_c_per_m=poisson_diag.specified_charge_c_per_m,
                dirichlet_reaction_charge_c_per_m=(
                    poisson_diag.dirichlet_reaction_charge_c_per_m),
                charge_balance_c_per_m=poisson_diag.charge_balance_c_per_m)
        field_history.append(field_diag); species_quadrature = {}
        for species in positive_species + negative_species:
            if adaptive_quadrature is None:
                estimate = adjoint_boundary_state_face_flux(
                    boundary_state, species.name, potential, solid, cells, normals_array,
                    proposal_species=proposals.get(species.name), n_face_position=n_face_position)
                normalized = estimate["per_face"]
                normalized_stderr = np.zeros_like(normalized)
                normalized_replicates = normalized[None, :]
                endpoint_stderr = np.zeros((len(cells), 2))
                endpoint_replicates = estimate["per_face_endpoint"][None, :, :]
                species_quadrature[species.name] = estimate
            else:
                options = dict(adaptive_quadrature)
                bidirectional = bool(options.pop("bidirectional", False))
                forward_options = dict(options.pop("forward_options", {}))
                switch_factor = float(options.pop("method_switch_factor", 2.0))
                consistency_sigma = float(options.pop("consistency_sigma", 5.0))
                support_sigma = float(options.pop("support_sigma", 2.0))
                support_ratio = float(options.pop("support_ratio", 0.5))
                freeze_method_hint = bool(options.pop("freeze_method_hint", False))
                warm_start_backoff = int(options.pop("warm_start_backoff", 0))
                if warm_start_backoff < 0:
                    raise ValueError("warm_start_backoff must be nonnegative")
                options.setdefault("n_face_position", n_face_position)
                if species.name in adaptive_levels:
                    base_level = int(options.get("base_log2", 6))
                    options.setdefault(
                        "initial_log2_samples",
                        np.maximum(adaptive_levels[species.name] - warm_start_backoff,
                                   base_level))
                if bidirectional:
                    if species.name in forward_adaptive_levels:
                        forward_base = int(forward_options.get("base_log2", 8))
                        forward_options.setdefault(
                            "initial_log2_samples",
                            np.maximum(
                                forward_adaptive_levels[species.name] - warm_start_backoff,
                                forward_base))
                    hybrid = bidirectional_boundary_state_cell_flux(
                        boundary_state, species.name, potential, solid, cells, normals_array,
                        proposal_species=proposals.get(species.name),
                        adjoint_options=options, forward_options=forward_options,
                        element_absolute_tolerance=options.get("element_absolute_tolerance", 1e-3),
                        element_relative_tolerance=options.get("element_relative_tolerance", 0.05),
                        method_hint=hybrid_hint.get(species.name), switch_factor=switch_factor,
                        consistency_sigma=consistency_sigma, support_sigma=support_sigma,
                        support_ratio=support_ratio,
                        freeze_method_hint=freeze_method_hint)
                    hybrid_hint[species.name] = hybrid["method"].copy()
                    adaptive_levels[species.name] = hybrid["adjoint"].log2_samples.copy()
                    forward_adaptive_levels[species.name] = (
                        hybrid["forward"].log2_samples.copy())
                    species_quadrature[species.name] = hybrid
                    if not hybrid["converged"]:
                        failed = np.where(~hybrid["cell_converged"])[0]
                        worst = int(failed[np.argmax(
                            hybrid["estimator_discrepancy_sigma"][failed])])
                        raise AdaptiveQuadratureConvergenceError(
                            f"nodal charging quadrature failed for {species.name!r}, cell="
                            f"{tuple(hybrid['unique_cells'][worst])}",
                            iteration=iteration + 1, species=species.name, quadrature=hybrid,
                            surface_voltage=surface_readout(), potential=potential,
                            cells=cells, normals=normals,
                            accepted_state=accepted_checkpoint(),
                            rejected_state=rejected_quadrature_checkpoint(potential))
                    normalized = hybrid["selected_face_mean"]
                    normalized_stderr = hybrid["selected_face_stderr"]
                    normalized_replicates = hybrid["selected_face_replicates"]
                    endpoint_stderr = hybrid["selected_endpoint_stderr"]
                    endpoint_replicates = hybrid["selected_endpoint_replicates"]
                else:
                    estimate = adaptive_adjoint_boundary_state_face_flux(
                        boundary_state, species.name, potential, solid, cells, normals_array,
                        proposal_species=proposals.get(species.name), **options)
                    species_quadrature[species.name] = estimate
                    if not estimate.converged:
                        raise AdaptiveQuadratureConvergenceError(
                            f"nodal charging adjoint quadrature failed for {species.name!r}",
                            iteration=iteration + 1, species=species.name,
                            quadrature=estimate, surface_voltage=surface_readout(),
                            potential=potential, cells=cells, normals=normals,
                            accepted_state=accepted_checkpoint(),
                            rejected_state=rejected_quadrature_checkpoint(potential))
                    normalized = estimate.element_mean
                    normalized_stderr = estimate.element_stderr
                    normalized_replicates = estimate.element_replicates
                    endpoint_replicates = estimate.auxiliary_replicates
                    endpoint_stderr = endpoint_replicates.std(axis=0, ddof=1) / np.sqrt(
                        endpoint_replicates.shape[0])
            scale = species.flux_m2_s * abs(species.charge_number)
            species_face_current[species.name] = normalized * scale
            species_face_stderr[species.name] = normalized_stderr * scale
            species_face_replicates[species.name] = normalized_replicates * scale
            species_endpoint_stderr[species.name] = endpoint_stderr * scale
            species_endpoint_replicates[species.name] = endpoint_replicates * scale
        quadrature_history.append(species_quadrature)

        def assemble(items):
            replicate_counts = {species_face_replicates[item.name].shape[0] for item in items}
            if len(replicate_counts) != 1:
                raise ValueError("charging species require a common replicate count")
            replicate_count = replicate_counts.pop()
            replicates = np.zeros((replicate_count, dof_count))
            floor_variance = np.zeros(dof_count)
            for species in items:
                face_replicates = species_face_replicates[species.name]
                endpoint_error = species_endpoint_stderr[species.name]
                endpoint_replicates = species_endpoint_replicates[species.name]
                for face_index, (component, endpoints) in enumerate(zip(face_components, face_nodes)):
                    if component > 0:
                        index = component_index[int(component)]
                        replicates[:, index] += face_replicates[:, face_index]
                        floor_variance[index] += np.sum(endpoint_error[face_index] ** 2)
                    else:
                        for endpoint_index, node in enumerate(endpoints):
                            node_comp = node_component.get(node, 0)
                            index = (component_index[node_comp] if node_comp > 0
                                     else dielectric_index[node])
                            replicates[:, index] += endpoint_replicates[:, face_index, endpoint_index]
                            floor_variance[index] += endpoint_error[face_index, endpoint_index] ** 2
            value = replicates.mean(axis=0)
            replicate_stderr = (replicates.std(axis=0, ddof=1) / np.sqrt(replicate_count)
                                if replicate_count > 1 else np.zeros(dof_count))
            stderr = np.maximum(replicate_stderr, np.sqrt(floor_variance))
            return value / current_scale, stderr / current_scale

        ion_current, ion_stderr = assemble(positive_species)
        electron_current, electron_stderr = assemble(negative_species)
        raw = np.log((ion_current + 1e-12) / (electron_current + 1e-12))
        certified, ilo, ihi, elo, ehi = _confidence_separated_log_ratio(
            ion_current, electron_current, ion_stderr, electron_stderr,
            current_confidence_sigma)
        # A zero-hit estimate can have zero mean but finite support. Use the confidence upper bounds
        # for activity so an unresolved rare current is not silently declared physically absent.
        activity = (ihi + ehi) >= active_flux
        raw_active = np.abs(raw[activity]); cert_active = np.abs(certified[activity])
        balance = dict(
            log_ratio=raw.copy(), active=activity.copy(),
            max_abs_log_ratio=float(raw_active.max()) if raw_active.size else 0.0,
            rms_log_ratio=float(np.sqrt(np.mean(raw_active ** 2))) if raw_active.size else 0.0)
        interval_balance = dict(
            log_ratio=certified.copy(), active=activity.copy(),
            max_abs_log_ratio=float(cert_active.max()) if cert_active.size else 0.0,
            rms_log_ratio=float(np.sqrt(np.mean(cert_active ** 2))) if cert_active.size else 0.0)
        log_lower = np.full(dof_count, -np.inf)
        log_upper = np.full(dof_count, np.inf)
        positive_ion_lower = ilo > 0.0
        positive_electron_lower = elo > 0.0
        log_lower[positive_ion_lower] = np.log(
            ilo[positive_ion_lower] / np.maximum(ehi[positive_ion_lower], 1e-300))
        log_upper[positive_electron_lower] = np.log(
            np.maximum(ihi[positive_electron_lower], 1e-300)
            / elo[positive_electron_lower])
        confidence_envelope = np.maximum(
            np.abs(log_lower[activity]), np.abs(log_upper[activity]))
        interval_balance.update(
            log_ratio_interval_lower=log_lower,
            log_ratio_interval_upper=log_upper,
            confidence_envelope_max_abs_log_ratio=(
                float(confidence_envelope.max()) if confidence_envelope.size else 0.0),
            confidence_envelope_rms_log_ratio=(
                float(np.sqrt(np.mean(confidence_envelope ** 2)))
                if confidence_envelope.size else 0.0),
            confidence_intervals_finite=bool(np.all(np.isfinite(confidence_envelope))))
        # Trust-region acceptance follows the physical mean-current equation under common deterministic
        # samples. The confidence-separated distance can shrink merely because uncertainty widened and
        # therefore cannot be used as an optimization merit.
        merit = balance["rms_log_ratio"]
        maximum_merit = balance["max_abs_log_ratio"]
        trial_merit_history.append(merit); trial_max_history.append(maximum_merit)
        # The tolerance is a fixed numerical-noise tube around the best accepted component values,
        # not a per-step allowance that may ratchet upward over a long continuation.
        if (trust_region and pending_step is not None
                and _trust_merit_worsened(
                    merit, maximum_merit, trust_best_rms, trust_best_max,
                    trust_growth_tolerance, trust_merit)):
            boundary_voltage[:] = pending_step["boundary_voltage"]
            conductor_voltage[:] = pending_step["conductor_voltage"]
            surface_charge_node[:] = pending_step["surface_charge_node"]
            hybrid_hint = {
                name: value.copy() for name, value in pending_step["hybrid_hint"].items()}
            adaptive_levels = {
                name: value.copy() for name, value in pending_step["adaptive_levels"].items()}
            forward_adaptive_levels = {
                name: value.copy()
                for name, value in pending_step["forward_adaptive_levels"].items()}
            beta_current *= 0.5; rejected_steps += 1; pending_step = None
            anderson_x.clear(); anderson_residual.clear()
            if beta_current < minimum_beta:
                raise RuntimeError(
                    f"nodal charging trust region collapsed below minimum_beta={minimum_beta:g}")
            continue
        if (trust_region and pending_step is not None
                and _trust_merit_strongly_improved(
                    merit, maximum_merit, pending_step["merit"],
                    pending_step["maximum_merit"], trust_merit)):
            beta_current = min(float(beta), 1.2 * beta_current)
        pending_step = None
        trust_best_rms = min(trust_best_rms, merit)
        trust_best_max = min(trust_best_max, maximum_merit)
        history.append(balance); interval_history.append(interval_balance)
        accepted_beta_history.append(beta_current)
        gain_iteration = int(initial_accepted_iterations) + len(history)
        iteration_gain = beta_current * (1.0 + gain_iteration / gain_offset) ** (-gain_decay)
        accepted_gain_history.append(iteration_gain)
        last_accepted_state = dict(
            boundary_voltage=boundary_voltage.copy(),
            conductor_voltage=conductor_voltage.copy(), potential=potential.copy(),
            surface_charge_node=surface_charge_node.copy(),
            surface_voltage=surface_readout(),
            ion_current=ion_current.copy(), electron_current=electron_current.copy(),
            ion_stderr=ion_stderr.copy(), electron_stderr=electron_stderr.copy(),
            species_face_current={
                name: value.copy() for name, value in species_face_current.items()},
            species_face_stderr={
                name: value.copy() for name, value in species_face_stderr.items()},
            species_face_replicates={
                name: value.copy() for name, value in species_face_replicates.items()},
            species_endpoint_stderr={
                name: value.copy() for name, value in species_endpoint_stderr.items()},
            species_endpoint_replicates={
                name: value.copy() for name, value in species_endpoint_replicates.items()},
            quadrature=species_quadrature,
            method_hint={name: value.copy() for name, value in hybrid_hint.items()},
            adaptive_levels={name: value.copy() for name, value in adaptive_levels.items()},
            forward_adaptive_levels={
                name: value.copy() for name, value in forward_adaptive_levels.items()},
            beta_current=float(beta_current),
            raw_max_abs_log_ratio=balance["max_abs_log_ratio"],
            raw_rms_log_ratio=balance["rms_log_ratio"],
            confidence_envelope_max_abs_log_ratio=(
                interval_balance["confidence_envelope_max_abs_log_ratio"]),
            confidence_envelope_rms_log_ratio=(
                interval_balance["confidence_envelope_rms_log_ratio"]),
            trust_best_rms=float(trust_best_rms), trust_best_max=float(trust_best_max),
            anderson_x=(np.stack(anderson_x) if anderson_x
                        else np.empty((0, dof_count))),
            anderson_residual=(np.stack(anderson_residual) if anderson_residual
                               else np.empty((0, dof_count))),
            accepted_iterations_total=int(initial_accepted_iterations) + len(history))
        if (balance_tol is not None and len(history) >= int(min_iter)
                and interval_history[-1]["confidence_envelope_max_abs_log_ratio"] <= balance_tol):
            break
        if trust_region:
            pending_step = dict(
                boundary_voltage=boundary_voltage.copy(),
                conductor_voltage=conductor_voltage.copy(), merit=merit,
                maximum_merit=maximum_merit,
                surface_charge_node=surface_charge_node.copy(),
                hybrid_hint={name: value.copy() for name, value in hybrid_hint.items()},
                adaptive_levels={
                    name: value.copy() for name, value in adaptive_levels.items()},
                forward_adaptive_levels={
                    name: value.copy() for name, value in forward_adaptive_levels.items()})
        # Confidence intervals gate the direction; the residual magnitude remains the physical mean
        # equation. Using an interval edge as the residual would move the fixed point with sampling
        # tolerance. A sample ladder must still narrow the final confidence envelope below the claim.
        direction_resolved = certified != 0.0
        update_residual = np.where(activity & direction_resolved, raw, 0.0)
        preconditioned_state = np.empty(dof_count)
        for index, node in enumerate(dielectric_nodes):
            preconditioned_state[index] = (
                boundary_voltage[node] if poisson_system is None
                else surface_charge_node[node] / surface_capacitance[index])
        for component in components:
            preconditioned_state[component_index[component]] = conductor_voltage[component]
        fixed_point_residual = response_energy_eV * update_residual
        if nonlinear_update == "anderson":
            step = _anderson_step(
                preconditioned_state, fixed_point_residual,
                anderson_x, anderson_residual, iteration_gain, anderson_depth)
        else:
            step = iteration_gain * fixed_point_residual
        step = np.clip(step, -dVmax, dVmax)
        for index, node in enumerate(dielectric_nodes):
            if activity[index]:
                if poisson_system is None:
                    boundary_voltage[node] += step[index]
                else:
                    surface_charge_node[node] += surface_capacitance[index] * step[index]
        for component in components:
            index = component_index[component]
            if activity[index]:
                conductor_voltage[component] += step[index]

    # A final fixed-point proposal has not had its currents evaluated. Return the last accepted state
    # so voltage, potential, residual, and checkpoint restarts all describe one physical iterate.
    if last_accepted_state is None:
        raise RuntimeError("nodal charging solve produced no accepted state")
    converged = bool(
        balance_tol is not None and len(history) >= int(min_iter)
        and interval_history[-1]["confidence_envelope_max_abs_log_ratio"]
        <= float(balance_tol))
    boundary_voltage[:] = last_accepted_state["boundary_voltage"]
    conductor_voltage[:] = last_accepted_state["conductor_voltage"]
    surface_charge_node[:] = last_accepted_state["surface_charge_node"]
    ion_current = last_accepted_state["ion_current"]
    electron_current = last_accepted_state["electron_current"]
    ion_stderr = last_accepted_state["ion_stderr"]
    electron_stderr = last_accepted_state["electron_stderr"]
    species_face_current = last_accepted_state["species_face_current"]
    species_face_stderr = last_accepted_state["species_face_stderr"]
    species_face_replicates = last_accepted_state["species_face_replicates"]
    species_endpoint_stderr = last_accepted_state["species_endpoint_stderr"]
    species_endpoint_replicates = last_accepted_state["species_endpoint_replicates"]
    impose_conductors()
    if poisson_system is None:
        potential, field_final = solve_nodal_laplace(
            solid, boundary_nodal_voltage=boundary_voltage,
            sweeps=field_sweeps, omega=1.7, tolerance=field_tolerance)
    else:
        potential, poisson_diag = poisson_system.solve(surface_charge_node)
        boundary_voltage[:] = potential
        field_final = dict(
            sweeps=1, max_abs=poisson_diag.max_abs_residual_v,
            rms=poisson_diag.rms_residual_v,
            active_nodes=int(np.prod(potential.shape)), free_nodes=poisson_diag.free_nodes,
            electrostatic_energy_j_per_m=poisson_diag.electrostatic_energy_j_per_m,
            specified_charge_c_per_m=poisson_diag.specified_charge_c_per_m,
            dirichlet_reaction_charge_c_per_m=(
                poisson_diag.dirichlet_reaction_charge_c_per_m),
            charge_balance_c_per_m=poisson_diag.charge_balance_c_per_m)
    return dict(
        solid=solid.copy(), conductor_ids=conductor_ids.copy(),
        surface_voltage=surface_readout(), boundary_nodal_voltage=boundary_voltage.copy(),
        surface_charge_node_c_per_m=surface_charge_node.copy(),
        surface_charge_density_c_per_m2=np.divide(
            surface_charge_node, node_surface_length_m, out=np.zeros_like(surface_charge_node),
            where=node_surface_length_m > 0.0),
        node_surface_length_m=node_surface_length_m,
        potential=potential, cells=np.asarray(cells, dtype=int), normals=normals_array,
        face_components=face_components, dielectric_nodes=np.asarray(dielectric_nodes, dtype=int),
        conductor_voltage=conductor_voltage, iterations=len(history),
        converged=converged,
        termination_reason=(
            "balance_tolerance" if converged else
            "fixed_iteration_budget" if balance_tol is None else "iteration_limit"),
        requested_balance_tolerance=(None if balance_tol is None else float(balance_tol)),
        balance_history=history, balance_final=history[-1],
        interval_balance_history=interval_history, interval_balance_final=interval_history[-1],
        field_history=field_history, field_final=field_final,
        quadrature_history=quadrature_history, quadrature_final=last_accepted_state["quadrature"],
        species_face_current=species_face_current,
        species_face_stderr=species_face_stderr,
        species_face_replicates=species_face_replicates, current_scale_m2_s=current_scale,
        species_endpoint_stderr=species_endpoint_stderr,
        species_endpoint_replicates=species_endpoint_replicates,
        ion_current=ion_current, electron_current=electron_current,
        ion_current_stderr=ion_stderr, electron_current_stderr=electron_stderr,
        active_flux=active_flux, surface_discretization="boundary_nodal",
        electrostatic_state=("surface_charge_poisson" if charge_mode else "dirichlet_voltage"),
        update_residual="mean_log_current_ratio_when_direction_resolved",
        trust_merit=("mean_rms_and_max_log_current_ratio" if trust_merit == "pareto"
                     else f"mean_{trust_merit}_log_current_ratio"),
        trust_merit_mode=trust_merit, rejected_steps=rejected_steps,
        trust_best_rms=float(last_accepted_state["trust_best_rms"]),
        trust_best_max=float(last_accepted_state["trust_best_max"]),
        beta_final=beta_current, trial_merit_history=np.asarray(trial_merit_history),
        trial_max_history=np.asarray(trial_max_history),
        accepted_beta_history=np.asarray(accepted_beta_history),
        accepted_gain_history=np.asarray(accepted_gain_history),
        gain_decay=float(gain_decay), gain_offset=float(gain_offset),
        nonlinear_update=nonlinear_update, anderson_depth=int(anderson_depth),
        accepted_iterations_total=int(initial_accepted_iterations) + len(history),
        restart_accepted_iterations=max(
            int(initial_accepted_iterations) + len(history) - 1, 0),
        restart_beta=float(last_accepted_state["beta_current"]),
        anderson_x_history=last_accepted_state["anderson_x"].copy(),
        anderson_residual_history=last_accepted_state["anderson_residual"].copy(),
        method_hint={
            name: value.copy() for name, value in last_accepted_state["method_hint"].items()},
        adaptive_levels={
            name: value.copy()
            for name, value in last_accepted_state["adaptive_levels"].items()},
        forward_adaptive_levels={
            name: value.copy()
            for name, value in last_accepted_state["forward_adaptive_levels"].items()})
