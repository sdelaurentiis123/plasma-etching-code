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


def solve_boundary_state_charging_nodal(
        solid, conductor_ids, boundary_state, *, ion_species=None, electron_species=None,
        initial_surface_voltage=None, initial_boundary_nodal_voltage=None,
        n_iter=40, beta=0.5, response_energy_eV=4.0,
        dVmax=8.0, balance_tol=1e-3, min_iter=2, field_sweeps=500,
        field_tolerance=1e-9, boundary_proposals=None, n_face_position=8,
        adaptive_quadrature=None, active_flux=1e-4, current_confidence_sigma=2.0,
        trust_region=True, trust_growth_tolerance=0.02, minimum_beta=1e-4):
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
    conductor_voltage = np.zeros(int(conductor_ids.max()) + 1)
    for component in components:
        values = (np.asarray([boundary_voltage[node] for node, owner in node_component.items()
                             if owner == component]) if initial_boundary_nodal_voltage is not None
                  else initial_cell[conductor_ids == component])
        conductor_voltage[component] = float(values.mean()) if values.size else 0.0

    def impose_conductors():
        for node, component in node_component.items():
            boundary_voltage[node] = conductor_voltage[component]
        boundary_voltage[:, 0] = 0.0

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

    hybrid_hint = {}
    history = []; interval_history = []; field_history = []; quadrature_history = []
    species_face_current = {}; species_face_stderr = {}; species_face_replicates = {}
    species_endpoint_stderr = {}; species_endpoint_replicates = {}
    beta_current = float(beta); pending_step = None; rejected_steps = 0
    last_accepted_state = None
    trial_merit_history = []; accepted_beta_history = []
    for iteration in range(int(n_iter)):
        impose_conductors()
        potential, field_diag = solve_nodal_laplace(
            solid, boundary_nodal_voltage=boundary_voltage,
            sweeps=field_sweeps, omega=1.7, tolerance=field_tolerance)
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
                options.pop("warm_start_backoff", None)
                options.setdefault("n_face_position", n_face_position)
                if bidirectional:
                    hybrid = bidirectional_boundary_state_cell_flux(
                        boundary_state, species.name, potential, solid, cells, normals_array,
                        proposal_species=proposals.get(species.name),
                        adjoint_options=options, forward_options=forward_options,
                        element_absolute_tolerance=options.get("element_absolute_tolerance", 1e-3),
                        element_relative_tolerance=options.get("element_relative_tolerance", 0.05),
                        method_hint=hybrid_hint.get(species.name), switch_factor=switch_factor,
                        consistency_sigma=consistency_sigma)
                    hybrid_hint[species.name] = hybrid["method"].copy()
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
                            cells=cells, normals=normals)
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
                            potential=potential, cells=cells, normals=normals)
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
        activity = (ion_current + electron_current) >= active_flux
        ilo = np.maximum(ion_current - current_confidence_sigma * ion_stderr, 0.0)
        ihi = ion_current + current_confidence_sigma * ion_stderr
        elo = np.maximum(electron_current - current_confidence_sigma * electron_stderr, 0.0)
        ehi = electron_current + current_confidence_sigma * electron_stderr
        certified = np.zeros(dof_count)
        ion_high = ilo > ehi; electron_high = elo > ihi
        certified[ion_high] = np.log(ilo[ion_high] / np.maximum(ehi[ion_high], 1e-300))
        certified[electron_high] = np.log(
            np.maximum(ihi[electron_high], 1e-300) / elo[electron_high])
        raw_active = np.abs(raw[activity]); cert_active = np.abs(certified[activity])
        balance = dict(
            log_ratio=raw.copy(), active=activity.copy(),
            max_abs_log_ratio=float(raw_active.max()) if raw_active.size else 0.0,
            rms_log_ratio=float(np.sqrt(np.mean(raw_active ** 2))) if raw_active.size else 0.0)
        interval_balance = dict(
            log_ratio=certified.copy(), active=activity.copy(),
            max_abs_log_ratio=float(cert_active.max()) if cert_active.size else 0.0,
            rms_log_ratio=float(np.sqrt(np.mean(cert_active ** 2))) if cert_active.size else 0.0)
        merit = interval_balance["rms_log_ratio"]
        trial_merit_history.append(merit)
        if (trust_region and pending_step is not None
                and merit > pending_step["merit"] * (1.0 + trust_growth_tolerance)):
            boundary_voltage[:] = pending_step["boundary_voltage"]
            conductor_voltage[:] = pending_step["conductor_voltage"]
            hybrid_hint = {
                name: value.copy() for name, value in pending_step["hybrid_hint"].items()}
            beta_current *= 0.5; rejected_steps += 1; pending_step = None
            if beta_current < minimum_beta:
                raise RuntimeError(
                    f"nodal charging trust region collapsed below minimum_beta={minimum_beta:g}")
            continue
        if trust_region and pending_step is not None and merit < 0.8 * pending_step["merit"]:
            beta_current = min(float(beta), 1.2 * beta_current)
        pending_step = None
        history.append(balance); interval_history.append(interval_balance)
        accepted_beta_history.append(beta_current)
        last_accepted_state = dict(
            boundary_voltage=boundary_voltage.copy(),
            conductor_voltage=conductor_voltage.copy(), potential=potential.copy(),
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
            quadrature=species_quadrature)
        if (balance_tol is not None and len(history) >= int(min_iter)
                and interval_history[-1]["max_abs_log_ratio"] <= balance_tol):
            break
        if trust_region:
            pending_step = dict(
                boundary_voltage=boundary_voltage.copy(),
                conductor_voltage=conductor_voltage.copy(), merit=merit,
                hybrid_hint={name: value.copy() for name, value in hybrid_hint.items()})
        # Iterate on the physical mean-current root. Confidence intervals certify that root; their
        # nearest edge is not a residual and would make the fixed point depend on sampling tolerance.
        update_residual = np.where(activity, raw, 0.0)
        step = np.clip(beta_current * response_energy_eV * update_residual, -dVmax, dVmax)
        for index, node in enumerate(dielectric_nodes):
            if activity[index]:
                boundary_voltage[node] += step[index]
        for component in components:
            index = component_index[component]
            if activity[index]:
                conductor_voltage[component] += step[index]

    # A final fixed-point proposal has not had its currents evaluated. Return the last accepted state
    # so voltage, potential, residual, and checkpoint restarts all describe one physical iterate.
    if last_accepted_state is None:
        raise RuntimeError("nodal charging solve produced no accepted state")
    boundary_voltage[:] = last_accepted_state["boundary_voltage"]
    conductor_voltage[:] = last_accepted_state["conductor_voltage"]
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
    potential, field_final = solve_nodal_laplace(
        solid, boundary_nodal_voltage=boundary_voltage,
        sweeps=field_sweeps, omega=1.7, tolerance=field_tolerance)
    return dict(
        solid=solid.copy(), conductor_ids=conductor_ids.copy(),
        surface_voltage=surface_readout(), boundary_nodal_voltage=boundary_voltage.copy(),
        potential=potential, cells=np.asarray(cells, dtype=int), normals=normals_array,
        face_components=face_components, dielectric_nodes=np.asarray(dielectric_nodes, dtype=int),
        conductor_voltage=conductor_voltage, iterations=len(history),
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
        update_residual="mean_log_current_ratio", rejected_steps=rejected_steps,
        beta_final=beta_current, trial_merit_history=np.asarray(trial_merit_history),
        accepted_beta_history=np.asarray(accepted_beta_history))
