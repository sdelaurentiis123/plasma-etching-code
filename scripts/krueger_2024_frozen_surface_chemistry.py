#!/usr/bin/env python3
"""Frozen-geometry finite-chemistry diagnostic for the Krueger base checkpoint.

The production feature step couples transport, stateful surface chemistry, and
level-set motion.  A zero-duration feature audit cannot expose SiO2 removal:
the reduced mechanism defines oxide velocity from *integrated* formula-unit
removal and therefore returns zero at exactly dt=0.  This diagnostic evaluates
transport once, freezes those face fluxes, and advances only the exact routed
surface mechanism over a bounded horizon.  It never receives a level set in
its chemistry kernel, so advection, redistancing, topology cleanup, and surface
state remapping are impossible by construction.

This is a diagnostic screen, not an endpoint prediction or validation claim.
Every usable horizon must pass chemistry-step refinement, material-ledger,
gross-displacement, and neutral-reaction-probability drift gates.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from krueger_2024_endpoint_operator_audit import _evaluate  # noqa: E402
from krueger_2024_frozen_checkpoint_2x2 import (  # noqa: E402
    EvaluationDeadlineExceeded,
    _hard_deadline,
    _hash_manifest,
    _inputs_unchanged,
    _jsonable,
    _load_source,
    _operator_config,
    _sha256,
    _snapshot_inputs,
    _write_json_atomic,
)
from krueger_2024_trench_pilot import _maximum_ledger_residual  # noqa: E402
from petch.amorphous_carbon_mask import (  # noqa: E402
    build_krueger_2024_material_router_3d,
)
from petch.feature_step_3d import _select_surface_fluxes  # noqa: E402
from petch.surface_kinetics import (  # noqa: E402
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)


SCHEMA = "petch.krueger-2024.frozen-surface-chemistry.v1"
PARAMETERS = {
    "r17": {
        "effective_mask_crosslinked_growth_fraction": 0.8934059741411972,
        "oxide_etch_yield_scale": 0.5667632723491973,
    },
    "r19": {
        "effective_mask_crosslinked_growth_fraction": 0.9004722559883319,
        "oxide_etch_yield_scale": 0.5586489665864749,
    },
}
TRANSPORT_OPERATOR = {
    "fidelity": "diagnostic_q3_screen",
    "authority": False,
    "boundary_case": "base",
    "boundary_mode": "legacy_compressed_tensor",
    "ion_energy_bin_eV": 500.0,
    "ion_angle_bin_deg": 0.5,
    "ballistic_transport": "face_gather",
    "ballistic_face_quadrature_points": 3,
    "n_position": 16,
    "neutral_radiosity_rays_per_face": 8,
    "transport_device": "cpu",
    "duration_s": 0.0,
}
GATES = {
    "maximum_refinement_relative_error": 0.01,
    "maximum_gross_displacement_dx": 0.05,
    "maximum_neutral_reaction_probability_absolute_drift": 0.01,
    "maximum_flux_weighted_absorption_drift": 0.01,
    "maximum_radiosity_relative_balance_error": 5.0e-12,
    # Coverage = 1-exp(-N/Nmono).  Below one part per million, the film changes
    # neutral access by less than the diagnostic's one-percent flux gate.
    "effective_film_coverage_threshold": 1.0e-6,
}
DIAGNOSTIC_SOURCE_PATHS = (
    "scripts/krueger_2024_frozen_surface_chemistry.py",
    "scripts/krueger_2024_frozen_checkpoint_2x2.py",
    "scripts/krueger_2024_endpoint_operator_audit.py",
    "scripts/krueger_2024_trench_pilot.py",
)
RUNTIME_SOURCE_PATHS = tuple(sorted(
    str(path.relative_to(ROOT))
    for path in (ROOT / "src" / "petch").glob("*.py")
    if not path.name.startswith("._")
))
SOURCE_PATHS = DIAGNOSTIC_SOURCE_PATHS + RUNTIME_SOURCE_PATHS
BASE_INPUT_PATHS = (
    "data/experimental/krueger_2024/base_case_boundary_fluxes.csv",
    "data/experimental/krueger_2024/digitized_figure4_iead.csv",
    "data/experimental/krueger_2024/digitized_figure4_iead_metadata.json",
)


def _array_hash(digest, label, value):
    array = np.ascontiguousarray(np.asarray(value))
    digest.update(str(label).encode("utf-8"))
    digest.update(b"\0")
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(array.shape).encode("ascii"))
    digest.update(b"\0")
    digest.update(array.tobytes())
    digest.update(b"\n")


def surface_flux_sha256(fluxes):
    """Hash the exact immutable neutral and energetic face measures."""
    if not isinstance(fluxes, SurfaceFluxes):
        raise TypeError("fluxes must be SurfaceFluxes")
    digest = sha256()
    for name, values in sorted(fluxes.neutral_flux_m2_s.items()):
        _array_hash(digest, f"neutral.{name}", values)
    for index, population in enumerate(fluxes.energetic_fluxes):
        prefix = f"energetic.{index}.{population.name}"
        digest.update(type(population).__name__.encode("ascii") + b"\n")
        if isinstance(population, FaceResolvedEnergeticFlux):
            digest.update(f"{prefix}.face_count={population.face_count}\n".encode("ascii"))
            for field in (
                    "event_face", "event_flux_m2_s", "event_energy_eV",
                    "event_cosine_incidence", "event_position",
                    "event_incident_direction"):
                value = getattr(population, field)
                if value is not None:
                    _array_hash(digest, f"{prefix}.{field}", value)
        elif isinstance(population, EnergeticFlux):
            for field in ("flux_m2_s", "energy_eV", "cosine_incidence", "weight"):
                _array_hash(digest, f"{prefix}.{field}", getattr(population, field))
        else:  # pragma: no cover - SurfaceFluxes validates this.
            raise TypeError(type(population).__name__)
    return digest.hexdigest()


def _inventory_add(target, source):
    for name, values in source.items():
        values = np.asarray(values, dtype=float)
        if name not in target:
            target[name] = np.zeros(values.shape, dtype=float)
        target[name] += values


def _integrated_inventory(inventory, area_m2):
    return {
        name: float(np.sum(np.asarray(values, dtype=float) * area_m2))
        for name, values in sorted(inventory.items())
    }


def _reaction_probability_drift(
        mechanism, initial_state, final_state, active_material, fixed_fluxes,
        active_area_m2):
    initial = mechanism.neutral_reaction_probability_by_material(
        initial_state, active_material)
    final = mechanism.neutral_reaction_probability_by_material(
        final_state, active_material)
    species = sorted(set(initial) | set(final))
    by_species = {}
    maximum_absolute = 0.0
    maximum_weighted = 0.0
    for name in species:
        before = np.asarray(initial.get(name, 0.0), dtype=float)
        after = np.asarray(final.get(name, 0.0), dtype=float)
        before, after = np.broadcast_arrays(before, after)
        absolute = np.abs(after - before)
        local_maximum = float(np.max(absolute)) if absolute.size else 0.0
        flux = np.asarray(fixed_fluxes.neutral_flux_m2_s.get(name, 0.0), dtype=float)
        flux = np.broadcast_to(flux, active_material.shape)
        incident = float(np.sum(active_area_m2 * flux))
        weighted = (
            float(np.sum(active_area_m2 * flux * absolute)) / incident
            if incident > 0.0 else 0.0
        )
        by_species[name] = {
            "maximum_absolute_probability_drift": local_maximum,
            "flux_weighted_absorption_drift": weighted,
            "incident_rate_s-1": incident,
        }
        maximum_absolute = max(maximum_absolute, local_maximum)
        maximum_weighted = max(maximum_weighted, weighted)
    return {
        "maximum_absolute_probability_drift": maximum_absolute,
        "maximum_flux_weighted_absorption_drift": maximum_weighted,
        "by_species": by_species,
    }


def advance_frozen_surface_chemistry_3d(
        state, fixed_surface_fluxes, active_face_area_m2, active_face_material_id,
        mechanism, *, horizon_s, substeps, dx_m):
    """Advance exact routed chemistry under fixed flux without accepting geometry.

    The input state and flux objects are immutable contracts.  A new state and
    accumulated integrated exchange are returned; no mesh, level set, or remap
    entry point is present.
    """
    horizon_s = float(horizon_s)
    substeps = int(substeps)
    dx_m = float(dx_m)
    material = np.asarray(active_face_material_id, dtype=int)
    area = np.asarray(active_face_area_m2, dtype=float)
    if (not np.isfinite(horizon_s) or horizon_s <= 0.0
            or substeps <= 0 or not np.isfinite(dx_m) or dx_m <= 0.0
            or area.shape != material.shape or np.any(area <= 0.0)):
        raise ValueError("invalid frozen-chemistry horizon, mesh measure, or substeps")
    state_before = {
        name: np.asarray(values).copy() for name, values in state.fields.items()
    }
    flux_hash_before = surface_flux_sha256(fixed_surface_fluxes)
    initial_state = state
    current = state
    removed = {}; outgoing = {}; unresolved = {}; deposited = {}
    gross_displacement = np.zeros(material.shape)
    net_displacement = np.zeros(material.shape)
    maximum_ledger_residual = 0.0
    all_valid = True
    validity_reasons = []
    step_s = horizon_s / substeps
    for _ in range(substeps):
        result = mechanism.advance_by_material(
            current, fixed_surface_fluxes, step_s, material)
        exchange = result.material_exchange
        _inventory_add(removed, exchange.removed_units_m2)
        _inventory_add(outgoing, exchange.outgoing_units_m2)
        _inventory_add(unresolved, exchange.unresolved_units_m2)
        _inventory_add(deposited, exchange.deposited_units_m2)
        maximum_ledger_residual = max(
            maximum_ledger_residual, float(_maximum_ledger_residual(exchange)))
        recession = np.asarray(result.etch_velocity_m_s, dtype=float)
        growth = np.asarray(result.normal_growth_velocity_m_s, dtype=float)
        gross_displacement += (recession + growth) * step_s
        net_displacement += (recession - growth) * step_s
        all_valid &= bool(result.validity.within_declared_scope)
        validity_reasons.extend(result.validity.reasons)
        current = result.state
    if any(
            not np.array_equal(values, state.fields[name])
            for name, values in state_before.items()):
        raise RuntimeError("frozen chemistry mutated its input surface state")
    if surface_flux_sha256(fixed_surface_fluxes) != flux_hash_before:
        raise RuntimeError("frozen chemistry mutated its fixed transport flux")

    cumulative_residual = 0.0
    for name, values in removed.items():
        residual = (
            values - np.asarray(outgoing.get(name, 0.0))
            - np.asarray(unresolved.get(name, 0.0))
        )
        cumulative_residual = max(
            cumulative_residual,
            float(np.max(np.abs(residual))) if residual.size else 0.0,
        )
    probability = _reaction_probability_drift(
        mechanism, initial_state, current, material, fixed_surface_fluxes, area)
    oxide = material == 1
    polymer_name = "m1__polymer_units_m2"
    polymer_initial = np.asarray(initial_state.fields[polymer_name], dtype=float)
    polymer_final = np.asarray(current.fields[polymer_name], dtype=float)
    monolayer = float(
        mechanism.mechanisms[1].parameters.polymer_monolayer_density_m2)
    coverage_initial = 1.0 - np.exp(-polymer_initial / monolayer)
    coverage_final = 1.0 - np.exp(-polymer_final / monolayer)
    exact_depleted = oxide & (polymer_initial > 0.0) & (polymer_final == 0.0)
    effective_depleted = (
        oxide
        & (coverage_initial > GATES["effective_film_coverage_threshold"])
        & (coverage_final <= GATES["effective_film_coverage_threshold"])
    )
    density = float(mechanism.mechanisms[1].parameters.bulk_formula_density_m3)
    integrated_removed = _integrated_inventory(removed, area)
    oxide_units = float(integrated_removed.get("SiO2_formula_unit", 0.0))
    oxide_volume = oxide_units / density
    return {
        "horizon_s": horizon_s,
        "substeps": substeps,
        "chemistry_step_s": step_s,
        "state": current,
        "fixed_surface_flux_sha256": flux_hash_before,
        "integrated_exchange": {
            "removed_units": integrated_removed,
            "outgoing_units": _integrated_inventory(outgoing, area),
            "unresolved_units": _integrated_inventory(unresolved, area),
            "deposited_units": _integrated_inventory(deposited, area),
            "maximum_step_ledger_residual_units_m2": maximum_ledger_residual,
            "maximum_cumulative_ledger_residual_units_m2": cumulative_residual,
        },
        "oxide_removal": {
            "integrated_formula_units": oxide_units,
            "integrated_volume_m3": oxide_volume,
            "mean_normal_thickness_over_exposed_oxide_m": (
                oxide_volume / float(np.sum(area[oxide])) if np.any(oxide) else 0.0
            ),
        },
        "displacement": {
            "maximum_gross_displacement_m": float(np.max(gross_displacement)),
            "maximum_gross_displacement_dx": float(
                np.max(gross_displacement) / dx_m),
            "maximum_absolute_net_displacement_m": float(
                np.max(np.abs(net_displacement))),
        },
        "neutral_reaction_probability_drift": probability,
        "film_depletion": {
            "oxide_face_count": int(np.sum(oxide)),
            "initially_film_covered_oxide_face_count": int(
                np.sum(oxide & (polymer_initial > 0.0))),
            "exactly_depleted_oxide_face_count": int(np.sum(exact_depleted)),
            "effectively_depleted_oxide_face_count": int(np.sum(effective_depleted)),
            "effectively_depleted_oxide_face_fraction": float(
                np.sum(effective_depleted) / max(np.sum(oxide), 1)),
            "effective_film_coverage_threshold": GATES[
                "effective_film_coverage_threshold"],
        },
        "validity": {
            "within_declared_scope": bool(all_valid),
            "reasons": sorted(set(str(item) for item in validity_reasons)),
        },
    }


def _refinement_relative_error(coarse, fine, initial_state):
    by_exchange = {}
    maximum = 0.0
    for exchange_kind in ("removed_units", "deposited_units"):
        fields = set(coarse["integrated_exchange"][exchange_kind])
        fields |= set(fine["integrated_exchange"][exchange_kind])
        local = {}
        for name in sorted(fields):
            left = float(coarse["integrated_exchange"][exchange_kind].get(name, 0.0))
            right = float(fine["integrated_exchange"][exchange_kind].get(name, 0.0))
            relative = abs(left - right) / max(abs(left), abs(right), 1.0)
            local[name] = relative
            maximum = max(maximum, relative)
        by_exchange[exchange_kind] = local

    by_state_increment = {}
    for name in sorted(initial_state.fields):
        initial = np.asarray(initial_state.fields[name], dtype=float)
        coarse_increment = np.asarray(coarse["state"].fields[name]) - initial
        fine_increment = np.asarray(fine["state"].fields[name]) - initial
        difference = float(np.max(np.abs(coarse_increment - fine_increment)))
        scale = max(
            float(np.max(np.abs(coarse_increment))),
            float(np.max(np.abs(fine_increment))),
            64.0 * np.finfo(float).eps
            * max(float(np.max(np.abs(initial))), 1.0),
        )
        relative = difference / scale
        by_state_increment[name] = relative
        maximum = max(maximum, relative)
    return {
        "maximum_relative_error": maximum,
        "by_exchange_inventory": by_exchange,
        "by_final_state_increment": by_state_increment,
    }


def evaluate_frozen_horizon(
        state, fixed_surface_fluxes, active_face_area_m2, active_face_material_id,
        mechanism, *, horizon_s, coarse_substeps, dx_m):
    """Evaluate one horizon with N/2N chemistry refinement and all refusal gates."""
    coarse = advance_frozen_surface_chemistry_3d(
        state, fixed_surface_fluxes, active_face_area_m2, active_face_material_id,
        mechanism, horizon_s=horizon_s, substeps=coarse_substeps, dx_m=dx_m)
    fine = advance_frozen_surface_chemistry_3d(
        state, fixed_surface_fluxes, active_face_area_m2, active_face_material_id,
        mechanism, horizon_s=horizon_s, substeps=2 * coarse_substeps, dx_m=dx_m)
    refinement = _refinement_relative_error(coarse, fine, state)
    probability = fine["neutral_reaction_probability_drift"]
    exchange = fine["integrated_exchange"]
    gates = {
        "chemistry_step_refinement": bool(
            refinement["maximum_relative_error"]
            <= GATES["maximum_refinement_relative_error"]),
        "exact_material_ledger": bool(
            exchange["maximum_step_ledger_residual_units_m2"] == 0.0
            and exchange["maximum_cumulative_ledger_residual_units_m2"] == 0.0),
        "gross_displacement_within_frozen_geometry_limit": bool(
            fine["displacement"]["maximum_gross_displacement_dx"]
            <= GATES["maximum_gross_displacement_dx"]),
        "neutral_reaction_probability_stable": bool(
            probability["maximum_absolute_probability_drift"]
            <= GATES["maximum_neutral_reaction_probability_absolute_drift"]),
        "flux_weighted_absorption_stable": bool(
            probability["maximum_flux_weighted_absorption_drift"]
            <= GATES["maximum_flux_weighted_absorption_drift"]),
        "within_declared_scope": bool(
            coarse["validity"]["within_declared_scope"]
            and fine["validity"]["within_declared_scope"]),
    }
    return {
        "horizon_s": float(horizon_s),
        "coarse_substeps": int(coarse_substeps),
        "fine_substeps": int(2 * coarse_substeps),
        "refinement": refinement,
        "fine": {key: value for key, value in fine.items() if key != "state"},
        "gates": gates,
        "all_gates_pass": bool(all(gates.values())),
    }


def _maximum_radiosity_error(result):
    return max(
        (float(item["relative_balance_error"])
         for item in result.diagnostics.get("neutral_radiosity", {}).values()),
        default=0.0,
    )


def _parameter_transport(source, parameter_label, *, seed, deadline_s):
    parameters = PARAMETERS[parameter_label]
    config = _operator_config(source["config"], parameters)
    snapshot = _snapshot_inputs(source["geometry"], source["state"])
    started = perf_counter()
    with _hard_deadline(deadline_s):
        result, boundary, reported_wall = _evaluate(
            source["geometry"], source["state"], source["fingerprint"],
            boundary_mode=TRANSPORT_OPERATOR["boundary_mode"],
            ion_bins=(TRANSPORT_OPERATOR["ion_energy_bin_eV"],
                      TRANSPORT_OPERATOR["ion_angle_bin_deg"]),
            face_points=TRANSPORT_OPERATOR["ballistic_face_quadrature_points"],
            pilot_config=config,
            radiosity_rays=TRANSPORT_OPERATOR["neutral_radiosity_rays_per_face"],
            seed=int(seed), ballistic_transport="face_gather",
            n_position=TRANSPORT_OPERATOR["n_position"], transport_device="cpu")
    wall = perf_counter() - started
    if not _inputs_unchanged(snapshot, source["geometry"], source["state"]):
        raise RuntimeError("zero-duration transport mutated checkpoint inputs")
    role = {
        species.name: (
            "energetic_bombardment" if species.charge_number != 0
            else "neutral_reactant")
        for species in boundary.species
    }
    active_flux = _select_surface_fluxes(
        result.transport.surface_fluxes, result.active_face_index,
        len(result.face_material_id), role)
    active_material = np.asarray(result.face_material_id, dtype=int)[
        np.asarray(result.active_face_index, dtype=int)]
    active_area_m2 = (
        np.asarray(result.active_face_area, dtype=float)
        * source["geometry"].mesh_length_unit_m ** 2)
    mechanism = build_krueger_2024_material_router_3d(**parameters)
    radiosity_error = _maximum_radiosity_error(result)
    transport_gates = {
        "input_checkpoint_unchanged": True,
        "material_ledger_exact": bool(
            _maximum_ledger_residual(result.surface.material_exchange) == 0.0),
        "radiosity_balance_within_tolerance": bool(
            radiosity_error <= GATES["maximum_radiosity_relative_balance_error"]),
        "within_declared_scope": bool(result.validity.within_declared_scope),
    }
    return {
        "parameter_label": parameter_label,
        "parameters": parameters,
        "wall_time_s": wall,
        "reported_evaluator_wall_time_s": float(reported_wall),
        "fixed_surface_flux_sha256": surface_flux_sha256(active_flux),
        "maximum_radiosity_relative_balance_error": radiosity_error,
        "gates": transport_gates,
        "all_gates_pass": bool(all(transport_gates.values())),
        "boundary_provenance": _jsonable(boundary.provenance),
        "active_flux": active_flux,
        "active_material": active_material,
        "active_area_m2": active_area_m2,
        "mechanism": mechanism,
    }


def _public_transport(item):
    return {
        key: value for key, value in item.items()
        if key not in {"active_flux", "active_material", "active_area_m2", "mechanism"}
    }


def run(args):
    started = perf_counter()
    source = _load_source("r19", args.r19_source)
    r17_audit_path = Path(args.r17_source) / "audit.json"
    r17_audit = json.loads(r17_audit_path.read_text(encoding="utf-8"))
    if r17_audit.get("status") != "complete":
        raise ValueError("R17 parameter provenance requires a completed base audit")
    if r17_audit["configuration"].get("boundary_case") != "base":
        raise ValueError("R17 parameter provenance is not a base-boundary run")
    for label, config in (("r17", r17_audit["configuration"]),
                          ("r19", source["config"])):
        observed = {name: float(config[name]) for name in PARAMETERS[label]}
        if observed != PARAMETERS[label]:
            raise ValueError(f"{label.upper()} mechanism parameters disagree with provenance")
    if source["checkpoint_metadata"].get("physical_time_s") != 60.0:
        raise ValueError("frozen chemistry requires the completed 60 s R19 checkpoint")
    dt_next = float(source["checkpoint_metadata"]["next_step_duration_s"])
    horizon_fractions = (1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0)
    destination = Path(args.output)
    payload = {
        "schema": SCHEMA,
        "status": "running",
        "scientific_scope": (
            "diagnostic-only finite surface chemistry under frozen q3 transport and frozen "
            "geometry; no profile evolution, remap, calibration, or validation claim"),
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
        },
        "transport_operator": TRANSPORT_OPERATOR,
        "gates": GATES,
        "checkpoint": {
            "audit_sha256": _sha256(source["audit_path"]),
            "checkpoint_sha256": _sha256(source["checkpoint_path"]),
            "metadata": _jsonable(source["checkpoint_metadata"]),
            "metrics": _jsonable(source["metrics"]),
        },
        "parameter_provenance": {
            "r17_audit_sha256": _sha256(r17_audit_path),
            "r19_audit_sha256": _sha256(source["audit_path"]),
            "parameter_pairs": PARAMETERS,
        },
        "execution_budget": {
            "maximum_transport_cell_wall_s": float(args.maximum_transport_cell_wall_s),
            "maximum_total_wall_s": float(args.maximum_total_wall_s),
            "next_profile_step_s": dt_next,
            "horizon_fractions": list(horizon_fractions),
        },
        "transport": {},
        "horizons": [],
        "largest_common_passing_horizon": None,
        "provenance": {
            "source": _hash_manifest(SOURCE_PATHS),
            "base_inputs": _hash_manifest(BASE_INPUT_PATHS),
        },
    }
    runtime = {}
    for label in ("r17", "r19"):
        remaining = float(args.maximum_total_wall_s) - (perf_counter() - started)
        deadline = min(float(args.maximum_transport_cell_wall_s), remaining)
        try:
            item = _parameter_transport(
                source, label, seed=int(args.seed), deadline_s=deadline)
        except EvaluationDeadlineExceeded as error:
            payload["status"] = "bounded_transport_timeout"
            payload["execution_budget"]["timeout"] = {
                "parameter_label": label, "deadline_s": deadline, "reason": str(error)}
            payload["total_wall_time_s"] = float(perf_counter() - started)
            _write_json_atomic(destination, payload)
            return payload
        runtime[label] = item
        payload["transport"][label] = _public_transport(item)
        _write_json_atomic(destination, payload)
        if not item["all_gates_pass"]:
            payload["status"] = "transport_gate_failure"
            payload["total_wall_time_s"] = float(perf_counter() - started)
            _write_json_atomic(destination, payload)
            return payload

    dx_m = float(source["geometry"].dx * source["geometry"].mesh_length_unit_m)
    largest = None
    for index, fraction in enumerate(horizon_fractions):
        if perf_counter() - started >= float(args.maximum_total_wall_s):
            payload["status"] = "bounded_total_timeout"
            break
        horizon = dt_next * fraction
        entry = {
            "fraction_of_next_profile_step": fraction,
            "horizon_s": horizon,
            "parameter_results": {},
            "common_pass": False,
        }
        for label in ("r17", "r19"):
            item = runtime[label]
            result = evaluate_frozen_horizon(
                source["state"], item["active_flux"], item["active_area_m2"],
                item["active_material"], item["mechanism"],
                horizon_s=horizon, coarse_substeps=2 ** index, dx_m=dx_m)
            entry["parameter_results"][label] = result
            if not result["all_gates_pass"]:
                entry["first_failure"] = {"parameter_label": label, "gates": result["gates"]}
                break
        if len(entry["parameter_results"]) == 2 and all(
                item["all_gates_pass"] for item in entry["parameter_results"].values()):
            entry["common_pass"] = True
            r17_removed = entry["parameter_results"]["r17"]["fine"][
                "oxide_removal"]["integrated_formula_units"]
            r19_removed = entry["parameter_results"]["r19"]["fine"][
                "oxide_removal"]["integrated_formula_units"]
            entry["paired_oxide_removal_direction"] = {
                "r19_minus_r17_integrated_formula_units": r19_removed - r17_removed,
                "r19_to_r17_ratio": (
                    r19_removed / r17_removed if r17_removed > 0.0 else None),
                "direction": (
                    "r19_lower" if r19_removed < r17_removed
                    else "r19_higher" if r19_removed > r17_removed else "equal"),
            }
            largest = entry
        payload["horizons"].append(entry)
        _write_json_atomic(destination, payload)
        if not entry["common_pass"]:
            payload["status"] = "stopped_at_first_failed_horizon"
            break
    else:
        payload["status"] = "pass"
    if largest is not None:
        payload["largest_common_passing_horizon"] = {
            "fraction_of_next_profile_step": largest["fraction_of_next_profile_step"],
            "horizon_s": largest["horizon_s"],
            "paired_oxide_removal_direction": largest[
                "paired_oxide_removal_direction"],
        }
    payload["total_wall_time_s"] = float(perf_counter() - started)
    _write_json_atomic(destination, payload)
    print(json.dumps({
        "status": payload["status"],
        "largest_common_passing_horizon": payload["largest_common_passing_horizon"],
        "output": str(destination),
        "total_wall_time_s": payload["total_wall_time_s"],
    }, indent=2, sort_keys=True))
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--r17-source",
        default=(ROOT / "results" / "krueger_2024_base_calibration_r17"
                 / "axisym_candidate"))
    parser.add_argument(
        "--r19-source",
        default=ROOT / "results" / "krueger_2024_r19_response_check" / "remote_artifacts")
    parser.add_argument(
        "--output",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "frozen_surface_chemistry" / "audit.json"))
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument("--maximum-transport-cell-wall-s", type=float, default=180.0)
    parser.add_argument("--maximum-total-wall-s", type=float, default=480.0)
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
