#!/usr/bin/env python3
"""Frozen-geometry co-integration of Krueger surface state and neutral radiosity.

The diagnostic evaluates direct q3 transport once, estimates one periodic
diffuse form-factor operator on the identical checkpoint mesh, and then
subcycles the exact material-router chemistry while resolving neutral
radiosity whenever the accepted surface state changes.  Geometry is never an
input to the chemistry integrator: no advection, redistance, topology cleanup,
or remap can occur.

This remains a non-authority diagnostic.  Its cached form-factor implementation
must first reproduce the production one-shot radiosity facewise and the
archived q3 integrated audit.  Each adaptive candidate is evaluated as one
full step and as two half steps with a midpoint radiosity re-solve.  Only the
fine path can be accepted, and only when every state increment, exchange
inventory, and per-face recession/growth measure agrees within tolerance.
Local reaction-probability change remains a separate hard safety cap.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
from unittest.mock import patch

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
from krueger_2024_frozen_surface_chemistry import (  # noqa: E402
    BASE_INPUT_PATHS,
    PARAMETERS,
    _integrated_inventory,
    _inventory_add,
    surface_flux_sha256,
)
from krueger_2024_trench_pilot import _maximum_ledger_residual  # noqa: E402
from petch.amorphous_carbon_mask import (  # noqa: E402
    build_krueger_2024_material_router_3d,
)
import petch.feature_step_3d as feature_step_module  # noqa: E402
from petch.feature_step_3d import (  # noqa: E402
    _select_surface_fluxes,
)
from petch.neutral_radiosity_3d import (  # noqa: E402
    DiffuseFormFactors3D,
    DiffuseNeutralNoSinkError,
    solve_diffuse_neutral_radiosity_3d,
)
from petch.surface_kinetics import SurfaceFluxes  # noqa: E402
from petch.threed import extract_mesh_3d  # noqa: E402


SCHEMA = "petch.krueger-2024.frozen-radiosity-chemistry.v1"
TRANSPORT_OPERATOR = {
    "fidelity": "diagnostic_q3_direct_then_cached_radiosity",
    "authority": False,
    "boundary_case": "base",
    "boundary_mode": "legacy_compressed_tensor",
    "ion_energy_bin_eV": 500.0,
    "ion_angle_bin_deg": 0.5,
    "ballistic_transport": "face_gather",
    "ballistic_face_quadrature_points": 3,
    "n_position": 16,
    "transport_device": "cpu",
    "duration_s": 0.0,
}
RADIOSITY = {
    "rays_per_face": 8,
    "maximum_rays_per_face": 64,
    "seed_offset": 10000,
    "periodic_lateral": True,
    "ray_offset_dx": 1.0e-3,
    "relative_tolerance": 1.0e-12,
    "maximum_iterations": 2000,
}
GATES = {
    "maximum_embedded_relative_error": 0.01,
    "maximum_local_reaction_probability_change": 0.01,
    "tight_maximum_local_reaction_probability_change": 0.005,
    "minimum_substep_s": 1.0e-4,
    "maximum_radiosity_relative_balance_error": 5.0e-12,
    "maximum_cached_reference_relative_error": 1.0e-12,
    "maximum_tolerance_halving_relative_oxide_error": 0.01,
    "maximum_gross_displacement_dx": 0.05,
}
HORIZON_FRACTIONS = (1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0)
MAXIMUM_CLI_RAYS_PER_FACE = 1024
DIAGNOSTIC_SOURCES = (
    "scripts/krueger_2024_frozen_radiosity_chemistry.py",
    "scripts/krueger_2024_frozen_surface_chemistry.py",
    "scripts/krueger_2024_frozen_checkpoint_2x2.py",
    "scripts/krueger_2024_endpoint_operator_audit.py",
    "scripts/krueger_2024_trench_pilot.py",
)
RUNTIME_SOURCES = tuple(sorted(
    str(path.relative_to(ROOT))
    for path in (ROOT / "src" / "petch").glob("*.py")
    if not path.name.startswith("._")
))
SOURCE_PATHS = DIAGNOSTIC_SOURCES + RUNTIME_SOURCES


class CoIntegrationRefusal(RuntimeError):
    """A chemistry/radiosity substep could not pass above the minimum dt."""


def _power_of_two_rays(value):
    try:
        rays = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("rays per face must be an integer") from error
    if (str(value).strip() != str(rays) or rays <= 0
            or rays & (rays - 1) or rays > MAXIMUM_CLI_RAYS_PER_FACE):
        raise argparse.ArgumentTypeError(
            "rays per face must be a positive power of two no greater than "
            f"{MAXIMUM_CLI_RAYS_PER_FACE}")
    return rays


def _unit_horizon_fraction(value):
    try:
        fraction = float(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("horizon fraction must be numeric") from error
    if not np.isfinite(fraction) or not 0.0 < fraction <= 1.0:
        raise argparse.ArgumentTypeError("horizon fraction must lie in (0, 1]")
    return fraction


def _select_horizon_fractions(requested, maximum):
    fractions = tuple(float(value) for value in requested)
    maximum = float(maximum)
    if (not fractions or any(not np.isfinite(value) or not 0.0 < value <= 1.0
                             for value in fractions)
            or not np.isfinite(maximum) or not 0.0 < maximum <= 1.0):
        raise ValueError("horizon fractions and maximum must lie in (0, 1]")
    if any(right <= left for left, right in zip(fractions, fractions[1:])):
        raise ValueError("horizon fractions must be unique and strictly increasing")
    selected = tuple(value for value in fractions if value <= maximum)
    if not selected:
        raise ValueError(
            "maximum horizon fraction excludes every requested horizon")
    return selected


def _certification_mode_for_rays(rays_per_face):
    return (
        "archived_q3_8ray_plus_production_facewise_cache_exact"
        if int(rays_per_face) == RADIOSITY["rays_per_face"]
        else "nondefault_rays_production_facewise_cache_exact")


def form_factor_sha256(factors):
    if not isinstance(factors, DiffuseFormFactors3D):
        raise TypeError("factors must be DiffuseFormFactors3D")
    digest = sha256()
    digest.update(f"face_count={factors.face_count}\n".encode("ascii"))
    digest.update(f"rays_per_face={factors.rays_per_face}\n".encode("ascii"))
    for name in ("source_face", "target_face", "transfer_fraction", "escape_fraction"):
        array = np.ascontiguousarray(np.asarray(getattr(factors, name)))
        digest.update(name.encode("ascii") + b"\0")
        digest.update(str(array.dtype).encode("ascii") + b"\0")
        digest.update(json.dumps(array.shape).encode("ascii") + b"\0")
        digest.update(array.tobytes())
    return digest.hexdigest()


def _active_probability(mechanism, state, active_material):
    return {
        name: np.asarray(values, dtype=float)
        for name, values in mechanism.neutral_reaction_probability_by_material(
            state, active_material).items()
    }


def _maximum_probability_change(before, after):
    maximum = 0.0
    by_species = {}
    for name in sorted(set(before) | set(after)):
        left, right = np.broadcast_arrays(
            np.asarray(before.get(name, 0.0), dtype=float),
            np.asarray(after.get(name, 0.0), dtype=float))
        value = float(np.max(np.abs(right - left))) if left.size else 0.0
        by_species[name] = value
        maximum = max(maximum, value)
    return maximum, by_species


def solve_cached_radiosity(
        direct_surface_fluxes, full_face_area_m2, factors, state, mechanism,
        active_face_index, active_material, *, relative_tolerance,
        maximum_iterations):
    """Solve the production radiosity equation on one cached geometric operator."""
    area = np.asarray(full_face_area_m2, dtype=float)
    active = np.asarray(active_face_index, dtype=int)
    if factors.face_count != area.size:
        raise ValueError("cached form factors do not match the direct transport mesh")
    if active.size != area.size or not np.array_equal(active, np.arange(area.size)):
        raise ValueError(
            "diagnostic requires every extracted interface face to use the routed mechanism")
    active_probability = _active_probability(mechanism, state, active_material)
    neutral = {}
    diagnostics = {}
    full_probability = {}
    for name, direct in direct_surface_fluxes.neutral_flux_m2_s.items():
        probability = np.asarray(active_probability.get(name, np.zeros(active.size)), dtype=float)
        full_probability[name] = probability
        direct = np.asarray(direct, dtype=float)
        if not np.any(probability > 0.0):
            source_rate = float(np.sum(area * direct))
            neutral[name] = direct
            diagnostics[name] = {
                "source_rate_s": source_rate,
                "reacted_rate_s": 0.0,
                "escaped_rate_s": source_rate,
                "relative_balance_error": 0.0,
                "relative_linear_residual": 0.0,
                "solver_method": "analytic_zero_reaction_elision",
            }
            continue
        solution = solve_diffuse_neutral_radiosity_3d(
            direct, area, factors.source_face, factors.target_face,
            factors.transfer_fraction, factors.escape_fraction, probability,
            relative_tolerance=float(relative_tolerance),
            maximum_iterations=int(maximum_iterations))
        neutral[name] = solution.incident_flux_m2_s
        diagnostics[name] = {
            "source_rate_s": solution.source_rate_s,
            "reacted_rate_s": solution.reacted_rate_s,
            "escaped_rate_s": solution.escaped_rate_s,
            "relative_balance_error": solution.relative_balance_error,
            "relative_linear_residual": solution.relative_linear_residual,
            "solver_method": solution.solver_method,
            "iteration_count": solution.iteration_count,
            "inactive_face_count": solution.inactive_face_count,
        }
    return (
        SurfaceFluxes(neutral, direct_surface_fluxes.energetic_fluxes),
        diagnostics,
        full_probability,
    )


def _maximum_balance(diagnostics):
    return max(
        (float(item["relative_balance_error"]) for item in diagnostics.values()),
        default=0.0,
    )


def _exchange_accumulator():
    return {"removed": {}, "outgoing": {}, "unresolved": {}, "deposited": {}}


def _accumulate_exchange(target, exchange):
    _inventory_add(target["removed"], exchange.removed_units_m2)
    _inventory_add(target["outgoing"], exchange.outgoing_units_m2)
    _inventory_add(target["unresolved"], exchange.unresolved_units_m2)
    _inventory_add(target["deposited"], exchange.deposited_units_m2)


def _cumulative_ledger_residual(exchange):
    maximum = 0.0
    for name, removed in exchange["removed"].items():
        residual = (
            removed - np.asarray(exchange["outgoing"].get(name, 0.0))
            - np.asarray(exchange["unresolved"].get(name, 0.0))
        )
        maximum = max(
            maximum,
            float(np.max(np.abs(residual))) if residual.size else 0.0)
    return maximum


def _exchange_from_steps(*steps):
    exchange = _exchange_accumulator()
    for step in steps:
        _accumulate_exchange(exchange, step.material_exchange)
    return exchange


def _relative_array_error(left, right, *, scale_floor):
    left, right = np.broadcast_arrays(
        np.asarray(left, dtype=float), np.asarray(right, dtype=float))
    difference = float(np.max(np.abs(left - right))) if left.size else 0.0
    scale = max(
        float(np.max(np.abs(left))) if left.size else 0.0,
        float(np.max(np.abs(right))) if right.size else 0.0,
        float(scale_floor),
    )
    return difference / scale


def _embedded_step_error(
        initial_state, coarse_step, fine_step, coarse_exchange, fine_exchange,
        coarse_recession_m, fine_recession_m, coarse_growth_m, fine_growth_m,
        *, dx_m):
    """Compare one full step to two half steps on every accepted quantity."""
    maximum = 0.0
    state_fields = {}
    for name in sorted(initial_state.fields):
        initial = np.asarray(initial_state.fields[name], dtype=float)
        coarse_increment = np.asarray(coarse_step.state.fields[name]) - initial
        fine_increment = np.asarray(fine_step.state.fields[name]) - initial
        floor = 64.0 * np.finfo(float).eps * max(
            float(np.max(np.abs(initial))) if initial.size else 0.0, 1.0)
        relative = _relative_array_error(
            coarse_increment, fine_increment, scale_floor=floor)
        state_fields[name] = relative
        maximum = max(maximum, relative)

    inventories = {}
    for kind in ("removed", "outgoing", "unresolved", "deposited"):
        local = {}
        names = sorted(set(coarse_exchange[kind]) | set(fine_exchange[kind]))
        for name in names:
            relative = _relative_array_error(
                coarse_exchange[kind].get(name, 0.0),
                fine_exchange[kind].get(name, 0.0), scale_floor=1.0)
            local[name] = relative
            maximum = max(maximum, relative)
        inventories[kind] = local

    displacement_floor = 64.0 * np.finfo(float).eps * float(dx_m)
    recession = _relative_array_error(
        coarse_recession_m, fine_recession_m,
        scale_floor=displacement_floor)
    growth = _relative_array_error(
        coarse_growth_m, fine_growth_m,
        scale_floor=displacement_floor)
    maximum = max(maximum, recession, growth)
    return {
        "maximum_relative_error": maximum,
        "by_state_field_increment": state_fields,
        "by_exchange_inventory": inventories,
        "per_face_integrated_recession_relative_error": recession,
        "per_face_integrated_growth_relative_error": growth,
    }


def _update_embedded_error_maxima(target, error):
    target["maximum_relative_error"] = max(
        target["maximum_relative_error"], error["maximum_relative_error"])
    for name, value in error["by_state_field_increment"].items():
        target["by_state_field_increment"][name] = max(
            target["by_state_field_increment"].get(name, 0.0), value)
    for kind, values in error["by_exchange_inventory"].items():
        for name, value in values.items():
            target["by_exchange_inventory"][kind][name] = max(
                target["by_exchange_inventory"][kind].get(name, 0.0), value)
    for name in (
            "per_face_integrated_recession_relative_error",
            "per_face_integrated_growth_relative_error"):
        target[name] = max(target[name], error[name])


def co_integrate_frozen_radiosity_chemistry(
        state, direct_surface_fluxes, full_face_area_m2, factors,
        active_face_index, active_material, species_role, mechanism, *,
        horizon_s, maximum_local_dp, minimum_substep_s,
        relative_tolerance, maximum_iterations, dx_m):
    """Adaptively co-integrate chemistry and radiosity without a geometry path."""
    horizon_s = float(horizon_s)
    minimum_substep_s = float(minimum_substep_s)
    dx_m = float(dx_m)
    if (not np.isfinite(horizon_s) or horizon_s <= 0.0
            or not np.isfinite(minimum_substep_s) or minimum_substep_s <= 0.0
            or 2.0 * minimum_substep_s > horizon_s
            or not np.isfinite(dx_m) or dx_m <= 0.0):
        raise ValueError("invalid co-integration horizon, minimum substep, or dx")
    state_before = {
        name: np.asarray(values).copy() for name, values in state.fields.items()}
    flux_hash = surface_flux_sha256(direct_surface_fluxes)
    current = state
    elapsed = 0.0
    trial_dt = horizon_s
    accepted = 0
    rejected = 0
    minimum_accepted = np.inf
    maximum_accepted = 0.0
    maximum_dp_seen = 0.0
    maximum_balance_seen = 0.0
    maximum_embedded_error_seen = 0.0
    accepted_embedded_error_maxima = {
        "maximum_relative_error": 0.0,
        "by_state_field_increment": {},
        "by_exchange_inventory": {
            kind: {} for kind in (
                "removed", "outgoing", "unresolved", "deposited")},
        "per_face_integrated_recession_relative_error": 0.0,
        "per_face_integrated_growth_relative_error": 0.0,
    }
    rejection_history = []
    exchange_total = _exchange_accumulator()
    recession_displacement = np.zeros(np.asarray(active_material).shape)
    growth_displacement = np.zeros(np.asarray(active_material).shape)
    maximum_step_ledger_residual = 0.0
    roundoff_s = max(np.finfo(float).eps * horizon_s * 8.0, 1.0e-18)
    while True:
        remaining = horizon_s - elapsed
        if remaining <= roundoff_s:
            elapsed = horizon_s
            break
        minimum_macrostep_s = 2.0 * minimum_substep_s
        if remaining < minimum_macrostep_s:
            raise CoIntegrationRefusal(
                "final embedded remainder falls below two minimum substeps: "
                f"remaining={remaining:.9g}, minimum={minimum_substep_s:.9g}")
        trial_dt = min(trial_dt, remaining)
        if 0.0 < remaining - trial_dt < minimum_macrostep_s:
            trial_dt = remaining
        failure = None
        try:
            start_flux, start_diagnostics, probability_start = solve_cached_radiosity(
                direct_surface_fluxes, full_face_area_m2, factors, current, mechanism,
                active_face_index, active_material,
                relative_tolerance=relative_tolerance,
                maximum_iterations=maximum_iterations)
            active_start_flux = _select_surface_fluxes(
                start_flux, active_face_index, len(full_face_area_m2), species_role)
            half_dt = 0.5 * trial_dt
            coarse = mechanism.advance_by_material(
                current, active_start_flux, trial_dt, active_material)
            fine_first = mechanism.advance_by_material(
                current, active_start_flux, half_dt, active_material)
            midpoint_flux, midpoint_diagnostics, probability_midpoint = (
                solve_cached_radiosity(
                    direct_surface_fluxes, full_face_area_m2, factors,
                    fine_first.state, mechanism, active_face_index, active_material,
                    relative_tolerance=relative_tolerance,
                    maximum_iterations=maximum_iterations))
            active_midpoint_flux = _select_surface_fluxes(
                midpoint_flux, active_face_index, len(full_face_area_m2), species_role)
            fine_second = mechanism.advance_by_material(
                fine_first.state, active_midpoint_flux, half_dt, active_material)
            probability_final = _active_probability(
                mechanism, fine_second.state, active_material)
            dp_start_mid, by_species_start_mid = _maximum_probability_change(
                probability_start, probability_midpoint)
            dp_mid_final, by_species_mid_final = _maximum_probability_change(
                probability_midpoint, probability_final)
            dp_start_final, by_species_start_final = _maximum_probability_change(
                probability_start, probability_final)
            maximum_dp = max(dp_start_mid, dp_mid_final, dp_start_final)
            maximum_balance = max(
                _maximum_balance(start_diagnostics),
                _maximum_balance(midpoint_diagnostics))
            ledgers = [
                float(_maximum_ledger_residual(item.material_exchange))
                for item in (coarse, fine_first, fine_second)]
            ledger = max(ledgers)
            coarse_exchange = _exchange_from_steps(coarse)
            fine_exchange = _exchange_from_steps(fine_first, fine_second)
            coarse_recession = (
                np.asarray(coarse.etch_velocity_m_s, dtype=float) * trial_dt)
            coarse_growth = (
                np.asarray(coarse.normal_growth_velocity_m_s, dtype=float) * trial_dt)
            fine_recession = half_dt * (
                np.asarray(fine_first.etch_velocity_m_s, dtype=float)
                + np.asarray(fine_second.etch_velocity_m_s, dtype=float))
            fine_growth = half_dt * (
                np.asarray(fine_first.normal_growth_velocity_m_s, dtype=float)
                + np.asarray(fine_second.normal_growth_velocity_m_s, dtype=float))
            embedded = _embedded_step_error(
                current, coarse, fine_second, coarse_exchange, fine_exchange,
                coarse_recession, fine_recession, coarse_growth, fine_growth,
                dx_m=dx_m)
            if maximum_dp > float(maximum_local_dp):
                failure = {
                    "kind": "reaction_probability_change",
                    "maximum_local_dp": maximum_dp,
                    "start_to_midpoint": by_species_start_mid,
                    "midpoint_to_final": by_species_mid_final,
                    "start_to_final": by_species_start_final,
                }
            elif (embedded["maximum_relative_error"]
                  > GATES["maximum_embedded_relative_error"]):
                failure = {
                    "kind": "embedded_step_error",
                    "embedded_error": embedded,
                }
            elif maximum_balance > GATES["maximum_radiosity_relative_balance_error"]:
                failure = {
                    "kind": "radiosity_balance",
                    "maximum_relative_balance_error": maximum_balance,
                }
            elif ledger != 0.0 or _cumulative_ledger_residual(fine_exchange) != 0.0:
                failure = {
                    "kind": "material_ledger",
                    "maximum_step_residual": ledger,
                    "fine_cumulative_residual": _cumulative_ledger_residual(
                        fine_exchange),
                }
            elif not all(item.validity.within_declared_scope for item in (
                    coarse, fine_first, fine_second)):
                failure = {
                    "kind": "mechanism_scope",
                    "reasons": sorted(set(
                        str(reason) for item in (coarse, fine_first, fine_second)
                        for reason in item.validity.reasons)),
                }
        except (DiffuseNeutralNoSinkError, RuntimeError, ValueError) as error:
            failure = {"kind": "mechanism_exception", "reason": str(error)}
            maximum_dp = np.inf
            ledger = np.inf
            embedded = {"maximum_relative_error": np.inf}
        if failure is not None:
            rejected += 1
            rejection_history.append({
                "elapsed_s": elapsed, "rejected_dt_s": trial_dt, **failure})
            reduced = 0.5 * trial_dt
            if reduced < minimum_macrostep_s:
                raise CoIntegrationRefusal(
                    f"substep failed above minimum dt: elapsed={elapsed:.9g}, "
                    f"dt={trial_dt:.9g}, minimum={minimum_substep_s:.9g}, "
                    f"failure={failure['kind']}")
            trial_dt = reduced
            continue
        _accumulate_exchange(exchange_total, fine_first.material_exchange)
        _accumulate_exchange(exchange_total, fine_second.material_exchange)
        recession_displacement += fine_recession
        growth_displacement += fine_growth
        maximum_step_ledger_residual = max(maximum_step_ledger_residual, ledger)
        maximum_dp_seen = max(maximum_dp_seen, maximum_dp)
        maximum_balance_seen = max(maximum_balance_seen, maximum_balance)
        maximum_embedded_error_seen = max(
            maximum_embedded_error_seen, embedded["maximum_relative_error"])
        _update_embedded_error_maxima(
            accepted_embedded_error_maxima, embedded)
        current = fine_second.state
        elapsed += trial_dt
        accepted += 1
        minimum_accepted = min(minimum_accepted, half_dt)
        maximum_accepted = max(maximum_accepted, half_dt)
        if (elapsed < horizon_s - roundoff_s
                and maximum_dp < 0.25 * float(maximum_local_dp)
                and embedded["maximum_relative_error"]
                < 0.25 * GATES["maximum_embedded_relative_error"]):
            trial_dt = min(2.0 * trial_dt, horizon_s - elapsed)
    if any(
            not np.array_equal(values, state.fields[name])
            for name, values in state_before.items()):
        raise RuntimeError("co-integration mutated its input surface state")
    if surface_flux_sha256(direct_surface_fluxes) != flux_hash:
        raise RuntimeError("co-integration mutated its direct transport flux")
    area = np.asarray(full_face_area_m2, dtype=float)[active_face_index]
    integrated = {
        kind: _integrated_inventory(values, area)
        for kind, values in exchange_total.items()
    }
    oxide_units = float(integrated["removed"].get("SiO2_formula_unit", 0.0))
    density = float(mechanism.mechanisms[1].parameters.bulk_formula_density_m3)
    oxide_area = float(np.sum(area[np.asarray(active_material) == 1]))
    oxide_volume = oxide_units / density
    return {
        "state": current,
        "horizon_s": horizon_s,
        "accepted_substeps": accepted,
        "rejected_substeps": rejected,
        "minimum_accepted_substep_s": float(minimum_accepted),
        "maximum_accepted_substep_s": maximum_accepted,
        "maximum_local_dp_seen": maximum_dp_seen,
        "maximum_embedded_relative_error_seen": maximum_embedded_error_seen,
        "accepted_embedded_error_maxima": accepted_embedded_error_maxima,
        "maximum_radiosity_balance_seen": maximum_balance_seen,
        "rejection_history": rejection_history,
        "fixed_direct_flux_sha256": flux_hash,
        "integrated_exchange": integrated,
        "per_face_integrated_exchange_units_m2": {
            kind: {
                name: np.asarray(values).copy()
                for name, values in sorted(inventory.items())
            }
            for kind, inventory in sorted(exchange_total.items())
        },
        "per_face_integrated_exchange_sha256": per_face_exchange_sha256(
            exchange_total),
        "maximum_step_ledger_residual_units_m2": maximum_step_ledger_residual,
        "maximum_cumulative_ledger_residual_units_m2": _cumulative_ledger_residual(
            exchange_total),
        "oxide_removal": {
            "integrated_formula_units": oxide_units,
            "integrated_volume_m3": oxide_volume,
            "effective_mean_normal_velocity_m_s": (
                oxide_volume / oxide_area / horizon_s if oxide_area > 0.0 else 0.0),
        },
        "effective_velocity": {
            "maximum_gross_velocity_m_s": float(
                np.max(recession_displacement + growth_displacement) / horizon_s),
            "maximum_absolute_net_velocity_m_s": float(
                np.max(np.abs(recession_displacement - growth_displacement)) / horizon_s),
        },
        "displacement": {
            "maximum_gross_displacement_m": float(
                np.max(recession_displacement + growth_displacement)),
            "maximum_gross_displacement_dx": float(
                np.max(recession_displacement + growth_displacement) / dx_m),
            "per_face_integrated_recession_m": recession_displacement,
            "per_face_integrated_growth_m": growth_displacement,
        },
    }


def _evaluate_with_captured_radiosity(
        source, config, *, seed, rays_per_face=RADIOSITY["rays_per_face"]):
    """Run one production q3 operator and retain its exact pre-radiosity transport.

    Periodic first-hit transport is enabled by the production radiosity option.  Running a
    nominally "direct-only" feature step would silently disable periodic ray wrapping, so the
    direct transport must be captured at the production wrapper boundary instead.
    """
    captured = []
    captured_base_transport = []
    estimator = feature_step_module.estimate_diffuse_form_factors_3d
    production_apply = feature_step_module._apply_diffuse_neutral_transport

    def capture_factors(*args, **kwargs):
        factors = estimator(*args, **kwargs)
        captured.append(factors)
        return factors

    def capture_base_transport(transport, *args, **kwargs):
        captured_base_transport.append(transport)
        return production_apply(transport, *args, **kwargs)

    with patch.object(
            feature_step_module, "estimate_diffuse_form_factors_3d",
            side_effect=capture_factors), patch.object(
                feature_step_module, "_apply_diffuse_neutral_transport",
                side_effect=capture_base_transport):
        result, boundary, reported_wall = _evaluate(
            source["geometry"], source["state"], source["fingerprint"],
            boundary_mode=TRANSPORT_OPERATOR["boundary_mode"],
            ion_bins=(TRANSPORT_OPERATOR["ion_energy_bin_eV"],
                      TRANSPORT_OPERATOR["ion_angle_bin_deg"]),
            face_points=TRANSPORT_OPERATOR["ballistic_face_quadrature_points"],
            pilot_config=config, radiosity_rays=int(rays_per_face),
            seed=int(seed), ballistic_transport="face_gather",
            n_position=TRANSPORT_OPERATOR["n_position"], transport_device="cpu")
    if not captured or len(captured_base_transport) != 1:
        raise RuntimeError(
            "production q3 evaluation did not expose exactly one radiosity boundary")
    _verts, faces, _centroids, areas = extract_mesh_3d(
        source["geometry"].phi, source["geometry"].dx)
    active = np.asarray(result.active_face_index, dtype=int)
    if (not np.array_equal(active, np.arange(len(faces)))
            or not np.array_equal(np.asarray(result.active_face_area), areas)):
        raise RuntimeError("production q3 result and extracted form-factor mesh differ")
    return (
        result,
        boundary,
        reported_wall,
        captured_base_transport[0],
        captured[-1],
        areas * source["geometry"].mesh_length_unit_m ** 2,
        len(captured),
    )


def _facewise_flux_comparison(reference, cached, face_area_m2):
    area = np.asarray(face_area_m2, dtype=float)
    names = sorted(
        set(reference.neutral_flux_m2_s) | set(cached.neutral_flux_m2_s))
    by_species = {}
    maximum_linf = 0.0
    maximum_area_l1 = 0.0
    for name in names:
        expected = np.asarray(reference.neutral_flux_m2_s.get(name, 0.0), dtype=float)
        observed = np.asarray(cached.neutral_flux_m2_s.get(name, 0.0), dtype=float)
        expected, observed = np.broadcast_arrays(expected, observed)
        difference = np.abs(observed - expected)
        linf = float(np.max(difference)) / max(
            float(np.max(np.abs(expected))), 1.0)
        area_l1 = float(np.sum(area * difference)) / max(
            float(np.sum(area * np.abs(expected))), 1.0)
        by_species[name] = {
            "relative_linf_error": linf,
            "area_weighted_relative_l1_error": area_l1,
        }
        maximum_linf = max(maximum_linf, linf)
        maximum_area_l1 = max(maximum_area_l1, area_l1)
    reference_hash = surface_flux_sha256(reference)
    cached_hash = surface_flux_sha256(cached)
    return {
        "declared_norms": (
            "facewise neutral-flux relative Linf and physical-area-weighted relative L1; "
            "denominators floored at one SI rate unit"),
        "reference_surface_flux_sha256": reference_hash,
        "cached_surface_flux_sha256": cached_hash,
        "hash_equal": bool(reference_hash == cached_hash),
        "maximum_relative_linf_error": maximum_linf,
        "maximum_area_weighted_relative_l1_error": maximum_area_l1,
        "by_species": by_species,
    }


def _integrated_species_rates(fluxes, area, active_material):
    output = {}
    all_flux = dict(fluxes.neutral_flux_m2_s)
    all_flux.update({item.name: item.flux_m2_s for item in fluxes.energetic_fluxes})
    for name, values in sorted(all_flux.items()):
        values = np.asarray(values, dtype=float)
        by_material = {}
        for material_id, label in ((1, "sio2"), (2, "amorphous_carbon_mask")):
            selected = np.asarray(active_material) == material_id
            by_material[label] = float(np.sum(area[selected] * values[selected]))
        output[name] = {
            "total": float(np.sum(area * values)),
            "by_material": by_material,
        }
    return output


def certify_cached_initial_radiosity(
        reference_path, checkpoint_path, production_fluxes, cached_fluxes,
        active_area_m2, active_material, production_diagnostics,
        cached_diagnostics, *, expected_transport_seed):
    reference_path = Path(reference_path)
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    if reference.get("status") != "single_cell_complete":
        raise ValueError("one-shot q3 reference is not a completed single-cell audit")
    reference_operator = reference.get("current_operator", {})
    if reference_operator.get("transport_seed") != int(expected_transport_seed):
        raise ValueError(
            "one-shot q3 reference lacks the exact declared transport seed")
    if reference_operator.get("neutral_radiosity_seed") != (
            int(expected_transport_seed) + RADIOSITY["seed_offset"]):
        raise ValueError(
            "one-shot q3 reference lacks the exact declared radiosity seed")
    if reference["checkpoints"]["r19"]["checkpoint_sha256"] != _sha256(checkpoint_path):
        raise ValueError("one-shot q3 reference uses a different checkpoint")
    cell = reference["evaluations"].get("r19_checkpoint__r17_parameters")
    if cell is None:
        raise ValueError("one-shot q3 reference lacks the R17-parameter checkpoint cell")
    observed = _integrated_species_rates(cached_fluxes, active_area_m2, active_material)
    comparisons = {}
    maximum = 0.0
    for name, summary in cell["summary"]["incident_flux_by_species"].items():
        expected = {
            "total": float(summary["total"]["incident_rate_s-1"]),
            "sio2": float(summary["by_material"]["sio2"]["incident_rate_s-1"]),
            "amorphous_carbon_mask": float(
                summary["by_material"]["amorphous_carbon_mask"]["incident_rate_s-1"]),
        }
        actual = {
            "total": observed[name]["total"],
            "sio2": observed[name]["by_material"]["sio2"],
            "amorphous_carbon_mask": observed[name]["by_material"][
                "amorphous_carbon_mask"],
        }
        local = {}
        for scope in expected:
            relative = abs(actual[scope] - expected[scope]) / max(abs(expected[scope]), 1.0)
            local[scope] = {
                "expected_rate_s-1": expected[scope],
                "cached_rate_s-1": actual[scope],
                "relative_error": relative,
            }
            maximum = max(maximum, relative)
        comparisons[name] = local
    production = certify_production_facewise_cache_exact(
        production_fluxes, cached_fluxes, active_area_m2,
        production_diagnostics, cached_diagnostics)
    gates = dict(production["gates"])
    gates["integrated_flux_reproduction"] = bool(
        maximum <= GATES["maximum_cached_reference_relative_error"])
    return {
        "reference_sha256": _sha256(reference_path),
        "transport_seed": int(expected_transport_seed),
        "neutral_radiosity_seed": (
            int(expected_transport_seed) + RADIOSITY["seed_offset"]),
        "maximum_integrated_flux_relative_error": maximum,
        "maximum_radiosity_relative_balance_error": production[
            "maximum_radiosity_relative_balance_error"],
        "facewise_production_replay": production[
            "facewise_production_replay"],
        "comparisons": comparisons,
        "gates": gates,
        "all_gates_pass": bool(all(gates.values())),
    }


def certify_production_facewise_cache_exact(
        production_fluxes, cached_fluxes, active_area_m2,
        production_diagnostics, cached_diagnostics):
    """Certify cache replay without comparing incompatible archived ray levels.

    This is the only certification used for a nondefault form-factor ray level.
    Production and cached solves share that run's exact captured operator, so
    facewise flux hashes must match exactly; toleranced norms and both projectile
    balances remain independent defensive gates.
    """
    facewise = _facewise_flux_comparison(
        production_fluxes, cached_fluxes, active_area_m2)
    maximum_balance = max(
        _maximum_balance(production_diagnostics),
        _maximum_balance(cached_diagnostics))
    gates = {
        "facewise_hash_reproduction": bool(facewise["hash_equal"]),
        "facewise_norm_reproduction": bool(
            facewise["maximum_relative_linf_error"]
            <= GATES["maximum_cached_reference_relative_error"]
            and facewise["maximum_area_weighted_relative_l1_error"]
            <= GATES["maximum_cached_reference_relative_error"]),
        "radiosity_balance": bool(
            maximum_balance <= GATES["maximum_radiosity_relative_balance_error"]),
    }
    return {
        "mode": "nondefault_rays_production_facewise_cache_exact",
        "archived_q3_comparison_performed": False,
        "maximum_radiosity_relative_balance_error": maximum_balance,
        "facewise_production_replay": facewise,
        "gates": gates,
        "all_gates_pass": bool(all(gates.values())),
    }


def _array_manifest_sha256(items):
    """Hash a named array manifest without depending on mapping insertion order."""
    digest = sha256(b"petch.named-array-manifest.v1\0")
    for name, values in sorted(items, key=lambda item: item[0]):
        array = np.ascontiguousarray(np.asarray(values))
        metadata = json.dumps(
            {
                "name": str(name),
                "dtype": array.dtype.str,
                "shape": list(array.shape),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        digest.update(len(metadata).to_bytes(8, "big"))
        digest.update(metadata)
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def surface_state_fields_sha256(state):
    """Return an order-independent, byte-exact digest of all surface-state fields."""
    return _array_manifest_sha256(state.fields.items())


def per_face_exchange_sha256(exchange):
    """Return an exact digest of every named per-face material-exchange array."""
    return _array_manifest_sha256(
        (f"{kind}/{name}", values)
        for kind, inventory in exchange.items()
        for name, values in inventory.items()
    )


def _public_integration(result):
    public = {key: value for key, value in result.items() if key != "state"}
    state = result["state"]
    public["final_state_fields"] = {
        name: np.asarray(values).copy()
        for name, values in sorted(state.fields.items())
    }
    public["final_state_fields_sha256"] = surface_state_fields_sha256(state)
    return public


def evaluate_co_integrated_horizon(
        state, direct_surface_fluxes, full_face_area_m2, factors,
        active_face_index, active_material, species_role, mechanism, *,
        horizon_s, minimum_substep_s, relative_tolerance,
        maximum_iterations, dx_m):
    """Run the declared probability-tolerance pair and apply frozen-state gates."""
    common = dict(
        state=state,
        direct_surface_fluxes=direct_surface_fluxes,
        full_face_area_m2=full_face_area_m2,
        factors=factors,
        active_face_index=active_face_index,
        active_material=active_material,
        species_role=species_role,
        mechanism=mechanism,
        horizon_s=horizon_s,
        minimum_substep_s=minimum_substep_s,
        relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations,
        dx_m=dx_m,
    )
    nominal = co_integrate_frozen_radiosity_chemistry(
        maximum_local_dp=GATES["maximum_local_reaction_probability_change"],
        **common)
    tight = co_integrate_frozen_radiosity_chemistry(
        maximum_local_dp=GATES[
            "tight_maximum_local_reaction_probability_change"],
        **common)
    left = nominal["oxide_removal"]["integrated_formula_units"]
    right = tight["oxide_removal"]["integrated_formula_units"]
    relative = abs(left - right) / max(abs(left), abs(right), 1.0)
    exact_ledger = all(
        item["maximum_step_ledger_residual_units_m2"] == 0.0
        and item["maximum_cumulative_ledger_residual_units_m2"] == 0.0
        for item in (nominal, tight))
    maximum_displacement_dx = max(
        item["displacement"]["maximum_gross_displacement_dx"]
        for item in (nominal, tight))
    maximum_balance = max(
        item["maximum_radiosity_balance_seen"] for item in (nominal, tight))
    gates = {
        "embedded_step_doubling": bool(all(
            item["maximum_embedded_relative_error_seen"]
            <= GATES["maximum_embedded_relative_error"]
            for item in (nominal, tight))),
        "local_probability_safety": bool(
            nominal["maximum_local_dp_seen"]
            <= GATES["maximum_local_reaction_probability_change"]
            and tight["maximum_local_dp_seen"]
            <= GATES["tight_maximum_local_reaction_probability_change"]),
        "tolerance_halving_agreement": bool(
            relative <= GATES[
                "maximum_tolerance_halving_relative_oxide_error"]),
        "exact_material_ledger": bool(exact_ledger),
        "gross_displacement_within_frozen_geometry_limit": bool(
            maximum_displacement_dx <= GATES["maximum_gross_displacement_dx"]),
        "radiosity_balance_within_tolerance": bool(
            maximum_balance <= GATES[
                "maximum_radiosity_relative_balance_error"]),
    }
    return {
        "horizon_s": float(horizon_s),
        "nominal": _public_integration(nominal),
        "tight": _public_integration(tight),
        "tolerance_halving_relative_oxide_error": relative,
        "maximum_gross_displacement_dx": maximum_displacement_dx,
        "maximum_radiosity_relative_balance_error": maximum_balance,
        "gates": gates,
        "all_gates_pass": bool(all(gates.values())),
    }


def run(args):
    started = perf_counter()
    destination = Path(args.output)
    requested_rays_per_face = int(args.rays_per_face)
    requested_horizon_fractions = tuple(
        float(value) for value in args.requested_horizon_fractions)
    selected_horizon_fractions = tuple(
        float(value) for value in args.horizon_fractions)
    maximum_horizon_fraction = float(args.maximum_horizon_fraction)
    certification_mode = _certification_mode_for_rays(requested_rays_per_face)
    radiosity_operator = {
        **RADIOSITY,
        "rays_per_face": requested_rays_per_face,
        # This is the maximum used by the production helper when its sampled
        # graph forces nested no-sink refinement.
        "maximum_rays_per_face": 8 * requested_rays_per_face,
    }
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
        raise ValueError("co-integration requires the completed 60 s R19 checkpoint")
    dt_next = float(source["checkpoint_metadata"]["next_step_duration_s"])
    dx_m = float(source["geometry"].dx * source["geometry"].mesh_length_unit_m)
    payload = {
        "schema": SCHEMA,
        "status": "running",
        "scientific_scope": (
            "diagnostic-only frozen-geometry co-integration of exact routed surface chemistry "
            "and cached diffuse-neutral radiosity over one checkpoint next-step horizon"),
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
        },
        "transport_operator": {
            **TRANSPORT_OPERATOR,
            "transport_seed": int(args.seed),
            "neutral_radiosity_seed": (
                int(args.seed) + radiosity_operator["seed_offset"]),
        },
        "radiosity_operator": radiosity_operator,
        "form_factor_certification_mode": certification_mode,
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
            "maximum_direct_transport_wall_s": float(args.maximum_transport_wall_s),
            "maximum_total_wall_s": float(args.maximum_total_wall_s),
            "next_profile_step_s": dt_next,
            # Retain the original artifact key.  At default CLI settings its
            # value is byte-for-byte unchanged.
            "horizon_fractions": list(selected_horizon_fractions),
            "requested_horizon_fractions": list(requested_horizon_fractions),
            "maximum_horizon_fraction": maximum_horizon_fraction,
            "selected_horizon_fractions": list(selected_horizon_fractions),
        },
        "form_factor_certification": None,
        "horizons": [],
        "largest_common_passing_horizon": None,
        "provenance": {
            "source": _hash_manifest(SOURCE_PATHS),
            "base_inputs": _hash_manifest(BASE_INPUT_PATHS),
            "runtime_selection": {
                "requested_rays_per_face": requested_rays_per_face,
                "certification_mode": certification_mode,
                "requested_horizon_fractions": list(requested_horizon_fractions),
                "maximum_horizon_fraction": maximum_horizon_fraction,
                "selected_horizon_fractions": list(selected_horizon_fractions),
            },
        },
    }
    config = _operator_config(source["config"], PARAMETERS["r17"])
    config["radiosity_enabled"] = True
    snapshot = _snapshot_inputs(source["geometry"], source["state"])
    direct_started = perf_counter()
    try:
        with _hard_deadline(min(
                float(args.maximum_transport_wall_s),
                float(args.maximum_total_wall_s) - (perf_counter() - started))):
            (production_result, boundary, direct_reported_wall, base_transport,
             factors, full_area, form_factor_estimate_count) = (
                _evaluate_with_captured_radiosity(
                    source, config, seed=int(args.seed),
                    rays_per_face=requested_rays_per_face))
    except EvaluationDeadlineExceeded as error:
        payload["status"] = "bounded_direct_transport_timeout"
        payload["execution_budget"]["timeout_reason"] = str(error)
        payload["execution_budget"]["timeout_interpretation"] = {
            "classification": "implementation_or_controller_evidence_only",
            "physics_conclusion_permitted": False,
        }
        payload["total_wall_time_s"] = float(perf_counter() - started)
        _write_json_atomic(destination, payload)
        return payload
    direct_wall = perf_counter() - direct_started
    if not _inputs_unchanged(snapshot, source["geometry"], source["state"]):
        raise RuntimeError("direct transport mutated checkpoint inputs")
    role = {
        species.name: (
            "energetic_bombardment" if species.charge_number != 0
            else "neutral_reactant")
        for species in boundary.species
    }
    active = np.asarray(production_result.active_face_index, dtype=int)
    active_material = np.asarray(
        production_result.face_material_id, dtype=int)[active]
    direct_flux = base_transport.surface_fluxes
    production_initial = production_result.transport.surface_fluxes
    production_diagnostics = production_result.diagnostics.get(
        "neutral_radiosity", {})
    payload["direct_transport"] = {
        "wall_time_s": float(direct_wall),
        "reported_evaluator_wall_time_s": float(direct_reported_wall),
        "direct_surface_flux_sha256": surface_flux_sha256(direct_flux),
        "production_surface_flux_sha256": surface_flux_sha256(production_initial),
        "boundary_provenance": _jsonable(boundary.provenance),
        "input_checkpoint_unchanged": True,
        "periodic_first_hit_captured_before_radiosity": True,
    }
    _write_json_atomic(destination, payload)
    try:
        remaining = float(args.maximum_total_wall_s) - (perf_counter() - started)
        with _hard_deadline(remaining):
            r17_mechanism = build_krueger_2024_material_router_3d(
                **PARAMETERS["r17"])
            cached_initial, cached_diagnostics, _ = solve_cached_radiosity(
                direct_flux, full_area, factors, source["state"], r17_mechanism,
                active, active_material,
                relative_tolerance=radiosity_operator["relative_tolerance"],
                maximum_iterations=radiosity_operator["maximum_iterations"])
            if certification_mode == (
                    "archived_q3_8ray_plus_production_facewise_cache_exact"):
                if int(factors.rays_per_face) != RADIOSITY["rays_per_face"]:
                    raise ValueError(
                        "the archived q3 gate requires an actual 8-ray production "
                        f"operator, not {factors.rays_per_face} rays per face")
                certification = certify_cached_initial_radiosity(
                    args.q3_reference, source["checkpoint_path"], production_initial,
                    cached_initial, full_area[active], active_material,
                    production_diagnostics, cached_diagnostics,
                    expected_transport_seed=int(args.seed))
            else:
                certification = certify_production_facewise_cache_exact(
                    production_initial, cached_initial, full_area[active],
                    production_diagnostics, cached_diagnostics)
    except EvaluationDeadlineExceeded as error:
        payload["status"] = "bounded_form_factor_timeout"
        payload["execution_budget"]["timeout_reason"] = str(error)
        payload["execution_budget"]["timeout_interpretation"] = {
            "classification": "implementation_or_controller_evidence_only",
            "physics_conclusion_permitted": False,
        }
        payload["total_wall_time_s"] = float(perf_counter() - started)
        _write_json_atomic(destination, payload)
        return payload
    except (DiffuseNeutralNoSinkError, RuntimeError, ValueError) as error:
        payload["status"] = "form_factor_cache_certification_failure"
        payload["form_factor_certification"] = {
            "all_gates_pass": False,
            "reason": str(error),
            "exception_type": type(error).__name__,
        }
        payload["total_wall_time_s"] = float(perf_counter() - started)
        _write_json_atomic(destination, payload)
        return payload
    payload["form_factors"] = {
        "sha256": form_factor_sha256(factors),
        "face_count": factors.face_count,
        "requested_rays_per_face": requested_rays_per_face,
        "rays_per_face": factors.rays_per_face,
        "production_estimate_count": form_factor_estimate_count,
    }
    payload["form_factor_certification"] = certification
    if not _inputs_unchanged(snapshot, source["geometry"], source["state"]):
        raise RuntimeError("radiosity certification mutated checkpoint inputs")
    if (surface_flux_sha256(direct_flux)
            != payload["direct_transport"]["direct_surface_flux_sha256"]):
        raise RuntimeError("radiosity certification mutated direct transport flux")
    payload["form_factor_certification"]["input_checkpoint_unchanged"] = True
    payload["form_factor_certification"]["direct_flux_unchanged"] = True
    _write_json_atomic(destination, payload)
    if not certification["all_gates_pass"]:
        payload["status"] = "form_factor_cache_certification_failure"
        payload["total_wall_time_s"] = float(perf_counter() - started)
        _write_json_atomic(destination, payload)
        return payload

    largest = None
    for fraction in selected_horizon_fractions:
        horizon = dt_next * fraction
        entry = {
            "fraction_of_next_profile_step": fraction,
            "horizon_s": horizon,
            "parameter_results": {},
            "common_pass": False,
        }
        for label in ("r17", "r19"):
            mechanism = build_krueger_2024_material_router_3d(**PARAMETERS[label])
            try:
                remaining = float(args.maximum_total_wall_s) - (perf_counter() - started)
                with _hard_deadline(remaining):
                    result = evaluate_co_integrated_horizon(
                        source["state"], direct_flux, full_area, factors, active,
                        active_material, role, mechanism, horizon_s=horizon,
                        minimum_substep_s=GATES["minimum_substep_s"],
                        relative_tolerance=radiosity_operator["relative_tolerance"],
                        maximum_iterations=radiosity_operator["maximum_iterations"],
                        dx_m=dx_m)
            except EvaluationDeadlineExceeded as error:
                entry["first_failure"] = {
                    "parameter_label": label,
                    "kind": "bounded_total_timeout",
                    "reason": str(error),
                    "classification": "implementation_or_controller_evidence_only",
                    "physics_conclusion_permitted": False,
                }
                payload["status"] = "bounded_total_timeout"
                break
            except CoIntegrationRefusal as error:
                entry["parameter_results"][label] = {
                    "parameters": PARAMETERS[label],
                    "status": "refused",
                    "reason": str(error),
                }
                entry["first_failure"] = {
                    "parameter_label": label,
                    "kind": "co_integration_refusal",
                    "reason": str(error),
                }
                payload["status"] = "stopped_at_first_failed_horizon"
                break
            result["parameters"] = PARAMETERS[label]
            entry["parameter_results"][label] = result
            if not result["all_gates_pass"]:
                entry["first_failure"] = {
                    "parameter_label": label,
                    "kind": "co_integration_gate_failure",
                    "gates": result["gates"],
                }
                payload["status"] = "stopped_at_first_failed_horizon"
                break
        if len(entry["parameter_results"]) == 2 and all(
                item.get("all_gates_pass", False)
                for item in entry["parameter_results"].values()):
            entry["common_pass"] = True
            r17 = entry["parameter_results"]["r17"]["tight"]["oxide_removal"]
            r19 = entry["parameter_results"]["r19"]["tight"]["oxide_removal"]
            difference = (
                r19["integrated_formula_units"] - r17["integrated_formula_units"])
            entry["paired_oxide_removal_direction"] = {
                "r19_minus_r17_integrated_formula_units": difference,
                "r19_to_r17_ratio": (
                    r19["integrated_formula_units"] / r17["integrated_formula_units"]
                    if r17["integrated_formula_units"] > 0.0 else None),
                "direction": (
                    "r19_lower" if difference < 0.0
                    else "r19_higher" if difference > 0.0 else "equal"),
            }
            largest = entry
        payload["horizons"].append(entry)
        _write_json_atomic(destination, payload)
        if not entry["common_pass"]:
            if payload["status"] == "running":
                payload["status"] = "stopped_at_first_failed_horizon"
            break
    else:
        payload["status"] = "pass"
    if largest is not None:
        payload["largest_common_passing_horizon"] = {
            "fraction_of_next_profile_step": largest[
                "fraction_of_next_profile_step"],
            "horizon_s": largest["horizon_s"],
            "paired_oxide_removal_direction": largest[
                "paired_oxide_removal_direction"],
        }
    payload["total_wall_time_s"] = float(perf_counter() - started)
    _write_json_atomic(destination, payload)
    print(json.dumps({
        "status": payload["status"],
        "largest_common_passing_horizon": payload[
            "largest_common_passing_horizon"],
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
        "--q3-reference",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "frozen_checkpoint_q3_seed241_current" / "audit.json"))
    parser.add_argument(
        "--output",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "frozen_radiosity_chemistry" / "audit.json"))
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument(
        "--rays-per-face", type=_power_of_two_rays,
        default=RADIOSITY["rays_per_face"],
        help=(
            "initial production/cached diffuse form-factor rays per face; "
            "nondefault levels use production-facewise cache certification only"))
    parser.add_argument(
        "--horizon-fractions", type=_unit_horizon_fraction, nargs="+",
        default=HORIZON_FRACTIONS,
        metavar="FRACTION",
        help="strictly increasing fractions of the checkpoint's next profile step")
    parser.add_argument(
        "--maximum-horizon-fraction", type=_unit_horizon_fraction,
        default=1.0,
        help="truncate the requested horizon ladder at this fraction")
    parser.add_argument("--maximum-transport-wall-s", type=float, default=180.0)
    parser.add_argument("--maximum-total-wall-s", type=float, default=480.0)
    args = parser.parse_args(argv)
    requested = tuple(float(value) for value in args.horizon_fractions)
    try:
        selected = _select_horizon_fractions(
            requested, args.maximum_horizon_fraction)
    except ValueError as error:
        parser.error(str(error))
    args.requested_horizon_fractions = requested
    args.horizon_fractions = selected
    return args


if __name__ == "__main__":
    run(parse_args())
