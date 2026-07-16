#!/usr/bin/env python3
"""Run the fixed-duration Jeong 2023 transfer matrix through the unified 3-D engine.

This is an experiment adapter, not a second feature solver.  It consumes the same boundary,
surface-kinetics, transport, and profile-evolution APIs as the Jeon 2022 campaign.  Exactly one
experimental marker is eligible for magnitude calibration; every other marker is scored held out.
The first campaign stage is deliberately charge-off.  A miss confined to the 60 nm etch-stop rows
is the preregistered causality gate for the physical-time charging path, not permission to tune the
chemistry against those rows.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np

from petch.charged_surface_response_3d import GrazingSpecularIonReflection3D
from petch.experimental_boundary import (
    Jeong2023IonBoundaryClosure,
    build_jeong_2023_boundary_state,
)
from petch.experimental_data import (
    load_jeong_2023_etch_depths,
    load_jeong_2023_radical_densities,
)
from petch.feature_step_3d import make_rectangular_trench_geometry_3d, solve_feature_3d
from petch.fluorocarbon_lamagna import (
    LaMagnaFluorocarbonParameters, LaMagnaGarozzoFluorocarbonMechanism,
)
from petch.surface_kinetics import (
    ReducedSiO2FluorocarbonMechanism, ReducedSiO2FluorocarbonParameters,
)

from jeon_unified_baseline import _floor_height, _write_json_atomic, baseline_mechanism


ROOT = Path(__file__).parents[1]
DATA = ROOT / "data" / "experimental" / "jeong_2023"


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git_revision():
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _git_worktree_dirty():
    try:
        return bool(subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        return None


def _canonical_input_deck(args):
    """Return the replay deck and a content hash independent of output location."""
    deck = {}
    for key, value in sorted(vars(args).items()):
        if key in {"output", "quiet"}:
            continue
        if isinstance(value, Path):
            value = str(value)
        deck[key] = value
    encoded = json.dumps(
        deck, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return deck, hashlib.sha256(encoded).hexdigest()


def _campaign_provenance(args):
    deck, config_hash = _canonical_input_deck(args)
    evidence_files = [
        DATA / "README.md",
        DATA / "digitized_figure6_radicals.csv",
        DATA / "digitized_figure7_depths.csv",
    ]
    if args.ion_closure_json is not None:
        evidence_files.append(args.ion_closure_json)
    implementation_files = (
        Path(__file__),
        ROOT / "src" / "petch" / "experimental_boundary.py",
        ROOT / "src" / "petch" / "surface_kinetics.py",
        ROOT / "src" / "petch" / "feature_step_3d.py",
    )
    return {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_hash_sha256": config_hash,
        "input_deck": deck,
        "evidence_checksums_sha256": {
            str(path.relative_to(ROOT)): _sha256(path) for path in evidence_files},
        "implementation_checksums_sha256": {
            str(path.resolve().relative_to(ROOT.resolve())): _sha256(path)
            for path in implementation_files},
        "git_revision": _git_revision(),
        "git_worktree_dirty_at_start": _git_worktree_dirty(),
        "runtime": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        },
    }


def _campaign_payload(args, runs, provenance, expected_runs):
    complete = len(runs) == expected_runs
    transport_only = bool(args.transport_only)
    return {
        "campaign": "jeong_2023_fixed_duration_transfer",
        "stage": (
            "fixed_geometry_transport_audit"
            if transport_only else "charge_off_chemistry_transport"),
        "status": "complete" if complete else "partial_resumable_evidence",
        "provenance": provenance,
        "progress": {
            "completed_runs": len(runs),
            "expected_runs": expected_runs,
        },
        "calibration_contract": {
            "eligible_anchor": "ion_energy/200nm/890V only",
            "held_out_points": 17,
            "tuning_on_held_out_forbidden": True,
            "score_eligible": not transport_only,
            "transport_only_runs_are_diagnostic": transport_only,
        },
        "chemistry": {
            "model": args.chemistry_model,
            "lamagna_neutral_transport_mode": args.lamagna_neutral_transport_mode,
            "complex_probability": args.complex_probability,
            "deposition_probability": args.deposition_probability,
            "complex_removal_reaction_order": args.complex_removal_reaction_order,
            "energetic_response_scale": args.energetic_response_scale,
            "bare_reference_yield_at_900eV": args.bare_reference_yield,
            "complex_reference_yield_at_900eV": args.complex_reference_yield,
            "polymer_reference_yield_at_900eV": args.polymer_reference_yield,
            "parameter_source": (
                "Huang_Kushner_2019_reaction_table_plus_one_Jeong_2023_anchor"
                if args.chemistry_model == "huang_2019_reduced"
                else "Jeon_2022_development_screen_then_one_Jeong_2023_anchor"),
        },
        "transport": {
            "ion_reflection_enabled": bool(args.ion_reflection),
            "ion_reflection_material_ids": [1, 2] if args.ion_reflection else [],
            "ballistic_transport": args.ballistic_transport,
        },
        "runs": list(runs),
    }


def _select_control(rows, control_mode, width_nm, control_value):
    if control_mode == "ion_energy":
        matches = [item for item in rows if item.control_mode == control_mode
                   and item.trench_width_nm == width_nm
                   and np.isclose(item.self_bias_magnitude_v, control_value)]
    else:
        matches = [item for item in rows if item.control_mode == control_mode
                   and item.trench_width_nm == width_nm
                   and np.isclose(item.electron_density_m3, control_value)]
    if len(matches) != 1:
        raise ValueError("requested Jeong 2023 control must match exactly one frozen target")
    return matches[0]


def _load_ion_closure(path):
    if path is None:
        return Jeong2023IonBoundaryClosure.all_argon_development()
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {
        "species_mass_amu",
        "normal_energy_fraction_of_self_bias",
        "provenance",
        "supports_prediction_within_declared_domain",
    }
    optional = {
        "species_density_fraction",
        "explicit_species_flux_m2_s",
        "positive_ion_density_over_electron_density",
    }
    if not isinstance(payload, dict) or not required <= set(payload) <= required | optional:
        raise ValueError(
            "ion closure JSON must contain the declared Jeong ion-boundary schema")
    return Jeong2023IonBoundaryClosure(
        species_mass_amu=payload["species_mass_amu"],
        normal_energy_fraction_of_self_bias=payload[
            "normal_energy_fraction_of_self_bias"],
        species_density_fraction=payload.get("species_density_fraction", {}),
        explicit_species_flux_m2_s=payload.get("explicit_species_flux_m2_s", {}),
        positive_ion_density_over_electron_density=payload.get(
            "positive_ion_density_over_electron_density", 1.0),
        provenance=payload["provenance"],
        supports_prediction_within_declared_domain=payload[
            "supports_prediction_within_declared_domain"],
    )


def _weighted_quantile(values, weights, quantiles):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)
    if (values.ndim != 1 or weights.shape != values.shape
            or quantiles.ndim != 1 or np.any(~np.isfinite(values))
            or np.any(~np.isfinite(weights)) or np.any(weights < 0.0)
            or np.any((quantiles < 0.0) | (quantiles > 1.0))):
        raise ValueError("invalid weighted-quantile inputs")
    if not values.size or weights.sum() <= 0.0:
        return [None for _ in quantiles]
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    indices = np.searchsorted(
        cumulative, quantiles * cumulative[-1], side="left")
    indices = np.minimum(indices, len(sorted_values) - 1)
    return [float(sorted_values[index]) for index in indices]


def _reimpact_energy_diagnostics(
        final_step, boundary, realized_domain_um, geometry, ion_species_name):
    """Audit the energy measure generated by the certified reflection cascade.

    The Huang--Kushner reaction table contains direct fluorocarbon-ion polymer deposition in the
    5--70 eV band.  Jeong does not report the positive-ion composition, so this diagnostic measures
    only whether the common transport operator creates that energy support.  It does not assign
    fluorocarbon identity to the engine's current effective ``Ar+`` population.

    Re-impact particle rates use the cascade's authoritative emitted-event weights rather than
    summing face flux densities.  This keeps the diagnostic invariant to triangle area and exactly
    consistent with the cascade particle ledger.
    """
    cascade = final_step.charged_surface_cascade
    if cascade is None:
        return {
            "available": False,
            "reason": "ion reflection disabled",
            "transport_species": ion_species_name,
            "physical_interpretation_limitation": (
                "no reflected energy measure exists when the reflection channel is disabled"),
        }
    source_area_m2 = float(
        realized_domain_um[0] * realized_domain_um[1]
        * geometry.mesh_length_unit_m ** 2)
    boundary_rate_s = float(
        boundary.get(ion_species_name).flux_m2_s * source_area_m2)
    primary_landed_rate_s = float(
        boundary_rate_s * final_step.transport.hit_probability[ion_species_name])
    floor_threshold_mesh = float(
        np.min(final_step.active_face_centroid[:, 2]) + geometry.dx)
    energy_chunks = []
    rate_chunks = []
    floor_chunks = []
    bounce_rows = []
    for bounce_index, flights in enumerate(cascade.flights_by_bounce, start=1):
        bounce_energy = []
        bounce_rate = []
        bounce_floor = []
        emitted_rate_s = 0.0
        landed_rate_s = 0.0
        escaped_rate_s = 0.0
        for flight in flights:
            landed = flight.termination == 1
            landed_weights = flight.emitted.event_rate_s[landed]
            if landed_weights.shape != flight.incident.event_energy_eV.shape:
                raise RuntimeError(
                    "charged-cascade landed lineage does not align with incident energies")
            emitted_rate_s += flight.emitted_rate_s
            landed_rate_s += flight.landed_rate_s
            escaped_rate_s += flight.escaped_rate_s
            if not landed_weights.size:
                continue
            positions = flight.incident.event_position
            if positions is None:
                raise RuntimeError(
                    "charged-cascade energy audit requires landed impact positions")
            floor = positions[:, 2] <= floor_threshold_mesh
            bounce_energy.append(flight.incident.event_energy_eV)
            bounce_rate.append(landed_weights)
            bounce_floor.append(floor)
        if bounce_energy:
            energy = np.concatenate(bounce_energy)
            rate = np.concatenate(bounce_rate)
            floor = np.concatenate(bounce_floor)
        else:
            energy = np.empty(0)
            rate = np.empty(0)
            floor = np.empty(0, dtype=bool)
        energy_chunks.append(energy)
        rate_chunks.append(rate)
        floor_chunks.append(floor)

        def band_rate(low, high, *, floor_only=False):
            selected = (energy >= low) & (energy <= high)
            if floor_only:
                selected &= floor
            return float(rate[selected].sum())

        bounce_rows.append({
            "bounce": bounce_index,
            "emitted_rate_s": float(emitted_rate_s),
            "landed_rate_s": float(landed_rate_s),
            "escaped_rate_s": float(escaped_rate_s),
            "landed_rate_5_to_30_eV_s": band_rate(5.0, 30.0),
            "landed_rate_5_to_70_eV_s": band_rate(5.0, 70.0),
            "floor_landed_rate_5_to_70_eV_s": band_rate(
                5.0, 70.0, floor_only=True),
            "landed_energy_quantiles_eV": dict(zip(
                ("q05", "q50", "q95"),
                _weighted_quantile(energy, rate, (0.05, 0.50, 0.95)))),
        })
    energy = np.concatenate(energy_chunks) if energy_chunks else np.empty(0)
    rate = np.concatenate(rate_chunks) if rate_chunks else np.empty(0)
    floor = np.concatenate(floor_chunks) if floor_chunks else np.empty(0, dtype=bool)
    total_reimpact_rate_s = float(rate.sum())

    def aggregate_band(low, high, *, floor_only=False):
        selected = (energy >= low) & (energy <= high)
        if floor_only:
            selected &= floor
        return float(rate[selected].sum())

    low_30 = aggregate_band(5.0, 30.0)
    low_70 = aggregate_band(5.0, 70.0)
    floor_low_70 = aggregate_band(5.0, 70.0, floor_only=True)
    return {
        "available": True,
        "transport_species": ion_species_name,
        "source_boundary_particle_rate_s": boundary_rate_s,
        "primary_landed_particle_rate_s": primary_landed_rate_s,
        "all_reimpact_landed_particle_rate_s": total_reimpact_rate_s,
        "reimpact_landed_over_primary_landed": (
            total_reimpact_rate_s / primary_landed_rate_s
            if primary_landed_rate_s > 0.0 else None),
        "reimpact_landed_rate_5_to_30_eV_s": low_30,
        "reimpact_landed_rate_5_to_70_eV_s": low_70,
        "floor_reimpact_landed_rate_5_to_70_eV_s": floor_low_70,
        "fraction_of_reimpact_rate_5_to_30_eV": (
            low_30 / total_reimpact_rate_s if total_reimpact_rate_s > 0.0 else 0.0),
        "fraction_of_reimpact_rate_5_to_70_eV": (
            low_70 / total_reimpact_rate_s if total_reimpact_rate_s > 0.0 else 0.0),
        "rate_5_to_70_eV_over_primary_landed": (
            low_70 / primary_landed_rate_s if primary_landed_rate_s > 0.0 else None),
        "landed_energy_quantiles_eV": dict(zip(
            ("q05", "q50", "q95"),
            _weighted_quantile(energy, rate, (0.05, 0.50, 0.95)))),
        "floor_definition": {
            "coordinate_system": "input_step_mesh_units",
            "maximum_z": floor_threshold_mesh,
            "band_thickness_mesh_units": float(geometry.dx),
        },
        "by_bounce": bounce_rows,
        "physical_interpretation_limitation": (
            "the v1 cascade propagates the reflected effective positive ion as charged; "
            "Huang--Kushner neutralize a surface-striking ion and propagate its hot-neutral "
            "partner, and Jeong does not report the fluorocarbon-ion fraction"),
    }


def run_target(target, radicals, args):
    width_um = target.trench_width_nm * 1e-3
    pitch_um = args.mask_line_width_um + width_um
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=pitch_um, cell_length=args.cell_length_um,
        domain_height=args.domain_height_um, dx=args.dx_um, opening_width=width_um,
        mask_thickness=args.mask_thickness_um, substrate_top=args.substrate_top_um,
        etched_depth=args.initial_etched_depth_um)
    realized_domain_um = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(realized_domain_um[2])
    radical_support_density = (
        args.energy_radical_support_density_m3
        if target.control_mode == "ion_energy" else target.electron_density_m3)
    selected_radicals = tuple(
        item for item in radicals
        if np.isclose(item.electron_density_m3, radical_support_density))
    radical_channel_mode = {
        "reduced_si_o2": "aggregate",
        "huang_2019_reduced": "huang_2019_reduced",
        "lamagna_garozzo": "heavy_light",
    }[args.chemistry_model]
    boundary = build_jeong_2023_boundary_state(
        target, selected_radicals,
        reference_plane_m=source_z * geometry.mesh_length_unit_m,
        neutral_temperature_K=args.neutral_temperature_k,
        electron_temperature_eV=args.electron_temperature_ev,
        ion_tangential_temperature_eV=args.ion_tangential_temperature_ev,
        n_transverse_ion=args.ion_transverse_quadrature,
        n_transverse_neutral=args.neutral_transverse_quadrature,
        n_normal_neutral=args.neutral_normal_quadrature,
        radical_channel_mode=radical_channel_mode,
        ion_closure=_load_ion_closure(args.ion_closure_json))
    positive_ions = tuple(
        item for item in boundary.species if item.charge_number > 0)
    positive_ion_names = tuple(item.name for item in positive_ions)
    if not positive_ions:
        raise RuntimeError("Jeong boundary contains no positive ions")
    if args.ion_reflection and len(positive_ions) != 1:
        raise ValueError(
            "the v1 certified reflection law is species-specific; a multispecies "
            "Jeong ion boundary must use --no-ion-reflection until each ion has "
            "a declared material-response law")
    if args.chemistry_model == "reduced_si_o2":
        mechanism = baseline_mechanism(
            args.complex_probability, args.deposition_probability,
            args.complex_removal_reaction_order,
            bare_reference_yield=args.bare_reference_yield,
            complex_reference_yield=args.complex_reference_yield,
            polymer_reference_yield=args.polymer_reference_yield)
        species_role = {
            **{name: "energetic_bombardment" for name in positive_ion_names},
            "FC_total": "neutral_reactant",
        }
        mask_reaction_probability = {
            "FC_total": args.mask_reaction_probability}
        surface_fixed_point_tolerance = None
    elif args.chemistry_model == "huang_2019_reduced":
        mechanism = ReducedSiO2FluorocarbonMechanism(
            ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection(
                energetic_response_scale=args.energetic_response_scale))
        species_role = {
            name: "energetic_bombardment" for name in positive_ion_names}
        species_role.update({
            name: "neutral_reactant"
            for name in ("CF", "CF2", "FC_complex_02", "FC_polymer_heavy")})
        mask_reaction_probability = {
            name: args.mask_reaction_probability
            for name in ("CF", "CF2", "FC_complex_02", "FC_polymer_heavy")}
        surface_fixed_point_tolerance = None
    else:
        mechanism = LaMagnaGarozzoFluorocarbonMechanism(
            LaMagnaFluorocarbonParameters.viennaps_4_6_1_reference(
                reference_etchant_flux_m2_s=boundary.get("FC_etchant").flux_m2_s,
                neutral_transport_mode=args.lamagna_neutral_transport_mode))
        species_role = {
            **{name: "energetic_bombardment" for name in positive_ion_names},
            "FC_etchant": "neutral_reactant",
            "FC_polymer": "neutral_reactant",
        }
        mask_reaction_probability = {
            "FC_etchant": args.mask_reaction_probability,
            "FC_polymer": args.mask_reaction_probability,
        }
        surface_fixed_point_tolerance = args.surface_fixed_point_tolerance
    reflection = (
        GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
            (1, 2), ion_species_name=positive_ion_names[0])
        if args.ion_reflection else None)
    reflection_options = None if reflection is None else {
        "fixed_dt": args.reflection_fixed_dt,
        "max_steps": args.reflection_max_steps,
        "trajectory_adaptive_horizon": True,
        "trajectory_emergency_max_steps": args.reflection_emergency_max_steps,
        "max_bounces": args.reflection_max_bounces,
        "relative_tail_tolerance": args.reflection_tail_tolerance,
        "adaptive_bounce_extension": True,
        "emergency_max_bounces": args.reflection_emergency_max_bounces,
        "periodic_lateral": True,
    }
    initial_floor = _floor_height(geometry.phi, geometry.dx)
    simulated_duration_s = 0.0 if args.transport_only else target.etch_duration_s
    simulated_steps = 1 if args.transport_only else args.steps
    started = perf_counter()
    result = solve_feature_3d(
        geometry, boundary, species_role,
        mechanism, etchable_material_ids=(1,), duration_s=simulated_duration_s,
        n_steps=simulated_steps,
        source_bounds=(0.0, float(realized_domain_um[0]),
                       0.0, float(realized_domain_um[1])),
        source_z=source_z, n_position=args.source_positions, seed=args.seed,
        cfl_number=0.3, reinitialize=args.reinitialize, transport_device="cpu",
        neutral_radiosity_options={
            "rays_per_face": args.form_factor_rays, "seed": args.seed + 1000,
            "periodic_lateral": True,
            "domain_size": realized_domain_um,
            "nonetchable_reaction_probability_by_material": {
                2: mask_reaction_probability},
        },
        neutral_surface_fixed_point_tolerance=surface_fixed_point_tolerance,
        neutral_surface_fixed_point_max_iterations=(
            args.surface_fixed_point_max_iterations),
        charged_surface_response=reflection,
        charged_surface_response_options=reflection_options,
        ballistic_transport=args.ballistic_transport,
        ballistic_face_quadrature_points=args.ballistic_face_quadrature_points,
        reinitialization_method=args.reinitialization_method)
    wall_time_s = perf_counter() - started
    final_floor = _floor_height(result.geometry.phi, geometry.dx)
    prediction_nm = float((initial_floor - final_floor) * 1000.0)
    residual_nm = (
        None if args.transport_only else prediction_nm - target.etch_depth_nm)
    final_step = result.steps[-1]
    final_centroid = final_step.next_active_face_centroid
    final_area = final_step.next_active_face_area
    floor_band = final_centroid[:, 2] <= float(np.min(final_centroid[:, 2]) + geometry.dx)

    def weighted_mean(value):
        value = np.asarray(value, dtype=float)
        return float(np.dot(value[floor_band], final_area[floor_band])
                     / final_area[floor_band].sum())

    transport_floor = (
        final_step.active_face_centroid[:, 2]
        <= float(np.min(final_step.active_face_centroid[:, 2]) + geometry.dx))
    transport_area = final_step.active_face_area
    ion_populations = tuple(
        item for item in final_step.transport.surface_fluxes.energetic_fluxes
        if item.name in positive_ion_names)
    if {item.name for item in ion_populations} != set(positive_ion_names):
        raise RuntimeError("transport did not return every declared positive-ion population")
    active_ion_flux = sum(
        (item.flux_m2_s[final_step.active_face_index]
         for item in ion_populations),
        start=np.zeros(len(final_step.active_face_index)))
    boundary_ion_flux = float(sum(item.flux_m2_s for item in positive_ions))
    final_etch_velocity = np.broadcast_to(
        np.asarray(final_step.surface.etch_velocity_m_s, dtype=float),
        (len(final_step.active_face_index),))
    final_growth_velocity = np.broadcast_to(
        np.asarray(getattr(
            final_step.surface, "normal_growth_velocity_m_s", 0.0), dtype=float),
        (len(final_step.active_face_index),))

    def transport_floor_mean(value):
        value = np.asarray(value, dtype=float)
        return float(np.dot(value[transport_floor], transport_area[transport_floor])
                     / transport_area[transport_floor].sum())

    if args.chemistry_model in {"reduced_si_o2", "huang_2019_reduced"}:
        surface_state = {
            "floor_complex_fraction": weighted_mean(
                result.surface_state.complex_fraction),
            "floor_polymer_units_m2": weighted_mean(
                result.surface_state.polymer_units_m2),
            "maximum_polymer_units_m2": float(np.max(
                result.surface_state.polymer_units_m2)),
        }
    else:
        surface_state = {
            "floor_etchant_coverage": weighted_mean(
                result.surface_state.etchant_coverage),
            "floor_polymer_coverage": weighted_mean(
                result.surface_state.polymer_coverage),
            "floor_etchant_on_polymer_coverage": weighted_mean(
                result.surface_state.etchant_on_polymer_coverage),
            "floor_polymer_film_units_m2": weighted_mean(
                result.surface_state.polymer_film_units_m2),
            "maximum_polymer_film_units_m2": float(np.max(
                result.surface_state.polymer_film_units_m2)),
        }
    radical_channel_flux = {
        species.name: species.flux_m2_s for species in boundary.species
        if species.charge_number == 0}

    return {
        "target": {
            "source_figure": target.source_figure,
            "control_mode": target.control_mode,
            "trench_width_nm": target.trench_width_nm,
            "self_bias_magnitude_v": target.self_bias_magnitude_v,
            "electron_density_m3": target.electron_density_m3,
            "etch_duration_s": target.etch_duration_s,
            "etch_depth_nm": target.etch_depth_nm,
            "digitization_interval_nm": [
                target.etch_depth_nm - target.digitization_uncertainty_nm,
                target.etch_depth_nm + target.digitization_uncertainty_nm],
            "measurement_uncertainty_semantics": target.measurement_uncertainty_semantics,
            "split": target.split,
            "role": target.role,
        },
        "prediction": {
            "etch_depth_nm": prediction_nm,
            "residual_nm": residual_nm,
            "within_digitization_interval": (
                None if args.transport_only else bool(
                    abs(residual_nm) <= target.digitization_uncertainty_nm)),
            "score_eligible": not args.transport_only,
            "charging_mode": "off_preregistered_chemistry_transport_stage",
        },
        "boundary": {
            "ion_flux_m2_s": boundary_ion_flux,
            "ion_species_flux_m2_s": {
                item.name: item.flux_m2_s for item in positive_ions},
            "ion_species_mean_energy_eV": {
                item.name: item.mean_energy_eV for item in positive_ions},
            "aggregate_radical_flux_m2_s": float(sum(radical_channel_flux.values())),
            "radical_channel_flux_m2_s": radical_channel_flux,
            "radical_channel_mode": radical_channel_mode,
            "radical_support_electron_density_m3": radical_support_density,
            "provenance": dict(boundary.provenance),
        },
        "chemistry": {
            "model": args.chemistry_model,
            "provenance": dict(mechanism.provenance),
        },
        "surface_state": surface_state,
        "final_floor_transport": {
            "ion_flux_m2_s": transport_floor_mean(active_ion_flux),
            "ion_flux_over_boundary": (
                transport_floor_mean(active_ion_flux) / boundary_ion_flux),
            "etch_velocity_nm_s": 1.0e9 * transport_floor_mean(
                final_etch_velocity),
            "growth_velocity_nm_s": 1.0e9 * transport_floor_mean(
                final_growth_velocity),
            "net_recession_velocity_nm_s": 1.0e9 * transport_floor_mean(
                final_etch_velocity - final_growth_velocity),
            "primary_hit_probability": float(sum(
                item.flux_m2_s
                * final_step.transport.hit_probability[item.name]
                for item in positive_ions) / boundary_ion_flux),
            "primary_escape_probability": float(sum(
                item.flux_m2_s
                * final_step.transport.escape_probability[item.name]
                for item in positive_ions) / boundary_ion_flux),
            "primary_hit_probability_by_species": {
                item.name: final_step.transport.hit_probability[item.name]
                for item in positive_ions},
            "primary_escape_probability_by_species": {
                item.name: final_step.transport.escape_probability[item.name]
                for item in positive_ions},
            "lineage_replay_count": final_step.transport.lineage_replay_count,
            "lineage_replay_fraction": final_step.transport.lineage_replay_fraction,
        },
        "validity": {
            "within_declared_scope": result.validity.within_declared_scope,
            "parameter_evidence_supports_prediction": (
                result.validity.parameter_evidence_supports_prediction),
            "nonpredictive_parameters": list(result.validity.nonpredictive_parameters),
            "known_limitations": list(result.validity.known_limitations),
        },
        "numerics": {
            "dx_um": args.dx_um,
            "requested_steps": args.steps,
            "simulated_steps": simulated_steps,
            "simulated_duration_s": simulated_duration_s,
            "transport_only": bool(args.transport_only),
            "requested_pitch_um": pitch_um,
            "requested_cell_length_um": args.cell_length_um,
            "initial_etched_depth_um": args.initial_etched_depth_um,
            "realized_domain_um": realized_domain_um.tolist(),
            "source_positions": args.source_positions,
            "form_factor_rays": args.form_factor_rays,
            "neutral_surface_fixed_point_tolerance": surface_fixed_point_tolerance,
            "neutral_surface_fixed_point_max_iterations": (
                args.surface_fixed_point_max_iterations),
            "final_neutral_surface_fixed_point_iterations": final_step.diagnostics[
                "neutral_surface_fixed_point_iterations"],
            "final_neutral_surface_fixed_point_residual": final_step.diagnostics[
                "neutral_surface_fixed_point_residual"],
            "ballistic_transport": args.ballistic_transport,
            "ballistic_face_quadrature_points": args.ballistic_face_quadrature_points,
            "reinitialize": args.reinitialize,
            "reinitialization_method": args.reinitialization_method,
            "seed": args.seed, "wall_time_s": wall_time_s,
        },
        "ion_reflection": {
            "enabled": bool(args.ion_reflection),
            "model": None if reflection is None else dict(reflection.provenance),
            "options": reflection_options,
            "final_step_bounces": final_step.diagnostics[
                "charged_surface_response_bounces"],
            "final_step_reimpact_events": final_step.diagnostics[
                "charged_surface_response_reimpact_events"],
            "final_step_relative_charge_error": final_step.diagnostics[
                "charged_surface_response_relative_charge_error"],
            "final_step_maximum_energy_error": final_step.diagnostics[
                "charged_surface_response_maximum_energy_error"],
            "final_step_tail_l1_error_bound": final_step.diagnostics[
                "charged_surface_response_tail_l1_error_bound"],
            "final_step_energy_diagnostics": _reimpact_energy_diagnostics(
                final_step, boundary, realized_domain_um, geometry,
                positive_ion_names[0]),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--control-mode", choices=("ion_energy", "ion_flux"), required=True)
    parser.add_argument(
        "--control-value", type=float, required=True,
        help="self-bias V for ion_energy; electron density m^-3 for ion_flux")
    parser.add_argument("--widths-nm", type=float, nargs="+", default=[200.0])
    parser.add_argument("--dx-um", type=float, default=0.02)
    parser.add_argument("--steps", type=int, default=40)
    parser.add_argument(
        "--transport-only", action="store_true",
        help=(
            "evaluate one exact fixed-geometry transport/response step at zero physical "
            "duration; emit cascade energy diagnostics but do not score etch depth"))
    parser.add_argument("--mask-line-width-um", type=float, default=1.0)
    parser.add_argument("--cell-length-um", type=float, default=0.06)
    parser.add_argument("--mask-thickness-um", type=float, default=1.45)
    parser.add_argument("--substrate-top-um", type=float, default=2.40)
    parser.add_argument("--domain-height-um", type=float, default=4.10)
    parser.add_argument(
        "--initial-etched-depth-um", type=float, default=0.0,
        help=(
            "diagnostic fixed-geometry trench depth; nonzero values require "
            "--transport-only and are never scored as experimental predictions"))
    parser.add_argument(
        "--chemistry-model",
        choices=("huang_2019_reduced", "reduced_si_o2", "lamagna_garozzo"),
        default="huang_2019_reduced",
        help=("select the source-faithful Huang reduced projection, the older aggregate "
              "development law, or the ViennaPS/La Magna three-coverage law"))
    parser.add_argument(
        "--lamagna-neutral-transport-mode",
        choices=("viennaps_4_6_1", "species_specific"), default="viennaps_4_6_1")
    parser.add_argument("--surface-fixed-point-tolerance", type=float, default=1e-4)
    parser.add_argument("--surface-fixed-point-max-iterations", type=int, default=20)
    parser.add_argument("--complex-probability", type=float, default=2.818e-4)
    parser.add_argument("--deposition-probability", type=float, default=1.145e-4)
    parser.add_argument("--complex-removal-reaction-order", type=int, choices=(1, 2), default=1)
    parser.add_argument(
        "--energetic-response-scale", type=float, default=1.0,
        help=("single Jeong-anchor scale multiplying all published Huang energetic p0 values; "
              "must be frozen before held-out scoring"))
    parser.add_argument("--bare-reference-yield", type=float, default=0.0)
    parser.add_argument("--complex-reference-yield", type=float, default=2.217)
    parser.add_argument("--polymer-reference-yield", type=float, default=0.50)
    parser.add_argument("--mask-reaction-probability", type=float, default=1e-3)
    parser.add_argument("--neutral-temperature-k", type=float, default=300.0)
    parser.add_argument("--electron-temperature-ev", type=float, default=3.0)
    parser.add_argument("--ion-tangential-temperature-ev", type=float, default=0.026)
    parser.add_argument("--energy-radical-support-density-m3", type=float, default=1.9e15)
    parser.add_argument("--ion-transverse-quadrature", type=int, default=3)
    parser.add_argument(
        "--ion-closure-json", type=Path,
        help=(
            "explicit Jeong positive-ion mixture/flux closure; omitted preserves "
            "the historical all-Ar development boundary"))
    parser.add_argument("--neutral-transverse-quadrature", type=int, default=3)
    parser.add_argument("--neutral-normal-quadrature", type=int, default=4)
    parser.add_argument("--source-positions", type=int, default=16)
    parser.add_argument("--form-factor-rays", type=int, default=16)
    parser.add_argument(
        "--ballistic-transport", choices=("forward", "face_gather"), default="forward")
    parser.add_argument("--ballistic-face-quadrature-points", type=int, choices=(1, 3), default=3)
    parser.add_argument(
        "--ion-reflection", action=argparse.BooleanOptionalAction, default=True,
        help="use the common certified ACL+SiO2 grazing-ion cascade")
    parser.add_argument("--reflection-fixed-dt", type=float, default=0.005)
    parser.add_argument("--reflection-max-steps", type=int, default=512)
    parser.add_argument("--reflection-emergency-max-steps", type=int, default=8192)
    parser.add_argument("--reflection-max-bounces", type=int, default=64)
    parser.add_argument("--reflection-emergency-max-bounces", type=int, default=512)
    parser.add_argument("--reflection-tail-tolerance", type=float, default=1e-8)
    parser.add_argument("--reinitialize", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--reinitialization-method", choices=("skfmm", "fsm", "cr2"), default="cr2")
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    if (args.dx_um <= 0.0 or args.steps <= 0 or args.source_positions <= 0
            or args.form_factor_rays <= 0 or any(
                not np.isfinite(value) or value < 0.0 or value > 5.0 for value in (
                    args.bare_reference_yield, args.complex_reference_yield,
                    args.polymer_reference_yield))
            or not 0.0 < args.energetic_response_scale <= 4.0
            or not 0.0 < args.surface_fixed_point_tolerance < 1.0
            or args.surface_fixed_point_max_iterations <= 0
            or args.reflection_fixed_dt <= 0.0
            or args.reflection_max_steps <= 0
            or args.reflection_emergency_max_steps < args.reflection_max_steps
            or args.reflection_max_bounces <= 0
            or args.reflection_emergency_max_bounces < args.reflection_max_bounces
            or not 0.0 < args.reflection_tail_tolerance < 1.0
            or not np.isfinite(args.initial_etched_depth_um)
            or args.initial_etched_depth_um < 0.0):
        raise ValueError("invalid Jeong 2023 campaign controls")
    if args.initial_etched_depth_um > 0.0 and not args.transport_only:
        raise ValueError(
            "a nonzero initial etched depth is diagnostic-only and requires --transport-only")
    if args.ion_reflection and args.ballistic_transport != "forward":
        raise ValueError("certified ion reflection requires forward impact-position lineage")
    if args.ion_closure_json is not None and not args.ion_closure_json.is_file():
        raise ValueError("ion closure JSON does not exist")
    # Parse before any geometry/transport work so malformed or unsupported closures
    # fail immediately rather than after an expensive profile allocation.
    ion_closure = _load_ion_closure(args.ion_closure_json)
    if args.ion_reflection and len(ion_closure.species_mass_amu) != 1:
        raise ValueError(
            "multispecies ion closures currently require --no-ion-reflection")

    depths = load_jeong_2023_etch_depths(DATA / "digitized_figure7_depths.csv")
    radicals = load_jeong_2023_radical_densities(DATA / "digitized_figure6_radicals.csv")
    provenance = _campaign_provenance(args)
    runs = []
    for width in args.widths_nm:
        target = _select_control(depths, args.control_mode, width, args.control_value)
        run = run_target(target, radicals, args)
        runs.append(run)
        print(
            f"completed {args.control_mode}={args.control_value:g}, width={width:g} nm: "
            f"prediction={run['prediction']['etch_depth_nm']:.6g} nm, "
            f"target={target.etch_depth_nm:.6g} nm, wall={run['numerics']['wall_time_s']:.2f} s",
            flush=True)
        if args.output is not None:
            _write_json_atomic(
                args.output,
                _campaign_payload(args, runs, provenance, len(args.widths_nm)))
    output = _campaign_payload(args, runs, provenance, len(args.widths_nm))
    if args.output is not None:
        _write_json_atomic(args.output, output)
    if not args.quiet:
        print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
