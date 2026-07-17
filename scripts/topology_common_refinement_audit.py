#!/usr/bin/env python3
"""Bounded close/seal/reopen audit for the conservative moving-surface engine.

The fixture is a periodic, extruded keyhole with a narrow neck above a wider
cell-resolved chamber.  A declared conformal normal-growth law closes the neck;
an externally visible directional strip law reopens it.  This is a numerical
conformance audit, not a chemistry validation or parameter calibration.

Every grid uses the same physical geometry and a timestep proportional to dx.
Runs have hard physical-time ceilings and write a compact JSON receipt plus a
PNG diagnostic.  The audit never reads experimental or held-out data.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter
from types import SimpleNamespace

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/petch-matplotlib")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from petch.boundary_state import (  # noqa: E402
    PlasmaBoundaryState,
    SpeciesBoundaryState,
)
from petch.feature_step_3d import (  # noqa: E402
    FeatureGeometry3D,
    _periodic_physical_volume_topology_signature,
    advance_feature_step_3d,
)
from petch.surface_exchange import SurfaceMaterialExchange  # noqa: E402
from petch.surface_kinetics import (  # noqa: E402
    MechanismValidity,
    SiO2SurfaceState,
)
from petch.threed import reinit_narrow  # noqa: E402


PHYSICAL_GEOMETRY = {
    "cell_width_um": 0.60,
    "cell_length_um": 0.10,
    "domain_height_um": 1.00,
    "center_um": 0.30,
    "floor_um": 0.15,
    "shoulder_um": 0.55,
    "top_um": 0.75,
    "chamber_half_width_um": 0.20,
    "neck_half_width_um": 0.075,
}


def _jsonable(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return str(value)


def _write_json(path, payload):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(target)


def _source_sha256():
    return sha256(Path(__file__).read_bytes()).hexdigest()


def _git_revision():
    return subprocess.run(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def build_keyhole_geometry(dx_um):
    """Build the same physical keyhole at one endpoint-compatible spacing."""
    dx = float(dx_um)
    extents = np.asarray((0.60, 0.10, 1.00), dtype=float)
    intervals = np.rint(extents / dx).astype(int)
    if (np.any(intervals < 2)
            or not np.allclose(intervals * dx, extents, rtol=0.0, atol=1e-13)):
        raise ValueError("dx must divide every physical extent with at least two intervals")
    shape = tuple(int(value) + 1 for value in intervals)
    x, y, z = (np.arange(size, dtype=float) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    floor = PHYSICAL_GEOMETRY["floor_um"]
    top = PHYSICAL_GEOMETRY["top_um"]
    center = PHYSICAL_GEOMETRY["center_um"]
    wide_wall = np.minimum.reduce((
        Z - floor,
        top - Z,
        np.abs(X - center) - PHYSICAL_GEOMETRY["chamber_half_width_um"],
    ))
    narrow_neck = np.minimum.reduce((
        Z - PHYSICAL_GEOMETRY["shoulder_um"],
        top - Z,
        np.abs(X - center) - PHYSICAL_GEOMETRY["neck_half_width_um"],
    ))
    phi = reinit_narrow(
        np.maximum.reduce((floor - Z, wide_wall, narrow_neck)), dx, 2.0)
    material = np.where(phi > 0.0, 1, 0)
    return FeatureGeometry3D(
        phi, material, dx, 1.0e-6, material_levelsets={1: phi})


def boundary(mode):
    if mode == "coat":
        species = (
            SpeciesBoundaryState(
                "coat", 0, 50.0, 1.0e20, [[0.0, 0.0, 1.0]], [1.0]),
            SpeciesBoundaryState(
                "probe+", 1, 40.0, 1.0e20, [[0.0, 0.0, 10.0]], [1.0]),
        )
    elif mode == "etch":
        species = (SpeciesBoundaryState(
            "etch+", 1, 40.0, 1.0e20, [[0.0, 0.0, 10.0]], [1.0]),)
    else:
        raise ValueError(mode)
    return PlasmaBoundaryState(species, reference_plane_m=1.0e-6)


class ManufacturedReversibleMotion:
    """Declared geometry law used only to certify close/seal/reopen numerics."""

    density_m3 = 1.0e28
    coat_velocity_m_s = 2.5e-8
    etch_velocity_m_s = 5.0e-8

    @staticmethod
    def initial_state(shape=()):
        return SiO2SurfaceState.bare(shape)

    @staticmethod
    def _validity():
        return MechanismValidity(
            within_declared_scope=True,
            reasons=(),
            unsupported_neutral_species=(),
            known_model_form_omissions=(
                "manufactured conformal-coat/directional-strip certification law",),
            parameter_evidence_supports_prediction=False,
            nonpredictive_parameters=("manufactured_normal_velocity",),
        )

    def advance(self, state, fluxes, duration_s):
        shape = state.polymer_units_m2.shape
        if "coat" in fluxes.neutral_flux_m2_s:
            growth = np.full(shape, self.coat_velocity_m_s)
            etch = np.zeros(shape)
        else:
            population = next(
                item for item in fluxes.energetic_fluxes if item.name == "etch+")
            incident = np.asarray(population.flux_m2_s, dtype=float)
            maximum = float(np.max(incident)) if incident.size else 0.0
            etch = (
                np.zeros(shape) if maximum == 0.0
                else self.etch_velocity_m_s * incident / maximum)
            growth = np.zeros(shape)
        removed = etch * self.density_m3 * float(duration_s)
        deposited = growth * self.density_m3 * float(duration_s)
        return SimpleNamespace(
            state=state,
            etch_velocity_m_s=etch,
            normal_growth_velocity_m_s=growth,
            material_exchange=SurfaceMaterialExchange(
                removed_units_m2={"solid_unit": removed},
                outgoing_units_m2={},
                unresolved_units_m2={"solid_unit": removed},
                deposited_units_m2={"solid_unit": deposited},
                known_limitations=("manufactured topology-conformance law",),
            ),
            product_populations=(),
            validity=self._validity(),
        )


def _center_opening_width(geometry, z_um):
    """Linearly interpolate the center-connected gas interval at one z plane."""
    field = np.asarray(geometry.phi[:, 0, :], dtype=float)
    k = int(np.clip(np.rint(float(z_um) / geometry.dx), 0, field.shape[1] - 1))
    line = field[:, k]
    center = int(np.rint(PHYSICAL_GEOMETRY["center_um"] / geometry.dx))
    if line[center] > 0.0:
        return 0.0
    left = center
    while left > 0 and line[left - 1] <= 0.0:
        left -= 1
    right = center
    while right + 1 < line.size and line[right + 1] <= 0.0:
        right += 1

    def crossing(solid_index, gas_index):
        solid_value = line[solid_index]
        gas_value = line[gas_index]
        fraction = solid_value / (solid_value - gas_value)
        return (solid_index + fraction * (gas_index - solid_index)) * geometry.dx

    left_x = 0.0 if left == 0 else crossing(left - 1, left)
    right_x = ((line.size - 1) * geometry.dx if right == line.size - 1
               else crossing(right + 1, right))
    return float(max(right_x - left_x, 0.0))


def _record(step, phase, phase_time_s, physical_time_s):
    receipt = step.state_remap_diagnostics.get("geometry_receipt", {})
    material = step.state_remap_diagnostics["materials"][1]
    event = step.diagnostics["topology_event"]
    return {
        "phase": phase,
        "phase_time_s": float(phase_time_s),
        "physical_time_s": float(physical_time_s),
        "topology": _periodic_physical_volume_topology_signature(
            step.geometry, (1,)),
        "topology_event": None if event is None else event["kind"],
        "neck_width_um": _center_opening_width(step.geometry, 0.65),
        "active_face_count": int(step.next_active_face_area.size),
        "filled_zero_volume_gas_nodes": int(
            step.diagnostics["filled_unresolved_gas_cavity_cells"]),
        "maximum_relative_conservation_residual": float(
            material["max_relative_conservation_residual"]),
        "old_matched_area_fraction": receipt.get("old_matched_area_fraction"),
        "new_matched_area_fraction": receipt.get("new_matched_area_fraction"),
        "capacity_projection_area_reduction": receipt.get(
            "capacity_projection_area_reduction"),
        "raw_overlap_area": receipt.get("raw_overlap_area"),
        "capacity_projection_iterations": receipt.get(
            "capacity_projection_iterations"),
    }


def run_level(dx_um, *, maximum_coat_time_s, maximum_etch_time_s):
    geometry = build_keyhole_geometry(dx_um)
    mechanism = ManufacturedReversibleMotion()
    state = None
    fingerprint = None
    physical_time = 0.0
    phase_time = 0.0
    dt = 5.0 * float(dx_um)
    history = []
    maximum_conservation_residual = 0.0
    minimum_matched_area_fraction = 1.0
    maximum_projection_reduction_fraction = 0.0
    started = perf_counter()

    def take_step(mode, duration_s):
        nonlocal geometry, state, fingerprint, physical_time
        plasma = boundary(mode)
        role = {
            species.name: (
                "neutral_reactant" if species.charge_number == 0
                else "energetic_bombardment")
            for species in plasma.species
        }
        result = advance_feature_step_3d(
            geometry,
            plasma,
            role,
            mechanism,
            etchable_material_ids=(1,),
            duration_s=float(duration_s),
            source_bounds=(-0.005, 0.605, -0.005, 0.105),
            source_z=1.0,
            surface_state=state,
            surface_state_mesh_fingerprint=fingerprint,
            n_position=4,
            seed=29,
            cfl_number=0.25,
            reinitialize=True,
            reinitialization_method="cr2",
            profile_periodic_lateral=True,
            transport_device="cpu",
            ballistic_transport="face_gather",
            ballistic_face_quadrature_points=1,
            surface_state_remap_backend="common_refinement",
            topology_change_policy="continue_gas_cavity",
        )
        geometry = result.geometry
        state = result.next_surface_state
        fingerprint = result.next_surface_state_mesh_fingerprint
        physical_time += float(duration_s)
        return result

    closure_time = None
    coat_steps = int(np.ceil(float(maximum_coat_time_s) / dt))
    for _ in range(coat_steps):
        step = take_step("coat", dt)
        phase_time += dt
        history.append(_record(step, "coat", phase_time, physical_time))
        if step.diagnostics["topology_event"] is not None:
            closure_time = phase_time
            break

    sealed_flux_is_zero = None
    reopening_time = None
    if closure_time is not None:
        sealed = take_step("coat", 0.4 * dt)
        centroid = sealed.active_face_centroid
        cavity = (
            (np.abs(centroid[:, 0] - PHYSICAL_GEOMETRY["center_um"]) < 0.12)
            & (centroid[:, 2] < 0.68))
        probe = next(
            item for item in sealed.transport.surface_fluxes.energetic_fluxes
            if item.name == "probe+")
        sealed_flux_is_zero = bool(
            np.any(cavity)
            and np.all(probe.flux_m2_s[cavity] == 0.0)
            and np.all(sealed.transport.surface_fluxes.neutral_flux_m2_s[
                "coat"][cavity] == 0.0))
        history.append(_record(
            sealed, "sealed_hold", 0.4 * dt, physical_time))

        phase_time = 0.0
        etch_steps = int(np.ceil(float(maximum_etch_time_s) / dt))
        for _ in range(etch_steps):
            step = take_step("etch", dt)
            phase_time += dt
            history.append(_record(step, "etch", phase_time, physical_time))
            if step.diagnostics["topology_event"] is not None:
                reopening_time = phase_time
                break

    for item in history:
        maximum_conservation_residual = max(
            maximum_conservation_residual,
            float(item["maximum_relative_conservation_residual"]))
        matched = [item["old_matched_area_fraction"], item["new_matched_area_fraction"]]
        matched = [float(value) for value in matched if value is not None]
        if matched:
            minimum_matched_area_fraction = min(
                minimum_matched_area_fraction, *matched)
        reduction = item["capacity_projection_area_reduction"]
        raw_overlap = item["raw_overlap_area"]
        if reduction is not None and raw_overlap is not None:
            maximum_projection_reduction_fraction = max(
                maximum_projection_reduction_fraction,
                float(reduction) / max(
                    float(raw_overlap), np.finfo(float).tiny))

    return {
        "dx_um": float(dx_um),
        "dt_s": dt,
        "closure_time_s": closure_time,
        "reopening_time_s": reopening_time,
        "sealed_external_flux_exactly_zero": sealed_flux_is_zero,
        "maximum_relative_conservation_residual": maximum_conservation_residual,
        "minimum_matched_area_fraction": minimum_matched_area_fraction,
        "maximum_capacity_projection_reduction_fraction": (
            maximum_projection_reduction_fraction),
        "wall_time_s": perf_counter() - started,
        "history": history,
        "passed": bool(
            closure_time is not None
            and reopening_time is not None
            and sealed_flux_is_zero
            and maximum_conservation_residual <= 5e-13),
    }


def _plot(audit, output):
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    for level in audit["levels"]:
        for phase, linestyle in (("coat", "-"), ("etch", "--")):
            selected = [item for item in level["history"] if item["phase"] == phase]
            if selected:
                axes[0].plot(
                    [item["phase_time_s"] for item in selected],
                    [item["neck_width_um"] for item in selected],
                    linestyle=linestyle,
                    marker="o",
                    markersize=2.5,
                    label=f"dx={level['dx_um'] * 1000:g} nm {phase}",
                )
    axes[0].set_xlabel("phase time (s)")
    axes[0].set_ylabel("neck opening (µm)")
    axes[0].set_title("Resolved keyhole closure and reopening")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=7, ncol=2)

    dx_nm = np.asarray([item["dx_um"] * 1000.0 for item in audit["levels"]])
    closure = np.asarray([item["closure_time_s"] for item in audit["levels"]])
    reopening = np.asarray([item["reopening_time_s"] for item in audit["levels"]])
    axes[1].plot(dx_nm, closure, "o-", label="enclosure")
    axes[1].plot(dx_nm, reopening, "s-", label="reopening after sealed hold")
    axes[1].invert_xaxis()
    axes[1].set_xlabel("grid spacing (nm; finer →)")
    axes[1].set_ylabel("event time (s)")
    axes[1].set_title("Event-time refinement")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def _refinement_summary(levels):
    """Score event-time agreement against one coarser cell-crossing time.

    A topology event inferred from a nodal sign field cannot be localized more
    tightly than a grid cell without a separate subcell event reconstruction.
    The declared normal velocities turn that spatial resolution into a physical
    time resolution; this avoids an arbitrary percentage threshold.
    """
    ordered = sorted(levels, key=lambda item: float(item["dx_um"]), reverse=True)
    if len(ordered) < 2:
        return {
            "passed": False,
            "reason": "at least two grid levels are required",
            "adjacent_pairs": [],
        }
    coat_speed_um_s = ManufacturedReversibleMotion.coat_velocity_m_s / 1.0e-6
    etch_speed_um_s = ManufacturedReversibleMotion.etch_velocity_m_s / 1.0e-6
    pairs = []
    for coarser, finer in zip(ordered[:-1], ordered[1:]):
        closure_difference = abs(
            float(coarser["closure_time_s"]) - float(finer["closure_time_s"]))
        reopening_difference = abs(
            float(coarser["reopening_time_s"]) - float(finer["reopening_time_s"]))
        closure_resolution = float(coarser["dx_um"]) / coat_speed_um_s
        reopening_resolution = float(coarser["dx_um"]) / etch_speed_um_s
        pair_passed = bool(
            closure_difference <= closure_resolution
            and reopening_difference <= reopening_resolution)
        pairs.append({
            "coarser_dx_um": float(coarser["dx_um"]),
            "finer_dx_um": float(finer["dx_um"]),
            "closure_time_difference_s": closure_difference,
            "closure_one_coarse_cell_crossing_time_s": closure_resolution,
            "reopening_time_difference_s": reopening_difference,
            "reopening_one_coarse_cell_crossing_time_s": reopening_resolution,
            "passed": pair_passed,
        })
    authoritative = pairs[-1]
    return {
        "criterion": (
            "adjacent event times agree within one coarser-grid cell crossing "
            "at the declared maximum phase normal velocity"),
        "coat_velocity_um_s": coat_speed_um_s,
        "etch_velocity_um_s": etch_speed_um_s,
        "adjacent_pairs": pairs,
        "all_adjacent_pairs_passed": bool(all(item["passed"] for item in pairs)),
        "authoritative_pair": authoritative,
        "coarsest_level_is_diagnostic_only": bool(
            len(pairs) > 1 and not pairs[0]["passed"]),
        "passed": bool(authoritative["passed"]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--levels-um", type=float, nargs="+", default=(0.05, 0.025, 0.0125))
    parser.add_argument("--maximum-coat-time-s", type=float, default=8.0)
    parser.add_argument("--maximum-etch-time-s", type=float, default=12.0)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "results" / "topology_common_refinement_audit")
    parser.add_argument(
        "--summarize-existing", action="store_true",
        help="recompute derived refinement/plot fields without rerunning the engine")
    args = parser.parse_args()
    if args.maximum_coat_time_s <= 0.0 or args.maximum_etch_time_s <= 0.0:
        raise ValueError("physical-time ceilings must be positive")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.summarize_existing:
        audit_path = output / "audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit["refinement"] = _refinement_summary(audit["levels"])
        audit["postprocessor_script_sha256"] = _source_sha256()
        audit["passed"] = bool(
            all(item["passed"] for item in audit["levels"])
            and audit["refinement"]["passed"])
        _write_json(audit_path, audit)
        _plot(audit, output / "topology_refinement.png")
        print(json.dumps({
            "passed": audit["passed"],
            "refinement": audit["refinement"],
        }, indent=2))
        return 0 if audit["passed"] else 2
    audit = {
        "schema": "petch.topology_common_refinement_audit.v1",
        "scope": "manufactured_numerical_conformance_not_chemistry_validation",
        "git_revision": _git_revision(),
        "script_sha256": _source_sha256(),
        "physical_geometry": PHYSICAL_GEOMETRY,
        "operator": {
            "surface_state_remap_backend": "common_refinement",
            "topology_change_policy": "continue_gas_cavity",
            "profile_periodic_lateral": True,
            "reinitialization_method": "cr2",
            "ballistic_transport": "face_gather",
            "ballistic_face_quadrature_points": 1,
            "n_position": 4,
            "seed": 29,
            "timestep_rule": "dt_s = 5 * dx_um",
        },
        "ceilings": {
            "maximum_coat_time_s": args.maximum_coat_time_s,
            "maximum_etch_time_s": args.maximum_etch_time_s,
        },
        "levels": [],
    }
    for dx in args.levels_um:
        level = run_level(
            dx,
            maximum_coat_time_s=args.maximum_coat_time_s,
            maximum_etch_time_s=args.maximum_etch_time_s,
        )
        audit["levels"].append(level)
        _write_json(output / "audit.json", audit)
    audit["refinement"] = _refinement_summary(audit["levels"])
    audit["passed"] = bool(
        all(item["passed"] for item in audit["levels"])
        and audit["refinement"]["passed"])
    audit["total_wall_time_s"] = float(sum(
        item["wall_time_s"] for item in audit["levels"]))
    _write_json(output / "audit.json", audit)
    _plot(audit, output / "topology_refinement.png")
    print(json.dumps({
        "passed": audit["passed"],
        "total_wall_time_s": audit["total_wall_time_s"],
        "levels": [{
            "dx_um": item["dx_um"],
            "closure_time_s": item["closure_time_s"],
            "reopening_time_s": item["reopening_time_s"],
            "wall_time_s": item["wall_time_s"],
        } for item in audit["levels"]],
    }, indent=2))
    return 0 if audit["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
