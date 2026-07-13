"""Replay one charging state and expose the incident energy/angle evidence by surface region.

This is deliberately not a surface-response model.  It evaluates the unchanged absorbing,
hard-visibility kinetic operator and reports the event measure that a future material response would
consume.  Energy-threshold fractions are diagnostics only; they are not secondary-electron yields.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from petch.charging_coupled_3d import advance_dielectric_charging_3d

import charging_task1_physical_time_3d as task1


DEFAULT_STATE = (
    ROOT / "results/charging_task1_3d_refined_15us_level13_audit/final_states.npz")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _weighted_quantile(values, weights, probabilities):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    probabilities = np.asarray(probabilities, dtype=float)
    if (values.shape != weights.shape or values.ndim != 1 or values.size == 0
            or np.any(~np.isfinite(values)) or np.any(~np.isfinite(weights))
            or np.any(weights < 0.0) or not np.sum(weights) > 0.0
            or np.any((probabilities < 0.0) | (probabilities > 1.0))):
        raise ValueError("weighted quantiles require a nonempty finite positive event measure")
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    targets = probabilities * cumulative[-1]
    indices = np.searchsorted(cumulative, targets, side="left")
    return sorted_values[np.minimum(indices, sorted_values.size - 1)]


def _surface_regions(centroids):
    """Classify the fixed Task-1 trench without using noisy triangle-normal thresholds."""
    centroids = np.asarray(centroids, dtype=float)
    distance = np.column_stack((
        np.abs(centroids[:, 2] - 0.5),
        np.abs(centroids[:, 0] - 0.25),
        np.abs(centroids[:, 0] - 0.75),
        np.abs(centroids[:, 2] - 1.5),
    ))
    base = np.asarray(("floor", "left_wall", "right_wall", "top"), dtype=object)[
        np.argmin(distance, axis=1)]
    region = base.copy()
    wall = np.isin(base, ("left_wall", "right_wall"))
    region[wall & (centroids[:, 2] >= 0.875)] = "upper_wall"
    region[wall & (centroids[:, 2] < 0.875)] = "lower_wall"
    return region.astype("U16")


def _sobolewski_2021_ar_kinetic_yield(energy_eV):
    """SiO2 Ar+ kinetic-emission fit, Sobolewski PSST 30 025004 (2021), Eq. 8."""
    energy = np.asarray(energy_eV, dtype=float)
    return 0.030 * energy ** 2 / (200.0 + energy) ** 1.5


def _event_rows(populations, regions, areas, area_unit_m2, thresholds):
    rows = []
    for population in populations:
        event_region = regions[population.event_face]
        # event_flux is a face flux density contribution.  Multiplication by physical face area
        # gives the particle-rate measure required for region aggregation.
        event_rate = (
            population.event_flux_m2_s
            * areas[population.event_face]
            * float(area_unit_m2))
        for region in ("top", "upper_wall", "lower_wall", "floor", "all"):
            selected = (np.ones(event_rate.size, dtype=bool) if region == "all"
                        else event_region == region)
            if not np.any(selected) or not np.sum(event_rate[selected]) > 0.0:
                continue
            energy = population.event_energy_eV[selected]
            cosine = population.event_cosine_incidence[selected]
            weight = event_rate[selected]
            quantile = _weighted_quantile(energy, weight, (0.1, 0.5, 0.9, 0.99))
            row = dict(
                species=population.name,
                region=region,
                incident_rate_s=float(np.sum(weight)),
                mean_energy_eV=float(np.average(energy, weights=weight)),
                energy_p10_eV=float(quantile[0]),
                energy_p50_eV=float(quantile[1]),
                energy_p90_eV=float(quantile[2]),
                energy_p99_eV=float(quantile[3]),
                mean_cosine_incidence=float(np.average(cosine, weights=weight)),
                sobolewski2021_kinetic_see_mean_yield=None,
            )
            for threshold in thresholds:
                row[f"rate_fraction_energy_ge_{threshold:g}_eV"] = float(
                    np.sum(weight[energy >= threshold]) / np.sum(weight))
            if population.name == "Ar+":
                row["sobolewski2021_kinetic_see_mean_yield"] = float(np.average(
                    _sobolewski_2021_ar_kinetic_yield(energy), weights=weight))
            rows.append(row)
    return rows


def _write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot(path, rows, populations, regions, areas, area_unit_m2, thresholds):
    import matplotlib.pyplot as plt

    display_regions = ("top", "upper_wall", "lower_wall", "floor")
    colors = dict(top="#4c78a8", upper_wall="#e45756", lower_wall="#f2cf5b", floor="#54a24b")
    by_key = {(row["species"], row["region"]): row for row in rows}
    ion_name = next(item.name for item in populations if item.name != "electron")
    signed = []
    for region in display_regions:
        ion = by_key[(ion_name, region)]["incident_rate_s"]
        electron = by_key[("electron", region)]["incident_rate_s"]
        signed.append((ion - electron) / (ion + electron))

    figure, axes = plt.subplots(1, 3, figsize=(14.4, 4.4), constrained_layout=True)
    axes[0].bar(display_regions, signed, color=[colors[item] for item in display_regions])
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].axhline(0.08, color="0.35", ls="--", lw=1)
    axes[0].axhline(-0.08, color="0.35", ls="--", lw=1, label="±0.08 contract")
    axes[0].set(title="Integrated signed current", ylabel="(ion − electron) / total")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].legend(fontsize=8)

    electron = next(item for item in populations if item.name == "electron")
    event_region = regions[electron.event_face]
    event_rate = (
        electron.event_flux_m2_s * areas[electron.event_face] * float(area_unit_m2))
    for region in display_regions:
        selected = event_region == region
        order = np.argsort(electron.event_energy_eV[selected], kind="stable")
        energy = electron.event_energy_eV[selected][order]
        cumulative = np.cumsum(event_rate[selected][order])
        axes[1].step(energy, cumulative / cumulative[-1], where="post",
                     color=colors[region], label=region)
    axes[1].set_xscale("log")
    axes[1].set(title="Incident-electron energy CDF", xlabel="impact energy (eV)", ylabel="rate CDF")
    axes[1].legend(fontsize=8)

    width = 0.19
    positions = np.arange(len(thresholds), dtype=float)
    for index, region in enumerate(display_regions):
        row = by_key[("electron", region)]
        fraction = [row[f"rate_fraction_energy_ge_{value:g}_eV"] for value in thresholds]
        axes[2].bar(positions + (index - 1.5) * width, fraction, width,
                    color=colors[region], label=region)
    axes[2].set_xticks(positions, [f"≥{value:g}" for value in thresholds])
    axes[2].set_ylim(0.0, 1.0)
    axes[2].set(title="Electron rate above energy", xlabel="threshold (eV)", ylabel="fraction")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(axis="y", alpha=0.22)
    figure.suptitle("Surface-response preflight — unchanged absorbing hard-visibility operator")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--charge-key", default="refined_charge_node_c")
    parser.add_argument("--method-key", default="refined_method_hint_Ar+")
    parser.add_argument("--forward-level", type=int, default=11)
    parser.add_argument("--thresholds-eV", type=float, nargs="+", default=(10.0, 25.0, 50.0, 100.0))
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "results/charging_surface_response_preflight_3d")
    args = parser.parse_args()
    thresholds = tuple(sorted({float(item) for item in args.thresholds_eV}))
    if not thresholds or thresholds[0] < 0.0:
        parser.error("energy thresholds must be nonnegative")

    geometry, poisson = task1._geometry_and_poisson(0.125)
    verts, faces, centroids, areas = task1.extract_mesh_3d(geometry.phi, geometry.dx)
    normals = task1._surface_gas_normals(verts, faces, centroids, geometry)
    with np.load(args.state) as archive:
        charge = np.asarray(archive[args.charge_key], dtype=float).copy()
        method = np.asarray(archive[args.method_key]).astype("U7")
    if charge.shape != poisson.shape or method.shape != (len(faces),):
        parser.error("checkpoint charge or frozen method map does not match the Task-1 grid")

    source_z = 2.0
    boundary = task1._boundary(source_z * geometry.mesh_length_unit_m)
    proposals = {
        "Ar+": task1._ion_proposal(boundary, 9, 79),
        "electron": task1._electron_proposal(boundary, 9, 83),
    }
    bidirectional = dict(
        forward_log2_samples=args.forward_level,
        adjoint_log2_samples=9,
        n_replicates=4,
        max_forward_log2_samples=args.forward_level,
        max_adjoint_log2_samples=args.forward_level,
        max_face_quadrature_points=7,
        element_absolute_tolerance=0.02,
        element_relative_tolerance=0.1,
        face_quadrature_points=3,
        method_hint={"Ar+": method},
        require_certification=False)
    result = advance_dielectric_charging_3d(
        poisson, charge, boundary, verts, faces, areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        duration_s=1.0e-15, mesh_length_unit_m=geometry.mesh_length_unit_m,
        n_position=256, seed=79, trajectory_fixed_dt=0.005,
        trajectory_max_steps=50000, phase_space_log2_samples=args.forward_level,
        periodic_lateral=True,
        transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1.0e-4,
        adjoint_proposals=proposals,
        adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
        bidirectional_options=bidirectional, transport_device="cpu")
    populations = result.transport.surface_fluxes.energetic_fluxes
    regions = _surface_regions(centroids)
    area_unit_m2 = geometry.mesh_length_unit_m ** 2
    rows = _event_rows(populations, regions, areas, area_unit_m2, thresholds)

    config = dict(
        schema="petch.charging.surface_response_preflight.v1",
        state=args.state.name,
        state_sha256=_sha256(args.state),
        charge_key=args.charge_key,
        method_key=args.method_key,
        forward_level=args.forward_level,
        thresholds_eV=list(thresholds),
        diagnostic_yield_model=(
            "Ar+ kinetic SEE only: Sobolewski, PSST 30 025004 (2021), Eq. 8; "
            "not applied to charge or transport"),
        exact_operator="absorbing first-hit, hard visibility; no surface response applied")
    config_hash = hashlib.sha256(
        json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    summary = dict(
        config_hash=config_hash,
        config=config,
        rows=rows,
        limitations=(
            "energy-threshold fractions are not emission or reflection yields",
            "the reported Ar+ kinetic-SEE yield is a diagnostic fit, not an applied response",
            "material state, roughness, fluorocarbon coverage, and surface conditioning are unresolved",
            "the replay uses the archived frozen ion estimator map and unchanged absorbing current law",
        ))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _write_csv(args.output_dir / "incident_event_summary.csv", rows)
    _plot(
        args.output_dir / "surface_response_preflight.png", rows, populations, regions, areas,
        area_unit_m2, thresholds)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
