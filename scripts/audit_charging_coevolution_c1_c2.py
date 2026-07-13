#!/usr/bin/env python3
"""Reproduce the bounded C1 remap and C2 grazing-reflection certification figures."""
from __future__ import annotations

import csv
from hashlib import sha256
import json
import platform
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from petch.boundary_transport_3d import trace_charged_surface_events_field_3d
from petch.charged_surface_response_3d import (
    ChargedSurfaceContext3D,
    GrazingSpecularIonReflection3D,
)
from petch.surface_charge_remap_3d import remap_surface_charge_3d
from petch.surface_kinetics import FaceResolvedEnergeticFlux


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "charging_coevolution_c1_c2"


def plane_mesh(cells, z=0.0):
    coordinate = np.linspace(0.0, 1.0, int(cells) + 1)
    x, y = np.meshgrid(coordinate, coordinate, indexing="ij")
    vertices = np.column_stack((x.ravel(), y.ravel(), np.full(x.size, float(z))))
    faces = []
    stride = int(cells) + 1
    for i in range(int(cells)):
        for j in range(int(cells)):
            lower = i * stride + j
            faces.extend(((lower, lower + stride, lower + stride + 1),
                          (lower, lower + stride + 1, lower + 1)))
    return vertices, np.asarray(faces, dtype=int)


def centroids(vertices, faces):
    return np.asarray(vertices)[np.asarray(faces)].mean(axis=1)


def sigma(point):
    point = np.asarray(point)
    return 2.0e-6 + 0.4e-6 * np.sin(2.0 * np.pi * point[:, 0]) * np.cos(
        2.0 * np.pi * point[:, 1])


def remap_refinement():
    rows = []
    previous = None
    for cells in (4, 8, 16, 32):
        old_vertices, old_faces = plane_mesh(cells, z=0.0)
        new_vertices, new_faces = plane_mesh(2 * cells, z=0.1)
        old_sigma = sigma(centroids(old_vertices, old_faces))
        result = remap_surface_charge_3d(
            old_vertices, old_faces, old_sigma, np.ones(len(old_faces), dtype=int),
            np.full(len(old_faces), 0.1), new_vertices, new_faces,
            np.ones(len(new_faces), dtype=int), mesh_length_unit_m=1.0,
            neighbor_count=4, maximum_distance=0.3)
        exact = sigma(centroids(new_vertices, new_faces))
        error = float(np.sqrt(np.mean((result.sigma_c_per_m2 - exact) ** 2)))
        order = None if previous is None else float(np.log2(previous / error))
        rows.append(dict(
            cells=cells, old_faces=len(old_faces), new_faces=len(new_faces),
            rms_sigma_error_c_m2=error, observed_order=order,
            relative_charge_balance_error=result.relative_charge_balance_error))
        previous = error
    return rows


def reflection_model():
    return GrazingSpecularIonReflection3D.literature_bounded_sensitivity(
        "Si", ion_species_name="Ar+")


def corner_reimpact(model):
    vertices = np.array([
        [0.25, 0.0, 0.25], [0.25, 1.0, 0.25], [0.25, 0.0, 1.0],
        [0.25, 0.0, 0.25], [1.0, 0.0, 0.25], [0.25, 1.0, 0.25]])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    areas = np.full(2, 0.375)
    normals = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    context = ChargedSurfaceContext3D(
        areas * 1e-12, normals, np.array(["Si", "Si"]))
    count = 128
    height = np.linspace(0.35, 0.90, count)
    cosine = 0.1
    direction = np.tile([-cosine, 0.0, -np.sqrt(1.0 - cosine ** 2)], (count, 1))
    incident_rate = np.full(count, 1.0e6)
    incident = FaceResolvedEnergeticFlux(
        "Ar+", 2, np.zeros(count, dtype=int),
        incident_rate / context.face_area_m2[0], np.full(count, 100.0),
        np.full(count, cosine),
        event_position=np.column_stack((np.full(count, 0.25), np.full(count, 0.02), height)),
        event_incident_direction=direction)
    transfer = model.evaluate((incident,), {"Ar+": 1}, context)
    flight, = trace_charged_surface_events_field_3d(
        transfer.outgoing, vertices, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6,
        launch_offset=1e-5, fixed_dt=0.005, max_steps=1000, device="cpu")
    floor = flight.hit_face == 1
    distance = flight.incident.event_position[:, 0] - 0.25
    weights = (
        flight.incident.event_flux_m2_s
        * context.face_area_m2[flight.incident.event_face])
    bins = np.linspace(0.0, 0.08, 17)
    histogram, edges = np.histogram(distance, bins=bins, weights=weights)
    incident_total_rate = float(np.sum(incident_rate))
    reflected_total_rate = float(np.sum(transfer.outgoing[0].event_rate_s))
    deposited_particle_rate = incident_total_rate - reflected_total_rate
    return dict(
        incident_event_count=count,
        reflected_probability=float(model.reflection_probability([cosine])[0]),
        landed_floor_event_count=int(np.count_nonzero(floor)),
        escaped_event_count=int(np.count_nonzero(flight.termination == 2)),
        minimum_corner_distance=float(np.min(distance)),
        maximum_corner_distance=float(np.max(distance)),
        incident_particle_rate_s=incident_total_rate,
        reflected_particle_rate_s=reflected_total_rate,
        deposited_particle_rate_s=deposited_particle_rate,
        relative_particle_balance_error=float(abs(
            incident_total_rate - reflected_total_rate - deposited_particle_rate)
            / incident_total_rate),
        relative_charge_balance_error=transfer.relative_charge_balance_error,
        relative_kinetic_energy_balance_error=(
            transfer.relative_kinetic_energy_balance_error),
        reflected_flight_relative_particle_balance_error=(
            flight.relative_particle_balance_error),
        bin_left=edges[:-1].tolist(), bin_right=edges[1:].tolist(),
        reflected_rate_by_bin_s=histogram.tolist())


def make_figure(remap, corner, model):
    plt.style.use("default")
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.4))
    cells = np.array([row["cells"] for row in remap])
    error = np.array([row["rms_sigma_error_c_m2"] for row in remap])
    axes[0].loglog(1.0 / cells, error, "o-", color="#285f9e", label="measured")
    reference = error[0] * ((1.0 / cells) / (1.0 / cells[0]))
    axes[0].loglog(1.0 / cells, reference, "--", color="#666666", label="first order")
    axes[0].set_xlabel("surface spacing h")
    axes[0].set_ylabel("RMS sheet-charge error [C m$^{-2}$]")
    axes[0].set_title("C1: remap refinement")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)

    angle = np.linspace(0.0, 89.0, 300)
    cosine = np.cos(np.deg2rad(angle))
    central = model.reflection_probability(cosine)
    low = 0.80 * (1.0 - cosine ** 2.0)
    high = np.minimum(1.0, 1.0 * (1.0 - cosine ** 8.0))
    axes[1].fill_between(angle, np.minimum(low, high), np.maximum(low, high),
                         color="#d0a23c", alpha=0.25, label="declared bounds")
    axes[1].plot(angle, central, color="#a26700", label="v1 central")
    axes[1].axvline(75.0, color="#777777", linewidth=1, linestyle=":")
    axes[1].set(xlabel="incidence angle from normal [deg]",
                ylabel="reflection probability", ylim=(0.0, 1.03),
                title="C2: bounded grazing law")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, loc="upper left")

    left = np.asarray(corner["bin_left"])
    right = np.asarray(corner["bin_right"])
    rate = np.asarray(corner["reflected_rate_by_bin_s"])
    axes[2].bar(0.5 * (left + right), rate, width=0.9 * (right - left),
                color="#3f8f65", label="reflection on")
    axes[2].plot([left[0], right[-1]], [0.0, 0.0], color="#8a3b3b", linewidth=2,
                 label="reflection off")
    axes[2].set(xlabel="floor distance from wall corner [mesh units]",
                ylabel="reflected ion rate [s$^{-1}$]",
                title="C2: corner-focused floor impacts")
    axes[2].grid(True, axis="y", alpha=0.25)
    axes[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT / "c1_c2_certification.png", dpi=180)
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    remap = remap_refinement()
    model = reflection_model()
    corner = corner_reimpact(model)
    with (OUTPUT / "c1_remap_refinement.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(remap[0]))
        writer.writeheader(); writer.writerows(remap)
    config = dict(
        c1_cells=[4, 8, 16, 32], c1_new_mesh_refinement=2,
        c1_neighbor_count=4, c1_translation_mesh_units=0.1,
        c2_geometry="interior straight wall/floor corner",
        c2_incident_energy_eV=100.0, c2_incidence_cosine=0.1,
        c2_event_count=128, c2_fixed_dt=0.005, c2_launch_offset=1e-5,
        c2_max_steps=1000, c2_field="zero potential", device="cpu",
        sampling_mode="deterministic manufactured events", seed=None,
        visibility="exact hard triangle intersections")
    config_hash = sha256(
        json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    source_paths = [
        ROOT / "src/petch/surface_charge_remap_3d.py",
        ROOT / "src/petch/charged_surface_response_3d.py",
        ROOT / "src/petch/charged_surface_cascade_3d.py",
        Path(__file__).resolve()]
    source_checksums = {
        path.relative_to(ROOT).as_posix(): sha256(path.read_bytes()).hexdigest()
        for path in source_paths}
    summary = dict(
        run_manifest=dict(
            config_hash=config_hash,
            config=config,
            engine_git_revision=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            source_sha256=source_checksums,
            hardware=platform.platform(),
            python=platform.python_version()),
        closure=("charge rides advancing material; etched charged layer is removed; "
                 "newly exposed surface is uncharged"),
        c1=remap,
        c2=dict(
            parameter_values={
                "grazing_reflection_probability": model.grazing_reflection_probability,
                "angular_exponent": model.angular_exponent,
                "energy_retention_fraction": model.energy_retention_fraction},
            parameter_bounds=dict(model.parameter_bounds),
            provenance={name: value.source for name, value in model.parameter_evidence.items()},
            corner_reimpact=corner))
    (OUTPUT / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figure(remap, corner, model)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
