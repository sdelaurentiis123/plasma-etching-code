"""Angular-convergence harness (roadmap P0 / IADF research S0) — the ruler.

Measures how much of the deep-feature flux split is decided by *angular
representation* rather than physics, on ideal static geometry with the
production deterministic face gather.  No chemistry, no surface evolution,
no box: every number here is pure transport.

Two experiments:

EXP A  quadrature ablation of the current collisionless boundary — Gauss-Hermite
       transverse order (the virtual-sheath path) and axisymmetric azimuthal
       order / polar bin width (the digitized-IEAD path) at AR 30/100/200.

EXP B  the mouth question at Krueger AR ~9: wall-flux profiles under the narrow
       collisionless beam, the measured two-component (core+tail) beam, and the
       digitized Krueger IEAD actually consumed by the feature validation.

Self-consistency gate: a single collimated direction into a straight trench must
reproduce the analytic view factor of the unshadowed floor, 1 - AR*tan(theta).
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.boundary_transport_3d import gather_boundary_state_ballistic_3d
from petch.reactor_boundary import load_krueger_2024_digitized_iead

_DEFAULT_DATA = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "data", "experimental", "krueger_2024")

KRUEGER_MEAN_ION_ENERGY_EV = 3465.0
# Kim et al., JJAP 64, 05SP15 (2025), 0.1-deg-resolved MCP measurement:
# thermal core 0.044 eV / collisional tail 0.57 eV, projected to 1-D sigma at
# the Krueger mean ion energy (sigma = sqrt(T_perp / 2E)); tail fraction 0.65
# from the sheath survival law 1 - exp(-s/lambda) at Krueger conditions.
MEASURED_CORE_SIGMA_DEG = np.degrees(np.sqrt(0.044 / (2.0 * KRUEGER_MEAN_ION_ENERGY_EV)))
MEASURED_TAIL_SIGMA_DEG = np.degrees(np.sqrt(0.57 / (2.0 * KRUEGER_MEAN_ION_ENERGY_EV)))
MEASURED_TAIL_FRACTION = 0.65
# petch virtual-sheath width measured in RESEARCH_IADF_SUBDEGREE_AND_REACTOR.
PETCH_VIRTUAL_SHEATH_SIGMA_DEG = 0.1478


# ----------------------------------------------------------------- geometry


def trench_mesh(aspect, *, width=1.0, mask_pad=1.0, y_width=1.0, n_wall=60,
                n_floor=200):
    """Straight trench: floor, two sidewalls in ``n_wall`` depth bands, mask top.

    Returns (verts, faces, areas, centroids, normals, meta).  Gas normals point
    into the trench volume.  Face count is independent of aspect ratio so the
    ablation isolates angular resolution from spatial resolution.
    """
    depth = float(aspect) * float(width)
    x0, x1 = float(mask_pad), float(mask_pad) + float(width)
    lx = float(width) + 2.0 * float(mask_pad)
    ly = float(y_width)
    verts, faces, normals, tags = [], [], [], []

    def quad(p0, p1, p2, p3, normal, tag):
        base = len(verts)
        verts.extend([p0, p1, p2, p3])
        faces.append([base, base + 1, base + 2])
        faces.append([base, base + 2, base + 3])
        normals.extend([normal, normal])
        tags.extend([tag, tag])

    # The floor is banded across x so that partial shadowing is resolved
    # spatially; otherwise the floor/wall split is a step function of where the
    # centroid of a single wide face happens to land.
    floor_edges = np.linspace(x0, x1, int(n_floor) + 1)
    for lo, hi in zip(floor_edges[:-1], floor_edges[1:]):
        quad([lo, 0.0, 0.0], [hi, 0.0, 0.0], [hi, ly, 0.0], [lo, ly, 0.0],
             [0.0, 0.0, 1.0], "floor")
    edges = np.linspace(0.0, depth, int(n_wall) + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        quad([x0, 0.0, lo], [x0, 0.0, hi], [x0, ly, hi], [x0, ly, lo],
             [1.0, 0.0, 0.0], f"wall_lo:{0.5 * (lo + hi):.6f}")
        quad([x1, 0.0, lo], [x1, ly, lo], [x1, ly, hi], [x1, 0.0, hi],
             [-1.0, 0.0, 0.0], f"wall_hi:{0.5 * (lo + hi):.6f}")
    quad([0.0, 0.0, depth], [x0, 0.0, depth], [x0, ly, depth], [0.0, ly, depth],
         [0.0, 0.0, 1.0], "mask")
    quad([x1, 0.0, depth], [lx, 0.0, depth], [lx, ly, depth], [x1, ly, depth],
         [0.0, 0.0, 1.0], "mask")

    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=int)
    normals = np.asarray(normals, dtype=float)
    tri = verts[faces]
    areas = 0.5 * np.linalg.norm(np.cross(
        tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    centroids = tri.mean(axis=1)
    meta = {
        "aspect": float(aspect), "width": float(width), "depth": depth,
        "lx": lx, "ly": ly, "source_z": depth + float(width),
        "tags": np.asarray(tags), "n_wall": int(n_wall),
        "n_floor": int(n_floor),
    }
    return verts, faces, areas, centroids, normals, meta


# ----------------------------------------------------------------- beams


def _axisymmetric_species(polar_deg, weight, *, azimuthal_order,
                          energy_eV=KRUEGER_MEAN_ION_ENERGY_EV, flux=1.0e20,
                          name="ions"):
    """Lift a polar-angle measure to 3-D by uniform azimuth (production closure)."""
    polar = np.deg2rad(np.abs(np.asarray(polar_deg, dtype=float)))
    weight = np.asarray(weight, dtype=float)
    weight = weight / weight.sum()
    energy = np.broadcast_to(np.asarray(energy_eV, dtype=float), polar.shape)
    speed = np.sqrt(energy)
    order = int(azimuthal_order)
    azimuth = 2.0 * np.pi * (np.arange(order, dtype=float) + 0.5) / order
    transverse = (speed * np.sin(polar))[:, None]
    velocity = np.column_stack((
        (transverse * np.cos(azimuth)[None, :]).ravel(),
        (transverse * np.sin(azimuth)[None, :]).ravel(),
        np.repeat(speed * np.cos(polar), order),
    ))
    return SpeciesBoundaryState(
        name, 1, 39.948, float(flux), velocity_sqrt_eV=velocity,
        weight=np.repeat(weight / order, order))


def rayleigh_polar_measure(sigma_deg, *, n_polar=400, max_sigma=6.0):
    """Polar measure of an axisymmetric beam whose 1-D projected width is sigma.

    Transverse velocity components are independent Gaussians, so the polar angle
    carries the Rayleigh measure p(theta) d(theta) ~ theta exp(-theta^2/2 sigma^2).
    """
    sigma = float(sigma_deg)
    edges = np.linspace(0.0, max_sigma * sigma, int(n_polar) + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    weight = np.exp(-edges[:-1] ** 2 / (2.0 * sigma ** 2)) - np.exp(
        -edges[1:] ** 2 / (2.0 * sigma ** 2))
    return mid, weight


def gaussian_beam(sigma_deg, *, azimuthal_order=16, n_polar=400, **kw):
    mid, weight = rayleigh_polar_measure(sigma_deg, n_polar=n_polar)
    return _axisymmetric_species(mid, weight, azimuthal_order=azimuthal_order, **kw)


def two_component_beam(core_sigma_deg, tail_sigma_deg, tail_fraction, *,
                       azimuthal_order=16, n_polar=400, **kw):
    """Measured bi-Gaussian: thermal core plus sheath-collision tail."""
    span = max(core_sigma_deg, tail_sigma_deg)
    edges = np.linspace(0.0, 6.0 * span, int(n_polar) + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])

    def band(sigma):
        return (np.exp(-edges[:-1] ** 2 / (2.0 * sigma ** 2))
                - np.exp(-edges[1:] ** 2 / (2.0 * sigma ** 2)))

    weight = ((1.0 - float(tail_fraction)) * band(float(core_sigma_deg))
              + float(tail_fraction) * band(float(tail_sigma_deg)))
    return _axisymmetric_species(mid, weight, azimuthal_order=azimuthal_order, **kw)


def virtual_sheath_ion_species(order, *, azimuthal_order=None):
    """The current production virtual-sheath ion beam at Gauss-Hermite ``order``."""
    from petch.reactor_boundary import (
        PlasmaDiagnosticState, build_diagnostic_virtual_sheath_boundary)
    from petch.sheath import PeriodicSheathVoltage

    diagnostic = PlasmaDiagnosticState(
        electron_density_m3=3e16, electron_temperature_eV=3.6, ion_name="ions",
        ion_mass_amu=40.0, source="probe", ion_flux_m2_s=1e20,
        ion_flux_evidence_kind="assumed")
    waveform = PeriodicSheathVoltage.sinusoidal(
        dc_v=2500.0, amplitude_v=2300.0, frequency_hz=1e6, source="probe")
    boundary = build_diagnostic_virtual_sheath_boundary(
        diagnostic, waveform, reference_plane_m=0.0,
        collisionless_justification="probe", n_transverse_ion=int(order))
    ion = boundary.get("ions")
    return SpeciesBoundaryState(
        "ions", 1, float(ion.mass_amu), float(ion.flux_m2_s),
        velocity_sqrt_eV=np.asarray(ion.velocity_sqrt_eV, dtype=float),
        weight=np.asarray(ion.weight, dtype=float))


def krueger_iead_polar_scaled_species(scale, *, azimuthal_order=16,
                                      angle_bin_deg=0.25, energy_bin_eV=250.0,
                                      data_directory=None):
    """Digitized IEAD lifted with polar angles scaled by ``scale``.

    Diagnostic for the axisymmetric closure: the production lift uses the
    published *planar signed angle* directly as the polar angle, which makes the
    lifted beam's planar marginal narrower than the published one by exactly
    sqrt(2) for a Gaussian-transverse beam.  scale=sqrt(2) restores the
    published planar width.  This is a measurement of the closure error, not a
    proposed physical model.
    """
    iead = load_krueger_2024_digitized_iead(data_directory or _DEFAULT_DATA)
    energy, signed_angle, weight, _ = iead.development_quadrature(
        energy_bin_eV=energy_bin_eV, angle_bin_deg=angle_bin_deg)
    return _axisymmetric_species(
        np.abs(np.asarray(signed_angle, dtype=float)) * float(scale),
        np.asarray(weight, dtype=float), azimuthal_order=azimuthal_order,
        energy_eV=np.asarray(energy, dtype=float))


def krueger_iead_species(*, azimuthal_order=16, angle_bin_deg=None,
                         energy_bin_eV=None, data_directory=None):
    """The digitized Krueger Fig-4 IEAD as the feature validation consumes it."""
    iead = load_krueger_2024_digitized_iead(data_directory or _DEFAULT_DATA)
    return iead.development_species(
        1.0e20, effective_mass_amu=39.948,
        mixture_closure="harness diagnostic; aggregate singly charged",
        name="ions", energy_bin_eV=energy_bin_eV, angle_bin_deg=angle_bin_deg,
        azimuthal_closure="axisymmetric_uniform", azimuthal_order=int(azimuthal_order))


def beam_statistics(species):
    """Projected 1-D width, polar support and acceptance fractions of a beam."""
    velocity = np.asarray(species.velocity_sqrt_eV, dtype=float)
    weight = np.asarray(species.weight, dtype=float)
    weight = weight / weight.sum()
    theta_x = np.degrees(np.arctan2(velocity[:, 0], velocity[:, 2]))
    polar = np.degrees(np.arctan2(
        np.hypot(velocity[:, 0], velocity[:, 1]), velocity[:, 2]))
    stats = {
        "sigma_theta_1d_deg": float(np.sqrt(np.average(theta_x ** 2, weights=weight))),
        "max_polar_deg": float(polar.max()),
        "nodes": int(velocity.shape[0]),
    }
    for aspect in (9.0, 30.0, 100.0, 200.0):
        cone = np.degrees(np.arctan(1.0 / aspect))
        stats[f"accept_AR{int(aspect)}"] = float(weight[polar < cone].sum())
    return stats


# ----------------------------------------------------------------- transport


def gather(species, mesh, *, face_quadrature_points=1):
    verts, faces, areas, centroids, normals, meta = mesh
    boundary = PlasmaBoundaryState(
        (species,), reference_plane_m=meta["source_z"] * 1e-6)
    result = gather_boundary_state_ballistic_3d(
        boundary, {species.name: "energetic_bombardment"},
        verts, faces, areas, centroids, normals,
        source_bounds=(0.0, meta["lx"], 0.0, meta["ly"]),
        source_z=meta["source_z"], mesh_length_unit_m=1e-6,
        face_quadrature_points=int(face_quadrature_points),
        periodic_lateral=True,
        domain_size=(meta["lx"], meta["ly"], meta["source_z"]), device="cpu")
    energetic = result.surface_fluxes.energetic_fluxes
    rate = np.zeros(len(faces), dtype=float)
    for population in energetic:
        np.add.at(rate, np.asarray(population.event_face, dtype=int),
                  np.asarray(population.event_flux_m2_s, dtype=float)
                  * areas[np.asarray(population.event_face, dtype=int)])
    return rate


def split_observables(rate, mesh, species):
    verts, faces, areas, centroids, normals, meta = mesh
    tags = meta["tags"]
    mouth_rate = float(species.flux_m2_s) * meta["width"] * meta["ly"]
    is_floor = tags == "floor"
    is_wall = np.char.startswith(tags.astype(str), "wall")
    depth = np.array([
        float(tag.split(":")[1]) if tag.startswith("wall") else np.nan
        for tag in tags.astype(str)])
    edges = np.linspace(0.0, meta["depth"], meta["n_wall"] + 1)
    profile = np.zeros(meta["n_wall"], dtype=float)
    index = np.clip(np.searchsorted(edges, depth[is_wall], side="right") - 1,
                    0, meta["n_wall"] - 1)
    np.add.at(profile, index, rate[is_wall])
    return {
        "floor_fraction": float(rate[is_floor].sum() / mouth_rate),
        "wall_fraction": float(rate[is_wall].sum() / mouth_rate),
        "mask_fraction": float(rate[tags == "mask"].sum()
                               / (float(species.flux_m2_s)
                                  * (meta["lx"] - meta["width"]) * meta["ly"])),
        "wall_profile": profile / mouth_rate,
        "profile_depth": 0.5 * (edges[:-1] + edges[1:]),
        "mouth_rate": mouth_rate,
    }


# ----------------------------------------------------------------- gate


def view_factor_gate(aspect=30.0, theta_deg=1.0):
    """Collimated direction into a straight trench: floor fraction = 1 - AR*tan(theta)."""
    mesh = trench_mesh(aspect, n_wall=40)
    speed = np.sqrt(KRUEGER_MEAN_ION_ENERGY_EV)
    theta = np.deg2rad(theta_deg)
    species = SpeciesBoundaryState(
        "ions", 1, 39.948, 1.0e20,
        velocity_sqrt_eV=[[speed * np.sin(theta), 0.0, speed * np.cos(theta)]],
        weight=[1.0])
    observed = split_observables(gather(species, mesh), mesh, species)["floor_fraction"]
    expected = max(0.0, 1.0 - aspect * np.tan(theta))
    return observed, expected


# ----------------------------------------------------------------- experiments


def experiment_a(aspects=(30.0, 100.0, 200.0)):
    rows = []
    for aspect in aspects:
        mesh = trench_mesh(aspect, n_wall=60)
        for order in (3, 5, 9, 17):
            species = virtual_sheath_ion_species(order)
            stats = beam_statistics(species)
            observed = split_observables(gather(species, mesh), mesh, species)
            rows.append({
                "experiment": "A1_gauss_hermite", "aspect": aspect,
                "control": f"n_transverse={order}", **stats,
                "floor_fraction": observed["floor_fraction"],
                "wall_fraction": observed["wall_fraction"]})
            print(f"  A1 AR{aspect:.0f} GH={order:2d}  sigma={stats['sigma_theta_1d_deg']:.4f}"
                  f"  maxpolar={stats['max_polar_deg']:.3f}"
                  f"  floor={observed['floor_fraction']:.4f}"
                  f"  wall={observed['wall_fraction']:.4f}", flush=True)
        for order in (4, 8, 16, 32, 64):
            species = krueger_iead_species(azimuthal_order=order, angle_bin_deg=0.25,
                                           energy_bin_eV=250.0)
            stats = beam_statistics(species)
            observed = split_observables(gather(species, mesh), mesh, species)
            rows.append({
                "experiment": "A2_azimuth", "aspect": aspect,
                "control": f"azimuth={order}", **stats,
                "floor_fraction": observed["floor_fraction"],
                "wall_fraction": observed["wall_fraction"]})
            print(f"  A2 AR{aspect:.0f} az={order:3d}  sigma={stats['sigma_theta_1d_deg']:.4f}"
                  f"  floor={observed['floor_fraction']:.4f}"
                  f"  wall={observed['wall_fraction']:.4f}", flush=True)
        for label, kwargs in (("exact_digitized", {}),
                              ("bin_0.25deg", {"angle_bin_deg": 0.25, "energy_bin_eV": 250.0}),
                              ("bin_0.50deg", {"angle_bin_deg": 0.50, "energy_bin_eV": 250.0}),
                              ("bin_1.00deg", {"angle_bin_deg": 1.00, "energy_bin_eV": 250.0})):
            species = krueger_iead_species(azimuthal_order=16, **kwargs)
            stats = beam_statistics(species)
            observed = split_observables(gather(species, mesh), mesh, species)
            rows.append({
                "experiment": "A3_polar_bin", "aspect": aspect, "control": label,
                **stats, "floor_fraction": observed["floor_fraction"],
                "wall_fraction": observed["wall_fraction"]})
            print(f"  A3 AR{aspect:.0f} {label:16s} sigma={stats['sigma_theta_1d_deg']:.4f}"
                  f"  floor={observed['floor_fraction']:.4f}"
                  f"  wall={observed['wall_fraction']:.4f}", flush=True)
    return rows


def experiment_b(aspect=9.0, n_wall=120):
    mesh = trench_mesh(aspect, n_wall=n_wall)
    beams = {
        "narrow_collisionless": gaussian_beam(PETCH_VIRTUAL_SHEATH_SIGMA_DEG),
        "measured_two_component": two_component_beam(
            MEASURED_CORE_SIGMA_DEG, MEASURED_TAIL_SIGMA_DEG, MEASURED_TAIL_FRACTION),
        "krueger_digitized_iead": krueger_iead_species(
            azimuthal_order=16, angle_bin_deg=0.25, energy_bin_eV=250.0),
        "krueger_iead_sqrt2_corrected": krueger_iead_polar_scaled_species(np.sqrt(2.0)),
    }
    out = {}
    for label, species in beams.items():
        stats = beam_statistics(species)
        observed = split_observables(gather(species, mesh), mesh, species)
        out[label] = {**stats,
                      "floor_fraction": observed["floor_fraction"],
                      "wall_fraction": observed["wall_fraction"],
                      "wall_profile": observed["wall_profile"].tolist(),
                      "profile_depth": observed["profile_depth"].tolist()}
        print(f"  B {label:24s} sigma={stats['sigma_theta_1d_deg']:.4f}"
              f"  floor={observed['floor_fraction']:.4f}"
              f"  wall={observed['wall_fraction']:.4f}", flush=True)
    depth = np.asarray(out["narrow_collisionless"]["profile_depth"])
    mouth_band = depth > 0.75 * depth.max()          # upper quarter = mouth region
    ratios = {}
    for label in ("measured_two_component", "krueger_digitized_iead",
                  "krueger_iead_sqrt2_corrected"):
        num = np.asarray(out[label]["wall_profile"])[mouth_band].sum()
        den = np.asarray(out["narrow_collisionless"]["wall_profile"])[mouth_band].sum()
        ratios[label] = float(num / den) if den > 0 else float("inf")
    out["mouth_region_ratio_vs_narrow"] = ratios
    return out


def experiment_c(aspects=(9.0, 30.0, 100.0, 200.0)):
    """AR-resolved size of the axisymmetric-closure (sqrt2) width loss."""
    rows = []
    for aspect in aspects:
        mesh = trench_mesh(aspect, n_wall=60)
        production = krueger_iead_species(
            azimuthal_order=16, angle_bin_deg=0.25, energy_bin_eV=250.0)
        corrected = krueger_iead_polar_scaled_species(np.sqrt(2.0))
        a = split_observables(gather(production, mesh), mesh, production)
        b = split_observables(gather(corrected, mesh), mesh, corrected)
        rows.append({
            "aspect": aspect,
            "production_wall": a["wall_fraction"], "corrected_wall": b["wall_fraction"],
            "wall_ratio": (b["wall_fraction"] / a["wall_fraction"]
                           if a["wall_fraction"] > 0 else float("inf")),
            "production_floor": a["floor_fraction"], "corrected_floor": b["floor_fraction"],
            "floor_ratio": (b["floor_fraction"] / a["floor_fraction"]
                            if a["floor_fraction"] > 0 else float("inf"))})
        print(f"  C AR{aspect:6.0f}  wall {a['wall_fraction']:.4f} -> {b['wall_fraction']:.4f}"
              f"  (x{rows[-1]['wall_ratio']:.3f})   floor {a['floor_fraction']:.4f}"
              f" -> {b['floor_fraction']:.4f} (x{rows[-1]['floor_ratio']:.3f})", flush=True)
    return rows


def closure_regression():
    """P1a regression: the production lift must carry the published width.

    EXP C measured the production lift discarding exactly sqrt(2) of the
    published planar width (0.5893 deg lifted vs 0.8334 deg published).  This
    is the cheap standing check that the fix stays in — beam statistics only,
    no transport gather.
    """
    iead = load_krueger_2024_digitized_iead(_DEFAULT_DATA)
    tangent = np.tan(np.deg2rad(iead.signed_angle_deg))
    published = np.degrees(np.arctan(np.sqrt(np.average(
        tangent ** 2, weights=iead.probability_weight))))
    production = krueger_iead_species(
        azimuthal_order=16, angle_bin_deg=0.25, energy_bin_eV=250.0)
    velocity = np.asarray(production.velocity_sqrt_eV, dtype=float)
    weight = np.asarray(production.weight, dtype=float)
    weight = weight / weight.sum()
    lifted = np.degrees(np.arctan(np.sqrt(np.average(
        (velocity[:, 0] / velocity[:, 2]) ** 2, weights=weight))))
    row = {"published_planar_deg": float(published),
           "lifted_planar_deg": float(lifted),
           "deficit_ratio": float(published / lifted),
           "pre_fix_deficit_ratio": 1.4141}
    print(f"  published planar {published:.4f} deg; lifted planar {lifted:.4f} deg;"
          f" deficit x{row['deficit_ratio']:.4f} (was x1.4141)", flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="results/curated/angular_convergence_p0")
    parser.add_argument("--skip-a", action="store_true")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("GATE view-factor (collimated, AR 30):", flush=True)
    payload = {"gate": {}}
    for theta in (0.5, 1.0, 1.5):
        observed, expected = view_factor_gate(30.0, theta)
        error = abs(observed - expected)
        print(f"  theta={theta}deg observed={observed:.6f} expected={expected:.6f}"
              f" |err|={error:.2e}", flush=True)
        payload["gate"][f"theta_{theta}"] = {
            "observed": observed, "expected": expected, "error": error}

    if not args.skip_a:
        print("EXP A quadrature ablation:", flush=True)
        payload["experiment_a"] = experiment_a()
    print("EXP B mouth question (AR 9):", flush=True)
    payload["experiment_b"] = experiment_b()
    print("EXP C azimuthal-closure error vs AR:", flush=True)
    payload["experiment_c"] = experiment_c()
    print("P1a closure regression:", flush=True)
    payload["closure_regression"] = closure_regression()

    with open(os.path.join(args.out, "angular_convergence.json"), "w") as handle:
        json.dump(payload, handle, indent=2, default=float)
    print(f"wrote {args.out}/angular_convergence.json", flush=True)
    return payload


if __name__ == "__main__":
    main()
