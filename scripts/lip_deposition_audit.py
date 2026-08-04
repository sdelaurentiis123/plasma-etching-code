#!/usr/bin/env python3
"""Deposition-side audit of the mask-lip film budget (Track A, 2026-08-04).

`RESULTS_LIP_REMOVAL_AUDIT_2026-08-04.md` falsified the angular ion-removal law
as the cause of the top-band closure excess: no angular function can balance a
near-vertical lip, so the residual ~0.8 x deposition must live on the
*deposition* side -- either the depositor flux delivered to near-vertical faces
or the effective sticking applied there.

This audit measures both, per depth band, against analytic references that need
no free parameters:

* **Delivered flux vs unobstructed arrival.**  The neutral boundary carries the
  exact half-Maxwellian *flux* angular marginal ``p(mu, phi) = 2 mu / 2 pi``
  (`reactor_boundary._direction_marginalized_thermal_flux_species`).  For a face
  of normal ``n`` with no obstruction, the transport must deliver

      A(n) = sum_s w_s max(-d_s . n, 0) / |d_s,z|

  times the source flux, where ``(d_s, w_s)`` are that same quadrature's nodes
  and weights (`boundary_transport_3d`: ``projection = cos_inc / -d_z``).  A(n)
  is 1 for a fully exposed horizontal face and 1/2 for a fully exposed vertical
  wall (analytically ``1/pi * int sqrt(1-mu^2) dmu * int cos phi dphi``).  The
  ratio ``delivered / (source * A(n))`` is therefore a pure *visibility*
  fraction in [0, 1]: greater than one would mean the gather over-delivers.
* **Isotropy consistency.**  Oxygen and the depositors share one angular law, so
  ``(O_delivered / O_source) / (dep_delivered / dep_source)`` must be 1 face by
  face.  Any departure is a transport bug, and it would break the geometry-free
  0.1953 O share that `RESULTS_O_CHANNEL_2026-08-04.md` gates.
* **Effective sticking.**  ``deposited_polymer_units / delivered_depositor_cells``
  recovers which published row actually fires: 0.1 on fresh polymer, 0.02 on
  crosslinked polymer, 0.094 on bare mask, or a coverage/crosslink blend of
  them (`mixed_layer.py` deposition blend).

Bands follow the neck regrade: 0-50 nm below the mask top (the 10.8x excess),
50-100, 100-150, and 200-270 nm (the band that already matches Krueger to
0.88-1.11x -- the control).

Usage:
    python scripts/lip_deposition_audit.py --neck-nm 45 --dx 0.01
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from petch.feature_step_3d import (  # noqa: E402
    _extract_uniform_surface_arrays,
    _surface_gas_normals,
)
from mouth_equilibrium_probe import (  # noqa: E402
    CELL_WIDTH_UM,
    MASK_THICKNESS_UM,
    SUBSTRATE_TOP_UM,
    _boundary_and_mechanism,
    gather_transport,
    make_necked_trench_geometry_3d,
    relax_chemistry,
)

# Depth bands below the mask top, in nm (RESULTS_NECK_REGRADE_2026-08-04).
BANDS = ((0.0, 50.0), (50.0, 100.0), (100.0, 150.0), (200.0, 270.0))

# Published deposition rows (mixed_layer_mechanism.KRUEGER_2024_DEPOSITION_*).
STICKING_ROWS = {
    "on_fresh_polymer": 0.1,
    "on_crosslinked_polymer": 0.02,
    "on_bare_mask": 0.094,
}


def measure_local_aperture(geometry, z_um, mask_top_um):
    """Gas-phase aperture (um) at each face's height, from the level set itself."""
    phi = np.asarray(geometry.phi, dtype=float)
    dx = float(geometry.dx)
    ny = phi.shape[1] // 2
    open_cells = (phi[:, ny, :] <= 0.0).sum(axis=0)          # per z index
    aperture = open_cells * dx
    index = np.clip(np.rint(z_um / dx).astype(int), 0, aperture.size - 1)
    return aperture[index]


def unobstructed_arrival(boundary, normals, species_name):
    """A(n): delivered/source for an unobstructed face, from the boundary's own nodes."""
    species = next(s for s in boundary.species if s.name == species_name)
    velocity = np.asarray(species.velocity_sqrt_eV, dtype=float).copy()
    velocity[:, 2] *= -1.0                      # engine convention: toward the surface
    speed = np.linalg.norm(velocity, axis=1)
    direction = velocity / speed[:, None]
    weight = np.asarray(species.weight, dtype=float)
    cos_inc = np.clip(-np.einsum("sd,fd->sf", direction, normals), 0.0, 1.0)
    projection = cos_inc / (-direction[:, 2, None])
    return np.einsum("s,sf->f", weight, projection)


def analytic_slot_visibility(boundary, species_name, normals, points_um,
                             aperture_um, mask_top_um, cell_width_um):
    """Visibility of an infinite (y-invariant) parallel-walled slot.

    A ray leaving the face point and clearing the mask-top plane inside the
    local aperture escapes to the source; otherwise the opposing lip blocks it.
    This is the geometric reference the ray-traced gather should reproduce on
    near-vertical faces, and it is independent of the mesh's facet structure.
    """
    species = next(s for s in boundary.species if s.name == species_name)
    velocity = np.asarray(species.velocity_sqrt_eV, dtype=float).copy()
    velocity[:, 2] *= -1.0
    speed = np.linalg.norm(velocity, axis=1)
    direction = velocity / speed[:, None]
    weight = np.asarray(species.weight, dtype=float)
    reverse = -direction                       # face -> source
    cos_inc = np.clip(-np.einsum("sd,fd->sf", direction, normals), 0.0, 1.0)
    projection = cos_inc / (-direction[:, 2, None])
    rise = np.maximum(mask_top_um - points_um[:, 2], 0.0)          # (f,)
    with np.errstate(divide="ignore", invalid="ignore"):
        travel = rise[None, :] / np.maximum(reverse[:, 2, None], 1e-12)
    exit_x = points_um[None, :, 0] + travel * reverse[:, 0, None]
    half = 0.5 * aperture_um
    centre = 0.5 * cell_width_um
    visible = (np.abs(exit_x - centre) <= half) & (reverse[:, 2, None] > 0.0)
    return np.einsum("s,sf->f", weight, projection * visible)


def band_report(args):
    geometry = make_necked_trench_geometry_3d(
        neck_width_um=float(args.neck_nm) * 1e-3, dx=float(args.dx))
    result, mechanism, role, elapsed = gather_transport(
        geometry, transport_device=str(args.transport_device))
    step, flux, active, material, trace = relax_chemistry(
        result, mechanism, role, relax_s=float(args.relax_s))

    boundary, _, _, _, _ = _boundary_and_mechanism(geometry)
    source_flux = {s.name: float(s.flux_m2_s) for s in boundary.species}
    depositors = [s.name for s in boundary.species
                  if s.charge_number == 0 and s.name.startswith("C")]
    oxygen = [s.name for s in boundary.species if s.name == "O"]

    # The gas-facing normals the transport itself uses, on the same mesh and
    # face ordering that `active_face_index` indexes.
    verts, faces, centroids_all, _, _ = _extract_uniform_surface_arrays(geometry)
    normals = np.asarray(
        _surface_gas_normals(verts, faces, centroids_all, geometry),
        dtype=float)[active]
    centroid = np.asarray(result.active_face_centroid, dtype=float)
    area = np.asarray(result.active_face_area, dtype=float)
    scale = geometry.dx if centroid[:, 2].max() > 10.0 else 1.0
    z_um = centroid[:, 2] * scale
    mask_top = SUBSTRATE_TOP_UM + MASK_THICKNESS_UM
    depth_nm = (mask_top - z_um) * 1e3

    neutral = result.transport.surface_fluxes.neutral_flux_m2_s
    delivered_dep = np.zeros(active.size)
    source_dep = 0.0
    for name in depositors:
        delivered_dep += np.asarray(neutral[name], dtype=float)[active]
        source_dep += source_flux[name]
    delivered_o = np.zeros(active.size)
    source_o = 0.0
    for name in oxygen:
        delivered_o += np.asarray(neutral[name], dtype=float)[active]
        source_o += source_flux[name]

    # Unobstructed arrival factor A(n) per face (identical for every neutral,
    # they share the angular law; use the first depositor's nodes).
    arrival = unobstructed_arrival(boundary, normals, depositors[0])

    mask_result = step.material_results.get(2)
    mask_faces = np.where(material == 2)[0]
    order = {int(face): index for index, face in enumerate(mask_faces)}

    def mask_field(name):
        value = np.asarray(getattr(mask_result, name), dtype=float)
        return np.broadcast_to(value, (mask_faces.size,))

    deposited = mask_field("deposited_polymer_units_m2")
    removed_polymer = mask_field("removed_polymer_units_m2")
    state = mask_result.state
    film_total = (np.asarray(state.n_c_film, dtype=float)
                  + np.asarray(state.n_f_film, dtype=float))
    xl_fraction = np.where(film_total > 0.0,
                           np.asarray(state.n_xl_film, dtype=float)
                           / np.maximum(film_total, 1e-300), 0.0)
    net = (np.asarray(step.etch_velocity_m_s, dtype=float)
           - np.asarray(step.normal_growth_velocity_m_s, dtype=float))

    # Lateral (aperture-setting) mask faces: normal dominated by x.
    lateral = (np.abs(normals[:, 0]) > np.abs(normals[:, 2])) & (material == 2)
    tilt_deg = np.degrees(np.arcsin(np.clip(np.abs(normals[:, 2]), 0.0, 1.0)))

    # Controls: the fully exposed mask top (must read delivered/source ~ 1) and
    # the trench floor (deeply shadowed).  These bound the transport's
    # normalisation independently of any wall geometry.
    horizontal = np.abs(normals[:, 2]) > np.abs(normals[:, 0])
    top_face = horizontal & (np.abs(depth_nm) <= 15.0) & (normals[:, 2] > 0.0)
    floor_face = horizontal & (depth_nm > 1000.0) & (normals[:, 2] > 0.0)
    controls = {}
    for label, select in (("mask_top", top_face), ("trench_floor", floor_face)):
        if not np.any(select):
            controls[label] = None
            continue
        w = area[select]
        controls[label] = {
            "face_count": int(np.count_nonzero(select)),
            "delivered_over_source_depositor": float(
                np.average(delivered_dep[select] / source_dep, weights=w)),
            "unobstructed_arrival_A": float(np.average(arrival[select], weights=w)),
            "visibility_fraction": float(np.average(
                (delivered_dep[select] / source_dep)
                / np.maximum(arrival[select], 1e-300), weights=w)),
        }

    # Analytic slot reference for the lateral faces.
    aperture_um = measure_local_aperture(geometry, z_um, mask_top)
    slot = analytic_slot_visibility(
        boundary, depositors[0], normals,
        np.column_stack((centroid[:, 0] * scale, centroid[:, 1] * scale, z_um)),
        aperture_um, mask_top, CELL_WIDTH_UM)

    bands = []
    for low, high in BANDS:
        select = lateral & (depth_nm >= low) & (depth_nm < high)
        if not np.any(select):
            bands.append({"band_nm": [low, high], "face_count": 0})
            continue
        w = area[select]
        local = np.asarray([order[int(f)] for f in np.where(select)[0]], dtype=int)
        dep_ratio = delivered_dep[select] / source_dep
        o_ratio = delivered_o[select] / source_o
        visibility = dep_ratio / np.maximum(arrival[select], 1e-300)
        with np.errstate(invalid="ignore", divide="ignore"):
            isotropy = o_ratio / np.maximum(dep_ratio, 1e-300)
            # Depositor arrival expressed in film cells: the mechanism counts
            # C-carrying precursor cells, so compare like with like via the
            # deposited units actually produced per delivered depositor flux.
            effective_sticking = np.where(
                delivered_dep[select] > 0.0,
                deposited[local] / np.maximum(delivered_dep[select], 1e-300),
                np.nan)
        bands.append({
            "band_nm": [low, high],
            "face_count": int(np.count_nonzero(select)),
            "wall_tilt_deg_mean": float(np.average(tilt_deg[select], weights=w)),
            "delivered_over_source_depositor": float(np.average(dep_ratio, weights=w)),
            "unobstructed_arrival_A": float(np.average(arrival[select], weights=w)),
            "visibility_fraction": float(np.average(visibility, weights=w)),
            "visibility_max": float(np.max(visibility)),
            "analytic_slot_delivered_over_source": float(
                np.average(slot[select], weights=w)),
            "measured_over_analytic_slot": float(
                np.average(dep_ratio, weights=w)
                / max(np.average(slot[select], weights=w), 1e-300)),
            "isotropy_ratio_O_over_depositor": float(np.average(isotropy, weights=w)),
            "effective_sticking": float(np.average(effective_sticking, weights=w)),
            "crosslinked_fraction": float(np.average(xl_fraction[local], weights=w)),
            "deposited_units_m2_s": float(np.average(deposited[local], weights=w)),
            "removed_polymer_units_m2_s": float(np.average(removed_polymer[local], weights=w)),
            "removal_over_deposition": float(
                np.average(removed_polymer[local], weights=w)
                / max(np.average(deposited[local], weights=w), 1e-300)),
            "net_velocity_nm_s": float(np.average(net[select], weights=w)) * 1e9,
        })

    return {
        "neck_nm": float(args.neck_nm),
        "dx_um": float(args.dx),
        "gather_seconds": float(elapsed),
        "relax_trace": [float(v) for v in trace],
        "source_flux_m2_s": {"depositors": source_dep, "oxygen": source_o},
        "depositor_species": depositors,
        "sticking_rows": STICKING_ROWS,
        "controls": controls,
        "bands": bands,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neck-nm", type=float, default=45.0)
    parser.add_argument("--dx", type=float, default=0.01)
    parser.add_argument("--relax-s", type=float, default=0.4)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "curated" / "lip_deposition_audit")
    args = parser.parse_args(argv)

    payload = band_report(args)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / f"audit_neck{int(args.neck_nm)}_dx{args.dx}.json").write_text(
        json.dumps(payload, indent=2) + "\n")

    header = (f"{'band(nm)':>12} {'tilt':>6} {'dep/src':>9} {'A(n)':>7} "
              f"{'vis':>7} {'O/dep':>7} {'s_eff':>8} {'x_xl':>6} "
              f"{'meas/slot':>9} {'rem/dep':>8} {'net(nm/s)':>10}")
    print(header)
    for label, value in (payload.get("controls") or {}).items():
        if value:
            print(f"  control {label:>12}: delivered/source="
                  f"{value['delivered_over_source_depositor']:.4f} "
                  f"A(n)={value['unobstructed_arrival_A']:.4f} "
                  f"visibility={value['visibility_fraction']:.4f}")
    for band in payload["bands"]:
        if not band.get("face_count"):
            print(f"{str(band['band_nm']):>12} (no faces)")
            continue
        print(f"{band['band_nm'][0]:5.0f}-{band['band_nm'][1]:<6.0f} "
              f"{band['wall_tilt_deg_mean']:6.2f} "
              f"{band['delivered_over_source_depositor']:9.4f} "
              f"{band['unobstructed_arrival_A']:7.4f} "
              f"{band['visibility_fraction']:7.4f} "
              f"{band['isotropy_ratio_O_over_depositor']:7.4f} "
              f"{band['effective_sticking']:8.4f} "
              f"{band['crosslinked_fraction']:6.3f} "
              f"{band['measured_over_analytic_slot']:8.3f} "
              f"{band['removal_over_deposition']:8.4f} "
              f"{band['net_velocity_nm_s']:10.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
