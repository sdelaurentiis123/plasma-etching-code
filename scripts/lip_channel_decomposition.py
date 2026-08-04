#!/usr/bin/env python3
"""Per-channel decomposition of the mask-lip film budget (Track A, 2026-08-04).

`RESULTS_LIP_DEPOSITION_AUDIT_2026-08-04.md` established, with the transport
defect fixed, that the deposition side is faithful (visibility <= 1 and within
7-15% of an analytic slot reference, published sticking rows firing, O/depositor
isotropy exact) while the top band still sits at removal/deposition = 0.216
against 0.88-0.99 in the bands that match Krueger.

This script attributes that removal to its channels by *ablation* rather than by
re-deriving the laws: one frozen-geometry gather is reused, and the surface
mechanism is relaxed four times on modified flux fields --

    full            all channels
    no_oxygen       O flux zeroed
    no_hot_neutral  the reflection cascade's secondary population dropped
    no_energetic    every energetic population dropped

so that

    O           = removal(full) - removal(no_oxygen)
    hot neutral = removal(full) - removal(no_hot_neutral)
    primary ion = removal(no_hot_neutral) - removal(no_energetic)

Ablation uses the module itself, so no published law is reimplemented here and
the attribution cannot drift from what the engine actually integrates.

It also reports, per band, the energetic flux and mean incidence split by
population (does the cascade deliver hot neutrals to the mask top at all?) and
a face-level listing of the convex mask corner (do chamfer facets exist there
to receive the peak of the angular yield?).

Usage:
    python scripts/lip_channel_decomposition.py --neck-nm 45 --dx 0.01
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
    _select_surface_fluxes,
    _surface_gas_normals,
)
from petch.surface_kinetics import SurfaceFluxes  # noqa: E402

from mouth_equilibrium_probe import (  # noqa: E402
    MASK_THICKNESS_UM,
    SUBSTRATE_TOP_UM,
    _boundary_and_mechanism,
    gather_transport,
    make_necked_trench_geometry_3d,
)

BANDS = ((0.0, 50.0), (50.0, 100.0), (100.0, 150.0), (200.0, 270.0))
HOT_SUFFIX = ":hot_neutral"


def variant_fluxes(surface_fluxes, *, drop_oxygen=False, drop_hot=False,
                   drop_energetic=False):
    """A copy of the gathered field with whole channels removed."""
    neutral = {name: np.array(value, dtype=float, copy=True)
               for name, value in surface_fluxes.neutral_flux_m2_s.items()}
    if drop_oxygen:
        for name in list(neutral):
            if name == "O":
                neutral[name][:] = 0.0
    populations = tuple(surface_fluxes.energetic_fluxes)
    if drop_energetic:
        populations = ()
    elif drop_hot:
        populations = tuple(p for p in populations
                            if not p.name.endswith(HOT_SUFFIX))
    return SurfaceFluxes(neutral, populations)


def relax(mechanism, flux, material, *, relax_s=0.4, rounds=5, tolerance=5e-3):
    state = mechanism.initial_state_by_material(material)
    previous = None
    step = None
    for _ in range(int(rounds)):
        step = mechanism.advance_by_material(state, flux, float(relax_s), material)
        net = (np.asarray(step.etch_velocity_m_s, dtype=float)
               - np.asarray(step.normal_growth_velocity_m_s, dtype=float))
        if previous is not None:
            scale = max(float(np.max(np.abs(net))), 1e-30)
            if float(np.max(np.abs(net - previous))) / scale <= float(tolerance):
                state = step.state
                break
        previous = net
        state = step.state
    return step


def mask_arrays(step, material):
    """Deposition/removal per mask face, in polymer units, aligned to mask faces."""
    result = step.material_results.get(2)
    count = int((material == 2).sum())
    def field(name):
        value = np.asarray(getattr(result, name), dtype=float)
        return np.broadcast_to(value, (count,))
    return field("deposited_polymer_units_m2"), field("removed_polymer_units_m2")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--neck-nm", type=float, default=45.0)
    parser.add_argument("--dx", type=float, default=0.01)
    parser.add_argument("--relax-s", type=float, default=0.4)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--output", default="results/curated/lip_channel_decomposition")
    args = parser.parse_args()

    geometry = make_necked_trench_geometry_3d(
        neck_width_um=float(args.neck_nm) * 1e-3, dx=float(args.dx))
    result, mechanism, role, elapsed = gather_transport(
        geometry, transport_device=str(args.transport_device))
    active = np.asarray(result.active_face_index, dtype=int)
    material = np.asarray(result.face_material_id, dtype=int)[active]
    n_faces = len(result.face_material_id)

    verts, faces, centroids_all, _, _ = _extract_uniform_surface_arrays(geometry)
    normals = np.asarray(
        _surface_gas_normals(verts, faces, centroids_all, geometry),
        dtype=float)[active]
    centroid = np.asarray(result.active_face_centroid, dtype=float)
    scale = geometry.dx if centroid[:, 2].max() > 10.0 else 1.0
    z_um = centroid[:, 2] * scale
    mask_top = SUBSTRATE_TOP_UM + MASK_THICKNESS_UM
    depth_nm = (mask_top - z_um) * 1e3
    tilt_deg = np.degrees(np.arcsin(np.clip(np.abs(normals[:, 2]), 0.0, 1.0)))
    lateral = (np.abs(normals[:, 0]) > np.abs(normals[:, 2])) & (material == 2)

    variants = {
        "full": dict(),
        "no_oxygen": dict(drop_oxygen=True),
        "no_hot_neutral": dict(drop_hot=True),
        "no_energetic": dict(drop_energetic=True),
    }
    per_variant = {}
    for name, kwargs in variants.items():
        selected = _select_surface_fluxes(
            variant_fluxes(result.transport.surface_fluxes, **kwargs),
            active, n_faces, role)
        step = relax(mechanism, selected, material, relax_s=float(args.relax_s))
        deposited, removed = mask_arrays(step, material)
        net = (np.asarray(step.etch_velocity_m_s, dtype=float)
               - np.asarray(step.normal_growth_velocity_m_s, dtype=float))
        per_variant[name] = dict(deposited=deposited, removed=removed,
                                 net=np.broadcast_to(net, (active.size,)))

    mask_faces = np.where(material == 2)[0]
    mask_lateral = lateral[mask_faces]
    mask_depth = depth_nm[mask_faces]
    area = np.asarray(result.active_face_area, dtype=float)[mask_faces]

    # Energetic flux per population, per face.
    populations = {}
    for population in result.transport.surface_fluxes.energetic_fluxes:
        flux = np.zeros(n_faces)
        cosw = np.zeros(n_faces)
        np.add.at(flux, population.event_face,
                  population.event_flux_m2_s)
        np.add.at(cosw, population.event_face,
                  population.event_flux_m2_s * population.event_cosine_incidence)
        populations[population.name] = (flux[active][mask_faces],
                                        np.divide(cosw[active][mask_faces],
                                                  np.maximum(flux[active][mask_faces], 1e-300)))

    rows = []
    for low, high in BANDS:
        band = mask_lateral & (mask_depth >= low) & (mask_depth < high)
        if not band.any():
            continue
        weight = area[band]
        total = float(weight.sum())

        def wsum(values):
            return float((np.asarray(values)[band] * weight).sum())

        dep = wsum(per_variant["full"]["deposited"])
        rem_full = wsum(per_variant["full"]["removed"])
        rem_no_o = wsum(per_variant["no_oxygen"]["removed"])
        rem_no_hot = wsum(per_variant["no_hot_neutral"]["removed"])
        rem_no_e = wsum(per_variant["no_energetic"]["removed"])
        net_nm_s = float((np.asarray(per_variant["full"]["net"])[mask_faces][band]
                          * weight).sum() / total) * 1e9
        row = dict(
            band_nm=[low, high], faces=int(band.sum()),
            mean_tilt_deg=float((tilt_deg[mask_faces][band] * weight).sum() / total),
            deposition_units=dep,
            removal_units=rem_full,
            removal_over_deposition=rem_full / dep if dep else float("nan"),
            oxygen_share=(rem_full - rem_no_o) / dep if dep else float("nan"),
            hot_neutral_share=(rem_full - rem_no_hot) / dep if dep else float("nan"),
            primary_ion_share=(rem_no_hot - rem_no_e) / dep if dep else float("nan"),
            residual_share=rem_no_e / dep if dep else float("nan"),
            net_nm_s=net_nm_s,
        )
        for name, (flux, cosine) in populations.items():
            row[f"flux[{name}]"] = float((flux[band] * weight).sum() / total)
            row[f"cos[{name}]"] = float((cosine[band] * flux[band] * weight).sum()
                                        / max(float((flux[band] * weight).sum()), 1e-300))
        rows.append(row)

    corner = np.argsort(mask_depth)[:14]
    corner_rows = [dict(depth_nm=float(mask_depth[i]),
                        tilt_deg=float(tilt_deg[mask_faces][i]),
                        lateral=bool(mask_lateral[i]),
                        net_nm_s=float(np.asarray(per_variant["full"]["net"])[mask_faces][i]) * 1e9,
                        removal_over_deposition=float(
                            per_variant["full"]["removed"][i]
                            / max(per_variant["full"]["deposited"][i], 1e-300)))
                   for i in corner]

    report = dict(neck_nm=float(args.neck_nm), dx=float(args.dx),
                  gather_s=float(elapsed), bands=rows, corner_faces=corner_rows)
    out = ROOT / str(args.output)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"channels_neck{int(args.neck_nm)}_dx{args.dx}.json").write_text(
        json.dumps(report, indent=2))

    print(f"gather {elapsed:.1f}s  mask lateral faces {int(mask_lateral.sum())}")
    header = ("band        tilt   rem/dep    O     hot   prim   resid   net nm/s"
              "   ion flux   hot flux")
    print(header)
    for row in rows:
        print("%3.0f-%-3.0f  %6.2f  %7.3f %6.3f %6.3f %6.3f %6.3f %9.3f %10.2e %10.2e"
              % (row["band_nm"][0], row["band_nm"][1], row["mean_tilt_deg"],
                 row["removal_over_deposition"], row["oxygen_share"],
                 row["hot_neutral_share"], row["primary_ion_share"],
                 row["residual_share"], row["net_nm_s"],
                 row.get("flux[ions]", 0.0), row.get("flux[ions:hot_neutral]", 0.0)))
    print("\ncorner faces (shallowest first)")
    for row in corner_rows:
        print("  depth %6.1f nm  tilt %5.1f  lateral %-5s  rem/dep %7.3f  net %8.3f nm/s"
              % (row["depth_nm"], row["tilt_deg"], row["lateral"],
                 row["removal_over_deposition"], row["net_nm_s"]))


if __name__ == "__main__":
    main()
