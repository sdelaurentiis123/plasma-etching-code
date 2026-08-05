"""Grade the reflection cascade's etch-front delivery against Huang's funnelling curve.

Huang's thesis reports the observable directly (tmp/pdfs/huang_thesis.txt):

  L5405-5407: "This shadowing contributes to a decrease in ion flux to the etch
  front from 2.0 x 10^15 to 0.3 x 10^15 cm-2s-1."  (oxide AR 0 -> 40)

  L5408-5413: "The flux of hot neutrals to the etch front increases from
  3.1 x 10^15 to 8.0 x 10^15 cm-2s-1 as the etch depth increases from 0 to
  480 nm (AR = 4). ... As the etch depth increases from 480 to 4,800 nm
  (AR = 40), the flux of hot neutrals to the etch front decreases to
  1.1 x 10^15 cm-2s-1, which is due to diffusive scattering from the sidewalls
  and thermalization of the hot particles following several collisions."

  L5399-5402: "Each strike of the etch front by a neutral particle (either hot
  or thermal) increments the flux-count."  (re-arrival counting)

  L5578-5580: "There is a small initial increase in etch rate which is due to
  tapering of the feature that funnels hot neutrals to the etch front.  The
  etch rate then decreases by 80% by the time the AR reaches 40."

His base case has a photoresist AR of 13 above the oxide, so his "oxide AR 0"
point already sits under a deep mask -- the same situation as the Krueger cell
(850 nm mask over a 90 nm opening, mask AR 9.4).  The comparison is therefore
made on ratios normalised at oxide AR 0, which cancels the differing cell
widths and source magnitudes.

The scan builds straight (unnecked) trenches at a series of etched depths, runs
one frozen-geometry transport gather each with the cascade active, and reports
the delivery to the floor band split into direct ions and cascade hot neutrals.
No fitting: the published ratios are the grading reference.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from petch.threed import reinit_narrow                       # noqa: E402
from petch.feature_step_3d import FeatureGeometry3D          # noqa: E402
from petch.surface_kinetics import FaceResolvedEnergeticFlux  # noqa: E402
from mouth_equilibrium_probe import (                        # noqa: E402
    CELL_WIDTH_UM, CELL_LENGTH_UM, MASK_THICKNESS_UM, OPEN_WIDTH_UM,
    HEADSPACE_UM, gather_transport, _boundary_and_mechanism,
)
from petch.feature_step_3d import advance_feature_step_3d     # noqa: E402


def gather_transport_pilot(geometry, *, transport_device="cpu",
                           ion_azimuthal_order=16, face_quadrature_points=3):
    """One frozen gather under the PRODUCTION pilot transport configuration.

    The pilot routes thermal neutrals through the deterministic extruded
    radiosity exchange and leaves the ballistic (energetic) gather to its own
    default, i.e. it does not pass ``ballistic_periodic_lateral``.  Running the
    identical geometry through both configurations isolates backend
    differences in etch-front delivery from geometry and chemistry.
    """
    from time import perf_counter
    boundary, mechanism, role, realized, source_z = _boundary_and_mechanism(
        geometry, ion_azimuthal_order=ion_azimuthal_order)
    radiosity = {
        "form_factor_backend": "deterministic_extruded_2d",
        "periodic_lateral": True,
        "domain_size": tuple(float(v) for v in realized),
        "relative_tolerance": 1e-4,
        "maximum_iterations": 200,
        "deterministic_extruded_options": {
            "exchange_method": "analytic_occlusion",
            "exchange_relative_tolerance": 1e-4,
            "maximum_refinement_level": 24,
        },
    }
    started = perf_counter()
    result = advance_feature_step_3d(
        geometry, boundary, role, mechanism, etchable_material_ids=(1, 2),
        duration_s=0.0, source_bounds=(0.0, realized[0], 0.0, realized[1]),
        source_z=source_z, n_position=1, seed=1, cfl_number=0.25,
        reinitialize=True, reinitialization_method="cr2",
        profile_periodic_lateral=True,
        transport_device=str(transport_device),
        neutral_radiosity_options=radiosity,
        ballistic_transport="face_gather", grazing_ion_reflection={},
        ballistic_face_quadrature_points=int(face_quadrature_points),
        topology_change_policy="continue_gas_cavity",
        surface_state_remap_backend="common_refinement")
    return result, mechanism, role, perf_counter() - started

# Huang thesis Fig. 6.6(a) anchors, normalised at oxide AR 0 (see docstring).
HUANG = {
    "ion_ar0_cm2s": 2.0e15,
    "ion_ar40_cm2s": 0.3e15,
    "hot_ar0_cm2s": 3.1e15,
    "hot_ar4_cm2s": 8.0e15,
    "hot_ar40_cm2s": 1.1e15,
}
HUANG_RATIOS = {
    "hot_peak_over_ar0": HUANG["hot_ar4_cm2s"] / HUANG["hot_ar0_cm2s"],   # 2.58
    "ion_ar40_over_ar0": HUANG["ion_ar40_cm2s"] / HUANG["ion_ar0_cm2s"],  # 0.15
    "hot_ar40_over_ar0": HUANG["hot_ar40_cm2s"] / HUANG["hot_ar0_cm2s"],  # 0.355
    "total_ar40_over_ar0": ((HUANG["ion_ar40_cm2s"] + HUANG["hot_ar40_cm2s"])
                            / (HUANG["ion_ar0_cm2s"] + HUANG["hot_ar0_cm2s"])),
    "hot_over_ion_ar0": HUANG["hot_ar0_cm2s"] / HUANG["ion_ar0_cm2s"],    # 1.55
    "hot_over_ion_ar40": HUANG["hot_ar40_cm2s"] / HUANG["ion_ar40_cm2s"],  # 3.67
}
FLOOR_BAND_UM = 0.03      # floor faces within this height of the etch front


def make_straight_trench_geometry_3d(
        *, etched_depth, dx, cell_width=CELL_WIDTH_UM,
        cell_length=CELL_LENGTH_UM, opening_width=OPEN_WIDTH_UM,
        mask_thickness=MASK_THICKNESS_UM, floor_pad=0.10,
        headspace=HEADSPACE_UM, mesh_length_unit_m=1e-6,
        substrate_material_id=1, mask_material_id=2):
    """Straight-walled trench etched ``etched_depth`` into the substrate.

    Same cell as the Krueger validation (0.09 um opening under a 0.85 um mask);
    only the etch depth varies, so the sweep isolates aspect-ratio transport.
    """
    substrate_top = float(etched_depth) + float(floor_pad)
    domain_height = substrate_top + mask_thickness + headspace
    shape = tuple(max(3, int(round(length / dx)) + 1)
                  for length in (cell_width, cell_length, domain_height))
    x, y, z = (np.arange(size) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    radius = np.abs(X - 0.5 * cell_width)
    floor = substrate_top - float(etched_depth)
    mask_top = substrate_top + mask_thickness
    half_open = 0.5 * opening_width

    base = floor - Z
    substrate_wall = np.minimum(np.minimum(Z - floor, substrate_top - Z),
                                radius - half_open)
    substrate_levelset = np.maximum(base, substrate_wall)
    mask_levelset = np.minimum(np.minimum(Z - substrate_top, mask_top - Z),
                               radius - half_open)

    substrate_phi = reinit_narrow(substrate_levelset, dx, domain_height + cell_width)
    mask_phi = reinit_narrow(mask_levelset, dx, domain_height + cell_width)
    phi = reinit_narrow(np.maximum(substrate_phi, mask_phi), dx,
                        domain_height + cell_width)

    substrate_solid = (Z < substrate_top) & ~((Z > floor) & (radius < half_open))
    mask_solid = (Z >= substrate_top) & (Z < mask_top) & (radius >= half_open)
    material = np.zeros(shape, dtype=int)
    material[substrate_solid] = int(substrate_material_id)
    material[mask_solid] = int(mask_material_id)
    unlabeled = (phi > 0.0) & (material == 0)
    owner = substrate_levelset >= mask_levelset
    material[unlabeled] = np.where(owner[unlabeled], int(substrate_material_id),
                                   int(mask_material_id))
    return FeatureGeometry3D(
        phi, material, dx, mesh_length_unit_m,
        material_levelsets={int(substrate_material_id): substrate_phi,
                            int(mask_material_id): mask_phi}), floor


def floor_face_mask(result, *, floor_z_um, band=FLOOR_BAND_UM,
                    opening=OPEN_WIDTH_UM, cell_width=CELL_WIDTH_UM):
    """Active faces lying on the etch front: inside the opening, at the floor."""
    centroids = np.asarray(result.active_face_centroid, dtype=float)
    radius = np.abs(centroids[:, 0] - 0.5 * cell_width)
    return ((np.abs(centroids[:, 2] - floor_z_um) <= band)
            & (radius < 0.5 * opening))


def delivery_split(result, keep_mask):
    """Area-weighted flux to the kept faces, split into direct ions / hot neutrals."""
    active = np.asarray(result.active_face_index, dtype=int)
    old_to_new = np.full(len(result.face_material_id), -1, dtype=int)
    old_to_new[active] = np.arange(active.size)
    areas = np.asarray(result.active_face_area, dtype=float)
    kept_area = float(areas[keep_mask].sum())
    out = {"direct_ion": 0.0, "hot_neutral": 0.0,
           "direct_ion_power": 0.0, "hot_neutral_power": 0.0}
    neutral = result.transport.surface_fluxes.neutral_flux_m2_s
    for name, field in dict(neutral).items():
        value = np.asarray(field, dtype=float)
        if value.ndim == 0:
            out["neutral_" + name] = float(value)
            continue
        per_face = value[np.asarray(result.active_face_index, dtype=int)] \
            if value.size == len(result.face_material_id) else value
        out["neutral_" + name] = float(
            (per_face[keep_mask] * areas[keep_mask]).sum())
    for population in result.transport.surface_fluxes.energetic_fluxes:
        is_hot = population.name.endswith(":hot_neutral")
        key = "hot_neutral" if is_hot else "direct_ion"
        if not isinstance(population, FaceResolvedEnergeticFlux):
            continue
        face = old_to_new[np.asarray(population.event_face, dtype=int)]
        keep = (face >= 0)
        face = face[keep]
        keep &= True
        flux = np.asarray(population.event_flux_m2_s, dtype=float)[keep]
        energy = np.asarray(population.event_energy_eV, dtype=float)[keep]
        on_floor = keep_mask[face]
        rate = flux[on_floor] * areas[face[on_floor]]
        out[key] += float(rate.sum())
        out[key + "_power"] += float((rate * energy[on_floor]).sum())
    if kept_area > 0.0:
        for key in list(out):
            out[key] /= kept_area
    out["floor_area_um2"] = kept_area
    return out


def scan(aspect_ratios, *, dx, device, opening=OPEN_WIDTH_UM, pipeline="probe"):
    gather = gather_transport if pipeline == "probe" else gather_transport_pilot
    records = []
    for ar in aspect_ratios:
        etched = float(ar) * opening
        geometry, floor_z = make_straight_trench_geometry_3d(
            etched_depth=etched, dx=dx)
        result, _, _, elapsed = gather(
            geometry, transport_device=device)
        keep = floor_face_mask(result, floor_z_um=floor_z)
        record = {"pipeline": pipeline,
                  "oxide_aspect_ratio": float(ar),
                  "etched_depth_um": etched,
                  "total_aspect_ratio": (etched + MASK_THICKNESS_UM) / opening,
                  "n_floor_faces": int(keep.sum()),
                  "gather_s": elapsed}
        record.update(delivery_split(result, keep))
        record["total_energetic"] = record["direct_ion"] + record["hot_neutral"]
        records.append(record)
        print(f"AR {ar:5.1f}  ion {record['direct_ion']:.4e}  "
              f"hot {record['hot_neutral']:.4e}  "
              f"faces {record['n_floor_faces']:3d}  {elapsed:5.1f}s", flush=True)
    return records


def grade(records):
    base = next(r for r in records if r["oxide_aspect_ratio"] == 0.0)
    for record in records:
        for key in ("direct_ion", "hot_neutral", "total_energetic"):
            denom = base[key]
            record[key + "_over_ar0"] = (record[key] / denom if denom > 0
                                         else float("nan"))
        record["hot_over_ion"] = (record["hot_neutral"] / record["direct_ion"]
                                  if record["direct_ion"] > 0 else float("nan"))
    peak = max(records, key=lambda r: r["hot_neutral_over_ar0"])
    return {
        "huang_reference": HUANG,
        "huang_ratios": HUANG_RATIOS,
        "measured_hot_peak_over_ar0": peak["hot_neutral_over_ar0"],
        "measured_hot_peak_at_ar": peak["oxide_aspect_ratio"],
        "measured_hot_over_ion_ar0": base["hot_neutral"] / base["direct_ion"],
        "records": records,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aspect-ratios", type=float, nargs="+",
                        default=[0.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0])
    parser.add_argument("--dx-um", type=float, default=0.01)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--pipeline", choices=("probe", "pilot"), default="probe")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "curated"
                        / "cascade_funnelling" / "scan.json")
    args = parser.parse_args(argv)

    records = scan(args.aspect_ratios, dx=args.dx_um, device=args.device,
                   pipeline=args.pipeline)
    payload = grade(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2))
    print(json.dumps({k: v for k, v in payload.items() if k != "records"},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
