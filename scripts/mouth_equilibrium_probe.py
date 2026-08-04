#!/usr/bin/env python3
"""Frozen-geometry interrogation of the Krüger trench neck (mouth) equilibrium.

The 60 s evolution runs answer "where does the mouth end up" only after five
hours of wall clock, and they conflate transport, chemistry and level-set
kinematics.  This probe asks the equilibrium question directly:

    for a trench whose neck aperture is *prescribed* at w, does the neck film
    grow (aperture closes) or erode (aperture opens)?

The aperture where the net normal velocity changes sign is the equilibrium
aperture the evolution run would converge to.  Bisecting that sign change costs
one transport gather per geometry (~75 s on CPU) instead of a full campaign.

Three experiments:

1. ``sweep``      net neck velocity over prescribed apertures, bisected.
2. ``budget``     the lip removal/deposition decomposition at the experimental
                  neck (39 nm), including removal versus ion incidence angle
                  against a cosine reference (You et al., *Coatings* 13, 1452
                  (2023), Fig. 6b: measured normalised etch rates lie *above*
                  cosine out to 50-60 deg).
3. ``resolution`` the same 39 nm evaluation at a second dx, which separates a
                  physical neck balance from a discretisation floor (one cell
                  of film per side is 2*dx of aperture).

Chemistry is relaxed to its local steady state on the frozen mesh: the transport
gather runs once at ``duration_s=0`` (no geometry motion, no topology risk) and
the surface mechanism is then integrated on that fixed flux field until the
face velocities stop changing.  Geometry never moves, so the reported velocity
is the equilibrium balance at exactly the prescribed aperture.

Geometry follows the digitised Fig. 7 profile (RESEARCH_MOUTH_MECHANISM_KRUEGER
_2026-08-02.md): 90 nm at the mask top, a smooth constriction to the prescribed
aperture ~250 nm into the 850 nm mask, re-opening to 90 nm below it, over a
300 nm etched oxide trench.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.feature_step_3d import (
    FeatureGeometry3D,
    advance_feature_step_3d,
    _select_surface_fluxes,
)
from petch.threed import reinit_narrow
from petch.reactor_boundary import build_krueger_2024_development_boundary
from petch.surface_kinetics import FaceResolvedEnergeticFlux

DATA = ROOT / "data" / "experimental" / "krueger_2024"

# Digitised Fig. 7 geometry (RESEARCH_MOUTH_MECHANISM_KRUEGER_2026-08-02 §1.5).
CELL_WIDTH_UM = 0.13
CELL_LENGTH_UM = 0.02
MASK_THICKNESS_UM = 0.85
OPEN_WIDTH_UM = 0.09
NECK_DEPTH_UM = 0.25          # below the mask top; SEM minimum 0.20, MCFPM 0.271
NECK_AXIAL_SIGMA_UM = 0.10    # simulated neck is ~240 nm axially at 1.5x minimum
SUBSTRATE_TOP_UM = 0.40
ETCHED_DEPTH_UM = 0.30
HEADSPACE_UM = 0.15

EXPERIMENTAL_NECK_NM = 39.0   # SEM minimum aperture
TARGET_MOUTH_NM = 45.0        # Krüger's fitted w_m target


def make_necked_trench_geometry_3d(
        *, neck_width_um, dx, cell_width=CELL_WIDTH_UM, cell_length=CELL_LENGTH_UM,
        opening_width=OPEN_WIDTH_UM, mask_thickness=MASK_THICKNESS_UM,
        substrate_top=SUBSTRATE_TOP_UM, etched_depth=ETCHED_DEPTH_UM,
        neck_depth=NECK_DEPTH_UM, neck_sigma=NECK_AXIAL_SIGMA_UM,
        headspace=HEADSPACE_UM, mesh_length_unit_m=1e-6,
        substrate_material_id=1, mask_material_id=2):
    """Rectangular trench whose mask aperture is constricted to ``neck_width_um``.

    Identical in construction to ``make_rectangular_trench_geometry_3d`` except
    that the mask half-width varies axially: the aperture is ``opening_width``
    at the mask top and bottom and narrows smoothly to ``neck_width_um`` at
    ``neck_depth`` below the mask top.  The constriction is labelled mask
    material, which routes it to the a-C mechanism whose fluorocarbon film is
    the physical neck.
    """
    domain_height = substrate_top + mask_thickness + headspace
    if not (0.0 < neck_width_um <= opening_width < cell_width):
        raise ValueError("neck width must be positive and no wider than the opening")
    shape = tuple(max(3, int(round(length / dx)) + 1)
                  for length in (cell_width, cell_length, domain_height))
    x, y, z = (np.arange(size) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    radius = np.abs(X - 0.5 * cell_width)
    floor = substrate_top - etched_depth
    mask_top = substrate_top + mask_thickness

    base = floor - Z
    substrate_wall_slab = np.minimum(Z - floor, substrate_top - Z)
    substrate_wall = np.minimum(substrate_wall_slab, radius - 0.5 * opening_width)
    substrate_levelset = np.maximum(base, substrate_wall)

    # Axial aperture profile: Gaussian constriction centred NECK_DEPTH below the
    # mask top, floored at the prescribed neck width.
    z_neck = mask_top - neck_depth
    bump = np.exp(-(((Z - z_neck) / neck_sigma) ** 2))
    half_width = 0.5 * (opening_width - (opening_width - neck_width_um) * bump)
    mask_slab = np.minimum(Z - substrate_top, mask_top - Z)
    mask_levelset = np.minimum(mask_slab, radius - half_width)

    substrate_phi = reinit_narrow(substrate_levelset, dx, domain_height + cell_width)
    mask_phi = reinit_narrow(mask_levelset, dx, domain_height + cell_width)
    phi = reinit_narrow(np.maximum(substrate_phi, mask_phi), dx,
                        domain_height + cell_width)

    substrate_solid = (Z < substrate_top) & ~(
        (Z > floor) & (radius < 0.5 * opening_width))
    mask_solid = (Z >= substrate_top) & (Z < mask_top) & (radius >= half_width)
    material = np.zeros(shape, dtype=int)
    material[substrate_solid] = int(substrate_material_id)
    material[mask_solid] = int(mask_material_id)
    unlabeled_solid = (phi > 0.0) & (material == 0)
    substrate_owner = substrate_levelset >= mask_levelset
    material[unlabeled_solid] = np.where(
        substrate_owner[unlabeled_solid],
        int(substrate_material_id), int(mask_material_id))
    return FeatureGeometry3D(
        phi, material, dx, mesh_length_unit_m,
        material_levelsets={
            int(substrate_material_id): substrate_phi,
            int(mask_material_id): mask_phi,
        })


ANGLE_BAND_LO_UM = 0.015      # measurement band above the apex (excludes it)
ANGLE_BAND_HI_UM = 0.055
ANGLE_BAND_MID_UM = 0.035     # aperture is pinned to the reference width here


def make_sloped_wall_geometry_3d(
        *, wall_angle_deg, aperture_nm, dx, cell_width=CELL_WIDTH_UM,
        cell_length=CELL_LENGTH_UM, opening_width=OPEN_WIDTH_UM,
        mask_thickness=MASK_THICKNESS_UM, substrate_top=SUBSTRATE_TOP_UM,
        etched_depth=ETCHED_DEPTH_UM, neck_depth=NECK_DEPTH_UM,
        headspace=HEADSPACE_UM, mesh_length_unit_m=1e-6,
        substrate_material_id=1, mask_material_id=2):
    """Trench whose mask wall is a STRAIGHT taper of prescribed angle.

    The Gaussian constriction of :func:`make_necked_trench_geometry_3d` pins the
    wall vertical at its own minimum -- any smooth minimum has zero slope, so a
    sweep over aperture cannot vary the wall angle independently.  Here the wall
    is a straight wedge of half-angle ``wall_angle_deg`` from vertical, and the
    apex aperture is chosen so that the aperture at the band midpoint equals
    ``aperture_nm`` for every angle.  Angle is then a controlled variable at
    fixed local aperture.
    """
    tan_alpha = float(np.tan(np.deg2rad(float(wall_angle_deg))))
    apex_half = 0.5 * float(aperture_nm) * 1e-3 - ANGLE_BAND_MID_UM * tan_alpha
    if apex_half <= 0.5 * dx:
        raise ValueError("wall angle too steep for the reference aperture")
    domain_height = substrate_top + mask_thickness + headspace
    shape = tuple(max(3, int(round(length / dx)) + 1)
                  for length in (cell_width, cell_length, domain_height))
    x, y, z = (np.arange(size) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    radius = np.abs(X - 0.5 * cell_width)
    floor = substrate_top - etched_depth
    mask_top = substrate_top + mask_thickness
    z_apex = mask_top - neck_depth

    base = floor - Z
    substrate_wall_slab = np.minimum(Z - floor, substrate_top - Z)
    substrate_wall = np.minimum(substrate_wall_slab, radius - 0.5 * opening_width)
    substrate_levelset = np.maximum(base, substrate_wall)

    half_width = np.minimum(apex_half + np.abs(Z - z_apex) * tan_alpha,
                            0.5 * opening_width)
    mask_slab = np.minimum(Z - substrate_top, mask_top - Z)
    mask_levelset = np.minimum(mask_slab, radius - half_width)

    substrate_phi = reinit_narrow(substrate_levelset, dx, domain_height + cell_width)
    mask_phi = reinit_narrow(mask_levelset, dx, domain_height + cell_width)
    phi = reinit_narrow(np.maximum(substrate_phi, mask_phi), dx,
                        domain_height + cell_width)

    substrate_solid = (Z < substrate_top) & ~(
        (Z > floor) & (radius < 0.5 * opening_width))
    mask_solid = (Z >= substrate_top) & (Z < mask_top) & (radius >= half_width)
    material = np.zeros(shape, dtype=int)
    material[substrate_solid] = int(substrate_material_id)
    material[mask_solid] = int(mask_material_id)
    unlabeled_solid = (phi > 0.0) & (material == 0)
    substrate_owner = substrate_levelset >= mask_levelset
    material[unlabeled_solid] = np.where(
        substrate_owner[unlabeled_solid],
        int(substrate_material_id), int(mask_material_id))
    return FeatureGeometry3D(
        phi, material, dx, mesh_length_unit_m,
        material_levelsets={
            int(substrate_material_id): substrate_phi,
            int(mask_material_id): mask_phi,
        })


def measure_aperture_profile(geometry, *, substrate_top=SUBSTRATE_TOP_UM,
                             mask_thickness=MASK_THICKNESS_UM):
    """Aperture (nm) versus depth into the mask, from the nodal level set."""
    x, _, z = geometry.coordinate_arrays
    row = geometry.phi.shape[1] // 2
    center = int(np.argmin(np.abs(x - 0.5 * (x[0] + x[-1]))))
    unit_nm = geometry.mesh_length_unit_m * 1e9
    mask_top = substrate_top + mask_thickness
    depths, apertures = [], []
    for k, z_value in enumerate(z):
        if not (substrate_top + 0.25 * geometry.dx
                <= z_value <= mask_top - 0.25 * geometry.dx):
            continue
        line = geometry.phi[:, row, k]
        if line[center] > 0.0:
            depths.append((mask_top - z_value) * unit_nm)
            apertures.append(0.0)
            continue
        left = center
        while left > 0 and line[left] <= 0.0:
            left -= 1
        right = center
        while right < line.size - 1 and line[right] <= 0.0:
            right += 1
        # Linear interface crossings on each side.
        def cross(i, j):
            fi, fj = line[i], line[j]
            if fi == fj:
                return x[i]
            return x[i] + (x[j] - x[i]) * (0.0 - fi) / (fj - fi)
        x_left = cross(left, left + 1) if line[left] > 0.0 else x[0]
        x_right = cross(right, right - 1) if line[right] > 0.0 else x[-1]
        depths.append((mask_top - z_value) * unit_nm)
        apertures.append((x_right - x_left) * unit_nm)
    return np.asarray(depths), np.asarray(apertures)


def _boundary_and_mechanism(geometry, *, ion_azimuthal_order=16):
    realized = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    source_z = float(realized[2])
    boundary = build_krueger_2024_development_boundary(
        DATA, n_transverse_neutral=5, n_normal_neutral=8,
        reference_plane_m=source_z * geometry.mesh_length_unit_m,
        neutral_direction_polar_order=12, neutral_direction_azimuthal_order=24,
        ion_energy_bin_eV=None, ion_angle_bin_deg=None,
        ion_azimuthal_closure="axisymmetric_uniform",
        ion_azimuthal_order=int(ion_azimuthal_order))
    mechanism = build_krueger_2024_material_router_3d(
        surface_model="mixed_layer", mixed_layer_volatilization_yield=1.0)
    role = {
        species.name: ("energetic_bombardment" if species.charge_number != 0
                       else "neutral_reactant")
        for species in boundary.species}
    return boundary, mechanism, role, realized, source_z


def gather_transport(geometry, *, transport_device="cpu", ion_azimuthal_order=16,
                     face_quadrature_points=3):
    """One frozen-geometry transport gather (duration 0: no motion, no topology risk)."""
    boundary, mechanism, role, realized, source_z = _boundary_and_mechanism(
        geometry, ion_azimuthal_order=ion_azimuthal_order)
    started = perf_counter()
    result = advance_feature_step_3d(
        geometry, boundary, role, mechanism, etchable_material_ids=(1, 2),
        duration_s=0.0, source_bounds=(0.0, realized[0], 0.0, realized[1]),
        source_z=source_z, n_position=1, seed=1, cfl_number=0.25,
        reinitialize=True, reinitialization_method="cr2",
        profile_periodic_lateral=True,
        # The feature cell is periodic in x and y, so the ballistic first-hit
        # gather must be too.  Without this the gather keeps only the rays whose
        # back-projection lands inside the 0.13 x 0.02 um source rectangle,
        # which discards broad-angle thermal neutrals ~36x harder than the
        # near-vertical ion beam and silently distorts every removal/deposition
        # ratio (RESULTS_LIP_DEPOSITION_AUDIT_2026-08-04.md).  The production
        # pilot never had this defect: --radiosity-backend
        # deterministic_extruded_2d sets periodic neutral transport itself.
        ballistic_periodic_lateral=True,
        transport_device=str(transport_device),
        ballistic_transport="face_gather", grazing_ion_reflection={},
        ballistic_face_quadrature_points=int(face_quadrature_points),
        topology_change_policy="continue_gas_cavity",
        surface_state_remap_backend="common_refinement")
    return result, mechanism, role, perf_counter() - started


def relax_chemistry(result, mechanism, role, *, relax_s=0.4, rounds=5,
                    tolerance=5e-3):
    """Integrate the surface mechanism on the frozen flux field to steady state.

    Returns the final ``MaterialSurfaceStepResult3D`` plus the convergence trace
    of the maximum absolute net velocity change between rounds.
    """
    active = np.asarray(result.active_face_index, dtype=int)
    material = np.asarray(result.face_material_id, dtype=int)[active]
    flux = _select_surface_fluxes(
        result.transport.surface_fluxes, active,
        len(result.face_material_id), role)
    state = mechanism.initial_state_by_material(material)
    previous = None
    trace = []
    step = None
    for _ in range(int(rounds)):
        step = mechanism.advance_by_material(state, flux, float(relax_s), material)
        net = (np.asarray(step.etch_velocity_m_s, dtype=float)
               - np.asarray(step.normal_growth_velocity_m_s, dtype=float))
        if previous is not None:
            scale = max(float(np.max(np.abs(net))), 1e-30)
            change = float(np.max(np.abs(net - previous))) / scale
            trace.append(change)
            if change <= float(tolerance):
                state = step.state
                break
        previous = net
        state = step.state
    return step, flux, active, material, trace


def neck_face_mask(result, geometry, *, substrate_top=SUBSTRATE_TOP_UM,
                   mask_thickness=MASK_THICKNESS_UM, neck_depth=NECK_DEPTH_UM,
                   axial_half_window_um=0.08):
    """Select the mask sidewall faces inside the constriction band.

    Neck faces are (a) mask material, (b) within +/- the axial window of the
    neck plane, and (c) laterally facing (the aperture-setting faces, whose
    normal is dominated by x).
    """
    active = np.asarray(result.active_face_index, dtype=int)
    centroid = np.asarray(result.active_face_centroid, dtype=float)
    material = np.asarray(result.face_material_id, dtype=int)[active]
    mask_top = substrate_top + mask_thickness
    z_neck = mask_top - neck_depth
    z = centroid[:, 2] * geometry.dx if centroid.max() > 10.0 else centroid[:, 2]
    in_band = np.abs(z - z_neck) <= float(axial_half_window_um)
    return (material == 2) & in_band, z


def face_incidence(result, *, name_filter=("ions",)):
    """Flux-weighted mean incidence cosine and total energetic flux per active face."""
    n_active = np.asarray(result.active_face_index, dtype=int).size
    weight = np.zeros(n_active)
    weighted_cos = np.zeros(n_active)
    energy = np.zeros(n_active)
    active = np.asarray(result.active_face_index, dtype=int)
    old_to_new = np.full(len(result.face_material_id), -1, dtype=int)
    old_to_new[active] = np.arange(active.size)
    for population in result.transport.surface_fluxes.energetic_fluxes:
        base = population.name.rsplit(":hot_neutral", 1)[0]
        if base not in name_filter:
            continue
        if isinstance(population, FaceResolvedEnergeticFlux):
            face = old_to_new[np.asarray(population.event_face, dtype=int)]
            keep = face >= 0
            face = face[keep]
            flux = np.asarray(population.event_flux_m2_s, dtype=float)[keep]
            cosine = np.asarray(population.event_cosine_incidence, dtype=float)[keep]
            e_value = np.asarray(population.event_energy_eV, dtype=float)[keep]
            np.add.at(weight, face, flux)
            np.add.at(weighted_cos, face, flux * cosine)
            np.add.at(energy, face, flux * e_value)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_cos = np.where(weight > 0.0, weighted_cos / np.maximum(weight, 1e-300), np.nan)
        mean_energy = np.where(weight > 0.0, energy / np.maximum(weight, 1e-300), np.nan)
    return weight, mean_cos, mean_energy


def evaluate_geometry(neck_nm, *, dx, transport_device="cpu", relax_s=0.4,
                      rounds=5):
    """Full frozen-geometry evaluation at a prescribed neck aperture."""
    geometry = make_necked_trench_geometry_3d(
        neck_width_um=float(neck_nm) * 1e-3, dx=float(dx))
    depths, apertures = measure_aperture_profile(geometry)
    interior = apertures[(depths > 30.0) & (depths < 700.0)]
    realised_neck = float(np.min(interior)) if interior.size else float("nan")
    result, mechanism, role, gather_s = gather_transport(
        geometry, transport_device=transport_device)
    step, flux, active, material, trace = relax_chemistry(
        result, mechanism, role, relax_s=relax_s, rounds=rounds)
    etch = np.asarray(step.etch_velocity_m_s, dtype=float)
    growth = np.asarray(step.normal_growth_velocity_m_s, dtype=float)
    net = etch - growth
    is_neck, z_face = neck_face_mask(result, geometry)
    area = np.asarray(result.active_face_area, dtype=float)
    if not np.any(is_neck):
        raise RuntimeError(f"no neck faces selected at {neck_nm} nm")
    weights = area[is_neck]
    net_neck = float(np.average(net[is_neck], weights=weights))
    ion_flux, mean_cos, mean_energy = face_incidence(result)
    return {
        "prescribed_neck_nm": float(neck_nm),
        "realised_neck_nm": realised_neck,
        "dx_um": float(dx),
        "net_neck_velocity_nm_s": net_neck * 1e9,
        "etch_neck_nm_s": float(np.average(etch[is_neck], weights=weights)) * 1e9,
        "growth_neck_nm_s": float(np.average(growth[is_neck], weights=weights)) * 1e9,
        "neck_face_count": int(np.count_nonzero(is_neck)),
        "neck_mean_incidence_cos": float(
            np.nanmean(mean_cos[is_neck])) if np.any(is_neck) else float("nan"),
        "neck_ion_flux_m2_s": float(np.average(ion_flux[is_neck], weights=weights)),
        "gather_wall_s": float(gather_s),
        "relaxation_trace": [float(item) for item in trace],
        "active_faces": int(active.size),
    }, dict(result=result, step=step, geometry=geometry, is_neck=is_neck,
            net=net, etch=etch, growth=growth, area=area, z_face=z_face,
            ion_flux=ion_flux, mean_cos=mean_cos, mean_energy=mean_energy,
            depths=depths, apertures=apertures, material=material)


def sloped_band_mask(result, geometry, *, substrate_top=SUBSTRATE_TOP_UM,
                     mask_thickness=MASK_THICKNESS_UM, neck_depth=NECK_DEPTH_UM):
    """Mask sidewall faces inside the straight-taper measurement band."""
    active = np.asarray(result.active_face_index, dtype=int)
    centroid = np.asarray(result.active_face_centroid, dtype=float)
    material = np.asarray(result.face_material_id, dtype=int)[active]
    z_apex = substrate_top + mask_thickness - neck_depth
    z = centroid[:, 2] * geometry.dx if centroid.max() > 10.0 else centroid[:, 2]
    above = z - z_apex
    in_band = (above >= ANGLE_BAND_LO_UM) & (above <= ANGLE_BAND_HI_UM)
    return (material == 2) & in_band, z


def evaluate_wall_angle(wall_angle_deg, *, aperture_nm, dx,
                        transport_device="cpu", relax_s=0.4, rounds=5):
    """Frozen-geometry evaluation at a prescribed WALL ANGLE, fixed aperture."""
    geometry = make_sloped_wall_geometry_3d(
        wall_angle_deg=float(wall_angle_deg), aperture_nm=float(aperture_nm),
        dx=float(dx))
    depths, apertures = measure_aperture_profile(geometry)
    mask_top = SUBSTRATE_TOP_UM + MASK_THICKNESS_UM
    apex_depth_nm = NECK_DEPTH_UM * 1e3
    band = ((depths >= apex_depth_nm - ANGLE_BAND_HI_UM * 1e3)
            & (depths <= apex_depth_nm - ANGLE_BAND_LO_UM * 1e3))
    band_aperture = float(np.mean(apertures[band])) if band.any() else float("nan")
    if band.sum() >= 2:
        fit = np.polyfit(depths[band], apertures[band] / 2.0, 1)
        realised_angle = float(np.degrees(np.arctan(abs(fit[0]))))
    else:
        realised_angle = float("nan")
    result, mechanism, role, gather_s = gather_transport(
        geometry, transport_device=transport_device)
    step, flux, active, material, trace = relax_chemistry(
        result, mechanism, role, relax_s=relax_s, rounds=rounds)
    etch = np.asarray(step.etch_velocity_m_s, dtype=float)
    growth = np.asarray(step.normal_growth_velocity_m_s, dtype=float)
    net = etch - growth
    is_band, z_face = sloped_band_mask(result, geometry)
    area = np.asarray(result.active_face_area, dtype=float)
    if not np.any(is_band):
        raise RuntimeError(f"no band faces at alpha={wall_angle_deg} deg")
    weights = area[is_band]
    ion_flux, mean_cos, mean_energy = face_incidence(result)
    return {
        "wall_angle_deg": float(wall_angle_deg),
        "realised_wall_angle_deg": realised_angle,
        "reference_aperture_nm": float(aperture_nm),
        "band_aperture_nm": band_aperture,
        "dx_um": float(dx),
        "net_band_velocity_nm_s": float(
            np.average(net[is_band], weights=weights)) * 1e9,
        "etch_band_nm_s": float(
            np.average(etch[is_band], weights=weights)) * 1e9,
        "growth_band_nm_s": float(
            np.average(growth[is_band], weights=weights)) * 1e9,
        "band_face_count": int(np.count_nonzero(is_band)),
        "band_mean_incidence_cos": float(np.nanmean(mean_cos[is_band])),
        "band_ion_flux_m2_s": float(np.average(ion_flux[is_band], weights=weights)),
        "gather_wall_s": float(gather_s),
        "relaxation_trace": [float(item) for item in trace],
    }


def angle_sweep(args, output):
    """Check (c): net velocity versus WALL ANGLE at fixed local aperture."""
    records = []
    for angle in args.angle_deg:
        try:
            record = evaluate_wall_angle(
                angle, aperture_nm=args.angle_aperture_nm, dx=args.dx_um,
                transport_device=args.transport_device, relax_s=args.relax_s,
                rounds=args.rounds)
        except ValueError as error:
            print(f"alpha {angle:5.2f} deg: SKIPPED ({error})", flush=True)
            continue
        records.append(record)
        print(f"alpha {angle:5.2f} deg (realised {record['realised_wall_angle_deg']:5.2f}, "
              f"aperture {record['band_aperture_nm']:6.2f} nm) net "
              f"{record['net_band_velocity_nm_s']:+10.5f} nm/s "
              f"(etch {record['etch_band_nm_s']:.5f} / growth "
              f"{record['growth_band_nm_s']:.5f}) cos {record['band_mean_incidence_cos']:.4f} "
              f"[{record['gather_wall_s']:.0f}s]", flush=True)
    zero_cross = None
    ordered = sorted(records, key=lambda item: item["wall_angle_deg"])
    for lower, upper in zip(ordered, ordered[1:]):
        a, b = lower["net_band_velocity_nm_s"], upper["net_band_velocity_nm_s"]
        if a * b < 0.0:
            span = upper["wall_angle_deg"] - lower["wall_angle_deg"]
            zero_cross = lower["wall_angle_deg"] + span * (0.0 - a) / (b - a)
            break
    payload = {"angle_sweep": ordered, "zero_cross_deg": zero_cross,
               "reference_aperture_nm": args.angle_aperture_nm}
    (output / "angle_sweep.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"\nnet-velocity zero crossing: {zero_cross} deg "
          f"(None = no sign change in the swept range)")
    return payload


def bisect_equilibrium(records):
    """Linear interpolation of the aperture where net neck velocity changes sign."""
    ordered = sorted(records, key=lambda item: item["realised_neck_nm"])
    for lower, upper in zip(ordered, ordered[1:]):
        a, b = lower["net_neck_velocity_nm_s"], upper["net_neck_velocity_nm_s"]
        if a == 0.0:
            return lower["realised_neck_nm"]
        if a * b < 0.0:
            span = upper["realised_neck_nm"] - lower["realised_neck_nm"]
            return lower["realised_neck_nm"] + span * (0.0 - a) / (b - a)
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/curated/mouth_equilibrium_probe")
    parser.add_argument("--dx-um", type=float, default=0.01)
    parser.add_argument("--fine-dx-um", type=float, default=0.005)
    parser.add_argument("--neck-nm", type=float, nargs="*",
                        default=[45.0, 39.0, 33.0, 27.0, 21.0])
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--relax-s", type=float, default=0.4)
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--skip-resolution", action="store_true")
    parser.add_argument("--angle-deg", type=float, nargs="*", default=None,
                        help="check (c): sweep WALL ANGLE at fixed aperture")
    parser.add_argument("--angle-aperture-nm", type=float, default=45.0)
    parser.add_argument("--depth-profile-nm", type=float, default=None,
                        help="run a single geometry and dump net velocity versus "
                             "depth for every mask face (locates the closure)")
    args = parser.parse_args(argv)

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    if args.angle_deg:
        return angle_sweep(args, output)
    if args.depth_profile_nm is not None:
        return depth_profile(args, output)
    payload = {"sweep": [], "budget": None, "resolution": None}
    detail_by_neck = {}

    for neck in args.neck_nm:
        record, detail = evaluate_geometry(
            neck, dx=args.dx_um, transport_device=args.transport_device,
            relax_s=args.relax_s, rounds=args.rounds)
        payload["sweep"].append(record)
        detail_by_neck[float(neck)] = detail
        print(f"neck {neck:5.1f} nm (realised {record['realised_neck_nm']:6.2f}) "
              f"net {record['net_neck_velocity_nm_s']:+9.4f} nm/s "
              f"(etch {record['etch_neck_nm_s']:.4f} / growth "
              f"{record['growth_neck_nm_s']:.4f})  [{record['gather_wall_s']:.0f}s]",
              flush=True)

    payload["equilibrium_aperture_nm"] = bisect_equilibrium(payload["sweep"])

    # Budget decomposition at the experimental neck.
    reference = min(args.neck_nm, key=lambda value: abs(value - EXPERIMENTAL_NECK_NM))
    detail = detail_by_neck[float(reference)]
    payload["budget"] = budget_decomposition(detail, reference)

    if not args.skip_resolution:
        fine_record, fine_detail = evaluate_geometry(
            reference, dx=args.fine_dx_um,
            transport_device=args.transport_device, relax_s=args.relax_s,
            rounds=args.rounds)
        coarse = next(item for item in payload["sweep"]
                      if item["prescribed_neck_nm"] == float(reference))
        payload["resolution"] = {
            "coarse": coarse, "fine": fine_record,
            "net_velocity_ratio": (
                fine_record["net_neck_velocity_nm_s"]
                / coarse["net_neck_velocity_nm_s"]
                if coarse["net_neck_velocity_nm_s"] != 0.0 else None),
        }
        print(f"resolution dx={args.fine_dx_um}: net "
              f"{fine_record['net_neck_velocity_nm_s']:+9.4f} nm/s "
              f"vs {coarse['net_neck_velocity_nm_s']:+9.4f} at dx={args.dx_um}",
              flush=True)

    (output / "probe.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_plots(output, payload, detail_by_neck, reference)
    print(f"equilibrium aperture: {payload['equilibrium_aperture_nm']}")
    return payload


def depth_profile(args, output):
    """Net normal velocity versus depth for every mask face at one aperture.

    The sweep shows the neck band never reverses sign; this locates where the
    closure actually happens, which the single ``mask_opening`` metric hides.
    """
    neck = float(args.depth_profile_nm)
    record, detail = evaluate_geometry(
        neck, dx=args.dx_um, transport_device=args.transport_device,
        relax_s=args.relax_s, rounds=args.rounds)
    material = detail["material"]
    z = np.asarray(detail["z_face"], dtype=float)
    mask_top = SUBSTRATE_TOP_UM + MASK_THICKNESS_UM
    depth_nm = (mask_top - z) * 1e3
    is_mask = material == 2
    net_nm_s = detail["net"] * 1e9
    growth_nm_s = detail["growth"] * 1e9
    ion = detail["ion_flux"]
    cosine = detail["mean_cos"]
    order = np.argsort(depth_nm[is_mask])
    rows = []
    for index in order:
        selection = np.where(is_mask)[0][index]
        rows.append({
            "depth_nm": float(depth_nm[selection]),
            "net_nm_s": float(net_nm_s[selection]),
            "growth_nm_s": float(growth_nm_s[selection]),
            "ion_flux_m2_s": float(ion[selection]),
            "incidence_cos": float(cosine[selection]),
        })
    # Bin by depth for a readable profile.
    edges = np.arange(0.0, 860.0, 50.0)
    binned = []
    depths = np.asarray([row["depth_nm"] for row in rows])
    values = np.asarray([row["net_nm_s"] for row in rows])
    fluxes = np.asarray([row["ion_flux_m2_s"] for row in rows])
    for low, high in zip(edges, edges[1:]):
        keep = (depths >= low) & (depths < high)
        if not np.any(keep):
            continue
        binned.append({
            "depth_lo_nm": float(low), "depth_hi_nm": float(high),
            "faces": int(np.count_nonzero(keep)),
            "mean_net_nm_s": float(np.mean(values[keep])),
            "min_net_nm_s": float(np.min(values[keep])),
            "mean_ion_flux_m2_s": float(np.mean(fluxes[keep])),
        })
    payload = {"aperture_nm": neck, "record": record, "binned": binned,
               "faces": rows}
    (output / f"depth_profile_{int(neck)}nm.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    for item in binned:
        print(f"depth {item['depth_lo_nm']:5.0f}-{item['depth_hi_nm']:5.0f} nm  "
              f"faces {item['faces']:3d}  mean net {item['mean_net_nm_s']:+.5f} "
              f"min {item['min_net_nm_s']:+.5f} nm/s  "
              f"ion {item['mean_ion_flux_m2_s']:.3e}", flush=True)
    return payload


def budget_decomposition(detail, neck_nm):
    """Per-face lip budget at the reference neck, plus removal versus angle."""
    step = detail["step"]
    is_neck = detail["is_neck"]
    area = detail["area"][is_neck]
    mask_result = step.material_results.get(2)
    if mask_result is None:
        raise RuntimeError("mask mechanism result missing")
    material = detail["material"]
    local = np.where(material[is_neck[material == material]] == 2)[0] if False else None
    # Mask-local index of the selected neck faces.
    mask_faces = np.where(material == 2)[0]
    order = {int(face): index for index, face in enumerate(mask_faces)}
    selected = np.where(is_neck)[0]
    local_index = np.asarray([order[int(face)] for face in selected], dtype=int)

    def local_field(name):
        value = np.asarray(getattr(mask_result, name), dtype=float)
        return np.broadcast_to(value, (mask_faces.size,))[local_index]

    deposited = local_field("deposited_polymer_units_m2")
    removed_polymer = local_field("removed_polymer_units_m2")
    removed_bare = local_field("removed_bare_formula_units_m2")
    removed_complex = local_field("removed_complex_units_m2")
    cosine = detail["mean_cos"][is_neck]
    ion_flux = detail["ion_flux"][is_neck]
    net = detail["net"][is_neck]
    total_removal = removed_polymer + removed_bare + removed_complex
    with np.errstate(invalid="ignore", divide="ignore"):
        removal_per_ion = np.where(ion_flux > 0.0, total_removal / ion_flux, np.nan)
    finite = np.isfinite(cosine) & np.isfinite(removal_per_ion) & (removal_per_ion > 0.0)
    normalised = (removal_per_ion[finite] / np.nanmax(removal_per_ion[finite])
                  if np.any(finite) else np.asarray([]))
    angle_deg = np.degrees(np.arccos(np.clip(cosine[finite], -1.0, 1.0)))
    cosine_reference = cosine[finite] / np.nanmax(cosine[finite]) if np.any(finite) else []
    above = (np.asarray(normalised) > np.asarray(cosine_reference)) if np.any(finite) else []
    band = (angle_deg >= 40.0) & (angle_deg <= 70.0) if np.any(finite) else []
    return {
        "reference_neck_nm": float(neck_nm),
        "neck_face_count": int(np.count_nonzero(is_neck)),
        "area_weighted": {
            "deposited_polymer_units_m2": float(np.average(deposited, weights=area)),
            "removed_polymer_units_m2": float(np.average(removed_polymer, weights=area)),
            "removed_bare_units_m2": float(np.average(removed_bare, weights=area)),
            "removed_complex_units_m2": float(np.average(removed_complex, weights=area)),
            "net_velocity_nm_s": float(np.average(net, weights=area)) * 1e9,
        },
        "removal_over_deposition": (
            float(np.average(removed_polymer + removed_bare + removed_complex,
                             weights=area)
                  / np.average(deposited, weights=area))
            if np.average(deposited, weights=area) > 0.0 else None),
        "angle_deg": [float(value) for value in angle_deg],
        "normalised_removal": [float(value) for value in np.asarray(normalised)],
        "cosine_reference": [float(value) for value in np.asarray(cosine_reference)],
        "fraction_above_cosine": (
            float(np.mean(above)) if np.size(above) else None),
        "fraction_above_cosine_40_70deg": (
            float(np.mean(np.asarray(above)[band])) if np.any(band) else None),
    }


def write_plots(output, payload, detail_by_neck, reference):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - plotting is optional
        return
    sweep = payload["sweep"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    aperture = [item["realised_neck_nm"] for item in sweep]
    velocity = [item["net_neck_velocity_nm_s"] for item in sweep]
    axes[0].axhline(0.0, color="0.6", lw=1.0)
    axes[0].plot(aperture, velocity, "o-", color="#1f77b4")
    equilibrium = payload.get("equilibrium_aperture_nm")
    if equilibrium is not None:
        axes[0].axvline(equilibrium, color="#d62728", ls="--",
                        label=f"equilibrium {equilibrium:.1f} nm")
    axes[0].axvline(EXPERIMENTAL_NECK_NM, color="#2ca02c", ls=":",
                    label=f"SEM neck {EXPERIMENTAL_NECK_NM:.0f} nm")
    axes[0].set_xlabel("prescribed neck aperture (nm)")
    axes[0].set_ylabel("net normal velocity at neck (nm/s)\n>0 opens, <0 closes")
    axes[0].set_title("Frozen-geometry neck balance")
    axes[0].legend(fontsize=8)

    budget = payload["budget"]
    angle = np.asarray(budget["angle_deg"])
    removal = np.asarray(budget["normalised_removal"])
    if angle.size:
        order = np.argsort(angle)
        axes[1].plot(angle[order], removal[order], "o", ms=4, color="#1f77b4",
                     label="petch lip removal / ion")
        grid = np.linspace(0.0, 89.0, 200)
        axes[1].plot(grid, np.cos(np.radians(grid)), "-", color="0.4",
                     label="cosine reference")
        axes[1].axvspan(50.0, 60.0, color="#ffcc99", alpha=0.4,
                        label="You 2023: above cosine to 50-60 deg")
    axes[1].set_xlabel("ion incidence angle (deg)")
    axes[1].set_ylabel("normalised removal per incident ion")
    axes[1].set_title(f"Lip removal vs angle at {reference:.0f} nm neck")
    axes[1].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(output / "neck_balance_and_angle.png", dpi=150)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(5.4, 4.4))
    for neck, detail in sorted(detail_by_neck.items()):
        axis.plot(detail["apertures"], detail["depths"], lw=1.2,
                  label=f"{neck:.0f} nm")
    axis.invert_yaxis()
    axis.set_xlabel("aperture (nm)")
    axis.set_ylabel("depth below mask top (nm)")
    axis.set_title("Probe geometries (digitised Fig. 7 profile)")
    axis.legend(fontsize=8, title="prescribed neck")
    figure.tight_layout()
    figure.savefig(output / "probe_geometries.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main()
