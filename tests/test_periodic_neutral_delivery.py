"""The periodic ballistic gather must deliver the full thermal neutral flux.

`RESULTS_LIP_DEPOSITION_AUDIT_2026-08-04.md`: the mouth-equilibrium probe ran
`advance_feature_step_3d` with `profile_periodic_lateral=True` but without
`ballistic_periodic_lateral=True`, so the first-hit gather kept only the rays
whose back-projection landed inside the finite source rectangle.  That cut is
strongly angle selective -- it discards the broad cosine-law thermal neutrals far
harder than the near-vertical ion beam -- so every removal/deposition ratio the
probe reported was distorted without any refusal firing.

These gates pin the contract from both sides: under periodic lateral transport an
exposed face receives exactly the source flux, and the non-periodic finite-source
path is confirmed to be the (legitimate, but very different) alternative so the
two can never be silently swapped again.
"""

from pathlib import Path

import numpy as np
import pytest

from petch.boundary_transport_3d import gather_boundary_state_ballistic_3d
from petch.reactor_boundary import build_krueger_2024_development_boundary

KRUEGER_2024_DATA = (
    Path(__file__).parents[1] / "data" / "experimental" / "krueger_2024")


def _boundary():
    return build_krueger_2024_development_boundary(
        KRUEGER_2024_DATA, n_transverse_neutral=5, n_normal_neutral=8,
        reference_plane_m=1.0e-6,
        neutral_direction_polar_order=12, neutral_direction_azimuthal_order=24,
        ion_azimuthal_closure="axisymmetric_uniform", ion_azimuthal_order=8)


def _flat_plane(size=1.0):
    verts = np.array([[0.0, 0.0, 0.0], [size, 0.0, 0.0],
                      [size, size, 0.0], [0.0, size, 0.0]], dtype=float)
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    tri = verts[faces]
    areas = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    return verts, faces, areas, tri.mean(axis=1)


def _gather(periodic, size=1.0):
    boundary = _boundary()
    verts, faces, areas, centroids = _flat_plane(size)
    role = {species.name: ("energetic_bombardment" if species.charge_number
                           else "neutral_reactant")
            for species in boundary.species}
    result = gather_boundary_state_ballistic_3d(
        boundary, role, verts, faces, areas, centroids,
        np.tile([0.0, 0.0, 1.0], (faces.shape[0], 1)),
        source_bounds=(0.0, size, 0.0, size), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3, device="cpu",
        periodic_lateral=periodic)
    return boundary, result


def test_periodic_gather_delivers_full_neutral_flux():
    """An exposed face under periodic transport receives the source flux."""
    boundary, result = _gather(periodic=True)
    delivered = result.surface_fluxes.neutral_flux_m2_s
    for species in boundary.species:
        if species.charge_number:
            continue
        ratio = np.asarray(delivered[species.name]) / species.flux_m2_s
        assert np.allclose(ratio, 1.0, rtol=1e-6), (species.name, ratio)


def test_periodic_delivery_preserves_species_ratios():
    """Oxygen and the depositors share one angular law, so their delivered
    ratio must equal their source ratio face by face -- the invariant the
    geometry-free O-share gate (RESULTS_O_CHANNEL) rests on."""
    boundary, result = _gather(periodic=True)
    delivered = result.surface_fluxes.neutral_flux_m2_s
    source = {s.name: s.flux_m2_s for s in boundary.species}
    oxygen = np.asarray(delivered["O"]) / source["O"]
    for name in ("CF", "CF2", "CF3", "C2F3"):
        assert np.allclose(np.asarray(delivered[name]) / source[name], oxygen,
                           rtol=1e-9)


def test_nonperiodic_finite_source_truncates_thermal_neutrals():
    """The non-periodic path keeps only rays back-projecting into the source
    rectangle.  That is correct for a genuinely finite source and catastrophic
    when used for a periodic feature cell: it removes most of the thermal
    neutral flux while leaving the near-vertical ion beam nearly intact."""
    boundary, result = _gather(periodic=False)
    delivered = result.surface_fluxes.neutral_flux_m2_s
    neutral_ratio = float(np.mean(
        np.asarray(delivered["CF2"]) / next(
            s.flux_m2_s for s in boundary.species if s.name == "CF2")))
    assert neutral_ratio < 0.5

    ion = next(p for p in result.surface_fluxes.energetic_fluxes
               if p.name == "ions")
    ion_source = next(s.flux_m2_s for s in boundary.species
                      if s.charge_number != 0)
    _, _, areas, _ = _flat_plane()
    ion_landed = float(
        (np.asarray(ion.event_flux_m2_s)
         * areas[np.asarray(ion.event_face)]).sum() / areas.sum())
    ion_ratio = ion_landed / ion_source
    # The angular selectivity is the whole point: the near-vertical beam keeps
    # far more of its flux than the cosine-law neutrals do.
    assert ion_ratio > 3.0 * neutral_ratio


def test_periodic_flux_balance_closes():
    """Opaque periodic cell: every downward source ray lands."""
    boundary, result = _gather(periodic=True)
    _, _, areas, _ = _flat_plane()
    delivered = result.surface_fluxes.neutral_flux_m2_s
    for species in boundary.species:
        if species.charge_number:
            continue
        landed = float((np.asarray(delivered[species.name]) * areas).sum())
        assert landed == pytest.approx(species.flux_m2_s * areas.sum(), rel=1e-6)
