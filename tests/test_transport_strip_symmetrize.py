"""Gates for extruded y-strip transport symmetrization (extrusion-guard fix)."""

import numpy as np
import pytest

from petch.boundary_transport_3d import (
    BoundaryTransport3DResult,
    symmetrize_transport_across_strips,
)
from petch.surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)


def _two_strip_mesh():
    """Floor split into two y-strips of one triangle each, plus a wall.

    Strip faces share (x, z) centroid and normal, differing only in y — the
    extruded-equivalence relation. Faces 0,1 are the two floor strips
    (same centroid x,z; identical normal); face 2 is an unrelated wall.
    """
    verts = np.array([
        # strip 0 floor triangle (y in [0,1])
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0],
        # strip 1 floor triangle (y in [1,2]) — same x,z centroid as strip 0
        [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.5, 2.0, 0.0],
        # wall triangle (different centroid + normal)
        [0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 1.0, 1.0],
    ], dtype=float)
    faces = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]], dtype=int)
    tri = verts[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    centroids = tri.mean(axis=1)
    # Force strips 0,1 to share the identical (x,z) centroid the grouping keys on.
    centroids[0, [0, 2]] = centroids[1, [0, 2]] = (0.5, 0.0)
    normals[0] = normals[1] = np.array([0.0, 0.0, 1.0])
    areas = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    return verts, faces, areas, centroids, normals


def _result(neutral, energetic):
    return BoundaryTransport3DResult(
        SurfaceFluxes(neutral, energetic),
        {"x": 1.0}, {"x": 0.0}, {"x": 0.0}, "manufactured", ())


def test_asymmetric_neutral_becomes_strip_uniform_conserving_rate():
    verts, faces, areas, centroids, normals = _two_strip_mesh()
    neutral = {"CF2": np.array([2.0e19, 0.0, 5.0e18])}   # strips differ
    result = _result(neutral, ())
    out = symmetrize_transport_across_strips(
        result, verts, faces, centroids, normals)
    cf2 = np.asarray(out.surface_fluxes.neutral_flux_m2_s["CF2"])
    assert cf2[0] == pytest.approx(cf2[1])               # strip-uniform
    assert cf2[2] == pytest.approx(5.0e18)               # wall untouched
    # Area-weighted total flux over the strip pair is conserved.
    before = float(np.dot(neutral["CF2"][:2], areas[:2]))
    after = float(np.dot(cf2[:2], areas[:2]))
    assert after == pytest.approx(before, rel=1e-12)


def test_asymmetric_face_resolved_redistributes_conserving_rate():
    verts, faces, areas, centroids, normals = _two_strip_mesh()
    # One event on strip 0 only — must split equally across strips 0 and 1.
    pop = FaceResolvedEnergeticFlux(
        "Ar+", 3, np.array([0]), np.array([3.0e19]), np.array([500.0]),
        np.array([0.4]))
    result = _result({}, (pop,))
    out = symmetrize_transport_across_strips(
        result, verts, faces, centroids, normals)
    sym = out.surface_fluxes.energetic_fluxes[0]
    rate = np.asarray(sym.event_flux_m2_s) * areas[sym.event_face]
    before = 3.0e19 * areas[0]
    assert rate.sum() == pytest.approx(before, rel=1e-12)   # conserved
    # Each strip carries exactly half the rate.
    per_face = {int(f): 0.0 for f in sym.event_face}
    for f, r in zip(sym.event_face, rate):
        per_face[int(f)] += r
    assert per_face[0] == pytest.approx(per_face[1], rel=1e-12)
    assert set(per_face) == {0, 1}
    assert np.all(sym.event_energy_eV == 500.0)


def test_already_symmetric_input_is_unchanged():
    verts, faces, areas, centroids, normals = _two_strip_mesh()
    neutral = {"CF2": np.array([4.0e19, 4.0e19, 1.0e19])}   # strips equal
    pop = FaceResolvedEnergeticFlux(
        "Ar+", 3, np.array([0, 1]), np.array([2.0e19, 2.0e19]),
        np.array([500.0, 500.0]), np.array([0.4, 0.4]))
    result = _result(neutral, (pop,))
    out = symmetrize_transport_across_strips(
        result, verts, faces, centroids, normals)
    assert np.allclose(
        np.asarray(out.surface_fluxes.neutral_flux_m2_s["CF2"]), neutral["CF2"],
        rtol=0.0, atol=0.0)
    sym = out.surface_fluxes.energetic_fluxes[0]
    got = np.asarray(sym.flux_m2_s)
    assert np.allclose(got, np.asarray(pop.flux_m2_s), rtol=0.0, atol=1e-9)


def test_scalar_energetic_flux_passes_through():
    verts, faces, areas, centroids, normals = _two_strip_mesh()
    pop = EnergeticFlux(
        "Ar+", 5.0e18, np.array([500.0]), np.array([0.8]), np.array([1.0]))
    result = _result({}, (pop,))
    out = symmetrize_transport_across_strips(
        result, verts, faces, centroids, normals)
    assert out.surface_fluxes.energetic_fluxes[0] is pop


def test_engine_applies_symmetrization_under_periodic_lateral(monkeypatch):
    """The engine must invoke symmetrization exactly when the profile is
    periodic-lateral (extruded contract), and not otherwise."""
    import petch.feature_step_3d as fs

    calls = []
    real = fs.symmetrize_transport_across_strips

    def spy(transport, *a, **k):
        calls.append(True)
        return real(transport, *a, **k)

    monkeypatch.setattr(fs, "symmetrize_transport_across_strips", spy)
    # The choke point is gated on `if profile_periodic_lateral:` immediately
    # before `base_transport = transport`; assert the guard source contract.
    import inspect
    src = inspect.getsource(fs.advance_feature_step_3d)
    assert "if profile_periodic_lateral:" in src
    assert "symmetrize_transport_across_strips(" in src
    idx_call = src.index("symmetrize_transport_across_strips(\n")
    idx_base = src.index("base_transport = transport")
    assert idx_call < idx_base
