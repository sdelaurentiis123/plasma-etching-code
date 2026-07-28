"""Gates for the grazing-ion specular reflection splitter (audit P0.2)."""

import numpy as np
import pytest

from petch.boundary_transport_3d import split_grazing_ion_reflection
from petch.surface_kinetics import FaceResolvedEnergeticFlux


def _box_trench_mesh():
    """Two facing vertical walls plus a floor, 1x1 lateral, walls at x=0.2/0.8."""
    quads = [
        # left wall (gas normal +x), from z=0 to z=1 at x=0.2
        ([0.2, 0.0, 0.0], [0.2, 1.0, 0.0], [0.2, 1.0, 1.0], [0.2, 0.0, 1.0], [1.0, 0.0, 0.0]),
        # right wall (gas normal -x) at x=0.8
        ([0.8, 0.0, 0.0], [0.8, 0.0, 1.0], [0.8, 1.0, 1.0], [0.8, 1.0, 0.0], [-1.0, 0.0, 0.0]),
        # floor (gas normal +z) at z=0
        ([0.2, 0.0, 0.0], [0.8, 0.0, 0.0], [0.8, 1.0, 0.0], [0.2, 1.0, 0.0], [0.0, 0.0, 1.0]),
    ]
    verts = []
    faces = []
    normals = []
    for a, b, c, d, n in quads:
        base = len(verts)
        verts.extend([a, b, c, d])
        faces.append([base, base + 1, base + 2])
        faces.append([base, base + 2, base + 3])
        normals.extend([n, n])
    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=int)
    normals = np.asarray(normals, dtype=float)
    tri = verts[faces]
    areas = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    centroids = tri.mean(axis=1)
    return verts, faces, areas, centroids, normals


def _population(face, flux, energy, cosine, direction):
    return FaceResolvedEnergeticFlux(
        "Ar+", 6, np.asarray(face), np.asarray(flux, dtype=float),
        np.asarray(energy, dtype=float), np.asarray(cosine, dtype=float),
        event_incident_direction=np.asarray(direction, dtype=float))


def test_normal_incidence_is_a_bitwise_noop():
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    population = _population([4], [1.0e19], [1000.0], [1.0], [[0.0, 0.0, -1.0]])
    primary, secondary, diag = split_grazing_ion_reflection(
        population, verts, faces, areas, centroids, normals,
        domain_size=(1.0, 1.0, 1.0), periodic_lateral=False)
    assert secondary is None
    assert diag["reflected_rate"] == 0.0
    assert primary is population  # collision weight untouched


def test_grazing_wall_event_reflects_to_facing_wall_with_conservation():
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    # Ion nearly parallel to the left wall: mostly -z, slight -x INTO the wall
    # (cos=0.02: truly grazing, where the expectation form keeps continuation)
    direction = np.array([[-0.02, 0.0, -0.9998]])
    direction /= np.linalg.norm(direction)
    cosine = float(abs(direction[0] @ np.array([1.0, 0.0, 0.0])))
    population = _population([0], [2.0e19], [1500.0], [cosine], direction)
    primary, secondary, diag = split_grazing_ion_reflection(
        population, verts, faces, areas, centroids, normals,
        domain_size=(1.0, 1.0, 1.0), periodic_lateral=False)
    assert secondary is not None
    # ADDITIVE: the collision keeps full weight; the continuing hot neutral
    # carries the reflected measure exactly (spawned + escaped == reflected).
    assert primary is population
    original = float(population.event_flux_m2_s[0] * areas[0])
    cos0 = float(population.event_cosine_incidence[0])
    kress = max((1.0 + 9.3 * (1.0 - cos0 ** 2)) * cos0, 0.0)
    y_react = min(max(0.9 * (1500.0 ** 0.5 - 20.0 ** 0.5)
                      / (500.0 ** 0.5 - 20.0 ** 0.5) * kress, 0.0), 1.0)
    expected_reflected = (0.95 * (1.0 - cos0 ** 3)
                          * (1.0 - y_react) * original)
    landed = float((secondary.event_flux_m2_s
                    * areas[secondary.event_face]).sum())
    assert landed + diag["escaped_rate"] == pytest.approx(
        expected_reflected, rel=1e-12)
    assert diag["reflected_rate"] == pytest.approx(expected_reflected, rel=1e-12)
    # Energy retention is exactly the declared fraction
    assert secondary.event_energy_eV[0] == pytest.approx(0.90 * 1500.0)
    # The specular partner of a wall-grazing downward ray continues downward:
    # it must land on the floor or the facing wall, never back on wall 0/1
    assert int(secondary.event_face[0]) in (2, 3, 4, 5)
    assert secondary.name == "Ar+:hot_neutral"


def test_specular_direction_reverses_normal_component_only():
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    direction = np.array([[-0.08, 0.0, -0.9968]])
    direction /= np.linalg.norm(direction)
    population = _population(
        [0], [1.0e19], [120.0],
        [float(abs(direction[0] @ np.array([1.0, 0.0, 0.0])))], direction)
    _, secondary, _ = split_grazing_ion_reflection(
        population, verts, faces, areas, centroids, normals,
        domain_size=(1.0, 1.0, 1.0), periodic_lateral=False)
    reflected = np.asarray(secondary.event_incident_direction[0])
    expected = direction[0] - 2.0 * (direction[0] @ np.array([1.0, 0.0, 0.0])) \
        * np.array([1.0, 0.0, 0.0])
    assert np.allclose(reflected, expected / np.linalg.norm(expected), atol=1e-6)


def test_reflected_weight_follows_declared_angular_law():
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    for cosine in (0.05, 0.3, 0.7):
        direction = np.array([[np.sqrt(1 - 0.0), 0.0, 0.0]])  # placeholder dir
        d = np.array([[cosine, 0.0, -np.sqrt(max(1 - cosine**2, 0.0))]])
        population = _population([0], [1.0e19], [1000.0], [cosine], d)
        primary, secondary, diag = split_grazing_ion_reflection(
            population, verts, faces, areas, centroids, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=False)
        kress = max((1.0 + 9.3 * (1.0 - cosine ** 2)) * cosine, 0.0)
        y_react = min(max(0.9 * (1000.0 ** 0.5 - 20.0 ** 0.5)
                          / (500.0 ** 0.5 - 20.0 ** 0.5) * kress, 0.0), 1.0)
        expected_weight = 0.95 * (1.0 - cosine ** 3) * (1.0 - y_react)
        original = float(population.event_flux_m2_s[0] * areas[0])
        reflected = diag["reflected_rate"]
        assert reflected == pytest.approx(expected_weight * original, rel=1e-12)


def test_empty_options_dict_activates_reflection_in_gather():
    """Regression: the pilot passes {} for default options; {} is falsy, and
    an `if options:` gate silently disabled reflection for a whole campaign
    run. Empty options MUST mean 'on with defaults'."""
    from tests.test_boundary_transport_3d import _boundary, _flat_unit_plane
    from petch.boundary_transport_3d import gather_boundary_state_ballistic_3d

    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    result = gather_boundary_state_ballistic_3d(
        _boundary(), {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"},
        verts, faces, areas, centroids, np.tile([0.0, 0.0, 1.0], (2, 1)),
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3, device="cpu",
        grazing_ion_reflection={})
    # Normal incidence on a flat plane: reflection runs but reflects ~nothing;
    # the diagnostic entry proves the code path executed.
    assert "Ar+:hot_neutral" in result.hit_probability
    diag = result.hit_probability["Ar+:hot_neutral"]
    assert diag["reflected_rate"] >= 0.0


def test_role_selection_keeps_hot_neutrals():
    """Regression: 'ions:hot_neutral' is not a declared boundary species, so
    the role filter silently dropped the entire reflected population (ml10
    ran bitwise-identical to no-reflection despite 2e6 secondary events per
    step in the receipts). Hot neutrals inherit the parent ion's role."""
    from petch.feature_step_3d import _select_surface_fluxes
    from petch.surface_kinetics import (
        FaceResolvedEnergeticFlux, SurfaceFluxes as EngineFluxes)

    ions = FaceResolvedEnergeticFlux(
        "ions", 4, np.array([0, 2]), np.array([1e19, 2e19]),
        np.array([1000.0, 1200.0]), np.array([1.0, 0.5]))
    hot = FaceResolvedEnergeticFlux(
        "ions:hot_neutral", 4, np.array([1]), np.array([5e18]),
        np.array([900.0]), np.array([0.3]))
    fluxes = EngineFluxes({"CF2": np.full(4, 1e19)}, (ions, hot))
    selected = _select_surface_fluxes(
        fluxes, np.arange(4), 4,
        species_role={"ions": "energetic_bombardment",
                      "CF2": "neutral_reactant"})
    names = {p.name for p in selected.energetic_fluxes}
    assert names == {"ions", "ions:hot_neutral"}


def test_reacting_fraction_suppresses_continuation_at_high_yield():
    """Expectation form: where the published film law reacts strongly the
    continuation collapses (Krueger MC either-or, in expectation). A grazing
    low-energy event continues near-fully; a mid-angle high-energy one
    (kress peak, yield ~1) continues barely."""
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    def continued(cosine, energy):
        d = np.array([[-cosine, 0.0, -np.sqrt(max(1 - cosine ** 2, 1e-9))]])
        d /= np.linalg.norm(d)
        pop = _population([0], [1.0e19], [energy], [cosine], d)
        _, _, diag = split_grazing_ion_reflection(
            pop, verts, faces, areas, centroids, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=False)
        return diag["reflected_rate"] / (1.0e19 * areas[0])
    assert continued(0.05, 100.0) > 0.5      # grazing, low E: funnel survives
    assert continued(0.6, 1500.0) < 0.05     # kress peak, high E: consumed
