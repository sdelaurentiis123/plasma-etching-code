"""Gates for the exact MCFPM reflection cascade (Eq. 2.34 + leftover rule)."""

import numpy as np
import pytest

from petch.boundary_transport_3d import split_grazing_ion_reflection
from petch.surface_kinetics import FaceResolvedEnergeticFlux


def _box_trench_mesh():
    quads = [
        ([0.2, 0.0, 0.0], [0.2, 1.0, 0.0], [0.2, 1.0, 1.0], [0.2, 0.0, 1.0], [1.0, 0.0, 0.0]),
        ([0.8, 0.0, 0.0], [0.8, 0.0, 1.0], [0.8, 1.0, 1.0], [0.8, 1.0, 0.0], [-1.0, 0.0, 0.0]),
        ([0.2, 0.0, 0.0], [0.8, 0.0, 0.0], [0.8, 1.0, 0.0], [0.2, 1.0, 0.0], [0.0, 0.0, 1.0]),
    ]
    verts, faces, normals = [], [], []
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


def _split(pop, **kw):
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    return split_grazing_ion_reflection(
        pop, verts, faces, areas, centroids, normals,
        domain_size=(1.0, 1.0, 1.0), periodic_lateral=False, **kw)


def test_normal_incidence_fully_consumed():
    """cos=1: kress=1, react=0.9 -> weight 0.1 continues but arrives with
    retained energy 0 (theta<70) -> thermalized, no secondary population."""
    pop = _population([4], [1.0e19], [1000.0], [1.0], [[0.0, 0.0, -1.0]])
    primary, secondary, diag = _split(pop)
    assert primary is pop
    assert secondary is None
    assert diag["thermalized_rate"] > 0.0


def test_grazing_fast_ray_cascades_specular_full_energy():
    """cos=0.02 @1500 eV: react=0.9*kress(0.02)~0.17 -> ~0.83 continues at
    FULL energy (specular branch), multi-bounce down the trench."""
    d = np.array([[-0.02, 0.0, -0.9998]])
    d /= np.linalg.norm(d)
    pop = _population([0], [2.0e19], [1500.0], [0.02], d)
    primary, secondary, diag = _split(pop)
    assert secondary is not None
    assert np.all(secondary.event_energy_eV == 1500.0)  # full retention
    assert diag["bounce_generations"] >= 1
    first_weight = 1.0 - min(0.9 * (1 + 9.3 * (1 - 0.02**2)) * 0.02, 1.0)
    total_first = float((secondary.event_flux_m2_s
                         * np.asarray(_box_trench_mesh()[2])[secondary.event_face])[
                             :1].sum())
    # First-generation landed rate equals leftover weight exactly (no escape
    # for a downward ray in a closed trench).
    original = 2.0e19 * _box_trench_mesh()[2][0]
    assert total_first == pytest.approx(first_weight * original, rel=1e-9)


def test_leftover_rule_zero_when_reactions_saturate():
    """kress peak (cos~0.6): react capped at 1 -> nothing continues."""
    d = np.array([[-0.6, 0.0, -0.8]])
    d /= np.linalg.norm(d)
    pop = _population([0], [1.0e19], [1500.0], [0.6], d)
    _, secondary, diag = _split(pop)
    assert secondary is None
    assert diag["reflected_rate"] == 0.0


def test_slow_particle_thermalizes():
    """E<10 eV: below diffusive cutoff, dropped from the cascade."""
    d = np.array([[-0.02, 0.0, -0.9998]])
    d /= np.linalg.norm(d)
    pop = _population([0], [1.0e19], [8.0], [0.02], d)
    _, secondary, diag = _split(pop)
    assert secondary is None
    assert diag["thermalized_rate"] > 0.0


def test_role_selection_keeps_hot_neutrals():
    from petch.feature_step_3d import _select_surface_fluxes
    from petch.surface_kinetics import SurfaceFluxes as EngineFluxes

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


def test_empty_options_dict_activates_reflection_in_gather():
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
    assert "Ar+:hot_neutral" in result.hit_probability
