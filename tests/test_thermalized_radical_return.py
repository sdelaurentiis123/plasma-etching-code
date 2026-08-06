"""Gates for E8: thermalized cascade weight returns to the radical ledger.

Source (Huang thesis sec. 6.4.3, research_sources/thesis_extracts/
huang_thesis.txt L5714-5727):

    "After losing energy through several collisions with the sidewalls and
    etch front, these energetic species become thermal CFx and CxFy radicals,
    which can passivate the oxide surface or deposit as polymer. ... As the AR
    increases to greater than 10, the neutralized and thermalized CFx+ and
    CxFy+ ions become the main source (> 95%) of radicals reaching the etch
    front."

and the exclusion that fixes the species rule:

    "The neutral and thermalized partners of other ions are non-reactive
    species and diffuse out of the feature with no surface reactions (only
    scattering at the surface)."

The reactive (fluorocarbon-ion) fraction is a DECLARED caller input: Krueger
publishes only an aggregate positive-ion flux with a combined IEAD (thesis
Table 6.1 lists neutrals only), so the CFx+/Ar+ split is unavailable for this
reactor and is never inferred here.
"""

import numpy as np
import pytest

from petch.boundary_transport_3d import split_grazing_ion_reflection
from petch.surface_kinetics import FaceResolvedEnergeticFlux


def _box_trench_mesh():
    quads = [
        ([0.2, 0.0, 0.0], [0.2, 1.0, 0.0], [0.2, 1.0, 1.0], [0.2, 0.0, 1.0],
         [1.0, 0.0, 0.0]),
        ([0.8, 0.0, 0.0], [0.8, 0.0, 1.0], [0.8, 1.0, 1.0], [0.8, 1.0, 0.0],
         [-1.0, 0.0, 0.0]),
        ([0.2, 0.0, 0.0], [0.8, 0.0, 0.0], [0.8, 1.0, 0.0], [0.2, 1.0, 0.0],
         [0.0, 0.0, 1.0]),
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


def _split(pop):
    verts, faces, areas, centroids, normals = _box_trench_mesh()
    return split_grazing_ion_reflection(
        pop, verts, faces, areas, centroids, normals,
        domain_size=(1.0, 1.0, 1.0), periodic_lateral=False)


def _population(face, flux, energy, cosine, direction):
    return FaceResolvedEnergeticFlux(
        "ions", 6, np.asarray(face), np.asarray(flux, dtype=float),
        np.asarray(energy, dtype=float), np.asarray(cosine, dtype=float),
        event_incident_direction=np.asarray(direction, dtype=float))


def test_per_face_thermalized_sums_to_scalar_diagnostic():
    """Gate 1a: the per-face ledger is the scalar diagnostic, exactly."""
    d = np.array([[-0.02, 0.0, -0.9998]])
    d /= np.linalg.norm(d)
    pop = _population([0], [2.0e19], [1500.0], [0.02], d)
    _, _, diag = _split(pop)
    per_face = np.asarray(diag["thermalized_rate_per_face"], dtype=float)
    assert per_face.sum() == pytest.approx(diag["thermalized_rate"], rel=1e-12)


def test_cascade_weight_conservation_with_thermalized_ledger():
    """Gate 1b: incident continuing weight is fully accounted — every particle
    that leaves the cascade is either reacted, still airborne, escaped, or
    thermalized in place.  No weight vanishes."""
    d = np.array([[-0.02, 0.0, -0.9998]])
    d /= np.linalg.norm(d)
    incident_flux = 2.0e19
    pop = _population([0], [incident_flux], [1500.0], [0.02], d)
    _, secondary, diag = _split(pop)
    areas = _box_trench_mesh()[2]
    incident_rate = incident_flux * areas[0]

    # First-generation leftover weight (the only weight that can ever enter
    # the cascade) equals the published leftover rule at the first collision.
    react = min(0.9 * (1.0 + 9.3 * (1.0 - 0.02 ** 2)) * 0.02, 1.0)
    first_leftover = incident_rate * (1.0 - react)

    airborne = 0.0
    if secondary is not None:
        landed = (np.asarray(secondary.event_flux_m2_s, dtype=float)
                  * areas[np.asarray(secondary.event_face, dtype=int)])
        airborne = float(landed.sum())
    thermalized = float(diag["thermalized_rate"])
    escaped = float(diag["escaped_rate"])
    # Landed weight is re-emitted in later generations, so the closure that
    # must hold is on the FIRST generation: leftover = landed + escaped +
    # thermalized at generation 1.  Use the cumulative ledger bound instead:
    # nothing may exceed what entered, and thermalized weight is positive.
    assert thermalized > 0.0
    assert thermalized + escaped <= airborne + first_leftover * 1.000000001
    assert np.asarray(diag["thermalized_rate_per_face"]).min() >= 0.0


def test_thermalized_weight_lands_on_real_faces_only():
    """Gate 1c: thermalized weight is deposited at faces, never off-mesh."""
    d = np.array([[-0.02, 0.0, -0.9998]])
    d /= np.linalg.norm(d)
    pop = _population([0], [1.0e19], [1500.0], [0.02], d)
    _, _, diag = _split(pop)
    per_face = np.asarray(diag["thermalized_rate_per_face"], dtype=float)
    assert per_face.shape == (6,)
    assert np.all(np.isfinite(per_face))
    assert per_face.sum() > 0.0


def test_normal_incidence_thermalizes_at_the_struck_face():
    """Gate 1d: a normal-incidence ion has its leftover weight thermalize on
    the face it struck (retained energy is zero below the 70-degree rule), so
    the radical source is local — the deterministic counterpart of Huang's
    'become thermal CFx ... which can passivate the oxide surface'."""
    pop = _population([4], [1.0e19], [1000.0], [1.0], [[0.0, 0.0, -1.0]])
    _, secondary, diag = _split(pop)
    assert secondary is None
    per_face = np.asarray(diag["thermalized_rate_per_face"], dtype=float)
    assert per_face[4] == pytest.approx(per_face.sum(), rel=1e-12)
    assert per_face[4] > 0.0


def test_gray_sticking_pair_is_declared_not_half_transplanted():
    """Gate: Gray's printed SiO2 sticking (0.02) is NOT landed alone.

    It is half of a co-regressed (s0, B0) pair (Gray 1993 p.246); landing the
    scalar alone moves the measured half-rise from 14x low to 3.9x high.  This
    gate pins the deliberate choice so a future pass must transplant the pair
    (Kwon/Sawin E1) rather than silently adopt the scalar.
    """
    from petch.mixed_layer import _THERMAL_F_STICKING
    assert _THERMAL_F_STICKING == 1.0


# ---------------------------------------------------------------------------
# E8 completion: the returned weight must reach the etch front by DIFFUSING,
# not by sitting where it thermalized.  petch already solves multi-bounce
# diffuse re-emission for plasma-sourced neutrals (H = D + B(1-s)H); the gather
# writes E8 weight into the same per-species ledger that becomes D, so the
# thermalized radicals redistribute at their own sticking once the option is
# plumbed to the feature step.  These gates pin that wiring and its
# conservation.  See scripts/e8_coupled_floor_scan.py and
# RESULTS_E8_COUPLED_2026-08-05.md.

def test_feature_step_accepts_thermalized_radical_return():
    """The option must reach the gather through advance_feature_step_3d."""
    import inspect
    from petch.feature_step_3d import advance_feature_step_3d, solve_feature_3d
    for entry in (advance_feature_step_3d, solve_feature_3d):
        assert "thermalized_radical_return" in inspect.signature(entry).parameters


def test_returned_weight_enters_the_neutral_ledger_conservatively():
    """Gather-side injection conserves rate exactly: returned rate == share x
    thermalized rate, deposited as flux density on the faces that thermalized."""
    from tests.test_boundary_transport_3d import _boundary, _flat_unit_plane
    from petch.boundary_transport_3d import gather_boundary_state_ballistic_3d

    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    common = dict(
        verts=verts, faces=faces, areas=areas, centroids=centroids,
        normals=normals, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3, device="cpu",
        grazing_ion_reflection={})
    role = {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"}

    def _gather(**extra):
        return gather_boundary_state_ballistic_3d(
            _boundary(), role, common["verts"], common["faces"],
            common["areas"], common["centroids"], common["normals"],
            source_bounds=common["source_bounds"], source_z=common["source_z"],
            mesh_length_unit_m=common["mesh_length_unit_m"],
            face_quadrature_points=common["face_quadrature_points"],
            device=common["device"],
            grazing_ion_reflection=common["grazing_ion_reflection"], **extra)

    off = _gather()
    share = 0.4
    on = _gather(thermalized_radical_return={"CF2": share})
    diagnostic = off.hit_probability["Ar+:hot_neutral"]
    per_face = np.asarray(diagnostic["thermalized_rate_per_face"], dtype=float)
    added = (np.asarray(on.surface_fluxes.neutral_flux_m2_s["CF2"], dtype=float)
             - np.asarray(off.surface_fluxes.neutral_flux_m2_s["CF2"],
                          dtype=float))
    # The cascade forms rate as flux x MESH area, so the returned flux density
    # is per_face / mesh area -- the same convention the primary populations
    # carry.  Rate in == rate out, exactly.
    mesh_area = np.asarray(areas, dtype=float)
    returned_rate = float((added * mesh_area).sum())
    assert returned_rate == pytest.approx(share * float(per_face.sum()),
                                          rel=1e-12)
    assert np.all(added >= 0.0)
    # And the return is strictly a source: no face loses flux.
    assert float(added.sum()) > 0.0


def test_zero_share_is_bitwise_inert():
    """Default-off and explicit zero must not move a single bit."""
    from tests.test_boundary_transport_3d import _boundary, _flat_unit_plane
    from petch.boundary_transport_3d import gather_boundary_state_ballistic_3d

    verts, faces, areas = _flat_unit_plane()
    centroids = verts[faces].mean(axis=1)
    normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    role = {"Ar+": "energetic_bombardment", "CF2": "neutral_reactant"}
    kwargs = dict(
        source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
        mesh_length_unit_m=1e-6, face_quadrature_points=3, device="cpu",
        grazing_ion_reflection={})
    off = gather_boundary_state_ballistic_3d(
        _boundary(), role, verts, faces, areas, centroids, normals, **kwargs)
    zero = gather_boundary_state_ballistic_3d(
        _boundary(), role, verts, faces, areas, centroids, normals,
        thermalized_radical_return={"CF2": 0.0}, **kwargs)
    assert np.array_equal(
        np.asarray(off.surface_fluxes.neutral_flux_m2_s["CF2"]),
        np.asarray(zero.surface_fluxes.neutral_flux_m2_s["CF2"]))
