"""Gates for the axisymmetric hole evolution driver (hole study phase 2).

The driver introduces no new physics: it moves a profile under the transport
that phase 1 characterised frozen, coupled to the mixed-layer chemistry the
trench campaign validated.  These gates hold it to exactly that claim --
the transport it computes must reproduce the gated benchmarks, its ledgers must
close, and its degenerate limit must be the 0-D blanket.
"""

import importlib.util
from math import pi
from pathlib import Path

import numpy as np
import pytest

from petch.axisymmetric_evolution import (
    HoleEvolutionState,
    HoleGeometry,
    advance_hole_step,
    build_hole_enclosure,
    cascade_hole_delivery,
    evolve_hole,
    solve_diffuse_hole_delivery,
)
from petch.axisymmetric_exchange_3d import cylinder_clausing_transmission
from petch.iadf_two_component import kim_2025_reference_iadf
from petch.mixed_layer_mechanism import (
    MixedLayerSurfaceState,
    build_krueger_2024_mixed_layer_mechanisms,
)
from petch.reactor_boundary import build_krueger_2024_development_boundary
from petch.surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "krueger_2024"
REFERENCE_ENERGY_EV = 1500.0
#: Coarse cascade quadrature: these gates test bookkeeping and limits, not the
#: converged delivery numbers (those are phase 1's, regressed below at the
#: production quadrature).
FAST_CASCADE = dict(n_polar=48, n_azimuth=16, n_radial=8)


def _boundary_fluxes():
    boundary = build_krueger_2024_development_boundary(
        DATA, reference_plane_m=2.0e-6)
    mouth = {item.name: float(item.flux_m2_s)
             for item in boundary.species if item.charge_number == 0}
    return mouth, float(boundary.get("ions").flux_m2_s)


def _bare(faces):
    return MixedLayerSurfaceState.bare((faces,))


# --- transport regression ------------------------------------------------------


@pytest.mark.parametrize("aspect", [10.0, 50.0, 100.0])
def test_enclosure_reproduces_gated_clausing_transmission(aspect):
    """Perfectly reflecting walls + absorbing floor IS the Clausing problem.

    The driver assembles its enclosure from the gated band-exchange algebra plus
    disk factors; solved in that limit it must return the benchmark value the
    phase-1 study quotes, not merely something close.
    """
    radius, depth = 0.5, float(aspect)
    edges = np.linspace(0.0, depth, max(24, int(24 * aspect)) + 1)
    enclosure = build_hole_enclosure(radius, depth, edges)
    sticking = np.zeros(enclosure.area.size)
    sticking[-1] = 1.0
    _arrival, absorbed, escaped = solve_diffuse_hole_delivery(
        enclosure, sticking, 1.0)
    reference = cylinder_clausing_transmission(aspect)
    assert absorbed[-1] == pytest.approx(reference, rel=1e-9)
    assert absorbed.sum() + escaped == pytest.approx(1.0, abs=1e-12)


def test_enclosure_geometry_factors_close():
    """Every face's exchange row plus its escape must sum to one."""
    edges = np.linspace(0.0, 20.0, 241)
    enclosure = build_hole_enclosure(0.5, 20.0, edges)
    rows = enclosure.face_to_face.sum(axis=1) + enclosure.face_to_mouth
    assert np.allclose(rows, 1.0, atol=1e-9)
    assert enclosure.mouth_to_face.sum() == pytest.approx(1.0, abs=1e-9)
    assert enclosure.closure_residual < 1e-9


@pytest.mark.parametrize("tail", [0.0, 0.65])
@pytest.mark.parametrize("aspect", [10.0, 50.0])
def test_cascade_reproduces_phase1_bottom_delivery(tail, aspect):
    """The band-resolved cascade must agree with the phase-1 characterisation.

    Same rule, same quadrature, different bookkeeping (arrivals per band rather
    than only the reacting share), so bottom delivery must match to the last
    bit -- any drift means the two implementations have diverged.
    """
    spec = importlib.util.spec_from_file_location(
        "hole_study_phase1", ROOT / "scripts" / "hole_study_phase1.py")
    phase1 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(phase1)
    iadf = kim_2025_reference_iadf(tail_fraction=tail)
    reference = phase1.cascade_delivery(aspect, iadf, REFERENCE_ENERGY_EV)
    edges = np.linspace(0.0, aspect, max(24, int(24 * aspect)) + 1)
    mine = cascade_hole_delivery(0.5, aspect, edges, iadf, REFERENCE_ENERGY_EV)
    assert mine["total_bottom"] == pytest.approx(reference["total_bottom"], rel=1e-12)
    assert mine["direct_bottom"] == pytest.approx(reference["direct_bottom"], rel=1e-12)


def test_cascade_weight_budget_closes():
    """Bottom + reacted-at-wall + thermalised = every particle that entered."""
    iadf = kim_2025_reference_iadf(tail_fraction=0.65)
    edges = np.linspace(0.0, 40.0, 481)
    record = cascade_hole_delivery(0.5, 40.0, edges, iadf, REFERENCE_ENERGY_EV)
    assert abs(record["closure_residual"]) < 1e-12
    # Arrivals exceed unity at depth: one particle strikes many bands, and the
    # surface sees each strike.  That is re-arrival counting, not a leak.
    assert record["wall_hist"].sum() > record["absorbed_per_band"].sum()


# --- conservation --------------------------------------------------------------


def test_diffuse_solve_conserves_with_born_source():
    """The E8 birth term must be conserved like any other supply."""
    edges = np.linspace(0.0, 8.0, 97)
    enclosure = build_hole_enclosure(0.5, 8.0, edges)
    faces = enclosure.area.size
    sticking = np.full(faces, 0.2)
    born = np.zeros(faces)
    born[10] = 3.0
    arrival, absorbed, escaped = solve_diffuse_hole_delivery(
        enclosure, sticking, 5.0, born_rate=born)
    assert absorbed.sum() + escaped == pytest.approx(5.0 + born.sum(), rel=1e-12)
    assert np.all(arrival >= 0.0)


def test_step_ledgers_close():
    mouth, ion = _boundary_fluxes()
    oxide, _mask = build_krueger_2024_mixed_layer_mechanisms()
    geometry = HoleGeometry.straight(45e-9, 4.5e-7, 3e-8)
    state = HoleEvolutionState(geometry, _bare(geometry.band_count + 1))
    _next, record = advance_hole_step(
        state, oxide, mouth_flux_m2_s=mouth, ion_flux_m2_s=ion,
        iadf=kim_2025_reference_iadf(tail_fraction=0.65),
        energy_eV=REFERENCE_ENERGY_EV, dt_s=1e-3, cascade_kwargs=FAST_CASCADE)
    assert record["enclosure_closure_residual"] < 1e-9
    assert record["neutral_balance_residual"] < 1e-9
    assert abs(record["cascade_closure_residual"]) < 1e-9


# --- degenerate limit ----------------------------------------------------------


def test_wide_short_hole_reproduces_blanket_rate():
    """As the hole opens out, the floor must see the unobstructed plasma.

    A hole 200x wider than it is deep has no shadowing and no cascade to speak
    of, so its floor chemistry must reduce to the 0-D blanket fed by the same
    mouth fluxes and the same beam -- the driver's transport layer contributing
    nothing but unity.
    """
    mouth, ion = _boundary_fluxes()
    oxide, _mask = build_krueger_2024_mixed_layer_mechanisms()
    iadf = kim_2025_reference_iadf(tail_fraction=0.65)
    radius, depth = 1e-5, 1e-7
    geometry = HoleGeometry.straight(radius, depth, depth)
    state = HoleEvolutionState(geometry, _bare(geometry.band_count + 1))
    _next, record = advance_hole_step(
        state, oxide, mouth_flux_m2_s=mouth, ion_flux_m2_s=ion, iadf=iadf,
        energy_eV=REFERENCE_ENERGY_EV, dt_s=1e-3, cascade_kwargs=FAST_CASCADE)

    # 0-D reference: identical construction with delivery forced to unity.
    polar_deg, polar_w = iadf.polar_quadrature(REFERENCE_ENERGY_EV, n_polar=48)
    cosine = np.cos(np.deg2rad(polar_deg))
    edges = np.linspace(0.0, 1.0, 9)
    index = np.clip(np.digitize(cosine, edges) - 1, 0, 7)
    histogram = np.zeros(8)
    np.add.at(histogram, index, polar_w)
    centres = 0.5 * (edges[:-1] + edges[1:])
    live = histogram > 0.0
    events = FaceResolvedEnergeticFlux(
        "ions", 1, np.zeros(int(live.sum()), dtype=int),
        ion * histogram[live], np.full(int(live.sum()), REFERENCE_ENERGY_EV),
        centres[live])
    blanket = oxide.advance(
        _bare(1), SurfaceFluxes({k: np.full(1, v) for k, v in mouth.items()},
                                (events,)),
        1e-3, strict=False)
    reference = float(np.asarray(blanket.etch_velocity_m_s)[0])
    assert reference > 0.0
    assert record["floor_etch_velocity_m_s"] == pytest.approx(reference, rel=0.01)


# --- evolution -----------------------------------------------------------------


def test_smoke_evolution_advances_floor_and_passivates_walls():
    """An AR-10 hole must etch downward while its walls take on film."""
    mouth, ion = _boundary_fluxes()
    oxide, _mask = build_krueger_2024_mixed_layer_mechanisms()
    geometry = HoleGeometry.straight(45e-9, 9e-7, 6e-8)
    state = HoleEvolutionState(geometry, _bare(geometry.band_count + 1))
    final, records, reason = evolve_hole(
        state, oxide, mouth_flux_m2_s=mouth, ion_flux_m2_s=ion,
        iadf=kim_2025_reference_iadf(tail_fraction=0.65),
        energy_eV=REFERENCE_ENERGY_EV, duration_s=0.05, dt_s=1e-2,
        cascade_kwargs=FAST_CASCADE)

    assert reason == "duration"
    assert len(records) == 5
    depths = [row["floor_depth_m"] for row in records]
    assert np.all(np.diff(depths) > 0.0)                  # the floor advances
    assert final.geometry.floor_depth > geometry.floor_depth
    assert records[-1]["max_wall_growth_velocity_m_s"] > 0.0   # walls passivate
    assert records[-1]["floor_etch_velocity_m_s"] > 0.0
    for row in records:
        assert row["neutral_balance_residual"] < 1e-9
        assert abs(row["cascade_closure_residual"]) < 1e-9
    # Film on the wall must not exceed the etched volume budget: the hole grows.
    assert records[-1]["solid_volume_removed_m3"] > 0.0


def test_evolution_stops_on_the_declared_straightness_envelope():
    """Past the declared tolerance the driver stops rather than extrapolating."""
    mouth, ion = _boundary_fluxes()
    oxide, _mask = build_krueger_2024_mixed_layer_mechanisms()
    geometry = HoleGeometry.straight(45e-9, 4.5e-7, 3e-8)
    state = HoleEvolutionState(geometry, _bare(geometry.band_count + 1))
    _final, records, reason = evolve_hole(
        state, oxide, mouth_flux_m2_s=mouth, ion_flux_m2_s=ion,
        iadf=kim_2025_reference_iadf(tail_fraction=0.65),
        energy_eV=REFERENCE_ENERGY_EV, duration_s=10.0, dt_s=1e-2,
        straightness_tolerance=1e-4, cascade_kwargs=FAST_CASCADE)
    assert reason == "profile_left_straight_wall_envelope"
    assert records[-1]["straightness_deviation"] > 1e-4


def test_geometry_exposes_new_bands_as_the_floor_descends():
    geometry = HoleGeometry.straight(1.0, 2.0, 1.0)
    assert geometry.band_count == 2
    deeper = HoleGeometry(1.0, np.full(3, 1.0), 2.5, 1.0)
    assert deeper.band_count == 3
    assert deeper.band_exposed_height[-1] == pytest.approx(0.5)
    assert deeper.band_area[-1] == pytest.approx(2.0 * pi * 1.0 * 0.5)
    with pytest.raises(ValueError):
        HoleGeometry(1.0, np.full(2, 1.0), 2.5, 1.0)
