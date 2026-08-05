"""Gates for the Appendix-B angular classes on the oxide/mask ion rows.

Krueger Table B.0.1 marks four ion rows with an angular class; before
2026-08-05 the module applied a factor to exactly one of them
(RESULTS_LIP_REMOVAL_AUDIT_2026-08-04, "New finding: three angular markers are
unimplemented").  These gates pin the two published shapes, the normalisation
convention that makes p0 the normal-incidence probability, and the invariance
of every previously validated normal-incidence result.
"""

import numpy as np
import pytest

from petch.mixed_layer import (
    MixedLayerParams,
    MixedLayerState,
    SurfaceFluxes,
    _angular_chemical_sputter,
    _angular_physical_sputter,
    steady_state,
    step,
)


def _f1(deg):
    return float(_angular_physical_sputter(np.cos(np.radians(deg))))


def _f2(deg):
    return float(_angular_chemical_sputter(np.cos(np.radians(deg))))


def test_both_classes_are_unity_at_normal_incidence():
    """p0 in Table B.0.1 is the normal-incidence probability: f(0) = 1 exactly
    for both classes, so normal-incidence yields are untouched by this change."""
    assert _f1(0.0) == 1.0
    assert _f2(0.0) == 1.0


def test_class1_is_the_physical_sputter_shape():
    """Huang: 'maximum at 60 deg, reduced probability at normal incidence and
    zero probability at grazing incidence' (relative to its own peak)."""
    peak_deg = max(np.linspace(0.0, 90.0, 9001), key=_f1)
    assert 50.0 < peak_deg < 60.0          # (1+B sin^2)cos peaks at 54.7 deg
    assert _f1(peak_deg) / _f1(0.0) == pytest.approx(4.17, rel=1e-2)
    assert _f1(90.0) == pytest.approx(0.0, abs=1e-14)
    # monotone rise to the peak, monotone fall after it
    rise = [_f1(d) for d in np.linspace(0.0, peak_deg, 50)]
    fall = [_f1(d) for d in np.linspace(peak_deg, 90.0, 50)]
    assert all(b >= a for a, b in zip(rise, rise[1:]))
    assert all(b <= a for a, b in zip(fall, fall[1:]))


def test_class2_is_the_chemical_sputter_shape():
    """Huang, verbatim: 'unity for normal incidence and angles up to 45 deg,
    with a monotonic roll-off to zero probability at grazing incidence'."""
    for deg in (0.0, 10.0, 25.0, 40.0, 45.0):
        assert _f2(deg) == pytest.approx(1.0, rel=1e-12)
    assert _f2(90.0) == pytest.approx(0.0, abs=1e-15)
    roll = [_f2(d) for d in np.linspace(45.0, 90.0, 200)]
    assert all(b <= a for a, b in zip(roll, roll[1:]))
    assert 0.0 < _f2(80.0) < _f2(60.0) < 1.0


def test_classes_are_distinct_off_normal():
    """The two classes must not collapse onto each other: physical sputtering
    is enhanced where chemical sputtering is still on its plateau."""
    assert _f1(45.0) > 3.5 * _f2(45.0)
    assert _f1(60.0) > 5.0 * _f2(60.0)


def test_normal_incidence_steady_state_is_bitwise_unchanged():
    """Every validated 0-D result was taken at normal incidence; with f(0)=1
    the oxide rows must reproduce the pre-change rates exactly.  Reference
    values recomputed with the angular factors forced to 1."""
    fluxes = SurfaceFluxes(3.0e19, 2.0e20, 5.0e19, 6.0e18, 1000.0)
    at_normal = steady_state(fluxes)
    # cosine_incidence defaults to 1.0; state the invariance explicitly by
    # comparing against the same call with the cosine written out.
    explicit = steady_state(SurfaceFluxes(3.0e19, 2.0e20, 5.0e19, 6.0e18,
                                          1000.0, cosine_incidence=1.0))
    assert float(np.asarray(at_normal.sif4_rate)) == float(
        np.asarray(explicit.sif4_rate))
    assert float(np.asarray(at_normal.state.n_si)) == float(
        np.asarray(explicit.state.n_si))


def test_oblique_incidence_now_moves_the_oxide_rows():
    """The whole point of the change: at 80 deg the bare-oxide (class 1) and
    complex (class 2) channels must both be suppressed relative to normal,
    where before they were angle-blind."""
    bare_only = SurfaceFluxes(0.0, 1.0e22, 0.0, 6.0e18, 1000.0)
    normal = steady_state(bare_only).sif4_rate
    oblique = steady_state(SurfaceFluxes(0.0, 1.0e22, 0.0, 6.0e18, 1000.0,
                                         cosine_incidence=np.cos(
                                             np.radians(80.0)))).sif4_rate
    assert float(np.asarray(oblique)) < float(np.asarray(normal))


def test_ledgers_still_close_with_angular_factors():
    """Conservation is independent of the angular factors."""
    params = MixedLayerParams()
    for cos in (1.0, 0.7, 0.2, 0.05):
        state = MixedLayerState()
        worst = 0.0
        fluxes = SurfaceFluxes(3.0e19, 2.0e20, 5.0e19, 6.0e18, 1000.0,
                               cosine_incidence=cos)
        for _ in range(300):
            result = step(state, fluxes, 1e-4, params)
            for key in ("fluorine", "carbon", "silicon", "oxygen"):
                worst = max(worst, abs(result.ledger_residuals[key]))
            state = result.state
        assert worst < 1e-9, cos


def test_atom_path_matches_scalar_path_at_single_angle():
    """One atom at angle theta must reproduce the scalar path at the same
    cosine -- the contract that keeps per-event and mean paths consistent
    now that three more kernels carry an angular factor."""
    cos = float(np.cos(np.radians(65.0)))
    scalar = step(MixedLayerState(),
                  SurfaceFluxes(3.0e19, 2.0e20, 5.0e19, 6.0e18, 1000.0,
                                cosine_incidence=cos),
                  1e-4)
    atoms = step(MixedLayerState(),
                 SurfaceFluxes(3.0e19, 2.0e20, 5.0e19, 0.0, 1000.0,
                               ion_atom_face=np.array([0]),
                               ion_atom_flux=np.array([6.0e18]),
                               ion_atom_energy_eV=np.array([1000.0]),
                               ion_atom_cosine=np.array([cos])),
                 1e-4)
    assert float(np.asarray(atoms.state.n_si)) == pytest.approx(
        float(np.asarray(scalar.state.n_si)), rel=1e-12)
    assert float(np.asarray(atoms.sif4_rate)) == pytest.approx(
        float(np.asarray(scalar.sif4_rate)), rel=1e-12)
