"""Gates for the O-radical polymer-etch channel budget.

Krueger's published row is ``O(g) + P(s) -> products`` at ``p_ox`` per
collision with EXPOSED POLYMER: a per-cell probability that already subsumes
which atom in the polymer unit the oxygen lands on.  The film-oxidation term
must therefore gate on film coverage alone, and each reaction must remove one
whole polymer unit (one carbon plus the local film F/C ratio).

Scaling the carbon removal by the film composition ``x_c`` in addition
throttled the channel by ``1/x_c`` -- 2.69x at the F/C = 1.69 composition the
evolution runs actually reach -- which is why raising ``p_ox`` by 48% between
ml16a and ml16b moved the neck only 11%.  See RESULTS_O_CHANNEL_2026-08-04.md
and RESULTS_WALL_SLOPE_FALSIFICATION_2026-08-04.md.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from petch.mixed_layer import (
    MixedLayerParams,
    MixedLayerState,
    SurfaceFluxes,
    _FC_FILM_ATOM_DENSITY_M3,
    _MONOLAYER_AREAL_M2,
    step,
)
from petch.mixed_layer_mechanism import KRUEGER_2024_DEPOSITION_ON_POLYMER

# Krueger 2024 Table-I base fluxes (m^-2 s^-1), the deck the pilot runs.
_TABLE_I_FLUX = {
    "CF": 4.4e20, "CF2": 9.4e20, "CF3": 8.4e19, "C2F3": 6.8e20,
    "C3F4": 9.5e20, "O": 7.7e20,
}
# Deposition events per second on polymer: sum_i p_i J_i.  Both this and the
# oxidation rate are per-cell quantities, so their ratio is the geometry-free
# budget fraction Krueger's mechanism implies.
_DEPOSITION_EVENTS = sum(
    KRUEGER_2024_DEPOSITION_ON_POLYMER.get(name, 0.0) * flux
    for name, flux in _TABLE_I_FLUX.items())
_TARGET_RATIO = 0.1953  # p_ox * J_O / sum_i p_i J_i at Table-I fluxes

# Film composition measured on the ml16a evolution checkpoint (F/C = 1.690,
# areal total 1.173e21 m^-2 -> 15.6 nm, i.e. theta_film = 1 to machine zero).
_THICK_FILM_C = 4.36e20
_THICK_FILM_F = 7.37e20


def _thick_film_state():
    return MixedLayerState(n_c_film=_THICK_FILM_C, n_f_film=_THICK_FILM_F)


def _theta_film(state):
    thickness_nm = float(np.asarray(state.film_thickness_nm()))
    return 1.0 - np.exp(-thickness_nm * 1e-9 * _FC_FILM_ATOM_DENSITY_M3
                        / _MONOLAYER_AREAL_M2)


def test_thick_film_is_fully_covered():
    """The pinch regime the gate below describes really is theta_film = 1."""
    assert _theta_film(_thick_film_state()) == pytest.approx(1.0, abs=1e-12)


def test_o_removal_over_deposition_matches_krueger_budget():
    """O-removal / deposition must reproduce the published 0.195 budget.

    Both sides are per-cell rates at the same (thermal, geometry-free) fluxes,
    so this ratio is a property of the mechanism alone.  Gate: within 10%.
    """
    params = MixedLayerParams(substrate="carbon")
    ox_c = params.oxidation_probability * _TABLE_I_FLUX["O"] * _theta_film(
        _thick_film_state())
    ratio = ox_c / _DEPOSITION_EVENTS
    assert ratio == pytest.approx(_TARGET_RATIO, rel=0.10)


def test_oxidation_removes_one_polymer_unit_per_collision():
    """Each reactive O collision removes 1 C plus the film's F/C, not a
    composition-weighted fraction of one atom.

    Measured through the module: with only an oxygen flux acting on a thick
    film, the carbon lost per unit time must equal p_ox * J_O (one carbon per
    reactive collision), and the fluorine lost must equal that times F/C.
    """
    params = MixedLayerParams(substrate="carbon")
    state = _thick_film_state()
    fluxes = SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=_TABLE_I_FLUX["O"],
        ion_flux=0.0, ion_energy_eV=0.0)
    dt = 1e-6
    result = step(state, fluxes, dt, params)
    lost_c = (state.n_c_film - float(np.asarray(result.state.n_c_film))) / dt
    lost_f = (state.n_f_film - float(np.asarray(result.state.n_f_film))) / dt
    expected_c = params.oxidation_probability * _TABLE_I_FLUX["O"]
    f_over_c = _THICK_FILM_F / _THICK_FILM_C
    # rel=1e-6: the rate is a difference of ~1e21 reservoirs over dt=1e-6,
    # which leaves ~9 significant digits after float64 cancellation.
    assert lost_c == pytest.approx(expected_c, rel=1e-6)
    assert lost_f == pytest.approx(expected_c * f_over_c, rel=1e-6)
    # And the composition-throttled form is what we are NOT doing: it would
    # remove 2.69x less carbon at this composition.
    throttled = expected_c * _THICK_FILM_C / (_THICK_FILM_C + _THICK_FILM_F)
    assert lost_c / throttled == pytest.approx(1.0 + f_over_c, rel=1e-6)


def test_o_channel_responds_proportionally_to_p_ox():
    """The channel must not be inert: a +48% p_ox change (the ml16a->ml16b
    step) must move film carbon removal by +48%, not by a throttled fraction.
    """
    state = _thick_film_state()
    fluxes = SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=_TABLE_I_FLUX["O"],
        ion_flux=0.0, ion_energy_eV=0.0)
    dt = 1e-6

    def carbon_loss(p_ox):
        params = MixedLayerParams(substrate="carbon", oxidation_probability=p_ox)
        result = step(state, fluxes, dt, params)
        return (state.n_c_film - float(np.asarray(result.state.n_c_film))) / dt

    low, high = carbon_loss(0.0423), carbon_loss(0.0628)
    assert high / low == pytest.approx(0.0628 / 0.0423, rel=1e-6)


def test_ledger_still_closes_with_the_stronger_channel():
    """Element ledgers must close to machine precision at the new rate."""
    params = MixedLayerParams(substrate="carbon")
    fluxes = SurfaceFluxes(
        precursor_flux=3.0e19, fluorine_flux=2.0e20,
        oxygen_flux=_TABLE_I_FLUX["O"], ion_flux=6.0e18, ion_energy_eV=1000.0)
    state = _thick_film_state()
    worst = 0.0
    # dt = 1e-4 is the repo's standard ledger step: the storage term is a
    # difference of ~1e21 reservoirs, so a smaller dt divides float64
    # cancellation noise by dt and reports it as a spurious residual (~3e-9 at
    # dt = 1e-6, identically so before this change).
    for _ in range(200):
        result = step(state, fluxes, 1e-4, params)
        for key in ("fluorine", "carbon", "silicon", "oxygen"):
            worst = max(worst, abs(result.ledger_residuals[key]))
        state = result.state
    assert worst < 1e-9
