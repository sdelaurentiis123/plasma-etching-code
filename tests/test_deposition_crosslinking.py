"""Gates for deposition-driven crosslinking (Krueger Table 6.2 / sec. 2.2.3).

His crosslink module creates PC *during deposition* ("Crosslinking occurs
during the deposition of eligible materials"); ions only BREAK crosslinks
(`CF(xs) + M -> CF(s) + M`, 0.3 @ 8 eV).  Creation therefore does not collapse
on a near-vertical wall, while ion-driven creation does (~200x through the
double cosine).  See RESULTS_LIP_CROSSLINK_2026-08-04.md.
"""

import numpy as np
import pytest

from petch.mixed_layer import MixedLayerParams, MixedLayerState, step

_SRC = dict(precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=0.0)


def _relax(fluxes, params=None, dt=2.0, steps=3000):
    from petch.mixed_layer import SurfaceFluxes
    params = params or MixedLayerParams(substrate="carbon")
    state = MixedLayerState()
    result = None
    for _ in range(steps):
        result = step(state, fluxes, dt, params)
        state = result.state
    total = float(np.asarray(state.n_c_film)) + float(np.asarray(state.n_f_film))
    return state, result, total


def _fluxes(dep_c, dep_f, ion_flux, cosine, oxygen=0.0):
    from petch.mixed_layer import SurfaceFluxes
    return SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=oxygen,
        ion_flux=ion_flux, ion_energy_eV=1500.0, cosine_incidence=cosine,
        film_deposition_carbon_flux=dep_c,
        film_deposition_fluorine_flux=dep_f,
        substrate_deposition_carbon_flux=dep_c,
        substrate_deposition_fluorine_flux=dep_f)


def test_ion_free_steady_state_is_the_row_stoichiometry():
    """`P(s)+P(s) -> PC(s)+PC(s)` converts two units per deposition event, so
    with no ions to break bonds the crosslinked fraction relaxes to the exact
    analytic balance 2(1-x) = x, i.e. x = 2/3.  No constant is fitted: the 2 is
    the published row's own stoichiometry."""
    state, _, total = _relax(_fluxes(1.0e20, 1.7e20, 0.0, 1.0))
    x = float(np.asarray(state.n_xl_film)) / total
    assert x == pytest.approx(2.0 / 3.0, rel=0.05)


def test_crosslinking_survives_grazing_incidence():
    """The mechanism's signature: on a near-vertical wall the ion channel
    collapses but deposition does not, so the film still crosslinks.  Ion-dose
    creation alone leaves the lip fresh (measured 0.163 in the feature runs)."""
    cos_lip = float(np.sin(np.deg2rad(0.472)))          # audited top-band tilt
    state, _, total = _relax(
        _fluxes(1.0e20, 1.7e20, 9.6e19 * cos_lip * 0.742, cos_lip))
    x_lip = float(np.asarray(state.n_xl_film)) / total
    assert x_lip > 0.6
    # ...and the same fluxes at normal incidence stay crosslinked too, so the
    # channel is not a grazing-only artefact.
    state_n, _, total_n = _relax(_fluxes(1.0e20, 1.7e20, 9.6e19, 1.0))
    assert float(np.asarray(state_n.n_xl_film)) / total_n > 0.6


def test_crosslinked_lip_slows_film_growth():
    """Crosslinked film attaches at the published 0.02 row instead of 0.1, so a
    crosslinked lip must grow markedly slower than a fresh one at identical
    delivered flux -- the lever that sets mouth closure."""
    from petch.mixed_layer import SurfaceFluxes
    cos_lip = float(np.sin(np.deg2rad(0.472)))
    fresh = _fluxes(1.0e20, 1.7e20, 9.6e19 * cos_lip * 0.742, cos_lip)
    # Same fluxes, but the crosslinked attachment row supplied explicitly.
    crosslinked = SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=0.0,
        ion_flux=9.6e19 * cos_lip * 0.742, ion_energy_eV=1500.0,
        cosine_incidence=cos_lip,
        film_deposition_carbon_flux=1.0e20,
        film_deposition_fluorine_flux=1.7e20,
        crosslinked_deposition_carbon_flux=2.0e19,
        crosslinked_deposition_fluorine_flux=3.4e19,
        substrate_deposition_carbon_flux=1.0e20,
        substrate_deposition_fluorine_flux=1.7e20)
    _, _, total_blend = _relax(crosslinked, steps=1500)
    _, _, total_fresh = _relax(fresh, steps=1500)
    assert total_blend < total_fresh


def test_crosslinked_atoms_never_exceed_the_film():
    """State invariant: n_xl is a subset of the film inventory at all times,
    including the first steps when the new channel is at its most aggressive."""
    from petch.mixed_layer import SurfaceFluxes
    params = MixedLayerParams(substrate="carbon")
    state = MixedLayerState()
    fluxes = _fluxes(5.0e20, 8.5e20, 1.0e18, 0.05)
    for _ in range(200):
        result = step(state, fluxes, 1e-3, params)
        state = result.state
        total = (float(np.asarray(state.n_c_film))
                 + float(np.asarray(state.n_f_film)))
        assert float(np.asarray(state.n_xl_film)) <= total + 1e-6 * max(total, 1.0)


def test_no_crosslinking_without_deposition_or_ions():
    """Eligibility gate: with no film and no deposition there is nothing to
    bond to, so the channel is exactly inert (no spontaneous creation)."""
    state, _, _ = _relax(_fluxes(0.0, 0.0, 0.0, 1.0), steps=50)
    assert float(np.asarray(state.n_xl_film)) == 0.0
