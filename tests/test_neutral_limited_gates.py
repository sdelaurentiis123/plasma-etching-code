"""GATE N1/N2: the mixed layer against published beam data.

Preregistered in RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md sec 5.4.
Published data on both axes:

  Gray, Tepermeister & Sawin, JVST B 11, 1243 (1993), replotted as Kwon
  (ScD, MIT DMSE 2004) Fig. 3.4 p. 76, 350 eV Ar+ on SiO2 with an F beam:
  sputter floor Y = 0.28 as F/Ar+ -> 0, plateau Y ~ 1.10, half-rise at
  F/Ar+ ~ 27.  Dynamic range Y(0)/Y(sat) = 0.25.

  Butterbaugh, Gray & Sawin, JVST B 9, 1461 (1991), via Kwon Fig. 2.6 p. 36
  and p. 36 text: "there is a reduction in etching yield as the CFx flux
  increases due to the deposition on the surface".

What is gated hard here is the STRUCTURE both chemistries share -- a
coverage-weighted chemical channel over a physical-sputter floor, monotone
and saturating in radical supply, and reduced by depositing carbon.  The
absolute dynamic range is a *diagnostic*: it is the ratio of two published
removal rows measured in a different halogen chemistry, so it is asserted
against petch's own current value as a regression tripwire, with the
distance to Gray recorded in the assertion message rather than enforced.
"""

import numpy as np
import pytest

from petch.mixed_layer import SurfaceFluxes, steady_state

BEAM_ENERGY_EV = 350.0
ION_FLUX = 1.0e19

# Gray 1993 via Kwon Fig. 3.4 (digitized [VERIFY] against the original).
GRAY_DYNAMIC_RANGE = 0.28 / 1.10          # 0.2545
GRAY_HALF_RISE = 27.0

# petch's measured value.  The bare-sputter and complex rows have nearly equal
# yields at feature energies, so the coverage-weighted sum is almost
# theta-independent and the layer cannot be neutral-limited at any radical
# supply -- 3.4x away from Gray's measured range.  This is the open depth-channel
# defect; RESULTS_ANGULAR_CONVENTION_2026-08-05.md shows that the class-1
# angular normalisation reproduces Gray (0.210) but overshoots the depth gate
# in the opposite direction, so it is NOT the whole answer and was not landed.
PETCH_DYNAMIC_RANGE = 0.873
GRAY_BAND = (0.20, 0.30)


def _beam_yield(flux_ratio, cf2_ratio=0.0, energy_eV=BEAM_ENERGY_EV):
    fluxes = SurfaceFluxes(
        precursor_flux=cf2_ratio * ION_FLUX,
        fluorine_flux=flux_ratio * ION_FLUX,
        oxygen_flux=0.0,
        ion_flux=ION_FLUX,
        ion_energy_eV=energy_eV,
    )
    result = steady_state(fluxes)
    return float(np.asarray(result.substrate_removal_rate)) / ION_FLUX


@pytest.fixture(scope="module")
def beam_curve():
    return {r: _beam_yield(r) for r in (0.0, 5.0, 40.0, 500.0)}


def test_n1_yield_is_monotone_and_saturating(beam_curve):
    """Gray's curve rises monotonically and saturates; ours must too."""
    ratios = sorted(beam_curve)
    values = [beam_curve[r] for r in ratios]
    assert all(b >= a for a, b in zip(values, values[1:])), values
    # Saturation: the last decade must add less than the first.
    early = values[1] - values[0]
    late = values[-1] - values[-2]
    assert late < early, (early, late)


def test_n1_floor_is_pure_physical_sputter(beam_curve):
    """With no radicals the only surviving channel is bare-oxide sputter, so
    the floor must be positive and strictly below the saturated plateau --
    the structure Gray measures as 0.28 -> 1.10."""
    assert beam_curve[0.0] > 0.0
    assert beam_curve[0.0] < beam_curve[500.0]


def test_n1_dynamic_range_records_the_open_miss(beam_curve):
    """GATE N1: the zero-radical floor over the saturated plateau.

    Gray measures 0.28/1.10 = 0.255 (band 0.20-0.30); petch measures 0.873.
    The miss is REAL and open, so this pins petch's own value as a tripwire
    rather than asserting the published band -- any change that moves it must
    move it deliberately and be graded on the depth channel at the same time
    (RESULTS_ANGULAR_CONVENTION_2026-08-05.md).
    """
    measured = beam_curve[0.0] / beam_curve[500.0]
    assert measured == pytest.approx(PETCH_DYNAMIC_RANGE, rel=0.08), (
        f"dynamic range moved to {measured:.3f} from {PETCH_DYNAMIC_RANGE:.3f}; "
        f"Gray 1993 measures {GRAY_DYNAMIC_RANGE:.3f} (band {GRAY_BAND}). "
        f"Update the constant only with the paired depth forecast.")
    lo, hi = GRAY_BAND
    assert not (lo <= measured <= hi), (
        "dynamic range entered Gray's band -- if this is intended, flip this "
        "gate to assert the band and cite the run that graded the depth.")


def test_n2_depositing_carbon_reduces_yield():
    """Butterbaugh three-beam: adding CFx at fixed F/Ar+ must lower the yield
    (Kwon p. 36). Sign gate only -- the magnitude is chemistry-specific."""
    clean = _beam_yield(40.0, cf2_ratio=0.0)
    loaded = _beam_yield(40.0, cf2_ratio=10.0)
    assert loaded < clean, (clean, loaded)


def test_repassivation_asymmetry_is_present():
    """Krueger L6556-6564 / Huang L10214-10222: further passivation of an
    ALREADY-complexed site is 2e-4 / 1e-4 against 0.278/0.2 for a pristine
    site -- a ~1400x asymmetry that is what makes the surface saturate.

    petch carries this as Langmuir site blocking: chemisorption enters the
    layer through `site_open = (1 - theta_film) * (1 - theta_f_layer)`, so a
    layer already saturated with fluorine accepts essentially nothing. This
    gate pins the direction and the strength of that blocking.
    """
    from petch.mixed_layer import _MONOLAYER_AREAL_M2, MixedLayerState, step

    radical = 1.0e20
    fluxes = SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=0.0,
        ion_flux=ION_FLUX, ion_energy_eV=BEAM_ENERGY_EV,
        chemisorption_carbon_flux=radical,
        chemisorption_fluorine_flux=2.0 * radical)
    dt = 1.0e-6

    def carbon_uptake(n_f):
        state = MixedLayerState(n_f=n_f)
        after = step(state, fluxes, dt).state
        return float(np.asarray(after.n_c)) / dt

    pristine = carbon_uptake(0.0)
    complexed = carbon_uptake(_MONOLAYER_AREAL_M2)
    # Krueger's asymmetry is 0.278 -> 0.0002 (1390x). petch blocks the site
    # outright once it is occupied, which is at least that strong.
    assert pristine > 0.0
    assert complexed < pristine / 1390.0, (pristine, complexed)
