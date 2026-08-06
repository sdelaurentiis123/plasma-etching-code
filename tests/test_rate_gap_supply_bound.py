"""The rate gap is not a channel-magnitude error.

Two independent gates, both against measurements archived in
``research_sources/``:

1.  CROSS-EXPERIMENT YIELD CHECK.  The two SiO2 ion channels carry Gray's
    beam-measured laws (MIT thesis 1993, Ar+/F on SiO2, QCM).  Karahashi
    measured the same quantities on a completely different apparatus --
    mass-analyzed single-species CFx+ ion beam, Osaka, UHV, no radical flux --
    and published absolute yields at 1000 eV
    (``research_sources/thesis_extracts/karahashi_2007_sio2_cfx_ionbeam.txt``
    L118-127, the full-text open-access review of JVST A 22, 1166 (2004)):

        "F+ イオンに関しては 0.3 molecules／ion と低く ... CF3+ イオンの場合
         エッチングイールドの値が 1.5 molecules／ion となる"
        = F+ gives 0.3 molecules/ion (physical-sputter-like, close to Ne+);
          CF3+ gives 1.5 molecules/ion (the chemically enhanced maximum).

    petch's laws must reproduce that bracket without ever having been fitted
    to it.  This is the strongest available check on the absolute magnitudes.

2.  SUPPLY BOUND.  At the floor of a HAR feature the complex channel is
    carbon/fluorine-supply limited, so its magnitude has no authority over the
    etch rate there.  Scaling the magnitude 8x must move the floor rate by
    less than 25%.  This is what refutes "calibrate the chemical-channel
    magnitude to close the depth gate": no value of that constant can.

See RESULTS_RATE_GAP_CLOSURE_2026-08-06.md.
"""
import sys

import numpy as np
import pytest

sys.path.insert(0, "src")

import petch.mixed_layer as ml
from petch.chemistry_deck import build_mixed_layer_mechanisms_from_deck
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes

# Krueger et al., JVST A 42, 043008 (2024), Table I -- "Ions 1.2 x 10^16 cm-2 s-1".
_ION_PUBLISHED_M2_S = 1.2e20
_TABLE_I_NEUTRALS = {
    "CF": 4.4e20, "CF2": 9.4e20, "C2F3": 6.8e20, "CF3": 8.4e19,
    "O": 7.7e20, "C3F4": 9.5e20,
}

# Karahashi 2007 (full text archived), Fig. 3, SiO2 at 1000 eV.
_KARAHASHI_F_PLUS = 0.3          # molecules/ion, physical-sputter-like
_KARAHASHI_CF3_PLUS = 1.5        # molecules/ion, chemically enhanced maximum


def _floor_rate_nm_s(complex_scale, *, neutral_delivery=0.10,
                     ion_delivery=0.70, energy_eV=3406.0, steps=120, dt=2.0):
    """Steady oxide recession at HAR-floor delivery, complex yield scaled."""
    original = ml._GRAY_BETA_A
    ml._GRAY_BETA_A = original * float(complex_scale)
    try:
        oxide, _ = build_mixed_layer_mechanisms_from_deck()
        neutrals = {k: v * neutral_delivery for k, v in _TABLE_I_NEUTRALS.items()}
        ion = EnergeticFlux(
            name="Ar+", flux_m2_s=_ION_PUBLISHED_M2_S * ion_delivery,
            energy_eV=np.array([energy_eV]), cosine_incidence=np.array([1.0]),
            weight=np.array([1.0]))
        fluxes = SurfaceFluxes(neutral_flux_m2_s=neutrals,
                              energetic_fluxes=(ion,))
        state = oxide.initial_state(())
        result = None
        for _ in range(steps):
            result = oxide.advance(state, fluxes, dt)
            state = result.state
        return float(np.asarray(result.etch_velocity_m_s)) * 1e9
    finally:
        ml._GRAY_BETA_A = original


def test_physical_channel_matches_karahashi_f_plus():
    """Gray's sputter law vs Karahashi's F+ measurement at 1000 eV."""
    petch = float(ml._bare_sputter_yield(np.array([1000.0])))
    # 0.381 vs 0.3 measured: within 30%, two different beams, no fitting.
    assert petch == pytest.approx(_KARAHASHI_F_PLUS, rel=0.30)


def test_chemical_channel_matches_karahashi_cf3_plus():
    """Gray's beta_e vs Karahashi's CF3+ measurement at 1000 eV.

    Two independent apparatus, eleven years apart, never cross-fitted:
    petch reads 1.570 against 1.5 molecules/ion measured -- 4.7%.  The channel's
    ABSOLUTE per-ion magnitude is therefore corroborated to within experimental
    scatter.  (Not in tension with the earlier "2.8x too weak" reading, which
    compared against Gray's saturated plateau at 350 eV -- a coverage-curve
    quantity rather than an absolute yield.)
    """
    petch = float(ml._complex_yield(np.array([1000.0])))
    assert petch == pytest.approx(_KARAHASHI_CF3_PLUS, rel=0.10)


def test_complex_magnitude_cannot_close_the_floor_rate():
    """The floor is supply-bounded: an 8x magnitude buys under 25%.

    This is the gate that refutes a chemical-channel magnitude calibration as
    the depth fix.  Closing the measured depth gate needs 2.38x on the feature
    average; the magnitude cannot deliver it at any value.
    """
    base = _floor_rate_nm_s(1.0)
    scaled = _floor_rate_nm_s(8.0)
    assert base > 0.0
    assert scaled / base < 1.25, (
        f"floor rate moved {scaled / base:.2f}x under an 8x complex-yield "
        "scale; the supply bound this gate asserts has been broken")


def test_floor_rate_saturates_in_the_magnitude():
    """4x and 8x agree: the channel is at its supply ceiling, not its law."""
    four = _floor_rate_nm_s(4.0)
    eight = _floor_rate_nm_s(8.0)
    assert four == pytest.approx(eight, rel=0.02)
