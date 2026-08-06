"""Sensitivity gates for the incumbent neutral-assisted complex channel.

These tests establish one narrow result: once the existing mixed-layer
complex channel becomes neutral-supply limited at a HAR floor, multiplying
its Gray-derived magnitude cannot close the feature-depth gap.

They do *not* bound missing channels. Karahashi's CFx+ measurements are
species-resolved reactive-ion totals, whereas the incumbent kernels are
Ar-like physical removal plus neutral-F-assisted removal. Comparing one
formula value to one CF3+ marker was not a cross-experiment validation and has
been removed. The actual species ladder is gated end-to-end in
``tests/test_reactive_ion_beam.py``.
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


def test_complex_magnitude_cannot_close_the_floor_rate():
    """The floor is supply-bounded: an 8x magnitude buys under 25%.

    This rejects the existing complex-channel magnitude as the sole depth fix.
    It says nothing about missing reactive-ion or parent-molecule channels.
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
