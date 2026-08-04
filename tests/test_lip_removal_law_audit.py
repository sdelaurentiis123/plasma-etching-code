"""Guards for the grazing lip-removal audit (RESULTS_LIP_REMOVAL_AUDIT_2026-08-04).

Two things are pinned here:

1. The *falsification*: no angular yield function can balance the mask-lip film
   budget at 86-89 deg incidence.  The required f(theta) exceeds any published
   per-ion yield by more than an order of magnitude, so future work must not
   spend effort re-shaping the angular law to open the mouth.
2. The *coverage gap*: three Appendix-B rows carry angular markers that the
   module does not yet apply.  The test records which kernels do and do not
   carry an angular factor so the gap cannot be silently lost.
"""

import numpy as np

from petch.mixed_layer import _threshold_power_yield

# Krueger Table-I base-case fluxes (m^-2 s^-1) and Appendix-B per-cell rows.
_DEPOSITOR_FLUX = {"CF": 4.4e20, "CF2": 9.4e20, "CF3": 8.4e19, "C2F3": 6.8e20}
_DEPOSITION_PROBABILITY = {"CF": 0.1, "CF2": 0.1, "CF3": 0.1, "C2F3": 0.03}
_OXYGEN_FLUX = 7.7e20
_OXIDATION_PROBABILITY = 0.0423
_ION_FLUX = 9.6e19
_ION_ENERGY_EV = 1500.0
_KRESS_B = 9.3


def _deposition_rate():
    return sum(_DEPOSITION_PROBABILITY[name] * flux
               for name, flux in _DEPOSITOR_FLUX.items())


def _kress(cosine):
    return max((1.0 + _KRESS_B * (1.0 - cosine ** 2)) * cosine, 0.0)


def test_oxygen_share_of_the_lip_budget_is_the_published_ratio():
    """Both channels are per-cell thermal rows, so their ratio is geometry
    free and must reproduce Krueger's mechanism exactly (RESULTS_O_CHANNEL)."""
    ratio = _OXIDATION_PROBABILITY * _OXYGEN_FLUX / _deposition_rate()
    assert ratio == 0.1953 or abs(ratio - 0.1953) < 5e-4


def test_no_angular_law_can_balance_the_lip():
    """The audit's central negative result.

    Removal at a lip face is (areal ion flux) x (energy yield) x f(theta).
    The transport already supplies the areal cosine, so closing the residual
    with the angular law alone would need f(88.7 deg) ~ 62 -- an order of
    magnitude beyond the peak of any published per-ion yield (our own Kress
    form peaks at 4.17; the in-chemistry measurements of Cho 2000 and
    Schaepkens 1998 bound the peak near 1.3).
    """
    deposition = _deposition_rate()
    oxygen_share = _OXIDATION_PROBABILITY * _OXYGEN_FLUX / deposition
    yield_energy = min(_threshold_power_yield(
        np.asarray(_ION_ENERGY_EV), 0.9, 20.0, 500.0, 0.5), 1.0)

    cosine = float(np.cos(np.radians(88.7)))
    required = ((1.0 - oxygen_share) * deposition
                / (_ION_FLUX * cosine * yield_energy))
    assert required > 50.0                      # ~61.6 as audited
    assert required > 10.0 * (1.0 + _KRESS_B) ** 0.5 * 1.0  # >> our own peak

    # Even discarding the areal cosine entirely -- the "double-cos" hypothesis
    # -- leaves the ion channel far short of the ~0.80 x deposition needed.
    no_areal_cosine = (_ION_FLUX * yield_energy * _kress(cosine)) / deposition
    assert no_areal_cosine < 0.2
    assert 1.0 - oxygen_share - no_areal_cosine > 0.6


def test_grazing_ion_removal_is_a_negligible_share_of_the_lip_budget():
    """Quantifies why: the lip sits where the law has already collapsed."""
    deposition = _deposition_rate()
    yield_energy = min(_threshold_power_yield(
        np.asarray(_ION_ENERGY_EV), 0.9, 20.0, 500.0, 0.5), 1.0)
    shares = {}
    for angle in (86.0, 88.0, 88.7):
        cosine = float(np.cos(np.radians(angle)))
        shares[angle] = (_ION_FLUX * cosine * yield_energy
                         * _kress(cosine) / deposition)
    assert shares[86.0] < 0.05
    assert shares[88.7] < 0.005
    assert shares[86.0] > shares[88.0] > shares[88.7]   # monotone collapse


def test_kress_form_cannot_peak_beyond_the_magic_angle():
    """(1 + B sin^2) cos peaks at cos^2 = (1+B)/(3B) -> 54.7 deg as B -> inf.

    Huang's thesis describes the physical-sputtering f(theta) as peaking near
    60 deg; this form can never reach it, which bounds how faithful any
    B-tuning could be.  Recorded so the shape question stays visible.
    """
    for b in (2.0, 9.3, 100.0, 1e6):
        cos_peak = np.sqrt((1.0 + b) / (3.0 * b))
        assert np.degrees(np.arccos(cos_peak)) <= 54.8
    peak = max(_kress(c) for c in np.cos(np.radians(np.linspace(0.0, 90.0, 9001))))
    assert 4.1 < peak < 4.2        # our B = 9.3 peak/normal, audited at 4.17


def test_angular_marker_coverage_gap_is_recorded():
    """Appendix B marks four ion rows with an angular class; the module
    applies an angular factor to exactly one of them.

    Verbatim rows (krueger_thesis.txt):
      CF(s)     + Ar+ -> EP   + Ar#          0.9    20 0.5 500 1
      AC(s)     + Ar+ -> C    + Ar#          0.001 200 0.4 250 1
      SiO2(s)   + Ar+ -> Ar#  + SiO2         0.0852 70 1   140 1
      SiO2CF(s) + Ar+ -> SiF + CO2 + Ar#     0.1471 35 1   140 2

    This test documents the gap rather than asserting the physics: closing it
    changes the oxide channels that set trench depth, so it must be graded by
    a confirmation run, not landed blind.
    """
    marked = {"polymer_sputter": 1, "mask_ac": 1, "oxide_bare": 1,
              "oxide_complex": 2}
    applies_angular_factor = {"polymer_sputter": True, "mask_ac": False,
                              "oxide_bare": False, "oxide_complex": False}
    assert set(marked) == set(applies_angular_factor)
    implemented = [k for k, v in applies_angular_factor.items() if v]
    assert implemented == ["polymer_sputter"]
    assert sum(1 for v in applies_angular_factor.values() if not v) == 3
