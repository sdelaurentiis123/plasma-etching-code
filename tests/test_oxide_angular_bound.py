"""The oxide/mask class-1 angular shape is bounded by in-chemistry measurement.

Krueger's cited class-1 source (Kress 1999) is a Cu/Ar molecular-dynamics study
whose shape gives peak/normal 4.17.  The only angular sputter measurements on
SiO2 in fluorocarbon -- Cho 2000 (JVST A 18, 2705) and Schaepkens 1998 (JVST A
16, 3281) -- bound peak/normal at 1.30-1.33.  The oxide/mask rows use the
measurement-bounded shape; the polymer row keeps Krueger's cited value.
"""
import numpy as np
import pytest

from petch.mixed_layer import (_angular_oxide_sputter, _angular_physical_sputter,
                               _KRUEGER_CLASS1_B, _OXIDE_CLASS1_B)

_TH = np.linspace(0.0, 90.0, 901)
_COS = np.cos(np.deg2rad(_TH))


def test_normal_incidence_is_unity_for_both_shapes():
    """f(0) = 1 for any B -- every blanket/normal-incidence result is unchanged."""
    assert float(_angular_oxide_sputter(1.0)) == 1.0
    assert float(_angular_physical_sputter(1.0)) == 1.0


def test_oxide_peak_sits_inside_the_measured_band():
    peak = float(np.max(_angular_oxide_sputter(_COS)))
    assert 1.30 <= peak <= 1.35, peak


def test_polymer_row_keeps_kruegers_cited_value():
    assert _KRUEGER_CLASS1_B == 9.3
    assert float(np.max(_angular_physical_sputter(_COS))) == pytest.approx(4.17, abs=0.02)


def test_oxide_shape_is_3x_below_the_cu_md_form_at_its_peak():
    """The quantity the swap is about: 4.17 -> 1.31 peak/normal ratio."""
    assert (float(np.max(_angular_physical_sputter(_COS)))
            / float(np.max(_angular_oxide_sputter(_COS)))) == pytest.approx(3.2, abs=0.15)


def test_both_shapes_vanish_at_grazing_and_are_non_negative():
    for f in (_angular_oxide_sputter, _angular_physical_sputter):
        vals = f(_COS)
        assert np.all(vals >= 0.0)
        assert float(f(0.0)) == 0.0


def test_oxide_shape_bounds_the_timestep_amplification():
    """The ml21 dt collapse: removal amplification on ~50 deg faces.

    With B=9.3 the oxide rows amplify off-normal removal 2.5x (measured dt
    collapse 2.63x); the measured bound keeps it under 1.2x.
    """
    amp_kruger = float(np.max(_angular_physical_sputter(_COS)))
    amp_oxide = float(np.max(_angular_oxide_sputter(_COS)))
    assert amp_kruger > 2.5
    assert amp_oxide < 1.4
