import numpy as np
import pytest

from petch.surface_kinetics import EnergeticYield


def test_yield_knee_reduces_to_legacy_below_and_saturates_above():
    legacy = EnergeticYield(
        reference_yield=1.0, threshold_energy_eV=60.0, reference_energy_eV=1000.0,
        energy_exponent=0.5)
    kneed = EnergeticYield(
        reference_yield=1.0, threshold_energy_eV=60.0, reference_energy_eV=1000.0,
        energy_exponent=0.5, knee_energy_eV=3000.0)
    low = np.array([100.0, 500.0, 2999.0])
    assert np.allclose(kneed.evaluate(low, 1.0), legacy.evaluate(low, 1.0))
    high = np.array([3000.0, 3600.0, 8000.0])
    values = kneed.evaluate(high, 1.0)
    assert np.allclose(values, values[0])
    assert values[0] < legacy.evaluate(np.array([3600.0]), 1.0)[0]


def test_yield_knee_must_exceed_threshold():
    with pytest.raises(ValueError, match="knee"):
        EnergeticYield(
            reference_yield=1.0, threshold_energy_eV=60.0,
            reference_energy_eV=1000.0, knee_energy_eV=50.0)


def test_oxygen_langmuir_form_limits():
    # rate form phi/(1 + phi/phi_half): linear at low flux, asymptote phi_half at high
    phi_half = 4.0e20
    low = 1.0e18
    assert low / (1 + low / phi_half) == pytest.approx(low, rel=3e-3)
    high = 4.0e22
    assert high / (1 + high / phi_half) == pytest.approx(phi_half, rel=2e-2)
