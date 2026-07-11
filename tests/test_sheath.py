import numpy as np

from petch.sheath import CollisionlessRFSheath, child_langmuir_sheath_thickness


def test_child_langmuir_thickness_is_positive_and_density_scaling_is_exact():
    s1 = child_langmuir_sheath_thickness(100.0, 4.0, 40.0, 1e16)
    s4 = child_langmuir_sheath_thickness(100.0, 4.0, 40.0, 4e16)
    assert s1 > 0.0
    assert np.isclose(s4 / s1, 0.5, rtol=1e-12)


def test_static_sheath_recovers_bohm_plus_voltage_energy():
    sheath = CollisionlessRFSheath(
        V_dc=80.0, V_rf=0.0, frequency_hz=1e6, Te_eV=4.0,
        ion_mass_amu=40.0, thickness_m=1e-3)
    energies = sheath.ion_impact_energies(np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False))
    assert np.allclose(energies, 82.0, atol=0.08)


def test_finite_transit_time_phase_mixes_high_frequency_iedf():
    phases = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)
    common = dict(V_dc=80.0, V_rf=40.0, Te_eV=4.0, ion_mass_amu=40.0, thickness_m=1e-3)
    low = CollisionlessRFSheath(frequency_hz=1e5, **common).ion_impact_energies(phases)
    high = CollisionlessRFSheath(frequency_hz=2e7, **common).ion_impact_energies(phases)
    assert np.ptp(high) < 0.25 * np.ptp(low)
    assert abs(high.mean() - 82.0) < 1.0
