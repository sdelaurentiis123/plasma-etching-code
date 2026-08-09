import numpy as np
import pytest

from petch.chang_sawin_chlorine_sio2 import (
    ChangSawinArClSiO2Mechanism,
    chang_sawin_sio2_angular_factor,
)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


def _flux(energy, ratio, angle_deg=0.0, ion_flux=2.0e18):
    return SurfaceFluxes(
        {"Cl": ratio * ion_flux},
        (
            EnergeticFlux(
                "Ar+",
                ion_flux,
                np.asarray([energy]),
                np.asarray([np.cos(np.deg2rad(angle_deg))]),
                np.asarray([1.0]),
            ),
        ),
    )


@pytest.mark.parametrize(
    "energy,y0,beta,s",
    [(70.0, 0.01, 0.04, 0.001), (100.0, 0.02, 0.08, 0.005)],
)
def test_replays_table_4_2_cards(energy, y0, beta, s):
    mechanism = ChangSawinArClSiO2Mechanism()
    ratio = 90.0
    result = mechanism.advance(
        mechanism.initial_state(), _flux(energy, ratio), 1.0
    )
    theta = s * ratio / (s * ratio + beta)
    expected = y0 + (beta - y0) * theta
    assert float(result.chlorination_fraction) == pytest.approx(theta)
    assert float(result.mean_yield_sio2_formula_per_ion) == pytest.approx(expected)
    assert result.validity.within_declared_scope


def test_zero_chlorine_recovers_physical_yield():
    mechanism = ChangSawinArClSiO2Mechanism()
    result = mechanism.advance(
        mechanism.initial_state(), _flux(100.0, 0.0), 2.5
    )
    assert float(result.chlorination_fraction) == 0.0
    assert float(result.mean_yield_sio2_formula_per_ion) == pytest.approx(0.02)
    assert float(result.removed_sio2_formula_units_m2) == pytest.approx(
        2.0e18 * 0.02 * 2.5
    )
    assert not result.material_exchange.product_routing_complete


def test_measured_oxide_angle_peaks_near_sixty_degrees():
    angles = np.asarray([0.0, 30.0, 60.0, 75.0, 90.0])
    factor = chang_sawin_sio2_angular_factor(
        np.cos(np.deg2rad(angles))
    )
    assert factor[0] == pytest.approx(1.0)
    assert factor[2] > factor[1]
    assert factor[2] > factor[3]
    assert factor[-1] == pytest.approx(0.0, abs=1.0e-12)


def test_energy_extrapolation_is_explicitly_evidence_gated():
    mechanism = ChangSawinArClSiO2Mechanism()
    fluxes = _flux(50.0, 90.0)
    with pytest.raises(ValueError, match="outside declared scope"):
        mechanism.advance(mechanism.initial_state(), fluxes, 1.0)
    result = mechanism.advance(
        mechanism.initial_state(), fluxes, 1.0, strict=False
    )
    assert not result.validity.within_declared_scope
    assert "70--100 eV" in result.validity.reasons[0]


def test_absolute_velocity_uses_formula_unit_density():
    mechanism = ChangSawinArClSiO2Mechanism()
    result = mechanism.advance(
        mechanism.initial_state(), _flux(100.0, 90.0), 1.0
    )
    expected = (
        float(result.removal_rate_sio2_formula_m2_s)
        / mechanism.parameters.bulk_sio2_formula_density_m3
    )
    assert float(result.etch_velocity_m_s) == pytest.approx(expected)
