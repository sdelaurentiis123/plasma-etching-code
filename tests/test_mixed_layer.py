"""Gate tests for the element-resolved mixed-layer v1 (standalone module)."""

import math

import pytest

from petch.mixed_layer import (
    MixedLayerParams,
    MixedLayerState,
    SurfaceFluxes,
    interface_energy_eV,
    steady_state,
    step,
)

_BASE = SurfaceFluxes(
    precursor_flux=3.0e19,
    fluorine_flux=2.0e20,
    oxygen_flux=5.0e19,
    ion_flux=6.0e18,
    ion_energy_eV=1000.0,
)


def _integrate(fluxes, n_steps=2000, dt=1e-4, params=MixedLayerParams()):
    state = MixedLayerState()
    worst = 0.0
    result = None
    for _ in range(n_steps):
        result = step(state, fluxes, dt, params)
        for key in ("fluorine", "carbon", "silicon", "oxygen"):
            worst = max(worst, abs(result.ledger_residuals[key]))
        state = result.state
    return result, worst


def test_ledger_closes_to_machine_precision():
    _, worst = _integrate(_BASE)
    assert worst < 1e-9


def test_ledger_closes_under_starved_and_rich_conditions():
    for flux in (
        SurfaceFluxes(1e20, 1e19, 0.0, 1e18, 500.0),      # polymer rich
        SurfaceFluxes(1e17, 5e20, 2e20, 2e19, 2000.0),    # ion/F rich
        SurfaceFluxes(0.0, 0.0, 0.0, 1e19, 1500.0),       # ions only
    ):
        _, worst = _integrate(flux, n_steps=500)
        assert worst < 1e-9, flux


def test_degenerate_no_precursor_keeps_film_and_carbon_empty():
    fluxes = SurfaceFluxes(0.0, 2.0e20, 0.0, 6.0e18, 1000.0)
    result, _ = _integrate(fluxes, n_steps=1000)
    assert result.state.n_c_film == 0.0
    assert result.state.n_c == 0.0
    assert result.sif4_rate > 0.0


def test_interface_energy_monotone_in_film_thickness():
    params = MixedLayerParams()
    energies = [interface_energy_eV(1000.0, d, params) for d in (0.0, 1.0, 3.0, 8.0)]
    assert energies[0] == pytest.approx(1000.0)
    assert all(a > b for a, b in zip(energies, energies[1:]))


def test_recession_is_supply_capacity_minimum():
    """Raising F supply saturates at the ion capacity (the ceiling)."""
    rates = []
    for j_f in (1.0e18, 1.0e19, 1.0e20, 1.0e21, 4.0e21):
        fluxes = SurfaceFluxes(0.0, j_f, 0.0, 6.0e18, 1000.0)
        result = steady_state(fluxes)
        rates.append(result.sif4_rate)
    assert rates[1] > rates[0]
    assert rates[-1] == pytest.approx(rates[-2], rel=5e-2)  # capacity plateau
    assert rates[-1] < 0.25 * 4.0e21  # not supply-limited at the top


def test_clog_boundary_emerges_and_moves_with_ion_energy():
    def clogs(precursor, energy):
        fluxes = SurfaceFluxes(precursor, 5.0e19, 0.0, 3.0e18, energy)
        return steady_state(fluxes).ledger_residuals.get("clogged", False)

    assert not clogs(1.0e18, 1000.0)
    assert clogs(3.0e21, 1000.0)
    # Higher ion energy keeps a flux etching that clogs at lower energy.
    boundary_low = clogs(2.0e20, 200.0)
    boundary_high = clogs(2.0e20, 3000.0)
    assert boundary_low and not boundary_high


def test_oxygen_thins_film_and_rate_saturates():
    thicknesses = []
    rates = []
    for j_o in (0.0, 5.0e19, 5.0e20, 2.0e21, 4.0e21):
        fluxes = SurfaceFluxes(3.0e19, 2.0e20, j_o, 6.0e18, 1000.0)
        result = steady_state(fluxes)
        thicknesses.append(result.state.film_thickness_nm())
        rates.append(result.sif4_rate)
    assert all(a >= b for a, b in zip(thicknesses, thicknesses[1:]))
    assert rates[-2] > rates[0]  # oxygen enables etching through the film
    # Once the film is consumed, additional oxygen buys almost nothing:
    # the saturation point is C availability, not a fitted constant.
    assert thicknesses[-1] < 0.1
    assert rates[-1] == pytest.approx(rates[-2], rel=0.1)


def test_bitwise_determinism():
    a, _ = _integrate(_BASE, n_steps=500)
    b, _ = _integrate(_BASE, n_steps=500)
    assert a.state == b.state
    assert a.sif4_rate == b.sif4_rate
