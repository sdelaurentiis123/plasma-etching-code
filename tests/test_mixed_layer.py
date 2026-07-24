"""Gate tests for the element-resolved mixed-layer v1 (standalone module)."""

import math

import numpy as np
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
    # Higher ion energy keeps a flux etching that clogs at lower energy. The
    # flux sits between the two energies' total removal capacities — at much
    # higher precursor flux deposition beats removal at ANY energy (that
    # regime spuriously passed before the accelerated integrator existed).
    boundary_low = clogs(3.0e19, 200.0)
    boundary_high = clogs(3.0e19, 3000.0)
    assert boundary_low and not boundary_high


def test_oxygen_thins_film_and_moves_clog_boundary():
    """Gas oxygen's role: thin the film monotonically, and rescue a flux
    combination that clogs without it (the Krueger O2-sweep mechanism).
    In a healthy etching state the demand is already saturated by lattice
    oxygen, so extra gas O buys almost nothing — C availability, not a
    fitted constant, sets the saturation point."""
    # The film balance is a cliff (removal saturates with coverage, deposition
    # does not): without oxygen this flux clogs; enough oxygen rescues it, and
    # within the etching branch the residual thin film keeps thinning with O.
    outcomes = []
    for j_o in (0.0, 5.0e20, 2.0e21):
        fluxes = SurfaceFluxes(6.0e19, 2.0e20, j_o, 6.0e18, 1000.0)
        result = steady_state(fluxes)
        outcomes.append((result.ledger_residuals.get("clogged", False),
                         float(np.asarray(result.state.film_thickness_nm())),
                         float(np.asarray(result.sif4_rate))))
    assert outcomes[0][0] and outcomes[0][2] == 0.0          # clogged, no etch
    assert not outcomes[1][0] and not outcomes[2][0]         # rescued
    assert outcomes[1][2] > 0.0 and outcomes[2][2] > 0.0
    assert outcomes[2][1] < outcomes[1][1]                   # residual film thins

    # Healthy etching regime: lattice oxygen already saturates the demand, so
    # gas oxygen buys nothing — the saturation point is C availability.
    lean = [steady_state(SurfaceFluxes(3.0e19, 2.0e20, j_o, 6.0e18, 1000.0)).sif4_rate
            for j_o in (0.0, 2.0e21)]
    assert lean[1] == pytest.approx(lean[0], rel=0.1)

    def clogs(j_o):
        fluxes = SurfaceFluxes(2.0e20, 5.0e19, j_o, 3.0e18, 600.0)
        return steady_state(fluxes).ledger_residuals.get("clogged", False)

    assert clogs(0.0) and not clogs(4.0e21)


def test_selectivity_emerges_from_lattice_oxygen():
    """Same plasma over SiO2 vs an a-C mask: the mask must grow the thicker
    film, see less interface energy, and etch much slower — with no
    selectivity parameter anywhere (Standaert/Oehrlein mechanism)."""
    oxide = steady_state(_BASE, MixedLayerParams(substrate="sio2"))
    mask = steady_state(_BASE, MixedLayerParams(substrate="carbon"))
    assert mask.state.film_thickness_nm() > 2.0 * oxide.state.film_thickness_nm()
    assert mask.interface_energy_eV < oxide.interface_energy_eV
    assert oxide.recession_velocity_m_s > 3.0 * mask.recession_velocity_m_s


def test_carbon_substrate_ledger_closes():
    params = MixedLayerParams(substrate="carbon")
    _, worst = _integrate(_BASE, n_steps=1000, params=params)
    assert worst < 1e-9


def test_rate_follows_derived_energy_law():
    """In the capacity-limited film-free regime, the substrate removal rate
    must scale exactly as the ZBL deposited-in-layer factor — the same law
    the K24-DEKNOB-1 study validated. No fitted energy dependence."""
    from petch.mixed_layer import _deposited_energy

    params = MixedLayerParams()
    rates = {}
    for energy in (300.0, 1000.0, 3000.0):
        fluxes = SurfaceFluxes(0.0, 1.0e22, 0.0, 6.0e18, energy)
        rates[energy] = steady_state(fluxes, params).sif4_rate
    for energy in (300.0, 3000.0):
        expected = (_deposited_energy(energy, 1.0, params)[0]
                    / _deposited_energy(1000.0, 1.0, params)[0])
        assert rates[energy] / rates[1000.0] == pytest.approx(expected, rel=5e-3)


def test_rung0_degenerate_matches_langmuir_closed_form():
    """Rung 0 (design doc 5.5): with no carbon anywhere, the layer must
    reduce to the Belen/ViennaPS coverage structure exactly — steady state
    theta_F = J/(J + 4*capacity), rate = capacity * theta_F, with capacity
    from the derived deposited-energy factor (nothing fitted)."""
    from petch.mixed_layer import _deposited_energy

    params = MixedLayerParams()
    for j_f, energy in ((5.0e19, 400.0), (2.0e20, 1000.0), (1.0e21, 2000.0)):
        fluxes = SurfaceFluxes(0.0, j_f, 0.0, 6.0e18, energy)
        result = steady_state(fluxes, params)
        eps_dep = _deposited_energy(energy, 1.0, params)[0]
        capacity = (params.volatilization_yield * fluxes.ion_flux
                    * eps_dep / params.reference_energy_eV)
        theta = j_f / (j_f + 4.0 * capacity)
        assert result.sif4_rate == pytest.approx(capacity * theta, rel=1e-5)


def test_vectorized_step_matches_scalar_bitwise():
    """One array call over N faces must equal N scalar calls exactly — the
    contract that lets the feature engine (and later the GPU port) batch
    per-face chemistry without changing a single bit."""
    import numpy as np

    conditions = [
        (3.0e19, 2.0e20, 5.0e19, 6.0e18, 1000.0, 1.0),
        (1.5e20, 5.0e19, 0.0, 3.0e18, 600.0, 0.7),
        (0.0, 1.0e21, 2.0e20, 1.0e19, 2000.0, 0.3),
        (8.0e19, 0.0, 1.0e20, 2.0e18, 150.0, 1.0),
    ]
    params = MixedLayerParams()
    dt = 1e-4
    # Scalar path, advanced two steps from bare.
    scalar_states = []
    for cond in conditions:
        state = MixedLayerState()
        for _ in range(2):
            state = step(state, SurfaceFluxes(*cond), dt, params).state
        scalar_states.append(state)
    # Vector path: all faces at once.
    arrays = [np.array(col) for col in zip(*conditions)]
    vec_fluxes = SurfaceFluxes(*arrays)
    vec_state = MixedLayerState(*[np.zeros(len(conditions))] * 6)
    for _ in range(2):
        vec_state = step(vec_state, vec_fluxes, dt, params).state
    for i, expected in enumerate(scalar_states):
        for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f"):
            assert float(np.asarray(getattr(vec_state, name))[i]) == float(
                getattr(expected, name)), (i, name)


def test_bitwise_determinism():
    a, _ = _integrate(_BASE, n_steps=500)
    b, _ = _integrate(_BASE, n_steps=500)
    assert a.state == b.state
    assert a.sif4_rate == b.sif4_rate
