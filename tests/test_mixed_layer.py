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


def test_finite_range_interface_transport_stops_and_obeys_slant_path():
    from petch.ion_energy_deposition import csda_path_nm, FLUOROCARBON_FILM

    params = MixedLayerParams(film_energy_transport="csda_finite_range")
    incident = 200.0
    path_nm = csda_path_nm(
        incident,
        params.ion_atomic_number,
        params.ion_mass_amu,
        FLUOROCARBON_FILM,
    )
    identity = interface_energy_eV(incident, 0.0, params, 0.0)
    normal = interface_energy_eV(incident, 0.3 * path_nm, params, 1.0)
    slanted = interface_energy_eV(incident, 0.3 * path_nm, params, 0.5)
    stopped = interface_energy_eV(incident, 1.01 * path_nm, params, 1.0)
    assert identity == incident
    assert 0.0 <= slanted < normal < incident
    assert stopped == 0.0


def test_interface_transport_refuses_unknown_model():
    params = MixedLayerParams(film_energy_transport="not-a-physics-model")
    with pytest.raises(ValueError, match="film_energy_transport"):
        interface_energy_eV(200.0, 0.5, params)


def test_recession_is_supply_capacity_minimum():
    """Raising F supply saturates at the ion capacity (the ceiling).

    RESTATED 2026-08-06: fluxes scaled by 1/s0 (Gray Table 5-10, s0 = 0.02)
    so the sweep still spans starved -> saturated now that the supply term
    carries its measured sticking coefficient.  Assertions unchanged.
    """
    from petch.mixed_layer import _THERMAL_F_STICKING

    rates = []
    for j_f in (1.0e18 / _THERMAL_F_STICKING, 1.0e19 / _THERMAL_F_STICKING,
                1.0e20 / _THERMAL_F_STICKING, 1.0e21 / _THERMAL_F_STICKING,
                4.0e21 / _THERMAL_F_STICKING):
        fluxes = SurfaceFluxes(0.0, j_f, 0.0, 6.0e18, 1000.0)
        result = steady_state(fluxes)
        rates.append(result.sif4_rate)
    assert rates[1] > rates[0]
    assert rates[-1] == pytest.approx(rates[-2], rel=5e-2)  # capacity plateau
    # not supply-limited at the top (compare against the DELIVERED supply)
    assert rates[-1] < 0.25 * _THERMAL_F_STICKING * 4.0e21 / _THERMAL_F_STICKING


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
    # Flux points sit at the published-law clog boundary (film sputter
    # p0=0.9 @ 20 eV): p=2e20 clogs dry, oxygen rescues it, and more oxygen
    # thins the surviving film further.
    outcomes = []
    for j_o in (0.0, 2.0e21, 8.0e21):
        fluxes = SurfaceFluxes(2.0e20, 2.0e20, j_o, 6.0e18, 1000.0)
        result = steady_state(fluxes)
        outcomes.append((result.ledger_residuals.get("clogged", False),
                         float(np.asarray(result.state.film_thickness_nm())),
                         float(np.asarray(result.sif4_rate))))
    assert outcomes[0][0] and outcomes[0][2] == 0.0          # clogged, no etch
    assert not outcomes[1][0] and not outcomes[2][0]         # rescued
    assert outcomes[1][2] > 0.0 and outcomes[2][2] > 0.0
    assert outcomes[2][1] <= outcomes[1][1] + 1e-12          # film thins with O

    # Healthy etching regime: lattice oxygen already saturates the demand, so
    # gas oxygen buys nothing — the saturation point is C availability.
    # RESTATED 2026-08-06: the healthy-etching probe moves 3.0e19 -> 1.0e19
    # precursor with Gray's s0 = 0.02 (Table 5-10).  Delivered F drops 50x, so
    # the film wins at a lower precursor flux and the "healthy branch" the
    # assertion is about now sits below 3e19.  Assertion unchanged: on that
    # branch lattice oxygen already saturates demand and gas O buys nothing.
    lean = [steady_state(SurfaceFluxes(1.0e19, 2.0e20, j_o, 6.0e18, 1000.0)).sif4_rate
            for j_o in (0.0, 2.0e21)]
    assert lean[1] == pytest.approx(lean[0], rel=0.1)

    # Probe flux raised 1.0e21 -> 2.0e21 when the O channel was un-throttled:
    # Krueger's p_ox row is per collision with exposed polymer, so removing the
    # spurious composition factor strengthened O etching by 1/x_c = 2.69x and
    # moved the clog boundary up by the same factor in precursor flux. The
    # asserted physics is unchanged (a boundary exists; more O rescues it) —
    # only the probe point tracks the corrected channel. See
    # RESULTS_O_CHANNEL_2026-08-04.md and tests/test_o_channel_budget.py.
    def clogs(j_o):
        fluxes = SurfaceFluxes(2.0e21, 2.0e20, j_o, 6.0e18, 1000.0)
        return steady_state(fluxes).ledger_residuals.get("clogged", False)

    assert clogs(2.0e21) and not clogs(8.0e21)


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


def test_rate_follows_measured_sqrt_energy_law():
    """In the capacity-limited film-free regime the removal rate must follow
    the energy law of the PRIMARY SOURCE, not a reference-normalised shape.

    Gray (MIT thesis 1993) measured this system over 20-2000 eV and fitted
    beta_e = 0.053(sqrt(E) - sqrt(4)) (Eq. 5-35), having explicitly tested and
    rejected the Sigmund form that is linear in E (Fig. 5-2, p.161).  petch
    previously carried a ZBL-shaped factor that is linear to within 2 percent
    above 700 eV, i.e. the rejected shape wearing a different label; against
    Gray's six measured points it scores R2 = -1.84 where the sqrt law scores
    0.994 (RESEARCH_ENERGY_SCALING_2026-08-05).

    At F-saturated coverage the composite is dominated by the complex channel,
    so the rate ratio must track Gray's sqrt law directly.
    """
    from petch.mixed_layer import _complex_yield

    params = MixedLayerParams()
    rates = {}
    # F flux scaled by 1/s0 (Gray Table 5-10) to hold the same F-saturated
    # coverage regime the assertion was written for.  RESTATED 2026-08-06.
    from petch.mixed_layer import _THERMAL_F_STICKING
    for energy in (300.0, 1000.0, 3000.0):
        fluxes = SurfaceFluxes(0.0, 1.0e22 / _THERMAL_F_STICKING, 0.0,
                               6.0e18, energy)
        rates[energy] = steady_state(fluxes, params).sif4_rate
    for energy in (300.0, 3000.0):
        expected = float(_complex_yield(energy) / _complex_yield(1000.0))
        # The bare-sputter channel adds a small sqrt component of its own on
        # top; the composite tracks the complex law within 3 percent.
        assert rates[energy] / rates[1000.0] == pytest.approx(expected, rel=3e-2)


def test_reproduces_gray_measured_yield_points():
    """The complex channel must reproduce Gray's Table 5-10 measurements.

    beta_e = "the number of SiF4 molecules removed from fluorine saturated
    surface regions per incoming ion", measured at six energies.  Tolerance is
    Gray's own fit residual (Eq. 5-35 scores R2 = 0.9927 against these points);
    the 250 eV point sits below his own trend, so it carries the loosest bound.
    """
    from petch.mixed_layer import _complex_yield, _bare_sputter_yield

    measured = {20.0: 0.13, 150.0: 0.55, 250.0: 0.60,
                350.0: 0.85, 500.0: 1.10, 2000.0: 2.25}
    for energy, value in measured.items():
        tol = 0.25 if energy == 250.0 else 0.06
        assert float(_complex_yield(energy)) == pytest.approx(value, rel=tol)
    # Table 5-1, the physical channel, exactly.
    assert float(_bare_sputter_yield(350.0)) == pytest.approx(
        0.0139 * (350.0 ** 0.5 - 18.0 ** 0.5), rel=1e-12)


def test_rung0_degenerate_matches_langmuir_closed_form():
    """Rung 0 (design doc 5.5): with no carbon anywhere, the layer must
    reduce to the Belen/ViennaPS coverage structure exactly — steady state
    theta_F = s0*J/(s0*J + 4*capacity), rate = capacity * theta_F, with
    capacity from the measured energy law (nothing fitted).

    RESTATED 2026-08-06, not weakened: the adsorption coefficient s0 = 0.02
    is now carried explicitly (Gray, MIT thesis 1993, Table 5-10 and p.246 --
    the co-regressed partner of the beta_e law this same closed form uses).
    The Langmuir structure being asserted is identical; only the supply term
    now carries its measured sticking coefficient instead of an implicit 1.0.
    Probe fluxes are scaled by 1/s0 so the comparison stays in the same
    coverage regime it was written for."""
    from petch.mixed_layer import _THERMAL_F_STICKING, _complex_yield

    params = MixedLayerParams()
    for j_f, energy in ((5.0e19, 400.0), (2.0e20, 1000.0), (1.0e21, 2000.0)):
        j_f = j_f / _THERMAL_F_STICKING
        fluxes = SurfaceFluxes(0.0, j_f, 0.0, 6.0e18, energy)
        result = steady_state(fluxes, params)
        capacity = (params.volatilization_yield * fluxes.ion_flux
                    * float(_complex_yield(energy)))
        supply = _THERMAL_F_STICKING * j_f
        theta = supply / (supply + 4.0 * capacity)
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


def test_activated_site_balance_is_exact():
    """theta_act equilibrium must satisfy formation == activated-channel
    consumption + removal share (exact site bookkeeping, not the blend)."""
    fluxes = SurfaceFluxes(
        precursor_flux=0.0, fluorine_flux=0.0, oxygen_flux=0.0,
        ion_flux=6.0e18, ion_energy_eV=1000.0,
        chemisorption_carbon_flux=5.0e19, chemisorption_fluorine_flux=1.0e20,
        chemisorption_activated_carbon_flux=1.5e20,
        chemisorption_activated_fluorine_flux=3.0e20)
    result = steady_state(fluxes)
    state = result.state
    theta = float(np.asarray(state.n_act)) / 1.0e19
    assert 0.0 < theta < 1.0
    # Reconstruct the code's own site balance from the converged state:
    # formation on open oxide == activated-channel consumption + removal.
    from petch.mixed_layer import _MONOLAYER_AREAL_M2
    theta_f_layer = min(float(np.asarray(state.n_f)) / _MONOLAYER_AREAL_M2, 1.0)
    site_open = 1.0 - theta_f_layer          # film-free case
    formation = 0.9 * 6.0e18 * (1.0 - theta)
    consumption = (theta * 1.5e20 * site_open
                   + float(np.asarray(result.substrate_removal_rate)) * theta)
    assert formation == pytest.approx(consumption, rel=0.15)
