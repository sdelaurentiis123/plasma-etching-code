"""Gates for atom-resolved (per-event) ion chemistry in the mixed layer."""

import numpy as np
import pytest

from petch.mixed_layer import (
    MixedLayerParams,
    MixedLayerState,
    SurfaceFluxes,
    step,
)


def _scalar_fluxes(**kw):
    base = dict(precursor_flux=3.0e19, fluorine_flux=2.0e20,
                oxygen_flux=5.0e19, ion_flux=6.0e18, ion_energy_eV=1000.0,
                cosine_incidence=0.8)
    base.update(kw)
    return SurfaceFluxes(**base)


def _atom_fluxes(faces, fluxes_list, energies, cosines, n_faces=1, **kw):
    base = dict(
        precursor_flux=np.full(n_faces, 3.0e19),
        fluorine_flux=np.full(n_faces, 2.0e20),
        oxygen_flux=np.full(n_faces, 5.0e19),
        ion_flux=np.zeros(n_faces),
        ion_energy_eV=np.zeros(n_faces),
        cosine_incidence=np.ones(n_faces),
        ion_atom_face=np.asarray(faces, dtype=int),
        ion_atom_flux=np.asarray(fluxes_list, dtype=float),
        ion_atom_energy_eV=np.asarray(energies, dtype=float),
        ion_atom_cosine=np.asarray(cosines, dtype=float))
    for key, value in kw.items():
        base[key] = value
    return SurfaceFluxes(**base)


def _states(n=1):
    zero = np.zeros(n)
    return MixedLayerState(zero, zero, zero, zero, zero, zero, zero)


def test_single_atom_matches_scalar_path_exactly():
    """One atom with the scalar (flux, E, cos) must reproduce the scalar
    path to machine precision at every state field."""
    params = MixedLayerParams()
    state_s = MixedLayerState()
    state_a = _states(1)
    dt = 1e-4
    for _ in range(200):
        rs = step(state_s, _scalar_fluxes(), dt, params)
        ra = step(state_a, _atom_fluxes([0], [6.0e18], [1000.0], [0.8]),
                  dt, params)
        state_s, state_a = rs.state, ra.state
        for name in ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f",
                     "n_xl_film"):
            scalar_value = float(np.asarray(getattr(state_s, name)))
            atom_value = float(np.asarray(getattr(state_a, name))[0])
            assert atom_value == pytest.approx(scalar_value, rel=1e-12, abs=1e-6)


def test_two_atoms_differ_from_compressed_mean_with_concave_sign():
    """Splitting the same total flux into a low/high energy pair must give a
    LOWER sputter-driven removal than the flux-weighted mean energy (the
    published sqrt-E film law is concave: Jensen), proving atoms bypass the
    compression."""
    params = MixedLayerParams()
    dt = 1e-4
    total = 6.0e18
    low, high = 200.0, 1800.0
    mean = 0.5 * (low + high)
    state_mean = _states(1)
    state_pair = _states(1)
    for _ in range(400):
        r_mean = step(state_mean, _atom_fluxes(
            [0], [total], [mean], [1.0]), dt, params)
        r_pair = step(state_pair, _atom_fluxes(
            [0, 0], [0.5 * total, 0.5 * total], [low, high], [1.0, 1.0]),
            dt, params)
        state_mean, state_pair = r_mean.state, r_pair.state
    film_mean = float(np.asarray(state_mean.n_c_film + state_mean.n_f_film))
    film_pair = float(np.asarray(state_pair.n_c_film + state_pair.n_f_film))
    # Concave sputter law: the pair sputters LESS than the mean -> thicker film.
    assert film_pair > film_mean * 1.001


def test_atom_ledgers_close():
    params = MixedLayerParams()
    state = _states(2)
    fluxes = _atom_fluxes(
        [0, 0, 1], [3.0e18, 3.0e18, 6.0e18], [400.0, 1600.0, 1000.0],
        [0.9, 0.3, 1.0], n_faces=2)
    worst = 0.0
    dt = 1e-4
    for _ in range(500):
        result = step(state, fluxes, dt, params)
        for key in ("fluorine", "carbon", "silicon", "oxygen"):
            worst = max(worst, float(np.max(np.abs(
                np.asarray(result.ledger_residuals[key])))))
        state = result.state
    assert worst < 1e-9


def test_adapter_builds_atoms_and_drops_compression_omission():
    from petch.mixed_layer_mechanism import (
        build_krueger_2024_mixed_layer_mechanisms,
    )
    from petch.surface_kinetics import (
        EnergeticFlux,
        FaceResolvedEnergeticFlux,
        SurfaceFluxes as EngineFluxes,
    )

    oxide, _ = build_krueger_2024_mixed_layer_mechanisms()
    events = FaceResolvedEnergeticFlux(
        "Ar+", 3, np.array([0, 1, 1]), np.array([1.0e19, 5.0e18, 5.0e18]),
        np.array([1500.0, 400.0, 1600.0]), np.array([1.0, 0.2, 0.9]))
    hot = FaceResolvedEnergeticFlux(
        "Ar+:hot_neutral", 3, np.array([2]), np.array([2.0e18]),
        np.array([1350.0]), np.array([0.6]))
    fluxes = EngineFluxes(
        neutral_flux_m2_s={"CF2": np.full(3, 1e19), "O": np.zeros(3),
                           "C3F4": np.zeros(3)},
        energetic_fluxes=(events, hot))
    module = oxide._module_fluxes(fluxes, (3,))
    assert module.ion_atom_face is not None
    assert len(module.ion_atom_face) == 4  # 3 primary + 1 hot neutral
    assert float(np.asarray(module.ion_atom_energy_eV)[-1]) == 1350.0
    omissions = oxide.validity(fluxes).known_model_form_omissions
    assert not any("compressed" in item for item in omissions)
    result = oxide.advance(oxide.initial_state((3,)), fluxes, 0.5)
    assert result.validity.within_declared_scope
    assert np.all(np.asarray(result.etch_velocity_m_s) >= 0.0)
