"""Gate tests for the mixed-layer feature-engine adapter."""

import numpy as np
import pytest

from petch.mixed_layer import MixedLayerParams
from petch.mixed_layer_mechanism import (
    MixedLayerMechanism,
    MixedLayerSurfaceState,
)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


def _ion(flux, energy, cosine=1.0):
    return EnergeticFlux(
        name="Ar+", flux_m2_s=flux,
        energy_eV=np.array([energy]),
        cosine_incidence=np.array([cosine]),
        weight=np.array([1.0]))


def _fluxes(shape=(), cfx=3.0e19, f=2.0e20, o=5.0e19, ion=6.0e18, e=1000.0):
    def field(v):
        return np.full(shape, v) if shape else v
    return SurfaceFluxes(
        neutral_flux_m2_s={"CFx": field(cfx), "F": field(f), "O": field(o)},
        energetic_fluxes=(_ion(field(ion), e),))


def test_advance_produces_positive_oxide_velocity():
    mech = MixedLayerMechanism(MixedLayerParams())
    state = mech.initial_state((4,))
    result = mech.advance(state, _fluxes((4,)), 5.0)
    assert result.etch_velocity_m_s.shape == (4,)
    assert np.all(result.etch_velocity_m_s > 0.0)
    assert np.all(np.asarray(result.state.removed_formula_units_m2) > 0.0)
    assert result.validity.within_declared_scope


def test_mask_advance_is_slower_than_oxide():
    oxide = MixedLayerMechanism(MixedLayerParams(substrate="sio2"))
    mask = MixedLayerMechanism(MixedLayerParams(substrate="carbon"))
    fluxes = _fluxes()
    v_oxide = oxide.advance(oxide.initial_state(()), fluxes, 5.0).etch_velocity_m_s
    v_mask = mask.advance(mask.initial_state(()), fluxes, 5.0).etch_velocity_m_s
    assert float(v_oxide) > 3.0 * float(v_mask)


def test_unmapped_species_refused():
    mech = MixedLayerMechanism(MixedLayerParams())
    fluxes = SurfaceFluxes(
        neutral_flux_m2_s={"CFx": 1e19, "F": 1e19, "O": 0.0, "SiF4": 1e18},
        energetic_fluxes=(_ion(1e18, 500.0),))
    validity = mech.validity(fluxes)
    assert not validity.within_declared_scope
    assert "SiF4" in validity.unsupported_neutral_species
    with pytest.raises(ValueError):
        mech.advance(mech.initial_state(()), fluxes, 1.0)


def test_remap_contract_roundtrip():
    state = MixedLayerSurfaceState(1e18, 2e18, 3e17, 4e16, 5e16, 6e17, 7e19)
    fields = state.conservative_surface_fields()
    rebuilt = state.with_conservative_surface_fields(fields)
    for name, value in fields.items():
        assert np.array_equal(getattr(rebuilt, name), value)
    with pytest.raises(ValueError):
        state.with_conservative_surface_fields({"n_f": 1.0})


def test_material_exchange_closes():
    mech = MixedLayerMechanism(MixedLayerParams())
    result = mech.advance(mech.initial_state(()), _fluxes(), 2.0)
    exchange = result.material_exchange
    assert not exchange.product_routing_complete  # declared unresolved, honestly
    for name in exchange.removed_units_m2:
        assert np.all(np.abs(exchange.residual_units_m2(name)) == 0.0)
    assert float(exchange.removed_units_m2["sio2_formula"]) > 0.0


def test_reaction_probabilities_shape_and_range():
    mech = MixedLayerMechanism(MixedLayerParams())
    state = mech.initial_state((3,))
    probability = mech.neutral_reaction_probability(state)
    assert set(probability) == {"CFx", "F", "O"}
    for values in probability.values():
        arr = np.asarray(values)
        assert arr.shape == (3,)
        assert np.all((arr >= 0.0) & (arr <= 1.0))


def test_spectrum_compression_uses_flux_weighted_mean():
    mech = MixedLayerMechanism(MixedLayerParams())
    two_pop = SurfaceFluxes(
        neutral_flux_m2_s={"CFx": 0.0, "F": 1e21, "O": 0.0},
        energetic_fluxes=(_ion(3e18, 400.0), _ion(1e18, 2000.0)))
    module = mech._module_fluxes(two_pop, ())
    assert float(module.ion_flux) == pytest.approx(4e18)
    assert float(module.ion_energy_eV) == pytest.approx(
        (3e18 * 400.0 + 1e18 * 2000.0) / 4e18)


def test_krueger_species_stoichiometry_carries_carbon_and_fluorine():
    from petch.mixed_layer_mechanism import (
        KRUEGER_2024_PRECURSOR_STOICHIOMETRY,
        build_krueger_2024_mixed_layer_mechanisms,
    )
    oxide, mask = build_krueger_2024_mixed_layer_mechanisms()
    assert oxide.parameters.substrate == "sio2"
    assert mask.parameters.substrate == "carbon"
    fluxes = SurfaceFluxes(
        neutral_flux_m2_s={"CF2": 1.0e19, "C2F3": 2.0e19, "O": 0.0,
                           "C3F4": 5.0e18},
        energetic_fluxes=(_ion(1e18, 800.0),))
    module = oxide._module_fluxes(fluxes, ())
    # carbon flux = 1*1e19 + 2*2e19; bound F = 2*1e19 + 3*2e19
    assert float(module.precursor_flux) == pytest.approx(5.0e19)
    assert float(module.precursor_fc_ratio) == pytest.approx(8.0e19 / 5.0e19)
    assert oxide.validity(fluxes).within_declared_scope  # C3F4 declared inert
    assert set(KRUEGER_2024_PRECURSOR_STOICHIOMETRY) == {
        "CF", "CF2", "CF3", "C2F3", "C2F4", "C3F5", "C3F6"}


def test_router_mixed_layer_option():
    from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
    from petch.mixed_layer_mechanism import MixedLayerMechanism

    router = build_krueger_2024_material_router_3d(surface_model="mixed_layer")
    mechanisms = getattr(router, "mechanisms", None) or getattr(
        router, "_mechanisms", None)
    if mechanisms is None:
        pytest.skip("router does not expose its mechanism map")
    assert all(isinstance(m, MixedLayerMechanism) for m in dict(mechanisms).values())
    with pytest.raises(ValueError):
        build_krueger_2024_material_router_3d(
            surface_model="mixed_layer", oxide_etch_yield_scale=0.5)


def test_duration_zero_is_identity():
    mech = MixedLayerMechanism(MixedLayerParams())
    state = mech.initial_state((2,))
    result = mech.advance(state, _fluxes((2,)), 0.0)
    assert result.state is state
    assert np.all(result.etch_velocity_m_s == 0.0)
