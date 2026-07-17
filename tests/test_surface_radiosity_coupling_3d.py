from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pytest

from petch.material_mechanism_3d import (
    MaterialMechanismRouter3D,
    MaterialSurfaceState3D,
)
from petch.neutral_radiosity_3d import (
    DiffuseFormFactors3D,
    solve_diffuse_neutral_radiosity_3d,
)
from petch.surface_exchange import SurfaceMaterialExchange, SurfaceProductPopulation
from petch.surface_kinetics import MechanismValidity, SurfaceFluxes
from petch.surface_radiosity_coupling_3d import (
    SurfaceRadiosityCacheIdentityError,
    SurfaceRadiosityCouplingRefusal,
    SurfaceRadiosityOperatorCache3D,
    integrate_surface_radiosity_chemistry_3d,
)


VALID = MechanismValidity(True, (), (), (), True, ())


@dataclass(frozen=True)
class _ScalarSurfaceState:
    coverage: np.ndarray
    integrated_incident_units_m2: np.ndarray

    def __post_init__(self):
        coverage, integrated = np.broadcast_arrays(
            np.asarray(self.coverage, dtype=float),
            np.asarray(self.integrated_incident_units_m2, dtype=float))
        coverage = np.array(coverage, copy=True)
        integrated = np.array(integrated, copy=True)
        if (np.any(~np.isfinite(coverage)) or np.any(coverage < 0.0)
                or np.any(coverage > 1.0 + 32.0 * np.finfo(float).eps)
                or np.any(~np.isfinite(integrated)) or np.any(integrated < 0.0)):
            raise ValueError("invalid manufactured scalar state")
        coverage = np.minimum(coverage, 1.0)
        coverage.setflags(write=False)
        integrated.setflags(write=False)
        object.__setattr__(self, "coverage", coverage)
        object.__setattr__(self, "integrated_incident_units_m2", integrated)

    def conservative_surface_fields(self):
        return {
            "coverage": self.coverage,
            "integrated_incident_units_m2": self.integrated_incident_units_m2,
        }

    def conservative_surface_upper_bounds(self):
        return {"coverage": 1.0, "integrated_incident_units_m2": None}

    def surface_field_remap_modes(self):
        return {"coverage": "intensive", "integrated_incident_units_m2": "conservative"}

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("manufactured state contract changed")
        return type(self)(fields["coverage"], fields["integrated_incident_units_m2"])


class _ScalarSurfaceMechanism:
    def __init__(
            self, probability, *, coverage_drift=0.0, flux_coupling=0.0,
            emit_product=False):
        self.probability = probability
        self.coverage_drift = float(coverage_drift)
        self.flux_coupling = float(flux_coupling)
        self.emit_product = bool(emit_product)

    def initial_state(self, shape=()):
        return _ScalarSurfaceState(np.zeros(shape), np.zeros(shape))

    def neutral_reaction_probability(self, state):
        value = np.asarray(self.probability(state.coverage), dtype=float)
        return {"A": np.broadcast_to(value, state.coverage.shape)}

    def advance(self, state, fluxes, duration_s):
        duration = float(duration_s)
        incident = np.asarray(fluxes.neutral_flux_m2_s["A"], dtype=float)
        coverage_rate = self.coverage_drift + self.flux_coupling * incident
        coverage = np.minimum(state.coverage + duration * coverage_rate, 1.0)
        integrated = state.integrated_incident_units_m2 + duration * incident
        removed = duration * incident
        deposited = 0.125 * duration * incident * (1.0 - state.coverage)
        if self.emit_product:
            exchange = SurfaceMaterialExchange(
                {"removed_units": removed}, {"removed_units": removed}, {},
                {"deposited_units": deposited})
            products = (SurfaceProductPopulation(
                "P", "removed_units", removed, 1.0, 1.0,
                provenance={"source": "manufactured emitted-product refusal"}),)
        else:
            exchange = SurfaceMaterialExchange(
                {"removed_units": removed}, {}, {"removed_units": removed},
                {"deposited_units": deposited})
            products = ()
        return SimpleNamespace(
            state=_ScalarSurfaceState(coverage, integrated),
            etch_velocity_m_s=1.0e-9 * incident,
            normal_growth_velocity_m_s=(
                2.0e-10 * incident * (1.0 - state.coverage)),
            material_exchange=exchange,
            product_populations=products,
            validity=VALID)


def _router(probability, **kwargs):
    return MaterialMechanismRouter3D(
        {1: _ScalarSurfaceMechanism(probability, **kwargs)},
        provenance={1: {"source": "manufactured SC1 gate"}})


def _cache(
        direct, factors, *, area=None, active=None,
        identity=None):
    direct = np.asarray(direct, dtype=float)
    face_count = direct.size
    return SurfaceRadiosityOperatorCache3D(
        SurfaceFluxes({"A": direct}), factors,
        np.ones(face_count) if area is None else np.asarray(area, dtype=float),
        np.arange(face_count) if active is None else np.asarray(active, dtype=int),
        np.ones(face_count, dtype=int), {"A": "neutral_reactant"},
        {"geometry_sha256": "g", "operator_epoch": 7}
        if identity is None else identity,
        relative_tolerance=1.0e-12, maximum_iterations=2000)


def _open_cache(direct=2.0):
    factors = DiffuseFormFactors3D(
        1, np.asarray([], dtype=int), np.asarray([], dtype=int),
        np.asarray([], dtype=float), np.asarray([1.0]), 1)
    return _cache([direct], factors)


def _cavity_cache():
    factors = DiffuseFormFactors3D(
        2, np.asarray([0, 1]), np.asarray([1, 0]),
        np.asarray([0.8, 0.8]), np.asarray([0.2, 0.2]), 10)
    return _cache([1.0, 0.2], factors)


def _identity():
    return {"geometry_sha256": "g", "operator_epoch": 7}


def _state(router, face_count):
    return router.initial_state_by_material(np.ones(face_count, dtype=int))


def _integrate(cache, router, state, duration, *, error=0.02, dp=0.05, minimum=1e-6):
    return integrate_surface_radiosity_chemistry_3d(
        cache, router, state, duration,
        operator_identity_payload=_identity(),
        maximum_embedded_relative_error=error,
        maximum_local_reaction_probability_change=dp,
        minimum_chemistry_substep_s=minimum)


def test_state_independent_probability_reproduces_production_radiosity_equation():
    cache = _cavity_cache()
    router = _router(lambda coverage: np.full_like(coverage, 0.3))
    state = _state(router, 2)

    observed = cache.solve(
        state, router, operator_identity_payload=_identity())
    expected = solve_diffuse_neutral_radiosity_3d(
        [1.0, 0.2], [1.0, 1.0], [0, 1], [1, 0], [0.8, 0.8],
        [0.2, 0.2], [0.3, 0.3], relative_tolerance=1e-12,
        maximum_iterations=2000)

    assert np.array_equal(
        observed.surface_fluxes.neutral_flux_m2_s["A"],
        expected.incident_flux_m2_s)
    assert observed.maximum_relative_balance_error <= 2e-12
    assert not cache.face_area_m2.flags.writeable
    assert not cache.form_factors.transfer_fraction.flags.writeable


def test_one_face_open_limit_is_exact_and_inputs_remain_immutable():
    cache = _open_cache()
    router = _router(lambda coverage: np.full_like(coverage, 0.4))
    state = _state(router, 1)
    before = {name: value.copy() for name, value in state.fields.items()}
    cache_sha = cache.cache_sha256

    result = _integrate(cache, router, state, 0.25, error=1e-12, dp=0.01)

    assert result.diagnostics.rejected_trial_count == 0
    assert result.diagnostics.accepted_trial_count == 1
    assert result.diagnostics.maximum_accepted_embedded_relative_error == 0.0
    assert np.array_equal(result.material_exchange.removed_units_m2["removed_units"], [0.5])
    assert np.array_equal(result.material_exchange.unresolved_units_m2["removed_units"], [0.5])
    assert np.array_equal(result.integrated_recession_m, [0.5e-9])
    assert all(np.array_equal(state.fields[name], value) for name, value in before.items())
    assert cache.cache_sha256 == cache_sha
    assert result.diagnostics.rejected_trial_exchange_contribution_is_zero
    trial = result.diagnostics.trials[0]
    assert set(trial.embedded_error.exchange_inventory) == {
        "removed", "outgoing", "unresolved", "deposited"}
    assert set(trial.embedded_error.state_increment) == set(state.fields)


def test_two_face_state_dependent_cavity_agrees_with_tight_reference():
    cache = _cavity_cache()
    probability = lambda coverage: 0.12 + 0.7 * coverage
    router = _router(probability, coverage_drift=0.03, flux_coupling=0.22)
    state = _state(router, 2)

    candidate = _integrate(
        cache, router, state, 0.5, error=0.02, dp=0.05, minimum=1e-5)
    reference = _integrate(
        cache, router, state, 0.5, error=2e-4, dp=0.005, minimum=1e-7)

    for name in state.fields:
        assert np.allclose(
            candidate.state.fields[name], reference.state.fields[name],
            rtol=7e-3, atol=1e-12)
    assert np.allclose(
        candidate.material_exchange.removed_units_m2["removed_units"],
        reference.material_exchange.removed_units_m2["removed_units"],
        rtol=7e-3, atol=1e-12)
    assert candidate.diagnostics.maximum_radiosity_relative_balance_error <= 2e-12


def test_nonmonotone_probability_path_rejects_equal_endpoints_without_ledger_leak():
    cache = _open_cache(1.0)
    router = _router(
        lambda coverage: 4.0 * coverage * (1.0 - coverage),
        coverage_drift=1.0)
    state = _state(router, 1)

    result = _integrate(
        cache, router, state, 1.0, error=2.0, dp=0.3, minimum=1e-4)
    first = result.diagnostics.trials[0]

    assert not first.accepted
    assert first.rejection_reason == "reaction_probability_safety_cap"
    assert first.start_to_final_probability_change.maximum == 0.0
    assert first.start_to_midpoint_probability_change.maximum == 1.0
    assert result.diagnostics.rejected_trial_count > 0

    # Re-run exactly the accepted 1/16 path without the rejected proposals.  If a
    # rejected exchange leaked into the accumulator, these ledgers would disagree.
    manual_state = state
    manual_removed = np.zeros(1)
    for _ in range(16):
        step = _integrate(
            cache, router, manual_state, 1.0 / 16.0,
            error=2.0, dp=1.1, minimum=1e-4)
        manual_state = step.state
        manual_removed += step.material_exchange.removed_units_m2["removed_units"]
    assert np.array_equal(result.state.fields["m1__coverage"], manual_state.fields["m1__coverage"])
    assert np.array_equal(
        result.material_exchange.removed_units_m2["removed_units"], manual_removed)


def test_timestep_tightening_converges_toward_same_cavity_solution():
    cache = _cavity_cache()
    router = _router(
        lambda coverage: 0.1 + 0.75 * coverage,
        coverage_drift=0.02, flux_coupling=0.25)
    state = _state(router, 2)
    loose = _integrate(cache, router, state, 0.5, error=0.08, dp=0.12, minimum=1e-6)
    medium = _integrate(cache, router, state, 0.5, error=0.02, dp=0.04, minimum=1e-6)
    tight = _integrate(cache, router, state, 0.5, error=0.005, dp=0.015, minimum=1e-7)

    target = tight.state.fields["m1__integrated_incident_units_m2"]
    loose_error = np.linalg.norm(
        loose.state.fields["m1__integrated_incident_units_m2"] - target)
    medium_error = np.linalg.norm(
        medium.state.fields["m1__integrated_incident_units_m2"] - target)
    assert medium_error < loose_error
    refinement = (
        loose.diagnostics.minimum_accepted_chemistry_substep_s
        / medium.diagnostics.minimum_accepted_chemistry_substep_s)
    contraction = loose_error / medium_error
    assert np.isclose(refinement, 4.0)
    # The embedded method is first order globally.  Against the still-finer
    # reference, this manufactured cavity contracts slightly faster than the
    # nominal 4x step refinement; keep a quantitative, non-overfitted bound.
    assert 3.0 <= contraction <= 8.0


def test_deterministic_replay_is_exact_including_rejection_diagnostics():
    cache = _open_cache(1.0)
    router = _router(
        lambda coverage: 0.15 + 0.7 * coverage,
        coverage_drift=0.4)
    state = _state(router, 1)

    first = _integrate(cache, router, state, 0.25, error=0.01, dp=0.02, minimum=1e-5)
    second = _integrate(cache, router, state, 0.25, error=0.01, dp=0.02, minimum=1e-5)

    assert first.diagnostics.rejected_trial_count > 0
    assert first.diagnostics == second.diagnostics
    assert set(first.state.fields) == set(second.state.fields)
    assert all(np.array_equal(first.state.fields[name], second.state.fields[name])
               for name in first.state.fields)
    for attribute in (
            "removed_units_m2", "outgoing_units_m2",
            "unresolved_units_m2", "deposited_units_m2"):
        left = getattr(first.material_exchange, attribute)
        right = getattr(second.material_exchange, attribute)
        assert set(left) == set(right)
        assert all(np.array_equal(left[name], right[name]) for name in left)
    assert np.array_equal(first.integrated_recession_m, second.integrated_recession_m)
    assert np.array_equal(first.integrated_growth_m, second.integrated_growth_m)


def test_terminal_minimum_step_refusal_leaves_state_and_cache_unchanged():
    cache = _open_cache(1.0)
    router = _router(
        lambda coverage: 4.0 * coverage * (1.0 - coverage),
        coverage_drift=1.0)
    state = _state(router, 1)
    state_before = {name: value.copy() for name, value in state.fields.items()}
    direct_before = cache.direct_surface_fluxes.neutral_flux_m2_s["A"].copy()
    cache_sha = cache.cache_sha256

    with pytest.raises(SurfaceRadiosityCouplingRefusal, match="minimum"):
        _integrate(cache, router, state, 1.0, error=2.0, dp=0.01, minimum=0.1)

    assert all(np.array_equal(state.fields[name], value)
               for name, value in state_before.items())
    assert np.array_equal(
        cache.direct_surface_fluxes.neutral_flux_m2_s["A"], direct_before)
    assert cache.cache_sha256 == cache_sha


def test_cache_identity_and_incomplete_active_routing_refuse():
    cache = _open_cache()
    router = _router(lambda coverage: np.full_like(coverage, 0.5))
    state = _state(router, 1)
    with pytest.raises(SurfaceRadiosityCacheIdentityError, match="identity"):
        cache.solve(
            state, router,
            operator_identity_payload={"geometry_sha256": "changed", "operator_epoch": 7})

    factors = DiffuseFormFactors3D(
        2, np.asarray([0, 1]), np.asarray([1, 0]),
        np.asarray([0.5, 0.5]), np.asarray([0.5, 0.5]), 2)
    with pytest.raises(SurfaceRadiosityCouplingRefusal, match="incomplete active-face"):
        _cache([1.0, 1.0], factors, active=[0])

    class MissingProbabilityRouter:
        def neutral_reaction_probability_by_material(self, _state, _material):
            return {}

    with pytest.raises(
            SurfaceRadiosityCouplingRefusal, match="cover every cached neutral"):
        cache.solve(
            state, MissingProbabilityRouter(),
            operator_identity_payload=_identity())

    class SupersetProbabilityRouter:
        def neutral_reaction_probability_by_material(self, _state, _material):
            return {"A": np.asarray([0.5]), "unused": np.asarray([1.0])}

    solved = cache.solve(
        state, SupersetProbabilityRouter(),
        operator_identity_payload=_identity())
    assert set(solved.reaction_probability) == {"A"}


def test_product_populations_refuse_and_zero_duration_is_bitwise_noop():
    cache = _open_cache()
    router = _router(lambda coverage: np.full_like(coverage, 0.5))
    state = _state(router, 1)

    zero = _integrate(cache, router, state, 0.0, minimum=1e-6)
    assert zero.state is state
    assert zero.diagnostics.radiosity_solve_count == 0
    assert zero.diagnostics.trials == ()
    assert np.array_equal(zero.integrated_recession_m, np.zeros(1))
    assert not zero.integrated_recession_m.flags.writeable

    product_router = _router(
        lambda coverage: np.full_like(coverage, 0.5), emit_product=True)
    product_state = _state(product_router, 1)
    with pytest.raises(SurfaceRadiosityCouplingRefusal, match="product populations"):
        _integrate(cache, product_router, product_state, 0.1)
