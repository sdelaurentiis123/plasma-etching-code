import json
from pathlib import Path

import numpy as np
import pytest

from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)
from petch.reactor_global.electron_kinetics import (
    ELECTRON_SPEED_PER_SQRT_EV_M_S,
    ElectronCollisionMomentKernel,
    ElectronEnergyDistribution,
    ElectronEnergyGrid,
    DeterministicTwoTermBoltzmannSolver,
    ScharfetterGummelEnergyOperator,
    TwoTermBoltzmannCondition,
    normalize_eepf,
    normalize_eepf_jvp,
    normalize_eepf_vjp,
)


ROOT = Path(__file__).resolve().parents[1]


def _constant_elastic(maximum_energy_eV=200.0, cross_section_m2=2.0e-20):
    return ElectronCollisionProcess(
        kind="ELASTIC",
        target="manufactured",
        product=None,
        electron_energy_eV=(0.0, maximum_energy_eV),
        cross_section_m2=(cross_section_m2, cross_section_m2),
        mass_ratio=1.0e-5,
    )


def _deck(*processes):
    return ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256="c" * 64,
        source_database="manufactured-test-deck",
        retrieved_at="2026-08-08",
        source_reference="tests only; not physical evidence",
    )


def test_fixed_grid_carries_exact_piecewise_constant_eepf_weights():
    grid = ElectronEnergyGrid((0.0, 1.0, 3.0, 7.0))
    assert grid.cell_count == 3
    assert np.sum(grid.normalization_weights) == pytest.approx(
        (2.0 / 3.0) * 7.0 ** 1.5)
    assert np.sum(grid.mean_energy_weights_eV) == pytest.approx(
        (2.0 / 5.0) * 7.0 ** 2.5)

    constant = ElectronEnergyDistribution.from_unnormalized(
        grid, np.ones(grid.cell_count))
    assert constant.normalization == pytest.approx(1.0)
    assert constant.mean_energy_eV == pytest.approx(3.0 * 7.0 / 5.0)

    with pytest.raises(ValueError, match="invalid electron-energy grid"):
        ElectronEnergyGrid((0.0, 1.0, 1.0))


def test_maxwellian_eepf_recovers_mean_energy_and_batch_first_shape():
    grid = ElectronEnergyGrid.linear(120.0, 4800)
    temperatures = np.array([1.0, 2.5, 5.0])
    distribution = ElectronEnergyDistribution.maxwellian(
        grid, temperatures)
    assert distribution.batch_shape == (3,)
    assert distribution.normalization == pytest.approx(np.ones(3))
    assert distribution.mean_energy_eV == pytest.approx(
        1.5 * temperatures, rel=2.0e-4)
    receipt = distribution.convergence_receipt(
        tail_cell_count=4, maximum_tail_population_fraction=1.0e-8)
    assert receipt["normalization_passed"]
    assert receipt["positivity_passed"]
    assert receipt["tail_passed"]
    assert not receipt["supports_reactor_state_prediction"]
    assert not receipt["supports_feature_depth"]


def test_normalization_jvp_and_vjp_are_exact_duals():
    rng = np.random.default_rng(1945)
    grid = ElectronEnergyGrid.linear(20.0, 80)
    raw = np.exp(-grid.cell_centers_eV[None, :] / np.array([[1.2], [3.4]]))
    # Keep both finite-difference samples inside the nonnegative EEPF domain.
    tangent = raw * rng.normal(size=raw.shape)
    cotangent = rng.normal(size=raw.shape)

    jvp = normalize_eepf_jvp(grid, raw, tangent)
    vjp = normalize_eepf_vjp(grid, raw, cotangent)
    assert np.sum(cotangent * jvp) == pytest.approx(
        np.sum(vjp * tangent), rel=2.0e-13, abs=2.0e-13)
    assert np.sum(jvp * grid.normalization_weights, axis=-1) == pytest.approx(
        np.zeros(2), abs=2.0e-13)

    step = 1.0e-7
    finite_difference = (
        normalize_eepf(grid, raw + step * tangent)
        - normalize_eepf(grid, raw - step * tangent)
    ) / (2.0 * step)
    assert jvp == pytest.approx(finite_difference, rel=2.0e-7, abs=2.0e-9)


def test_collision_moments_match_constant_cross_section_maxwellian_limit():
    temperature_eV = 2.5
    cross_section_m2 = 2.0e-20
    grid = ElectronEnergyGrid.linear(100.0, 4000)
    distribution = ElectronEnergyDistribution.maxwellian(
        grid, temperature_eV)
    kernel = ElectronCollisionMomentKernel.from_process(
        grid, _constant_elastic(200.0, cross_section_m2))
    moments = kernel.evaluate(distribution)

    expected_rate = (
        cross_section_m2
        * ELECTRON_SPEED_PER_SQRT_EV_M_S
        * 2.0 / np.sqrt(np.pi)
        * np.sqrt(temperature_eV)
    )
    assert moments.rate_coefficient_m3_s == pytest.approx(
        expected_rate, rel=2.0e-4)
    assert moments.collision_weighted_mean_incident_energy_eV == pytest.approx(
        2.0 * temperature_eV, rel=3.0e-4)
    assert moments.unresolved_population_fraction == pytest.approx(0.0)
    assert not moments.supports_swarm_validation
    assert not moments.supports_feature_depth


def test_collision_kernel_fails_closed_on_unmeasured_eepf_support():
    grid = ElectronEnergyGrid.linear(20.0, 400)
    distribution = ElectronEnergyDistribution.maxwellian(grid, 5.0)
    kernel = ElectronCollisionMomentKernel.from_process(
        grid, _constant_elastic(10.0))
    with pytest.raises(ValueError, match="outside collision cross-section"):
        kernel.evaluate(
            distribution,
            maximum_unresolved_population_fraction=1.0e-4,
        )


def test_collision_rate_and_energy_moment_jvp_vjp_dot_products_close():
    rng = np.random.default_rng(2005)
    grid = ElectronEnergyGrid.linear(40.0, 240)
    kernel = ElectronCollisionMomentKernel.from_process(
        grid, _constant_elastic(80.0))
    tangent = rng.normal(size=(3, grid.cell_count))
    rate_cotangent = rng.normal(size=3)
    energy_cotangent = rng.normal(size=3)

    rate_left = np.sum(
        rate_cotangent * np.asarray(kernel.rate_jvp(tangent)))
    rate_right = np.sum(kernel.rate_vjp(rate_cotangent) * tangent)
    assert rate_left == pytest.approx(rate_right, rel=2.0e-15)

    energy_left = np.sum(
        energy_cotangent
        * np.asarray(kernel.incident_energy_jvp(tangent)))
    energy_right = np.sum(
        kernel.incident_energy_vjp(energy_cotangent) * tangent)
    assert energy_left == pytest.approx(energy_right, rel=2.0e-15)


def test_scharfetter_gummel_is_exact_for_local_exponential_equilibrium():
    grid = ElectronEnergyGrid.linear(20.0, 160)
    drift = np.full(grid.cell_count - 1, -0.7)
    diffusion = np.full(grid.cell_count - 1, 1.3)
    operator = ScharfetterGummelEnergyOperator(grid, drift, diffusion)
    solution = operator.solve(np.zeros(grid.cell_count))

    expected = normalize_eepf(
        grid,
        np.exp((drift[0] / diffusion[0]) * grid.cell_centers_eV),
    )
    assert solution.distribution.eepf_eV_minus_3_over_2 == pytest.approx(
        expected, rel=5.0e-10, abs=3.0e-13)
    assert solution.energy_flux_faces == pytest.approx(
        np.zeros(grid.cell_count + 1), abs=3.0e-14)
    assert solution.compatibility_multiplier == pytest.approx(0.0, abs=1.0e-13)
    assert solution.maximum_augmented_residual < 1.0e-12
    assert solution.maximum_physical_residual < 1.0e-12
    assert not solution.supports_collision_boltzmann_solve
    assert not solution.supports_feature_depth


def test_scharfetter_gummel_batches_and_implicit_adjoint_are_exact_duals():
    rng = np.random.default_rng(44)
    grid = ElectronEnergyGrid.linear(12.0, 60)
    drift = np.stack((
        np.full(grid.cell_count - 1, -0.4),
        np.full(grid.cell_count - 1, -0.9),
    ))
    diffusion = np.stack((
        np.full(grid.cell_count - 1, 1.1),
        np.full(grid.cell_count - 1, 1.7),
    ))
    operator = ScharfetterGummelEnergyOperator(grid, drift, diffusion)
    solution = operator.solve(np.zeros((2, grid.cell_count)))
    assert solution.distribution.batch_shape == (2,)
    assert solution.distribution.normalization == pytest.approx(np.ones(2))

    source_tangent = rng.normal(size=(2, grid.cell_count))
    distribution_cotangent = rng.normal(size=(2, grid.cell_count))
    multiplier_cotangent = rng.normal(size=2)
    state_tangent, compatibility_tangent = operator.implicit_source_jvp(
        source_tangent)
    source_cotangent = operator.implicit_source_vjp(
        distribution_cotangent, multiplier_cotangent)
    left = (
        np.sum(distribution_cotangent * state_tangent)
        + np.sum(multiplier_cotangent * compatibility_tangent)
    )
    right = np.sum(source_cotangent * source_tangent)
    assert left == pytest.approx(right, rel=5.0e-12, abs=5.0e-12)
    assert np.sum(
        state_tangent * grid.normalization_weights, axis=-1
    ) == pytest.approx(np.zeros(2), abs=3.0e-12)


def test_scharfetter_gummel_refuses_nonconservative_closed_domain_source():
    grid = ElectronEnergyGrid.linear(10.0, 40)
    operator = ScharfetterGummelEnergyOperator(
        grid,
        np.full(grid.cell_count - 1, -0.5),
        np.full(grid.cell_count - 1, 1.0),
    )
    with pytest.raises(ValueError, match="incompatible with zero-flux"):
        operator.solve(np.ones(grid.cell_count))


def test_two_term_zero_field_elastic_limit_is_gas_maxwellian():
    grid = ElectronEnergyGrid.linear(1.0, 240)
    solver = DeterministicTwoTermBoltzmannSolver(
        grid, _deck(_constant_elastic(10.0)))
    condition = TwoTermBoltzmannCondition(
        reduced_electric_field_Td=0.0,
        gas_temperature_K=300.0,
        target_mole_fractions={"manufactured": 1.0},
        initial_electron_temperature_eV=0.5,
    )
    solution = solver.solve(condition, damping=1.0, relative_tolerance=1.0e-11)
    gas_thermal_energy_eV = 1.380649e-23 * 300.0 / 1.602176634e-19
    expected = ElectronEnergyDistribution.maxwellian(
        grid, gas_thermal_energy_eV)
    assert solution.distribution.eepf_eV_minus_3_over_2 == pytest.approx(
        expected.eepf_eV_minus_3_over_2, rel=8.0e-10, abs=1.0e-13)
    assert solution.distribution.mean_energy_eV == pytest.approx(
        1.5 * gas_thermal_energy_eV, rel=2.0e-3)
    assert solution.iteration_count == 1
    assert solution.maximum_equation_residual_m3_s < 2.0e-25
    assert (
        solution.transport_moments.reduced_field_power_gain_eV_m3_s
        == 0.0
    )
    expected_mean_speed = (
        ELECTRON_SPEED_PER_SQRT_EV_M_S
        * 2.0 / np.sqrt(np.pi)
        * np.sqrt(gas_thermal_energy_eV)
    )
    assert solution.transport_moments.mean_electron_speed_m_s == (
        pytest.approx(expected_mean_speed, rel=2.0e-3))
    assert solution.transport_moments.isotropic_wall_flux_coefficient_m_s == (
        pytest.approx(0.25 * expected_mean_speed, rel=2.0e-3))
    assert solution.transport_moments.mean_wall_loss_electron_energy_eV == (
        pytest.approx(2.0 * gas_thermal_energy_eV, rel=3.0e-3))
    assert solution.transport_moments.supports_flux_transport_moments
    assert not solution.transport_moments.supports_direct_swarm_grade
    assert not solution.transport_moments.supports_reactor_state_prediction
    assert solution.particle_growth_closure_error_m3_s == pytest.approx(0.0)
    assert solution.supports_collision_boltzmann_solve
    assert not solution.supports_direct_swarm_grade
    assert not solution.supports_reactor_state_prediction
    assert not solution.supports_feature_depth


def test_two_term_excitation_source_conserves_particles_and_cools_field_tail():
    elastic = _constant_elastic(300.0, 2.0e-20)
    excitation = ElectronCollisionProcess(
        kind="EXCITATION",
        target="manufactured",
        product="manufactured(v=1)",
        electron_energy_eV=(0.0, 2.0, 5.0, 300.0),
        cross_section_m2=(0.0, 0.0, 1.2e-20, 1.2e-20),
        energy_loss_eV=2.0,
    )
    grid = ElectronEnergyGrid.linear(200.0, 600)
    condition = TwoTermBoltzmannCondition(
        reduced_electric_field_Td=5.0,
        gas_temperature_K=300.0,
        target_mole_fractions={"manufactured": 1.0},
    )
    elastic_only = DeterministicTwoTermBoltzmannSolver(
        grid, _deck(elastic)).solve(
            condition,
            damping=1.0,
            maximum_tail_population_fraction=1.0e-6,
        )
    with_excitation = DeterministicTwoTermBoltzmannSolver(
        grid, _deck(elastic, excitation)).solve(
            condition,
            damping=1.0,
            maximum_tail_population_fraction=1.0e-6,
        )
    assert with_excitation.particle_source_from_collision_operator_m3_s == (
        pytest.approx(0.0, abs=2.0e-27))
    assert with_excitation.particle_growth_closure_error_m3_s == pytest.approx(
        0.0, abs=2.0e-27)
    assert with_excitation.roundoff_negative_population_fraction < 1.0e-15
    assert with_excitation.distribution.mean_energy_eV < (
        elastic_only.distribution.mean_energy_eV)
    excitation_moment = with_excitation.collision_moments[1]
    assert excitation_moment.rate_coefficient_m3_s > 0.0
    assert excitation_moment.collision_weighted_mean_incident_energy_eV > 2.0


def test_two_term_requires_growth_model_for_attachment_or_ionization():
    attachment = ElectronCollisionProcess(
        kind="ATTACHMENT",
        target="manufactured",
        product="negative + fragment",
        electron_energy_eV=(0.0, 1.0, 20.0),
        cross_section_m2=(1.0e-21, 1.0e-21, 0.0),
        energy_loss_eV=0.0,
    )
    solver = DeterministicTwoTermBoltzmannSolver(
        ElectronEnergyGrid.linear(10.0, 100),
        _deck(_constant_elastic(20.0), attachment),
    )
    with pytest.raises(ValueError, match="require temporal_growth"):
        solver.solve(TwoTermBoltzmannCondition(
            reduced_electric_field_Td=20.0,
            gas_temperature_K=300.0,
            target_mole_fractions={"manufactured": 1.0},
        ))


def test_temporal_growth_ionization_closes_particle_source_and_rate_moment():
    elastic = _constant_elastic(200.0, 2.0e-20)
    excitation = ElectronCollisionProcess(
        kind="EXCITATION",
        target="manufactured",
        product="manufactured(v=1)",
        electron_energy_eV=(0.0, 2.0, 5.0, 200.0),
        cross_section_m2=(0.0, 0.0, 1.0e-20, 1.0e-20),
        energy_loss_eV=2.0,
    )
    ionization = ElectronCollisionProcess(
        kind="IONIZATION",
        target="manufactured",
        product="manufactured+",
        electron_energy_eV=(0.0, 10.0, 15.0, 200.0),
        cross_section_m2=(0.0, 0.0, 4.0e-21, 4.0e-21),
        energy_loss_eV=10.0,
    )
    solution = DeterministicTwoTermBoltzmannSolver(
        ElectronEnergyGrid.linear(100.0, 600),
        _deck(elastic, excitation, ionization),
    ).solve(
        TwoTermBoltzmannCondition(
            reduced_electric_field_Td=50.0,
            gas_temperature_K=300.0,
            target_mole_fractions={"manufactured": 1.0},
            growth_model="temporal_growth",
        ),
        damping=0.5,
        relative_tolerance=1.0e-8,
        maximum_tail_population_fraction=1.0e-5,
    )
    ionization_rate = solution.collision_moments[2].rate_coefficient_m3_s
    assert solution.net_growth_rate_coefficient_m3_s == pytest.approx(
        ionization_rate, rel=3.0e-13)
    assert solution.particle_source_from_collision_operator_m3_s == (
        pytest.approx(ionization_rate, rel=3.0e-13))
    assert solution.particle_growth_closure_error_m3_s == pytest.approx(
        0.0, abs=1.0e-28)
    assert solution.maximum_equation_residual_m3_s < 1.0e-23
    assert solution.iteration_count < 100


def test_multiple_ionization_growth_uses_declared_electron_multiplicity():
    double_ionization = ElectronCollisionProcess(
        kind="IONIZATION",
        target="manufactured",
        product="manufactured++",
        electron_energy_eV=(0.0, 15.0, 20.0, 200.0),
        cross_section_m2=(0.0, 0.0, 4.0e-21, 4.0e-21),
        energy_loss_eV=15.0,
        electron_number_change=2,
    )
    solution = DeterministicTwoTermBoltzmannSolver(
        ElectronEnergyGrid.linear(100.0, 600),
        _deck(_constant_elastic(200.0), double_ionization),
    ).solve(
        TwoTermBoltzmannCondition(
            reduced_electric_field_Td=80.0,
            gas_temperature_K=300.0,
            target_mole_fractions={"manufactured": 1.0},
            growth_model="temporal_growth",
        ),
        relative_tolerance=1.0e-8,
        maximum_tail_population_fraction=1.0e-5,
    )
    ionization_rate = solution.collision_moments[1].rate_coefficient_m3_s
    assert solution.net_growth_rate_coefficient_m3_s == pytest.approx(
        2.0 * ionization_rate, rel=3.0e-13)
    assert solution.particle_growth_closure_error_m3_s == 0.0


def test_temporal_growth_eigenroot_crosses_attachment_ionization_balance():
    attachment = ElectronCollisionProcess(
        kind="ATTACHMENT",
        target="manufactured",
        product="negative + fragment",
        electron_energy_eV=(0.0, 0.2, 1.0, 100.0),
        cross_section_m2=(2.0e-20, 2.0e-20, 0.0, 0.0),
        energy_loss_eV=0.0,
    )
    excitation = ElectronCollisionProcess(
        kind="EXCITATION",
        target="manufactured",
        product="manufactured(v=1)",
        electron_energy_eV=(0.0, 2.0, 5.0, 100.0),
        cross_section_m2=(0.0, 0.0, 1.0e-20, 1.0e-20),
        energy_loss_eV=2.0,
    )
    ionization = ElectronCollisionProcess(
        kind="IONIZATION",
        target="manufactured",
        product="positive + electron",
        electron_energy_eV=(0.0, 10.0, 15.0, 100.0),
        cross_section_m2=(0.0, 0.0, 2.0e-20, 2.0e-20),
        energy_loss_eV=10.0,
    )
    solver = DeterministicTwoTermBoltzmannSolver(
        ElectronEnergyGrid.linear(80.0, 480),
        _deck(
            _constant_elastic(100.0, 2.0e-20),
            excitation,
            attachment,
            ionization,
        ),
    )
    solutions = tuple(
        solver.solve(
            TwoTermBoltzmannCondition(
                reduced_electric_field_Td=field_Td,
                gas_temperature_K=300.0,
                target_mole_fractions={"manufactured": 1.0},
                growth_model="temporal_growth",
            ),
            relative_tolerance=1.0e-8,
            maximum_tail_population_fraction=1.0e-5,
        )
        for field_Td in (30.0, 100.0)
    )
    assert solutions[0].net_growth_rate_coefficient_m3_s < 0.0
    assert solutions[1].net_growth_rate_coefficient_m3_s > 0.0
    for solution in solutions:
        assert solution.particle_growth_closure_error_m3_s == 0.0
        assert solution.maximum_equation_residual_m3_s < 1.0e-23
        assert solution.weighted_iteration_residual < 1.0e-8
        assert solution.iteration_count < 100


def test_independent_bolos_oracle_receipt_converges_without_claim_inflation():
    receipt = json.loads((
        ROOT / "results/curated/reactor_global_kinetics/"
        "two_term_bolos_oracle_v1.json"
    ).read_text())
    assert receipt["schema"] == "petch.two_term_bolos_oracle.v1"
    assert [row["cell_count"] for row in receipt["rows"]] == [600, 1200, 2400]
    assert receipt["finest_grid_pass"]
    finest = receipt["rows"][-1]
    assert finest["mean_energy_relative_error"] < 0.01
    assert finest["excitation_rate_relative_error"] < 0.01
    assert finest["flux_reduced_mobility_relative_error"] < 0.01
    assert finest["scalar_reduced_diffusion_relative_error"] < 0.01
    assert finest["eepf_weighted_l1"] < 0.002
    assert not receipt["supports_direct_swarm_grade"]
    assert not receipt["supports_reactor_state_prediction"]
    assert not receipt["supports_wafer_flux"]
    assert not receipt["supports_feature_depth"]
