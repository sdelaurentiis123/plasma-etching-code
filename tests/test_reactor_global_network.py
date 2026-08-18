import numpy as np
import pytest

from petch.reactor_global import (
    CM3_TO_M3,
    E_CHARGE_C,
    ELECTRON_MASS_KG,
    INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT,
    KEMANECI_ELASTIC_ENERGY_APPROXIMATION,
    STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT,
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronInverseTemperaturePolynomialRateCoefficient,
    ElectronLogTemperatureInversePolynomialRateCoefficient,
    ElectronBase10LogPolynomialRateCoefficient,
    ElectronMaxwellianCrossSectionRateCoefficient,
    ElectronTabulatedCrossSectionSupport,
    ElectronTemperatureTabulatedRateCoefficient,
    ElectronAnalyticRateTerm,
    ElectronCompositeRateCoefficient,
    ElectronDetailedBalanceRateCoefficient,
    GasTemperatureArrheniusRateCoefficient,
    RateContext,
    Reaction,
    ReactionNetwork,
    Species,
)


def _toy_species():
    return (
        Species(
            name="A", mass_amu=10.0, charge_number=0,
            composition={"X": 1}, role="neutral",
            source="manufactured"),
        Species(
            name="B", mass_amu=10.0, charge_number=0,
            composition={"X": 1}, role="neutral",
            source="manufactured"),
        Species(
            name="C", mass_amu=20.0, charge_number=0,
            composition={"X": 2}, role="neutral",
            source="manufactured"),
    )


def test_cm3_per_second_conversion_is_exact_at_declared_scale():
    coefficient = ConstantRateCoefficient.from_cm3_per_s(
        2.5e-7, source="manufactured")
    assert coefficient.value_si == 2.5e-7 * CM3_TO_M3
    assert coefficient.coefficient_si(RateContext(3.0)) == 2.5e-13


def test_first_order_decay_source_matches_analytic_derivative():
    species = _toy_species()[:2]
    reaction = Reaction(
        name="A_to_B",
        reactants={"A": 1},
        products={"B": 1},
        kinetic_orders={"A": 1},
        rate_coefficient=ConstantRateCoefficient.from_per_second(
            2.0, source="manufactured"),
        electron_energy_loss_eV=0.0,
        source="manufactured",
    )
    network = ReactionNetwork(species=species, reactions=(reaction,))
    source = network.source_vector_m3_s(
        {"A": 3.0, "B": 7.0}, RateContext(1.0))
    np.testing.assert_allclose(source, [-6.0, 6.0], rtol=0.0, atol=0.0)


def test_bimolecular_reaction_preserves_analytic_invariants():
    reaction = Reaction(
        name="A_plus_B_to_C",
        reactants={"A": 1, "B": 1},
        products={"C": 1},
        kinetic_orders={"A": 1, "B": 1},
        rate_coefficient=ConstantRateCoefficient(
            value_si=0.25, density_order=2, source="manufactured",
            source_units="m^3 s^-1", evidence_kind="derived"),
        electron_energy_loss_eV=0.0,
        source="manufactured",
    )
    network = ReactionNetwork(species=_toy_species(), reactions=(reaction,))
    source = network.source_vector_m3_s(
        {"A": 2.0, "B": 3.0, "C": 0.0}, RateContext(1.0))
    np.testing.assert_allclose(source, [-1.5, -1.5, 1.5], rtol=0.0, atol=0.0)
    assert source[0] - source[1] == 0.0
    assert source[0] + source[1] + 2.0 * source[2] == 0.0
    assert network.source_conservation_report(
        {"A": 2.0, "B": 3.0, "C": 0.0},
        RateContext(1.0),
    )["normalized_maximum_residual"] == 0.0
    assert network.stoichiometric_matrix is network.stoichiometric_matrix


def test_electron_impact_ionization_conserves_charge_with_explicit_electrons():
    electron = Species(
        name="e", mass_amu=5.48579909065e-4, charge_number=-1,
        composition={}, role="electron", source="CODATA",
        evidence_kind="measured")
    argon = Species(
        name="Ar", mass_amu=39.948, charge_number=0,
        composition={"Ar": 1}, role="neutral", source="NIST",
        evidence_kind="measured")
    ion = Species(
        name="Ar+", mass_amu=39.948, charge_number=1,
        composition={"Ar": 1}, role="positive_ion", source="NIST",
        evidence_kind="measured")
    reaction = Reaction(
        name="argon_ionization",
        reactants={"e": 1, "Ar": 1},
        products={"Ar+": 1, "e": 2},
        kinetic_orders={"e": 1, "Ar": 1},
        rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
            1.23e-7, activation_eV=18.68,
            source="lee-lieberman-1994-global Table 3"),
        electron_energy_loss_eV=15.76,
        source="lee-lieberman-1994-global Table 3",
    )
    network = ReactionNetwork(
        species=(electron, argon, ion), reactions=(reaction,))
    residual = network.reaction_conservation_residuals()["argon_ionization"]
    assert residual == {"elements": {"Ar": 0.0}, "charge_number": 0.0}


def test_nonconserving_plasma_shorthand_is_rejected():
    electron = Species(
        name="e", mass_amu=5.48579909065e-4, charge_number=-1,
        composition={}, role="electron", source="CODATA")
    metastable = Species(
        name="Ar*", mass_amu=39.948, charge_number=0,
        composition={"Ar": 1}, role="excited_neutral", source="NIST")
    argon = Species(
        name="Ar", mass_amu=39.948, charge_number=0,
        composition={"Ar": 1}, role="neutral", source="NIST")
    ion = Species(
        name="Ar+", mass_amu=39.948, charge_number=1,
        composition={"Ar": 1}, role="positive_ion", source="NIST")
    printed_shorthand = Reaction(
        name="pooling_without_emitted_electron",
        reactants={"Ar*": 2},
        products={"Ar": 1, "Ar+": 1},
        kinetic_orders={"Ar*": 2},
        rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
            6.2e-10, source="lee-lieberman-1994-global Table 3"),
        electron_energy_loss_eV=0.0,
        source="lee-lieberman-1994-global Table 3 printed shorthand",
    )
    with pytest.raises(ValueError, match="does not conserve atoms and charge"):
        ReactionNetwork(
            species=(electron, metastable, argon, ion),
            reactions=(printed_shorthand,),
        )


def test_rate_evaluation_rejects_negative_or_nonfinite_density():
    species = _toy_species()[:2]
    reaction = Reaction(
        name="A_to_B",
        reactants={"A": 1}, products={"B": 1},
        kinetic_orders={"A": 1},
        rate_coefficient=ConstantRateCoefficient.from_per_second(
            1.0, source="manufactured"),
        electron_energy_loss_eV=0.0, source="manufactured")
    network = ReactionNetwork(species=species, reactions=(reaction,))
    for invalid in (-1.0, np.nan, np.inf):
        with pytest.raises(ValueError):
            network.event_rates_m3_s(
                {"A": invalid, "B": 0.0}, RateContext(1.0))


def test_inverse_temperature_polynomial_rate_replays_printed_expression():
    coefficient = (
        ElectronInverseTemperaturePolynomialRateCoefficient.from_cm3_per_s(
            3.69e-10,
            inverse_temperature_coefficients=(
                -1.68, 1.457, -0.44, 0.0572, -0.0026),
            source="lee-lieberman-1994-global Table 2 k3",
        )
    )
    temperature = 3.25
    expected = 3.69e-10 * np.exp(
        -1.68 / temperature
        + 1.457 / temperature ** 2
        - 0.44 / temperature ** 3
        + 0.0572 / temperature ** 4
        - 0.0026 / temperature ** 5
    ) * CM3_TO_M3
    assert np.isclose(
        coefficient.coefficient_si(RateContext(temperature)),
        expected,
        rtol=2.0e-16,
        atol=0.0,
    )


def test_log_temperature_inverse_polynomial_replays_kokkoris_equation2():
    printed = (-34.92, 1.487, -2.377, -29.71, -0.1449)
    coefficient = (
        ElectronLogTemperatureInversePolynomialRateCoefficient.from_log_si(
            printed[0],
            temperature_power=printed[1],
            inverse_temperature_coefficients=printed[2:],
            source="Kokkoris et al. 2009 Table 1 G11 Druyvesteyn",
        )
    )
    temperature = 3.25
    expected = np.exp(
        printed[0]
        + printed[1] * np.log(temperature)
        + printed[2] / temperature
        + printed[3] / temperature ** 2
        + printed[4] / temperature ** 3
    )
    assert np.isclose(
        coefficient.coefficient_si(RateContext(temperature)),
        expected,
        rtol=1.0e-14,
        atol=0.0,
    )


def test_base10_log_polynomial_rate_replays_lennon_equation6():
    printed = (
        1.419e-7,
        -1.864e-8,
        -5.439e-8,
        3.306e-8,
        -3.54e-9,
        -2.915e-8,
    )
    coefficient = ElectronBase10LogPolynomialRateCoefficient.from_cm3_per_s(
        printed,
        reference_temperature_eV=12.96,
        activation_eV=12.96,
        temperature_power=0.5,
        source="lee-lieberman-1994-global Table 2 k4",
    )
    temperature = 4.0
    ratio = temperature / 12.96
    expected = (
        ratio ** 0.5
        * np.exp(-12.96 / temperature)
        * sum(
            value * np.log10(ratio) ** order
            for order, value in enumerate(printed)
        )
        * CM3_TO_M3
    )
    assert np.isclose(
        coefficient.coefficient_si(RateContext(temperature)),
        expected,
        rtol=2.0e-16,
        atol=0.0,
    )


@pytest.mark.parametrize("temperature", [1.29, 130.0])
def test_lennon_equation6_rejects_outside_published_temperature_domain(
        temperature):
    coefficient = ElectronBase10LogPolynomialRateCoefficient.from_cm3_per_s(
        (1.419e-7,),
        reference_temperature_eV=12.96,
        activation_eV=12.96,
        temperature_power=0.5,
        source="lennon-1988-ionization Eq. 6",
    )
    with pytest.raises(ValueError, match="outside the Lennon Eq.-6"):
        coefficient.coefficient_si(RateContext(temperature))


def test_tabulated_cross_section_maxwellian_integral_is_analytic():
    cross_section = 2.5e-20
    maximum_energy = 100.0
    temperature = 2.0
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, maximum_energy),
        cross_section_m2=(cross_section, cross_section),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured constant cross section",
        evidence_kind="derived",
    )
    support_ratio = maximum_energy / temperature
    retained_kernel = 1.0 - (
        support_ratio + 1.0) * np.exp(-support_ratio)
    expected = (
        cross_section
        * np.sqrt(8.0 * E_CHARGE_C * temperature
                  / (np.pi * ELECTRON_MASS_KG))
        * retained_kernel
    )
    assert coefficient.coefficient_si(
        RateContext(temperature)) == pytest.approx(
            expected, rel=3.0e-15, abs=0.0)


def test_tabulated_cross_section_incident_energy_moment_is_analytic():
    cross_section = 2.5e-20
    maximum_energy = 100.0
    temperature = 2.0
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, maximum_energy),
        cross_section_m2=(cross_section, cross_section),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured constant cross section",
        evidence_kind="derived",
    )
    support_ratio = maximum_energy / temperature
    retained_kernel = 1.0 - 0.5 * (
        support_ratio ** 2 + 2.0 * support_ratio + 2.0
    ) * np.exp(-support_ratio)
    expected = (
        2.0
        * temperature
        * cross_section
        * np.sqrt(8.0 * E_CHARGE_C * temperature
                  / (np.pi * ELECTRON_MASS_KG))
        * retained_kernel
    )
    assert coefficient.incident_energy_moment_eV_m3_s(
        RateContext(temperature)) == pytest.approx(
            expected, rel=3.0e-15, abs=0.0)


def test_energy_moment_applies_stricter_tail_gate_than_particle_rate():
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, 34.0),
        cross_section_m2=(1.0e-20, 1.0e-20),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured energy-tail gate",
        evidence_kind="derived",
    )
    context = RateContext(2.0)
    assert coefficient.coefficient_si(context) > 0.0
    with pytest.raises(ValueError, match="incident-energy Maxwellian"):
        coefficient.incident_energy_moment_eV_m3_s(context)


def test_tabulated_cross_section_uses_physical_threshold():
    common = dict(
        electron_energy_eV=(0.0, 1.0, 2.0, 100.0),
        threshold_eV=1.5,
        relative_uncertainty=0.2,
        source="manufactured threshold test",
        evidence_kind="measured",
    )
    first = ElectronMaxwellianCrossSectionRateCoefficient(
        cross_section_m2=(0.0, 9.0e-20, 1.0e-20, 1.0e-20),
        **common,
    )
    second = ElectronMaxwellianCrossSectionRateCoefficient(
        cross_section_m2=(8.0e-20, 7.0e-20, 1.0e-20, 1.0e-20),
        **common,
    )
    assert first.coefficient_si(RateContext(2.0)) == pytest.approx(
        second.coefficient_si(RateContext(2.0)),
        rel=0.0,
        abs=0.0,
    )


def test_tabulated_cross_section_allows_threshold_just_below_first_sample():
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(11.5, 12.0, 100.0),
        cross_section_m2=(0.03e-20, 0.11e-20, 6.19e-20),
        threshold_eV=11.481,
        relative_uncertainty=None,
        source="manufactured sampled-above-threshold test",
        evidence_kind="derived",
    )
    assert coefficient.coefficient_si(RateContext(2.0)) > 0.0


def test_tabulated_cross_section_fails_on_unmeasured_maxwellian_tail():
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, 100.0),
        cross_section_m2=(1.0e-20, 1.0e-20),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured short support",
        evidence_kind="derived",
    )
    with pytest.raises(ValueError, match="unmeasured cross-section tail"):
        coefficient.coefficient_si(RateContext(50.0))


def test_tabulated_cross_section_support_integrates_rate_and_energy_moments():
    cross_section = 2.5e-20
    maximum_energy = 100.0
    temperature = 2.0
    support = ElectronTabulatedCrossSectionSupport(
        electron_energy_eV=(0.0, maximum_energy),
        cross_section_m2=(cross_section, cross_section),
        relative_uncertainty=None,
        source="manufactured constant support",
        evidence_kind="derived",
    )
    ratio = maximum_energy / temperature
    rate_retained = 1.0 - (ratio + 1.0) * np.exp(-ratio)
    energy_retained = 1.0 - 0.5 * (
        ratio ** 2 + 2.0 * ratio + 2.0) * np.exp(-ratio)
    speed_scale = np.sqrt(
        8.0 * E_CHARGE_C * temperature / (np.pi * ELECTRON_MASS_KG)
    )
    expected_rate = cross_section * speed_scale * rate_retained
    expected_energy = (
        2.0 * temperature * cross_section * speed_scale * energy_retained
    )
    context = RateContext(temperature)
    assert support.tabulated_rate_coefficient_si(context) == pytest.approx(
        expected_rate, rel=3.0e-15, abs=0.0)
    assert (
        support.tabulated_incident_energy_moment_eV_m3_s(context)
        == pytest.approx(expected_energy, rel=3.0e-15, abs=0.0)
    )
    assert support.rate_kernel_missing_fractions(temperature) == (
        pytest.approx(0.0, abs=0.0),
        pytest.approx(
            (ratio + 1.0) * np.exp(-ratio), rel=2.0e-15, abs=0.0),
    )
    assert support.incident_energy_kernel_missing_fractions(temperature) == (
        pytest.approx(0.0, abs=0.0),
        pytest.approx(
            0.5 * (ratio ** 2 + 2.0 * ratio + 2.0) * np.exp(-ratio),
            rel=2.0e-15,
            abs=0.0,
        ),
    )


def test_tabulated_cross_section_support_does_not_hide_missing_tails():
    support = ElectronTabulatedCrossSectionSupport(
        electron_energy_eV=(0.05, 11.8),
        cross_section_m2=(1.83e-20, 0.0043e-20),
        relative_uncertainty=None,
        source="manufactured finite support",
        evidence_kind="derived",
    )
    rate_low, rate_high = support.rate_kernel_missing_fractions(3.0)
    energy_low, energy_high = (
        support.incident_energy_kernel_missing_fractions(3.0)
    )
    assert 0.0 < rate_low < 0.001
    assert 0.09 < rate_high < 0.10
    assert 0.0 < energy_low < 1.0e-5
    assert 0.24 < energy_high < 0.25


def test_tabulated_temperature_rate_is_exact_at_nodes_and_positive_between():
    coefficient = ElectronTemperatureTabulatedRateCoefficient(
        electron_temperature_eV=(1.0, 2.0, 4.0),
        coefficient_m3_s=(1.0e-18, 1.0e-15, 1.0e-13),
        source="manufactured bounded table",
        evidence_kind="derived",
    )
    assert coefficient.coefficient_si(RateContext(1.0)) == pytest.approx(
        1.0e-18, rel=2.0e-15, abs=0.0)
    assert coefficient.coefficient_si(RateContext(2.0)) == pytest.approx(
        1.0e-15, rel=2.0e-15, abs=0.0)
    between = coefficient.coefficient_si(RateContext(3.0))
    assert 1.0e-15 < between < 1.0e-13


@pytest.mark.parametrize("temperature", [0.99, 4.01])
def test_tabulated_temperature_rate_refuses_extrapolation(temperature):
    coefficient = ElectronTemperatureTabulatedRateCoefficient(
        electron_temperature_eV=(1.0, 4.0),
        coefficient_m3_s=(1.0e-18, 1.0e-13),
        source="manufactured bounded table",
        evidence_kind="derived",
    )
    with pytest.raises(ValueError, match="outside the tabulated rate domain"):
        coefficient.coefficient_si(RateContext(temperature))


def test_incomplete_electron_energy_ledger_fails_closed():
    species = _toy_species()[:2]
    reaction = Reaction(
        name="unknown_energy_A_to_B",
        reactants={"A": 1},
        products={"B": 1},
        kinetic_orders={"A": 1},
        rate_coefficient=ConstantRateCoefficient.from_per_second(
            2.0, source="manufactured"),
        electron_energy_loss_eV=None,
        source="manufactured unresolved energy",
    )
    network = ReactionNetwork(species=species, reactions=(reaction,))
    assert not network.has_complete_electron_energy_ledger
    with pytest.raises(ValueError, match="electron-energy ledger is incomplete"):
        network.electron_power_loss_density_W_m3(
            {"A": 3.0, "B": 7.0},
            RateContext(1.0),
        )


def test_particle_removing_collision_uses_same_cross_section_energy_moment():
    electron = Species(
        name="e", mass_amu=5.48579909065e-4, charge_number=-1,
        composition={}, role="electron", source="CODATA")
    molecule = Species(
        name="Cl2", mass_amu=70.906, charge_number=0,
        composition={"Cl": 2}, role="neutral", source="manufactured")
    atom = Species(
        name="Cl", mass_amu=35.453, charge_number=0,
        composition={"Cl": 1}, role="neutral", source="manufactured")
    anion = Species(
        name="Cl-", mass_amu=35.453, charge_number=-1,
        composition={"Cl": 1}, role="negative_ion", source="manufactured")
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, 100.0),
        cross_section_m2=(2.5e-20, 2.5e-20),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured complete attachment cross section",
        evidence_kind="derived",
    )
    reaction = Reaction(
        name="manufactured_dissociative_attachment",
        reactants={"e": 1, "Cl2": 1},
        products={"Cl": 1, "Cl-": 1},
        kinetic_orders={"e": 1, "Cl2": 1},
        rate_coefficient=coefficient,
        electron_energy_loss_eV=None,
        electron_energy_loss_moment=(
            INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT),
        source="manufactured Boltzmann moment",
    )
    network = ReactionNetwork(
        species=(electron, molecule, atom, anion),
        reactions=(reaction,),
    )
    densities = {
        "e": 2.0e16,
        "Cl2": 3.0e20,
        "Cl": 0.0,
        "Cl-": 0.0,
    }
    context = RateContext(2.0)
    expected = (
        E_CHARGE_C
        * coefficient.incident_energy_moment_eV_m3_s(context)
        * densities["e"]
        * densities["Cl2"]
    )
    assert network.has_complete_electron_energy_ledger
    assert network.electron_power_loss_density_W_m3(
        densities, context) == pytest.approx(
            expected, rel=3.0e-15, abs=0.0)


def test_incident_energy_moment_rejects_unrelated_rate_or_scalar_override():
    common = dict(
        name="invalid_moment",
        reactants={"e": 1, "Cl2": 1},
        products={"Cl": 1, "Cl-": 1},
        kinetic_orders={"e": 1, "Cl2": 1},
        electron_energy_loss_moment=(
            INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT),
        source="manufactured invalid moment",
    )
    with pytest.raises(ValueError, match="same Maxwellian cross section"):
        Reaction(
            rate_coefficient=(
                ElectronArrheniusRateCoefficient.from_cm3_per_s(
                    1.0e-8,
                    activation_eV=1.0,
                    source="unrelated fit",
                )
            ),
            electron_energy_loss_eV=None,
            **common,
        )
    with pytest.raises(ValueError, match="invalid reactor reaction"):
        Reaction(
            rate_coefficient=ElectronMaxwellianCrossSectionRateCoefficient(
                electron_energy_eV=(0.0, 100.0),
                cross_section_m2=(1.0e-20, 1.0e-20),
                threshold_eV=0.0,
                relative_uncertainty=None,
                source="manufactured complete table",
                evidence_kind="derived",
            ),
            electron_energy_loss_eV=2.0,
            **common,
        )


def test_elastic_energy_moment_uses_same_momentum_transfer_cross_section():
    electron = Species(
        name="e", mass_amu=5.48579909065e-4, charge_number=-1,
        composition={}, role="electron", source="CODATA")
    molecule = Species(
        name="Cl2", mass_amu=70.906, charge_number=0,
        composition={"Cl": 2}, role="neutral", source="manufactured")
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, 100.0),
        cross_section_m2=(2.5e-20, 2.5e-20),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured complete momentum-transfer cross section",
        evidence_kind="derived",
    )
    target_mass_kg = 70.906 * 1.66053906892e-27
    common = dict(
        reactants={"e": 1, "Cl2": 1},
        products={"e": 1, "Cl2": 1},
        kinetic_orders={"e": 1, "Cl2": 1},
        rate_coefficient=coefficient,
        electron_energy_loss_eV=None,
        elastic_target_mass_kg=target_mass_kg,
        source="manufactured elastic transfer",
    )
    exact = Reaction(
        name="manufactured_exact_elastic",
        electron_energy_loss_moment=(
            STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT),
        **common,
    )
    kemaneci = Reaction(
        name="manufactured_kemaneci_elastic",
        electron_energy_loss_moment=(
            KEMANECI_ELASTIC_ENERGY_APPROXIMATION),
        **common,
    )
    densities = {"e": 2.0e16, "Cl2": 3.0e20}
    context = RateContext(2.0)
    transfer_fraction = (
        2.0 * ELECTRON_MASS_KG * target_mass_kg
        / (ELECTRON_MASS_KG + target_mass_kg) ** 2
    )
    expected_exact = (
        transfer_fraction
        * coefficient.incident_energy_moment_eV_m3_s(context)
        * densities["e"]
        * densities["Cl2"]
    )
    assert exact.electron_energy_loss_rate_eV_m3_s(
        densities, context) == pytest.approx(
            expected_exact, rel=3.0e-15, abs=0.0)

    expected_kemaneci = (
        3.0 * context.electron_temperature_eV
        * ELECTRON_MASS_KG / target_mass_kg
        * coefficient.coefficient_si(context)
        * densities["e"]
        * densities["Cl2"]
    )
    assert kemaneci.electron_energy_loss_rate_eV_m3_s(
        densities, context) == pytest.approx(
            expected_kemaneci, rel=3.0e-15, abs=0.0)
    assert expected_kemaneci / expected_exact == pytest.approx(
        0.7500116051257543, rel=2.0e-15, abs=0.0)

    network = ReactionNetwork(
        species=(electron, molecule), reactions=(exact, kemaneci))
    assert network.has_complete_electron_energy_ledger


def test_elastic_energy_moment_rejects_missing_mass_and_changed_target():
    coefficient = ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(0.0, 100.0),
        cross_section_m2=(1.0e-20, 1.0e-20),
        threshold_eV=0.0,
        relative_uncertainty=None,
        source="manufactured momentum-transfer table",
        evidence_kind="derived",
    )
    common = dict(
        name="invalid_elastic_transfer",
        reactants={"e": 1, "Cl2": 1},
        kinetic_orders={"e": 1, "Cl2": 1},
        rate_coefficient=coefficient,
        electron_energy_loss_eV=None,
        electron_energy_loss_moment=(
            STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT),
        source="manufactured invalid elastic transfer",
    )
    with pytest.raises(ValueError, match="elastic electron-energy transfer"):
        Reaction(products={"e": 1, "Cl2": 1}, **common)
    with pytest.raises(ValueError, match="elastic electron-energy transfer"):
        Reaction(
            products={"e": 1, "Cl": 2},
            elastic_target_mass_kg=70.906 * 1.66053906892e-27,
            **common,
        )


def test_electron_power_ledger_rejects_unknown_density_species():
    reaction = Reaction(
        name="A_to_B_with_energy",
        reactants={"A": 1},
        products={"B": 1},
        kinetic_orders={"A": 1},
        rate_coefficient=ConstantRateCoefficient.from_per_second(
            2.0, source="manufactured"),
        electron_energy_loss_eV=1.0,
        source="manufactured",
    )
    network = ReactionNetwork(
        species=_toy_species()[:2], reactions=(reaction,))
    with pytest.raises(KeyError, match="unknown density species"):
        network.electron_power_loss_density_W_m3(
            {"A": 3.0, "B": 7.0, "ghost": 1.0},
            RateContext(1.0),
        )


def test_composite_electron_fit_replays_terms_and_enforces_domain():
    coefficient = ElectronCompositeRateCoefficient(
        terms=(
            ElectronAnalyticRateTerm(
                prefactor_si=3.43e-15,
                temperature_power=-1.18,
                inverse_temperature_coefficients=(-3.98,),
            ),
            ElectronAnalyticRateTerm(
                prefactor_si=3.05e-16,
                temperature_power=-1.33,
                shifted_inverse_coefficient_eV=-0.11,
                shifted_inverse_offset_eV=0.014,
            ),
        ),
        minimum_temperature_eV=0.5,
        maximum_temperature_eV=10.0,
        density_order=2.0,
        source="manufactured Kemaneci-form fit",
        source_units="m^3 s^-1; Te in eV",
        evidence_kind="published_compilation",
    )
    temperature = 2.5
    expected = (
        3.43e-15 * temperature ** -1.18 * np.exp(-3.98 / temperature)
        + 3.05e-16 * temperature ** -1.33
        * np.exp(-0.11 / (temperature + 0.014))
    )
    assert coefficient.coefficient_si(
        RateContext(temperature)) == pytest.approx(
            expected, rel=3.0e-15, abs=0.0)
    for invalid in (0.499, 10.001):
        with pytest.raises(ValueError, match="published fit domain"):
            coefficient.coefficient_si(RateContext(invalid))


def test_composite_electron_fit_supports_inverse_and_log_gaussian_terms():
    coefficient = ElectronCompositeRateCoefficient(
        terms=(
            ElectronAnalyticRateTerm(
                prefactor_si=3.28e-17,
                temperature_power=-1.12,
                inverse_temperature_coefficients=(-0.37,),
            ),
            ElectronAnalyticRateTerm(
                prefactor_si=2.86e-17,
                log_temperature_shift=0.99,
                log_temperature_width=1.06,
            ),
        ),
        minimum_temperature_eV=0.5,
        maximum_temperature_eV=10.0,
        density_order=2.0,
        source="manufactured Kemaneci-form log fit",
        source_units="m^3 s^-1; Te in eV",
        evidence_kind="published_compilation",
    )
    temperature = 3.0
    expected = (
        3.28e-17 * temperature ** -1.12 * np.exp(-0.37 / temperature)
        + 2.86e-17 * np.exp(
            -(np.log(temperature) + 0.99) ** 2 / (2.0 * 1.06 ** 2))
    )
    assert coefficient.coefficient_si(
        RateContext(temperature)) == pytest.approx(
            expected, rel=3.0e-15, abs=0.0)


def test_electron_detailed_balance_recovers_boltzmann_equilibrium():
    forward = ElectronArrheniusRateCoefficient.from_cm3_per_s(
        3.2e-8,
        activation_eV=0.4,
        temperature_power=-0.25,
        source="manufactured excitation coefficient",
    )
    coefficient = ElectronDetailedBalanceRateCoefficient(
        forward_rate_coefficient=forward,
        excitation_energy_eV=0.109,
        lower_statistical_weight=4.0,
        upper_statistical_weight=2.0,
        source="microscopic reversibility regression",
    )
    context = RateContext(2.3)
    # Independent 50-digit evaluation of the complete reverse coefficient.
    assert coefficient.coefficient_si(context) == pytest.approx(
        4.579312925020376e-14, rel=3.0e-15, abs=0.0)

    upper_to_lower_population = 0.5 * np.exp(-0.109 / 2.3)
    forward_event_rate = forward.coefficient_si(context)
    reverse_event_rate = (
        upper_to_lower_population * coefficient.coefficient_si(context))
    assert reverse_event_rate == pytest.approx(
        forward_event_rate, rel=3.0e-15, abs=0.0)
    assert coefficient.density_order == forward.density_order
    assert coefficient.source_units == forward.source_units


@pytest.mark.parametrize(
    "gap,lower_weight,upper_weight",
    [(0.0, 4.0, 2.0), (0.109, 0.0, 2.0), (0.109, 4.0, -2.0)],
)
def test_electron_detailed_balance_rejects_nonphysical_levels(
        gap, lower_weight, upper_weight):
    forward = ElectronArrheniusRateCoefficient.from_cm3_per_s(
        1.0e-8, activation_eV=0.1, source="manufactured")
    with pytest.raises(ValueError, match="invalid detailed-balance"):
        ElectronDetailedBalanceRateCoefficient(
            forward_rate_coefficient=forward,
            excitation_energy_eV=gap,
            lower_statistical_weight=lower_weight,
            upper_statistical_weight=upper_weight,
            source="manufactured",
        )


def test_gas_temperature_rate_requires_and_replays_declared_temperature():
    coefficient = GasTemperatureArrheniusRateCoefficient(
        prefactor_si=5.0e-14,
        temperature_power=-0.5,
        activation_temperature_K=0.0,
        reference_temperature_K=300.0,
        density_order=2.0,
        source="manufactured heavy-particle fit",
        source_units="m^3 s^-1; Tg in K",
        evidence_kind="published_compilation",
    )
    with pytest.raises(ValueError, match="gas temperature is required"):
        coefficient.coefficient_si(RateContext(2.0))
    assert coefficient.coefficient_si(
        RateContext(2.0, gas_temperature_K=1200.0)
    ) == pytest.approx(2.5e-14, rel=2.0e-15, abs=0.0)
