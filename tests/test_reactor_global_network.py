import numpy as np
import pytest

from petch.reactor_global import (
    CM3_TO_M3,
    E_CHARGE_C,
    ELECTRON_MASS_KG,
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronInverseTemperaturePolynomialRateCoefficient,
    ElectronBase10LogPolynomialRateCoefficient,
    ElectronMaxwellianCrossSectionRateCoefficient,
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
        RateContext(temperature)) == pytest.approx(expected, rel=3.0e-15)


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
