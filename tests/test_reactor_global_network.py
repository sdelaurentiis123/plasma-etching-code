import numpy as np
import pytest

from petch.reactor_global import (
    CM3_TO_M3,
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
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
