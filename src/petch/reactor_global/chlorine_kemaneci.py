"""Kemaneci 2014 detailed-chlorine source-reproduction networks.

This module reproduces the 36 non-elastic *forward* volume reactions printed
in Table 4 of ``kemaneci-2014-chlorine-global``.  It deliberately omits the
two elastic cross-section channels and the eight reverse excitation reactions
that the paper says are obtained by detailed balance but does not print.

The source's printed charge-exchange range ``(28)--(32)`` is inconsistent:
``Cl2(v=0--3)`` expands to four reactions.  The official COMSOL reproduction
also implements four.  Its raw model contains 38 forward rows (including two
elastic channels) plus eight detailed-balance reverses, for 46 volume
reaction features.  This replay therefore retains Table-4 labels 28--31 and
records label 32 as a source defect instead of inventing a fifth molecular
state.

All electron fits fail outside the paper's declared ``0.5--10 eV`` domain.
No electron-event energy is inferred from a fit exponent.  The forward-only
paper replay and the separate COMSOL nonelastic replay are particle/chemistry
verification infrastructure and intentionally fail the electron-power ledger.
"""
from __future__ import annotations

from .argon import ELECTRON_MASS_AMU
from .chlorine import (
    CHLORINE_ATOM_MASS_AMU,
    CHLORINE_MOLECULE_MASS_AMU,
)
from .network import (
    ConstantRateCoefficient,
    ElectronAnalyticRateTerm,
    ElectronCompositeRateCoefficient,
    ElectronDetailedBalanceRateCoefficient,
    GasTemperatureArrheniusRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)

_SOURCE = "kemaneci-2014-chlorine-global Table 4"
_EVIDENCE = "published_compilation"
_ELECTRON_UNITS = "m^3 s^-1; Te in eV; 0.5 <= Te <= 10"
_COMSOL_SOURCE = (
    "comsol-6.4 chlorine_global_model raw dmodel.xml; "
    "Kemaneci 2014 implementation reproduction"
)


def _term(
        prefactor: float, *, power: float = 0.0,
        inverse: tuple[float, ...] = (),
        shifted: float = 0.0, offset: float | None = None,
        log_shift: float | None = None,
        log_width: float | None = None) -> ElectronAnalyticRateTerm:
    return ElectronAnalyticRateTerm(
        prefactor_si=prefactor,
        temperature_power=power,
        inverse_temperature_coefficients=inverse,
        shifted_inverse_coefficient_eV=shifted,
        shifted_inverse_offset_eV=offset,
        log_temperature_shift=log_shift,
        log_temperature_width=log_width,
    )


def _electron_rate(
        row: str, *terms: ElectronAnalyticRateTerm,
        source: str | None = None,
) -> ElectronCompositeRateCoefficient:
    return ElectronCompositeRateCoefficient(
        terms=terms,
        minimum_temperature_eV=0.5,
        maximum_temperature_eV=10.0,
        density_order=2.0,
        source=source or f"{_SOURCE} reaction {row}",
        source_units=_ELECTRON_UNITS,
        evidence_kind=_EVIDENCE,
    )


def _gas_rate(
        row: str, prefactor: float, *, power: float = 0.0,
        activation_K: float = 0.0, density_order: float = 2.0
) -> GasTemperatureArrheniusRateCoefficient:
    units = "m^3 s^-1; Tg in K"
    if density_order == 3.0:
        units = "m^6 s^-1; Tg in K"
    return GasTemperatureArrheniusRateCoefficient(
        prefactor_si=prefactor,
        temperature_power=power,
        activation_temperature_K=activation_K,
        reference_temperature_K=300.0,
        density_order=density_order,
        source=f"{_SOURCE} reaction {row}",
        source_units=units,
        evidence_kind=_EVIDENCE,
    )


def kemaneci_2014_chlorine_species() -> tuple[Species, ...]:
    """Return the ten unique heavy states plus explicit electrons."""
    source = "kemaneci-2014-chlorine-global Table 2"
    return (
        Species(
            name="e", mass_amu=ELECTRON_MASS_AMU, charge_number=-1,
            composition={}, role="electron", source="CODATA electron mass",
            evidence_kind="measured"),
        Species(
            name="Cl2", mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=0, composition={"Cl": 2}, role="neutral",
            source=f"{source}; unqualified Cl2 is v=0",
            evidence_kind=_EVIDENCE),
        Species(
            name="Cl2(v=1)", mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=0, composition={"Cl": 2}, role="excited_neutral",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl2(v=2)", mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=0, composition={"Cl": 2}, role="excited_neutral",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl2(v=3)", mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=0, composition={"Cl": 2}, role="excited_neutral",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl", mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=0, composition={"Cl": 1}, role="neutral",
            source=f"{source}; Cl(2P3/2) ground state",
            evidence_kind=_EVIDENCE),
        Species(
            name="Cl(2P1/2)", mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=0, composition={"Cl": 1}, role="excited_neutral",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl(1P5/2)", mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=0, composition={"Cl": 1}, role="excited_neutral",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl2+", mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=1, composition={"Cl": 2}, role="positive_ion",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl+", mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=1, composition={"Cl": 1}, role="positive_ion",
            source=source, evidence_kind=_EVIDENCE),
        Species(
            name="Cl-", mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=-1, composition={"Cl": 1}, role="negative_ion",
            source=source, evidence_kind=_EVIDENCE),
    )


def build_kemaneci_2014_forward_chlorine_network() -> ReactionNetwork:
    """Return the source's 36 printed non-elastic forward reactions."""
    reactions: list[Reaction] = []

    def electron(
            label: int, name: str, reactants: dict[str, float],
            products: dict[str, float], rate: ElectronCompositeRateCoefficient
    ) -> None:
        reactions.append(Reaction(
            name=f"k{label:02d}_{name}",
            reactants=reactants,
            products=products,
            kinetic_orders=reactants,
            rate_coefficient=rate,
            electron_energy_loss_eV=None,
            source=f"{_SOURCE} reaction {label}",
        ))

    electron(1, "Cl2_dissociation", {"e": 1, "Cl2": 1},
             {"Cl": 2, "e": 1},
             _electron_rate("1", _term(
                 1.04e-13, power=-0.29, inverse=(-8.84,))))
    electron(2, "Cl2_ionization", {"e": 1, "Cl2": 1},
             {"Cl2+": 1, "e": 2},
             _electron_rate("2", _term(
                 5.12e-14, power=0.48, inverse=(-12.34,))))
    electron(3, "Cl2_dissociative_ionization", {"e": 1, "Cl2": 1},
             {"Cl": 1, "Cl+": 1, "e": 2},
             _electron_rate("3", _term(
                 2.14e-13, power=-0.07, inverse=(-25.26,))))
    electron(4, "Cl2_double_dissociative_ionization",
             {"e": 1, "Cl2": 1}, {"Cl+": 2, "e": 3},
             _electron_rate("4", _term(
                 2.27e-16, power=1.92, inverse=(-21.26,))))

    attachment_scales = (
        (5, "Cl2", 3.43e-15, 3.05e-16),
        (6, "Cl2(v=1)", 14.06e-15, 12.51e-16),
        (7, "Cl2(v=2)", 30.18e-15, 26.84e-16),
        (8, "Cl2(v=3)", 46.31e-15, 41.18e-16),
    )
    for label, molecule, first, second in attachment_scales:
        electron(
            label,
            f"{molecule}_dissociative_attachment",
            {"e": 1, molecule: 1},
            {"Cl": 1, "Cl-": 1},
            _electron_rate(
                str(label),
                _term(first, power=-1.18, inverse=(-3.98,)),
                _term(
                    second, power=-1.33, shifted=-0.11, offset=0.014),
            ),
        )

    electron(9, "Cl2_ion_pair", {"e": 1, "Cl2": 1},
             {"Cl+": 1, "Cl-": 1, "e": 1},
             _electron_rate("9", _term(
                 2.94e-16, power=0.19, inverse=(-18.79,))))
    electron(10, "Cl2_v0_to_v1", {"e": 1, "Cl2": 1},
             {"Cl2(v=1)": 1, "e": 1},
             _electron_rate("10", _term(
                 3.99e-12, power=-1.5, inverse=(-7.51, -0.0001))))
    electron(11, "Cl2_v0_to_v2", {"e": 1, "Cl2": 1},
             {"Cl2(v=2)": 1, "e": 1},
             _electron_rate(
                 "11",
                 _term(3.28e-17, power=-1.12, inverse=(-0.37,)),
                 _term(2.86e-17, log_shift=0.99, log_width=1.06)))
    electron(12, "Cl2_v0_to_v3", {"e": 1, "Cl2": 1},
             {"Cl2(v=3)": 1, "e": 1},
             _electron_rate(
                 "12",
                 _term(1.30e-17, power=-1.24, inverse=(-0.41,)),
                 _term(6.08e-18, log_shift=0.94, log_width=1.02)))
    for label, lower, upper in (
        (13, "Cl2(v=1)", "Cl2(v=2)"),
        (14, "Cl2(v=2)", "Cl2(v=3)"),
    ):
        electron(
            label, f"{lower}_to_{upper}", {"e": 1, lower: 1},
            {upper: 1, "e": 1},
            _electron_rate(
                str(label),
                _term(3.00e-16, power=-1.00, inverse=(-0.37,)),
                _term(4.61e-16, log_shift=1.04, log_width=1.10)))
    electron(15, "Cl2_v1_to_v3", {"e": 1, "Cl2(v=1)": 1},
             {"Cl2(v=3)": 1, "e": 1},
             _electron_rate(
                 "15",
                 _term(1.25e-16, power=-1.13, inverse=(-0.36,)),
                 _term(1.06e-16, log_shift=1.01, log_width=1.06)))
    electron(16, "Cl2plus_dissociative_recombination",
             {"e": 1, "Cl2+": 1}, {"Cl": 2},
             _electron_rate("16", _term(9.00e-14, power=-0.50)))
    electron(17, "Cl_excitation_2P1_2", {"e": 1, "Cl": 1},
             {"Cl(2P1/2)": 1, "e": 1},
             _electron_rate("17", _term(
                 4.55e-14, power=-0.46, inverse=(-2.01, -0.001))))
    electron(18, "Cl_excitation_1P5_2", {"e": 1, "Cl": 1},
             {"Cl(1P5/2)": 1, "e": 1},
             _electron_rate("18", _term(
                 7.03e-17, power=0.55,
                 inverse=(-2.15, -1.5, -2.05))))
    electron(19, "Cl_ionization", {"e": 1, "Cl": 1},
             {"Cl+": 1, "e": 2},
             _electron_rate("19", _term(
                 3.17e-14, power=0.53, inverse=(-13.29,))))
    electron(20, "Cl_2P1_2_ionization", {"e": 1, "Cl(2P1/2)": 1},
             {"Cl+": 1, "e": 2},
             _electron_rate("20", _term(
                 3.17e-14, power=0.53, inverse=(-13.19,))))
    electron(21, "Cl_1P5_2_ionization", {"e": 1, "Cl(1P5/2)": 1},
             {"Cl+": 1, "e": 2},
             _electron_rate("21", _term(
                 4.33e-14, power=0.55, inverse=(-0.15, -0.85))))

    reactions.append(Reaction(
        name="k22_Cl_1P5_2_radiative_decay",
        reactants={"Cl(1P5/2)": 1}, products={"Cl": 1},
        kinetic_orders={"Cl(1P5/2)": 1},
        rate_coefficient=ConstantRateCoefficient.from_per_second(
            1.0e5, source=f"{_SOURCE} reaction 22",
            evidence_kind=_EVIDENCE),
        electron_energy_loss_eV=0.0,
        source=f"{_SOURCE} reaction 22",
    ))
    electron(23, "Clminus_detachment", {"e": 1, "Cl-": 1},
             {"Cl": 1, "e": 2},
             _electron_rate("23", _term(
                 9.02e-15, power=0.92, inverse=(-4.88,))))
    electron(24, "Clminus_double_detachment", {"e": 1, "Cl-": 1},
             {"Cl+": 1, "e": 3},
             _electron_rate("24", _term(
                 3.62e-15, power=0.72, inverse=(-25.38,))))

    def heavy(
            label: int, name: str, reactants: dict[str, float],
            products: dict[str, float], rate
    ) -> None:
        reactions.append(Reaction(
            name=f"k{label:02d}_{name}", reactants=reactants,
            products=products, kinetic_orders=reactants,
            rate_coefficient=rate, electron_energy_loss_eV=0.0,
            source=f"{_SOURCE} reaction {label}",
        ))

    heavy(25, "Cl2plus_Clminus_to_3Cl", {"Cl2+": 1, "Cl-": 1},
          {"Cl": 3}, _gas_rate("25", 5.0e-14, power=-0.5))
    heavy(26, "Cl2plus_Clminus_to_Cl_Cl2", {"Cl2+": 1, "Cl-": 1},
          {"Cl": 1, "Cl2": 1}, ConstantRateCoefficient(
              value_si=5.0e-14, density_order=2.0,
              source=f"{_SOURCE} reaction 26", source_units="m^3 s^-1",
              evidence_kind=_EVIDENCE))
    heavy(27, "Clplus_Clminus_to_2Cl", {"Cl+": 1, "Cl-": 1},
          {"Cl": 2}, _gas_rate("27", 5.0e-14, power=-0.5))

    for label, molecule in enumerate(
            ("Cl2", "Cl2(v=1)", "Cl2(v=2)", "Cl2(v=3)"), start=28):
        heavy(
            label, f"{molecule}_Clplus_charge_exchange",
            {molecule: 1, "Cl+": 1}, {"Cl": 1, "Cl2+": 1},
            ConstantRateCoefficient(
                value_si=5.40e-16, density_order=2.0,
                source=f"{_SOURCE} reactions 28--32 printed range anomaly",
                source_units="m^3 s^-1", evidence_kind=_EVIDENCE))

    heavy(33, "three_body_Cl2_association", {"Cl": 2, "Cl2": 1},
          {"Cl2": 2}, _gas_rate(
              "33", 3.50e-45, activation_K=810.0, density_order=3.0))
    heavy(34, "three_body_Cl_association", {"Cl": 3},
          {"Cl2": 1, "Cl": 1}, _gas_rate(
              "34", 8.75e-46, activation_K=810.0, density_order=3.0))
    for label, upper, lower in (
        (35, "Cl2(v=1)", "Cl2"),
        (36, "Cl2(v=2)", "Cl2(v=1)"),
        (37, "Cl2(v=3)", "Cl2(v=2)"),
    ):
        heavy(
            label, f"Cl_quench_{upper}_to_{lower}",
            {"Cl": 1, upper: 1}, {"Cl": 1, lower: 1},
            _gas_rate(str(label), 1.30e-17, power=0.5))

    return ReactionNetwork(
        species=kemaneci_2014_chlorine_species(),
        reactions=tuple(reactions),
    )


def build_kemaneci_2014_comsol_nonelastic_chlorine_network(
) -> ReactionNetwork:
    """Return the official COMSOL 6.4 model's 44 nonelastic reactions.

    This is an implementation-reproduction mode, not an evaluated-physics
    deck.  It intentionally preserves COMSOL's unit statistical-weight ratio,
    its absolute Figure-10 atomic energies used as de-excitation gaps, and its
    reaction-20 ``13.29 eV`` fit parameter.  Those choices are explicit so a
    successful replay cannot be mistaken for independent physical validation.
    """
    forward_network = build_kemaneci_2014_forward_chlorine_network()
    reactions = list(forward_network.reactions)

    reaction20_index = next(
        index for index, reaction in enumerate(reactions)
        if reaction.name.startswith("k20_"))
    paper_reaction20 = reactions[reaction20_index]
    comsol_reaction20_rate = _electron_rate(
        "20",
        _term(3.17e-14, power=0.53, inverse=(-13.29,)),
        source=f"{_COMSOL_SOURCE}; forward reaction 20",
    )
    reactions[reaction20_index] = Reaction(
        name=paper_reaction20.name,
        reactants=paper_reaction20.reactants,
        products=paper_reaction20.products,
        kinetic_orders=paper_reaction20.kinetic_orders,
        rate_coefficient=comsol_reaction20_rate,
        electron_energy_loss_eV=None,
        source=f"{_COMSOL_SOURCE}; forward reaction 20",
    )

    forward_by_label = {
        int(reaction.name[1:3]): reaction for reaction in reactions
    }
    reverse_rows = (
        (10, "Cl2_v1_to_v0", 0.07),
        (11, "Cl2_v2_to_v0", 0.14),
        (12, "Cl2_v3_to_v0", 0.21),
        (13, "Cl2_v2_to_v1", 0.07),
        (14, "Cl2_v3_to_v2", 0.07),
        (15, "Cl2_v3_to_v1", 0.14),
        (17, "Cl_2P1_2_to_ground", 1.35),
        (18, "Cl_1P5_2_to_ground", 10.17),
    )
    for label, name, gap_eV in reverse_rows:
        forward = forward_by_label[label]
        reactions.append(Reaction(
            name=f"k{label:02d}r_{name}_comsol",
            reactants=forward.products,
            products=forward.reactants,
            kinetic_orders=forward.products,
            rate_coefficient=ElectronDetailedBalanceRateCoefficient(
                forward_rate_coefficient=forward.rate_coefficient,
                excitation_energy_eV=gap_eV,
                lower_statistical_weight=1.0,
                upper_statistical_weight=1.0,
                source=(
                    f"{_COMSOL_SOURCE}; reverse of reaction {label}; "
                    "unit statistical-weight ratio"
                ),
                evidence_kind=_EVIDENCE,
            ),
            electron_energy_loss_eV=-gap_eV,
            source=f"{_COMSOL_SOURCE}; reverse of reaction {label}",
        ))

    return ReactionNetwork(
        species=forward_network.species,
        reactions=tuple(reactions),
    )
