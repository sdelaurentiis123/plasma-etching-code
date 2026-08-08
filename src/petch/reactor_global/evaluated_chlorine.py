"""Evaluated electron-collision data for predictive chlorine reactors.

This module is intentionally separate from the Lee--Lieberman reproduction
deck.  The values below are measured/evaluated cross sections, not parameters
tuned to a reactor state or an etched profile.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from .chlorine import build_lee_lieberman_chlorine_particle_network
from .network import (
    ElectronMaxwellianCrossSectionRateCoefficient,
    ElectronTabulatedCrossSectionSupport,
    ElectronTemperatureTabulatedRateCoefficient,
    Reaction,
    ReactionNetwork,
)

ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV = 12.967633
MOLECULAR_CHLORINE_TOTAL_IONIZATION_THRESHOLD_EV = 11.481

_HAYES_ATOMIC_CHLORINE_IONIZATION_ENERGY_EV = (
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0,
    45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0,
    95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0,
    135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 170.0, 180.0,
    190.0, 200.0,
)

_HAYES_ATOMIC_CHLORINE_IONIZATION_CROSS_SECTION_1E20_M2 = (
    0.00, 0.01, 0.02, 0.24, 0.52, 0.74, 1.01, 1.27, 1.50, 1.65,
    1.99, 2.34, 2.59, 2.80, 2.96, 3.16, 3.20, 3.27, 3.35, 3.35,
    3.43, 3.44, 3.47, 3.49, 3.49, 3.47, 3.44, 3.43, 3.43, 3.37,
    3.34, 3.31, 3.23, 3.20, 3.21, 3.15, 3.13, 3.07, 3.05, 3.01,
    2.97, 2.96, 2.91, 2.85, 2.81, 2.72, 2.68, 2.63,
)

_NIST_MOLECULAR_CHLORINE_TOTAL_IONIZATION_ENERGY_EV = (
    11.5, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
    60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0,
)

_NIST_MOLECULAR_CHLORINE_TOTAL_IONIZATION_CROSS_SECTION_1E20_M2 = (
    0.03, 0.11, 0.25, 0.43, 0.69, 0.99, 1.32, 1.67, 2.06, 2.47,
    3.25, 3.79, 4.17, 4.51, 4.80, 5.26, 5.49, 5.68, 5.87, 6.03,
    6.15, 6.25, 6.32, 6.33, 6.31, 6.28, 6.25, 6.22, 6.19,
)

_NIST_CL2_DISSOCIATIVE_ATTACHMENT_ENERGY_EV = (
    0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.0, 1.2,
    1.6, 2.0, 2.2, 2.6, 3.0, 3.2, 3.6, 4.0, 4.2, 4.6, 5.0,
    5.2, 5.6, 6.0, 6.2, 6.6, 7.0, 7.2, 7.6, 8.0, 8.2, 8.6,
    9.0, 9.2, 9.6, 10.0, 10.2, 10.6, 11.0, 11.2, 11.6, 11.8,
)

_NIST_CL2_DISSOCIATIVE_ATTACHMENT_CROSS_SECTION_1E20_M2 = (
    1.83, 1.04, 0.32, 0.081, 0.026, 0.013, 0.0088, 0.0065,
    0.0055, 0.0062, 0.011, 0.024, 0.032, 0.036, 0.025, 0.018,
    0.012, 0.017, 0.022, 0.033, 0.047, 0.053, 0.062, 0.062,
    0.060, 0.052, 0.039, 0.030, 0.018, 0.0091, 0.0066, 0.0053,
    0.0051, 0.0049, 0.0051, 0.0049, 0.0048, 0.0046, 0.0045,
    0.0042, 0.0041, 0.0043,
)

HAMILTON_2018_CL2_DISSOCIATION_STATES = (
    ("a_3Pi_u", 3.252),
    ("A_1Pi_u", 4.348),
    ("b_3Pi_g", 6.498),
    ("B_1Pi_g", 7.537),
    ("C_1Delta_g", 7.790),
    ("c_3Sigma_g_minus", 7.257),
    ("D_1Sigma_g_plus", 8.228),
    ("e_3Sigma_u_plus", 9.219),
)
_HAMILTON_2018_RATE_TABLE = (
    Path(__file__).with_name("data")
    / "hamilton_2018_cl2_state_maxwellian_rates.csv"
)


def nist_hayes_atomic_chlorine_ionization_rate(
        *, maximum_kernel_tail_fraction: float = 1.0e-6
) -> ElectronMaxwellianCrossSectionRateCoefficient:
    """Return the NIST-recommended measured Cl-I ionization rate provider.

    Christophorou and Olthoff reproduce the selected Hayes et al. measurements
    in their Table 25 and quote a +/-14% absolute uncertainty.  NIST ASD's
    ground-state ionization limit supplies the physical threshold.  The raw
    table's tiny sub-threshold values are retained in the source transcription
    but excluded from the physical Maxwellian integral.
    """
    return ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(
            _HAYES_ATOMIC_CHLORINE_IONIZATION_ENERGY_EV),
        cross_section_m2=tuple(
            value * 1.0e-20
            for value in (
                _HAYES_ATOMIC_CHLORINE_IONIZATION_CROSS_SECTION_1E20_M2)
        ),
        threshold_eV=ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV,
        relative_uncertainty=0.14,
        source=(
            "christophorou-olthoff-1999-cl2 Table 25; Hayes et al. "
            "measurements; NIST ASD Cl I threshold"
        ),
        evidence_kind="measured",
        maximum_kernel_tail_fraction=maximum_kernel_tail_fraction,
    )


def nist_molecular_chlorine_total_ionization_rate(
        *, maximum_kernel_tail_fraction: float = 1.0e-6
) -> ElectronMaxwellianCrossSectionRateCoefficient:
    """Return the NIST-suggested total Cl2 ionization rate provider.

    Christophorou and Olthoff's Table 12 is the average of two incompatible
    experimental datasets.  The review says their difference exceeds their
    combined quoted uncertainties, so this provider deliberately carries no
    scalar ``relative_uncertainty``.  It is evaluated measurement evidence,
    not a precision claim.

    This is a *total* positive-ion source only.  The same review says the
    relative electron-impact production of Cl2+ and Cl+ is unknown.  The
    provider therefore must not be substituted for either species-resolved
    Lee--Lieberman reaction without a separately validated branching model.
    """
    return ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(
            _NIST_MOLECULAR_CHLORINE_TOTAL_IONIZATION_ENERGY_EV),
        cross_section_m2=tuple(
            value * 1.0e-20
            for value in (
                _NIST_MOLECULAR_CHLORINE_TOTAL_IONIZATION_CROSS_SECTION_1E20_M2)
        ),
        threshold_eV=MOLECULAR_CHLORINE_TOTAL_IONIZATION_THRESHOLD_EV,
        relative_uncertainty=None,
        source=(
            "christophorou-olthoff-1999-cl2 Table 12 total ionization; "
            "evaluated average of Kurepa--Belic and Stevie--Vasile "
            "measurements"
        ),
        evidence_kind="published_compilation",
        maximum_kernel_tail_fraction=maximum_kernel_tail_fraction,
    )


def nist_cl2_dissociative_attachment_cross_section_support(
) -> ElectronTabulatedCrossSectionSupport:
    """Return NIST's evaluated Cl2 attachment table without tail closure.

    Table 16 begins at 0.05 eV and ends at 11.8 eV.  Dissociative
    attachment removes the incident electron, so its electron-fluid power
    sink requires the energy-weighted moment ``<sigma v E>`` rather than a
    guessed event threshold.  The returned object evaluates both that moment
    and the particle-rate moment only on the printed support and exposes the
    missing low/high Maxwellian kernel fractions.  It is therefore evidence
    for an energy audit, not by itself a complete reactor rate provider.

    The review adjusted the Kurepa--Belic cross section upward by 30 percent
    to reconcile swarm measurements but did not assign the resulting table a
    scalar uncertainty.  No uncertainty is invented here.
    """
    return ElectronTabulatedCrossSectionSupport(
        electron_energy_eV=(
            _NIST_CL2_DISSOCIATIVE_ATTACHMENT_ENERGY_EV),
        cross_section_m2=tuple(
            value * 1.0e-20
            for value in (
                _NIST_CL2_DISSOCIATIVE_ATTACHMENT_CROSS_SECTION_1E20_M2)
        ),
        relative_uncertainty=None,
        source=(
            "christophorou-olthoff-1999-cl2 Table 16 suggested total "
            "dissociative attachment cross section"
        ),
        evidence_kind="published_compilation",
    )


@lru_cache(maxsize=1)
def _hamilton_2018_rate_rows():
    with _HAMILTON_2018_RATE_TABLE.open(
            newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 200:
        raise RuntimeError("Hamilton Cl2 state-rate table is incomplete")
    return tuple(rows)


def hamilton_2018_cl2_state_dissociation_rates():
    """Return eight fast, state-resolved Cl2 dissociation rate providers.

    The compact table was generated from all 50,000 points of each source
    cross section and reproduces the authors' independently supplied
    Maxwellian total-rate array within 0.445% over 0.3--5 eV.  Its tabulated
    providers refuse extrapolation outside that source-stated industrial
    temperature range.
    """
    rows = _hamilton_2018_rate_rows()
    temperatures = tuple(
        float(row["electron_temperature_eV"]) for row in rows)
    return tuple(
        (
            state,
            excitation_eV,
            ElectronTemperatureTabulatedRateCoefficient(
                electron_temperature_eV=temperatures,
                coefficient_m3_s=tuple(
                    float(row[f"{state}_m3_s"]) for row in rows),
                source=(
                    "hamilton-2018-cl2-dissociation exact OPJ "
                    f"{state} cross section; analytic Maxwellian reduction"
                ),
                evidence_kind="semi_empirical",
                relative_uncertainty=None,
            ),
        )
        for state, excitation_eV
        in HAMILTON_2018_CL2_DISSOCIATION_STATES
    )


def hamilton_2018_cl2_state_dissociation_reactions():
    """Return energy-resolved ``e + Cl2 -> e + 2Cl`` reactions.

    Hamilton et al. identify all eight retained excited states as
    dissociating to two ground-state Cl atoms. Their calculated vertical
    excitation energies are the event-specific electron-energy losses.
    """
    return tuple(
        Reaction(
            name=f"e_Cl2_dissociation_{state}",
            reactants={"e": 1, "Cl2": 1},
            products={"e": 1, "Cl": 2},
            kinetic_orders={"e": 1, "Cl2": 1},
            rate_coefficient=rate,
            electron_energy_loss_eV=excitation_eV,
            source=(
                "hamilton-2018-cl2-dissociation state-resolved R-matrix "
                f"cross section and Table 2 VEE: {state}"
            ),
        )
        for state, excitation_eV, rate
        in hamilton_2018_cl2_state_dissociation_rates()
    )


def build_hamilton_dissociation_chlorine_particle_network(
) -> ReactionNetwork:
    """Replace only Lee's lumped neutral-dissociation row with Hamilton.

    All other Lee--Lieberman particle reactions are retained verbatim so this
    deck isolates one evidence upgrade. In particular, total molecular
    ionization cannot replace the two species-resolved Lee channels because
    the evaluated source does not resolve the ``Cl2+``/``Cl+`` branch.
    """
    legacy = build_lee_lieberman_chlorine_particle_network()
    retained = tuple(
        reaction for reaction in legacy.reactions
        if reaction.name != "e_Cl2_dissociation"
    )
    if len(retained) != len(legacy.reactions) - 1:
        raise RuntimeError(
            "legacy chlorine deck does not contain exactly one neutral "
            "dissociation row")
    network = ReactionNetwork(
        species=legacy.species,
        reactions=(
            *retained,
            *hamilton_2018_cl2_state_dissociation_reactions(),
        ),
    )
    network.assert_closed_conservation()
    return network
