import csv
from pathlib import Path

import numpy as np
import pytest

from petch.gray_argon_fluorine_sio2 import (
    GrayArFSiO2Mechanism,
    GrayArFSiO2Parameters,
    GrayArFSiO2State,
)
from petch.surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)


ROOT = Path(__file__).parents[1]
TABLE = (
    ROOT / "data/surface_interactions/gray_1993/"
    "table5_10_ar_f_sio2_model_parameters.csv"
)


def _ions(flux=1e20, *, energy=150.0, cosine=1.0):
    return EnergeticFlux(
        "Ar+", flux, np.atleast_1d(energy), np.atleast_1d(cosine),
        np.ones(np.atleast_1d(energy).shape),
    )


def _source_law(energy, ratio):
    p0 = 0.0139 * max(np.sqrt(energy) - np.sqrt(18.0), 0.0)
    beta2 = 0.053 * max(np.sqrt(energy) - np.sqrt(4.0), 0.0)
    branch = 0.007 * np.sqrt(energy)
    theta = 0.02 * ratio / (
        0.02 * ratio + 2.0 * beta2 * (1.0 + branch))
    total = p0 * (1.0 - theta) + beta2 * (1.0 + branch) * theta
    return p0, beta2, branch, theta, total


def test_equations_5_30_5_31_and_sio2_product_branches_replay_exactly():
    mechanism = GrayArFSiO2Mechanism()
    ion_flux = 2e20
    ratios = np.array([20.0, 50.0, 200.0])
    result = mechanism.advance(
        mechanism.initial_state(ratios.shape),
        SurfaceFluxes(
            {"F": ion_flux * ratios},
            (_ions(np.full(ratios.shape, ion_flux)),),
        ),
        2.0,
    )
    expected = np.array([_source_law(150.0, ratio) for ratio in ratios])

    assert np.allclose(result.fluorinated_fraction, expected[:, 3])
    observed_yield = (
        result.etch_velocity_m_s * 2.2e28 / ion_flux)
    assert np.allclose(observed_yield, expected[:, 4])
    assert np.allclose(
        result.sif2_removed_formula_units_m2
        / result.sif4_removed_formula_units_m2,
        expected[:, 2],
    )
    assert np.max(
        result.steady_site_balance_abs_residual_f_atoms_m2_s
    ) / ion_flux < 2e-15


def test_sio2_and_f_atoms_close_with_source_product_stoichiometry():
    result = GrayArFSiO2Mechanism().advance(
        GrayArFSiO2State.bare(),
        SurfaceFluxes({"F": 5e22}, (_ions(1e20, energy=350.0),)),
        3.0,
    )
    removed = (
        result.physical_removed_formula_units_m2
        + result.sif2_removed_formula_units_m2
        + result.sif4_removed_formula_units_m2
    )
    assert np.array_equal(
        result.material_exchange.removed_units_m2["SiO2_formula"], removed)
    assert np.array_equal(
        result.material_exchange.outgoing_units_m2["SiO2_formula"], removed)
    assert np.allclose(
        result.consumed_f_atoms_m2,
        2.0 * result.sif2_removed_formula_units_m2
        + 4.0 * result.sif4_removed_formula_units_m2,
    )
    assert result.consumed_f_atoms_m2 <= result.incident_f_atoms_m2
    assert result.material_exchange.product_routing_complete
    assert {
        item.name for item in result.product_populations
    } == {
        "SiO2_physical",
        "SiF2_plus_O2_bundle",
        "SiF4_plus_O2_bundle",
    }
    assert all(not item.transport_ready for item in result.product_populations)


def test_joint_energy_distribution_is_integrated_before_sio2_site_balance():
    mechanism = GrayArFSiO2Mechanism()
    population = EnergeticFlux(
        "Ar+", 1e20,
        np.array([20.0, 500.0]),
        np.ones(2),
        np.array([0.25, 0.75]),
    )
    result = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes({"F": 5e22}, (population,)),
        1.0,
    )
    p = np.array([_source_law(energy, 0.0)[0]
                  for energy in [20.0, 500.0]])
    beta = np.array([_source_law(energy, 0.0)[1]
                     for energy in [20.0, 500.0]])
    branch = np.array([_source_law(energy, 0.0)[2]
                       for energy in [20.0, 500.0]])
    physical_capacity = 1e20 * np.dot([0.25, 0.75], p)
    chemical_capacity = 1e20 * np.dot(
        [0.25, 0.75], beta * (1.0 + branch))
    theta = 0.02 * 5e22 / (0.02 * 5e22 + 2.0 * chemical_capacity)
    expected_rate = (
        physical_capacity * (1.0 - theta) + chemical_capacity * theta)

    assert result.fluorinated_fraction == pytest.approx(theta)
    assert result.etch_velocity_m_s * 2.2e28 == pytest.approx(expected_rate)


def test_face_event_measure_preserves_sio2_energy_without_mean_bias():
    mechanism = GrayArFSiO2Mechanism()
    events = FaceResolvedEnergeticFlux(
        "Ar+", 2,
        event_face=np.array([0, 0, 1]),
        event_flux_m2_s=np.array([1e20, 2e20, 3e20]),
        event_energy_eV=np.array([20.0, 500.0, 150.0]),
        event_cosine_incidence=np.ones(3),
    )
    result = mechanism.advance(
        mechanism.initial_state((2,)),
        SurfaceFluxes({"F": np.array([1e23, 1e23])}, (events,)),
        0.0,
    )
    assert result.etch_velocity_m_s.shape == (2,)
    assert result.etch_velocity_m_s[0] != result.etch_velocity_m_s[1]
    assert np.all(result.fluorinated_fraction > 0.0)


def test_table5_10_and_continuous_sio2_laws_are_not_conflated():
    with TABLE.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 6
    energies = np.array([float(row["ion_energy_eV"]) for row in rows])
    printed_branch = np.array([
        float(row["fragment_branch_b"]) for row in rows])
    printed_beta = np.array([
        float(row["constant_s0_beta2"]) for row in rows])
    law_branch = 0.007 * np.sqrt(energies)
    law_beta = 0.053 * np.maximum(np.sqrt(energies) - 2.0, 0.0)

    assert np.max(np.abs(printed_branch - law_branch)) < 8e-4
    # The continuous regression is close but not identical to the individual
    # parenthesized Table-5.10 values; the 250 eV point is the largest
    # residual from Gray's global square-root fit.
    assert np.max(np.abs(printed_beta - law_beta) / printed_beta) == (
        pytest.approx(0.22000596657436755, abs=1e-12))


def test_sio2_scope_refuses_angular_species_energy_and_low_ratio_transfer():
    mechanism = GrayArFSiO2Mechanism()
    state = mechanism.initial_state()
    with pytest.raises(ValueError, match="normal incidence"):
        mechanism.advance(
            state,
            SurfaceFluxes(
                {"F": 5e22},
                (_ions(cosine=np.cos(np.deg2rad(1.0))),)),
            1.0,
        )
    with pytest.raises(ValueError, match="20--2000"):
        mechanism.advance(
            state,
            SurfaceFluxes({"F": 5e22}, (_ions(energy=2500.0),)),
            1.0,
        )
    with pytest.raises(ValueError, match="R>15"):
        mechanism.advance(
            state,
            SurfaceFluxes({"F": 1e21}, (_ions(1e20),)),
            1.0,
        )
    with pytest.raises(ValueError, match="no declared"):
        mechanism.advance(
            state,
            SurfaceFluxes(
                {"F": 5e22},
                (EnergeticFlux(
                    "CF3+", 1e20, [150.0], [1.0], [1.0]),),
            ),
            1.0,
        )


def test_sio2_provenance_does_not_call_beam_regression_first_principles():
    mechanism = GrayArFSiO2Mechanism()
    provenance = mechanism.provenance
    assert "not an elementary first-principles potential" in provenance["claim"]
    assert provenance["declared_domain"]["incidence"] == "normal only"
    assert all(
        item.supports_prediction_within_declared_domain
        for item in GrayArFSiO2Parameters.beam_energy_law().evidence.values()
    )
    with pytest.raises(TypeError):
        mechanism.parameters.evidence["fragment_branch"] = None
