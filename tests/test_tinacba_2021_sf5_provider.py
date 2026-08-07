from pathlib import Path

import numpy as np
import pytest

from petch.interaction_data import load_tinacba_2021_sf5_tables
from petch.material_mechanism_3d import MaterialMechanismRouter3D
from petch.surface_interaction_table import SurfaceInteractionDomainError
from petch.surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    ParameterEvidence,
    SurfaceFluxes,
)
from petch.tabulated_chemistry import TabulatedNormalIonRemovalMechanism


ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "data"
    / "experimental"
    / "tinacba_2021"
    / "figure8_sf5_md_experiment.csv"
)
TARGET_DENSITY_M3 = {
    "Si": 5.13992e28,
    "SiO2": 2.461176470588235e28,
}


def _tables():
    return load_tinacba_2021_sf5_tables(DATA)


def _evidence(material):
    return ParameterEvidence(
        (
            "Tinacba 2021 Figures 8/10: 2000 eV MD yield divided by "
            f"printed depth-per-dose slope for {material}"
        ),
        "source_derived_model_film_density",
        supports_prediction_within_declared_domain=True,
    )


def _mechanisms():
    tables = _tables()
    return {
        1: TabulatedNormalIonRemovalMechanism(
            tables.silicon, TARGET_DENSITY_M3["Si"], _evidence("Si")
        ),
        2: TabulatedNormalIonRemovalMechanism(
            tables.silicon_dioxide,
            TARGET_DENSITY_M3["SiO2"],
            _evidence("SiO2"),
        ),
    }


def test_sf5_tables_are_checksum_bound_and_replay_all_md_nodes():
    tables = _tables()
    energy = np.asarray([150, 300, 500, 1000, 1500, 2000], dtype=float)
    assert np.array_equal(
        tables.silicon.evaluate({"ion_energy": energy}).values[
            "target_removal_yield"
        ],
        [2.0207, 2.3523, 2.8497, 4.2487, 5.8238, 6.4249],
    )
    assert np.array_equal(
        tables.silicon_dioxide.evaluate({"ion_energy": energy}).values[
            "target_removal_yield"
        ],
        [0.8187, 1.3575, 1.6995, 2.1451, 2.6839, 3.3472],
    )
    assert tables.silicon.provenance["validation"][
        "beam_overlap_energies_eV"
    ] == [150, 2000]


def test_sf5_atomistic_provider_runs_through_common_material_router():
    mechanisms = _mechanisms()
    router = MaterialMechanismRouter3D(
        mechanisms,
        provenance={
            1: {"source": "Tinacba 2021 SF5+ MD", "material": "Si"},
            2: {"source": "Tinacba 2021 SF5+ MD", "material": "SiO2"},
        },
    )
    material = np.asarray([1, 1, 2, 2])
    state = router.initial_state_by_material(material)
    events = FaceResolvedEnergeticFlux(
        "SF5+",
        4,
        event_face=np.arange(4),
        event_flux_m2_s=np.full(4, 1.0e20),
        event_energy_eV=np.asarray([150.0, 2000.0, 150.0, 2000.0]),
        event_cosine_incidence=np.ones(4),
    )
    result = router.advance_by_material(
        state, SurfaceFluxes({}, (events,)), 1.0, material
    )

    expected_yield = np.asarray([2.0207, 6.4249, 0.8187, 3.3472])
    expected_density = np.asarray([
        TARGET_DENSITY_M3["Si"],
        TARGET_DENSITY_M3["Si"],
        TARGET_DENSITY_M3["SiO2"],
        TARGET_DENSITY_M3["SiO2"],
    ])
    assert np.allclose(
        result.etch_velocity_m_s, 1.0e20 * expected_yield / expected_density
    )
    assert np.allclose(
        result.material_exchange.unresolved_units_m2["Si_atom"],
        [2.0207e20, 6.4249e20, 0.0, 0.0],
        rtol=1.0e-15,
        atol=0.0,
    )
    assert np.allclose(
        result.material_exchange.unresolved_units_m2["SiO2_formula"],
        [0.0, 0.0, 0.8187e20, 3.3472e20],
        rtol=1.0e-15,
        atol=0.0,
    )
    assert not result.material_exchange.product_routing_complete
    assert result.product_populations == ()
    assert result.validity.parameter_evidence_supports_prediction


def test_sf5_provider_depth_per_dose_matches_atomistic_board_conversion():
    mechanisms = _mechanisms()
    reference_flux = 1.0e20
    expected_depth_nm = {
        (1, 150.0): 3.931384146056748,
        (1, 2000.0): 12.5,
        (2, 150.0): 3.32645793499044,
        (2, 2000.0): 13.6,
    }
    for material_id, mechanism in mechanisms.items():
        for energy in (150.0, 2000.0):
            ions = EnergeticFlux(
                "SF5+",
                reference_flux,
                np.asarray([energy]),
                np.asarray([1.0]),
                np.asarray([1.0]),
            )
            result = mechanism.advance(
                mechanism.initial_state(),
                SurfaceFluxes({}, (ions,)),
                1.0,
            )
            assert np.isclose(
                float(result.etch_velocity_m_s) * 1.0e9,
                expected_depth_nm[(material_id, energy)],
                rtol=0.0,
                atol=1.0e-13,
            )
            assert np.array_equal(
                result.material_exchange.residual_units_m2(
                    mechanism.target_material_unit
                ),
                0.0,
            )


def test_sf5_provider_refuses_unmeasured_dimensions_and_extrapolation():
    mechanism = _mechanisms()[1]
    state = mechanism.initial_state()

    def ions(*, energy=150.0, cosine=1.0, name="SF5+"):
        return EnergeticFlux(
            name,
            1.0e20,
            np.asarray([energy]),
            np.asarray([cosine]),
            np.asarray([1.0]),
        )

    with pytest.raises(ValueError, match="normal-incidence atomistic board"):
        mechanism.advance(
            state, SurfaceFluxes({}, (ions(cosine=0.99),)), 1.0
        )
    with pytest.raises(ValueError, match="no atomistic target-removal"):
        mechanism.advance(
            state, SurfaceFluxes({"F": 1.0e20}, (ions(),)), 1.0
        )
    with pytest.raises(ValueError, match="no atomistic target-removal"):
        mechanism.advance(
            state, SurfaceFluxes({}, (ions(name="SF4+"),)), 1.0
        )
    with pytest.raises(SurfaceInteractionDomainError, match="ion_energy"):
        mechanism.advance(
            state, SurfaceFluxes({}, (ions(energy=100.0),)), 1.0
        )
    omissions = mechanism.validity(SurfaceFluxes({})).known_model_form_omissions
    assert any("S-Si" in item for item in omissions)
    assert any("product identities" in item for item in omissions)
