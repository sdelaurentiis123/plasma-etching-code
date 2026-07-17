import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.physical_api import PhysicalProcess
from petch.reactor_coupling import (
    BoundReactorBoundaryProvider,
    ReactorBoundaryPrediction,
    ReactorBoundaryQuery,
    SurfaceFeedbackState,
    resolve_reactor_boundary,
)
from petch.feature_step_3d import make_rectangular_trench_geometry_3d


def _query(**overrides):
    values = dict(
        condition_id="base",
        tool_id="tool-A",
        recipe={"pressure_Pa": 1.333, "bias_power_W": 4000.0},
        recipe_units={"pressure_Pa": "Pa", "bias_power_W": "W"},
        wafer_position_m=(0.0, 0.075),
        process_time_s=12.0,
        substrate_temperature_K=300.0,
        provenance={"source": "manufactured coupling gate"},
    )
    values.update(overrides)
    return ReactorBoundaryQuery(**values)


def _boundary(*, predictive):
    ion = SpeciesBoundaryState(
        name="Ar+",
        charge_number=1,
        mass_amu=39.948,
        flux_m2_s=1.0e20,
        velocity_sqrt_eV=np.array([[0.0, 0.0, 10.0]]),
        weight=np.array([1.0]),
        provenance={"source": "manufactured IEAD"},
    )
    return PlasmaBoundaryState(
        species=(ion,),
        reference_plane_m=2.0e-6,
        provenance={
            "provider": "manufactured reactor",
            "supports_prediction": predictive,
        },
    )


def test_reactor_query_hash_is_canonical_and_changes_with_the_operating_point():
    first = _query()
    reordered = _query(
        recipe={"bias_power_W": 4000.0, "pressure_Pa": 1.333},
        recipe_units={"bias_power_W": "W", "pressure_Pa": "Pa"},
    )
    changed = _query(process_time_s=13.0)

    assert first.query_sha256 == reordered.query_sha256
    assert first.query_sha256 != changed.query_sha256


def test_bound_provider_refuses_recipe_or_feedback_reuse():
    query = _query()
    feedback = SurfaceFeedbackState(
        species_loss_probability={"CF2": 0.3},
        product_flux_m2_s={"COF2": 2.0e18},
        net_current_density_A_m2=0.0,
        surface_temperature_K=300.0,
        reference_area_m2=1.0e-8,
        provenance={"source": "manufactured feature feedback"},
    )
    provider = BoundReactorBoundaryProvider(
        query_sha256=query.query_sha256,
        feedback_sha256=feedback.feedback_sha256,
        boundary=_boundary(predictive=False),
        provider_name="published-deck",
        provider_version="v1",
        supports_prediction=False,
    )

    assert provider.predict_boundary(query, feedback).boundary is provider.boundary
    with pytest.raises(ValueError, match="operating point"):
        provider.predict_boundary(_query(process_time_s=13.0), feedback)
    with pytest.raises(ValueError, match="surface feedback"):
        provider.predict_boundary(query)


def test_predictive_resolution_cannot_launder_a_development_boundary():
    query = _query()
    with pytest.raises(ValueError, match="cannot promote"):
        ReactorBoundaryPrediction(
            boundary=_boundary(predictive=False),
            query_sha256=query.query_sha256,
            provider_name="bad provider",
            provider_version="v1",
            supports_prediction=True,
        )

    provider = BoundReactorBoundaryProvider(
        query.query_sha256,
        _boundary(predictive=False),
        "development provider",
        "v1",
        False,
    )
    assert resolve_reactor_boundary(
        provider, query, claim_mode="development").supports_prediction is False
    with pytest.raises(ValueError, match="predictive reactor boundary"):
        resolve_reactor_boundary(provider, query, claim_mode="predictive")


def test_resolved_boundary_feeds_the_existing_common_engine_without_an_adapter_solver():
    query = _query()
    provider = BoundReactorBoundaryProvider(
        query.query_sha256,
        _boundary(predictive=True),
        "validated reactor provider",
        "v2",
        True,
        relative_standard_uncertainty={"Ar+.flux": 0.05},
    )
    prediction = resolve_reactor_boundary(
        provider, query, claim_mode="predictive")
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=7.0,
        cell_length=7.0,
        domain_height=7.0,
        dx=1.0,
        opening_width=3.0,
        mask_thickness=2.0,
        substrate_top=4.0,
        etched_depth=2.0,
        mesh_length_unit_m=1.0e-7,
    )

    process = PhysicalProcess(
        geometry=geometry,
        boundary=prediction.boundary,
        species_role={"Ar+": "ion"},
        mechanism=object(),
        etchable_material_ids=(1,),
        duration_s=1.0e-6,
        n_steps=1,
        source_bounds=(0.0, 7.0, 0.0, 7.0),
        source_z=7.0,
    )

    assert process.boundary is prediction.boundary
    assert prediction.relative_standard_uncertainty["Ar+.flux"] == 0.05
    assert (
        process.boundary.provenance["reactor_coupling"]["query_sha256"]
        == query.query_sha256
    )
