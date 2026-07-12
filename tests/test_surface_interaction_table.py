import numpy as np
import pytest

from petch.surface_interaction_table import (
    InteractionAxis,
    SurfaceInteractionDomainError,
    SurfaceInteractionTable,
)


def _table():
    axes = (
        InteractionAxis("ion_energy", np.array([10.0, 30.0]), "eV"),
        InteractionAxis("cosine_incidence", np.array([0.5, 1.0]), "1"),
        InteractionAxis("neutral_to_ion_flux_ratio", np.array([1.0, 100.0]), "1", "log"),
    )
    energy, cosine, ratio = np.meshgrid(
        axes[0].values, axes[1].values, axes[2].values, indexing="ij")
    log_ratio = np.log(ratio)
    etch_yield = 0.01 * energy + 0.2 * cosine + 0.03 * log_ratio
    chloride_branch = np.broadcast_to(0.25 + 0.1 * cosine, etch_yield.shape)
    silicon_branch = 1.0 - chloride_branch
    return SurfaceInteractionTable(
        material="Si(100)", incident_species=("Ar+", "Cl"), axes=axes,
        outputs={
            "etch_yield": etch_yield,
            "chloride_product_fraction": chloride_branch,
            "silicon_product_fraction": silicon_branch,
        },
        output_units={
            "etch_yield": "Si/ion",
            "chloride_product_fraction": "1",
            "silicon_product_fraction": "1",
        },
        provenance={
            "source": "manufactured interpolation gate",
            "evidence_type": "analytic",
            "intended_source": "OSTI 2589032",
        },
        standard_uncertainty={"etch_yield": 0.02},
        bounds={
            "etch_yield": (0.0, None),
            "chloride_product_fraction": (0.0, 1.0),
            "silicon_product_fraction": (0.0, 1.0),
        },
        conservation_groups={
            "emitted_product_fraction": (
                "chloride_product_fraction", "silicon_product_fraction")},
    )


def test_surface_interaction_table_replays_every_source_node_exactly():
    table = _table()
    coordinates = np.meshgrid(*[axis.values for axis in table.axes], indexing="ij")
    evaluated = table.evaluate({
        axis.name: value for axis, value in zip(table.axes, coordinates)})

    for name, expected in table.outputs.items():
        assert np.array_equal(evaluated.values[name], expected)
    assert evaluated.extrapolated_fraction == 0.0
    assert evaluated.table_fingerprint == table.fingerprint


def test_surface_interaction_table_interpolates_in_declared_physical_coordinates():
    table = _table()
    energy = np.array([15.0, 25.0]); cosine = np.array([0.6, 0.9])
    ratio = np.array([10.0, 10.0])
    evaluated = table.evaluate({
        "ion_energy": energy,
        "cosine_incidence": cosine,
        "neutral_to_ion_flux_ratio": ratio,
    })
    expected = 0.01 * energy + 0.2 * cosine + 0.03 * np.log(ratio)

    assert np.allclose(evaluated.values["etch_yield"], expected, rtol=1e-14)
    assert np.allclose(evaluated.standard_uncertainty["etch_yield"], 0.02)
    assert np.allclose(
        evaluated.values["chloride_product_fraction"]
        + evaluated.values["silicon_product_fraction"], 1.0)


def test_surface_interaction_table_refuses_silent_extrapolation_and_reports_explicit_use():
    table = _table()
    coordinates = {
        "ion_energy": np.array([20.0, 40.0]),
        "cosine_incidence": np.array([0.75, 0.75]),
        "neutral_to_ion_flux_ratio": np.array([10.0, 10.0]),
    }
    with pytest.raises(SurfaceInteractionDomainError, match="ion_energy"):
        table.evaluate(coordinates)

    evaluated = table.evaluate(coordinates, extrapolation="linear")
    assert evaluated.extrapolated_fraction == 0.5
    assert evaluated.outside_axes == ("ion_energy",)


def test_interaction_axis_snaps_only_machine_roundoff_at_validated_endpoints():
    table = _table()
    coordinates = {
        "ion_energy": 10.0,
        "cosine_incidence": 0.5,
        "neutral_to_ion_flux_ratio": np.nextafter(1.0, 0.0),
    }
    evaluated = table.evaluate(coordinates)
    exact = table.evaluate({**coordinates, "neutral_to_ion_flux_ratio": 1.0})
    assert evaluated.extrapolated_fraction == 0.0
    assert np.array_equal(evaluated.values["etch_yield"], exact.values["etch_yield"])
    with pytest.raises(SurfaceInteractionDomainError):
        table.evaluate({**coordinates, "neutral_to_ion_flux_ratio": 0.999})


def test_surface_interaction_payload_round_trip_is_bitwise_replayable():
    table = _table()
    replay = SurfaceInteractionTable.from_payload(table.to_payload())

    assert replay.fingerprint == table.fingerprint
    assert replay.to_payload() == table.to_payload()


def test_surface_interaction_table_enforces_branch_conservation_and_output_bounds():
    table = _table(); payload = table.to_payload()
    payload["outputs"]["chloride_product_fraction"][0][0][0] = 1.2
    with pytest.raises(ValueError, match="upper physical bound"):
        SurfaceInteractionTable.from_payload(payload)

    payload = table.to_payload()
    payload["outputs"]["silicon_product_fraction"][0][0][0] = 0.5
    with pytest.raises(ValueError, match="sum to one"):
        SurfaceInteractionTable.from_payload(payload)


def test_one_axis_leave_one_out_audit_reports_interpolation_error_separately():
    axis = InteractionAxis("energy", np.array([10.0, 20.0, 30.0, 40.0]), "eV")
    table = SurfaceInteractionTable(
        material="Si", incident_species=("Ar+",), axes=(axis,),
        outputs={"yield": np.array([1.0, 4.0, 9.0, 16.0])},
        output_units={"yield": "Si/ion"},
        standard_uncertainty={"yield": np.array([0.1, 0.2, 0.3, 0.4])},
        bounds={"yield": (0.0, None)},
        provenance={"source": "manufactured quadratic", "evidence_type": "analytic"})

    audit = table.leave_one_out_interpolation_audit("yield")

    assert np.array_equal(audit.coordinates, [20.0, 30.0])
    assert np.array_equal(audit.observed, [4.0, 9.0])
    assert np.array_equal(audit.predicted, [5.0, 10.0])
    assert np.array_equal(audit.absolute_error, [1.0, 1.0])
    assert audit.table_fingerprint == table.fingerprint
    assert not hasattr(audit, "combined_uncertainty")


def test_leave_one_out_audit_refuses_multiaxis_and_unknown_outputs():
    with pytest.raises(ValueError, match="exactly one"):
        _table().leave_one_out_interpolation_audit("etch_yield")
    axis = InteractionAxis("energy", np.array([10.0, 20.0, 30.0]), "eV")
    table = SurfaceInteractionTable(
        material="Si", incident_species=("Ar+",), axes=(axis,),
        outputs={"yield": np.ones(3)}, output_units={"yield": "Si/ion"},
        provenance={"source": "manufactured", "evidence_type": "analytic"})
    with pytest.raises(ValueError, match="unknown"):
        table.leave_one_out_interpolation_audit("damage")
