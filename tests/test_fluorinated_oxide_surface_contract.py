from dataclasses import replace

from petch import (
    EnergeticFlux,
    FluorinatedOxideSurfaceState,
    ParameterEvidence,
    ReducedFluorinatedOxideMechanism,
    ReducedFluorinatedOxideParameters,
    SurfaceFluxes,
)


def test_generic_oxide_contract_relabels_material_without_transferring_coefficients():
    source = ReducedFluorinatedOxideParameters.huang_kushner_2019_reduced_projection()
    parameters = replace(
        source,
        material_name="TiO2 manufactured plumbing test",
        material_inventory_name="TiO2_formula_unit",
    )
    mechanism = ReducedFluorinatedOxideMechanism(parameters)
    fluxes = SurfaceFluxes({}, (
        EnergeticFlux("Ar+", 1.0e18, [140.0], [1.0], [1.0]),
    ))
    result = mechanism.advance(
        FluorinatedOxideSurfaceState.bare(), fluxes, 1.0, strict=False
    )

    assert mechanism.provenance["material"] == {
        "name": "TiO2 manufactured plumbing test",
        "inventory_name": "TiO2_formula_unit",
    }
    assert "TiO2_formula_unit" in result.material_exchange.removed_units_m2
    assert "SiO2_formula_unit" not in result.material_exchange.removed_units_m2
    # Relabeling changes bookkeeping only. The source evidence remains
    # nonpredictive for TiO2 and must not silently become a target deck.
    assert result.validity.parameter_evidence_supports_prediction is False
    assert "polymer_bulk_unit_density_m3" not in mechanism.provenance["parameters"]


def test_default_sio2_inventory_contract_is_backward_compatible():
    parameters = ReducedFluorinatedOxideParameters.huang_kushner_2019_reduced_projection()

    assert parameters.material_name == "SiO2"
    assert parameters.material_inventory_name == "SiO2_formula_unit"
    assert parameters.polymer_bulk_unit_density_m3 is None


def test_evidence_bearing_passivation_inventory_drives_physical_growth():
    source = ReducedFluorinatedOxideParameters.huang_kushner_2019_reduced_projection()
    evidence = dict(source.evidence)
    evidence["polymer_bulk_unit_density_m3"] = ParameterEvidence(
        "manufactured passivation-growth conservation test",
        "manufactured",
        supports_prediction_within_declared_domain=False,
    )
    density = 7.5e28
    parameters = replace(
        source,
        material_name="TiO2 manufactured passivation test",
        material_inventory_name="TiO2_formula_unit",
        polymer_bulk_unit_density_m3=density,
        evidence=evidence,
    )
    mechanism = ReducedFluorinatedOxideMechanism(parameters)
    initial = FluorinatedOxideSurfaceState(
        complex_fraction=1.0,
        polymer_units_m2=0.0,
        activated_complex_fraction=1.0,
    )
    duration_s = 2.0
    result = mechanism.advance(
        initial,
        SurfaceFluxes({"CF2": 1.0e21}),
        duration_s,
    )

    assert result.deposited_polymer_units_m2 > 0.0
    assert mechanism.provenance["parameters"][
        "polymer_bulk_unit_density_m3"
    ] == density
    assert result.removed_polymer_units_m2 == 0.0
    assert result.etch_velocity_m_s == 0.0
    assert result.normal_growth_velocity_m_s == (
        result.deposited_polymer_units_m2 / density / duration_s
    )
    assert result.state.polymer_units_m2 == result.deposited_polymer_units_m2
    assert result.validity.parameter_evidence_supports_prediction is False
