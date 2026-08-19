from dataclasses import replace

from petch import (
    EnergeticFlux,
    FluorinatedOxideSurfaceState,
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


def test_default_sio2_inventory_contract_is_backward_compatible():
    parameters = ReducedFluorinatedOxideParameters.huang_kushner_2019_reduced_projection()

    assert parameters.material_name == "SiO2"
    assert parameters.material_inventory_name == "SiO2_formula_unit"
