import numpy as np
import pytest

from petch.surface_exchange import (
    SurfaceMaterialExchange, SurfaceProductPopulation, unresolved_surface_exchange,
    validate_surface_product_routing,
)


def test_unresolved_exchange_closes_removed_material_without_inventing_products():
    removed = np.array([0.0, 2.0, 5.0])
    exchange = unresolved_surface_exchange(
        removed_units_m2={"SiO2_formula_unit": removed},
        deposited_units_m2={"fluorocarbon_film_unit": np.array([1.0, 0.0, 3.0])})

    assert not exchange.product_routing_complete
    assert np.array_equal(exchange.unresolved_units_m2["SiO2_formula_unit"], removed)
    assert np.array_equal(exchange.residual_units_m2("SiO2_formula_unit"), np.zeros(3))
    assert not exchange.removed_units_m2["SiO2_formula_unit"].flags.writeable


def test_routed_physical_sputter_exchange_closes_exactly():
    removed = np.array([2.0, 5.0])
    exchange = SurfaceMaterialExchange(
        removed_units_m2={"SiO2_formula_unit": removed},
        outgoing_units_m2={"SiO2_formula_unit": removed},
        unresolved_units_m2={}, deposited_units_m2={})

    assert exchange.product_routing_complete
    assert np.array_equal(exchange.residual_units_m2("SiO2_formula_unit"), np.zeros(2))


def test_exchange_rejects_material_creation_or_untracked_outgoing_inventory():
    with pytest.raises(ValueError, match="does not close"):
        SurfaceMaterialExchange(
            removed_units_m2={"Si": np.array([1.0])},
            outgoing_units_m2={"Si": np.array([1.1])},
            unresolved_units_m2={}, deposited_units_m2={})
    with pytest.raises(ValueError, match="originate"):
        SurfaceMaterialExchange(
            removed_units_m2={"Si": np.array([1.0])},
            outgoing_units_m2={"C": np.array([1.0])},
            unresolved_units_m2={"Si": np.array([1.0])}, deposited_units_m2={})


def test_explicit_product_population_closes_outgoing_material_ledger():
    removed = np.array([2.0, 5.0])
    exchange = SurfaceMaterialExchange(
        removed_units_m2={"SiO2_formula_unit": removed},
        outgoing_units_m2={"SiO2_formula_unit": removed},
        unresolved_units_m2={}, deposited_units_m2={})
    product = SurfaceProductPopulation(
        "redeposited_SiO2", "SiO2_formula_unit", removed / 2.0,
        material_units_per_particle=2.0, mass_amu=120.168,
        angular_model="diffuse_cosine", energy_model="monoenergetic",
        energy_parameters={"energy_eV": 1.0},
        provenance={"source": "manufactured routing gate"})

    assert validate_surface_product_routing(exchange, (product,)) == (product,)
    assert product.transport_ready


def test_compound_product_closes_each_elemental_inventory():
    count = np.array([2.0, 5.0])
    exchange = SurfaceMaterialExchange(
        removed_units_m2={
            "Si_atom": count,
            "Cl_atom": 2.0 * count,
        },
        outgoing_units_m2={
            "Si_atom": count,
            "Cl_atom": 2.0 * count,
        },
        unresolved_units_m2={},
        deposited_units_m2={},
    )
    product = SurfaceProductPopulation(
        "SiCl2", "Si_atom", count, 1.0, 98.991,
        additional_source_inventories_per_particle={"Cl_atom": 2.0},
        provenance={"source": "manufactured compound-routing gate"})

    assert validate_surface_product_routing(exchange, (product,)) == (product,)
    assert np.array_equal(
        product.integrated_material_inventories_m2["Si_atom"], count)
    assert np.array_equal(
        product.integrated_material_inventories_m2["Cl_atom"], 2.0 * count)
    with pytest.raises(TypeError):
        product.additional_source_inventories_per_particle["Cl_atom"] = 1.0


def test_compound_product_rejects_duplicate_or_nonphysical_stoichiometry():
    with pytest.raises(ValueError, match="stoichiometry"):
        SurfaceProductPopulation(
            "SiCl", "Si_atom", [1.0], 1.0, 63.538,
            additional_source_inventories_per_particle={"Si_atom": 1.0})
    with pytest.raises(ValueError, match="stoichiometry"):
        SurfaceProductPopulation(
            "SiCl", "Si_atom", [1.0], 1.0, 63.538,
            additional_source_inventories_per_particle={"Cl_atom": 0.0})


def test_product_population_cannot_hide_outgoing_material_error():
    exchange = SurfaceMaterialExchange(
        removed_units_m2={"Si": np.array([2.0])},
        outgoing_units_m2={"Si": np.array([2.0])},
        unresolved_units_m2={}, deposited_units_m2={})
    product = SurfaceProductPopulation(
        "Si_fragment", "Si", np.array([1.9]), 1.0, 28.085,
        provenance={"source": "manufactured routing gate"})

    with pytest.raises(ValueError, match="does not close"):
        validate_surface_product_routing(exchange, (product,))


def test_materially_identified_product_can_expose_missing_launch_distribution():
    product = SurfaceProductPopulation(
        "Si", "Si_atom", np.array([2.0]), 1.0, 28.085,
        provenance={"source": "sourced yield without differential emission data"})

    assert not product.transport_ready
    assert product.angular_model is None
    assert product.energy_model is None


def test_product_energy_law_requires_physical_parameters():
    with pytest.raises(ValueError, match="binding and maximum"):
        SurfaceProductPopulation(
            "Si", "Si_atom", [1.0], 1.0, 28.085,
            angular_model="diffuse_cosine", energy_model="thompson")
    with pytest.raises(ValueError, match="invalid Thompson"):
        SurfaceProductPopulation(
            "Si", "Si_atom", [1.0], 1.0, 28.085,
            angular_model="diffuse_cosine", energy_model="thompson",
            energy_parameters={
                "surface_binding_energy_eV": 5.0, "maximum_energy_eV": 4.0})
