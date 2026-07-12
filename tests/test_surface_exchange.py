import numpy as np
import pytest

from petch.surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange


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
