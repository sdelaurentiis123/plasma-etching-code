import numpy as np
import pytest

from petch.neutral_radiosity_3d import DiffuseFormFactors3D
from petch.surface_exchange import SurfaceProductPopulation
from petch.surface_product_redeposition_3d import (
    SurfaceProductRedepositionContract3D,
    SurfaceProductRedepositionLaw3D,
    transport_surface_product_redeposition_3d,
)


def _law(sticking, deposited_material_id=1):
    return SurfaceProductRedepositionLaw3D(
        population_name="Si", deposited_material_id=deposited_material_id,
        sticking_probability_by_material=sticking,
        bulk_material_unit_density_m3=5e28,
        parameter_sources={
            "sticking_probability_by_material": "manufactured sticking gate",
            "bulk_material_unit_density_m3": "manufactured density gate",
        },
        parameter_bounds={
            "sticking_probability_by_material": (0.0, 1.0),
            "bulk_material_unit_density_m3": (4e28, 6e28),
        })


def _population():
    return SurfaceProductPopulation(
        "Si", "Si_atom", integrated_particle_count_m2=np.array([2.0, 0.0]),
        material_units_per_particle=1.0, mass_amu=28.085,
        angular_model="diffuse_cosine", energy_model="thompson",
        energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0},
        provenance={"source": "manufactured redeposition gate"})


def _factors():
    return DiffuseFormFactors3D(
        face_count=2, source_face=np.array([0, 1]), target_face=np.array([1, 0]),
        transfer_fraction=np.array([0.4, 0.2]), escape_fraction=np.array([0.6, 0.8]),
        rays_per_face=8)


def test_same_material_redeposition_closes_emitted_deposited_escaped_ledger():
    result = transport_surface_product_redeposition_3d(
        (_population(),), duration_s=2.0, face_area_m2=np.array([1.0, 2.0]),
        form_factors=_factors(), face_material_id=np.array([1, 1]),
        evolving_face_mask=np.array([True, True]),
        contract=SurfaceProductRedepositionContract3D((_law({1: 1.0}),)))

    assert np.allclose(result.deposited_units_m2["Si_atom"], [0.0, 0.4])
    assert np.allclose(result.normal_growth_velocity_m_s, [0.0, 4e-30])
    assert np.isclose(
        result.emitted_material_units_s["Si"],
        result.deposited_material_units_s["Si"]
        + result.escaped_material_units_s["Si"])
    assert result.maximum_relative_balance_error < 1e-13


def test_redeposition_refuses_cross_material_film_and_unbounded_parameters():
    contract = SurfaceProductRedepositionContract3D((_law({1: 1.0, 2: 1.0}),))
    with pytest.raises(ValueError, match="cross-material coating"):
        transport_surface_product_redeposition_3d(
            (_population(),), 2.0, np.array([1.0, 2.0]), _factors(),
            np.array([1, 2]), np.array([True, True]), contract)

    with pytest.raises(ValueError, match="inside their declared bounds"):
        SurfaceProductRedepositionLaw3D(
            "Si", 1, {1: 0.9}, 5e28,
            parameter_sources={
                "sticking_probability_by_material": "source",
                "bulk_material_unit_density_m3": "source"},
            parameter_bounds={
                "sticking_probability_by_material": (0.0, 0.5),
                "bulk_material_unit_density_m3": (4e28, 6e28)})
