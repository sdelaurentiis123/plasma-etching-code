import pytest

from scripts.audit_zhu_npg80_moving_cr_profiles import (
    CHROMIUM_REFERENCE_DENSITY_KG_M3,
    _router,
    chromium_atom_density_m3,
)


def test_chromium_atom_density_is_dimensional_and_validated():
    assert chromium_atom_density_m3() == pytest.approx(8.328e28, rel=2.0e-4)
    with pytest.raises(ValueError):
        chromium_atom_density_m3(0.0)


def test_moving_mask_router_preserves_conditional_selectivity():
    router = _router(
        scenario_name="ions",
        tio2_rate_nm_min=42.0,
        selectivity=14.0,
        density_kg_m3=3700.0,
    )
    tio2 = router.mechanisms[1].parameters
    chromium = router.mechanisms[2].parameters

    assert (
        tio2.blanket_removal_velocity_m_s
        / chromium.blanket_removal_velocity_m_s
    ) == pytest.approx(14.0)
    assert chromium.bulk_material_unit_density_m3 == pytest.approx(
        chromium_atom_density_m3(CHROMIUM_REFERENCE_DENSITY_KG_M3)
    )
    assert router.provenance["materials"]["2"]["evidence"]["source"] == (
        "janissen-2016-tio2-rie"
    )
