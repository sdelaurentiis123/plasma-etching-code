import pytest
import numpy as np
from types import SimpleNamespace

from scripts.audit_zhu_npg80_moving_cr_profiles import (
    CHROMIUM_REFERENCE_DENSITY_KG_M3,
    _mask_metrics,
    _process_pool_options,
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


def test_subcell_mask_is_unresolved_but_not_relabelled_physically_exhausted():
    dx = 0.01
    mask = -np.ones((41, 41, 101), dtype=float)
    # Retain a real center interval thinner than one vertical cell.  The
    # production driver must stop its quantitative claim, but must not call
    # the remaining material physically exhausted.
    x = np.arange(mask.shape[0]) * dx
    z = np.arange(mask.shape[2]) * dx
    middle = int(np.argmin(np.abs(x - 0.2)))
    mask[middle, middle, :] = np.abs(z - 0.82) - 0.004
    geometry = SimpleNamespace(material_levelsets={2: mask}, dx=dx)

    metrics = _mask_metrics(geometry, pitch_nm=400.0)
    assert 0.0 < metrics["center_remaining_thickness_nm"] < 10.0
    assert metrics["mask_below_vertical_resolution_at_center"] is True
    assert metrics["mask_exhausted_at_center"] is False


def test_cuda_multiworker_campaign_uses_spawn_context():
    # Keep this as a structural regression: actually exercising CUDA belongs
    # to the hardware parity sentinel, while this test prevents a future
    # refactor from reintroducing the known forked-primary-context failure.
    assert _process_pool_options("cpu") == {}
    assert _process_pool_options("cuda:0")[
        "mp_context"
    ].get_start_method() == "spawn"
