import pytest
import numpy as np
from types import SimpleNamespace

from scripts.audit_zhu_npg80_moving_cr_profiles import (
    ANALOG_BOARD,
    CHROMIUM_REFERENCE_DENSITY_KG_M3,
    PREREGISTRATION,
    REACTOR_DOSE,
    _load,
    _mask_metrics,
    _process_pool_options,
    _run_trajectory,
    _router,
    _scenario_inputs,
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


def test_oxford_thin_cr_first_steps_do_not_create_subcell_material_island():
    """Pin the exact geometry path that failed the first production board.

    At 10 nm resolution, the first Cr-mask redistance repair used to expose a
    second one-node component.  That is a discretization artifact, not mask
    exhaustion and not a permissible topology-change endpoint.
    """
    preregistration = _load(PREREGISTRATION)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_DOSE)
    scenario = _scenario_inputs(preregistration, reactor)[0]
    rates = (
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"]),
    )

    profiles = _run_trajectory(
        width_nm=80.0,
        scenario=scenario,
        rates_nm_min=rates,
        selectivity=14.0,
        duration_s=6.0,
        dx_nm=10.0,
        preregistration=preregistration,
    )

    assert len(profiles) == 2
    assert [item["accepted_profile_steps"] for item in profiles] == [1, 2]
    assert all(
        item["terminal_reason"] == "requested_duration"
        for item in profiles
    )


def test_interior_minimum_thickness_ends_trajectory_before_mesh_degeneracy():
    """The v3 mask-exhaustion guard: the trajectory must end as a declared
    physical event when the Cr mask thins below one vertical cell anywhere
    inside the footprint interior, never by feeding marching cubes a
    degenerate sliver (the second production-board failure, w80/low/tail0)."""
    from scripts.audit_zhu_npg80_moving_cr_profiles import (
        MASK_MATERIAL,
        _geometry,
        _mask_interior_minimum_thickness_nm,
    )
    preregistration = _load(PREREGISTRATION)
    geometry, pitch_nm = _geometry(
        width_nm=80.0, dx_nm=10.0, preregistration=preregistration)

    initial = _mask_interior_minimum_thickness_nm(
        geometry, pitch_nm=pitch_nm, width_nm=80.0)
    # full 45 nm mask everywhere in the interior at t=0
    assert initial == pytest.approx(45.0, abs=2.0)

    # carve one interior column of the Cr level set to a sub-cell sliver
    mask = np.array(geometry.material_levelsets[MASK_MATERIAL], dtype=float)
    dx = float(geometry.dx)
    coordinate = np.arange(mask.shape[0]) * dx
    center = 0.5 * pitch_nm * 1.0e-3
    middle = int(np.argmin(np.abs(coordinate - center)))
    column = mask[middle + 1, middle + 1, :]
    z = np.arange(len(column)) * dx
    bottom = float(np.min(z[column >= 0.0]))
    mask[middle + 1, middle + 1, :] = np.minimum(
        column, (bottom + 0.5 * dx) - z)
    carved_geometry = SimpleNamespace(
        material_levelsets={MASK_MATERIAL: mask}, dx=dx)

    carved = _mask_interior_minimum_thickness_nm(
        carved_geometry, pitch_nm=pitch_nm, width_nm=80.0)
    assert carved < 10.0
    assert carved < initial
