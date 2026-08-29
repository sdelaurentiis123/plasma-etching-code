from scripts.audit_zhu_npg80_conditional_profiles import _hash, _load
from scripts.audit_zhu_npg80_gds_square_profiles import (
    MODEL_REVISION,
    PREREGISTRATION,
    _job_spec,
    _jobs,
    _physics_preregistration,
)


def test_exact_gds_profile_board_is_bound_to_received_layout():
    document = _load(PREREGISTRATION)

    assert document["geometry_source"]["sha256"] == (
        "1378b31c6b206a5b62c0254979c0fce36219dc625f716cbd00f858b53f36832b"
    )
    assert document["exact_layout_geometry"]["pitch_nm"] == 350.0
    assert document["exact_layout_geometry"]["square_width_nm"] == [
        105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0,
        145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 185.0, 190.0,
        195.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0,
    ]
    assert document["exact_layout_geometry"][
        "mask_polarity_confirmed_by_operator"] is False


def test_exact_layout_adapts_without_mutating_preregistration():
    document = _load(PREREGISTRATION)
    adapted = _physics_preregistration(document)

    assert "inferred_geometry_board" not in document
    assert adapted["inferred_geometry_board"]["pitch_nm"] == 350.0
    assert adapted["inferred_geometry_board"]["target_layout_confirmed"] is True
    assert adapted["inferred_geometry_board"]["evidence_class"] == (
        "operator_supplied_exact_GDSII"
    )


def test_exact_gds_board_has_all_184_frozen_square_trajectories():
    jobs, rates, selectivities, scenarios = _jobs(smoke=False)

    assert len(rates) == 2
    assert len(selectivities) == 2
    assert len(scenarios) == 4
    assert len(jobs) == 23 * 2 * 4
    assert {job[0] for job in jobs} == set(
        _load(PREREGISTRATION)["exact_layout_geometry"]["square_width_nm"]
    )
    assert {job[5] for job in jobs} == {10.0}


def test_exact_gds_cache_spec_binds_revision_preregistration_and_gds():
    job = _jobs(smoke=False)[0][0]
    spec = _job_spec(job)

    assert spec["model_revision"] == MODEL_REVISION
    assert spec["preregistration_sha256"] == _hash(PREREGISTRATION)
    assert spec["gds_sha256"] == _load(PREREGISTRATION)["geometry_source"][
        "sha256"]
    assert spec["mesh_spacing_nm"] == 10.0
