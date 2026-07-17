import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.feature_step_3d import (
    _face_material_ids,
    _surface_mesh_fingerprint,
    make_rectangular_trench_geometry_3d,
)
from petch.threed import extract_mesh_3d


SCRIPT = Path(__file__).parents[1] / "scripts" / "krueger_2024_multiresolution_audit.py"
SPEC = importlib.util.spec_from_file_location("krueger_2024_multiresolution_audit", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _manufactured_checkpoint_state():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.04,
        cell_length=0.02,
        domain_height=0.10,
        dx=0.005,
        opening_width=0.02,
        mask_thickness=0.02,
        substrate_top=0.05,
        etched_depth=0.01,
    )
    verts, faces, centroids, _ = extract_mesh_3d(geometry.phi, geometry.dx)
    material = _face_material_ids(centroids, geometry)
    active = np.flatnonzero(np.isin(material, AUDIT.ETCHABLE))
    mechanism = build_krueger_2024_material_router_3d(**AUDIT.CALIBRATION)
    state = mechanism.initial_state_by_material(material[active])
    fingerprint = _surface_mesh_fingerprint(
        verts, faces, active, material, geometry
    )
    return geometry, state, fingerprint


def test_grid_contract_refuses_nominal_twenty_nm_for_published_cell():
    twenty = AUDIT.grid_contract(0.02)
    ten = AUDIT.grid_contract(0.01)
    five = AUDIT.grid_contract(0.005)

    assert not twenty["compatible"]
    assert np.allclose(twenty["realized_extent_um"], (0.12, 0.04, 2.8))
    assert "does not divide" in twenty["reason"]
    assert "two intervals" in twenty["reason"]
    assert ten["compatible"]
    assert five["compatible"]


def test_plan_is_calibration_only_and_blocks_changed_domain(tmp_path):
    plan = AUDIT.build_plan(tmp_path / "missing-checkpoint.npz")

    assert "no held-out outcomes read" in plan["scientific_status"]
    by_dx = {item["dx_nm"]: item for item in plan["levels"]}
    assert by_dx[20.0]["initial_case"] == "blocked"
    assert by_dx[20.0]["late_case"] == "blocked"
    assert by_dx[10.0]["initial_case"] == "eligible"
    assert plan["pairing_contract"]["schedule_owner"] == "5 nm"


def test_late_source_provenance_refuses_a_heldout_boundary(tmp_path):
    checkpoint = tmp_path / "checkpoint.npz"
    checkpoint.write_bytes(b"provenance-name-only")
    configuration = {
        "boundary_case": "oxygen_ratio",
        "oxygen_to_fluorocarbon_ratio": 1.5,
        "low_frequency_power_kw": 6.0,
        "effective_mask_crosslinked_growth_fraction": (
            AUDIT.CALIBRATION["effective_mask_crosslinked_growth_fraction"]
        ),
        "oxide_etch_yield_scale": AUDIT.CALIBRATION["oxide_etch_yield_scale"],
        "dx_um": 0.005,
        "seed": 241,
    }
    (tmp_path / "audit.json").write_text(
        json.dumps({
            "configuration": configuration,
            "history": [{"physical_time_s": 56.0}],
        }),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="not the declared base-calibration source"):
        AUDIT._validate_source_provenance(
            checkpoint, {"dx": 0.005, "physical_time_s": 56.0}
        )


def test_aligned_checkpoint_restriction_preserves_topology_and_conservative_state():
    geometry, state, fingerprint = _manufactured_checkpoint_state()
    restricted, remapped, next_fingerprint, provenance = (
        AUDIT.restrict_checkpoint_state(geometry, state, fingerprint, 0.01)
    )

    assert restricted.dx == 0.01
    assert restricted.phi.shape == (5, 3, 11)
    assert provenance["restriction_factor"] == 2
    assert provenance["method"].startswith("aligned_nodal_material_levelset_restriction")
    assert next_fingerprint != fingerprint
    assert set(remapped.fields) == set(state.fields)
    residual = max(
        item["max_relative_conservation_residual"]
        for item in provenance["surface_state_remap"]["materials"].values()
    )
    assert residual < 1e-12


def test_identity_restriction_is_exact():
    geometry, state, fingerprint = _manufactured_checkpoint_state()
    restricted, remapped, next_fingerprint, provenance = (
        AUDIT.restrict_checkpoint_state(geometry, state, fingerprint, geometry.dx)
    )

    assert restricted is geometry
    assert remapped is state
    assert next_fingerprint == fingerprint
    assert provenance["method"] == "identity"


def test_exact_zero_depth_widths_use_declared_initial_csg_not_experiment():
    completed, provenance = AUDIT._complete_analytic_initial_metrics({
        "mask_opening_nm": 90.0,
        "top_feature_width_nm": np.nan,
        "maximum_feature_width_nm": np.nan,
    })

    assert completed["top_feature_width_nm"] == 90.0
    assert completed["maximum_feature_width_nm"] == 90.0
    assert set(provenance["completed_fields"]) == {
        "top_feature_width_nm", "maximum_feature_width_nm"
    }
    assert "CSG" in provenance["method"]


def test_worker_command_replays_fine_schedule(tmp_path):
    args = AUDIT.parse_args([
        "--output", str(tmp_path / "out"),
        "--source-checkpoint", str(tmp_path / "source.npz"),
    ])
    schedule = tmp_path / "schedule.json"
    command = AUDIT._worker_command(
        args, "late", 10.0, tmp_path / "late_10nm", schedule=schedule
    )

    assert command[command.index("--worker-schedule") + 1] == str(schedule)
    assert command[command.index("--worker-dx-nm") + 1] == "10"
    assert command[command.index("--topology-policy") + 1] == "refuse"


def test_worker_environment_unifies_transport_and_level_set_device():
    environment = AUDIT._worker_environment("cuda:0")

    assert environment["PETCH_DEVICE"] == "cuda:0"
    assert environment["OMP_NUM_THREADS"] == "1"
