import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[1] / "scripts"
    / "block_levelset_manufactured_audit_3d.py")
SPEC = importlib.util.spec_from_file_location(
    "block_levelset_manufactured_audit_3d", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_candidate_shapes_respect_thin_periodic_axis_and_partial_x():
    shapes = AUDIT.candidate_block_shapes((27, 5, 561))

    assert len(shapes) == 4
    assert all(shape[1] == 4 for shape in shapes)
    assert (26, 4, 32) in shapes
    assert all(all(value > 0 for value in shape) for shape in shapes)


def test_decision_refuses_optimistic_deep_storage_below_gates():
    records = []
    for value in (1.2, 1.8):
        records.append({
            "dx_nm": 5.0,
            "etched_depth_um": 0.9,
            "band_width_cells": 8,
            "integrity_pass": True,
            "core_memory_reduction": value,
            "indexed_halo_only_memory_reduction": value - 0.1,
            "unique_node_work_upper_bound": value - 0.2,
        })

    decision = AUDIT.decide(records)

    assert decision["status"] == "fixed_dx_sparse_no_go_for_krueger"
    assert not decision["pass"]
    assert "surface/transport" in decision["scientific_action"]
