from petch.stl_import import StlMesh, write_stl
from scripts.audit_stl_geometry import _planar_extrusion_payload, build

import numpy as np


def _box_with_point_facet():
    vertices = np.asarray([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
        [0, 0, 2], [1, 0, 2], [0, 1, 2], [1, 1, 2],
        [1, 1, 2],
    ], dtype=float)
    faces = np.asarray([
        [0, 2, 3], [0, 3, 1], [4, 7, 6], [4, 5, 7],
        [0, 1, 5], [0, 5, 4], [2, 6, 7], [2, 7, 3],
        [0, 4, 6], [0, 6, 2], [1, 3, 7], [1, 7, 5],
        [8, 8, 8],
    ], dtype=int)
    return StlMesh(vertices, faces)


def test_partner_stl_audit_repairs_only_point_facets_and_refuses_units(tmp_path):
    source = write_stl(tmp_path / "input.stl", _box_with_point_facet())
    repaired, audit = build(source)

    assert len(repaired.faces) == 12
    assert audit["repair"]["removed_face_count"] == 1
    assert audit["clean_diagnostics"]["failure_reason"] is None
    assert audit["connected_component_count"] == 1
    # An axis-aligned manufactured box is extruded along all three axes; the
    # audit must report that ambiguity instead of choosing a convenient axis.
    assert audit["planar_extrusion"]["unique_candidate"] is False
    assert audit["planar_extrusion"]["candidate_axis"] is None
    assert audit["planar_extrusion"]["candidate_axes"] == [0, 1, 2]
    assert audit["planar_extrusion"]["projected_fill_fraction"] is None
    assert audit["units"]["physical_scale_identified"] is False
    assert audit["simulation_readiness"]["geometry_topology_ready"] is True
    assert audit["simulation_readiness"]["physical_scale_ready"] is False


def test_unique_planar_extrusion_reports_scale_free_projected_metrics():
    payload = _planar_extrusion_payload(
        extent=np.asarray([2.0, 5.0, 4.0]),
        components=[
            {"signed_volume_file_units_cubed": 12.0},
            {"signed_volume_file_units_cubed": 4.0},
        ],
        candidates=[0],
    )

    assert payload["candidate_axis"] == 0
    assert payload["extrusion_thickness_file_units"] == 2.0
    assert payload["projected_solid_area_file_units_squared"] == 8.0
    assert payload["projected_bounding_box_area_file_units_squared"] == 20.0
    assert payload["projected_fill_fraction"] == 0.4
    assert payload["component_projected_area_file_units_squared"] == [6.0, 2.0]
