import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = runpy.run_path(
    str(ROOT / "scripts" / "export_lan_jeon_bolsig_reference.py"))


def test_exported_collision_topology_and_units_are_explicit():
    payload = MODULE["collision_text"]()

    assert payload.count("\nELASTIC\n") == 1
    assert payload.count("\nATTACHMENT\n") == 1
    assert payload.count("\nIONIZATION\n") == 1
    assert payload.count("\nEXCITATION\n") == 7
    assert "0\t8.3900000000000002e-19" in payload
    assert "362\t1.1e-20" in payload
    assert "doi:10.3938/jkps.64.1320" in payload


def test_reference_run_requests_bulk_transport_at_all_source_fields():
    payload = MODULE["instruction_text"]()

    assert "4 / Density-gradient expansion" in payload
    assert "800 / Number of grid points" in payload
    assert payload.count("/ File") == 2
    fields = MODULE["_figure7_fields"]()
    assert len(fields) == 18
    assert fields[0] > 100.0
    assert fields[-1] < 1200.0
    for field in fields:
        assert f"{field:.17g}\n" in payload
