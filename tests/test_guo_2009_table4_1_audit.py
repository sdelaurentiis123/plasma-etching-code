import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = (
    ROOT / "results" / "curated" / "guo_2009_table4_1" / "audit.json")


def _load_auditor():
    path = ROOT / "scripts" / "audit_guo_2009_table4_1.py"
    spec = importlib.util.spec_from_file_location("guo_table4_1_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_guo_table_audit_matches_rebuild_and_closes_atoms():
    frozen = json.loads(AUDIT_PATH.read_text())
    rebuilt = _load_auditor().build_audit(ROOT)

    assert frozen == rebuilt
    assert rebuilt["status"] == "passed"
    assert rebuilt["row_count"] == 20
    assert rebuilt["sequential_reaction_ids"]
    assert rebuilt["table_csv_hash_matches_manifest"]
    assert rebuilt["thresholds_match_visual_transcription"]
    assert rebuilt["all_atom_counted_reactions_close"]
    assert len(rebuilt["atom_counted_reactions"]) == 15


def test_guo_coefficients_keep_rate_semantics_and_evidence_ceiling():
    audit = json.loads(AUDIT_PATH.read_text())

    assert "not bounded sticking probabilities" in audit[
        "coefficient_semantics"]
    assert audit["coefficients_greater_than_one"] == {
        "S_F": 20.0,
        "S_N_on_C": 3.5,
        "S_N_on_O": 1.8,
        "beta_C-V": 10000.0,
        "beta_CF2": 2.0,
        "beta_SiF2": 6.75,
        "beta_d_V_by_I": 1.66,
    }
    assert audit["evidence_ceiling"] == "L1_yield_regressed"
    assert "not independent validation" in audit["claim"]
    assert "Krueger absolute-depth fit" in audit["claim"]


def test_visual_manifest_pins_source_and_preserves_source_contradictions():
    audit = json.loads(AUDIT_PATH.read_text())

    assert audit["source_pdf_sha256"] == (
        "f5c78c0089fe4104019435c6fd34e9b8f284358dda1df0101ecec54c592750d2")
    pages = audit["visual_audit"]["pages"]
    assert [page["pdf_page"] for page in pages] == [
        100, 101, 102, 103, 35, 45,
    ]
    assert all(page["render_mode"] == "RGB" for page in pages)
    assert all(len(page["render_sha256"]) == 64 for page in pages)
    assert len(audit["source_inconsistencies_preserved"]) == 3
    assert len(audit["visual_audit"]["formalism_pages"]) == 3
