from pathlib import Path

import pytest

from scripts.audit_guo_karahashi_transfer import build_audit


ROOT = Path(__file__).parents[1]


@pytest.fixture(scope="module")
def audit():
    return build_audit(ROOT)


def test_guo_karahashi_board_is_complete_and_has_no_fit(audit):
    assert len(audit["point_board"]) == 21
    assert {
        row["species"] for row in audit["point_board"]
    } == {"F+", "CF+", "CF2+", "CF3+"}
    assert audit["calibration_firewall"]["fitted_parameters"] == []
    assert not audit["calibration_firewall"][
        "karahashi_yields_used_by_surface_solver"]


def test_guo_karahashi_receipts_are_checksum_bound(audit):
    receipts = audit["source_receipts"]
    assert receipts["karahashi_data_sha256"] == (
        receipts["loader_pinned_karahashi_data_sha256"])
    assert receipts["karahashi_data_sha256"] == (
        "5d3b58a6d23e3fa77b5e1484407dbe0f19f23ec89687f184cd11ceaebb017c26")
    assert receipts["guo_reaction_deck_sha256"] == (
        "2f93b31c1862095d133b765cddc1957ec27f0e0d05b7b2dfbafc0908f54b9600")


def test_guo_karahashi_atomic_and_numerical_gates_pass(audit):
    gates = audit["numerical_and_atomic_gates"]
    assert gates["all_steady_states_pass"]
    assert gates["all_atom_ledgers_pass"]


def test_guo_karahashi_independent_transfer_exposes_species_defects(audit):
    board = {
        (row["species"], row["energy_eV"]): row
        for row in audit["point_board"]
    }
    assert board[("CF2+", 1000.0)]["predicted_yield_sio2_per_ion"] == (
        pytest.approx(1.4183200338533348, rel=2e-9))
    assert board[("CF3+", 1000.0)]["predicted_yield_sio2_per_ion"] == (
        pytest.approx(1.7608577618293395, rel=2e-9))
    assert board[("CF3+", 250.0)]["predicted_yield_sio2_per_ion"] < 0.0
    assert board[("CF3+", 500.0)]["predicted_yield_sio2_per_ion"] < 0.0
    assert not audit["verdict"]["species_resolved_transfer_validated"]


def test_committed_guo_karahashi_audit_is_current(audit):
    path = (
        ROOT / "results" / "curated" / "guo_karahashi_transfer"
        / "audit.json")
    assert path.exists()
    import json
    assert json.loads(path.read_text(encoding="utf-8")) == audit
