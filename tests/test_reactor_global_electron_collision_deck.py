import hashlib

import pytest

from petch.reactor_global.electron_collision_deck import (
    load_bolsig_lxcat_file,
    parse_bolsig_lxcat_bytes,
)


DECK = b"""Example preamble that must not be interpreted as data.

COMMENT
Cl2
References: manufactured parser fixture only
------------------------------------------------------------

ELASTIC
Cl2
7.680e-6 / electron-to-target mass ratio
COMMENT: manufactured constant momentum transfer
------------------------------------------------------------
0.0 2.0e-20
1.0 2.0e-20
100.0 2.0e-20
------------------------------------------------------------

EXCITATION
Cl2 -> Cl2(v=1)
0.069 / threshold energy
COMMENT: manufactured vibrational excitation
------------------------------------------------------------
0.0 0.0
0.069 0.0
1.0 1.0e-21
------------------------------------------------------------

IONIZATION
Cl2 -> Cl2+
11.49 / threshold energy
------------------------------------------------------------
0.0 0.0
11.49 0.0
100.0 4.0e-20
------------------------------------------------------------

ATTACHMENT
Cl2 -> Cl- + Cl
COMMENT: attachment has no third parameter line
------------------------------------------------------------
0.0 0.0
0.1 1.0e-20
10.0 1.0e-23
------------------------------------------------------------
"""


def _parse(payload=DECK, **changes):
    metadata = {
        "source_database": "manufactured-test-deck",
        "retrieved_at": "2026-08-08",
        "source_reference": "tests only; not physical evidence",
        "target": "Cl2",
    }
    metadata.update(changes)
    return parse_bolsig_lxcat_bytes(payload, **metadata)


def test_parser_hash_locks_user_supplied_deck_and_preserves_process_metadata():
    deck = _parse()
    assert deck.payload_sha256 == hashlib.sha256(DECK).hexdigest()
    assert deck.targets == ("Cl2",)
    assert [item.kind for item in deck.processes] == [
        "ELASTIC", "EXCITATION", "IONIZATION", "ATTACHMENT",
    ]

    elastic, excitation, ionization, attachment = deck.processes
    assert elastic.mass_ratio == pytest.approx(7.680e-6)
    assert elastic.energy_loss_eV is None
    assert excitation.product == "Cl2(v=1)"
    assert excitation.energy_loss_eV == pytest.approx(0.069)
    assert ionization.energy_loss_eV == pytest.approx(11.49)
    assert elastic.electron_number_change == 0
    assert excitation.electron_number_change == 0
    assert ionization.electron_number_change == 1
    assert attachment.electron_number_change == -1
    assert attachment.energy_loss_eV == 0.0
    assert attachment.product == "Cl- + Cl"
    assert attachment.comments == (
        "COMMENT: attachment has no third parameter line",
    )


def test_explicit_multiple_ionization_multiplicity_is_never_inferred_from_label():
    base = _parse().processes[2]
    from dataclasses import replace

    double = replace(
        base,
        product="Cl2++",
        electron_number_change=2,
    )
    assert double.electron_number_change == 2
    with pytest.raises(ValueError, match="electron-number change"):
        replace(
            _parse().processes[1],
            electron_number_change=1,
        )


def test_structural_readiness_is_not_promoted_to_reactor_or_depth_evidence():
    readiness = _parse().structural_kinetic_readiness("Cl2")
    assert not readiness["structurally_ready_for_kinetic_input"]
    assert readiness["issues"] == (
        "multiterm_requires_explicit_elastic_angular_closure",
    )
    assert readiness["process_count"] == 4
    assert readiness["momentum_process_count"] == 1
    assert readiness["ionization_process_count"] == 1
    assert readiness["attachment_process_count"] == 1
    assert not readiness["supports_swarm_validation"]
    assert not readiness["supports_reactor_state_prediction"]
    assert not readiness["supports_wafer_flux"]
    assert not readiness["supports_feature_depth"]

    replay = _parse().structural_kinetic_readiness(
        "Cl2", elastic_angular_closure="isotropic_source_reproduction")
    assert replay["structurally_ready_for_kinetic_input"]
    assert replay["angular_evidence_class"] == (
        "source_reproduction_assumption")
    assert not replay["contains_differential_elastic_cross_sections"]

    with pytest.raises(ValueError, match="unsupported elastic angular"):
        _parse().structural_kinetic_readiness(
            "Cl2", elastic_angular_closure="silent_isotropic_default")


def test_parser_refuses_hash_mismatch_and_nonmonotone_cross_section_grid():
    with pytest.raises(RuntimeError, match="hash mismatch"):
        _parse(expected_sha256="0" * 64)

    damaged = DECK.replace(
        b"0.0 2.0e-20\n1.0 2.0e-20\n100.0 2.0e-20",
        b"0.0 2.0e-20\n1.0 2.0e-20\n0.5 2.0e-20",
    )
    with pytest.raises(ValueError, match="invalid electron collision process"):
        _parse(damaged)


def test_file_loader_never_packages_or_redistributes_source_bytes(tmp_path):
    path = tmp_path / "user-supplied-lxcat.txt"
    path.write_bytes(DECK)
    deck = load_bolsig_lxcat_file(
        path,
        source_database="manufactured-test-deck",
        retrieved_at="2026-08-08",
        source_reference="tests only; not physical evidence",
        target="Cl2",
    )
    assert not deck.packaged_or_redistributed
    assert not hasattr(deck, "source_payload")
    assert deck.payload_sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
