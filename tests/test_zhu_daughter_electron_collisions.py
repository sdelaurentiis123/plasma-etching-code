from dataclasses import replace

import numpy as np
import pytest

from petch.reactor_global.argon import ELECTRON_MASS_AMU
from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)
from petch.reactor_global.zhu_daughter_electron_collisions import (
    HF_DISSOCIATION_THRESHOLD_EV,
    HF_IONIZATION_THRESHOLD_EV,
    HF_MASS_AMU,
    deconvolve_siglo_f2_effective_momentum,
    derive_huang_2020_partial_hf_replay,
)


def _deck(processes, *, target):
    return ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256=("a" if target == "HCl" else "b") * 64,
        source_database="manufactured test deck",
        retrieved_at="2026-08-18",
        source_reference="tests only; not physical evidence",
    )


def _hcl_deck():
    common = {
        "target": "HCl",
        "electron_energy_eV": (0.0, 20.0, 120.0),
        "cross_section_m2": (2.0e-20, 2.0e-20, 2.0e-20),
    }
    return _deck((
        ElectronCollisionProcess(
            kind="ELASTIC", product=None, mass_ratio=1.5e-5, **common
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            product="HCl(v1)",
            energy_loss_eV=0.349,
            **common,
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            product="HCl(v2)",
            energy_loss_eV=0.691,
            **common,
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            product="HCl*(5.29eV)",
            energy_loss_eV=5.29,
            **common,
        ),
        ElectronCollisionProcess(
            kind="IONIZATION",
            product="HCl+",
            energy_loss_eV=12.5,
            **common,
        ),
        ElectronCollisionProcess(
            kind="ATTACHMENT",
            product="HCl-",
            energy_loss_eV=0.0,
            **common,
        ),
    ), target="HCl")


def _f2_deck(*, effective_scale=1.0):
    return _deck((
        ElectronCollisionProcess(
            kind="EFFECTIVE",
            target="F2",
            product=None,
            electron_energy_eV=(0.0, 10.0, 120.0),
            cross_section_m2=tuple(
                effective_scale * value
                for value in (8.0e-20, 8.0e-20, 8.0e-20)
            ),
            mass_ratio=1.43e-5,
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            target="F2",
            product="F2(v)",
            electron_energy_eV=(0.1, 5.0),
            cross_section_m2=(1.0e-20, 0.0),
            energy_loss_eV=0.1,
        ),
        ElectronCollisionProcess(
            kind="IONIZATION",
            target="F2",
            product="F2+",
            electron_energy_eV=(10.0, 120.0),
            cross_section_m2=(0.0, 2.0e-20),
            energy_loss_eV=10.0,
        ),
    ), target="F2")


def test_huang_partial_hf_replay_transfers_only_declared_hcl_channels():
    replay = derive_huang_2020_partial_hf_replay(_hcl_deck())
    assert [item.kind for item in replay.derived_deck.processes] == [
        "ELASTIC", "EXCITATION", "IONIZATION",
    ]
    elastic, dissociation, ionization = replay.derived_deck.processes
    assert elastic.mass_ratio == pytest.approx(ELECTRON_MASS_AMU / HF_MASS_AMU)
    assert dissociation.energy_loss_eV == HF_DISSOCIATION_THRESHOLD_EV
    assert dissociation.electron_energy_eV[0] == pytest.approx(0.58)
    assert dissociation.product == "H + F"
    assert ionization.energy_loss_eV == HF_IONIZATION_THRESHOLD_EV
    assert ionization.electron_energy_eV[0] == pytest.approx(3.507)
    assert ionization.product == "HF+"
    assert not replay.supports_complete_hf_eedf
    assert not replay.supports_unique_reactor_state
    assert not replay.supports_feature_depth
    assert all("vibrational" not in (item.product or "") for item in (
        dissociation, ionization
    ))


def test_f2_effective_deconvolution_reconstructs_total_momentum_exactly():
    source = _f2_deck()
    replay = deconvolve_siglo_f2_effective_momentum(source)
    assert replay.derived_deck.processes[0].kind == "ELASTIC"
    assert all(
        process.kind != "EFFECTIVE"
        for process in replay.derived_deck.processes
    )
    energy = np.asarray(replay.derived_deck.processes[0].electron_energy_eV)
    reconstructed = np.zeros_like(energy)
    for process in replay.derived_deck.processes:
        reconstructed += np.interp(
            energy,
            process.electron_energy_eV,
            process.cross_section_m2,
            left=0.0,
            right=process.cross_section_m2[-1],
        )
    effective = source.processes[0]
    expected = np.interp(
        energy, effective.electron_energy_eV, effective.cross_section_m2
    )
    assert reconstructed == pytest.approx(expected, rel=2.0e-15)
    assert replay.minimum_elastic_cross_section_m2 > 0.0
    assert replay.derived_deck.processes[1].electron_energy_eV[-1] == 120.0


def test_f2_deconvolution_fails_closed_on_inconsistent_effective_set():
    with pytest.raises(ValueError, match="became nonpositive"):
        deconvolve_siglo_f2_effective_momentum(
            _f2_deck(effective_scale=0.05)
        )

    source = _f2_deck()
    nonzero_short_tail = replace(
        source.processes[1], cross_section_m2=(1.0e-20, 1.0e-21)
    )
    damaged = replace(
        source,
        processes=(source.processes[0], nonzero_short_tail, source.processes[2]),
    )
    with pytest.raises(ValueError, match="nonzero tail"):
        deconvolve_siglo_f2_effective_momentum(damaged)
