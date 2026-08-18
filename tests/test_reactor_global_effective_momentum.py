from dataclasses import replace

import numpy as np
import pytest

from petch.reactor_global.effective_momentum import (
    deconvolve_effective_momentum,
)
from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)


def _deck(effective=(5.0e-20, 6.0e-20, 5.0e-20)):
    processes = (
        ElectronCollisionProcess(
            kind="EFFECTIVE",
            target="H2",
            product=None,
            mass_ratio=2.72e-4,
            electron_energy_eV=(0.0, 2.0, 10.0),
            cross_section_m2=effective,
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            target="H2",
            product="H2(v=1)",
            energy_loss_eV=1.0,
            electron_energy_eV=(1.0, 4.0, 10.0),
            cross_section_m2=(0.0, 2.0e-20, 1.0e-20),
        ),
    )
    return ElectronCollisionDeck(
        processes=processes,
        payload_sha256="a" * 64,
        source_database="manufactured effective set",
        retrieved_at="2026-08-18",
        source_reference="tests only",
    )


def _derive(deck=None):
    return deconvolve_effective_momentum(
        _deck() if deck is None else deck,
        "H2",
        retrieved_at="2026-08-18",
        source_reference="manufactured deconvolution gate",
    )


def test_deconvolution_uses_union_of_knots_and_recomposes_effective_row():
    result = _derive()
    momentum, excitation = result.derived_deck.processes
    assert momentum.kind == "MOMENTUM"
    assert momentum.electron_energy_eV == (0.0, 1.0, 2.0, 4.0, 10.0)
    energy = np.asarray(momentum.electron_energy_eV)
    effective = np.interp(energy, (0.0, 2.0, 10.0), (5e-20, 6e-20, 5e-20))
    inelastic = np.interp(
        energy, (1.0, 4.0, 10.0), (0.0, 2e-20, 1e-20),
        left=0.0, right=0.0,
    )
    assert np.asarray(momentum.cross_section_m2) + inelastic == pytest.approx(
        effective)
    assert result.maximum_recomposition_relative_residual < 1.0e-15
    assert result.inelastic_process_count == 1
    assert not result.supports_swarm_validation
    assert not result.supports_feature_depth
    assert result.derived_deck.payload_sha256 != result.source_deck.payload_sha256


def test_deconvolution_fails_closed_for_negative_or_truncated_inelastic_set():
    negative = _deck(effective=(1.0e-20, 1.0e-20, 1.0e-20))
    with pytest.raises(ValueError, match="negative elastic momentum"):
        _derive(negative)

    source = _deck()
    truncated = replace(
        source.processes[1],
        electron_energy_eV=(1.0, 4.0),
        cross_section_m2=(0.0, 2.0e-20),
    )
    damaged = replace(source, processes=(source.processes[0], truncated))
    with pytest.raises(ValueError, match="ends below effective support"):
        _derive(damaged)
