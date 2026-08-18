"""Explicit deconvolution of LXCat/BOLSIG effective momentum sets.

An ``EFFECTIVE`` row is the elastic momentum-transfer cross section plus the
sum of the set's inelastic cross sections.  Petch's collision operator already
adds every inelastic row to momentum relaxation, so consuming an effective row
directly would count those processes twice.  This module performs the inverse
operation on the complete union of source knots and records a derivation hash.
It packages no source cross-section bytes.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import json
import math

import numpy as np

from .electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)


_MOMENTUM_KINDS = frozenset({"EFFECTIVE", "ELASTIC", "MOMENTUM"})


@dataclass(frozen=True)
class EffectiveMomentumDeconvolution:
    """One source set and its solver-ready elastic-momentum derivation."""

    source_deck: ElectronCollisionDeck
    derived_deck: ElectronCollisionDeck
    target: str
    effective_process_index: int
    inelastic_process_count: int
    minimum_elastic_cross_section_m2: float
    maximum_recomposition_relative_residual: float
    algorithm: str = "source_linear_union_knots_subtract_all_inelastic_v1"
    supports_collision_input_replay: bool = True
    supports_swarm_validation: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            self.derived_deck is self.source_deck
            or self.target not in self.source_deck.targets
            or self.derived_deck.targets != (self.target,)
            or self.effective_process_index < 0
            or self.inelastic_process_count < 1
            or not math.isfinite(self.minimum_elastic_cross_section_m2)
            or self.minimum_elastic_cross_section_m2 <= 0.0
            or not math.isfinite(
                self.maximum_recomposition_relative_residual)
            or self.maximum_recomposition_relative_residual < 0.0
            or self.algorithm
            != "source_linear_union_knots_subtract_all_inelastic_v1"
            or not self.supports_collision_input_replay
            or self.supports_swarm_validation
            or self.supports_reactor_state_prediction
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid effective-momentum deconvolution")


def _cross_section_on_grid(
    process: ElectronCollisionProcess,
    energy_eV: np.ndarray,
) -> np.ndarray:
    source_energy = np.asarray(process.electron_energy_eV, dtype=float)
    source_sigma = np.asarray(process.cross_section_m2, dtype=float)
    if energy_eV[-1] > source_energy[-1] and source_sigma[-1] != 0.0:
        raise ValueError(
            "nonzero inelastic cross section ends below effective support"
        )
    return np.interp(
        energy_eV,
        source_energy,
        source_sigma,
        left=0.0,
        right=0.0,
    )


def deconvolve_effective_momentum(
    deck: ElectronCollisionDeck,
    target: str,
    *,
    retrieved_at: str,
    source_reference: str,
) -> EffectiveMomentumDeconvolution:
    """Replace one target's ``EFFECTIVE`` row by elastic momentum transfer.

    The source deck must contain only the selected target, exactly one
    effective row, no independently supplied elastic row, and at least one
    inelastic row.  Subtraction is evaluated on the union of every knot so a
    negative segment cannot hide between the effective row's native knots.
    """

    if not isinstance(deck, ElectronCollisionDeck):
        raise TypeError("an ElectronCollisionDeck is required")
    target = str(target).strip()
    date = str(retrieved_at).strip()
    reference = str(source_reference).strip()
    if not target or not date or not reference:
        raise ValueError("target and derivation metadata must be non-empty")
    if deck.targets != (target,):
        raise ValueError("source deck must contain exactly the selected target")
    momentum = [
        (index, process)
        for index, process in enumerate(deck.processes)
        if process.kind in _MOMENTUM_KINDS
    ]
    if len(momentum) != 1 or momentum[0][1].kind != "EFFECTIVE":
        raise ValueError("source target requires exactly one EFFECTIVE row")
    effective_index, effective = momentum[0]
    inelastic = tuple(
        process for process in deck.processes
        if process.kind not in _MOMENTUM_KINDS
    )
    if not inelastic:
        raise ValueError("effective deconvolution requires inelastic rows")

    effective_energy = np.asarray(effective.electron_energy_eV, dtype=float)
    union = {float(value) for value in effective_energy}
    lower = float(effective_energy[0])
    upper = float(effective_energy[-1])
    for process in inelastic:
        union.update(
            float(value) for value in process.electron_energy_eV
            if lower <= value <= upper
        )
    energy = np.asarray(sorted(union), dtype=float)
    effective_sigma = np.interp(
        energy,
        effective_energy,
        np.asarray(effective.cross_section_m2, dtype=float),
    )
    inelastic_sum = np.zeros_like(energy)
    for process in inelastic:
        inelastic_sum += _cross_section_on_grid(process, energy)
    elastic_sigma = effective_sigma - inelastic_sum
    scale = max(float(np.max(effective_sigma)), np.finfo(float).tiny)
    roundoff_tolerance = 128.0 * np.finfo(float).eps * scale
    if np.any(elastic_sigma < -roundoff_tolerance):
        worst = int(np.argmin(elastic_sigma))
        raise ValueError(
            "effective set implies negative elastic momentum transfer at "
            f"{energy[worst]:.12g} eV: {elastic_sigma[worst]:.12g} m2"
        )
    elastic_sigma = np.maximum(elastic_sigma, 0.0)
    if np.any(elastic_sigma <= 0.0):
        raise ValueError("effective set implies nonpositive elastic momentum")
    recomposed = elastic_sigma + inelastic_sum
    relative = np.abs(recomposed - effective_sigma) / np.maximum(
        effective_sigma, np.finfo(float).tiny)

    derivation = {
        "schema": "petch.effective_momentum_deconvolution.v1",
        "source_payload_sha256": deck.payload_sha256,
        "target": target,
        "effective_process_index": effective_index,
        "inelastic_process_count": len(inelastic),
        "algorithm": "source_linear_union_knots_subtract_all_inelastic_v1",
        "source_cross_sections_modified": True,
        "source_bytes_redistributed": False,
    }
    digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    derived_momentum = replace(
        effective,
        kind="MOMENTUM",
        electron_energy_eV=tuple(float(value) for value in energy),
        cross_section_m2=tuple(float(value) for value in elastic_sigma),
        comments=effective.comments + (
            "Derived elastic momentum = source EFFECTIVE - sum(source "
            "inelastic), evaluated on the complete union of knots.",
        ),
    )
    derived_processes = list(deck.processes)
    derived_processes[effective_index] = derived_momentum
    derived = ElectronCollisionDeck(
        processes=tuple(derived_processes),
        payload_sha256=digest,
        source_database=f"derived from {deck.source_database}",
        retrieved_at=date,
        source_reference=(
            reference + "; derivation="
            + json.dumps(derivation, sort_keys=True)
        ),
    )
    return EffectiveMomentumDeconvolution(
        source_deck=deck,
        derived_deck=derived,
        target=target,
        effective_process_index=effective_index,
        inelastic_process_count=len(inelastic),
        minimum_elastic_cross_section_m2=float(np.min(elastic_sigma)),
        maximum_recomposition_relative_residual=float(np.max(relative)),
    )
