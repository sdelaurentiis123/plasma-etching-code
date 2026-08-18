"""Deterministic composition of independently provenance-locked decks."""
from __future__ import annotations

from hashlib import sha256
import json
from typing import Iterable

from .electron_collision_deck import ElectronCollisionDeck


def compose_electron_collision_decks(
    decks: Iterable[ElectronCollisionDeck],
    *,
    retrieved_at: str,
    mixture_name: str,
) -> ElectronCollisionDeck:
    """Combine disjoint target decks without changing any cross section.

    Mole fractions remain part of the Boltzmann condition, not the collision
    data.  This function only provides one immutable process collection and a
    derivation hash that is independent of mixture composition.
    """

    components = tuple(decks)
    if not components or any(
        not isinstance(item, ElectronCollisionDeck) for item in components
    ):
        raise ValueError("at least one electron collision deck is required")
    name = str(mixture_name).strip()
    date = str(retrieved_at).strip()
    if not name or not date:
        raise ValueError("mixture name and retrieval date must be non-empty")
    owner: dict[str, int] = {}
    for index, deck in enumerate(components):
        for target in deck.targets:
            if target in owner:
                raise ValueError(
                    f"electron target {target!r} appears in multiple decks"
                )
            owner[target] = index
    derivation = {
        "schema": "petch.electron_collision_mixture.v1",
        "mixture_name": name,
        "components": [
            {
                "payload_sha256": deck.payload_sha256,
                "targets": list(deck.targets),
            }
            for deck in components
        ],
        "cross_sections_modified": False,
        "mole_fractions_embedded": False,
    }
    digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    return ElectronCollisionDeck(
        processes=tuple(
            process for deck in components for process in deck.processes
        ),
        payload_sha256=digest,
        source_database=f"composed electron collision deck: {name}",
        retrieved_at=date,
        source_reference=(
            "No source curve modified; component references: "
            + " | ".join(deck.source_reference for deck in components)
            + "; derivation="
            + json.dumps(derivation, sort_keys=True)
        ),
    )
