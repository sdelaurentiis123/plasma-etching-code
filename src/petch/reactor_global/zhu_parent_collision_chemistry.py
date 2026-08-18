"""Conserved parent-feed electron chemistry for the Zhu NPG80 condition.

This module maps every non-momentum row in the independently audited
CHF3/SF6/O2 collision decks to explicit heavy products.  It is the bridge
between an EEPF and particle balances: product labels alone are not chemistry.

The mapping is deliberately honest about two remaining branch closures.  The
NIST SF6 neutral-dissociation total is mapped to its literature-dominant
SF5 + F channel until its partial curves are installed, and the Song O2 total
positive ionization is mapped to O2+ until its O+/O2+ split is installed.
Neither closure can be selected from the withheld TiO2 SEM.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .argon import ELECTRON_MASS_AMU
from .chf3_electron_collisions import (
    NISTEvaluatedCHF3Replay,
    derive_nist_evaluated_chf3_replay,
    load_kushner_zhang_2000_chf3_replay,
)
from .electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionHeavyMapping,
)
from .electron_collision_deck import ElectronCollisionDeck
from .electron_collision_mixture import compose_electron_collision_decks
from .network import Species
from .o2_electron_collisions import Song2026O2Replay, load_song_2026_o2_replay
from .sf6_electron_collisions import (
    NISTProductResolvedSF6Replay,
    derive_nist_product_resolved_sf6_replay,
)


_ATOMIC_MASS_AMU = {
    "C": 12.011,
    "H": 1.00784,
    "F": 18.998403163,
    "O": 15.9994,
    "S": 32.065,
}


@dataclass(frozen=True)
class ZhuParentCollisionChemistry:
    """One immutable collision deck plus its fully covered product mapping."""

    chf3_replay: NISTEvaluatedCHF3Replay
    sf6_replay: NISTProductResolvedSF6Replay
    o2_replay: Song2026O2Replay
    mixed_deck: ElectronCollisionDeck
    species: tuple[Species, ...]
    collision_chemistry: ElectronCollisionChemistry
    sf6_neutral_dissociation_closure: str = "dominant_SF5_plus_F"
    o2_positive_ionization_closure: str = "all_O2plus"
    supports_parent_collision_sources: bool = True
    supports_complete_daughter_eedf: bool = False
    supports_unique_reactor_state: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        nonmomentum = sum(
            process.kind not in {"MOMENTUM", "ELASTIC", "EFFECTIVE"}
            for process in self.mixed_deck.processes
        )
        if (
            self.collision_chemistry.collision_deck is not self.mixed_deck
            or self.collision_chemistry.species != self.species
            or len(self.collision_chemistry.mappings) != nonmomentum
            or self.sf6_neutral_dissociation_closure
            != "dominant_SF5_plus_F"
            or self.o2_positive_ionization_closure != "all_O2plus"
            or not self.supports_parent_collision_sources
            or self.supports_complete_daughter_eedf
            or self.supports_unique_reactor_state
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid Zhu parent collision chemistry")


def _mass(composition: dict[str, int]) -> float:
    return sum(
        _ATOMIC_MASS_AMU[element] * count
        for element, count in composition.items()
    )


def zhu_parent_collision_species() -> tuple[Species, ...]:
    """Return the exact union needed by all parent collision products."""

    neutral = (
        ("CHF3", {"C": 1, "H": 1, "F": 3}),
        ("CF3", {"C": 1, "F": 3}),
        ("CHF2", {"C": 1, "H": 1, "F": 2}),
        ("CF2", {"C": 1, "F": 2}),
        ("CHF", {"C": 1, "H": 1, "F": 1}),
        ("CF", {"C": 1, "F": 1}),
        ("H", {"H": 1}),
        ("HF", {"H": 1, "F": 1}),
        ("F", {"F": 1}),
        ("F2", {"F": 2}),
        ("SF6", {"S": 1, "F": 6}),
        ("SF5", {"S": 1, "F": 5}),
        ("SF4", {"S": 1, "F": 4}),
        ("SF3", {"S": 1, "F": 3}),
        ("SF2", {"S": 1, "F": 2}),
        ("SF", {"S": 1, "F": 1}),
        ("S", {"S": 1}),
        ("O2", {"O": 2}),
        ("O", {"O": 1}),
    )
    positive = (
        ("CF3+", {"C": 1, "F": 3}, 1),
        ("CHF2+", {"C": 1, "H": 1, "F": 2}, 1),
        ("CF2+", {"C": 1, "F": 2}, 1),
        ("CHF+", {"C": 1, "H": 1, "F": 1}, 1),
        ("CF+", {"C": 1, "F": 1}, 1),
        ("CH+", {"C": 1, "H": 1}, 1),
        ("F+", {"F": 1}, 1),
        ("F2+", {"F": 2}, 1),
        ("SF5+", {"S": 1, "F": 5}, 1),
        ("SF4+", {"S": 1, "F": 4}, 1),
        ("SF3+", {"S": 1, "F": 3}, 1),
        ("SF2+", {"S": 1, "F": 2}, 1),
        ("SF+", {"S": 1, "F": 1}, 1),
        ("S+", {"S": 1}, 1),
        ("SF4++", {"S": 1, "F": 4}, 2),
        ("SF2++", {"S": 1, "F": 2}, 2),
        ("O2+", {"O": 2}, 1),
    )
    negative = (
        ("F-", {"F": 1}),
        ("SF6-", {"S": 1, "F": 6}),
        ("SF5-", {"S": 1, "F": 5}),
        ("SF4-", {"S": 1, "F": 4}),
        ("SF3-", {"S": 1, "F": 3}),
        ("SF2-", {"S": 1, "F": 2}),
        ("F2-", {"F": 2}),
        ("O-", {"O": 1}),
    )
    source = "parent collision products; component primary sources in deck"
    return (
        Species(
            name="e",
            mass_amu=ELECTRON_MASS_AMU,
            charge_number=-1,
            composition={},
            role="electron",
            source="CODATA electron mass",
        ),
        *(Species(
            name=name,
            mass_amu=_mass(composition),
            charge_number=0,
            composition=composition,
            role="neutral",
            source=source,
            evidence_kind="published_compilation",
        ) for name, composition in neutral),
        *(Species(
            name=name,
            mass_amu=_mass(composition),
            charge_number=charge,
            composition=composition,
            role="positive_ion",
            source=source,
            evidence_kind="published_compilation",
        ) for name, composition, charge in positive),
        *(Species(
            name=name,
            mass_amu=_mass(composition),
            charge_number=-1,
            composition=composition,
            role="negative_ion",
            source=source,
            evidence_kind="published_compilation",
        ) for name, composition in negative),
    )


_CHF3_PRODUCTS = {
    "CF3 + H": {"CF3": 1, "H": 1},
    "CHF2 + F": {"CHF2": 1, "F": 1},
    "CF2 + H + F": {"CF2": 1, "H": 1, "F": 1},
    "CHF + F + F": {"CHF": 1, "F": 2},
    "CF + H + F + F": {"CF": 1, "H": 1, "F": 2},
    "CF + H + F2": {"CF": 1, "H": 1, "F2": 1},
    "CF3+ + H": {"CF3+": 1, "H": 1},
    "CHF2+ + F": {"CHF2+": 1, "F": 1},
    "CF2+ + HF": {"CF2+": 1, "HF": 1},
    "CHF+ + F + F": {"CHF+": 1, "F": 2},
    "CF+ + HF + F": {"CF+": 1, "HF": 1, "F": 1},
    "CH+ + F2 + F": {"CH+": 1, "F2": 1, "F": 1},
    "F+ + CHF2": {"F+": 1, "CHF2": 1},
    "CHF2 + F-": {"CHF2": 1, "F-": 1},
    "CHF2+ + F-": {"CHF2+": 1, "F-": 1},
}


_SF6_ION_PRODUCTS = {
    "SF5+": {"SF5+": 1, "F": 1},
    "SF4+": {"SF4+": 1, "F": 2},
    "SF3+": {"SF3+": 1, "F": 3},
    "SF2+": {"SF2+": 1, "F2": 1, "F": 2},
    "SF+": {"SF+": 1, "F2": 1, "F": 3},
    "S+": {"S+": 1, "F2": 3},
    "F+": {"F+": 1, "SF4": 1, "F": 1},
    "SF4++": {"SF4++": 1, "F": 2},
    "SF2++": {"SF2++": 1, "F2": 1, "F": 2},
}


_SF6_ATTACHMENT_PRODUCTS = {
    "SF6-": {"SF6-": 1},
    "SF5-": {"SF5-": 1, "F": 1},
    "SF4-": {"SF4-": 1, "F": 2},
    "SF3-": {"SF3-": 1, "F": 3},
    "SF2-": {"SF2-": 1, "F": 4},
    "F2-": {"F2-": 1, "SF4": 1},
    "F-": {"F-": 1, "SF5": 1},
}


def _heavy_products(target: str, kind: str, product: str | None):
    if target == "CHF3":
        if kind == "EXCITATION" and product and product.startswith("CHF3(v"):
            return {"CHF3": 1}, "measured_energy_loss_state_collapsed"
        try:
            return _CHF3_PRODUCTS[product], "swarm_regressed_working_set"
        except KeyError as exc:
            raise RuntimeError(f"unmapped CHF3 collision product {product!r}") from exc
    if target == "SF6":
        if product == "SF6(v) aggregate":
            return {"SF6": 1}, "nist_evaluated_vibrational_aggregate"
        if product == "SFx + F aggregate":
            return {"SF5": 1, "F": 1}, "dominant_branch_sensitivity"
        if kind == "IONIZATION":
            try:
                return _SF6_ION_PRODUCTS[product], "nist_anchor_threshold_closure"
            except KeyError as exc:
                raise RuntimeError(f"unmapped SF6 positive ion {product!r}") from exc
        if kind == "ATTACHMENT":
            try:
                return _SF6_ATTACHMENT_PRODUCTS[product], "nist_product_resolved"
            except KeyError as exc:
                raise RuntimeError(f"unmapped SF6 negative ion {product!r}") from exc
    if target == "O2":
        if product == "O + O aggregate":
            return {"O": 2}, "song_recommended_total_dissociation"
        if product == "O2+/O+ aggregate":
            return {"O2+": 1}, "unresolved_all_O2plus_sensitivity"
        if product == "O- + O":
            return {"O-": 1, "O": 1}, "song_recommended_attachment"
        if kind == "EXCITATION":
            return {"O2": 1}, "excited_state_collapsed_to_ground_inventory"
    raise RuntimeError(
        f"unmapped parent collision target={target!r} kind={kind!r} product={product!r}"
    )


def _build_mappings(deck: ElectronCollisionDeck):
    mappings = []
    for index, process in enumerate(deck.processes):
        if process.kind in {"MOMENTUM", "ELASTIC", "EFFECTIVE"}:
            continue
        products, evidence = _heavy_products(
            process.target, process.kind, process.product)
        mappings.append(ElectronCollisionHeavyMapping(
            process_index=index,
            reaction_name=(
                f"parent_electron_{index:02d}_{process.target}_"
                f"{process.kind.lower()}_{str(process.product).replace(' ', '_')}"
            ),
            heavy_reactants={process.target: 1},
            heavy_products=products,
            source=(
                "component collision source and explicit heavy mapping; "
                f"target={process.target}; product={process.product}"
            ),
            evidence_kind=evidence,
        ))
    return tuple(mappings)


def build_zhu_parent_collision_chemistry(
    song_o2_workbook: str | Path,
) -> ZhuParentCollisionChemistry:
    """Build the exact 55/5/1 parent-feed collision chemistry provider."""

    chf3 = derive_nist_evaluated_chf3_replay(
        load_kushner_zhang_2000_chf3_replay())
    sf6 = derive_nist_product_resolved_sf6_replay()
    o2 = load_song_2026_o2_replay(song_o2_workbook)
    mixed = compose_electron_collision_decks(
        (chf3.derived_deck, sf6.derived_deck, o2.derived_deck),
        retrieved_at="2026-08-18",
        mixture_name=(
            "Zhu NPG80 parent chemistry: 55 CHF3 / 5 SF6 / 1 O2 feed"
        ),
    )
    species = zhu_parent_collision_species()
    chemistry = ElectronCollisionChemistry(
        mixed, species, _build_mappings(mixed))
    return ZhuParentCollisionChemistry(
        chf3_replay=chf3,
        sf6_replay=sf6,
        o2_replay=o2,
        mixed_deck=mixed,
        species=species,
        collision_chemistry=chemistry,
    )
