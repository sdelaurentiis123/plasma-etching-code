"""Rights-safe daughter-electron collision sensitivities for Zhu NPG80.

The Oxford reactor state is dominated by HF after feed dissociation, but the
published parent-only calculation treats every daughter as transparent to the
two-term electron operator.  This module implements two literature-declared
closures without packaging any LXCat source bytes:

* Huang et al. (JVST A 38, 023007, 2020) use the HCl momentum-transfer curve
  for HF and threshold-shift the corresponding HCl dissociation and ionization
  curves.  HF vibration and attachment come from different sources and are
  deliberately *not* invented here.
* The SIGLO F2 set tabulates an EFFECTIVE momentum curve.  LXCat defines that
  curve as elastic momentum transfer plus the total inelastic cross sections.
  The petch operator assembles all inelastic momentum explicitly, so the
  elastic residual must be recovered before composition or those channels are
  counted twice.

Both transformations are deterministic sensitivities, not new measurements,
swarm validation, a closed nonlinear reactor state, or a depth prediction.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import json
import math

import numpy as np

from .argon import ELECTRON_MASS_AMU
from .electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)
from .electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionHeavyMapping,
)
from .electron_collision_mixture import compose_electron_collision_decks
from .network import Species
from .zhu_parent_collision_chemistry import ZhuParentCollisionChemistry


HF_MASS_AMU = 20.006243163
HF_DISSOCIATION_THRESHOLD_EV = 5.87
HF_IONIZATION_THRESHOLD_EV = 16.007
HUANG_2020_DOI = "10.1116/1.5125568"
LXCAT_EFFECTIVE_DEFINITION_DOI = "10.1002/ppap.201600098"


@dataclass(frozen=True)
class HFDaughterCollisionSensitivity:
    derived_deck: ElectronCollisionDeck
    source_hcl_payload_sha256: str
    omitted_hf_channels: tuple[str, ...] = (
        "vibrational excitation: requires Rohr and Linder HF data",
        "dissociative attachment: requires Xu/Gallup/Fabrikant HF data",
        "rotational excitation and superelastic channels",
    )
    evidence_class: str = "huang_2020_hcl_surrogate_partial_hf"
    supports_complete_hf_eedf: bool = False
    supports_unique_reactor_state: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            self.derived_deck.targets != ("HF",)
            or len(self.derived_deck.processes) != 3
            or not self.omitted_hf_channels
            or self.evidence_class != "huang_2020_hcl_surrogate_partial_hf"
            or self.supports_complete_hf_eedf
            or self.supports_unique_reactor_state
            or self.supports_feature_depth
        ):
            raise ValueError("invalid partial HF collision sensitivity")


@dataclass(frozen=True)
class F2EffectiveDeconvolution:
    derived_deck: ElectronCollisionDeck
    source_f2_payload_sha256: str
    maximum_energy_eV: float
    minimum_elastic_cross_section_m2: float
    evidence_class: str = "siglo_effective_momentum_deconvolution"
    supports_direct_f2_swarm_grade: bool = False
    supports_unique_reactor_state: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            self.derived_deck.targets != ("F2",)
            or self.maximum_energy_eV <= 0.0
            or self.minimum_elastic_cross_section_m2 <= 0.0
            or self.evidence_class != "siglo_effective_momentum_deconvolution"
            or self.supports_direct_f2_swarm_grade
            or self.supports_unique_reactor_state
            or self.supports_feature_depth
        ):
            raise ValueError("invalid F2 effective-momentum deconvolution")


@dataclass(frozen=True)
class ZhuAugmentedCollisionChemistry:
    """One fully mapped parent+HF+F2 collision provider for reactor closure."""

    parent_chemistry: ZhuParentCollisionChemistry
    hf_replay: HFDaughterCollisionSensitivity
    f2_replay: F2EffectiveDeconvolution
    mixed_deck: ElectronCollisionDeck
    species: tuple[Species, ...]
    collision_chemistry: ElectronCollisionChemistry
    supplemental_reactions_replaced: tuple[str, ...]
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
            or not self.supplemental_reactions_replaced
            or not self.supports_parent_collision_sources
            or self.supports_complete_daughter_eedf
            or self.supports_unique_reactor_state
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid augmented Zhu collision chemistry")


def _derivation_hash(payload: dict[str, object]) -> str:
    return sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def _one_process(
    deck: ElectronCollisionDeck,
    *,
    kind: str,
    energy_loss_eV: float | None = None,
) -> ElectronCollisionProcess:
    candidates = [
        process for process in deck.processes
        if process.kind == kind
        and (
            energy_loss_eV is None
            or math.isclose(
                float(process.energy_loss_eV or 0.0),
                energy_loss_eV,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        )
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"expected one {kind} process at loss {energy_loss_eV!r}"
        )
    return candidates[0]


def _threshold_shift(
    process: ElectronCollisionProcess,
    *,
    threshold_eV: float,
    product: str,
) -> ElectronCollisionProcess:
    source_threshold = float(process.energy_loss_eV or 0.0)
    shift = float(threshold_eV) - source_threshold
    energies = tuple(value + shift for value in process.electron_energy_eV)
    if energies[0] < 0.0:
        raise ValueError("threshold shift produced negative electron energy")
    return replace(
        process,
        target="HF",
        product=product,
        electron_energy_eV=energies,
        energy_loss_eV=float(threshold_eV),
        comments=process.comments + (
            f"Huang 2020 HCl-to-HF threshold shift: {shift:.12g} eV",
        ),
    )


def derive_huang_2020_partial_hf_replay(
    hcl_deck: ElectronCollisionDeck,
) -> HFDaughterCollisionSensitivity:
    """Apply only the HCl-to-HF transfers declared by Huang et al.

    The Hayashi HCl 5.29 eV aggregate excitation is the corresponding
    dissociative-loss curve used for the threshold transfer.  The separate HCl
    vibration and attachment rows are intentionally excluded because Huang's
    HF mechanism cites different HF-specific sources for those processes.
    """

    if hcl_deck.targets != ("HCl",):
        raise ValueError("HF replay requires one target-pure HCl deck")
    elastic = _one_process(hcl_deck, kind="ELASTIC")
    dissociation = _one_process(
        hcl_deck, kind="EXCITATION", energy_loss_eV=5.29
    )
    ionization = _one_process(hcl_deck, kind="IONIZATION")
    hf_elastic = replace(
        elastic,
        target="HF",
        product=None,
        mass_ratio=ELECTRON_MASS_AMU / HF_MASS_AMU,
        comments=elastic.comments + (
            "Huang 2020: HF momentum transfer assumed equal to HCl; "
            "mass ratio changed to HF",
        ),
    )
    hf_dissociation = _threshold_shift(
        dissociation,
        threshold_eV=HF_DISSOCIATION_THRESHOLD_EV,
        product="H + F",
    )
    hf_ionization = _threshold_shift(
        ionization,
        threshold_eV=HF_IONIZATION_THRESHOLD_EV,
        product="HF+",
    )
    derivation = {
        "schema": "petch.huang-2020-partial-hf-collision-replay.v1",
        "source_hcl_payload_sha256": hcl_deck.payload_sha256,
        "source_hcl_database": hcl_deck.source_database,
        "momentum_rule": "same cross section as HCl; HF mass ratio",
        "dissociation_rule": {
            "source_hcl_threshold_eV": dissociation.energy_loss_eV,
            "hf_threshold_eV": HF_DISSOCIATION_THRESHOLD_EV,
        },
        "ionization_rule": {
            "source_hcl_threshold_eV": ionization.energy_loss_eV,
            "hf_threshold_eV": HF_IONIZATION_THRESHOLD_EV,
        },
        "hf_vibration_included": False,
        "hf_attachment_included": False,
    }
    deck = ElectronCollisionDeck(
        processes=(hf_elastic, hf_dissociation, hf_ionization),
        payload_sha256=_derivation_hash(derivation),
        source_database="Huang-2020 partial HF replay from user HCl deck",
        retrieved_at=hcl_deck.retrieved_at,
        source_reference=(
            f"Huang et al., JVST A 38 023007 ({HUANG_2020_DOI}); "
            f"source deck: {hcl_deck.source_reference}; derivation="
            + json.dumps(derivation, sort_keys=True)
        ),
    )
    return HFDaughterCollisionSensitivity(
        derived_deck=deck,
        source_hcl_payload_sha256=hcl_deck.payload_sha256,
    )


def _interpolate_process(
    process: ElectronCollisionProcess,
    energies_eV: np.ndarray,
) -> np.ndarray:
    source_energy = np.asarray(process.electron_energy_eV, dtype=float)
    source_sigma = np.asarray(process.cross_section_m2, dtype=float)
    return np.interp(
        energies_eV,
        source_energy,
        source_sigma,
        left=0.0,
        right=float(source_sigma[-1]),
    )


def _extend_zero_tail(
    process: ElectronCollisionProcess,
    *,
    maximum_energy_eV: float,
) -> ElectronCollisionProcess:
    final_energy = process.electron_energy_eV[-1]
    if final_energy >= maximum_energy_eV:
        return process
    if process.cross_section_m2[-1] != 0.0:
        raise ValueError(
            "F2 process ends below the working domain with a nonzero tail"
        )
    return replace(
        process,
        electron_energy_eV=(
            *process.electron_energy_eV, float(maximum_energy_eV)
        ),
        cross_section_m2=(*process.cross_section_m2, 0.0),
        comments=process.comments + (
            f"zero tail extended to {maximum_energy_eV:g} eV",
        ),
    )


def deconvolve_siglo_f2_effective_momentum(
    f2_deck: ElectronCollisionDeck,
    *,
    maximum_energy_eV: float = 120.0,
) -> F2EffectiveDeconvolution:
    """Recover the F2 elastic residual on one declared working domain."""

    maximum = float(maximum_energy_eV)
    if f2_deck.targets != ("F2",) or not math.isfinite(maximum) or maximum <= 0:
        raise ValueError("F2 replay requires one target-pure deck and domain")
    effective = _one_process(f2_deck, kind="EFFECTIVE")
    if maximum > effective.electron_energy_eV[-1]:
        raise ValueError("F2 working domain exceeds effective-curve support")
    inelastic = tuple(
        _extend_zero_tail(process, maximum_energy_eV=maximum)
        for process in f2_deck.processes
        if process is not effective
    )
    energy = np.asarray([
        value for value in effective.electron_energy_eV if value < maximum
    ] + [maximum], dtype=float)
    effective_sigma = _interpolate_process(effective, energy)
    inelastic_sum = np.sum([
        _interpolate_process(process, energy) for process in inelastic
    ], axis=0)
    elastic_sigma = effective_sigma - inelastic_sum
    if np.any(~np.isfinite(elastic_sigma)) or np.any(elastic_sigma <= 0.0):
        index = int(np.argmin(elastic_sigma))
        raise ValueError(
            "F2 effective deconvolution became nonpositive at "
            f"{energy[index]:g} eV"
        )
    elastic = ElectronCollisionProcess(
        kind="ELASTIC",
        target="F2",
        product=None,
        electron_energy_eV=tuple(energy.tolist()),
        cross_section_m2=tuple(elastic_sigma.tolist()),
        mass_ratio=float(effective.mass_ratio),
        comments=effective.comments + (
            "elastic residual = SIGLO effective minus every tabulated "
            "inelastic cross section",
            f"deconvolution working domain: 0--{maximum:g} eV",
        ),
    )
    derivation = {
        "schema": "petch.siglo-f2-effective-deconvolution.v1",
        "source_f2_payload_sha256": f2_deck.payload_sha256,
        "maximum_energy_eV": maximum,
        "inelastic_process_count": len(inelastic),
        "minimum_elastic_cross_section_m2": float(np.min(elastic_sigma)),
        "effective_definition_doi": LXCAT_EFFECTIVE_DEFINITION_DOI,
        "zero_tail_extension_only_when_source_endpoint_is_zero": True,
    }
    deck = ElectronCollisionDeck(
        processes=(elastic, *inelastic),
        payload_sha256=_derivation_hash(derivation),
        source_database="SIGLO F2 effective-momentum deconvolution",
        retrieved_at=f2_deck.retrieved_at,
        source_reference=(
            f"LXCat effective definition {LXCAT_EFFECTIVE_DEFINITION_DOI}; "
            f"source deck: {f2_deck.source_reference}; derivation="
            + json.dumps(derivation, sort_keys=True)
        ),
    )
    return F2EffectiveDeconvolution(
        derived_deck=deck,
        source_f2_payload_sha256=f2_deck.payload_sha256,
        maximum_energy_eV=maximum,
        minimum_elastic_cross_section_m2=float(np.min(elastic_sigma)),
    )


def zhu_hf_f2_replaced_supplemental_reactions(
    *,
    kokkoris_eedf_shape: str = "druyvesteyn",
) -> tuple[str, ...]:
    """Scalar rows superseded by energy-resolved HF/F2 collision moments."""

    if kokkoris_eedf_shape not in {"druyvesteyn", "maxwellian"}:
        raise ValueError("unsupported Kokkoris EEDF shape")
    return (
        f"kokkoris_2009_G6_{kokkoris_eedf_shape}",
        f"kokkoris_2009_G7_{kokkoris_eedf_shape}",
        f"kokkoris_2009_G16_{kokkoris_eedf_shape}",
        f"kokkoris_2009_G19_{kokkoris_eedf_shape}",
        "lim_2014_R14",
    )


def _daughter_heavy_products(
    process: ElectronCollisionProcess,
) -> tuple[dict[str, int], str]:
    if process.target == "HF":
        if process.kind == "EXCITATION" and process.product == "H + F":
            return {"H": 1, "F": 1}, "huang_hcl_threshold_shift"
        if process.kind == "IONIZATION" and process.product == "HF+":
            return {"HF+": 1}, "huang_hcl_threshold_shift"
    if process.target == "F2":
        if process.kind == "ATTACHMENT":
            return {"F-": 1, "F": 1}, "siglo_product_resolved"
        if process.kind == "IONIZATION":
            return {"F2+": 1}, "siglo_product_resolved"
        if process.kind == "EXCITATION":
            if process.energy_loss_eV in {3.16, 4.34}:
                return {"F": 2}, "kokkoris_dissociation_interpretation"
            return {"F2": 1}, "internal_state_collapsed_to_inventory"
    raise RuntimeError(
        "unmapped daughter collision "
        f"{process.target}/{process.kind}/{process.product}"
    )


def build_zhu_augmented_collision_chemistry(
    parent: ZhuParentCollisionChemistry,
    hf_replay: HFDaughterCollisionSensitivity,
    f2_replay: F2EffectiveDeconvolution,
    *,
    reactor_species: tuple[Species, ...],
    kokkoris_eedf_shape: str = "druyvesteyn",
) -> ZhuAugmentedCollisionChemistry:
    """Compose and atom/charge-map the enlarged reactor collision basis."""

    if not isinstance(parent, ZhuParentCollisionChemistry):
        raise TypeError("a Zhu parent collision provider is required")
    if not isinstance(hf_replay, HFDaughterCollisionSensitivity):
        raise TypeError("an HF daughter replay is required")
    if not isinstance(f2_replay, F2EffectiveDeconvolution):
        raise TypeError("an F2 daughter replay is required")
    species = tuple(reactor_species)
    mixed = compose_electron_collision_decks(
        (parent.mixed_deck, hf_replay.derived_deck, f2_replay.derived_deck),
        retrieved_at="2026-08-18",
        mixture_name="Zhu NPG80 parent plus partial HF and SIGLO F2",
    )
    mappings = [
        ElectronCollisionHeavyMapping(
            process_index=mapping.process_index,
            reaction_name=mapping.reaction_name,
            heavy_reactants=mapping.heavy_reactants,
            heavy_products=mapping.heavy_products,
            source=mapping.source,
            evidence_kind=mapping.evidence_kind,
        )
        for mapping in parent.collision_chemistry.mappings
    ]
    parent_count = len(parent.mixed_deck.processes)
    for index, process in enumerate(
        mixed.processes[parent_count:], start=parent_count
    ):
        if process.kind in {"MOMENTUM", "ELASTIC", "EFFECTIVE"}:
            continue
        products, evidence = _daughter_heavy_products(process)
        mappings.append(ElectronCollisionHeavyMapping(
            process_index=index,
            reaction_name=(
                f"daughter_electron_{index:02d}_{process.target}_"
                f"{process.kind.lower()}_"
                f"{str(process.product).replace(' ', '_')}"
            ),
            heavy_reactants={process.target: 1},
            heavy_products=products,
            source=(
                "derived daughter collision source and explicit heavy "
                f"mapping; target={process.target}; product={process.product}"
            ),
            evidence_kind=evidence,
        ))
    chemistry = ElectronCollisionChemistry(
        mixed, species, tuple(mappings)
    )
    return ZhuAugmentedCollisionChemistry(
        parent_chemistry=parent,
        hf_replay=hf_replay,
        f2_replay=f2_replay,
        mixed_deck=mixed,
        species=species,
        collision_chemistry=chemistry,
        supplemental_reactions_replaced=(
            zhu_hf_f2_replaced_supplemental_reactions(
                kokkoris_eedf_shape=kokkoris_eedf_shape
            )
        ),
    )
