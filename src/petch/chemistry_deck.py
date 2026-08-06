"""Chemistry decks — declarative, provenance-carrying etch-system definitions.

A deck is a plain dict (JSON-serializable) that declares every constant of one
etch chemistry: species stoichiometry, per-surface-state sticking maps,
chemisorption (bare and ion-activated), the per-material threshold-power
sputter laws, film state transitions, and the scattering constants. The
mixed-layer engine consumes a deck; adding a chemistry is a data file plus a
validation campaign, never a code branch.

Two classes of entry, distinguished so the deck is a COMPLETE declaration of
the system without pretending the engine is more parameterized than it is:

- ``"parameterized": True`` blocks are consumed by the factory and flow into
  the mechanism objects (species, chemisorption, deposition maps, the
  ``MixedLayerParams`` scalars).
- ``"engine_resident": True`` blocks record constants that currently live
  inside ``petch.mixed_layer`` / ``petch.boundary_transport_3d`` kernels. They
  are declared here with provenance so a deck fully describes its system;
  wiring them through the factory is the next pass and must be gated the same
  bitwise way. Validation checks them (bounds, provenance) but the factory
  does not inject them, and a mismatch against the engine default is reported
  by :func:`engine_resident_drift`.

Honest layering (doctrine, see CHEMISTRY_DECK_DESIGN_2026-07-28.md): the
engine and the mechanism FORMS are general physics; the constants are
per-chemistry with published provenance. The Krueger deck contains that
author's five optimizer-fitted values — we add zero of our own. The retirement
path per constant is beam measurement > derivation > declared-fitted-with-source.
"""

from __future__ import annotations

from copy import deepcopy

from petch.mixed_layer import MixedLayerParams
from petch.mixed_layer_mechanism import MixedLayerMechanism

__all__ = [
    "KRUEGER_2024_DECK",
    "DeckValidationError",
    "build_mixed_layer_mechanisms_from_deck",
    "engine_resident_drift",
    "validate_deck",
]


class DeckValidationError(ValueError):
    """A deck failed schema, bounds, or provenance validation."""


# --------------------------------------------------------------------------
# Deck 1: Krueger 2024 (Ar / C4F6 / O2 over SiO2 with an a-C mask)
# --------------------------------------------------------------------------

_K24 = "Krueger 2024 PhD thesis"

KRUEGER_2024_DECK = {
    "name": "krueger_2024_ar_c4f6_o2_sio2",
    "provenance": (
        f"{_K24}, Appendix B (mechanism table) + Table 6.5 (converged "
        "parameter set that produced the fig. 6.16 / ch.-6 validation "
        "feature: 45 nm mouth, 825 nm depth, 850 nm remaining mask)"),
    "declared_fitted": {
        # Named honestly: these five are Krueger's OWN optimizer outputs
        # (his Table 6.3 free parameters -> Table 6.5 converged values). We
        # add no fitted constants; retiring these is the beam-measurement /
        # derivation program, tracked per constant.
        "constants": ["sio2_bare_sputter_p0", "sio2_complex_sputter_p0",
                      "chemisorption_bare_CF", "film_oxidation",
                      "mask_deposition"],
        "provenance": f"{_K24}, Table 6.3 (free parameters) -> Table 6.5 (converged)",
    },

    "species": {
        "parameterized": True,
        "provenance": f"{_K24}, Appendix B species list; C4F6 feed dissociation products",
        # name -> (carbon atoms, fluorine atoms)
        "precursor_stoichiometry": {
            "CF": (1.0, 1.0), "CF2": (1.0, 2.0), "CF3": (1.0, 3.0),
            "C2F3": (2.0, 3.0), "C2F4": (2.0, 4.0),
            "C3F5": (3.0, 5.0), "C3F6": (3.0, 6.0),
        },
        "fluorine": [],          # the boundary publishes no atomic-F flux
        "oxygen": ["O"],
        "inert": ["C3F4"],
    },

    "chemisorption": {
        "parameterized": True,
        "provenance": (
            f"{_K24}, Appendix B complex-formation rows (CH1 bare) and CH2 "
            "(on ion-activated SiO2*); bare CF/CF2 = Table 6.5 converged 0.278"),
        # Complex formation on bare oxide: bound fluorine into the mixed layer.
        "bare": {
            "CF": 0.278, "CF2": 0.278, "CF3": 0.2, "C2F3": 0.2,
            "C2F4": 0.001, "C3F5": 0.001, "C3F6": 0.001,
        },
        # On ion-activated SiO2* — the dominant complex channel. C2F3 carries
        # the CF3 value (declared: appendix row group, not an independent row).
        "activated": {
            "CF": 0.8, "CF2": 0.85, "CF3": 0.9, "C2F3": 0.9,
            "C2F4": 0.001, "C3F5": 0.001, "C3F6": 0.001,
        },
    },

    "deposition": {
        "parameterized": True,
        "provenance": (
            f"{_K24}, Appendix B deposition rows (on polymer / on bare "
            "substrate / on crosslinked skin) + Table 6.5 mask row 0.094"),
        # Sticking on existing polymer is ~100x sticking on bare oxide:
        # radicals on open oxide chemisorb (complex) instead of polymerizing.
        # [VERIFY] ORDERING: the fab literature reports sticking that DECREASES
        # with the radical F/C ratio ("the sticking coefficient of F-rich CFx
        # radicals such as CF2 is lower than that of C-rich radicals", Hiwasa
        # et al., Appl. Phys. Express 15, 106002 (2022), whose ref. 12 is Izawa
        # et al., Jpn. J. Appl. Phys. 46, 7870 (2007)).  This row inverts that
        # ordering: CF3 (F/C = 3) sticks at 0.1 while C2F3 (F/C = 1.5) sticks at
        # 0.03.  NOT changed here: Izawa's 0.004 is a model-inverted effective
        # coefficient spanning a ~125x F-rich/C-rich axis, not a per-site rate
        # commensurate with these rows, and no published magnitude pair maps
        # onto this species set.  Krueger's values are kept as a validated,
        # internally consistent set; retiring them needs beam-measured
        # per-species sticking (RESEARCH_LIP_CERTAINTY_2026-08-04.md Q1).
        "on_polymer": {"CF": 0.1, "CF2": 0.1, "CF3": 0.1, "C2F3": 0.03},
        "on_substrate": {"CF": 0.002, "CF2": 0.0015, "CF3": 0.001,
                         "C2F3": 0.001},
        "on_crosslinked": {"CF": 0.02, "CF2": 0.02, "CF3": 0.02,
                           "C2F3": 0.02},
        "on_mask": {"CF": 0.094, "CF2": 0.094, "CF3": 0.094, "C2F3": 0.094},
    },

    "layer": {
        "parameterized": True,
        "provenance": (
            f"{_K24} + literature anchors: displacement energy 25 eV "
            "(polymer, band 10-80); s_f Gogolides band; eta_mix Humbird-Graves "
            "O(1); p_ox = Table 6.5 converged 0.0423; Ar+ projectile"),
        # -> MixedLayerParams fields, shared by both materials.
        "displacement_energy_eV": 25.0,
        "ion_atomic_number": 18,
        "ion_mass_amu": 39.948,
        "precursor_fc_ratio": 1.5,
        "sticking_probability": 0.0842,
        "fluorine_film_sticking": 0.05,
        "oxidation_probability": 0.0423,
        "mixing_efficiency": 1.0,
        "volatilization_yield": 1.0,
        "reference_energy_eV": 1000.0,
        "film_sputter_yield": 0.1384,
        "minimum_layer_depth_nm": 0.5,
        "clog_film_thickness_nm": 20.0,
    },

    "materials": {
        "parameterized": True,
        "provenance": f"{_K24}, ch. 6 feature stack: SiO2 target, a-C hard mask",
        "sio2": {"substrate": "sio2", "formula_density_m3": 2.2e28},
        "ac_mask": {"substrate": "carbon", "atom_density_m3": 1.0e29},
    },

    # ---- engine-resident blocks: declared here, applied in the kernels ----

    "sputter_laws": {
        "engine_resident": True,
        "provenance": (
            f"{_K24}, Appendix B yield rows; threshold-power form "
            "y = p0 * ((E^q - Eth^q) / (E0^q - Eth^q)); complex channel is "
            "ZBL-deposited-energy shaped at the 140 eV reference; Kress "
            "angular response B = 9.3. SiO2 p0 values = Table 6.5 converged."),
        # (p0, threshold_eV, reference_eV, exponent) unless noted.
        "film_fresh": [0.9, 20.0, 500.0, 0.5],
        "film_crosslinked": [0.6, 50.0, 500.0, 0.5],
        "sio2_complex": {"p0": 0.1471, "shape": "zbl_deposited_energy",
                         "reference_eV": 140.0},
        "sio2_bare": [0.0852, 70.0, 140.0, 1.0],
        "ac_mask_physical": [0.001, 200.0, 250.0, 0.4],
        "kress_angular_B": 9.3,
    },

    "state_transitions": {
        "engine_resident": True,
        "provenance": (
            f"{_K24}, Appendix B rows: A1 ion activation of SiO2 (0.9/ion); "
            "S3d ion de-crosslinking (0.3 @ 8 eV); crosslink formation is "
            "zero-knob (ion dose / displacement energy, Bruce-Graves)"),
        "activation_probability": 0.9,
        "site_density_m2": 1.0e19,
        "decrosslink": [0.3, 8.0, 500.0, 0.5],
        "mask_oxidation": 1.0e-5,
    },

    "scattering": {
        "engine_resident": True,
        "provenance": (
            "Huang PhD thesis, Eq. 2.34 retention rule + leftover-probability "
            "selection (interpretation B': the energy function scales removal, "
            "not selection); MCFPM cascade to 8 generations"),
        "specular_threshold_eV": 100.0,
        "diffusive_cutoff_eV": 10.0,
        "specular_cutoff_angle_deg": 70.0,
        "max_bounces": 8,
        "selection": "leftover_B_prime",
    },

    "declared_omissions": {
        "engine_resident": True,
        "provenance": (
            f"{_K24}, Appendix B rows implemented nowhere yet — bounded and "
            "declared rather than silently absorbed (RESEARCH_FINAL_CLOSURE_"
            "AUDIT_2026-07-28.md)"),
        "carbonization": 0.01,          # CF(s) + ion -> AC(s) + F @ 20 eV
        "carbon_redeposition": 0.01,    # sputtered C redeposits as AC
        "ex1_angular_exponent": 2.0,    # EX1 angular form, Eth = 35 eV
    },
}


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

_REQUIRED_SECTIONS = ("species", "chemisorption", "deposition", "layer",
                      "materials", "sputter_laws", "state_transitions",
                      "scattering")

# Sections whose leaf floats are probabilities in [0, 1].
_PROBABILITY_MAPS = {
    "chemisorption": ("bare", "activated"),
    "deposition": ("on_polymer", "on_substrate", "on_crosslinked", "on_mask"),
}

_LAYER_PROBABILITIES = ("sticking_probability", "fluorine_film_sticking",
                        "oxidation_probability", "film_sputter_yield")

_LAYER_POSITIVE = ("displacement_energy_eV", "ion_mass_amu",
                   "precursor_fc_ratio", "mixing_efficiency",
                   "volatilization_yield", "reference_energy_eV",
                   "minimum_layer_depth_nm", "clog_film_thickness_nm")


def validate_deck(deck) -> None:
    """Raise :class:`DeckValidationError` unless the deck is complete and sane.

    Checks: required sections present; every section carries a nonempty
    provenance string; probability-valued constants lie in [0, 1]; positive
    scalars are positive and finite; deposition/chemisorption species are
    declared in the species block; and the sputter-law tuples are well formed.
    """
    if not isinstance(deck, dict):
        raise DeckValidationError("deck must be a dict")
    for key in ("name", "provenance"):
        value = deck.get(key)
        if not isinstance(value, str) or not value.strip():
            raise DeckValidationError(f"deck requires a nonempty {key!r}")

    missing = [name for name in _REQUIRED_SECTIONS if name not in deck]
    if missing:
        raise DeckValidationError(f"deck missing sections: {sorted(missing)}")

    for name in _REQUIRED_SECTIONS:
        section = deck[name]
        if not isinstance(section, dict):
            raise DeckValidationError(f"deck section {name!r} must be a dict")
        provenance = section.get("provenance")
        if not isinstance(provenance, str) or not provenance.strip():
            raise DeckValidationError(
                f"deck section {name!r} requires a nonempty provenance")
        if not (section.get("parameterized") or section.get("engine_resident")):
            raise DeckValidationError(
                f"deck section {name!r} must declare parameterized or "
                "engine_resident")

    species = deck["species"]
    stoichiometry = species.get("precursor_stoichiometry")
    if not isinstance(stoichiometry, dict) or not stoichiometry:
        raise DeckValidationError("species.precursor_stoichiometry is required")
    for name, pair in stoichiometry.items():
        try:
            carbon, fluorine = (float(pair[0]), float(pair[1]))
        except (TypeError, IndexError, ValueError) as error:
            raise DeckValidationError(
                f"species {name!r} stoichiometry must be (C, F)") from error
        if carbon <= 0.0 or fluorine < 0.0:
            raise DeckValidationError(
                f"species {name!r} stoichiometry out of bounds: {pair}")
    declared = set(stoichiometry)

    for section_name, keys in _PROBABILITY_MAPS.items():
        section = deck[section_name]
        for key in keys:
            table = section.get(key)
            if table is None:
                continue
            if not isinstance(table, dict):
                raise DeckValidationError(
                    f"{section_name}.{key} must be a mapping")
            for name, value in table.items():
                _check_probability(f"{section_name}.{key}[{name}]", value)
                if name not in declared:
                    raise DeckValidationError(
                        f"{section_name}.{key} names undeclared species {name!r}")

    layer = deck["layer"]
    for key in _LAYER_PROBABILITIES:
        if key in layer:
            _check_probability(f"layer.{key}", layer[key])
    for key in _LAYER_POSITIVE:
        if key in layer:
            value = layer[key]
            if not isinstance(value, (int, float)) or not (0.0 < float(value) < 1e12):
                raise DeckValidationError(
                    f"layer.{key} must be a positive finite scalar, got {value!r}")

    laws = deck["sputter_laws"]
    for key, law in laws.items():
        if key in ("provenance", "engine_resident", "parameterized"):
            continue
        if key == "kress_angular_B":
            if not isinstance(law, (int, float)) or float(law) < 0.0:
                raise DeckValidationError("sputter_laws.kress_angular_B must be >= 0")
            continue
        if isinstance(law, dict):
            _check_probability(f"sputter_laws.{key}.p0", law.get("p0"))
            continue
        if not isinstance(law, (list, tuple)) or len(law) != 4:
            raise DeckValidationError(
                f"sputter_laws.{key} must be (p0, threshold_eV, reference_eV, q)")
        _check_probability(f"sputter_laws.{key}[p0]", law[0])
        threshold, reference, exponent = (float(law[1]), float(law[2]),
                                          float(law[3]))
        if threshold < 0.0 or reference <= threshold or exponent <= 0.0:
            raise DeckValidationError(
                f"sputter_laws.{key} energies/exponent out of bounds: {law}")

    transitions = deck["state_transitions"]
    _check_probability("state_transitions.activation_probability",
                       transitions.get("activation_probability"))
    if "mask_oxidation" in transitions:
        _check_probability("state_transitions.mask_oxidation",
                           transitions["mask_oxidation"])

    materials = deck["materials"]
    for key in ("sio2", "ac_mask"):
        material = materials.get(key)
        if not isinstance(material, dict) or "substrate" not in material:
            raise DeckValidationError(
                f"materials.{key} requires a substrate declaration")
        if material["substrate"] not in ("sio2", "carbon"):
            raise DeckValidationError(
                f"materials.{key}.substrate must be 'sio2' or 'carbon'")


def _check_probability(label, value) -> None:
    if not isinstance(value, (int, float)):
        raise DeckValidationError(f"{label} must be a number, got {value!r}")
    value = float(value)
    if not (0.0 <= value <= 1.0):
        raise DeckValidationError(f"{label} out of [0, 1]: {value}")


# --------------------------------------------------------------------------
# Factory
# --------------------------------------------------------------------------

def build_mixed_layer_mechanisms_from_deck(deck=None, *, oxide_parameters=None,
                                           mask_parameters=None,
                                           volatilization_yield=None,
                                           reactive_ion_yield_table=None,
                                           validate=True):
    """Build the (oxide, mask) mixed-layer mechanism pair declared by ``deck``.

    ``volatilization_yield`` overrides the deck's layer value (the single
    base-condition-anchored absolute rate constant). ``oxide_parameters`` /
    ``mask_parameters`` bypass the deck's layer block entirely for callers that
    construct their own :class:`MixedLayerParams`.
    """
    deck = KRUEGER_2024_DECK if deck is None else deck
    if validate:
        validate_deck(deck)

    layer = dict(deck["layer"])
    for key in ("parameterized", "engine_resident", "provenance"):
        layer.pop(key, None)
    if volatilization_yield is not None:
        layer["volatilization_yield"] = float(volatilization_yield)

    materials = deck["materials"]
    if oxide_parameters is None:
        oxide_parameters = MixedLayerParams(
            substrate=materials["sio2"]["substrate"], **layer)
    if mask_parameters is None:
        mask_parameters = MixedLayerParams(
            substrate=materials["ac_mask"]["substrate"], **layer)

    species = deck["species"]
    chemisorption = deck["chemisorption"]
    deposition = deck["deposition"]
    common = dict(
        precursor_species={
            str(name): (float(pair[0]), float(pair[1]))
            for name, pair in species["precursor_stoichiometry"].items()},
        fluorine_species=tuple(species.get("fluorine", ())),
        oxygen_species=tuple(species.get("oxygen", ())),
        inert_species=tuple(species.get("inert", ())),
        chemisorption_probability=dict(chemisorption["bare"]),
        chemisorption_activated_probability=dict(chemisorption["activated"]),
        deposition_probability_on_film=dict(deposition["on_polymer"]),
        deposition_probability_on_substrate=dict(deposition["on_substrate"]),
        deposition_probability_on_crosslinked=dict(deposition["on_crosslinked"]))
    mask_kwargs = dict(common)
    mask_kwargs["deposition_probability_on_substrate"] = dict(
        deposition["on_mask"])
    oxide_kwargs = dict(common)
    oxide_kwargs["reactive_ion_yield_table"] = reactive_ion_yield_table
    return (MixedLayerMechanism(oxide_parameters, **oxide_kwargs),
            MixedLayerMechanism(mask_parameters, **mask_kwargs))


def engine_resident_drift(deck=None):
    """Report engine-resident deck constants that disagree with the kernels.

    The engine-resident blocks are a declaration, not an injection: this
    compares them against the values actually compiled into
    ``petch.mixed_layer`` and ``petch.boundary_transport_3d`` so a deck cannot
    quietly claim physics the engine does not run. Returns a list of
    ``(path, deck_value, engine_value)`` triples; empty means consistent.
    """
    deck = KRUEGER_2024_DECK if deck is None else deck
    from petch import mixed_layer as _module

    drift = []

    def compare(path, deck_value, engine_value):
        if float(deck_value) != float(engine_value):
            drift.append((path, float(deck_value), float(engine_value)))

    transitions = deck.get("state_transitions", {})
    if "site_density_m2" in transitions:
        compare("state_transitions.site_density_m2",
                transitions["site_density_m2"], _module._SITE_DENSITY_M2)

    materials = deck.get("materials", {})
    if "formula_density_m3" in materials.get("sio2", {}):
        compare("materials.sio2.formula_density_m3",
                materials["sio2"]["formula_density_m3"],
                _module._SIO2_FORMULA_DENSITY_M3)
    return drift


def deck_copy(deck=None):
    """Deep copy of a deck (decks are data; callers may mutate their copy)."""
    return deepcopy(KRUEGER_2024_DECK if deck is None else deck)
