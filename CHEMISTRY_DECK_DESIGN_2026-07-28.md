# Chemistry Deck — the generalizable etch-system format (design v1)

Goal: a new etch chemistry is a DATA FILE with provenance, never code. The
mixed-layer engine consumes a deck; all system-specific constants live in it.

## Deck schema (Python dict / JSON; every constant carries provenance)

```
deck = {
  "name": "krueger_2024_ar_c4f6_o2_sio2",
  "provenance": "Krueger PhD thesis 2024, Appendix B + Table 6.5 (converged)",
  "species": {                      # gas species -> (C_atoms, F_atoms, role)
    "CF": (1, 1, "precursor"), ..., "O": (0, 0, "oxidant"),
    "C3F4": (3, 4, "inert"), "ions": ("aggregate", "energetic")},
  "scattering": {"E_ts_eV": 100.0, "E_c_eV": 10.0, "theta_c_deg": 70.0,
                  "selection": "leftover_B_prime"},
  "materials": {
    "sio2": {
      "densities": {"formula_m3": 2.2e28},
      "products": {"complex": {"formula": "SiF4", "f_cost": 4},
                    "bare_sputter": {"f_cost": 0}},
      "sputter": {"bare": [0.0852, 70, 1.0, 140, "kress"],
                   "complex": [0.1471, 35, 1.0, 140, "chang_sawin"]},
      "activation": {"probability": 0.9,
                      "chemisorption_activated": {"CF": 0.8, "CF2": 0.85,
                                                    "CF3": 0.9, "C2F3": 0.9}},
      "chemisorption_bare": {"CF": 0.278, "CF2": 0.278, "CF3": 0.2, ...},
      "deposition_on_substrate": {"CF": 0.002, ...}},
    "ac_mask": {
      "densities": {"atom_m3": 1.0e29},
      "sputter": {"physical": [0.001, 200, 0.4, 250, "kress"]},
      "deposition_on_substrate": {"CF": 0.094, ...},
      "oxidation": 1.0e-5}},
  "film": {
    "density_m3": 7.5e28,
    "deposition_on_fresh": {"CF": 0.1, "CF2": 0.1, "CF3": 0.1, "C2F3": 0.03},
    "deposition_on_crosslinked": {"CF": 0.02, ...},
    "sputter_fresh": [0.9, 20, 0.5, 500, "kress"],
    "sputter_crosslinked": [0.6, 50, 0.5, 500, "kress"],
    "decrosslink": [0.3, 8, 0.5, 500],
    "oxidation": 0.0423,
    "displacement_energy_eV": 25.0},
}
```

## Factory
`build_mechanisms_from_deck(deck) -> {material_id: MixedLayerMechanism}` —
replaces every hardcoded KRUEGER_2024_* constant; the router takes the deck.
Validation at load: schema, bounds, provenance-nonempty, F-ledger
consistency of product stoichiometry.

## Decks roadmap
1. krueger_2024 (extract the constants now in code — verbatim move, gated
   bitwise vs current behavior).
2. belen_sf6_o2_si (the existing validated Si arm — proves generality:
   two chemistries, one engine, zero code branches).
3. sawin-derived decks (beam-measured constants; blind campaigns).
4. Partner/fab decks (their measured constants — the commercial surface).

## Honest layering (doctrine)
Engine + mechanism FORMS are general physics. Constants are per-chemistry
with published provenance; the Krueger deck contains his five
optimizer-fitted values (we add zero). Retirement path per constant:
beam measurement > DEKNOB-style derivation > declared-fitted-with-source.
