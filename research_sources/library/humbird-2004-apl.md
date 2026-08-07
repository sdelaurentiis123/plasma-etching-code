# humbird-2004-apl

**Humbird, Graves, Hua & Oehrlein (2004), Ar+-driven F transport through FC films**

- **DOI:** https://doi.org/10.1063/1.1644338
- **Author full text:** https://www.researchgate.net/publication/224427642_Molecular_dynamics_simulations_of_Ar-induced_transport_of_fluorine_through_fluorocarbon_films
- **Status:** PRIMARY FULL TEXT ONLINE + VERIFIED EXCERPT
- **Extract:** `research_sources/thesis_extracts/humbird_2004_apl_primary_excerpt.txt`
- **Topic:** atomistic ion mixing, defluorination, Si-F bond formation

## Verified claims

| claim | use |
|---|---|
| thermal CF2 + normal Ar+, 9:1, compared at 20 and 200 eV | declares the MD boundary |
| Si-C-F depth increases from about 15 to 30 Å across 20-200 eV | energy-dependent mixed-layer depth |
| destroyed C-F bonds correspond nearly one-to-one with formed Si-F bonds | atom-balanced transfer channel |
| pure-Si physical sputter is only 0.11 Si/ion at 200 eV in the same potential | chemical enhancement is required |

These are atomistic model constraints, not reactor flux measurements.

## Related primary-author steady-state evidence

- **Primary-author seminar:** David Humbird and David Graves, UC Berkeley,
  FLCC Research Seminar, 23 February 2004,
  [Molecular Dynamics Simulations of Plasma-Surface Interactions and
  Etching](https://www.slideserve.com/moeshe/molecular-dynamics-simulations-of-plasma-surface-interactions-and-etching).
  Evidence grade: primary-author seminar, not peer reviewed.
- **Digitized board:**
  `data/experimental/humbird_graves_2004/seminar_surface_state_curves.csv`
- **Replay:** `scripts/digitize_humbird_graves_2004_seminar.py`
- **Verified source raster hashes:** slide 20
  `20ee4c5c2b9b0f9e522343073d10cbca800a33806e694c93a11607b2985ace1c`;
  slide 21
  `ff8d8610a8785c37046c86b69702314725be72a784635a5e57cdd75b1a425bd5`.

The seminar extends the CF2/Ar+ trajectories far beyond the short APL
comparison.  At 200 eV the C inventory continues to grow while the
instantaneous Si yield decays.  Adding thermal F sustains Si removal and
reduces the blocking inventory.  The authors' mechanism slide places a
Si--F reaction front beneath a finite-residence Si--C transport layer, with
fluorine recycling.  That topology directly rejects a same-step
``Si in == SiF4 out`` closure.

## Longer peer-reviewed mechanism papers

- Humbird and Graves, *Fluorocarbon plasma etching of silicon: Factors
  controlling etch rate*, **J. Appl. Phys. 96, 65 (2004)**,
  DOI [10.1063/1.1736321](https://doi.org/10.1063/1.1736321).  The abstract
  identifies the amorphous Si:C mixed layer, rather than only the overlying
  fluorocarbon film, as the primary rate limiter; small thermal-F addition
  increases its permeability.
- Humbird and Graves, *Mechanism of silicon etching in the presence of
  CF2, F, and Ar+*, **J. Appl. Phys. 96, 2466 (2004)**,
  DOI [10.1063/1.1769602](https://doi.org/10.1063/1.1769602).  The abstract
  reports a leading SiFx front followed by a Si--C layer and identifies ion
  penetration to that front as rate critical.

The seminar curves are used only as a surface-state validation board.  They
must not be relabeled as beam measurements, oxide data, quantum-accurate
dynamics, or a source of wafer flux.

## Implemented topology contract

`src/petch/stratified_fluorocarbon_si.py` implements the source-constrained
reservoir topology without embedding a fitted yield:

- separate FC film, Si--C transport layer, and Si--F reaction front;
- explicit C--F, C--C-crosslink, and Si--F bond inventories;
- paired C--F scission / F transfer / Si--F formation;
- nonzero Si residence (only start-of-step transported Si can leave);
- number-density-times-depth capacities rather than one-monolayer caps; and
- finite CSDA ion transmission through both live layers, including the
  incidence-angle slant path.

The module is deliberately an event-ledger core.  It does not infer event
probabilities from this paper.  Those belong to the separately preregistered
surface calibration in
`data/experimental/humbird_graves_2004/model_protocol.json`.
