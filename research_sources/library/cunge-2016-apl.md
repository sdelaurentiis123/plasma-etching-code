# cunge-2016-apl

**Cunge et al., APL 108, 093109 (2016)**

- **DOI/URL:** APL 108, 093109 (2016)
- **Retrieval route:** publisher
- **Status:** ABSTRACT-ONLY
- **Topic:** iadf-sheath — Ion angular/energy distributions and sheath collisions

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | > "A series of capillary plates were placed on the RFEA surface, each with a specified AR > through which the IED was measured. The ion energy distributions were measured at the > wafer surface in an inductive coupled plasma (ICP) reactor from Applied Materials. A > helium plasma was used at a pressure of 10 mTorr. ICP power of 750 W and bias power set > for a DC self-bias voltage of 100 V." | `RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`:338 |
| Q2 | 1. **No in-feature measurement of anything at AR > ~50.** The best in-feature transport measurement is Cunge 2016's capillary-plate IED, run at modest AR in a He ICP. Ishikawa's 2018 review says it outright: "The problem of absolute radical fluxes inside HAR holes has not been completely solved yet" [Q-relay]. **Nobody has measured ion flux, radical flux, energy spectrum, or potential inside a 100:1 feature.** Every deep model in the field, including ours, is unconstrained there. *Fillable by us + a hardware partner:* a capillary/AR-plate transmission sweep to AR 100–200 is a one-fixture experiment that directly gates our exact cone-acceptance table with zero chemistry. | `RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`:671 |
| U3 | [unquoted — verify on next use] Impedans' Vertex generalisation replaces the physical capillary with an **electrostatically synthesised effective AR** (potential difference between grids 2 and 3), giving AR-resolved IEDs without venting [Q]. | `RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`:349 |
| U4 | [unquoted — verify on next use] **[INF] This is a direct experimental analogue of our cone-acceptance table.** Our study predicts, from an exact geometric cone plus the two-component beam, the fraction of beam transmitted at each AR (0.9669 → 0.5736 for core-only from AR 100 → 200; 0.4888 → 0.2421 at tail 0.65). Cunge's method measures exactly that transmitted fraction against AR, and inverts it for T_i. **A capillary-plate transmission-vs-AR curve is a preregistrable, zero- chemistry, zero-charging blind gate for our transport operator + beam model, and I found no evidence anyone has run it above AR ~50.** This is arguably the cheapest experiment a hardware partner could run that would validate the exact part of our engine that is strongest. | `RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`:353 |

_Harvested 2 quoted + 2 unquoted mentions across the repo's docs._
