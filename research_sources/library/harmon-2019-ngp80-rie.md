# harmon-2019-ngp80-rie

**Exact Oxford PlasmaPro NGP80 CHF3 process and within-run self-bias drift**

- **Citation:** Jeffrey Lee Harmon, *4H-SiC Trench-Gate MOSFET: Practical
  Surface-Channel Mobility Extraction*, PhD dissertation, North Carolina State
  University (2019).
- **Official record:**
  `https://repository.lib.ncsu.edu/items/73f6b7b4-34ad-4cfb-9bf3-a465df95eb77`
- **Official full text:**
  `https://repository.lib.ncsu.edu/server/api/core/bitstreams/619c77c4-da13-4040-9c28-92420dc4b447/content`
- **User-provided PDF SHA256:**
  `1f7c3530e062f045c809f1e407e14c4647f0d2a1671e7e868c33c9426747b634`
- **Status:** PRIMARY FULL THESIS; PDF 148 / PRINTED 126 VISUALLY AUDITED

The PDF metadata says author “Skyler Bunn,” but the title page and official
repository record identify Jeffrey Lee Harmon. The title page and repository
record control.

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| H1 | The NCSU-NNF tool is an Oxford PlasmaPro NGP80 RIE without an independent ICP source; its single RF generator produces both plasma and DC self-bias. | Exact equipment-class evidence that source density and bombardment cannot be varied independently in this topology. It does not make forward power a voltage measurement. |
| H2 | The etch step is `40 sccm CHF3`, `25 mTorr`, `200 W`, and `20 C`; the procedure says the DC-bias magnitude typically starts above `300 V` and ends below about `200 V`. | Two censored, time-resolved observations. They are preserved as `>300` and `<~200`, never collapsed into a measured `[200,300]` interval. |
| H3 | An empty-chamber oxygen clean is repeated until a pale-yellow plasma color is recovered, and post-clean duration scales with fluorocarbon etch time. | Direct evidence that chamber conditioning/wall inventory is a state variable for this exact model. |
| H4 | The developed condition reports `27.4 nm/min` SiO2 and `13.5 nm/min` photoresist removal, with repeatable `2:1` selectivity. | Adjacent CHF3 process validation only; neither rates nor selectivity transfer to TiO2/Cr. |

## Executable decision

`oxford80_self_bias.py` encodes H2 as separate censored start/end rows. The
threshold-to-threshold history is a deterministic sensitivity witness, not a
claim that either endpoint equals the true target-machine voltage.
