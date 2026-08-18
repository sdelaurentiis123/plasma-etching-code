# janissen-2016-tio2-rie

**Janissen et al., tunable CHF3/O2 RIE of single-crystal TiO2 with Cr masks**

- **Citation:** R. Janissen et al., “Tunable top-down fabrication and
  functional surface coating of single-crystal titanium dioxide
  nanostructures and nanoparticles,” *Nanoscale* **8**, 14010–14017 (2016).
- **DOI:** `10.1039/C6NR00898D`
- **Primary full text:**
  `https://pubs.rsc.org/en/content/articlehtml/2016/nr/c6nr00898d`
- **Author-hosted PDF:**
  `https://www.bionanotech.eu/publications/Janissen_Nanoscale_2016.pdf`
- **Supplement reproduced in primary-author thesis:** S. Ha, *Single-Crystal
  Titanium Dioxide Nanostructures*, TU Delft (2018), Chapter 3.
- **Thesis full text:**
  `https://pure.tudelft.nl/ws/portalfiles/portal/47058891/PhD_Thesis_S.Ha.pdf`
- **Thesis PDF SHA256:**
  `d4e8afedd1349a91b14ee10f589e9cae01aab48eb31d2be44433cd26e6fce912`
- **Status:** PRIMARY FULL TEXT + SUPPLEMENTARY TABLES S3.1--S3.3
  CHECKSUM-PINNED AND VISUALLY AUDITED
- **Extraction receipt:**
  `data/experimental/janissen_2016_tio2/extraction_manifest.json`
- **Topic:** TiO2/Cr feature etching, fluorocarbon passivation, selectivity,
  and machine-specific oxygen sensitivity

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The reported optimized CHF3/O2 RIE process removes single-crystal TiO2 at approximately 40 nm/min with approximately 14:1 TiO2:Cr selectivity and produces near-vertical nanocylinders around 89 degrees. | Adjacent material/mask comparison for the Zhu NPG80 condition. It is not an absolute prediction for a different tool or for a CHF3/SF6/O2 feed. |
| Q2 | CHF3 supplies both TiO2-etching species and fluorocarbon passivation species; small O2 changes tune the balance between lateral passivation and removal. | Supports treating O2 as a surface-state/profile control rather than equating its feed fraction with a wafer reaction flux. |
| Q3 | In the reported machine, low O2 flows from zero to about 1 sccm span positive, near-vertical, and negative sidewall responses, while substantially higher O2 can produce hourglass profiles. | The supplied 1 sccm condition makes a near-vertical regime plausible, but the added SF6 and machine transfer prevent a quantitative profile claim. |
| Q4 | A second nominally similar machine required a substantially different O2-flow range to obtain comparable profiles. | Direct evidence that a recipe flow cannot be transferred across chambers without machine-specific boundary calibration. |
| Q5 | The within-batch `Ra` power sweep at `50 ubar`, `50/5/30 sccm CHF3/O2/Ar` reports TiO2 rates of `30`, `58`, and `68 nm/min` at `100`, `165`, and `200 W`; Cr rates are `1.7`, `3.2`, and `3.7 nm/min`, leaving the printed selectivity near `18:1`. | A quantitative, out-of-target process-response board. It can grade a reactor-plus-surface model without using the Zhu SEM. It does not isolate ion flux, energy, or radical flux. |
| Q6 | The closest printed stack witness (Fig. 3.2a) used a `45 nm` Cr mask of `175 nm` diameter, `50/0.5 sccm CHF3/O2`, `200 W`, `50 ubar` (`37.5 mTorr`), measured `-950 V` DC bias, and `11 min`; the main text reports a `430 nm`-tall, `110 nm`-diameter, `89 degree` cylinder and about `40 nm/min`, `14:1`. | Directly demonstrates that comparable RIE forward-power/pressure can coexist with a voltage far above the Oxford-family witnesses. The machine is a Fluor Z401S and the substrate is single-crystal rutile, so the voltage and yield are not target coefficients. |
| Q7 | On a second nominally identical Fluor Z401S at `-1100 V`, `200 W`, `50 ubar`, and `50/4 sccm CHF3/O2`, the source reports `652 nm` height after `15 min` and `273 nm` after `8 min`, with global height RSD `1.4%` and `3.1%`. | Two scale-bearing feature-depth outcomes for model validation. The non-proportional implied rates expose pattern/loading/startup dependence and forbid treating one profile height divided by time as a universal surface yield. |

## Use in petch

The exact Tables S3.1--S3.3 are now committed as replayable CSVs after a
300-dpi original-resolution visual audit. They form a process-response and
feature-depth validation board. No row is installed directly as an Oxford
NPG80 reactor input or a TiO2 surface coefficient. At the comparison
selectivity, a 45 nm Cr mask protects only 630 nm of TiO2 removal, exposing the
mask-survival gate before the Zhu SEM is revealed. Any surface-law inference
must preserve the source's single-crystal-rutile, machine, mask, and loading
boundaries; the Zhu film is ALD TiO2 with phase and density still unreported.
