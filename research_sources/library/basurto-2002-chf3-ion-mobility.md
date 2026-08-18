# basurto-2002-chf3-ion-mobility

**Mass-resolved molecular-ion mobility in CHF3**

- **Citation:** E. Basurto and J. de Urquijo, "Mobility of CF3+ in CF4,
  CHF2+ in CHF3, and C+ in Ar," *Journal of Applied Physics* **91**, 36--39
  (2002).
- **DOI:** `10.1063/1.1421034`
- **Full text retrieved from author-hosted PDF:**
  `https://tpc_gas.cfnssbu.physics.sunysb.edu/tpc_gas/Papers/cf3_mobility.pdf`
- **PDF SHA256:**
  `a959267e28687a1497ddb96669bdc096d096446a3c08765e1a5a5c3ecdd53b48`
- **Local text extraction:**
  `research_sources/thesis_extracts/basurto-2002-chf3-ion-mobility.txt`
- **Status:** PRIMARY FULL TEXT + FIGURE 1 CHECKSUM-PINNED PIL DIGITIZATION
- **Digitization receipt:**
  `data/experimental/basurto_2002_chf3/digitization_manifest.json`
- **Replay script:**
  `scripts/digitize_basurto_2002_chf2_chf3_mobility.py`
- **Digitized CSV SHA256:**
  `96dc38fc862888ae2252be123a2209166ac3a6d35c2e2c0718ab33aead227775`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| B1 | The experiment measures mass-resolved CHF2+ mobility in CHF3 over about `30--750 Td`, at `5--100 mTorr` and `293--310 K`, with final mobility uncertainties of `2--4%`. | Direct swarm-scale evidence for low-energy molecular-ion transport in the target parent gas. It is not a CF3+ measurement. |
| B2 | The reported low-field reduced mobility for CHF2+ in CHF3 is `0.53 +/- 0.01 cm2 V-1 s-1`, versus a `1.33 cm2 V-1 s-1` polarization-limit estimate. | Rejects a naive Langevin/polarization closure for the molecular ion; clustering or strong valence interactions are plausible source explanations. |
| B3 | CHF2+ is the dominant ion created in the CHF3 source, followed by smaller CF+ and CF3+ populations. | Species-topology evidence for a future CHF3 reactor deck; it does not establish the Zhu tool's wafer fractions. |
| B4 | The source does not invert its mobility curve into a unique elastic differential cross section. | Figure 1 can grade a candidate molecular collision provider, but cannot by itself close the feature-facing angular kernel. |

## Executable decision

The open-circle Figure 1 curve is installed as a no-extrapolation, C1 PCHIP
mobility closure in `petch.reactor_global.chf3_ion_mobility`.  It can close
measured CHF2+ drift and effective momentum-relaxation scales over the
digitized `45.18--450.35 Td` support.  It cannot be inverted into a unique
elastic differential cross section and does not close a target sheath IEAD.
Its low-field discrepancy forbids silently using the Langevin value.
