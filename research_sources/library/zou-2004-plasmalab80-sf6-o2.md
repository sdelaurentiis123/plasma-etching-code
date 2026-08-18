# zou-2004-plasmalab80-sf6-o2

**Oxford PlasmaLab 80+ SF6/O2 pressure, self-bias, loading, and profile study**

- **Citation:** H. Zou, “Anisotropic Si deep beam etching with profile control
  using SF6/O2 plasma,” *Microsystem Technologies* **10**, 603--607 (2004).
- **DOI:** `10.1007/s00542-003-0338-3`
- **User-provided PDF SHA256:**
  `b23135afe677a4816ff460e20949c8173a9b71f123cbe8702f3414ad570f8421`
- **Status:** PRIMARY FULL ARTICLE; FIGURES 1--4 AND TABLE 1 VISUALLY AUDITED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Z1 | Experiments use an Oxford PlasmaLab 80+ parallel-plate RIE at `13.56 MHz`, with a `20 cm` water-cooled base electrode held at `20 C`. | Same equipment family and electrode-temperature boundary as the Zhu condition; not the same serial tool. |
| Z2 | Figure 2 uses `12 sccm SF6`, `160 W`, and an O2 sweep; DC bias spans `360--387 V` at `30 mTorr` and `117--175 V` at `200 mTorr`. | Direct 30 mTorr voltage range for machine-family sensitivity. The oxygen sweep means it is not one fixed-composition point. |
| Z3 | Samples are about `3.5 cm2` on a dummy four-inch wafer for consistent loading, with a `200 nm` Cr mask. | Loading and powered-area context that blocks careless transfer to an unknown chip/carrier geometry. |
| Z4 | The paper states that DC bias rises with dark-space voltage, which rises with increasing RF input power and decreasing pressure, with gas-composition dependence. | Qualitative reduced-drive prior only; no universal `P/p` law is installed. |

## Executable decision

The measured `360--387 V` range supplies two fixed sensitivity witnesses. The
code preserves the range, chemistry sweep, and loading difference rather than
fitting a voltage-power curve through it.
