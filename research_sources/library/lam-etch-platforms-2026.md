# lam-etch-platforms-2026

**Lam Research, public etch-platform architecture and application map**

- **Primary sources:** current official Lam product pages for Akara, Flex,
  Kiyo, Vantex, Syndion, and DSiE; Lam’s official “Etch Essentials” technical
  overview; official 2016 Flex CCP release; official 2009 annual report; and
  Lam’s 16 April 2025 UC Berkeley equipment announcement.
- **Retrieval date:** 2026-08-07
- **Status:** PRIMARY MANUFACTURER PRODUCT/TECHNICAL RECORDS READ
- **Topic:** public reactor class, market/application segmentation, and
  validation-access routes

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Lam describes conductor etch as commonly ICP/high-density and dielectric etch as commonly CCP/lower-density with high-energy ions. | The reactor provider must branch by source class; one 0-D closure cannot be relabeled across all Lam tools. |
| Q2 | Kiyo is a conductor-etch family with a transformer-coupled symmetric chamber, advanced pulsing, independent process tuning, and wafer-edge sheath control. | Prioritize separate source-power, bias/sheath, and spatial-uniformity interfaces for a Kiyo-class model. Public pages do not disclose chamber dimensions or recipes. |
| Q3 | Akara builds on Kiyo and introduces DirectDrive solid-state plasma sourcing, TEMPO species/power pulsing, and SNAP ion-energy control; Lam reports more than 30,000 Kiyo chambers in its installed-base lineage. | Highest-value conductor-etch integration target. Marketing performance statements are not validation measurements. |
| Q4 | Flex is dielectric etch; Lam’s 2016 release identifies its small-volume design as a CCP reactor. Current products add multi-frequency confined plasma, RF pulsing, ALE, and cryogenic operation for HARC holes/trenches. | Closest public Lam product class to fluorocarbon dielectric/HARC work, but not proof of Krüger apparatus equivalence. |
| Q5 | Vantex is a Sense.i dielectric/HARC platform with advanced RF and cryogenic capability. | Market-critical HARC target; the public product page does not publish enough electrical geometry to instantiate a device-specific digital twin. |
| Q6 | Syndion uses a transformer-coupled symmetric source and fast alternating etch/deposition for deep-Si TSV, HBM, trench, and power-device applications. DSiE supports Bosch and steady-state modes with fast gas switching. | Direct fit to a future SF6/C4F8 cyclic reactor+feature board and absolute-depth validation. |
| Q7 | Selis is described as a high-pressure ICP with high radical density for isotropic selective etch. | Requires radical residence-time/wall-recombination accuracy rather than a HARC ion-only emphasis. |
| Q8 | Sense.i evolved from Kiyo/Flex and emphasizes sensors, autonomous calibration, repeatability, and compact productivity. | Equipment telemetry is a boundary-data opportunity, not itself a plasma-physics validation claim. |
| Q9 | Lam’s 2025 Berkeley donation combines Kiyo conductor/metal etch, Flex dielectric etch, and Syndion GP DRIE chambers on a 2300 platform. | Concrete academic access path for cross-platform experiments, subject to Berkeley/Lam data permissions; no access is assumed. |
| Q10 | Lam's 2024 Cryo 3.0 release reports more than 7,500 Lam HARC dielectric chambers in NAND production, nearly 1,000 using cryogenic etch, and five million wafers processed with Lam cryogenic etch. | Manufacturer-reported installed-base evidence makes HARC/cryogenic transfer and fleet calibration a large product target. It is not independent market share or a profile-accuracy validation. |
| Q11 | The original Lam TCP validation leads are Ra, Bradley & Chen (1994), DOI `10.1116/1.579316`, and Patrick et al. (1997), *JVST A* 15, 1250–1256. | Their abstract/bibliographic records identify ion-current and V/I/impedance measurements in Lam-class ICP etching. Quantitative import remains blocked until full primary text and original pixels are recovered. |

## Use decision

The public record supports a **platform-class emulator**, not a claim that the
project has reconstructed proprietary Lam devices. The defensible branches
are:

1. transformer-coupled ICP conductor/deep-Si (`Kiyo`, `Akara`, `Syndion`);
2. small-volume multi-frequency CCP dielectric/HARC (`Flex`, likely the
   Vantex application class without asserting undisclosed internals); and
3. high-pressure radical-rich ICP selective etch (`Selis`).

A Lam-specific digital twin requires chamber dimensions, gas distribution,
RF measurement nodes and calibration, wall thermal state, pumping
conductance, and time-resolved recipe data from an actual tool.
