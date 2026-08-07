# Lam reactor architectures: what the public record actually lets us build

## Verdict

Lam is highly relevant to the reactor program, but the useful target is not
“clone a Lam chamber from marketing pages.” The public record is sufficient
to define three equipment classes, the measurements each class must expose,
and the chemistry/depth boards that matter commercially. It is not sufficient
for a device-specific digital twin.

| public Lam class | public source architecture | market/application | first defensible petch target |
|---|---|---|---|
| Kiyo / Akara | transformer-coupled ICP lineage; Akara adds solid-state DirectDrive, fast plasma/species pulsing, and shaped ion-energy control | conductor etch for logic, DRAM, NAND | source/bias-separated ICP boundary; Ar then Cl2/Ar plasma state; time-resolved IEAD |
| Flex | small-volume, multi-frequency confined CCP | dielectric/HARC, fluorocarbon SiO2/SiN, ALE, cryogenic | dual/multi-frequency CCP boundary; C4F6/Ar/O2 species and ion-energy board |
| Vantex | public page identifies advanced-RF dielectric/HARC platform; internal source details remain undisclosed | deepest memory holes/trenches | HARC transport + cryogenic surface state after a measured equipment boundary is available |
| Syndion / DSiE | transformer-coupled deep-Si source with fast gas switching; Bosch and steady-state modes | TSV/HBM, MEMS, power trench, image sensor | pulsed SF6/C4F8 or SF6/O2 global balances coupled to time-resolved feature depth |
| Selis | high-pressure, radical-rich ICP | isotropic selective etch | residence-time, gas heating, wall recombination, and radical-flux validation |

## Why DirectDrive matters to the current failure

The independent Mahoney argon board failed because its `100 W` observable is
forward-minus-reflected generator power, while the model balance needs power
absorbed in the plasma. The source explicitly leaves matching-network and
coil loss unmeasured.

Hopwood measured 70–90% coupling in a different planar ICP by measuring the
matching/inductor dissipation. That proves the measurement topology, not a
portable efficiency constant.

Lam’s public DirectDrive patent family improves the upstream boundary:

1. voltage and current are sensed at the RF output;
2. phase is tracked for fast frequency control;
3. a reference network can connect the output to an RF power meter and dummy
   load for calibration; and
4. traditional slow motorized matching is removed from the direct-drive
   path.

Lam's newer `WO2024015694A1` record makes the measurement topology even more
concrete. It places V/I or phase-magnitude sensors at the matchless-source
output, at the output of a conventional matching network, or at upstream
nodes combined with a lumped matching-network model. It also describes
per-coil current sensing for multi-coil sources. Those signals can identify a
coil load and plasma transition, but load impedance still includes the
current-carrying source hardware. A state-matched resistive-loss measurement
or model is still required before `Re(VI*)` becomes plasma absorption.

The last physics step remains: calibrated source-output power minus
coil/window/electrode and other downstream hardware dissipation equals plasma
absorption. DirectDrive makes that step cleaner to instrument; it does not
delete it.

This distinction is now executable in `petch.reactor_global.power`. The code
will not promote either forward-minus-reflected power or DirectDrive output
power to predictive absorbed power without an independent loss closure.

## Market-relevant validation ladder

### 1. Universal kernel and RF boundary

- Conserve every element and charge in volume and wall reactions.
- Carry pumping, feed, wall return, and surface exchange explicitly.
- Accept direct calorimetry or a calibrated RF node plus independently
  measured hardware loss.
- Retain pulse waveforms/states rather than hiding them in a scalar average.

This layer is equipment-agnostic and is now under test in pure argon.

### 2. ICP conductor/deep-Si branch

- Verify Ar electron state and ion flux on independent planar ICP data.
- Add the Lee–Lieberman Cl2 network only after every rate and wall channel is
  provenance graded.
- Validate species-resolved boundary fluxes before any Lam-class profile
  claim.
- For Syndion/DSiE class work, couple fast gas switching to separate
  SF6-etch and C4F8-passivation surface states; cycle-averaged chemistry is
  insufficient for scallop and absolute-depth prediction.

The locally archived Zhang thesis has measured He/Cl2 profiles from a
commercial Lam ICP, but its reported plasma inputs are HPEM predictions and
its own feature comparison notes missing measured IEAD/flux ratios. It is a
useful later feature board, not independent reactor validation.

The highest-value public Lam-class plasma-state target is now the original
TCP paper by Ra, Bradley, and Chen, “Etching of aluminum alloys in the
transformer-coupled plasma etcher,” *JVST A* **12**, 1328–1333 (1994), DOI
`10.1116/1.579316`. Lee and Lieberman's Figure 8 identifies Ra's symbols as an
argon ion-density/current comparison, and the publisher abstract identifies
the Cl2/BCl3 TCP apparatus and source/bias separation. A second high-value
target is Patrick et al., *JVST A* **15**, 1250–1256 (1997), which reports
discharge-impedance and V/I-sensor characterization of aluminum etching.
Neither paper's quantitative values are admitted yet: only bibliographic and
abstract records have been recovered, so full primary text and original-pixel
figure audit are the retrieval gate.

### 3. CCP dielectric/HARC branch

- Build a two-frequency or multi-frequency sheath/source interface rather
  than reusing an electropositive ICP sheath.
- Validate C4F6/Ar mixture and pressure trends against Benck’s absolute,
  mass-resolved ion-current/IED board.
- Add O2 only with atom-balanced C/O/F gas and wall reactions.
- Couple the species/IED boundary to the already measured surface-yield and
  extreme-AR transport gates.

This is the direct route from the reactor program back to honest Krüger-class
depth. It will not reconstruct Krüger’s missing condition by naming the
reactor “Flex-like.”

### 4. Hardware partnership boundary

A device-specific prediction needs:

- chamber and plasma-volume geometry;
- gas-injection and pumping conductance;
- pressure and neutral-temperature measurements;
- source and bias voltage/current waveforms at declared nodes;
- plasma-off or calorimetric RF hardware-loss maps;
- wall/chuck/window temperatures and material state;
- species-resolved ion/radical flux or at least independent plasma-state
  diagnostics; and
- fixed-time blanket and feature outcomes held out from boundary
  calibration.

Lam’s 2025 announcement that UC Berkeley received a 2300 platform combining
Kiyo, Flex, and Syndion GP chambers is a concrete potential validation-access
route. No access or unpublished data is assumed.

## Commercial implication

The product need is not another profile fitter. It is the auditable bridge
from equipment telemetry and recipe states to species-resolved wafer
boundaries, then to measured surface physics and depth. Lam’s public
segmentation shows that the same kernel can serve several valuable markets,
but only through architecture-specific source, sheath, wall, and pulsing
providers:

- Flex/Vantex: dielectric HARC and cryogenic depth;
- Kiyo/Akara: conductor etch, atomic ion-energy/profile control;
- Syndion/DSiE: absolute deep-Si depth and cycle morphology; and
- Selis: radical-selective isotropic removal.

The installed bases make the prioritization economically concrete. Lam says
Akara inherits experience from more than **30,000 Kiyo chambers** in
production. Its Cryo 3.0 release says more than **7,500 Lam HARC dielectric
etch chambers** were then used in NAND production, nearly **1,000** with
cryogenic etch, and that cryogenic tools had processed five million wafers.
Those are manufacturer-reported installed-base facts, not independent market
share or performance measurements. They nevertheless show that an
architecture-aware boundary model can address both a very large conductor
fleet and a large, rapidly differentiating HARC fleet; the reusable value is
recipe transfer, diagnostics, and chamber-to-chamber calibration, not only
one simulated profile.

That is the reactor roadmap. “Lam digital twin” becomes an earned claim only
after one real chamber supplies the missing hardware boundary.

## Primary sources

- Lam, Akara:
  `https://www.lamresearch.com/product/akara/`
- Lam, Kiyo:
  `https://www.lamresearch.com/product/kiyo-product-family/`
- Lam, Flex:
  `https://www.lamresearch.com/product/flex-product-family/`
- Lam, Vantex:
  `https://www.lamresearch.com/product/vantex/`
- Lam, Syndion:
  `https://www.lamresearch.com/product/syndion-product-family/`
- Lam, DSiE:
  `https://www.lamresearch.com/product/dsie-product-family/`
- Lam, etch architecture overview:
  `https://newsroom.lamresearch.com/etch-essentials-semiconductor-manufacturing`
- Lam, 2016 Flex CCP release:
  `https://newsroom.lamresearch.com/2016-09-06-Lam-Research-Introduces-Dielectric-Atomic-Layer-Etching-Capability-for-Advanced-Logic`
- Lam, 2024 Cryo 3.0 release:
  `https://newsroom.lamresearch.com/2024-07-31-Lam-Research-Introduces-Lam-Cryo-TM-3-0-Cryogenic-Etch-Technology-to-Accelerate-Scaling-of-3D-NAND-for-the-AI-Era`
- Lam, 2025 Berkeley 2300 platform:
  `https://investor.lamresearch.com/2025-04-16-Lam-Research-Donates-Leading-Edge-Etch-System-to-Accelerate-Nanofabrication-R-D-at-UC-Berkeley`
- Lam patents: `US10515781B1`, `US12165841B2`, `US12106947B2`,
  `WO2024015694A1`.
- Hopwood 1994:
  `https://doi.org/10.1088/0963-0252/3/4/002`
- Mahoney et al. 1994:
  `https://doi.org/10.1063/1.357672`
- Ra, Bradley, and Chen 1994:
  `https://doi.org/10.1116/1.579316` (**abstract/bibliographic record only;
  full-text retrieval required before quantitative use**).
- Patrick et al. 1997:
  *JVST A* **15**, 1250–1256 (**bibliographic record only; full-text
  retrieval required before quantitative use**).
