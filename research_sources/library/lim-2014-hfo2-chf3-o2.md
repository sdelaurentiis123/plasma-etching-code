# lim-2014-hfo2-chf3-o2

**CHF3/O2 daughter chemistry, diagnosed plasma, and HfO2 etch rates**

- **Citation:** N. Lim, A. Efremov, G. Y. Yeom, and K.-H. Kwon, "On the
  Etching Characteristics and Mechanisms of HfO2 Thin Films in CF4/O2/Ar and
  CHF3/O2/Ar Plasma for Nano-Devices," *Journal of Nanoscience and
  Nanotechnology* **14** (2014), 9670--9679.
- **DOI:** `10.1166/jnn.2014.10171`
- **Author/institution-hosted full text read online:**
  `https://spl.skku.ac.kr/_res/pnpl/etc/2014-14.pdf`
- **Status:** PRIMARY FULL TEXT; TABLE I TEXT AND PAGE IMAGE VERIFIED ONLINE;
  CONSERVATION-CHECKED EXECUTABLE SUBSET

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| L1 | Table I prints daughter-electron reactions R6, R10--R14, R17, R20--R21, and R23--R25 and neutral reactions R26--R75. | Landed as 62 conserved reactions in `lim_2014_chf3_oxygen_chemistry.py`; parent-feed and duplicate F2/H2/O2/O electron rows are excluded. |
| L2 | Neutral coefficients R26--R75 were compiled for an assumed gas temperature of 700 K; gas temperature was not measured. | Coefficients remain labeled `published_compilation`, and `supports_target_temperature_transfer` is false for the Zhu 350 K central state. |
| L3 | R26 gives `CHF3 + F -> HF + CF3` as `1.58e-13 cm3/s`, while Voloshin's 350 K mechanism gives `1.82e-12 cm3/s`. | Exposed as `lim_700K` and `voloshin_350K` branches; never averaged or selected by SEM/depth fit. |
| L4 | The surface scheme uses loss probabilities 0.05 for F, H, and CF3; 0.1 for CF2, CF, and O; and 1 for C and O(1d). | These are central development defaults for the radical-wall sensitivity, not measurements of the target Oxford chamber. |
| L5 | The plasma state was constrained using double-Langmuir-probe Te, ion current, floating potential, and total positive-ion density; Te and n+ entered the 0-D chemistry model. | The chemistry validation does not prove that forward generator power uniquely predicts an Oxford wafer flux. |
| L6 | HfO2 etch rate is non-monotonic with mixture; radical fluxes and polymer formation/removal matter in addition to F and ion energy flux. | A depth model must couple F, CFx, O, H, and ion energy/flux rather than collapse the reactor output to a single fluorine flux. |

## Executable decision

The conservation-checked Table-I subset closes the main C/H/F/O daughter
chain without duplicating the measured parent EEPF operator.  It is used as a
source branch and sensitivity layer.  It does not identify the Zhu machine's
absorbed power, surface conditioning, or gas temperature, and it cannot be
used as an etch-depth calibration because its target is HfO2 in a different
ICP reactor.
