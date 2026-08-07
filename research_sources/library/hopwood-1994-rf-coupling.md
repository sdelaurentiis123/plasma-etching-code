# hopwood-1994-rf-coupling

**Hopwood, measured planar-ICP RF coupling efficiency**

- **Citation:** J. Hopwood, “Planar RF induction plasma coupling
  efficiency,” *Plasma Sources Science and Technology* **3**, 460–464
  (1994).
- **DOI:** `10.1088/0963-0252/3/4/002`
- **Primary author/institution record:**
  `https://research.ibm.com/publications/planar-rf-induction-plasma-coupling-efficiency`
- **Status:** PRIMARY AUTHOR ABSTRACT READ; FULL TEXT NOT LOCALLY ARCHIVED
- **Topic:** measured matching-network/inductor loss and the distinction
  between RF delivery and plasma absorption

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Coupling efficiency was determined by measuring dissipation in the matching network and inductive coupler. | Direct precedent for a separate hardware-loss measurement in the power boundary. |
| Q2 | The abstract reports 20–60 A RMS RF current and an equivalent coupling-network resistance of 0.09 ohm. | Shows why even small RF resistance creates material `I²R` loss. These numbers are apparatus-specific. |
| Q3 | Reported plasma-coupling efficiency was 70–90%, with the remainder assigned to circuit ohmic heating. | Contextual range only; it cannot be transplanted into Mahoney, Kiyo, Akara, or Krüger. |
| Q4 | Efficiency was nearly constant from 200 to 2000 W in that apparatus, but lower near 1 mTorr and without magnetic confinement. | Coupling is state and hardware dependent even when a broad power plateau exists. |

## Use decision

This primary abstract validates the topology of the new power contract:
measure RF power at a declared node and subtract independently measured
hardware dissipation. It does not land a universal 70–90% efficiency
constant. Full text and apparatus matching are required before any numerical
loss law is used in a predictive deck.
